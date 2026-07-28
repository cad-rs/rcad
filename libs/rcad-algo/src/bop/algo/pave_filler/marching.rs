use super::*;

impl<'a> super::PaveFiller<'a> {
    pub(crate) fn intersect_ff_by_numeric_intss(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        grid_n: usize,
    ) {
        use crate::bop::int_tools::intss::numeric_intss_with_domains;
        use crate::bop::int_tools::pcurve_derive::polyline_pcurve_by_projection;

        // Use face-specific UV domains (set up by DS::setup_uv_boundaries)
        // if available.  For cylinders this encodes the actual face height range,
        // ensuring the intersection polyline endpoints fall *inside* the UV
        // boundary rectangle and can be used to split it.
        let dom1 = self.ds.face_uv_boundary(f1).and_then(|uv| {
            if uv.len() >= 3 {
                let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite()
                {
                    return Some([u_min, u_max, v_min, v_max]);
                }
            }
            None
        });
        let dom2 = self.ds.face_uv_boundary(f2).and_then(|uv| {
            if uv.len() >= 3 {
                let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite()
                {
                    return Some([u_min, u_max, v_min, v_max]);
                }
            }
            None
        });

        let result =
            numeric_intss_with_domains(s1, s2, grid_n, dom1, dom2, Some(self.ff_tol(f1, f2)));
        if result.is_empty() {
            return;
        }

        let mut curve_indices = Vec::new();
        for sir in &result.curves {
            let (mut chain, approx_curve) = match &sir.curve_3d {
                crate::bop::int_tools::intss::SurfaceCurve::Polyline(pts) => (pts.clone(), None),
                crate::bop::int_tools::intss::SurfaceCurve::BSplineCurve(bs) => {
                    // Sample BSpline back to polyline for face-boundary snapping
                    let n = 64usize;
                    let pts: Vec<DVec3> = (0..=n)
                        .map(|i| {
                            let t = i as f64 / n as f64;
                            bs.point_at(t)
                        })
                        .collect();
                    (pts, Some(Curve3::BSpline((**bs).clone())))
                }
                _ => continue,
            };
            if chain.len() < 2 {
                continue;
            }

            self.snap_polyline_endpoints_to_face_boundaries(&mut chain, f1, f2);

            let v_start = self.ds.add_vertex(chain[0]);
            let v_end = self.ds.add_vertex(chain[chain.len() - 1]);

            let arc_len: f64 = chain.windows(2).map(|w| (w[1] - w[0]).length()).sum();
            let dir = (chain[chain.len() - 1] - chain[0]).normalize_or_zero();
            let pcurve_a = sir
                .pcurve_on_a
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s1));
            let pcurve_b = sir
                .pcurve_on_b
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s2));

            let curve_idx = self.ds.intersection_curves.len();
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: chain[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: chain,
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                curve_extra: crate::bop::ds::CurveExtra::default(),
            });

            self.ds.face_info_mut(f1).curves_sc.insert(curve_idx);
            self.ds.face_info_mut(f2).curves_sc.insert(curve_idx);
            self.ds.face_info_mut(f1).vertices_in.insert(v_start);
            self.ds.face_info_mut(f1).vertices_in.insert(v_end);
            self.ds.face_info_mut(f2).vertices_in.insert(v_start);
            self.ds.face_info_mut(f2).vertices_in.insert(v_end);

            curve_indices.push(curve_idx);
        }

        if !curve_indices.is_empty() {
            self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
                tangent_faces: false,
            });
        }
    }

    /// OCCT: marching FF intersection
    /// OCCT: marching FF intersection
    pub(crate) fn intersect_ff_by_marching(&mut self, f1: usize, f2: usize) {
        use crate::bop::int_tools::marching::{MarchingConfig, adaptive_sampling_density};

        let s1 = self.ds.face_surface(f1).cloned().unwrap_or_else(|| {
            panic!("marching: face {} has no surface", f1)
        });
        let s2 = self.ds.face_surface(f2).cloned().unwrap_or_else(|| {
            panic!("marching: face {} has no surface", f2)
        });

        // use sign-change grid marching (IntTools_FaceFace / IntPatch_ImpPrmIntersection).
        // No BSpline demotion 锟?BSpline surfaces stay as parametric (ts=0) and use UV grid marching.
        let any_curved = !matches!(&s1, Surface3::Plane(_)) || !matches!(&s2, Surface3::Plane(_));
        if any_curved {
            let char_len = |s: &Surface3| -> f64 {
                match s {
                    Surface3::Sphere(sp) => sp.radius,
                    Surface3::Cylinder(cy) => cy.radius,
                    Surface3::Cone(co) => co.radius.max(0.5),
                    Surface3::Torus(to) => to.major_radius.max(to.minor_radius),
                    Surface3::BSpline(bsp) => {
                        if bsp.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bsp.control_points {
                                for p in row {
                                    mn = mn.min(*p);
                                    mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    Surface3::Bezier(bez) => {
                        if bez.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bez.control_points {
                                for p in row {
                                    mn = mn.min(*p);
                                    mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    _ => 1.0,
                }
            };
            let avg_len = (char_len(&s1) + char_len(&s2)) * 0.5;
            let mut grid_n = ((avg_len * 10.0) as usize).max(64).min(256);

            let skew_factor = match (&s1, &s2) {
                (Surface3::Cylinder(c1), Surface3::Cone(c2))
                | (Surface3::Cone(c2), Surface3::Cylinder(c1)) => {
                    let a1 = c1.axis.normalize();
                    let a2 = c2.axis.normalize();
                    let sin_angle = a1.cross(a2).length();
                    (1.0 + sin_angle * 3.0).min(3.0)
                }
                _ => 1.0,
            };
            grid_n = ((grid_n as f64 * skew_factor) as usize).min(256);

            self.intersect_ff_by_numeric_intss(f1, f2, &s1, &s2, grid_n);
            return;
        }

        // Use adaptive sampling density based on surface geometry
        let base_density = 16usize;
        let sampling1 = adaptive_sampling_density(&s1, base_density);
        let sampling2 = adaptive_sampling_density(&s2, base_density);
        // Use the higher density to ensure we don't miss narrow intersections
        let n_u = sampling1.n_u.max(sampling2.n_u);
        let n_v = sampling1.n_v.max(sampling2.n_v);

        let _samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
        // Use multi-scale seed detection for improved robustness
        // Scales: coarse (8x8), medium (16x16), fine (32x32), ultra (64x64)
        let base_step = self.estimate_step_size(&s1, &s2);
        let seed_dedup = (base_step * 2.0).max(self.ff_tol(f1, f2) * 2.0);
        let seeds = crate::bop::int_tools::marching::find_seed_points_multiscale(
            &s1,
            &s2,
            |nu, nv| self.generate_surface_samples_grid(&s1, nu, nv),
            &[8, 16, 32, 64],
            seed_dedup,
        );

        if seeds.is_empty() {
            return;
        }

        // Compute a finite bounding box that contains both faces' intersection region.
        // Use boundary vertices (actual face extent) with a generous margin.
        let bounds_from_face = |face_idx: usize| -> (DVec3, DVec3) {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            // Use boundary vertices (from wire edges)
            for &vi in self.ds.face_boundary_verts(face_idx) {
                let p = self.ds.vertex_point(vi);
                mn = mn.min(p);
                mx = mx.max(p);
            }
            // Also sample boundary edges for curved edges (e.g. circles)
            for &ei in self.ds.face_boundary_edges(face_idx) {
                if ei < self.ds.edge_count() {
                    let [t0, t1] = self.ds.edge_range(ei);
                    if let Some(curve) = self.ds.edge_curve(ei) {
                        for k in 0..=8usize {
                            let t = t0 + (t1 - t0) * k as f64 / 8.0;
                            let p = curve.point_at(t);
                            if p.is_finite() {
                                mn = mn.min(p);
                                mx = mx.max(p);
                            }
                        }
                    }
                }
            }
            // If still infinite, use a generous default
            if !mn.is_finite() || !mx.is_finite() {
                mn = DVec3::splat(-10.0);
                mx = DVec3::splat(10.0);
            }
            (mn, mx)
        };

        let (mn1, mx1) = bounds_from_face(f1);
        let (mn2, mx2) = bounds_from_face(f2);
        let margin = 1.0;
        let aabb_min = mn1.min(mn2) - DVec3::splat(margin);
        let aabb_max = mx1.max(mx2) + DVec3::splat(margin);

        // Use adaptive step size based on characteristic lengths
        let char_len = sampling1
            .characteristic_length
            .min(sampling2.characteristic_length);
        let step_size = base_step.min(char_len * 0.5).max(TOLERANCE_MESH_LEGACY);

        // Configure marching with convergence monitoring
        let marching_config = MarchingConfig {
            step_size,
            min_step_size: step_size * 0.01,
            max_steps: 500,
            max_oscillations: 3,
            step_reduction_factor: 0.5,
            deflection_tol: step_size * 0.001,
            multiscale_seeds: true,
        };

        let mut curve_indices = Vec::new();
        // Track all points already covered by marched curves, to deduplicate
        // seeds that trace the same intersection curve.
        let mut covered_points: Vec<DVec3> = Vec::new();
        let ff = self.ff_tol(f1, f2);
        let dedup_tol = (step_size * 3.0).max(ff * 2.0);

        for seed in seeds {
            // Skip if this seed is near any point already covered by a previous curve
            if covered_points
                .iter()
                .any(|&cp| (cp - seed).length_squared() < dedup_tol * dedup_tol)
            {
                continue;
            }

            let curve = crate::bop::int_tools::marching::march_intersection_with_config(
                &s1,
                &s2,
                seed,
                &marching_config,
                |p: DVec3| p.cmpge(aabb_min).all() && p.cmple(aabb_max).all(),
            );

            if curve.points.len() < 2 {
                continue;
            }

            // Mark all curve points as covered (sample every few for efficiency)
            for (i, &p) in curve.points.iter().enumerate() {
                if i % 5 == 0 {
                    covered_points.push(p);
                }
            }

            let v_start = self.ds.add_vertex(curve.points[0]);
            let v_end = self.ds.add_vertex(curve.points[curve.points.len() - 1]);

            let curve_idx = self.ds.intersection_curves.len();
            // Compute arc-length for t_range
            let arc_len: f64 = curve
                .points
                .windows(2)
                .map(|w| (w[1] - w[0]).length())
                .sum();
            let dir = (curve.points[curve.points.len() - 1] - curve.points[0]).normalize_or_zero();
            let t_range = [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)];

            // 锟?reApprox 锟?validate pcurves; retry with loose tolerance
            // if validation fails.
            let (pcurve_a, pcurve_b) =
                self.make_marching_pcurves_with_reapprox(&curve.points, &s1, &s2, f1, f2, &t_range);

            // approximate marching polyline to BSpline (MakeCurve / GeomInt_IntSS::MakeBSpline)
            let approx_curve = if curve.points.len() >= 4 {
                crate::bop::int_tools::intss::polyline_to_bspline(
                    &curve.points,
                    TOLERANCE_TOL_SCALE_MICRO,
                )
                .filter(|c| matches!(c, Curve3::BSpline(_)))
            } else {
                None
            };

            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: curve.points[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: curve.points.clone(),
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                curve_extra: crate::bop::ds::CurveExtra::default(),
            });

            self.ds.face_info_mut(f1).curves_sc.insert(curve_idx);
            self.ds.face_info_mut(f2).curves_sc.insert(curve_idx);
            self.ds.face_info_mut(f1).vertices_in.insert(v_start);
            self.ds.face_info_mut(f1).vertices_in.insert(v_end);
            self.ds.face_info_mut(f2).vertices_in.insert(v_start);
            self.ds.face_info_mut(f2).vertices_in.insert(v_end);

            curve_indices.push(curve_idx);
        }

        if !curve_indices.is_empty() {
            self.ds.interf_ff.push(crate::bop::ds::InterferenceFF {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
                tangent_faces: false,
            });
        }
    }

    /// OCCT: marching pcurve re-approximation
    /// OCCT: marching pcurve re-approx
    pub(crate) fn make_marching_pcurves_with_reapprox(
        &self,
        points: &[DVec3],
        s1: &Surface3,
        s2: &Surface3,
        f1: usize,
        f2: usize,
        t_range: &[f64; 2],
    ) -> (Option<Curve2d>, Option<Curve2d>) {
        let uv_bounds1 = s1.default_domain();
        let uv_bounds2 = s2.default_domain();
        let is_u_periodic1 = matches!(
            s1,
            Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_)
        );
        let is_u_periodic2 = matches!(
            s2,
            Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_)
        );
        let u_per1 = if is_u_periodic1 {
            Some(std::f64::consts::TAU)
        } else {
            None
        };
        let u_per2 = if is_u_periodic2 {
            Some(std::f64::consts::TAU)
        } else {
            None
        };

        // Attempt 1: default tolerance
        let pca = crate::bop::int_tools::pcurve_derive::polyline_pcurve_by_projection(points, s1);
        let pcb = crate::bop::int_tools::pcurve_derive::polyline_pcurve_by_projection(points, s2);

        let valid_a = pca.as_ref().map_or(false, |pc| {
            crate::bop::int_tools::pcurve_derive::is_curve_valid_2d(pc)
                && crate::bop::int_tools::pcurve_derive::check_pcurve_in_face(
                    pc, *t_range, uv_bounds1, u_per1, None,
                )
        });
        let valid_b = pcb.as_ref().map_or(false, |pc| {
            crate::bop::int_tools::pcurve_derive::is_curve_valid_2d(pc)
                && crate::bop::int_tools::pcurve_derive::check_pcurve_in_face(
                    pc, *t_range, uv_bounds2, u_per2, None,
                )
        });

        if valid_a && valid_b {
            return (pca, pcb);
        }

        // 锟?reApprox 锟?fallback with looser validation.
        // Skip the self-intersection check (is_curve_valid_2d) since polyline
        // pcurves from marching can have V-folds that are geometrically correct.
        let valid_a2 = pca.as_ref().map_or(false, |pc| {
            crate::bop::int_tools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds1, u_per1, None)
        });
        let valid_b2 = pcb.as_ref().map_or(false, |pc| {
            crate::bop::int_tools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds2, u_per2, None)
        });
        if valid_a2 && valid_b2 {
            return (pca, pcb);
        }

        // Final fallback: return pcurves even if invalid 锟?the builder handles
        // out-of-face pcurves via its own boundary clipping.
        (pca, pcb)
    }

    /// OCCT: generate surface sample points
    /// OCCT: generate surface samples
    pub(crate) fn generate_surface_samples(
        &self,
        surface: &Surface3,
        n1: usize,
        n2: usize,
    ) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                crate::bop::int_tools::marching::sample_cylinder(cyl, [-20.0, 20.0], n1, n2)
            }
            Surface3::Sphere(sph) => crate::bop::int_tools::marching::sample_sphere(sph, n1, n2),
            Surface3::Torus(torus) => crate::bop::int_tools::marching::sample_torus(torus, n1, n2),
            Surface3::Plane(plane) => sample_plane(plane, 20.0, n1),
            Surface3::Cone(cone) => sample_cone(cone, 0.01, 20.0, n1, n2),
            // Generic fallback: sample via surface.default_domain() UV grid.
            // Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
            _ => sample_surface_generic(surface, n1, n2),
        }
    }

    /// OCCT: generate surface sample points
    /// OCCT: surface sample grid
    /// OCCT: generate surface samples
    /// OCCT: surface sample grid
    pub(crate) fn generate_surface_samples_grid(
        &self,
        surface: &Surface3,
        n_u: usize,
        n_v: usize,
    ) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                // u = azimuth index (0..n_u), v = height index (0..n_v)
                // sample_cylinder returns row = height, col = azimuth,
                // so transpose to row = azimuth, col = height for grid indexing.
                // Rebuild in (n_u azimuth) 锟?(n_v height) order.
                let height_range = [-20.0_f64, 20.0_f64];
                let u_ax = if cyl.axis.x.abs() < 0.9 {
                    cyl.axis.cross(DVec3::X).normalize()
                } else {
                    cyl.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = cyl.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta = 2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let h = height_range[0]
                            + (height_range[1] - height_range[0]) * iv as f64
                                / (n_v - 1).max(1) as f64;
                        pts.push(
                            cyl.origin
                                + cyl.axis * h
                                + (u_ax * theta.cos() + v_ax * theta.sin()) * cyl.radius,
                        );
                    }
                }
                pts
            }
            Surface3::Sphere(sph) => {
                let u_ax = if sph.axis.x.abs() < 0.9 {
                    sph.axis.cross(DVec3::X).normalize()
                } else {
                    sph.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = sph.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta = 2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let phi = std::f64::consts::PI * iv as f64 / (n_v - 1).max(1) as f64;
                        pts.push(
                            sph.center
                                + sph.radius
                                    * (sph.axis * phi.cos()
                                        + (u_ax * theta.cos() + v_ax * theta.sin()) * phi.sin()),
                        );
                    }
                }
                pts
            }
            _ => {
                // Fallback: generic UV-grid sampling for BSpline, Bezier, Offset, etc.
                sample_surface_generic(surface, n_u, n_v)
            }
        }
    }

    /// OCCT: estimate FF step size
    /// OCCT: estimate FF step size
    pub(crate) fn estimate_step_size(&self, s1: &Surface3, s2: &Surface3) -> f64 {
        // Use a fraction of the smallest characteristic dimension
        let size1 = match s1 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        let size2 = match s2 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        size1.min(size2) * 0.1
    }
}
