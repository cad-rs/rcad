fn numeric_intss_impl(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
    dom1_override: Option<[f64; 4]>,
    dom2_override: Option<[f64; 4]>,
    geom_tol_floor: Option<f64>,
) -> SurfaceSurfaceIntersection {
    let dom1 = s1.default_domain();
    let dom2 = s2.default_domain();

    // Clamp infinite domains. For cylinders the v-domain is [-鈭? +鈭瀅; we use
    // a range large enough to cover any practical intersection geometry.
    // 500 units covers large mechanical parts; callers should pass explicit
    // domain overrides for parts exceeding this range.
    const DOMAIN_CLAMP: f64 = 500.0;
    let clamp_dom = |[u0, u1, v0, v1]: [f64; 4]| -> [f64; 4] {
        [
            if u0.is_finite() { u0 } else { -DOMAIN_CLAMP },
            if u1.is_finite() { u1 } else { DOMAIN_CLAMP },
            if v0.is_finite() { v0 } else { -DOMAIN_CLAMP },
            if v1.is_finite() { v1 } else { DOMAIN_CLAMP },
        ]
    };
    let [u1_0, u1_1, v1_0, v1_1] = dom1_override.unwrap_or_else(|| clamp_dom(dom1));
    let [u2_0, u2_1, v2_0, v2_1] = dom2_override.unwrap_or_else(|| clamp_dom(dom2));

    // Pre-sample s2 on a grid for fast approximate distance computation.
    // For pairs where both surfaces are curved (e.g. cone 脳 cylinder) a
    // denser s2 grid captures the surface shape more accurately, improving
    // sign-change detection near the intersection band.  Cap at 80 to keep
    // the O(n1虏 脳 n2虏) distance evaluation tractable.
    let n2 = n.min(80);
    let mut s2_pts: Vec<DVec3> = Vec::with_capacity(n2 * n2);
    for i in 0..n2 {
        for j in 0..n2 {
            let u = u2_0 + (u2_1 - u2_0) * i as f64 / (n2 - 1).max(1) as f64;
            let v = v2_0 + (v2_1 - v2_0) * j as f64 / (n2 - 1).max(1) as f64;
            let p = s2.point_at(u, v);
            if p.is_finite() {
                s2_pts.push(p);
            }
        }
    }

    if s2_pts.is_empty() {
        return SurfaceSurfaceIntersection::default();
    }

    // Approximate distance from 3D point to s2 surface via closest sample,
    // optionally refined with analytic distance when the projection falls within
    // the face UV domain (only when bounded domain overrides are provided by
    // the caller 鈥?the boolean Pave-Filler path).  For unbounded surfaces
    // (direct `intersect_surfaces` calls without domain overrides) the
    // original nn-distance is used to avoid altering sign-change patterns for
    // large default domains.
    //
    // 鉁?OCCT-aligned: Implicit signed distance for Plane (IntSurf_Quadric / gp_Pln),
    // matching IntPatch_ImpPrmIntersection which uses F(P) = n路(P - origin) as the
    // signed-distance function.  Signed values enable zero-crossing detection when
    // the grid passes through the plane, rather than relying on threshold proximity.
    let approx_dist_to_s2 = |p: DVec3| -> f64 {
        if !p.is_finite() {
            return f64::INFINITY;
        }

        // Explicit signed distance for Plane.  Returned directly 鈥?does not fall
        // through to nn_d or proj distance so the sign is preserved for the grid.
        if let Surface3::Plane(pl) = s2 {
            return pl.normal.dot(p - pl.origin);
        }

        // Analytic unsigned distance for Cone and Cylinder surfaces (O(1), exact).
        // This catches narrow intersection bands that the pre-sampled s2 grid
        // misses, without requiring domain overrides.  The nn_d below still
        // provides bounded-domain awareness.
        let analytic_d = match s2 {
            Surface3::Cylinder(cyl) => {
                let axis = cyl.axis.normalize_or_zero();
                let to_pt = p - cyl.origin;
                let radial_len = (to_pt - axis * to_pt.dot(axis)).length();
                (radial_len - cyl.radius).abs()
            }
            Surface3::Cone(cone) => {
                let axis = cone.axis_dir();
                let apex = cone.apex;
                let to_pt = p - apex;
                let along = to_pt.dot(axis);
                let radial_len = (to_pt - axis * along).length();
                let (sin_h, cos_h) = cone.half_angle_rad.sin_cos();
                // Signed distance to the single nappe:
                //   F(P) = (r - R)路cos(尾) - z路sin(尾)   where z = along, r = radial_len
                //   F(P) = 0 on the cone surface
                //   |F(P)| = Euclidean distance to the cone lateral surface
                let d = (radial_len - cone.radius) * cos_h - along * sin_h;
                d.abs()
            }
            _ => f64::INFINITY,
        };

        let nn_d = s2_pts
            .iter()
            .map(|q| (*q - p).length())
            .fold(f64::INFINITY, f64::min);
        if let Some([u_min, u_max, v_min, v_max]) = dom2_override {
            let proj = closest_point_on_surface(s2, p, 16);
            if proj.distance.is_finite() {
                let (u, v) = proj.params;
                // Accept analytic distance when the projection lands within
                // the face UV boundary (with a small numerical slop to avoid
                // rejecting borderline projections at discrete grid edges).
                let uv_tol = TOLERANCE_ABS
                    .max(TOLERANCE_COORD_SUB * (u_max - u_min + v_max - v_min) * 0.5);
                if u >= u_min - uv_tol
                    && u <= u_max + uv_tol
                    && v >= v_min - uv_tol
                    && v <= v_max + uv_tol
                {
                    return proj.distance.min(nn_d).min(analytic_d);
                }
            }
        }
        nn_d.min(analytic_d)
    };

    // Threshold: treated as "on the surface" if distance < this.
    // Use the average cell size on s1 as a reference scale.
    let du = (u1_1 - u1_0) / n as f64;
    let dv = (v1_1 - v1_0) / n as f64;
    let p00 = s1.point_at(u1_0, v1_0);
    let p10 = s1.point_at(u1_0 + du, v1_0);
    let p01 = s1.point_at(u1_0, v1_0 + dv);
    let cell_size = (p10 - p00).length().max((p01 - p00).length()).max(TOLERANCE_MESH_LEGACY);
    let floor = geom_tol_floor
        .filter(|t| t.is_finite() && *t > 0.0)
        .unwrap_or(TOLERANCE_ABS)
        .max(TOLERANCE_ABS);
    let threshold = (cell_size * 2.0).max(floor);

    // Compute distance at each grid point
    let nn = n + 1; // grid has (n+1) 脳 (n+1) nodes
    let mut dist: Vec<f64> = Vec::with_capacity(nn * nn);
    let mut pts: Vec<DVec3> = Vec::with_capacity(nn * nn);
    let mut uvs: Vec<DVec2> = Vec::with_capacity(nn * nn); // UV for UV-space Newton refinement
    for i in 0..nn {
        for j in 0..nn {
            let u = u1_0 + (u1_1 - u1_0) * i as f64 / n as f64;
            let v = v1_0 + (v1_1 - v1_0) * j as f64 / n as f64;
            let p = s1.point_at(u, v);
            if !p.is_finite() {
                pts.push(DVec3::ZERO);
                dist.push(f64::INFINITY);
                uvs.push(DVec2::ZERO);
                continue;
            }
            pts.push(p);
            dist.push(approx_dist_to_s2(p));
            uvs.push(DVec2::new(u, v));
        }
    }

    let idx = |i: usize, j: usize| i * nn + j;

    // Find sign-change edges and interpolate crossing points
    let mut crossing_pts: Vec<DVec3> = Vec::new();

    // Helper: refine a 3D crossing point using UV-space Newton (OCCT IntPatch_TheSearchInside).
    // For grid-edge crossings, UV is interpolated from the edge's UV coordinates.
    let refine_crossing = |p3: DVec3, uv: DVec2| -> DVec3 {
        let refined = refine_uv_intersection(s1, s2, uv.x, uv.y, threshold);
        refined.map_or(p3, |r| r.0)
    };

    // Extract crossing point from grid edge indices + UV interpolation + UV Newton
    let edge_crossing = |a: usize, b: usize, t: f64| -> DVec3 {
        let t = t.clamp(0.0, 1.0);
        let p3 = pts[a].lerp(pts[b], t);
        let uv = uvs[a].lerp(uvs[b], t);
        refine_crossing(p3, uv)
    };

    // Horizontal edges: (i,j) 鈥?(i, j+1)
    for i in 0..nn {
        for j in 0..n {
            let a = idx(i, j);
            let b = idx(i, j + 1);
            let da = dist[a];
            let db = dist[b];
            // Zero crossing: signed distance changes sign.
            // Without signed distance, fall back to near/far threshold change.
            let zero_cross = da.is_finite() && db.is_finite() && da * db < 0.0;
            if zero_cross || (da < threshold) != (db < threshold) {
                let t = if (da - db).abs() < TOLERANCE_FLOAT_DEDUP {
                    0.5
                } else if zero_cross {
                    // Zero crossing: interpolate to d=0 (the true intersection)
                    da / (da - db)
                } else {
                    // Near/far boundary: interpolate to threshold
                    (threshold - da) / (db - da)
                };
                crossing_pts.push(edge_crossing(a, b, t));
            } else if da.abs() > threshold && db.abs() > threshold {
                // Grazing band: trough along the chord, both endpoints "outside".
                let (t, dmin) =
                    min_approx_dist_on_segment(pts[a], pts[b], &approx_dist_to_s2);
                if dmin <= threshold {
                    crossing_pts.push(edge_crossing(a, b, t));
                }
            } else if da.abs() <= threshold && db.abs() <= threshold {
                // Both near the surface.  If signs differ, zero_cross already caught it.
                if da * db >= 0.0 {
                    let pm = pts[a].lerp(pts[b], 0.5);
                    let dm = approx_dist_to_s2(pm);
                    if dm.abs() > threshold {
                        let (t1, d1) = min_approx_dist_on_segment(pts[a], pm, &approx_dist_to_s2);
                        if d1 <= threshold {
                            crossing_pts.push(refine_crossing(pts[a].lerp(pm, t1), uvs[a].lerp(uvs[b], 0.5 * t1)));
                        }
                        let (t2, d2) = min_approx_dist_on_segment(pm, pts[b], &approx_dist_to_s2);
                        if d2 <= threshold {
                            crossing_pts.push(refine_crossing(pm.lerp(pts[b], t2), uvs[a].lerp(uvs[b], 0.5 + 0.5 * t2)));
                        }
                    }
                }
            }
        }
    }

    // Vertical edges: (i,j) 鈥?(i+1, j)
    for i in 0..n {
        for j in 0..nn {
            let a = idx(i, j);
            let b = idx(i + 1, j);
            let da = dist[a];
            let db = dist[b];
            let zero_cross = da.is_finite() && db.is_finite() && da * db < 0.0;
            if zero_cross || (da < threshold) != (db < threshold) {
                let t = if (da - db).abs() < TOLERANCE_FLOAT_DEDUP {
                    0.5
                } else if zero_cross {
                    da / (da - db)
                } else {
                    (threshold - da) / (db - da)
                };
                crossing_pts.push(edge_crossing(a, b, t));
            } else if da.abs() > threshold && db.abs() > threshold {
                let (t, dmin) =
                    min_approx_dist_on_segment(pts[a], pts[b], &approx_dist_to_s2);
                if dmin <= threshold {
                    crossing_pts.push(edge_crossing(a, b, t));
                }
            } else if da.abs() <= threshold && db.abs() <= threshold {
                if da * db >= 0.0 {
                    let pm = pts[a].lerp(pts[b], 0.5);
                    let dm = approx_dist_to_s2(pm);
                    if dm.abs() > threshold {
                        let (t1, d1) = min_approx_dist_on_segment(pts[a], pm, &approx_dist_to_s2);
                        if d1 <= threshold {
                            crossing_pts.push(refine_crossing(pts[a].lerp(pm, t1), uvs[a].lerp(uvs[b], 0.5 * t1)));
                        }
                        let (t2, d2) = min_approx_dist_on_segment(pm, pts[b], &approx_dist_to_s2);
                        if d2 <= threshold {
                            crossing_pts.push(refine_crossing(pm.lerp(pts[b], t2), uvs[a].lerp(uvs[b], 0.5 + 0.5 * t2)));
                        }
                    }
                }
            }
        }
    }

    // OCCT-aligned: IntStart_SearchInside / IntPatch_TheSearchInside 鈥?
    // scan all grid points, run UV-constrained Newton-Raphson on F(u,v)=0.
    // Catches intersection points that grid-edge crossing misses (e.g. when
    // ALL signed distances are on one side of the surface).
    //
    // OCCT reference: IntStart_SearchInside.gxx lines 117-262.
    // For each sample point, creates a UV search box [uv卤du/2, uv卤dv/2]
    // (~Binf/Bsup in OCCT), runs math_FunctionSetRoot, checks |F| <= Tol.
    {
        let mut near_zero_pts: Vec<DVec3> = Vec::new();
        // UV grid cell half-size for the OCCT-aligned search box
        let half_du = du * 0.5;
        let half_dv = dv * 0.5;
        for i in 0..nn {
            for j in 0..nn {
                let k = idx(i, j);
                let d = dist[k];
                if !d.is_finite() {
                    continue;
                }
                let uv = uvs[k];
                // OCCT-aligned: search box = [uv 卤 du/2, uv 卤 dv/2] clamped to domain
                let u0 = (uv.x - half_du).max(u1_0);
                let u1_b = (uv.x + half_du).min(u1_1);
                let v0 = (uv.y - half_dv).max(v1_0);
                let v1_b = (uv.y + half_dv).min(v1_1);
                if u0 >= u1_b || v0 >= v1_b {
                    continue;
                }
                if let Some((p, _, _)) =
                    refine_uv_intersection_bounded(s1, s2, uv.x, uv.y, u0, u1_b, v0, v1_b, threshold)
                {
                    near_zero_pts.push(p);
                }
            }
        }
        // Dedup near-zero points against existing crossing points
        if !near_zero_pts.is_empty() {
            let dedup_sq = threshold * threshold;
            for p in near_zero_pts {
                if !crossing_pts.iter().any(|cp| (cp - p).length_squared() < dedup_sq) {
                    crossing_pts.push(p);
                }
            }
        }
    }

    let refine_tol = threshold.max(floor).max(TOLERANCE_ABS);

    // Localized closest approach: avoids returning nothing when crossings stay in a tangential pocket
    // but grid edges fail the XOR topology test because of coarse sampling noise.
    if crossing_pts.is_empty() {
        let mut min_d = f64::INFINITY;
        let mut max_d = f64::NEG_INFINITY;
        let mut min_p = DVec3::ZERO;
        for i in 0..dist.len() {
            let di = dist[i];
            if !di.is_finite() {
                continue;
            }
            if di < min_d {
                min_d = di;
                min_p = pts[i];
            }
            if di > max_d {
                max_d = di;
            }
        }
        let spread = max_d - min_d;
        let spread_tol = cell_size.max(threshold * 0.5).max(floor);
        if spread.is_finite() && min_d <= threshold && spread >= spread_tol {
            crossing_pts.push(min_p);
        }
    }

    crossing_pts = crossing_pts
        .into_iter()
        .map(|p| project_onto_intersection_tol(s1, s2, p, refine_tol))
        .collect();

    // Keep dedup **much** tighter than `refine_tol` / `threshold` so polylines are not collapsed
    // to a few clusters (e.g. Steinmetz cylinder-cylinder numeric paths).
    let dedup_tol = (cell_size * 0.08)
        .max(floor)
        .max(TOLERANCE_MESH_LEGACY);
    crossing_pts = dedup_points_spatial(crossing_pts, dedup_tol);

    let mut out = SurfaceSurfaceIntersection::default();

    if crossing_pts.is_empty() {
        return out;
    }

    if crossing_pts.len() == 1 {
        let p = crossing_pts[0];
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Point(p),
            pcurve_on_a: None,
            pcurve_on_b: None,
        });
        return out;
    }

    // BFS-greedy ordering: connect nearest unvisited neighbors into chains.
    // This works well for smooth curves; for self-intersecting surfaces it may
    // produce slightly wrong orderings near the crossing, which is acceptable
    // for topological boolean operations.
    let ordered = greedy_order_points(crossing_pts, floor);

    for chain in ordered {
        if chain.len() < 2 {
            continue;
        }

        // 鉁?OCCT-aligned: Try BSpline conversion for compact / evaluable
        //    representation (GeomInt_IntSS::MakeBSpline).
        let curve_3d = if chain.len() >= 4 {
            polyline_to_bspline(&chain, TOLERANCE_TOL_SCALE_MICRO)
                .map(|c| match c {
                    Curve3::BSpline(b) => SurfaceCurve::BSplineCurve(Box::new(b)),
                    _ => SurfaceCurve::Polyline(chain.clone()),
                })
                .unwrap_or_else(|| SurfaceCurve::Polyline(chain.clone()))
        } else {
            SurfaceCurve::Polyline(chain.clone())
        };

        let pca = polyline_pcurve_by_projection(&chain, s1);
        let pcb = polyline_pcurve_by_projection(&chain, s2);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d,
            pcurve_on_a: pca,
            pcurve_on_b: pcb,
        });
    }

    out
}

