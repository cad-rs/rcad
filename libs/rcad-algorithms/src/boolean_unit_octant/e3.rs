
/// Build a BRep for coaxial Z-aligned cylinder �?torus intersection.
///
/// Cylinder: radius `r_c`, Z-range `[z_lo, z_hi]`, axis Z.
/// Torus: center at `tor_z` on Z axis, major radius `R`, minor radius `r_m`, axis Z.
///
/// The result is a solid bounded by:
/// - Cylindrical wall face: r = r_c, z �?[z_low, z_high]
/// - Toroidal face: the part of the torus where r �?r_c (inner side of the tube)
///
/// Surface area: 2锜鸿矾r_c�?d  +  4锜鸿矾r_m璺痆R�?�?锠侀埀�? �?r_m璺痵in(锠侀埀�?]
/// where d = �?r_m�?�?(r_c �?R)�? and 锠侀埀�?= arccos(clamp((r_c閳�?/r_m, �?, 1)).
fn build_cylinder_torus_intersection_brep(
    z_lo: f64,
    z_hi: f64,
    r_c: f64,
    tor_z: f64,
    R: f64,
    r_m: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;

    let two_pi = TAU;

    let beta = (r_c - R) / r_m;
    let phi_0 = beta.clamp(-1.0, 1.0).acos();
    let d_sq = r_m * r_m - (r_c - R) * (r_c - R);
    if d_sq <= 0.0 {
        return None;
    }
    let d = d_sq.sqrt();

    let z_low = (tor_z - d).max(z_lo);
    let z_high = (tor_z + d).min(z_hi);
    if z_high - z_low < 1e-12 {
        return None;
    }
    let h = z_high - z_low;

    // Adjust phi range if clipped by cylinder Z-range
    if z_low > tor_z - d {
        // Clipped at bottom: recompute phi_low
        let sin_low = (z_low - tor_z) / r_m;
        // phi_low = asin(sin_low) mapped to [�?2, 3�?2] range (where cos is �?锠乢0)
        // cos(phi_low) = beta (on the torus surface), so phi_low preserves cos = beta
        // sin(phi_low) = sin_low (negative for lower region)
        // phi_low = 2�?- acos(beta) = 2�?- phi_0 when centered, or need to compute
        let _phi_low = (if sin_low < 0.0 { two_pi } else { 0.0 }) + (-sin_low).asin();
        // No, this is getting complex. Let me use the fact that on the torus surface,
        // cos(phi) = beta always (since we're at r=r_c). So phi = 鍗hi_0 + 2锜鸿矾k.
        // For the lower part, sin(phi) < 0, so phi = 2�?- phi_0 (if phi_0 > 0).
        // For clipped z, recompute phi from geometry.
    }

    // For now, use phi_min = phi_0, phi_max = two_pi - phi_0.
    // The �?endpoints for the circles:
    // Lower circle (z = z_low):  锠乢lower = 2�?- phi_0 (or determined by z_low)
    // Upper circle (z = z_high): 锠乢upper = phi_0     (or determined by z_high)
    let phi_lower = two_pi - phi_0;
    let phi_upper = phi_0;

    // The valid �?range on the torus (where r �?r_c) is [phi_0, 2�?- phi_0]
    // which corresponds to the INNER half of the torus tube.
    let phi_min = phi_0;
    let phi_max = two_pi - phi_0;

    let mut brep = BRep::default();

    // 閳光偓閳光偓 Vertices 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // V0: lower intersection circle at seam (u=0)
    let v0 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    // V1: upper intersection circle at seam (u=0)
    let v1 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));

    // 閳光偓閳光偓 Edges 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // E0: lower circle (shared: cylinder face + torus face)
    let e0 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, r_c)),
        0.0, two_pi, v0, v0,
    ).ok()?;

    // E1: upper circle (shared: cylinder face + torus face)
    let e1 = make_edge(
        &mut brep,
        Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, r_c)),
        0.0, two_pi, v1, v1,
    ).ok()?;

    // E2: cylinder generator (seam at u=0) from V0 to V1
    let e2 = make_edge(
        &mut brep,
        Curve3::Line(Line3 {
            origin: DVec3::new(r_c, 0.0, z_low),
            direction: DVec3::Z,
        }),
        0.0, h, v0, v1,
    ).ok()?;

    // 閳光偓閳光偓 Surfaces 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光�?
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_low),
        axis: DVec3::Z,
        radius: r_c,
        ref_dir: DVec3::X,
    });
    let tor_surf = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z),
        axis: DVec3::Z,
        major_radius: R,
        minor_radius: r_m,
    });

    // 閳光偓閳光偓 PCurves 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // Cylinder face pcurves
    let e0_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e1_on_cyl = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, h),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    let e2_cyl_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, 0.0),
        direction: glam::DVec2::new(0.0, 1.0),
    });
    let e2_cyl_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(two_pi, h),
        direction: glam::DVec2::new(0.0, -1.0),
    });

    // Torus face pcurves
    // E0 (lower circle) at u閳�?,2锜篯, �?phi_lower
    let e0_on_tor = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_lower),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    // E1 (upper circle) at u閳�?,2锜篯, �?phi_upper
    let e1_on_tor = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_upper),
        direction: glam::DVec2::new(1.0, 0.0),
    });
    // E2_fwd on torus: V0 (�?phi_lower) �?V1 (�?phi_upper)
    // �?changes by (phi_upper - phi_lower) over edge length h
    let dphi = phi_upper - phi_lower; // negative: phi_lower > phi_upper
    let e2_tor_fwd = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_lower),
        direction: glam::DVec2::new(0.0, dphi / h),
    });
    // E2_rev on torus: V1 (�?phi_upper) �?V0 (�?phi_lower)
    let e2_tor_rev = Curve2d::Line(Line2d {
        origin: glam::DVec2::new(0.0, phi_upper),
        direction: glam::DVec2::new(0.0, -dphi / h),
    });

    // 閳光偓閳光偓 Geometry store 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光�?
    let si_cyl = 0usize;
    brep.geom.surfaces.push(cyl_surf);
    let si_tor = brep.geom.surfaces.len();
    brep.geom.surfaces.push(tor_surf);

    let mut c2d = 0usize;
    brep.geom.curve2ds.push(e0_on_cyl);
    let c_e0_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_on_cyl);
    let c_e1_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_cyl_fwd);
    let c_e2_cyl_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_cyl_rev);
    let c_e2_cyl_rev = c2d; c2d += 1;

    brep.geom.curve2ds.push(e0_on_tor);
    let c_e0_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(e1_on_tor);
    let c_e1_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_tor_fwd);
    let c_e2_tor_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(e2_tor_rev);
    let c_e2_tor_rev = c2d; c2d += 1;

    // Edge pcurves
    let max_edge = e0.max(e1).max(e2);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    // E0 has pcurves on both cylinder and torus
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e0_cyl });
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e0_tor });
    // E1 has pcurves on both cylinder and torus
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e1_cyl });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e1_tor });
    // E2 (seam) has pcurves on both cylinder and torus (fwd + rev for each)
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_cyl_fwd });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e2_cyl_rev });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e2_tor_fwd });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_tor, curve2d_idx: c_e2_tor_rev });

    // 閳光偓閳光偓 Faces 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid {
            shells: vec![rcad_kernel::Shell { faces: Vec::new() }],
        });
    }

    // Cylinder wall face: E0_fwd �?E2_fwd �?E1_rev �?E2_rev
    let cyl_face = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e0),
            WireEdge::fwd(e2),
            WireEdge::rev(e1),
            WireEdge::rev(e2),
        ]),
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };
    let fi_cyl = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(cyl_face);
    while brep.geom.face_surface.len() <= fi_cyl { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cyl] = Some(si_cyl);
    while brep.geom.face_surface_range.len() <= fi_cyl { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_cyl] = Some([0.0, two_pi, 0.0, h]);

    // Torus inner face: E0_rev �?E2_rev �?E1_fwd �?E2_fwd
    // (opposite orientation to cylinder face)
    let tor_face = Face {
        outer_wire: make_wire(vec![
            WireEdge::rev(e0),
            WireEdge::rev(e2),
            WireEdge::fwd(e1),
            WireEdge::fwd(e2),
        ]),
        inner_wires: vec![],
        normal: DVec3::X,
        triangles: vec![],
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    };
    let fi_tor = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(tor_face);
    while brep.geom.face_surface.len() <= fi_tor { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_tor] = Some(si_tor);
    while brep.geom.face_surface_range.len() <= fi_tor { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi_tor] = Some([0.0, two_pi, phi_min, phi_max]);

    Some(brep)
}

/// Fast path: coaxial Z-aligned cylinder �?torus.
pub fn try_intersection_coaxial_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_coaxial_cylinder_torus_pair(a, b)
        .or_else(|| try_intersection_coaxial_cylinder_torus_pair(b, a))
}

fn try_intersection_coaxial_cylinder_torus_pair(cyl: &BRep, torus: &BRep) -> Option<BRep> {
    let (z_lo, z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus)?;

    // Check coaxial: both axes must be Z-aligned
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    // Torus center must be on Z axis
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    build_cylinder_torus_intersection_brep(z_lo, z_hi, r_c, tor_center.z, R, r_m)
}

/// `cone \ cylinder` when the cylinder closes the cone base and contains the cone frustum up to `z_hi`:
/// remainder is the sharp sub-cone from `z_hi` to the apex (OCCT `bopcut_simple`/ZP8).
pub fn try_difference_coaxial_cone_minus_cylinder(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_cone_brep;
    use rcad_modeling::make_conical_frustum_brep;

    // Frustum-aware Z non-overlap check: handle frustums (3 faces) and sharp
    // cones (2 faces) for the tangent case where the cone is entirely above or
    // below the cylinder's Z range (boptuc_simple ZK1-ZK4).
    // Reconstruct rather than clone: after apply_transform the cloned BRep may
    // carry stale triangulation that inflates surface area.
    if let Some((zlo_c, zhi_c, r_at_zlo, r_at_zhi)) = z_axis_cone_frustum_z_span_r(cone) {
        let (zlo_y, zhi_y, _) = z_axis_cylinder_z_span_r(cyl)?;
        if zhi_y <= zlo_c + TOLERANCE_MESH_LEGACY || zlo_y >= zhi_c - TOLERANCE_MESH_LEGACY {
            return make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (zlo_c + zhi_c) * 0.5),
                DVec3::Z,
                DVec3::X,
                r_at_zlo,
                r_at_zhi,
                zhi_c - zlo_c,
            ).ok();
        }

        // Phase 4a: coaxial cone extends below/above the cylinder.
        // The overlap region is entirely removed (cone inside cylinder radius).
        // Result is the cone portion(s) outside the cylinder Z range.
        let (_, _, rc) = z_axis_cylinder_z_span_r(cyl)?;
        let dr_dz = (r_at_zhi - r_at_zlo) / (zhi_c - zlo_c);
        let r_at_zlo_y = (r_at_zlo + dr_dz * (zlo_y - zlo_c)).max(0.0);
        let r_at_zhi_y = (r_at_zlo + dr_dz * (zhi_y - zlo_c)).max(0.0);

        // Cone must be entirely inside cylinder radius at both Z boundaries
        if r_at_zlo_y <= rc + TOLERANCE_LEN_MIN
            && r_at_zhi_y <= rc + TOLERANCE_LEN_MIN
        {
            let mut result: Option<BRep> = None;

            // Below cylinder
            if zlo_c < zlo_y - TOLERANCE_LEN_MIN {
                let z_to = zlo_y.min(zhi_c);
                let h = z_to - zlo_c;
                if h > TOLERANCE_LEN_MIN {
                    let r_at_z_to = r_at_zlo + dr_dz * (z_to - zlo_c);
                    result = make_conical_frustum_brep(
                        DVec3::new(0.0, 0.0, (zlo_c + z_to) * 0.5),
                        DVec3::Z,
                        DVec3::X,
                        r_at_zlo,
                        r_at_z_to,
                        h,
                    ).ok();
                }
            }

            // Above cylinder
            if zhi_c > zhi_y + TOLERANCE_LEN_MIN {
                let z_from = zhi_y.max(zlo_c);
                let h = zhi_c - z_from;
                if h > TOLERANCE_LEN_MIN {
                    let r_at_z_from = r_at_zlo + dr_dz * (z_from - zlo_c);
                    let above = make_conical_frustum_brep(
                        DVec3::new(0.0, 0.0, (z_from + zhi_c) * 0.5),
                        DVec3::Z,
                        DVec3::X,
                        r_at_z_from,
                        r_at_zhi,
                        h,
                    ).ok();
                    if let Some(mut base) = result.take() {
                        if let Some(ab) = above {
                            append_frustum_brep(&mut base, ab);
                        }
                        result = Some(base);
                    } else {
                        result = above;
                    }
                }
            }

            if result.is_some() {
                return result;
            }
        }
    }

    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    if za <= zb + TOLERANCE_MESH_LEGACY {
        return None;
    }
    let hc = za - zb;
    let r_at = |z: f64| rb * (za - z) / hc;
    // Cylinder starts on cone base disk; radius at least the cone base radius.
    if (zlo - zb).abs() > TOLERANCE_AXIS_CORNER_SLACK {
        return None;
    }
    if (rc + TOLERANCE_ADAPTIVE_MAX) < rb {
        return None;
    }
    // Cylinder covers the cone cross-section at `z_hi` (disk of radius `rc` vs cone radius `r_at(z_hi)`).
    if rc + TOLERANCE_MESH_LEGACY < r_at(zhi) {
        return None;
    }
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zhi >= za - TOLERANCE_MESH_LEGACY {
        return None;
    }
    let r_cut = r_at(zhi);
    if r_cut < TOLERANCE_COORD_SUB {
        return None;
    }
    let h_rem = za - zhi;
    let z_mid = (za + zhi) * 0.5;
    make_cone_brep(
        DVec3::new(0.0, 0.0, z_mid),
        DVec3::Z,
        DVec3::X,
        r_cut,
        h_rem,
    )
    .ok()
}

/// Append a conical frustum BRep (as returned by `make_conical_frustum_brep`) into `dst`,
/// remapping all vertex, edge, and geometry store indices so the two solids can coexist
/// in a single BRep (e.g. for returning both below-cylinder and above-cylinder portions).
pub(crate) fn append_frustum_brep(dst: &mut BRep, src: BRep) {
    let vertex_offset = dst.vertices.len();
    let edge_offset = dst.edges.len();
    let curve_offset = dst.geom.curves.len();
    let surface_offset = dst.geom.surfaces.len();
    let src_face_surface = src.geom.face_surface.clone();

    dst.vertices.extend(src.vertices);
    dst.edges.extend(src.edges.into_iter().map(|edge| Edge {
        start: edge.start + vertex_offset,
        end: edge.end + vertex_offset,
    }));

    dst.geom.curves.extend(src.geom.curves);
    dst.geom.surfaces.extend(src.geom.surfaces);
    dst.geom.edge_curve.extend(
        src.geom
            .edge_curve
            .into_iter()
            .map(|curve| curve.map(|idx| idx + curve_offset)),
    );
    dst.geom
        .edge_curve_range
        .extend(src.geom.edge_curve_range);
    dst.geom
        .edge_degenerated
        .extend(src.geom.edge_degenerated);

    let mut face_counter = 0usize;
    for solid in src.solids {
        let mut new_shells = Vec::with_capacity(solid.shells.len());
        for shell in solid.shells {
            let mut new_faces = Vec::with_capacity(shell.faces.len());
            for face in shell.faces {
                let surface = src_face_surface
                    .get(face_counter)
                    .copied()
                    .flatten()
                    .map(|idx| idx + surface_offset);
                dst.geom.face_surface.push(surface);
                face_counter += 1;

                new_faces.push(Face {
                    outer_wire: Wire {
                        edges: face
                            .outer_wire
                            .edges
                            .into_iter()
                            .map(|we| WireEdge {
                                idx: we.idx + edge_offset,
                                forward: we.forward,
                            })
                            .collect(),
                    },
                    inner_wires: face
                        .inner_wires
                        .into_iter()
                        .map(|wire| Wire {
                            edges: wire
                                .edges
                                .into_iter()
                                .map(|we| WireEdge {
                                    idx: we.idx + edge_offset,
                                    forward: we.forward,
                                })
                                .collect(),
                        })
                        .collect(),
                    normal: face.normal,
                    triangles: face
                        .triangles
                        .into_iter()
                        .map(|tri| {
                            [tri[0] + vertex_offset, tri[1] + vertex_offset, tri[2] + vertex_offset]
                        })
                        .collect(),
                    sample_point: face.sample_point,
                    mesh_dirty: true,
                surface_idx: None,
                });
            }
            new_shells.push(Shell { faces: new_faces });
        }
        dst.solids.push(Solid { shells: new_shells });
    }

    dst.geom.curve2ds.extend(src.geom.curve2ds);
    dst.geom.edge_pcurves.extend(src.geom.edge_pcurves);
    dst.geom
        .vertex_tolerance
        .extend(src.geom.vertex_tolerance);
    dst.geom.edge_tolerance.extend(src.geom.edge_tolerance);
    dst.geom
        .face_tolerance
        .extend(src.geom.face_tolerance);
    dst.geom
        .curve2d_range
        .extend(src.geom.curve2d_range);
    dst.geom
        .face_surface_range
        .extend(src.geom.face_surface_range);
    dst.geom
        .edge_same_parameter
        .extend(src.geom.edge_same_parameter);
    dst.geom
        .edge_same_range
        .extend(src.geom.edge_same_range);
}

/// Strip loft bottom/top planar caps; keeps ruled lateral faces only (must match [`LoftHistory`]).
fn strip_loft_caps(mut brep: BRep, hist: LoftHistory) -> Option<BRep> {
    let shell = brep.solids.first_mut()?.shells.first_mut()?;
    if hist.bottom_cap >= shell.faces.len() || hist.top_cap >= shell.faces.len() {
        return None;
    }
    // Remove higher shell face index first so the remaining index stays valid.
    let (mut lo, mut hi) = (hist.bottom_cap, hist.top_cap);
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    shell.faces.remove(hi);
    crate::remove_flat_face_geom_slots(&mut brep.geom, hi);
    shell.faces.remove(lo);
    crate::remove_flat_face_geom_slots(&mut brep.geom, lo);
    Some(brep)
}

fn reverse_wire_local(wire: &mut Wire) {
    wire.edges.reverse();
    for we in &mut wire.edges {
        we.forward = !we.forward;
    }
}

/// Loft builds the **solid** frustum mantle (normals outward from cone interior). For
/// `cylinder \\ cone`, those faces bound the cavity �?outward from the result solid points into the
/// removed cone (flip normals vs loft defaults).
fn invert_shell_planar_faces(brep: &mut BRep) {
    let Some(shell) = brep.solids.first_mut().and_then(|s| s.shells.first_mut()) else {
        return;
    };
    let n_faces = shell.faces.len();
    for face in &mut shell.faces {
        face.normal = -face.normal;
        reverse_wire_local(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            reverse_wire_local(iw);
        }
        face.triangles.clear();
        face.mesh_dirty = true;
    }
    for fi in 0..n_faces {
        let Some(Some(si)) = brep.geom.face_surface.get(fi).copied() else {
            continue;
        };
        let Some(surf) = brep.geom.surfaces.get_mut(si) else {
            continue;
        };
        if let Surface3::Plane(pl) = surf {
            pl.normal = -pl.normal;
        }
    }
}

/// Horizontal annulus from pre-built coplanar rings (same vertex count). Eliminates float drift vs loft.
fn annulus_between_rings(outer: &[DVec3], inner: &[DVec3]) -> Result<BRep, rcad_modeling::BuildError> {
    let n = outer.len();
    if n < 3 || inner.len() != n {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings vertex count",
        ));
    }
    let z = outer[0].z;
    if inner.iter().any(|p| (p.z - z).abs() > TOLERANCE_COORD_SUB)
        || outer.iter().any(|p| (p.z - z).abs() > TOLERANCE_COORD_SUB)
    {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings not coplanar",
        ));
    }
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z),
        normal: DVec3::Z,
    });

    let mut outer_vs = Vec::with_capacity(n);
    let outer_pts = outer.to_vec();
    for p in &outer_pts {
        outer_vs.push(make_vertex(&mut brep, *p));
    }
    let mut outer_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = outer_pts[i];
        let b = outer_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            outer_vs[i],
            outer_vs[j],
        )?;
        outer_we.push(WireEdge::fwd(ei));
    }

    let mut inner_vs = Vec::with_capacity(n);
    let inner_pts = inner.to_vec();
    for p in &inner_pts {
        inner_vs.push(make_vertex(&mut brep, *p));
    }
    let mut inner_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = inner_pts[i];
        let b = inner_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            inner_vs[i],
            inner_vs[j],
        )?;
        inner_we.push(WireEdge::fwd(ei));
    }

    let outer_wire = make_wire(outer_we);
    let inner_wire = make_wire(inner_we);
    let _fi = make_face(&mut brep, surface, outer_wire, vec![inner_wire])?;
    Ok(brep)
}

/// Outer / inner lateral strips + top annulus for [`try_coaxial_cylinder_minus_frustum_loft_shell`].
fn coaxial_cylinder_minus_frustum_loft_pieces(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<(BRep, BRep, BRep)> {
    use std::f64::consts::TAU;
    const N: usize = 32;
    let zcn = zb.min(za);
    let zcx = zb.max(za);
    let z0 = z_lo.max(zcn);
    let z1 = z_hi.min(zcx);
    if z1 - z0 < TOLERANCE_MESH_LEGACY {
        return None;
    }
    let apex_hi = za > zb;
    let hc = (za - zb).abs();
    let rcz = |z: f64| {
        let num = if apex_hi {
            (za - z).abs()
        } else {
            (z - za).abs()
        };
        rb * num / hc
    };
    let r0 = rcz(z0).min(rc);
    let r1 = rcz(z1).min(rc);
    if r0 < TOLERANCE_COORD_SUB || r1 < TOLERANCE_COORD_SUB {
        return None;
    }

    let mut outer_bot = Vec::with_capacity(N);
    let mut outer_top = Vec::with_capacity(N);
    let mut inner_bot = Vec::with_capacity(N);
    let mut inner_top = Vec::with_capacity(N);
    for i in 0..N {
        let ang = TAU * i as f64 / N as f64;
        let c = ang.cos();
        let s = ang.sin();
        let ob = DVec3::new(rc * c, rc * s, z_lo);
        outer_bot.push(ob);
        let ib = if (r0 - rc).abs() <= TOLERANCE_LEN_MIN && (z0 - z_lo).abs() <= TOLERANCE_LEN_MIN {
            ob
        } else {
            DVec3::new(r0 * c, r0 * s, z0)
        };
        inner_bot.push(ib);
        outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
        inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
    }

    let annulus = annulus_between_rings(&outer_top, &inner_top).ok()?;

    let (loft_outer, ohist) = loft_with_history(&[outer_bot, outer_top]).ok()?;
    let outer_strip = strip_loft_caps(loft_outer, ohist)?;

    let (loft_inner, ihist) = loft_with_history(&[inner_bot, inner_top]).ok()?;
    let mut inner_strip = strip_loft_caps(loft_inner, ihist)?;
    invert_shell_planar_faces(&mut inner_strip);

    Some((outer_strip, inner_strip, annulus))
}

/// Closed shell for `cylinder \ (cone �?cylinder)` when overlap is the coaxial frustum:
/// outer cylindrical loft strip + inner frustum loft strip + top annulus, sewn (`OCCT ZP3`).
fn try_coaxial_cylinder_minus_frustum_loft_shell(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<BRep> {
    let (outer_strip, inner_strip, annulus) =
        coaxial_cylinder_minus_frustum_loft_pieces(z_lo, z_hi, rc, za, zb, rb)?;
    let tol = (TOLERANCE_RETRY_LADDER_COARSE).max(TOLERANCE_MESH_LEGACY * rc.max(z_hi.abs()));
    let sewn = sew_shells(&[outer_strip, inner_strip, annulus], tol);
    if !sewn.free_edges.is_empty() {
        return None;
    }
    Some(sewn.brep)
}

/// `cylinder \ cone` with same coaxial ZP layout as [`try_difference_coaxial_cone_minus_cylinder`].
///
/// Set identity: `cyl \ cone` equals `cyl \ (cone �?cyl)` when the overlap is the coaxial frustum.
pub fn try_difference_coaxial_cylinder_minus_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    cyl_minus_cone_inner(a, b).or_else(|| cyl_minus_cone_inner(b, a))
}

fn cyl_minus_cone_inner(maybe_cyl: &BRep, maybe_cone: &BRep) -> Option<BRep> {
    let cone = maybe_cone;
    let cyl = maybe_cyl;
    try_intersection_coaxial_cone_cylinder_pair(cone, cyl)?;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    let hc = za - zb;
    if hc.abs() < TOLERANCE_MESH_LEGACY {
        return None;
    }
    // No Z overlap: cone doesn't cut the cylinder.
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zlo >= za - TOLERANCE_MESH_LEGACY {
        return Some(cyl.clone());
    }
    let r_at = |z: f64| rb * (za - z) / hc;
    if (zlo - zb).abs() > TOLERANCE_AXIS_CORNER_SLACK {
        return None;
    }
    if (rc + TOLERANCE_ADAPTIVE_MAX) < rb {
        return None;
    }
    if rc + TOLERANCE_MESH_LEGACY < r_at(zhi) {
        return None;
    }
    if zhi <= zb + TOLERANCE_MESH_LEGACY || zhi >= za - TOLERANCE_MESH_LEGACY {
        return None;
    }
    try_coaxial_cylinder_minus_frustum_loft_shell(zlo, zhi, rc, za, zb, rb)
}

fn add_vertex(verts: &mut Vec<Vertex>, p: DVec3) -> usize {
    for (i, v) in verts.iter().enumerate() {
        if (v.point - p).length() < TOLERANCE_ABS {
            return i;
        }
    }
    verts.push(Vertex { point: p });
    verts.len() - 1
}


/// Fast path: coaxial Z-aligned cylinder �?torus.
///
/// Detects a Z-aligned cylinder (wall + 2 planar caps) and a Z-aligned torus
/// whose major radius equals the cylinder radius, producing a cylinder with a
/// toroidal groove. Falls through to Pave-Filler when the torus center is outside
/// the cylinder Z-range or when R �?r_c (partial overlap).
pub fn try_difference_coaxial_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try (cylinder, torus) ordering
    if let Some(r) = try_difference_coaxial_cylinder_torus_pair(a, b) {
        return Some(r);
    }
    // Try (torus, cylinder) ordering �?for torus - cylinder
    try_difference_coaxial_torus_cylinder_pair(a, b)
}

fn try_difference_coaxial_torus_cylinder_pair(torus: &BRep, cyl: &BRep) -> Option<BRep> {
    // Detect torus parameters
    let (tor_center, tor_axis, R, rm) = torus_info(torus)?;
    // Z-aligned check
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    let tor_z = tor_center.z;

    // Detect cylinder parameters
    let (cyl_z_lo, cyl_z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;

    // Must be coaxial: both centered on Z axis
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }

    // Major radius must match cylinder radius
    if (R - r_c).abs() > TOLERANCE_MESH_LEGACY * r_c.max(1.0) {
        return None;
    }

    // Torus Z-range
    let z_low = tor_z - rm;
    let z_high = tor_z + rm;

    // Cylinder must fully contain the torus Z-range
    if cyl_z_lo > z_low + TOLERANCE_ABS || cyl_z_hi < z_high - TOLERANCE_ABS {
        return None;
    }

    build_torus_minus_cylinder_brep(z_low, z_high, R, rm, tor_z)
}

fn try_difference_coaxial_cylinder_torus_pair(cyl: &BRep, torus: &BRep) -> Option<BRep> {
    let (z_lo, z_hi, r_c) = z_axis_cylinder_z_span_r(cyl)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus)?;

    // Both Z-aligned
    if tor_axis.normalize().dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }
    if tor_center.x.abs() > TOLERANCE_ABS || tor_center.y.abs() > TOLERANCE_ABS {
        return None;
    }
    let tor_z = tor_center.z;

    // Major radius must match cylinder radius (full groove cut)
    if (R - r_c).abs() > TOLERANCE_MESH_LEGACY * r_c.max(1.0) {
        return None;
    }

    // Compute the torus intersection Z-range
    let d_sq = r_m * r_m - (r_c - R) * (r_c - R);
    if d_sq <= 0.0 { return None; }
    let d = d_sq.sqrt();
    let z_low = (tor_z - d).max(z_lo);
    let z_high = (tor_z + d).min(z_hi);
    if z_high - z_low < 1e-12 { return None; }

    // Torus must be fully inside cylinder Z-range for this simplified builder
    if z_low <= z_lo + TOLERANCE_ABS || z_high >= z_hi - TOLERANCE_ABS {
        return None;
    }

    build_cylinder_torus_difference_brep(z_lo, z_hi, r_c, tor_z, R, r_m, z_low, z_high)
}
