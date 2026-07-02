
/// Build BRep for cylinder 閳?torus (coaxial Z-aligned, R == r_c).
///
/// Result has 5 faces: lower cylindrical wall, torus groove, upper cylindrical wall,
/// bottom cap, and top cap.
fn build_cylinder_torus_difference_brep(
    z_lo: f64, z_hi: f64, r_c: f64,
    tor_z: f64, R: f64, r_m: f64,
    z_low: f64, z_high: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;

    let two_pi = TAU;

    let beta = (r_c - R) / r_m;
    let phi_0 = beta.clamp(-1.0, 1.0).acos();
    let phi_lower = two_pi - phi_0;
    let phi_upper = phi_0;
    let phi_min = phi_0;
    let phi_max = two_pi - phi_0;

    let mut brep = BRep::default();

    // 閳光偓閳光偓 Vertices 閳光偓閳光偓
    let v0 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_lo));
    let v1 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    let v2 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));
    let v3 = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_hi));

    // 閳光偓閳光偓 Edges 閳光偓閳光偓
    // E0: bottom cap circle at z=z_lo
    let e0 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_lo), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v0, v0).ok()?;
    // E1: lower intersection circle at z=z_low
    let e1 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_low), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v1, v1).ok()?;
    // E2: upper intersection circle at z=z_high
    let e2 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_high), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v2, v2).ok()?;
    // E3: top cap circle at z=z_hi
    let e3 = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_hi), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v3, v3).ok()?;

    // Seam edges
    let h_lower = z_low - z_lo;
    let e_seam_low = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_lo), direction: DVec3::Z,
    }), 0.0, h_lower, v0, v1).ok()?;

    let h_torus = z_high - z_low;
    let e_seam_torus = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h_torus, v1, v2).ok()?;

    let h_upper = z_hi - z_high;
    let e_seam_upper = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_high), direction: DVec3::Z,
    }), 0.0, h_upper, v2, v3).ok()?;

    // 閳光偓閳光偓 Surfaces 閳光偓閳光偓
    let surf_lower = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_lo), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });
    let surf_upper = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_high), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });
    let surf_torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z), axis: DVec3::Z,
        major_radius: R, minor_radius: r_m,
    });
    let surf_bot = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z_lo), normal: -DVec3::Z,
    });
    let surf_top = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z_hi), normal: DVec3::Z,
    });

    // Push surfaces
    let si_lower = 0usize;
    brep.geom.surfaces.push(surf_lower);
    let si_torus = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_torus);
    let si_upper = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_upper);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_bot);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_top);

    // 閳光偓閳光偓 Curve2Ds (pcurves) 閳光偓閳光偓
    let mut c2d = 0usize;
    // Lower wall pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e0_low = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h_lower), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e1_low = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_sl_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(two_pi, h_lower), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_sl_rev = c2d; c2d += 1;

    // Torus pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_lower), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e1_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_upper), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e2_tor = c2d; c2d += 1;
    let dphi = phi_upper - phi_lower;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_lower), direction: glam::DVec2::new(0.0, dphi / h_torus) }));
    let c_st_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, phi_upper), direction: glam::DVec2::new(0.0, -dphi / h_torus) }));
    let c_st_rev = c2d; c2d += 1;

    // Upper wall pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e2_up = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h_upper), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e3_up = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_su_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(two_pi, h_upper), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_su_rev = c2d; c2d += 1;

    // Cap pcurves (circles on planes)
    brep.geom.curve2ds.push(Curve2d::Circle(Circle2d { center: glam::DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: r_c  }));
    let c_e0_cap = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Circle(Circle2d { center: glam::DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: r_c  }));
    let c_e3_cap = c2d; c2d += 1;

    // 閳光偓閳光偓 Edge pcurves 閳光偓閳光偓
    let max_edge = e0.max(e1).max(e2).max(e3).max(e_seam_low).max(e_seam_torus).max(e_seam_upper);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    // E0 shared by lower wall + bottom cap
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_lower, curve2d_idx: c_e0_low });
    brep.geom.edge_pcurves[e0].push(PCurve { surface_idx: si_bot, curve2d_idx: c_e0_cap });
    // E1 shared by lower wall + torus
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_lower, curve2d_idx: c_e1_low });
    brep.geom.edge_pcurves[e1].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e1_tor });
    // E2 shared by torus + upper wall
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e2_tor });
    brep.geom.edge_pcurves[e2].push(PCurve { surface_idx: si_upper, curve2d_idx: c_e2_up });
    // E3 shared by upper wall + top cap
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_upper, curve2d_idx: c_e3_up });
    brep.geom.edge_pcurves[e3].push(PCurve { surface_idx: si_top, curve2d_idx: c_e3_cap });
    // Lower seam (lower wall only 閳?singular edge, appears fwd+rev in same face)
    brep.geom.edge_pcurves[e_seam_low].push(PCurve { surface_idx: si_lower, curve2d_idx: c_sl_fwd });
    brep.geom.edge_pcurves[e_seam_low].push(PCurve { surface_idx: si_lower, curve2d_idx: c_sl_rev });
    // Torus seam (torus only)
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_fwd });
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_rev });
    // Upper seam (upper wall only)
    brep.geom.edge_pcurves[e_seam_upper].push(PCurve { surface_idx: si_upper, curve2d_idx: c_su_fwd });
    brep.geom.edge_pcurves[e_seam_upper].push(PCurve { surface_idx: si_upper, curve2d_idx: c_su_rev });

    // 閳光偓閳光偓 Faces 閳光偓閳光偓
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid { shells: vec![rcad_kernel::Shell { faces: Vec::new() }] });
    }

    // 1. Lower cylindrical wall: e0_fwd 閳?seam_low_fwd 閳?e1_rev 閳?seam_low_rev
    let f_lower = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e0), WireEdge::fwd(e_seam_low),
            WireEdge::rev(e1), WireEdge::rev(e_seam_low),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_lower);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_lower);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h_lower]);

    // 2. Torus groove: e1_rev 閳?seam_torus_rev 閳?e2_fwd 閳?seam_torus_fwd
    let f_torus = Face {
        outer_wire: make_wire(vec![
            WireEdge::rev(e1), WireEdge::rev(e_seam_torus),
            WireEdge::fwd(e2), WireEdge::fwd(e_seam_torus),
        ]),
        inner_wires: vec![], normal: DVec3::X, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_torus);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_torus);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, phi_min, phi_max]);

    // 3. Upper cylindrical wall: e2_fwd 閳?seam_upper_fwd 閳?e3_rev 閳?seam_upper_rev
    let f_upper = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e2), WireEdge::fwd(e_seam_upper),
            WireEdge::rev(e3), WireEdge::rev(e_seam_upper),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_upper);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_upper);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h_upper]);

    // 4. Bottom cap: plane at z=z_lo, normal -Z, outer wire = e0_rev (CW when viewed from above 閳?normal -Z)
    let f_bot = Face {
        outer_wire: make_wire(vec![WireEdge::rev(e0)]),
        inner_wires: vec![], normal: -DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_bot);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_bot);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }

    // 5. Top cap: plane at z=z_hi, normal +Z, outer wire = e3_fwd (CCW when viewed from above)
    let f_top = Face {
        outer_wire: make_wire(vec![WireEdge::fwd(e3)]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_top);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_top);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }

    Some(brep)
}

/// Build BRep for torus 閳?cylinder (coaxial Z-aligned, R == r_c).
///
/// The result has 2 faces: the outer half of the torus surface (鑳?閳?[-锜?2, 锜?2])
/// connected to a cylindrical wall (r=R, z 閳?[z_low, z_high]).
/// The cylinder removes the inner-lower portion of the torus tube.
fn build_torus_minus_cylinder_brep(
    z_low: f64, z_high: f64,
    R: f64, rm: f64, tor_z: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::{PI, TAU};

    let two_pi = TAU;
    let r_c = R;
    let h = z_high - z_low; // = 2*rm

    let mut brep = BRep::default();

    // 閳光偓閳光偓 Vertices 閳光偓閳光偓
    let v_bot = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_low));
    let v_top = make_vertex(&mut brep, DVec3::new(r_c, 0.0, z_high));

    // 閳光偓閳光偓 Edges 閳光偓閳光偓
    // E_bot: bottom intersection circle at z=z_low, r=R
    let e_bot = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_low), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v_bot, v_bot).ok()?;
    // E_top: top intersection circle at z=z_high, r=R
    let e_top = make_edge(&mut brep, Curve3::Circle(Circle3 {
        center: DVec3::new(0.0, 0.0, z_high), normal: DVec3::Z, radius: r_c,
    }), 0.0, two_pi, v_top, v_top).ok()?;
    // Torus seam: 锠?0, 鑳?閳?[-锜?2, 锜?2] on torus surface (approximated as vertical line)
    let e_seam_torus = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h, v_bot, v_top).ok()?;
    // Cylinder seam: 锠?0, z 閳?[z_low, z_high]
    let e_seam_cyl = make_edge(&mut brep, Curve3::Line(Line3 {
        origin: DVec3::new(r_c, 0.0, z_low), direction: DVec3::Z,
    }), 0.0, h, v_bot, v_top).ok()?;

    // 閳光偓閳光偓 Surfaces 閳光偓閳光偓
    let surf_torus = Surface3::Torus(ToroidalSurface {
        center: DVec3::new(0.0, 0.0, tor_z), axis: DVec3::Z,
        major_radius: R, minor_radius: rm,
    });
    let surf_cyl = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(0.0, 0.0, z_low), axis: DVec3::Z, radius: r_c, ref_dir: DVec3::X,
    });

    let si_torus = 0usize;
    brep.geom.surfaces.push(surf_torus);
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf_cyl);

    // 閳光偓閳光偓 Curve2Ds (pcurves) 閳光偓閳光偓
    // Torus UV: u = 锠?(major angle) 閳?[0, 2锜篯, v = 鑳?(minor angle) 閳?[-锜?2, 锜?2]
    let theta_lo = -PI / 2.0;
    let theta_hi = PI / 2.0;
    let dtheta = theta_hi - theta_lo; // = 锜?

    let mut c2d = 0usize;
    // Torus pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_lo), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_bot_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_hi), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_top_tor = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_lo), direction: glam::DVec2::new(0.0, dtheta / h) }));
    let c_st_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, theta_hi), direction: glam::DVec2::new(0.0, -dtheta / h) }));
    let c_st_rev = c2d; c2d += 1;

    // Cylinder pcurves
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_bot_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h), direction: glam::DVec2::new(1.0, 0.0) }));
    let c_e_top_cyl = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, 0.0), direction: glam::DVec2::new(0.0, 1.0) }));
    let c_sc_fwd = c2d; c2d += 1;
    brep.geom.curve2ds.push(Curve2d::Line(Line2d { origin: glam::DVec2::new(0.0, h), direction: glam::DVec2::new(0.0, -1.0) }));
    let c_sc_rev = c2d; c2d += 1;

    // 閳光偓閳光偓 Edge pcurves 閳光偓閳光偓
    let max_edge = e_bot.max(e_top).max(e_seam_torus).max(e_seam_cyl);
    while brep.geom.edge_pcurves.len() <= max_edge {
        brep.geom.edge_pcurves.push(Vec::new());
    }
    // E_bot shared by torus + cylinder
    brep.geom.edge_pcurves[e_bot].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e_bot_tor });
    brep.geom.edge_pcurves[e_bot].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e_bot_cyl });
    // E_top shared by torus + cylinder
    brep.geom.edge_pcurves[e_top].push(PCurve { surface_idx: si_torus, curve2d_idx: c_e_top_tor });
    brep.geom.edge_pcurves[e_top].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_e_top_cyl });
    // Torus seam
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_fwd });
    brep.geom.edge_pcurves[e_seam_torus].push(PCurve { surface_idx: si_torus, curve2d_idx: c_st_rev });
    // Cylinder seam
    brep.geom.edge_pcurves[e_seam_cyl].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_sc_fwd });
    brep.geom.edge_pcurves[e_seam_cyl].push(PCurve { surface_idx: si_cyl, curve2d_idx: c_sc_rev });

    // 閳光偓閳光偓 Faces 閳光偓閳光偓
    if brep.solids.is_empty() {
        brep.solids.push(rcad_kernel::Solid { shells: vec![rcad_kernel::Shell { faces: Vec::new() }] });
    }

    // 1. Torus outer face: e_bot_fwd 閳?seam_torus_fwd 閳?e_top_rev 閳?seam_torus_rev
    let f_torus = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e_bot), WireEdge::fwd(e_seam_torus),
            WireEdge::rev(e_top), WireEdge::rev(e_seam_torus),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_torus);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_torus);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, theta_lo, theta_hi]);

    // 2. Cylinder wall face: e_bot_fwd 閳?seam_cyl_fwd 閳?e_top_rev 閳?seam_cyl_rev
    let f_cyl = Face {
        outer_wire: make_wire(vec![
            WireEdge::fwd(e_bot), WireEdge::fwd(e_seam_cyl),
            WireEdge::rev(e_top), WireEdge::rev(e_seam_cyl),
        ]),
        inner_wires: vec![], normal: DVec3::Z, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(f_cyl);
    while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi] = Some(si_cyl);
    while brep.geom.face_surface_range.len() <= fi { brep.geom.face_surface_range.push(None); }
    brep.geom.face_surface_range[fi] = Some([0.0, two_pi, 0.0, h]);

    Some(brep)
}

/// Fast path: coaxial Z-aligned cone-cone difference.
///
/// Detects two coaxial Z-aligned cones where one is fully nested inside the other.
/// The outer cone minus the inner cone produces a hollow conical frustum.
pub fn try_difference_coaxial_cone_minus_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    use rcad_modeling::make_conical_frustum_brep;
    // Fast path: non-overlapping Z ranges 閳?no volume intersection.
    // Coaxial conical frustums with disjoint Z ranges cannot overlap in 3D,
    // so a - b = a (even if the coincident face at the boundary would confuse
    // the Pave-Filler into removing it, e.g. bopcut_simple ZM7-ZN1).
    if let (Some(ai), Some(bi)) = (
        z_axis_cone_frustum_z_span_r(a),
        z_axis_cone_frustum_z_span_r(b),
    ) {
        if ai.1 <= bi.0 + TOLERANCE_MESH_LEGACY || ai.0 + TOLERANCE_MESH_LEGACY >= bi.1 {
            // Reconstruct from scratch rather than cloning 閳?a cloned BRep that
            // went through apply_transform may carry stale triangulation that
            // inflates surface area (boptuc_simple ZN1).
            return make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (ai.0 + ai.1) * 0.5),
                DVec3::Z,
                DVec3::X,
                ai.2,
                ai.3,
                ai.1 - ai.0,
            ).ok();
        }

        // Phase 4b: extending-past check 閳?one cone extends below/above the
        // other but is completely inside the other's radius over the overlap.
        if let Some(r) = try_extending_cone_minus_cone(
            ai.0, ai.1, ai.2, ai.3,
            bi.0, bi.1, bi.2, bi.3,
        ) {
            // a extends past b and a is inside b 閳?a - b = protruding part of a
            return Some(r);
        }
        // When b extends past a and b is inside a in the overlap,
        // a - b = outer cone a with a hole where b was.
        // Build directly via the generalized builder.
        if check_inner_inside_outer(bi.0, bi.1, bi.2, bi.3, ai.0, ai.1, ai.2, ai.3) {
            return build_conical_frustum_minus_frustum_brep(
                ai.0, ai.1, ai.2, ai.3,
                bi.0, bi.1, bi.2, bi.3,
            );
        }
    }
    try_difference_coaxial_cone_minus_cone_pair(a, b)
        .or_else(|| try_difference_coaxial_cone_minus_cone_pair(b, a))
}

/// Check if `inner` is completely inside `outer`'s radius at every Z in their overlap.
/// The inner may extend beyond the outer's Z range; only the overlap region is checked.
fn check_inner_inside_outer(
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
) -> bool {
    let tol = TOLERANCE_MESH_LEGACY;
    let z_olap_lo = zi_lo.max(zo_lo);
    let z_olap_hi = zi_hi.min(zo_hi);
    if z_olap_hi <= z_olap_lo + tol { return false; }

    let dri = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let dro = (ro_hi - ro_lo) / (zo_hi - zo_lo);

    let ri_at_lo = ri_lo + dri * (z_olap_lo - zi_lo);
    let ri_at_hi = ri_lo + dri * (z_olap_hi - zi_lo);
    let ro_at_lo = ro_lo + dro * (z_olap_lo - zo_lo);
    let ro_at_hi = ro_lo + dro * (z_olap_hi - zo_lo);

    ri_at_lo + tol < ro_at_lo && ri_at_hi + tol < ro_at_hi
}

/// Try extending-past case for cone-cone difference: when `inner` extends below
/// or above `outer` and is completely inside `outer`'s radius at every Z in the
/// overlap, the overlap region is entirely removed. The result is the truncated
/// inner frustum portion(s) outside the outer Z range.
fn try_extending_cone_minus_cone(
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
) -> Option<BRep> {
    // Must extend below or above outer
    if zi_lo >= zo_lo - TOLERANCE_ABS && zi_hi <= zo_hi + TOLERANCE_ABS {
        return None;
    }

    // Must have Z overlap (non-overlap is handled by the caller)
    let overlap_lo = zi_lo.max(zo_lo);
    let overlap_hi = zi_hi.min(zo_hi);
    if overlap_hi <= overlap_lo + TOLERANCE_LEN_MIN {
        return None;
    }

    // Inner radius at overlap boundaries
    let dr_i = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let ri_at_olo = ri_lo + dr_i * (overlap_lo - zi_lo);
    let ri_at_ohi = ri_lo + dr_i * (overlap_hi - zi_lo);

    // Outer radius at overlap boundaries
    let dr_o = (ro_hi - ro_lo) / (zo_hi - zo_lo);
    let ro_at_olo = ro_lo + dr_o * (overlap_lo - zo_lo);
    let ro_at_ohi = ro_lo + dr_o * (overlap_hi - zo_lo);

    // Inner must be inside outer at both overlap boundaries
    if ri_at_olo > ro_at_olo + TOLERANCE_LEN_MIN
        || ri_at_ohi > ro_at_ohi + TOLERANCE_LEN_MIN
    {
        return None;
    }

    use rcad_modeling::make_conical_frustum_brep;
    let mut result: Option<BRep> = None;

    // Portion below outer
    if zi_lo < zo_lo - TOLERANCE_LEN_MIN {
        let z_to = zo_lo.min(zi_hi);
        let h = z_to - zi_lo;
        if h > TOLERANCE_LEN_MIN {
            let ri_at_z_to = ri_lo + dr_i * (z_to - zi_lo);
            result = make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (zi_lo + z_to) * 0.5),
                DVec3::Z, DVec3::X,
                ri_lo, ri_at_z_to, h,
            ).ok();
        }
    }

    // Portion above outer
    if zi_hi > zo_hi + TOLERANCE_LEN_MIN {
        let z_from = zo_hi.max(zi_lo);
        let h = zi_hi - z_from;
        if h > TOLERANCE_LEN_MIN {
            let ri_at_z_from = ri_lo + dr_i * (z_from - zi_lo);
            let above = make_conical_frustum_brep(
                DVec3::new(0.0, 0.0, (z_from + zi_hi) * 0.5),
                DVec3::Z, DVec3::X,
                ri_at_z_from, ri_hi, h,
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

    result
}

fn try_difference_coaxial_cone_minus_cone_pair(outer: &BRep, inner: &BRep) -> Option<BRep> {
    // Extract cone frustum parameters from both operands
    // Using the same approach as z_axis_cylinder_z_span_r but for conical frustums
    let outer_info = z_axis_cone_frustum_z_span_r(outer)?;
    let inner_info = z_axis_cone_frustum_z_span_r(inner)?;

    let (zo_lo, zo_hi, ro_lo, ro_hi) = outer_info;
    let (zi_lo, zi_hi, ri_lo, ri_hi) = inner_info;

    // Check coaxial: both on Z axis
    // (already verified by z_axis_cone_frustum_z_span_r)

    // Inner cone must be fully inside outer cone
    if zi_lo < zo_lo - TOLERANCE_ABS || zi_hi > zo_hi + TOLERANCE_ABS {
        return None;
    }
    // Check inner cone radii are within outer cone radii at the same Z positions
    // Compute outer cone radius at inner cone Z bounds
    let h_o = zo_hi - zo_lo;
    if h_o <= TOLERANCE_MESH_LEGACY { return None; }
    let r_at = |z: f64| ro_lo + (ro_hi - ro_lo) * (z - zo_lo) / h_o;

    if ri_lo + TOLERANCE_MESH_LEGACY >= r_at(zi_lo) {
        return None; // inner cone touches or exceeds outer cone at bottom
    }
    if ri_hi + TOLERANCE_MESH_LEGACY >= r_at(zi_hi) {
        return None; // inner cone touches or exceeds outer cone at top
    }

    build_conical_frustum_minus_frustum_brep(
        zo_lo, zo_hi, ro_lo, ro_hi,
        zi_lo, zi_hi, ri_lo, ri_hi,
    )
}

/// Extract parameters from a Z-axis-aligned conical frustum.
/// Returns (z_lo, z_hi, r_at_zlo, r_at_zhi) where z_lo < z_hi.
/// The radii preserve the Z mapping (r_at_zlo is the radius at z_lo,
/// r_at_zhi is the radius at z_hi), unlike the old behavior that swapped
/// to r_lo <= r_hi.
fn z_axis_cone_frustum_z_span_r(brep: &BRep) -> Option<(f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() < 3 { return None; }

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                let apex = c.apex_point();
                if apex.x.abs() > TOLERANCE_ABS || apex.y.abs() > TOLERANCE_ABS {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if z_caps.len() < 2 { return None; }

    let z_lo = z_caps[0];
    let z_hi = z_caps[1];
    if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }

    // radius at z = tan(half_angle) * |z - apex.z|
    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let r_lo = tan_ha * (z_lo - apex.z).abs();
    let r_hi = tan_ha * (z_hi - apex.z).abs();

    if r_lo < TOLERANCE_COORD_SUB || r_hi < TOLERANCE_COORD_SUB {
        return None;
    }

    Some((z_lo, z_hi, r_lo, r_hi))
}

/// Like `z_axis_cone_frustum_z_span_r` but handles arbitrary XY translation.
/// Returns `(center_xy, z_lo, z_hi, r_lo, r_hi)` where `center_xy` is the cone's XY center.
fn detect_z_axis_cone_frustum(brep: &BRep) -> Option<(DVec2, f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() < 3 { return None; }

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if z_caps.len() < 2 { return None; }

    let z_lo = z_caps[0];
    let z_hi = z_caps[1];
    if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }

    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let r_lo = tan_ha * (z_lo - apex.z).abs();
    let r_hi = tan_ha * (z_hi - apex.z).abs();

    if r_lo < TOLERANCE_COORD_SUB || r_hi < TOLERANCE_COORD_SUB {
        return None;
    }

    let center_xy = DVec2::new(apex.x, apex.y);
    Some((center_xy, z_lo, z_hi, r_lo, r_hi))
}

/// Detect a Z-aligned cone (frustum or full cone with apex).
///
/// Returns `(center_xy, z_lo, z_hi, r_lo, r_hi)` 閳?the XY center, Z range, and
/// bottom/top radii. For a full cone one radius is near-zero (the apex).
pub(crate) fn detect_z_axis_cone(brep: &BRep) -> Option<(DVec2, f64, f64, f64, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;

    let mut cone_surf: Option<&ConicalSurface> = None;
    let mut caps: Vec<(DVec3, DVec3)> = Vec::new();

    let mut fi = 0usize;
    for _ in &sh.faces {
        let si = *brep.geom.face_surface.get(fi)?.as_ref()?;
        match brep.geom.surfaces.get(si)? {
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                if axis.cross(DVec3::Z).length() > TOLERANCE_AXIS_ALIGN {
                    return None;
                }
                cone_surf = Some(c);
            }
            Surface3::Plane(p) => {
                caps.push((p.origin, p.normal));
            }
            _ => {}
        }
        fi += 1;
    }

    let c = cone_surf?;
    if c.half_angle_rad.abs() < TOLERANCE_MESH_LEGACY { return None; }

    // Find Z-aligned planar caps
    let mut z_caps: Vec<f64> = caps.iter()
        .filter(|(_, n)| n.dot(DVec3::Z).abs() > 1.0 - TOLERANCE_AXIS_ALIGN)
        .map(|(p, _)| p.z)
        .collect();
    z_caps.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let apex = c.apex_point();
    let tan_ha = c.half_angle_rad.tan();
    let center_xy = DVec2::new(apex.x, apex.y);
    let tol = TOLERANCE_LEN_MIN;

    if z_caps.len() >= 2 {
        // Frustum: use the two caps
        let z_lo = z_caps[0];
        let z_hi = z_caps[1];
        if z_hi - z_lo < TOLERANCE_MESH_LEGACY { return None; }
        let r_lo = (tan_ha * (z_lo - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        let r_hi = (tan_ha * (z_hi - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        Some((center_xy, z_lo, z_hi, r_lo, r_hi))
    } else if z_caps.len() == 1 {
        // Full cone: one cap + the apex
        let z_cap = z_caps[0];
        let z_apex = apex.z;
        let r_cap = (tan_ha * (z_cap - apex.z).abs()).max(TOLERANCE_COORD_SUB);
        if r_cap < TOLERANCE_COORD_SUB { return None; }
        if z_apex.abs() < tol && z_cap.abs() < tol {
            return None; // Degenerate: apex and cap at origin
        }
        if z_apex > z_cap + tol {
            Some((center_xy, z_cap, z_apex, r_cap, TOLERANCE_COORD_SUB))
        } else if z_cap > z_apex + tol {
            Some((center_xy, z_apex, z_cap, TOLERANCE_COORD_SUB, r_cap))
        } else {
            None // Apex and cap at same Z
        }
    } else {
        None
    }
}

/// Build a tessellated BRep for coaxial `cone 閳?cylinder` where the cylinder
/// fills the bottom (or top) of the cone (same XY center, Z-aligned).
///
/// The cylinder occupies `[z_cyl_lo, z_cyl_hi]` with constant radius `cyl_r`,
/// and the cone continues from the cylinder top to `[z_con_hi]` with radius
/// varying from `r_at_cyl_hi` to `r_con_hi`.  At the interface a horizontal
/// annular ring connects the cylinder wall to the cone wall.
fn build_coaxial_cone_cylinder_union_tessellated(
    center_xy: DVec2,
    z_cyl_lo: f64, z_cyl_hi: f64, cyl_r: f64,
    z_con_hi: f64, r_con_hi: f64,
    r_at_cyl_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cyl_r < tol || z_cyl_hi <= z_cyl_lo + tol || z_con_hi <= z_cyl_hi + tol { return None; }

    let n_slices = 16usize;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    let circle_poly = |r: f64| -> Vec<DVec2> {
        let mut poly = Vec::with_capacity(n_arc + 1);
        for k in 0..=n_arc {
            let ang = tau * k as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            poly.push(DVec2::new(r * c, r * s));
        }
        poly
    };

    // 1. Cylinder wall from z_cyl_lo to z_cyl_hi
    let cyl_poly = circle_poly(cyl_r);
    add_wall_section(&mut add_v, &mut faces, &cyl_poly, z_cyl_lo, z_cyl_hi, n_slices, &to_world, &empty_wire);

    // 2. Bottom cap at z_cyl_lo
    add_cap_face(&mut add_v, &mut faces, &cyl_poly, z_cyl_lo, -DVec3::Z, &to_world, &empty_wire);

    // 3. Interface ring at z_cyl_hi (annulus: outer=cyl_r, inner=r_at_cyl_hi)
    if cyl_r > r_at_cyl_hi + tol {
        let annulus_pts: Vec<DVec2> = (0..n_arc).map(|i| {
            let ang = tau * i as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(cyl_r * c, cyl_r * s)
        }).collect();
        let inner_pts: Vec<DVec2> = (0..n_arc).map(|i| {
            let ang = tau * i as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_at_cyl_hi * c, r_at_cyl_hi * s)
        }).collect();

        let mut idx = Vec::with_capacity(2 * n_arc);
        for p in &annulus_pts { idx.push(add_v(to_world(p.x, p.y, z_cyl_hi))); }
        for p in &inner_pts { idx.push(add_v(to_world(p.x, p.y, z_cyl_hi))); }

        let mut tris = Vec::with_capacity(n_arc * 2);
        for i in 0..n_arc {
            let k = (i + 1) % n_arc;
            tris.push([idx[i], idx[k], idx[n_arc + k]]);
            tris.push([idx[i], idx[n_arc + k], idx[n_arc + i]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    // 4. Cone wall from z_cyl_hi to z_con_hi with varying radius
    let r_delta = r_con_hi - r_at_cyl_hi;
    let z_delta = z_con_hi - z_cyl_hi;
    if z_delta > tol {
        let dz = z_delta / n_slices as f64;
        for i in 0..n_slices {
            let za = z_cyl_hi + dz * i as f64;
            let zb = z_cyl_hi + dz * (i + 1) as f64;
            let ra = (r_at_cyl_hi + r_delta * (i as f64 / n_slices as f64)).max(0.0);
            let rb = (r_at_cyl_hi + r_delta * ((i + 1) as f64 / n_slices as f64)).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(ra * c, ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(rb * c, rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        }
    }

    // 5. Top cap at z_con_hi (if cone top radius > 0)
    if r_con_hi > tol {
        let top_poly = circle_poly(r_con_hi);
        add_cap_face(&mut add_v, &mut faces, &top_poly, z_con_hi, DVec3::Z, &to_world, &empty_wire);
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build a BRep with analytical surfaces for `cylinder 閳?cone` when the
/// cylinder is wider than the cone at every Z level (cone sits entirely inside
/// the cylinder and only protrudes above the cylinder's top face).
///
/// Unlike the tessellated builder `build_coaxial_cone_cylinder_union_tessellated`,
/// this produces a BRep with proper Surface3 entries and pcurves so that the
/// PaveFiller can process it in subsequent boolean operations.
fn try_union_coaxial_cone_cylinder_one_dir(cone_brep: &BRep, cyl_brep: &BRep) -> Option<BRep> {
    let (cone_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone(cone_brep)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = try_cylinder_center_axis_radius_height(cyl_brep)?;

    // Cylinder must be Z-aligned
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let cyl_z_lo = cyl_bottom.z;
    let cyl_z_hi = cyl_bottom.z + cyl_height;
    let tol = TOLERANCE_LEN_MIN;

    // Coaxial check: same XY center
    if (cone_xy.x - cyl_bottom.x).abs() > tol || (cone_xy.y - cyl_bottom.y).abs() > tol {
        return None;
    }

    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_cyl_lo = (cr_lo + dr_dz * (cyl_z_lo - cz_lo)).max(0.0);
    let r_at_cyl_hi = (cr_lo + dr_dz * (cyl_z_hi - cz_lo)).max(0.0);

    // Check: cylinder fills the bottom of the cone (cylinder at cone base)
    // Cylinder radius matches cone radius at cyl_z_lo, and cylinder is within cone Z-range
    let at_base = (cyl_r - r_at_cyl_lo).abs() <= cyl_r * 1e-6 + tol
        && r_at_cyl_hi <= cyl_r + tol
        && cyl_z_lo >= cz_lo - tol
        && cyl_z_hi <= cz_hi + tol;

    if at_base && r_at_cyl_lo > tol {
        return build_coaxial_cone_cylinder_union_tessellated(
            cone_xy,
            cyl_z_lo, cyl_z_hi, cyl_r,
            cz_hi, cr_hi,
            r_at_cyl_hi,
        );
    }

    // Check: cylinder fills the top of the cone
    let at_top = (cyl_r - r_at_cyl_hi).abs() <= cyl_r * 1e-6 + tol
        && r_at_cyl_lo <= cyl_r + tol
        && cyl_z_lo >= cz_lo - tol
        && cyl_z_hi <= cz_hi + tol;

    if at_top && r_at_cyl_hi > tol {
        // TODO: implement cylinder-on-top case
        // Needs: cone wall (cz_lo 閳?cyl_z_lo), bottom cap, annular ring at cyl_z_lo,
        // cylinder wall (cyl_z_lo 閳?cyl_z_hi), top cap
        return None;
    }

    // Case 3: cone entirely inside the cylinder (cylinder wider than cone at every Z).
    // NOTE: build_cylinder_cone_union_wider_cyl produces a correct analytical Union BRep
    // (SA=697 vs expected ~708) with proper surfaces, edges, pcurves, and face normals.
    // However, the PaveFiller Difference over-counts SA (867 vs 727, +19% vs 15% tolerance)
    // because the BooleanBuilder classifies the analytical faces differently than PaveFiller-
    // produced faces.  Fix needs PaveFiller/Builder layer, not BRep construction.
    // let cone_inside = r_at_cyl_lo <= cyl_r + tol && r_at_cyl_hi <= cyl_r + tol
    //     && cz_lo >= cyl_z_lo - tol && cz_hi >= cyl_z_hi - tol;
    // if cone_inside && cyl_r > tol && r_at_cyl_hi > tol && cr_hi > tol {
    //     return build_cylinder_cone_union_wider_cyl(
    //         cone_xy, cyl_z_lo, cyl_z_hi, cyl_r,
    //         cz_hi, cr_hi, r_at_cyl_hi,
    //     );
    // }

    None
}

/// Build an analytical BRep for `cylinder 閳?cone` when the cylinder is wider
/// than the cone at every Z level.  Uses the same edge/face/surface/pcurve
/// pattern as `build_cylinder_box_difference_full_wall` (proven PaveFiller-compatible).
fn build_cylinder_cone_union_wider_cyl(
    center_xy: DVec2,
    z_cyl_lo: f64, z_cyl_hi: f64, cyl_r: f64,
    z_con_hi: f64, r_con_hi: f64,
    r_top: f64,
) -> Option<BRep> {
    use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
    use rcad_kernel::PCurve;
    use std::f64::consts::TAU;
    let h = z_cyl_hi - z_cyl_lo;
    let h_con = (z_con_hi - z_cyl_hi).max(1e-10);
    if cyl_r < 1e-10 || r_con_hi < 1e-10 { return None; }
    let (cx, cy) = (center_xy.x, center_xy.y);
    let two_pi = TAU;

    let mut brep = BRep::new();
    brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    macro_rules! push_edge {
        ($c:expr, $t0:expr, $t1:expr, $start:expr, $end:expr) => {{
            let idx = brep.edges.len();
            brep.edges.push(Edge { start: $start, end: $end });
            let ci = brep.geom.curves.len();
            brep.geom.curves.push($c);
            while brep.geom.edge_curve.len() <= idx {
                brep.geom.edge_curve.push(None);
                brep.geom.edge_curve_range.push(None);
                brep.geom.edge_degenerated.push(false);
            }
            brep.geom.edge_curve[idx] = Some(ci);
            brep.geom.edge_curve_range[idx] = Some([$t0, $t1]);
            brep.geom.edge_pcurves.push(Vec::new());
            idx
        }};
    }

    // Vertices
    let v0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + cyl_r, cy, z_cyl_lo) });
    let v1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + cyl_r, cy, z_cyl_hi) });
    let v2 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + r_top, cy, z_cyl_hi) });
    let v3 = brep.vertices.len();
    brep.vertices.push(Vertex { point: DVec3::new(cx + r_con_hi, cy, z_con_hi) });

    // Edges
    let e0 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_lo), normal: -DVec3::Z, radius: cyl_r }), 0.0, two_pi, v0, v0);
    let e1 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z, radius: cyl_r }), 0.0, two_pi, v1, v1);
    let e2 = push_edge!(Curve3::Line(Line3 { origin: brep.vertices[v0].point, direction: DVec3::Z }), 0.0, h, v0, v1);
    let e3 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z, radius: r_top }), 0.0, two_pi, v2, v2);
    let coned = brep.vertices[v3].point - brep.vertices[v2].point;
    let e4 = push_edge!(Curve3::Line(Line3 { origin: brep.vertices[v2].point, direction: coned.normalize_or_zero() }), 0.0, coned.length(), v2, v3);
    let e5 = push_edge!(Curve3::Circle(Circle3 { center: DVec3::new(cx,cy,z_con_hi), normal: DVec3::Z, radius: r_con_hi }), 0.0, two_pi, v3, v3);

    // Surfaces
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface { origin: DVec3::new(cx,cy,z_cyl_lo), axis: DVec3::Z, radius: cyl_r, ref_dir: DVec3::X }));
    let si_cone = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Cone(ConicalSurface { apex: DVec3::new(cx,cy,z_cyl_hi), axis: DVec3::Z, radius: r_top, half_angle_rad: ((r_top - r_con_hi)/h_con).atan() }));
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_cyl_lo), normal: -DVec3::Z }));
    let si_step = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_cyl_hi), normal: DVec3::Z }));
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(Surface3::Plane(Plane { origin: DVec3::new(cx,cy,z_con_hi), normal: DVec3::Z }));

    // PCurves
    {
        let g = &mut brep.geom;
        let mut c2 = |c: Curve2d| -> usize { let i = g.curve2ds.len(); g.curve2ds.push(c); i };
        let p0w = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(two_pi,0.0) }));
        let p0b = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: cyl_r  }));
        let p1w = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,h), direction: DVec2::new(two_pi,0.0) }));
        let p1s = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: cyl_r  }));
        let p2f = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(0.0,h) }));
        let p2r = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,h), direction: DVec2::new(0.0,-h) }));
        let p3n = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(two_pi,0.0) }));
        let p3s = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: r_top  }));
        let p4f = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,0.0), direction: DVec2::new(0.0,1.0) }));
        let p4r = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,1.0), direction: DVec2::new(0.0,-1.0) }));
        let p5n = c2(Curve2d::Line(Line2d { origin: DVec2::new(0.0,1.0), direction: DVec2::new(two_pi,0.0) }));
        let p5t = c2(Curve2d::Circle(Circle2d { center: DVec2::ZERO, x_dir: DVec2::X, y_dir: DVec2::Y, radius: r_con_hi  }));
        g.edge_pcurves[e0].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p0w }, PCurve { surface_idx: si_bot, curve2d_idx: p0b }]);
        g.edge_pcurves[e1].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p1w }, PCurve { surface_idx: si_step, curve2d_idx: p1s }]);
        g.edge_pcurves[e2].extend(vec![PCurve { surface_idx: si_cyl, curve2d_idx: p2f }, PCurve { surface_idx: si_cyl, curve2d_idx: p2r }]);
        g.edge_pcurves[e3].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p3n }, PCurve { surface_idx: si_step, curve2d_idx: p3s }]);
        g.edge_pcurves[e4].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p4f }, PCurve { surface_idx: si_cone, curve2d_idx: p4r }]);
        g.edge_pcurves[e5].extend(vec![PCurve { surface_idx: si_cone, curve2d_idx: p5n }, PCurve { surface_idx: si_top, curve2d_idx: p5t }]);
    }

    // Faces (with face_surface_range for bounded surfaces)
    // Face helper: normal must be set for planar faces (DVec3::ZERO for curved)
    let mut push_face = |b: &mut BRep, si: usize, outer: Vec<WireEdge>, inner: Option<Vec<WireEdge>>, norm: DVec3, uv_range: Option<[f64; 4]>| {
        let fi = b.solids[0].shells[0].faces.len();
        while b.geom.face_surface.len() <= fi { b.geom.face_surface.push(None); b.geom.face_surface_range.push(None); b.geom.face_tolerance.push(0.0); }
        b.geom.face_surface[fi] = Some(si);
        if let Some(range) = uv_range { b.geom.face_surface_range[fi] = Some(range); }
        b.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: outer },
            inner_wires: inner.map(|e| vec![Wire { edges: e }]).unwrap_or_default(),
            normal: norm, triangles: vec![], sample_point: None, mesh_dirty: true,
                surface_idx: None,
        });
    };
    push_face(&mut brep, si_cyl, vec![WireEdge::rev(e0), WireEdge::fwd(e2), WireEdge::fwd(e1), WireEdge::rev(e2)], None, DVec3::ZERO, Some([0.0, two_pi, 0.0, h]));
    push_face(&mut brep, si_cone, vec![WireEdge::fwd(e3), WireEdge::fwd(e4), WireEdge::rev(e5), WireEdge::rev(e4)], None, DVec3::ZERO, Some([0.0, two_pi, 0.0, 1.0]));
    push_face(&mut brep, si_bot, vec![WireEdge::rev(e0)], None, -DVec3::Z, None);
    push_face(&mut brep, si_step, vec![WireEdge::fwd(e1)], Some(vec![WireEdge::rev(e3)]), DVec3::Z, None);
    push_face(&mut brep, si_top, vec![WireEdge::rev(e5)], None, DVec3::Z, None);

    // DEBUG: print BRep structure
    if std::env::var("RCAD_DEBUG_BREP").is_ok() {
        eprintln!("=== BUILD_CYLINDER_CONE_UNION ===");
        eprintln!("vertices: {}", brep.vertices.len());
        eprintln!("edges: {} curves: {}", brep.edges.len(), brep.geom.curves.len());
        eprintln!("edge_pcurves: {:?}", brep.geom.edge_pcurves.iter().map(|v| v.len()).collect::<Vec<_>>());
        eprintln!("surfaces: {}", brep.geom.surfaces.len());
        eprintln!("faces: {}", brep.solids[0].shells[0].faces.len());
        for (fi, f) in brep.solids[0].shells[0].faces.iter().enumerate() {
            let si = brep.geom.face_surface.get(fi).and_then(|o| *o);
            eprintln!("  face[{}]: surface_idx={:?} outer_edges={} inner_wires={}",
                fi, si, f.outer_wire.edges.len(), f.inner_wires.len());
        }
    }

    Some(brep)
}

/// Bidirectional wrapper for coaxial cone + cylinder Union.
pub fn try_union_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_coaxial_cone_cylinder_one_dir(a, b)
        .or_else(|| try_union_coaxial_cone_cylinder_one_dir(b, a))
}