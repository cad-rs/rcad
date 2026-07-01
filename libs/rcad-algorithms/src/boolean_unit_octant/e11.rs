
/// bounded by a vertical plane parallel to the cylinder axis.
///
/// The result is a portion of the cylinder cut lengthwise by a plane parallel to
/// its axis. Only the side where `(P - center)璺痗lip_n 閳?-cut_dist` is kept.
/// `clip_n` must be a horizontal unit vector (z=0) pointing into the kept half.
/// `cut_dist` is the distance from the cylinder center to the cut plane (measured
/// in the direction opposite to `clip_n`, i.e., into the kept half).
///
/// Build a full Z-aligned cylinder BRep with selective inclusion of end caps.
///
/// When `include_bottom_cap` is false, the face at `center.z - height/2` is omitted.
/// When `include_top_cap` is false, the face at `center.z + height/2` is omitted.
/// The lateral face and the included cap faces are always created.
fn build_cylinder_brep_caps(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
    include_bottom_cap: bool,
    include_top_cap: bool,
) -> Result<BRep, rcad_modeling::BuildError> {
    let mut brep = make_cylinder_brep(center, axis, ref_dir, radius, height)?;
    if include_bottom_cap && include_top_cap {
        return Ok(brep);
    }
    let half_h = height * 0.5;
    let z_lo = center.z - half_h;
    let z_hi = center.z + half_h;
    let shell = &mut brep.solids[0].shells[0];
    shell.faces.retain(|face| {
        // Determine if this face is a planar cap by checking its outer_wire
        // vertex Z coordinates. A cap at z_lo has all vertices near z_lo.
        let zs: Vec<f64> = face.outer_wire.edges.iter()
            .filter_map(|we| brep.edges.get(we.idx))
            .flat_map(|e| [brep.vertices.get(e.start).map(|v| v.point.z),
                           brep.vertices.get(e.end).map(|v| v.point.z)])
            .flatten()
            .collect();
        if zs.is_empty() {
            return true;
        }
        let z_avg = zs.iter().sum::<f64>() / zs.len() as f64;
        let is_bottom = (z_avg - z_lo).abs() < 0.1 * height;
        let is_top = (z_avg - z_hi).abs() < 0.1 * height;
        if is_bottom && !include_bottom_cap { return false; }
        if is_top && !include_top_cap { return false; }
        true
    });
    Ok(brep)
}

/// bounded by a vertical plane parallel to the cylinder axis.
///
/// The result is a portion of the cylinder cut lengthwise by a plane parallel to
/// its axis. Only the side where `(P - center)璺痗lip_n 閳?-cut_dist` is kept.
/// `clip_n` must be a horizontal unit vector (z=0) pointing into the kept half.
/// `cut_dist` is the distance from the cylinder center to the cut plane (measured
/// in the direction opposite to `clip_n`, i.e., into the kept half).
///
/// When cut_dist=0, the cut plane passes through the cylinder axis and the result
/// is a clean half-cylinder. When cut_dist > 0, the cut plane is offset outward
/// from the axis, and more than half the cylinder is kept.
///
/// This is used when a box fully contains the cylinder in one XY axis but only
/// partially contains it in the other.
fn build_half_cylinder_intersection_brep(
    center: DVec3,   // cylinder center (cx, cy, cz)
    r: f64,          // radius
    h: f64,          // height
    clip_n: DVec3,   // horizontal unit normal pointing into the kept half
    cut_dist: f64,   // distance from center to cut plane (閳?, into kept half)
) -> BRep {
    debug_assert!(cut_dist >= 0.0 && cut_dist <= r + 1e-12,
        "cut_dist must be in [0, r]. Got {cut_dist} for r={r}");

    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    // Azimuth angle of clip_n in XY plane.
    let phi = clip_n.y.atan2(clip_n.x);

    // Half-angle of the kept arc: 浼?= arccos(-cut_dist/r).
    // For cut_dist=0 (center cut): 浼?= 锜?2 閳?铻杣 = 锜?(half-cylinder).
    // For cut_dist=r (full cylinder): 浼?= arccos(-1) = 锜?閳?铻杣 = 2锜?(full cylinder).
    let alpha = (-cut_dist / r).acos();

    // Vertices at the intersection of the cut plane with the cylinder surface.
    // V0/V3 at u = 锠?- 浼?(left generator), V1/V2 at u = 锠?+ 浼?(right generator).
    let (sa, ca) = alpha.sin_cos();
    let (sp, cp) = phi.sin_cos();

    // (cos(锠佸崵浼?, sin(锠佸崵浼?) = (cp*ca 閳?sp*sa, sp*ca 鍗?cp*sa)
    let cos_phi_minus_alpha = cp * ca + sp * sa;
    let sin_phi_minus_alpha = sp * ca - cp * sa;
    let cos_phi_plus_alpha = cp * ca - sp * sa;
    let sin_phi_plus_alpha = sp * ca + cp * sa;

    let v0_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_lo);
    let v1_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_lo);
    let v2_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_hi);
    let v3_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_hi);

    // --- Build BRep directly ---
    let mut brep = BRep::new();

    // Vertices
    let v0 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v0_p });
    let v1 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v1_p });
    let v2 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v2_p });
    let v3 = brep.vertices.len();
    brep.vertices.push(Vertex { point: v3_p });

    // Edge index helpers 閳?push an edge and return its index.
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
        // Ensure parallel vecs
        while brep.geom.edge_curve.len() <= idx {
            brep.geom.edge_curve.push(None);
            brep.geom.edge_curve_range.push(None);
            brep.geom.edge_degenerated.push(false);
        }
        brep.geom.edge_curve[idx] = Some(ci);
        brep.geom.edge_curve_range[idx] = Some([t0, t1]);
        brep.geom.edge_pcurves.push(Vec::new());
        idx
    };

    // E0: bottom arc (V1閳壐0 along kept side of bottom circle).
    // Circle3 with normal=-Z maps CylSurf u 閳?Circle3 -u, so V1(CylSurf u=锠?浼?
    // 閳?Circle3 u=-(锠?浼? and V0(CylSurf u=锠?浼? 閳?Circle3 u=-(锠?浼?.
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
        radius: r,
    });
    let e0 = next_curve(circle_bot, -phi - alpha, -phi + alpha, v1, v0);

    // E1: right generator (V1閳壐2)
    let line_r = Curve3::Line(Line3 {
        origin: v1_p,
        direction: DVec3::Z,
    });
    let e1 = next_curve(line_r, 0.0, h, v1, v2);

    // E2: top arc (V3閳壐2 along kept side of top circle).
    // Circle3 with normal=Z uses the same param as CylSurf, so V3(CylSurf u=锠?浼?
    // 閳?Circle3 u=锠?浼?and V2(CylSurf u=锠?浼? 閳?Circle3 u=锠?浼?
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
        radius: r,
    });
    let e2 = next_curve(circle_top, phi - alpha, phi + alpha, v3, v2);

    // E3: left generator (V0閳壐3)
    let line_l = Curve3::Line(Line3 {
        origin: v0_p,
        direction: DVec3::Z,
    });
    let e3 = next_curve(line_l, 0.0, h, v0, v3);

    // E4: cut bottom (V0閳壐1)
    let line_cb = Curve3::Line(Line3 {
        origin: v0_p,
        direction: v1_p - v0_p,
    });
    let e4 = next_curve(line_cb, 0.0, 1.0, v0, v1);

    // E5: cut top (V2閳壐3)
    let line_ct = Curve3::Line(Line3 {
        origin: v2_p,
        direction: v3_p - v2_p,
    });
    let e5 = next_curve(line_ct, 0.0, 1.0, v2, v3);

    // --- Surfaces ---
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(center.x, center.y, cz_lo),
        axis: DVec3::Z,
        radius: r,
        ref_dir: DVec3::X,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
    });
    let bot_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
    });
    let cut_plane = Surface3::Plane(Plane {
        origin: center - clip_n * cut_dist,
        normal: -clip_n, // outward from the solid
    });

    // Push surfaces and register face-surface mapping.
    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cyl_surf);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(top_plane);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(bot_plane);
    let si_cut = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cut_plane);

    // Helper to push a face.
    let mut push_face = |outer: Wire, surf_idx: usize, normal: DVec3| -> usize {
        let fi = if brep.solids.is_empty() {
            brep.solids.push(Solid {
                shells: vec![Shell { faces: Vec::new() }],
            });
            0
        } else {
            brep.solids[0].shells[0].faces.len()
        };
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: outer,
            inner_wires: Vec::new(),
            normal,
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(surf_idx);
        fi
    };

    // F0: Cylindrical face 閳?wire V0閳壐1閳壐2閳壐3閳壐0
    // E0 stored as V1閳壐0 閳?E0_rev = V0閳壐1
    // E1 stored as V1閳壐2 閳?E1_fwd = V1閳壐2
    // E2 stored as V3閳壐2 閳?E2_rev = V2閳壐3
    // E3 stored as V0閳壐3 閳?E3_rev = V3閳壐0
    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(e0),
            WireEdge::fwd(e1),
            WireEdge::rev(e2),
            WireEdge::rev(e3),
        ],
    };
    let _f0 = push_face(cyl_wire, si_cyl, clip_n);
    // Set face_surface_range for analytic SA: u 閳?[锠?浼? 锠?浼猐, v 閳?[0, h]
    while brep.geom.face_surface_range.len() <= _f0 {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[_f0] = Some([phi - alpha, phi + alpha, 0.0, h]);

    // F1: Top half-disk (normal=+Z). Wire: V2閳壐3 (cut top) 閳?V3閳壐2 (top arc)
    // E5_fwd: V2閳壐3, E2_fwd: V3閳壐2
    let top_wire = Wire {
        edges: vec![WireEdge::fwd(e5), WireEdge::fwd(e2)],
    };
    let _f1 = push_face(top_wire, si_top, DVec3::Z);

    // F2: Bottom half-disk (normal=-Z). Wire: V0閳壐1 (cut bottom) 閳?V1閳壐0 (bottom arc)
    // E4_fwd: V0閳壐1, E0_fwd: V1閳壐0
    let bot_wire = Wire {
        edges: vec![WireEdge::fwd(e4), WireEdge::fwd(e0)],
    };
    let _f2 = push_face(bot_wire, si_bot, -DVec3::Z);

    // F3: Cut face (normal=-clip_n). Wire: V0閳壐1閳壐2閳壐3閳壐0
    // E4_fwd: V0閳壐1, E1_fwd: V1閳壐2, E5_fwd: V2閳壐3, E3_rev: V3閳壐0
    let cut_wire = Wire {
        edges: vec![
            WireEdge::fwd(e4),
            WireEdge::fwd(e1),
            WireEdge::fwd(e5),
            WireEdge::rev(e3),
        ],
    };
    let _f3 = push_face(cut_wire, si_cut, -clip_n);

    brep
}

/// Shared helper for intersection: cyl_brep is the cylinder, box_brep is the box.
fn try_intersect_cylinder_box_one_dir(cyl_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let ca = try_cylinder_center_axis_radius_height(cyl_brep)?;
    let (cyl_bottom, cyl_axis, cyl_r, cyl_height) = ca;
    // try_cylinder_center_axis_radius_height returns the CylindricalSurface origin
    // (bottom of cylinder), not the geometric center. Compute the actual center.
    let cyl_center = cyl_bottom + cyl_axis * (cyl_height / 2.0);
    let bx = try_as_box(box_brep)?;

    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN {
        return None;
    }

    let z_idx = find_z_axis_index(&bx)?;
    let (u_idx, v_idx) = match z_idx {
        0 => (1, 2),
        1 => (2, 0),
        _ => (0, 1),
    };
    let u_ax = bx.axes[u_idx];
    let v_ax = bx.axes[v_idx];
    let eu = bx.extents[u_idx];
    let ev = bx.extents[v_idx];
    let ew = bx.extents[z_idx];
    let bc = bx.center;

    let cyl_z_lo = cyl_center.z - cyl_height / 2.0;
    let cyl_z_hi = cyl_center.z + cyl_height / 2.0;
    let box_z_lo = bc.z - ew;
    let box_z_hi = bc.z + ew;

    let tol = TOL * 10.0;

    // Intersection Z range.
    let inter_lo = cyl_z_lo.max(box_z_lo);
    let inter_hi = cyl_z_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol {
        return Some(BRep::default());
    }

    // Check XY containment.
    let cu = (cyl_center - bc).dot(u_ax);
    let cv = (cyl_center - bc).dot(v_ax);
    let full_u = cu - cyl_r >= -eu - tol && cu + cyl_r <= eu + tol;
    let full_v = cv - cyl_r >= -ev - tol && cv + cyl_r <= ev + tol;
    let tangent_splits = cylinder_box_tangent_split_thetas(
        full_u, full_v, cu, cv, u_ax, v_ax, eu, ev, cyl_r, tol,
    );

    if full_u && full_v {
        let h = inter_hi - inter_lo;
        let cz = inter_lo + h / 2.0;
        let center = DVec3::new(cyl_center.x, cyl_center.y, cz);
        return if tangent_splits.is_empty() {
            Some(make_cylinder_brep(center, cyl_axis, u_ax, cyl_r, h).ok()?)
        } else {
            Some(build_cylinder_box_intersection_brep_with_splits(
                center,
                cyl_r,
                h,
                &[],
                &tangent_splits,
            ))
        };
    }

    // Collect clip planes from all partially-contained axes.
    // Each partial axis may have 0, 1, or 2 clipped sides.
    let mut clip_planes: Vec<(DVec3, f64)> = Vec::new();
    for &(full_axis, cp, ax, ext) in &[(full_u, cu, u_ax, eu), (full_v, cv, v_ax, ev)] {
        if !full_axis {
            if cp - cyl_r < -ext - tol {
                clip_planes.push((ax, cp + ext));
            }
            if cp + cyl_r > ext + tol {
                clip_planes.push((-ax, ext - cp));
            }
        }
    }

    if clip_planes.is_empty() {
        return None;
    }

    let h = inter_hi - inter_lo;
    let cz = inter_lo + h / 2.0;
    let adj_center = DVec3::new(cyl_center.x, cyl_center.y, cz);

    // Single clip plane 閳?existing half-cylinder builder (backward compat, efficient).
    if clip_planes.len() == 1 {
        let (clip_dir, cut_dist) = clip_planes[0];
        // When cut_dist < -r the cylinder center is so far outside the box
        // that the cylinder doesn't reach the box face 閳?empty intersection.
        // (The multi-plane builder below handles this via zero-range valid 鑳?
        // intervals, but the half-cylinder builder asserts cut_dist 閳?0.)
        if cut_dist < -cyl_r + 1e-12 {
            return Some(BRep::default());
        }
        if tangent_splits.is_empty() {
            return Some(build_half_cylinder_intersection_brep(adj_center, cyl_r, h, clip_dir, cut_dist));
        }
    }

    // Multiple clip planes 閳?general multi-plane builder.
    Some(build_cylinder_box_intersection_brep_with_splits(
        adj_center,
        cyl_r,
        h,
        &clip_planes,
        &tangent_splits,
    ))
}