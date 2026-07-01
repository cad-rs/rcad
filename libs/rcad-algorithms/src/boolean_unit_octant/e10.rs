
/// Full-wall difference case: box entirely within the cylinder cross-section.
/// The wall is the full cylinder. Caps are donuts (full circle with box polygon
/// as inner hole). Side faces on each clip plane are created separately.
fn build_cylinder_box_difference_full_wall(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let mut brep = BRep::new();
    brep.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // Helper macro to push a curve and edge (avoids closure borrow issues).
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

    let canon = |t: f64| if t >= two_pi - 1e-12 { 0.0 } else { t };

    // ---- 1. Pre-compute plane info & theta intersection points ----
    let info: Vec<(f64, f64, DVec3, f64)> = clip_planes.iter().map(|&(n, d)| {
        let alpha = (-d / r).clamp(-1.0, 1.0).acos();
        let phi = n.y.atan2(n.x);
        (phi, alpha, n, d)
    }).collect();

    // Collect unique theta endpoints where clip planes meet the circle.
    struct VEntry { theta: f64, lo_idx: usize, hi_idx: usize }
    let mut vtab: Vec<VEntry> = Vec::new();
    for &(phi, alpha, _n, _d) in &info {
        if alpha >= pi - 1e-12 { continue; }
        for raw_t in [phi - alpha, phi + alpha] {
            let t = canon(raw_t.rem_euclid(two_pi));
            if !vtab.iter().any(|v| (v.theta - t).abs() < 1e-12) {
                let (c, s) = (t.cos(), t.sin());
                let lo = brep.vertices.len();
                brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, cz_lo) });
                let hi = brep.vertices.len();
                brep.vertices.push(Vertex { point: DVec3::new(center.x + r * c, center.y + r * s, cz_hi) });
                vtab.push(VEntry { theta: t, lo_idx: lo, hi_idx: hi });
            }
        }
    }
    vtab.sort_by(|a, b| a.theta.partial_cmp(&b.theta).unwrap());
    if vtab.len() < 2 { return brep; }

    // ---- 2. Full cylinder wall ----
    let cyl_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(center.x, center.y, cz_lo),
            axis: DVec3::Z, radius: r, ref_dir: DVec3::X,
        }));
        si
    };
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z, radius: r,
    });
    let ba = push_edge!(circle_bot, -pi / 2.0 - two_pi, -pi / 2.0, 0, 0);
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z, radius: r,
    });
    let ta = push_edge!(circle_top, -pi / 2.0, two_pi - pi / 2.0, 0, 0);

    let v0_lo = vtab[0].lo_idx;
    let v0_hi = vtab[0].hi_idx;
    let seam_curve = Curve3::Line(Line3 { origin: brep.vertices[v0_lo].point, direction: DVec3::Z });
    let seam_gen = push_edge!(seam_curve, 0.0, h, v0_lo, v0_hi);

    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(ba),     // bottom arc V0_lo閳壐0_lo
            WireEdge::fwd(seam_gen), // up
            WireEdge::fwd(ta),     // top arc V0_lo閳壐0_hi
            WireEdge::rev(seam_gen), // down
        ],
    };
    let fi_cyl = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: cyl_wire, inner_wires: Vec::new(),
        normal: DVec3::ZERO, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
                surface_idx: None,
    });
    while brep.geom.face_surface.len() <= fi_cyl { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_cyl] = Some(cyl_surf_idx);
    brep.geom.face_surface_range.push(Some([0.0, two_pi, 0.0, h]));

    // ---- 3. Side faces on clip planes ----
    for pi_idx in 0..clip_planes.len() {
        let (n, _d) = clip_planes[pi_idx];
        let (_phi, alpha, _n2, _d2) = info[pi_idx];
        if alpha >= pi - 1e-12 { continue; }

        let raw_lo = (info[pi_idx].0 - info[pi_idx].1).rem_euclid(two_pi);
        let raw_hi = (info[pi_idx].0 + info[pi_idx].1).rem_euclid(two_pi);
        let t_lo = canon(raw_lo);
        let t_hi = canon(raw_hi);

        if let (Some(ilo), Some(ihi)) = (
            vtab.iter().position(|v| (v.theta - t_lo).abs() < 1e-12),
            vtab.iter().position(|v| (v.theta - t_hi).abs() < 1e-12),
        ) {
            let (vlo_lo, vlo_hi) = (vtab[ilo].lo_idx, vtab[ilo].hi_idx);
            let (vhi_lo, vhi_hi) = (vtab[ihi].lo_idx, vtab[ihi].hi_idx);

            let chord_bot = Curve3::Line(Line3 {
                origin: brep.vertices[vlo_lo].point,
                direction: brep.vertices[vhi_lo].point - brep.vertices[vlo_lo].point,
            });
            let eb = push_edge!(chord_bot, 0.0, 1.0, vlo_lo, vhi_lo);
            let chord_top = Curve3::Line(Line3 {
                origin: brep.vertices[vlo_hi].point,
                direction: brep.vertices[vhi_hi].point - brep.vertices[vlo_hi].point,
            });
            let et = push_edge!(chord_top, 0.0, 1.0, vlo_hi, vhi_hi);
            let gen_lo_curve = Curve3::Line(Line3 { origin: brep.vertices[vlo_lo].point, direction: DVec3::Z });
            let gen_lo = push_edge!(gen_lo_curve, 0.0, h, vlo_lo, vlo_hi);
            let gen_hi_curve = Curve3::Line(Line3 { origin: brep.vertices[vhi_lo].point, direction: DVec3::Z });
            let gen_hi = push_edge!(gen_hi_curve, 0.0, h, vhi_lo, vhi_hi);

            let side_plane_idx = {
                let si = brep.geom.surfaces.len();
                brep.geom.surfaces.push(Surface3::Plane(Plane {
                    origin: center - n * info[pi_idx].3,
                    normal: -n,
                }));
                si
            };
            let fi = brep.solids[0].shells[0].faces.len();
            brep.solids[0].shells[0].faces.push(Face {
                outer_wire: Wire {
                    edges: vec![
                        WireEdge::fwd(eb), WireEdge::fwd(gen_hi),
                        WireEdge::rev(et), WireEdge::rev(gen_lo),
                    ],
                },
                inner_wires: Vec::new(), normal: -n,
                triangles: Vec::new(), sample_point: None, mesh_dirty: true,
                surface_idx: None,
            });
            while brep.geom.face_surface.len() <= fi { brep.geom.face_surface.push(None); }
            brep.geom.face_surface[fi] = Some(side_plane_idx);
        }
    }

    // ---- 4. Cap faces with inner wire (box polygon) ----
    let inner_v: Vec<usize> = vtab.iter().map(|v| v.lo_idx).collect();
    let n_inner = inner_v.len();
    let mut inner_edges: Vec<WireEdge> = Vec::with_capacity(n_inner);
    for i in 0..n_inner {
        let j = (i + 1) % n_inner;
        let p_a = brep.vertices[inner_v[i]].point;
        let p_b = brep.vertices[inner_v[j]].point;
        let chord = Curve3::Line(Line3 { origin: p_a, direction: p_b - p_a });
        let e = push_edge!(chord, 0.0, 1.0, inner_v[i], inner_v[j]);
        inner_edges.push(WireEdge::fwd(e));
    }

    let bot_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
        }));
        si
    };
    let top_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
        }));
        si
    };

    // Bottom cap: outer=full circle (CCW), inner=box polygon (CW)
    let fi_bot = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: Wire { edges: vec![WireEdge::fwd(ba)] },
        inner_wires: vec![Wire { edges: inner_edges.clone() }],
        normal: -DVec3::Z, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
                surface_idx: None,
    });
    while brep.geom.face_surface.len() <= fi_bot { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_bot] = Some(bot_surf_idx);

    // Top cap: same
    let fi_top = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: Wire { edges: vec![WireEdge::fwd(ta)] },
        inner_wires: vec![Wire { edges: inner_edges }],
        normal: DVec3::Z, triangles: Vec::new(),
        sample_point: None, mesh_dirty: true,
                surface_idx: None,
    });
    while brep.geom.face_surface.len() <= fi_top { brep.geom.face_surface.push(None); }
    brep.geom.face_surface[fi_top] = Some(top_surf_idx);

    brep
}

/// Build the middle Z-slice of C \ B for the parallel-only clip plane case.
///
/// When all clip planes have parallel normals (e.g. cylinder center inside box
/// in one XY axis), the gap routing in `build_cylinder_box_clipped_brep` fails
/// because `build_plane_chain` cannot construct valid corners between parallel
/// planes.  Instead, we build each outside-slab arc independently and compound.
fn build_cylinder_box_difference_parallel_only(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    build_cylinder_box_difference_parallel_only_skip(center, r, h, clip_planes, false, false)
}

/// Same as [`build_cylinder_box_difference_parallel_only`] but with optional cap skipping.
fn build_cylinder_box_difference_parallel_only_skip(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    // Collect arcs and merge into a single non-compound BRep (avoids nested
    // compound when the caller compounds this with top/bottom pieces).
    //
    // Each clip_plane (n, d) has n pointing INTO the box interior.  The
    // outside-slab (difference) arc on the opposite side uses clip direction
    // = -n (from center toward the clip face) with the same cut_dist d.
    let mut arcs: Vec<BRep> = Vec::new();
    for &(n, d) in clip_planes {
        let clip_dir = -n; // from center toward the clip plane (outward from box)
        if d >= r - 1e-12 {
            continue;
        }
        arcs.push(build_cylinder_arc_for_difference_skip(
            center, r, h, clip_dir, d, skip_bottom_cap, skip_top_cap,
        ));
    }
    if arcs.is_empty() {
        return BRep::default();
    }
    if arcs.len() == 1 {
        let arc = arcs.into_iter().next().unwrap();
        return arc;
    }
    let mut merged = BRep::new();
    for arc in &arcs {
        merged.append_disjoint_brep(arc);
    }
    merged
}

pub(crate) fn build_cylinder_box_clipped_brep(
    center: DVec3,
    r: f64,
    h: f64,
    intervals: &[(f64, f64)],
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
    use_chain_routing: bool,
) -> BRep {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let mut brep = BRep::new();
    brep.solids.push(Solid {
        shells: vec![Shell { faces: Vec::new() }],
    });

    // ---- 1. use provided intervals ----
    let mut intervals = intervals.to_vec();
    intervals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    if intervals.is_empty() {
        return brep;
    }

    // Pre-compute (phi, alpha, normal, cut_dist) for each plane
    let info: Vec<(f64, f64, DVec3, f64)> = clip_planes.iter().map(|&(n, d)| {
        let alpha = (-d / r).clamp(-1.0, 1.0).acos();
        let phi = n.y.atan2(n.x);
        (phi, alpha, n, d)
    }).collect();

    // ---- 2. vertices 閳?one pair (lo, hi) per unique 鑳?----
    // Canonicalize: 鑳?near 0 or 2锜?閳?0
    let canon = |t: f64| if t >= two_pi - 1e-12 { 0.0 } else { t };

    fn push_vertex(brep: &mut BRep, p: DVec3) -> usize {
        let idx = brep.vertices.len();
        brep.vertices.push(Vertex { point: p });
        idx
    }

    struct VEntry { theta: f64, lo: usize, hi: usize }
    let mut vtab: Vec<VEntry> = Vec::new();
    let mut interval_verts: Vec<(usize, usize)> = Vec::new(); // (s_entry, e_entry)

    for &(s_raw, e_raw) in &intervals {
        let s = canon(s_raw);
        let e = canon(e_raw);
        for t in [s, e] {
            if !vtab.iter().any(|ve| (ve.theta - t).abs() < 1e-12) {
                let (c, sn) = (t.cos(), t.sin());
                let lo = push_vertex(&mut brep, DVec3::new(center.x + r * c, center.y + r * sn, cz_lo));
                let hi = push_vertex(&mut brep, DVec3::new(center.x + r * c, center.y + r * sn, cz_hi));
                vtab.push(VEntry { theta: t, lo, hi });
            }
        }
        let sidx = vtab.iter().position(|ve| (ve.theta - s).abs() < 1e-12).unwrap();
        let eidx = vtab.iter().position(|ve| (ve.theta - e).abs() < 1e-12).unwrap();
        interval_verts.push((sidx, eidx));
    }

    // Quick helpers
    let v_lo = |entry: usize| vtab[entry].lo;
    let v_hi = |entry: usize| vtab[entry].hi;

    // ---- 3. edge helper (same pattern as build_half_cylinder_intersection_brep) ----
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
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

    // ---- 4. cylindrical wall faces (one per interval) ----
    // We also collect generator edges & chord edges needed for cap / side faces.
    let n_int = intervals.len();

    // Store per-interval edges: (bottom_arc, right_gen, top_arc, left_gen)
    struct IntervalEdges { ba: usize, rg: usize, ta: usize, lg: usize }
    let mut interval_edges: Vec<IntervalEdges> = Vec::with_capacity(n_int);

    for (i, &(s_raw, e_raw)) in intervals.iter().enumerate() {
        let s = canon(s_raw);
        let e = canon(e_raw);
        let si = interval_verts[i].0;
        let ei = interval_verts[i].1;

        // Bottom circle (normal = 閳妬)
        // Circle3(normal=-Z): C(t) = center + r*(-sin(t), -cos(t), 0)
        // For vertex at standard CCW angle 鑳? P = center + r*(cos(鑳?, sin(鑳?, 0)
        // Mapping: (-sin(t), -cos(t)) = (cos(鑳?, sin(鑳?) 閳?t = -锜?2 - 鑳?
        let circle_bot = Curve3::Circle(Circle3 {
            center: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
            radius: r,
        });
        // Stored as V_e 閳?V_s, used as rev in wires 閳?effective V_s 閳?V_e (CCW)
        let ba = next_curve(circle_bot, -pi / 2.0 - e_raw, -pi / 2.0 - s_raw, v_lo(ei), v_lo(si));

        // Right generator at 鑳?= e
        let p_e_lo = DVec3::new(center.x + r * e.cos(), center.y + r * e.sin(), cz_lo);
        let rg = next_curve(
            Curve3::Line(Line3 { origin: p_e_lo, direction: DVec3::Z }),
            0.0, h, v_lo(ei), v_hi(ei),
        );

        // Top circle (normal = +Z)
        // Circle3(normal=+Z): C(t) = center + r*(-sin(t), cos(t), 0)
        // Mapping: (-sin(t), cos(t)) = (cos(鑳?, sin(鑳?) 閳?t = 鑳?- 锜?2
        let circle_top = Curve3::Circle(Circle3 {
            center: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
            radius: r,
        });
        // Stored as V_s 閳?V_e (fwd in cap wire = 鑳僟s閳崹绔塭 CCW; rev in cyl wall = 鑳僟e閳崹绔塻)
        let ta = next_curve(circle_top, s_raw - pi / 2.0, e_raw - pi / 2.0, v_hi(si), v_hi(ei));

        // Left generator at 鑳?= s
        let p_s_lo = DVec3::new(center.x + r * s.cos(), center.y + r * s.sin(), cz_lo);
        let lg = next_curve(
            Curve3::Line(Line3 { origin: p_s_lo, direction: DVec3::Z }),
            0.0, h, v_lo(si), v_hi(si),
        );

        // Cylindrical wall face
        let cyl_wire = Wire {
            edges: vec![
                WireEdge::rev(ba),  // V_s_lo 閳?V_e_lo
                WireEdge::fwd(rg),  // V_e_lo 閳?V_e_hi
                WireEdge::rev(ta),  // V_e_hi 閳?V_s_hi
                WireEdge::rev(lg),  // V_s_hi 閳?V_s_lo
            ],
        };

        let si_cyl = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(center.x, center.y, cz_lo),
            axis: DVec3::Z,
            radius: r,
            ref_dir: DVec3::X,
        }));

        let mid_theta = 0.5 * (s_raw + e_raw);
        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: cyl_wire,
            inner_wires: Vec::new(),
            normal: DVec3::new(mid_theta.cos(), mid_theta.sin(), 0.0),
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
                surface_idx: None,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(si_cyl);
        while brep.geom.face_surface_range.len() <= fi {
            brep.geom.face_surface_range.push(None);
        }
        // UV range: u 閳?[鑳僟s, 鑳僟e], v 閳?[0, h]
        brep.geom.face_surface_range[fi] = Some([s_raw, e_raw, 0.0, h]);

        interval_edges.push(IntervalEdges { ba, rg, ta, lg });
    }

    // ---- 5. chord edges between interval endpoints (on clip planes) ----
    // Each gap from interval[i].e 閳?interval[i+1].s (mod n_int) may have 1+
    // segments (one per clip plane traversed).  We record gaps per interval as
    // a list of segments, each with bottom/top chord edges and the generator
    // edges at their endpoints (needed later for side faces).

    struct GapSeg {
        bot_chord: usize,
        top_chord: usize,
        plane: usize,
        gen_from: usize, // V_from_lo 閳?V_from_hi (stored direction)
        gen_to: usize,   // V_to_lo 閳?V_to_hi (stored direction)
    }

    let mut gap_segs: Vec<Vec<GapSeg>> = Vec::with_capacity(n_int);

    for gi in 0..n_int {
        let i0 = gi;
        let i1 = if gi + 1 < n_int { gi + 1 } else { 0 };
        let theta_from = intervals[i0].1; // e of interval i0
        let theta_to = intervals[i1].0;   // s of interval i1

        if (canon(theta_from) - canon(theta_to)).abs() < 1e-12 {
            gap_segs.push(Vec::new());
            continue;
        }

        let p_from = find_plane_for_theta(theta_from, &info, r);
        let p_to = find_plane_for_theta(theta_to, &info, r);

        match (p_from, p_to) {
            (Some(pidx), Some(pidx2)) if pidx == pidx2 => {
                // Same plane 閳?single chord
                let vi_fr = interval_verts[i0].1;
                let vi_to = interval_verts[i1].0;

                let chord_bot = Curve3::Line(Line3 {
                    origin: brep.vertices[v_lo(vi_fr)].point,
                    direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[v_lo(vi_fr)].point,
                });
                let eb = next_curve(chord_bot, 0.0, 1.0, v_lo(vi_fr), v_lo(vi_to));

                let chord_top = Curve3::Line(Line3 {
                    origin: brep.vertices[v_hi(vi_fr)].point,
                    direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[v_hi(vi_fr)].point,
                });
                let et = next_curve(chord_top, 0.0, 1.0, v_hi(vi_fr), v_hi(vi_to));

                gap_segs.push(vec![GapSeg {
                    bot_chord: eb,
                    top_chord: et,
                    plane: pidx,
                    gen_from: interval_edges[i0].rg,
                    gen_to: interval_edges[i1].lg,
                }]);
            }
            (Some(p1), Some(p2)) if p1 != p2 => {
                let vi_fr = interval_verts[i0].1;
                let vi_to = interval_verts[i1].0;

                if use_chain_routing {
                    // Chain routing: route through intermediate clip planes when
                    // the direct corner is degenerate (parallel planes).  Used for
                    // the intersection case where corners are inside the cylinder.
                    let chain = build_plane_chain(p1, p2, clip_planes, center);
                    let n_segs = chain.len();

                    // Step 1: Create corner vertices and generator edges.
                    let n_corners = n_segs.saturating_sub(1);
                    struct CornerData { lo: usize, hi: usize, gen_edge: usize }
                    let mut corner_data: Vec<Option<CornerData>> = Vec::new();
                    corner_data.resize_with(n_corners, || None);
                    for j in 0..n_corners {
                        let (pa, pb) = (chain[j], chain[j + 1]);
                        let (na, da) = (clip_planes[pa].0, clip_planes[pa].1);
                        let (nb, db) = (clip_planes[pb].0, clip_planes[pb].1);
                        let corner_xy = corner_of_planes(na, da, nb, db, center);
                        let lo = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_lo) });
                        let hi = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_hi) });
                        let gen_edge = next_curve(
                            Curve3::Line(Line3 { origin: brep.vertices[lo].point, direction: DVec3::Z }),
                            0.0, h, lo, hi,
                        );
                        corner_data[j] = Some(CornerData { lo, hi, gen_edge });
                    }

                    // Step 2: Build one segment per plane in the chain.
                    let mut segs: Vec<GapSeg> = Vec::with_capacity(n_segs);
                    for j in 0..n_segs {
                        let plane = chain[j];
                        let is_first = j == 0;
                        let is_last = j == n_segs - 1;

                        let (start_lo, start_hi) = if is_first {
                            (v_lo(vi_fr), v_hi(vi_fr))
                        } else {
                            let c = corner_data[j - 1].as_ref().unwrap();
                            (c.lo, c.hi)
                        };
                        let (end_lo, end_hi) = if is_last {
                            (v_lo(vi_to), v_hi(vi_to))
                        } else {
                            let c = corner_data[j].as_ref().unwrap();
                            (c.lo, c.hi)
                        };
                        let gen_from = if is_first {
                            interval_edges[i0].rg
                        } else {
                            corner_data[j - 1].as_ref().unwrap().gen_edge
                        };
                        let gen_to = if is_last {
                            interval_edges[i1].lg
                        } else {
                            corner_data[j].as_ref().unwrap().gen_edge
                        };

                        let bot_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[start_lo].point,
                            direction: brep.vertices[end_lo].point - brep.vertices[start_lo].point,
                        });
                        let eb = next_curve(bot_chord, 0.0, 1.0, start_lo, end_lo);

                        let top_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[start_hi].point,
                            direction: brep.vertices[end_hi].point - brep.vertices[start_hi].point,
                        });
                        let et = next_curve(top_chord, 0.0, 1.0, start_hi, end_hi);

                        segs.push(GapSeg {
                            bot_chord: eb,
                            top_chord: et,
                            plane,
                            gen_from,
                            gen_to,
                        });
                    }
                    gap_segs.push(segs);
                } else {
                    // Direct corner routing for difference case: connect gap
                    // endpoints through the direct intersection corner of the
                    // two bounding clip planes.  If the corner is outside the
                    // cylinder or the planes are parallel, use a single chord
                    // between cylinder wall points to avoid "fin" faces that
                    // extend outside the cylinder cross-section.
                    let (n1, d1) = (clip_planes[p1].0, clip_planes[p1].1);
                    let (n2, d2) = (clip_planes[p2].0, clip_planes[p2].1);
                    let corner_xy = corner_of_planes(n1, d1, n2, d2, center);
                    let corner_dist_sq = (corner_xy.x - center.x).powi(2)
                        + (corner_xy.y - center.y).powi(2);

                    // Only use 2-segment corner routing when the corner is
                    // inside the cylinder AND the planes are non-parallel.
                    let non_parallel = (n1.x * n2.y - n1.y * n2.x).abs() > 1e-12;
                    if non_parallel && corner_dist_sq <= r * r + 1e-9 {
                        // 2 segments through the corner
                        let lo_corner = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_lo) });
                        let hi_corner = brep.vertices.len();
                        brep.vertices.push(Vertex { point: DVec3::new(corner_xy.x, corner_xy.y, cz_hi) });
                        let gen_corner = next_curve(
                            Curve3::Line(Line3 { origin: brep.vertices[lo_corner].point, direction: DVec3::Z }),
                            0.0, h, lo_corner, hi_corner,
                        );

                        // Segment 1: p1 閳?corner
                        let bot_chord1 = Curve3::Line(Line3 {
                            origin: brep.vertices[v_lo(vi_fr)].point,
                            direction: brep.vertices[lo_corner].point - brep.vertices[v_lo(vi_fr)].point,
                        });
                        let eb1 = next_curve(bot_chord1, 0.0, 1.0, v_lo(vi_fr), lo_corner);
                        let top_chord1 = Curve3::Line(Line3 {
                            origin: brep.vertices[v_hi(vi_fr)].point,
                            direction: brep.vertices[hi_corner].point - brep.vertices[v_hi(vi_fr)].point,
                        });
                        let et1 = next_curve(top_chord1, 0.0, 1.0, v_hi(vi_fr), hi_corner);

                        // Segment 2: corner 閳?p2
                        let bot_chord2 = Curve3::Line(Line3 {
                            origin: brep.vertices[lo_corner].point,
                            direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[lo_corner].point,
                        });
                        let eb2 = next_curve(bot_chord2, 0.0, 1.0, lo_corner, v_lo(vi_to));
                        let top_chord2 = Curve3::Line(Line3 {
                            origin: brep.vertices[hi_corner].point,
                            direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[hi_corner].point,
                        });
                        let et2 = next_curve(top_chord2, 0.0, 1.0, hi_corner, v_hi(vi_to));

                        gap_segs.push(vec![
                            GapSeg {
                                bot_chord: eb1, top_chord: et1, plane: p1,
                                gen_from: interval_edges[i0].rg,
                                gen_to: gen_corner,
                            },
                            GapSeg {
                                bot_chord: eb2, top_chord: et2, plane: p2,
                                gen_from: gen_corner,
                                gen_to: interval_edges[i1].lg,
                            },
                        ]);
                    } else {
                        // Single direct chord between cylinder wall points
                        // (avoids fins for corners outside the cylinder, and
                        // handles parallel planes).
                        let bot_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[v_lo(vi_fr)].point,
                            direction: brep.vertices[v_lo(vi_to)].point - brep.vertices[v_lo(vi_fr)].point,
                        });
                        let eb = next_curve(bot_chord, 0.0, 1.0, v_lo(vi_fr), v_lo(vi_to));
                        let top_chord = Curve3::Line(Line3 {
                            origin: brep.vertices[v_hi(vi_fr)].point,
                            direction: brep.vertices[v_hi(vi_to)].point - brep.vertices[v_hi(vi_fr)].point,
                        });
                        let et = next_curve(top_chord, 0.0, 1.0, v_hi(vi_fr), v_hi(vi_to));
                        gap_segs.push(vec![GapSeg {
                            bot_chord: eb, top_chord: et, plane: p1,
                            gen_from: interval_edges[i0].rg,
                            gen_to: interval_edges[i1].lg,
                        }]);
                    }
                }
            }
            _ => {
                gap_segs.push(Vec::new());
            }
        }
    }

    // ---- 6. bottom & top planar cap faces ----
    let bot_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_lo),
            normal: -DVec3::Z,
        }));
        si
    };
    let top_surf_idx = {
        let si = brep.geom.surfaces.len();
        brep.geom.surfaces.push(Surface3::Plane(Plane {
            origin: DVec3::new(center.x, center.y, cz_hi),
            normal: DVec3::Z,
        }));
        si
    };

    // Bottom cap: V_s_lo 閳?V_e_lo (arc) 閳?V_s_{i+1}_lo (chords) 閳?...
    let mut bot_wire_edges: Vec<WireEdge> = Vec::new();
    let mut top_wire_edges: Vec<WireEdge> = Vec::new();
    for i in 0..n_int {
        // Arc from start-vertex to end-vertex (bottom: V_s_lo閳壐_e_lo, top: V_s_hi閳壐_e_hi)
        bot_wire_edges.push(WireEdge::rev(interval_edges[i].ba));
        top_wire_edges.push(WireEdge::fwd(interval_edges[i].ta));
        // Gap chord segments from end-i to start-(i+1)
        for seg in &gap_segs[i] {
            bot_wire_edges.push(WireEdge::fwd(seg.bot_chord));
            top_wire_edges.push(WireEdge::fwd(seg.top_chord));
        }
    }

    let push_face = |brep: &mut BRep, wire: Wire, surf_idx: usize, normal: DVec3| {
        let fi = brep.solids[0].shells[0].faces.len();
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: wire,
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
    };
    if !skip_bottom_cap {
        push_face(&mut brep, Wire { edges: bot_wire_edges }, bot_surf_idx, -DVec3::Z);
    }
    if !skip_top_cap {
        push_face(&mut brep, Wire { edges: top_wire_edges }, top_surf_idx, DVec3::Z);
    }

    // ---- 7. side faces (one per gap segment, on the segment's clip plane) ----
    // The face wire runs: bot_chord_fwd 閳?gen_to_fwd 閳?top_chord_rev 閳?gen_from_rev
    for segs in &gap_segs {
        for seg in segs {
            // Determine the plane for this side face.
            // Most gap segments lie on a clip-plane, but the "single direct chord"
            // fallback for parallel clip planes (bopcut_simple t6/v4) creates a
            // chord on the cylinder wall that does NOT lie on the assigned plane.
            // In that case we infer the correct vertical plane from the chord.
            let (n, d) = {
                let (n0, d0) = clip_planes[seg.plane];
                let bot_e = &brep.edges[seg.bot_chord];
                let s_pt = brep.vertices[bot_e.start].point;
                let e_pt = brep.vertices[bot_e.end].point;
                // Check whether both endpoints lie on the assigned clip plane.
                let tol = 1e-8;
                let on_p0 = (n0.dot(s_pt - center) - (-d0)).abs() < tol
                         && (n0.dot(e_pt - center) - (-d0)).abs() < tol;
                if on_p0 {
                    (n0, d0)
                } else {
                    // Chord endpoints don't lie on the assigned clip plane.
                    // Compute a vertical plane through the chord: the face
                    // winding normal is cross(chord_dir, Z), and we need -n
                    // to match it (since push_face stores outward = -n).
                    let chord_dir = e_pt - s_pt;
                    let cross_z = chord_dir.cross(DVec3::Z);
                    let len = cross_z.length();
                    if len < 1e-15 {
                        // Degenerate chord 閳?fall back to original plane.
                        (n0, d0)
                    } else {
                        let n = -cross_z / len; // so that -n = winding normal
                        let d = -n.dot(s_pt - center);
                        (n, d)
                    }
                }
            };
            let si_side = {
                let si = brep.geom.surfaces.len();
                brep.geom.surfaces.push(Surface3::Plane(Plane {
                    origin: center - n * d,
                    normal: -n,
                }));
                si
            };
            let side_wire = Wire {
                edges: vec![
                    WireEdge::fwd(seg.bot_chord),
                    WireEdge::fwd(seg.gen_to),
                    WireEdge::rev(seg.top_chord),
                    WireEdge::rev(seg.gen_from),
                ],
            };
            push_face(&mut brep, side_wire, si_side, -n);
        }
    }

    brep
}
/// Build an arc BRep for the parallel-only cylinder-box difference case.
///
/// This builds the portion of a Z-aligned cylinder satisfying
/// `(P - center)璺痗lip_n 閳?cut_dist` (the outside-slab region). The arc is
/// centered on `clip_n` with half-angle `浼?= acos(cut_dist/r)`. The clip face
/// is on the plane `center + clip_n璺痗ut_dist` with outward normal `-clip_n`.
///
/// `clip_n`: horizontal unit normal from center toward the clip plane.
/// `cut_dist`: distance from center to clip plane (閳?, 閳槝).
fn build_cylinder_arc_for_difference(
    center: DVec3,
    r: f64,
    h: f64,
    clip_n: DVec3,
    cut_dist: f64,
) -> BRep {
    build_cylinder_arc_for_difference_skip(center, r, h, clip_n, cut_dist, false, false)
}

/// Same as [`build_cylinder_arc_for_difference`] but with optional cap skipping.
fn build_cylinder_arc_for_difference_skip(
    center: DVec3,
    r: f64,
    h: f64,
    clip_n: DVec3,
    cut_dist: f64,
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    let alpha = (cut_dist / r).clamp(-1.0, 1.0).acos();
    let phi = clip_n.y.atan2(clip_n.x);
    let half_h = h * 0.5;
    let cz_lo = center.z - half_h;
    let cz_hi = center.z + half_h;

    let (sa, ca) = alpha.sin_cos();
    let (sp, cp) = phi.sin_cos();

    // (cos(锠佸崵浼?, sin(锠佸崵浼?)
    let cos_phi_minus_alpha = cp * ca + sp * sa;
    let sin_phi_minus_alpha = sp * ca - cp * sa;
    let cos_phi_plus_alpha = cp * ca - sp * sa;
    let sin_phi_plus_alpha = sp * ca + cp * sa;

    let v0_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_lo);
    let v1_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_lo);
    let v2_p = DVec3::new(center.x + r * cos_phi_plus_alpha, center.y + r * sin_phi_plus_alpha, cz_hi);
    let v3_p = DVec3::new(center.x + r * cos_phi_minus_alpha, center.y + r * sin_phi_minus_alpha, cz_hi);

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

    // Edge helper
    let mut next_curve = |c: Curve3, t0: f64, t1: f64, start: usize, end: usize| -> usize {
        let idx = brep.edges.len();
        brep.edges.push(Edge { start, end });
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(c);
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

    // E0: bottom arc (V1閳壐0), same convention as build_half_cylinder_intersection_brep
    let circle_bot = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
        radius: r,
    });
    let e0 = next_curve(circle_bot, -phi - alpha, -phi + alpha, v1, v0);

    // E1: right generator (V1閳壐2)
    let line_r = Curve3::Line(Line3 { origin: v1_p, direction: DVec3::Z });
    let e1 = next_curve(line_r, 0.0, h, v1, v2);

    // E2: top arc (V3閳壐2)
    let circle_top = Curve3::Circle(Circle3 {
        center: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
        radius: r,
    });
    let e2 = next_curve(circle_top, phi - alpha, phi + alpha, v3, v2);

    // E3: left generator (V0閳壐3)
    let line_l = Curve3::Line(Line3 { origin: v0_p, direction: DVec3::Z });
    let e3 = next_curve(line_l, 0.0, h, v0, v3);

    // E4: bottom chord on clip plane (V0閳壐1)
    let line_cb = Curve3::Line(Line3 { origin: v0_p, direction: v1_p - v0_p });
    let e4 = next_curve(line_cb, 0.0, 1.0, v0, v1);

    // E5: top chord on clip plane (V2閳壐3)
    let line_ct = Curve3::Line(Line3 { origin: v2_p, direction: v3_p - v2_p });
    let e5 = next_curve(line_ct, 0.0, 1.0, v2, v3);

    // --- Surfaces ---
    let cyl_surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(center.x, center.y, cz_lo),
        axis: DVec3::Z, radius: r, ref_dir: DVec3::X,
    });
    let top_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_hi),
        normal: DVec3::Z,
    });
    let bot_plane = Surface3::Plane(Plane {
        origin: DVec3::new(center.x, center.y, cz_lo),
        normal: -DVec3::Z,
    });
    // Clip plane at center + clip_n*cut_dist, outward normal = -clip_n
    let clip_surf = Surface3::Plane(Plane {
        origin: center + clip_n * cut_dist,
        normal: -clip_n,
    });

    let si_cyl = brep.geom.surfaces.len();
    brep.geom.surfaces.push(cyl_surf);
    let si_top = brep.geom.surfaces.len();
    brep.geom.surfaces.push(top_plane);
    let si_bot = brep.geom.surfaces.len();
    brep.geom.surfaces.push(bot_plane);
    let si_clip = brep.geom.surfaces.len();
    brep.geom.surfaces.push(clip_surf);

    // Face helper
    let mut push_face = |outer: Wire, surf_idx: usize, normal: DVec3| -> usize {
        let fi = if brep.solids.is_empty() {
            brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });
            0
        } else {
            brep.solids[0].shells[0].faces.len()
        };
        brep.solids[0].shells[0].faces.push(Face {
            outer_wire: outer, inner_wires: Vec::new(), normal,
            triangles: Vec::new(), sample_point: None, mesh_dirty: true,
                surface_idx: None,
        });
        while brep.geom.face_surface.len() <= fi {
            brep.geom.face_surface.push(None);
        }
        brep.geom.face_surface[fi] = Some(surf_idx);
        fi
    };

    // F0: Cylindrical wall 閳?wire V0閳壐1閳壐2閳壐3閳壐0
    let cyl_wire = Wire {
        edges: vec![
            WireEdge::rev(e0), WireEdge::fwd(e1),
            WireEdge::rev(e2), WireEdge::rev(e3),
        ],
    };
    let _f0 = push_face(cyl_wire, si_cyl, clip_n);
    while brep.geom.face_surface_range.len() <= _f0 {
        brep.geom.face_surface_range.push(None);
    }
    brep.geom.face_surface_range[_f0] = Some([phi - alpha, phi + alpha, 0.0, h]);

    // F1: Top cap (normal=+Z)
    if !skip_top_cap {
        let _f1 = push_face(Wire { edges: vec![WireEdge::fwd(e5), WireEdge::fwd(e2)] }, si_top, DVec3::Z);
    }

    // F2: Bottom cap (normal=-Z)
    if !skip_bottom_cap {
        let _f2 = push_face(Wire { edges: vec![WireEdge::fwd(e4), WireEdge::fwd(e0)] }, si_bot, -DVec3::Z);
    }

    // F3: Clip face (bounding the arc on the box face, outward normal = -clip_n)
    let clip_wire = Wire {
        edges: vec![
            WireEdge::fwd(e4), WireEdge::fwd(e1),
            WireEdge::fwd(e5), WireEdge::rev(e3),
        ],
    };
    let _f3 = push_face(clip_wire, si_clip, -clip_n);

    brep
}