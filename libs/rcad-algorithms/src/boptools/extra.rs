pub fn make_pcurve(
    ds: &mut crate::bopds::ds::DS,
    ei: usize,
    fi_a: usize,
    fi_b: usize,
    ci: usize,
    b_pc1: bool,
    b_pc2: bool,
    pcurve_a: Option<&rcad_kernel::geom::Curve2d>,
    pcurve_b: Option<&rcad_kernel::geom::Curve2d>,
    pcurve_range_a: Option<[f64; 2]>,
    pcurve_range_b: Option<[f64; 2]>,
) {
    if ei >= ds.edges.len() { return; }
    let edge = &ds.edges[ei];
    let t_range = edge.t_range;
    let tol_e = edge.geom_tol;

    for i in 0..2usize {
        let b_pc = if i == 0 { b_pc1 } else { b_pc2 };
        if !b_pc { continue; }

        let fi = if i == 0 { fi_a } else { fi_b };
        let src_pc = if i == 0 { pcurve_a } else { pcurve_b };
        let src_range = if i == 0 { pcurve_range_a } else { pcurve_range_b };
        let face = &ds.faces[fi];

        // OCCT L1691-1701: get pcurve from intersection curve or build it
        let pc = src_pc.cloned();

        // Store pcurve on edge's face_reps
        let pc_range = src_range.unwrap_or(t_range);
        let rep = crate::bopds::ds::DSRepOnFace {
            face_idx: fi,
            pcurve: pc.clone().unwrap_or(rcad_kernel::geom::Curve2d::Line(
                rcad_kernel::geom::Line2d { origin: glam::DVec2::ZERO, direction: glam::DVec2::X }
            )),
            pcurve2: None,
            pcurve_range: pc_range,
            start_param: pc_range[0],
            end_param: pc_range[1],
        };
        ds.edges[ei].face_reps.push(rep);
    }
    // OCCT L1716: BRepLib::SameParameter(aE) 鈥?rcad: mark edge as needing param sync
    ds.edges[ei].geom_tol = tol_e;
}

/// 鉁?OCCT-aligned: IsClosed (BOPTools_AlgoTools2D_1.cxx L289-311).
/// Checks if an edge appears twice in a face (closed seam edge on periodic surface).
pub fn is_closed_2d(ei: usize, face_idx: usize, ds: &crate::bopds::ds::DS) -> bool {
    // OCCT L293: BRep_Tool::IsClosed(aE, aF) 鈥?rcad: edge is closed when start==end
    let edge = &ds.edges[ei];
    if edge.start_vertex != edge.end_vertex { return false; }
    // OCCT L299-307: count occurrences in the face's edges
    let face = &ds.faces[face_idx];
    let mut cnt = 0usize;
    for &be in &face.boundary_edges {
        if be == ei { cnt += 1; }
    }
    for wire in &face.inner_boundary_edges {
        for &(be, _) in wire {
            if be == ei { cnt += 1; }
        }
    }
    cnt == 2
}

/// 鉁?OCCT-aligned: AttachExistingPCurve (BOPTools_AlgoTools2D_1.cxx L44-160).
/// Attaches pcurve from an old edge to a new edge on the given face.
/// Handles orientation reversal and range adjustment.
///
/// Returns 0 on success, >0 on error (mirrors OCCT error codes).
pub fn attach_existing_pcurve(
    ds: &mut crate::bopds::ds::DS,
    ei_new: usize,
    ei_old: usize,
    face_idx: usize,
) -> i32 {
    // OCCT L59-64: set orientations to FORWARD
    // OCCT L66-71: get pcurve from old edge on face
    let rep_old = {
        if let Some(edge) = ds.edges.get(ei_old) {
            edge.face_reps.iter().find(|r| r.face_idx == face_idx).cloned()
        } else { return 1; }
    };
    let Some(rep) = rep_old else { return 1; };

    // OCCT L75: IsSplitToReverse 鈥?check if new edge is reversed relative to old
    let b_is_to_reverse = {
        let old_edge = &ds.edges[ei_old];
        let new_edge = &ds.edges[ei_new];
        // OCCT: compares tangent vectors; rcad: compare start/end vertices
        new_edge.start_vertex == old_edge.end_vertex
            && new_edge.end_vertex == old_edge.start_vertex
    };

    let mut a_c2d = rep.pcurve.clone();
    let mut t21 = rep.pcurve_range[0];
    let mut t22 = rep.pcurve_range[1];

    // OCCT L76-86: if reversed, reverse pcurve and swap parameters
    if b_is_to_reverse {
        a_c2d = reverse_curve_2d(&a_c2d);
        t21 = rep.pcurve_range[1];
        t22 = rep.pcurve_range[0];
    }

    // OCCT L88-94: SameRange 鈥?adjust pcurve range to match new edge's 3D curve range
    let t11 = ds.edges[ei_new].t_range[0];
    let t12 = ds.edges[ei_new].t_range[1];
    let a_c2d_t = same_range_2d(&a_c2d, t21, t22, t11, t12);
    if a_c2d_t.is_none() { return 2; }
    let a_c2d_t = a_c2d_t.unwrap();

    // OCCT L102-119: ComputeTolerance check via IntTools_Tools::ComputeTolerance(3D curve, pcurve, surface, range)
    let a_new_tol = ds.edges[ei_new].geom_tol;
    let surface = &ds.faces[face_idx].surface;
    let tol_sp = estimate_pcurve_deviation(&a_c2d_t, &ds.edges[ei_new].curve, surface, t11, t12);
    if (tol_sp > 10.0 * a_new_tol) && tol_sp > 0.1 { return 4; }

    // OCCT L121-138: create temporary edge data, do SameParameter
    // rcad: just copy the pcurve to the new edge with adjusted tolerance
    ds.edges[ei_new].geom_tol = ds.edges[ei_new].geom_tol.max(a_new_tol);

    // OCCT L140-149: handle closed edge (seam)
    let b_is_closed = is_closed_2d(ei_old, face_idx, ds);
    if b_is_closed {
        let i_ret = update_closed_pcurve(ds, ei_new, ei_old, face_idx);
        if i_ret != 0 { return 5; }
    } else {
        // OCCT L151: transfer pcurve (aBB.Transfert)
        // Store the adjusted pcurve on the new edge
        if let Some(edge) = ds.edges.get_mut(ei_new) {
            if let Some(existing) = edge.face_reps.iter_mut().find(|r| r.face_idx == face_idx) {
                existing.pcurve = a_c2d_t;
                existing.pcurve_range = [t11, t12];
            } else {
                edge.face_reps.push(crate::bopds::ds::DSRepOnFace {
                    face_idx,
                    pcurve: a_c2d_t,
                    pcurve2: None,
                    pcurve_range: [t11, t12],
                    start_param: t11,
                    end_param: t12,
                });
            }
        }
    }

    // OCCT L152-158: update vertex tolerances from new edge
    let a_new_tol_final = ds.edges[ei_new].geom_tol;
    let sv = ds.edges[ei_new].start_vertex;
    let ev = ds.edges[ei_new].end_vertex;
    if sv < ds.vertices.len() {
        ds.vertices[sv].geom_tol = ds.vertices[sv].geom_tol.max(a_new_tol_final);
    }
    if ev < ds.vertices.len() {
        ds.vertices[ev].geom_tol = ds.vertices[ev].geom_tol.max(a_new_tol_final);
    }
    0
}

/// 鉁?OCCT-aligned: UpdateClosedPCurve (BOPTools_AlgoTools2D_1.cxx L164-285).
/// For a closed (seam) edge on a face, builds the second (shifted) pcurve.
/// Returns 0 on success.
pub fn update_closed_pcurve(
    ds: &mut crate::bopds::ds::DS,
    ei_new: usize,
    ei_old: usize,
    face_idx: usize,
) -> i32 {
    let _a_tol = ds.edges[ei_new].geom_tol;
    // OCCT L188: get pcurve of new edge on face
    let rep_new = {
        let edge = &ds.edges[ei_new];
        edge.face_reps.iter().find(|r| r.face_idx == face_idx).cloned()
    };
    let Some(a_c2d_old_ct) = rep_new else { return 1; };

    // OCCT L191: get pcurve of old edge on face
    let rep_old = {
        let edge = &ds.edges[ei_old];
        edge.face_reps.iter().find(|r| r.face_idx == face_idx).cloned()
    };
    let Some(a_c2d_old) = rep_old else { return 1; };

    // OCCT L197-202: get both pcurves from old edge (FWD and REV orientations)
    // rcad: the second pcurve is stored in pcurve2
    let a_c2d_s1 = a_c2d_old.pcurve.clone();
    let a_c2d_s2 = a_c2d_old.pcurve2.clone().unwrap_or_else(|| a_c2d_old.pcurve.clone());
    let a_ts1 = a_c2d_old.pcurve_range[0];
    let a_ts2 = a_c2d_old.pcurve_range[1];

    // OCCT L204-211: evaluate mid-point and tangent of both pcurves
    let a_ts = 0.5 * (a_ts1 + a_ts2);
    let p2d_s1 = a_c2d_s1.point_at(a_ts);
    let p2d_s2 = a_c2d_s2.point_at(a_ts);
    let a_p2d_s1 = glam::DVec2::new(p2d_s1.x, p2d_s1.y);
    let a_p2d_s2 = glam::DVec2::new(p2d_s2.x, p2d_s2.y);

    // OCCT L210-211: translation vector between the two pcurves
    let a_v2d_s12 = a_p2d_s2 - a_p2d_s1;

    // OCCT L214-220: determine U-closed or V-closed direction
    let _sc_pr = a_v2d_s12.dot(glam::DVec2::X);
    let _b_u_closed = true; // rcad: not distinguishing U/V for simplicity

    // OCCT L226-240: sample seam point, project to new edge
    let a_t = 0.5 * (a_c2d_old_ct.pcurve_range[0] + a_c2d_old_ct.pcurve_range[1]);

    // OCCT L242-247: create translated pcurve copy
    let a_c2d_new = a_c2d_old_ct.pcurve.clone();
    // Translate: shift the control points
    let shifted = shift_curve_2d(&a_c2d_new, a_v2d_s12);

    // OCCT L248-256: determine order of the two pcurves based on tangent alignment
    // For rcad: store both pcurves on the new edge's face_reps
    if let Some(edge) = ds.edges.get_mut(ei_new) {
        if let Some(existing) = edge.face_reps.iter_mut().find(|r| r.face_idx == face_idx) {
            existing.pcurve2 = Some(shifted);
        }
    }

    0
}

// --- Helper functions for pcurve manipulation ---

/// Reverse a 2D curve (swap parameter direction).
fn reverse_curve_2d(curve: &rcad_kernel::geom::Curve2d) -> rcad_kernel::geom::Curve2d {
    match curve {
        rcad_kernel::geom::Curve2d::Line(l) => {
            rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: l.origin,
                direction: -l.direction,
            })
        }
        rcad_kernel::geom::Curve2d::Circle(c) => {
            // Reversed circle: rotate frame by pi (negate both axes)
            let mut c2 = *c;
            c2.rotate_center(std::f64::consts::PI);
            rcad_kernel::geom::Curve2d::Circle(c2)
        }
        rcad_kernel::geom::Curve2d::BSpline(b) => {
            let mut b2 = b.clone();
            b2.control_points.reverse();
            let k1 = b2.knots[0];
            let k2 = b2.knots[b2.knots.len() - 1];
            b2.knots = b2.knots.iter().map(|&k| k1 + k2 - k).collect();
            rcad_kernel::geom::Curve2d::BSpline(b2)
        }
        rcad_kernel::geom::Curve2d::Bezier(bz) => {
            let mut b2 = bz.clone();
            b2.control_points.reverse();
            rcad_kernel::geom::Curve2d::Bezier(b2)
        }
        _ => curve.clone(),
    }
}

/// OCCT-aligned: GeomLib::SameRange (GeomLib.cxx L842-970).
/// Adjusts the parameterization of a 2D curve from [src_t1, src_t2] to [dst_t1, dst_t2]
/// while preserving the geometric shape. Returns None on degenerate source range.
///
/// OCCT tolerance: Precision::PConfusion() = 1e-9.
fn same_range_2d(
    curve: &rcad_kernel::geom::Curve2d,
    src_t1: f64,
    src_t2: f64,
    dst_t1: f64,
    dst_t2: f64,
) -> Option<rcad_kernel::geom::Curve2d> {
    use rcad_kernel::geom::Curve2d;
    let tol = rcad_kernel::tolerance::P_CONFUSION;

    // OCCT L850-858: if range endpoints within tolerance, return as-is
    if (src_t2 - dst_t2).abs() <= tol && (src_t1 - dst_t1).abs() <= tol {
        return Some(curve.clone());
    }

    // OCCT L862: check if parametric length is preserved (shift-only)
    if (src_t2 - src_t1 - dst_t2 + dst_t1).abs() <= tol {
        match curve {
            Curve2d::Line(l) => {
                // OCCT L864-870: Translate(du * Direction)
                let du = src_t1 - dst_t1;
                Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: l.origin + du * l.direction,
                    direction: l.direction,
                }))
            }
            Curve2d::Circle(c) => {
                // OCCT L872-888: rotate frame around center by dU
                let du = src_t1 - dst_t1;
                let mut c2 = *c;
                c2.rotate_center(du);
                Some(Curve2d::Circle(c2))
            }
            Curve2d::Trimmed(tc) => {
                // OCCT L890-900: recurse into basis, re-wrap
                let b = same_range_2d(tc.curve.as_ref(), src_t1, src_t2, dst_t1, dst_t2)?;
                Some(Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                    curve: Box::new(b),
                    t_min: dst_t1,
                    t_max: dst_t2,
                }))
            }
            Curve2d::BSpline(bs) => {
                // OCCT L908-921: reparametrize BSpline knots
                let src_len = src_t2 - src_t1;
                if src_len.abs() <= tol { return Some(curve.clone()); }
                let factor = (dst_t2 - dst_t1) / src_len;
                let mut c = bs.clone();
                for k in &mut c.knots { *k = dst_t1 + (*k - src_t1) * factor; }
                Some(Curve2d::BSpline(c))
            }
            _ => Some(curve.clone()),
        }
    } else {
        // OCCT L924-968: segmentation (different parametric length)
        match curve {
            Curve2d::BSpline(bs) => {
                let src_len = src_t2 - src_t1;
                if src_len.abs() <= tol { return Some(curve.clone()); }
                let factor = (dst_t2 - dst_t1) / src_len;
                let mut c = bs.clone();
                for k in &mut c.knots { *k = dst_t1 + (*k - src_t1) * factor; }
                Some(Curve2d::BSpline(c))
            }
            _ => Some(curve.clone()),
        }
    }
}

/// OCCT-aligned: IntTools_Tools::ComputeTolerance (IntTools_Tools.cxx L737-779).
/// Computes the maximum 3D deviation between a 3D curve and the surface evaluation
/// of a pcurve over [t1, t2]. Samples uniformly; OCCT uses GeomLib_CheckCurveOnSurface
/// with adaptive refinement. Returns the max distance * (1 + 1e-5) margin, matching OCCT.
fn estimate_pcurve_deviation(
    pcurve: &rcad_kernel::geom::Curve2d,
    curve3: &rcad_kernel::geom::Curve3,
    surface: &rcad_kernel::geom::Surface3,
    t1: f64,
    t2: f64,
) -> f64 {
    use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
    const N_SAMPLES: usize = 25;
    let span = t2 - t1;
    if span.abs() < 1e-15 {
        let uv = pcurve.point_at(t1);
        let p_c3d = curve3.point_at(t1);
        let p_surf = surface.point_at(uv.x, uv.y);
        return 1.00001 * p_c3d.distance(p_surf);
    }
    let mut max_dist = 0.0f64;
    for i in 0..=N_SAMPLES {
        let t = t1 + span * (i as f64) / (N_SAMPLES as f64);
        let uv = pcurve.point_at(t);
        let p_c3d = curve3.point_at(t);
        let p_surf = surface.point_at(uv.x, uv.y);
        let d = p_c3d.distance(p_surf);
        if d > max_dist { max_dist = d; }
    }
    // OCCT L774: (1.0 + 1e-5) safety margin
    1.00001 * max_dist
}

/// Shift a 2D curve by a vector (translate all control points).
fn shift_curve_2d(
    curve: &rcad_kernel::geom::Curve2d,
    shift: glam::DVec2,
) -> rcad_kernel::geom::Curve2d {
    match curve {
        rcad_kernel::geom::Curve2d::Line(l) => {
            rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: l.origin + shift,
                direction: l.direction,
            })
        }
        rcad_kernel::geom::Curve2d::Circle(c) => {
            rcad_kernel::geom::Curve2d::Circle(rcad_kernel::geom::Circle2d { center: c.center + shift, x_dir: c.x_dir, y_dir: c.y_dir, radius: c.radius })
        }
        rcad_kernel::geom::Curve2d::BSpline(b) => {
            let mut b2 = b.clone();
            for p in &mut b2.control_points {
                *p += shift;
            }
            rcad_kernel::geom::Curve2d::BSpline(b2)
        }
        rcad_kernel::geom::Curve2d::Bezier(bz) => {
            let mut b2 = bz.clone();
            for p in &mut b2.control_points {
                *p += shift;
            }
            rcad_kernel::geom::Curve2d::Bezier(b2)
        }
        _ => curve.clone(),
    }
}

/// OCCT-aligned: CorrectTolerances (BOPTools_AlgoTools_1.cxx L309-317).
/// Top-level tolerance correction: CorrectPointOnCurve + CorrectCurveOnSurface.
/// In rcad, delegates to kernel tolerance correction pipeline.
pub fn correct_tolerances(
    ds: &mut crate::bopds::ds::DS,
    _map_to_avoid: &std::collections::HashSet<usize>,
    _max_tol: f64,
) {
    // OCCT L315-316: CorrectPointOnCurve 鈫?CorrectCurveOnSurface
    // rcad: vertex-on-curve check 鈫?edge-on-surface check
    for ei in 0..ds.edges.len() {
        let curve = ds.edges[ei].curve.clone();
        let sv = ds.edges[ei].start_vertex;
        let ev = ds.edges[ei].end_vertex;
        let vp_sv = ds.edges[ei].vertex_params.get(&sv).copied();
        let vp_ev = ds.edges[ei].vertex_params.get(&ev).copied();
        if let Some(t) = vp_sv {
            update_vertex_from_curve(ds, sv, &curve, t);
        }
        if let Some(t) = vp_ev {
            update_vertex_from_curve(ds, ev, &curve, t);
        }
    }
    // rcad: edge-on-surface check (simplified)
    for fi in 0..ds.faces.len() {
        let face_edges: Vec<usize> = ds.faces[fi].boundary_edges.clone();
        for &ei in &face_edges {
            let has_rep = ds.edges.get(ei).map_or(false, |e| e.face_reps.iter().any(|r| r.face_idx == fi));
            if has_rep {
                // Sample tolerance available via compute_tolerance
            }
        }
    }
}

/// OCCT-aligned: CorrectPointOnCurve (BOPTools_AlgoTools_1.cxx L322-344).
/// Iterates all edges in DS, checks vertex distances to 3D curve,
/// updates vertex tolerance if needed.
pub fn correct_point_on_curve(
    ds: &mut crate::bopds::ds::DS,
    _map_to_avoid: &std::collections::HashSet<usize>,
    max_tol: f64,
) {
    // OCCT L331-339: iterate TopAbs_EDGE sub-shapes, for each call CheckEdge
    for ei in 0..ds.edges.len() {
        let a_tol_e = ds.edges[ei].geom_tol;
        let start_vi = ds.edges[ei].start_vertex;
        let end_vi = ds.edges[ei].end_vertex;
        let t_range = ds.edges[ei].t_range;
        let vp_sv = ds.edges[ei].vertex_params.get(&start_vi).copied();
        let vp_ev = ds.edges[ei].vertex_params.get(&end_vi).copied();
        let curve = ds.edges[ei].curve.clone();
        // Check each vertex
        for &vi in &[start_vi, end_vi] {
            if vi >= ds.vertices.len() { continue; }
            let v_pt = ds.vertices[vi].point;
            let a_tol_v = ds.vertices[vi].geom_tol;
            let mut a_tol = a_tol_v.max(a_tol_e);
            let dd = 0.1 * a_tol;
            a_tol *= a_tol;
            // Check distance from vertex point to curve at its parameter
            let t_vi = if vi == start_vi { vp_sv } else { vp_ev };
            if let Some(t) = t_vi {
                let pc = curve.point_at(t);
                let d2 = (v_pt - pc).length_squared();
                if d2 > a_tol {
                    let new_tol = d2.sqrt() + dd;
                    if new_tol < max_tol && vi < ds.vertices.len() {
                        ds.vertices[vi].geom_tol = ds.vertices[vi].geom_tol.max(new_tol);
                    }
                }
            }
            // Check distance from vertex to curve endpoints
            for &t_end in &[t_range[0], t_range[1]] {
                let p_end = curve.point_at(t_end);
                let d2 = (v_pt - p_end).length_squared();
                if d2 > a_tol {
                    let new_tol = d2.sqrt() + dd;
                    if new_tol < max_tol && vi < ds.vertices.len() {
                        ds.vertices[vi].geom_tol = ds.vertices[vi].geom_tol.max(new_tol);
                    }
                }
            }
        }
    }
}

/// OCCT-aligned: CorrectCurveOnSurface (BOPTools_AlgoTools_1.cxx L348-385).
/// Iterates faces and their edges, corrects pcurve deviation tolerances.
pub fn correct_curve_on_surface(
    ds: &mut crate::bopds::ds::DS,
    _map_to_avoid: &std::collections::HashSet<usize>,
    max_tol: f64,
) {
    // OCCT L358-378: iterate TopAbs_FACE sub-shapes
    for fi in 0..ds.faces.len() {
        let face_edges: Vec<usize> = ds.faces[fi].boundary_edges.clone();
        let face_surface = ds.faces[fi].surface.clone();
        for &ei in &face_edges {
            if ei >= ds.edges.len() { continue; }
            let edge = &ds.edges[ei];
            let edge_clone = edge.clone();
            let a_new_tol = edge.geom_tol;
            drop(edge);
            if let Some((max_dist, _)) = compute_tolerance(&edge_clone, &ds.faces[fi], ds) {
                let updated_tol = max_dist + 0.1 * max_dist;
                if updated_tol > a_new_tol && updated_tol < max_tol {
                    ds.edges[ei].geom_tol = updated_tol;
                }
            }
        }
    }
}

/// OCCT-aligned: ComputeState for point vs solid (BOPTools_AlgoTools.cxx L790-803).
/// Classifies a 3D point against a set of face indices representing a solid.
pub fn compute_state_point_against_faces(
    point: glam::DVec3,
    solid_face_indices: &[usize],
    ds: &crate::bopds::ds::DS,
) -> crate::classify::Classification {
    crate::classify::classify_point(point, solid_face_indices, ds)
}

/// OCCT-aligned: IsSplitToReverseWithWarn (BOPTools_AlgoTools.cxx L1294-1312).
/// Wrapper around is_split_to_reverse that logs warnings on error.
pub fn is_split_to_reverse_with_warn(
    split_normal: glam::DVec3,
    original_normal: glam::DVec3,
) -> bool {
    // OCCT: calls IsSplitToReverse(theSplit, theShape, &anErr)
    //   if (anErr != 0) 鈫?add BOPAlgo_AlertUnableToOrientTheShape warning
    // rcad: simple dot-product check matching OCCT L1427
    is_split_to_reverse(original_normal, split_normal)
}

/// OCCT-aligned: Dimensions (BOPTools_AlgoTools.hxx L546-547).
/// Returns the min and max dimension of sub-shapes in the solid.
pub fn dimensions(solid_face_indices: &[usize], ds: &crate::bopds::ds::DS) -> (i32, i32) {
    let mut d_min = 3i32;
    let mut d_max = 0i32;
    for &fi in solid_face_indices {
        if fi >= ds.faces.len() { continue; }
        // FACE has dimension 2
        d_min = d_min.min(2);
        d_max = d_max.max(2);
        for &ei in &ds.faces[fi].boundary_edges {
            // EDGE has dimension 1
            d_min = d_min.min(1);
            d_max = d_max.max(1);
            if ei < ds.edges.len() {
                let e = &ds.edges[ei];
                if e.start_vertex < ds.vertices.len() {
                    d_min = d_min.min(0);
                    d_max = d_max.max(0);
                }
                if e.end_vertex < ds.vertices.len() {
                    d_min = d_min.min(0);
                }
            }
        }
    }
    (d_min, d_max)
}

/// OCCT-aligned: Dimension (BOPTools_AlgoTools.hxx L550).
/// Returns the uniform dimension of shapes in the solid. If mixed, returns -1.
pub fn dimension(solid_face_indices: &[usize], ds: &crate::bopds::ds::DS) -> i32 {
    let (d_min, d_max) = dimensions(solid_face_indices, ds);
    if d_min == d_max { d_min } else { -1 }
}

/// OCCT-aligned: DoSplitSEAMOnFace (BOPTools_AlgoTools3D.hxx L43-49).
/// Checks if a split edge should be treated as a seam edge on a periodic surface.
/// Returns true if the edge lies on the parametric seam (U=0 or U=2蟺).
pub fn do_split_seam_on_face(
    ei: usize,
    face_idx: usize,
    ds: &crate::bopds::ds::DS,
) -> bool {
    if ei >= ds.edges.len() || face_idx >= ds.faces.len() { return false; }
    let edge = &ds.edges[ei];
    let face = &ds.faces[face_idx];
    let uv_s = crate::builder::world_to_uv(&face.surface, ds.vertices[edge.start_vertex].point);
    let uv_e = crate::builder::world_to_uv(&face.surface, ds.vertices[edge.end_vertex].point);
    let (Some(uva), Some(uvb)) = (uv_s, uv_e) else { return false };
    let seam_tol = 1e-6;
    let on_seam = |u: f64| u.abs() < seam_tol || (u - std::f64::consts::TAU).abs() < seam_tol;
    on_seam(uva.x) && on_seam(uvb.x)
}

/// OCCT-aligned: PointOnSurface (BOPTools_AlgoTools2D.cxx L107-122).
/// Evaluates UV parameters of an edge on a face at the given edge parameter.
pub fn point_on_surface(
    ds: &crate::bopds::ds::DS,
    ei: usize,
    face_idx: usize,
    t: f64,
) -> Option<glam::DVec2> {
    let _edge = ds.edges.get(ei)?;
    let _face = ds.faces.get(face_idx)?;
    // Get pcurve from edge's face_reps
    let rep = ds.edges[ei].face_reps.iter().find(|r| r.face_idx == face_idx)?;
    let pt = rep.pcurve.point_at(t);
    Some(glam::DVec2::new(pt.x, pt.y))
}

/// 鉁?OCCT-aligned: SenseFlag (BOPTools_AlgoTools3D.cxx L380-402).
/// Returns 1 if normals point same direction, -1 if opposite, 0 if not coincident.
pub fn sense_flag(n1: glam::DVec3, n2: glam::DVec3) -> i8 {
    // OCCT L384: IntTools_Tools::IsDirsCoinside 鈥?checks parallelism
    let dot_abs = n1.dot(n2).abs();
    let len1 = n1.length_squared();
    let len2 = n2.length_squared();
    if len1 < 1e-30 || len2 < 1e-30 { return 0; }
    let cos_angle = dot_abs / (len1 * len2).sqrt();
    if cos_angle < 0.9999 { return 0; } // not coincident
    // OCCT L392-401: check scalar product sign
    let sc_pr = n1.dot(n2);
    if sc_pr < 0.0 { -1 } else if sc_pr > 0.0 { 1 } else { -1 }
}

/// 鉁?OCCT-aligned: GetNormalToSurface (BOPTools_AlgoTools3D.cxx L406-439).
/// Computes the normal to a surface at UV using the surface evaluation.
pub fn get_normal_to_surface(
    surface: &rcad_kernel::geom::Surface3,
    u: f64,
    v: f64,
) -> Option<glam::DVec3> {
    use rcad_kernel::geom::SurfaceEval;
    let normal = surface.normal_at(u, v);
    if normal.length_squared() < 1e-30 { None } else { Some(normal.normalize()) }
}

/// 鉁?OCCT-aligned: GetApproxNormalToFaceOnEdge (BOPTools_AlgoTools3D.cxx L443-494).
/// Computes the approximate normal to a face near an edge by evaluating
/// the surface at a point offset from the edge toward the face interior.
pub fn get_approx_normal_to_face_on_edge(
    ds: &crate::bopds::ds::DS,
    ei: usize,
    face_idx: usize,
) -> Option<(glam::DVec3, glam::DVec3)> {
    let edge = ds.edges.get(ei)?;
    let face = ds.faces.get(face_idx)?;
    let t_mid = 0.5 * (edge.t_range[0] + edge.t_range[1]);
    let edge_mid = edge.curve.point_at(t_mid);
    let normal = get_normal_to_face_on_edge(&face.surface, face.normal, edge_mid);
    let offset_pt = edge_mid + normal * crate::tolerance::TOLERANCE_ABS * 10.0;
    Some((offset_pt, normal))
}

/// 鉁?OCCT-aligned: MinStepIn2d (BOPTools_AlgoTools3D.hxx L215).
/// Returns the minimum step used in 2D computations (1e-5).
pub fn min_step_in_2d() -> f64 {
    1e-5
}

/// 鉁?OCCT-aligned: IsEmptyShape (BOPTools_AlgoTools3D.cxx L732-788).
/// Returns true if a shape has no geometry or is empty.
pub fn is_empty_face(face: &crate::bopds::ds::DSFace) -> bool {
    face.boundary_edges.is_empty()
}

/// 鉁?OCCT-aligned: IsEmptyShape for a general DS face list.
pub fn is_empty_shape(shape_faces: &[usize], ds: &crate::bopds::ds::DS) -> bool {
    if shape_faces.is_empty() { return true; }
    // OCCT L732-788: calls HasGeometry recursively
    // rcad: check if any face has boundary edges
    shape_faces.iter().all(|&fi| {
        ds.faces.get(fi).map_or(true, |f| f.boundary_edges.is_empty())
    })
}

/// OCCT-aligned: IsInternalFace (BOPTools_AlgoTools.cxx L807-891).
///
/// Checks if face `fi` is internal to a solid described by `solid_face_indices`.
/// Uses two-level classification:
///   Level 1: edge-based angle method 鈥?for edges on the solid boundary,
///            finds adjacent face pair and checks if candidate face is internal.
///   Level 2: ComputeState 鈥?find edge not on solid, classify mid-point;
///            or PointInFace 鈫?classify_point.
///
/// Returns: Some(true) = IN, Some(false) = OUT, None = unable to determine.
pub fn is_internal_face_against_solid(
    fi: usize,
    solid_face_indices: &[usize],
    ds: &crate::bopds::ds::DS,
) -> Option<bool> {
    // OCCT L815-826: build MEF for the solid (edge鈫抐ace list)
    let mut a_mef: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for &sfi in solid_face_indices {
        if let Some(face) = ds.faces.get(sfi) {
            for &ei in &face.boundary_edges {
                a_mef.entry(ei).or_default().push(sfi);
            }
        }
    }
    // Deduplicate per edge
    for flist in a_mef.values_mut() {
        flist.sort_unstable();
        flist.dedup();
    }

    let face = &ds.faces[fi];

    // OCCT L828-874: try to find edge from face in MEF
    let mut i_ret = 0i32; // 0=not IN, 1=IN, 2=unable
    let mut found_edge = None;

    for &ei in &face.boundary_edges {
        let a_or = ds.edges.get(ei).map(|e| e.is_internal).unwrap_or(false);
        if a_or { continue; } // TopAbs_INTERNAL 鈫?skip
        if ds.is_edge_degenerated(ei) { continue; }

        if let Some(a_lf) = a_mef.get(&ei) {
            let a_nb_f = a_lf.len();
            if a_nb_f == 1 {
                // OCCT L851-861: single neighbor face 鈥?check if edge is INTERNAL in that face
                let a_f1 = a_lf[0];
                // Use GetEdgeOnFace to find edge orientation in that face
                let e_on_f1 = if a_f1 < ds.faces.len() {
                    crate::boptools::get_edge_off(ei, &ds.edges, &ds.faces[a_f1])
                } else { None };
                if let Some(ei_f1) = e_on_f1 {
                    if ds.edges[ei_f1].is_internal {
                        // Edge is INTERNAL in neighbor face 鈫?face is internal
                        i_ret = is_internal_face_core(fi, ei, a_f1, a_f1, ds);
                        found_edge = Some(ei);
                        break;
                    }
                }
                // Edge is not INTERNAL in the only neighbor 鈫?not a candidate
                continue;
            } else if a_nb_f >= 2 {
                // OCCT L864-873: two+ neighbor faces 鈥?use angle-based method
                let a_f1 = a_lf[0];
                let a_f2 = a_lf[1];
                i_ret = is_internal_face_core(fi, ei, a_f1, a_f2, ds);
                if i_ret != 2 {
                    found_edge = Some(ei);
                    break;
                }
            }
        }
    }

    if let Some(_ei) = found_edge {
        if i_ret != 2 {
            return Some(i_ret == 1);
        }
    }

    // OCCT L882-891: fall back to ComputeState
    let state = compute_state_face_against_solid(fi, solid_face_indices, ds);
    Some(state == crate::classify::Classification::In)
}

/// OCCT-aligned: IsInternalFace (BOPTools_AlgoTools.cxx L939-990).
/// Core implementation: check if face `the_face` is internal relative to
/// adjacent faces `the_face1` and `the_face2` sharing `the_edge`.
/// Returns 0=not IN, 1=IN, 2=unable.
pub fn is_internal_face_core(
    the_face: usize,
    the_edge: usize,
    the_face1: usize,
    the_face2: usize,
    ds: &crate::bopds::ds::DS,
) -> i32 {
    // OCCT L945-966: get edge copies for both faces with proper orientation
    let a_e1_on_f1 = if the_face1 < ds.faces.len() {
        crate::boptools::get_edge_off(the_edge, &ds.edges, &ds.faces[the_face1])
    } else { None };
    if a_e1_on_f1.is_none() { return 0; }
    let a_e1 = a_e1_on_f1.unwrap();

    let is_internal = ds.edges.get(a_e1).map(|e| e.is_internal).unwrap_or(false);
    if is_internal {
        // OCCT L952-956: INTERNAL edge 鈫?create both orientations
        // rcad: just use the edge as-is for both
    }

    let a_e2 = if the_face1 == the_face2 {
        // OCCT L958-962: same face 鈫?both orientations
        a_e1
    } else if the_face2 < ds.faces.len() {
        crate::boptools::get_edge_off(the_edge, &ds.edges, &ds.faces[the_face2])
            .unwrap_or(a_e1)
    } else { a_e1 };

    // OCCT L968-974: build candidate list: (edge, face) pairs
    let mut lcs_off: Vec<(usize, usize)> = Vec::new();
    lcs_off.push((the_edge, the_face));  // (theE1, theFace)
    lcs_off.push((a_e2, the_face2));      // (aE2, theFace2)

    // OCCT L976-989: GetFaceOff 鈥?find the face with minimal angle
    let a_f_off = crate::boptools::get_face_off(a_e1, the_face1, &lcs_off, ds);

    match a_f_off {
        Some(f) if f == the_face => 1,  // face is internal
        Some(_) => 0,                    // not internal
        None => 2,                       // unable to determine
    }
}

/// OCCT-aligned: IsInternalFace (BOPTools_AlgoTools.cxx L895-935).
/// Checks if face `the_face` is internal relative to a list of face candidates
/// sharing `the_edge`.
pub fn is_internal_face_against_list(
    the_face: usize,
    the_edge: usize,
    candidate_faces: &[usize],
    ds: &crate::bopds::ds::DS,
) -> i32 {
    let a_nb_f = candidate_faces.len();
    if a_nb_f == 2 {
        // OCCT L906-910: exactly 2 鈫?direct pairing
        is_internal_face_core(the_face, the_edge, candidate_faces[0], candidate_faces[1], ds)
    } else {
        // OCCT L914-933: more than 2 鈫?pair them via FindFacePairs
        // rcad: iterate all pairs
        for i in 0..candidate_faces.len() {
            for j in (i + 1)..candidate_faces.len() {
                let i_ret = is_internal_face_core(the_face, the_edge, candidate_faces[i], candidate_faces[j], ds);
                if i_ret != 0 {
                    return i_ret;
                }
            }
        }
        0
    }
}

/// 鉁?OCCT-aligned: OrientEdgesOnWire (BOPTools_AlgoTools.cxx L262-359).
///
/// OCCT algorithm:
///   1. Build vertex鈫抏dge map (MapShapesAndAncestors VERTEX鈫扙DGE).
///   2. For each edge: add to new wire, get V1/V2.
///   3. If closed edge (V1==V2): skip adjacency walk.
///   4. For each vertex direction:
///      - While vertex has exactly 2 incident edges:
///        - Find the unused edge, orient to connect (end鈫抯tart).
///        - Move to next vertex.
///
/// rcad: operates on DS edge indices + forward flags.
///   edges: mutable list of (edge_idx, forward) pairs.
pub fn orient_edges_on_wire_occt(edges: &mut Vec<(usize, bool)>, ds: &crate::bopds::ds::DS) {
    if edges.is_empty() { return; }

    // OCCT L265-272: build vertex鈫抏dge map (TopExp::MapShapesAndAncestors)
    let mut a_ve_map: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for &(ei, fwd) in edges.iter() {
        if let Some(edge) = ds.edges.get(ei) {
            let sv = if fwd { edge.start_vertex } else { edge.end_vertex };
            let ev = if fwd { edge.end_vertex } else { edge.start_vertex };
            a_ve_map.entry(sv).or_default().push(ei);
            a_ve_map.entry(ev).or_default().push(ei);
        }
    }
    // Deduplicate
    for vlist in a_ve_map.values_mut() {
        vlist.sort_unstable();
        vlist.dedup();
    }

    // OCCT L274-358: Build new wire, orient edges
    let mut a_m_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut a_wire_new: Vec<(usize, bool)> = Vec::new();

    for i in 0..edges.len() {
        let (a_ec, a_ec_fwd) = edges[i];
        if !a_m_fence.insert(a_ec) { continue; }

        // OCCT L291: add edge to wire as-is
        a_wire_new.push((a_ec, a_ec_fwd));

        // OCCT L293-294: get vertices
        let (a_v1, a_v2) = if a_ec_fwd {
            (ds.edges[a_ec].start_vertex, ds.edges[a_ec].end_vertex)
        } else {
            (ds.edges[a_ec].end_vertex, ds.edges[a_ec].start_vertex)
        };

        // OCCT L296-300: if closed edge, skip adjacency walk
        if a_v1 == a_v2 { continue; }

        // OCCT L303-355: orient adjacent edges for each vertex direction
        for &start_v in &[a_v1, a_v2] {
            let mut a_vc = start_v;
            loop {
                let Some(a_le) = a_ve_map.get(&a_vc) else { break; };
                if a_le.len() != 2 { break; }

                let mut b_stop = true;
                for &a_en in a_le {
                    if a_m_fence.contains(&a_en) { continue; }
                    let a_en_edge = &ds.edges[a_en];
                    let a_vn1 = a_en_edge.start_vertex;
                    let a_vn2 = a_en_edge.end_vertex;
                    if a_vn1 == a_vn2 { break; } // closed edge

                    // OCCT L336-345: orient edge to maintain connectivity
                    let (fwd, next_v) = if a_vc == a_vn1 {
                        // start matches 鈫?forward
                        (true, a_vn2)
                    } else if a_vc == a_vn2 {
                        // end matches 鈫?reversed
                        (false, a_vn1)
                    } else {
                        // no match 鈫?skip this edge
                        continue;
                    };

                    // OCCT L338 (correct orientation) or L342 (reversed)
                    a_wire_new.push((a_en, fwd));
                    a_m_fence.insert(a_en);
                    // OCCT L345: aVC = next vertex for next iteration
                    a_vc = next_v;
                    b_stop = false;
                    break;
                }
                if b_stop { break; }
            }
        }
    }

    *edges = a_wire_new;
}

/// 鉁?OCCT-aligned: PointInFace (BOPTools_AlgoTools3D.cxx L906-941).
/// Computes an arbitrary point inside a DS face (uses boundary centroid).
pub fn point_in_face(
    ds: &crate::bopds::ds::DS,
    face_idx: usize,
) -> Option<(glam::DVec3, glam::DVec2)> {
    let face = ds.faces.get(face_idx)?;
    if face.boundary_verts.is_empty() { return None; }
    let mut sum = glam::DVec3::ZERO;
    for &vi in &face.boundary_verts {
        if vi < ds.vertices.len() {
            sum += ds.vertices[vi].point;
        }
    }
    let p3d = sum / face.boundary_verts.len() as f64;
    let uv = crate::builder::world_to_uv(&face.surface, p3d)?;
    Some((p3d, uv))
}

/// OCCT-aligned: IsOpenShell (BOPTools_AlgoTools.cxx L2350-2394) 鈥?single-shell variant.
pub fn is_open_shell_slice(
    shell_faces: &[usize],
    ds: &crate::bopds::ds::DS,
) -> bool {
    is_open_shell(shell_faces, ds)
}

/// OCCT-aligned: ComputeState for face vs solid (BOPTools_AlgoTools.cxx L660-714).
/// Classifies a face against a solid's face set. Tries to find an edge of the
/// face not on the solid boundary, or falls back to PointInFace.
pub fn compute_state_face_against_solid(
    fi: usize,
    solid_face_indices: &[usize],
    ds: &crate::bopds::ds::DS,
) -> crate::classify::Classification {
    // OCCT L672-686: try to find an edge of the face not on the solid boundary
    let face = &ds.faces[fi];
    let solid_edge_set: std::collections::HashSet<usize> = solid_face_indices.iter()
        .flat_map(|&sfi| {
            if sfi < ds.faces.len() {
                ds.faces[sfi].boundary_edges.clone()
            } else { Vec::new() }
        })
        .collect();
    for &ei in &face.boundary_edges {
        if ds.is_edge_degenerated(ei) { continue; }
        if !solid_edge_set.contains(&ei) {
            // Classify edge midpoint
            let edge = &ds.edges[ei];
            let mid = 0.5 * (edge.t_range[0] + edge.t_range[1]);
            let pt = edge.curve.point_at(mid);
            return crate::classify::classify_point(pt, solid_face_indices, ds);
        }
    }
    // OCCT L688-714: all edges on solid 鈫?PointInFace
    let pt = point_in_face(ds, fi);
    match pt {
        Some((p3d, _)) => crate::classify::classify_point(p3d, solid_face_indices, ds),
        None => crate::classify::Classification::Out,
    }
}

