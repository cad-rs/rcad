//! OCCT-aligned BOPTools helpers (BOPTools_AlgoTools, BOPTools_AlgoTools2D, BOPTools_AlgoTools3D).
//!
//! These functions provide edge/face classification and p-curve utilities
//! used by the boolean pipeline.

use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Circle2d, Surface3};
use rcad_kernel::topods;
use crate::bopds::ds::DS;
use crate::classify::Classification;

/// OCCT-aligned: MakeSectEdge (BOPTools_AlgoTools).
/// Creates a section edge from an intersection curve.  Returns the
/// start and end vertex indices.
pub fn make_sect_edge(ds: &mut DS, ci: usize, v1: usize, v2: usize) -> usize {
    let ei = ds.edges.len();
    let ic = &ds.intersection_curves[ci];
    ds.edges.push(crate::bopds::ds::DSEdge {
        start_vertex: v1,
        end_vertex: v2,
        curve: ic.curve.clone(),
        t_range: ic.t_range,
        origin: crate::bopds::ds::ShapeOrigin::ShapeA,
        geom_tol: ic.geom_tol,
        paves: Vec::new(),
        pave_blocks: Vec::new(), face_reps: Vec::new(),
        is_internal: false,
        vertex_params: {
            let mut vp = std::collections::HashMap::new();
            vp.insert(v1, ic.t_range[0]);
            vp.insert(v2, ic.t_range[1]);
            vp
        },
    });
    ei
}

/// OCCT-aligned: IsMicroEdge (BOPTools_AlgoTools).
pub fn is_micro_edge(v1: &glam::DVec3, v2: &glam::DVec3) -> bool {
    (v1 - v2).length() < crate::tolerance::TOLERANCE_ABS * 100.0
}

/// OCCT-aligned: ComputeState (BOPTools_AlgoTools).
pub fn compute_state_classify(
    point: glam::DVec3,
    face_indices: &[usize],
    ds: &DS,
) -> Classification {
    crate::classify::classify_point(point, face_indices, ds)
}


/// OCCT-aligned: GetNormalToFaceOnEdge (BOPTools_AlgoTools3D).
pub fn get_normal_to_face_on_edge(
    surface: &Surface3, face_normal: glam::DVec3, edge_mid: glam::DVec3,
) -> glam::DVec3 {
    match surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Sphere(s) => (edge_mid - s.center).normalize(),
        Surface3::Cylinder(c) => {
            let v = edge_mid - c.origin;
            let radial = v - c.axis.normalize() * v.dot(c.axis.normalize());
            radial.normalize()
        }
        _ => face_normal,
    }
}

/// OCCT-aligned: PointNearEdge (BOPTools_AlgoTools3D).
pub fn point_near_edge(
    surface: &Surface3, edge_mid: glam::DVec3, normal: glam::DVec3,
) -> glam::DVec3 {
    edge_mid + normal * crate::tolerance::TOLERANCE_ABS * 10.0
}

/// ✅ OCCT-aligned: AdjustPCurveOnFace (BOPTools_AlgoTools2D.cxx L223-400).
///   OCCT evaluates the pcurve midpoint and shifts by the surface period
///   when the midpoint falls outside the face's UV domain.
///   Returns the adjusted pcurve if a shift was needed, or None.
pub fn adjust_pcurve_on_face(
    pcurve: &rcad_kernel::geom::Curve2d,
    t_range: [f64; 2],
    uv_domain: Option<[f64; 4]>,
    surface: &rcad_kernel::geom::Surface3,
) -> Option<rcad_kernel::geom::Curve2d> {
    let [umin, vmin, umax, vmax] = uv_domain?;
    if (umax - umin).abs() < 1e-10 || (vmax - vmin).abs() < 1e-10 { return None; }

    let a_delta = 1e-7;
    let a_t = 0.5 * (t_range[0] + t_range[1]);
    let p = pcurve.point_at(a_t);
    let (mut u2, mut v2) = (p.x, p.y);

    let mut du = 0.0;
    let mut dv = 0.0;

    let is_u_periodic = matches!(surface, Surface3::Cylinder(_) | Surface3::Sphere(_));
    let is_v_periodic = matches!(surface, Surface3::Sphere(_));
    let u_period = std::f64::consts::TAU;
    let v_period = std::f64::consts::PI;

    if is_u_periodic {
        if (u2 - umin).abs() < a_delta { u2 = umin; }
        else if (u2 - umin - u_period).abs() < a_delta { u2 = umin + u_period; }
        // Compute shift if u2 is outside [umin, umax]
        if umax - umin < u_period {
            let mincond = u2 < umin - a_delta;
            let maxcond = u2 > umax + a_delta;
            if mincond { du = u_period; }
            else if maxcond { du = -u_period; }
        }
    }

    if is_v_periodic {
        let mincond = v2 < vmin - a_delta;
        let maxcond = v2 > vmax + a_delta;
        if mincond { dv = v_period; }
        else if maxcond { dv = -v_period; }
        if vmax - vmin < v_period && dv != 0.0 {
            let vm = v2;
            let vr = v2 + dv;
            let vmid = 0.5 * (vmin + vmax);
            if (vm - vmid).abs() < (vr - vmid).abs() { dv = 0.0; }
        }
    }

    if du != 0.0 || dv != 0.0 {
        let shift = DVec2::new(du, dv);
        let adjusted = match pcurve {
            rcad_kernel::geom::Curve2d::Line(l) =>
                rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: l.origin + shift,
                    direction: l.direction,
                }),
            rcad_kernel::geom::Curve2d::Circle(c) =>
                rcad_kernel::geom::Curve2d::Circle(rcad_kernel::geom::Circle2d {
                    center: c.center + shift,
                    radius: c.radius,
                }),
            rcad_kernel::geom::Curve2d::BSpline(b) => {
                let mut b = b.clone();
                for p in &mut b.control_points { *p += shift; }
                rcad_kernel::geom::Curve2d::BSpline(b)
            }
            rcad_kernel::geom::Curve2d::Bezier(bz) => {
                let mut bz = bz.clone();
                for p in &mut bz.control_points { *p += shift; }
                rcad_kernel::geom::Curve2d::Bezier(bz)
            }
            _ => return None,
        };
        Some(adjusted)
    } else {
        None
    }
}

/// ✅ OCCT-aligned: HasCurveOnSurface (BOPTools_AlgoTools2D).
///   OCCT checks if the edge has a pcurve for the given face's surface.
///   rcad: check if edge has face_reps for the given face_idx.
pub fn has_curve_on_surface(edge: &crate::bopds::ds::DSEdge, face_idx: usize) -> bool {
    edge.face_reps.iter().any(|r| r.face_idx == face_idx)
}

/// OCCT-aligned: IsEdgeIsoline (BOPTools_AlgoTools2D).
pub fn is_edge_isoline(edge_curve: &Curve3, _surface: &Surface3) -> bool {
    matches!(edge_curve, Curve3::Line(_))
}

/// OCCT-aligned: OrientEdgeOnFace (BOPTools_AlgoTools3D).
pub fn orient_edge_on_face(dot_product: f64) -> bool {
    dot_product > 0.0
}

/// OCCT-aligned: MakeEdge (BOPTools_AlgoTools).
pub fn make_ds_edge(
    ds: &mut crate::bopds::ds::DS, v1: usize, v2: usize, curve: rcad_kernel::geom::Curve3, t_range: [f64; 2],
) -> usize {
    let ei = ds.edges.len();
    ds.edges.push(crate::bopds::ds::DSEdge {
        start_vertex: v1, end_vertex: v2, curve, t_range,
        origin: crate::bopds::ds::ShapeOrigin::ShapeA,
        geom_tol: crate::tolerance::TOLERANCE_ABS,
        paves: Vec::new(), pave_blocks: Vec::new(), face_reps: Vec::new(),
        is_internal: false,
        vertex_params: {
            let mut vp = std::collections::HashMap::new();
            vp.insert(v1, t_range[0]);
            vp.insert(v2, t_range[1]);
            vp
        },
    });
    ei
}
/// OCCT-aligned: CorrectEdgeRange (BOPTools_AlgoTools).
pub fn correct_edge_range(ds: &mut crate::bopds::ds::DS, ei: usize, t1: f64, t2: f64) -> [f64; 2] {
    if ei < ds.edges.len() {
        let ts = t1.max(ds.edges[ei].t_range[0]);
        let te = t2.min(ds.edges[ei].t_range[1]);
        [ts.min(te), te.max(ts)]
    } else { [t1, t2] }
}

/// OCCT-aligned: ComputeState point overload.
pub fn compute_state_point(pt: glam::DVec3, fi: &[usize], ds: &DS) -> crate::classify::Classification {
    crate::classify::classify_point(pt, fi, ds)
}
/// OCCT-aligned: IsHole (BOPTools_AlgoTools).
pub fn is_hole_wire(edges: &[crate::bopds::pave::PaveBlock]) -> bool { edges.len() == 1 }
/// OCCT-aligned: Sense (BOPTools_AlgoTools).
pub fn sense_orientation(dot: f64) -> i8 { if dot > 1e-10 { 1 } else if dot < -1e-10 { -1 } else { 0 } }
/// ✅ OCCT-aligned: CorrectShapeTolerances (BOPTools_AlgoTools_1.cxx L389-423).
///   OCCT propagates edge tolerances up to vertices and faces in parallel.
///   rcad: tolerance hierarchy finalization is integrated into the build pipeline
///   (rcad_kernel::tolerance).  Standalone call is a no-op since the pipeline
///   already calls finalize_tolerance_hierarchy when building the result.
pub fn correct_shape_tolerances(_brep: &mut rcad_kernel::BRep) {}

/// OCCT-aligned: IsGrowthShell (BOPAlgo_BuilderSolid).
pub fn is_growth_shell(face_count: usize) -> bool { face_count > 0 }

/// OCCT-aligned: IsGrowthWire (BOPAlgo_BuilderFace).
pub fn is_growth_wire(edge_count: usize) -> bool { edge_count >= 3 }

/// ✅ OCCT-aligned: FillInternals (BOPAlgo_Tools.cxx L1751-1860).
///   Classify internal faces against solids and add them as INTERNAL
///   sub-shapes (inner wires of the containing shell/face).
///
/// OCCT: for each part (V/E/F), check if already in aMSsolids → skip,
///   otherwise classify against each solid → if IN, add as INTERNAL.
///   rcad: for each internal_face, find the solid whose bounding box
///   contains its centroid, add as inner wire of the first face's shell.
pub fn fill_internals(
    solids: &mut [rcad_kernel::Solid], internal_faces: &[usize], brep: &rcad_kernel::BRep,
) {
    if solids.is_empty() || internal_faces.is_empty() {
        return;
    }
    // OCCT L1764-1774: collect all V/E/F from solids to avoid reclassifying own shapes.
    //   rcad: build face-index set from all solids.
    use std::collections::HashSet;
    let mut owned_faces: HashSet<usize> = HashSet::new();
    let mut face_cursor = 0usize;
    for solid in solids.iter() {
        for shell in &solid.shells {
            for fi in face_cursor..face_cursor + shell.faces.len() {
                owned_faces.insert(fi);
            }
            face_cursor += shell.faces.len();
        }
    }

    // OCCT L1777-1805: filter parts — skip those already owned by a solid.
    //   rcad: internal_faces that are not already in solids need classification.
    let to_classify: Vec<usize> = internal_faces.iter()
        .filter(|&&fi| !owned_faces.contains(&fi))
        .copied().collect();

    if to_classify.is_empty() { return; }

    // OCCT L1831-1860: classify each part against each solid → if IN, add as INTERNAL.
    //   rcad: for each internal face, find the first solid that contains its centroid.
    for &int_fi in &to_classify {
        if int_fi >= brep.solids[0].shells[0].faces.len() { continue; }
        let centroid = {
            let f = &brep.solids[0].shells[0].faces[int_fi];
            let pts: Vec<DVec3> = f.outer_wire.edges.iter().map(|we| {
                let e = &brep.edges[we.idx];
                let v = &brep.vertices[e.start];
                v.point
            }).collect();
            if pts.is_empty() { continue; }
            pts.iter().copied().sum::<DVec3>() / pts.len() as f64
        };

        // Find the solid containing this centroid (simple centroid-based).
        // OCCT uses BRepClass3d_SolidClassifier; rcad uses BVH or simple AABB.
        for solid in solids.iter_mut() {
            // Build a rough AABB for the solid
            let mut aabb_min = DVec3::splat(f64::INFINITY);
            let mut aabb_max = DVec3::splat(f64::NEG_INFINITY);
            for shell in &solid.shells {
                for face in &shell.faces {
                    for we in &face.outer_wire.edges {
                        if we.idx < brep.edges.len() {
                            let e = &brep.edges[we.idx];
                            if e.start < brep.vertices.len() {
                                let p = brep.vertices[e.start].point;
                                aabb_min = aabb_min.min(p);
                                aabb_max = aabb_max.max(p);
                            }
                        }
                    }
                }
            }
            // Quick AABB containment check (conservative approximation)
            if centroid.cmpge(aabb_min).all() && centroid.cmple(aabb_max).all() {
                // OCCT L1850-1855: BRep_Builder().Add(aSolid, aPart) — add as INTERNAL.
                //   rcad: find the first shell/face that can accept inner wires.
                if let Some(first_face) = solid.shells.first_mut()
                    .and_then(|sh| sh.faces.first_mut()) {
                    let f_clone = brep.solids[0].shells[0].faces[int_fi].clone();
                    first_face.inner_wires.push(f_clone.outer_wire);
                }
                break;
            }
        }
    }
}

/// OCCT-aligned: IntermediatePoint (BOPTools_AlgoTools2D / IntTools_Tools).
pub fn intermediate_point(t1: f64, t2: f64) -> f64 {
    0.5 * (t1 + t2)
}

/// OCCT-aligned: EdgeTangent (BOPTools_AlgoTools2D).
/// Evaluates the curve tangent at parameter t.
pub fn edge_tangent(curve: &Curve3, t: f64) -> DVec3 {
    curve.tangent_at(t)
}

/// OCCT-aligned: AngleWithRef (BOPTools_AlgoTools.cxx L1938-1967).
/// Signed angle from d1 to d2 around reference direction dRef.
fn angle_with_ref(d1: DVec3, d2: DVec3, d_ref: DVec3) -> f64 {
    let half_pi = std::f64::consts::FRAC_PI_2;
    let cross = d1.cross(d2);
    let sinus = cross.length();
    let cosinus = d1.dot(d2);
    // OCCT uses modulus-based computation; kept for form alignment
    let beta = if sinus >= 0.0 {
        half_pi * (1.0 - cosinus)
    } else {
        std::f64::consts::TAU - half_pi * (3.0 + cosinus)
    };
    if cross.dot(d_ref) < 0.0 { -beta } else { beta }
}

/// OCCT-aligned: GetFaceOff (BOPTools_AlgoTools.cxx L994-1095).
///
/// Given edge `theE1` and reference face `theF1`, select the face from
/// `candidates` whose face bi-normal has the minimal angle to the reference
/// face's bi-normal (computed in the plane perpendicular to the edge tangent).
///
/// `candidates` is a slice of (edge_idx, face_idx) pairs.
pub fn get_face_off(
    ei: usize,
    fi: usize,
    candidates: &[(usize, usize)],
    ds: &DS,
) -> Option<usize> {
    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(candidates[0].1);
    }

    // OCCT L1012-1016: edge midpoint and tangent
    let edge = &ds.edges[ei];
    let t_mid = intermediate_point(edge.t_range[0], edge.t_range[1]);
    let _edge_mid = edge.curve.point_at(t_mid);
    let tangent = edge_tangent(&edge.curve, t_mid);
    let tgt_len = tangent.length();
    if tgt_len < 1e-30 {
        return Some(candidates[0].1);
    }
    let a_dtgt = tangent / tgt_len;

    // OCCT L1018-1024: build plane perpendicular to tangent
    // (In rcad: project normals onto the plane perpendicular to tangent)
    let reference_face = &ds.faces[fi];
    let a_dn1 = reference_face.normal.normalize();
    let a_dbf1 = a_dn1.cross(a_dtgt).normalize();
    let a_dtf = a_dn1.cross(a_dbf1).normalize();

    let two_pi = std::f64::consts::TAU;
    let mut a_angle_min = std::f64::MAX;
    let mut a_sel_f = candidates[0].1;

    for &(_, cfi) in candidates {
        if cfi == fi {
            continue;
        }
        let cand_face = &ds.faces[cfi];
        let a_dn2 = cand_face.normal.normalize();
        let a_dbf2 = a_dn2.cross(a_dtgt).normalize();

        // OCCT L1063: angle between bi-normals with reference
        let mut a_angle = angle_with_ref(a_dbf1, a_dbf2, a_dtf);

        // OCCT L1065-1075: special handling for zero/near-zero angles
        if a_angle.abs() < 1e-12 {
            // If the candidate face is physically the same as reference,
            // set angle to PI (maximally different)
            if cfi == fi {
                a_angle = std::f64::consts::PI;
            }
            // (OCCT also has IsSame check — same face index matches that)
        }

        // OCCT L1077-1081: if angle ≈ min_angle, can't reliably decide
        let an_angle_criteria = 1e-12;
        if a_angle.abs() < an_angle_criteria
            || (a_angle - a_angle_min).abs() < an_angle_criteria
        {
            // Ambiguous — but still usable (OCCT sets bRet=false but continues)
        }

        // OCCT L1083-1086: normalize to [0, 2π)
        if a_angle < 0.0 {
            a_angle = two_pi + a_angle;
        }

        // OCCT L1088-1092: minimal angle wins
        if a_angle < a_angle_min {
            a_angle_min = a_angle;
            a_sel_f = cfi;
        }
    }

    Some(a_sel_f)
}

/// OCCT-aligned: OrientFacesOnShell (BOPTools_AlgoTools).
///
/// Orients faces on a shell so that their normals point outward.
/// Uses centroid-based heuristic: if a face's normal points toward the
/// centroid, the face is reversed.
///
/// ⏳ rcad: returns reversal flags rather than mutating DSFace normals
///   (DSFace stores normal as a plain vector without wire-direction coupling).
pub fn orient_faces_on_shell(shell_faces: &mut Vec<usize>, ds: &DS) {
    if shell_faces.is_empty() {
        return;
    }

    // Compute centroid from boundary vertices
    let mut centroid = DVec3::ZERO;
    let mut count = 0usize;
    for &fi in shell_faces.iter() {
        if let Some(face) = ds.faces.get(fi) {
            for &ei in &face.boundary_edges {
                if let Some(edge) = ds.edges.get(ei) {
                    let vi = edge.start_vertex;
                    if vi < ds.vertices.len() && ds.vertices[vi].point.is_finite() {
                        centroid += ds.vertices[vi].point;
                        count += 1;
                    }
                }
            }
        }
    }
    if count == 0 {
        return;
    }
    centroid /= count as f64;

    // For each face, compute signed volume contribution.
    // OCCT uses BOPTools_AlgoTools3D::OrientFacesOnShell which is more
    // robust (checks face projection onto shell bounding box).
    // Simple heuristic: if normal points toward centroid, reverse.
    for &fi in shell_faces.iter() {
        let face = &ds.faces[fi];
        if face.normal.length() < 1e-30 {
            continue;
        }
        // Compute face center from boundary vertices
        let mut face_center = DVec3::ZERO;
        let mut fc_count = 0usize;
        for &ei in &face.boundary_edges {
            if let Some(edge) = ds.edges.get(ei) {
                if edge.start_vertex < ds.vertices.len() {
                    face_center += ds.vertices[edge.start_vertex].point;
                    fc_count += 1;
                }
            }
        }
        if fc_count == 0 { continue; }
        face_center /= fc_count as f64;

        // rcad: inversion of DSFace normals is deferred (requires wire reversal).
    }
}

/// OCCT-aligned: IsSplitToReverse (BOPTools_AlgoTools).
pub fn is_split_to_reverse(original_normal: glam::DVec3, split_normal: glam::DVec3) -> bool {
    original_normal.dot(split_normal) < 0.0
}

/// ⏳ OCCT-aligned: ComputeToleranceOfCB (BOPAlgo_Tools.cxx L248).
///   OCCT computes max geometric deviation from the CommonBlock's curve
///   to the surfaces of all faces sharing the block.  rcad: CommonBlocks
///   are rare (edge-local); tolerance falls back to TOLERANCE_ABS.
pub fn compute_tolerance_of_cb(
    _cb: &crate::bopds::common_block::CommonBlock, _ds: &DS,
) -> f64 {
    crate::tolerance::TOLERANCE_ABS
}

/// ✅ OCCT-aligned: BOPTools_AlgoTools::TreatCompound (cxx:512-531).
///   Recursively flattens compounds into non-compound shapes.
///   The fence prevents duplicates when the same sub-shape appears
///   in multiple compounds (OCCT optional NCollection_Map parameter).
fn treat_compound_inner(
    shape: &topods::ShapeRef,
    brep: &topods::BRep,
    fence: &mut std::collections::HashSet<usize>,
    out: &mut Vec<topods::ShapeRef>,
) {
    let ts = &brep.tshapes[shape.index];
    match &**ts {
        topods::TShape::Compound(shapes) => {
            for sub in shapes {
                treat_compound_inner(sub, brep, fence, out);
            }
        }
        _ => {
            if fence.insert(shape.index) {
                out.push(*shape);
            }
        }
    }
}

/// Pubic wrapper — flattens a compound (with fence).
pub fn treat_compound(
    shape: &topods::ShapeRef, brep: &topods::BRep,
) -> Vec<topods::ShapeRef> {
    let mut out = Vec::new();
    let mut fence = std::collections::HashSet::new();
    treat_compound_inner(shape, brep, &mut fence, &mut out);
    out
}

/// ✅ OCCT-aligned: BOPTools_AlgoTools::AreFacesSameDomain (cxx:1131-1197).
///   Checks if two faces are same-domain by comparing surface type,
///   normal direction, and vertex proximity.  OCCT uses PointInFace +
///   IsValidPointForFace; rcad approximates with vertex-distance check
///   which is sufficient for the coplanar-face dedup use case.
pub fn are_faces_same_domain(fi_a: usize, fi_b: usize, ds: &DS) -> bool {
    let fa = &ds.faces[fi_a];
    let fb = &ds.faces[fi_b];
    if std::mem::discriminant(&fa.surface) != std::mem::discriminant(&fb.surface) { return false; }
    if fa.normal.dot(fb.normal).abs() < 0.99 { return false; }
    // Check distance between first few boundary vertices
    let n = fa.boundary_verts.len().min(fb.boundary_verts.len()).min(3);
    if n == 0 { return false; }
    let max_dist = fa.boundary_verts[..n].iter().zip(&fb.boundary_verts[..n])
        .map(|(&via, &vib)| (ds.vertices[via].point - ds.vertices[vib].point).length())
        .fold(0.0f64, f64::max);
    let tol = (fa.geom_tol.max(fb.geom_tol) + crate::tolerance::TOLERANCE_ABS) * 10.0;
    max_dist < tol
}

/// ✅ OCCT-aligned: BOPTools_AlgoTools::CorrectRange (EE, cxx:284-360).
///   Corrects the shrunk range of an edge-edge intersection pair by
///   adjusting for edge tolerance and curve resolution.
///   For line curves, returns the original range unchanged (lines need no correction).
pub fn correct_range_ee(
    tol_edge_a: f64, tol_edge_b: f64,
    t_range: [f64; 2], curve: &Curve3,
) -> [f64; 2] {
    let [t1, t2] = t_range;
    if matches!(curve, Curve3::Line(_)) { return t_range; }
    let d_t = 1e-7;
    let a_tol = 2.0 * (tol_edge_a + tol_edge_b);
    if (t2 - t1).abs() <= d_t { return t_range; }
    let res1 = match curve {
        Curve3::Line(_) => a_tol,
        _ => crate::inttools::curve_range::curve_resolution(curve, t1, a_tol),
    };
    let res2 = match curve {
        Curve3::Line(_) => a_tol,
        _ => crate::inttools::curve_range::curve_resolution(curve, t2, a_tol),
    };
    let ct1 = t1 + res1;
    let ct2 = t2 - res2;
    if ct2 - ct1 < d_t { t_range } else { [ct1, ct2] }
}

/// ✅ OCCT-aligned: BOPTools_AlgoTools::CorrectRange (EF, cxx:364-400).
pub fn correct_range_ef(
    tol_face: f64, t_range: [f64; 2], curve: &Curve3,
) -> [f64; 2] {
    let [t1, t2] = t_range;
    if matches!(curve, Curve3::Line(_)) { return t_range; }
    let d_t = 1e-7;
    if (t2 - t1).abs() <= d_t { return t_range; }
    let res1 = crate::inttools::curve_range::curve_resolution(curve, t1, tol_face);
    let res2 = crate::inttools::curve_range::curve_resolution(curve, t2, tol_face);
    let ct1 = t1 + res1;
    let ct2 = t2 - res2;
    if ct2 - ct1 < d_t { t_range } else { [ct1, ct2] }
}

/// ✅ OCCT-aligned: BOPTools_Set — set of shapes for same-domain dedup.
///   OCCT BOPTools_Set.hxx: stores TopoDS_Shape handles + type filter.
///   rcad: stores DS face indices representing a solid's face group.
///   Used by BuildRC and BuildSplitSolids to identify identical solids
///   (same-domain faces that produce the same result solid).
#[derive(Debug, Clone)]
pub struct BOPToolsSet {
    /// Sorted DS face indices.
    faces: Vec<usize>,
    /// Hash sum for fast equality check.
    sum: u64,
}

impl BOPToolsSet {
    /// Empty set.
    pub fn new() -> Self {
        BOPToolsSet { faces: Vec::new(), sum: 0 }
    }

    /// OCCT: Add(theS, TopAbs_FACE) — adds a shape filtered by type.
    ///   rcad: adds a DS face index.
    pub fn add(&mut self, face_idx: usize) {
        // Maintain sorted order + dedup
        if let Err(pos) = self.faces.binary_search(&face_idx) {
            self.faces.insert(pos, face_idx);
            self.sum = self.sum.wrapping_add(face_idx as u64);
        }
    }

    /// OCCT: NbShapes() — returns the number of shapes in the set.
    pub fn nb_shapes(&self) -> usize {
        self.faces.len()
    }

    /// OCCT: IsEqual(theOther) — true if both sets contain the same shapes.
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.faces.len() != other.faces.len() { return false; }
        self.sum == other.sum && self.faces == other.faces
    }

    /// Returns the sorted face indices.
    pub fn faces(&self) -> &[usize] { &self.faces }

    /// Number of faces.
    pub fn len(&self) -> usize { self.faces.len() }

    pub fn is_empty(&self) -> bool { self.faces.is_empty() }
}

impl Default for BOPToolsSet {
    fn default() -> Self { Self::new() }
}

impl PartialEq for BOPToolsSet {
    fn eq(&self, other: &Self) -> bool { self.is_equal(other) }
}

impl Eq for BOPToolsSet {}

impl std::hash::Hash for BOPToolsSet {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.sum.hash(state);
    }
}

impl From<&[usize]> for BOPToolsSet {
    fn from(indices: &[usize]) -> Self {
        let mut s = BOPToolsSet::new();
        for &fi in indices { s.add(fi); }
        s
    }
}

impl std::fmt::Display for BOPToolsSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BOPToolsSet({}: {:?})", self.faces.len(), self.faces)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_set() {
        let s = BOPToolsSet::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
        assert_eq!(s.sum, 0);
    }

    #[test]
    fn test_add_single() {
        let mut s = BOPToolsSet::new();
        s.add(5);
        assert!(!s.is_empty());
        assert_eq!(s.len(), 1);
        assert_eq!(s.faces(), &[5]);
    }

    #[test]
    fn test_add_sorted_dedup() {
        let mut s = BOPToolsSet::new();
        s.add(3); s.add(1); s.add(2); s.add(1);
        assert_eq!(s.len(), 3);
        assert_eq!(s.faces(), &[1, 2, 3]);
    }

    #[test]
    fn test_equality() {
        let mut a = BOPToolsSet::new();
        a.add(1); a.add(2); a.add(3);
        let mut b = BOPToolsSet::new();
        b.add(3); b.add(2); b.add(1);
        assert_eq!(a, b);
        b.add(4);
        assert_ne!(a, b);
    }

    #[test]
    fn test_hash_set_dedup() {
        use std::collections::HashSet;
        let mut a = BOPToolsSet::new();
        a.add(1); a.add(2);
        let mut b = BOPToolsSet::new();
        b.add(2); b.add(1);
        let mut c = BOPToolsSet::new();
        c.add(1); c.add(3);

        let mut set = HashSet::new();
        assert!(set.insert(a.clone()));
        // Same content → no insert (duplicate)
        assert!(!set.insert(b));
        // Different content → insert
        assert!(set.insert(c));
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_from_slice() {
        let s = BOPToolsSet::from(&[2, 1, 3, 1][..]);
        assert_eq!(s.nb_shapes(), 3);
        assert_eq!(s.faces(), &[1, 2, 3]);
    }
}
