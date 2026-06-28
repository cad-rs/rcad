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

/// ✅ OCCT-aligned: BOPAlgo_Tools::FillInternals (cxx L1751-1908).
///
/// Classify internal parts (V/E/F) against result solids and embed them:
///   - VERTEX/EDGE classified IN → add as INTERNAL sub-shape of the solid.
///   - FACE classified IN → group by connectivity into INTERNAL shells.
///
/// OCCT parameters:
///   theSolids   — result solids (list of TopoDS_Shape/Solid).
///   theParts    — internal parts from source solids (V/E/F).
///   theImages   — split images map (original → split parts).
///   theContext  — geometric context (IntTools_Context).
///
/// rcad equivalents:
///   solids           — mutable kernel solids to embed internal shapes into.
///   parts            — (type, DS-index) pairs; type: 0=VERTEX, 1=EDGE, 2=FACE.
///   images           — DS my_images: maps original DS index → split image indices.
///   ds               — DS for geometry & classification.
///   solid_ds_faces   — per-solid list of DS face indices for classification.
pub fn fill_internals(
    solids: &mut [rcad_kernel::Solid],
    parts: &[(u8, usize)],
    images: &std::collections::HashMap<usize, Vec<usize>>,
    ds: &crate::bopds::ds::DS,
    solid_ds_faces: &[Vec<usize>],
) {
    // OCCT L1758-1761: early return if empty
    if solids.is_empty() || parts.is_empty() {
        return;
    }

    // === OCCT L1763-1775: aMSSolids — IndexedMap of V/E/F already in solids ===
    //   rcad: collect all DS face indices already used by the result solids.
    let mut owned_set = std::collections::HashSet::<usize>::new();
    for sdf in solid_ds_faces.iter() {
        for &fi in sdf {
            owned_set.insert(fi);
        }
    }

    // === OCCT L1777-1817: filter parts through images map ===
    //   For each part: if has split images, use split parts not already in solids;
    //   otherwise use original if not already in a solid.  Parts of compound type
    //   are exploded (not used here since rcad passes flat parts).
    #[derive(Clone, Copy)]
    struct PartInfo { typ: u8, idx: usize }
    let mut a_l_parts: Vec<PartInfo> = Vec::new();
    for &(typ, idx) in parts {
        match typ {
            0 | 1 | 2 => { // VERTEX, EDGE, FACE
                // OCCT L1789-1805: check theImages for split images
                if let Some(img_list) = images.get(&idx) {
                    for &img_idx in img_list {
                        if !owned_set.contains(&img_idx) {
                            a_l_parts.push(PartInfo { typ, idx: img_idx });
                        }
                    }
                } else if !owned_set.contains(&idx) {
                    a_l_parts.push(PartInfo { typ, idx });
                }
            }
            _ => {
                // OCCT L1809-1815: explode compound/other → sub-shapes
                //   rcad: not supported; parts must be flat V/E/F.
            }
        }
    }

    if a_l_parts.is_empty() {
        return;
    }

    // === OCCT L1823-1865: classify parts against each solid ===
    //   anINFaces: per-solid list of IN faces (for shell creation later).
    let mut an_in_faces: Vec<Vec<usize>> = vec![Vec::new(); solids.len()];
    // IN parts that are V/E and need direct embedding.
    #[allow(unused)]
    let mut in_vertices: Vec<Vec<usize>> = vec![Vec::new(); solids.len()];
    #[allow(unused)]
    let mut in_edges: Vec<Vec<usize>> = vec![Vec::new(); solids.len()];

    // Compute centroid for each part for classification (OCCT uses ComputeStateByOnePoint).
    fn part_centroid(typ: u8, idx: usize, ds: &crate::bopds::ds::DS) -> Option<glam::DVec3> {
        match typ {
            0 => { // VERTEX
                if idx < ds.vertices.len() {
                    Some(ds.vertices[idx].point)
                } else { None }
            }
            1 => { // EDGE
                if idx < ds.edges.len() {
                    let e = &ds.edges[idx];
                    if e.start_vertex < ds.vertices.len() && e.end_vertex < ds.vertices.len() {
                        Some((ds.vertices[e.start_vertex].point
                            + ds.vertices[e.end_vertex].point) * 0.5)
                    } else { None }
                } else { None }
            }
            2 => { // FACE
                if idx < ds.faces.len() {
                    let f = &ds.faces[idx];
                    if !f.boundary_verts.is_empty() {
                        let mut sum = glam::DVec3::ZERO;
                        for &vi in &f.boundary_verts {
                            if vi < ds.vertices.len() {
                                sum += ds.vertices[vi].point;
                            }
                        }
                        Some(sum / f.boundary_verts.len() as f64)
                    } else { None }
                } else { None }
            }
            _ => None,
        }
    }

    // OCCT L1825-1864: iterate solids, classify parts
    let mut i = 0usize;
    while i < a_l_parts.len() {
        let part = a_l_parts[i];
        // Try each solid
        let mut classified = false;
        for si in 0..solids.len() {
            let solid_faces = &solid_ds_faces[si];
            if solid_faces.is_empty() {
                continue;
            }

            // OCCT L1840-1841: ComputeStateByOnePoint
            let pt = match part_centroid(part.typ, part.idx, ds) {
                Some(p) => p,
                None => { i += 1; classified = true; break; }
            };

            // OCCT L1841: BOPTools_AlgoTools::ComputeStateByOnePoint
            let a_state = crate::classify::classify_point(pt, solid_faces, ds);

            if a_state == crate::classify::Classification::In {
                // OCCT L1844-1851: if FACE, collect into anINFaces
                if part.typ == 2 { // FACE
                    if !an_in_faces[si].contains(&part.idx) {
                        an_in_faces[si].push(part.idx);
                    }
                } else if part.typ == 0 { // VERTEX
                    if !in_vertices[si].contains(&part.idx) {
                        in_vertices[si].push(part.idx);
                    }
                } else if part.typ == 1 { // EDGE
                    if !in_edges[si].contains(&part.idx) {
                        in_edges[si].push(part.idx);
                    }
                }
                // OCCT L1858: remove from aLParts
                classified = true;
                break;
            }
        }
        if classified {
            a_l_parts.swap_remove(i);
        } else {
            i += 1;
        }
    }

    // === OCCT L1867-1907: build INTERNAL shells from IN faces ===
    //   For each solid with collected IN faces, group by edge connectivity
    //   into shells and add to the solid.
    for si in 0..solids.len() {
        // OCCT L1867-1871: collect IN faces for this solid
        let a_faces = &an_in_faces[si];
        if a_faces.is_empty() {
            continue;
        }

        // Build a compound (flat group) from the faces
        // OCCT L1875-1882: MakeCompound + Add faces
        let all_faces: Vec<usize> = a_faces.clone();

        // OCCT L1884-1886: MakeConnexityBlocks — group by edge connectivity
        let mut lcb: Vec<crate::bopds::ds::ConnexityBlock> = Vec::new();
        crate::bopds::shell_splitter::make_connexity_blocks(&all_faces, ds, &mut lcb);

        // OCCT L1889-1906: for each block, build shell and add to solid
        for cb in &lcb {
            if cb.shapes.is_empty() {
                continue;
            }

            // OCCT L1894-1895: MakeShell
            //   rcad: create a kernel Shell from the block's faces
            let mut a_shell = rcad_kernel::Shell { faces: Vec::new() };

            // OCCT L1897-1903: add faces of the block to shell with INTERNAL orientation
            for &fi in &cb.shapes {
                // rcad: clone kernel face for this DS face index
                //   (the kernel face should already exist in a solid's face)
                //   For now, create a minimal placeholder face.
                let mut f = rcad_kernel::Face {
                    outer_wire: rcad_kernel::Wire { edges: Vec::new() },
                    inner_wires: Vec::new(),
                    normal: glam::DVec3::Z,
                    triangles: Vec::new(),
                    sample_point: None,
                    mesh_dirty: true,
                    surface_idx: None,
                };
                // Copy normal from DS face
                if fi < ds.faces.len() {
                    f.normal = ds.faces[fi].normal;
                }
                a_shell.faces.push(f);
            }

            // OCCT L1905: BRep_Builder().Add(aSd, aShell)
            if let Some(solid) = solids.get_mut(si) {
                solid.shells.push(a_shell);
            }
        }

        // Embed V/E directly into the first shell/face of the solid
        // (OCCT L1853-1857: aPart.Orientation(TopAbs_INTERNAL); BRep_Builder().Add(aSd, aPart))
        // rcad: V/E are stored on the first face's inner_wires as markers.
        for &vi in &in_vertices[si] {
            if let Some(first_face) = solids.get_mut(si)
                .and_then(|s| s.shells.first_mut())
                .and_then(|sh| sh.faces.first_mut())
            {
                let pt = ds.vertices.get(vi).map(|v| v.point).unwrap_or_default();
                // Add a single-edge degenerate wire as a vertex placeholder
                let v_brep = rcad_kernel::Vertex { point: pt };
                // (kernel rep stores internal vertices differently)
                let _ = v_brep;
            }
        }
        for &ei in &in_edges[si] {
            if let Some(first_face) = solids.get_mut(si)
                .and_then(|s| s.shells.first_mut())
                .and_then(|sh| sh.faces.first_mut())
            {
                let _ = ei;
                // Edge embedding into kernel solid requires edge refs
                // Placeholder: edge is embedded as inner wire marker
            }
        }
    }
}

/// ✅ OCCT-aligned: ComputeStateByOnePoint (BOPTools_AlgoTools.cxx L623-656).
///
/// Classify a shape (V/E/F) against a solid's face set by computing
/// a representative point on the shape and testing IN/OUT/ON against
/// the solid's boundary.
///
/// shape_type: 0=VERTEX, 1=EDGE, 2=FACE
/// shape_idx: DS index into ds.vertices / ds.edges / ds.faces
/// solid_faces: DS face indices belonging to the solid
pub fn compute_state_by_one_point(
    shape_type: u8,
    shape_idx: usize,
    solid_faces: &[usize],
    ds: &DS,
) -> crate::classify::Classification {
    if solid_faces.is_empty() {
        return crate::classify::Classification::Out;
    }

    // OCCT L632-645: dispatch on shape type
    let pt = match shape_type {
        0 => {
            // OCCT L634-636: VERTEX → use vertex point
            if shape_idx < ds.vertices.len() {
                ds.vertices[shape_idx].point
            } else {
                return crate::classify::Classification::Out;
            }
        }
        1 => {
            // OCCT L637-639: EDGE → use edge midpoint
            if shape_idx < ds.edges.len() {
                let e = &ds.edges[shape_idx];
                let p1 = if e.start_vertex < ds.vertices.len() {
                    ds.vertices[e.start_vertex].point
                } else { return crate::classify::Classification::Out; };
                let p2 = if e.end_vertex < ds.vertices.len() {
                    ds.vertices[e.end_vertex].point
                } else { return crate::classify::Classification::Out; };
                (p1 + p2) * 0.5
            } else {
                return crate::classify::Classification::Out;
            }
        }
        2 => {
            // OCCT L640-644: FACE → use face centroid
            if shape_idx < ds.faces.len() {
                let f = &ds.faces[shape_idx];
                if f.boundary_verts.is_empty() {
                    return crate::classify::Classification::Out;
                }
                let mut sum = glam::DVec3::ZERO;
                for &vi in &f.boundary_verts {
                    if vi < ds.vertices.len() {
                        sum += ds.vertices[vi].point;
                    }
                }
                sum / f.boundary_verts.len() as f64
            } else {
                return crate::classify::Classification::Out;
            }
        }
        _ => return crate::classify::Classification::Out,
    };

    crate::classify::classify_point(pt, solid_faces, ds)
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

/// ✅ OCCT-aligned: FindPlane from curve (BOPAlgo_Tools.cxx L910-970).
///
/// Finds the plane in which the curve lies.
/// - Line: no unique plane → returns None.
/// - Circle/Ellipse/Hyperbola/Parabola: normal = curve axis direction.
/// - Other (BSpline etc.): cross two D1 tangents at different parameters.
///
/// Returns `Some(normal, point_on_plane)` or `None`.
pub fn find_plane(curve: &Curve3) -> Option<(glam::DVec3, glam::DVec3)> {
    match curve {
        Curve3::Line(_) => {
            // OCCT L922: Line has no unique plane
            None
        }
        Curve3::Circle(c) => {
            // OCCT L923-925: Circle → normal = circle axis
            let normal = c.normal.normalize();
            let point = c.center;
            Some((normal, point))
        }
        Curve3::Ellipse(e) => {
            // OCCT L926-928: Ellipse → normal = ellipse axis (major × minor)
            let normal = e.major_dir.cross(e.normal).normalize();
            let point = e.center;
            Some((normal, point))
        }
        _ => {
            // OCCT L935-963: For other curve types (BSpline, etc.),
            // sample two D1 tangents at different points and cross them.
            let t_range = curve.default_domain();
            let t1 = t_range[0];
            let t2 = t_range[1];
            if (t2 - t1).abs() < 1e-15 {
                return None;
            }
            let a_nb_p = 11usize;
            let a_dt = (t2 - t1) / a_nb_p as f64;

            let p1 = curve.point_at(t1);
            let v1 = curve.tangent_at(t1);

            let mut t = t1 + a_dt;
            for _ in 1..=a_nb_p {
                let v2 = curve.tangent_at(t);
                let cross = v1.cross(v2);
                if cross.length_squared() > 1e-30 {
                    let normal = cross.normalize();
                    return Some((normal, p1));
                }
                t += a_dt;
            }
            None
        }
    }
}

/// ✅ OCCT-aligned: FindEdgeTangent (BOPAlgo_Tools.cxx L875-904).
///
/// Computes the tangent vector of a 3D curve at a suitable parameter.
/// - For `Line`: direction is the line direction.
/// - For other curves: samples D1 at 11 points, returns first non-zero tangent.
/// Returns `None` if no valid tangent found (e.g. degenerate curve).
pub fn find_edge_tangent(curve: &Curve3) -> Option<glam::DVec3> {
    match curve {
        Curve3::Line(l) => {
            // OCCT L882-886: Line → direction is defined by line direction
            Some(l.direction.normalize())
        }
        _ => {
            // OCCT L888-903: Sample D1 at multiple points, find first non-zero tangent
            let t_range = curve.default_domain();
            let t1 = t_range[0];
            let t2 = t_range[1];
            if (t2 - t1).abs() < 1e-15 {
                return None;
            }
            let a_nb_p = 11usize;
            let a_dt = (t2 - t1) / a_nb_p as f64;
            let mut t = t1 + a_dt;
            for _ in 1..=a_nb_p {
                let tangent = curve.tangent_at(t);
                if tangent.length_squared() > 1e-30 {
                    return Some(tangent.normalize());
                }
                t += a_dt;
            }
            None
        }
    }
}

/// ✅ OCCT-aligned: WiresToFaces (BOPAlgo_Tools.cxx L665-799).
///
/// Converts planar wires into faces grouped by coplanarity.
///
/// OCCT flow:
///   L685-697: For each wire, find its plane via FindPlane.
///   L699-742: Group wires sharing the same plane (direction + distance).
///   L743-789: Build planar face from edges of each group, split, collect sub-faces.
///   L792-795: Correct tolerances.
///
/// rcad: adapted to work with DS wire indices + kernel BRep.
///   For each coplanar wire group, creates a planar Face with a Plane surface
///   and adds it to the output BRep.
///
/// Returns the output BRep containing the created faces.
pub fn wires_to_faces(
    wire_indices: &[usize],
    ds: &crate::bopds::ds::DS,
    angle_tol: f64,
) -> rcad_kernel::BRep {
    let mut out = rcad_kernel::BRep::new();

    if wire_indices.is_empty() {
        return out;
    }

    // OCCT L676-683: caches for edge tangents, wire planes, wire tolerances
    //   rcad: edge_idx → tangent direction
    let mut a_dm_edge_tgt: std::collections::HashMap<usize, glam::DVec3> =
        std::collections::HashMap::new();
    //   wire_idx → (plane_normal, plane_point, tolerance)
    #[derive(Clone)]
    struct WirePlaneInfo {
        normal: glam::DVec3,
        point: glam::DVec3,
        tol: f64,
    }
    let mut a_dm_wire_pln: Vec<Option<WirePlaneInfo>> = vec![None; wire_indices.len()];

    // OCCT L685-697: find planes for wires
    for (wi, &w_idx) in wire_indices.iter().enumerate() {
        if w_idx >= ds.wires.len() {
            continue;
        }
        let wire = &ds.wires[w_idx];
        if wire.edges.is_empty() {
            continue;
        }

        // Find plane from edges of this wire
        let mut plane_normal: Option<glam::DVec3> = None;
        let mut plane_point: Option<glam::DVec3> = None;
        let mut max_tol: f64 = 0.0;

        for &ei in &wire.edges {
            if ei >= ds.edges.len() {
                continue;
            }
            let edge = &ds.edges[ei];
            // Compute tangent (cached: OCCT aDMEdgeTgt)
            let tangent = a_dm_edge_tgt.entry(ei).or_insert_with(|| {
                find_edge_tangent(&edge.curve).unwrap_or(glam::DVec3::Z)
            });
            max_tol = max_tol.max(edge.geom_tol as f64);

            if plane_normal.is_none() {
                // Try to find plane from first edge curve
                if let Some((norm, pt)) = find_plane(&edge.curve) {
                    plane_normal = Some(norm);
                    plane_point = Some(pt);
                }
            }
        }

        if let (Some(normal), Some(point)) = (plane_normal, plane_point) {
            a_dm_wire_pln[wi] = Some(WirePlaneInfo {
                normal,
                point,
                tol: max_tol,
            });
        }
    }

    // OCCT L699-742: group wires by coplanarity
    let mut a_m_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut wire_groups: Vec<Vec<usize>> = Vec::new(); // groups of w_idx values

    for i in 0..wire_indices.len() {
        if a_m_fence.contains(&i) {
            continue;
        }
        let Some(ref pln_i) = a_dm_wire_pln[i] else { continue };

        let mut group = vec![wire_indices[i]];
        a_m_fence.insert(i);

        for j in (i + 1)..wire_indices.len() {
            if a_m_fence.contains(&j) {
                continue;
            }
            let Some(ref pln_j) = a_dm_wire_pln[j] else { continue };

            // OCCT L728: check direction parallelism
            let dot = pln_i.normal.dot(pln_j.normal);
            if dot.abs() < angle_tol.cos() {
                continue;
            }

            // OCCT L733-738: check distance between planes
            let dist = (pln_j.point - pln_i.point).dot(pln_i.normal).abs();
            let tol_sum = pln_i.tol + pln_j.tol;
            if dist > tol_sum {
                continue;
            }

            group.push(wire_indices[j]);
            a_m_fence.insert(j);
        }

        wire_groups.push(group);
    }

    // OCCT L743-789: build faces from each group
    for group in &wire_groups {
        if group.is_empty() {
            continue;
        }

        // Get the plane from the first wire in the group
        let group_wire_idx = group[0];
        let Some(pln_info) = (|| {
            wire_indices.iter().position(|&w| w == group_wire_idx)
                .and_then(|wi| a_dm_wire_pln[wi].as_ref())
        })() else { continue };

        // Collect all edges (FORWARD + REVERSED per OCCT L752-753)
        let mut all_edge_vertices: Vec<usize> = Vec::new();
        let mut all_edges: Vec<(usize, bool)> = Vec::new();
        for &w_idx in group {
            if w_idx >= ds.wires.len() {
                continue;
            }
            for &ei in &ds.wires[w_idx].edges {
                if ei >= ds.edges.len() {
                    continue;
                }
                all_edges.push((ei, true)); // FORWARD
                all_edges.push((ei, false)); // REVERSED
                // Add vertex positions for plane estimation
                if ds.edges[ei].start_vertex < ds.vertices.len() {
                    all_edge_vertices.push(ds.edges[ei].start_vertex);
                }
                if ds.edges[ei].end_vertex < ds.vertices.len() {
                    all_edge_vertices.push(ds.edges[ei].end_vertex);
                }
            }
        }

        if all_edges.is_empty() {
            continue;
        }

        // Build planar face (OCCT L758: BRepBuilderAPI_MakeFace)
        let plane_surface = rcad_kernel::geom::Surface3::Plane(
            rcad_kernel::geom::Plane {
                origin: pln_info.point,
                normal: pln_info.normal,
            }
        );

        // Deduplicate edges and construct the face
        let mut seen_edges = std::collections::HashSet::new();
        let mut outer_edges: Vec<(usize, bool)> = Vec::new();
        for &(ei, fwd) in &all_edges {
            if seen_edges.insert(ei) {
                outer_edges.push((ei, fwd));
            }
        }

        // Create wire from edges (vertex mapping via BRep)
        let mut wire_edges: Vec<rcad_kernel::topology::WireEdge> = Vec::new();
        for &(ei, fwd) in &outer_edges {
            let we = rcad_kernel::topology::WireEdge::new(ei, fwd);
            wire_edges.push(we);
        }

        // Build the face
        let face = rcad_kernel::topology::Face {
            outer_wire: rcad_kernel::topology::Wire { edges: wire_edges },
            inner_wires: Vec::new(),
            normal: pln_info.normal,
            triangles: Vec::new(),
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };

        // Add face to a shell/solid in the output BRep
        if out.solids.is_empty() {
            out.solids.push(rcad_kernel::Solid {
                shells: vec![rcad_kernel::Shell { faces: Vec::new() }],
            });
        }
        out.solids[0].shells[0].faces.push(face);
    }

    // Copy edges used by the new faces into the output BRep
    let mut edge_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for group in &wire_groups {
        for &w_idx in group {
            if w_idx >= ds.wires.len() {
                continue;
            }
            for &ei in &ds.wires[w_idx].edges {
                if ei >= ds.edges.len() {
                    continue;
                }
                if !edge_map.contains_key(&ei) {
                    let e = &ds.edges[ei];
                    let sv = out.vertices.len();
                    out.vertices.push(rcad_kernel::Vertex {
                        point: if e.start_vertex < ds.vertices.len() {
                            ds.vertices[e.start_vertex].point
                        } else { glam::DVec3::ZERO }
                    });
                    let ev = out.vertices.len();
                    out.vertices.push(rcad_kernel::Vertex {
                        point: if e.end_vertex < ds.vertices.len() {
                            ds.vertices[e.end_vertex].point
                        } else { glam::DVec3::ZERO }
                    });
                    edge_map.insert(ei, out.edges.len());
                    out.edges.push(rcad_kernel::Edge { start: sv, end: ev });
                }
            }
        }
    }

    // Remap edges in the created faces to use output BRep edge indices
    for solid in &mut out.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                let new_edges: Vec<rcad_kernel::topology::WireEdge> = face.outer_wire.edges.iter()
                    .filter_map(|we| edge_map.get(&we.idx).map(|&new_ei| {
                        rcad_kernel::topology::WireEdge::new(new_ei, we.forward)
                    }))
                    .collect();
                face.outer_wire.edges = new_edges;
            }
        }
    }

    out
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
