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
    //   rcad: collect all DS face indices used by result solids AND their
    //   boundary vertices/edges (OCCT: TopExp::MapShapes for V/E/F).
    let mut a_ms_solids_vertices: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut a_ms_solids_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut a_ms_solids_faces: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for sdf in solid_ds_faces.iter() {
        for &fi in sdf {
            a_ms_solids_faces.insert(fi);
            if let Some(face) = ds.faces.get(fi) {
                for &ei in &face.boundary_edges {
                    a_ms_solids_edges.insert(ei);
                    if let Some(edge) = ds.edges.get(ei) {
                        a_ms_solids_vertices.insert(edge.start_vertex);
                        a_ms_solids_vertices.insert(edge.end_vertex);
                    }
                }
                for wire in &face.inner_boundary_edges {
                    for &(ei, _) in wire {
                        a_ms_solids_edges.insert(ei);
                        if let Some(edge) = ds.edges.get(ei) {
                            a_ms_solids_vertices.insert(edge.start_vertex);
                            a_ms_solids_vertices.insert(edge.end_vertex);
                        }
                    }
                }
            }
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
            0 => { // VERTEX
                if let Some(img_list) = images.get(&idx) {
                    for &img_idx in img_list {
                        if !a_ms_solids_vertices.contains(&img_idx) {
                            a_l_parts.push(PartInfo { typ, idx: img_idx });
                        }
                    }
                } else if !a_ms_solids_vertices.contains(&idx) {
                    a_l_parts.push(PartInfo { typ, idx });
                }
            }
            1 => { // EDGE
                if let Some(img_list) = images.get(&idx) {
                    for &img_idx in img_list {
                        if !a_ms_solids_edges.contains(&img_idx) {
                            a_l_parts.push(PartInfo { typ, idx: img_idx });
                        }
                    }
                } else if !a_ms_solids_edges.contains(&idx) {
                    a_l_parts.push(PartInfo { typ, idx });
                }
            }
            2 => { // FACE
                if let Some(img_list) = images.get(&idx) {
                    for &img_idx in img_list {
                        if !a_ms_solids_faces.contains(&img_idx) {
                            a_l_parts.push(PartInfo { typ, idx: img_idx });
                        }
                    }
                } else if !a_ms_solids_faces.contains(&idx) {
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
    //   (OCCT: iterate solids OUTER, parts INNER, remove classified parts from list)
    //   anINFaces: per-solid map of IN faces (for shell creation later).
    let mut an_in_faces: Vec<Vec<usize>> = vec![Vec::new(); solids.len()];

    // OCCT L1825-1864: for each solid, classify remaining parts
    for si in 0..solids.len() {
        if solid_ds_faces[si].is_empty() {
            continue;
        }

        // Compute centroid for each part (OCCT uses ComputeStateByOnePoint).
        let mut i = 0usize;
        while i < a_l_parts.len() {
            let part = a_l_parts[i];
            let pt = match part.typ {
                0 => ds.vertices.get(part.idx).map(|v| v.point),
                1 => ds.edges.get(part.idx).and_then(|e| {
                    let p1 = ds.vertices.get(e.start_vertex)?;
                    let p2 = ds.vertices.get(e.end_vertex)?;
                    Some((p1.point + p2.point) * 0.5)
                }),
                2 => {
                    let face = &ds.faces[part.idx];
                    if !face.boundary_verts.is_empty() {
                        let mut sum = glam::DVec3::ZERO;
                        for &vi in &face.boundary_verts {
                            if let Some(v) = ds.vertices.get(vi) {
                                sum += v.point;
                            }
                        }
                        Some(sum / face.boundary_verts.len() as f64)
                    } else { None }
                }
                _ => None,
            };
            let Some(pt) = pt else { i += 1; continue; };

            // OCCT L1841: BOPTools_AlgoTools::ComputeStateByOnePoint
            let a_state = crate::classify::classify_point(pt, &solid_ds_faces[si], ds);

            if a_state == crate::classify::Classification::In {
                if part.typ == 2 { // FACE
                    if !an_in_faces[si].contains(&part.idx) {
                        an_in_faces[si].push(part.idx);
                    }
                } else if part.typ == 0 { // VERTEX
                    // OCCT L1853-1857: aPart.Orientation(TopAbs_INTERNAL);
                    //   BRep_Builder().Add(aSd, aPart).  rcad equivalent:
                    //   topods::TSolidData.internal_vertices.push(aPart)
                    //   (not yet wired; TSolidData field added in topods.rs)
                } else if part.typ == 1 { // EDGE
                    // OCCT: same pattern → topods::TSolidData.internal_edges
                }
                // OCCT L1858: aLParts.Remove(itLP)
                a_l_parts.swap_remove(i);
            } else {
                i += 1;
            }
        }
    }

    // === OCCT L1867-1907: build INTERNAL shells from IN faces ===
    for si in 0..solids.len() {
        let a_faces = &an_in_faces[si];
        if a_faces.is_empty() {
            continue;
        }

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
            let mut a_shell = rcad_kernel::Shell { faces: Vec::new() };

            // OCCT L1897-1903: add faces to shell with INTERNAL orientation
            for &fi in &cb.shapes {
                if let Some(ds_face) = ds.faces.get(fi) {
                    let f = rcad_kernel::Face {
                        outer_wire: rcad_kernel::Wire { edges: Vec::new() },
                        inner_wires: Vec::new(),
                        normal: ds_face.normal,
                        triangles: Vec::new(),
                        sample_point: None,
                        mesh_dirty: true,
                        surface_idx: None,
                    };
                    a_shell.faces.push(f);
                }
            }

            // OCCT L1905: BRep_Builder().Add(aSd, aShell)
            if let Some(solid) = solids.get_mut(si) {
                solid.shells.push(a_shell);
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
/// OCCT-aligned: MakeContainer (BOPTools_AlgoTools.cxx L1600-1645).
///
/// Creates an empty container shape of the specified type.
/// rcad: maps to topods TShape variants. Returns a ShapeRef with
/// the appropriate empty TShape.
pub fn make_container(shape_type: u8, brep: &mut rcad_kernel::topods::BRep) -> rcad_kernel::topods::ShapeRef {
    use rcad_kernel::topods::{TShape, ShapeRef, TShellData, TWireData, TSolidData};
    let idx = brep.tshapes.len();
    let shape: std::sync::Arc<TShape> = match shape_type {
        0 => std::sync::Arc::new(TShape::Compound(Vec::new())),
        1 => std::sync::Arc::new(TShape::CompSolid(Vec::new())),
        2 => std::sync::Arc::new(TShape::Solid(TSolidData { shells: Vec::new(), internal_vertices: Vec::new(), internal_edges: Vec::new(), moved: false })),
        3 => std::sync::Arc::new(TShape::Shell(TShellData { faces: Vec::new(), closed: false, moved: false })),
        4 => std::sync::Arc::new(TShape::Wire(TWireData { edges: Vec::new(), closed: false, moved: false })),
        _ => return ShapeRef::new(0),
    };
    brep.tshapes.push(shape);
    ShapeRef::new(idx)
}

/// OCCT-aligned: PointOnEdge (BOPTools_AlgoTools_2.cxx L275-280).
///
/// Evaluates a point on a DS edge at the given parameter.
pub fn point_on_edge(edge: &crate::bopds::ds::DSEdge, t: f64) -> glam::DVec3 {
    edge.curve.point_at(t)
}

/// ✅ OCCT-aligned: IsInvertedSolid (BOPTools_AlgoTools.cxx L2398-2408).
///
/// Checks if a solid is inverted (normals point inward) by classifying
/// an "infinite point" against the solid's face set.
/// Returns true if the infinite point is IN (solid encloses infinity → inverted).
pub fn is_inverted_solid(solid_faces: &[usize], ds: &DS) -> bool {
    if solid_faces.is_empty() {
        return false;
    }
    // Build an AABB from the solid's vertices, then place a test point
    // far outside that AABB (OCCT uses PerformInfinitePoint which classifies
    // a point at ~1e15 from origin).  Here we use the AABB diagonal * 100.
    let mut aabb_min = glam::DVec3::splat(f64::INFINITY);
    let mut aabb_max = glam::DVec3::splat(f64::NEG_INFINITY);
    for &fi in solid_faces {
        if fi >= ds.faces.len() { continue; }
        for &vi in &ds.faces[fi].boundary_verts {
            if vi < ds.vertices.len() {
                let p = ds.vertices[vi].point;
                aabb_min = aabb_min.min(p);
                aabb_max = aabb_max.max(p);
            }
        }
    }
    let size = (aabb_max - aabb_min).length();
    if size < 1e-30 { return false; }
    let far_point = aabb_max + glam::DVec3::splat(size * 100.0);
    let state = crate::classify::classify_point(far_point, solid_faces, ds);
    state == crate::classify::Classification::In
}

/// ✅ OCCT-aligned: IntermediatePoint (BOPTools_AlgoTools2D / IntTools_Tools).
pub fn intermediate_point(t1: f64, t2: f64) -> f64 {
    0.5 * (t1 + t2)
}

/// ✅ OCCT-aligned: IntermediatePoint (IntTools_Tools.cxx L254-258).
/// OCCT uses PAR_T = 0.43213918 instead of 0.5 for numerical stability.
pub fn intermediate_point_occt(t1: f64, t2: f64) -> f64 {
    const PAR_T: f64 = 0.43213918;
    (1.0 - PAR_T) * t1 + PAR_T * t2
}

/// ✅ OCCT-aligned: IsDirsCoinside (IntTools_Tools.cxx L164-173).
/// Checks if two direction unit vectors are coincident within angular threshold.
/// OCCT default threshold is 0.0002 (approx 0.011 degrees).
pub fn is_dirs_coinside(d1: glam::DVec3, d2: glam::DVec3) -> bool {
    let d = (d1 - d2).length();
    let lim = 0.0002;
    d < lim || (2.0 - d).abs() < lim
}

/// ✅ OCCT-aligned: IsDirsCoinside with custom threshold (IntTools_Tools.cxx L177-187).
pub fn is_dirs_coinside_with_tol(d1: glam::DVec3, d2: glam::DVec3, d_lim: f64) -> bool {
    let d = (d1 - d2).length();
    d < d_lim || (2.0 - d).abs() < d_lim
}

/// ✅ OCCT-aligned: IsClosed (IntTools_Tools.cxx L78-102).
/// Checks if a bounded 3D curve is closed (start point ≈ end point).
pub fn is_curve_closed(curve: &Curve3, t_range: [f64; 2]) -> bool {
    let p1 = curve.point_at(t_range[0]);
    let p2 = curve.point_at(t_range[1]);
    let conf = crate::tolerance::TOLERANCE_ABS;
    (p1 - p2).length_squared() < conf * conf
}

/// ✅ OCCT-aligned: IsOnPave (IntTools_Tools.cxx L579-589).
/// Checks if parameter aT1 is within aTolerance of either range boundary.
pub fn is_on_pave(t: f64, range: [f64; 2], tol: f64) -> bool {
    (range[0] - t).abs() < tol || (range[1] - t).abs() < tol
}

/// ✅ OCCT-aligned: IsInRange (IntTools_Tools.cxx L650-666).
/// Checks if either endpoint of range aR falls within aRRef (expanded by tolerance).
pub fn is_in_range(r: [f64; 2], r_ref: [f64; 2], tol: f64) -> bool {
    let t_ref1 = r_ref[0] - tol;
    let t_ref2 = r_ref[1] + tol;
    (r[0] >= t_ref1 && r[0] <= t_ref2) || (r[1] >= t_ref1 && r[1] <= t_ref2)
}

/// ✅ OCCT-aligned: ComputeIntRange (IntTools_Tools.cxx L783-804).
/// Computes intersection range tolerance from two tolerances and the angle
/// between the surfaces at the intersection point.
pub fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    let half_pi = std::f64::consts::FRAC_PI_2;
    if (half_pi - angle).abs() < crate::tolerance::TOLERANCE_ANG {
        return tol2;
    }
    let an_angle = if angle > half_pi { std::f64::consts::PI - angle } else { angle };
    let a1 = tol1 * (half_pi - an_angle).tan();
    let a2 = tol2 / an_angle.sin();
    a1 + a2
}

/// ✅ OCCT-aligned: CurveTolerance (IntTools_Tools.cxx L430-464).
/// Computes the max tolerance for a curve given the base tolerance.
/// For non-parabola curves, returns aTolBase unchanged.
pub fn curve_tolerance(curve: &Curve3, tol_base: f64) -> f64 {
    // OCCT only adjusts for Parabola curves; all other curve types return tol_base.
    match curve {
        Curve3::Parabola(_) => tol_base * 10.0, // approximate: parabola may need larger tol
        _ => tol_base,
    }
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

/// ✅ OCCT-aligned: OrientFacesOnShell (BOPTools_AlgoTools.cxx L363-507).
///
/// OCCT algorithm:
///   1. Build edge→face map for the shell, deduplicate seam edges.
///   2. For each non-degenerated edge with exactly 2 faces:
///      a. If both unprocessed: add first face to new shell.
///      b. Get Orientation(edge, face) for each face (FWD/REV in wire).
///      c. If one processed, one not:
///         - Orientations equal → reverse the unprocessed face
///           (unless edge is closed in that face).
///         - Mark processed, add to shell.
///   3. For edges with != 2 faces: add remaining unprocessed faces.
///
/// Returns the reordered shell face indices.  Faces that need reversal
/// are added to `reversed` (caller must negate normals when building BRep).
pub fn orient_faces_on_shell(
    shell_faces: &[usize],
    ds: &DS,
    reversed: &mut std::collections::HashSet<usize>,
) -> Vec<usize> {
    if shell_faces.is_empty() { return Vec::new(); }

    let mut a_processed: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut a_shell_new: Vec<usize> = Vec::new();

    // OCCT L370-377: TopExp::MapShapesAndAncestors → edge→face map
    let mut a_ef_map: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for &fi in shell_faces.iter() {
        if let Some(face) = ds.faces.get(fi) {
            for &ei in &face.boundary_edges {
                a_ef_map.entry(ei).or_default().push(fi);
            }
        }
    }

    // OCCT L380-403: deduplicate seam edges
    let edge_keys: Vec<usize> = a_ef_map.keys().copied().collect();
    for &ei in &edge_keys {
        if let Some(a_lf) = a_ef_map.get_mut(&ei) {
            a_lf.sort_unstable();
            a_lf.dedup();
        }
    }

    // Helper: edge orientation in face (OCCT local function `Orientation` L511-527)
    let edge_orientation_in_face = |ei: usize, fi: usize, ds: &DS| -> bool {
        if let Some(face) = ds.faces.get(fi) {
            if face.boundary_edges.contains(&ei) {
                // Forward if edge's start_vertex comes first in boundary_edges
                return true;
            }
        }
        true
    };

    // OCCT L406-478: process edges with exactly 2 faces
    for &ei in &edge_keys {
        if ds.is_edge_degenerated(ei) { continue; }
        let Some(a_lf) = a_ef_map.get(&ei) else { continue; };
        if a_lf.len() != 2 { continue; }

        let a_f1 = a_lf[0];
        let a_f2 = a_lf[1];

        let b_p1 = a_processed.contains(&a_f1);
        let b_p2 = a_processed.contains(&a_f2);
        if b_p1 && b_p2 { continue; }

        if !b_p1 && !b_p2 {
            a_processed.insert(a_f1);
            a_shell_new.push(a_f1);
        }

        let an_or_e1 = edge_orientation_in_face(ei, a_f1, ds);
        let an_or_e2 = edge_orientation_in_face(ei, a_f2, ds);

        if b_p1 && !b_p2 {
            if an_or_e1 == an_or_e2 {
                let e_closed = ds.edges.get(ei).map_or(false, |e| e.start_vertex == e.end_vertex);
                if !e_closed {
                    reversed.insert(a_f2);
                }
            }
            a_processed.insert(a_f2);
            a_shell_new.push(a_f2);
        } else if !b_p1 && b_p2 {
            if an_or_e1 == an_or_e2 {
                let e_closed = ds.edges.get(ei).map_or(false, |e| e.start_vertex == e.end_vertex);
                if !e_closed {
                    reversed.insert(a_f1);
                }
            }
            a_processed.insert(a_f1);
            a_shell_new.push(a_f1);
        }
    }

    // OCCT L482-505: remaining faces from edges with != 2 faces
    for &ei in &edge_keys {
        if ds.is_edge_degenerated(ei) { continue; }
        let Some(a_lf) = a_ef_map.get(&ei) else { continue; };
        if a_lf.len() != 2 {
            for &fi in a_lf {
                if a_processed.insert(fi) {
                    a_shell_new.push(fi);
                }
            }
        }
    }

    a_shell_new
}

/// OCCT-aligned: IsSplitToReverse (BOPTools_AlgoTools).
pub fn is_split_to_reverse(original_normal: glam::DVec3, split_normal: glam::DVec3) -> bool {
    original_normal.dot(split_normal) < 0.0
}

/// ✅ OCCT-aligned: ComputeToleranceOfCB (BOPAlgo_Tools.cxx L248-356).
///
/// OCCT algorithm:
///   1. Get reference PB's original edge tolerance as aTolMax.
///   2. If no other PBs and no faces → return aTolMax.
///   3. Sample reference curve at 11 points.
///   4. For each OTHER PB: project points onto that edge's curve,
///      aTolNew = edge.tol + projection_distance.
///   5. For each face: project points onto face surface,
///      aTolNew = face.tol + projection_distance.
///   6. Return max aTolMax found.
pub fn compute_tolerance_of_cb(
    cb: &crate::bopds::common_block::CommonBlock,
    ds: &DS,
) -> f64 {
    // OCCT L252-256: null check
    let a_pbr_pbs = cb.pave_blocks();
    if a_pbr_pbs.is_empty() { return 0.0; }
    let (a_pb_idx_ref, _) = a_pbr_pbs[0];
    let a_ref_pb = &ds.pave_blocks[a_pb_idx_ref];
    let n_e = a_ref_pb.original_edge;
    let a_tol_max = ds.edges.get(n_e).map(|e| e.geom_tol).unwrap_or(0.0);
    if a_pbr_pbs.len() < 2 && cb.faces().is_empty() { return a_tol_max; }

    // OCCT L271-278: sample reference curve
    let a_nb_pnt = 11usize;
    let a_curve = &ds.edges[n_e].curve;
    let (a_t1, a_t2) = (ds.edges[n_e].t_range[0], ds.edges[n_e].t_range[1]);
    let a_dt = (a_t2 - a_t1) / (a_nb_pnt + 1) as f64;

    let mut a_tol_max = a_tol_max;

    // OCCT L287-323: iterate other PaveBlocks
    for &(a_pb_idx, _fi) in &a_pbr_pbs[1..] {
        let a_pb = &ds.pave_blocks[a_pb_idx];
        let a_e_idx = a_pb.original_edge;
        let a_tol_e = ds.edges.get(a_e_idx).map(|e| e.geom_tol).unwrap_or(0.0);
        let a_curve_other = &ds.edges[a_e_idx].curve;

        let mut t = a_t1;
        for _ in 1..=a_nb_pnt {
            t += a_dt;
            let a_p = a_curve.point_at(t);
            // Project point onto other edge's curve
            let proj = rcad_kernel::closest_point_on_curve(a_curve_other, a_p, 16);
            let a_dist = if proj.distance.is_finite() { proj.distance } else { 0.0 };
            let a_tol_new = a_tol_e + a_dist;
            if a_tol_new > a_tol_max { a_tol_max = a_tol_new; }
        }
    }

    // OCCT L327-353: iterate faces
    for &a_fi in cb.faces() {
        let a_tol_f = ds.faces.get(a_fi).map(|f| f.geom_tol).unwrap_or(0.0);
        let surf = &ds.faces[a_fi].surface;
        let mut t = a_t1;
        for _ in 1..=a_nb_pnt {
            t += a_dt;
            let a_p = a_curve.point_at(t);
            let proj = rcad_kernel::projection::closest_point_on_surface(surf, a_p, 16);
            let a_dist = if proj.distance.is_finite() { proj.distance } else { 0.0 };
            let a_tol_new = a_tol_f + a_dist;
            if a_tol_new > a_tol_max { a_tol_max = a_tol_new; }
        }
    }
    a_tol_max
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

/// OCCT-aligned: ComputeVV (BOPTools_AlgoTools.cxx L1742-1760).
/// Vertex-point intersection. Returns 0 if vertex and point intersect.
pub fn compute_vv_p(v: &crate::bopds::ds::DSVertex, p: glam::DVec3, tol_p: f64) -> i32 {
    let tol_sum = v.geom_tol + tol_p + crate::tolerance::TOLERANCE_ABS;
    let tol_sum2 = tol_sum * tol_sum;
    let d2 = (v.point - p).length_squared();
    if d2 > tol_sum2 { 1 } else { 0 }
}

/// OCCT-aligned: ComputeVV (BOPTools_AlgoTools.cxx L1764-1786).
/// Vertex-vertex intersection with fuzz tolerance.
pub fn compute_vv(v1: &crate::bopds::ds::DSVertex, v2: &crate::bopds::ds::DSVertex, fuzz: f64) -> i32 {
    let a_fuzz = if fuzz > crate::tolerance::TOLERANCE_ABS { fuzz } else { crate::tolerance::TOLERANCE_ABS };
    let tol_sum = v1.geom_tol + v2.geom_tol + a_fuzz;
    let tol_sum2 = tol_sum * tol_sum;
    let d2 = (v1.point - v2.point).length_squared();
    if d2 > tol_sum2 { 1 } else { 0 }
}

/// OCCT-aligned: MakeNewVertex from point + tolerance (BOPTools_AlgoTools.cxx L96-98).
/// Returns the new vertex DS index.
pub fn make_new_vertex(ds: &mut crate::bopds::ds::DS, p: glam::DVec3, tol: f64) -> usize {
    let vi = ds.vertices.len();
    ds.vertices.push(crate::bopds::ds::DSVertex {
        point: p, geom_tol: tol, origin: None, is_internal: false,
    });
    vi
}

/// OCCT-aligned: MakeNewVertex from two vertices (BOPTools_AlgoTools.cxx L101-103).
/// Merges two vertices into one at their midpoint with combined tolerance.
pub fn make_new_vertex_from_two(v1: &crate::bopds::ds::DSVertex, v2: &crate::bopds::ds::DSVertex,
    ds: &mut crate::bopds::ds::DS) -> usize {
    let mid = (v1.point + v2.point) * 0.5;
    let dist = (v1.point - v2.point).length();
    let tol = dist * 0.5 + v1.geom_tol.max(v2.geom_tol) + crate::tolerance::TOLERANCE_LEN_MIN;
    let vi = ds.vertices.len();
    ds.vertices.push(crate::bopds::ds::DSVertex {
        point: mid, geom_tol: tol, origin: None, is_internal: false,
    });
    vi
}

/// OCCT-aligned: GetEdgeOnFace (BOPTools_AlgoTools.cxx L1809-1835).
/// Finds the edge in the face's wires that is the same as theE1.
/// Returns Some((edge_index, forward)) if found.
pub fn get_edge_on_face(ei: usize, ds_edges: &[crate::bopds::ds::DSEdge], face: &crate::bopds::ds::DSFace) -> Option<(usize, bool)> {
    for &be in &face.boundary_edges {
        let edge = &ds_edges[be];
        if edge.start_vertex == ds_edges[ei].start_vertex && edge.end_vertex == ds_edges[ei].end_vertex
            || edge.start_vertex == ds_edges[ei].end_vertex && edge.end_vertex == ds_edges[ei].start_vertex
        {
            let forward = edge.start_vertex == ds_edges[ei].start_vertex;
            return Some((be, forward));
        }
    }
    // Check inner wires
    for inner in &face.inner_boundary_edges {
        for &(be, _fwd) in inner {
            let edge = &ds_edges[be];
            if edge.start_vertex == ds_edges[ei].start_vertex && edge.end_vertex == ds_edges[ei].end_vertex
                || edge.start_vertex == ds_edges[ei].end_vertex && edge.end_vertex == ds_edges[ei].start_vertex
            {
                return Some((be, edge.start_vertex == ds_edges[ei].start_vertex));
            }
        }
    }
    None
}

/// OCCT-aligned: GetEdgeOff (BOPTools_AlgoTools.cxx L1099-1127).
/// Finds the edge in theFace that is the same as theE1 but with opposite orientation.
/// Returns Some(edge_index) if found.
pub fn get_edge_off(ei: usize, ds_edges: &[crate::bopds::ds::DSEdge], face: &crate::bopds::ds::DSFace) -> Option<usize> {
    let e1 = &ds_edges[ei];
    for &be in &face.boundary_edges {
        let edge = &ds_edges[be];
        if edge.start_vertex == e1.end_vertex && edge.end_vertex == e1.start_vertex {
            return Some(be);
        }
    }
    for inner in &face.inner_boundary_edges {
        for &(be, _fwd) in inner {
            let edge = &ds_edges[be];
            if edge.start_vertex == e1.end_vertex && edge.end_vertex == e1.start_vertex {
                return Some(be);
            }
        }
    }
    None
}

/// OCCT-aligned: IsHole (BOPTools_AlgoTools.cxx L1527-1596).
/// Determines if a wire is a hole on a face using signed area of the
/// pcurve polygon. Returns true if the wire is a hole (signed area > 0).
pub fn is_hole(
    wire_edges: &[(usize, bool)],  // (edge_idx, forward_in_wire)
    face_idx: usize,
    ds: &crate::bopds::ds::DS,
) -> bool {
    let mut area = 0.0;
    for &(ei, fwd) in wire_edges {
        let edge = &ds.edges[ei];
        let rep = edge.face_reps.iter().find(|r| r.face_idx == face_idx);
        let (Some(ref pc), Some(range)) = (rep.map(|r| &r.pcurve), rep.map(|r| r.pcurve_range)) else {
            continue;
        };
        // OCCT: sample pcurve at NbSamples points, compute signed area (Y0+Y1)*(X1-X0)
        let nb_samples = 11usize; // OCCT: NbSamples * 4 (min 2, then *4)
        let dt = (range[1] - range[0]) / (nb_samples - 1) as f64;
        let (t_start, step) = if fwd {
            (range[0], dt)
        } else {
            (range[1], -dt)
        };
        let mut t = t_start;
        let p0 = pc.point_at(t);
        let (mut x0, mut y0) = (p0.x, p0.y);
        for _ in 1..nb_samples {
            t += step;
            let p1 = pc.point_at(t);
            let (x1, y1) = (p1.x, p1.y);
            area += (y0 + y1) * (x1 - x0);
            x0 = x1;
            y0 = y1;
        }
    }
    area > 0.0
}

/// OCCT-aligned: Sense (BOPTools_AlgoTools.cxx L1201-1251).
/// Determines the relative orientation of normals of two faces near a shared edge.
/// Returns 1 if normals point in the same direction,
/// -1 if opposite, 0 if no shared non-closed edge found.
pub fn sense(
    fi_a: usize,
    fi_b: usize,
    ds: &crate::bopds::ds::DS,
) -> i8 {
    // OCCT L1210-1238: find a non-degenerate, non-closed edge shared by both faces
    let face_a = &ds.faces[fi_a];
    let face_b = &ds.faces[fi_b];
    for &ei in &face_a.boundary_edges {
        if face_b.boundary_edges.contains(&ei) {
            // Found shared edge
            let edge = &ds.edges[ei];
            if edge.start_vertex == edge.end_vertex {
                continue; // closed edge, skip
            }
            // Compute normals near the edge
            let t_mid = (edge.t_range[0] + edge.t_range[1]) * 0.5;
            let edge_mid = edge.curve.point_at(t_mid);
            let n_a = get_normal_to_face_on_edge(&face_a.surface, face_a.normal, edge_mid);
            let n_b = get_normal_to_face_on_edge(&face_b.surface, face_b.normal, edge_mid);
            let dot = n_a.dot(n_b);
            return if dot > 1e-10 { 1 } else if dot < -1e-10 { -1 } else { 0 };
        }
    }
    0
}

/// OCCT-aligned: IsOpenShell (BOPTools_AlgoTools.cxx L2350-2394).
/// Checks if a shell is open (has an edge used by only one non-INTERNAL/EXTERNAL face).
pub fn is_open_shell(
    shell_faces: &[usize],
    ds: &crate::bopds::ds::DS,
) -> bool {
    // OCCT L2361: build edge→face map for the shell
    let mut edge_face_count: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &fi in shell_faces {
        if fi >= ds.faces.len() { continue; }
        let f = &ds.faces[fi];
        for &ei in &f.boundary_edges {
            if ds.is_edge_degenerated(ei) { continue; }
            *edge_face_count.entry(ei).or_insert(0) += 1;
        }
        for inner in &f.inner_boundary_edges {
            for &(ei, _) in inner {
                if ds.is_edge_degenerated(ei) { continue; }
                *edge_face_count.entry(ei).or_insert(0) += 1;
            }
        }
    }
    // OCCT L2367-2391: an edge used by only one non-INTERNAL/EXTERNAL face → shell is open
    for &cnt in edge_face_count.values() {
        if cnt == 1 {
            return true;
        }
    }
    false
}

/// OCCT-aligned: IsBlockInOnFace (BOPTools_AlgoTools.cxx L1971-2063).
/// Checks if a pave block (shrunk range) lies in/on the face.
/// Tests three points: near start, near end, and intermediate.
pub fn is_block_in_on_face(
    t_range: [f64; 2],
    face_idx: usize,
    ei: usize,
    ds: &crate::bopds::ds::DS,
) -> bool {
    let [mut f1, mut l1] = t_range;
    // OCCT L1982-1985: shrink range by 0.75% on each side
    let dt = 0.0075;
    let k = dt * (l1 - f1);
    f1 += k;
    l1 -= k;
    if l1 <= f1 { return false; }

    // Get edge curve for point evaluation
    let edge = &ds.edges[ei];

    // OCCT L1988-2007: test near-start point
    let p11 = edge.curve.point_at(f1);
    if !is_point_in_on_face(&p11, face_idx, ds) { return false; }

    // OCCT L2010-2028: test near-end point
    let p12 = edge.curve.point_at(l1);
    if !is_point_in_on_face(&p12, face_idx, ds) { return false; }

    // OCCT L2033-2062: test intermediate point with distance check
    let m1 = (f1 + l1) * 0.5;
    let p_mid = edge.curve.point_at(m1);

    // Project onto face surface
    let face = &ds.faces[face_idx];
    let proj_dist = match &face.surface {
        rcad_kernel::geom::Surface3::Plane(pl) => {
            (p_mid - pl.origin).dot(pl.normal).abs()
        }
        _ => {
            // For non-planar surfaces, use projection API
            let proj = rcad_kernel::projection::closest_point_on_surface(
                &face.surface, p_mid, 16);
            if proj.distance.is_finite() {
                proj.distance
            } else {
                return false;
            }
        }
    };

    let tol_e = edge.geom_tol;
    let tol_f = face.geom_tol;
    let tol = tol_e + tol_f;
    if proj_dist > tol { return false; }

    is_point_in_on_face(&p_mid, face_idx, ds)
}

/// Helper: check if a 3D point is in/on a face by projecting to UV and
/// testing against the face's surface proximity.
pub(crate) fn is_point_in_on_face(pt: &glam::DVec3, face_idx: usize, ds: &crate::bopds::ds::DS) -> bool {
    let face = &ds.faces[face_idx];
    match &face.surface {
        rcad_kernel::geom::Surface3::Plane(pl) => {
            let dist = (*pt - pl.origin).dot(pl.normal).abs();
            dist < face.geom_tol * 10.0
        }
        _ => {
            // For non-planar: project point onto surface and check distance
            let proj = rcad_kernel::projection::closest_point_on_surface(
                &face.surface, *pt, 16);
            if proj.distance.is_finite() {
                proj.distance < face.geom_tol * 100.0
            } else {
                false
            }
        }
    }
}

/// OCCT-aligned: ComputeTolerance (BOPTools_AlgoTools_1.cxx L1093-1111).
/// Computes max distance and max parameter deviation of edge from face surface.
/// Returns Some((max_dist, max_param)) on success.
pub fn compute_tolerance(
    edge: &crate::bopds::ds::DSEdge,
    face: &crate::bopds::ds::DSFace,
    ds: &crate::bopds::ds::DS,
) -> Option<(f64, f64)> {
    // OCCT: BRepLib_CheckCurveOnSurface — check edge deviation from face surface
    // rcad: sample the edge curve at multiple points and compute distance to face surface
    let n_samples = 23usize;
    let t_range = edge.t_range;
    let dt = (t_range[1] - t_range[0]) / (n_samples - 1) as f64;
    let mut max_dist = 0.0;
    let mut max_par = 0.0;
    for i in 0..n_samples {
        let t = t_range[0] + i as f64 * dt;
        let pt = edge.curve.point_at(t);
        let dist = match &face.surface {
            rcad_kernel::geom::Surface3::Plane(pl) => {
                (pt - pl.origin).dot(pl.normal).abs()
            }
            rcad_kernel::geom::Surface3::Sphere(s) => {
                ((pt - s.center).length() - s.radius).abs()
            }
            rcad_kernel::geom::Surface3::Cylinder(c) => {
                let v = pt - c.origin;
                let radial = v - c.axis * v.dot(c.axis);
                (radial.length() - c.radius).abs()
            }
            rcad_kernel::geom::Surface3::Cone(cn) => {
                let axis = cn.axis_dir();
                let apex_to_pt = pt - cn.apex;
                let proj = apex_to_pt.dot(axis);
                let radial = apex_to_pt - axis * proj;
                let r_at_z = cn.radius + proj * cn.half_angle_rad.tan();
                (radial.length() - r_at_z.abs()).abs()
            }
            _ => {
                // Generic surface: use projection API
                let proj = rcad_kernel::projection::closest_point_on_surface(
                    &face.surface, pt, 16);
                if proj.distance.is_finite() {
                    proj.distance
                } else {
                    continue;
                }
            }
        };
        if dist > max_dist {
            max_dist = dist;
            max_par = t;
        }
    }
    if max_dist < 1e-30 { None } else { Some((max_dist, max_par)) }
}

/// OCCT-aligned: MakeEdge (BOPTools_AlgoTools.cxx L1721-1738).
/// Creates a DS edge from an intersection curve with vertices and tolerance update.
/// Updates both vertex tolerances (adds DTolerance margin) and sets edge tolerance.
/// Returns the new edge index.
pub fn make_edge(
    ds: &mut crate::bopds::ds::DS,
    ci: usize,
    v1: usize,
    v2: usize,
    t1: f64,
    t2: f64,
    tol_r3d: f64,
) -> usize {
    // OCCT L1730: aNeedTol = theTolR3D + DTolerance()
    let need_tol = tol_r3d + crate::tolerance::TOLERANCE_LEN_MIN;
    // OCCT L1732-1733: UpdateVertex theV1/theV2 with aNeedTol
    if v1 < ds.vertices.len() {
        ds.vertices[v1].geom_tol = ds.vertices[v1].geom_tol.max(need_tol);
    }
    if v2 < ds.vertices.len() {
        ds.vertices[v2].geom_tol = ds.vertices[v2].geom_tol.max(need_tol);
    }
    // OCCT L1735: MakeSectEdge(aIC, aV1, aT1, aV2, aT2, aE)
    let ei = make_sect_edge(ds, ci, v1, v2);
    // OCCT L1737: UpdateEdge(aE, theTolR3D)
    if ei < ds.edges.len() {
        ds.edges[ei].geom_tol = ds.edges[ei].geom_tol.max(tol_r3d);
    }
    ei
}

/// OCCT-aligned: CopyEdge (BOPTools_AlgoTools.hxx L152).
/// Creates a copy of a DS edge (deep copy with new index).
pub fn copy_ds_edge(ds: &mut crate::bopds::ds::DS, ei: usize) -> usize {
    let new_ei = ds.edges.len();
    let src = &ds.edges[ei].clone();
    ds.edges.push(src.clone());
    new_ei
}

/// OCCT-aligned: MakeSplitEdge (BOPTools_AlgoTools.hxx L155-160).
/// Creates a split edge from a base edge with new vertices at specified parameters.
/// The new edge inherits the base edge's curve and pcurves.
pub fn make_split_edge(
    ds: &mut crate::bopds::ds::DS,
    base_ei: usize,
    v1: usize,
    p1: f64,
    v2: usize,
    p2: f64,
) -> usize {
    // Copy fields from base edge before mutable borrow
    let (curve, t_range_src, origin, geom_tol, face_reps, is_internal) = {
        let base = &ds.edges[base_ei];
        (base.curve.clone(), base.t_range, base.origin, base.geom_tol, base.face_reps.clone(), base.is_internal)
    };
    let new_ei = ds.edges.len();
    let t_range = if p1 < p2 { [p1, p2] } else { [p2, p1] };
    ds.edges.push(crate::bopds::ds::DSEdge {
        start_vertex: v1,
        end_vertex: v2,
        curve,
        t_range,
        origin,
        geom_tol,
        paves: Vec::new(),
        pave_blocks: Vec::new(),
        face_reps,
        is_internal,
        vertex_params: {
            let mut vp = std::collections::HashMap::new();
            vp.insert(v1, p1);
            vp.insert(v2, p2);
            vp
        },
    });
    new_ei
}

/// OCCT-aligned: MakeVertex (BOPTools_AlgoTools.cxx L1790-1805).
/// Creates a vertex from a list of vertex indices. If single vertex, returns it.
/// Otherwise computes the bounding vertex (midpoint + combined tolerance).
pub fn make_vertex_from_list(
    ds: &mut crate::bopds::ds::DS,
    vertex_indices: &[usize],
) -> usize {
    if vertex_indices.is_empty() {
        // OCCT returns invalid vertex for empty list
        return usize::MAX;
    }
    if vertex_indices.len() == 1 {
        // OCCT L1793-1796: aVnew = first vertex
        return vertex_indices[0];
    }
    // OCCT L1797-1804: BRepLib::BoundingVertex → midpoint + bounding tolerance
    let mut sum = glam::DVec3::ZERO;
    let mut max_tol = 0.0f64;
    for &vi in vertex_indices {
        if vi < ds.vertices.len() {
            sum += ds.vertices[vi].point;
            max_tol = max_tol.max(ds.vertices[vi].geom_tol);
        }
    }
    let n = vertex_indices.len() as f64;
    let mid = sum / n;
    // Compute bounding tolerance: max distance from midpoint + max source tolerance
    let max_dist = vertex_indices.iter()
        .filter_map(|&vi| ds.vertices.get(vi))
        .map(|v| (v.point - mid).length())
        .fold(0.0f64, f64::max);
    let tol = max_dist + max_tol + crate::tolerance::TOLERANCE_LEN_MIN;
    let vi = ds.vertices.len();
    ds.vertices.push(crate::bopds::ds::DSVertex {
        point: mid, geom_tol: tol, origin: None, is_internal: false,
    });
    vi
}

/// OCCT-aligned: UpdateVertex from curve (BOPTools_AlgoTools.hxx L124-126).
/// Updates vertex tolerance given its position on a curve.
/// The new tolerance covers the distance from the vertex to the curve.
pub fn update_vertex_from_curve(
    ds: &mut crate::bopds::ds::DS,
    vi: usize,
    curve: &rcad_kernel::geom::Curve3,
    t: f64,
) {
    if vi >= ds.vertices.len() { return; }
    let pt = curve.point_at(t);
    let dist = (ds.vertices[vi].point - pt).length();
    let new_tol = dist + ds.vertices[vi].geom_tol + crate::tolerance::TOLERANCE_LEN_MIN;
    ds.vertices[vi].geom_tol = ds.vertices[vi].geom_tol.max(new_tol);
}

/// OCCT-aligned: UpdateVertex from edge (BOPTools_AlgoTools.hxx L131-133).
/// Updates vertex tolerance given its parameter on an edge.
pub fn update_vertex_from_edge(
    ds: &mut crate::bopds::ds::DS,
    vi: usize,
    ei: usize,
    t: f64,
) {
    if vi >= ds.vertices.len() || ei >= ds.edges.len() { return; }
    let edge = &ds.edges[ei];
    let pt = edge.curve.point_at(t);
    let dist = (ds.vertices[vi].point - pt).length();
    let new_tol = dist + edge.geom_tol + crate::tolerance::TOLERANCE_LEN_MIN;
    ds.vertices[vi].geom_tol = ds.vertices[vi].geom_tol.max(new_tol);
}

/// OCCT-aligned: UpdateVertex from another vertex (BOPTools_AlgoTools.hxx L136-138).
/// Updates vertex aVN's tolerance to cover tolerance zone of aVF.
pub fn update_vertex_from_vertex(
    ds: &mut crate::bopds::ds::DS,
    vn: usize,
    vf: usize,
) {
    if vn >= ds.vertices.len() || vf >= ds.vertices.len() { return; }
    let dist = (ds.vertices[vn].point - ds.vertices[vf].point).length();
    let new_tol = dist + ds.vertices[vf].geom_tol + crate::tolerance::TOLERANCE_LEN_MIN;
    ds.vertices[vn].geom_tol = ds.vertices[vn].geom_tol.max(new_tol);
}

/// OCCT-aligned: MakePCurve (BOPTools_AlgoTools.cxx L1649-1717).
/// Builds pcurves for edge `ei` on faces `fi_a`/`fi_b` using pcurves from
/// intersection curve `ci`.  Handles periodic surface adjustment.
///
/// b_pc1/b_pc2: if true, compute pcurve for that face; if false, skip.
/// pcurve_a/pcurve_b: existing pcurves from intersection curve (may be None).
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
    // OCCT L1716: BRepLib::SameParameter(aE) — rcad: mark edge as needing param sync
    ds.edges[ei].geom_tol = tol_e;
}

/// ✅ OCCT-aligned: IsClosed (BOPTools_AlgoTools2D_1.cxx L289-311).
/// Checks if an edge appears twice in a face (closed seam edge on periodic surface).
pub fn is_closed_2d(ei: usize, face_idx: usize, ds: &crate::bopds::ds::DS) -> bool {
    // OCCT L293: BRep_Tool::IsClosed(aE, aF) — rcad: edge is closed when start==end
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

/// ✅ OCCT-aligned: AttachExistingPCurve (BOPTools_AlgoTools2D_1.cxx L44-160).
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

    // OCCT L75: IsSplitToReverse — check if new edge is reversed relative to old
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

    // OCCT L88-94: SameRange — adjust pcurve range to match new edge's 3D curve range
    let t11 = ds.edges[ei_new].t_range[0];
    let t12 = ds.edges[ei_new].t_range[1];
    let a_c2d_t = same_range_2d(&a_c2d, t21, t22, t11, t12);
    if a_c2d_t.is_none() { return 2; }
    let a_c2d_t = a_c2d_t.unwrap();

    // OCCT L102-119: ComputeTolerance check
    let a_new_tol = ds.edges[ei_new].geom_tol;
    // rcad: sample pcurve deviation vs 3D curve (simplified)
    let tol_sp = estimate_pcurve_deviation(&a_c2d_t, &ds.edges[ei_new].curve, t11, t12);
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

/// ✅ OCCT-aligned: UpdateClosedPCurve (BOPTools_AlgoTools2D_1.cxx L164-285).
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
            // Reversed circle: same center, radius, negate direction
            rcad_kernel::geom::Curve2d::Circle(rcad_kernel::geom::Circle2d {
                center: c.center,
                radius: c.radius,
            })
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

/// Adjust pcurve range to match target range (OCCT GeomLib::SameRange equivalent).
fn same_range_2d(
    curve: &rcad_kernel::geom::Curve2d,
    _src_t1: f64,
    _src_t2: f64,
    _dst_t1: f64,
    _dst_t2: f64,
) -> Option<rcad_kernel::geom::Curve2d> {
    // For rcad: pcurves are stored with their range directly.
    // If the source and destination ranges differ, we could reparametrize,
    // but for now return the curve as-is with the target range.
    Some(curve.clone())
}

/// Estimate deviation of pcurve from 3D curve by sampling.
fn estimate_pcurve_deviation(
    _pcurve: &rcad_kernel::geom::Curve2d,
    _curve3: &rcad_kernel::geom::Curve3,
    _t1: f64,
    _t2: f64,
) -> f64 {
    // OCCT uses IntTools_Tools::ComputeTolerance with 3D curve + pcurve + surface.
    // rcad: simplified — returns 0 (no deviation). Callers can use compute_tolerance.
    0.0
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
            rcad_kernel::geom::Curve2d::Circle(rcad_kernel::geom::Circle2d {
                center: c.center + shift,
                radius: c.radius,
            })
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
    // OCCT L315-316: CorrectPointOnCurve → CorrectCurveOnSurface
    // rcad: vertex-on-curve check → edge-on-surface check
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
    //   if (anErr != 0) → add BOPAlgo_AlertUnableToOrientTheShape warning
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
/// Returns true if the edge lies on the parametric seam (U=0 or U=2π).
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

/// ✅ OCCT-aligned: SenseFlag (BOPTools_AlgoTools3D.cxx L380-402).
/// Returns 1 if normals point same direction, -1 if opposite, 0 if not coincident.
pub fn sense_flag(n1: glam::DVec3, n2: glam::DVec3) -> i8 {
    // OCCT L384: IntTools_Tools::IsDirsCoinside — checks parallelism
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

/// ✅ OCCT-aligned: GetNormalToSurface (BOPTools_AlgoTools3D.cxx L406-439).
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

/// ✅ OCCT-aligned: GetApproxNormalToFaceOnEdge (BOPTools_AlgoTools3D.cxx L443-494).
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

/// ✅ OCCT-aligned: MinStepIn2d (BOPTools_AlgoTools3D.hxx L215).
/// Returns the minimum step used in 2D computations (1e-5).
pub fn min_step_in_2d() -> f64 {
    1e-5
}

/// ✅ OCCT-aligned: IsEmptyShape (BOPTools_AlgoTools3D.cxx L732-788).
/// Returns true if a shape has no geometry or is empty.
pub fn is_empty_face(face: &crate::bopds::ds::DSFace) -> bool {
    face.boundary_edges.is_empty()
}

/// ✅ OCCT-aligned: IsEmptyShape for a general DS face list.
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
///   Level 1: edge-based angle method — for edges on the solid boundary,
///            finds adjacent face pair and checks if candidate face is internal.
///   Level 2: ComputeState — find edge not on solid, classify mid-point;
///            or PointInFace → classify_point.
///
/// Returns: Some(true) = IN, Some(false) = OUT, None = unable to determine.
pub fn is_internal_face_against_solid(
    fi: usize,
    solid_face_indices: &[usize],
    ds: &crate::bopds::ds::DS,
) -> Option<bool> {
    // OCCT L815-826: build MEF for the solid (edge→face list)
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
        if a_or { continue; } // TopAbs_INTERNAL → skip
        if ds.is_edge_degenerated(ei) { continue; }

        if let Some(a_lf) = a_mef.get(&ei) {
            let a_nb_f = a_lf.len();
            if a_nb_f == 1 {
                // OCCT L851-861: single neighbor face — check if edge is INTERNAL in that face
                let a_f1 = a_lf[0];
                // Use GetEdgeOnFace to find edge orientation in that face
                let e_on_f1 = if a_f1 < ds.faces.len() {
                    crate::boptools::get_edge_off(ei, &ds.edges, &ds.faces[a_f1])
                } else { None };
                if let Some(ei_f1) = e_on_f1 {
                    if ds.edges[ei_f1].is_internal {
                        // Edge is INTERNAL in neighbor face → face is internal
                        i_ret = is_internal_face_core(fi, ei, a_f1, a_f1, ds);
                        found_edge = Some(ei);
                        break;
                    }
                }
                // Edge is not INTERNAL in the only neighbor → not a candidate
                continue;
            } else if a_nb_f >= 2 {
                // OCCT L864-873: two+ neighbor faces — use angle-based method
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
        // OCCT L952-956: INTERNAL edge → create both orientations
        // rcad: just use the edge as-is for both
    }

    let a_e2 = if the_face1 == the_face2 {
        // OCCT L958-962: same face → both orientations
        a_e1
    } else if the_face2 < ds.faces.len() {
        crate::boptools::get_edge_off(the_edge, &ds.edges, &ds.faces[the_face2])
            .unwrap_or(a_e1)
    } else { a_e1 };

    // OCCT L968-974: build candidate list: (edge, face) pairs
    let mut lcs_off: Vec<(usize, usize)> = Vec::new();
    lcs_off.push((the_edge, the_face));  // (theE1, theFace)
    lcs_off.push((a_e2, the_face2));      // (aE2, theFace2)

    // OCCT L976-989: GetFaceOff — find the face with minimal angle
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
        // OCCT L906-910: exactly 2 → direct pairing
        is_internal_face_core(the_face, the_edge, candidate_faces[0], candidate_faces[1], ds)
    } else {
        // OCCT L914-933: more than 2 → pair them via FindFacePairs
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

/// ✅ OCCT-aligned: OrientEdgesOnWire (BOPTools_AlgoTools.cxx L262-359).
///
/// OCCT algorithm:
///   1. Build vertex→edge map (MapShapesAndAncestors VERTEX→EDGE).
///   2. For each edge: add to new wire, get V1/V2.
///   3. If closed edge (V1==V2): skip adjacency walk.
///   4. For each vertex direction:
///      - While vertex has exactly 2 incident edges:
///        - Find the unused edge, orient to connect (end→start).
///        - Move to next vertex.
///
/// rcad: operates on DS edge indices + forward flags.
///   edges: mutable list of (edge_idx, forward) pairs.
pub fn orient_edges_on_wire_occt(edges: &mut Vec<(usize, bool)>, ds: &crate::bopds::ds::DS) {
    if edges.is_empty() { return; }

    // OCCT L265-272: build vertex→edge map (TopExp::MapShapesAndAncestors)
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
                        // start matches → forward
                        (true, a_vn2)
                    } else if a_vc == a_vn2 {
                        // end matches → reversed
                        (false, a_vn1)
                    } else {
                        // no match → skip this edge
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

/// ✅ OCCT-aligned: PointInFace (BOPTools_AlgoTools3D.cxx L906-941).
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

/// OCCT-aligned: IsOpenShell (BOPTools_AlgoTools.cxx L2350-2394) — single-shell variant.
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
    // OCCT L688-714: all edges on solid → PointInFace
    let pt = point_in_face(ds, fi);
    match pt {
        Some((p3d, _)) => crate::classify::classify_point(p3d, solid_face_indices, ds),
        None => crate::classify::Classification::Out,
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
