//! OCCT ChFi2d package-level shared functions — 1:1 translation.
//!
//! Sources:
//!   - ChFi2d.hxx (package docs, class shell)
//!   - ChFi2d.cxx (CommonVertex, FindConnectedEdges)
//!   - ChFi2d_ConstructionError.hxx (error enum)
//!
//! OCCT ChFi2d class has no data members (only static methods) — the Rust
//! translation keeps them as free functions prefixed `chfi2d_`.

use rcad_kernel::topo::topods::Shape;

// =========================================================================
// OCCT ChFi2d_ConstructionError.hxx L21-39 — error enum.
// The OCCT enum values keep their original order; rcad reuses it as the
// error/status carrier for the whole ChFi2d package.
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFi2dConstructionError {
    /// OCCT ChFi2d_NotPlanar — the face is not planar
    NotPlanar,
    /// OCCT ChFi2d_NoFace — the face is null
    NoFace,
    /// OCCT ChFi2d_InitialisationError — the two faces used for the
    /// initialisation are uncompatible
    InitialisationError,
    /// OCCT ChFi2d_ParametersError — the parameters as distances or angle
    /// for chamfer are less or equal to zero
    ParametersError,
    /// OCCT ChFi2d_Ready — the initialization has been successful
    Ready,
    IsDone,
    /// OCCT ChFi2d_ComputationError — the algorithm could not find a solution
    ComputationError,
    /// OCCT ChFi2d_ConnexionError — the vertex given to locate the fillet or
    /// the chamfer is not connected to 2 edges
    ConnexionError,
    /// OCCT ChFi2d_TangencyError — the two edges connected to the vertex are
    /// tangent
    TangencyError,
    /// OCCT ChFi2d_FirstEdgeDegenerated — the first edge is degenerated
    FirstEdgeDegenerated,
    /// OCCT ChFi2d_LastEdgeDegenerated — the last edge is degenerated
    LastEdgeDegenerated,
    /// OCCT ChFi2d_BothEdgesDegenerated — the two edges are degenerated
    BothEdgesDegenerated,
    /// OCCT ChFi2d_NotAuthorized — One or the two edges connected to the
    /// vertex is a fillet or a chamfer; One or the two edges connected to the
    /// vertex is not a line or a circle
    NotAuthorized,
}

/// OCCT TopTools_ShapeMapHasher identity: TShape + Location
/// (TopTools_ShapeMapHasher.hxx L37: `return S1.IsSame(S2);`).
#[inline]
fn shape_map_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopExp::MapShapesAndAncestors(S, T, A, M) restricted to the
/// (TopAbs_VERTEX, TopAbs_EDGE) instantiation used by ChFi2d.
///
/// Architecture difference: rcad stores the shape tree on the TShape handles
/// themselves, so the traversal reads the face's wires (outer first, then
/// inner) and each wire's edges (TWireData::edges) instead of walking the
/// OCCT handle graph. The result keeps the OCCT
/// NCollection_IndexedDataMap insertion order as a Vec of (vertex key,
/// ancestor edge list) pairs.
fn topexp_map_shapes_and_ancestors_face_edges(face: &Shape) -> Vec<((u64, u32), Vec<Shape>)> {
    let mut map: Vec<((u64, u32), Vec<Shape>)> = Vec::new();

    // Face sub-shape traversal order (OCCT TopExp_Explorer visits wires then
    // edges per wire); rcad TFaceData keeps outer_wire + inner_wires.
    let fd = match face.as_face() {
        Some(fd) => fd,
        None => return map,
    };
    let mut wires: Vec<Shape> = Vec::new();
    if !fd.outer_wire.is_null() {
        wires.push(fd.outer_wire.clone());
    }
    for w in &fd.inner_wires {
        wires.push(w.clone());
    }

    for w in &wires {
        let wd = match w.as_wire() {
            Some(wd) => wd,
            None => continue,
        };
        for e in &wd.edges {
            let ed = match e.as_edge() {
                Some(ed) => ed,
                None => continue,
            };
            for v in [&ed.first, &ed.last] {
                if v.as_vertex().is_none() {
                    continue;
                }
                let key = shape_map_key(v);
                if let Some(entry) = map.iter_mut().find(|(k, _)| *k == key) {
                    entry.1.push(e.clone());
                } else {
                    map.push((key, vec![e.clone()]));
                }
            }
        }
    }
    map
}

/// OCCT ChFi2d::CommonVertex (ChFi2d.cxx L30-47).
///
/// Returns the common vertex and `true` when found (OCCT out parameter `V`
/// becomes the first tuple element; `false` maps to OCCT `return false`).
pub fn chfi2d_common_vertex(e1: &Shape, e2: &Shape) -> (Shape, bool) {
    // OCCT L32-34: TopExp::Vertices(E1, firstVertex1, lastVertex1) etc.
    // Without CumOri the canonical edge ends are the TShape's first/last.
    let ed1 = e1.as_edge().expect("ChFi2d::CommonVertex: not an edge");
    let ed2 = e2.as_edge().expect("ChFi2d::CommonVertex: not an edge");
    let first_vertex1 = &ed1.first;
    let last_vertex1 = &ed1.last;
    let first_vertex2 = &ed2.first;
    let last_vertex2 = &ed2.last;

    // OCCT L36-40
    if first_vertex1.is_same(first_vertex2) || first_vertex1.is_same(last_vertex2) {
        return (first_vertex1.clone(), true);
    }
    // OCCT L41-45
    if last_vertex1.is_same(first_vertex2) || last_vertex1.is_same(last_vertex2) {
        return (last_vertex1.clone(), true);
    }
    // OCCT L46
    (Shape::null(), false)
}

/// OCCT ChFi2d::FindConnectedEdges (ChFi2d.cxx L51-92).
///
/// Locates the two edges of the face `f` connected at the vertex `v` and
/// writes them into `e1` / `e2` (OCCT out parameters).
pub fn chfi2d_find_connected_edges(
    f: &Shape,
    v: &Shape,
    e1: &mut Shape,
    e2: &mut Shape,
) -> ChFi2dConstructionError {
    // OCCT L56-58:
    // NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
    //   TopTools_ShapeMapHasher> vertexMap;
    // TopExp::MapShapesAndAncestors(F, TopAbs_VERTEX, TopAbs_EDGE, vertexMap);
    let vertex_map = topexp_map_shapes_and_ancestors_face_edges(f);

    // OCCT L60: if (vertexMap.Contains(V))
    if let Some((_, ancestors)) = vertex_map
        .iter()
        .find(|(key, _)| *key == shape_map_key(v))
    {
        // OCCT L62: list iterator over vertexMap.FindFromKey(V).
        // `it` is the cursor standing in for iterator.More()/Next().
        let list: &[Shape] = ancestors;
        let mut it = 0usize;

        // OCCT L63-71
        if it < list.len() {
            *e1 = list[it].clone();
            it += 1;
        } else {
            return ChFi2dConstructionError::ConnexionError;
        }
        // OCCT L72-80
        if it < list.len() {
            *e2 = list[it].clone();
            it += 1;
        } else {
            return ChFi2dConstructionError::ConnexionError;
        }
        // OCCT L82-85
        if it < list.len() {
            return ChFi2dConstructionError::ConnexionError;
        }
    }
    // OCCT L87-90
    else {
        return ChFi2dConstructionError::ConnexionError;
    }
    // OCCT L91
    ChFi2dConstructionError::IsDone
}
