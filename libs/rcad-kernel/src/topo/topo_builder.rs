// OCCT TopoDS_Builder.cxx L28-103 (MakeShape / Add) + TopoDS_Iterator.cxx
// L24-80 (the myShapes walk) + BRep_Builder.hxx L57
// (`class BRep_Builder : public TopoDS_Builder`) + TopoDS_TShape.hxx L185
// (`NCollection_List<TopoDS_Shape> myShapes;`).
//
// Architecture difference (the rcad BRep pool).  In OCCT a TopoDS_Shape holds
// its own TShape handle, so a shape handed out by any builder stays usable
// everywhere and no re-hosting statement exists anywhere in the library.  rcad
// splits the pair: `topo_shape::Shape` carries the TShape handle PLUS a flat
// `index`, and every BRepTool accessor reads `brep.tshapes[shape.index]`
// (topods.rs L1762-1802).  A shape is therefore only resolvable inside a pool
// that holds its TShape at that flat index.  The function below is that glue:
// it materializes the reachable TShape graph of a shape in a fresh BRep at its
// SOURCE flat index, which is the rcad equivalent of the OCCT statements:
//
//   TopoDS_Builder::MakeShape (L28-33) — bind the shape to its TShape.
//   TopoDS_Builder::Add (L37-103)      — the parent references its components
//                                        (`aTShape->myShapes.Append(aChild)`,
//                                        L91-92); it never copies them, so the
//                                        rcad pool shares the Arc.
//   TopoDS_Iterator::Initialize (L24-57) — `myIterator.Init(S.TShape()->myShapes)`
//                                        (L50) plus the updateCurrentShape
//                                        step (L72-80): the sub-shape walk.
//
// The pool index preservation has no OCCT counterpart (OCCT has no flat
// index), and the component walk is written against the rcad typed component
// fields rather than `TShape::my_shapes` (only the `add_t*` builders fill
// `my_shapes`; the typed fields are the maintained form) — the same
// reachability as `bop/algo/builder.rs` `push_shape_recursive`.

use std::collections::HashSet;
use std::sync::Arc;

use glam::DAffine3;

use crate::geom::Curve2d;
use crate::topo::topo_shape::Shape;
use crate::topo::topods::{BRep, CurveRepresentation, TShape, TVertexData};

/// The pool slot filler for source indices that are not reachable from the
/// adopted root (no OCCT counterpart: OCCT has no pool, hence no flat index
/// that could be looked up without a TShape).  Every `brep.tshapes[i]` below
/// the root index must exist because rcad consumers index the pool directly
/// (`brep_top_shapes`, `BRep::edge_wrapper_locations`, ...).
fn unused_tshape_slot() -> Arc<TShape> {
    Arc::new(TShape::Vertex(TVertexData {
        my_shapes: Vec::new(),
        flags: 0,
        point: glam::DVec3::ZERO,
        tolerance: 0.0,
        points: Vec::new(),
    }))
}

/// Adopt `shape` into a standalone BRep pool.
///
/// OCCT anchors: TopoDS_Builder::MakeShape (TopoDS_Builder.cxx L28-33),
/// TopoDS_Builder::Add (L37-103), TopoDS_Iterator::Initialize
/// (TopoDS_Iterator.cxx L24-57).  See the module header for the rcad
/// BRep-pool architecture difference.
///
/// Every TShape reachable from `shape` is copied to its original flat index
/// (the Arc is SHARED, OCCT TopoDS_Builder.cxx L91-92 references the component
/// TShape), so the resulting pool is directly usable by the pool-bound
/// consumers (`BRepTool` accessors, `brep_top_shapes_with_locations`,
/// `BRepTopolTool`).  `locations` is the locations table of the pool `shape`
/// was built against and is carried over verbatim — a Shape does not reference
/// its source pool, so a shape whose `location` is non-identity can only be
/// adopted correctly when the caller supplies that table (the same convention
/// as `algo_ext::topods_ext::extract_result_brep`).
pub fn brep_from_shape(shape: &Shape, locations: &[DAffine3]) -> BRep {
    let mut out = BRep::new();
    out.locations = locations.to_vec();
    let mut visited: HashSet<u64> = HashSet::new();
    place(&mut out, shape, &mut visited);
    out
}

/// One TopoDS_Iterator::updateCurrentShape step of the adoption walk
/// (TopoDS_Iterator.cxx L72-80): store the component and recurse into its own
/// components.
fn place(brep: &mut BRep, sr: &Shape, visited: &mut HashSet<u64>) {
    if !visited.insert(sr.ptr_id()) {
        return;
    }
    if sr.index == usize::MAX {
        // A pool-free child (the rcad free-standing encoding carries
        // index == usize::MAX): pad-to-index would grow the arena toward
        // usize::MAX and exhaust memory (observed: the 2^36-byte OOM on
        // the offset a3 path).  Skip the materialization — the original
        // handle stays pool-free and the engine dual-path reads it through
        // the pool-free accessors; engine WRITES to such children require
        // the pool-model migration (recorded).
        return;
    }
    if brep.tshapes.len() <= sr.index {
        let dummy = unused_tshape_slot();
        while brep.tshapes.len() <= sr.index {
            brep.tshapes.push(dummy.clone());
        }
    }
    brep.tshapes[sr.index] = sr.data.clone();
    match &*sr.data {
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                place(brep, sh, visited);
            }
            for v in &sd.internal_vertices {
                place(brep, v, visited);
            }
            for e in &sd.internal_edges {
                place(brep, e, visited);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                place(brep, f, visited);
            }
        }
        TShape::Face(fd) => {
            place(brep, &fd.outer_wire, visited);
            for w in &fd.inner_wires {
                place(brep, w, visited);
            }
            for v in &fd.internal_vertices {
                place(brep, v, visited);
            }
        }
        TShape::Wire(wd) => {
            for e in &wd.edges {
                place(brep, e, visited);
            }
        }
        TShape::Edge(ed) => {
            place(brep, &ed.first, visited);
            place(brep, &ed.last, visited);
        }
        TShape::CompSolid(cs) => {
            for s in cs {
                place(brep, s, visited);
            }
        }
        TShape::Compound(cd) => {
            for s in cd {
                place(brep, s, visited);
            }
        }
        // A VERTEX carries children only when a shape was Added into it (the
        // compatibility table TopoDS_Builder.cxx L47-71 admits COMPOUND /
        // SOLID / FACE / EDGE into a VERTEX) and rcad keeps those children in
        // the untyped TVertexData::my_shapes, which no rcad builder
        // populates — an empty walk.
        TShape::Vertex(_) => {}
    }
}

/// Materialize the OCCT invariant "a pcurve IS a BRep_CurveOnSurface
/// representation" over the rcad transition store.
///
/// OCCT anchor: a BRep_TEdge carries ONE curve-representation list
/// (`BRep_TEdge.cxx` L109-122 region — `myCurves`, copied whole by
/// `EmptyCopy`) and `BRep_Tool::CurveOnSurface` reads pcurves from it
/// (BRep_Tool.cxx L339-364).  The rcad `TEdgeData::pcurves` map is the extra
/// transition store still fed by map-only writers; this helper closes the
/// gap left by those writers by pushing, for every map entry whose face key
/// has NO matching `CurveOnSurface` / `CurveOnClosedSurface` representation,
/// the equivalent `CurveOnSurface` representation (pcurve and range taken
/// from the map entry).
///
/// Closed-surface pairs: an entry already backed by a
/// `CurveOnClosedSurface` representation is skipped — that representation
/// is authoritative and carries both pcurves.  A map-only entry cannot
/// reconstruct a second pcurve (the map store holds one pcurve per face
/// key), so it materializes as the single `CurveOnSurface` representation —
/// exactly the answer the map-era reader produced for it.
///
/// The walk mirrors [`place`]: visited set on `ptr_id`, and pool-free
/// children (`index == usize::MAX`) are skipped — their TShapes are not
/// slots of `brep`, so there is nothing to migrate.
pub fn migrate_map_to_representations(brep: &mut BRep, shape: &Shape) {
    let mut visited: HashSet<u64> = HashSet::new();
    let mut edge_indices: Vec<usize> = Vec::new();
    collect_edge_indices(shape, &mut visited, &mut edge_indices);
    for idx in edge_indices {
        let Some(slot) = brep.tshapes.get_mut(idx) else {
            continue;
        };
        if let TShape::Edge(ed) = Arc::make_mut(slot) {
            // Collect the map entries without a matching representation
            // first (the push below borrows the representations field).
            let missing: Vec<((u64, u32), Curve2d, [f64; 2])> = ed
                .pcurves
                .iter()
                .filter(|(face_key, _)| {
                    !ed.representations.iter().any(|a_cr| match a_cr {
                        CurveRepresentation::CurveOnSurface { face, .. }
                        | CurveRepresentation::CurveOnClosedSurface { face, .. } => {
                            face == *face_key
                        }
                        _ => false,
                    })
                })
                .map(|(face_key, (pcurve, a_f, a_l))| {
                    (*face_key, pcurve.clone(), [*a_f, *a_l])
                })
                .collect();
            for (face_key, pcurve, range) in missing {
                ed.representations
                    .push(CurveRepresentation::CurveOnSurface {
                        face: face_key,
                        pcurve,
                        range,
                    });
            }
        }
    }
}

/// The edge collector of [`migrate_map_to_representations`] — one
/// TopoDS_Iterator::updateCurrentShape step (TopoDS_Iterator.cxx L72-80)
/// over the same typed component fields `place` walks.
fn collect_edge_indices(sr: &Shape, visited: &mut HashSet<u64>, out: &mut Vec<usize>) {
    if !visited.insert(sr.ptr_id()) {
        return;
    }
    // The pool-free skip of `place` (see the recorded pool-model note there).
    if sr.index == usize::MAX {
        return;
    }
    match &*sr.data {
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                collect_edge_indices(sh, visited, out);
            }
            for v in &sd.internal_vertices {
                collect_edge_indices(v, visited, out);
            }
            for e in &sd.internal_edges {
                collect_edge_indices(e, visited, out);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                collect_edge_indices(f, visited, out);
            }
        }
        TShape::Face(fd) => {
            collect_edge_indices(&fd.outer_wire, visited, out);
            for w in &fd.inner_wires {
                collect_edge_indices(w, visited, out);
            }
            for v in &fd.internal_vertices {
                collect_edge_indices(v, visited, out);
            }
        }
        TShape::Wire(wd) => {
            for e in &wd.edges {
                collect_edge_indices(e, visited, out);
            }
        }
        TShape::Edge(ed) => {
            out.push(sr.index);
            collect_edge_indices(&ed.first, visited, out);
            collect_edge_indices(&ed.last, visited, out);
        }
        TShape::CompSolid(cs) => {
            for s in cs {
                collect_edge_indices(s, visited, out);
            }
        }
        TShape::Compound(cd) => {
            for s in cd {
                collect_edge_indices(s, visited, out);
            }
        }
        TShape::Vertex(_) => {}
    }
}

/// The OCCT contract behind the adoption: the adopted pool must resolve the
/// root shape's flat index, its location index and every accessor used by the
/// consumers.  Guards the round trip a caller performs after `brep_from_shape`.
pub fn adopted_shape(brep: &BRep, shape: &Shape) -> Shape {
    // OCCT TopoDS_Builder::MakeShape (TopoDS_Builder.cxx L28-33): the Shape
    // bound to this pool keeps its own orientation (S.Orientation is untouched
    // by MakeShape; only Location/Orientation of a NEW shape are reset).  The
    // rcad Shape carries (data, index, location, orientation) already, so the
    // re-hosting is the identity on these four fields.
    Shape {
        data: brep.tshapes[shape.index].clone(),
        index: shape.index,
        location: shape.location,
        orientation: shape.orientation,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topo::topods::{BRepBuilder, ShapeType};
    use glam::DVec3;

    /// A box-less minimal graph: one edge with two vertices, adopted into its
    /// own pool; the flat indices and the Arc identity survive.
    #[test]
    fn brep_from_shape_preserves_flat_index_and_arc() {
        let mut src = BRep::new();
        let mut b = BRepBuilder::new();
        let v0 = b.add_vertex(&mut src, DVec3::ZERO, 1e-7);
        let v1 = b.add_vertex(&mut src, DVec3::X, 1e-7);
        let e = b.add_edge(
            &mut src,
            Some(crate::geom::Curve3::Line(crate::geom::Line3::new(
                DVec3::ZERO,
                DVec3::X,
            ))),
            v0.clone(),
            v1.clone(),
            [0.0, 1.0],
        );

        let dst = brep_from_shape(&e, &[]);
        assert!(dst.tshapes.len() > e.index);
        assert_eq!(dst.tshapes[e.index].shape_type(), ShapeType::Edge);
        // OCCT TopoDS_Builder::Add references the component TShape: the Arc
        // is shared, not copied.
        assert_eq!(dst.tshapes[e.index].as_ref() as *const TShape, e.data.as_ref() as *const TShape);
        let back = adopted_shape(&dst, &e);
        assert!(back.is_same(&e));
        assert_eq!(back.index, e.index);
    }

    /// A vertex reachable from the root is adopted too, at ITS index.
    #[test]
    fn brep_from_shape_adopts_sub_shapes() {
        let mut src = BRep::new();
        let mut b = BRepBuilder::new();
        let v0 = b.add_vertex(&mut src, DVec3::ZERO, 1e-7);
        let v1 = b.add_vertex(&mut src, DVec3::X, 1e-7);
        let e = b.add_edge(
            &mut src,
            Some(crate::geom::Curve3::Line(crate::geom::Line3::new(
                DVec3::ZERO,
                DVec3::X,
            ))),
            v0.clone(),
            v1.clone(),
            [0.0, 1.0],
        );

        let dst = brep_from_shape(&e, &[]);
        assert_eq!(dst.tshapes[v0.index].shape_type(), ShapeType::Vertex);
        assert_eq!(dst.tshapes[v1.index].shape_type(), ShapeType::Vertex);
        // The adopted pool resolves its own accessors (the pool-bound
        // consumer contract: brep.edge(shape)).
        let ed = dst.edge(adopted_shape(&dst, &e));
        assert_eq!(ed.first.index, v0.index);
    }

    /// The locations table is carried over verbatim so location indices of the
    /// adopted shapes resolve in the destination pool.
    #[test]
    fn brep_from_shape_carries_locations() {
        let mut src = BRep::new();
        let mut b = BRepBuilder::new();
        let v = b.add_vertex(&mut src, DVec3::ZERO, 1e-7);
        let loc = src.add_location(glam::DAffine3::from_translation(DVec3::X));
        let mut moved = v.clone();
        moved.location = loc;

        let dst = brep_from_shape(&moved, &src.locations);
        assert_eq!(dst.locations, src.locations);
        assert_eq!(dst.get_location(loc), src.get_location(loc));
    }
}
