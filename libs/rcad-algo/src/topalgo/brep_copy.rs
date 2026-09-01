//! OCCT BRepBuilderAPI_Copy (TKTopAlgo) — deep copy of a shape.
//!
//! BRepBuilderAPI_Copy(S, copyGeom) builds a new TopoDS_Shape whose TShapes
//! are freshly allocated (BRepBuilderAPI_Copy.cxx Perform, via
//! BRepTools_Modifier with an empty modification): with copyGeom=true every
//! geometry TShape is duplicated, so the copy is NOT IsEqual to the source
//! (different TShape pointers); with copyGeom=false the TShapes are shared
//! with the source.

use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{self, BRep, CurveRepresentation, TShape};

/// OCCT BRepBuilderAPI_Copy(shape, copyGeom).Shape() — a new BRep whose
/// TShapes are deep-copied when `copy_geom` is true, shared otherwise.
/// The graph copy (crate::bop::ds::clone_shape_graph / clone_tshape) mirrors
/// the BRepTools_Modifier copy of every TShape with fresh handles; the
/// face/vertex ptr keys stored on the edges (pcurves, vertex_params,
/// representations) are remapped to the cloned handles, exactly like
/// crate::bop::ds::clone_arguments_private.
pub fn copy_brep(brep: &BRep, copy_geom: bool) -> BRep {
    let mut out = BRep::new();
    out.locations = brep.locations.clone();
    if copy_geom {
        let mut cache: HashMap<u64, Arc<TShape>> = HashMap::new();
        for (i, ts) in brep.tshapes.iter().enumerate() {
            // clone_shape_graph returns a Shape whose data Arc comes from the
            // cache — the SAME Arc the edge ptr-key remap below resolves to.
            let shape = Shape::from_parts(ts.clone(), i, 0, topods::Orientation::Forward);
            let cloned = crate::bop::ds::clone_shape_graph(&shape, &mut cache);
            out.tshapes.push(cloned.data);
        }
        // The deep copy rebuilds every Arc, so the identity keys (OCCT
        // TopoDS_Shape::IsSame semantics) stored on the edges are stale —
        // remap each from the source TShape ptr to the cloned Arc ptr.
        let remap: HashMap<u64, u64> = cache
            .iter()
            .map(|(k, v)| (*k, std::sync::Arc::as_ptr(v) as u64))
            .collect();
        for ts in &mut out.tshapes {
            // The cloned TShapes are shared through multiple references (the
            // cache dedups Arc clones), so Arc::make_mut would clone the Arc
            // and split the topology.  Mutate in place like
            // crate::bop::ds::clone_arguments_private — OCCT's TShape is a
            // shared handle, so the in-place mutation (single-threaded) is
            // safe and matches that model.
            let ptr = Arc::as_ptr(ts) as *mut TShape;
            let ts = unsafe { &mut *ptr };
            let TShape::Edge(ed) = ts else { continue };
            ed.pcurves = ed
                .pcurves
                .iter()
                .map(|(&(p, l), v)| ((remap.get(&p).copied().unwrap_or(p), l), v.clone()))
                .collect();
            if !ed.vertex_params.is_empty() {
                ed.vertex_params = ed
                    .vertex_params
                    .iter()
                    .map(|(&k, &v)| (remap.get(&k).copied().unwrap_or(k), v))
                    .collect();
            }
            ed.representations = ed
                .representations
                .iter()
                .map(|r| match r {
                    CurveRepresentation::CurveOnSurface {
                        face: (p, l),
                        pcurve,
                        range,
                    } => CurveRepresentation::CurveOnSurface {
                        face: (remap.get(p).copied().unwrap_or(*p), *l),
                        pcurve: pcurve.clone(),
                        range: *range,
                    },
                    CurveRepresentation::CurveOnClosedSurface {
                        face: (p, l),
                        pcurve1,
                        pcurve2,
                        range,
                    } => CurveRepresentation::CurveOnClosedSurface {
                        face: (remap.get(p).copied().unwrap_or(*p), *l),
                        pcurve1: pcurve1.clone(),
                        pcurve2: pcurve2.clone(),
                        range: *range,
                    },
                    other => other.clone(),
                })
                .collect();
        }
    } else {
        out.tshapes = brep.tshapes.clone();
    }
    out
}

/// Shape-level deep copy (BRepBuilderAPI_Copy of a single shape reference):
/// clones the referenced TShape graph with a shared cache, preserving the
/// index/location/orientation.
pub fn copy_shape(s: &Shape, cache: &mut HashMap<u64, Arc<TShape>>) -> Shape {
    crate::bop::ds::clone_shape_graph(s, cache)
}
