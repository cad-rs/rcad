//! Free-form BRep construction (OCCT `BRepBuilderAPI` equivalent).
//!
//! These functions provide a low-level API to incrementally build a BRep by
//! appending curves, edges, wires, faces, and solids. The caller is responsible
//! for topological consistency.

use std::sync::Arc;

use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Surface3, SurfaceEval};
use rcad_kernel::topods;
use rcad_kernel::topods::{Orientation, Shape};
use rcad_kernel::topology::{Face, Shell, Wire, WireEdge};

use crate::builder::BuildError;

/// Convert a tshape index to a `Shape` (valid for vertices/edges that exist as TShapes).
fn idx_to_shaperef(brep: &BRep, idx: usize) -> Option<Shape> {
    brep.tshapes.get(idx).map(|ts| Shape {
        data: ts.clone(),
        index: idx,
        orientation: Orientation::Forward,
        location: 0,
    })
}

/// Build a TShape::Wire from a topology `Wire`, adding it to the BRep.
fn birep_wire(brep: &mut BRep, wire: &Wire) -> Shape {
    let edge_refs: Vec<Shape> = wire
        .edges
        .iter()
        .map(|we| {
            brep.tshapes
                .get(we.idx)
                .map(|ts| Shape {
                    data: ts.clone(),
                    index: we.idx,
                    orientation: if we.forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    },
                    location: 0,
                })
                .unwrap_or(Shape::null())
        })
        .collect();
    brep.add_twire(edge_refs)
}

/// Adds a new vertex to the BRep and returns its tshape index.
pub fn make_vertex(brep: &mut BRep, point: glam::DVec3) -> usize {
    brep.add_tvertex(point).index
}

/// Adds a new edge (with associated curve and parameter range) to the BRep.
///
/// Returns the tshape index of the new edge.
pub fn make_edge(
    brep: &mut BRep,
    curve: Curve3,
    t1: f64,
    t2: f64,
    v0: usize,
    v1: usize,
) -> Result<usize, BuildError> {
    let v0_sr = brep
        .tshapes
        .get(v0)
        .map(|ts| Shape {
            data: ts.clone(),
            index: v0,
            orientation: Orientation::Forward,
            location: 0,
        })
        .ok_or(BuildError::InvalidIndex(v0))?;
    let v1_sr = brep
        .tshapes
        .get(v1)
        .map(|ts| Shape {
            data: ts.clone(),
            index: v1,
            orientation: Orientation::Forward,
            location: 0,
        })
        .ok_or(BuildError::InvalidIndex(v1))?;
    let sr = brep.add_tedge(Some(curve), v0_sr, v1_sr, [t1, t2]);
    Ok(sr.index)
}

/// Constructs a `Wire` from a list of `WireEdge`s without modifying the BRep.
pub fn make_wire(edges: Vec<WireEdge>) -> Wire {
    Wire { edges }
}

/// Adds a new face to the BRep and returns its tshape index.
///
/// The face normal is derived from `surface.normal_at(0.0, 0.0)`.
pub fn make_face(
    brep: &mut BRep,
    surface: Surface3,
    outer: Wire,
    inner_wires: Vec<Wire>,
) -> Result<usize, BuildError> {
    if outer.edges.is_empty() {
        return Err(BuildError::DegenerateGeometry("outer wire has no edges"));
    }
    let outer_sr = birep_wire(brep, &outer);
    let inner_srs: Vec<Shape> = inner_wires.iter().map(|w| birep_wire(brep, w)).collect();
    let sr = brep.add_tface(Some(surface), outer_sr, inner_srs, None, None, vec![], true);
    Ok(sr.index)
}

/// Appends a new solid (composed of the given shells) to the BRep and returns
/// its tshape index.
///
/// Each face in each shell should already have been created via [`make_face`].
/// The shells are searched by face wire structure to find the matching TShapes.
pub fn make_solid(brep: &mut BRep, shells: Vec<Shell>) -> usize {
    let shell_srs: Vec<Shape> = shells
        .iter()
        .map(|shell| {
            let face_srs: Vec<Shape> = shell
                .faces
                .iter()
                .filter_map(|_face| {
                    // Search for the face TShape by matching wire structure.
                    // This is a best-effort lookup; in practice callers should
                    // use higher-level operations that return Shape directly.
                    brep.tshapes.iter().enumerate().find_map(|(_i, ts)| {
                        if let topods::TShape::Face(fd) = &**ts {
                            // Match by wire edge count as a heuristic
                            if fd.outer_wire.index < brep.tshapes.len() {
                                if let topods::TShape::Wire(wd) =
                                    &*brep.tshapes[fd.outer_wire.index]
                                {
                                    if wd.edges.len() == _face.outer_wire.edges.len() {
                                        return idx_to_shaperef(brep, fd.outer_wire.index);
                                    }
                                }
                            }
                        }
                        None
                    })
                })
                .collect();
            brep.add_tshell(face_srs)
        })
        .collect();
    let sr = brep.add_tsolid(shell_srs);
    sr.index
}
