// OCCT BRepBndLib — bounding boxes for topological shapes.
//
// OCCT ref: TKTopAlgo/BRepBndLib/BRepBndLib.hxx / .cxx
//
// Provides static methods to compute bounding boxes of shapes.
// rcad: uses rcad-kernel BndBox internally.

use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{TShape, ShapeType};
use glam::DVec3;

/// OCCT BRepBndLib — static methods for computing bounding boxes of topological shapes.
pub struct BRepBndLib;

impl BRepBndLib {
    /// OCCT: Add(S, B, useTriangulation) — add shape S to bounding box B.
    ///
    /// Traverses the shape hierarchy and accumulates bounding boxes of
    /// all sub-shapes (faces → edges → vertices).
    pub fn add(shape: &Shape, box_: &mut BndBox, _use_triangulation: bool) {
        // OCCT: TopExp_Explorer for FACEs, EDGEs, VERTEXs
        // rcad: recursive traversal of shape sub-shapes
        Self::add_shape(shape, box_);
    }

    /// OCCT: AddClose(S, B) — quick bounding box (assumes polygonal faces).
    pub fn add_close(shape: &Shape, box_: &mut BndBox) {
        Self::add(shape, box_, false);
    }

    /// Recursively add a shape and its sub-shapes to the bounding box.
    fn add_shape(shape: &Shape, box_: &mut BndBox) {
        match &*shape.data {
            TShape::Vertex(vd) => {
                box_.add_point(vd.point);
            }
            TShape::Edge(ed) => {
                // Add edge endpoints
                if let TShape::Vertex(ref v1) = *ed.first.data {
                    box_.add_point(v1.point);
                }
                if let TShape::Vertex(ref v2) = *ed.last.data {
                    box_.add_point(v2.point);
                }
            }
            TShape::Face(fd) => {
                // Add face: collect all vertex points from outer + inner wires
                Self::add_wire(&fd.outer_wire, box_);
                for iw in &fd.inner_wires {
                    Self::add_wire(iw, box_);
                }
            }
            TShape::Shell(sd) => {
                for f in &sd.faces {
                    Self::add_shape(f, box_);
                }
            }
            TShape::Solid(sd) => {
                for sh in &sd.shells {
                    Self::add_shape(sh, box_);
                }
            }
            TShape::Compound(cd) => {
                for s in cd {
                    Self::add_shape(s, box_);
                }
            }
            _ => {}
        }
    }

    /// Add a wire's edge vertices to the bounding box.
    fn add_wire(wire: &Shape, box_: &mut BndBox) {
        if let TShape::Wire(ref wd) = *wire.data {
            for e in &wd.edges {
                Self::add_shape(e, box_);
            }
        }
    }
}
