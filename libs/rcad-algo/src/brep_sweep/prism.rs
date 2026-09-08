//! OCCT BRepSweep_Prism (TKPrim/BRepSweep) — natural constructors to build
//! BRepSweep translated swept Primitives.
//!
//! Sources:
//! - BRepSweep_Prism.hxx L34-104
//! - BRepSweep_Prism.cxx L30-145

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topods::ShapeType;

use super::num_linear_regular_sweep::NumLinearRegularSweepSlots;
use super::sweep_num_shape::SweepNumShape;
use super::translation::BRepSweepTranslation;

/// OCCT BRepSweep_Prism (BRepSweep_Prism.hxx L34-104).
pub struct BRepSweepPrism {
    /// OCCT: myTranslation (BRepSweep_Translation).
    pub my_translation: BRepSweepTranslation,
}

impl BRepSweepPrism {
    /// OCCT BRepSweep_Prism::BRepSweep_Prism(S, V, C, Canonize) (cxx L30-38)
    /// — builds the prism of base S and vector V.  If C is true, S is
    /// copied.  If Canonize is true then generated surfaces are attempted to
    /// be canonized in simple types.
    pub fn with_vec(s: &Shape, v: glam::DVec3, c: bool, canonize: bool) -> Self {
        // OCCT L34: myTranslation(S, NumShape(), Location(V), V, C, Canonize)
        // — NumShape() = Sweep_NumShape(2, TopAbs_EDGE) (cxx L102-105);
        // Location(V) = the translation location (cxx L125-131).
        let my_translation = BRepSweepTranslation::new(
            s,
            &Self::num_shape(),
            Self::location_vec(v),
            v,
            c,
            canonize,
        );
        // OCCT L36-37: Standard_ConstructionError_Raise_if(V.Magnitude() <=
        // Precision::Confusion(), "BRepSweep_Prism::Constructor").
        if v.length() <= CONFUSION {
            panic!("Standard_ConstructionError: BRepSweep_Prism::Constructor");
        }
        BRepSweepPrism { my_translation }
    }

    /// OCCT BRepSweep_Prism::BRepSweep_Prism(S, D, Inf, C, Canonize)
    /// (cxx L42-49) — builds a semi-infinite or an infinite prism of base S.
    pub fn with_dir(s: &Shape, d: glam::DVec3, inf: bool, c: bool, canonize: bool) -> Self {
        // OCCT L47: myTranslation(S, NumShape(Inf), Location(D), D, C,
        // Canonize).
        let my_translation = BRepSweepTranslation::new(
            s,
            &Self::num_shape_inf(inf),
            Self::location_vec(d),
            d,
            c,
            canonize,
        );
        BRepSweepPrism { my_translation }
    }

    /// OCCT BRepSweep_Prism::Shape() (cxx L53-56) — the TopoDS Shape attached
    /// to the prism.
    pub fn shape(&mut self) -> Shape {
        self.my_translation.shape_full()
    }

    /// OCCT BRepSweep_Prism::Shape(aGenS) (cxx L60-63).
    pub fn shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_translation.shape_of_gen(a_gen_s)
    }

    /// OCCT BRepSweep_Prism::FirstShape() (cxx L67-70) — the bottom of the
    /// prism.
    pub fn first_shape(&mut self) -> Shape {
        self.my_translation.first_shape_full()
    }

    /// OCCT BRepSweep_Prism::FirstShape(aGenS) (cxx L74-77).
    pub fn first_shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_translation.first_shape(a_gen_s)
    }

    /// OCCT BRepSweep_Prism::LastShape() (cxx L81-84) — the top of the prism.
    pub fn last_shape(&mut self) -> Shape {
        self.my_translation.last_shape_full()
    }

    /// OCCT BRepSweep_Prism::LastShape(aGenS) (cxx L88-91).
    pub fn last_shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_translation.last_shape(a_gen_s)
    }

    /// OCCT BRepSweep_Prism::Vec() (cxx L95-98) — the vector of the Prism; if
    /// it is an infinite prism the Vec is unitar.
    pub fn vec(&self) -> glam::DVec3 {
        self.my_translation.vec()
    }

    /// OCCT BRepSweep_Prism::IsUsed(aGenS) (cxx L135-138).
    pub fn is_used(&self, a_gen_s: &Shape) -> bool {
        self.my_translation.is_used(a_gen_s)
    }

    /// OCCT BRepSweep_Prism::GenIsUsed(theS) (cxx L142-145).
    pub fn gen_is_used(&self, the_s: &Shape) -> bool {
        self.my_translation.gen_is_used(the_s)
    }

    /// OCCT BRepSweep_Prism::NumShape() (cxx L102-105) — the NumShape of a
    /// limited prism (private).
    fn num_shape() -> SweepNumShape {
        SweepNumShape::new_indexed(2, ShapeType::Edge)
    }

    /// OCCT BRepSweep_Prism::NumShape(Inf) (cxx L109-121) — the NumShape of
    /// an infinite prism (private).
    fn num_shape_inf(inf: bool) -> SweepNumShape {
        let mut n = SweepNumShape::new();
        if inf {
            n.init(0, ShapeType::Edge, false, true, true);
        } else {
            n.init(1, ShapeType::Edge, false, false, true);
        }
        n
    }

    /// OCCT BRepSweep_Prism::Location(V) (cxx L125-131) — the translation
    /// location (private).  Architecture difference: the TopLoc_Location
    /// travels as a DAffine3.
    fn location_vec(v: glam::DVec3) -> glam::DAffine3 {
        glam::DAffine3::from_translation(v)
    }
}


#[cfg(test)]
#[cfg(test)]
mod prism_tests {
    //! Phase-2 anchor assets (the translation-period acceptance is the form
    //! audit; these smoke anchors guard the end-to-end engine once the
    //! phase-2 debugging opens).

    use super::*;
    use crate::brep_sweep::BRepSweepBuilder;
    use crate::brep_algo::tool::explorer;
    use rcad_kernel::geom::{Plane, Surface3};
    use rcad_kernel::topods::ShapeType;
    use std::collections::HashSet;

    /// A unit square face in the z = 0 plane (the classic prism test
    /// generatrix), built with the sweep's own BRep_Builder re-host stack
    /// (MakeVertex / MakeEdge(C, Tol) / Add(E, V) / MakeWire / Add(W, E) /
    /// MakeFace(S, Tol) / Add(F, W)).
    fn unit_square_face() -> Shape {
        let b = BRepSweepBuilder::new().my_builder;
        let corners = [
            glam::DVec3::new(0.0, 0.0, 0.0),
            glam::DVec3::new(1.0, 0.0, 0.0),
            glam::DVec3::new(1.0, 1.0, 0.0),
            glam::DVec3::new(0.0, 1.0, 0.0),
        ];
        let verts: Vec<Shape> = corners
            .iter()
            .map(|p| b.make_vertex(*p, CONFUSION))
            .collect();
        let mut edges = Vec::new();
        for i in 0..4 {
            let p0 = corners[i];
            let p1 = corners[(i + 1) % 4];
            let e = b.make_edge_curve(
                &rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize(),
                }),
                CONFUSION,
            );
            // OCCT BRep_Builder::Range(E, First, Last) on the 3D curve, plus
            // the vertex parameters BRepLib_MakeEdge stores at construction.
            let len = (p1 - p0).length();
            if let rcad_kernel::topods::TShape::Edge(ed) = unsafe {
                &mut *(std::sync::Arc::as_ptr(&e.data) as *mut rcad_kernel::topods::TShape)
            } {
                ed.range = [0.0, len];
                ed.vertex_params.insert(verts[i].ptr_id(), 0.0);
                ed.vertex_params.insert(verts[(i + 1) % 4].ptr_id(), len);
            }
            // Add(E, V): FORWARD -> first, REVERSED -> last.
            b.add(&e, &verts[i]);
            let mut vlast = verts[(i + 1) % 4].clone();
            vlast.orientation = rcad_kernel::topods::Orientation::Reversed;
            b.add(&e, &vlast);
            edges.push(e);
        }
        let wire = b.make_wire();
        for e in &edges {
            b.add(&wire, e);
        }
        let face = b.make_face(
            &Surface3::Plane(Plane::new(glam::DVec3::ZERO, glam::DVec3::Z)),
            CONFUSION,
        );
        b.add(&face, &wire);
        face
    }

    #[test]
    fn prism_of_unit_square_face_builds_a_solid() {
        let base = unit_square_face();
        let mut prism =
            BRepSweepPrism::with_vec(&base, glam::DVec3::new(0.0, 0.0, 2.0), false, true);
        let solid = prism.shape();
        // NOTE: Shape::is_null() keys on index == usize::MAX, which pool-free
        // builder shapes also carry — assert on the type instead (the
        // consumer-side switch notes carry the same caveat).
        assert_eq!(solid.shape_type(), ShapeType::Solid);
        // The classic prism topology: 8 unique vertices / 12 unique edges /
        // 6 faces (2 caps + 4 lateral walls).  The TopExp_Explorer walk
        // counts per-occurrence (OCCT does not dedup), so the unique-TShape
        // counts are the assertion.
        let faces = explorer(&solid, ShapeType::Face, ShapeType::Shape);
        let edges = explorer(&solid, ShapeType::Edge, ShapeType::Shape);
        let verts = explorer(&solid, ShapeType::Vertex, ShapeType::Shape);
        assert_eq!(faces.len(), 6, "6 faces: 2 caps + 4 lateral");
        for (fi, f) in faces.iter().enumerate() {
            if let rcad_kernel::topods::TShape::Face(fd) = f.data.as_ref() {
                let ws: Vec<u64> = explorer(&fd.outer_wire, ShapeType::Edge, ShapeType::Shape).iter().map(|e| e.ptr_id()).collect();
                eprintln!("[FACE {}] surf={:?} edges={:?}", fi, fd.surface.as_ref().map(|s| std::mem::discriminant(s)), ws);
            }
        }
        let uniq = |v: &[Shape]| -> usize {
            v.iter().map(|s| s.ptr_id()).collect::<HashSet<u64>>().len()
        };
        
        for (i, e) in edges.iter().enumerate() {
            if let rcad_kernel::topods::TShape::Edge(ed) = e.data.as_ref() {
                if let Some(c) = &ed.curve {
                    eprintln!("[EDGE {}] ptr={:x} origin={:?} range={:?}", i, e.ptr_id(), c, ed.range);
                }
            }
        }
        assert_eq!(uniq(&edges), 12, "12 unique prism edges");
        assert_eq!(uniq(&verts), 8, "8 unique prism vertices");
    }
}
