//! OCCT BRepOffsetAPI_MakeThickSolid — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeThickSolid.cxx (L30-119) +
//!         BRepOffsetAPI_MakeThickSolid.hxx (L47-158).
//!
//! OCCT inheritance chain (hxx L47): BRepOffsetAPI_MakeThickSolid ->
//! BRepOffsetAPI_MakeOffsetShape.  Rust has no inheritance: the subclass
//! composes the base struct; the C++ `protected` members keep the pub(crate)
//! visibility and are reached through `self.base.` (the Rust-no-inheritance
//! architecture difference, annotated at each use).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepOffset_MakeOffset (base) is the parallel Stage 2
//!    batch translation in offset/brep_offset_make_offset.rs; the engine
//!    BRepOffset_MakeOffset::MakeThickSolid / AddFace /
//!    OffsetFacesFromShapes / ClosingFaces calls keep the OCCT call form.
//! 3. `TopoDS::Face(it.Value())` and `it.ChangeValue().Reverse()` map to the
//!    Shape orientation update in place (NCollection_List::Iterator
//!    ChangeValue is an in-place handle edit).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;

use crate::brep_fill::offset_wire::GeomAbsJoinType;

use super::brep_offset_api_make_offset_shape::{
    BRepOffsetAPIMakeOffsetShape, OffsetAlgoType,
};
use super::brep_offset_make_offset::BRepOffset_Mode;

/// OCCT BRepOffsetAPI_MakeThickSolid (hxx L47-158).
pub struct BRepOffsetAPIMakeThickSolid {
    // OCCT base BRepOffsetAPI_MakeOffsetShape (hxx L47 inheritance).
    pub(crate) base: BRepOffsetAPIMakeOffsetShape,
}

impl BRepOffsetAPIMakeThickSolid {
    /// OCCT BRepOffsetAPI_MakeThickSolid::BRepOffsetAPI_MakeThickSolid()
    /// (cxx L32-36).
    pub fn new() -> Self {
        let mut r = BRepOffsetAPIMakeThickSolid {
            base: BRepOffsetAPIMakeOffsetShape::new(),
        };
        // OCCT L34-35: Build only solids.
        r.base.my_simple_offset_shape.set_build_solid_flag(true);
        r
    }

    /// OCCT BRepOffsetAPI_MakeThickSolid::MakeThickSolidByJoin(...) (cxx
    /// L38-67).
    pub fn make_thick_solid_by_join(
        &mut self,
        s: &Shape,
        closing_faces: &[Shape],
        offset: f64,
        tol: f64,
        mode: BRepOffset_Mode,
        intersection: bool,
        self_inter: bool,
        join: GeomAbsJoinType,
        remove_int_edges: bool,
    ) {
        // OCCT L41: NotDone().
        self.base.my_done = false;
        // OCCT L42: myLastUsedAlgo = OffsetAlgo_JOIN.
        self.base.my_last_used_algo = OffsetAlgoType::OffsetAlgoJoin;

        // OCCT L44-46: myOffsetShape.Initialize(S, Offset, Tol, Mode,
        // Intersection, SelfInter, Join, false, RemoveIntEdges).
        self.base.my_offset_shape.initialize(
            s,
            offset,
            tol,
            mode,
            intersection,
            self_inter,
            join,
            false,
            remove_int_edges,
        );
        // OCCT L47-50: for each closing face, myOffsetShape.AddFace(face).
        for it in closing_faces {
            self.base.my_offset_shape.add_face(it);
        }

        // OCCT L52: myOffsetShape.MakeThickSolid(theRange).
        self.base.my_offset_shape.make_thick_solid();
        // OCCT L53-56.
        if !self.base.my_offset_shape.is_done() {
            return;
        }

        // OCCT L58-59.
        self.base.my_shape = self.base.my_offset_shape.shape().clone();
        self.base.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeThickSolid::MakeThickSolidBySimple(theS,
    /// theOffsetValue) (cxx L69-84).
    pub fn make_thick_solid_by_simple(&mut self, the_s: &Shape, the_offset_value: f64) {
        // OCCT L71: NotDone().
        self.base.my_done = false;
        // OCCT L72: myLastUsedAlgo = OffsetAlgo_SIMPLE.
        self.base.my_last_used_algo = OffsetAlgoType::OffsetAlgoSimple;

        // OCCT L74-76.
        self.base
            .my_simple_offset_shape
            .initialize(the_s, the_offset_value);
        self.base.my_simple_offset_shape.perform();

        // OCCT L78-81.
        if !self.base.my_simple_offset_shape.is_done() {
            return;
        }

        // OCCT L83-84.
        self.base.my_shape = self.base.my_simple_offset_shape.get_result_shape();
        self.base.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeThickSolid::Build(...) (cxx L86-87) — does
    /// nothing.
    pub fn build(&mut self) {}

    /// OCCT BRepOffsetAPI_MakeThickSolid::Modified(S) (cxx L89-119).
    pub fn modified(&mut self, f: &Shape) -> Vec<Shape> {
        // OCCT L90: myGenerated.Clear().
        self.base.my_generated.clear();

        if self.base.my_last_used_algo == OffsetAlgoType::OffsetAlgoJoin
            && self
                .base
                .my_offset_shape
                .offset_faces_from_shapes()
                .has_image(f)
        {
            // OCCT L94-97.
            if self.base.my_offset_shape.closing_faces().contains(f) {
                // OCCT L96: LastImage(F, myGenerated).
                self.base
                    .my_offset_shape
                    .offset_faces_from_shapes()
                    .last_image(f, &mut self.base.my_generated);

                // OCCT L99-109: Reverse generated shapes in case of small
                // solids — the iterator in-place reverse (TopoDS_Shape::
                // Reverse through ChangeValue; the rcad Shape orientation
                // field is the public carrier).
                for it in self.base.my_generated.iter_mut() {
                    it.orientation = match it.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        o => o,
                    };
                }
            }
        } else if self.base.my_last_used_algo == OffsetAlgoType::OffsetAlgoSimple {
            // OCCT L112-118.
            let a_mod_shape = self.base.my_simple_offset_shape.modified(f);
            if !a_mod_shape.is_null() {
                self.base.my_generated.push(a_mod_shape);
            }
        }

        self.base.my_generated.clone()
    }

    /// The composed base (the OCCT BRepOffsetAPI_MakeOffsetShape subobject)
    /// — exposes the inherited accessors (Generated/IsDeleted/... carry the
    /// base behavior).
    pub fn base(&self) -> &BRepOffsetAPIMakeOffsetShape {
        &self.base
    }

    /// The mutable composed base (inherited PerformByJoin/PerformBySimple/
    /// Generated/IsDeleted/GetJoinType/MakeOffset are reachable through it).
    pub fn base_mut(&mut self) -> &mut BRepOffsetAPIMakeOffsetShape {
        &mut self.base
    }
}

impl Default for BRepOffsetAPIMakeThickSolid {
    fn default() -> Self {
        Self::new()
    }
}
