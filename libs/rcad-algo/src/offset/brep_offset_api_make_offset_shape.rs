//! OCCT BRepOffsetAPI_MakeOffsetShape — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeOffsetShape.cxx (L24-163) +
//!         BRepOffsetAPI_MakeOffsetShape.hxx (L44-189).
//!
//! OCCT inheritance chain (hxx L52): BRepOffsetAPI_MakeOffsetShape ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated, the Done/NotDone flag) are kept as plain
//! fields of the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepOffset_MakeOffset myOffsetShape is the parallel
//!    Stage 2 batch translation in offset/brep_offset_make_offset.rs (the
//!    import follows the OCCT hxx member form); the engine member
//!    BRepOffset_MakeSimpleOffset mySimpleOffsetShape is the landed
//!    offset/brep_offset_make_simple_offset.rs translation.
//! 3. OCCT enum OffsetAlgo_Type (hxx L156-160) is the local enum below.

use rcad_kernel::topo_shape::Shape;

use crate::brep_fill::offset_wire::GeomAbsJoinType;

use super::brep_offset_make_offset::BRepOffsetMakeOffset;
use super::brep_offset_make_offset::BRepOffset_Mode;
use super::brep_offset_make_simple_offset::BRepOffsetMakeSimpleOffset;

/// OCCT enum OffsetAlgo_Type (hxx L156-160).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OffsetAlgoType {
    OffsetAlgoNone,
    OffsetAlgoJoin,
    OffsetAlgoSimple,
}

/// OCCT BRepOffsetAPI_MakeOffsetShape (hxx L52-189).
pub struct BRepOffsetAPIMakeOffsetShape {
    // OCCT BRepBuilderAPI base members.
    pub(crate) my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    pub(crate) my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    pub(crate) my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT protected members (hxx L162-168); pub(crate) carries the C++
    // `protected` visibility for the BRepOffsetAPI_MakeThickSolid subclass
    // composition (Rust has no inheritance).
    pub(crate) my_last_used_algo: OffsetAlgoType,       // OCCT: myLastUsedAlgo
    pub(crate) my_offset_shape: BRepOffsetMakeOffset,   // OCCT: myOffsetShape
    pub(crate) my_simple_offset_shape: BRepOffsetMakeSimpleOffset, // OCCT: mySimpleOffsetShape
}

impl BRepOffsetAPIMakeOffsetShape {
    /// OCCT BRepOffsetAPI_MakeOffsetShape::BRepOffsetAPI_MakeOffsetShape()
    /// (cxx L26-29).
    pub fn new() -> Self {
        BRepOffsetAPIMakeOffsetShape {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_last_used_algo: OffsetAlgoType::OffsetAlgoNone,
            my_offset_shape: BRepOffsetMakeOffset::new(),
            my_simple_offset_shape: BRepOffsetMakeSimpleOffset::new(),
        }
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::PerformByJoin(...) (cxx L31-58).
    pub fn perform_by_join(
        &mut self,
        s: &Shape,
        offset: f64,
        tol: f64,
        mode: BRepOffset_Mode,
        intersection: bool,
        self_inter: bool,
        join: GeomAbsJoinType,
        remove_int_edges: bool,
    ) {
        // OCCT L33: NotDone().
        self.my_done = false;
        // OCCT L34: myLastUsedAlgo = OffsetAlgo_JOIN.
        self.my_last_used_algo = OffsetAlgoType::OffsetAlgoJoin;

        // OCCT L36-38: myOffsetShape.Initialize(S, Offset, Tol, Mode,
        // Intersection, SelfInter, Join, false, RemoveIntEdges)
        // + myOffsetShape.MakeOffsetShape(theRange).
        self.my_offset_shape.initialize(
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
        self.my_offset_shape.make_offset_shape();

        // OCCT L40-43.
        if !self.my_offset_shape.is_done() {
            return;
        }

        // OCCT L45-46.
        self.my_shape = self.my_offset_shape.shape().clone();
        self.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::PerformBySimple(theS,
    /// theOffsetValue) (cxx L60-75).
    pub fn perform_by_simple(&mut self, the_s: &Shape, the_offset_value: f64) {
        // OCCT L62: NotDone().
        self.my_done = false;
        // OCCT L63: myLastUsedAlgo = OffsetAlgo_SIMPLE.
        self.my_last_used_algo = OffsetAlgoType::OffsetAlgoSimple;

        // OCCT L65-67.
        self.my_simple_offset_shape
            .initialize(the_s, the_offset_value);
        self.my_simple_offset_shape.perform();

        // OCCT L69-72.
        if !self.my_simple_offset_shape.is_done() {
            return;
        }

        // OCCT L74-75.
        self.my_shape = self.my_simple_offset_shape.get_result_shape();
        self.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::MakeOffset() (cxx L77-80).
    pub fn make_offset(&self) -> &BRepOffsetMakeOffset {
        &self.my_offset_shape
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::Build(...) (cxx L82-83) — does
    /// nothing.
    pub fn build(&mut self) {}

    /// OCCT BRepOffsetAPI_MakeOffsetShape::Generated(S) (cxx L85-104).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L86: myGenerated.Clear().
        self.my_generated.clear();
        if self.my_last_used_algo == OffsetAlgoType::OffsetAlgoJoin {
            // OCCT L89.
            self.my_generated = self.my_offset_shape.generated(s).to_vec();
        } else if self.my_last_used_algo == OffsetAlgoType::OffsetAlgoSimple {
            // OCCT L92-99.
            let a_gen_shape = self.my_simple_offset_shape.generated(s);
            if !a_gen_shape.is_null() && !a_gen_shape.is_same(s) {
                self.my_generated.push(a_gen_shape);
            }
        }

        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::Modified(S) (cxx L106-125).
    pub fn modified(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L107: myGenerated.Clear().
        self.my_generated.clear();
        if self.my_last_used_algo == OffsetAlgoType::OffsetAlgoJoin {
            // OCCT L110.
            self.my_generated = self.my_offset_shape.modified(s).to_vec();
        } else if self.my_last_used_algo == OffsetAlgoType::OffsetAlgoSimple {
            // OCCT L113-120.
            let a_gen_shape = self.my_simple_offset_shape.modified(s);
            if !a_gen_shape.is_null() && !a_gen_shape.is_same(s) {
                self.my_generated.push(a_gen_shape);
            }
        }

        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::IsDeleted(S) (cxx L127-135).
    pub fn is_deleted(&mut self, s: &Shape) -> bool {
        if self.my_last_used_algo == OffsetAlgoType::OffsetAlgoJoin {
            return self.my_offset_shape.is_deleted(s);
        }
        false
    }

    /// OCCT BRepOffsetAPI_MakeOffsetShape::GetJoinType() (cxx L137-140).
    pub fn get_join_type(&self) -> GeomAbsJoinType {
        self.my_offset_shape.get_join_type()
    }

    /// OCCT BRepBuilderAPI_Command::IsDone() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_Command.hxx: Standard_Boolean IsDone() const);
    /// the rcad field keeps the crate visibility, the accessor restores the
    /// OCCT-public surface (form restoration, no behavior change).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_MakeShape.hxx: Standard_EXPORT const TopoDS_Shape&
    /// Shape() const; raises StdFail_NotDone when not done).
    pub fn shape(&self) -> Shape {
        assert!(
            self.my_done,
            "StdFail_NotDone: BRepOffsetAPI_MakeOffsetShape::Shape()"
        );
        self.my_shape.clone()
    }

    /// Test-world extraction bridge (the FilletResult.brep pattern,
    /// algo_ext::topods_ext::extract_result_brep): flattens the result root
    /// shape into the self-contained BRep the test world consumes
    /// (StepWriter / total_surface_area).  OCCT has no equivalent (the
    /// TopoDS_Shape carries its arena implicitly) — the rcad BRep-pool
    /// architecture difference #4 glue.
    pub fn result_brep(&mut self) -> Option<rcad_kernel::topo::topods::BRep> {
        if !self.my_done || self.my_shape.is_null() {
            return None;
        }
        let locations = self.my_offset_shape.my_brep.locations.clone();
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            locations,
        ))
    }
}

impl Default for BRepOffsetAPIMakeOffsetShape {
    fn default() -> Self {
        Self::new()
    }
}
