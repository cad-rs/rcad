//! OCCT BRepOffsetAPI_MakeDraft — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeDraft.cxx (L21-85) +
//!         BRepOffsetAPI_MakeDraft.hxx (L121-182).
//!
//! OCCT inheritance chain (hxx L121): BRepOffsetAPI_MakeDraft ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated, the Done/NotDone flag) are kept as plain
//! fields of the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepFill_Draft myDraft (TKBool/BRepFill) is the
//!    brep_fill/brep_fill_draft.rs translation (imported below following
//!    the OCCT hxx member form).
//! 3. OCCT enum BRepBuilderAPI_TransitionMode (TKBRep/BRepBuilderAPI) is the
//!    local enum below (no rcad BRepBuilderAPI enum module yet);
//!    BRepFill_TransitionStyle (TKBool/BRepFill) is the engine enum,
//!    re-exported under its OCCT-named path (the sibling facade
//!    brep_offset_api_make_pipe_shell.rs imports it from here).

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo_shape::Shape;

use glam::DVec3;

use crate::brep_fill::brep_fill_draft::BRepFillDraft;

// OCCT BRepFill_TransitionStyle (TKBool/BRepFill_TransitionStyle.hxx L23-27)
// — the engine enum, re-exported for the sibling facades.
pub use crate::brep_fill::brep_fill_pipe::BRepFillTransitionStyle;

/// OCCT BRepBuilderAPI_TransitionMode (TKBRep/BRepBuilderAPI_TransitionMode.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepBuilderAPITransitionMode {
    TransitionMode_Transformed,
    RightCorner,
    RoundCorner,
}

/// OCCT BRepOffsetAPI_MakeDraft (hxx L121-182).
pub struct BRepOffsetAPIMakeDraft {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private member (hxx L181).
    my_draft: BRepFillDraft, // OCCT: myDraft
}

impl BRepOffsetAPIMakeDraft {
    /// OCCT BRepOffsetAPI_MakeDraft::BRepOffsetAPI_MakeDraft(Shape, Dir,
    /// Angle) (cxx L21-27).
    pub fn new(the_shape: &Shape, the_dir: &DVec3, the_angle: f64) -> Self {
        let mut r = BRepOffsetAPIMakeDraft {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            // OCCT L24: myDraft(Shape, Dir, Angle).
            my_draft: BRepFillDraft::new(the_shape, *the_dir, the_angle),
        };
        // OCCT L26: NotDone().
        r.my_done = false;
        r
    }

    /// OCCT BRepOffsetAPI_MakeDraft::SetOptions(Style, AngleMin, AngleMax)
    /// (cxx L29-39).
    pub fn set_options(
        &mut self,
        style: BRepBuilderAPITransitionMode,
        angle_min: f64,
        angle_max: f64,
    ) {
        // OCCT L33-37.
        let mut style_bf = BRepFillTransitionStyle::TransitionStyle_Right;
        if style == BRepBuilderAPITransitionMode::RoundCorner {
            style_bf = BRepFillTransitionStyle::TransitionStyle_Round;
        }
        self.my_draft.set_options(style_bf, angle_min, angle_max);
    }

    /// OCCT BRepOffsetAPI_MakeDraft::SetDraft(IsInternal) (cxx L41-44).
    pub fn set_draft(&mut self, is_internal: bool) {
        self.my_draft.set_draft(is_internal);
    }

    /// OCCT BRepOffsetAPI_MakeDraft::Perform(LengthMax) (cxx L46-54).
    pub fn perform(&mut self, length_max: f64) {
        self.my_draft.perform(length_max);
        if self.my_draft.is_done() {
            // OCCT L51: Done().
            self.my_done = true;
            self.my_shape = self.my_draft.shape();
        }
    }

    /// OCCT BRepOffsetAPI_MakeDraft::Perform(Surface, KeepInsideSurface)
    /// (cxx L56-65).
    pub fn perform_with_surface(&mut self, surface: &Surface3, keep_inside_surface: bool) {
        self.my_draft
            .perform_with_surface(surface, keep_inside_surface);
        if self.my_draft.is_done() {
            // OCCT L61: Done().
            self.my_done = true;
            self.my_shape = self.my_draft.shape();
        }
    }

    /// OCCT BRepOffsetAPI_MakeDraft::Perform(StopShape, KeepOutSide)
    /// (cxx L67-75).
    pub fn perform_with_shape(&mut self, stop_shape: &Shape, keep_out_side: bool) {
        self.my_draft.perform_with_shape(stop_shape, keep_out_side);
        if self.my_draft.is_done() {
            // OCCT L72: Done().
            self.my_done = true;
            self.my_shape = self.my_draft.shape();
        }
    }

    /// OCCT BRepOffsetAPI_MakeDraft::Shell() (cxx L77-80).
    pub fn shell(&self) -> Shape {
        self.my_draft.shell()
    }

    /// OCCT BRepOffsetAPI_MakeDraft::Generated(S) (cxx L82-85).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        self.my_draft.generated(s)
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
            "StdFail_NotDone: BRepOffsetAPI_MakeDraft::Shape()"
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
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            Vec::new(),
        ))
    }
}
