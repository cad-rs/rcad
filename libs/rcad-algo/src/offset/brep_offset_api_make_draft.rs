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
//! 2. The engine member BRepFill_Draft myDraft (TKBool/BRepFill) has no rcad
//!    translation yet — the BRepFillDraft carrier below keeps the OCCT
//!    constructor/accessor surface with GAP panics (port plan section 0.6);
//!    every facade body around the engine calls is translated 1:1.
//! 3. OCCT enums BRepBuilderAPI_TransitionMode (TKBRep/BRepBuilderAPI) and
//!    BRepFill_TransitionStyle (TKBool/BRepFill) are the local enums below
//!    (no rcad BRepBuilderAPI/BRepFill enum modules yet).

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo_shape::Shape;

use glam::DVec3;

/// OCCT BRepBuilderAPI_TransitionMode (TKBRep/BRepBuilderAPI_TransitionMode.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepBuilderAPITransitionMode {
    TransitionMode_Transformed,
    RightCorner,
    RoundCorner,
}

/// OCCT BRepFill_TransitionStyle (TKBool/BRepFill_TransitionStyle.hxx L23-27).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum BRepFillTransitionStyle {
    TransitionStyle_Modified,
    TransitionStyle_Right,
    TransitionStyle_Round,
}

// ---------------------------------------------------------------------------
// GAP carrier (architecture difference #2).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Draft (TKBool/BRepFill, BRepFill_Draft.hxx L42-99) — the
/// draft-surface engine of MakeDraft (architecture difference #2; GAP: no
/// rcad translation yet — the GAP panics are the section 0.6 annotation; the
/// constructor and field storage keep the OCCT form).
pub struct BRepFillDraft {
    my_shape: Shape,         // OCCT: myShape
    my_shell: Shape,         // OCCT: myShell
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated
    my_done: bool,           // OCCT: myDone
}

impl BRepFillDraft {
    /// OCCT BRepFill_Draft::BRepFill_Draft(Shape, Dir, Angle)
    /// (BRepFill_Draft.hxx L42) — GAP: the engine computation is not
    /// translated; the storage form is kept.
    pub fn new(_the_shape: &Shape, _the_dir: &DVec3, _the_angle: f64) -> Self {
        BRepFillDraft {
            my_shape: Shape::null(),
            my_shell: Shape::null(),
            my_generated: Vec::new(),
            my_done: false,
        }
    }

    /// OCCT BRepFill_Draft::SetOptions(Style, AngleMin, AngleMax)
    /// (BRepFill_Draft.hxx L44) — GAP.
    pub fn set_options(
        &mut self,
        _the_style: BRepFillTransitionStyle,
        _the_angle_min: f64,
        _the_angle_max: f64,
    ) {
        panic!("GAP: BRepFill_Draft::SetOptions (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Draft::SetDraft(IsInternal) (BRepFill_Draft.hxx L48) —
    /// GAP.
    pub fn set_draft(&mut self, _the_is_internal: bool) {
        panic!("GAP: BRepFill_Draft::SetDraft (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Draft::Perform(LengthMax) (BRepFill_Draft.hxx L50) — GAP.
    pub fn perform(&mut self, _the_length_max: f64) {
        panic!("GAP: BRepFill_Draft::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Draft::Perform(Surface, KeepInsideSurface)
    /// (BRepFill_Draft.hxx L52) — GAP.
    pub fn perform_with_surface(&mut self, _the_surface: &Surface3, _the_keep_inside: bool) {
        panic!("GAP: BRepFill_Draft::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Draft::Perform(StopShape, KeepOutSide)
    /// (BRepFill_Draft.hxx L55) — GAP.
    pub fn perform_with_shape(&mut self, _the_stop_shape: &Shape, _the_keep_outside: bool) {
        panic!("GAP: BRepFill_Draft::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Draft::IsDone() (BRepFill_Draft.hxx L57).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepFill_Draft::Shape() (BRepFill_Draft.hxx L68).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_Draft::Shell() (BRepFill_Draft.hxx L62).
    pub fn shell(&self) -> Shape {
        self.my_shell.clone()
    }

    /// OCCT BRepFill_Draft::Generated(S) (BRepFill_Draft.hxx L66).
    pub fn generated(&self, _the_s: &Shape) -> &Vec<Shape> {
        panic!("GAP: BRepFill_Draft::Generated (TKBool/BRepFill not translated)");
    }
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
            my_draft: BRepFillDraft::new(the_shape, the_dir, the_angle),
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
        self.my_draft.generated(s).clone()
    }
}
