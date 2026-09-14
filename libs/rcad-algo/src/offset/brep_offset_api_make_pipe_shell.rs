//! OCCT BRepOffsetAPI_MakePipeShell — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakePipeShell.cxx (L26-267) +
//!         BRepOffsetAPI_MakePipeShell.hxx (L91-395).
//!
//! OCCT inheritance chain (hxx L91): BRepOffsetAPI_MakePipeShell ->
//! BRepPrimAPI_MakeSweep -> BRepBuilderAPI_MakeShape.  Rust has no
//! inheritance: the base-class members (myShape, myGenerated, the Done
//! flag) are kept as plain fields of the struct (the Stage 2e facade
//! precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member Handle(BRepFill_PipeShell) myPipe (TKBool/BRepFill)
//!    is the landed translation
//!    crate::brep_fill::brep_fill_pipe_shell::BRepFillPipeShell; every
//!    facade body forwards to it following the OCCT facade statements,
//!    including the GeomFill_PipeError -> BRepBuilderAPI_PipeError switch
//!    (cxx L138-156).
//! 3. The rcad BRep-pool architecture difference #4: OCCT shapes carry
//!    their arena implicitly while the rcad engine methods take the pool
//!    as an explicit argument — the facade owns its pool (`my_brep`, the
//!    BRepOffsetAPI_ThruSections facade precedent).
//! 4. OCCT enums BRepBuilderAPI_PipeError (BRepBuilderAPI_PipeError.hxx)
//!    is the local enum below; GeomFill_PipeError is the landed
//!    crate::geomalgo::geomfill::trihedron_law::PipeError and
//!    BRepFill_TypeOfContact / BRepFill_TransitionStyle are the landed
//!    crate::brep_fill::brep_fill_pipe_shell_b enums.
//! 5. const occ::handle<Law_Function>& L -> &LawFunctionHandle (the landed
//!    geomalgo::law::law_function handle type).

use rcad_kernel::math::gp::Ax2;
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;

use glam::DVec3;

use crate::brep_fill::brep_fill_pipe_shell::BRepFillPipeShell;
use crate::brep_fill::brep_fill_pipe_shell_b::{BRepFillTransitionStyle, BRepFillTypeOfContact};
use crate::geomalgo::geomfill::trihedron_law::PipeError;
use crate::geomalgo::law::law_function::LawFunctionHandle;

use super::brep_offset_api_make_draft::BRepBuilderAPITransitionMode;

// ---------------------------------------------------------------------------
// OCCT enum (architecture difference #5).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_PipeError (BRepBuilderAPI_PipeError.hxx L25-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepBuilderAPIPipeError {
    PipeDone,
    PipeNotDone,
    PlaneNotIntersectGuide,
    ImpossibleContact,
}

/// OCCT BRepOffsetAPI_MakePipeShell (hxx L91-395).
pub struct BRepOffsetAPIMakePipeShell {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // The rcad arena stand-in (architecture difference #3; the
    // BRepOffsetAPI_ThruSections facade precedent).
    my_brep: BRep, // rcad pool (arch. diff. #4)
    // OCCT private member (hxx L393): Handle(BRepFill_PipeShell) myPipe.
    my_pipe: BRepFillPipeShell, // OCCT: myPipe
}

impl BRepOffsetAPIMakePipeShell {
    /// OCCT BRepOffsetAPI_MakePipeShell::BRepOffsetAPI_MakePipeShell(Spine)
    /// (cxx L28-34).
    pub fn new(the_spine: &Shape) -> Self {
        // OCCT L30: myPipe = new (BRepFill_PipeShell)(Spine) — the rcad pool
        // (arch. diff. #4) is the arena the engine mutates.
        let mut my_brep = BRep::new();
        let my_pipe = BRepFillPipeShell::new(&mut my_brep, the_spine);
        let mut r = BRepOffsetAPIMakePipeShell {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_brep,
            my_pipe,
        };
        // OCCT L31: SetTolerance().
        r.set_tolerance(1.0e-4, 1.0e-4, 1.0e-2);
        // OCCT L32: SetTransitionMode().
        r.set_transition_mode(BRepBuilderAPITransitionMode::TransitionMode_Transformed);
        // OCCT L33: NotDone().
        r.my_done = false;
        r
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(IsFrenet) (cxx L36-39):
    /// myPipe->Set(IsFrenet).
    pub fn set_mode(&mut self, is_frenet: bool) {
        self.my_pipe.set(&self.my_brep, is_frenet);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetDiscreteMode() (cxx L41-44):
    /// myPipe->SetDiscrete().
    pub fn set_discrete_mode(&mut self) {
        self.my_pipe.set_discrete(&self.my_brep);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(Axe) (cxx L46-49):
    /// myPipe->Set(Axe).
    pub fn set_mode_with_axe(&mut self, axe: &Ax2) {
        self.my_pipe.set_with_axe(&self.my_brep, *axe);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(BiNormal) (cxx L51-54):
    /// myPipe->Set(BiNormal).
    pub fn set_mode_with_binormal(&mut self, binormal: &DVec3) {
        self.my_pipe.set_with_binormal(&self.my_brep, *binormal);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(SpineSupport)
    /// (cxx L56-59): return myPipe->Set(SpineSupport).
    pub fn set_mode_with_support(&mut self, spine_support: &Shape) -> bool {
        self.my_pipe.set_spine_support(&self.my_brep, spine_support)
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(AuxiliarySpine,
    /// CurvilinearEquivalence, KeepContact) (cxx L61-67):
    /// myPipe->Set(AuxiliarySpine, CurvilinearEquivalence, KeepContact).
    pub fn set_mode_with_auxiliary_spine(
        &mut self,
        auxiliary_spine: &Shape,
        curvilinear_equivalence: bool,
        keep_contact: BRepFillTypeOfContact,
    ) {
        self.my_pipe.set_with_auxiliary_spine(
            &mut self.my_brep,
            auxiliary_spine,
            curvilinear_equivalence,
            keep_contact,
        );
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Add(Profile, WithContact,
    /// WithCorrection) (cxx L69-74): myPipe->Add(Profile, WithContact,
    /// WithCorrection).
    pub fn add(&mut self, profile: &Shape, with_contact: bool, with_correction: bool) {
        self.my_pipe
            .add(&mut self.my_brep, profile, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Add(Profile, Location, WithContact,
    /// WithCorrection) (cxx L76-82): myPipe->Add(Profile, Location,
    /// WithContact, WithCorrection).
    pub fn add_with_location(
        &mut self,
        profile: &Shape,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe.add_with_location(
            &mut self.my_brep,
            profile,
            location,
            with_contact,
            with_correction,
        );
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetLaw(Profile, L, WithContact,
    /// WithCorrection) (cxx L84-90): myPipe->SetLaw(Profile, L, WithContact,
    /// WithCorrection).
    pub fn set_law(
        &mut self,
        profile: &Shape,
        l: &LawFunctionHandle,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe.set_law(profile, l, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetLaw(Profile, L, Location,
    /// WithContact, WithCorrection) (cxx L92-99).
    pub fn set_law_with_location(
        &mut self,
        profile: &Shape,
        l: &LawFunctionHandle,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe
            .set_law_with_location(profile, l, location, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Delete(Profile) (cxx L101-105):
    /// myPipe->DeleteProfile(Profile).
    pub fn delete(&mut self, profile: &Shape) {
        self.my_pipe.delete_profile(profile);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::IsReady() (cxx L107-111).
    pub fn is_ready(&self) -> bool {
        self.my_pipe.is_ready()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::GetStatus() (cxx L113-134).
    pub fn get_status(&self) -> BRepBuilderAPIPipeError {
        // OCCT L116-117: GeomFill_PipeError stat = myPipe->GetStatus().
        let stat = self.my_pipe.get_status();
        // OCCT L118-133: the switch.
        match stat {
            PipeError::PipeOk => BRepBuilderAPIPipeError::PipeDone,
            PipeError::PlaneNotIntersectGuide => BRepBuilderAPIPipeError::PlaneNotIntersectGuide,
            PipeError::ImpossibleContact => BRepBuilderAPIPipeError::ImpossibleContact,
            _ => BRepBuilderAPIPipeError::PipeNotDone,
        }
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetTolerance(Tol3d, BoundTol,
    /// TolAngular) (cxx L136-141).
    pub fn set_tolerance(&mut self, tol3d: f64, bound_tol: f64, tol_angular: f64) {
        self.my_pipe.set_tolerance(tol3d, bound_tol, tol_angular);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMaxDegree(NewMaxDegree)
    /// (cxx L143-147).
    pub fn set_max_degree(&mut self, new_max_degree: i32) {
        self.my_pipe.set_max_degree(new_max_degree);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMaxSegments(NewMaxSegments)
    /// (cxx L149-153).
    pub fn set_max_segments(&mut self, new_max_segments: i32) {
        self.my_pipe.set_max_segments(new_max_segments);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetForceApproxC1(ForceApproxC1)
    /// (cxx L155-165).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_pipe.set_force_approx_c1(force_approx_c1);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetTransitionMode(Mode)
    /// (cxx L167-171): myPipe->SetTransition((BRepFill_TransitionStyle)Mode).
    /// The C++ cast is by ordinal: BRepBuilderAPI_Transformed(0) ->
    /// BRepFill_Modified(0); BRepBuilderAPI_RightCorner(1) -> BRepFill_Right
    /// (1) — the rcad engine enum spells BRepFill_Right "RightCorner"
    /// (BRepFill_Sweep.cxx L3914 `Transition == BRepFill_Right` maps to
    /// brep_fill_sweep.rs L1404 `transition == RightCorner`);
    /// BRepBuilderAPI_RoundCorner(2) -> BRepFill_Round(2).  The engine
    /// SetTransition default arguments (BRepFill_PipeShell.hxx L176-179:
    /// Angmin = 1.0e-2, Angmax = 6.0) are stated explicitly.
    pub fn set_transition_mode(&mut self, mode: BRepBuilderAPITransitionMode) {
        let style = match mode {
            BRepBuilderAPITransitionMode::TransitionMode_Transformed => {
                BRepFillTransitionStyle::Modified
            }
            BRepBuilderAPITransitionMode::RightCorner => BRepFillTransitionStyle::RightCorner,
            BRepBuilderAPITransitionMode::RoundCorner => BRepFillTransitionStyle::Round,
        };
        self.my_pipe.set_transition(style, 1.0e-2, 6.0);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Simulate(N, R) (cxx L173-177):
    /// myPipe->Simulate(N, R).
    pub fn simulate(&mut self, n: i32, r: &mut Vec<Shape>) {
        self.my_pipe.simulate(&mut self.my_brep, n, r);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Build(...) (cxx L179-193).
    pub fn build(&mut self) {
        // OCCT L181-182: bool Ok; Ok = myPipe->Build().
        let ok = self.my_pipe.build(&mut self.my_brep);
        if ok {
            // OCCT L184-186: myShape = myPipe->Shape(); Done().
            self.my_shape = self.my_pipe.shape();
            self.my_done = true;
        } else {
            // OCCT L188-191: NotDone().
            self.my_done = false;
        }
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::MakeSolid() (cxx L195-208).
    pub fn make_solid(&mut self) -> bool {
        // OCCT L197-199: StdFail_NotDone_Raise_if(!IsDone(), ...).
        assert!(
            self.my_done,
            "BRepOffsetAPI_MakePipeShell::MakeSolid"
        );
        // OCCT L200-201: bool Ok; Ok = myPipe->MakeSolid().
        let ok = self.my_pipe.make_solid(&mut self.my_brep);
        if ok {
            // OCCT L203: myShape = myPipe->Shape().
            self.my_shape = self.my_pipe.shape();
        }
        // OCCT L207: return Ok.
        ok
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::FirstShape() (cxx L210-213).
    pub fn first_shape(&mut self) -> Shape {
        self.my_pipe.first_shape()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::LastShape() (cxx L215-218).
    pub fn last_shape(&mut self) -> Shape {
        self.my_pipe.last_shape()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Generated(S) (cxx L220-224).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L222: myPipe->Generated(S, myGenerated).
        self.my_generated = self.my_pipe.generated(s);
        // OCCT L223: return myGenerated.
        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::ErrorOnSurface() (cxx L226-230).
    pub fn error_on_surface(&self) -> f64 {
        self.my_pipe.error_on_surface()
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
    /// Shape() const; raises StdFail_NotDone when not done).  The engine-level
    /// BRepFillPipeShell::shape() accessor is the OCCT
    /// BRepFill_PipeShell::Shape() member.
    pub fn shape(&self) -> Shape {
        assert!(
            self.my_done,
            "StdFail_NotDone: BRepOffsetAPI_MakePipeShell::Shape()"
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
