//! OCCT BRepOffsetAPI_MakeEvolved — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeEvolved.cxx (L28-149) +
//!         BRepOffsetAPI_MakeEvolved.hxx (L76-151).
//!
//! OCCT inheritance chain (hxx L76): BRepOffsetAPI_MakeEvolved ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated, the Done flag) are kept as plain fields of
//! the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>; OCCT's file-static
//!    anEmptyList (cxx L28) maps to the `Vec::new()` of the not-done branch
//!    (the OCCT reference-to-static-list form cannot be carried).
//! 2. The engine member BRepFill_Evolved myEvolved (TKBool/BRepFill) is the
//!    real crate::brep_fill::brep_fill_evolved translation (its single 1:1
//!    body; the local GAP carrier was retired). The member
//!    BRepFill_AdvancedEvolved myVolume (TKBool/BRepFill) has no rcad
//!    translation yet — the carrier below keeps the OCCT
//!    constructor/method surface with GAP panics (port plan section 0.6);
//!    every facade body around the engine calls is translated 1:1.
//! 3. The static BRepFill::Axe(Spine, Profil, Axis, POS, Tol)
//!    (TKBool/BRepFill, BRepFill.cxx) is the GAP fn below.
//! 4. TopoDS_Iterator(Spine).Value() maps to the single-child wire read of
//!    the face (the tool.rs sub_shapes carrier).

use rcad_kernel::math::gp::Ax3;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;

use crate::brep_fill::offset_wire::GeomAbsJoinType;
use crate::brep_algo::tool::sub_shapes;
// OCCT BRepFill_Evolved (TKBool/BRepFill/BRepFill_Evolved.hxx / .cxx) — the
// single 1:1 body lives in crate::brep_fill::brep_fill_evolved (its correct
// location per the module map); the local GAP carrier of this file was
// retired.
use crate::brep_fill::brep_fill_evolved::BRepFillEvolved;

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #2/#3).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_AdvancedEvolved (TKBool/BRepFill) — the BOPAlgo_
/// MakerVolume-driven engine of the theIsVolume branch (architecture
/// difference #2; GAP: no rcad translation yet — the GAP panics are the
/// section 0.6 annotation).
pub struct BRepFillAdvancedEvolved {
    my_is_done: bool, // OCCT: myIsDone
    my_shape: Shape,  // OCCT: myShape
}

impl BRepFillAdvancedEvolved {
    /// OCCT BRepFill_AdvancedEvolved::BRepFill_AdvancedEvolved().
    pub fn new() -> Self {
        BRepFillAdvancedEvolved {
            my_is_done: false,
            my_shape: Shape::null(),
        }
    }

    /// OCCT BRepFill_AdvancedEvolved::SetParallelMode(theRunInParallel) —
    /// GAP.
    pub fn set_parallel_mode(&mut self, _the_run_in_parallel: bool) {
        panic!("GAP: BRepFill_AdvancedEvolved::SetParallelMode (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_AdvancedEvolved::Perform(Spine, Profil, Tol, Solid) —
    /// GAP.
    pub fn perform(&mut self, _the_spine: &Shape, _the_profil: &Shape, _the_tol: f64, _the_solid: bool) {
        panic!("GAP: BRepFill_AdvancedEvolved::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_AdvancedEvolved::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_AdvancedEvolved::Shape().
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }
}

impl Default for BRepFillAdvancedEvolved {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT static BRepFill::Axe(Spine, Profil, Axe, POS, Tol)
/// (TKBool/BRepFill, BRepFill.cxx) — computes the profile axis; GAP.
pub fn brep_fill_axe(
    _the_spine: &Shape,
    _the_profil: &Shape,
    _the_axe: &mut Ax3,
    _the_pos: &mut bool,
    _the_tol: f64,
) {
    panic!("GAP: BRepFill::Axe (TKBool/BRepFill not translated)");
}

/// OCCT BRepOffsetAPI_MakeEvolved (hxx L76-151).
pub struct BRepOffsetAPIMakeEvolved {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private members (hxx L147-150).
    my_evolved: BRepFillEvolved,           // OCCT: myEvolved
    my_volume: BRepFillAdvancedEvolved,    // OCCT: myVolume
    my_is_volume: bool,                    // OCCT: myIsVolume
}

impl Default for BRepOffsetAPIMakeEvolved {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetAPIMakeEvolved {
    /// OCCT BRepOffsetAPI_MakeEvolved::BRepOffsetAPI_MakeEvolved() (cxx L30).
    pub fn new() -> Self {
        BRepOffsetAPIMakeEvolved {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_evolved: BRepFillEvolved::new(),
            my_volume: BRepFillAdvancedEvolved::new(),
            my_is_volume: false,
        }
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::BRepOffsetAPI_MakeEvolved(theSpine,
    /// theProfile, theJoinType, theIsAxeProf, theIsSolid, theIsProfOnSpine,
    /// theTol, theIsVolume, theRunInParallel) (cxx L32-88).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_args(
        the_spine: &Shape,
        the_profil: &Shape,
        the_join_type: GeomAbsJoinType,
        the_is_axe_prof: bool,
        the_is_solid: bool,
        the_is_prof_on_spine: bool,
        the_tol: f64,
        the_is_volume: bool,
        the_run_in_parallel: bool,
    ) -> Self {
        // OCCT L41: myIsVolume(theIsVolume).
        let mut r = BRepOffsetAPIMakeEvolved {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_evolved: BRepFillEvolved::new(),
            my_volume: BRepFillAdvancedEvolved::new(),
            my_is_volume: the_is_volume,
        };
        // OCCT L43-46: the spine type guard.
        if the_spine.shape_type() != ShapeType::Wire && the_spine.shape_type() != ShapeType::Face {
            panic!("BRepOffsetAPI_MakeEvolved: face or wire is expected as a spine");
        }
        // OCCT L47-60: the volume branch.
        if r.my_is_volume {
            // OCCT L49.
            r.my_volume.set_parallel_mode(the_run_in_parallel);
            // OCCT L50-56: TopoDS_Wire aSpine — the wire itself, or the first
            // child of the face (TopoDS_Iterator(Spine).Value()).
            let a_spine = if the_spine.shape_type() == ShapeType::Wire {
                the_spine.clone()
            } else {
                let children = sub_shapes(the_spine);
                children
                    .first()
                    .cloned()
                    .unwrap_or_else(Shape::null)
            };
            // OCCT L57.
            r.my_volume.perform(&a_spine, the_profil, the_tol, the_is_solid);
            // OCCT L58-61.
            if !r.my_volume.is_done() {
                return r;
            }
        } else {
            // OCCT L64: gp_Ax3 Axis(gp::Origin(), gp::DZ(), gp::DX()).
            let mut axis = Ax3::from_pnt_n_vx(glam::DVec3::ZERO, glam::DVec3::Z, glam::DVec3::X);

            // OCCT L66-75.
            if !the_is_axe_prof {
                let mut pos = false;
                // OCCT L69: BRepFill::Axe(Spine, Profil, Axis, POS,
                // max(Tol, Precision::Confusion())).
                brep_fill_axe(
                    the_spine,
                    the_profil,
                    &mut axis,
                    &mut pos,
                    the_tol.max(CONFUSION),
                );
                // OCCT L70-74.
                if the_is_prof_on_spine && !pos {
                    return r;
                }
            }
            // OCCT L76-84.
            if the_spine.shape_type() == ShapeType::Wire {
                r.my_evolved.perform_with_wire_spine(
                    the_spine,
                    the_profil,
                    &axis,
                    the_join_type,
                    the_is_solid,
                );
            } else {
                r.my_evolved.perform_with_face_spine(
                    the_spine,
                    the_profil,
                    &axis,
                    the_join_type,
                    the_is_solid,
                );
            }
        }

        // OCCT L86: Build().
        r.build();
        r
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::Evolved() (cxx L90-100).
    pub fn evolved(&self) -> &BRepFillEvolved {
        if self.my_is_volume {
            panic!(
                "BRepOffsetAPI_MakeEvolved: myEvolved is accessed while in volume mode"
            );
        }
        &self.my_evolved
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::Build(...) (cxx L102-113).
    pub fn build(&mut self) {
        // OCCT L104-107.
        if self.my_evolved.is_done() {
            self.my_shape = self.my_evolved.shape().clone();
        } else if self.my_volume.is_done() {
            self.my_shape = self.my_volume.shape();
        }

        // OCCT L111: Done().
        self.my_done = true;
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::Top() (cxx L115-119).
    pub fn top(&self) -> &Shape {
        self.my_evolved.top()
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::Bottom() (cxx L121-125).
    pub fn bottom(&self) -> &Shape {
        self.my_evolved.bottom()
    }

    /// OCCT BRepOffsetAPI_MakeEvolved::GeneratedShapes(SpineShape, ProfShape)
    /// (cxx L127-140).
    pub fn generated_shapes(&self, spine_shape: &Shape, prof_shape: &Shape) -> Vec<Shape> {
        // OCCT L129-132: if (!myEvolved.IsDone()) return anEmptyList.
        if !self.my_evolved.is_done() {
            return Vec::new();
        }

        // OCCT L134: return myEvolved.GeneratedShapes(SpineShape, ProfShape).
        self.my_evolved.generated_shapes(spine_shape, prof_shape)
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
            "StdFail_NotDone: BRepOffsetAPI_MakeEvolved::Shape()"
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
