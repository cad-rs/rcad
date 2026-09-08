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
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>; the file-static
//!    anEmptyList (cxx L28) is the static empty Vec below (GeneratedShapes
//!    returns a slice; the OCCT reference-to-static-list form cannot be
//!    carried).
//! 2. The engine members BRepFill_Evolved myEvolved and
//!    BRepFill_AdvancedEvolved myVolume (TKBool/BRepFill) have no rcad
//!    translation yet — the two carriers below keep the OCCT
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

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #2/#3).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Evolved (TKBool/BRepFill, BRepFill_Evolved.hxx) — the
/// evolved-sweep engine of MakeEvolved (architecture difference #2; GAP: no
/// rcad translation yet — the GAP panics are the section 0.6 annotation; the
/// storage form is kept).
pub struct BRepFillEvolved {
    my_is_done: bool, // OCCT: myIsDone
    my_top: Shape,    // OCCT: myTop
    my_bottom: Shape, // OCCT: myBottom
    my_shape: Shape,  // OCCT: myShape
}

impl BRepFillEvolved {
    /// OCCT BRepFill_Evolved::BRepFill_Evolved().
    pub fn new() -> Self {
        BRepFillEvolved {
            my_is_done: false,
            my_top: Shape::null(),
            my_bottom: Shape::null(),
            my_shape: Shape::null(),
        }
    }

    /// OCCT BRepFill_Evolved::Perform(Spine, Profil, Axis, Join, Solid) —
    /// GAP.
    pub fn perform_with_wire(
        &mut self,
        _the_spine: &Shape,
        _the_profil: &Shape,
        _the_axis: &Ax3,
        _the_join: GeomAbsJoinType,
        _the_solid: bool,
    ) {
        panic!("GAP: BRepFill_Evolved::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Evolved::Perform(Spine, Profil, Axis, Join, Solid) —
    /// the face-spine form; GAP.
    pub fn perform_with_face(
        &mut self,
        _the_spine: &Shape,
        _the_profil: &Shape,
        _the_axis: &Ax3,
        _the_join: GeomAbsJoinType,
        _the_solid: bool,
    ) {
        panic!("GAP: BRepFill_Evolved::Perform (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Evolved::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_Evolved::Shape().
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_Evolved::Top().
    pub fn top(&self) -> &Shape {
        &self.my_top
    }

    /// OCCT BRepFill_Evolved::Bottom().
    pub fn bottom(&self) -> &Shape {
        &self.my_bottom
    }

    /// OCCT BRepFill_Evolved::GeneratedShapes(SpineShape, ProfShape) — GAP.
    pub fn generated_shapes(
        &self,
        _the_spine_shape: &Shape,
        _the_prof_shape: &Shape,
    ) -> &Vec<Shape> {
        panic!("GAP: BRepFill_Evolved::GeneratedShapes (TKBool/BRepFill not translated)")
    }
}

impl Default for BRepFillEvolved {
    fn default() -> Self {
        Self::new()
    }
}

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

/// OCCT file-static anEmptyList (cxx L28) — the empty GeneratedShapes result.
static AN_EMPTY_LIST: std::sync::OnceLock<Vec<Shape>> = std::sync::OnceLock::new();

fn an_empty_list() -> &'static Vec<Shape> {
    AN_EMPTY_LIST.get_or_init(Vec::new)
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
                r.my_evolved
                    .perform_with_wire(the_spine, the_profil, &axis, the_join_type, the_is_solid);
            } else {
                r.my_evolved
                    .perform_with_face(the_spine, the_profil, &axis, the_join_type, the_is_solid);
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
            self.my_shape = self.my_evolved.shape();
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
    pub fn generated_shapes(&self, spine_shape: &Shape, prof_shape: &Shape) -> &Vec<Shape> {
        // OCCT L129-132: if (!myEvolved.IsDone()) return anEmptyList.
        if !self.my_evolved.is_done() {
            return an_empty_list();
        }

        // OCCT L134: return myEvolved.GeneratedShapes(SpineShape, ProfShape).
        self.my_evolved.generated_shapes(spine_shape, prof_shape)
    }
}
