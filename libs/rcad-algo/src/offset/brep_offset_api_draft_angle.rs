//! OCCT BRepOffsetAPI_DraftAngle — 1:1 translation (class + build methods).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_DraftAngle.cxx (L55-229) +
//!         BRepOffsetAPI_DraftAngle.hxx (L60-205).
//! CorrectWires / CorrectVertexTol (cxx L231-992) live in
//! brep_offset_api_draft_angle_b.rs.
//!
//! OCCT inheritance chain (hxx L60): BRepOffsetAPI_DraftAngle ->
//! BRepBuilderAPI_ModifyShape -> BRepBuilderAPI_MakeShape.  Rust has no
//! inheritance: the base members (myShape, myGenerated, myDone,
//! myInitialShape, myModifier) are plain fields; the two DoModif drivers
//! (BRepBuilderAPI_ModifyShape.cxx L31-66) are the `do_modif` /
//! `do_modif_shape` fns below.
//!
//! Architecture differences:
//! 1. NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> myVtxToReplace ->
//!    HashMap<u64, Shape> keyed by the TShape identity (the
//!    TopTools_ShapeMapHasher key).
//! 2. BRepTools_ReShape mySubs -> the landed
//!    shhealing::shape_build::reshape::ShapeBuildReShape (it consumes the
//!    rcad BRep pool — the class holds my_brep, architecture difference #4;
//!    the OCCT const accessors take &mut self, the
//!    brep_offset_make_simple_offset.rs precedent).
//! 3. Handle(BRepTools_Modification) myModification holds the
//!    Draft_Modification down-cast target -> Option<DraftModification> (the
//!    landed offset/draft_modification.rs translation); the OCCT
//!    down_cast sites are direct field accesses.
//! 4. BRepTools_Modifier myModifier (TKTopAlgo/BRepTools) has no rcad
//!    translation of its Perform/ModifiedShape engine yet — the
//!    BRepToolsModifierForDraft carrier below keeps the OCCT method surface
//!    with GAP panics (port plan section 0.6); every facade body around the
//!    engine calls is translated 1:1.
//! 5. gp_Dir -> DVec3; gp_Pln -> rcad_kernel::geom::Plane.

use std::collections::HashMap;

use rcad_kernel::geom::Plane;
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;

use glam::DVec3;

use crate::offset::draft_error_status::DraftErrorStatus;
use crate::offset::draft_modification::DraftModification;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

// ---------------------------------------------------------------------------
// GAP carrier (architecture difference #4).
// ---------------------------------------------------------------------------

/// OCCT BRepTools_Modifier (TKTopAlgo/BRepTools, BRepTools_Modifier.hxx L44+)
/// — the shape rebuild vehicle driven by a BRepTools_Modification
/// (architecture difference #4; GAP: the Perform/ModifiedShape engine has no
/// rcad translation yet — the GAP panics are the section 0.6 annotation;
/// Init keeps the OCCT storage form; the BRepToolsModifier carrier of
/// brep_offset_make_simple_offset.rs is the same-shape precedent typed to
/// BRepOffset_SimpleOffset).
pub struct BRepToolsModifierForDraft {
    my_shape: Shape,  // OCCT: myShape
    my_is_done: bool, // OCCT: myIsDone
}

impl BRepToolsModifierForDraft {
    /// OCCT BRepTools_Modifier::BRepTools_Modifier(S) (hxx L50).
    pub fn new(the_s: &Shape) -> Self {
        BRepToolsModifierForDraft {
            my_shape: the_s.clone(),
            my_is_done: false,
        }
    }

    /// OCCT BRepTools_Modifier::Init(S).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_is_done = false;
    }

    /// OCCT BRepTools_Modifier::Perform(M) — GAP.
    pub fn perform(&mut self, the_m: &DraftModification) {
        let _ = the_m;
        panic!("GAP: BRepTools_Modifier::Perform (TKTopAlgo/BRepTools not translated)")
    }

    /// OCCT BRepTools_Modifier::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepTools_Modifier::ModifiedShape(S) — GAP.
    pub fn modified_shape(&self, the_s: &Shape) -> Shape {
        let _ = the_s;
        panic!("GAP: BRepTools_Modifier::ModifiedShape (TKTopAlgo/BRepTools not translated)")
    }
}

/// OCCT BRepOffsetAPI_DraftAngle (hxx L60-205).
pub struct BRepOffsetAPIDraftAngle {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    pub(crate) my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT BRepBuilderAPI_ModifyShape members.
    pub(crate) my_initial_shape: Shape, // OCCT: myInitialShape
    my_modifier: BRepToolsModifierForDraft, // OCCT: myModifier
    // The rcad arena stand-in consumed by the ReShape engine (architecture
    // difference #2).
    pub(crate) my_brep: BRep, // rcad pool (arch. diff. #4)
    // OCCT private members (hxx L199-203).
    my_modification: Option<DraftModification>, // OCCT: myModification (the Draft_Modification handle)
    pub(crate) my_vtx_to_replace: HashMap<u64, (Shape, Shape)>, // OCCT: myVtxToReplace (old shape carried for the DataMap Key()/Value() reads)
    pub(crate) my_subs: ShapeBuildReShape,     // OCCT: mySubs (BRepTools_ReShape)
}

impl Default for BRepOffsetAPIDraftAngle {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetAPIDraftAngle {
    /// OCCT BRepOffsetAPI_DraftAngle::BRepOffsetAPI_DraftAngle() (cxx L55).
    pub fn new() -> Self {
        BRepOffsetAPIDraftAngle {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_initial_shape: Shape::null(),
            my_modifier: BRepToolsModifierForDraft::new(&Shape::null()),
            my_brep: BRep::new(),
            my_modification: None,
            my_vtx_to_replace: HashMap::new(),
            my_subs: ShapeBuildReShape::new(),
        }
    }

    /// OCCT BRepOffsetAPI_DraftAngle::BRepOffsetAPI_DraftAngle(S) (cxx
    /// L57-61).
    pub fn new_with_shape(s: &Shape) -> Self {
        let mut r = Self::new();
        // OCCT L59-60.
        r.my_initial_shape = s.clone();
        r.my_modification = Some(DraftModification::new(s));
        r
    }

    /// OCCT BRepBuilderAPI_ModifyShape::DoModif()
    /// (BRepBuilderAPI_ModifyShape.cxx L31-45).
    fn do_modif(&mut self) {
        // OCCT L33-36.
        if self.my_initial_shape.is_null() || self.my_modification.is_none() {
            panic!("Standard_NullObject");
        }
        // OCCT L37: myModifier.Perform(myModification).
        self.my_modifier
            .perform(self.my_modification.as_ref().unwrap());
        // OCCT L38-45.
        if self.my_modifier.is_done() {
            // OCCT L40-41: Done(); myShape = myModifier.ModifiedShape(...).
            self.my_done = true;
            let initial = self.my_initial_shape.clone();
            self.my_shape = self.my_modifier.modified_shape(&initial);
        } else {
            // OCCT L43-44.
            self.my_done = false;
        }
    }

    /// OCCT BRepBuilderAPI_ModifyShape::DoModif(S)
    /// (BRepBuilderAPI_ModifyShape.cxx L47-53).
    fn do_modif_shape(&mut self, s: &Shape) {
        // OCCT L49-52.
        if !s.is_equal(&self.my_initial_shape) || !self.my_done {
            self.my_initial_shape = s.clone();
            self.my_modifier.init(s);
            self.do_modif();
        }
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Clear() (cxx L63-68).
    pub fn clear(&mut self) {
        if let Some(dmod) = self.my_modification.as_mut() {
            dmod.clear();
        }
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Init(S) (cxx L70-79).
    pub fn init(&mut self, s: &Shape) {
        // OCCT L72-73.
        self.my_initial_shape = s.clone();
        self.my_done = false;
        if let Some(dmod) = self.my_modification.as_mut() {
            dmod.init(s);
        } else {
            self.my_modification = Some(DraftModification::new(s));
        }
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Add(F, Direction, Angle, NeutralPlane,
    /// Flag) (cxx L81-93).
    pub fn add(
        &mut self,
        f: &Shape,
        direction: &DVec3,
        angle: f64,
        neutral_plane: &Plane,
        flag: bool,
    ) {
        // OCCT L84-87: POP-DPF : protection.
        if angle.abs() <= 1.0e-04 {
            return;
        }
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::Add() - initial shape is not set"
        );
        // OCCT L91-92.
        self.my_modification
            .as_mut()
            .unwrap()
            .add(f, *direction, angle, neutral_plane, flag);
    }

    /// OCCT BRepOffsetAPI_DraftAngle::AddDone() (cxx L95-100).
    pub fn add_done(&self) -> bool {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::AddDone() - initial shape is not set"
        );
        // OCCT L99: ProblematicShape().IsNull().
        self.my_modification
            .as_ref()
            .unwrap()
            .problematic_shape()
            .is_null()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Remove(F) (cxx L102-107).
    pub fn remove(&mut self, f: &Shape) {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::Remove() - initial shape is not set"
        );
        self.my_modification.as_mut().unwrap().remove(f);
    }

    /// OCCT BRepOffsetAPI_DraftAngle::ProblematicShape() (cxx L109-115).
    pub fn problematic_shape(&self) -> &Shape {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::ProblematicShape() - initial shape is not set"
        );
        self.my_modification.as_ref().unwrap().problematic_shape()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Status() (cxx L117-122).
    pub fn status(&self) -> DraftErrorStatus {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::Status() - initial shape is not set"
        );
        self.my_modification.as_ref().unwrap().error()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::ConnectedFaces(F) (cxx L124-130).
    pub fn connected_faces(&mut self, f: &Shape) -> &Vec<Shape> {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::ConnectedFaces() - initial shape is not set"
        );
        self.my_modification.as_mut().unwrap().connected_faces(f)
    }

    /// OCCT BRepOffsetAPI_DraftAngle::ModifiedFaces() (cxx L132-138).
    pub fn modified_faces(&mut self) -> &Vec<Shape> {
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::ModifiedFaces() - initial shape is not set"
        );
        self.my_modification.as_mut().unwrap().modified_faces()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Generated(S) (cxx L140-167).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L141: myGenerated.Clear().
        self.my_generated.clear();
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::Generated() - initial shape is not set"
        );
        // OCCT L147-153: if (S.ShapeType() == TopAbs_FACE).
        if s.shape_type() == ShapeType::Face {
            // OCCT L148-153: the NewSurface output carriers.
            let mut surf: Option<rcad_kernel::geom::Surface3> = None;
            let mut l: u32 = 0;
            let mut tol: f64 = 0.0;
            let mut rw: bool = false;
            let mut rf: bool = false;
            let dmod = self.my_modification.as_mut().unwrap();
            let has_new_surface =
                dmod.new_surface(s, &mut surf, &mut l, &mut tol, &mut rw, &mut rf);
            if has_new_surface {
                // OCCT L156-163.
                if self.my_vtx_to_replace.is_empty() {
                    let modified = self.modified_shape(s);
                    self.my_generated.push(modified);
                } else {
                    let modified = self.modified_shape(s);
                    let value = self.my_subs.value(&mut self.my_brep, &modified);
                    self.my_generated.push(value);
                }
            }
        }
        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Modified(S) (cxx L169-198).
    pub fn modified(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L170: myGenerated.Clear().
        self.my_generated.clear();
        assert!(
            !self.my_initial_shape.is_null(),
            "BRepOffsetAPI_DraftAngle::Modified() - initial shape is not set"
        );
        // OCCT L175-190: if (S.ShapeType() == TopAbs_FACE).
        if s.shape_type() == ShapeType::Face {
            let mut surf: Option<rcad_kernel::geom::Surface3> = None;
            let mut l: u32 = 0;
            let mut tol: f64 = 0.0;
            let mut rw: bool = false;
            let mut rf: bool = false;
            let dmod = self.my_modification.as_mut().unwrap();
            let has_new_surface =
                dmod.new_surface(s, &mut surf, &mut l, &mut tol, &mut rw, &mut rf);
            if !has_new_surface {
                // OCCT L181: Ce n est pas une generation => peut etre une
                // modif.
                if self.my_vtx_to_replace.is_empty() {
                    let modified = self.modified_shape(s);
                    self.my_generated.push(modified);
                } else {
                    let modified = self.modified_shape(s);
                    let value = self.my_subs.value(&mut self.my_brep, &modified);
                    self.my_generated.push(value);
                }
                // OCCT L191-194.
                if self.my_generated.len() == 1 && self.my_generated[0].is_same(s) {
                    self.my_generated.clear();
                }
            }
        }
        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_DraftAngle::ModifiedShape(S) (cxx L200-216).
    pub fn modified_shape(&mut self, s: &Shape) -> Shape {
        // OCCT L202-206.
        if s.shape_type() == ShapeType::Vertex {
            if let Some((_, replacement)) = self.my_vtx_to_replace.get(&s.ptr_id()) {
                return replacement.clone();
            }
        }
        if self.my_vtx_to_replace.is_empty() {
            // OCCT L210-211.
            self.my_modifier.modified_shape(s)
        } else {
            // OCCT L212-215: mySubs.Value(myModifier.ModifiedShape(S)) — the
            // ReShape::Value consumes the rcad pool (architecture difference
            // #2; the OCCT const accessor takes &mut self on the carrier).
            let a_ns = self.my_modifier.modified_shape(s);
            self.my_subs.value(&mut self.my_brep, &a_ns)
        }
    }

    /// OCCT BRepOffsetAPI_DraftAngle::Build(...) (cxx L218-230).
    pub fn build(&mut self) {
        // OCCT L220: Perform() on the Draft_Modification.
        self.my_modification.as_mut().unwrap().perform();
        // OCCT L221-229.
        if !self.my_modification.as_ref().unwrap().is_done() {
            self.my_done = false;
        } else {
            self.do_modif_shape(&self.my_initial_shape.clone());
            self.correct_wires();
            self.correct_vertex_tol();
        }
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
            "StdFail_NotDone: BRepOffsetAPI_DraftAngle::Shape()"
        );
        self.my_shape.clone()
    }

    /// Test-world extraction bridge (the FilletResult.brep pattern,
    /// algo_ext::topods_ext::extract_result_brep): flattens the
    /// (my_brep arena, my_shape root) pair into the self-contained BRep the
    /// test world consumes (StepWriter / total_surface_area).  OCCT has no
    /// equivalent (the TopoDS_Shape carries its arena implicitly) — the
    /// rcad BRep-pool architecture difference #4 glue.
    pub fn result_brep(&mut self) -> Option<BRep> {
        if !self.my_done || self.my_shape.is_null() {
            return None;
        }
        let locations = self.my_brep.locations.clone();
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            locations,
        ))
    }
}
