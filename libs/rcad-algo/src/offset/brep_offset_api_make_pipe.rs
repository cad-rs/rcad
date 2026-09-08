//! OCCT BRepOffsetAPI_MakePipe — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakePipe.cxx (L26-133) +
//!         BRepOffsetAPI_MakePipe.hxx (L54-100).
//!
//! OCCT inheritance chain (hxx L54): BRepOffsetAPI_MakePipe ->
//! BRepPrimAPI_MakeSweep -> BRepBuilderAPI_MakeShape.  Rust has no
//! inheritance: the base-class members (myShape, myGenerated, the Done flag)
//! are kept as plain fields of the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>;
//!    NCollection_IndexedMap<TopoDS_Shape> -> the occurrence Vec of the
//!    TopExp::MapShapes walk (only the extent is consumed, cxx L72-81).
//! 2. The engine member BRepFill_Pipe myPipe is the D1-approved translation
//!    feat/loc_ope_pipe.rs (TKBool/BRepFill).  Its 4-argument constructor
//!    (aMode, ForceApproxC1) and the Generated/ErrorOnSurface accessors are
//!    not exposed by that carrier yet — the two local carriers below keep
//!    the OCCT call form (reported gap).
//! 3. OCCT enum GeomFill_Trihedron (TKGeomAlgo/GeomFill_Trihedron.hxx
//!    L20-31) is the local enum below (no rcad GeomFill module yet).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;

use crate::brep_algo::tool::sub_shapes;
use crate::feat::loc_ope_pipe::BRepFillPipe;

/// OCCT enum GeomFill_Trihedron (GeomFill_Trihedron.hxx L20-31).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillTrihedron {
    IsCorrectedFrenet,
    IsFixed,
    IsFrenet,
    IsConstantNormal,
    IsDarboux,
    IsGuideAC,
    IsGuidePlan,
    IsGuideACWithContact,
    IsGuidePlanWithContact,
    IsDiscreteTrihedron,
}

// ---------------------------------------------------------------------------
// Engine accessors missing from the landed feat/loc_ope_pipe.rs carrier
// (architecture difference #2 — reported gap; the call sites keep the OCCT
// form).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Pipe::BRepFill_Pipe(Spine, Profile, aMode, ForceApproxC1)
/// (BRepFill_Pipe.hxx L60) — GAP.
fn brep_fill_pipe_new_with_mode(
    _the_spine: &Shape,
    _the_profile: &Shape,
    _the_mode: GeomFillTrihedron,
    _the_force_approx_c1: bool,
) -> BRepFillPipe {
    panic!("GAP: BRepFill_Pipe (aMode, ForceApproxC1) constructor (feat/loc_ope_pipe.rs carrier misses it)");
}

/// OCCT BRepFill_Pipe::Generated(S, List) (BRepFill_Pipe.cxx) — GAP.
fn brep_fill_pipe_generated(_the_pipe: &mut BRepFillPipe, _the_s: &Shape) -> Vec<Shape> {
    panic!("GAP: BRepFill_Pipe::Generated (feat/loc_ope_pipe.rs carrier misses it)");
}

/// OCCT BRepFill_Pipe::ErrorOnSurface() (BRepFill_Pipe.hxx) — GAP.
fn brep_fill_pipe_error_on_surface(_the_pipe: &BRepFillPipe) -> f64 {
    panic!("GAP: BRepFill_Pipe::ErrorOnSurface (feat/loc_ope_pipe.rs carrier misses it)");
}

/// OCCT BRepOffsetAPI_MakePipe (hxx L54-100).
pub struct BRepOffsetAPIMakePipe {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private member (hxx L98).
    my_pipe: BRepFillPipe, // OCCT: myPipe
}

impl BRepOffsetAPIMakePipe {
    /// OCCT BRepOffsetAPI_MakePipe::BRepOffsetAPI_MakePipe(Spine, Profile)
    /// (cxx L26-29).
    pub fn new(the_spine: &Shape, the_profile: &Shape) -> Self {
        // OCCT L28: myPipe(Spine, Profile).
        let my_pipe = BRepFillPipe::new(the_spine, the_profile);
        let mut r = BRepOffsetAPIMakePipe {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_pipe,
        };
        // OCCT L29: Build().
        r.build();
        r
    }

    /// OCCT BRepOffsetAPI_MakePipe::BRepOffsetAPI_MakePipe(Spine, Profile,
    /// aMode, ForceApproxC1) (cxx L44-48).
    pub fn new_with_mode(
        the_spine: &Shape,
        the_profile: &Shape,
        a_mode: GeomFillTrihedron,
        force_approx_c1: bool,
    ) -> Self {
        // OCCT L47: myPipe(Spine, Profile, aMode, ForceApproxC1).
        let my_pipe = brep_fill_pipe_new_with_mode(the_spine, the_profile, a_mode, force_approx_c1);
        let mut r = BRepOffsetAPIMakePipe {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_pipe,
        };
        // OCCT L48: Build().
        r.build();
        r
    }

    /// OCCT BRepOffsetAPI_MakePipe::Pipe() (cxx L50-53).
    pub fn pipe(&self) -> &BRepFillPipe {
        &self.my_pipe
    }

    /// OCCT BRepOffsetAPI_MakePipe::Build(...) (cxx L55-82).
    pub fn build(&mut self) {
        // OCCT L57: myShape = myPipe.Shape().
        self.my_shape = self.my_pipe.shape().clone().unwrap_or_else(Shape::null);
        // OCCT L58-59: Check for emptiness of result —
        // NCollection_IndexedMap theMap; TopExp::MapShapes(myShape, theMap).
        let the_map = sub_shapes(&self.my_shape);
        let the_map_extent = if self.my_shape.is_null() { 0 } else { 1 + the_map.len() };
        // OCCT L60-81.
        if the_map_extent == 1 {
            // OCCT L62: NotDone().
            self.my_done = false;
        } else {
            // OCCT L76: Done().
            self.my_done = true;
        }
    }

    /// OCCT BRepOffsetAPI_MakePipe::FirstShape() (cxx L84-87).
    pub fn first_shape(&self) -> Shape {
        self.my_pipe.first_shape().unwrap_or_else(Shape::null)
    }

    /// OCCT BRepOffsetAPI_MakePipe::LastShape() (cxx L89-92).
    pub fn last_shape(&self) -> Shape {
        self.my_pipe.last_shape().unwrap_or_else(Shape::null)
    }

    /// OCCT BRepOffsetAPI_MakePipe::Generated(S) (cxx L94-98).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L96: myPipe.Generated(S, myGenerated).
        self.my_generated = brep_fill_pipe_generated(&mut self.my_pipe, s);
        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_MakePipe::Generated(SSpine, SProfile)
    /// (cxx L100-131) — returns the generated elementary subshape.
    pub fn generated_elementary(&mut self, s_spine: &Shape, s_profile: &Shape) -> Shape {
        // OCCT L104-112.
        if s_profile.shape_type() == ShapeType::Edge {
            return self.my_pipe.face(s_spine, s_profile);
        } else if s_profile.shape_type() == ShapeType::Vertex {
            // OCCT L107-109: myPipe.Edge(TopoDS::Edge(SSpine),
            // TopoDS::Vertex(SProfile)).
            return self.my_pipe.edge(s_spine, s_profile);
        }

        // OCCT L114-116: POP pour NT — TopoDS_Shape bid; return bid.
        Shape::null()
    }

    /// OCCT BRepOffsetAPI_MakePipe::ErrorOnSurface() (cxx L133-136).
    pub fn error_on_surface(&self) -> f64 {
        brep_fill_pipe_error_on_surface(&self.my_pipe)
    }
}
