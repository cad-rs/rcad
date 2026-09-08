//! OCCT BRepOffsetAPI_NormalProjection — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_NormalProjection.cxx (L20-99) +
//!         BRepOffsetAPI_NormalProjection.hxx (L134-232).
//!
//! OCCT inheritance chain (hxx L134): BRepOffsetAPI_NormalProjection ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated) are kept as plain fields of the struct
//! (the BRepOffsetAPI facade precedent of this Stage 2e batch).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepAlgo_NormalProjection myNormalProjector is the
//!    landed brep_algo::normal_projection::BRepAlgoNormalProjection; the
//!    facade forwards every call 1:1 (cxx L20-99).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::GeomAbsShape;

use crate::brep_algo::normal_projection::BRepAlgoNormalProjection;

/// OCCT BRepOffsetAPI_NormalProjection (hxx L134-232).
pub struct BRepOffsetAPINormalProjection {
    // OCCT BRepBuilderAPI_Command base member (myDone) — the Done()/NotDone()
    // status flag.
    my_done: bool,
    // OCCT BRepBuilderAPI_MakeShape base members (BRepBuilderAPI_MakeShape.hxx).
    #[allow(dead_code)]
    my_shape: Shape,          // OCCT: myShape
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private member (hxx L231).
    my_normal_projector: BRepAlgoNormalProjection, // OCCT: myNormalProjector
}

impl Default for BRepOffsetAPINormalProjection {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetAPINormalProjection {
    /// OCCT BRepOffsetAPI_NormalProjection::BRepOffsetAPI_NormalProjection()
    /// (cxx L20).
    pub fn new() -> Self {
        BRepOffsetAPINormalProjection {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_normal_projector: BRepAlgoNormalProjection::new(),
        }
    }

    /// OCCT BRepOffsetAPI_NormalProjection::BRepOffsetAPI_NormalProjection(S)
    /// (cxx L22-25).
    pub fn new_with_shape(s: &Shape) -> Self {
        let mut r = Self::new();
        r.my_normal_projector.init(s);
        r
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Init(S) (cxx L27-30).
    pub fn init(&mut self, s: &Shape) {
        self.my_normal_projector.init(s);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Add(ToProj) (cxx L32-35).
    pub fn add(&mut self, to_proj: &Shape) {
        self.my_normal_projector.add(to_proj);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::SetParams(...) (cxx L37-44).
    pub fn set_params(
        &mut self,
        tol3_d: f64,
        tol2_d: f64,
        internal_continuity: GeomAbsShape,
        max_degree: usize,
        max_seg: usize,
    ) {
        self.my_normal_projector
            .set_params(tol3_d, tol2_d, internal_continuity, max_degree, max_seg);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::SetMaxDistance(MaxDist)
    /// (cxx L46-49).
    pub fn set_max_distance(&mut self, max_dist: f64) {
        self.my_normal_projector.set_max_distance(max_dist);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::SetLimit(FaceBounds)
    /// (cxx L51-54).
    pub fn set_limit(&mut self, face_bounds: bool) {
        self.my_normal_projector.set_limit(face_bounds);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Compute3d(With3d) (cxx L56-59).
    pub fn compute3d(&mut self, with3d: bool) {
        self.my_normal_projector.compute3d(with3d);
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Build(...) (cxx L61-66).
    pub fn build(&mut self) {
        self.my_normal_projector.build();
        self.my_shape = self.my_normal_projector.projection().clone();
        // OCCT L65: Done().
        self.my_done = true;
    }

    /// OCCT BRepOffsetAPI_NormalProjection::IsDone() (cxx L68-71).
    pub fn is_done(&self) -> bool {
        self.my_normal_projector.is_done()
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Couple(E) (cxx L73-76).
    pub fn couple(&self, e: &Shape) -> &Shape {
        self.my_normal_projector.couple(e)
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Generated(S) (cxx L78-82).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        self.my_normal_projector.generated(s)
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Projection() (cxx L84-87).
    pub fn projection(&self) -> &Shape {
        self.my_normal_projector.projection()
    }

    /// OCCT BRepOffsetAPI_NormalProjection::Ancestor(E) (cxx L89-92).
    pub fn ancestor(&self, e: &Shape) -> &Shape {
        self.my_normal_projector.ancestor(e)
    }

    /// OCCT BRepOffsetAPI_NormalProjection::BuildWire(ListOfWire)
    /// (cxx L96-99).
    pub fn build_wire(&self, list_of_wire: &mut Vec<Shape>) -> bool {
        self.my_normal_projector.build_wire(list_of_wire)
    }
}
