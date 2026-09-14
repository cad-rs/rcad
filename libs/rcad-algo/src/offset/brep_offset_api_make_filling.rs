//! OCCT BRepOffsetAPI_MakeFilling — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakeFilling.cxx (L22-212) +
//!         BRepOffsetAPI_MakeFilling.hxx (L75-343).
//!
//! OCCT inheritance chain (hxx L75): BRepOffsetAPI_MakeFilling ->
//! BRepBuilderAPI_MakeShape.  Rust has no inheritance: the base-class
//! members (myShape, myGenerated, the Done flag) are kept as plain fields of
//! the struct (the Stage 2e facade precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member BRepFill_Filling myFilling (TKBool/BRepFill) is the
//!    brep_fill/brep_fill_filling.rs translation (imported below following
//!    the OCCT hxx member form).
//! 3. The engine consumes the rcad BRep pool (the BRepAdaptor_Surface read
//!    in Add(U, V, Support, Order) and the BRepLib_MakeFace writes in
//!    Build): the class holds my_brep (the brep_offset_api_draft_angle.rs
//!    facade precedent, architecture difference #4).

use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::GeomAbsShape;

use glam::DVec3;

use crate::brep_fill::brep_fill_filling::BRepFillFilling;

/// OCCT BRepOffsetAPI_MakeFilling (hxx L75-343).
pub struct BRepOffsetAPIMakeFilling {
    // OCCT BRepBuilderAPI base members.
    #[allow(dead_code)]
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private member (hxx L341).
    my_filling: BRepFillFilling, // OCCT: myFilling
    // The rcad pool stand-in consumed by the engine (architecture
    // difference #4; the brep_offset_api_draft_angle.rs precedent).
    my_brep: BRep,
}

impl BRepOffsetAPIMakeFilling {
    /// OCCT BRepOffsetAPI_MakeFilling::BRepOffsetAPI_MakeFilling(Degree,
    /// NbPtsOnCur, NbIter, Anisotropie, Tol2d, Tol3d, TolAng, TolCurv,
    /// MaxDeg, MaxSegments) (cxx L24-46; hxx defaults Degree=3,
    /// NbPtsOnCur=15, NbIter=2, Anisotropie=false, Tol2d=0.00001,
    /// Tol3d=0.0001, TolAng=0.01, TolCurv=0.1, MaxDeg=8, MaxSegments=9).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
        tolcurv: f64,
        max_deg: i32,
        max_segments: i32,
    ) -> Self {
        // OCCT L40-44: myFilling(Degree, NbPtsOnCur, NbIter, Anisotropie,
        // Tol2d, Tol3d, TolAng, TolCurv, MaxDeg, MaxSegments).
        BRepOffsetAPIMakeFilling {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_filling: BRepFillFilling::new_with_args(
                degree, nb_pts_on_cur, nb_iter, anisotropie, tol2d, tol3d, tolang, tolcurv,
                max_deg, max_segments,
            ),
            my_brep: BRep::new(),
        }
    }

    /// OCCT BRepOffsetAPI_MakeFilling::SetConstrParam(Tol2d, Tol3d, TolAng,
    /// TolCurv) (cxx L48-55).
    pub fn set_constr_param(&mut self, tol2d: f64, tol3d: f64, tolang: f64, tolcurv: f64) {
        self.my_filling.set_constr_param(tol2d, tol3d, tolang, tolcurv);
    }

    /// OCCT BRepOffsetAPI_MakeFilling::SetResolParam(Degree, NbPtsOnCur,
    /// NbIter, Anisotropie) (cxx L57-64).
    pub fn set_resol_param(
        &mut self,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        anisotropie: bool,
    ) {
        self.my_filling
            .set_resol_param(degree, nb_pts_on_cur, nb_iter, anisotropie);
    }

    /// OCCT BRepOffsetAPI_MakeFilling::SetApproxParam(MaxDeg, MaxSegments)
    /// (cxx L66-71).
    pub fn set_approx_param(&mut self, max_deg: i32, max_segments: i32) {
        self.my_filling.set_approx_param(max_deg, max_segments);
    }

    /// OCCT BRepOffsetAPI_MakeFilling::LoadInitSurface(Surf) (cxx L73-76).
    pub fn load_init_surface(&mut self, surf: &Shape) {
        self.my_filling.load_init_surface(surf);
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Add(Constr, Order, IsBound)
    /// (cxx L78-84).
    pub fn add(&mut self, constr: &Shape, order: GeomAbsShape, is_bound: bool) -> i32 {
        self.my_filling.add(constr, order, is_bound)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Add(Constr, Support, Order, IsBound)
    /// (cxx L86-92) — adds an edge with supporting face as a constraint.
    pub fn add_with_support(
        &mut self,
        constr: &Shape,
        support: &Shape,
        order: GeomAbsShape,
        is_bound: bool,
    ) -> i32 {
        self.my_filling
            .add_with_support(constr, support, order, is_bound)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Add(Support, Order) (cxx L94-100) —
    /// adds a "free constraint": face without edge.
    pub fn add_free_constraint(&mut self, support: &Shape, order: GeomAbsShape) -> i32 {
        self.my_filling.add_free_constraint(support, order)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Add(Point) (cxx L102-107).
    pub fn add_point(&mut self, point: &DVec3) -> i32 {
        self.my_filling.add_point(*point)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Add(U, V, Support, Order)
    /// (cxx L109-115) — adds a point constraint on a face.
    pub fn add_point_on_support(
        &mut self,
        u: f64,
        v: f64,
        support: &Shape,
        order: GeomAbsShape,
    ) -> i32 {
        // The engine reads the support surface through the BRep pool
        // (BRepAdaptor_Surface; architecture difference #4).
        self.my_filling
            .add_point_on_support(&self.my_brep, u, v, support, order)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Build(...) (cxx L117-122).
    pub fn build(&mut self) {
        // OCCT L119: myFilling.Build().
        self.my_filling.build(&mut self.my_brep);
        // OCCT L120: myShape = myFilling.Face().
        self.my_shape = self.my_filling.face();
    }

    /// OCCT BRepOffsetAPI_MakeFilling::IsDone() (cxx L124-127).
    pub fn is_done(&self) -> bool {
        self.my_filling.is_done()
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Generated(S) (cxx L129-133) —
    /// returns the new edge (first in list) made from old edge "S".
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        self.my_filling.generated(s).clone()
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G0Error() (cxx L135-140) — returns
    /// maximum distance from boundary to the resulting surface.
    pub fn g0_error(&self) -> f64 {
        self.my_filling.g0_error()
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G1Error() (cxx L142-147) — returns
    /// maximum angle between the resulting surface and constraint surfaces
    /// at boundaries.
    pub fn g1_error(&self) -> f64 {
        self.my_filling.g1_error()
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G2Error() (cxx L149-154) — returns
    /// maximum difference of curvature between the resulting surface and
    /// constraint surfaces at boundaries.
    pub fn g2_error(&self) -> f64 {
        self.my_filling.g2_error()
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G0Error(Index) (cxx L156-161) —
    /// returns maximum distance between the constraint number Index and the
    /// resulting surface.
    pub fn g0_error_at(&mut self, index: i32) -> f64 {
        self.my_filling.g0_error_index(index)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G1Error(Index) (cxx L163-168) —
    /// returns maximum angle between the constraint number Index and the
    /// resulting surface.
    pub fn g1_error_at(&mut self, index: i32) -> f64 {
        self.my_filling.g1_error_index(index)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G2Error(Index) (cxx L170-175) —
    /// returns maximum difference of curvature between the constraint number
    /// Index and the resulting surface.
    pub fn g2_error_at(&mut self, index: i32) -> f64 {
        self.my_filling.g2_error_index(index)
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_MakeShape.hxx: Standard_EXPORT const TopoDS_Shape&
    /// Shape() const; raises StdFail_NotDone when not done).
    pub fn shape(&self) -> Shape {
        assert!(
            self.is_done(),
            "StdFail_NotDone: BRepOffsetAPI_MakeFilling::Shape()"
        );
        self.my_shape.clone()
    }

    /// Test-world extraction bridge (the FilletResult.brep pattern,
    /// algo_ext::topods_ext::extract_result_brep): flattens the
    /// (my_brep arena, my_shape root) pair into the self-contained BRep the
    /// test world consumes (StepWriter / total_surface_area).  OCCT has no
    /// equivalent (the TopoDS_Shape carries its arena implicitly) — the
    /// rcad BRep-pool architecture difference #4 glue.
    pub fn result_brep(&mut self) -> Option<rcad_kernel::topo::topods::BRep> {
        if !self.is_done() || self.my_shape.is_null() {
            return None;
        }
        let locations = self.my_brep.locations.clone();
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            locations,
        ))
    }
}
