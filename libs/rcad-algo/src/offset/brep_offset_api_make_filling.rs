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
//! 2. The engine member BRepFill_Filling myFilling (TKBool/BRepFill) has no
//!    rcad translation yet — the BRepFillFilling carrier below keeps the
//!    OCCT constructor/method surface with GAP panics (port plan section
//!    0.6); every facade body around the engine calls is translated 1:1.

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::GeomAbsShape;

use glam::DVec3;

// ---------------------------------------------------------------------------
// GAP carrier (architecture difference #2).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Filling (TKBool/BRepFill, BRepFill_Filling.hxx) — the
/// N-side filling engine of MakeFilling (architecture difference #2; GAP: no
/// rcad translation yet — the GAP panics are the section 0.6 annotation; the
/// constructor and field storage keep the OCCT form).
pub struct BRepFillFilling {
    my_is_done: bool, // OCCT: myIsDone
}

impl BRepFillFilling {
    /// OCCT BRepFill_Filling::BRepFill_Filling(Degree, NbPtsOnCur, NbIter,
    /// Anisotropie, Tol2d, Tol3d, TolAng, TolCurv, MaxDeg, MaxSegments)
    /// (BRepFill_Filling.cxx L43-63) — GAP: the parameter storage keeps the
    /// OCCT form, the engine computation is not translated.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        _degree: i32,
        _nb_pts_on_cur: i32,
        _nb_iter: i32,
        _anisotropie: bool,
        _tol2d: f64,
        _tol3d: f64,
        _tolang: f64,
        _tolcurv: f64,
        _max_deg: i32,
        _max_segments: i32,
    ) -> Self {
        BRepFillFilling { my_is_done: false }
    }

    /// OCCT BRepFill_Filling::SetConstrParam(Tol2d, Tol3d, TolAng, TolCurv)
    /// — GAP.
    pub fn set_constr_param(&mut self, _tol2d: f64, _tol3d: f64, _tolang: f64, _tolcurv: f64) {
        panic!("GAP: BRepFill_Filling::SetConstrParam (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::SetResolParam(Degree, NbPtsOnCur, NbIter,
    /// Anisotropie) — GAP.
    pub fn set_resol_param(
        &mut self,
        _degree: i32,
        _nb_pts_on_cur: i32,
        _nb_iter: i32,
        _anisotropie: bool,
    ) {
        panic!("GAP: BRepFill_Filling::SetResolParam (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::SetApproxParam(MaxDeg, MaxSegments) — GAP.
    pub fn set_approx_param(&mut self, _max_deg: i32, _max_segments: i32) {
        panic!("GAP: BRepFill_Filling::SetApproxParam (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::LoadInitSurface(Surf) — GAP.
    pub fn load_init_surface(&mut self, _surf: &Shape) {
        panic!("GAP: BRepFill_Filling::LoadInitSurface (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Add(Constr, Order, IsBound) — GAP.
    pub fn add(&mut self, _constr: &Shape, _order: GeomAbsShape, _is_bound: bool) -> i32 {
        panic!("GAP: BRepFill_Filling::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Add(Constr, Support, Order, IsBound) — GAP.
    pub fn add_with_support(
        &mut self,
        _constr: &Shape,
        _support: &Shape,
        _order: GeomAbsShape,
        _is_bound: bool,
    ) -> i32 {
        panic!("GAP: BRepFill_Filling::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Add(Support, Order) — GAP.
    pub fn add_free_constraint(&mut self, _support: &Shape, _order: GeomAbsShape) -> i32 {
        panic!("GAP: BRepFill_Filling::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Add(Point) — GAP.
    pub fn add_point(&mut self, _point: &DVec3) -> i32 {
        panic!("GAP: BRepFill_Filling::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Add(U, V, Support, Order) — GAP.
    pub fn add_point_on_support(
        &mut self,
        _u: f64,
        _v: f64,
        _support: &Shape,
        _order: GeomAbsShape,
    ) -> i32 {
        panic!("GAP: BRepFill_Filling::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Build() — GAP.
    pub fn build(&mut self) {
        panic!("GAP: BRepFill_Filling::Build (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_Filling::Generated(S) — GAP.
    pub fn generated(&mut self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepFill_Filling::Generated (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Face().
    pub fn face(&self) -> Shape {
        Shape::null()
    }

    /// OCCT BRepFill_Filling::G0Error() — GAP.
    pub fn g0_error(&self) -> f64 {
        panic!("GAP: BRepFill_Filling::G0Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::G1Error() — GAP.
    pub fn g1_error(&self) -> f64 {
        panic!("GAP: BRepFill_Filling::G1Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::G2Error() — GAP.
    pub fn g2_error(&self) -> f64 {
        panic!("GAP: BRepFill_Filling::G2Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::G0Error(Index) — GAP.
    pub fn g0_error_at(&mut self, _index: i32) -> f64 {
        panic!("GAP: BRepFill_Filling::G0Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::G1Error(Index) — GAP.
    pub fn g1_error_at(&mut self, _index: i32) -> f64 {
        panic!("GAP: BRepFill_Filling::G1Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::G2Error(Index) — GAP.
    pub fn g2_error_at(&mut self, _index: i32) -> f64 {
        panic!("GAP: BRepFill_Filling::G2Error (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_Filling::Surface() — the deformation carrier of
    /// MakeFilling (BRepFill_Filling.hxx Surface()); kept null under the
    /// same GAP.
    pub fn surface(&self) -> Option<Surface3> {
        None
    }
}

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
            my_filling: BRepFillFilling::new(
                degree, nb_pts_on_cur, nb_iter, anisotropie, tol2d, tol3d, tolang, tolcurv,
                max_deg, max_segments,
            ),
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
        self.my_filling.add_point(point)
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
        self.my_filling.add_point_on_support(u, v, support, order)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::Build(...) (cxx L117-122).
    pub fn build(&mut self) {
        // OCCT L119: myFilling.Build().
        self.my_filling.build();
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
        self.my_filling.generated(s)
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
        self.my_filling.g0_error_at(index)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G1Error(Index) (cxx L163-168) —
    /// returns maximum angle between the constraint number Index and the
    /// resulting surface.
    pub fn g1_error_at(&mut self, index: i32) -> f64 {
        self.my_filling.g1_error_at(index)
    }

    /// OCCT BRepOffsetAPI_MakeFilling::G2Error(Index) (cxx L170-175) —
    /// returns maximum difference of curvature between the constraint number
    /// Index and the resulting surface.
    pub fn g2_error_at(&mut self, index: i32) -> f64 {
        self.my_filling.g2_error_at(index)
    }
}
