//! OCCT Geom2dAPI_ProjectPointOnCurve (TKGeomAlgo/Geom2dAPI) — the 1:1
//! translation of Geom2dAPI_ProjectPointOnCurve.hxx / .cxx.
//!
//! OCCT members (hxx L129-133):
//!   L130: bool myIsDone;
//!   L131: int myIndex;
//!   L132: Extrema_ExtPC2d myExtPC;
//!   L133: Geom2dAdaptor_Curve myC;
//!
//! The engine is the kernel translation of `Extrema_ExtPC2d`
//! ([`ExtPC2d`], the `Extrema_GGExtPC` 2D instantiation).  The kernel engine
//! owns its results, so `my_ext_pc` is a direct member (the 3D twin
//! `geom_api_project_point_on_curve` needs a snapshot only because the 3D
//! engine borrows the curve-tool handle; the 2D one does not).
//!
//! Failure paths are the OCCT ones: `NbPoints()` answers 0 when not done and
//! `NearestPoint` / `Parameter` / `Distance` /
//! `LowerDistanceParameter` / `LowerDistance` raise `StdFail_NotDone`.

use glam::DVec2;
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::base::proj_lib::adaptor::Geom2dCurveAdaptor;
use rcad_kernel::geom::{Curve2d, Curve2dEval};

/// OCCT Geom2dAPI_ProjectPointOnCurve (hxx L32-134).
pub struct Geom2dAPIProjectPointOnCurve {
    /// hxx L130: bool myIsDone.
    my_is_done: bool,
    /// hxx L131: int myIndex (1-based, the OCCT index order; the default-ctor
    /// -1 cannot live in usize and every reader guards on myIsDone first).
    my_index: usize,
    /// hxx L132: Extrema_ExtPC2d myExtPC.
    my_ext_pc: Option<ExtPC2d>,
    /// hxx L133: Geom2dAdaptor_Curve myC.
    my_c: Option<Geom2dCurveAdaptor>,
}

impl Geom2dAPIProjectPointOnCurve {
    /// OCCT Geom2dAPI_ProjectPointOnCurve() (cxx L25-29).
    pub fn new() -> Self {
        Geom2dAPIProjectPointOnCurve {
            my_is_done: false,
            my_index: 0,
            my_ext_pc: None,
            my_c: None,
        }
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, Curve) (cxx L33-37) — the ctor
    /// body is `Init(P, Curve)`.
    pub fn new_point_curve(the_p: DVec2, the_curve: &Curve2d) -> Self {
        let mut a_this = Geom2dAPIProjectPointOnCurve::new();
        a_this.init_point_curve(the_p, the_curve);
        a_this
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, Curve, Umin, Usup) (cxx L41-47)
    /// — the ctor body is `Init(P, Curve, Umin, Usup)`.
    pub fn new_point_curve_ranged(
        the_p: DVec2,
        the_curve: &Curve2d,
        the_umin: f64,
        the_usup: f64,
    ) -> Self {
        let mut a_this = Geom2dAPIProjectPointOnCurve::new();
        a_this.init_point_curve_ranged(the_p, the_curve, the_umin, the_usup);
        a_this
    }

    /// OCCT Init(P, Curve) (cxx L51-54) — `Curve->FirstParameter()` /
    /// `Curve->LastParameter()` is the Geom2d_Curve natural domain (the
    /// rcad `default_domain()` per-type mapping).
    pub fn init_point_curve(&mut self, the_p: DVec2, the_curve: &Curve2d) {
        let [a_first, a_last] = the_curve.default_domain();
        self.init_point_curve_ranged(the_p, the_curve, a_first, a_last);
    }

    /// OCCT Init(P, Curve, Umin, Usup) (cxx L58-88).
    pub fn init_point_curve_ranged(
        &mut self,
        the_p: DVec2,
        the_curve: &Curve2d,
        the_umin: f64,
        the_usup: f64,
    ) {
        // cxx L63: myC.Load(Curve, Umin, Usup);
        self.my_c = Some(Geom2dCurveAdaptor::with_range(
            the_curve.clone(),
            the_umin,
            the_usup,
        ));
        let a_c = self.my_c.as_ref().unwrap();
        // cxx L65: Extrema_ExtPC2d theExtPC2d(P, myC); — the GGExtPC
        // (P, C, FirstParameter, LastParameter, theTolF = 1.0e-10) ctor
        // (Extrema_GGExtPC.hxx L99-106).
        let mut a_ext_pc_2d = ExtPC2d::new(
            the_p,
            &a_c.curve,
            1.0e-10,
            a_c.first,
            a_c.last,
        );
        // cxx L67: myExtPC = theExtPC2d;
        // cxx L69: myIsDone = myExtPC.IsDone() && (myExtPC.NbExt() > 0).
        self.my_is_done = a_ext_pc_2d.is_done() && a_ext_pc_2d.nb_ext() > 0;
        // cxx L71-87: evaluate the lower distance and its index.
        if self.my_is_done {
            let mut dist2;
            let mut dist2_min = a_ext_pc_2d.square_distance(1);
            self.my_index = 1;
            for i in 2..=a_ext_pc_2d.nb_ext() {
                dist2 = a_ext_pc_2d.square_distance(i);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    self.my_index = i;
                }
            }
        }
        self.my_ext_pc = Some(a_ext_pc_2d);
    }

    /// OCCT NbPoints() (cxx L92-102).
    pub fn nb_points(&self) -> usize {
        match (&self.my_is_done, &self.my_ext_pc) {
            (true, Some(a_ext_pc)) => a_ext_pc.nb_ext(),
            _ => 0,
        }
    }

    /// OCCT Point(Index) (cxx L106-111).
    pub fn point(&self, index: usize) -> DVec2 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: Geom2dAPI_ProjectPointOnCurve::Point");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .point(index)
            .point
    }

    /// OCCT Parameter(Index) (cxx L115-120).
    pub fn parameter(&self, index: usize) -> f64 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: Geom2dAPI_ProjectPointOnCurve::Parameter");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .point(index)
            .param
    }

    /// OCCT Parameter(Index, U&) (cxx L124-129).
    pub fn parameter_into(&self, index: usize, the_u: &mut f64) {
        *the_u = self.parameter(index);
    }

    /// OCCT Distance(Index) (cxx L133-138).
    pub fn distance(&self, index: usize) -> f64 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: Geom2dAPI_ProjectPointOnCurve::Distance");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .square_distance(index)
            .sqrt()
    }

    /// OCCT NearestPoint() (cxx L142-147).
    pub fn nearest_point(&self) -> DVec2 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: Geom2dAPI_ProjectPointOnCurve::NearestPoint");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .point(self.my_index)
            .point
    }

    /// OCCT LowerDistanceParameter() (cxx L165-170).
    pub fn lower_distance_parameter(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: Geom2dAPI_ProjectPointOnCurve::LowerDistanceParameter");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .point(self.my_index)
            .param
    }

    /// OCCT LowerDistance() (cxx L174-179).
    pub fn lower_distance(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: Geom2dAPI_ProjectPointOnCurve::LowerDistance");
        }
        self.my_ext_pc
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve: myExtPC is not loaded")
            .square_distance(self.my_index)
            .sqrt()
    }
}

impl Default for Geom2dAPIProjectPointOnCurve {
    fn default() -> Self {
        Self::new()
    }
}

// The `Geom2dAPI_ProjectPointOnCurve::operator int()` (cxx L151-154), `operator
// gp_Pnt2d()` (cxx L158-161) and `operator double()` (cxx L183-186) are C++
// implicit-conversion operators; Rust has no equivalent, so the callers state
// the queried method explicitly (`nb_points()` / `nearest_point()` /
// `lower_distance()`), which is exactly what the operators forward to.
