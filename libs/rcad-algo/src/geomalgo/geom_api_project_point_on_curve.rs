//! OCCT GeomAPI_ProjectPointOnCurve (TKGeomAlgo/GeomAPI) — the 1:1
//! translation of GeomAPI_ProjectPointOnCurve.hxx / .cxx / .lxx.
//!
//! OCCT members (hxx L129-133):
//!   L130: bool myIsDone;
//!   L131: int myIndex;
//!   L132: Extrema_ExtPC myExtPC;
//!   L133: GeomAdaptor_Curve myC;
//!
//! The engine is the kernel translation of `Extrema_ExtPC`
//! ([`ExtremaExtPC`], the `Extrema_GGExtPC` body).  Architecture difference:
//! `myExtPC` cannot be a struct member in rcad because the kernel engine
//! borrows the `Extrema_CurveTool` handle (`CurveToolHandle`) for its whole
//! lifetime, while the OCCT member copies the adaptor.  `my_ext_pc` therefore
//! holds exactly the engine queries the OCCT methods below read (`IsDone`,
//! `NbExt`, `SquareDistance(i)`, `Point(i)`), produced by the same
//! `Initialize` + `Perform` statements the OCCT Init/Perform run.  The
//! engine construction is fused into `ext_pc_perform` (cxx L58-59 / L95-96 /
//! L137) for that reason; `my_c` keeps the `GeomAdaptor_Curve` window so the
//! fused construction is the OCCT one.
//!
//! Failure paths are the OCCT ones: `NbPoints()` answers 0 when not done and
//! `NearestPoint` / `Parameter` / `Distance` /
//! `LowerDistanceParameter` / `LowerDistance` raise `StdFail_NotDone`.

use glam::DVec3;
use rcad_kernel::base::extrema::POnCurve;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::Curve3;

/// The `Extrema_ExtPC myExtPC` query surface (hxx L132) as read by
/// GeomAPI_ProjectPointOnCurve — see the module header for why it is a
/// snapshot instead of a member.
pub struct GeomAPIExtPCSnapshot {
    /// OCCT myExtPC.IsDone().
    pub is_done: bool,
    /// OCCT myExtPC.NbExt().
    pub nb_ext: usize,
    /// OCCT myExtPC.SquareDistance(i), in the OCCT 1-based index order.
    pub square_distance: Vec<f64>,
    /// OCCT myExtPC.Point(i), in the OCCT 1-based index order.
    pub point: Vec<POnCurve>,
}

impl GeomAPIExtPCSnapshot {
    fn empty() -> Self {
        GeomAPIExtPCSnapshot {
            is_done: false,
            nb_ext: 0,
            square_distance: Vec::new(),
            point: Vec::new(),
        }
    }
}

/// OCCT GeomAPI_ProjectPointOnCurve (hxx L32-134).
pub struct GeomAPIProjectPointOnCurve {
    /// hxx L130: bool myIsDone.
    my_is_done: bool,
    /// hxx L131: int myIndex (1-based, the OCCT index order).
    my_index: usize,
    /// hxx L132: Extrema_ExtPC myExtPC (see the module header).
    my_ext_pc: GeomAPIExtPCSnapshot,
    /// hxx L133: GeomAdaptor_Curve myC.
    my_c: Option<GeomCurveAdaptor>,
}

impl GeomAPIProjectPointOnCurve {
    /// OCCT GeomAPI_ProjectPointOnCurve() (cxx L25-29): myIsDone = false,
    /// myIndex = 0.
    pub fn new() -> Self {
        GeomAPIProjectPointOnCurve {
            my_is_done: false,
            my_index: 0,
            my_ext_pc: GeomAPIExtPCSnapshot::empty(),
            my_c: None,
        }
    }

    /// OCCT GeomAPI_ProjectPointOnCurve(P, Curve) (cxx L33-37) — the ctor
    /// body is `Init(P, Curve)`.
    pub fn new_point_curve(the_p: DVec3, the_curve: &Curve3) -> Self {
        let mut a_this = GeomAPIProjectPointOnCurve::new();
        a_this.init_point_curve(the_p, the_curve);
        a_this
    }

    /// OCCT GeomAPI_ProjectPointOnCurve(P, Curve, Umin, Usup) (cxx L41-47) —
    /// the ctor body is `Init(P, Curve, Umin, Usup)`.
    pub fn new_point_curve_ranged(
        the_p: DVec3,
        the_curve: &Curve3,
        the_umin: f64,
        the_usup: f64,
    ) -> Self {
        let mut a_this = GeomAPIProjectPointOnCurve::new();
        a_this.init_point_curve_ranged(the_p, the_curve, the_umin, the_usup);
        a_this
    }

    /// OCCT Init(P, Curve) (cxx L51-81).
    pub fn init_point_curve(&mut self, the_p: DVec3, the_curve: &Curve3) {
        // cxx L53: myC.Load(Curve);
        self.my_c = Some(GeomCurveAdaptor::new(the_curve.clone()));
        // cxx L54-57: (the commented-out `Extrema_ExtPC theExtPC(P, myC);
        // myExtPC = theExtPC;` block).
        // cxx L58-59: myExtPC.Initialize(myC, myC.FirstParameter(),
        // myC.LastParameter()); myExtPC.Perform(P).
        self.ext_pc_perform(the_p);
        // cxx L61-80.
        self.update_lower_distance();
    }

    /// OCCT Init(P, Curve, Umin, Usup) (cxx L85-118).
    pub fn init_point_curve_ranged(
        &mut self,
        the_p: DVec3,
        the_curve: &Curve3,
        the_umin: f64,
        the_usup: f64,
    ) {
        // cxx L90: myC.Load(Curve, Umin, Usup);
        self.my_c = Some(GeomCurveAdaptor::with_range(
            the_curve.clone(),
            the_umin,
            the_usup,
        ));
        // cxx L91-94: (the commented-out `Extrema_ExtPC theExtPC(P, myC);`
        // block).
        // cxx L95-96: myExtPC.Initialize(myC, myC.FirstParameter(),
        // myC.LastParameter()); myExtPC.Perform(P) — the loaded window makes
        // FirstParameter/LastParameter the (Umin, Usup) pair.
        self.ext_pc_perform(the_p);
        // cxx L98-117.
        self.update_lower_distance();
    }

    /// OCCT Init(Curve, Umin, Usup) (cxx L123-131).
    pub fn init_curve(&mut self, the_curve: &Curve3, the_umin: f64, the_usup: f64) {
        // cxx L127: myC.Load(Curve, Umin, Usup);
        self.my_c = Some(GeomCurveAdaptor::with_range(
            the_curve.clone(),
            the_umin,
            the_usup,
        ));
        // cxx L128: (the commented-out `myExtPC = Extrema_ExtPC(P, myC);`).
        // cxx L129: myExtPC.Initialize(myC, Umin, Usup); — the engine state is
        // rebuilt by the later Perform (cxx L137) with the same window, which
        // is what the fused construction in `ext_pc_perform` does.
        self.my_ext_pc = GeomAPIExtPCSnapshot::empty();
        // cxx L130: myIsDone = false;
        self.my_is_done = false;
    }

    /// OCCT Perform(aP3D) (cxx L135-156).
    pub fn perform(&mut self, the_p: DVec3) {
        // cxx L137: myExtPC.Perform(aP3D);
        self.ext_pc_perform(the_p);
        // cxx L139-155.
        self.update_lower_distance();
    }

    /// OCCT cxx L58-59 (the L95-96 and L137 twins) — the engine construction
    /// and query snapshot.  See the module header for the member-mapping
    /// reason (the kernel engine borrows `myC` for its lifetime).
    fn ext_pc_perform(&mut self, the_p: DVec3) {
        let a_c = self
            .my_c
            .as_ref()
            .expect("GeomAPI_ProjectPointOnCurve: myC is not loaded");
        // cxx L58: myExtPC.Initialize(myC, myC.FirstParameter(),
        // myC.LastParameter()) — the OCCT default theTolF = 1.0e-10.
        let a_tool = CurveToolHandle::for_curve3(&a_c.curve, a_c, a_c);
        let mut a_ext_pc = ExtremaExtPC::new();
        a_ext_pc.initialize(&a_tool, a_c.first, a_c.last, 1.0e-10);
        // cxx L59: myExtPC.Perform(P).
        a_ext_pc.perform(the_p);
        // The engine queries the methods below read (cxx L61 / L68-78).
        let a_is_done = a_ext_pc.is_done();
        let a_nb_ext = if a_is_done { a_ext_pc.nb_ext() } else { 0 };
        self.my_ext_pc = GeomAPIExtPCSnapshot {
            is_done: a_is_done,
            nb_ext: a_nb_ext,
            square_distance: (1..=a_nb_ext)
                .map(|i| a_ext_pc.square_distance(i))
                .collect(),
            point: (1..=a_nb_ext).map(|i| a_ext_pc.point(i).clone()).collect(),
        };
    }

    /// OCCT cxx L61-80 — the Init tail shared by the three Init forms (the
    /// cxx L98-117 and L139-155 bodies are the same statements).
    fn update_lower_distance(&mut self) {
        self.my_is_done = self.my_ext_pc.is_done && (self.my_ext_pc.nb_ext > 0);
        if self.my_is_done {
            // evaluate the lower distance and its index;
            let mut dist2;
            let mut dist2_min = self.my_ext_pc.square_distance[0];
            self.my_index = 1;
            for i in 2..=self.my_ext_pc.nb_ext {
                dist2 = self.my_ext_pc.square_distance[i - 1];
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    self.my_index = i;
                }
            }
        }
    }

    /// OCCT NbPoints() (cxx L161-171).
    pub fn nb_points(&self) -> usize {
        if self.my_is_done {
            self.my_ext_pc.nb_ext
        } else {
            0
        }
    }

    /// OCCT Point(Index) (cxx L175-180).
    pub fn point(&self, index: usize) -> DVec3 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: GeomAPI_ProjectPointOnCurve::Point");
        }
        self.my_ext_pc.point[index - 1].point
    }

    /// OCCT Parameter(Index) (cxx L184-189).
    pub fn parameter(&self, index: usize) -> f64 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: GeomAPI_ProjectPointOnCurve::Parameter");
        }
        self.my_ext_pc.point[index - 1].param
    }

    /// OCCT Parameter(Index, U&) (cxx L193-198).
    pub fn parameter_into(&self, index: usize, the_u: &mut f64) {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: GeomAPI_ProjectPointOnCurve::Parameter");
        }
        *the_u = self.my_ext_pc.point[index - 1].param;
    }

    /// OCCT Distance(Index) (cxx L202-207).
    pub fn distance(&self, index: usize) -> f64 {
        if index < 1 || index > self.nb_points() {
            panic!("Standard_OutOfRange: GeomAPI_ProjectPointOnCurve::Distance");
        }
        self.my_ext_pc.square_distance[index - 1].sqrt()
    }

    /// OCCT NearestPoint() (cxx L211-216).
    pub fn nearest_point(&self) -> DVec3 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: GeomAPI_ProjectPointOnCurve::NearestPoint");
        }
        self.my_ext_pc.point[self.my_index - 1].point
    }

    /// OCCT LowerDistanceParameter() (cxx L234-239).
    pub fn lower_distance_parameter(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: GeomAPI_ProjectPointOnCurve::LowerDistanceParameter");
        }
        self.my_ext_pc.point[self.my_index - 1].param
    }

    /// OCCT LowerDistance() (cxx L243-248).
    pub fn lower_distance(&self) -> f64 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: GeomAPI_ProjectPointOnCurve::LowerDistance");
        }
        self.my_ext_pc.square_distance[self.my_index - 1].sqrt()
    }

    /// OCCT Extrema() (lxx L19-22) — `return myExtPC;`.  The rcad engine is
    /// borrow-bound, so the accessor hands out the query snapshot (module
    /// header).
    pub fn extrema(&self) -> &GeomAPIExtPCSnapshot {
        &self.my_ext_pc
    }
}

impl Default for GeomAPIProjectPointOnCurve {
    fn default() -> Self {
        Self::new()
    }
}

// The `GeomAPI_ProjectPointOnCurve::operator int()` (cxx L220-223), `operator
// gp_Pnt()` (cxx L227-230) and `operator double()` (cxx L252-255) are C++
// implicit-conversion operators; Rust has no equivalent, so the callers state
// the queried method explicitly (`nb_points()` / `nearest_point()` /
// `lower_distance()`), which is exactly what the operators forward to.
