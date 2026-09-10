//! OCCT BRepLib_ValidateEdge (TKTopAlgo/BRepLib — `BRepLib_ValidateEdge.hxx`
//! L26-103 + `BRepLib_ValidateEdge.cxx` L24-244) — computes the max distance
//! between a 3D curve and a curve on surface. Two methods: approximate using
//! a finite number of points (default) and exact (through
//! `GeomLib_CheckCurveOnSurface`).
//!
//! Architecture bridges (the W1 shape_analysis/edge.rs numbering style):
//! 1. OCCT `occ::handle<Adaptor3d_Curve>` (the reference curve; the concrete
//!    instance used by every OCCT consumer is `GeomAdaptor_Curve(C, f, l)`)
//!    -> [`GeomAdaptorCurve`] value (clone-on-handle semantics).
//! 2. OCCT `occ::handle<Adaptor3d_CurveOnSurface>` -> [`Adaptor3dCurveOnSurface`]
//!    over the `Geom2dAdaptor_Curve` + `GeomAdaptor_Surface` re-hosts below.
//!    Only the members `BRepLib_ValidateEdge` consumes are carried
//!    (FirstParameter/LastParameter/Value/Resolution/GetSurface).
//! 3. `BRepCheck::PrecCurve/PrecSurface` (BRepCheck.cxx L70-127) -> the
//!    re-hosts below; they match the rcad `Curve3`/`Surface3` enums directly
//!    (the OCCT `Adaptor3d_Curve::GetType()` dispatch).
//! 4. `Extrema_LocateExtPC` (TKGeomBase) -> the kernel
//!    `extrema_locate_ext_pc` (generic over the `CurveEval` adaptor trait)
//!    wrapped in the [`ExtremaLocateExtPC`] Initialize/Perform contract below.
//! 5. `GeomLib_CheckCurveOnSurface` (TKGeomBase) is not translated yet; the
//!    existing GAP carrier `geomalgo::geom_lib_check_curve_on_surface` is
//!    called so `processExact` preserves the OCCT failure path (IsDone()
//!    false from the panic-free construction; Perform panics — plan §0.6).
//!    GAP: closes with the TKGeomBase GeomLib batch.

use rcad_kernel::base::extrema::POnCurve;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};

use crate::geomalgo::geom_lib_check_curve_on_surface::GeomLibCheckCurveOnSurface;

// ---------------------------------------------------------------------------
// Standard_Real re-hosts (pure math)
// ---------------------------------------------------------------------------

/// OCCT Standard_Real.hxx L176-179: RealEpsilon() = DBL_EPSILON.
const REAL_EPSILON: f64 = f64::EPSILON;

/// OCCT Standard_Real.hxx L240-247: Epsilon(theValue) - the absolute value of
/// the difference between theValue and the nearest representable value chosen
/// in the direction of infinity with the same sign; for 0 the minimal
/// positive representable value.
fn standard_real_epsilon(the_value: f64) -> f64 {
    // std::nextafter is stable only through the standard library's float
    // methods; the equivalent next-representable step is one ULP.
    if the_value >= 0.0 {
        let a_next = f64::from_bits(the_value.to_bits() + 1);
        a_next - the_value
    } else {
        let a_prev = f64::from_bits(the_value.to_bits() - 1);
        the_value - a_prev
    }
}

// ---------------------------------------------------------------------------
// GeomAdaptor_Curve re-host (architecture bridge #1)
// ---------------------------------------------------------------------------

/// OCCT GeomAdaptor_Curve (TKG3d/GeomAdaptor) — the adaptor over a
/// `Geom_Curve` restricted to [the_u_first, the_u_last]. Only the members the
/// `BRepLib_ValidateEdge` / `BRepCheck::PrecCurve` consumers use are carried.
#[derive(Debug, Clone)]
pub struct GeomAdaptorCurve {
    /// OCCT myCurve (the underlying Geom_Curve value).
    my_curve: Curve3,
    /// OCCT myFirst.
    my_first: f64,
    /// OCCT myLast.
    my_last: f64,
}

impl GeomAdaptorCurve {
    /// OCCT GeomAdaptor_Curve::GeomAdaptor_Curve(C, U1, U2).
    pub fn new(the_curve: Curve3, the_u_first: f64, the_u_last: f64) -> Self {
        GeomAdaptorCurve {
            my_curve: the_curve,
            my_first: the_u_first,
            my_last: the_u_last,
        }
    }

    /// OCCT FirstParameter() — the restricted range start.
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT LastParameter() — the restricted range end.
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT Value(U) — the curve point (no range clamp, the OCCT adaptor
    /// evaluates the underlying curve directly).
    pub fn value(&self, the_u: f64) -> glam::DVec3 {
        CurveEval::point_at(&self.my_curve, the_u)
    }

    /// OCCT GeomAdaptor_Curve::Resolution(R3D) (GeomAdaptor_Curve.cxx
    /// L1116-1147) — the kernel `Curve3::resolution` re-host carries the
    /// full OCCT switch.
    pub fn resolution(&self, the_r3d: f64) -> f64 {
        self.my_curve.resolution(the_r3d)
    }

    /// OCCT GetType() — the rcad Curve3 enum value (bridge #3: the OCCT
    /// GeomAbs_CurveType dispatch is the enum match in the PrecCurve re-host).
    pub fn curve(&self) -> &Curve3 {
        &self.my_curve
    }
}

impl CurveEval for GeomAdaptorCurve {
    fn point_at(&self, t: f64) -> glam::DVec3 {
        CurveEval::point_at(&self.my_curve, t)
    }
    fn tangent_at(&self, t: f64) -> glam::DVec3 {
        CurveEval::tangent_at(&self.my_curve, t)
    }
    fn default_domain(&self) -> [f64; 2] {
        [self.my_first, self.my_last]
    }
}

// ---------------------------------------------------------------------------
// Geom2dAdaptor_Curve re-host (architecture bridge #2)
// ---------------------------------------------------------------------------

/// OCCT Geom2dAdaptor_Curve (TKG3d/GeomAdaptor) — the adaptor over a
/// `Geom2d_Curve` restricted to [the_u_first, the_u_last].
#[derive(Debug, Clone)]
pub struct Geom2dAdaptorCurve {
    /// OCCT myCurve (the underlying Geom2d_Curve value).
    my_curve: Curve2d,
    /// OCCT myFirst.
    my_first: f64,
    /// OCCT myLast.
    my_last: f64,
}

impl Geom2dAdaptorCurve {
    /// OCCT Geom2dAdaptor_Curve::Geom2dAdaptor_Curve(C, U1, U2).
    pub fn new(the_curve: Curve2d, the_u_first: f64, the_u_last: f64) -> Self {
        Geom2dAdaptorCurve {
            my_curve: the_curve,
            my_first: the_u_first,
            my_last: the_u_last,
        }
    }

    /// OCCT FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT Value(U).
    pub fn value(&self, the_u: f64) -> glam::DVec2 {
        Curve2dEval::point_at(&self.my_curve, the_u)
    }
}

// ---------------------------------------------------------------------------
// GeomAdaptor_Surface re-host (architecture bridge #2)
// ---------------------------------------------------------------------------

/// OCCT GeomAdaptor_Surface (TKG3d/GeomAdaptor) — the adaptor over a
/// `Geom_Surface`. Only the members the `Adaptor3d_CurveOnSurface`
/// re-host consumes are carried.
#[derive(Debug, Clone)]
pub struct GeomAdaptorSurface {
    /// OCCT mySurface (the underlying Geom_Surface value).
    my_surface: Surface3,
}

impl GeomAdaptorSurface {
    /// OCCT GeomAdaptor_Surface::GeomAdaptor_Surface(S).
    pub fn new(the_surface: Surface3) -> Self {
        GeomAdaptorSurface {
            my_surface: the_surface,
        }
    }

    /// OCCT Value(U, V).
    pub fn value(&self, the_u: f64, the_v: f64) -> glam::DVec3 {
        SurfaceEval::point_at(&self.my_surface, the_u, the_v)
    }

    /// OCCT UResolution(R3d) — the kernel `Surface3::u_resolution` re-host
    /// (GeomAdaptor_Surface.cxx L1818-1898).
    pub fn u_resolution(&self, the_r3d: f64) -> f64 {
        self.my_surface.u_resolution(the_r3d)
    }

    /// OCCT VResolution(R3d).
    pub fn v_resolution(&self, the_r3d: f64) -> f64 {
        self.my_surface.v_resolution(the_r3d)
    }

    /// OCCT GetType() — the rcad Surface3 enum value (bridge #3).
    pub fn surface(&self) -> &Surface3 {
        &self.my_surface
    }
}

// ---------------------------------------------------------------------------
// Adaptor3d_CurveOnSurface re-host (architecture bridge #2)
// ---------------------------------------------------------------------------

/// OCCT Adaptor3d_CurveOnSurface (TKG3d/Adaptor3d) — the 3D curve defined as
/// the image of a 2D curve on a surface. Only the members
/// `BRepLib_ValidateEdge` consumes are carried
/// (FirstParameter/LastParameter cxx L977-988, Value/D0, Resolution
/// cxx L1364-1368, GetSurface).
#[derive(Debug, Clone)]
pub struct Adaptor3dCurveOnSurface {
    /// OCCT myCurve (the Geom2dAdaptor_Curve handle).
    my_curve: Geom2dAdaptorCurve,
    /// OCCT mySurface (the GeomAdaptor_Surface handle).
    my_surface: GeomAdaptorSurface,
}

impl Adaptor3dCurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface::Adaptor3d_CurveOnSurface(HC, HS)
    /// (cxx L968-972).
    pub fn new(the_curve: Geom2dAdaptorCurve, the_surface: GeomAdaptorSurface) -> Self {
        Adaptor3dCurveOnSurface {
            my_curve: the_curve,
            my_surface: the_surface,
        }
    }

    /// OCCT FirstParameter() (cxx L977-980) — myCurve->FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.my_curve.first_parameter()
    }

    /// OCCT LastParameter() (cxx L982-985) — myCurve->LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.my_curve.last_parameter()
    }

    /// OCCT Value(U) — D0: myCurve->D0(U, P2d); mySurf->D0(P2d.X, P2d.Y, P).
    pub fn value(&self, the_u: f64) -> glam::DVec3 {
        let a_p2d = self.my_curve.value(the_u);
        self.my_surface.value(a_p2d.x, a_p2d.y)
    }

    /// OCCT Resolution(R3d) (cxx L1364-1368): the pcurve resolution of
    /// min(UResolution, VResolution).
    pub fn resolution(&self, the_r3d: f64) -> f64 {
        let ru = self.my_surface.u_resolution(the_r3d);
        let rv = self.my_surface.v_resolution(the_r3d);
        self.my_curve.resolution_compat(ru.min(rv))
    }

    /// OCCT GetSurface() — the surface adaptor handle.
    pub fn get_surface(&self) -> &GeomAdaptorSurface {
        &self.my_surface
    }

    /// The underlying (curve2d, surface) values (bridge #5: the GAP-carrier
    /// GeomLib_CheckCurveOnSurface consumer signature).
    pub fn curve_on_surface_values(&self) -> (Curve2d, Surface3) {
        (
            self.my_curve.curve2d_clone(),
            self.my_surface.surface_clone(),
        )
    }
}

impl CurveEval for Adaptor3dCurveOnSurface {
    fn point_at(&self, t: f64) -> glam::DVec3 {
        Adaptor3dCurveOnSurface::value(self, t)
    }
    fn tangent_at(&self, t: f64) -> glam::DVec3 {
        // The CurveEval default tangent (normalized derivative) over the
        // chain-ruled point evaluation.
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        [self.my_curve.first_parameter(), self.my_curve.last_parameter()]
    }
}

impl Geom2dAdaptorCurve {
    /// The underlying Geom2d_Curve value (bridge #5).
    fn curve2d_clone(&self) -> Curve2d {
        self.my_curve.clone()
    }

    /// OCCT Geom2dAdaptor_Curve::Resolution(Ruv)
    /// (Geom2dAdaptor_Curve.cxx L1186-1216).  The kernel carries no
    /// `Curve2d::resolution` re-host yet; the Line/Circle/Ellipse arms are
    /// exact, the BSpline/Bezier arms fall to the OCCT default arm
    /// `Precision::Parametric(Ruv)` because the BSplCLib 2D resolution
    /// (ArrayDimension = 2, BSplCLib.cxx L4316-4820) has no kernel port.
    /// GAP: closes with the TKMath BSplCLib 2D-arm batch.
    fn resolution_compat(&self, the_ruv: f64) -> f64 {
        use rcad_kernel::geom::Curve2d as C2;
        use rcad_kernel::precision::parametric_default;
        match &self.my_curve {
            C2::Line(_) => the_ruv,
            C2::Circle(c) => {
                let r = c.radius;
                if r > the_ruv / 2.0 {
                    2.0 * (the_ruv / (2.0 * r)).asin()
                } else {
                    2.0 * std::f64::consts::PI
                }
            }
            C2::Ellipse(e) => the_ruv / e.major_radius,
            C2::BSpline(_) | C2::Bezier(_) => {
                // OCCT: the Geom2d_BSplineCurve/BezierCurve Resolution arms
                // (BSplCLib 2D); GAP-annotated above.
                parametric_default(the_ruv)
            }
            _ => parametric_default(the_ruv),
        }
    }
}

impl GeomAdaptorSurface {
    /// The underlying Geom_Surface value (bridge #5).
    fn surface_clone(&self) -> Surface3 {
        self.my_surface.clone()
    }
}

// ---------------------------------------------------------------------------
// BRepCheck::PrecCurve / PrecSurface re-hosts (architecture bridge #3)
// ---------------------------------------------------------------------------

/// OCCT BRepCheck::PrecCurve(aAC3D) (BRepCheck.cxx L70-97): the maximal
/// relative machine epsilon over the ellipse parameters, RealEpsilon()
/// otherwise.
pub fn brep_check_prec_curve(the_curve: &GeomAdaptorCurve) -> f64 {
    let mut a_xe_max = REAL_EPSILON;
    // GeomAbs_CurveType aCT = aAC3D.GetType();
    if let Curve3::Ellipse(a_el3d) = the_curve.curve() {
        let mut a_x = [
            a_el3d.center.x,
            a_el3d.center.y,
            a_el3d.center.z,
            a_el3d.major_radius,
            a_el3d.minor_radius,
        ];
        a_xe_max = -1.0;
        for a_xi in a_x.iter_mut() {
            if *a_xi < 0.0 {
                *a_xi = -*a_xi;
            }
            let a_xe = standard_real_epsilon(*a_xi);
            if a_xe > a_xe_max {
                a_xe_max = a_xe;
            }
        }
    }
    a_xe_max
}

/// OCCT BRepCheck::PrecSurface(aAHSurf) (BRepCheck.cxx L103-127): the maximal
/// relative machine epsilon over the cone parameters, RealEpsilon()
/// otherwise.
pub fn brep_check_prec_surface(the_surface: &GeomAdaptorSurface) -> f64 {
    let mut a_xe_max = REAL_EPSILON;
    // GeomAbs_SurfaceType aST = aAHSurf->GetType();
    if let Surface3::Cone(a_cone) = the_surface.surface() {
        let mut a_x = [
            a_cone.apex.x,
            a_cone.apex.y,
            a_cone.apex.z,
            a_cone.radius,
        ];
        a_xe_max = -1.0;
        for a_xi in a_x.iter_mut() {
            if *a_xi < 0.0 {
                *a_xi = -*a_xi;
            }
            let a_xe = standard_real_epsilon(*a_xi);
            if a_xe > a_xe_max {
                a_xe_max = a_xe;
            }
        }
    }
    a_xe_max
}

// ---------------------------------------------------------------------------
// Extrema_LocateExtPC re-host (architecture bridge #4)
// ---------------------------------------------------------------------------

/// OCCT Extrema_LocateExtPC (TKGeomBase) — the local Point-to-Curve extremum
/// from a seed parameter, with the OCCT Initialize/Perform/IsDone/
/// SquareDistance contract over the kernel `extrema_locate_ext_pc` Newton
/// search (the CurveEval-generic form).
#[derive(Debug, Clone)]
pub struct ExtremaLocateExtPC {
    /// OCCT myInit (the Initialize call was made).
    my_init: bool,
    /// OCCT myIsDone (the Perform call succeeded).
    my_is_done: bool,
    /// OCCT mySqDist (the square distance at the solution).
    my_sq_dist: f64,
    /// OCCT myPoint (the solution point on the curve).
    my_point: glam::DVec3,
    /// OCCT myPoint.Parameter().
    my_param: f64,
}

impl ExtremaLocateExtPC {
    /// OCCT Extrema_LocateExtPC::Extrema_LocateExtPC() — the default state
    /// (not done).
    pub fn new() -> Self {
        ExtremaLocateExtPC {
            my_init: false,
            my_is_done: false,
            my_sq_dist: 0.0,
            my_point: glam::DVec3::ZERO,
            my_param: 0.0,
        }
    }

    /// OCCT Initialize(C, U1, U2, Tol) — stores the search domain; the
    /// kernel Newton search consumes it in Perform.
    pub fn initialize<C: CurveEval>(&mut self, _the_c: &C, _the_u1: f64, _the_u2: f64, _the_tol: f64) {
        self.my_init = true;
    }

    /// OCCT Perform(P, U0) — the local search from the seed U0.
    pub fn perform<C: CurveEval>(
        &mut self,
        the_p: glam::DVec3,
        the_curve: &C,
        the_u0: f64,
        the_u1: f64,
        the_u2: f64,
        the_tol: f64,
    ) {
        let _ = self.my_init;
        match rcad_kernel::base::extrema::extrema_locate_ext_pc(
            the_p,
            the_curve,
            the_u0,
            the_u1,
            the_u2,
            the_tol,
        ) {
            Some(POnCurve { param, point }) => {
                self.my_is_done = true;
                self.my_point = point;
                self.my_param = param;
                self.my_sq_dist = (point - the_p).length_squared();
            }
            None => {
                self.my_is_done = false;
            }
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT SquareDistance().
    pub fn square_distance(&self) -> f64 {
        self.my_sq_dist
    }

    /// OCCT Point().Parameter() — the parameter of the solution point.
    pub fn param_of_point(&self) -> f64 {
        self.my_param
    }

    /// OCCT Point().
    pub fn point(&self) -> glam::DVec3 {
        self.my_point
    }
}

impl Default for ExtremaLocateExtPC {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// BRepLib_ValidateEdge
// ---------------------------------------------------------------------------

/// OCCT BRepLib_ValidateEdge (hxx L26-103): computes the max distance between
/// the 3D-curve <myReferenceCurve> and the curve on surface <myOtherCurve>.
pub struct BRepLibValidateEdge {
    /// OCCT myReferenceCurve.
    my_reference_curve: GeomAdaptorCurve,
    /// OCCT myOtherCurve.
    my_other_curve: Adaptor3dCurveOnSurface,
    /// OCCT mySameParameter.
    my_same_parameter: bool,
    /// OCCT myControlPointsNumber (hxx L94; initialized to 22 in cxx L30).
    my_control_points_number: i32,
    /// OCCT myToleranceForChecking.
    my_tolerance_for_checking: f64,
    /// OCCT myCalculatedDistance.
    my_calculated_distance: f64,
    /// OCCT myExitIfToleranceExceeded.
    my_exit_if_tolerance_exceeded: bool,
    /// OCCT myIsDone.
    my_is_done: bool,
    /// OCCT myIsExactMethod.
    my_is_exact_method: bool,
    /// OCCT myIsMultiThread.
    my_is_multi_thread: bool,
}

impl BRepLibValidateEdge {
    /// OCCT BRepLib_ValidateEdge::BRepLib_ValidateEdge(theReferenceCurve,
    /// theOtherCurve, theSameParameter) (cxx L24-38).
    pub fn new(
        the_reference_curve: GeomAdaptorCurve,
        the_other_curve: Adaptor3dCurveOnSurface,
        the_same_parameter: bool,
    ) -> Self {
        BRepLibValidateEdge {
            my_reference_curve: the_reference_curve,
            my_other_curve: the_other_curve,
            my_same_parameter: the_same_parameter,
            my_control_points_number: 22,
            my_tolerance_for_checking: 0.0,
            my_calculated_distance: 0.0,
            my_exit_if_tolerance_exceeded: false,
            my_is_done: false,
            my_is_exact_method: false,
            my_is_multi_thread: false,
        }
    }

    /// OCCT SetExactMethod(theIsExact) (hxx L39).
    pub fn set_exact_method(&mut self, the_is_exact: bool) {
        self.my_is_exact_method = the_is_exact;
    }

    /// OCCT IsExactMethod() (hxx L42).
    pub fn is_exact_method(&self) -> bool {
        self.my_is_exact_method
    }

    /// OCCT SetParallel(theIsMultiThread) (hxx L45).
    pub fn set_parallel(&mut self, the_is_multi_thread: bool) {
        self.my_is_multi_thread = the_is_multi_thread;
    }

    /// OCCT IsParallel() (hxx L48).
    pub fn is_parallel(&self) -> bool {
        self.my_is_multi_thread
    }

    /// OCCT SetControlPointsNumber(theControlPointsNumber) (hxx L51-54).
    pub fn set_control_points_number(&mut self, the_control_points_number: i32) {
        self.my_control_points_number = the_control_points_number;
    }

    /// OCCT SetExitIfToleranceExceeded(theToleranceForChecking)
    /// (cxx L81-85).
    pub fn set_exit_if_tolerance_exceeded(&mut self, the_tolerance_for_checking: f64) {
        self.my_exit_if_tolerance_exceeded = true;
        self.my_tolerance_for_checking = self.correct_tolerance(the_tolerance_for_checking);
    }

    /// OCCT Process() (cxx L89-99).
    pub fn process(&mut self) {
        if self.my_is_exact_method && self.my_same_parameter {
            self.process_exact();
        } else {
            self.process_approx();
        }
    }

    /// OCCT IsDone() (hxx L69).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT CheckTolerance(theToleranceToCheck) (cxx L42-45).
    pub fn check_tolerance(&self, the_tolerance_to_check: f64) -> bool {
        self.correct_tolerance(the_tolerance_to_check) > self.my_calculated_distance
    }

    /// OCCT GetMaxDistance() (cxx L49-53).
    pub fn get_max_distance(&self) -> f64 {
        let a_corrected_tolerance = self.my_calculated_distance * 1.00001;
        a_corrected_tolerance
    }

    /// OCCT UpdateTolerance(theToleranceToUpdate) (cxx L57-64).
    pub fn update_tolerance(&self, the_tolerance_to_update: &mut f64) {
        let a_corrected_tolerance = self.my_calculated_distance * 1.00001;
        if a_corrected_tolerance > *the_tolerance_to_update {
            *the_tolerance_to_update = a_corrected_tolerance;
        }
    }

    /// OCCT correctTolerance(theTolerance) (cxx L68-77) — adds some margin
    /// for distance checking.
    fn correct_tolerance(&self, the_tolerance: f64) -> f64 {
        let a_surface = self.my_other_curve.get_surface();
        let a_curve_precision = brep_check_prec_curve(&self.my_reference_curve);
        let a_surface_precision = brep_check_prec_surface(a_surface);
        let a_tolerance_delta = if a_curve_precision > a_surface_precision {
            a_curve_precision
        } else {
            a_surface_precision
        };
        let a_corrected_tolerance = the_tolerance + a_tolerance_delta;
        a_corrected_tolerance
    }

    /// OCCT processApprox() (cxx L103-230) — calculating in a finite number
    /// of points.
    fn process_approx(&mut self) {
        self.my_is_done = true;
        let a_square_tolerance_for_checking =
            self.my_tolerance_for_checking * self.my_tolerance_for_checking;
        let a_reference_first_param = self.my_reference_curve.first_parameter();
        let a_reference_last_param = self.my_reference_curve.last_parameter();
        let an_other_first_param = self.my_other_curve.first_parameter();
        let an_other_last_param = self.my_other_curve.last_parameter();
        let mut a_max_square_distance = 0.0f64;

        let a_control_points_number = if self.my_control_points_number < 1 {
            1
        } else {
            self.my_control_points_number
        };
        let an_is_projection = !self.my_same_parameter
            || (an_other_first_param - a_reference_first_param).abs()
                > rcad_kernel::PCONFUSION
            || (an_other_last_param - a_reference_last_param).abs() > rcad_kernel::PCONFUSION;

        if !an_is_projection {
            for i in 0..=a_control_points_number {
                let a_control_point_param = ((a_control_points_number - i) as f64
                    * a_reference_first_param
                    + i as f64 * a_reference_last_param)
                    / a_control_points_number as f64;
                let a_reference_point = self.my_reference_curve.value(a_control_point_param);
                let an_other_point = self.my_other_curve.value(a_control_point_param);
                let a_square_distance = a_reference_point.distance_squared(an_other_point);
                if a_square_distance > a_max_square_distance {
                    a_max_square_distance = a_square_distance;
                }
                // Stop process for best performance
                if self.my_exit_if_tolerance_exceeded
                    && a_max_square_distance > a_square_tolerance_for_checking
                {
                    self.my_calculated_distance = a_max_square_distance.sqrt();
                    return;
                }
            }
        } else {
            let a_reference_point = self.my_reference_curve.value(a_reference_first_param);
            let an_other_point = self.my_other_curve.value(an_other_first_param);
            let mut a_square_distance = a_reference_point.distance_squared(an_other_point);
            if a_square_distance > a_max_square_distance {
                a_max_square_distance = a_square_distance;
            }
            if self.my_exit_if_tolerance_exceeded
                && a_max_square_distance > a_square_tolerance_for_checking
            {
                self.my_calculated_distance = a_max_square_distance.sqrt();
                return;
            }

            let a_reference_point = self.my_reference_curve.value(a_reference_last_param);
            let an_other_point = self.my_other_curve.value(an_other_last_param);
            a_square_distance = a_reference_point.distance_squared(an_other_point);
            if a_square_distance > a_max_square_distance {
                a_max_square_distance = a_square_distance;
            }
            if self.my_exit_if_tolerance_exceeded
                && a_max_square_distance > a_square_tolerance_for_checking
            {
                self.my_calculated_distance = a_max_square_distance.sqrt();
                return;
            }

            let a_reference_resolution = self.my_reference_curve.resolution(rcad_kernel::CONFUSION);
            let an_other_resolution = self.my_other_curve.resolution(rcad_kernel::CONFUSION);
            let mut a_reference_extrema = ExtremaLocateExtPC::new();
            let mut an_other_extrema = ExtremaLocateExtPC::new();
            a_reference_extrema.initialize(
                &self.my_reference_curve,
                a_reference_first_param,
                a_reference_last_param,
                a_reference_resolution,
            );
            an_other_extrema.initialize(
                &self.my_other_curve,
                an_other_first_param,
                an_other_last_param,
                an_other_resolution,
            );
            for i in 1..a_control_points_number {
                let a_reference_param = ((a_control_points_number - i) as f64
                    * a_reference_first_param
                    + i as f64 * a_reference_last_param)
                    / a_control_points_number as f64;
                let a_reference_extrema_point =
                    self.my_reference_curve.value(a_reference_param);
                let an_other_param = ((a_control_points_number - i) as f64
                    * an_other_first_param
                    + i as f64 * an_other_last_param)
                    / a_control_points_number as f64;
                let an_other_extrema_point = self.my_other_curve.value(an_other_param);

                a_reference_extrema.perform(
                    an_other_extrema_point,
                    &self.my_reference_curve,
                    a_reference_param,
                    a_reference_first_param,
                    a_reference_last_param,
                    a_reference_resolution,
                );
                if a_reference_extrema.is_done() {
                    if a_reference_extrema.square_distance() > a_max_square_distance {
                        a_max_square_distance = a_reference_extrema.square_distance();
                    }
                    if self.my_exit_if_tolerance_exceeded
                        && a_max_square_distance > a_square_tolerance_for_checking
                    {
                        self.my_calculated_distance = a_max_square_distance.sqrt();
                        return;
                    }
                } else {
                    self.my_is_done = false;
                    // Stop process for best performance
                    return;
                }

                an_other_extrema.perform(
                    a_reference_extrema_point,
                    &self.my_other_curve,
                    an_other_param,
                    an_other_first_param,
                    an_other_last_param,
                    an_other_resolution,
                );
                if an_other_extrema.is_done() {
                    if an_other_extrema.square_distance() > a_max_square_distance {
                        a_max_square_distance = an_other_extrema.square_distance();
                    }
                    if self.my_exit_if_tolerance_exceeded
                        && a_max_square_distance > a_square_tolerance_for_checking
                    {
                        self.my_calculated_distance = a_max_square_distance.sqrt();
                        return;
                    }
                } else {
                    self.my_is_done = false;
                    // Stop process for best performance
                    return;
                }
            }
        }
        self.my_calculated_distance = a_max_square_distance.sqrt();
    }

    /// OCCT processExact() (cxx L234-244) — calculating through
    /// GeomLib_CheckCurveOnSurface.
    fn process_exact(&mut self) {
        // Bridge #5: the GAP-carrier GeomLibCheckCurveOnSurface consumes the
        // underlying (Curve3) and (Curve2d, Surface3) values.
        let (a_curve2d, a_surface) = self.my_other_curve.curve_on_surface_values();
        let mut a_check_curve_on_surface =
            GeomLibCheckCurveOnSurface::new(&self.my_reference_curve.curve3_clone());
        a_check_curve_on_surface.set_parallel(self.my_is_multi_thread);
        a_check_curve_on_surface.perform(&(a_curve2d, a_surface));
        self.my_is_done = a_check_curve_on_surface.is_done();
        if self.my_is_done {
            self.my_calculated_distance = a_check_curve_on_surface.max_distance();
        }
    }
}

impl GeomAdaptorCurve {
    /// The underlying Geom_Curve value (bridge #5).
    fn curve3_clone(&self) -> Curve3 {
        self.my_curve.clone()
    }
}
