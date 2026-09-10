//! OCCT GeomAPI_PointsToBSpline (TKGeomAlgo GeomAPI/GeomAPI_PointsToBSpline.cxx
//! L17-357) — 1:1 translation, escalated into the W1-4 batch per the
//! TKShHealing Path-B docket (dependency row `GeomAPI_PointsToBSpline`:
//! "verify at batch start; escalate to a new front batch if missing"; the
//! rcad kernel only had fitting utilities, no 1:1 port).
//!
//! The single ProjectCurveOnSurface call site (approximatePCurve, OCCT
//! L2252) uses the `(Points, Params, DegMin, DegMax, Continuity, Tol3D)`
//! constructor, whose `Init` drives `AppDef_BSplineCompute` — available as
//! `geomalgo::bspl_compute_line::BSplineCompute`.  The fourth `Init`
//! (`W1, W2, W3`) drives `AppDef_Variational`, which rcad has not ported;
//! that method keeps OCCT's failure path (returns with `IsDone() == false`)
//! and is marked as a GAP below.

use glam::DVec3;
use rcad_kernel::math::bspl_lib::reparametrize;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::app_def::MultiLine;
use crate::geomalgo::approx_int::ApproxParamType;
use crate::geomalgo::bspl_compute_line::BSplineCompute;

use super::array1::Array1;

// OCCT L33-40: anonymous-namespace helper.
fn has_meaningful_span(the_span: f64) -> bool {
    let a_resolution = rcad_kernel::math::gp::GP_RESOLUTION;
    the_span * the_span > a_resolution * a_resolution
}

/// Architecture bridge: OCCT `Geom_BSplineCurve(Poles, Knots, Mults, Degree)`
/// carries (distinct knots, multiplicities); the rcad `BSplineCurve3` encoding
/// stores the flat expanded knot vector.  Expands `the_mults` over
/// `the_knots`.
fn expand_knots(the_knots: &[f64], the_mults: &[usize]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (i, &k) in the_knots.iter().enumerate() {
        let m = the_mults.get(i).copied().unwrap_or(1);
        for _ in 0..m {
            flat.push(k);
        }
    }
    flat
}

/// OCCT GeomAPI_PointsToBSpline class.
#[derive(Debug, Clone, Default)]
pub struct PointsToBSpline {
    /// OCCT myCurve.
    my_curve: Option<rcad_kernel::geom::BSplineCurve3>,
    /// OCCT myIsDone.
    my_is_done: bool,
}

impl PointsToBSpline {
    /// OCCT L44-47: default constructor.
    pub fn new() -> Self {
        PointsToBSpline {
            my_curve: None,
            my_is_done: false,
        }
    }

    /// OCCT L51-59: constructor (Points, DegMin, DegMax, Continuity, Tol3D).
    pub fn with_points(
        the_points: &Array1<DVec3>,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) -> Self {
        let mut res = PointsToBSpline::new();
        res.init(the_points, the_deg_min, the_deg_max, the_continuity, the_tol3d);
        res
    }

    /// OCCT L63-72: constructor (Points, ParType, DegMin, DegMax, Continuity,
    /// Tol3D).
    pub fn with_points_par_type(
        the_points: &Array1<DVec3>,
        the_par_type: ApproxParamType,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) -> Self {
        let mut res = PointsToBSpline::new();
        res.init_par_type(
            the_points,
            the_par_type,
            the_deg_min,
            the_deg_max,
            the_continuity,
            the_tol3d,
        );
        res
    }

    /// OCCT L76-85: constructor (Points, Params, DegMin, DegMax, Continuity,
    /// Tol3D).
    pub fn with_points_params(
        the_points: &Array1<DVec3>,
        the_params: &Array1<f64>,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) -> Self {
        let mut res = PointsToBSpline::new();
        res.init_params(
            the_points,
            the_params,
            the_deg_min,
            the_deg_max,
            the_continuity,
            the_tol3d,
        );
        res
    }

    /// OCCT L89-99: constructor (Points, W1, W2, W3, DegMax, Continuity, Tol3D).
    pub fn with_points_weights(
        the_points: &Array1<DVec3>,
        the_w1: f64,
        the_w2: f64,
        the_w3: f64,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) -> Self {
        let mut res = PointsToBSpline::new();
        res.init_weights(
            the_points, the_w1, the_w2, the_w3, the_deg_max, the_continuity, the_tol3d,
        );
        res
    }

    /// OCCT L103-111: Init(Points, DegMin, DegMax, Continuity, Tol3D).
    pub fn init(
        &mut self,
        the_points: &Array1<DVec3>,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) {
        self.my_is_done = false;
        self.init_par_type(
            the_points,
            ApproxParamType::ChordLength,
            the_deg_min,
            the_deg_max,
            the_continuity,
            the_tol3d,
        );
    }

    /// OCCT L115-165: Init(Points, ParType, DegMin, DegMax, Continuity, Tol3D).
    pub fn init_par_type(
        &mut self,
        the_points: &Array1<DVec3>,
        the_par_type: ApproxParamType,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) {
        self.my_is_done = false;
        let the_tol2d = 0.0; // dummy argument for BSplineCompute.

        let nbit = 2;
        let mut use_squares = false;
        if the_tol3d <= 1.0e-3 {
            use_squares = true;
        }

        let mut the_computer = BSplineCompute::new(
            the_deg_min,
            the_deg_max,
            the_tol3d,
            the_tol2d,
            nbit,
            true,
            the_par_type,
            use_squares,
        );

        the_computer.set_continuity(continuity_order(the_continuity));

        // OCCT TheComputer.Perform(Points) wraps the point array into an
        // AppDef_MultiLine; rcad wraps via the AppDef_MultiLine(tabP3d) ctor.
        let a_tab: Vec<DVec3> = (the_points.lower()..=the_points.upper())
            .map(|i| *the_points.value(i))
            .collect();
        the_computer.perform(&MultiLine::new_tab_p3d(&a_tab));

        let the_curve = the_computer.value();

        let mut poles: Vec<DVec3> = Vec::new();
        the_curve.curve(1, &mut poles);

        self.my_curve = Some(rcad_kernel::geom::BSplineCurve3 {
            degree: the_curve.degree,
            // Architecture bridge: OCCT passes (Knots, Multiplicities); rcad
            // stores the flat expanded knot vector.
            knots: expand_knots(&the_curve.knots, &the_curve.mults),
            control_points: poles,
            weights: vec![1.0; the_curve.nb_poles()],
            is_periodic: false,
        });
        self.my_is_done = true;
    }

    /// OCCT L169-238: Init(Points, Params, DegMin, DegMax, Continuity, Tol3D).
    pub fn init_params(
        &mut self,
        the_points: &Array1<DVec3>,
        the_params: &Array1<f64>,
        the_deg_min: i32,
        the_deg_max: i32,
        the_continuity: GeomAbsShape,
        the_tol3d: f64,
    ) {
        self.my_is_done = false;
        if the_params.length() != the_points.length() {
            // OCCT throws Standard_OutOfRange("GeomAPI_PointsToBSpline::Init()
            // - invalid input").
            panic!("GeomAPI_PointsToBSpline::Init() - invalid input");
        }

        let the_tol2d = 0.0; // dummy argument for BSplineCompute.
        let nbp = the_params.length() as i32;
        let mut the_params_vec = Vector::new(1, nbp);
        the_params_vec.set(1, 0.0);
        the_params_vec.set(nbp, 1.0);

        let uf = *the_params.value(the_params.lower());
        let ul = *the_params.value(the_params.upper()) - uf;
        if !has_meaningful_span(ul) {
            return;
        }

        for i in 2..nbp {
            the_params_vec.set(i, (*the_params.value(i as usize) - uf) / ul);
        }

        let mut the_computer = BSplineCompute::new(
            the_deg_min,
            the_deg_max,
            the_tol3d,
            the_tol2d,
            0,
            true,
            ApproxParamType::IsoParametric,
            true,
        );

        the_computer.set_parameters(&the_params_vec);

        the_computer.set_continuity(continuity_order(the_continuity));

        let a_tab: Vec<DVec3> = (the_points.lower()..=the_points.upper())
            .map(|i| *the_points.value(i))
            .collect();
        the_computer.perform(&MultiLine::new_tab_p3d(&a_tab));

        let the_curve = the_computer.value();

        let mut poles: Vec<DVec3> = Vec::new();
        the_curve.curve(1, &mut poles);

        // OCCT Knots(TheCurve.Knots().Lower(), TheCurve.Knots().Upper());
        // Knots = TheCurve.Knots();
        let mut knots = the_curve.knots.clone();
        reparametrize(
            *the_params.value(the_params.lower()),
            *the_params.value(the_params.upper()),
            &mut knots,
        );

        self.my_curve = Some(rcad_kernel::geom::BSplineCurve3 {
            degree: the_curve.degree,
            // Architecture bridge: OCCT passes (Knots, Multiplicities); rcad
            // stores the flat expanded knot vector.
            knots: expand_knots(&knots, &the_curve.mults),
            control_points: poles,
            weights: vec![1.0; the_curve.nb_poles()],
            is_periodic: false,
        });
        self.my_is_done = true;
    }

    /// OCCT L242-329: Init(Points, W1, W2, W3, DegMax, Continuity, Tol3D).
    ///
    /// GAP: AppDef_Variational (the variational-smoothing approximation
    /// AppDef/AppParCurves stack member) has no rcad port.  OCCT's failure
    /// path is kept: the method returns with IsDone() == false exactly like a
    /// non-created / over-constrained / failed / not-done Variational run.
    /// The OCCT call shape is recorded for the future 1:1 port:
    ///   AppDef_MultiLine multL(NbPoint) of MultiPointConstraint(1,0) set to
    ///   the points; HArray1<AppParCurves_ConstraintCouple> TABofCC filled
    ///   with (i, AppParCurves_NoConstraint);
    ///   AppDef_Variational Variation(multL, 1, NbPoint, TABofCC);
    ///   SetMaxDegree(DegMax); SetContinuity(Continuity);
    ///   SetMaxSegment(1000); SetTolerance(Tol3D); SetWithMinMax(false);
    ///   SetNbIterations(nbit)  (nbit = 2, or 0 when Tol3D <= 1.e-3);
    ///   SetCriteriumWeight(W1, W2, W3);
    ///   bail-outs on !IsCreated() / IsOverConstrained() / Approximate()
    ///   exception / !IsDone(); TheCurve = Variation.Value().
    pub fn init_weights(
        &mut self,
        _the_points: &Array1<DVec3>,
        _the_w1: f64,
        _the_w2: f64,
        _the_w3: f64,
        _the_deg_max: i32,
        _the_continuity: GeomAbsShape,
        _the_tol3d: f64,
    ) {
        self.my_is_done = false;
        // GAP: AppDef_Variational not ported — OCCT failure path (IsDone false).
    }

    /// OCCT L336-343: Curve() — raises StdFail_NotDone when not done; rcad
    /// mirrors the raise as a panic.
    pub fn curve(&self) -> &rcad_kernel::geom::BSplineCurve3 {
        if !self.my_is_done {
            panic!("StdFail_NotDone: GeomAPI_PointsToBSpline::Curve ");
        }
        self.my_curve.as_ref().expect("GeomAPI_PointsToBSpline::Curve ")
    }

    // OCCT L347-350: `operator occ::handle<Geom_BSplineCurve>()` — the C++
    // implicit conversion operator has no Rust equivalent; `curve()` covers
    // the same access.

    /// OCCT L354-357: IsDone().
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }
}

/// OCCT L134-152 / L205-223: the Continuity → SetContinuity switch.
/// OCCT distinguishes GeomAbs_G1/G2; the rcad `GeomAbsShape` enum has no G
/// arms, so only the C orders are representable.
fn continuity_order(the_continuity: GeomAbsShape) -> i32 {
    match the_continuity {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::C1 => 1,
        GeomAbsShape::C2 => 2,
        _ => 3,
    }
}
