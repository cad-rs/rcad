//! OCCT Geom2d_BSplineCurve (TKG2d/Geom2d) — faithful port of the B-spline
//! curve data model and operation set consumed by
//! `Geom2dConvert_BSplineCurveToBezierCurve`.
//!
//! 1:1 translation of `Geom2d_BSplineCurve.cxx` (operations) and
//! `Geom2d_BSplineCurve_1.cxx` (read accessors), with the BSplCLib kernels
//! from [`crate::math::bspl_lib`].
//!
//! Architecture decision (commit 0fe5fea): the legacy flat-knot
//! [`BSplineCurve2`] is left untouched; this struct mirrors the OCCT member
//! data one-to-one (degree / distinct knots + mults / poles / weights /
//! periodic / derived flat knot sequence).  OCCT members without an
//! equivalent in this unit's dependency closure are omitted on purpose:
//! `myEvalRep` (deferred evaluation representation) and `myMaxDerivInv` /
//! `myMaxDerivInvOk` (Resolution() cache).

use glam::DVec2;

use crate::geom::BSplineCurve2;
use crate::math::bspl_lib::{
    knot_analysis, knot_sequence, knot_sequence_length, nb_poles, GeomAbsKnotDistribution,
    GP_RESOLUTION, BSPLIB_MAX_DEGREE,
};

/// OCCT GeomAbs_Shape values reachable from `updateKnots` (CN/C0..C3).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SmoothShape {
    CN,
    C0,
    C1,
    C2,
    C3,
}

/// OCCT Geom2d_BSplineCurve.
#[derive(Clone, Debug)]
pub struct Geom2dBSplineCurve {
    my_poles: Vec<DVec2>,
    my_weights: Vec<f64>,
    my_knots: Vec<f64>,
    my_flat_knots: Vec<f64>,
    my_mults: Vec<i32>,
    my_deg: usize,
    my_periodic: bool,
    my_rational: bool,
    my_knot_set: GeomAbsKnotDistribution,
    my_smooth: SmoothShape,
}

/// OCCT Epsilon(Value) (Standard_Real.hxx L242-246) — distance to the nearest
/// double in the direction of infinity with the same sign.
#[inline]
fn epsilon_of(value: f64) -> f64 {
    if value >= 0.0 {
        value.next_up() - value
    } else {
        value - value.next_down()
    }
}

/// OCCT BSplCLib::UnitWeights(N).
fn unit_weights(n: usize) -> Vec<f64> {
    vec![1.0; n]
}

/// OCCT `CheckCurveData` (Geom2d_BSplineCurve.cxx L60-92, static).
fn check_curve_data(
    poles: &[DVec2],
    knots: &[f64],
    mults: &[i32],
    degree: usize,
    periodic: bool,
) {
    if degree < 1 || degree > BSPLIB_MAX_DEGREE {
        panic!("Standard_ConstructionError: BSpline curve: invalid degree");
    }

    if poles.len() < 2 {
        panic!("Standard_ConstructionError: BSpline curve: at least 2 poles required");
    }
    if knots.len() != mults.len() {
        panic!("Standard_ConstructionError: BSpline curve: Knot and Mult array size mismatch");
    }

    let n = knots.len() as i32;
    for i in 1..n {
        if at_knot(knots, i + 1) - at_knot(knots, i) <= epsilon_of(at_knot(knots, i).abs()) {
            panic!("Standard_ConstructionError: BSpline curve: Knots interval values too close");
        }
    }

    if poles.len() != nb_poles(degree, periodic, mults) {
        panic!("Standard_ConstructionError: BSpline curve: # Poles and degree mismatch");
    }
}

/// OCCT `Rational` (Geom2d_BSplineCurve.cxx L95-105, static) — checks
/// rationality of an array of weights by consecutive comparison.
fn is_rational_weights(weights: &[f64]) -> bool {
    let upper = weights.len() as i32;
    for i in 1..upper {
        if (at_weight(weights, i) - at_weight(weights, i + 1)).abs() > GP_RESOLUTION {
            return true;
        }
    }
    false
}

#[inline]
fn at_knot(knots: &[f64], i: i32) -> f64 {
    knots[(i - 1) as usize]
}

#[inline]
fn at_weight(weights: &[f64], i: i32) -> f64 {
    weights[(i - 1) as usize]
}

impl Geom2dBSplineCurve {
    /// OCCT Geom2d_BSplineCurve(Poles, Knots, Mults, Degree, Periodic)
    /// (Geom2d_BSplineCurve.cxx L136-166) — non-rational constructor.
    pub fn new(
        poles: Vec<DVec2>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: usize,
        periodic: bool,
    ) -> Self {
        check_curve_data(&poles, &knots, &mults, degree, periodic);

        let mut curve = Geom2dBSplineCurve {
            my_poles: poles,
            my_weights: unit_weights(1),
            my_knots: knots,
            my_flat_knots: Vec::new(),
            my_mults: mults,
            my_deg: degree,
            my_periodic: periodic,
            my_rational: false,
            my_knot_set: GeomAbsKnotDistribution::NonUniform,
            my_smooth: SmoothShape::CN,
        };
        curve.my_weights = unit_weights(curve.my_poles.len());
        curve.update_knots();
        curve
    }

    /// OCCT Geom2d_BSplineCurve(Poles, Weights, Knots, Mults, Degree, Periodic)
    /// (Geom2d_BSplineCurve.cxx L170-225) — rational constructor.
    pub fn new_rational(
        poles: Vec<DVec2>,
        weights: Vec<f64>,
        knots: Vec<f64>,
        mults: Vec<i32>,
        degree: usize,
        periodic: bool,
    ) -> Self {
        check_curve_data(&poles, &knots, &mults, degree, periodic);

        if weights.len() != poles.len() {
            panic!(
                "Standard_ConstructionError: Geom2d_BSplineCurve: Weights and Poles array size mismatch"
            );
        }

        for i in 0..weights.len() {
            if weights[i] <= GP_RESOLUTION {
                panic!("Standard_ConstructionError: Geom2d_BSplineCurve: Weights values too small");
            }
        }

        // check really rational
        let rational = is_rational_weights(&weights);

        let mut curve = Geom2dBSplineCurve {
            my_poles: poles,
            my_weights: if rational {
                weights
            } else {
                unit_weights(1)
            },
            my_knots: knots,
            my_flat_knots: Vec::new(),
            my_mults: mults,
            my_deg: degree,
            my_periodic: periodic,
            my_rational: rational,
            my_knot_set: GeomAbsKnotDistribution::NonUniform,
            my_smooth: SmoothShape::CN,
        };
        if !curve.my_rational {
            curve.my_weights = unit_weights(curve.my_poles.len());
        }
        curve.update_knots();
        curve
    }

    /// Rebuild a faithful [`Geom2dBSplineCurve`] from the legacy flat-knot
    /// [`BSplineCurve2`] representation (architecture adaptation, see module
    /// docs).  `periodic` must be supplied by the caller: the legacy struct
    /// carries no periodic flag (its flat knot sequence is the expanded
    /// `myKnots x myMults` grouping, so equal consecutive values merge back
    /// into knots/multiplicities).
    pub fn from_bspline2(curve: &BSplineCurve2, periodic: bool) -> Self {
        let mut knots: Vec<f64> = Vec::new();
        let mut mults: Vec<i32> = Vec::new();
        for &k in &curve.knots {
            if knots.last() == Some(&k) {
                *mults.last_mut().expect("non-empty") += 1;
            } else {
                knots.push(k);
                mults.push(1);
            }
        }

        let poles = curve.control_points.clone();
        let rational =
            curve.weights.len() == poles.len() && is_rational_weights(&curve.weights);
        if rational {
            Self::new_rational(
                poles,
                curve.weights.clone(),
                knots,
                mults,
                curve.degree,
                periodic,
            )
        } else {
            Self::new(poles, knots, mults, curve.degree, periodic)
        }
    }

    /// OCCT Geom2d_BSplineCurve::Copy (Geom2d_BSplineCurve.cxx L109-112).
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// OCCT `updateKnots` (Geom2d_BSplineCurve.cxx L1280-1324).
    fn update_knots(&mut self) {
        let mut max_knot_mult = 0i32;
        let mut knot_set = GeomAbsKnotDistribution::NonUniform;
        knot_analysis(
            self.my_deg,
            self.my_periodic,
            &self.my_knots,
            &self.my_mults,
            &mut knot_set,
            &mut max_knot_mult,
        );
        self.my_knot_set = knot_set;

        if knot_set == GeomAbsKnotDistribution::Uniform && !self.my_periodic {
            self.my_flat_knots = self.my_knots.clone();
        } else {
            let len = knot_sequence_length(&self.my_mults, self.my_deg, self.my_periodic);
            self.my_flat_knots = vec![0.0; len];
            knot_sequence(
                &self.my_knots,
                &self.my_mults,
                self.my_deg,
                self.my_periodic,
                &mut self.my_flat_knots,
            );
        }

        self.my_smooth = if max_knot_mult == 0 {
            SmoothShape::CN
        } else {
            match self.my_deg as i32 - max_knot_mult {
                0 => SmoothShape::C0,
                1 => SmoothShape::C1,
                2 => SmoothShape::C2,
                3 => SmoothShape::C3,
                _ => SmoothShape::C3,
            }
        };
    }

    // ------------------------------------------------------------------
    // Read accessors (Geom2d_BSplineCurve_1.cxx)
    // ------------------------------------------------------------------

    /// OCCT Geom2d_BSplineCurve::IsPeriodic (_1.cxx L154-157).
    pub fn is_periodic(&self) -> bool {
        self.my_periodic
    }

    /// OCCT Geom2d_BSplineCurve::Continuity (_1.cxx L161-164).
    pub fn continuity(&self) -> SmoothShape {
        self.my_smooth
    }

    /// OCCT Geom2d_BSplineCurve::Degree (_1.cxx L168-171).
    pub fn degree(&self) -> usize {
        self.my_deg
    }

    /// OCCT Geom2d_BSplineCurve::FirstUKnotIndex (_1.cxx L335-345).
    pub fn first_uknot_index(&self) -> i32 {
        if self.my_periodic {
            1
        } else {
            crate::math::bspl_lib::first_uknot_index_mults(self.my_deg, &self.my_mults)
        }
    }

    /// OCCT Geom2d_BSplineCurve::FirstParameter (_1.cxx L349-352).
    pub fn first_parameter(&self) -> f64 {
        self.my_flat_knots[self.my_deg]
    }

    /// OCCT Geom2d_BSplineCurve::Knot (_1.cxx L356-360).
    pub fn knot(&self, index: i32) -> f64 {
        if index < 1 || index > self.my_knots.len() as i32 {
            panic!("Standard_OutOfRange: Geom2d_BSplineCurve::Knot");
        }
        at_knot(&self.my_knots, index)
    }

    /// OCCT Geom2d_BSplineCurve::KnotDistribution (_1.cxx L364-367).
    pub fn knot_distribution(&self) -> GeomAbsKnotDistribution {
        self.my_knot_set
    }

    /// OCCT Geom2d_BSplineCurve::LastUKnotIndex (_1.cxx L405-415).
    pub fn last_uknot_index(&self) -> i32 {
        if self.my_periodic {
            self.my_knots.len() as i32
        } else {
            crate::math::bspl_lib::last_uknot_index_mults(self.my_deg, &self.my_mults)
        }
    }

    /// OCCT Geom2d_BSplineCurve::LastParameter (_1.cxx L419-422).
    pub fn last_parameter(&self) -> f64 {
        let upper = self.my_flat_knots.len() as i32;
        at_knot(&self.my_flat_knots, upper - self.my_deg as i32)
    }

    /// OCCT Geom2d_BSplineCurve::Multiplicity (_1.cxx L575-580).
    pub fn multiplicity(&self, index: i32) -> i32 {
        if index < 1 || index > self.my_mults.len() as i32 {
            panic!("Standard_OutOfRange: Geom2d_BSplineCurve::Multiplicity");
        }
        self.my_mults[(index - 1) as usize]
    }

    /// OCCT Geom2d_BSplineCurve::NbKnots (_1.cxx L598-601).
    pub fn nb_knots(&self) -> i32 {
        self.my_knots.len() as i32
    }

    /// OCCT Geom2d_BSplineCurve::NbPoles (_1.cxx L605-608).
    pub fn nb_poles_curve(&self) -> i32 {
        self.my_poles.len() as i32
    }

    /// OCCT Geom2d_BSplineCurve::Pole (_1.cxx L612-616).
    pub fn pole(&self, index: i32) -> DVec2 {
        if index < 1 || index > self.my_poles.len() as i32 {
            panic!("Standard_OutOfRange: Geom2d_BSplineCurve::Pole");
        }
        self.my_poles[(index - 1) as usize]
    }

    /// OCCT Geom2d_BSplineCurve::Weight (_1.cxx L647-659).
    pub fn weight(&self, index: i32) -> f64 {
        if index < 1 || index > self.my_poles.len() as i32 {
            panic!("Standard_OutOfRange: Geom2d_BSplineCurve::Weight");
        }
        if self.my_rational {
            at_weight(&self.my_weights, index)
        } else {
            1.0
        }
    }

    /// OCCT Geom2d_BSplineCurve::IsRational (_1.cxx L691-694).
    pub fn is_rational(&self) -> bool {
        self.my_rational
    }

    /// OCCT Geom2d_BSplineCurve::PeriodicNormalization
    /// (Geom2d_BSplineCurve.cxx L1328-1344).
    pub fn periodic_normalization(&self, parameter: &mut f64) {
        if self.my_periodic {
            let upper = self.my_flat_knots.len() as i32;
            let deg = self.my_deg as i32;
            let period =
                at_knot(&self.my_flat_knots, upper - deg) - at_knot(&self.my_flat_knots, deg + 1);
            while *parameter > at_knot(&self.my_flat_knots, upper - deg) {
                *parameter -= period;
            }
            while *parameter < at_knot(&self.my_flat_knots, deg + 1) {
                *parameter += period;
            }
        }
    }
}
