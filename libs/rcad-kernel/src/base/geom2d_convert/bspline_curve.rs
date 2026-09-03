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
fn set_at_knot(knots: &mut [f64], i: i32, v: f64) {
    knots[(i - 1) as usize] = v;
}

#[inline]
fn set_at_mult(mults: &mut [i32], i: i32, v: i32) {
    mults[(i - 1) as usize] = v;
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

    // ------------------------------------------------------------------
    // Modification operations (Geom2d_BSplineCurve.cxx)
    // ------------------------------------------------------------------

    /// OCCT Geom2d_BSplineCurve::IncreaseDegree (Geom2d_BSplineCurve.cxx
    /// L236-293).  The BSplCLib gp_Pnt2d overload (BSplCLib_1.cxx L153-178)
    /// dispatches through `BSplCLib_IncreaseDegree`
    /// (BSplCLib_CurveComputation.pxx L538-590): rational curves are
    /// homogenized to dim 3 (PLib::SetPoles), the flat Dimension kernel runs,
    /// and the result is dehomogenized (PLib::GetPoles).
    pub fn increase_degree(&mut self, degree: usize) {
        if degree == self.my_deg {
            return;
        }

        if degree < self.my_deg || degree > BSPLIB_MAX_DEGREE {
            panic!("Standard_ConstructionError: BSpline curve: IncreaseDegree: bad degree value");
        }

        let from_k1 = self.first_uknot_index();
        let to_k2 = self.last_uknot_index();

        let step = (degree - self.my_deg) as i32;

        // OCCT npoles(1, myPoles.Length() + Step * (ToK2 - FromK1)) in poles;
        // the flat buffer holds `dim` doubles per pole.
        let npoles_len = self.my_poles.len() + (step * (to_k2 - from_k1)) as usize;

        let nbknots = crate::math::bspl_lib::increase_degree_count_knots(
            self.my_deg,
            degree,
            self.my_periodic,
            &self.my_mults,
        );

        let weights_opt = if self.my_rational {
            Some(&self.my_weights[..])
        } else {
            None
        };
        let (poles_flat, dim) = flatten_poles(&self.my_poles, weights_opt);
        let mut npoles_flat = vec![0.0f64; npoles_len * dim];
        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];

        crate::math::bspl_lib::increase_degree(
            self.my_deg,
            degree,
            self.my_periodic,
            dim,
            &poles_flat,
            &self.my_knots,
            &self.my_mults,
            &mut npoles_flat,
            &mut nknots,
            &mut nmults,
        );

        self.my_deg = degree;
        let (poles, weights) = unflatten_poles(&npoles_flat, dim);
        self.my_poles = poles;
        if self.my_rational {
            self.my_weights = weights;
        } else {
            self.my_weights = unit_weights(self.my_poles.len());
        }
        self.my_knots = nknots;
        self.my_mults = nmults;
        self.update_knots();
    }

    /// OCCT Geom2d_BSplineCurve::IncreaseMultiplicity(Index, M)
    /// (Geom2d_BSplineCurve.cxx L297-304).
    pub fn increase_multiplicity(&mut self, index: i32, m: i32) {
        let k = [at_knot(&self.my_knots, index)];
        let m_arr = [m - self.my_mults[(index - 1) as usize]];
        // OCCT passes Epsilon(1.) — machine epsilon at 1.0 == f64::EPSILON.
        self.insert_knots(&k, &m_arr, f64::EPSILON, true);
    }

    /// OCCT Geom2d_BSplineCurve::IncreaseMultiplicity(I1, I2, M)
    /// (Geom2d_BSplineCurve.cxx L308-318).
    pub fn increase_multiplicity_range(&mut self, i1: i32, i2: i32, m: i32) {
        let mut k = Vec::with_capacity((i2 - i1 + 1) as usize);
        let mut m_arr = Vec::with_capacity((i2 - i1 + 1) as usize);
        for i in i1..=i2 {
            k.push(at_knot(&self.my_knots, i));
            m_arr.push(m - self.my_mults[(i - 1) as usize]);
        }
        self.insert_knots(&k, &m_arr, f64::EPSILON, true);
    }

    /// OCCT Geom2d_BSplineCurve::InsertKnots(Knots, Mults, Epsilon, Add)
    /// (Geom2d_BSplineCurve.cxx L343-406).  Dispatch through
    /// `BSplCLib_InsertKnots` (BSplCLib_CurveComputation.pxx L381-437):
    /// rational curves homogenized to dim 3.
    pub fn insert_knots(&mut self, knots: &[f64], mults: &[i32], epsilon: f64, add: bool) {
        // Check and compute new sizes
        let mut nbpoles = 0i32;
        let mut nbknots = 0i32;

        if !crate::math::bspl_lib::prepare_insert_knots(
            self.my_deg,
            self.my_periodic,
            &self.my_knots,
            &self.my_mults,
            knots,
            Some(mults),
            &mut nbpoles,
            &mut nbknots,
            epsilon,
            add,
        ) {
            panic!("Standard_ConstructionError: Geom2d_BSplineCurve::InsertKnots");
        }

        if nbpoles as usize == self.my_poles.len() {
            return;
        }

        let weights_opt = if self.my_rational {
            Some(&self.my_weights[..])
        } else {
            None
        };
        let (poles_flat, dim) = flatten_poles(&self.my_poles, weights_opt);
        let mut npoles_flat = vec![0.0f64; nbpoles as usize * dim];
        let mut nknots = vec![0.0f64; nbknots as usize];
        let mut nmults = vec![0i32; nbknots as usize];

        crate::math::bspl_lib::insert_knots(
            self.my_deg,
            self.my_periodic,
            dim,
            &poles_flat,
            &self.my_knots,
            &self.my_mults,
            knots,
            Some(mults),
            &mut npoles_flat,
            &mut nknots,
            &mut nmults,
            epsilon,
            add,
        );

        let (poles, weights) = unflatten_poles(&npoles_flat, dim);
        if self.my_rational {
            self.my_weights = weights;
        } else {
            self.my_weights = unit_weights(nbpoles as usize);
        }
        self.my_poles = poles;
        self.my_knots = nknots;
        self.my_mults = nmults;
        self.update_knots();
    }

    /// OCCT Geom2d_BSplineCurve::Segment(U1, U2, theTolerance)
    /// (Geom2d_BSplineCurve.cxx L707-888).
    pub fn segment(&mut self, au1: f64, au2: f64, tolerance: f64) {
        if au2 < au1 {
            panic!("Standard_DomainError: Geom2d_BSplineCurve::Segment");
        }

        let mut new_u1 = 0.0f64;
        let mut new_u2 = 0.0f64;
        let mut u = 0.0f64;
        let mut du = 0.0f64;
        let mut adddu = 0.0f64;
        let was_periodic = self.my_periodic;

        let u1 = au1;
        let u2 = au2;

        let mut knots2 = [0.0f64; 2];
        let mults2 = [self.my_deg as i32, self.my_deg as i32];

        // define param distance to keep (eap, Apr 18 2002, occ311)
        if self.my_periodic {
            let period = self.last_parameter() - self.first_parameter();
            du = u2 - u1;
            if du - period > crate::core::precision::PCONFUSION {
                panic!("Standard_DomainError: Geom2d_BSplineCurve::Segment");
            }
            if du > period {
                du = period;
            }
            adddu = du;
        }

        let knots_upper = self.my_knots.len() as i32;
        let mut index = 0i32;
        crate::math::bspl_lib::locate_parameter_knots_mults(
            self.my_deg,
            &self.my_knots,
            &self.my_mults,
            u1,
            self.my_periodic,
            1,
            knots_upper,
            &mut index,
            &mut new_u1,
        );
        let mut index = 0i32;
        crate::math::bspl_lib::locate_parameter_knots_mults(
            self.my_deg,
            &self.my_knots,
            &self.my_mults,
            u2,
            self.my_periodic,
            1,
            knots_upper,
            &mut index,
            &mut new_u2,
        );

        //-- DBB
        let a_nu2 = new_u2;
        //-- DBB

        let mut abs_umax = new_u1.abs().max(new_u2.abs());
        abs_umax = abs_umax.max(self.first_parameter().abs().max(self.last_parameter().abs()));
        let eps = epsilon_of(abs_umax).max(tolerance);

        knots2[0] = new_u1.min(new_u2);
        knots2[1] = new_u1.max(new_u2);
        self.insert_knots(&knots2, &mults2, eps, false);

        if self.my_periodic {
            // set the origine at NewU1
            let mut index = 0i32;
            let knots_upper = self.my_knots.len() as i32;
            crate::math::bspl_lib::locate_parameter_knots_mults(
                self.my_deg,
                &self.my_knots,
                &self.my_mults,
                u1,
                self.my_periodic,
                1,
                knots_upper,
                &mut index,
                &mut u,
            );
            if (at_knot(&self.my_knots, index + 1) - u).abs() <= eps {
                index += 1;
            }
            self.set_origin(index);
            self.set_not_periodic();
            new_u2 = new_u1 + du;
        }

        // compute index1 and index2 to set the new knots and mults
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        let from_u1 = 1i32;
        let to_u2 = self.my_knots.len() as i32;
        let mut u = 0.0f64;
        crate::math::bspl_lib::locate_parameter_knots_mults(
            self.my_deg,
            &self.my_knots,
            &self.my_mults,
            new_u1,
            self.my_periodic,
            from_u1,
            to_u2,
            &mut index1,
            &mut u,
        );
        if (at_knot(&self.my_knots, index1 + 1) - u).abs() <= eps {
            index1 += 1;
        }
        crate::math::bspl_lib::locate_parameter_knots_mults(
            self.my_deg,
            &self.my_knots,
            &self.my_mults,
            new_u2,
            self.my_periodic,
            from_u1,
            to_u2,
            &mut index2,
            &mut u,
        );
        if (at_knot(&self.my_knots, index2 + 1) - u).abs() <= eps || index2 == index1 {
            index2 += 1;
        }

        let nbknots = (index2 - index1 + 1) as usize;
        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];

        // to restore changed U1
        if du > 0.0 {
            // if was periodic
            du = new_u1 - u1;
        }

        let mut k = 1i32;
        for i in index1..=index2 {
            set_at_knot(&mut nknots, k, at_knot(&self.my_knots, i) - du);
            set_at_mult(&mut nmults, k, self.my_mults[(i - 1) as usize]);
            k += 1;
        }
        set_at_mult(&mut nmults, 1, self.my_deg as i32 + 1);
        set_at_mult(&mut nmults, nbknots as i32, self.my_deg as i32 + 1);

        // compute index1 and index2 to set the new poles and weights
        let mut pindex1 =
            crate::math::bspl_lib::pole_index(self.my_deg, index1, self.my_periodic, &self.my_mults);
        let mut pindex2 =
            crate::math::bspl_lib::pole_index(self.my_deg, index2, self.my_periodic, &self.my_mults);

        pindex1 += 1;
        pindex2 = (pindex2 + 1).min(self.my_poles.len() as i32);

        let nbpoles = (pindex2 - pindex1 + 1) as usize;
        let mut npoles: Vec<DVec2> = Vec::with_capacity(nbpoles);
        let mut nweights: Vec<f64> = Vec::new();

        if self.my_rational {
            nweights.reserve(nbpoles);
            for i in pindex1..=pindex2 {
                npoles.push(self.my_poles[(i - 1) as usize]);
                nweights.push(at_weight(&self.my_weights, i));
            }
        } else {
            for i in pindex1..=pindex2 {
                npoles.push(self.my_poles[(i - 1) as usize]);
            }
        }

        //-- DBB
        if was_periodic {
            set_at_knot(&mut nknots, 1, u1);
            if a_nu2 < u2 {
                set_at_knot(&mut nknots, nbknots as i32, u1 + adddu);
            }
        }
        //-- DBB

        self.my_knots = nknots;
        self.my_mults = nmults;
        self.my_poles = npoles;
        if self.my_rational {
            self.my_weights = nweights;
        } else {
            self.my_weights = unit_weights(self.my_poles.len());
        }
        self.update_knots();
    }

    /// OCCT Geom2d_BSplineCurve::SetOrigin (Geom2d_BSplineCurve.cxx
    /// L989-1083).
    pub fn set_origin(&mut self, index: i32) {
        if !self.my_periodic {
            panic!("Standard_NoSuchObject: Geom2d_BSplineCurve::SetOrigin");
        }
        let first = self.first_uknot_index();
        let last = self.last_uknot_index();

        if (index < first) || (index > last) {
            panic!("Standard_DomainError: Geom2d_BSplineCurve::SetOrigin");
        }

        let nbknots = self.my_knots.len();
        let nbpoles = self.my_poles.len();

        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];

        // set the knots and mults
        let period = at_knot(&self.my_knots, last) - at_knot(&self.my_knots, first);
        let mut k = 1i32;
        for i in index..=last {
            set_at_knot(&mut nknots, k, at_knot(&self.my_knots, i));
            set_at_mult(&mut nmults, k, self.my_mults[(i - 1) as usize]);
            k += 1;
        }
        for i in (first + 1)..=index {
            set_at_knot(&mut nknots, k, at_knot(&self.my_knots, i) + period);
            set_at_mult(&mut nmults, k, self.my_mults[(i - 1) as usize]);
            k += 1;
        }

        // OCCT: int index = 1;
        //       for (i = first + 1; i <= Index; i++) index += myMults.Value(i);
        let origin = index;
        let mut index = 1i32;
        for i in (first + 1)..=origin {
            index += self.my_mults[(i - 1) as usize];
        }
        // set the poles and weights
        let mut npoles: Vec<DVec2> = vec![DVec2::ZERO; nbpoles];
        let mut nweights: Vec<f64> = Vec::new();
        let first = 1i32;
        let last = nbpoles as i32;
        if self.my_rational {
            nweights.resize(nbpoles, 0.0);
            let mut k = 1i32;
            for i in index..=last {
                npoles[(k - 1) as usize] = self.my_poles[(i - 1) as usize];
                nweights[(k - 1) as usize] = at_weight(&self.my_weights, i);
                k += 1;
            }
            for i in first..index {
                npoles[(k - 1) as usize] = self.my_poles[(i - 1) as usize];
                nweights[(k - 1) as usize] = at_weight(&self.my_weights, i);
                k += 1;
            }
        } else {
            let mut k = 1i32;
            for i in index..=last {
                npoles[(k - 1) as usize] = self.my_poles[(i - 1) as usize];
                k += 1;
            }
            for i in first..index {
                npoles[(k - 1) as usize] = self.my_poles[(i - 1) as usize];
                k += 1;
            }
        }

        self.my_poles = npoles;
        self.my_knots = nknots;
        self.my_mults = nmults;
        if self.my_rational {
            self.my_weights = nweights;
        } else {
            self.my_weights = unit_weights(nbpoles);
        }
        self.update_knots();
    }

    /// OCCT Geom2d_BSplineCurve::SetNotPeriodic
    /// (Geom2d_BSplineCurve.cxx L1087-1130).  Dispatch through
    /// `BSplCLib_Unperiodize` (BSplCLib_CurveComputation.pxx L602-642):
    /// rational curves homogenized to dim 3.
    pub fn set_not_periodic(&mut self) {
        if self.my_periodic {
            let mut nb_knots = 0i32;
            let mut nb_poles = 0i32;
            crate::math::bspl_lib::prepare_unperiodize(
                self.my_deg,
                &self.my_mults,
                &mut nb_knots,
                &mut nb_poles,
            );

            let weights_opt = if self.is_rational() {
                Some(&self.my_weights[..])
            } else {
                None
            };
            let (poles_flat, dim) = flatten_poles(&self.my_poles, weights_opt);
            let mut npoles_flat = vec![0.0f64; nb_poles as usize * dim];
            let mut nknots = vec![0.0f64; nb_knots as usize];
            let mut nmults = vec![0i32; nb_knots as usize];

            crate::math::bspl_lib::unperiodize(
                self.my_deg,
                &self.my_mults,
                &self.my_knots,
                &poles_flat,
                &mut nmults,
                &mut nknots,
                &mut npoles_flat,
            );

            let (poles, weights) = unflatten_poles(&npoles_flat, dim);
            self.my_poles = poles;
            if self.my_rational {
                self.my_weights = weights;
            } else {
                self.my_weights = unit_weights(self.my_poles.len());
            }
            self.my_mults = nmults;
            self.my_knots = nknots;
            self.my_periodic = false;
            self.update_knots();
        }
    }

    /// OCCT Geom2d_BSplineCurve::LocateU (Geom2d_BSplineCurve_1.cxx L710-760).
    pub fn locate_u(
        &self,
        u: f64,
        parametric_tolerance: f64,
        i1: &mut i32,
        i2: &mut i32,
        with_knot_repetition: bool,
    ) {
        let mut new_u = u;

        let cknots: &Vec<f64> = if with_knot_repetition {
            &self.my_flat_knots
        } else {
            &self.my_knots
        };

        self.periodic_normalization(&mut new_u); // Attention a la periode
        let knots_len = cknots.len() as i32;
        let ufirst = at_knot(cknots, 1);
        let ulast = at_knot(cknots, knots_len);
        let p_parametric_tolerance = parametric_tolerance.abs();
        if (new_u - ufirst).abs() <= p_parametric_tolerance {
            *i1 = 1;
            *i2 = 1;
        } else if (new_u - ulast).abs() <= p_parametric_tolerance {
            *i1 = knots_len;
            *i2 = knots_len;
        } else if new_u < ufirst {
            *i2 = 1;
            *i1 = 0;
        } else if new_u > ulast {
            *i1 = knots_len;
            *i2 = *i1 + 1;
        } else {
            *i1 = 1;
            crate::math::bspl_lib::hunt(cknots, new_u, i1);
            *i1 = (*i1).min(knots_len).max(1);
            while *i1 + 1 <= knots_len
                && (at_knot(cknots, *i1 + 1) - new_u).abs() <= p_parametric_tolerance
            {
                *i1 += 1;
            }
            if (at_knot(cknots, *i1) - new_u).abs() <= p_parametric_tolerance {
                *i2 = *i1;
            } else {
                *i2 = *i1 + 1;
            }
        }
    }
}

/// OCCT PLib::SetPoles — flatten poles to the flat kernel representation.
/// Non-rational: (x, y) per pole, dim 2.  Rational (PLib::SetPoles(Poles,
/// Weights, FPoles)): homogeneous (x*w, y*w, w) per pole, dim 3.
fn flatten_poles(poles: &[DVec2], weights: Option<&[f64]>) -> (Vec<f64>, usize) {
    match weights {
        None => {
            let mut flat = Vec::with_capacity(poles.len() * 2);
            for p in poles {
                flat.push(p.x);
                flat.push(p.y);
            }
            (flat, 2)
        }
        Some(w) => {
            let mut flat = Vec::with_capacity(poles.len() * 3);
            for (p, &wi) in poles.iter().zip(w.iter()) {
                flat.push(p.x * wi);
                flat.push(p.y * wi);
                flat.push(wi);
            }
            (flat, 3)
        }
    }
}

/// OCCT PLib::GetPoles — restore poles (and weights when homogeneous) from
/// the flat kernel representation.
fn unflatten_poles(flat: &[f64], dim: usize) -> (Vec<DVec2>, Vec<f64>) {
    if dim == 3 {
        let n = flat.len() / 3;
        let mut poles = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        for c in flat.chunks_exact(3) {
            weights.push(c[2]);
            poles.push(DVec2::new(c[0] / c[2], c[1] / c[2]));
        }
        (poles, weights)
    } else {
        let poles = flat
            .chunks_exact(2)
            .map(|c| DVec2::new(c[0], c[1]))
            .collect();
        (poles, Vec::new())
    }
}
