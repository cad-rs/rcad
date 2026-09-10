//! OCCT Law_BSpline (TKGeomAlgo/Law) — 1:1 port of Law_BSpline.hxx
//! (L94-112 members) and Law_BSpline.cxx (whole file L26-1770), the 1D
//! B-spline evolution law curve.
//!
//! Architecture mapping: `NCollection_Array1<double>` poles/knots /
//! `Array1<int>` mults -> `Vec<f64>` / `Vec<i32>` with OCCT 1-based
//! indexing inside loops; the `weights` null handle -> `Option<Vec<f64>>`;
//! `GeomAbs_BSplKnotDistribution` and the `smooth` (GeomAbs_Shape) values
//! are local enums (`KnotDistribution` / `SmoothShape`), mirroring the
//! Geom2dBSplineCurve port.  The static SetPoles/GetPoles/CheckCurveData/
//! KnotAnalysis/Rational helpers are ported at module level.

use rcad_kernel::math::bspl_lib as bspl;
use rcad_kernel::math::bspl_lib::{at, ati, GeomAbsKnotDistribution};

/// OCCT GeomAbs_Shape values reachable from `updateKnots` / `IsCN`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmoothShape {
    C0,
    G1,
    C1,
    G2,
    C2,
    C3,
    CN,
}

/// OCCT Epsilon(V) — the distance from V to the next representable double.
fn epsilon_of(value: f64) -> f64 {
    if value >= 0.0 {
        value.next_up() - value
    } else {
        value - value.next_down()
    }
}

/// OCCT static SetPoles (Law_BSpline.cxx L38-50): flatten (poles, weights)
/// to the interleaved homogeneous representation (P*w, w) — dimension 2.
fn set_poles(poles: &[f64], weights: &[f64], fp: &mut [f64]) {
    let mut j = 0usize;
    for i in 0..poles.len() {
        let w = weights[i];
        fp[j] = poles[i] * w;
        j += 1;
        fp[j] = w;
        j += 1;
    }
}

/// OCCT static GetPoles (Law_BSpline.cxx L52-65).
fn get_poles(fp: &[f64], poles: &mut [f64], weights: &mut [f64]) {
    let mut j = 0usize;
    for i in 0..poles.len() {
        let w = fp[j + 1];
        weights[i] = w;
        poles[i] = fp[j] / w;
        j += 2;
    }
}

/// OCCT static CheckCurveData (Law_BSpline.cxx L67-96).
fn check_curve_data(
    cpoles: &[f64],
    cknots: &[f64],
    cmults: &[i32],
    degree: usize,
    periodic: bool,
) {
    if degree < 1 || degree > LAW_BSPLINE_MAX_DEGREE {
        panic!("Standard_ConstructionError");
    }
    if cpoles.len() < 2 {
        panic!("Standard_ConstructionError");
    }
    if cknots.len() != cmults.len() {
        panic!("Standard_ConstructionError");
    }
    for i in 1..cknots.len() {
        if cknots[i] - cknots[i - 1] <= epsilon_of(cknots[i - 1].abs()) {
            panic!("Standard_ConstructionError");
        }
    }
    if cpoles.len() != bspl::nb_poles(degree, periodic, cmults) {
        panic!("Standard_ConstructionError");
    }
}

/// OCCT static Rational (Law_BSpline.cxx L156-171) — check the rationality
/// of an array of weights.
fn rational(w: &[f64]) -> bool {
    let n = w.len();
    let mut rat = false;
    for i in 0..n - 1 {
        rat = (w[i] - w[i + 1]).abs() > 2.2250738585072014e-308; // gp::Resolution()
        if rat {
            break;
        }
    }
    rat
}

/// OCCT Law_BSpline::MaxDegree() == BSplCLib::MaxDegree() == 25.
pub const LAW_BSPLINE_MAX_DEGREE: usize = 25;

/// OCCT Law_BSpline — definition of the 1D B_spline curve.
#[derive(Debug, Clone)]
pub struct LawBSpline {
    rational: bool,
    periodic: bool,
    knot_set: GeomAbsKnotDistribution,
    smooth: SmoothShape,
    deg: usize,
    poles: Vec<f64>,
    weights: Option<Vec<f64>>,
    flatknots: Vec<f64>,
    knots: Vec<f64>,
    mults: Vec<i32>,
}

impl LawBSpline {
    /// OCCT Law_BSpline(Poles, Knots, Mults, Degree, Periodic) (L199-223) —
    /// the non-rational constructor.
    pub fn new(
        poles: &[f64],
        knots: &[f64],
        mults: &[i32],
        degree: usize,
        periodic: bool,
    ) -> Self {
        check_curve_data(poles, knots, mults, degree, periodic);
        let mut curve = LawBSpline {
            rational: false,
            periodic,
            knot_set: GeomAbsKnotDistribution::NonUniform,
            smooth: SmoothShape::CN,
            deg: degree,
            poles: poles.to_vec(),
            weights: None,
            flatknots: Vec::new(),
            knots: knots.to_vec(),
            mults: mults.to_vec(),
        };
        curve.update_knots();
        curve
    }

    /// OCCT Law_BSpline(Poles, Weights, Knots, Mults, Degree, Periodic)
    /// (L225-268) — the rational constructor.
    pub fn new_rational(
        poles: &[f64],
        weights: &[f64],
        knots: &[f64],
        mults: &[i32],
        degree: usize,
        periodic: bool,
    ) -> Self {
        check_curve_data(poles, knots, mults, degree, periodic);
        if weights.len() != poles.len() {
            panic!("Standard_ConstructionError: Law_BSpline");
        }
        for i in 0..weights.len() {
            if weights[i] <= 2.2250738585072014e-308 {
                panic!("Standard_ConstructionError: Law_BSpline");
            }
        }
        // check really rational
        let rational = rational(weights);
        let mut curve = LawBSpline {
            rational,
            periodic,
            knot_set: GeomAbsKnotDistribution::NonUniform,
            smooth: SmoothShape::CN,
            deg: degree,
            poles: poles.to_vec(),
            weights: if rational {
                Some(weights.to_vec())
            } else {
                None
            },
            flatknots: Vec::new(),
            knots: knots.to_vec(),
            mults: mults.to_vec(),
        };
        curve.update_knots();
        curve
    }

    /// OCCT Copy (L173-196).
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// OCCT MaxDegree (L271-274).
    pub fn max_degree() -> usize {
        LAW_BSPLINE_MAX_DEGREE
    }

    /// OCCT IncreaseDegree (L276-343).
    pub fn increase_degree(&mut self, degree: usize) {
        if degree == self.deg {
            return;
        }
        if degree < self.deg || degree > Self::max_degree() {
            panic!("Standard_ConstructionError");
        }
        let from_k1 = self.first_uknot_index() as i32;
        let to_k2 = self.last_uknot_index() as i32;
        let step = (degree - self.deg) as i32;
        let mut npoles = vec![0.0f64; self.poles.len() + (step * (to_k2 - from_k1)) as usize];
        let nbknots = bspl::increase_degree_count_knots(self.deg, degree, self.periodic, &self.mults);
        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];
        let mut nweights: Option<Vec<f64>> = None;
        if self.is_rational() {
            let mut nw = vec![0.0f64; npoles.len()];
            let mut adimpol = vec![0.0f64; 2 * self.poles.len()];
            set_poles(&self.poles, self.weights.as_ref().unwrap(), &mut adimpol);
            let mut adimnpol = vec![0.0f64; 2 * npoles.len()];
            bspl::increase_degree(
                self.deg,
                degree,
                self.periodic,
                2,
                &adimpol,
                &self.knots,
                &self.mults,
                &mut adimnpol,
                &mut nknots,
                &mut nmults,
            );
            get_poles(&adimnpol, &mut npoles, &mut nw);
            nweights = Some(nw);
        } else {
            bspl::increase_degree(
                self.deg,
                degree,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                &mut npoles,
                &mut nknots,
                &mut nmults,
            );
        }
        self.deg = degree;
        self.poles = npoles;
        self.weights = nweights;
        self.knots = nknots;
        self.mults = nmults;
        self.update_knots();
    }

    /// OCCT IncreaseMultiplicity(Index, M) (L345-353).
    pub fn increase_multiplicity(&mut self, index: usize, m: i32) {
        let k = [at(&self.knots, index as i32)];
        let m_arr = [m - self.mults[index - 1]];
        self.insert_knots(&k, &m_arr, epsilon_of(1.0), true);
    }

    /// OCCT IncreaseMultiplicity(I1, I2, M) (L355-368).
    pub fn increase_multiplicity_range(&mut self, i1: usize, i2: usize, m: i32) {
        let mut k = Vec::with_capacity(i2 - i1 + 1);
        let mut m_arr = Vec::with_capacity(i2 - i1 + 1);
        for i in i1..=i2 {
            k.push(at(&self.knots, i as i32));
            m_arr.push(m - self.mults[i - 1]);
        }
        self.insert_knots(&k, &m_arr, epsilon_of(1.0), true);
    }

    /// OCCT IncrementMultiplicity(I1, I2, Step) (L370-380).
    pub fn increment_multiplicity_range(&mut self, i1: usize, i2: usize, step: i32) {
        let mut k = Vec::with_capacity(i2 - i1 + 1);
        let m_arr = vec![step; i2 - i1 + 1];
        for i in i1..=i2 {
            k.push(at(&self.knots, i as i32));
        }
        self.insert_knots(&k, &m_arr, epsilon_of(1.0), true);
    }

    /// OCCT InsertKnot(U, M, ParametricTolerance, Add) (L382-392).
    pub fn insert_knot(&mut self, u: f64, m: i32, parametric_tolerance: f64, add: bool) {
        let k = [u];
        let m_arr = [m];
        self.insert_knots(&k, &m_arr, parametric_tolerance, add);
    }

    /// OCCT InsertKnots(Knots, Mults, Epsilon, Add) (L394-469).
    pub fn insert_knots(
        &mut self,
        knots: &[f64],
        mults: &[i32],
        epsilon: f64,
        add: bool,
    ) {
        // Check and compute new sizes
        let mut nbpoles = 0i32;
        let mut nbknots = 0i32;
        if !bspl::prepare_insert_knots(
            self.deg,
            self.periodic,
            &self.knots,
            &self.mults,
            knots,
            Some(mults),
            &mut nbpoles,
            &mut nbknots,
            epsilon,
            add,
        ) {
            panic!("Standard_ConstructionError: Law_BSpline::InsertKnots");
        }
        if nbpoles as usize == self.poles.len() {
            return;
        }
        let mut npoles = vec![0.0f64; nbpoles as usize];
        let mut nknots = self.knots.clone();
        let mut nmults = self.mults.clone();
        if nbknots as usize != self.knots.len() {
            nknots = vec![0.0f64; nbknots as usize];
            nmults = vec![0i32; nbknots as usize];
        }
        if self.rational {
            let mut nweights = vec![0.0f64; nbpoles as usize];
            let mut adimpol = vec![0.0f64; 2 * self.poles.len()];
            set_poles(&self.poles, self.weights.as_ref().unwrap(), &mut adimpol);
            let mut adimnpol = vec![0.0f64; 2 * npoles.len()];
            bspl::insert_knots(
                self.deg,
                self.periodic,
                2,
                &adimpol,
                &self.knots,
                &self.mults,
                knots,
                Some(mults),
                &mut adimnpol,
                &mut nknots,
                &mut nmults,
                epsilon,
                add,
            );
            get_poles(&adimnpol, &mut npoles, &mut nweights);
            self.weights = Some(nweights);
        } else {
            bspl::insert_knots(
                self.deg,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                knots,
                Some(mults),
                &mut npoles,
                &mut nknots,
                &mut nmults,
                epsilon,
                add,
            );
        }
        self.poles = npoles;
        self.knots = nknots;
        self.mults = nmults;
        self.update_knots();
    }

    /// OCCT RemoveKnot(Index, M, Tolerance) (L471-548).
    pub fn remove_knot(&mut self, index: usize, m: i32, tolerance: f64) -> bool {
        if m < 0 {
            return true;
        }
        let i1 = self.first_uknot_index();
        let i2 = self.last_uknot_index();
        if !self.periodic && (index <= i1 || index >= i2) {
            panic!("Standard_OutOfRange");
        } else if self.periodic && (index < i1 || index > i2) {
            panic!("Standard_OutOfRange");
        }
        let step = self.mults[index - 1] - m;
        if step <= 0 {
            return true;
        }
        let mut npoles = vec![0.0f64; self.poles.len() - step as usize];
        let mut nknots = self.knots.clone();
        let mut nmults = self.mults.clone();
        if m == 0 {
            nknots = vec![0.0f64; self.knots.len() - 1];
            nmults = vec![0i32; self.knots.len() - 1];
        }
        if self.is_rational() {
            let mut nweights = vec![0.0f64; npoles.len()];
            let mut adimpol = vec![0.0f64; 2 * self.poles.len()];
            set_poles(&self.poles, self.weights.as_ref().unwrap(), &mut adimpol);
            let mut adimnpol = vec![0.0f64; 2 * npoles.len()];
            if !bspl::remove_knot(
                index,
                m,
                self.deg,
                self.periodic,
                2,
                &adimpol,
                &self.knots,
                &self.mults,
                &mut adimnpol,
                &mut nknots,
                &mut nmults,
                tolerance,
            ) {
                return false;
            }
            get_poles(&adimnpol, &mut npoles, &mut nweights);
            self.weights = Some(nweights);
        } else {
            if !bspl::remove_knot(
                index,
                m,
                self.deg,
                self.periodic,
                1,
                &self.poles,
                &self.knots,
                &self.mults,
                &mut npoles,
                &mut nknots,
                &mut nmults,
                tolerance,
            ) {
                return false;
            }
        }
        self.poles = npoles;
        self.knots = nknots;
        self.mults = nmults;
        self.update_knots();
        true
    }

    /// OCCT Reverse (L556-577).
    pub fn reverse(&mut self) {
        self.knots.reverse();
        self.mults.reverse();
        let last = if self.periodic {
            self.flatknots.len() - self.deg - 1
        } else {
            self.poles.len()
        };
        self.poles[..last].reverse();
        if self.rational {
            let weights = self.weights.as_mut().unwrap();
            weights[..last].reverse();
        }
        self.update_knots();
    }

    /// OCCT ReversedParameter (L579-582).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.first_parameter() + self.last_parameter() - u
    }

    /// OCCT Segment(U1, U2) (L584-711).
    pub fn segment(&mut self, u1: f64, u2: f64) {
        assert!(u2 >= u1, "Standard_DomainError: Law_BSpline::Segment");
        let eps = epsilon_of(u1.abs().max(u2.abs()));
        let delta = u2 - u1;
        let mut new_u1 = 0.0;
        let mut new_u2 = 0.0;
        let mut u = 0.0;
        {
            let mut index = 0i32;
            bspl::locate_parameter_knots_mults(
    self.deg,
    &self.knots.clone(),
    &self.mults.clone(),
    u1,
                self.periodic,
                1,
                self.knots.len() as i32,
                &mut index,
                &mut new_u1,
            );
            let mut index = 0i32;
            bspl::locate_parameter_knots_mults(
    self.deg,
    &self.knots.clone(),
    &self.mults.clone(),
    u2,
                self.periodic,
                1,
                self.knots.len() as i32,
                &mut index,
                &mut new_u2,
            );
        }
        let knots_seg = [new_u1.min(new_u2), new_u1.max(new_u2)];
        let mults_seg = [self.deg as i32, self.deg as i32];
        self.insert_knots(&knots_seg, &mults_seg, eps, true);
        if self.periodic {
            // set the origine at NewU1
            let mut index0 = 0i32;
            bspl::locate_parameter_knots_mults(
    self.deg,
    &self.knots.clone(),
    &self.mults.clone(),
    u1,
                self.periodic,
                1,
                self.knots.len() as i32,
                &mut index0,
                &mut u,
            );
            if (at(&self.knots, index0 + 1) - u).abs() < eps {
                index0 += 1;
            }
            self.set_origin(index0 as usize);
            self.set_not_periodic();
        }
        // compute index1 and index2 to set the new knots and mults
        let from_u1 = 1i32;
        let to_u2 = self.knots.len() as i32;
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        bspl::locate_parameter_knots_mults(
    self.deg,
    &self.knots.clone(),
    &self.mults.clone(),
    new_u1,
            self.periodic,
            from_u1,
            to_u2,
            &mut index1,
            &mut u,
        );
        bspl::locate_parameter_knots_mults(
    self.deg,
    &self.knots.clone(),
    &self.mults.clone(),
    new_u1 + delta,
            self.periodic,
            from_u1,
            to_u2,
            &mut index2,
            &mut u,
        );
        if (at(&self.knots, index2 + 1) - u).abs() < eps {
            index2 += 1;
        }
        let nbknots = (index2 - index1 + 1) as usize;
        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];
        let mut k = 1i32;
        for i in index1..=index2 {
            nknots[(k - 1) as usize] = at(&self.knots, i);
            nmults[(k - 1) as usize] = ati(&self.mults, i);
            k += 1;
        }
        nmults[0] = self.deg as i32 + 1;
        nmults[nbknots - 1] = self.deg as i32 + 1;
        // compute index1 and index2 to set the new poles and weights
        let mut pindex1 = bspl::pole_index(self.deg, index1, self.periodic, &self.mults);
        let mut pindex2 = bspl::pole_index(self.deg, index2, self.periodic, &self.mults);
        pindex1 += 1;
        pindex2 = (pindex2 + 1).min(self.poles.len() as i32);
        let nbpoles = (pindex2 - pindex1 + 1) as usize;
        let mut npoles = vec![0.0f64; nbpoles];
        let mut nweights = vec![0.0f64; nbpoles];
        let mut k = 1i32;
        if self.rational {
            for i in pindex1..=pindex2 {
                npoles[(k - 1) as usize] = at(&self.poles, i);
                nweights[(k - 1) as usize] = at(self.weights.as_ref().unwrap(), i);
                k += 1;
            }
        } else {
            for i in pindex1..=pindex2 {
                npoles[(k - 1) as usize] = at(&self.poles, i);
                k += 1;
            }
        }
        self.knots = nknots;
        self.mults = nmults;
        self.poles = npoles;
        if self.rational {
            self.weights = Some(nweights);
        }
        self.update_knots();
    }

    /// OCCT SetKnot(Index, K) (L713-741).
    pub fn set_knot(&mut self, index: usize, k: f64) {
        if index < 1 || index > self.knots.len() {
            panic!("Standard_OutOfRange");
        }
        let dk = epsilon_of(k).abs();
        if index == 1 {
            if k >= self.knots[1] - dk {
                panic!("Standard_ConstructionError");
            }
        } else if index == self.knots.len() {
            if k <= self.knots[self.knots.len() - 2] + dk {
                panic!("Standard_ConstructionError");
            }
        } else if k <= self.knots[index - 2] + dk || k >= self.knots[index] - dk {
            panic!("Standard_ConstructionError");
        }
        if k != self.knots[index - 1] {
            self.knots[index - 1] = k;
            self.update_knots();
        }
    }

    /// OCCT SetKnots(K) (L743-748).
    pub fn set_knots(&mut self, k: &[f64]) {
        check_curve_data(&self.poles, k, &self.mults, self.deg, self.periodic);
        self.knots = k.to_vec();
        self.update_knots();
    }

    /// OCCT SetKnot(Index, K, M) (L750-754).
    pub fn set_knot_with_mult(&mut self, index: usize, k: f64, m: i32) {
        self.increase_multiplicity(index, m);
        self.set_knot(index, k);
    }

    /// OCCT SetPeriodic (L756-788).
    pub fn set_periodic(&mut self) {
        let first = self.first_uknot_index();
        let last = self.last_uknot_index();
        let cknots: Vec<f64> = self.knots[first - 1..last].to_vec();
        self.knots = cknots;
        let mut cmults: Vec<i32> = self.mults[first - 1..last].to_vec();
        let len = cmults.len();
        cmults[0] = cmults[0].max(cmults[len - 1]);
        cmults[len - 1] = cmults[0];
        self.mults = cmults.clone();
        // compute new number of poles
        let nbp = bspl::nb_poles(self.deg, true, &cmults);
        self.poles = self.poles[..nbp].to_vec();
        if self.rational {
            let weights = self.weights.as_ref().unwrap();
            self.weights = Some(weights[..nbp].to_vec());
        }
        self.periodic = true;
        self.update_knots();
    }

    /// OCCT SetOrigin(Index) (L790-861).
    pub fn set_origin(&mut self, index: usize) {
        assert!(self.periodic, "Standard_NoSuchObject: Law_BSpline::SetOrigin");
        let first = self.first_uknot_index();
        let last = self.last_uknot_index();
        assert!(
            index >= first && index <= last,
            "Standard_DomainError: Law_BSpline::SetOrigine"
        );
        let nbknots = self.knots.len();
        let nbpoles = self.poles.len();
        let mut nknots = vec![0.0f64; nbknots];
        let mut nmults = vec![0i32; nbknots];
        // set the knots and mults
        let period = self.knots[last - 1] - self.knots[first - 1];
        let mut k = 1usize;
        for i in index..=last {
            nknots[k - 1] = self.knots[i - 1];
            nmults[k - 1] = self.mults[i - 1];
            k += 1;
        }
        for i in (first + 1)..=index {
            nknots[k - 1] = self.knots[i - 1] + period;
            nmults[k - 1] = self.mults[i - 1];
            k += 1;
        }
        let mut pole_index = 1usize;
        for i in (first + 1)..=index {
            pole_index += self.mults[i - 1] as usize;
        }
        // set the poles and weights
        let mut npoles = vec![0.0f64; nbpoles];
        let mut nweights = vec![0.0f64; nbpoles];
        let first_pole = 1usize;
        let last_pole = self.poles.len();
        if self.rational {
            let weights = self.weights.as_ref().unwrap();
            let mut k = 1usize;
            for i in pole_index..=last_pole {
                npoles[k - 1] = self.poles[i - 1];
                nweights[k - 1] = weights[i - 1];
                k += 1;
            }
            for i in first_pole..pole_index {
                npoles[k - 1] = self.poles[i - 1];
                nweights[k - 1] = weights[i - 1];
                k += 1;
            }
        } else {
            let mut k = 1usize;
            for i in pole_index..=last_pole {
                npoles[k - 1] = self.poles[i - 1];
                k += 1;
            }
            for i in first_pole..pole_index {
                npoles[k - 1] = self.poles[i - 1];
                k += 1;
            }
        }
        self.poles = npoles;
        self.knots = nknots;
        self.mults = nmults;
        if self.rational {
            self.weights = Some(nweights);
        }
        self.update_knots();
    }

    /// OCCT SetNotPeriodic (L863-911).
    pub fn set_not_periodic(&mut self) {
        if self.periodic {
            let mut nb_knots = 0i32;
            let mut nb_poles = 0i32;
            bspl::prepare_unperiodize(self.deg, &self.mults, &mut nb_knots, &mut nb_poles);
            let mut npoles = vec![0.0f64; nb_poles as usize];
            let mut nknots = vec![0.0f64; nb_knots as usize];
            let mut nmults = vec![0i32; nb_knots as usize];
            let mut nweights: Option<Vec<f64>> = None;
            if self.is_rational() {
                let mut nw = vec![0.0f64; nb_poles as usize];
                let mut adimpol = vec![0.0f64; 2 * self.poles.len()];
                set_poles(&self.poles, self.weights.as_ref().unwrap(), &mut adimpol);
                let mut adimnpol = vec![0.0f64; 2 * npoles.len()];
                bspl::unperiodize(
                    self.deg,
                    &self.mults,
                    &self.knots,
                    &adimpol,
                    &mut nmults,
                    &mut nknots,
                    &mut adimnpol,
                );
                get_poles(&adimnpol, &mut npoles, &mut nw);
                nweights = Some(nw);
            } else {
                bspl::unperiodize(
                    self.deg,
                    &self.mults,
                    &self.knots,
                    &self.poles,
                    &mut nmults,
                    &mut nknots,
                    &mut npoles,
                );
            }
            self.poles = npoles;
            self.weights = nweights;
            self.mults = nmults;
            self.knots = nknots;
            self.periodic = false;
            self.update_knots();
        }
    }

    /// OCCT SetPole(Index, P) (L913-921).
    pub fn set_pole(&mut self, index: usize, p: f64) {
        if index < 1 || index > self.poles.len() {
            panic!("Standard_OutOfRange");
        }
        self.poles[index - 1] = p;
    }

    /// OCCT SetPole(Index, P, W) (L923-927).
    pub fn set_pole_with_weight(&mut self, index: usize, p: f64, w: f64) {
        self.set_pole(index, p);
        self.set_weight(index, w);
    }

    /// OCCT SetWeight(Index, W) (L929-954).
    pub fn set_weight(&mut self, index: usize, w: f64) {
        if index < 1 || index > self.poles.len() {
            panic!("Standard_OutOfRange");
        }
        if w <= 2.2250738585072014e-308 {
            panic!("Standard_ConstructionError");
        }
        let mut rat = self.is_rational() || (w - 1.0).abs() > 2.2250738585072014e-308;
        if rat {
            if !self.is_rational() {
                self.weights = Some(vec![1.0; self.poles.len()]);
            }
            self.weights.as_mut().unwrap()[index - 1] = w;
            if self.is_rational() {
                rat = rational(self.weights.as_ref().unwrap());
                if !rat {
                    self.weights = None;
                }
            }
            self.rational = self.weights.is_some();
        }
    }

    /// OCCT UpdateKnots (L956-1000).
    fn update_knots(&mut self) {
        self.rational = self.weights.is_some();
        let mut max_knot_mult = 0i32;
        let mut knot_set = GeomAbsKnotDistribution::NonUniform;
        bspl::knot_analysis(
            self.deg,
            self.periodic,
            &self.knots,
            &self.mults,
            &mut knot_set,
            &mut max_knot_mult,
        );
        self.knot_set = knot_set;
        if knot_set == GeomAbsKnotDistribution::Uniform && !self.periodic {
            self.flatknots = self.knots.clone();
        } else {
            let len = bspl::knot_sequence_length(&self.mults, self.deg, self.periodic);
            let mut flat = vec![0.0f64; len];
            bspl::knot_sequence(
                &self.knots,
                &self.mults,
                self.deg,
                self.periodic,
                &mut flat,
            );
            self.flatknots = flat;
        }
        self.smooth = if max_knot_mult == 0 {
            SmoothShape::CN
        } else {
            match self.deg as i32 - max_knot_mult {
                0 => SmoothShape::C0,
                1 => SmoothShape::C1,
                2 => SmoothShape::C2,
                3 => SmoothShape::C3,
                _ => SmoothShape::C3,
            }
        };
    }

    /// OCCT PeriodicNormalization (L1005-1020).
    pub fn periodic_normalization(&self, parameter: &mut f64) {
        if self.periodic {
            let period = self.flatknots[self.flatknots.len() - self.deg - 1]
                - self.flatknots[self.deg];
            while *parameter > self.flatknots[self.flatknots.len() - self.deg - 1] {
                *parameter -= period;
            }
            while *parameter < self.flatknots[self.deg] {
                *parameter += period;
            }
        }
    }

    /// OCCT IsCN(N) (L1022-1050).
    pub fn is_cn(&self, n: i32) -> bool {
        assert!(n >= 0, "Standard_RangeError");
        match self.smooth {
            SmoothShape::CN => true,
            SmoothShape::C0 => n <= 0,
            SmoothShape::G1 => n <= 0,
            SmoothShape::C1 => n <= 1,
            SmoothShape::G2 => n <= 1,
            SmoothShape::C2 => n <= 2,
            SmoothShape::C3 => {
                if n <= 3 {
                    true
                } else {
                    n <= self.deg as i32
                        - bspl::max_knot_mult(&self.mults, 2, self.mults.len() as i32 - 1)
                }
            }
        }
    }

    /// OCCT IsClosed (L1052-1055).
    pub fn is_closed(&self) -> bool {
        (self.start_point() - self.end_point()).abs() <= 2.2250738585072014e-308
    }

    /// OCCT IsPeriodic (L1057-1060).
    pub fn is_periodic(&self) -> bool {
        self.periodic
    }

    /// OCCT Continuity (L1062-1065).
    pub fn continuity(&self) -> SmoothShape {
        self.smooth
    }

    /// OCCT Degree (L1067-1070).
    pub fn degree(&self) -> usize {
        self.deg
    }

    /// OCCT Value (L1072-1077).
    pub fn value(&self, u: f64) -> f64 {
        let mut p = 0.0;
        self.d0(u, &mut p);
        p
    }

    /// OCCT D0 (L1079-1092).
    pub fn d0(&self, u: f64, p: &mut f64) {
        let mut new_u = u;
        self.periodic_normalization(&mut new_u);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let mut results = vec![0.0f64; 2];
                bspl::eval_flat(
                    new_u, self.periodic, 0, &mut extrap, self.deg, &self.flatknots, 2, &flat,
                    &mut results,
                );
                *p = results[0] / results[1];
            }
            None => {
                let mut results = vec![0.0f64; 1];
                bspl::eval_flat(
                    new_u, self.periodic, 0, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                *p = results[0];
            }
        }
    }

    /// OCCT D1 (L1094-1110).
    pub fn d1(&self, u: f64, p: &mut f64, v1: &mut f64) {
        let mut new_u = u;
        self.periodic_normalization(&mut new_u);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let mut results = vec![0.0f64; 4];
                bspl::eval_flat(
                    new_u, self.periodic, 1, &mut extrap, self.deg, &self.flatknots, 2, &flat,
                    &mut results,
                );
                let f = results[0] / results[1];
                let d = (results[2] - f * results[3]) / results[1];
                *p = f;
                *v1 = d;
            }
            None => {
                let mut results = vec![0.0f64; 2];
                bspl::eval_flat(
                    new_u, self.periodic, 1, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                *p = results[0];
                *v1 = results[1];
            }
        }
    }

    /// OCCT D2 (L1112-1128).
    pub fn d2(&self, u: f64, p: &mut f64, v1: &mut f64, v2: &mut f64) {
        let mut new_u = u;
        self.periodic_normalization(&mut new_u);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let mut results = vec![0.0f64; 6];
                bspl::eval_flat(
                    new_u, self.periodic, 2, &mut extrap, self.deg, &self.flatknots, 2, &flat,
                    &mut results,
                );
                // PLib::RationalDerivatives (dim 1): D0 = P0/W0,
                // D1 = (P1 - D0*W1)/W0, D2 = (P2 - D0*W2 - 2*D1*W1)/W0.
                let d0 = results[0] / results[1];
                let d1 = (results[2] - d0 * results[3]) / results[1];
                let d2 = (results[4] - d0 * results[5] - 2.0 * d1 * results[3]) / results[1];
                *p = d0;
                *v1 = d1;
                *v2 = d2;
            }
            None => {
                let mut results = vec![0.0f64; 3];
                bspl::eval_flat(
                    new_u, self.periodic, 2, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                *p = results[0];
                *v1 = results[1];
                *v2 = results[2];
            }
        }
    }

    /// OCCT D3 (L1130-1156).
    pub fn d3(&self, u: f64, p: &mut f64, v1: &mut f64, v2: &mut f64, v3: &mut f64) {
        let mut new_u = u;
        self.periodic_normalization(&mut new_u);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let mut results = vec![0.0f64; 8];
                bspl::eval_flat(
                    new_u, self.periodic, 3, &mut extrap, self.deg, &self.flatknots, 2, &flat,
                    &mut results,
                );
                // PLib::RationalDerivatives (dim 1), extended to D3.
                let d0 = results[0] / results[1];
                let d1 = (results[2] - d0 * results[3]) / results[1];
                let d2 =
                    (results[4] - d0 * results[5] - 2.0 * d1 * results[3]) / results[1];
                let d3 = (results[6] - d0 * results[7] - 3.0 * d1 * results[5]
                    - 3.0 * d2 * results[3])
                    / results[1];
                *p = d0;
                *v1 = d1;
                *v2 = d2;
                *v3 = d3;
            }
            None => {
                let mut results = vec![0.0f64; 4];
                bspl::eval_flat(
                    new_u, self.periodic, 3, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                *p = results[0];
                *v1 = results[1];
                *v2 = results[2];
                *v3 = results[3];
            }
        }
    }

    /// OCCT DN (L1158-1170).
    pub fn dn(&self, u: f64, n: usize) -> f64 {
        let mut new_u = u;
        self.periodic_normalization(&mut new_u);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let count = n + 1;
                let mut results = vec![0.0f64; count * 2];
                bspl::eval_flat(
                    new_u, self.periodic, n as i32, &mut extrap, self.deg, &self.flatknots, 2,
                    &flat, &mut results,
                );
                // PLib::RationalDerivatives general form (dim 1).
                let mut ders = vec![0.0f64; count];
                for i in 0..count {
                    ders[i] = results[i * 2];
                    let mut binom = 1.0f64;
                    for j in 1..=i {
                        binom *= (i - j + 1) as f64 / j as f64;
                        ders[i] -= binom * ders[i - j] * results[j * 2 + 1];
                    }
                    ders[i] /= results[1];
                }
                ders[n]
            }
            None => {
                let mut results = vec![0.0f64; n + 1];
                bspl::eval_flat(
                    new_u, self.periodic, n as i32, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                results[n]
            }
        }
    }

    /// OCCT EndPoint (L1172-1181).
    pub fn end_point(&self) -> f64 {
        if *self.mults.last().unwrap() == self.deg as i32 + 1 {
            *self.poles.last().unwrap()
        } else {
            self.value(self.last_parameter())
        }
    }

    /// OCCT FirstUKnotIndex (L1183-1192).
    pub fn first_uknot_index(&self) -> usize {
        if self.periodic {
            1
        } else {
            bspl::first_uknot_index_mults(self.deg, &self.mults) as usize
        }
    }

    /// OCCT FirstParameter (L1194-1197).
    pub fn first_parameter(&self) -> f64 {
        self.flatknots[self.deg]
    }

    /// OCCT Knot(Index) (L1199-1203).
    pub fn knot(&self, index: usize) -> f64 {
        assert!(
            index >= 1 && index <= self.knots.len(),
            "Standard_OutOfRange: Law_BSpline::Knot"
        );
        self.knots[index - 1]
    }

    /// OCCT KnotDistribution (L1205-1208).
    pub fn knot_distribution(&self) -> GeomAbsKnotDistribution {
        self.knot_set
    }

    /// OCCT Knots(K) (L1210-1213) — returns the knot array.
    pub fn knots(&self) -> &Vec<f64> {
        &self.knots
    }

    /// OCCT KnotSequence(K) (L1215-1218) — returns the flat knot sequence.
    pub fn knot_sequence(&self) -> &Vec<f64> {
        &self.flatknots
    }

    /// OCCT LastUKnotIndex (L1220-1229).
    pub fn last_uknot_index(&self) -> usize {
        if self.periodic {
            self.knots.len()
        } else {
            bspl::last_uknot_index_mults(self.deg, &self.mults) as usize
        }
    }

    /// OCCT LastParameter (L1231-1234).
    pub fn last_parameter(&self) -> f64 {
        self.flatknots[self.flatknots.len() - self.deg - 1]
    }

    /// OCCT LocalValue (L1236-1242).
    pub fn local_value(&self, u: f64, from_k1: i32, to_k2: i32) -> f64 {
        let mut p = 0.0;
        self.local_d0(u, from_k1, to_k2, &mut p);
        p
    }

    /// OCCT LocalD0 (L1244-1262).
    pub fn local_d0(&self, u: f64, from_k1: i32, to_k2: i32, p: &mut f64) {
        assert!(from_k1 != to_k2, "Standard_DomainError: Law_BSpline::LocalValue");
        let mut uu = u;
        let mut index = 0i32;
        bspl::locate_parameter_flat(
            self.deg,
            &self.flatknots.clone(),
            u,
            self.periodic,
            from_k1,
            to_k2,
            &mut index,
            &mut uu,
        );
        // OCCT feeds this span index to BSplCLib::D0 as a hint; the rcad
        // eval_flat locates the span itself, so the value is unused.
        let _index = bspl::flat_index(self.deg, index as usize, &self.mults, self.periodic);
        let mut extrap = [0i32; 2];
        match &self.weights {
            Some(w) => {
                let mut flat = Vec::with_capacity(self.poles.len() * 2);
                for (po, wi) in self.poles.iter().zip(w.iter()) {
                    flat.push(po * wi);
                    flat.push(*wi);
                }
                let mut results = vec![0.0f64; 2];
                bspl::eval_flat(
                    uu, self.periodic, 0, &mut extrap, self.deg, &self.flatknots, 2, &flat,
                    &mut results,
                );
                *p = results[0] / results[1];
            }
            None => {
                let mut results = vec![0.0f64; 1];
                bspl::eval_flat(
                    uu, self.periodic, 0, &mut extrap, self.deg, &self.flatknots, 1,
                    &self.poles, &mut results,
                );
                *p = results[0];
            }
        }
    }

    /// OCCT LocalD1 (L1264-1286).
    pub fn local_d1(&self, u: f64, from_k1: i32, to_k2: i32, p: &mut f64, v1: &mut f64) {
        assert!(from_k1 != to_k2, "Standard_DomainError: Law_BSpline::LocalD1");
        let mut uu = u;
        let mut index = 0i32;
        bspl::locate_parameter_flat(
            self.deg,
            &self.flatknots.clone(),
            u,
            self.periodic,
            from_k1,
            to_k2,
            &mut index,
            &mut uu,
        );
        let _index = bspl::flat_index(self.deg, index as usize, &self.mults, self.periodic);
        self.d1(uu, p, v1);
    }

    /// OCCT LocalD2 (L1288-1312).
    pub fn local_d2(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
        p: &mut f64,
        v1: &mut f64,
        v2: &mut f64,
    ) {
        assert!(from_k1 != to_k2, "Standard_DomainError: Law_BSpline::LocalD2");
        let mut uu = u;
        let mut index = 0i32;
        bspl::locate_parameter_flat(
            self.deg,
            &self.flatknots.clone(),
            u,
            self.periodic,
            from_k1,
            to_k2,
            &mut index,
            &mut uu,
        );
        let _index = bspl::flat_index(self.deg, index as usize, &self.mults, self.periodic);
        self.d2(uu, p, v1, v2);
    }

    /// OCCT LocalD3 (L1314-1340).
    pub fn local_d3(
        &self,
        u: f64,
        from_k1: i32,
        to_k2: i32,
        p: &mut f64,
        v1: &mut f64,
        v2: &mut f64,
        v3: &mut f64,
    ) {
        assert!(from_k1 != to_k2, "Standard_DomainError: Law_BSpline::LocalD3");
        let mut uu = u;
        let mut index = 0i32;
        bspl::locate_parameter_flat(
            self.deg,
            &self.flatknots.clone(),
            u,
            self.periodic,
            from_k1,
            to_k2,
            &mut index,
            &mut uu,
        );
        let _index = bspl::flat_index(self.deg, index as usize, &self.mults, self.periodic);
        self.d3(uu, p, v1, v2, v3);
    }

    /// OCCT LocalDN (L1342-1364).
    pub fn local_dn(&self, u: f64, from_k1: i32, to_k2: i32, n: usize) -> f64 {
        assert!(from_k1 != to_k2, "Standard_DomainError: Law_BSpline::LocalD3");
        let mut uu = u;
        let mut index = 0i32;
        bspl::locate_parameter_flat(
            self.deg,
            &self.flatknots.clone(),
            u,
            self.periodic,
            from_k1,
            to_k2,
            &mut index,
            &mut uu,
        );
        let _index = bspl::flat_index(self.deg, index as usize, &self.mults, self.periodic);
        self.dn(uu, n)
    }

    /// OCCT Multiplicity(Index) (L1366-1370).
    pub fn multiplicity(&self, index: usize) -> i32 {
        assert!(
            index >= 1 && index <= self.mults.len(),
            "Standard_OutOfRange: Law_BSpline::Multiplicity"
        );
        self.mults[index - 1]
    }

    /// OCCT Multiplicities(M) (L1372-1375).
    pub fn multiplicities(&self) -> &Vec<i32> {
        &self.mults
    }

    /// OCCT NbKnots (L1377-1380).
    pub fn nb_knots(&self) -> usize {
        self.knots.len()
    }

    /// OCCT NbPoles (L1382-1385).
    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }

    /// OCCT Pole(Index) (L1387-1391).
    pub fn pole(&self, index: usize) -> f64 {
        assert!(
            index >= 1 && index <= self.poles.len(),
            "Standard_OutOfRange: Law_BSpline::Pole"
        );
        self.poles[index - 1]
    }

    /// OCCT Poles(P) (L1393-1396).
    pub fn poles(&self) -> &Vec<f64> {
        &self.poles
    }

    /// OCCT StartPoint (L1398-1407).
    pub fn start_point(&self) -> f64 {
        if self.mults[0] == self.deg as i32 + 1 {
            self.poles[0]
        } else {
            self.value(self.first_parameter())
        }
    }

    /// OCCT Weight(Index) (L1409-1419).
    pub fn weight(&self, index: usize) -> f64 {
        assert!(
            index >= 1 && index <= self.poles.len(),
            "Standard_OutOfRange: Law_BSpline::Weight"
        );
        if self.is_rational() {
            self.weights.as_ref().unwrap()[index - 1]
        } else {
            1.0
        }
    }

    /// OCCT Weights(W) (L1421-1434).
    pub fn weights(&self) -> Vec<f64> {
        if self.is_rational() {
            self.weights.as_ref().unwrap().clone()
        } else {
            vec![1.0; self.poles.len()]
        }
    }

    /// OCCT IsRational (L1436-1439).
    pub fn is_rational(&self) -> bool {
        self.weights.is_some()
    }

    /// OCCT LocateU(U, ParametricTolerance, I1, I2, WithKnotRepetition)
    /// (L1441-1496).
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
            &self.flatknots
        } else {
            &self.knots
        };
        self.periodic_normalization(&mut new_u); // Attention a la periode
        let ufirst = cknots[0];
        let ulast = cknots[cknots.len() - 1];
        if (u - ufirst).abs() <= parametric_tolerance.abs() {
            *i1 = 1;
            *i2 = 1;
        } else if (u - ulast).abs() <= parametric_tolerance.abs() {
            *i1 = cknots.len() as i32;
            *i2 = *i1;
        } else if new_u < ufirst - parametric_tolerance.abs() {
            *i2 = 1;
            *i1 = 0;
        } else if new_u > ulast + parametric_tolerance.abs() {
            *i1 = cknots.len() as i32;
            *i2 = *i1 + 1;
        } else {
            *i1 = 1;
            bspl::hunt(cknots, new_u, i1);
            *i1 = (*i1).min(cknots.len() as i32).max(1);
            while *i1 + 1 <= cknots.len() as i32
                && (at(cknots, *i1 + 1) - new_u).abs() <= parametric_tolerance.abs()
            {
                *i1 += 1;
            }
            if (at(cknots, *i1) - new_u).abs() <= parametric_tolerance.abs() {
                *i2 = *i1;
            } else {
                *i2 = *i1 + 1;
            }
        }
    }

    /// OCCT MovePointAndTangent (L1498-1532).
    #[allow(clippy::too_many_arguments)]
    pub fn move_point_and_tangent(
        &mut self,
        u: f64,
        p: f64,
        tangent: f64,
        tolerance: f64,
        starting_condition: i32,
        ending_condition: i32,
    ) -> i32 {
        let mut new_poles = vec![0.0f64; self.poles.len()];
        let dimension = 1usize;
        let mut delta = 0.0;
        let mut delta_derivative = 0.0;
        self.d1(u, &mut delta, &mut delta_derivative);
        let delta = p - delta;
        let delta_derivative = tangent - delta_derivative;
        let err = bspl::move_point_and_tangent(
            u,
            dimension,
            &[delta],
            &[delta_derivative],
            tolerance,
            self.deg,
            starting_condition,
            ending_condition,
            &self.poles,
            self.weights.as_deref(),
            &self.flatknots,
            &mut new_poles,
        );
        if err == 0 {
            self.poles = new_poles;
        }
        err
    }

    /// OCCT Resolution(Tolerance3D, UTolerance) (L1534-1560).
    pub fn resolution(&self, tolerance3d: f64) -> f64 {
        bspl::resolution(
            1,
            &self.poles,
            self.weights.as_deref(),
            &self.flatknots,
            self.deg,
            tolerance3d,
        )
    }
}
