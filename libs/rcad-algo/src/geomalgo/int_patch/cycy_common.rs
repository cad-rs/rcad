//! Common data types and helpers for the IntCyCy (cylinder-cylinder) numeric
//! engine — 1:1 translation of the corresponding pieces of
//! `IntPatch_ImpImpIntersection.cxx` and `Bnd_Range`.
//!
//! rcad data-model notes:
//!   - OCCT `Bnd_Range` (a real-valued interval) -> `BndRange { f, l }`.
//!   - OCCT `math_Matrix` (1-based, here 3x3 / 3x5) -> `Mat3` / `Mat35`,
//!     column access `Col(i)` mapped to 0-based `col(i-1)`.
//!   - OCCT `Standard_Real::RealSmall()`/`IsEqual` are exact-equality helpers.

use rcad_kernel::precision::{INFINITE_VALUE, is_infinite_value};

/// OCCT Standard_Real::RealSmall() (Standard_Real.hxx L132-135) = DBL_MIN.
pub fn real_small() -> f64 {
    f64::MIN_POSITIVE
}

/// OCCT Standard_Real::RealEpsilon() (Standard_Real.hxx L161-164) = DBL_EPSILON.
pub fn real_epsilon() -> f64 {
    f64::EPSILON
}

/// OCCT Standard_Real::RealLast() (Standard_Real.hxx L179-182) = DBL_MAX.
pub fn real_last() -> f64 {
    f64::MAX
}

/// OCCT Standard_Real::RealFirst() (Standard_Real.hxx L167-170) = -DBL_MAX.
pub fn real_first() -> f64 {
    -f64::MAX
}

/// OCCT static const double aNulValue (IntPatch_ImpImpIntersection.cxx L3930):
/// if std::abs(a) <= aNulValue then it is considered that a = 0.
pub const A_NUL_VALUE: f64 = 1.0e-11;

/// OCCT Standard_Real::IsEqual(a, b) (Standard_Real.hxx L148-151):
/// `std::abs(a - b) < RealSmall()` — effectively exact equality.
pub fn is_equal(a: f64, b: f64) -> bool {
    (a - b).abs() < real_small()
}

/// OCCT Precision::Infinite() (Precision.hxx L371) = 2e100.
pub fn precision_infinite() -> f64 {
    INFINITE_VALUE
}

/// OCCT Precision::IsInfinite(R) (Precision.hxx L350-353).
pub fn precision_is_infinite(r: f64) -> bool {
    is_infinite_value(r)
}

/// OCCT ElCLib::InPeriod(Par, Ufirst, Ulast) — bring Par into [Ufirst, Ulast].
pub fn in_period(par: f64, u_first: f64, u_last: f64) -> f64 {
    let period = u_last - u_first;
    let mut x = par;
    while x < u_first {
        x += period;
    }
    while x > u_last {
        x -= period;
    }
    x
}

// ============================================================================
// Bnd_Range — 1D real interval (Bnd_Range.hxx / Bnd_Range.cxx, TKMath/Bnd)
// ============================================================================

/// OCCT Bnd_Range::IntersectStatus (Bnd_Range.hxx L84-89).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectStatus {
    Out = 0,
    In = 1,
    Boundary = 2,
}

/// OCCT Bnd_Range — a range in 1D space restricted by two real values.
/// A range is void when it contains no point (myLast < myFirst).
#[derive(Debug, Clone, Copy)]
pub struct BndRange {
    f: f64,
    l: f64,
}

impl BndRange {
    /// OCCT Bnd_Range() — default constructor creates a VOID range.
    pub fn new() -> Self {
        BndRange { f: 0.0, l: -1.0 }
    }

    /// OCCT Bnd_Range(theMin, theMax) — never creates a VOID range.
    pub fn with_bounds(min: f64, max: f64) -> Self {
        BndRange { f: min, l: max }
    }

    /// OCCT Bnd_Range::Common(theOther) (Bnd_Range.cxx L21-36).
    pub fn common(&mut self, other: &BndRange) {
        if other.is_void() {
            self.set_void();
            return;
        }
        if self.is_void() {
            return;
        }
        self.f = self.f.max(other.f);
        self.l = self.l.min(other.l);
    }

    /// OCCT Bnd_Range::Union(theOther) (Bnd_Range.cxx L40-61) — joins to one
    /// interval; returns false if the operation cannot be done (empty or
    /// separated ranges).
    pub fn union(&mut self, other: &BndRange) -> bool {
        if self.is_void() || other.is_void() {
            return false;
        }
        if self.l < other.f {
            return false;
        }
        if self.f > other.l {
            return false;
        }
        self.f = self.f.min(other.f);
        self.l = self.l.max(other.l);
        true
    }

    /// OCCT Bnd_Range::IsIntersected(theVal, thePeriod) (Bnd_Range.cxx L65-128).
    pub fn is_intersected(&self, val: f64, period: f64) -> IntersectStatus {
        if self.is_void() {
            return IntersectStatus::Out;
        }
        let a_period = period.abs();
        let a_df = self.f - val;
        let a_dl = self.l - val;
        if a_period <= real_small() {
            let a_delta = a_df * a_dl;
            if is_equal(a_delta, 0.0) {
                return IntersectStatus::Boundary;
            }
            if a_delta > 0.0 {
                return IntersectStatus::Out;
            }
            return IntersectStatus::In;
        }
        let a_val1 = a_df / a_period;
        let a_val2 = a_dl / a_period;
        let a_par1 = a_val1.floor() as i64;
        let a_par2 = a_val2.floor() as i64;
        if a_par1 != a_par2 {
            // Interval (myFirst, myLast] intersects seam-edge
            if is_equal(a_val2, a_par2 as f64) {
                // myLast lies ON the seam-edge
                return IntersectStatus::Boundary;
            }
            return IntersectStatus::In;
        }
        // Here, aPar1 == aPar2.
        if is_equal(a_val1, a_par1 as f64) {
            // myFirst lies ON the seam-edge
            return IntersectStatus::Boundary;
        }
        IntersectStatus::Out
    }

    /// OCCT Bnd_Range::Split(theVal, theList, thePeriod) (Bnd_Range.cxx L132-171).
    /// Splits <this> to several sub-ranges by theVal; new ranges are appended
    /// to `the_list`.
    pub fn split(&self, val: f64, the_list: &mut Vec<BndRange>, period: f64) {
        let a_period = period.abs();
        if self.is_intersected(val, a_period) != IntersectStatus::In {
            the_list.push(*self);
            return;
        }
        let is_periodic = a_period > 0.0;
        if !is_periodic {
            the_list.push(BndRange::with_bounds(self.f, val));
            the_list.push(BndRange::with_bounds(val, self.l));
            return;
        }
        let mut a_val_prev = val + a_period * ((self.f - val) / a_period).ceil();
        // Now, (myFirst <= aValPrev < myFirst+aPeriod).
        if a_val_prev > self.f {
            the_list.push(BndRange::with_bounds(self.f, a_val_prev));
        }
        let mut a_val = a_val_prev + a_period;
        while a_val <= self.l {
            the_list.push(BndRange::with_bounds(a_val_prev, a_val));
            a_val_prev = a_val;
            a_val += a_period;
        }
        if a_val_prev < self.l {
            the_list.push(BndRange::with_bounds(a_val_prev, self.l));
        }
    }

    /// OCCT Bnd_Range::Add(theParameter) (Bnd_Range.hxx L100-110).
    pub fn add(&mut self, param: f64) {
        if self.is_void() {
            self.f = param;
            self.l = param;
            return;
        }
        self.f = self.f.min(param);
        self.l = self.l.max(param);
    }

    /// OCCT Bnd_Range::Add(theRange) (Bnd_Range.hxx L114-127).
    pub fn add_range(&mut self, range: &BndRange) {
        if range.is_void() {
            return;
        }
        if self.is_void() {
            *self = *range;
            return;
        }
        self.f = self.f.min(range.f);
        self.l = self.l.max(range.l);
    }

    /// OCCT Bnd_Range::GetMin(thePar) (Bnd_Range.hxx L131-140).
    pub fn get_min(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.f)
    }

    /// OCCT Bnd_Range::GetMax(thePar) (Bnd_Range.hxx L144-153).
    pub fn get_max(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.l)
    }

    /// OCCT Bnd_Range::GetBounds(theFirstPar, theLastPar) (Bnd_Range.hxx L157-167).
    pub fn get_bounds(&self) -> Option<(f64, f64)> {
        if self.is_void() {
            return None;
        }
        Some((self.f, self.l))
    }

    /// OCCT Bnd_Range::Delta() (Bnd_Range.hxx L214) — MAX-MIN (negative for VOID).
    pub fn delta(&self) -> f64 {
        self.l - self.f
    }

    /// OCCT Bnd_Range::IsVoid() (Bnd_Range.hxx L217).
    pub fn is_void(&self) -> bool {
        self.l < self.f
    }

    /// OCCT Bnd_Range::SetVoid() (Bnd_Range.hxx L220-224).
    pub fn set_void(&mut self) {
        self.l = -1.0;
        self.f = 0.0;
    }

    /// OCCT Bnd_Range::Shifted(theVal) (Bnd_Range.hxx L239-242).
    pub fn shifted(&self, val: f64) -> BndRange {
        if !self.is_void() {
            BndRange::with_bounds(self.f + val, self.l + val)
        } else {
            BndRange::new()
        }
    }

    /// OCCT Bnd_Range::Shift(theVal) (Bnd_Range.hxx L245-252).
    pub fn shift(&mut self, val: f64) {
        if !self.is_void() {
            self.f += val;
            self.l += val;
        }
    }

    /// OCCT Bnd_Range::Enlarge(theDelta) (Bnd_Range.hxx L227-236) — extends
    /// this range to the given value (in both sides).
    pub fn enlarge(&mut self, delta: f64) {
        if self.is_void() {
            return;
        }
        self.f -= delta;
        self.l += delta;
    }
}

// ============================================================================
// WLine — IntSurf_LineOn2S under construction (1-based OCCT semantics)
// ============================================================================

/// OCCT IntSurf_LineOn2S — the walking-line point list (the WLine's Curve()).
/// The OCCT NCollection_Sequence is 1-based; this wrapper keeps the same
/// 1-based `value`/`insert_before` semantics so the IntCyCy engine reads 1:1.
#[derive(Debug, Clone, Default)]
pub struct WLine {
    points: Vec<crate::geomalgo::int_patch::WLinePnt>,
}

impl WLine {
    pub fn new() -> Self {
        WLine { points: Vec::new() }
    }
    /// OCCT IntSurf_LineOn2S::NbPoints().
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }
    /// OCCT IntSurf_LineOn2S::Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> &crate::geomalgo::int_patch::WLinePnt {
        &self.points[index - 1]
    }
    /// OCCT IntSurf_LineOn2S::Value(Index, P) — 1-based setter.
    pub fn set_value(&mut self, index: usize, p: crate::geomalgo::int_patch::WLinePnt) {
        self.points[index - 1] = p;
    }
    /// OCCT IntSurf_LineOn2S::Append(P).
    pub fn append(&mut self, p: crate::geomalgo::int_patch::WLinePnt) {
        self.points.push(p);
    }
    /// OCCT IntSurf_LineOn2S::InsertBefore(Index, P) — 1-based.
    pub fn insert_before(&mut self, index: usize, p: crate::geomalgo::int_patch::WLinePnt) {
        if index > self.points.len() {
            self.points.push(p);
        } else {
            self.points.insert(index - 1, p);
        }
    }
    /// OCCT IntSurf_LineOn2S::RemovePoint(Index) — 1-based.
    pub fn remove_point(&mut self, index: usize) {
        self.points.remove(index - 1);
    }
    /// Consume the wrapped points as the final `IntPatchLine::wline_pnts`.
    pub fn into_points(self) -> Vec<crate::geomalgo::int_patch::WLinePnt> {
        self.points
    }
    /// Expose as a slice (for the WLine point count / iteration).
    pub fn as_slice(&self) -> &[crate::geomalgo::int_patch::WLinePnt] {
        &self.points
    }
}

// ============================================================================
// math_Matrix — small 3x3 / 3x5 helpers (1-based OCCT -> 0-based here)
// ============================================================================

/// OCCT math_Matrix(1, 3, 1, 3) — used by VBoundaryPrecise / StepComputing.
#[derive(Debug, Clone)]
pub struct Mat3 {
    m: [[f64; 3]; 3],
}

impl Mat3 {
    pub fn new() -> Self {
        Mat3 { m: [[0.0; 3]; 3] }
    }

    /// OCCT math_Matrix::SetCol(i, v) — 1-based column index.
    pub fn set_col(&mut self, col: usize, v: [f64; 3]) {
        self.m[col - 1] = v;
    }

    /// OCCT math_Matrix::Col(i) — 1-based column index.
    pub fn col(&self, col: usize) -> [f64; 3] {
        self.m[col - 1]
    }

    /// OCCT math_Matrix::Determinant().
    pub fn determinant(&self) -> f64 {
        let m = self.m;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
}

/// OCCT math_Matrix(1, 3, 1, 5) — used by StepComputing (3 rows x 5 cols).
#[derive(Debug, Clone)]
pub struct Mat35 {
    m: [[f64; 5]; 3],
}

impl Mat35 {
    pub fn new() -> Self {
        Mat35 { m: [[0.0; 5]; 3] }
    }

    /// OCCT math_Matrix::SetCol(i, v) — 1-based column index.
    pub fn set_col(&mut self, col: usize, v: [f64; 3]) {
        for r in 0..3 {
            self.m[r][col - 1] = v[r];
        }
    }

    /// OCCT math_Matrix::Col(i) — 1-based column index.
    pub fn col(&self, col: usize) -> [f64; 3] {
        [self.m[0][col - 1], self.m[1][col - 1], self.m[2][col - 1]]
    }

    /// OCCT math_Matrix::operator()(i, j) — 1-based (row, column).
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.m[row - 1][col - 1]
    }
}

// ============================================================================
// ShortCosForm / InscribePoint / InscribeInterval / ExcludeNearElements
// ============================================================================

/// OCCT ShortCosForm (IntPatch_ImpImpIntersection.cxx L5197-5250).
/// Represents theCosFactor*cosA + theSinFactor*sinA as theCoeff*cos(A-theAngle)
/// when possible (all angles in radians).  Returns (theCoeff, theAngle).
#[allow(unused_assignments)]
pub fn short_cos_form(cos_factor: f64, sin_factor: f64) -> (f64, f64) {
    let coeff = (cos_factor * cos_factor + sin_factor * sin_factor).sqrt();
    let mut angle = 0.0;
    if is_equal(coeff, 0.0) {
        angle = 0.0;
        return (coeff, angle);
    }
    angle = (cos_factor.abs() / coeff).acos();
    if sin_factor > 0.0 {
        if is_equal(cos_factor, 0.0) {
            angle = std::f64::consts::PI / 2.0;
        } else if cos_factor < 0.0 {
            angle = std::f64::consts::PI - angle;
        }
    } else if is_equal(sin_factor, 0.0) {
        if cos_factor < 0.0 {
            angle = std::f64::consts::PI;
        }
    }
    if sin_factor < 0.0 {
        if cos_factor > 0.0 {
            angle = 2.0 * std::f64::consts::PI - angle;
        } else if is_equal(cos_factor, 0.0) {
            angle = 3.0 * std::f64::consts::PI / 2.0;
        } else if cos_factor < 0.0 {
            angle = std::f64::consts::PI + angle;
        }
    }
    (coeff, angle)
}

/// OCCT InscribePoint (IntPatch_ImpImpIntersection.cxx L5453-5504).
/// Brings theUGiven into [theUfTarget, theUlTarget] (periodically) if possible;
/// returns false when the point cannot be inscribed.
pub fn inscribe_point(
    uf_target: f64,
    ul_target: f64,
    u_given: &mut f64,
    tol_2d: f64,
    period: f64,
    fl_force: bool,
) -> bool {
    if precision_is_infinite(*u_given) {
        return false;
    }
    if (uf_target - *u_given <= tol_2d) && (*u_given - ul_target <= tol_2d) {
        // It has already been inscribed.
        if fl_force {
            let mut u_temp = *u_given + period;
            if (uf_target - u_temp <= tol_2d) && (u_temp - ul_target <= tol_2d) {
                *u_given = u_temp;
                return true;
            }
            u_temp = *u_given - period;
            if (uf_target - u_temp <= tol_2d) && (u_temp - ul_target <= tol_2d) {
                *u_given = u_temp;
            }
        }
        return true;
    }
    let a_uf = uf_target - tol_2d;
    let a_ul = a_uf + period;
    *u_given = in_period(*u_given, a_uf, a_ul);
    (uf_target - *u_given <= tol_2d) && (*u_given - ul_target <= tol_2d)
}

/// OCCT InscribeInterval (IntPatch_ImpImpIntersection.cxx L5505-5566).
/// Shifts theRange to make at least one of its boundaries in
/// [theUfTarget, theUlTarget].
#[allow(unused_assignments)]
pub fn inscribe_interval(
    uf_target: f64,
    ul_target: f64,
    the_range: &mut BndRange,
    tol_2d: f64,
    period: f64,
) -> bool {
    let mut u_par = 0.0;
    let Some(min_val) = the_range.get_min() else {
        return false;
    };
    u_par = min_val;
    let a_delta = the_range.delta();
    let fl_force = (ul_target - u_par).abs() < tol_2d;
    if inscribe_point(uf_target, ul_target, &mut u_par, tol_2d, period, fl_force) {
        the_range.set_void();
        the_range.add(u_par);
        the_range.add(u_par + a_delta);
    } else {
        let Some(max_val) = the_range.get_max() else {
            return false;
        };
        u_par = max_val;
        let fl_force = (uf_target - u_par).abs() < tol_2d;
        if inscribe_point(uf_target, ul_target, &mut u_par, tol_2d, period, fl_force) {
            the_range.set_void();
            the_range.add(u_par);
            the_range.add(u_par - a_delta);
        } else {
            return false;
        }
    }
    true
}

/// OCCT ExcludeNearElements (IntPatch_ImpImpIntersection.cxx L5567-5614).
/// Checks if theArr contains two almost equal elements; if so, one of the equal
/// elements is excluded (made infinite).  theArr must be sorted ascending and
/// every non-infinite element is in [0, T].  Returns true if any element changed.
pub fn exclude_near_elements(
    arr: &mut [f64],
    n_of_members: usize,
    usurf1f: f64,
    usurf1l: f64,
    tol: f64,
) -> bool {
    let mut ret_val = false;
    for i in 1..n_of_members {
        let mut an_a = arr[i];
        let an_b = arr[i - 1];
        // Here, anA >= anB
        if precision_is_infinite(an_a) {
            break;
        }
        if (an_a - an_b) < tol {
            if (an_b != 0.0) && (an_b != usurf1f) && (an_b != usurf1l) {
                an_a = (an_a + an_b) / 2.0;
            } else {
                an_a = an_b;
            }
            // Make this element infinite and forget it.
            arr[i - 1] = precision_infinite();
            arr[i] = an_a;
            ret_val = true;
        }
    }
    ret_val
}

// ============================================================================
// MathRoot::Brent — bracketed root finding (MathRoot_Brent.hxx, TKMath/MathRoot)
// ============================================================================

/// OCCT MathUtils::Status (MathUtils_Types.hxx L29-41) — the subset used by the
/// IntCyCy numeric engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverStatus {
    Ok,
    NotConverged,
    MaxIterations,
    NumericalError,
    InvalidInput,
}

/// OCCT MathUtils::ScalarResult (MathUtils_Types.hxx L45-59).
#[derive(Debug, Clone, Copy)]
pub struct ScalarResult {
    pub status: SolverStatus,
    pub root: Option<f64>,
    pub value: Option<f64>,
}

/// OCCT MathUtils::Config (MathUtils_Config.hxx L106-141) — the subset used by
/// MathRoot::Brent.
#[derive(Debug, Clone, Copy)]
pub struct BrentConfig {
    pub max_iterations: i32,
    pub x_tolerance: f64,
    pub f_tolerance: f64,
}

impl BrentConfig {
    /// OCCT Config() default: MaxIterations = 100, XTolerance = FTolerance = 1e-10.
    pub fn new() -> Self {
        BrentConfig {
            max_iterations: 100,
            x_tolerance: 1.0e-10,
            f_tolerance: 1.0e-10,
        }
    }
}

impl Default for BrentConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT MathUtils::THE_EPSILON (MathUtils_Core.hxx L28) = DBL_EPSILON.
const THE_EPSILON: f64 = f64::EPSILON;
/// OCCT MathUtils::THE_ZERO_TOL (MathUtils_Core.hxx L32) = 1.0e-15.
const THE_ZERO_TOL: f64 = 1.0e-15;

/// OCCT MathRoot::Brent (MathRoot_Brent.hxx L44-193) — Brent's method for root
/// finding, combining bisection, secant, and inverse quadratic interpolation.
/// `f(x)` returns the function value, or None if the evaluation fails.
#[allow(unused_assignments)]
pub fn brent_root(
    f: &mut impl FnMut(f64) -> Option<f64>,
    lower: f64,
    upper: f64,
    config: &BrentConfig,
) -> ScalarResult {
    let mut result = ScalarResult { status: SolverStatus::NotConverged, root: None, value: None };

    let mut a_a = lower;
    let mut a_b = upper;
    let mut a_fa = 0.0;
    let mut a_fb = 0.0;

    // Evaluate at endpoints
    match f(a_a) {
        None => {
            result.status = SolverStatus::NumericalError;
            return result;
        }
        Some(v) => a_fa = v,
    }
    match f(a_b) {
        None => {
            result.status = SolverStatus::NumericalError;
            return result;
        }
        Some(v) => a_fb = v,
    }

    // Check that bracket is valid (sign change)
    if a_fa * a_fb > 0.0 {
        result.status = SolverStatus::InvalidInput;
        return result;
    }

    // Ensure |f(a)| >= |f(b)| (b is the better approximation)
    if a_fa.abs() < a_fb.abs() {
        std::mem::swap(&mut a_a, &mut a_b);
        std::mem::swap(&mut a_fa, &mut a_fb);
    }

    let mut a_c = a_a; // Previous iterate
    let mut a_fc = a_fa;
    let mut a_d = a_b - a_a; // Step size
    let mut a_e = a_d; // Previous step size

    for _ in 0..config.max_iterations {
        let a_tol = 2.0 * THE_EPSILON * a_b.abs() + 0.5 * config.x_tolerance;
        let a_m = 0.5 * (a_c - a_b);

        if a_fb.abs() < config.f_tolerance || a_fb == 0.0 || a_m.abs() <= a_tol {
            result.status = SolverStatus::Ok;
            result.root = Some(a_b);
            result.value = Some(a_fb);
            return result;
        }

        let mut a_s = 0.0; // New approximation

        // Try inverse quadratic interpolation if we have three distinct points
        if (a_fa - a_fc).abs() > THE_ZERO_TOL && (a_fb - a_fc).abs() > THE_ZERO_TOL {
            // Inverse quadratic interpolation
            a_s = a_a * a_fb * a_fc / ((a_fa - a_fb) * (a_fa - a_fc))
                + a_b * a_fa * a_fc / ((a_fb - a_fa) * (a_fb - a_fc))
                + a_c * a_fa * a_fb / ((a_fc - a_fa) * (a_fc - a_fb));
        } else {
            // Secant method
            a_s = a_b - a_fb * (a_b - a_a) / (a_fb - a_fa);
        }

        // Decide whether to accept the interpolation step
        let mut use_interp = false;

        // Check if s is between (3a+b)/4 and b
        let a_bound1 = (3.0 * a_a + a_b) / 4.0;
        if a_s > a_bound1.min(a_b) && a_s < a_bound1.max(a_b) {
            // Accept interpolation if step is smaller than half the previous step.
            if (a_s - a_b).abs() < a_e.abs() / 2.0 {
                use_interp = true;
            }
        }

        if !use_interp {
            // Bisection step
            a_s = a_b + a_m;
            a_e = a_m;
            a_d = a_m;
        } else {
            a_e = a_d;
            a_d = a_s - a_b;
        }

        // Update previous values
        a_a = a_b;
        a_fa = a_fb;

        // Compute new point, ensuring minimum step
        if a_d.abs() > a_tol {
            a_b = a_s;
        } else {
            a_b += if a_m > 0.0 { a_tol } else { -a_tol };
        }

        // Evaluate function at new point
        match f(a_b) {
            None => {
                result.status = SolverStatus::NumericalError;
                result.root = Some(a_b);
                return result;
            }
            Some(v) => a_fb = v,
        }

        // Update bracket
        if a_fb * a_fc > 0.0 {
            a_c = a_a;
            a_fc = a_fa;
            a_d = a_b - a_a;
            a_e = a_d;
        } else if a_fc.abs() < a_fb.abs() {
            // Swap b and c if c is better.
            a_a = a_b;
            a_fa = a_fb;
            std::mem::swap(&mut a_b, &mut a_c);
            std::mem::swap(&mut a_fb, &mut a_fc);
        }
    }

    // Maximum iterations reached
    result.status = SolverStatus::MaxIterations;
    result.root = Some(a_b);
    result.value = Some(a_fb);
    result
}
