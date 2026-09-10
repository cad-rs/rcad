//! OCCT Bnd_Range (FoundationClasses TKMath Bnd package).
//!
//! 1:1 translation of `Bnd_Range.hxx` (L28-325) + `Bnd_Range.cxx`
//! (L21-181).  OCCT out-parameters with a bool success flag map to
//! `Option` returns; `NCollection_List<Bnd_Range>` in `Split` maps to a
//! `Vec<Bnd_Range>` appended in the same order.

/// OCCT Standard::RealSmall() = DBL_MIN.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

/// OCCT Standard::IsEqual(v1, v2) — `|v1 - v2| < RealSmall()`.
fn std_is_equal(v1: f64, v2: f64) -> bool {
    (v1 - v2).abs() < REAL_SMALL
}

/// OCCT enum Bnd_Range::IntersectStatus (hxx L84-89).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntersectStatus {
    /// No intersection with theVal+k*thePeriod.
    Out = 0,
    /// Range strictly contains theVal+k*thePeriod.
    In = 1,
    /// Range boundary coincides with theVal+k*thePeriod.
    Boundary = 2,
}

/// OCCT Bnd_Range — a range in 1D space restricted by two real values.  A
/// range can be VOID (no point included): `myLast < myFirst`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Range {
    my_first: f64,
    my_last: f64,
}

impl Default for Range {
    fn default() -> Self {
        Range::new()
    }
}

impl Range {
    /// OCCT Bnd_Range() — hxx L43-47: creates a VOID range.
    pub fn new() -> Self {
        Range {
            my_first: 0.0,
            my_last: -1.0,
        }
    }

    /// OCCT Bnd_Range(theMin, theMax) — hxx L50-56: never creates a VOID
    /// range (Standard_ConstructionError if Last < First).
    pub fn from_bounds(the_min: f64, the_max: f64) -> Self {
        if the_max < the_min {
            panic!("Last < First");
        }
        Range {
            my_first: the_min,
            my_last: the_max,
        }
    }

    /// OCCT Common(theOther) — cxx L21-36: replaces this with the
    /// common part of this and theOther.
    pub fn common(&mut self, the_other: &Range) {
        if the_other.is_void() {
            self.set_void();
            return;
        }
        if self.is_void() {
            return;
        }
        self.my_first = self.my_first.max(the_other.my_first);
        self.my_last = self.my_last.min(the_other.my_last);
    }

    /// OCCT Union(theOther) — cxx L40-61: joins this and theOther into one
    /// interval; returns false when impossible (void or separated).
    pub fn union_with(&mut self, the_other: &Range) -> bool {
        if self.is_void() || the_other.is_void() {
            return false;
        }
        if self.my_last < the_other.my_first {
            return false;
        }
        if self.my_first > the_other.my_last {
            return false;
        }
        self.my_first = self.my_first.min(the_other.my_first);
        self.my_last = self.my_last.max(the_other.my_last);
        true
    }

    /// OCCT IsIntersected(theVal, thePeriod) — cxx L65-128: status of the
    /// intersection with values theVal+k*thePeriod.
    pub fn is_intersected(&self, the_val: f64, the_period: f64) -> IntersectStatus {
        if self.is_void() {
            return IntersectStatus::Out;
        }

        let a_period = the_period.abs();
        let a_df = self.my_first - the_val;
        let a_dl = self.my_last - the_val;

        if a_period <= REAL_SMALL {
            let a_delta = a_df * a_dl;
            if std_is_equal(a_delta, 0.0) {
                return IntersectStatus::Boundary;
            }
            if a_delta > 0.0 {
                return IntersectStatus::Out;
            }
            return IntersectStatus::In;
        }

        // The interval [aDF/aPeriod, aDL/aPeriod] must contain an integer.
        let a_val1 = a_df / a_period;
        let a_val2 = a_dl / a_period;
        let a_par1 = a_val1.floor() as i32;
        let a_par2 = a_val2.floor() as i32;
        if a_par1 != a_par2 {
            // Interval (myFirst, myLast] intersects the seam value.
            if std_is_equal(a_val2, a_par2 as f64) {
                // aVal2 is an integer => myLast lies ON the seam.
                return IntersectStatus::Boundary;
            }
            return IntersectStatus::In;
        }

        // Here aPar1 == aPar2.
        if std_is_equal(a_val1, a_par1 as f64) {
            // aVal1 is an integer => myFirst lies ON the seam.
            return IntersectStatus::Boundary;
        }
        IntersectStatus::Out
    }

    /// OCCT Split(theVal, theList, thePeriod) — cxx L132-171: splits this
    /// range into sub-ranges at theVal (+ multiples of thePeriod), appended
    /// to theList in order.
    pub fn split(&self, the_val: f64, the_list: &mut Vec<Range>, the_period: f64) {
        let a_period = the_period.abs();
        if self.is_intersected(the_val, a_period) != IntersectStatus::In {
            the_list.push(*self);
            return;
        }

        let is_periodic = a_period > 0.0;

        if !is_periodic {
            the_list.push(Range::from_bounds(self.my_first, the_val));
            the_list.push(Range::from_bounds(the_val, self.my_last));
            return;
        }

        let mut a_val_prev = the_val + a_period * ((self.my_first - the_val) / a_period).ceil();

        // Now, (myFirst <= aValPrev < myFirst+aPeriod).
        if a_val_prev > self.my_first {
            the_list.push(Range::from_bounds(self.my_first, a_val_prev));
        }

        let mut a_val = a_val_prev + a_period;
        while a_val <= self.my_last {
            the_list.push(Range::from_bounds(a_val_prev, a_val));
            a_val_prev = a_val;
            a_val += a_period;
        }

        if a_val_prev < self.my_last {
            the_list.push(Range::from_bounds(a_val_prev, self.my_last));
        }
    }

    /// OCCT Add(theParameter) — hxx L100-110: extends this to include the
    /// parameter.
    pub fn add_parameter(&mut self, the_parameter: f64) {
        if self.is_void() {
            self.my_first = the_parameter;
            self.my_last = the_parameter;
            return;
        }
        self.my_first = self.my_first.min(the_parameter);
        self.my_last = self.my_last.max(the_parameter);
    }

    /// OCCT Add(theRange) — hxx L114-127: extends this to include both
    /// ranges.
    pub fn add_range(&mut self, the_range: &Range) {
        if the_range.is_void() {
            return;
        }
        if self.is_void() {
            *self = *the_range;
            return;
        }
        self.my_first = self.my_first.min(the_range.my_first);
        self.my_last = self.my_last.max(the_range.my_last);
    }

    /// OCCT GetMin() — hxx L131-140.
    pub fn get_min(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.my_first)
    }

    /// OCCT GetMax() — hxx L144-153.
    pub fn get_max(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.my_last)
    }

    /// OCCT GetBounds() — hxx L157-167.
    pub fn get_bounds(&self) -> Option<(f64, f64)> {
        if self.is_void() {
            return None;
        }
        Some((self.my_first, self.my_last))
    }

    /// OCCT GetIntermediatePoint(theLambda) — hxx L193-202.
    pub fn get_intermediate_point(&self, the_lambda: f64) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.my_first + the_lambda * (self.my_last - self.my_first))
    }

    /// OCCT Center() — hxx L206-211.
    pub fn center(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(0.5 * (self.my_first + self.my_last))
    }

    /// OCCT Delta() — hxx L214: negative for a VOID range.
    pub fn delta(&self) -> f64 {
        self.my_last - self.my_first
    }

    /// OCCT IsVoid() — hxx L217.
    pub fn is_void(&self) -> bool {
        self.my_last < self.my_first
    }

    /// OCCT SetVoid() — hxx L220-224.
    pub fn set_void(&mut self) {
        self.my_last = -1.0;
        self.my_first = 0.0;
    }

    /// OCCT Enlarge(theDelta) — hxx L227-236.
    pub fn enlarge(&mut self, the_delta: f64) {
        if self.is_void() {
            return;
        }
        self.my_first -= the_delta;
        self.my_last += the_delta;
    }

    /// OCCT Shifted(theVal) — hxx L239-242.
    pub fn shifted(&self, the_val: f64) -> Range {
        if !self.is_void() {
            Range::from_bounds(self.my_first + the_val, self.my_last + the_val)
        } else {
            Range::new()
        }
    }

    /// OCCT Shift(theVal) — hxx L245-252.
    pub fn shift(&mut self, the_val: f64) {
        if !self.is_void() {
            self.my_first += the_val;
            self.my_last += the_val;
        }
    }

    /// OCCT TrimFrom(theValLower) — hxx L256-262.
    pub fn trim_from(&mut self, the_val_lower: f64) {
        if !self.is_void() {
            self.my_first = self.my_first.max(the_val_lower);
        }
    }

    /// OCCT TrimTo(theValUpper) — hxx L266-272.
    pub fn trim_to(&mut self, the_val_upper: f64) {
        if !self.is_void() {
            self.my_last = self.my_last.min(the_val_upper);
        }
    }

    /// OCCT IsOut(theValue) — hxx L275-278.
    pub fn is_out(&self, the_value: f64) -> bool {
        self.is_void() || the_value < self.my_first || the_value > self.my_last
    }

    /// OCCT IsOut(theRange) — hxx L281-284.
    pub fn is_out_range(&self, the_range: &Range) -> bool {
        self.is_void()
            || the_range.is_void()
            || the_range.my_last < self.my_first
            || the_range.my_first > self.my_last
    }

    /// OCCT Contains(theValue) — hxx L287.
    pub fn contains(&self, the_value: f64) -> bool {
        !self.is_out(the_value)
    }

    /// OCCT Intersects(theRange) — hxx L290-293.
    pub fn intersects(&self, the_range: &Range) -> bool {
        !self.is_out_range(the_range)
    }

    /// OCCT Min() — hxx L297-302.
    pub fn min(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.my_first)
    }

    /// OCCT Max() — hxx L306-311.
    pub fn max(&self) -> Option<f64> {
        if self.is_void() {
            return None;
        }
        Some(self.my_last)
    }
}
