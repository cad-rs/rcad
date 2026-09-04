//! OCCT Intf_TangentZone (TKGeomAlgo Intf package).
//!
//! 1:1 translation of `Intf_TangentZone.hxx` (L31-104) +
//! `Intf_TangentZone.cxx` (L26-330) + `Intf_TangentZone.lxx` (L17-50).

use super::intf::IntfSectionPoint;

/// OCCT Intrv-style bounds helpers from Standard (Intf_TangentZone.cxx
/// L28-29 use RealLast/RealFirst).
const REAL_LAST: f64 = f64::MAX;
const REAL_FIRST: f64 = -f64::MAX;

/// OCCT Intf_TangentZone — a zone of tangence between polygons or
/// polyhedra, as a sequence of intersection points plus the parameter
/// ranges on the two interfering arguments.
#[derive(Debug, Clone)]
pub struct TangentZone {
    /// NCollection_Sequence<Intf_SectionPoint> (1-based in OCCT).
    result: Vec<IntfSectionPoint>,
    param_on_first_min: f64,
    param_on_first_max: f64,
    param_on_second_min: f64,
    param_on_second_max: f64,
}

impl Default for TangentZone {
    fn default() -> Self {
        Self::new()
    }
}

impl TangentZone {
    /// OCCT Intf_TangentZone() — cxx L26-30: an empty tangent zone with
    /// inverted parameter bounds.
    pub fn new() -> Self {
        TangentZone {
            result: Vec::new(),
            param_on_first_min: REAL_LAST,
            param_on_second_min: REAL_LAST,
            param_on_first_max: REAL_FIRST,
            param_on_second_max: REAL_FIRST,
        }
    }

    fn update_param_bounds(&mut self, pi: &IntfSectionPoint) {
        // The min/max update block shared by Append/InsertAfter/InsertBefore
        // (cxx L40-56, L178-194, L202-218).
        if self.param_on_first_min > pi.param_on_first() {
            self.param_on_first_min = pi.param_on_first();
        }
        if self.param_on_second_min > pi.param_on_second() {
            self.param_on_second_min = pi.param_on_second();
        }
        if self.param_on_first_max < pi.param_on_first() {
            self.param_on_first_max = pi.param_on_first();
        }
        if self.param_on_second_max < pi.param_on_second() {
            self.param_on_second_max = pi.param_on_second();
        }
    }

    /// OCCT NumberOfPoints() — lxx L17-25.
    pub fn number_of_points(&self) -> usize {
        self.result.len()
    }

    /// OCCT Append(Pi) — cxx L37-57.
    pub fn append(&mut self, pi: &IntfSectionPoint) {
        self.result.push(*pi);
        self.update_param_bounds(pi);
    }

    /// OCCT Append(Tzi) — cxx L64-71: merges the tangent zone via
    /// PolygonInsert of each point.
    pub fn append_zone(&mut self, tzi: &TangentZone) {
        for ipi in 1..=tzi.number_of_points() {
            self.polygon_insert(&tzi.get_point(ipi));
        }
    }

    /// OCCT Insert(Pi) — cxx L79-119: the OCCT body is commented out and
    /// always returns false; kept verbatim.
    pub fn insert(&mut self, _pi: &IntfSectionPoint) -> bool {
        let inserted = false;
        inserted
    }

    /// OCCT PolygonInsert(Pi) — cxx L126-171.
    pub fn polygon_insert(&mut self, pi: &IntfSectionPoint) {
        let nbp_tz = self.number_of_points();
        if nbp_tz == 0 {
            self.append(pi);
            return;
        }
        if pi.param_on_first() >= self.param_on_first_max {
            self.append(pi);
        } else if pi.param_on_first() >= self.param_on_first_min {
            self.insert_before(1, pi);
        } else {
            // OCCT L167-169: points are appended unsorted ("On met les
            // points sans les classer").
            self.append(pi);
        }
    }

    /// OCCT InsertAfter(Index, Pi) — cxx L175-195.
    pub fn insert_after(&mut self, index: usize, pi: &IntfSectionPoint) {
        self.result.insert(index, *pi); // OCCT InsertAfter(index) == insert at index
        self.update_param_bounds(pi);
    }

    /// OCCT InsertBefore(Index, Pi) — cxx L199-219.
    pub fn insert_before(&mut self, index: usize, pi: &IntfSectionPoint) {
        self.result.insert(index - 1, *pi);
        self.update_param_bounds(pi);
    }

    /// OCCT GetPoint(Index) — cxx L226-229 (1-based).
    pub fn get_point(&self, index: usize) -> IntfSectionPoint {
        self.result[index - 1]
    }

    /// OCCT IsEqual(Other) — cxx L233-248.
    pub fn is_equal(&self, other: &TangentZone) -> bool {
        if self.result.len() != other.result.len() {
            return false;
        }
        for i in 1..=self.result.len() {
            if !self.result[i - 1].is_equal(&other.result[i - 1]) {
                return false;
            }
        }
        true
    }

    /// OCCT Contains(ThePI) — cxx L252-263.
    pub fn contains(&self, the_pi: &IntfSectionPoint) -> bool {
        for i in 1..=self.result.len() {
            if the_pi.is_equal(&self.result[i - 1]) {
                return true;
            }
        }
        false
    }

    /// OCCT ParamOnFirst(paraMin, paraMax) — lxx L33-41.
    pub fn param_on_first(&self) -> (f64, f64) {
        (self.param_on_first_min, self.param_on_first_max)
    }

    /// OCCT ParamOnSecond(paraMin, paraMax) — lxx L47-55.
    pub fn param_on_second(&self) -> (f64, f64) {
        (self.param_on_second_min, self.param_on_second_max)
    }

    /// OCCT InfoFirst(segMin, paraMin, segMax, paraMax) — cxx L267-274.
    pub fn info_first(&self) -> (i32, f64, i32, f64) {
        let (para_min, para_max) = self.param_on_first();
        let seg_min = para_min.trunc() as i32;
        let para_min = para_min - seg_min as f64;
        let seg_max = para_max.trunc() as i32;
        let para_max = para_max - seg_max as f64;
        (seg_min, para_min, seg_max, para_max)
    }

    /// OCCT InfoSecond(segMin, paraMin, segMax, paraMax) — cxx L278-285.
    pub fn info_second(&self) -> (i32, f64, i32, f64) {
        let (para_min, para_max) = self.param_on_second();
        let seg_min = para_min.trunc() as i32;
        let para_min = para_min - seg_min as f64;
        let seg_max = para_max.trunc() as i32;
        let para_max = para_max - seg_max as f64;
        (seg_min, para_min, seg_max, para_max)
    }

    /// OCCT RangeContains(ThePI) — cxx L289-296.
    pub fn range_contains(&self, the_pi: &IntfSectionPoint) -> bool {
        let (a, b) = self.param_on_first();
        let (c, d) = self.param_on_second();
        a <= the_pi.param_on_first()
            && the_pi.param_on_first() <= b
            && c <= the_pi.param_on_second()
            && the_pi.param_on_second() <= d
    }

    /// OCCT HasCommonRange(Other) — cxx L300-311.
    pub fn has_common_range(&self, other: &TangentZone) -> bool {
        let (a1, b1) = self.param_on_first();
        let (a2, b2) = self.param_on_second();
        let (c1, d1) = other.param_on_first();
        let (c2, d2) = other.param_on_second();

        ((c1 <= a1 && a1 <= d1) || (c1 <= b1 && b1 <= d1) || (a1 <= c1 && c1 <= b1))
            && ((c2 <= a2 && a2 <= d2) || (c2 <= b2 && b2 <= d2) || (a2 <= c2 && c2 <= b2))
    }
}
