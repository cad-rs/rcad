//! OCCT Intf_SectionLine (TKGeomAlgo Intf package).
//!
//! 1:1 translation of `Intf_SectionLine.hxx` (L31-97) +
//! `Intf_SectionLine.cxx` (L22-165) + `Intf_SectionLine.lxx` (L17-35).

use super::intf::IntfSectionPoint;

/// OCCT Intf_SectionLine — a polyline of intersection between two polyhedra
/// as a sequence of intersection points.
#[derive(Debug, Clone)]
pub struct SectionLine {
    /// NCollection_Sequence<Intf_SectionPoint> (1-based in OCCT).
    my_points: Vec<IntfSectionPoint>,
    closed: bool,
}

impl Default for SectionLine {
    fn default() -> Self {
        Self::new()
    }
}

impl SectionLine {
    /// OCCT Intf_SectionLine() — cxx L22-25.
    pub fn new() -> Self {
        SectionLine {
            my_points: Vec::new(),
            closed: false,
        }
    }

    /// OCCT Intf_SectionLine(Other) — cxx L29-33: copies the points, NOT
    /// the closed flag (kept verbatim).
    pub fn from_other(other: &SectionLine) -> Self {
        SectionLine {
            my_points: other.my_points.clone(),
            closed: false,
        }
    }

    /// OCCT NumberOfPoints() — lxx L17-25.
    pub fn number_of_points(&self) -> usize {
        self.my_points.len()
    }

    /// OCCT Append(Pi) — cxx L37-40.
    pub fn append(&mut self, pi: &IntfSectionPoint) {
        self.my_points.push(*pi);
    }

    /// OCCT Append(LS) — cxx L44-47: concatenates LS at the end.
    pub fn append_line(&mut self, ls: &mut SectionLine) {
        self.my_points.append(&mut ls.my_points);
    }

    /// OCCT Prepend(Pi) — cxx L51-54.
    pub fn prepend(&mut self, pi: &IntfSectionPoint) {
        self.my_points.insert(0, *pi);
    }

    /// OCCT Prepend(LS) — cxx L58-61: concatenates LS at the beginning.
    pub fn prepend_line(&mut self, ls: &mut SectionLine) {
        let mut pts = std::mem::take(&mut ls.my_points);
        pts.append(&mut self.my_points);
        self.my_points = pts;
    }

    /// OCCT Reverse() — cxx L65-68.
    pub fn reverse(&mut self) {
        self.my_points.reverse();
    }

    /// OCCT Close() — cxx L72-75.
    pub fn close(&mut self) {
        self.closed = true;
    }

    /// OCCT GetPoint(index) — cxx L79-82 (1-based).
    pub fn get_point(&self, index: usize) -> IntfSectionPoint {
        self.my_points[index - 1]
    }

    /// OCCT IsClosed() — cxx L86-95: OCCT ignores the closed flag here and
    /// compares the first and last points (comment kept in spirit); the
    /// flag itself is dead state in this version.
    pub fn is_closed(&self) -> bool {
        match (self.my_points.first(), self.my_points.last()) {
            (Some(f), Some(l)) => f.is_equal(l),
            _ => false,
        }
    }

    /// OCCT Contains(ThePI) — cxx L99-109.
    pub fn contains(&self, the_pi: &IntfSectionPoint) -> bool {
        for i in 1..=self.my_points.len() {
            if the_pi.is_equal(&self.my_points[i - 1]) {
                return true;
            }
        }
        false
    }

    /// OCCT IsEnd(ThePI) — cxx L113-124: 1 for the beginning,
    /// `NumberOfPoints` for the end (never exactly 2 for a single-segment
    /// line — the callers test `> 1`), otherwise 0.
    pub fn is_end(&self, the_pi: &IntfSectionPoint) -> i32 {
        if self.my_points[0].is_equal(the_pi) {
            return 1;
        }
        let last = self.my_points.len() - 1;
        if self.my_points[last].is_equal(the_pi) {
            return self.my_points.len() as i32;
        }
        0
    }

    /// OCCT IsEqual(Other) — cxx L128-142.
    pub fn is_equal(&self, other: &SectionLine) -> bool {
        if self.my_points.len() != other.my_points.len() {
            return false;
        }
        for i in 1..=self.my_points.len() {
            if !self.my_points[i - 1].is_equal(&other.my_points[i - 1]) {
                return false;
            }
        }
        true
    }
}
