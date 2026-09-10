//! OCCT GeomFill_Line (TKGeomAlgo/GeomFill) — 1:1 port of GeomFill_Line.hxx
//! (L27-43) + GeomFill_Line.cxx (L24-34) + GeomFill_Line.lxx (accessors).
//!
//! Architecture mapping: the OCCT class derives from Standard_Transient
//! (reference-counted handle base); rcad has no transient base, so the class
//! maps to a plain struct passed by value / Arc by callers.

/// OCCT GeomFill_Line — class for instantiation of AppBlend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Line {
    /// OCCT myNbPoints.
    my_nb_points: i32,
}

impl Line {
    /// OCCT GeomFill_Line::GeomFill_Line() (GeomFill_Line.cxx L24-27).
    pub fn new() -> Self {
        Line { my_nb_points: 0 }
    }

    /// OCCT GeomFill_Line::GeomFill_Line(const int NbPoints) (L31-34).
    pub fn with_points(nb_points: i32) -> Self {
        Line {
            my_nb_points: nb_points,
        }
    }

    /// OCCT NbPoints (GeomFill_Line.lxx L24-28).
    pub fn nb_points(&self) -> i32 {
        self.my_nb_points
    }

    /// OCCT Point (GeomFill_Line.lxx L30-34).
    pub fn point(&self, index: i32) -> i32 {
        index
    }
}

impl Default for Line {
    fn default() -> Self {
        Self::new()
    }
}
