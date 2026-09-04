//! Adaptor2d_Line2d (TKG2d) — the 2D line adaptor used by the TopolTool
//! restriction iteration.
//!
//! 1:1 translation of OCCT `Adaptor2d_Line2d.hxx` (L17-142) + `.cxx`
//! (L35-303): an iso-line in UV space trimmed to [myUfirst, myUlast].
//! Implements [`Curve2dAdaptor`] (the OCCT Adaptor2d_Curve2d abstract
//! interface).

use glam::DVec2;
use rcad_kernel::geom::Line2d;

use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType};
use rcad_kernel::math::GeomAbsShape;

/// OCCT Adaptor2d_Line2d — "used by the TopolTool to trim a surface".
#[derive(Debug, Clone, Copy)]
pub struct Line2dAdaptor {
    my_ufirst: f64,
    my_ulast: f64,
    /// OCCT gp_Ax2d myAx2d (location + unit direction).
    my_ax2d: Line2d,
}

impl Line2dAdaptor {
    /// OCCT Adaptor2d_Line2d() (cxx L35-39) — empty.
    pub fn new() -> Self {
        Line2dAdaptor {
            my_ufirst: 0.0,
            my_ulast: 0.0,
            my_ax2d: Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::X,
            },
        }
    }

    /// OCCT Adaptor2d_Line2d(P, D, UFirst, ULast) (cxx L43-51).
    pub fn new_pnt_dir(p: DVec2, d: DVec2, ufirst: f64, ulast: f64) -> Self {
        Line2dAdaptor {
            my_ufirst: ufirst,
            my_ulast: ulast,
            my_ax2d: Line2d {
                origin: p,
                direction: d,
            },
        }
    }

    /// OCCT Load(L) (cxx L68-73) — infinite range.
    pub fn load(&mut self, l: &Line2d) {
        self.my_ax2d = *l;
        self.my_ufirst = -rcad_kernel::precision::INFINITE_VALUE;
        self.my_ulast = rcad_kernel::precision::INFINITE_VALUE;
    }

    /// OCCT Load(L, Fi, La) (cxx L77-82).
    pub fn load_range(&mut self, l: &Line2d, fi: f64, la: f64) {
        self.my_ax2d = *l;
        self.my_ufirst = fi;
        self.my_ulast = la;
    }

    /// OCCT Trim(First, Last, Tol) (cxx L125-132) — a fresh adaptor with the
    /// same axis over [First, Last].
    pub fn trim(&self, first: f64, last: f64) -> Line2dAdaptor {
        let mut hl = Line2dAdaptor::new();
        let l = self.line();
        hl.load_range(&l, first, last);
        hl
    }

    /// OCCT Line() (cxx L229-232) — the infinite supporting line.
    pub fn line(&self) -> Line2d {
        self.my_ax2d
    }
}

impl Default for Line2dAdaptor {
    fn default() -> Self {
        Line2dAdaptor::new()
    }
}

impl Curve2dAdaptor for Line2dAdaptor {
    /// OCCT FirstParameter() (cxx L86-89).
    fn first_parameter(&self) -> f64 {
        self.my_ufirst
    }

    /// OCCT LastParameter() (cxx L92-96).
    fn last_parameter(&self) -> f64 {
        self.my_ulast
    }

    /// OCCT Value(X) (cxx L157-161) — ElCLib::LineValue.
    fn value(&self, x: f64) -> DVec2 {
        self.my_ax2d.origin + x * self.my_ax2d.direction
    }

    /// OCCT D1(X, P, V) (cxx L171-174) — ElCLib::LineD1.
    fn d1(&self, x: f64) -> (DVec2, DVec2) {
        (self.value(x), self.my_ax2d.direction)
    }

    /// OCCT D2(X, P, V1, V2) (cxx L178-182) — V2 is null.
    fn d2(&self, x: f64) -> (DVec2, DVec2, DVec2) {
        (self.value(x), self.my_ax2d.direction, DVec2::ZERO)
    }

    /// OCCT D3(X, P, V1, V2, V3) (cxx L186-195) — V2/V3 are null.
    fn d3(&self, x: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        (self.value(x), self.my_ax2d.direction, DVec2::ZERO, DVec2::ZERO)
    }

    /// OCCT DN(U, N) (cxx L200-211).
    fn dn(&self, _u: f64, n: i32) -> DVec2 {
        if n <= 0 {
            panic!("Standard_OutOfRange");
        }
        if n == 1 {
            return self.my_ax2d.direction;
        }
        DVec2::ZERO
    }

    /// OCCT GetType() (cxx L222-225).
    fn get_type(&self) -> Curve2dType {
        Curve2dType::Line
    }

    /// OCCT Resolution(R3d) (cxx L215-218) — returns R3d.
    fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::NbSamples — the line default (2 for
    /// a GeomAbs_Line source, clamped by the tool minimum of 20 as the
    /// landed Curve2d implementation).
    fn nb_samples(&self) -> i32 {
        20
    }

    /// OCCT IsClosed() (cxx L136-139).
    fn is_closed(&self) -> bool {
        false
    }

    /// OCCT IsPeriodic() (cxx L143-146).
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT Period() (cxx L150-153).
    fn period(&self) -> f64 {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT NbKnots() (cxx L285-288).
    fn nb_knots(&self) -> i32 {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT Degree() (cxx L264-267).
    fn degree(&self) -> i32 {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT NbPoles() (cxx L278-281).
    fn nb_poles(&self) -> i32 {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT Circle() (cxx L236-239).
    fn circle(&self) -> rcad_kernel::geom::Circle2d {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT Line() (cxx L229-232).
    fn line(&self) -> Line2d {
        self.my_ax2d
    }

    /// OCCT Ellipse() (cxx L243-246).
    fn ellipse(&self) -> rcad_kernel::geom::Ellipse2d {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT Parabola() (cxx L257-260).
    fn parabola(&self) -> rcad_kernel::geom::Parabola2d {
        panic!("Standard_NoSuchObject");
    }

    /// OCCT Hyperbola() (cxx L250-253).
    fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola2d {
        panic!("Standard_NoSuchObject");
    }
}

/// OCCT Continuity() (cxx L100-103) — GeomAbs_CN.
impl Line2dAdaptor {
    pub fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::CN
    }

    /// OCCT NbIntervals(S) (cxx L108-111).
    pub fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Intervals(T, S) (cxx L115-121).
    pub fn intervals(&self, t: &mut [f64], _s: GeomAbsShape) {
        t[0] = self.my_ufirst;
        t[1] = self.my_ulast;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: the restriction line Value at its parameters — a line
    /// through (0, 2) along X trimmed to [1, 4] evaluates to (x, 2).
    #[test]
    fn line2d_adaptor_value_anchor() {
        let l = Line2dAdaptor::new_pnt_dir(DVec2::new(0.0, 2.0), DVec2::X, 1.0, 4.0);
        assert_eq!(l.first_parameter(), 1.0);
        assert_eq!(l.last_parameter(), 4.0);
        let (p, v) = l.d1(3.0);
        assert_eq!(p, DVec2::new(3.0, 2.0));
        assert_eq!(v, DVec2::X);
        let t = l.trim(2.0, 3.5);
        assert_eq!(t.first_parameter(), 2.0);
        assert_eq!(t.value(2.5), DVec2::new(2.5, 2.0));
    }
}
