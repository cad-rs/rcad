//! OCCT Geom2dConvert_BSplineCurveToBezierCurve
//! (TKGeomBase/Geom2dConvert/Geom2dConvert_BSplineCurveToBezierCurve.cxx,
//! whole file L1-155).
//!
//! 1:1 translation on top of the faithful [`Geom2dBSplineCurve`] port.  The
//! OCCT `Arc` returns a `Geom2d_BezierCurve`; the rcad [`BezierCurve2`]
//! carries a weights array in all cases, so the non-rational branch stores
//! unit weights (architecture adaptation).

use crate::base::geom2d_convert::bspline_curve::Geom2dBSplineCurve;
use crate::core::precision::PCONFUSION;
use crate::geom::BezierCurve2;

/// OCCT Geom2dConvert_BSplineCurveToBezierCurve.
pub struct BSplineCurveToBezierCurve {
    my_curve: Geom2dBSplineCurve,
}

impl BSplineCurveToBezierCurve {
    /// OCCT constructor from the whole basis curve
    /// (Geom2dConvert_BSplineCurveToBezierCurve.cxx L27-44).
    pub fn new(basis_curve: &Geom2dBSplineCurve) -> Self {
        let mut my_curve = basis_curve.copy();
        // periodic curve can't be converted correctly by two main reasons:
        // last pole (equal to first one) is missing;
        // poles recomputation using default boor scheme is fails.
        if my_curve.is_periodic() {
            my_curve.set_not_periodic();
        }
        let uf = my_curve.first_parameter();
        let ul = my_curve.last_parameter();
        // OCCT Segment(Uf, Ul) — default theTolerance = Precision::PConfusion().
        my_curve.segment(uf, ul, PCONFUSION);
        my_curve.increase_multiplicity_range(
            my_curve.first_uknot_index(),
            my_curve.last_uknot_index(),
            my_curve.degree() as i32,
        );
        BSplineCurveToBezierCurve { my_curve }
    }

    /// OCCT constructor on the [U1, U2] range
    /// (Geom2dConvert_BSplineCurveToBezierCurve.cxx L48-91).
    pub fn new_with_range(
        basis_curve: &Geom2dBSplineCurve,
        u1: f64,
        u2: f64,
        parametric_tolerance: f64,
    ) -> Self {
        if u2 - u1 < parametric_tolerance {
            panic!("Standard_DomainError: GeomConvert_BSplineCurveToBezierSurface");
        }

        let mut uf = u1;
        let mut ul = u2;
        let ptol = parametric_tolerance / 2.0;

        let mut my_curve = basis_curve.copy();
        if my_curve.is_periodic() {
            my_curve.set_not_periodic();
        }

        let mut i1 = 0i32;
        let mut i2 = 0i32;
        my_curve.locate_u(u1, ptol, &mut i1, &mut i2, false);
        if i1 == i2 {
            // We are on the knot
            if my_curve.knot(i1) > u1 {
                uf = my_curve.knot(i1);
            }
        }

        my_curve.locate_u(u2, ptol, &mut i1, &mut i2, false);
        if i1 == i2 {
            // We are on the knot
            if my_curve.knot(i1) < u2 {
                ul = my_curve.knot(i1);
            }
        }

        my_curve.segment(uf, ul, PCONFUSION);
        my_curve.increase_multiplicity_range(
            my_curve.first_uknot_index(),
            my_curve.last_uknot_index(),
            my_curve.degree() as i32,
        );
        BSplineCurveToBezierCurve { my_curve }
    }

    /// OCCT Arc(Index) (Geom2dConvert_BSplineCurveToBezierCurve.cxx L95-125).
    pub fn arc(&self, index: i32) -> BezierCurve2 {
        if index < 1 || index > self.my_curve.nb_knots() - 1 {
            panic!("Standard_OutOfRange: Geom2dConvert_BSplineCurveToBezierCurve");
        }
        let deg = self.my_curve.degree() as i32;

        if self.my_curve.is_rational() {
            let mut poles = Vec::with_capacity((deg + 1) as usize);
            let mut weights = Vec::with_capacity((deg + 1) as usize);
            for i in 1..=deg + 1 {
                poles.push(self.my_curve.pole(i + deg * (index - 1)));
                weights.push(self.my_curve.weight(i + deg * (index - 1)));
            }
            BezierCurve2 {
                control_points: poles,
                weights,
            }
        } else {
            let mut poles = Vec::with_capacity((deg + 1) as usize);
            for i in 1..=deg + 1 {
                poles.push(self.my_curve.pole(i + deg * (index - 1)));
            }
            let n = poles.len();
            BezierCurve2 {
                control_points: poles,
                weights: vec![1.0; n],
            }
        }
    }

    /// OCCT Arcs(Curves) (Geom2dConvert_BSplineCurveToBezierCurve.cxx
    /// L129-137) — fills `curves[i - 1]` with `Arc(i)` for i = 1..NbArcs().
    pub fn arcs(&self, curves: &mut [BezierCurve2]) {
        let n = self.nb_arcs();
        for i in 1..=n {
            curves[(i - 1) as usize] = self.arc(i);
        }
    }

    /// OCCT Knots(TKnots) (Geom2dConvert_BSplineCurveToBezierCurve.cxx
    /// L141-148).
    pub fn knots(&self, tknots: &mut [f64]) {
        let mut kk = 0usize;
        for ii in 1..=self.my_curve.nb_knots() {
            tknots[kk] = self.my_curve.knot(ii);
            kk += 1;
        }
    }

    /// OCCT NbArcs (Geom2dConvert_BSplineCurveToBezierCurve.cxx L152-155).
    pub fn nb_arcs(&self) -> i32 {
        self.my_curve.nb_knots() - 1
    }
}
