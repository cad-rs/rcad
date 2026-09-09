//! OCCT Approx_SameParameter (TKGeomBase/Approx — Approx_SameParameter.cxx,
//! 992 lines) — GAP carrier (plan §0.6) consumed by the BRepFill_Sweep
//! SameParameter static (part B, cxx L359-390).
//!
//! The real body (SameParameter::Compute — the curve-on-surface
//! approximation through Approx_CurveOnSurface and the tolerance walk) is
//! out of this dispatch; the OCCT failure path is preserved (IsDone() ==
//! false, IsSameParameter() == false — the OCCT "echec SameParameter" arm).

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};

/// OCCT Approx_SameParameter (Approx_SameParameter.hxx L40-96).
pub struct ApproxSameParameter {
    /// OCCT myDone.
    my_done: bool,
    /// OCCT mySameParameter (the tolerance was reached without a rebuild).
    my_same_parameter: bool,
    /// OCCT myTolReached.
    my_tol_reached: f64,
    /// OCCT myCurve2d (the rebuilt pcurve).
    my_curve2d: Option<Curve2d>,
    /// OCCT myCurve3d (the rebuilt 3d curve).
    my_curve3d: Option<Curve3>,
    /// OCCT myCurveOnSurface (the curve-on-surface of the rebuild).
    my_curve_on_surface: Option<(Curve2d, Surface3)>,
}

impl ApproxSameParameter {
    /// OCCT Approx_SameParameter(C3d, Pcurv, S, Tol3d)
    /// (Approx_SameParameter.cxx L44-58).
    pub fn new(
        _c3d: &Curve3,
        _c3d_first: f64,
        _c3d_last: f64,
        _pcurv: &Curve2d,
        _s: &Surface3,
        _tol3d: f64,
    ) -> Self {
        panic!(
            "GAP: Approx_SameParameter (TKGeomBase/Approx, \
             Approx_SameParameter.cxx L44-58 + Compute) is not translated — \
             see file header (plan section 0.6)"
        );
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsSameParameter().
    pub fn is_same_parameter(&self) -> bool {
        self.my_same_parameter
    }

    /// OCCT TolReached().
    pub fn tol_reached(&self) -> f64 {
        self.my_tol_reached
    }

    /// OCCT Curve2d().
    pub fn curve2d(&self) -> Curve2d {
        self.my_curve2d.clone().expect("Approx_SameParameter::Curve2d")
    }

    /// OCCT Curve3d().
    pub fn curve3d(&self) -> Curve3 {
        self.my_curve3d.clone().expect("Approx_SameParameter::Curve3d")
    }

    /// OCCT CurveOnSurface() — the rcad (pcurve, surface) value pair (the
    /// Adaptor3d_CurveOnSurface handle mapping).
    pub fn curve_on_surface(&self) -> (Curve2d, Surface3) {
        self.my_curve_on_surface
            .clone()
            .expect("Approx_SameParameter::CurveOnSurface")
    }
}
