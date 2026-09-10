//! OCCT GeomLib_CheckCurveOnSurface (TKGeomBase/GeomLib —
//! GeomLib_CheckCurveOnSurface.cxx, 774 lines) — GAP carrier (plan §0.6)
//! consumed by the BRepFill_Sweep CheckSameParameterExact static (part B,
//! cxx L294-297).
//!
//! The real body (CheckCurveOnSurface::Perform — the extrema search between
//! the 3d curve and the curve-on-surface through
//! math_FunctionSetRoot/GeomLib_CheckCurveOnSurface evalutors) is out of
//! this dispatch; the OCCT failure path is preserved.

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};

/// OCCT GeomLib_CheckCurveOnSurface (GeomLib_CheckCurveOnSurface.hxx
/// L29-110).
pub struct GeomLibCheckCurveOnSurface {
    /// OCCT myErrorStatus.
    my_error_status: i32,
    /// OCCT myMaxDistance.
    my_max_distance: f64,
    /// OCCT myMaxParameter.
    my_max_parameter: f64,
    /// OCCT myParallel.
    my_parallel: bool,
}

impl GeomLibCheckCurveOnSurface {
    /// OCCT GeomLib_CheckCurveOnSurface(Curve)
    /// (GeomLib_CheckCurveOnSurface.cxx L40-46).
    pub fn new(_curve: &Curve3) -> Self {
        GeomLibCheckCurveOnSurface {
            my_error_status: 0,
            my_max_distance: 0.0,
            my_max_parameter: 0.0,
            my_parallel: true,
        }
    }

    /// OCCT SetParallel(theParallel) (cxx L53-56).
    pub fn set_parallel(&mut self, the_parallel: bool) {
        self.my_parallel = the_parallel;
    }

    /// OCCT Perform(CurveOnSurface) (cxx L62-120) — the max-distance
    /// computation between the 3d curve and its representation.
    pub fn perform(&mut self, _curve_on_surface: &(Curve2d, Surface3)) {
        panic!(
            "GAP: GeomLib_CheckCurveOnSurface::Perform (TKGeomBase/GeomLib, \
             GeomLib_CheckCurveOnSurface.cxx L62-120) is not translated — \
             see file header (plan section 0.6)"
        );
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.my_error_status == 0
    }

    /// OCCT MaxDistance().
    pub fn max_distance(&self) -> f64 {
        self.my_max_distance
    }

    /// OCCT MaxParameter().
    pub fn max_parameter(&self) -> f64 {
        self.my_max_parameter
    }
}
