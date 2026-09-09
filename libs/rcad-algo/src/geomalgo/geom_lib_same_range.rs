//! OCCT GeomLib statics (TKGeomBase/GeomLib — GeomLib.cxx) — GAP carriers
//! (plan §0.6) for the two routines consumed by BRepFill_Sweep part B:
//! - SameRange (GeomLib.cxx L842-908) — rebuilds a 2d curve over a target
//!   parameter range (consumed by Filling cxx L1037 and UpdateEdge
//!   cxx L1767);
//! - ExtendSurfByLength (GeomLib.cxx L1485+) — extends a bounded surface
//!   along an iso (consumed by BuildShell cxx L2334/L2341).
//!
//! The real bodies (SameRange's BSplCLib::ReparamMap / curve segmentation
//! and ExtendSurfByLength's prolongation approximation) are out of this
//! dispatch; the OCCT failure path is preserved.

use rcad_kernel::geom::{Curve2d, Surface3};

/// OCCT GeomLib::SameRange(Tolerance, Curve2d, First, Last, NewFirst,
/// NewLast, NewCurve2d) (GeomLib.cxx L842-908).
pub fn same_range(
    _tolerance: f64,
    curve2d: &Curve2d,
    _first: f64,
    _last: f64,
    _new_first: f64,
    _new_last: f64,
) -> Curve2d {
    let _ = curve2d;
    panic!(
        "GAP: GeomLib::SameRange (TKGeomBase/GeomLib.cxx L842-908) is not \
         translated — see file header (plan section 0.6)"
    );
}

/// OCCT GeomLib::ExtendSurfByLength(BoundedSurface, Length, Contour,
/// InU, After) (GeomLib.cxx L1485+).
pub fn extend_surf_by_length(
    _surface: &Surface3,
    _length: f64,
    _contour: i32,
    _in_u: bool,
    _after: bool,
) -> Surface3 {
    panic!(
        "GAP: GeomLib::ExtendSurfByLength (TKGeomBase/GeomLib.cxx L1485+) is \
         not translated — see file header (plan section 0.6)"
    );
}
