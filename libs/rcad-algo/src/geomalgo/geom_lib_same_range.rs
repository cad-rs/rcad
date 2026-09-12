//! OCCT GeomLib statics (TKGeomBase/GeomLib — GeomLib.cxx):
//! - SameRange (GeomLib.cxx L842-970) — rebuilds a 2d curve over a target
//!   parameter range (consumed by Filling cxx L1037 and UpdateEdge
//!   cxx L1767).  The 1:1 body lives in the kernel as
//!   `rcad_kernel::geom::same_range_2d`; this module is the OCCT-signature
//!   entry point (Tolerance first, `NewCurvePtr` out-parameter) used by the
//!   BRepFill/offset callers, so a single implementation serves both.
//! - ExtendSurfByLength (GeomLib.cxx L1485+) — extends a bounded surface
//!   along an iso (consumed by BuildShell cxx L2334/L2341).

use rcad_kernel::geom::{Curve2d, Surface3};

/// OCCT GeomLib::SameRange(Tolerance, Curve2d, First, Last, NewFirst,
/// NewLast, NewCurve2d) (GeomLib.cxx L842-970).
///
/// Returns the OCCT `NewCurvePtr` result; the OCCT null-handle outcome (the
/// L902-922 guard leaving `NewCurvePtr` untouched, or a failed
/// `CurveToBSplineCurve`) is the failure the callers observe as a null
/// handle.
pub fn same_range(
    tolerance: f64,
    curve2d: &Curve2d,
    first: f64,
    last: f64,
    new_first: f64,
    new_last: f64,
) -> Curve2d {
    match rcad_kernel::geom::same_range_2d(
        tolerance,
        curve2d.clone(),
        first,
        last,
        new_first,
        new_last,
    ) {
        Some(c) => c,
        None => panic!(
            "Standard_NullObject: GeomLib::SameRange produced a null curve \
             (GeomLib.cxx L842-970)"
        ),
    }
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
