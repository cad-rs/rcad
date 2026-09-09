//! OCCT Bisector — package function `IsConvex`, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector.hxx (L30-36) / Bisector.cxx (L24-33).

use std::sync::Arc;

use glam::DVec2;

use super::bisector_curve::BisectorCurve;

/// OCCT GeomAbs_JoinType (GeomAbs_JoinType.hxx).
///
/// GAP note: kernel-level enum; a private copy already exists in
/// `crate::brep_fill::offset_wire` — this local definition stands until the
/// enum is upstreamed to rcad-kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsJoinType {
    Arc,
    Intersection,
    Round,
}

/// OCCT Bisector::IsConvex(Cu, Sign) (Bisector.cxx L24-33).
pub fn is_convex(cu: &Arc<dyn BisectorCurve>, sign: f64) -> bool {
    let u1 = (cu.last_parameter() + cu.first_parameter()) / 2.;
    let (_p1, v1, v2) = cu.d2(u1);
    let tol = 1.0e-5;
    // OCCT: Sign * (V1 ^ V2) < Tol — the 2D cross product V1^V2.
    sign * cross2(v1, v2) < tol
}

/// OCCT gp_Vec2d::Crossed (gp_Vec2d.hxx) — the `^` 2D cross product.
#[inline]
pub(crate) fn cross2(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// OCCT gp_Vec2d/gp_Pnt2d::IsEqual(P, Tolerance).
#[inline]
pub(crate) fn pnt2d_equal(a: DVec2, b: DVec2, tolerance: f64) -> bool {
    a.distance(b) <= tolerance
}

/// OCCT gp_Dir2d::IsParallel(Other, AngularTolerance) — |angle| within the
/// tolerance of 0 or PI.
#[inline]
pub(crate) fn dir2d_parallel(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    let d = a.normalize_or_zero().dot(b.normalize_or_zero());
    let limit = (std::f64::consts::PI - angular_tolerance).cos();
    d >= limit || d <= -limit
}

/// OCCT gp_Vec2d::Rotate(PI/2) — rotation by +PI/2 (ccw).
#[inline]
pub(crate) fn rotate_half_pi(v: DVec2) -> DVec2 {
    DVec2::new(-v.y, v.x)
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    // OCCT Bisector::IsConvex — direct-API regression.
    #[test]
    fn is_convex_on_unit_circle() {
        use super::super::bisector_curve::Geom2dCurveHandle;
        use rcad_kernel::geom::{Circle2d, Curve2d};
        let c = Circle2d::new(glam::DVec2::ZERO, 1.0);
        let cu: Arc<dyn super::super::bisector_curve::BisectorCurve> =
            Arc::new(Geom2dCurveHandle::new(Curve2d::Circle(c)));
        // OCCT Bisector.cxx L33: Sign * (V1 ^ V2) < Tol.  At the mid
        // parameter of a ccw unit circle the cross product is +1, so
        // Sign=1 gives 1 < 1e-5 = false; Sign=-1 gives -1 < 1e-5 = true.
        assert!(!super::is_convex(&cu, 1.0));
        assert!(super::is_convex(&cu, -1.0));
    }
}
