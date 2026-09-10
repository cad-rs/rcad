// OCCT HLRBRep_LineTool (TKHLR) — the static tool over a gp_Lin (the
// projected-contour support line consumed by the walking/visibility
// algorithms).
//
// HLRBRep_LineTool.hxx L47-177 + .lxx (the one-line bodies; RealFirst/
// RealLast bounds, CN continuity, ElCLib evaluations).  The rcad line is
// the kernel `Line3`; the Bezier/BSpline/NbPoles-style members of the hxx
// are Standard_NoSuchObject for a line and stay unrepresentable here (the
// gp_Lin type itself answers them).

use rcad_kernel::geom::{Line3, Point3, Vec3};

/// OCCT FirstParameter (lxx) — RealFirst().
pub fn first_parameter(_c: &Line3) -> f64 {
    f64::MIN
}

/// OCCT LastParameter (lxx) — RealLast().
pub fn last_parameter(_c: &Line3) -> f64 {
    f64::MAX
}

/// OCCT Continuity (lxx) — GeomAbs_CN (4).
pub fn continuity(_c: &Line3) -> i32 {
    4
}

/// OCCT NbIntervals (lxx) — 1.
pub fn nb_intervals(_c: &Line3, _s: i32) -> i32 {
    1
}

/// OCCT IntervalFirst (lxx).
pub fn interval_first(_c: &Line3) -> f64 {
    f64::MIN
}

/// OCCT IntervalLast (lxx).
pub fn interval_last(_c: &Line3) -> f64 {
    f64::MAX
}

/// OCCT IntervalContinuity (lxx) — GeomAbs_CN (4).
pub fn interval_continuity(_c: &Line3) -> i32 {
    4
}

/// OCCT IsClosed (lxx).
pub fn is_closed(_c: &Line3) -> bool {
    false
}

/// OCCT IsPeriodic (lxx).
pub fn is_periodic(_c: &Line3) -> bool {
    false
}

/// OCCT Period (lxx).
pub fn period(_c: &Line3) -> f64 {
    0.0
}

/// OCCT Value(C, U) — ElCLib::Value(U, C).
pub fn value(c: &Line3, u: f64) -> Point3 {
    rcad_kernel::math::el::elclib_line_value(u, c.origin, c.direction)
}

/// OCCT D0(C, U, P).
pub fn d0(c: &Line3, u: f64) -> Point3 {
    value(c, u)
}

/// OCCT D1(C, U, P, V) — ElCLib::D1.
pub fn d1(c: &Line3, u: f64) -> (Point3, Vec3) {
    rcad_kernel::math::el::elclib_line_d1(u, c.origin, c.direction)
}

/// OCCT D2(C, U, P, V1, V2) — ElCLib::D1 + V2 = 0 (lxx).
pub fn d2(c: &Line3, u: f64) -> (Point3, Vec3, Vec3) {
    let (p, v1) = d1(c, u);
    (p, v1, Vec3::ZERO)
}

/// OCCT D3(C, U, P, V1, V2, V3) — ElCLib::D1 + V2 = V3 = 0 (lxx).
pub fn d3(c: &Line3, u: f64) -> (Point3, Vec3, Vec3, Vec3) {
    let (p, v1) = d1(c, u);
    (p, v1, Vec3::ZERO, Vec3::ZERO)
}

/// OCCT DN(C, U, N) — ElCLib::DN: the direction for N == 1, zero beyond.
pub fn dn(c: &Line3, u: f64, n: i32) -> Vec3 {
    let _ = u;
    if n == 1 {
        c.direction
    } else {
        Vec3::ZERO
    }
}

/// OCCT Resolution(C, R3D) (lxx) — R3D.
pub fn resolution(_c: &Line3, r3d: f64) -> f64 {
    r3d
}

/// OCCT GetType — GeomAbs_Line.
pub fn get_type(_c: &Line3) -> rcad_kernel::base::proj_lib::CurveType {
    rcad_kernel::base::proj_lib::CurveType::Line
}

/// OCCT Line(C) — identity.
pub fn line(c: &Line3) -> Line3 {
    c.clone()
}

/// OCCT NbSamples(C, U0, U1) (lxx) — 2 for a line.
pub fn nb_samples(_c: &Line3, _u0: f64, _u1: f64) -> i32 {
    2
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: the line statics — infinite bounds, CN, ElCLib
    /// evaluations, and the zero higher derivatives (lxx bodies).
    #[test]
    fn line_tool_statics() {
        let l = Line3::new(Point3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(first_parameter(&l), f64::MIN);
        assert_eq!(last_parameter(&l), f64::MAX);
        assert_eq!(continuity(&l), 4);
        assert!(!is_closed(&l) && !is_periodic(&l));

        let p = value(&l, 5.0);
        assert!((p.x - 1.0).abs() < 1e-12 && (p.y - 2.0).abs() < 1e-12 && (p.z - 8.0).abs() < 1e-12);

        let (p1, v1) = d1(&l, 0.0);
        assert!((v1.z - 1.0).abs() < 1e-12);
        let (_, _, v2) = d2(&l, 0.0);
        assert_eq!(v2, Vec3::ZERO);
        let (_, _, _, v3) = d3(&l, 0.0);
        assert_eq!(v3, Vec3::ZERO);
        assert_eq!(dn(&l, 0.0, 1), l.direction);
        assert_eq!(dn(&l, 0.0, 2), Vec3::ZERO);
    }
}
