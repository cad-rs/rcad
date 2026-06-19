use rcad_kernel::geom::{Curve3, CurveEval};

/// Curve parameter step: parameter increment needed to move `tol` distance along curve.
///
/// OCCT: Adaptor3d_Curve::Resolution(tol) = tol / |dP/dt|
/// (BRepLib_1.cxx L61, IntTools_ShrunkRange.cxx L162)
///
/// For a curve parameterized by `t`, the step `resolution` satisfies:
/// `|P(t + resolution) - P(t)| ≈ tol` (first-order approximation using tangent speed).
///
/// When the tangent speed is nearly zero (singularity), the resolution is clamped to `tol`.
fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {
    let speed = curve.tangent_at(t).length();
    if speed < 1e-15 {
        tol
    } else {
        tol / speed
    }
}

/// Compute the shrunk (valid) range for a curve segment, excluding the tolerance
/// spheres around the endpoint vertices.
///
/// OCCT-aligned: IntTools_ShrunkRange::Perform() (IntTools_ShrunkRange.cxx L107-191)
///              + BRepLib::FindValidRange (BRepLib_1.cxx L173-258)
///
/// The shrunk range [t_start, t_end] is the portion of [t1, t2] where the curve point
/// is outside both vertex tolerance spheres. This is critical for:
/// - Edge-Face intersection: the shrunk range tells us where the edge is "truly away"
///   from its endpoint vertices, avoiding false intersections near vertices
///   (OCCT IntTools_ShrunkRange purpose).
/// - Micro-edge detection: if the shrunk range is empty, the PaveBlock is too short
///   and should be removed (OCCT `!IsSplittable()` check).
///
/// OCCT algorithm (IntTools_ShrunkRange::Perform):
/// 1. Guard: return None if (t2 - t1) < Precision::PConfusion() (L117-120)
/// 2. Get vertex tolerances; for each vertex:
///    aTolV = max(aTolV, aTolE) + Precision::Confusion() (L129-142)
/// 3. Call BRepLib::FindValidRange which:
///    a. Computes eps = max(curve.Resolution(aTolE) * 0.1, Epsilon(aMaxPar),
///                          Precision::PConfusion())
///    b. For each endpoint, calls findNearestValidPoint that steps along the curve
///       until outside the tolerance sphere, then binary-search refines.
///    c. Returns (theFirst, theLast) with theFirst < theLast.
/// 4. Guard: return None if (myTS2 - myTS1) < Precision::PConfusion() (L152-155)
/// 5. Compute edge length on shrunk range (L159-170)
/// 6. Guard: return None if length < Precision::Confusion() (L171-174)
/// 7. Set is_splittable if length > 2*aTolE + 2*Precision::Confusion() (L184-187)
///
/// This implementation provides the core shrunk range computation (steps 1-4).
/// It uses curve_resolution at the endpoint parameters as an O(1) approximation
/// of the full stepping + binary search approach, which is accurate when the curve
/// speed is nearly constant near the endpoints.
///
/// # Arguments
/// * `curve` - The 3D curve
/// * `t_range` - The full parameter range [t1, t2] (t1 < t2)
/// * `v1_tol` - Geometric tolerance at the start vertex
/// * `v2_tol` - Geometric tolerance at the end vertex
/// * `edge_tol` - Geometric tolerance of the edge
///
/// # Returns
/// * `Some([t_start, t_end])` - Valid shrunk range where t_start < t_end
/// * `None` - Micro-edge; the entire range is covered by tolerance spheres
pub fn shrunk_range(
    curve: &Curve3,
    t_range: [f64; 2],
    v1_tol: f64,
    v2_tol: f64,
    edge_tol: f64,
) -> Option<[f64; 2]> {
    let [t1, t2] = t_range;

    // Step 1: Check for degenerate range
    // OCCT IntTools_ShrunkRange L117-120:
    //   if (myT2 - myT1 < aPDTol) return;  // aPDTol = Precision::PConfusion() = 1e-12
    if (t2 - t1).abs() < 1e-12 {
        return None;
    }

    // Step 2: Adjust vertex tolerances
    // OCCT IntTools_ShrunkRange L129-142:
    //   if (aTolV < aTolE) aTolV = aTolE;
    //   aTolV += aDTol;  // aDTol = Precision::Confusion() = 1e-7
    let confusion = 1e-7; // Precision::Confusion()
    let a_tol_v1 = v1_tol.max(edge_tol) + confusion;
    let a_tol_v2 = v2_tol.max(edge_tol) + confusion;

    // Step 3: Compute curve resolution at each endpoint
    // OCCT BRepLib_1.cxx L61:
    //   aStep = theCurve.Resolution(theTol) * 1.01;
    // Use curve_resolution (tol / tangent_speed) at the endpoint parameters.
    // This approximates the minimal parameter step needed to move a_tol_V distance
    // along the curve at the vertex location.
    let res1 = curve_resolution(curve, t1, a_tol_v1);
    let res2 = curve_resolution(curve, t2, a_tol_v2);

    // Step 4: Compute shrunk range [ts1, ts2]
    // OCCT BRepLib::FindValidRange returns theFirst/theLast via stepping + bisection.
    // Here we use the first-order approximation:
    //   ts1 = t1 + resolution_at_t1
    //   ts2 = t2 - resolution_at_t2
    let ts1 = t1 + res1;
    let ts2 = t2 - res2;

    // Step 5: Validate the shrunk range
    // OCCT BRepLib::FindValidRange L251-254:
    //   if (theFirst > theLast) → not valid (overlapping tolerance spheres)
    // OCCT IntTools_ShrunkRange L152-155:
    //   if ((myTS2 - myTS1) < aPDTol) → micro edge
    if ts1 >= ts2 || (ts2 - ts1) < 1e-12 {
        return None;
    }

    Some([ts1, ts2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::Line3;

    // ── Helper: line curve with unit-speed parameterization ──

    /// A line from (0,0,0) to (len,0,0), direction = X (unit length).
    /// Parameter t ∈ [0, len] maps to P(t) = (t, 0, 0).
    fn unit_line(len: f64) -> Curve3 {
        Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        })
    }

    // ── Helper: circle in XY plane, radius R, centered at origin ──
    fn unit_circle() -> Curve3 {
        Curve3::Circle(rcad_kernel::geom::Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        })
    }

    // ── shrunk_range: line tests ──

    #[test]
    fn long_line_has_valid_shrunk_range() {
        // For a line of length 100 with unit speed:
        // a_tol_V = max(1e-6, 1e-7) + 1e-7 = 1.1e-6
        // res = 1.1e-6 / 1.0 = 1.1e-6
        // ts1 ≈ 1.1e-6, ts2 ≈ 100 - 1.1e-6
        let line = unit_line(100.0);
        let result = shrunk_range(&line, [0.0, 100.0], 1e-6, 1e-6, 1e-7);
        assert!(result.is_some(), "Long line should have a valid shrunk range");
        let [ts1, ts2] = result.unwrap();
        assert!(ts1 > 0.0, "Start should be shrunk inward, got {}", ts1);
        assert!(ts2 < 100.0, "End should be shrunk inward, got {}", ts2);
        assert!(
            (ts2 - ts1 - 100.0).abs() < 1e-5,
            "Shrunk length should be ~100, got {}",
            ts2 - ts1
        );
    }

    #[test]
    fn large_tolerance_collapses_range() {
        // Tolerance sphere covers entire edge: a_tol = max(50, 1e-7) + 1e-7 ≈ 50
        // res ≈ 50 / 1.0 = 50, so ts1 = 50, ts2 = 50 - 50 = 0 → None
        let line = unit_line(100.0);
        let result = shrunk_range(&line, [0.0, 100.0], 50.0, 50.0, 1e-7);
        assert!(result.is_none(), "Large tolerance should collapse the range");
    }

    #[test]
    fn micro_edge_returns_none() {
        // t_range is below PConfusion threshold (1e-12)
        let line = unit_line(1.0);
        let result = shrunk_range(&line, [0.0, 1e-13], 1e-7, 1e-7, 1e-7);
        assert!(result.is_none(), "Micro-edge should return None");
    }

    #[test]
    fn zero_length_range_returns_none() {
        let line = unit_line(1.0);
        let result = shrunk_range(&line, [0.5, 0.5], 1e-7, 1e-7, 1e-7);
        assert!(result.is_none(), "Zero-length range should return None");
    }

    #[test]
    fn slightly_too_short_edge_returns_none() {
        // Edge of length 1e-9 with small tolerances.
        // a_tol_V = max(1e-7, 1e-7) + 1e-7 = 2e-7
        // res = 2e-7 / 1.0 = 2e-7
        // ts1 = t1 + 2e-7, ts2 = t2 - 2e-7 = t1 + 1e-9 - 2e-7 < ts1 → None
        let line = unit_line(1.0);
        let result = shrunk_range(&line, [0.0, 1e-9], 1e-7, 1e-7, 1e-7);
        assert!(result.is_none(), "Very short edge should be a micro-edge");
    }

    // ── shrunk_range: circle tests ──

    #[test]
    fn circle_full_turn_has_valid_shrunk_range() {
        // Circle radius 1, full turn [0, 2*PI].
        // tangent_at is normalized → speed = 1.0
        // a_tol_V = max(1e-6, 1e-7) + 1e-7 = 1.1e-6
        // res = 1.1e-6 / 1.0 = 1.1e-6
        let circle = unit_circle();
        let result = shrunk_range(&circle, [0.0, 2.0 * std::f64::consts::PI], 1e-6, 1e-6, 1e-7);
        assert!(result.is_some(), "Full circle should have a valid shrunk range");
        let [ts1, ts2] = result.unwrap();
        assert!(ts1 > 0.0, "Circle start should be shrunk inward");
        assert!(ts2 < 2.0 * std::f64::consts::PI, "Circle end should be shrunk inward");
        assert!(
            ts2 - ts1 > 6.0,
            "Most of circle should remain in shrunk range"
        );
    }

    #[test]
    fn circle_tiny_arc_returns_none() {
        // Tiny arc that is shorter than the tolerance allowance
        let circle = unit_circle();
        let result = shrunk_range(&circle, [0.0, 1e-8], 1e-7, 1e-7, 1e-7);
        assert!(result.is_none(), "Tiny arc should be a micro-edge");
    }

    // ── shrunk_range: asymmetric tolerances ──

    #[test]
    fn asymmetric_tolerances() {
        // Start vertex has large tolerance, end has small tolerance.
        // Start side: a_tol_V1 = max(10.0, 1e-7) + 1e-7 ≈ 10.0000001
        //   res1 = 10.0000001 / 1.0 = 10.0000001, ts1 ≈ 10.0000001
        // End side: a_tol_V2 = max(1e-7, 1e-7) + 1e-7 = 2e-7
        //   res2 = 2e-7 / 1.0 = 2e-7, ts2 ≈ 100 - 2e-7 ≈ 99.9999998
        let line = unit_line(100.0);
        let result = shrunk_range(&line, [0.0, 100.0], 10.0, 1e-7, 1e-7);
        assert!(result.is_some(), "Asymmetric tolerances should still produce a range");
        let [ts1, ts2] = result.unwrap();
        let expected_ts1 = 10.0 + 1e-7; // max(10.0, 1e-7) + 1e-7 = 10.0 + 1e-7 (the extra Confusion)
        assert!(
            (ts1 - expected_ts1).abs() < 1e-10,
            "Start should shrink by ~{:.10}, got {:.10}",
            expected_ts1,
            ts1
        );
        assert!(
            (ts2 - 100.0).abs() < 1e-3,
            "End should barely shrink, got {}",
            ts2
        );
    }

    #[test]
    fn asymmetric_tolerances_collapse() {
        // Large start tolerance covers the entire range
        let line = unit_line(1.0);
        let result = shrunk_range(&line, [0.0, 1.0], 0.6, 0.6, 1e-7);
        // a_tol_V = max(0.6, 1e-7) + 1e-7 ≈ 0.6000001
        // res = 0.6000001 / 1.0 = 0.6000001
        // ts1 ≈ 0.6000001, ts2 ≈ 1.0 - 0.6000001 ≈ 0.3999999
        // ts1 > ts2 → None
        assert!(result.is_none(), "Large start tol should collapse range");
    }

    // ── curve_resolution tests ──

    #[test]
    fn resolution_unit_speed_line() {
        let line = unit_line(1.0);
        // Speed = direction.length() = 1.0, so resolution = tol / 1.0 = tol
        let res = curve_resolution(&line, 0.5, 1e-7);
        assert!((res - 1e-7).abs() < 1e-20, "Resolution should equal tol for unit-speed line");
    }

    #[test]
    fn resolution_half_speed_line() {
        // Direction length = 0.5, so speed = 0.5
        let line = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X * 0.5,
        });
        // resolution = tol / 0.5 = 2 * tol
        let res = curve_resolution(&line, 0.5, 1e-7);
        assert!((res - 2e-7).abs() < 1e-20, "Resolution should be tol/speed for half-speed line");
    }
}
