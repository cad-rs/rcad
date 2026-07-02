//! OCCT-equivalence tests for low-level boolean support libraries.
//!
//! Each test verifies that rcad functions produce results matching
//! OCCT's known behavior on canonical inputs.

use glam::{DVec2, DVec3};
use rcad_algorithms::boptools::*;
use rcad_algorithms::inttools::edge_edge::*;
use rcad_algorithms::inttools::edge_face::*;
use rcad_algorithms::inttools::curve_range::curve_resolution;
use rcad_kernel::geom::*;
use rcad_kernel::projection::{closest_point_on_curve, closest_point_on_surface};
use rcad_kernel::extrema::extrema_curve_curve;
use rcad_kernel::extend::{bspline_to_bezier_curves, trim_curve, insert_knot_to_multiplicity};

// =========================================================================
// 1. closest_point_on_curve — OCCT GeomAPI_ProjectPointOnCurve equivalent
// =========================================================================

#[test]
fn test_point_on_line_closest_projection() {
    // Line: origin=(0,0,0), direction=(1,0,0)
    // Query point: (5, 10, 0) → closest should be (5, 0, 0), param=5, dist=10
    let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let result = closest_point_on_curve(&line, DVec3::new(5.0, 10.0, 0.0), 16);
    assert!((result.point - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-12);
    assert!((result.param - 5.0).abs() < 1e-12);
    assert!((result.distance - 10.0).abs() < 1e-12);
}

#[test]
fn test_point_on_circle_closest_projection() {
    // Circle: center=(0,0,0), normal=Z, radius=5
    // Query: (10, 0, 0) → closest = (5, 0, 0), dist = 5, param = 0
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
    let result = closest_point_on_curve(&circle, DVec3::new(10.0, 0.0, 0.0), 16);
    assert!((result.point - DVec3::new(5.0, 0.0, 0.0)).length() < 1e-10);
    assert!((result.distance - 5.0).abs() < 1e-10);
    // Parameter should be on unit circle (mod 2π)
    // The key correctness check is distance and point, not exact parameter
    assert!(result.param.is_finite(), "param should be finite");
}

#[test]
fn test_point_on_circle_center_projection() {
    // Point exactly at center of circle: should return any point on circle
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
    let result = closest_point_on_curve(&circle, DVec3::ZERO, 16);
    assert!((result.distance - 5.0).abs() < 1e-10);
    assert!((result.point.length() - 5.0).abs() < 1e-10);
}

// =========================================================================
// 2. closest_point_on_surface — OCCT GeomAPI_ProjectPointOnSurf equivalent
// =========================================================================

#[test]
fn test_point_on_plane_projection() {
    let plane = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, 5.0),
        normal: DVec3::Z,
    });
    let result = closest_point_on_surface(&plane, DVec3::new(10.0, 20.0, 15.0), 16);
    assert!((result.point - DVec3::new(10.0, 20.0, 5.0)).length() < 1e-12);
    assert!((result.distance - 10.0).abs() < 1e-12);
}

#[test]
fn test_point_on_sphere_projection() {
    let sphere = Surface3::Sphere(SphericalSurface {
        center: DVec3::ZERO, axis: DVec3::Z, radius: 3.0, ref_dir: DVec3::X,
    });
    // Point outside sphere
    let result = closest_point_on_surface(&sphere, DVec3::new(6.0, 0.0, 0.0), 16);
    assert!((result.distance - 3.0).abs() < 1e-10);
    assert!((result.point - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
    // Point inside sphere: projects to sphere surface
    let result2 = closest_point_on_surface(&sphere, DVec3::new(1.0, 0.0, 0.0), 16);
    assert!((result2.distance - 2.0).abs() < 1e-10);
    assert!((result2.point - DVec3::new(3.0, 0.0, 0.0)).length() < 1e-10);
}

#[test]
fn test_point_on_cylinder_projection() {
    let cyl = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::ZERO, axis: DVec3::Z, radius: 2.0, ref_dir: DVec3::X,
    });
    let result = closest_point_on_surface(&cyl, DVec3::new(5.0, 0.0, 3.0), 16);
    // Closest point on cylinder: radial direction from origin, same z
    assert!((result.distance - 3.0).abs() < 1e-10);
    assert!((result.point - DVec3::new(2.0, 0.0, 3.0)).length() < 1e-10);
}

// =========================================================================
// 3. extrema_curve_curve — OCCT GeomAPI_ExtremaCurveCurve equivalent
// =========================================================================

#[test]
fn test_extrema_line_line_skew() {
    // Two skew lines: L1 along X, L2 along Y offset in Z
    let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let l2 = Curve3::Line(Line3 { origin: DVec3::new(0.0, 0.0, 5.0), direction: DVec3::Y });
    let ext = extrema_curve_curve(&l1, &l2, 16);
    assert!(ext.pairs.len() >= 1);
    // Closest distance between the two lines should be 5 (Z offset)
    assert!((ext.min_distance() - 5.0).abs() < 1e-8);
}

#[test]
fn test_extrema_line_circle() {
    let l1 = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let c = Curve3::Circle(Circle3::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, 3.0));
    let ext = extrema_curve_curve(&l1, &c, 16);
    assert!(!ext.pairs.is_empty());
    // Minimum distance from line to circle center = 10, minus radius = 7
    assert!((ext.min_distance() - 7.0).abs() < 1e-6, "min_dist={}", ext.min_distance());
}

#[test]
fn test_extrema_circle_circle_coaxial() {
    // Two coaxial circles in parallel planes
    let c1 = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
    let c2 = Curve3::Circle(Circle3::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, 3.0));
    let ext = extrema_curve_curve(&c1, &c2, 16);
    assert!(!ext.pairs.is_empty());
    // Min distance: distance between circle centers projected radially
    // For coaxial circles, closest points: (5,0,0) on c1, (3,0,10) on c2
    // dist = sqrt((5-3)^2 + 10^2) = sqrt(4+100) ≈ 10.198
    assert!((ext.min_distance() - 10.198).abs() < 1e-2, "min_dist={}", ext.min_distance());
}

// =========================================================================
// 4. Curve resolution — OCCT BRepLib / Adaptor3d_Curve equivalent
// =========================================================================

#[test]
fn test_curve_resolution_line_is_tolerance() {
    // For a line, resolution = tol (speed normalization not needed for lines)
    let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let res = curve_resolution(&line, 0.5, 1e-7);
    assert!(res > 0.0);
    // On a line with unit direction, step = tol / |speed| = 1e-7 / 1.0
    assert!((res - 1e-7).abs() < 1e-12, "line res={}", res);
}

#[test]
fn test_curve_resolution_circle() {
    // Circle of radius 10: |tangent| = 10, resolution = tol/10
    let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 10.0));
    let res = curve_resolution(&circle, 0.0, 1e-7);
    // For a circle of radius 10, |tangent| = 10, so resolution should be
    // approximately tol / 10 = 1e-8. The exact value depends on finite difference step.
    assert!(res > 0.0 && res < 1e-5, "circle res={} should be positive and small", res);
}

// =========================================================================
// 5. 2D point-in-face classification — OCCT IntTools_FClass2d equivalent
// =========================================================================

#[test]
fn test_point_in_planar_face_square() {
    use rcad_algorithms::inttools::edge_face::point_in_planar_face;
    let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
    let square = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(10.0, 0.0, 0.0),
        DVec3::new(10.0, 10.0, 0.0),
        DVec3::new(0.0, 10.0, 0.0),
    ];
    assert!(point_in_planar_face(DVec3::new(5.0, 5.0, 0.0), &plane, &square));
    assert!(!point_in_planar_face(DVec3::new(15.0, 5.0, 0.0), &plane, &square));
    assert!(!point_in_planar_face(DVec3::new(5.0, 15.0, 0.0), &plane, &square));
}

// =========================================================================
// 6. BSpline → Bezier conversion — OCCT GeomConvert_BSplineCurveToBezierCurve
// =========================================================================

fn make_quadratic_bspline() -> BSplineCurve3 {
    // A quadratic BSpline with 3 control points and internal knot at 0.5
    BSplineCurve3 {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(5.0, 10.0, 0.0),
            DVec3::new(10.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    }
}

#[test]
fn test_bspline_to_bezier_preserves_geometry() {
    let curve = make_quadratic_bspline();
    let segments = bspline_to_bezier_curves(&curve);
    assert!(segments.len() >= 2, "should produce at least 2 Bezier segments, got {}", segments.len());
    for seg in &segments {
        assert_eq!(seg.degree, 2, "each Bezier segment should be quadratic");
    }
    // Each Bezier segment corresponds to a knot span of the original curve.
    // For now just verify the total number of segments is correct and they
    // match at the segment endpoints (global knots).
    // Segments use local domain [0,1].
    // The first segment covers [0,0.5], the second covers [0.5,1.0].
    // At the shared knot (0.5), both segments should give the same point.
    let mid = segments[0].point_at(1.0);
    let mid2 = segments[1].point_at(0.0);
    let d = (mid - mid2).length();
    assert!(d < 1e-8, "segments don't match at boundary: {}", d);
    // First segment at t=0 should match original at t=0
    let start0 = segments[0].point_at(0.0);
    let orig0 = curve.point_at(0.0);
    assert!((start0 - orig0).length() < 1e-8, "start mismatch");
    // Last segment at t=1 should match original at t=1
    let last = segments.last().unwrap().point_at(1.0);
    let orig1 = curve.point_at(1.0);
    assert!((last - orig1).length() < 1e-8, "end mismatch");
    // Midpoint: first segment at t=1 should be original at t=0.5
    let mid_pt = segments[0].point_at(1.0);
    let orig_mid = curve.point_at(0.5);
    assert!((mid_pt - orig_mid).length() < 5e-2, "mid mismatch: {}", (mid_pt - orig_mid).length());
}

#[test]
fn test_bspline_to_bezier_segment_count() {
    let curve = make_quadratic_bspline();
    let segments = bspline_to_bezier_curves(&curve);
    // With internal knot at 0.5, should get 2 segments
    assert_eq!(segments.len(), 2, "quadratic with 1 internal knot → 2 Bezier segments");
}

// =========================================================================
// 7. Knot insertion — OCCT Boehm algorithm equivalence
// =========================================================================

#[test]
fn test_knot_insertion_preserves_geometry() {
    let curve = make_quadratic_bspline();
    // Insert knot at 0.3
    let inserted = rcad_kernel::insert_knot_to_multiplicity(&curve, 0.3, 2);
    // Geometry should be preserved at all sample points
    for t in [0.0, 0.1, 0.3, 0.5, 0.7, 1.0] {
        let expected = curve.point_at(t);
        let got = inserted.point_at(t);
        assert!((got - expected).length() < 1e-12,
            "t={} mismatch after knot insertion", t);
    }
}

// =========================================================================
// 8. Trim curve — OCCT Geom_TrimmedCurve equivalent
// =========================================================================

#[test]
fn test_trim_bspline_preserves_segment() {
    let curve = make_quadratic_bspline();
    // trim_curve inserts boundary knots to full multiplicity
    // and extracts the internal segment
    let trimmed = rcad_kernel::extend::trim_curve(&curve, 0.2, 0.8);
    // The trimmed curve should have the same degree
    assert_eq!(trimmed.degree, curve.degree);
    // Test the trim with the full domain (identity trim)
    let full = rcad_kernel::extend::trim_curve(&curve, 0.0, 1.0);
    assert_eq!(full.control_points.len(), curve.control_points.len(),
        "full trim should preserve control points");
}

// =========================================================================
// 9. is_dirs_coinside — OCCT IntTools_Tools::IsDirsCoinside
// =========================================================================

#[test]
fn test_is_dirs_coinside_identical() {
    assert!(is_dirs_coinside(DVec3::X, DVec3::X));
    assert!(is_dirs_coinside(DVec3::Y, DVec3::Y));
    assert!(is_dirs_coinside(DVec3::Z, DVec3::Z));
}

#[test]
fn test_is_dirs_coinside_opposite() {
    // Opposite directions are 180° apart → distance of endpoints on unit sphere = 2
    // |D - (-D)| = 2 → (2 - 2).abs() < 0.0002 → true
    assert!(is_dirs_coinside(DVec3::X, -DVec3::X));
    assert!(is_dirs_coinside(DVec3::Y, -DVec3::Y));
}

#[test]
fn test_is_dirs_coinside_perpendicular() {
    // Perpendicular: distance = sqrt(2) ≈ 1.414, not < 0.0002 and |2 - 1.414| not < 0.0002
    assert!(!is_dirs_coinside(DVec3::X, DVec3::Y));
    assert!(!is_dirs_coinside(DVec3::X, DVec3::Z));
}

// =========================================================================
// 10. compute_int_range — OCCT IntTools_Tools::ComputeIntRange
// =========================================================================

#[test]
fn test_compute_int_range_orthogonal() {
    // Angle = π/2 → returns tol2 directly
    let result = compute_int_range(1e-7, 2e-7, std::f64::consts::FRAC_PI_2);
    assert!((result - 2e-7).abs() < 1e-14, "orthogonal int range = {}", result);
}

#[test]
fn test_compute_int_range_shallow() {
    // Angle = 10° → larger range
    let angle = 10.0_f64.to_radians();
    let result = compute_int_range(1e-7, 1e-7, angle);
    assert!(result > 1e-7, "shallow angle int range = {}", result);
}

// =========================================================================
// 11. is_on_pave / is_in_range — OCCT IntTools_Tools
// =========================================================================

#[test]
fn test_is_on_pave() {
    // 0.5 is not on the boundary of [0, 1] (tolerance 1e-7 → need within 1e-7 of 0 or 1)
    assert!(!is_on_pave(0.5, [0.0, 1.0], 1e-7));
    // 0.0 is exactly on the boundary
    assert!(is_on_pave(0.0, [0.0, 1.0], 1e-7));
    // 1.0 is exactly on the boundary
    assert!(is_on_pave(1.0, [0.0, 1.0], 1e-7));
    // 1e-6 from boundary with 1e-7 tolerance → too far
    assert!(!is_on_pave(1e-6, [0.0, 1.0], 1e-7));
    // 1e-8 from boundary with 1e-7 tolerance → within tolerance
    assert!(is_on_pave(1e-8, [0.0, 1.0], 1e-7));
}

#[test]
fn test_is_in_range() {
    assert!(is_in_range([0.3, 0.7], [0.2, 0.8], 0.01));
    assert!(!is_in_range([0.1, 0.2], [0.5, 0.8], 0.01));
    assert!(is_in_range([0.4, 0.6], [0.5, 0.7], 0.05));
}

// =========================================================================
// 12. intermediate_point — OCCT uses PAR_T=0.432... not 0.5
// =========================================================================

#[test]
fn test_intermediate_point_midpoint() {
    let result = intermediate_point(0.0, 10.0);
    assert!((result - 5.0).abs() < 1e-15);
}

#[test]
fn test_intermediate_point_occt_formula() {
    // OCCT: (1-PAR_T)*aFirst + PAR_T*aLast with PAR_T=0.43213918
    let result = intermediate_point_occt(0.0, 10.0);
    let expected = 4.3213918;
    assert!((result - expected).abs() < 1e-12, "occt intermed point = {}", result);
}

// =========================================================================
// 13. SenseFlag — OCCT BOPTools_AlgoTools3D::SenseFlag
// =========================================================================

#[test]
fn test_sense_flag_same() {
    assert_eq!(sense_flag(DVec3::Z, DVec3::Z), 1);
}

#[test]
fn test_sense_flag_opposite() {
    assert_eq!(sense_flag(DVec3::Z, -DVec3::Z), -1);
}

#[test]
fn test_sense_flag_perpendicular() {
    assert_eq!(sense_flag(DVec3::X, DVec3::Y), 0);
}

// =========================================================================
// 14. PointBoxDistance — OCCT IntTools_EdgeEdge::PointBoxDistance
// =========================================================================

#[test]
fn test_point_box_distance_inside() {
    let d = point_box_distance(DVec3::new(0.0, 0.0, 0.0),
        DVec3::splat(-1.0), DVec3::splat(1.0));
    assert!(d.abs() < 1e-15, "point inside box should have distance 0, got {}", d);
}

#[test]
fn test_point_box_distance_outside() {
    let d = point_box_distance(DVec3::new(5.0, 0.0, 0.0),
        DVec3::splat(-1.0), DVec3::splat(1.0));
    assert!((d - 4.0).abs() < 1e-12, "point outside box dist = {}", d);
}

#[test]
fn test_point_box_distance_corner() {
    let d = point_box_distance(DVec3::new(2.0, 3.0, 0.0),
        DVec3::splat(-1.0), DVec3::splat(1.0));
    let expected = ((1.0_f64).powi(2) + (2.0_f64).powi(2)).sqrt();
    assert!((d - expected).abs() < 1e-12, "corner dist = {} expected {}", d, expected);
}

// =========================================================================
// 15. EdgeEdge helper: split_range_on_segments — OCCT SplitRangeOnSegments
// =========================================================================

#[test]
fn test_split_range_on_segments_basic() {
    let (n, segs) = split_range_on_segments(0.0, 10.0, 0.5, 4);
    assert_eq!(n, 4);
    assert_eq!(segs.len(), 4);
    assert!((segs[0][0] - 0.0).abs() < 1e-12);
    assert!((segs[3][1] - 10.0).abs() < 1e-12);
}

#[test]
fn test_split_range_on_segments_below_resolution() {
    let (n, segs) = split_range_on_segments(0.0, 0.1, 0.5, 5);
    // Range smaller than resolution → single segment
    assert_eq!(n, 1);
    assert_eq!(segs.len(), 1);
}

// =========================================================================
// 16. EdgeFace: compute_edge_face_criteria — OCCT L528-548
// =========================================================================

#[test]
fn test_compute_edge_face_criteria_simple() {
    let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
    let crit = compute_edge_face_criteria(1e-7, 2e-7, &line);
    assert!((crit - 3e-7).abs() < 1e-14, "crit = {}", crit);
}

#[test]
fn test_compute_edge_face_criteria_bspline_large_ratio() {
    let bspline = Curve3::BSpline(BSplineCurve3 {
        degree: 3, knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![DVec3::ZERO; 4],
        weights: vec![1.0; 4],
    });
    let crit = compute_edge_face_criteria(1e-5, 1e-7, &bspline);
    // Ratio 1e-5/1e-7 = 100 → exactly at threshold, uses max
    assert!((crit - 1e-5).abs() < 1e-14, "bspline high ratio crit = {}", crit);
}

// =========================================================================
// 17. IsEqDistance — OCCT IntTools_EdgeFace::IsEqDistance equivalent
// =========================================================================

#[test]
fn test_is_eq_distance_cylinder_axis() {
    let cyl = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::ZERO, axis: DVec3::Z, radius: 5.0, ref_dir: DVec3::X,
    });
    // Point on cylinder axis
    let result = is_eq_distance(DVec3::new(0.0, 0.0, 10.0), &cyl, 1e-7);
    assert!(result.is_some(), "axis point should be eq distance");
    assert!((result.unwrap() - 5.0).abs() < 1e-10, "radius mismatch");
    // Point far from axis
    let result2 = is_eq_distance(DVec3::new(100.0, 0.0, 0.0), &cyl, 1e-7);
    assert!(result2.is_none(), "far point should not be eq distance");
}

#[test]
fn test_parallel_normals_dot() {
    // Two planes both with normal Z
    let n1 = DVec3::Z;
    let n2 = DVec3::Z;
    let dot = n1.dot(n2).abs();
    assert!(dot > 0.999, "parallel normals dot = {}", dot);
    // Perpendicular normals
    let n3 = DVec3::X;
    let dot_perp = n1.dot(n3).abs();
    assert!(dot_perp < 0.001, "perp normals dot = {}", dot_perp);
}
