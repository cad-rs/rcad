//! Point-set equation fitting — GProp_PEquation equivalent.
//!
//! ✅ OCCT-aligned: given a set of 3D points, determines whether they are:
//! - Coincident (→ single point)
//! - Collinear (→ line)
//! - Coplanar (→ plane)
//! - Space-filling (→ full 3D box)
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/GProp_PEquation_Test.cxx

use glam::{DVec3, Vec3};

const TOL: f64 = 1e-7;

/// Result of point-set equation analysis.
///
/// OCCT: GProp_PEquation — classifies a set of points by fitting the
/// lowest-dimensional linear manifold that contains them within tolerance.
#[derive(Debug, Clone, PartialEq)]
pub enum PointSetKind {
    /// All points are coincident (within tolerance).
    Point(DVec3),
    /// Points lie on a line (within tolerance).
    Line(DVec3, DVec3), // origin, direction
    /// Points lie on a plane (within tolerance).
    Plane(DVec3, DVec3), // origin, normal
    /// Points span 3D space (or no valid lower-dim fit).
    Space,
}

/// Analyze a set of 3D points and determine their geometric relationship.
///
/// OCCT: GProp_PEquation (constructor). Uses PCA to classify the point set.
/// Returns the lowest-dimensional manifold that fits the points within `tolerance`.
pub fn analyze_point_set(points: &[DVec3], tolerance: f64) -> PointSetKind {
    if points.is_empty() {
        return PointSetKind::Space;
    }

    let n = points.len() as f64;

    // Centroid
    let mut centroid = DVec3::ZERO;
    for &p in points { centroid += p; }
    centroid /= n;

    // Check if all points are coincident
    let max_dist = points.iter().map(|p| (*p - centroid).length()).fold(0.0, f64::max);
    if max_dist < tolerance {
        return PointSetKind::Point(centroid);
    }

    // Covariance matrix (3x3)
    let mut cxx = 0.0; let mut cxy = 0.0; let mut cxz = 0.0;
    let mut cyy = 0.0; let mut cyz = 0.0; let mut czz = 0.0;
    for &p in points {
        let dx = p.x - centroid.x;
        let dy = p.y - centroid.y;
        let dz = p.z - centroid.z;
        cxx += dx * dx; cxy += dx * dy; cxz += dx * dz;
        cyy += dy * dy; cyz += dy * dz; czz += dz * dz;
    }
    cxx /= n; cxy /= n; cxz /= n;
    cyy /= n; cyz /= n; czz /= n;

    // Power iteration to find largest eigenvalue/eigenvector
    let mut v = DVec3::new(1.0, 0.0, 0.0);
    for _ in 0..20 {
        let v2 = DVec3::new(
            cxx * v.x + cxy * v.y + cxz * v.z,
            cxy * v.x + cyy * v.y + cyz * v.z,
            cxz * v.x + cyz * v.y + czz * v.z,
        );
        let len = v2.length();
        if len > 1e-30 { v = v2 / len; } else { break; }
    }
    let eval1 = DVec3::new(
        cxx * v.x + cxy * v.y + cxz * v.z,
        cxy * v.x + cyy * v.y + cyz * v.z,
        cxz * v.x + cyz * v.y + czz * v.z,
    ).dot(v);

    // Deflate largest eigenvalue
    cxx -= eval1 * v.x * v.x;
    cxy -= eval1 * v.x * v.y;
    cxz -= eval1 * v.x * v.z;
    cyy -= eval1 * v.y * v.y;
    cyz -= eval1 * v.y * v.z;
    czz -= eval1 * v.z * v.z;

    // Second eigenvector
    let mut v2 = if (v - DVec3::X).length() > 0.1 { DVec3::new(0.0, 1.0, 0.0) } else { DVec3::new(1.0, 0.0, 0.0) };
    for _ in 0..20 {
        let v2_new = DVec3::new(
            cxx * v2.x + cxy * v2.y + cxz * v2.z,
            cxy * v2.x + cyy * v2.y + cyz * v2.z,
            cxz * v2.x + cyz * v2.y + czz * v2.z,
        );
        let len = v2_new.length();
        if len > 1e-30 { v2 = v2_new / len; } else { break; }
    }
    let eval2 = DVec3::new(
        cxx * v2.x + cxy * v2.y + cxz * v2.z,
        cxy * v2.x + cyy * v2.y + cyz * v2.z,
        cxz * v2.x + cyz * v2.y + czz * v2.z,
    ).dot(v2);

    // Third eigenvector = cross of first two
    let v3 = v.cross(v2);
    let eval3 = {
        let mut cxx3 = cxx - eval2 * v2.x * v2.x;
        let mut cxy3 = cxy - eval2 * v2.x * v2.y;
        let mut cxz3 = cxz - eval2 * v2.x * v2.z;
        let mut cyy3 = cyy - eval2 * v2.y * v2.y;
        let mut cyz3 = cyz - eval2 * v2.y * v2.z;
        let mut czz3 = czz - eval2 * v2.z * v2.z;
        cxx3 * v3.x * v3.x + 2.0 * cxy3 * v3.x * v3.y + 2.0 * cxz3 * v3.x * v3.z
            + cyy3 * v3.y * v3.y + 2.0 * cyz3 * v3.y * v3.z + czz3 * v3.z * v3.z
    };

    let sqrt_eval1 = eval1.abs().sqrt();
    let sqrt_eval2 = eval2.abs().sqrt();
    let sqrt_eval3 = eval3.abs().sqrt();

    // Check residual perpendicular to v3 (plane check)
    let plane_residual = sqrt_eval3;
    // Check residual perpendicular to v (line check) — use v2 and v3
    let line_residual = (sqrt_eval2 * sqrt_eval2 + sqrt_eval3 * sqrt_eval3).sqrt();

    if line_residual < tolerance.max(1e-10) {
        // Collinear: direction = v (first principal direction)
        PointSetKind::Line(centroid, v.normalize_or_zero())
    } else if plane_residual < tolerance.max(1e-10) {
        // Coplanar: normal = v3 (cross of first two principal directions)
        PointSetKind::Plane(centroid, v3.normalize_or_zero())
    } else {
        PointSetKind::Space
    }
}

// =============================================================================
// Tests — OCCT GProp_PEquation_Test.cxx
// =============================================================================

#[cfg(test)]
mod pequation_tests {
    use super::*;

    #[test]
    fn pequation_coincident_points() {
        let pts = vec![
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(1.0, 2.0, 3.0),
            DVec3::new(1.0, 2.0, 3.0),
        ];
        let result = analyze_point_set(&pts, 1e-6);
        match result {
            PointSetKind::Point(p) => {
                assert!((p - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
            }
            other => panic!("Expected Point, got {:?}", other),
        }
    }

    #[test]
    fn pequation_collinear_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
        ];
        let result = analyze_point_set(&pts, 1e-6);
        match result {
            PointSetKind::Line(origin, dir) => {
                assert!(dir.x.abs() > 0.9, "line direction should be along X");
            }
            other => panic!("Expected Line, got {:?}", other),
        }
    }

    #[test]
    fn pequation_coplanar_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
        ];
        let result = analyze_point_set(&pts, 1e-6);
        match result {
            PointSetKind::Plane(origin, normal) => {
                assert!(normal.z.abs() > 0.9, "plane normal should be along Z, got {:?}", normal);
            }
            other => panic!("Expected Plane, got {:?}", other),
        }
    }

    #[test]
    fn pequation_space_points() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
        ];
        let result = analyze_point_set(&pts, 1e-6);
        match result {
            PointSetKind::Space => {} // expected
            other => panic!("Expected Space, got {:?}", other),
        }
    }
}
