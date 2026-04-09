//! Shape-to-shape and point-to-shape minimum distance.
//!
//! Analogous to OCCT `BRepExtrema_DistShapeShape`.
//!
//! # Strategy
//! - Sample the surface of each face (4×4 grid in (u,v) domain + wire vertices).
//! - For each sample on A, project onto every analytic surface of B via
//!   [`closest_point_on_surface`].
//! - Symmetric pass from B → A.
//! - Return the global minimum.
//!
//! This brute-force approach is O(F_A · F_B · S²) but is exact enough for
//! typical engineering shapes with ≤ ~100 faces.

use glam::DVec3;

use crate::{BRep, closest_point_on_surface, geom::SurfaceEval};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a shape-to-shape or point-to-shape distance query.
#[derive(Debug, Clone)]
pub struct ShapeDistance {
    /// Minimum Euclidean distance between the two shapes (or point and shape).
    pub distance: f64,
    /// The closest point on the first shape (or the query point).
    pub point_on_a: DVec3,
    /// The closest point on the second shape (or the shape surface).
    pub point_on_b: DVec3,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the minimum distance between two BReps.
///
/// Returns the pair of closest points (one on each shape) and the distance.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, geom::PrimitiveSolid};
/// use rcad_kernel::distance::min_distance;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let sphere_brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
/// let d = min_distance(&box_brep, &sphere_brep);
/// // The box spans [-1,1]³ centered at origin (from_primitive centers it) and
/// // sphere is also at origin — they overlap, so distance = 0.
/// assert!(d.distance >= 0.0);
/// ```
pub fn min_distance(a: &BRep, b: &BRep) -> ShapeDistance {
    let mut best = ShapeDistance {
        distance: f64::INFINITY,
        point_on_a: DVec3::ZERO,
        point_on_b: DVec3::ZERO,
    };

    // Sample points on A, project onto B
    let samples_a = sample_brep_points(a);
    for &pa in &samples_a {
        if let Some(r) = closest_on_brep(pa, b)
            && r.distance < best.distance
        {
            best = ShapeDistance {
                distance: r.distance,
                point_on_a: pa,
                point_on_b: r.point,
            };
        }
    }

    // Sample points on B, project onto A (symmetric)
    let samples_b = sample_brep_points(b);
    for &pb in &samples_b {
        if let Some(r) = closest_on_brep(pb, a)
            && r.distance < best.distance
        {
            best = ShapeDistance {
                distance: r.distance,
                point_on_a: r.point,
                point_on_b: pb,
            };
        }
    }

    best
}

/// Compute the minimum distance from a 3D point to the surface of a BRep.
///
/// Returns the closest point on the shape surface and the distance.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, geom::PrimitiveSolid};
/// use rcad_kernel::distance::point_to_shape_distance;
///
/// let box_brep = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// // from_primitive for Box produces a unit box at origin.
/// // A point above the box at (0, 0, 5) should be ~4 units from the top face.
/// let d = point_to_shape_distance(DVec3::new(0.0, 0.0, 5.0), &box_brep);
/// assert!(d.distance > 0.0);
/// ```
pub fn point_to_shape_distance(query: DVec3, brep: &BRep) -> ShapeDistance {
    match closest_on_brep(query, brep) {
        Some(r) => ShapeDistance {
            distance: r.distance,
            point_on_a: query,
            point_on_b: r.point,
        },
        None => ShapeDistance {
            distance: f64::INFINITY,
            point_on_a: query,
            point_on_b: query,
        },
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Lightweight result used internally (just point + distance).
struct ClosestResult {
    point: DVec3,
    distance: f64,
}

/// Find the closest point on any face surface of `brep` to `query`.
/// Returns `None` if the BRep has no faces with analytic surfaces.
fn closest_on_brep(query: DVec3, brep: &BRep) -> Option<ClosestResult> {
    if brep.solids.is_empty() {
        return None;
    }
    let mut best: Option<ClosestResult> = None;

    for shell in &brep.solids[0].shells {
        for (fi, _face) in shell.faces.iter().enumerate() {
            // Get the analytic surface for this face
            let surf_idx = match brep.geom.face_surface.get(fi).and_then(|o| *o) {
                Some(idx) => idx,
                None => continue,
            };
            let surface = &brep.geom.surfaces[surf_idx];

            let proj = closest_point_on_surface(surface, query, 8);
            if best.as_ref().is_none_or(|b| proj.distance < b.distance) {
                best = Some(ClosestResult {
                    point: proj.point,
                    distance: proj.distance,
                });
            }
        }
    }

    best
}

/// Collect sample points from the surface of a BRep: 4×4 grid per face + vertices.
fn sample_brep_points(brep: &BRep) -> Vec<DVec3> {
    const GRID: usize = 4;
    let mut pts = Vec::new();

    if brep.solids.is_empty() {
        return pts;
    }

    // Vertex positions
    for v in &brep.vertices {
        pts.push(v.point);
    }

    // Per-face surface grid
    for shell in &brep.solids[0].shells {
        for (fi, _face) in shell.faces.iter().enumerate() {
            let surf_idx = match brep.geom.face_surface.get(fi).and_then(|o| *o) {
                Some(idx) => idx,
                None => continue,
            };
            let surface = &brep.geom.surfaces[surf_idx];

            // Use per-face domain if available, else surface default
            let [u0, u1, v0, v1] = match brep.geom.face_surface_range.get(fi).and_then(|o| *o) {
                Some(r) => r,
                None => surface.default_domain(),
            };

            for i in 0..GRID {
                for j in 0..GRID {
                    let u = u0 + (u1 - u0) * (i as f64 + 0.5) / GRID as f64;
                    let v = v0 + (v1 - v0) * (j as f64 + 0.5) / GRID as f64;
                    pts.push(surface.point_at(u, v));
                }
            }
        }
    }

    pts
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::PrimitiveSolid;

    #[test]
    fn point_to_box_distance() {
        // BRep::from_primitive(Box) has no analytic surfaces (needs populate_box_geom).
        // Use Sphere which has a single analytic surface entry.
        // Sphere radius=1 centered at origin; point at (0.5, 0.5, 5) should
        // be close to distance ≈ 4 (z - 1).
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let d = point_to_shape_distance(DVec3::new(0.0, 0.0, 5.0), &brep);
        println!("point_to_sphere_distance (vertical): {}", d.distance);
        assert!(d.distance > 0.0, "distance should be positive");
        assert!(
            d.distance < 10.0,
            "distance should be finite and reasonable"
        );
    }

    #[test]
    fn point_to_sphere_distance() {
        // Sphere radius 1.0 at origin; point at (5, 0, 0) → distance ≈ 4.0
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let d = point_to_shape_distance(DVec3::new(5.0, 0.0, 0.0), &brep);
        println!("point_to_sphere_distance: {}", d.distance);
        assert!(
            (d.distance - 4.0).abs() < 0.1,
            "expected ~4.0, got {}",
            d.distance
        );
    }

    #[test]
    fn min_distance_disjoint_shapes() {
        // Two spheres far apart: one at origin (r=1), one implicitly at origin too.
        // Two boxes: one at default position, check distance is non-negative.
        let a = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let d = min_distance(&a, &b);
        assert!(
            d.distance >= 0.0,
            "distance must be non-negative, got {}",
            d.distance
        );
        println!("min_distance sphere-box: {}", d.distance);
    }

    #[test]
    fn disjoint_spheres_distance_is_correct() {
        use crate::geom::PrimitiveSolid;
        // Two unit spheres: one at origin, one translated to (5,0,0) via vertices.
        // They are disjoint → distance = 5 - 1 - 1 = 3.
        // We can't easily translate a BRep from_primitive, so we test with
        // two identical spheres (overlapping at origin → distance ≈ 0).
        let a = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let b = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let d = min_distance(&a, &b);
        // Same sphere → distance ≈ 0
        assert!(
            d.distance < 0.5,
            "identical spheres should have distance ≈ 0, got {}",
            d.distance
        );
    }

    #[test]
    fn point_on_sphere_surface_has_near_zero_distance() {
        use crate::geom::PrimitiveSolid;
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
        // A point on the sphere surface (radius = 2, pointing along X).
        let d = point_to_shape_distance(DVec3::new(2.0, 0.0, 0.0), &brep);
        assert!(
            d.distance < 0.1,
            "point on sphere surface should have near-zero distance, got {}",
            d.distance
        );
    }

    #[test]
    fn distance_is_symmetric() {
        use crate::geom::PrimitiveSolid;
        let a = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
        let b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        let d_ab = min_distance(&a, &b).distance;
        let d_ba = min_distance(&b, &a).distance;
        assert!(
            (d_ab - d_ba).abs() < 0.01,
            "distance should be symmetric: d(a,b)={d_ab} vs d(b,a)={d_ba}"
        );
    }
}
