//! Point-set properties — GProp_PGProps equivalent.
//!
//! ✅ OCCT-aligned: computes mass (number of points), centre of mass (barycentre),
//! and matrix of inertia from a set of 3D points.
//!
//! OCCT source: src/ModelingData/TKGeomBase/GTests/GProp_PGProps_Test.cxx

use glam::DVec3;

const TOL: f64 = 1e-7;

/// Point-set properties accumulator.
///
/// OCCT: GProp_PGProps — accumulates point masses, computes centre of mass
/// and inertia tensor.
#[derive(Debug, Clone)]
pub struct PointSetProps {
    mass: f64,
    centre: DVec3,
    inertia: [f64; 6], // Ixx, Iyy, Izz, Ixy, Ixz, Iyz
}

impl PointSetProps {
    pub fn new() -> Self {
        Self { mass: 0.0, centre: DVec3::ZERO, inertia: [0.0; 6] }
    }

    /// Total mass (= number of points with unit weight).
    pub fn mass(&self) -> f64 { self.mass }

    /// Centre of mass (barycentre).
    pub fn centre_of_mass(&self) -> DVec3 { self.centre }

    /// Matrix of inertia (3×3 symmetric) in the global frame.
    pub fn matrix_of_inertia(&self) -> [f64; 9] {
        let [ixx, iyy, izz, ixy, ixz, iyz] = self.inertia;
        [
            ixx, ixy, ixz,
            ixy, iyy, iyz,
            ixz, iyz, izz,
        ]
    }

    /// Add a single point (unit weight).
    pub fn add_point(&mut self, pt: DVec3) {
        self.add_point_weighted(pt, 1.0);
    }

    /// Add a single point with weight.
    pub fn add_point_weighted(&mut self, pt: DVec3, weight: f64) {
        // Incremental formula for centre of mass:
        // new_centre = (old_mass * old_centre + weight * pt) / (old_mass + weight)
        let new_mass = self.mass + weight;
        if new_mass > 0.0 {
            self.centre = (self.centre * self.mass + pt * weight) / new_mass;
        }
        self.mass = new_mass;

        // Accumulate inertia contribution (parallel axis theorem applied at end)
        // For a point mass at position (x,y,z):
        // Ixx = m*(y²+z²), Iyy = m*(x²+z²), Izz = m*(x²+y²)
        // Ixy = -m*x*y, Ixz = -m*x*z, Iyz = -m*y*z
        let x = pt.x;
        let y = pt.y;
        let z = pt.z;
        self.inertia[0] += weight * (y * y + z * z); // Ixx
        self.inertia[1] += weight * (x * x + z * z); // Iyy
        self.inertia[2] += weight * (x * x + y * y); // Izz
        self.inertia[3] -= weight * x * y; // Ixy
        self.inertia[4] -= weight * x * z; // Ixz
        self.inertia[5] -= weight * y * z; // Iyz
    }

    /// Compute barycentre (centre of mass) of a set of points.
    ///
    /// OCCT: GProp_PGProps::Barycentre()
    pub fn barycentre(points: &[DVec3]) -> DVec3 {
        let mut props = PointSetProps::new();
        for &pt in points {
            props.add_point(pt);
        }
        if props.mass > 0.0 { props.centre } else { DVec3::ZERO }
    }
}

// =============================================================================
// Tests — OCCT GProp_PGProps_Test.cxx
// =============================================================================

#[cfg(test)]
mod pgprops_tests {
    use super::*;

    #[test]
    fn pgprops_empty_set() {
        let props = PointSetProps::new();
        assert!((props.mass() - 0.0).abs() < TOL);
    }

    #[test]
    fn pgprops_single_point() {
        let mut props = PointSetProps::new();
        props.add_point(DVec3::new(1.0, 2.0, 3.0));
        assert!((props.mass() - 1.0).abs() < TOL);
        assert!((props.centre_of_mass() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn pgprops_two_points_barycentre() {
        let mut props = PointSetProps::new();
        props.add_point(DVec3::new(0.0, 0.0, 0.0));
        props.add_point(DVec3::new(2.0, 4.0, 6.0));
        assert!((props.mass() - 2.0).abs() < TOL);
        assert!((props.centre_of_mass() - DVec3::new(1.0, 2.0, 3.0)).length() < TOL);
    }

    #[test]
    fn pgprops_weighted_points() {
        let mut props = PointSetProps::new();
        props.add_point_weighted(DVec3::new(0.0, 0.0, 0.0), 1.0);
        props.add_point_weighted(DVec3::new(4.0, 0.0, 0.0), 3.0);
        assert!((props.mass() - 4.0).abs() < TOL);
        // Weighted centroid: (0*1 + 4*3) / 4 = 3
        assert!((props.centre_of_mass() - DVec3::new(3.0, 0.0, 0.0)).length() < TOL);
    }

    #[test]
    fn pgprops_array_constructor() {
        let pts = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, -1.0, 0.0),
        ];
        let mut props = PointSetProps::new();
        for &pt in &pts { props.add_point(pt); }
        assert!((props.mass() - 4.0).abs() < TOL);
        assert!((props.centre_of_mass()).length() < TOL);
    }

    #[test]
    fn pgprops_static_barycentre() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
            DVec3::new(0.0, 6.0, 0.0),
        ];
        let g = PointSetProps::barycentre(&pts);
        assert!((g - DVec3::new(1.0, 2.0, 0.0)).length() < TOL);
    }

    #[test]
    fn pgprops_matrix_of_inertia_symmetric_points() {
        let pts = vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, -1.0, 0.0),
        ];
        let mut props = PointSetProps::new();
        for &pt in &pts { props.add_point(pt); }
        let m = props.matrix_of_inertia();
        // Ixx = sum(y²+z²) = (0) + (0) + (1) + (1) = 2
        assert!((m[0] - 2.0).abs() < TOL, "Ixx={}", m[0]);
        // Iyy = sum(x²+z²) = 1 + 1 + 0 + 0 = 2
        assert!((m[4] - 2.0).abs() < TOL, "Iyy={}", m[4]);
        // Izz = sum(x²+y²) = 1 + 1 + 1 + 1 = 4
        assert!((m[8] - 4.0).abs() < TOL, "Izz={}", m[8]);
        // Off-diagonal should be zero (symmetric)
        assert!((m[1]).abs() < TOL, "Ixy={}", m[1]);
        assert!((m[2]).abs() < TOL, "Ixz={}", m[2]);
        assert!((m[5]).abs() < TOL, "Iyz={}", m[5]);
    }
}
