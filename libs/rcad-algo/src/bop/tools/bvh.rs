//! BVH — re-exports from rcad-kernel + Aabb helper for gtest migration.
pub use rcad_kernel::math::bvh::{Bvh, BvhStats};
pub use rcad_kernel::math::bnd::BndBox;

use glam::DVec3;

/// Simple axis-aligned bounding box. Matches the `Aabb` used by GTest translations.
#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: DVec3,
    pub max: DVec3,
    pub gap: f64,
}

impl Aabb {
    pub fn empty() -> Self { Self { min: DVec3::splat(f64::INFINITY), max: DVec3::splat(f64::NEG_INFINITY), gap: 0.0 } }
    pub fn from_points(pts: &[DVec3]) -> Self { let mut a = Self::empty(); for &p in pts { a.expand_point(p); } a }
    pub fn expand_point(&mut self, p: DVec3) { self.min = self.min.min(p); self.max = self.max.max(p); }
    pub fn expand_aabb(&mut self, other: &Aabb) { self.min = self.min.min(other.min); self.max = self.max.max(other.max); self.gap = self.gap.max(other.gap); }
    pub fn center(&self) -> DVec3 { (self.min + self.max) * 0.5 }
    pub fn surface_area(&self) -> f64 { let d = self.max - self.min; if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 { return 0.0; } 2.0 * (d.x * d.y + d.y * d.z + d.z * d.x) }
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x - self.gap <= other.max.x + other.gap && self.max.x + self.gap >= other.min.x - other.gap
            && self.min.y - self.gap <= other.max.y + other.gap && self.max.y + self.gap >= other.min.y - other.gap
            && self.min.z - self.gap <= other.max.z + other.gap && self.max.z + self.gap >= other.min.z - other.gap
    }
    pub fn contains_point(&self, p: DVec3) -> bool { p.cmpge(self.min).all() && p.cmple(self.max).all() }
    pub fn point_dist_sq(&self, p: DVec3) -> f64 { let clamped = p.clamp(self.min, self.max); (p - clamped).length_squared() }
}
