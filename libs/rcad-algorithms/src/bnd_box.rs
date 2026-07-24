//! Axis-aligned bounding box (AABB) — OCCT Bnd_Box equivalent.
//!
//! Stores min/max extents and a gap (tolerance).  `intersects()` applies
//! the gap on both sides, matching OCCT `Bnd_Box::IsOut`:
//! `self.min - self.gap <= other.max + other.gap  &&  self.max + self.gap >= other.min - other.gap`

use glam::DVec3;

#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub min: DVec3,
    pub max: DVec3,
    pub gap: f64,
}

impl Aabb {
    /// Empty AABB (`min > max`, contains no points).
    pub fn empty() -> Self {
        Self {
            min: DVec3::splat(f64::INFINITY),
            max: DVec3::splat(f64::NEG_INFINITY),
            gap: 0.0,
        }
    }

    /// Build an AABB that encloses all given points.
    pub fn from_points(pts: &[DVec3]) -> Self {
        let mut aabb = Self::empty();
        for &p in pts {
            aabb.expand_point(p);
        }
        aabb
    }

    /// Expand to include one point.
    pub fn expand_point(&mut self, p: DVec3) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// Expand to include another AABB.  OCCT Bnd_Box::Add.
    pub fn expand_aabb(&mut self, other: &Aabb) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
        self.gap = self.gap.max(other.gap);
    }

    /// AABB center.
    pub fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    /// Surface area (for SAH cost).
    pub fn surface_area(&self) -> f64 {
        let d = self.max - self.min;
        if d.x < 0.0 || d.y < 0.0 || d.z < 0.0 {
            return 0.0; // empty AABB
        }
        2.0 * (d.x * d.y + d.y * d.z + d.z * d.x)
    }

    /// Whether this AABB intersects another.
    /// OCCT Bnd_Box::IsOut inverted — returns true when boxes overlap.
    /// Gap is added on both sides: min-gap <= max+gap && max+gap >= min-gap.
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.min.x - self.gap <= other.max.x + other.gap
            && self.max.x + self.gap >= other.min.x - other.gap
            && self.min.y - self.gap <= other.max.y + other.gap
            && self.max.y + self.gap >= other.min.y - other.gap
            && self.min.z - self.gap <= other.max.z + other.gap
            && self.max.z + self.gap >= other.min.z - other.gap
    }

    /// Ray–AABB intersection; returns entry parameter `t` along the ray (forward hits only).
    pub fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        let t1 = (self.min - origin) * inv_dir;
        let t2 = (self.max - origin) * inv_dir;

        let t_min = t1.min(t2);
        let t_max = t1.max(t2);

        let t_enter = t_min.x.max(t_min.y).max(t_min.z);
        let t_exit = t_max.x.min(t_max.y).min(t_max.z);

        if t_exit >= t_enter.max(0.0) {
            Some(t_enter.max(0.0))
        } else {
            None
        }
    }

    /// Squared minimum distance from a point to this AABB.
    pub fn point_dist_sq(&self, p: DVec3) -> f64 {
        let clamped = p.clamp(self.min, self.max);
        (p - clamped).length_squared()
    }

    /// Whether a point lies inside the AABB (inclusive of the boundary).
    pub fn contains_point(&self, p: DVec3) -> bool {
        p.cmpge(self.min).all() && p.cmple(self.max).all()
    }
}
