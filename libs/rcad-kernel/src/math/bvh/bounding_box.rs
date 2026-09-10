//! OCCT BVH_Box<T, N> equivalent for the 3D case (TKMath BVH package).
//!
//! 1:1 translation of the `BVH_Box` semantics used by the OCCT gtests:
//! default-constructed box is empty (`myIsInited == false`), `Area()` is 0
//! until the box is initialized, `Center()` is the midpoint of the corners
//! and `Intersects()` is an inclusive corner-overlap test
//! (BVH_Box.hxx L89-239).

use glam::DVec3;

/// OCCT BVH_Box<Standard_Real, 3>.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Lower corner (OCCT `myMin`).
    pub min: DVec3,
    /// Upper corner (OCCT `myMax`).
    pub max: DVec3,
    /// OCCT `myIsInited`.
    pub is_inited: bool,
}

impl Default for Aabb {
    fn default() -> Self {
        Self::empty()
    }
}

impl Aabb {
    /// OCCT BVH_Box(): default constructor creates an uninitialized
    /// (empty) box.
    pub fn empty() -> Self {
        Aabb {
            min: DVec3::ZERO,
            max: DVec3::ZERO,
            is_inited: false,
        }
    }

    /// OCCT BVH_Box(thePoint): a box initialized to a single point.
    pub fn from_point(p: DVec3) -> Self {
        Aabb {
            min: p,
            max: p,
            is_inited: true,
        }
    }

    /// OCCT BVH_Box::Add(thePoint): enlarge the box to include the point
    /// (initializes the box on the first point).
    pub fn add_point(&mut self, p: DVec3) {
        if self.is_inited {
            self.min = self.min.min(p);
            self.max = self.max.max(p);
        } else {
            self.min = p;
            self.max = p;
            self.is_inited = true;
        }
    }

    /// OCCT BVH_Box(thePoints, theNbPnts) / the box over a point set.
    pub fn from_points(points: &[DVec3]) -> Self {
        let mut b = Aabb::empty();
        for p in points {
            b.add_point(*p);
        }
        b
    }

    /// OCCT BVH_Box::IsValid — the box was initialized with at least one
    /// point.
    pub fn is_valid(&self) -> bool {
        self.is_inited
    }

    /// OCCT BVH_Box::Center(axis) — midpoint of the axis range.
    pub fn center(&self) -> DVec3 {
        (self.min + self.max) * 0.5
    }

    /// OCCT BVH_Box::Intersects(theBox): true when the boxes overlap
    /// (inclusive comparison on every axis).
    pub fn intersects(&self, other: &Aabb) -> bool {
        self.is_inited
            && other.is_inited
            && self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
            && self.min.z <= other.max.z
            && self.max.z >= other.min.z
    }

    /// OCCT BVH_Box::Area — for the 3D case the product of the axis sizes
    /// (0.0 while the box is empty).
    pub fn surface_area(&self) -> f64 {
        if !self.is_inited {
            return 0.0;
        }
        let d = self.max - self.min;
        d.x * d.y * d.z
    }
}
