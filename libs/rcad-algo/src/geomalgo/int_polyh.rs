//! OCCT IntPolyh_Point (IntPolyh_Point.cxx / .hxx) — a 3D point with its UV
//! parameters on a surface, used by the IntPolyh triangle-triangle
//! intersection chain.
//!
//! 1:1 translation of IntPolyh_Point.cxx L23-148 (Add/Sub/Divide/
//! Multiplication/SquareModulus/SquareDistance/Dot/Cross) plus the inline
//! accessors of IntPolyh_Point.hxx.

/// OCCT IntPolyh_Point.
#[derive(Debug, Clone, Copy, Default)]
pub struct IntPolyhPoint {
    x: f64,
    y: f64,
    z: f64,
    u: f64,
    v: f64,
}

impl IntPolyhPoint {
    /// OCCT IntPolyh_Point() — all components zero.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT IntPolyh_Point(x, y, z, u, v).
    pub fn new_uv(x: f64, y: f64, z: f64, u: f64, v: f64) -> Self {
        IntPolyhPoint { x, y, z, u, v }
    }

    /// OCCT IntPolyh_Point::X() (IntPolyh_Point.hxx L54).
    pub fn x(&self) -> f64 {
        self.x
    }

    /// OCCT IntPolyh_Point::Y().
    pub fn y(&self) -> f64 {
        self.y
    }

    /// OCCT IntPolyh_Point::Z().
    pub fn z(&self) -> f64 {
        self.z
    }

    /// OCCT IntPolyh_Point::U().
    pub fn u(&self) -> f64 {
        self.u
    }

    /// OCCT IntPolyh_Point::V().
    pub fn v(&self) -> f64 {
        self.v
    }

    /// OCCT IntPolyh_Point::SetX.
    pub fn set_x(&mut self, v: f64) {
        self.x = v;
    }

    /// OCCT IntPolyh_Point::SetY.
    pub fn set_y(&mut self, v: f64) {
        self.y = v;
    }

    /// OCCT IntPolyh_Point::SetZ.
    pub fn set_z(&mut self, v: f64) {
        self.z = v;
    }

    /// OCCT IntPolyh_Point::SetU.
    pub fn set_u(&mut self, v: f64) {
        self.u = v;
    }

    /// OCCT IntPolyh_Point::SetV.
    pub fn set_v(&mut self, v: f64) {
        self.v = v;
    }

    /// OCCT IntPolyh_Point::Add (IntPolyh_Point.cxx L39-49) — component-wise
    /// addition of all five components.
    pub fn add(&self, p1: &IntPolyhPoint) -> IntPolyhPoint {
        let mut res = IntPolyhPoint::new();
        res.set_x(self.x + p1.x());
        res.set_y(self.y + p1.y());
        res.set_z(self.z + p1.z());
        res.set_u(self.u + p1.u());
        res.set_v(self.v + p1.v());
        res
    }

    /// OCCT IntPolyh_Point::Sub (L53-63).
    pub fn sub(&self, p1: &IntPolyhPoint) -> IntPolyhPoint {
        let mut res = IntPolyhPoint::new();
        res.set_x(self.x - p1.x());
        res.set_y(self.y - p1.y());
        res.set_z(self.z - p1.z());
        res.set_u(self.u - p1.u());
        res.set_v(self.v - p1.v());
        res
    }

    /// OCCT IntPolyh_Point::Divide (L67-79): divide all five components by
    /// RR; when |RR| <= Precision::Computational() (machine epsilon) the
    /// default (zero) point is returned.
    pub fn divide(&self, rr: f64) -> IntPolyhPoint {
        let mut res = IntPolyhPoint::new();
        if rr.abs() > rcad_kernel::core::precision::COMPUTATIONAL {
            res.set_x(self.x / rr);
            res.set_y(self.y / rr);
            res.set_z(self.z / rr);
            res.set_u(self.u / rr);
            res.set_v(self.v / rr);
        }
        res
    }

    /// OCCT IntPolyh_Point::Multiplication (L83-93).
    pub fn multiplication(&self, rr: f64) -> IntPolyhPoint {
        let mut res = IntPolyhPoint::new();
        res.set_x(self.x * rr);
        res.set_y(self.y * rr);
        res.set_z(self.z * rr);
        res.set_u(self.u * rr);
        res.set_v(self.v * rr);
        res
    }

    /// OCCT IntPolyh_Point::SquareModulus (L97-101) — X^2 + Y^2 + Z^2.
    pub fn square_modulus(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// OCCT IntPolyh_Point::SquareDistance (L105-110) — 3D squared distance.
    pub fn square_distance(&self, p2: &IntPolyhPoint) -> f64 {
        (self.x - p2.x) * (self.x - p2.x) + (self.y - p2.y) * (self.y - p2.y)
            + (self.z - p2.z) * (self.z - p2.z)
    }

    /// OCCT IntPolyh_Point::Dot (L114-118) — 3D dot product.
    pub fn dot(&self, b: &IntPolyhPoint) -> f64 {
        self.x * b.x + self.y * b.y + self.z * b.z
    }

    /// OCCT IntPolyh_Point::Cross (L122-127) — sets self to a x b (3D).
    pub fn cross(&mut self, a: &IntPolyhPoint, b: &IntPolyhPoint) {
        self.x = a.y * b.z - a.z * b.y;
        self.y = a.z * b.x - a.x * b.z;
        self.z = a.x * b.y - a.y * b.x;
    }
}
