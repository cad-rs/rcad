// OCCT IntSurf (TKGeomAlgo) — support classes for surface intersection.
//
// IntSurf_PntOn2S: a 3D point with its parameters on two surfaces.
// IntSurf_Quadric: an implicit (analytic) surface in a unified form.
// IntSurf_Transition / TypeTrans / Situation: transition of an intersection
// line relative to a restriction arc on a surface (mirrors IntRes2d; rcad
// reuses crate::topalgo::int_res2d::{Transition, TypeTrans, Situation}).

pub mod line_on_2s;
pub mod quadric;

pub use line_on_2s::LineOn2S;
pub use quadric::{Quadric, QuadricType};

use glam::{DVec2, DVec3};

/// OCCT IntSurf_PntOn2S (IntSurf_PntOn2S.cxx) — a point on two surfaces.
#[derive(Debug, Clone)]
pub struct PntOn2S {
    pt: DVec3,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
}

impl PntOn2S {
    pub fn new() -> Self {
        PntOn2S {
            pt: DVec3::ZERO,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        }
    }

    /// OCCT SetValue(Pt, OnFirst, U, V) — set the 3D point and one surface's UV.
    pub fn set_value(&mut self, pt: DVec3, on_first: bool, u: f64, v: f64) {
        self.pt = pt;
        self.set_value_uv(on_first, u, v);
    }

    /// OCCT SetValue(OnFirst, U, V) — set only one surface's UV.
    pub fn set_value_uv(&mut self, on_first: bool, u: f64, v: f64) {
        if on_first {
            self.u1 = u;
            self.v1 = v;
        } else {
            self.u2 = u;
            self.v2 = v;
        }
    }

    /// OCCT SetValue(Pt) — set only the 3D point.
    pub fn set_value_pt(&mut self, pt: DVec3) {
        self.pt = pt;
    }

    /// OCCT Value() — the 3D point.
    pub fn value(&self) -> DVec3 {
        self.pt
    }

    /// OCCT ValueOnSurface(OnFirst).
    pub fn value_on_surface(&self, on_first: bool) -> DVec2 {
        if on_first {
            DVec2::new(self.u1, self.v1)
        } else {
            DVec2::new(self.u2, self.v2)
        }
    }

    /// OCCT Parameters(u1, v1, u2, v2) — all four parameters.
    pub fn parameters(&self) -> (f64, f64, f64, f64) {
        (self.u1, self.v1, self.u2, self.v2)
    }

    /// OCCT ParametersOnSurface(OnFirst, U, V).
    pub fn parameters_on_surface(&self, on_first: bool) -> (f64, f64) {
        if on_first {
            (self.u1, self.v1)
        } else {
            (self.u2, self.v2)
        }
    }

    /// OCCT IsSame (IntSurf_PntOn2S.cxx L85-113) — 3D proximity + optional 2D proximity.
    pub fn is_same(&self, other: &PntOn2S, tol_3d: f64, tol_2d: f64) -> bool {
        if self.pt.distance_squared(other.value()) > tol_3d * tol_3d {
            return false;
        }
        if tol_2d < 0.0 {
            // We need not compare 2D-coordinates of the points
            return true;
        }
        let (a_u1, a_v1, a_u2, a_v2) = other.parameters();
        let mut a_p1 = DVec2::new(self.u1, self.v1);
        let a_p2 = DVec2::new(a_u1, a_v1);
        if !a_p1.abs_diff_eq(a_p2, tol_2d) {
            return false;
        }
        a_p1 = DVec2::new(self.u2, self.v2);
        let a_p2 = DVec2::new(a_u2, a_v2);
        a_p1.abs_diff_eq(a_p2, tol_2d)
    }
}

impl Default for PntOn2S {
    fn default() -> Self {
        Self::new()
    }
}
