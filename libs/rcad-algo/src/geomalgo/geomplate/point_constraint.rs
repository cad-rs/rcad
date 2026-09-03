//! OCCT GeomPlate_PointConstraint (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_PointConstraint.cxx (whole file) for the point-only path.
//!
//! The second (surface-bound) constructor and D1/D2/LPropSurf depend on
//! GeomLProp_SLProps; they are out of the point-constraint anchor scope and
//! are left unported (OCCT ctor 2, L53-73; D1 L86-91; D2 L96-109;
//! LPropSurf L168-174).
//! gp_Pnt -> DVec3, gp_Pnt2d -> DVec2 (architecture mapping).

use glam::{DVec2, DVec3};

/// OCCT GeomPlate_PointConstraint.
#[derive(Debug, Clone, PartialEq)]
pub struct PointConstraint {
    my_order: i32,
    my_point: DVec3,
    my_u: f64,
    my_v: f64,
    my_tol_dist: f64,
    my_tol_ang: f64,
    my_tol_curv: f64,
    has_pnt2d_on_surf: bool,
    my_pt2d: DVec2,
}

impl PointConstraint {
    /// OCCT ctor with a point (GeomPlate_PointConstraint.cxx L31-48).
    pub fn new(pt: DVec3, order: i32, tol_dist: f64) -> Self {
        if order > 1 || order < -1 {
            panic!("GeomPlate_PointConstraint : the constraint must 0 or -1 with a point");
        }
        PointConstraint {
            my_order: order,
            my_point: pt,
            my_u: 0.0,
            my_v: 0.0,
            my_tol_dist: tol_dist,
            my_tol_ang: 0.0,
            my_tol_curv: 0.0,
            has_pnt2d_on_surf: false,
            my_pt2d: DVec2::ZERO,
        }
    }

    /// OCCT D0 (.cxx L78-81).
    pub fn d0(&self) -> DVec3 {
        self.my_point
    }

    /// OCCT SetG0Criterion (.cxx L114-117).
    pub fn set_g0_criterion(&mut self, tol_dist: f64) {
        self.my_tol_dist = tol_dist;
    }

    /// OCCT SetG1Criterion (.cxx L122-125).
    pub fn set_g1_criterion(&mut self, tol_ang: f64) {
        self.my_tol_ang = tol_ang;
    }

    /// OCCT SetG2Criterion (.cxx L130-133).
    pub fn set_g2_criterion(&mut self, tol_curv: f64) {
        self.my_tol_curv = tol_curv;
    }

    /// OCCT G0Criterion (.cxx L138-141).
    pub fn g0_criterion(&self) -> f64 {
        self.my_tol_dist
    }

    /// OCCT G1Criterion (.cxx L146-149).
    pub fn g1_criterion(&self) -> f64 {
        self.my_tol_ang
    }

    /// OCCT G2Criterion (.cxx L154-157).
    pub fn g2_criterion(&self) -> f64 {
        self.my_tol_curv
    }

    /// OCCT Order (.cxx L179-182).
    pub fn order(&self) -> i32 {
        self.my_order
    }

    /// OCCT SetOrder (.cxx L187-190).
    pub fn set_order(&mut self, order: i32) {
        self.my_order = order;
    }

    /// OCCT HasPnt2dOnSurf (.cxx L195-198).
    pub fn has_pnt2d_on_surf(&self) -> bool {
        self.has_pnt2d_on_surf
    }

    /// OCCT SetPnt2dOnSurf (.cxx L203-207).
    pub fn set_pnt2d_on_surf(&mut self, pnt2d: DVec2) {
        self.my_pt2d = pnt2d;
        self.has_pnt2d_on_surf = true;
    }

    /// OCCT Pnt2dOnSurf (.cxx L212-215).
    pub fn pnt2d_on_surf(&self) -> DVec2 {
        self.my_pt2d
    }
}
