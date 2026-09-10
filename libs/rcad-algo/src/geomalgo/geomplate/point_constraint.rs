//! OCCT GeomPlate_PointConstraint (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_PointConstraint.cxx (whole file).
//!
//! gp_Pnt -> DVec3, gp_Pnt2d -> DVec2 (architecture mapping).
//! `GeomLProp_SLProps myLProp` is a GAP leaf (the package is untranslated);
//! [`PointConstraint::lprop_surf`] keeps the OCCT SetParameters flow.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Surface3, SurfaceEval};

use super::curve_constraint::GeomLPropSlProps;

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
    // OCCT members myD11/myD12/myD21/myD22/myD23 (hxx) — filled by ctor 2;
    // the point-only ctor-1 leaves them uninitialized in OCCT, rcad zeroes.
    my_d11: DVec3,
    my_d12: DVec3,
    my_d21: DVec3,
    my_d22: DVec3,
    my_d23: DVec3,
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
            my_d11: DVec3::ZERO,
            my_d12: DVec3::ZERO,
            my_d21: DVec3::ZERO,
            my_d22: DVec3::ZERO,
            my_d23: DVec3::ZERO,
        }
    }

    /// OCCT ctor with a point on surface (GeomPlate_PointConstraint.cxx
    /// L53-73).
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_surface(
        u: f64,
        v: f64,
        surf: &Surface3,
        order: i32,
        tol_dist: f64,
        tol_ang: f64,
        tol_curv: f64,
    ) -> Self {
        // Surf->D2(myU, myV, myPoint, myD11, myD12, myD21, myD22, myD23).
        let (my_point, my_d11, my_d12, my_d21, my_d22, my_d23) = surf.derivatives2(u, v);
        // myLProp.SetSurface(Surf) — GAP (GeomLProp_SLProps), see lprop_surf.
        PointConstraint {
            my_order: order,
            my_point,
            my_u: u,
            my_v: v,
            my_tol_dist: tol_dist,
            my_tol_ang: tol_ang,
            my_tol_curv: tol_curv,
            has_pnt2d_on_surf: false,
            my_pt2d: DVec2::ZERO,
            my_d11,
            my_d12,
            my_d21,
            my_d22,
            my_d23,
        }
    }

    /// OCCT D0 (.cxx L78-81).
    pub fn d0(&self) -> DVec3 {
        self.my_point
    }

    /// OCCT D1 (.cxx L86-91) — returns (P, V1, V2).
    pub fn d1(&self) -> (DVec3, DVec3, DVec3) {
        (self.my_point, self.my_d11, self.my_d12)
    }

    /// OCCT D2 (.cxx L96-109) — returns (P, V1, V2, V3, V4, V5).
    #[allow(clippy::type_complexity)]
    pub fn d2(&self) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        (
            self.my_point,
            self.my_d11,
            self.my_d12,
            self.my_d21,
            self.my_d22,
            self.my_d23,
        )
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

    /// OCCT LPropSurf (.cxx L168-174).
    ///
    /// GAP leaf: the OCCT GeomLProp_SLProps return is an untranslated
    /// package; the carrier keeps the SetParameters(myU, myV) step.
    pub fn lprop_surf(&self) -> GeomLPropSlProps {
        GeomLPropSlProps {
            u: self.my_u,
            v: self.my_v,
        }
    }
}
