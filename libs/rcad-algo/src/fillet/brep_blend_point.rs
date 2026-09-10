//! OCCT Blend_Point (TKFillet/Blend) — 1:1 port of Blend_Point.hxx
//! (L30-304), Blend_Point.cxx (whole file L20-483) and Blend_Point.lxx
//! (L18-199, the inline accessors).
//!
//! Architecture mappings: `gp_Pnt` / `gp_Vec` map to `DVec3`; `gp_Vec2d`
//! maps to `DVec2`.  The OCCT constructors and the corresponding `SetValue`
//! overloads map to `new_*` constructors and `set_value_*` methods.

use glam::{DVec2, DVec3};

/// OCCT Blend_Point — point used in the blending computation, carrying the
/// parameter, the points and tangents on the two supports (surfaces or
/// curves) and the flags telling which supports are defined.
#[derive(Debug, Clone)]
pub struct BlendPoint {
    pt1: DVec3,
    pt2: DVec3,
    tg1: DVec3,
    tg2: DVec3,
    prm: f64,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    pc1: f64,
    pc2: f64,
    utg12d: f64,
    vtg12d: f64,
    utg22d: f64,
    vtg22d: f64,
    hass1: bool,
    hass2: bool,
    hasc1: bool,
    hasc2: bool,
    istgt: bool,
}

impl BlendPoint {
    /// OCCT Blend_Point::Blend_Point() (Blend_Point.cxx L20-23).
    pub fn new() -> Self {
        BlendPoint {
            pt1: DVec3::ZERO,
            pt2: DVec3::ZERO,
            tg1: DVec3::ZERO,
            tg2: DVec3::ZERO,
            prm: 0.0,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            pc1: 0.0,
            pc2: 0.0,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: false,
            hass2: false,
            hasc1: false,
            hasc2: false,
            istgt: true,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2, Tg1, Tg2, Tg12d,
    /// Tg22d) (Blend_Point.cxx L25-55) — a point on 2 surfaces, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_2_surfaces_with_tangents(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1,
            tg2,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1: 0.0,
            pc2: 0.0,
            utg12d: tg12d.x,
            vtg12d: tg12d.y,
            utg22d: tg22d.x,
            vtg22d: tg22d.y,
            hass1: true,
            hass2: true,
            hasc1: false,
            hasc2: false,
            istgt: false,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2)
    /// (Blend_Point.cxx L57-77) — a point on 2 surfaces, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_2_surfaces(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1: DVec3::ZERO,
            tg2: DVec3::ZERO,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1: 0.0,
            pc2: 0.0,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: true,
            hass2: true,
            hasc1: false,
            hasc2: false,
            istgt: true,
        }
    }

    /// OCCT Blend_Point(Pts, Ptc, Param, U, V, W, Tgs, Tgc, Tg2d)
    /// (Blend_Point.cxx L133-158) — a point on a surface and a curve,
    /// with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_surface_curve_with_tangents(
        ps: DVec3,
        pc: DVec3,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        tgs: DVec3,
        tgc: DVec3,
        tg2d: DVec2,
    ) -> Self {
        BlendPoint {
            pt1: ps,
            pt2: pc,
            tg1: tgs,
            tg2: tgc,
            prm: param,
            u1: u,
            v1: v,
            u2: 0.0,
            v2: 0.0,
            pc1: 0.0,
            pc2: w,
            utg12d: tg2d.x,
            vtg12d: tg2d.y,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: true,
            hass2: false,
            hasc1: false,
            hasc2: true,
            istgt: false,
        }
    }

    /// OCCT Blend_Point(Pts, Ptc, Param, U, V, W) (Blend_Point.cxx L160-178)
    /// — a point on a surface and a curve, without tangents.
    pub fn new_on_surface_curve(ps: DVec3, pc: DVec3, param: f64, u: f64, v: f64, w: f64) -> Self {
        BlendPoint {
            pt1: ps,
            pt2: pc,
            tg1: DVec3::ZERO,
            tg2: DVec3::ZERO,
            prm: param,
            u1: u,
            v1: v,
            u2: 0.0,
            v2: 0.0,
            pc1: 0.0,
            pc2: w,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: true,
            hass2: false,
            hasc1: false,
            hasc2: true,
            istgt: true,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2, PC, Tg1, Tg2,
    /// Tg12d, Tg22d) (Blend_Point.cxx L227-259) — a point on a surface and
    /// a curve on surface, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_surface_curve_on_surface_with_tangents(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1,
            tg2,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1: 0.0,
            pc2: pc,
            utg12d: tg12d.x,
            vtg12d: tg12d.y,
            utg22d: tg22d.x,
            vtg22d: tg22d.y,
            hass1: true,
            hass2: true,
            hasc1: false,
            hasc2: true,
            istgt: false,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2, PC)
    /// (Blend_Point.cxx L261-283) — a point on a surface and a curve on
    /// surface, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_surface_curve_on_surface(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1: DVec3::ZERO,
            tg2: DVec3::ZERO,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1: 0.0,
            pc2: pc,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: true,
            hass2: true,
            hasc1: false,
            hasc2: true,
            istgt: true,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2, PC1, PC2, Tg1, Tg2,
    /// Tg12d, Tg22d) (Blend_Point.cxx L343-377) — a point on two curves on
    /// surfaces, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_2_curves_on_surfaces_with_tangents(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1,
            tg2,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1,
            pc2,
            utg12d: tg12d.x,
            vtg12d: tg12d.y,
            utg22d: tg22d.x,
            vtg22d: tg22d.y,
            hass1: true,
            hass2: true,
            hasc1: true,
            hasc2: true,
            istgt: false,
        }
    }

    /// OCCT Blend_Point(Pt1, Pt2, Param, U1, V1, U2, V2, PC1, PC2)
    /// (Blend_Point.cxx L379-403) — a point on two curves on surfaces,
    /// without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn new_on_2_curves_on_surfaces(
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
    ) -> Self {
        BlendPoint {
            pt1: p1,
            pt2: p2,
            tg1: DVec3::ZERO,
            tg2: DVec3::ZERO,
            prm: param,
            u1,
            v1,
            u2,
            v2,
            pc1,
            pc2,
            utg12d: 0.0,
            vtg12d: 0.0,
            utg22d: 0.0,
            vtg22d: 0.0,
            hass1: true,
            hass2: true,
            hasc1: true,
            hasc2: true,
            istgt: true,
        }
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2, Tg1, Tg2, Tg12d, Tg22d)
    /// (Blend_Point.cxx L79-109) — a point on 2 surfaces, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_2_surfaces_with_tangents(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.hasc2 = false;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.x;
        self.vtg12d = tg12d.y;
        self.utg22d = tg22d.x;
        self.vtg22d = tg22d.y;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2)
    /// (Blend_Point.cxx L111-131) — a point on 2 surfaces, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_2_surfaces(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.hasc2 = false;
        self.istgt = true;
    }

    /// OCCT SetValue(Pts, Ptc, Param, U, V, W, Tgs, Tgc, Tg2d)
    /// (Blend_Point.cxx L180-205) — a point on a surface and a curve,
    /// with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_surface_curve_with_tangents(
        &mut self,
        ps: DVec3,
        pc: DVec3,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        tgs: DVec3,
        tgc: DVec3,
        tg2d: DVec2,
    ) {
        self.pt1 = ps;
        self.pt2 = pc;
        self.prm = param;
        self.u1 = u;
        self.v1 = v;
        self.hass1 = true;
        self.hass2 = false;
        self.hasc1 = false;
        self.pc2 = w;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tgs;
        self.tg2 = tgc;
        self.utg12d = tg2d.x;
        self.vtg12d = tg2d.y;
    }

    /// OCCT SetValue(Pts, Ptc, Param, U, V, W) (Blend_Point.cxx L207-225)
    /// — a point on a surface and a curve, without tangents.
    pub fn set_value_on_surface_curve(
        &mut self,
        ps: DVec3,
        pc: DVec3,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
    ) {
        self.pt1 = ps;
        self.pt2 = pc;
        self.prm = param;
        self.u1 = u;
        self.v1 = v;
        self.hass1 = true;
        self.hass2 = false;
        self.hasc1 = false;
        self.pc2 = w;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2, PC, Tg1, Tg2, Tg12d,
    /// Tg22d) (Blend_Point.cxx L285-317) — a point on a surface and a curve
    /// on surface, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_surface_curve_on_surface_with_tangents(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.pc2 = pc;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.x;
        self.vtg12d = tg12d.y;
        self.utg22d = tg22d.x;
        self.vtg22d = tg22d.y;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2, PC)
    /// (Blend_Point.cxx L319-341) — a point on a surface and a curve on
    /// surface, without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_surface_curve_on_surface(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.hasc1 = false;
        self.pc2 = pc;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2, PC1, PC2, Tg1, Tg2,
    /// Tg12d, Tg22d) (Blend_Point.cxx L405-439) — a point on two curves on
    /// surfaces, with tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_2_curves_on_surfaces_with_tangents(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
        tg1: DVec3,
        tg2: DVec3,
        tg12d: DVec2,
        tg22d: DVec2,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = false;
        self.tg1 = tg1;
        self.tg2 = tg2;
        self.utg12d = tg12d.x;
        self.vtg12d = tg12d.y;
        self.utg22d = tg22d.x;
        self.vtg22d = tg22d.y;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, U1, V1, U2, V2, PC1, PC2)
    /// (Blend_Point.cxx L441-465) — a point on two curves on surfaces,
    /// without tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn set_value_on_2_curves_on_surfaces(
        &mut self,
        p1: DVec3,
        p2: DVec3,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pc1: f64,
        pc2: f64,
    ) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.u1 = u1;
        self.v1 = v1;
        self.hass1 = true;
        self.u2 = u2;
        self.v2 = v2;
        self.hass2 = true;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// OCCT SetValue(Pt1, Pt2, Param, PC1, PC2) (Blend_Point.cxx L467-483)
    /// — a point on two curves.
    pub fn set_value_on_2_curves(&mut self, p1: DVec3, p2: DVec3, param: f64, pc1: f64, pc2: f64) {
        self.pt1 = p1;
        self.pt2 = p2;
        self.prm = param;
        self.hass1 = false;
        self.hass2 = false;
        self.pc1 = pc1;
        self.hasc1 = true;
        self.pc2 = pc2;
        self.hasc2 = true;
        self.istgt = true;
    }

    /// OCCT SetParameter (Blend_Point.lxx L18-21).
    pub fn set_parameter(&mut self, param: f64) {
        self.prm = param;
    }

    /// OCCT Parameter (Blend_Point.lxx L33-36).
    pub fn parameter(&self) -> f64 {
        self.prm
    }

    /// OCCT IsTangencyPoint (Blend_Point.lxx L58-61).
    pub fn is_tangency_point(&self) -> bool {
        self.istgt
    }

    /// OCCT PointOnS1 (Blend_Point.lxx L23-26).
    pub fn point_on_s1(&self) -> DVec3 {
        self.pt1
    }

    /// OCCT PointOnS2 (Blend_Point.lxx L28-31).
    pub fn point_on_s2(&self) -> DVec3 {
        self.pt2
    }

    /// OCCT ParametersOnS1 (Blend_Point.lxx L38-46).
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        if !self.hass1 {
            panic!("Standard_DomainError");
        }
        (self.u1, self.v1)
    }

    /// OCCT ParametersOnS2 (Blend_Point.lxx L48-56).
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        if !self.hass2 {
            panic!("Standard_DomainError");
        }
        (self.u2, self.v2)
    }

    /// OCCT TangentOnS1 (Blend_Point.lxx L63-70).
    pub fn tangent_on_s1(&self) -> DVec3 {
        if self.istgt {
            panic!("Standard_DomainError");
        }
        self.tg1
    }

    /// OCCT TangentOnS2 (Blend_Point.lxx L72-79).
    pub fn tangent_on_s2(&self) -> DVec3 {
        if self.istgt {
            panic!("Standard_DomainError");
        }
        self.tg2
    }

    /// OCCT Tangent2dOnS1 (Blend_Point.lxx L81-88).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        if self.istgt || !self.hass1 {
            panic!("Standard_DomainError");
        }
        DVec2::new(self.utg12d, self.vtg12d)
    }

    /// OCCT Tangent2dOnS2 (Blend_Point.lxx L90-97).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        if self.istgt || !self.hass2 {
            panic!("Standard_DomainError");
        }
        DVec2::new(self.utg22d, self.vtg22d)
    }

    /// OCCT PointOnS (Blend_Point.lxx L99-102).
    pub fn point_on_s(&self) -> DVec3 {
        self.pt1
    }

    /// OCCT PointOnC (Blend_Point.lxx L104-107).
    pub fn point_on_c(&self) -> DVec3 {
        self.pt2
    }

    /// OCCT ParametersOnS (Blend_Point.lxx L109-117).
    pub fn parameters_on_s(&self) -> (f64, f64) {
        if !self.hass1 {
            panic!("Standard_DomainError");
        }
        (self.u1, self.v1)
    }

    /// OCCT ParameterOnC (Blend_Point.lxx L119-126).
    pub fn parameter_on_c(&self) -> f64 {
        if !self.hasc2 {
            panic!("Standard_DomainError");
        }
        self.pc2
    }

    /// OCCT TangentOnS (Blend_Point.lxx L128-135).
    pub fn tangent_on_s(&self) -> DVec3 {
        if self.istgt || !self.hass1 {
            panic!("Standard_DomainError");
        }
        self.tg1
    }

    /// OCCT TangentOnC (Blend_Point.lxx L137-144).
    pub fn tangent_on_c(&self) -> DVec3 {
        if self.istgt {
            panic!("Standard_DomainError");
        }
        self.tg2
    }

    /// OCCT Tangent2d (Blend_Point.lxx L146-153).
    pub fn tangent_2d(&self) -> DVec2 {
        if self.istgt || !self.hass1 {
            panic!("Standard_DomainError");
        }
        DVec2::new(self.utg12d, self.vtg12d)
    }

    /// OCCT PointOnC1 (Blend_Point.lxx L155-158).
    pub fn point_on_c1(&self) -> DVec3 {
        self.pt1
    }

    /// OCCT PointOnC2 (Blend_Point.lxx L160-163).
    pub fn point_on_c2(&self) -> DVec3 {
        self.pt2
    }

    /// OCCT ParameterOnC1 (Blend_Point.lxx L165-172).
    pub fn parameter_on_c1(&self) -> f64 {
        if !self.hasc1 {
            panic!("Standard_DomainError");
        }
        self.pc1
    }

    /// OCCT ParameterOnC2 (Blend_Point.lxx L174-181).
    pub fn parameter_on_c2(&self) -> f64 {
        if !self.hasc2 {
            panic!("Standard_DomainError");
        }
        self.pc2
    }

    /// OCCT TangentOnC1 (Blend_Point.lxx L183-190).
    pub fn tangent_on_c1(&self) -> DVec3 {
        if self.istgt || !self.hass1 {
            panic!("Standard_DomainError");
        }
        self.tg1
    }

    /// OCCT TangentOnC2 (Blend_Point.lxx L192-199).
    pub fn tangent_on_c2(&self) -> DVec3 {
        if self.istgt {
            panic!("Standard_DomainError");
        }
        self.tg2
    }
}
