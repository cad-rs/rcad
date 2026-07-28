//! BRepLProp — local properties of edges and faces in a BRep.
//!
//! OCCT TKBRep BRepLProp package: BRepLProp_CLProps (edge local properties),
//! BRepLProp_SLProps (surface local properties).

use glam::DVec3;
use rcad_kernel::{Curve3, CurveEval, Surface3, SurfaceEval};
use crate::adaptor::{EdgeAdaptor, FaceAdaptor};

const DEFAULT_RESOLUTION: f64 = 1e-7;

// =============================================================================
// BRepLProp_CLProps — Edge Local Properties
// =============================================================================

/// Computes local properties of a BRep edge (OCCT BRepLProp_CLProps).
///
/// Parameter values are normalised [0, 1] (EdgeAdaptor convention).
pub struct BRepEdgeLProps {
    curve: Curve3,
    range: [f64; 2],
    n: i32,
    resolution: f64,
    u: f64,
    pnt: Option<DVec3>,
    d1v: Option<DVec3>,
    d2v: Option<DVec3>,
    d3v: Option<DVec3>,
}

impl BRepEdgeLProps {
    /// Create from an EdgeAdaptor (OCCT: BRepLProp_CLProps(Curve, N, Resolution)).
    pub fn from_adaptor(adaptor: &EdgeAdaptor<'_>, n: i32) -> Self {
        let (curve, range) = match adaptor.curve() {
            Some(c) => (c.clone(), adaptor.curve_range()),
            None => {
                let p0 = adaptor.point_at(0.0);
                let p1 = adaptor.point_at(1.0);
                (Curve3::Line(rcad_kernel::geom::Line3::new(p0, p1 - p0)), [0.0, 1.0])
            }
        };
        BRepEdgeLProps { curve, range, n, resolution: DEFAULT_RESOLUTION, u: f64::NAN, pnt: None, d1v: None, d2v: None, d3v: None }
    }

    fn compute(&mut self) {
        if self.pnt.is_some() { return; }
        let [r0, r1] = self.range;
        let span = r1 - r0;
        let t = (if span.abs() < 1e-15 { r0 } else { r0 + self.u.clamp(0.0, 1.0) * span }).clamp(r0, r1);
        self.pnt = Some(self.curve.point_at(t));
        if self.n >= 1 { self.d1v = Some(self.curve.derivative_at(t)); }
        if self.n >= 2 { self.d2v = Some(self.curve.derivative2_at(t)); }
        if self.n >= 3 { self.d3v = Some(self.curve.derivative3_at(t)); }
    }

    pub fn set_resolution(&mut self, res: f64) { self.resolution = res; }
    pub fn set_parameter(&mut self, u: f64) { self.u = u; self.pnt = None; self.d1v = None; self.d2v = None; self.d3v = None; }
    pub fn value(&mut self) -> DVec3 { self.compute(); self.pnt.unwrap_or(DVec3::ZERO) }
    pub fn d1(&mut self) -> DVec3 { self.compute(); self.d1v.unwrap_or(DVec3::ZERO) }
    pub fn d2(&mut self) -> DVec3 { self.compute(); self.d2v.unwrap_or(DVec3::ZERO) }
    pub fn d3(&mut self) -> DVec3 { self.compute(); self.d3v.unwrap_or(DVec3::ZERO) }
    pub fn is_tangent_defined(&mut self) -> bool { self.compute(); self.d1v.map_or(false, |d| d.length_squared() > self.resolution * self.resolution) }
    pub fn tangent(&mut self) -> Option<DVec3> { self.compute(); self.d1v.and_then(|d| { let l2 = d.length_squared(); if l2 > self.resolution * self.resolution { Some(d / l2.sqrt()) } else { None } }) }
    pub fn curvature(&mut self) -> f64 { self.compute(); match (self.d1v, self.d2v) { (Some(d1), Some(d2)) => { let s = d1.length(); if s < self.resolution { 0.0 } else { d1.cross(d2).length() / (s * s * s) } } _ => 0.0 } }
    pub fn normal(&mut self) -> Option<DVec3> { let c = self.curvature(); if c < self.resolution { None } else { match (self.d1v, self.d2v) { (Some(d1), Some(d2)) => { let cr = d1.cross(d2); let n = cr.cross(d1); let l = n.length(); if l < self.resolution { None } else { Some(n / l) } } _ => None } } }
    pub fn centre_of_curvature(&mut self) -> Option<DVec3> { let c = self.curvature(); if c < self.resolution { None } else { self.normal().map(|n| self.pnt.unwrap_or(DVec3::ZERO) + n / c) } }
}

// =============================================================================
// BRepLProp_SLProps — Face Local Properties
// =============================================================================

/// Computes local properties of a BRep face (OCCT BRepLProp_SLProps).
pub struct BRepFaceLProps {
    surface: Surface3,
    n: i32, resolution: f64,
    u: f64, v: f64,
    pnt: Option<DVec3>,
    d1u_v: Option<DVec3>, d1v_v: Option<DVec3>,
    d2u_v: Option<DVec3>, d2v_v: Option<DVec3>, duv_v: Option<DVec3>,
    normal_v: Option<DVec3>, normal_status: i32,
    curv_computed: bool,
    min_curv: f64, max_curv: f64, mean_curv: f64, gaus_curv: f64,
    dir_min_curv: DVec3, dir_max_curv: DVec3,
}

impl BRepFaceLProps {
    pub fn new(n: i32) -> Self {
        BRepFaceLProps {
            surface: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
            n, resolution: DEFAULT_RESOLUTION, u: f64::NAN, v: f64::NAN,
            pnt: None, d1u_v: None, d1v_v: None, d2u_v: None, d2v_v: None, duv_v: None,
            normal_v: None, normal_status: 0, curv_computed: false,
            min_curv: 0.0, max_curv: 0.0, mean_curv: 0.0, gaus_curv: 0.0,
            dir_min_curv: DVec3::ZERO, dir_max_curv: DVec3::ZERO,
        }
    }

    /// Create from a FaceAdaptor (extracts surface).
    pub fn from_adaptor(adaptor: &FaceAdaptor<'_>, n: i32) -> Self {
        let mut props = BRepFaceLProps::new(n);
        if let Some(surface) = adaptor.surface() {
            props.set_surface(surface);
        }
        props
    }

    fn invalidate(&mut self) { self.pnt = None; self.d1u_v = None; self.d1v_v = None; self.d2u_v = None; self.d2v_v = None; self.duv_v = None; self.normal_v = None; self.normal_status = 0; self.curv_computed = false; }

    fn compute(&mut self) {
        if self.pnt.is_some() { return; }
        let (u, v) = (self.u, self.v);
        self.pnt = Some(self.surface.point_at(u, v));
        if self.n >= 1 { let (_, pu, pv) = self.surface.derivatives(u, v); self.d1u_v = Some(pu); self.d1v_v = Some(pv); }
        if self.n >= 2 { let (_, pu, pv, puu, puv, pvv) = self.surface.derivatives2(u, v); self.d1u_v = Some(pu); self.d1v_v = Some(pv); self.d2u_v = Some(puu); self.d2v_v = Some(pvv); self.duv_v = Some(puv); }
    }

    fn comp_normal(&mut self) {
        if self.normal_status != 0 { return; }
        self.compute();
        match (self.d1u_v, self.d1v_v) {
            (Some(du), Some(dv)) => {
                let n = du.cross(dv); let l = n.length();
                if l > self.resolution { self.normal_v = Some(n / l); self.normal_status = 1; } else { self.normal_status = -1; }
            }
            _ => { self.normal_status = -1; }
        }
    }

    fn comp_curvature(&mut self) {
        if self.curv_computed { return; }
        if self.n < 2 { self.curv_computed = true; return; }
        self.comp_normal();
        if self.normal_status != 1 { self.curv_computed = true; return; }
        let du = self.d1u_v.unwrap(); let dv = self.d1v_v.unwrap();
        let duu = self.d2u_v.unwrap(); let dvv = self.d2v_v.unwrap();
        let n = self.normal_v.unwrap();
        let e = du.dot(du).max(1e-7); let f = du.dot(dv); let g = dv.dot(dv).max(1e-7);
        let l = duu.dot(n); let m = self.duv_v.unwrap_or(DVec3::ZERO).dot(n); let nv = dvv.dot(n);
        let egmf = e * g - f * f;
        if egmf.abs() < 1e-30 { self.curv_computed = true; return; }
        let w11 = (l * g - m * f) / egmf; let w12 = (m * g - nv * f) / egmf;
        let w21 = (m * e - l * f) / egmf; let w22 = (nv * e - m * f) / egmf;
        let tr = w11 + w22; let det = w11 * w22 - w12 * w21;
        let sd = (tr * tr - 4.0 * det).max(0.0).sqrt();
        self.max_curv = (0.5 * (tr + sd)).max(0.5 * (tr - sd));
        self.min_curv = (0.5 * (tr + sd)).min(0.5 * (tr - sd));
        self.mean_curv = 0.5 * tr; self.gaus_curv = det;
        if sd > 1e-15 {
            self.dir_max_curv = (w12 * du + (0.5 * (tr + sd) - w11) * dv).normalize_or_zero();
            self.dir_min_curv = (w12 * du + (0.5 * (tr - sd) - w11) * dv).normalize_or_zero();
        }
        self.curv_computed = true;
    }

    pub fn set_resolution(&mut self, r: f64) { self.resolution = r; self.invalidate(); }
    pub fn set_surface(&mut self, s: &Surface3) { self.surface = s.clone(); self.invalidate(); }
    pub fn set_parameters(&mut self, u: f64, v: f64) { self.u = u; self.v = v; self.invalidate(); }
    pub fn value(&mut self) -> DVec3 { self.compute(); self.pnt.unwrap_or(DVec3::ZERO) }
    pub fn d1u(&mut self) -> DVec3 { self.compute(); self.d1u_v.unwrap_or(DVec3::ZERO) }
    pub fn d1v(&mut self) -> DVec3 { self.compute(); self.d1v_v.unwrap_or(DVec3::ZERO) }
    pub fn d2u(&mut self) -> DVec3 { self.compute(); self.d2u_v.unwrap_or(DVec3::ZERO) }
    pub fn d2v(&mut self) -> DVec3 { self.compute(); self.d2v_v.unwrap_or(DVec3::ZERO) }
    pub fn duv(&mut self) -> DVec3 { self.compute(); self.duv_v.unwrap_or(DVec3::ZERO) }
    pub fn is_tangent_u_defined(&mut self) -> bool { self.compute(); self.d1u_v.map_or(false, |d| d.length_squared() > self.resolution * self.resolution) }
    pub fn is_tangent_v_defined(&mut self) -> bool { self.compute(); self.d1v_v.map_or(false, |d| d.length_squared() > self.resolution * self.resolution) }
    pub fn tangent_u(&mut self) -> Option<DVec3> { if self.is_tangent_u_defined() { Some(self.d1u_v.unwrap().normalize_or_zero()) } else { None } }
    pub fn tangent_v(&mut self) -> Option<DVec3> { if self.is_tangent_v_defined() { Some(self.d1v_v.unwrap().normalize_or_zero()) } else { None } }
    pub fn is_normal_defined(&mut self) -> bool { self.comp_normal(); self.normal_status == 1 }
    pub fn normal(&mut self) -> Option<DVec3> { self.comp_normal(); self.normal_v }
    pub fn min_curvature(&mut self) -> f64 { self.comp_curvature(); self.min_curv }
    pub fn max_curvature(&mut self) -> f64 { self.comp_curvature(); self.max_curv }
    pub fn gaussian_curvature(&mut self) -> f64 { self.comp_curvature(); self.gaus_curv }
    pub fn mean_curvature(&mut self) -> f64 { self.comp_curvature(); self.mean_curv }
    pub fn dir_min_curvature(&mut self) -> Option<DVec3> { self.comp_curvature(); if self.curv_computed && self.min_curv.abs() > 1e-15 { Some(self.dir_min_curv) } else { None } }
    pub fn dir_max_curvature(&mut self) -> Option<DVec3> { self.comp_curvature(); if self.curv_computed && self.max_curv.abs() > 1e-15 { Some(self.dir_max_curv) } else { None } }
}

// =============================================================================
// Free convenience functions
// =============================================================================

pub fn edge_curvature_at(edge: &EdgeAdaptor<'_>, u: f64) -> f64 { let mut p = BRepEdgeLProps::from_adaptor(edge, 2); p.set_parameter(u); p.curvature() }
pub fn edge_tangent_at(edge: &EdgeAdaptor<'_>, u: f64) -> Option<DVec3> { let mut p = BRepEdgeLProps::from_adaptor(edge, 1); p.set_parameter(u); p.tangent() }
pub fn face_normal_uv(surface: &Surface3, u: f64, v: f64) -> Option<DVec3> { let mut p = BRepFaceLProps::new(1); p.set_surface(surface); p.set_parameters(u, v); p.normal() }

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::PrimitiveSolid;
    use rcad_kernel::BRep;

    #[test]
    fn test_edge_curvature_straight() {
        let b = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
        let a = EdgeAdaptor::new(&b, 0);
        assert!(edge_curvature_at(&a, 0.5) < 1e-10);
    }

    #[test]
    fn test_edge_tangent() {
        let b = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
        let a = EdgeAdaptor::new(&b, 0);
        let t = edge_tangent_at(&a, 0.5);
        assert!(t.is_some());
        assert!((t.unwrap().length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_edge_domain() {
        let b = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
        let a = EdgeAdaptor::new(&b, 0);
        let mut p = BRepEdgeLProps::from_adaptor(&a, 2);
        p.set_parameter(0.0); let p0 = p.value();
        p.set_parameter(1.0); let p1 = p.value();
        assert!(((p1 - p0).length() - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_face_normal_plane() {
        let pl = Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z));
        assert!((face_normal_uv(&pl, 0.5, 0.5).unwrap() - DVec3::Z).length() < 1e-10);
    }

    #[test]
    fn test_face_curvature_plane() {
        let pl = Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z));
        let mut p = BRepFaceLProps::new(2);
        p.set_surface(&pl); p.set_parameters(0.0, 0.0);
        assert!(p.gaussian_curvature().abs() < 1e-10);
        assert!(p.mean_curvature().abs() < 1e-10);
    }
}
