// OCCT BRepAdaptor — adaptors for curves, surfaces, edges, faces.
//
// OCCT ref: TKBRep/BRepAdaptor/
//
// Provides OCCT-style adaptors for working with shapes.
// rcad: wraps DS edge/face access in adaptor pattern,
// so that bop (TKBO) code uses the same interface as OCCT.
//
// Also provides projection utilities that bridge bop → topalgo → rcad-kernel.

use crate::bop::ds::DS;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3, CurveEval, Curve2dEval, SurfaceEval};

/// OCCT BRepAdaptor_Curve — adapts an edge for curve queries.
///
/// In OCCT: wraps TopoDS_Edge with curve type, range, tolerance.
/// rcad: wraps a DS edge index.
pub struct CurveAdaptor<'a> {
    ds: &'a DS,
    edge_idx: usize,
    curve: Option<Curve3>,
    range: [f64; 2],
}

impl<'a> CurveAdaptor<'a> {
    /// Create adaptor from DS edge index.
    /// OCCT: BRepAdaptor_Curve(aE) — takes TopoDS_Edge.
    pub fn new(ds: &'a DS, edge_idx: usize) -> Self {
        let curve = ds.edge_curve(edge_idx).cloned();
        let range = ds.edge_range(edge_idx);
        CurveAdaptor { ds, edge_idx, curve, range }
    }

    /// OCCT: Value(T) — point on curve at parameter T.
    pub fn value(&self, t: f64) -> DVec3 {
        self.curve.as_ref().map(|c| c.point_at(t)).unwrap_or(DVec3::ZERO)
    }

    /// OCCT: D0(T, P) — same as Value.
    pub fn d0(&self, t: f64) -> DVec3 { self.value(t) }

    /// OCCT: D1(T, P, V1) — point and first derivative.
    pub fn d1(&self, t: f64) -> (DVec3, DVec3) {
        let p = self.value(t);
        let dt = 1e-7;
        let p2 = self.value(t + dt);
        (p, (p2 - p) / dt)
    }

    /// OCCT: FirstParameter() / LastParameter() — parameter range.
    pub fn first_parameter(&self) -> f64 { self.range[0] }
    pub fn last_parameter(&self) -> f64 { self.range[1] }

    /// OCCT: GetType() — curve type.
    pub fn get_type(&self) -> GeomAbsCurveType {
        match &self.curve {
            Some(Curve3::Line(_)) => GeomAbsCurveType::Line,
            Some(Curve3::Circle(_)) => GeomAbsCurveType::Circle,
            Some(Curve3::Ellipse(_)) => GeomAbsCurveType::Ellipse,
            Some(Curve3::BSpline(_)) => GeomAbsCurveType::BSpline,
            _ => GeomAbsCurveType::Other,
        }
    }

    /// OCCT: Resolution(R3D) — parameter tolerance from 3D tolerance.
    pub fn resolution(&self, r3d: f64) -> f64 {
        let dt = 1e-7;
        let p1 = self.value(self.range[0]);
        let p2 = self.value(self.range[0] + dt);
        let der_mag = (p2 - p1).length().max(1e-12);
        r3d / der_mag
    }

    /// Access the underlying curve.
    pub fn curve(&self) -> Option<&Curve3> { self.curve.as_ref() }
    pub fn range(&self) -> [f64; 2] { self.range }
    pub fn edge_idx(&self) -> usize { self.edge_idx }
}

/// OCCT GeomAbs_CurveType — curve type enum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeomAbsCurveType {
    Line, Circle, Ellipse, Hyperbola, Parabola, BSpline, Bezier, Offset, Other,
}

/// OCCT BRepAdaptor_Surface — adapts a face for surface queries.
///
/// rcad: wraps a DS face index.
pub struct SurfaceAdaptor<'a> {
    ds: &'a DS,
    face_idx: usize,
    surf: Option<Surface3>,
}

impl<'a> SurfaceAdaptor<'a> {
    pub fn new(ds: &'a DS, face_idx: usize) -> Self {
        let surf = ds.face_surface(face_idx);
        SurfaceAdaptor { ds, face_idx, surf }
    }

    /// OCCT: Value(U,V) — point on surface.
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        self.surf.as_ref().map(|s| s.point_at(u, v)).unwrap_or(DVec3::ZERO)
    }

    /// OCCT: D1(U,V,P,dU,dV) — point and first partial derivatives.
    pub fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        self.surf.as_ref().map(|s| s.derivatives(u, v))
            .unwrap_or((DVec3::ZERO, DVec3::ZERO, DVec3::ZERO))
    }

    /// OCCT: GetType() — surface type.
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        match &self.surf {
            Some(Surface3::Plane(_)) => GeomAbsSurfaceType::Plane,
            Some(Surface3::Cylinder(_)) => GeomAbsSurfaceType::Cylinder,
            Some(Surface3::Cone(_)) => GeomAbsSurfaceType::Cone,
            Some(Surface3::Sphere(_)) => GeomAbsSurfaceType::Sphere,
            Some(Surface3::Torus(_)) => GeomAbsSurfaceType::Torus,
            Some(Surface3::BSpline(_)) => GeomAbsSurfaceType::BSpline,
            _ => GeomAbsSurfaceType::Other,
        }
    }

    /// OCCT: UResolution(R3D) — U tolerance from 3D tolerance.
    pub fn u_resolution(&self, r3d: f64) -> f64 {
        let eps = 1e-6;
        let p = self.value(0.0, 0.0);
        let pu = self.value(eps, 0.0);
        let du = (pu - p).length().max(1e-12);
        r3d / du
    }

    /// OCCT: VResolution(R3D) — V tolerance from 3D tolerance.
    pub fn v_resolution(&self, r3d: f64) -> f64 {
        let eps = 1e-6;
        let p = self.value(0.0, 0.0);
        let pv = self.value(0.0, eps);
        let dv = (pv - p).length().max(1e-12);
        r3d / dv
    }

    pub fn surface(&self) -> Option<&Surface3> { self.surf.as_ref() }
    pub fn face_idx(&self) -> usize { self.face_idx }
}

/// OCCT GeomAbs_SurfaceType — surface type enum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeomAbsSurfaceType {
    Plane, Cylinder, Cone, Sphere, Torus, BSpline, Bezier, Offset, Other,
}

/// OCCT GeomAPI_ProjectPointOnCurve — projects a 3D point onto a curve.
pub fn project_point_on_curve(curve: &Curve3, query: DVec3) -> (f64, DVec3) {
    let proj = rcad_kernel::base::extrema::closest_point_on_curve(curve, query, 128);
    (proj.param, proj.point)
}

/// OCCT GeomAPI_ProjectPointOnSurf — projects a 3D point onto a surface.
pub fn project_point_on_surface(surface: &Surface3, point: DVec3) -> (DVec2, DVec3) {
    let proj = rcad_kernel::base::geom_api::project::closest_point_on_surface_near(
        surface, point, 64.0, 1e-7);
    (DVec2::new(proj.params.0, proj.params.1), proj.point)
}
