//! Surface type classification for analytic surfaces.
use rcad_kernel::geom::Surface3;

/// OCCT GeomAbs_SurfaceType equivalent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeomAbsSurfaceType {
    Plane, Cylinder, Cone, Sphere, Torus, Bezier, BSpline,
    SurfaceOfRevolution, SurfaceOfExtrusion, OtherSurface,
}

/// Classify a surface into GeomAbs_SurfaceType.
pub fn classify_surface_type(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::BSpline(_) | Surface3::Bezier(_) => GeomAbsSurfaceType::BSpline,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}

/// Quadric surface.
#[derive(Debug, Clone)]
pub enum Quadric {
    Plane(rcad_kernel::geom::Plane),
    Cylinder(rcad_kernel::geom::CylindricalSurface),
    Cone(rcad_kernel::geom::ConicalSurface),
    Sphere(rcad_kernel::geom::SphericalSurface),
    Torus(rcad_kernel::geom::ToroidalSurface),
}
