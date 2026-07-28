//! OCCT GeomAbs_SurfaceType.hxx — surface type classification

use rcad_kernel::geom::Surface3;

/// OCCT GeomAbs_SurfaceType.hxx
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

/// Convert rcad Surface3 to OCCT GeomAbs_SurfaceType.
pub fn classify_surface_type(surf: &Surface3) -> GeomAbsSurfaceType {
    match surf {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::BSpline(_) | Surface3::Bezier(_) => GeomAbsSurfaceType::BSplineSurface,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}
