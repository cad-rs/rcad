// OCCT Adaptor3d_Surface (TKG3d) — the instance interface used by Contap.
//
// OCCT passes `occ::handle<Adaptor3d_Surface>` through the whole Contap
// package and reaches it via the static `Adaptor3d_HSurfaceTool` facade;
// every tool call is a one-liner to the handle's virtual method.  The rcad
// encoding keeps that shape: this object-safe trait is the
// `handle<Adaptor3d_Surface>` and Contap calls it directly.
//
// [`GeomSurfaceAdapter`] is the GeomAdaptor_Surface equivalent over a
// kernel surface, with an optional restricted UV window
// (BRepAdaptor_Surface::FirstUParameter semantics — the face's wire-bounds
// domain takes precedence over the natural surface domain; see the 2a-4
// notes in `docs/tkhlr-port-plan.md`).

use glam::DVec3;
use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, Plane, Point3, SphericalSurface, Surface3, SurfaceEval,
    ToroidalSurface, Vec3,
};
use rcad_kernel::precision::CONFUSION;

use crate::geomalgo::int_patch::{classify_surface_type, GeomAbsSurfaceType as SurfaceType};

/// OCCT Adaptor3d_Surface — the instance interface with the methods Contap
/// uses (a subset of Adaptor3d_HSurfaceTool.hxx L40-293).
pub trait SurfaceAdapter {
    /// OCCT FirstUParameter().
    fn first_u_parameter(&self) -> f64;
    /// OCCT LastUParameter().
    fn last_u_parameter(&self) -> f64;
    /// OCCT FirstVParameter().
    fn first_v_parameter(&self) -> f64;
    /// OCCT LastVParameter().
    fn last_v_parameter(&self) -> f64;
    /// OCCT Value(U, V).
    fn value(&self, u: f64, v: f64) -> Point3;
    /// OCCT D1(U, V, P, D1U, D1V) — P is the surface point itself.
    fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3);
    /// OCCT D2(U, V, P, D1U, D1V, D2U, D2V, D2UV).
    fn d2(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3);
    /// OCCT UResolution(R3d).
    fn u_resolution(&self, r3d: f64) -> f64;
    /// OCCT VResolution(R3d).
    fn v_resolution(&self, r3d: f64) -> f64;
    /// OCCT GetType().
    fn get_type(&self) -> SurfaceType;
    /// OCCT Plane() — valid when GetType() == Plane.
    fn plane(&self) -> Plane;
    /// OCCT Cylinder() — valid when GetType() == Cylinder.
    fn cylinder(&self) -> CylindricalSurface;
    /// OCCT Cone() — valid when GetType() == Cone.
    fn cone(&self) -> ConicalSurface;
    /// OCCT Sphere() — valid when GetType() == Sphere.
    fn sphere(&self) -> SphericalSurface;
    /// OCCT Torus() — valid when GetType() == Torus.
    fn torus(&self) -> ToroidalSurface;
    /// OCCT NbUPoles().
    fn nb_u_poles(&self) -> usize;
    /// OCCT NbVPoles().
    fn nb_v_poles(&self) -> usize;
    /// OCCT NbUKnots() — distinct knots.
    fn nb_u_knots(&self) -> usize;
    /// OCCT NbVKnots() — distinct knots.
    fn nb_v_knots(&self) -> usize;
    /// OCCT UDegree().
    fn u_degree(&self) -> usize;
    /// OCCT VDegree().
    fn v_degree(&self) -> usize;
    /// OCCT IsUPeriodic().
    fn is_u_periodic(&self) -> bool;
    /// OCCT UPeriod().
    fn u_period(&self) -> f64;
    /// OCCT IsVPeriodic().
    fn is_v_periodic(&self) -> bool;
    /// OCCT VPeriod().
    fn v_period(&self) -> f64;

    /// The [uinf, usup, vinf, vsup] window from the First/Last parameters —
    /// the bounds the IntWalk engine derives Um/Vm/UM/VM from.
    fn uv_window(&self) -> [f64; 4] {
        [
            self.first_u_parameter(),
            self.last_u_parameter(),
            self.first_v_parameter(),
            self.last_v_parameter(),
        ]
    }

    /// The kernel surface when the adapter wraps one (the
    /// Contap_TheIWalking instantiation runs the landed Surface3-based
    /// IntWalk engine; the OCCT engine is driven by the adaptor handle).
    fn kernel_surface(&self) -> Option<&Surface3>;
}

/// OCCT `occ::handle<Adaptor3d_Surface>` — shared ownership.
pub type SurfaceHandle = std::sync::Arc<dyn SurfaceAdapter>;

/// OCCT GeomAdaptor_Surface equivalent — an adapter over a kernel surface
/// with an optional restricted UV window (the face domain).
#[derive(Debug, Clone)]
pub struct GeomSurfaceAdapter {
    surf: Surface3,
    /// The restricted face UV window [u0, u1, v0, v1]; when None the natural
    /// surface domain is used.
    uv_domain: Option<[f64; 4]>,
}

impl GeomSurfaceAdapter {
    /// Adaptor over the natural surface domain.
    pub fn new(surf: Surface3) -> Self {
        GeomSurfaceAdapter {
            surf,
            uv_domain: None,
        }
    }

    /// Adaptor over a restricted face UV window (BRepAdaptor_Surface
    /// semantics: FirstU/LastU/FirstV/LastV return the face domain).
    pub fn with_domain(surf: Surface3, uv_domain: [f64; 4]) -> Self {
        GeomSurfaceAdapter {
            surf,
            uv_domain: Some(uv_domain),
        }
    }

    /// The UTrim/VTrim narrowed-window copy (OCCT GeomAdaptor_Surface
    /// keeps the same underlying surface).
    pub fn clone_with_window(&self, uv_domain: [f64; 4]) -> Self {
        GeomSurfaceAdapter {
            surf: self.surf.clone(),
            uv_domain: Some(uv_domain),
        }
    }

    /// The underlying kernel surface (the Geom_Surface of OCCT).
    pub fn surface3(&self) -> &Surface3 {
        &self.surf
    }

    fn u_domain(&self) -> [f64; 2] {
        match self.uv_domain {
            Some(d) => [d[0], d[1]],
            None => {
                let d = SurfaceEval::default_domain(&self.surf);
                // OCCT returns RealFirst()/RealLast() (±DBL_MAX) for the
                // infinite bounds; the kernel domains use ±infinity.
                [occt_real_first(d[0]), occt_real_last(d[1])]
            }
        }
    }

    fn v_domain(&self) -> [f64; 2] {
        match self.uv_domain {
            Some(d) => [d[2], d[3]],
            None => {
                let d = SurfaceEval::default_domain(&self.surf);
                [occt_real_first(d[2]), occt_real_last(d[3])]
            }
        }
    }

    /// OCCT Precision::Parametric(P) = P * 0.01 (Precision.hxx L328).
    fn parametric(p: f64) -> f64 {
        p * 0.01
    }
}

/// OCCT RealFirst() = -DBL_MAX (Standard_Real.hxx); ±infinite kernel bounds
/// are reported in that style so the exact `== RealFirst()` comparisons of
/// OCCT (Contap_HContTool.cxx L136-162) work unchanged.
fn occt_real_first(r: f64) -> f64 {
    if r == f64::NEG_INFINITY {
        f64::MIN
    } else {
        r
    }
}

/// OCCT RealLast() = DBL_MAX.
fn occt_real_last(r: f64) -> f64 {
    if r == f64::INFINITY {
        f64::MAX
    } else {
        r
    }
}

impl SurfaceAdapter for GeomSurfaceAdapter {
    fn first_u_parameter(&self) -> f64 {
        self.u_domain()[0]
    }
    fn last_u_parameter(&self) -> f64 {
        self.u_domain()[1]
    }
    fn first_v_parameter(&self) -> f64 {
        self.v_domain()[0]
    }
    fn last_v_parameter(&self) -> f64 {
        self.v_domain()[1]
    }
    fn value(&self, u: f64, v: f64) -> Point3 {
        SurfaceEval::point_at(&self.surf, u, v)
    }
    fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        SurfaceEval::derivatives(&self.surf, u, v)
    }
    fn d2(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        SurfaceEval::derivatives2(&self.surf, u, v)
    }
    /// OCCT GeomAdaptor_Surface::UResolution (cxx L1818-1892).
    fn u_resolution(&self, r3d: f64) -> f64 {
        let res;
        match self.get_type() {
            SurfaceType::Torus => {
                let t = &self.surf;
                let r = match t {
                    Surface3::Torus(t) => t.major_radius + t.minor_radius,
                    _ => 0.0,
                };
                res = if r > CONFUSION { r3d / (2. * r) } else { 0. };
            }
            SurfaceType::Sphere => {
                let r = match &self.surf {
                    Surface3::Sphere(s) => s.radius,
                    _ => 0.0,
                };
                res = if r > CONFUSION { r3d / (2. * r) } else { 0. };
            }
            SurfaceType::Cylinder => {
                let r = match &self.surf {
                    Surface3::Cylinder(c) => c.radius,
                    _ => 0.0,
                };
                res = if r > CONFUSION { r3d / (2. * r) } else { 0. };
            }
            SurfaceType::Cone => {
                let vd = self.v_domain();
                if vd[1] - vd[0] > 1.0e10 {
                    // Not truly bounded => unknown resolution
                    return Self::parametric(r3d);
                }
                // VIso(VFirst/Last) circle radii: r(V) = refR + V * tan(semi).
                let (ref_r, semi) = match &self.surf {
                    Surface3::Cone(c) => (c.radius, c.half_angle_rad),
                    _ => (0.0, 0.0),
                };
                let r1 = ref_r + vd[1] * semi.tan();
                let r2 = ref_r + vd[0] * semi.tan();
                let r = if r1 > r2 { r1 } else { r2 };
                return if r > CONFUSION { r3d / r } else { 0. };
            }
            SurfaceType::Plane => {
                return r3d;
            }
            _ => return Self::parametric(r3d),
        }
        if res <= 1. {
            2. * res.asin()
        } else {
            2. * std::f64::consts::PI
        }
    }
    /// OCCT GeomAdaptor_Surface::VResolution (cxx L1896-1959).
    fn v_resolution(&self, r3d: f64) -> f64 {
        let res;
        match self.get_type() {
            SurfaceType::Torus => {
                let r = match &self.surf {
                    Surface3::Torus(t) => t.minor_radius,
                    _ => 0.0,
                };
                res = if r > CONFUSION { r3d / (2. * r) } else { 0. };
            }
            SurfaceType::Sphere => {
                let r = match &self.surf {
                    Surface3::Sphere(s) => s.radius,
                    _ => 0.0,
                };
                res = if r > CONFUSION { r3d / (2. * r) } else { 0. };
            }
            SurfaceType::Cylinder
            | SurfaceType::Cone
            | SurfaceType::Plane
            | SurfaceType::SurfaceOfExtrusion => {
                return r3d;
            }
            _ => return Self::parametric(r3d),
        }
        if res <= 1. {
            2. * res.asin()
        } else {
            2. * std::f64::consts::PI
        }
    }
    fn get_type(&self) -> SurfaceType {
        classify_surface_type(&self.surf)
    }
    fn plane(&self) -> Plane {
        match &self.surf {
            Surface3::Plane(p) => p.clone(),
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Plane"),
        }
    }
    fn cylinder(&self) -> CylindricalSurface {
        match &self.surf {
            Surface3::Cylinder(c) => c.clone(),
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Cylinder"),
        }
    }
    fn cone(&self) -> ConicalSurface {
        match &self.surf {
            Surface3::Cone(c) => c.clone(),
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Cone"),
        }
    }
    fn sphere(&self) -> SphericalSurface {
        match &self.surf {
            Surface3::Sphere(s) => s.clone(),
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Sphere"),
        }
    }
    fn torus(&self) -> ToroidalSurface {
        match &self.surf {
            Surface3::Torus(t) => t.clone(),
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Torus"),
        }
    }
    fn nb_u_poles(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => b.control_points.len(),
            Surface3::Bezier(b) => b.control_points.len(),
            _ => 0,
        }
    }
    fn nb_v_poles(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => b.control_points.first().map_or(0, |r| r.len()),
            Surface3::Bezier(b) => b.control_points.first().map_or(0, |r| r.len()),
            _ => 0,
        }
    }
    fn nb_u_knots(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => {
                // Distinct knots (OCCT NbUKnots counts knot groups).
                let mut n = 0usize;
                for (i, k) in b.knots_u.iter().enumerate() {
                    if i == 0 || *k != b.knots_u[i - 1] {
                        n += 1;
                    }
                }
                n
            }
            _ => 0,
        }
    }
    fn nb_v_knots(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => {
                let mut n = 0usize;
                for (i, k) in b.knots_v.iter().enumerate() {
                    if i == 0 || *k != b.knots_v[i - 1] {
                        n += 1;
                    }
                }
                n
            }
            _ => 0,
        }
    }
    fn u_degree(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => b.degree_u,
            _ => 0,
        }
    }
    fn v_degree(&self) -> usize {
        match &self.surf {
            Surface3::BSpline(b) => b.degree_v,
            _ => 0,
        }
    }
    fn is_u_periodic(&self) -> bool {
        matches!(self.get_type(), SurfaceType::Cylinder | SurfaceType::Torus)
    }
    fn u_period(&self) -> f64 {
        match self.get_type() {
            SurfaceType::Cylinder | SurfaceType::Torus => 2. * std::f64::consts::PI,
            _ => 0.0,
        }
    }
    fn kernel_surface(&self) -> Option<&Surface3> {
        Some(&self.surf)
    }
    fn is_v_periodic(&self) -> bool {
        self.get_type() == SurfaceType::Torus
    }
    fn v_period(&self) -> f64 {
        match self.get_type() {
            SurfaceType::Torus => 2. * std::f64::consts::PI,
            _ => 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::CylindricalSurface;

    fn unit_cylinder(r: f64) -> Surface3 {
        Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: r,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        })
    }

    /// OCCT anchor: GeomAdaptor_Surface::UResolution for a cylinder of
    /// radius R: 2*asin(R3d / (2R)); plane: R3d; VResolution cylinder: R3d.
    #[test]
    fn geom_surface_adapter_resolution() {
        let cyl = GeomSurfaceAdapter::new(unit_cylinder(2.0));
        let expected = 2. * (CONFUSION / (2. * 2.0)).asin();
        assert!((cyl.u_resolution(CONFUSION) - expected).abs() < 1e-15);
        assert!((cyl.v_resolution(CONFUSION) - CONFUSION).abs() < 1e-15);

        let pln = GeomSurfaceAdapter::new(Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 0.0, 1.0),
            u_dir: DVec3::new(1.0, 0.0, 0.0),
            v_dir: DVec3::new(0.0, 1.0, 0.0),
        }));
        assert!((pln.u_resolution(CONFUSION) - CONFUSION).abs() < 1e-15);
        assert!((pln.v_resolution(CONFUSION) - CONFUSION).abs() < 1e-15);
    }

    /// OCCT anchor: the restricted UV window overrides the natural domain
    /// (BRepAdaptor_Surface::FirstUParameter face semantics).
    #[test]
    fn geom_surface_adapter_face_domain() {
        let cyl = GeomSurfaceAdapter::with_domain(unit_cylinder(1.0), [0.5, 1.5, -2.0, 3.0]);
        assert!((cyl.first_u_parameter() - 0.5).abs() < 1e-15);
        assert!((cyl.last_u_parameter() - 1.5).abs() < 1e-15);
        assert!((cyl.first_v_parameter() + 2.0).abs() < 1e-15);
        assert!((cyl.last_v_parameter() - 3.0).abs() < 1e-15);
    }
}
