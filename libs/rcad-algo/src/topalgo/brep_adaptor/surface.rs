// OCCT BRepAdaptor_Surface (TKTopAlgo/BRepAdaptor) — the face-level surface
// adaptor, in the BRepAdaptor_Surface::Initialize(F, Restriction = true)
// form used by HLRBRep_Surface::Surface (HLRBRep_Surface.cxx L36-39).
//
// The rcad encoding follows the BRepAdaptor_Curve2d precedent
// (`topalgo/brep_adaptor/curve2d.rs`): `(brep: &'a BRep, my_face: Shape)`
// carries the OCCT TopoDS_Face handle; the restricted GeomAdaptor_Surface
// role (the face UV window over the location-transformed world surface) is
// delegated to `contap::GeomSurfaceAdapter`.

use rcad_kernel::geom::{
    ConicalSurface, CylindricalSurface, Plane, Point3, SphericalSurface, Surface3, ToroidalSurface,
    Vec3,
};
use rcad_kernel::topo::topods::{BRepTool, Shape, TShape};

use crate::geomalgo::int_patch::{classify_surface_type, GeomAbsSurfaceType};
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};

/// OCCT BRepAdaptor_Surface (restricted form).
#[derive(Debug, Clone)]
pub struct BRepAdaptorSurface<'a> {
    /// OCCT myFace.
    my_face: Shape,
    /// OCCT mySurf — the restricted GeomAdaptor over the world (located)
    /// face surface.
    my_surf: GeomSurfaceAdapter,
    _brep: std::marker::PhantomData<&'a ()>,
}

impl<'a> BRepAdaptorSurface<'a> {
    /// OCCT BRepAdaptor_Surface() — an undefined adaptor.
    pub fn new() -> Self {
        BRepAdaptorSurface {
            my_face: Shape::null(),
            my_surf: GeomSurfaceAdapter::new(Surface3::Plane(Plane {
                origin: Point3::ZERO,
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_dir: Vec3::new(1.0, 0.0, 0.0),
                v_dir: Vec3::new(0.0, 1.0, 0.0),
            })),
            _brep: std::marker::PhantomData,
        }
    }

    /// OCCT Initialize(F, R = true) — loads the face; with the restriction
    /// flag the adaptor window is the face's UV bounds (TFaceData.uv_domain,
    /// the 2a-4 face-window semantics), otherwise the natural surface
    /// domain.
    pub fn initialize_face(brep: &'a rcad_kernel::BRep, f: &Shape, restriction: bool) -> Self {
        let world = brep
            .face_surface_world(f)
            .unwrap_or_else(|| Surface3::Plane(Plane {
                origin: Point3::ZERO,
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_dir: Vec3::new(1.0, 0.0, 0.0),
                v_dir: Vec3::new(0.0, 1.0, 0.0),
            }));
        let uv_domain = if restriction {
            match &*f.data {
                TShape::Face(fd) => fd.uv_domain,
                _ => None,
            }
        } else {
            None
        };
        let my_surf = match uv_domain {
            Some(d) => GeomSurfaceAdapter::with_domain(world, d),
            None => GeomSurfaceAdapter::new(world),
        };
        BRepAdaptorSurface {
            my_face: f.clone(),
            my_surf,
            _brep: std::marker::PhantomData,
        }
    }

    /// OCCT Face() — the loaded face.
    pub fn face(&self) -> &Shape {
        &self.my_face
    }

    /// The restricted GeomAdaptor surface (OCCT mySurf, the
    /// GeomAdaptor_Surface member).
    pub fn adaptor_surface(&self) -> &GeomSurfaceAdapter {
        &self.my_surf
    }

    /// OCCT Bezier/BSpline Poles() — the pole grid (rows = U index) of the
    /// world surface; empty for the non-pole surface types.
    pub fn poles_grid(&self) -> Vec<Vec<Point3>> {
        match self.my_surf.surface3() {
            Surface3::Bezier(b) => b.control_points.clone(),
            Surface3::BSpline(b) => b.control_points.clone(),
            _ => Vec::new(),
        }
    }
}

impl Default for BRepAdaptorSurface<'_> {
    fn default() -> Self {
        Self::new()
    }
}

// The BRepAdaptor_Surface surface queries delegate to the restricted
// GeomAdaptor (OCCT BRepAdaptor_Surface forwards to mySurf the same way).
impl<'a> BRepAdaptorSurface<'a> {
    pub fn first_u_parameter(&self) -> f64 {
        self.my_surf.first_u_parameter()
    }
    pub fn last_u_parameter(&self) -> f64 {
        self.my_surf.last_u_parameter()
    }
    pub fn first_v_parameter(&self) -> f64 {
        self.my_surf.first_v_parameter()
    }
    pub fn last_v_parameter(&self) -> f64 {
        self.my_surf.last_v_parameter()
    }
    pub fn value(&self, u: f64, v: f64) -> Point3 {
        self.my_surf.value(u, v)
    }
    pub fn d1(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        self.my_surf.d1(u, v)
    }
    pub fn d2(&self, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        self.my_surf.d2(u, v)
    }
    pub fn u_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.u_resolution(r3d)
    }
    pub fn v_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.v_resolution(r3d)
    }
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        self.my_surf.get_type()
    }
    pub fn plane(&self) -> Plane {
        self.my_surf.plane()
    }
    pub fn cylinder(&self) -> CylindricalSurface {
        self.my_surf.cylinder()
    }
    pub fn cone(&self) -> ConicalSurface {
        self.my_surf.cone()
    }
    pub fn sphere(&self) -> SphericalSurface {
        self.my_surf.sphere()
    }
    pub fn torus(&self) -> ToroidalSurface {
        self.my_surf.torus()
    }
    pub fn nb_u_poles(&self) -> usize {
        self.my_surf.nb_u_poles()
    }
    pub fn nb_v_poles(&self) -> usize {
        self.my_surf.nb_v_poles()
    }
    pub fn nb_u_knots(&self) -> usize {
        self.my_surf.nb_u_knots()
    }
    pub fn nb_v_knots(&self) -> usize {
        self.my_surf.nb_v_knots()
    }
    pub fn u_degree(&self) -> usize {
        self.my_surf.u_degree()
    }
    pub fn v_degree(&self) -> usize {
        self.my_surf.v_degree()
    }
    pub fn is_u_periodic(&self) -> bool {
        self.my_surf.is_u_periodic()
    }
    pub fn u_period(&self) -> f64 {
        self.my_surf.u_period()
    }
    pub fn is_v_periodic(&self) -> bool {
        self.my_surf.is_v_periodic()
    }
    pub fn v_period(&self) -> f64 {
        self.my_surf.v_period()
    }
    pub fn is_u_closed(&self) -> bool {
        matches!(
            self.my_surf.get_type(),
            GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Torus
        )
    }
    pub fn is_v_closed(&self) -> bool {
        self.my_surf.get_type() == GeomAbsSurfaceType::Torus
    }
    /// OCCT BRepAdaptor_Surface::UTrim — the narrowed window on the same
    /// surface (GeomAdaptor semantics).
    pub fn u_trim(&self, first: f64, last: f64, _tol: f64) -> BRepAdaptorSurface<'a> {
        let d = self.my_surf.uv_window();
        BRepAdaptorSurface {
            my_face: self.my_face.clone(),
            my_surf: self.my_surf.clone_with_window([first.min(last), last.max(first), d[2], d[3]]),
            _brep: std::marker::PhantomData,
        }
    }
    /// OCCT BRepAdaptor_Surface::VTrim.
    pub fn v_trim(&self, first: f64, last: f64, _tol: f64) -> BRepAdaptorSurface<'a> {
        let d = self.my_surf.uv_window();
        BRepAdaptorSurface {
            my_face: self.my_face.clone(),
            my_surf: self.my_surf.clone_with_window([d[0], d[1], first.min(last), last.max(first)]),
            _brep: std::marker::PhantomData,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{CylindricalSurface, Line3};
    use rcad_kernel::topo::topods::BRepBuilder;

    /// OCCT anchor: Initialize(F, true) restricts the cylinder window to
    /// the face UV bounds; GetType stays Cylinder; FirstU/LastU return the
    /// face window.
    #[test]
    fn brep_adaptor_surface_cylinder_window() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, Point3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, Point3::new(0.0, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin: Point3::ZERO,
                direction: Vec3::new(0.0, 0.0, 1.0),
            })),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep
            .add_tface(
                Some(Surface3::Cylinder(CylindricalSurface {
                    origin: Point3::ZERO,
                    axis: Vec3::new(0.0, 0.0, 1.0),
                    radius: 1.5,
                    ref_dir: Vec3::new(1.0, 0.0, 0.0),
                    y_dir: None,
                })),
                wire,
                Vec::new(),
                None,
                Some([0.0, std::f64::consts::TAU, 0.0, 2.0]),
                Vec::new(),
                true,
            )
;

        let s = BRepAdaptorSurface::initialize_face(&brep, &face, true);
        assert_eq!(s.get_type(), GeomAbsSurfaceType::Cylinder);
        assert!((s.first_u_parameter() - 0.0).abs() < 1e-12);
        assert!((s.last_u_parameter() - std::f64::consts::TAU).abs() < 1e-12);
        assert!((s.first_v_parameter() - 0.0).abs() < 1e-12);
        assert!((s.last_v_parameter() - 2.0).abs() < 1e-12);
        // Value on the located surface: P(0, 1) = (R, 0, 1).
        let p = s.value(0.0, 1.0);
        assert!((p.x - 1.5).abs() < 1e-12 && p.y.abs() < 1e-12 && (p.z - 1.0).abs() < 1e-12);
    }
}
