//! OCCT ProjLib::Project free functions (TKGeomBase/ProjLib, ProjLib.cxx).
//!
//! 1:1 translation of ProjLib.cxx L50-233: the Project overloads over the
//! analytic surfaces (plane, cylinder, cone, sphere, torus), MakePCurveOfType
//! and IsAnaSurf.
//!
//! Architecture difference: the OCCT overloads construct a
//! ProjLib_Plane/Cylinder/Cone/Sphere/Torus member and return its result
//! (`ProjLib_Plane Proj(Pl, L); return Proj.Line();`).  rcad's members of the
//! ProjLib package are the projector structs of [`super`] (ProjLib_Projector /
//! ProjLib_Plane / ... encoding); the overloads below route through them and
//! unpack the 2D-embedded result into the matching 2D geometry type.  The
//! per-type Project bodies themselves are the members' own translations
//! (staged for 1:1 rewrite where the legacy member deviates).

use glam::{DVec2, DVec3};

use super::adaptor::{Adaptor3dSurface, GeomAbsSurfaceType};
use super::{
    ConeProjector, CylinderProjector, CurveType, PlaneProjector, Projector, SphereProjector,
    TorusProjector,
};
use crate::geom::{
    Circle2d, Circle3, ConicalSurface, Curve2d, CylindricalSurface, Ellipse2d, Ellipse3,
    Hyperbola2d, Hyperbola3, Line2d, Line3, Parabola2d, Parabola3, Plane, SphericalSurface,
    ToroidalSurface,
};

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Pnt& P) (ProjLib.cxx
/// L50-55): ElSLib::Parameters(Pl, P, U, V).
pub fn project_pln_pnt(pl: &Plane, p: DVec3) -> DVec2 {
    // OCCT: ElSLib::Parameters(Pl, P, U, V) — U = (P - Loc) . XDirection,
    // V = (P - Loc) . YDirection.
    let d = p - pl.origin;
    DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir))
}

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Lin& L) (L59-63).
pub fn project_pln_lin(pl: &Plane, l: &Line3) -> Line2d {
    let mut proj = PlaneProjector::with_plane(pl);
    proj.project_line(l);
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Circ& C) (L67-71).
pub fn project_pln_circ(pl: &Plane, c: &Circle3) -> Circle2d {
    let mut proj = PlaneProjector::with_plane(pl);
    proj.project_circle(c);
    as_circle2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Elips& E) (L75-79).
pub fn project_pln_elips(pl: &Plane, e: &Ellipse3) -> Ellipse2d {
    let mut proj = PlaneProjector::with_plane(pl);
    proj.project_ellipse(e);
    as_ellipse2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Parab& P) (L83-87).
pub fn project_pln_parab(pl: &Plane, p: &Parabola3) -> Parabola2d {
    let mut proj = PlaneProjector::with_plane(pl);
    proj.project_parabola(p);
    // OCCT returns Proj.Parabola() (gp_Parab2d); the projector stores the
    // 2D-embedded parabola (vertex in plane UV, axis direction in plane UV).
    let parab = proj.projector().parabola();
    Parabola2d {
        origin: DVec2::new(parab.vertex.x, parab.vertex.y),
        axis_dir: DVec2::new(parab.axis_dir.x, parab.axis_dir.y).normalize_or_zero(),
        focal_param: parab.focal_param,
    }
}

/// OCCT ProjLib::Project(const gp_Pln& Pl, const gp_Hypr& H) (L91-95).
pub fn project_pln_hypr(pl: &Plane, h: &Hyperbola3) -> Hyperbola2d {
    let mut proj = PlaneProjector::with_plane(pl);
    proj.project_hyperbola(h);
    // OCCT returns Proj.Hyperbola() (gp_Hypr2d); the projector stores the
    // 2D-embedded hyperbola (center in plane UV, major axis in plane UV).
    let hypr = proj.projector().hyperbola();
    Hyperbola2d {
        center: DVec2::new(hypr.center.x, hypr.center.y),
        major_dir: DVec2::new(hypr.major_dir.x, hypr.major_dir.y).normalize_or_zero(),
        semi_major: hypr.semi_major,
        semi_minor: hypr.semi_minor,
    }
}

/// OCCT ProjLib::Project(const gp_Cylinder& Cy, const gp_Pnt& P) (L99-104):
/// ElSLib::Parameters(Cy, P, U, V).
pub fn project_cylinder_pnt(cy: &CylindricalSurface, p: DVec3) -> DVec2 {
    cy.world_to_uv(p)
}

/// OCCT ProjLib::Project(const gp_Cylinder& Cy, const gp_Lin& L) (L108-112).
pub fn project_cylinder_lin(cy: &CylindricalSurface, l: &Line3) -> Line2d {
    let mut proj = CylinderProjector::with_cylinder(cy);
    proj.project_line(l);
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Cylinder& Cy, const gp_Circ& Ci) (L116-120)
/// — OCCT returns Proj.Line() (the v-parameter line of the circle's plane).
pub fn project_cylinder_circ(cy: &CylindricalSurface, ci: &Circle3) -> Line2d {
    let mut proj = CylinderProjector::with_cylinder(cy);
    proj.project_circle(ci);
    // OCCT: ProjLib_Cylinder::Project(gp_Circ) stores a line result; the
    // legacy member encodes the general case, so fall back to its Line arm
    // (the OCCT analytic case) through the shared unpack.
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Cone& Co, const gp_Pnt& P) (L124-129):
/// ElSLib::Parameters(Co, P, U, V).
pub fn project_cone_pnt(co: &ConicalSurface, p: DVec3) -> DVec2 {
    co.world_to_uv(p)
}

/// OCCT ProjLib::Project(const gp_Cone& Co, const gp_Lin& L) (L133-137).
pub fn project_cone_lin(co: &ConicalSurface, l: &Line3) -> Line2d {
    let mut proj = ConeProjector::with_cone(co);
    proj.project_line(l);
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Cone& Co, const gp_Circ& Ci) (L141-145) —
/// OCCT returns Proj.Line().
pub fn project_cone_circ(co: &ConicalSurface, ci: &Circle3) -> Line2d {
    let mut proj = ConeProjector::with_cone(co);
    proj.project_circle(ci);
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Sphere& Sp, const gp_Pnt& P) (L149-154):
/// ElSLib::Parameters(Sp, P, U, V).
pub fn project_sphere_pnt(sp: &SphericalSurface, p: DVec3) -> DVec2 {
    sp.world_to_uv(p)
}

/// OCCT ProjLib::Project(const gp_Sphere& Sp, const gp_Circ& Ci) (L158-162) —
/// OCCT returns Proj.Line().
pub fn project_sphere_circ(sp: &SphericalSurface, ci: &Circle3) -> Line2d {
    let mut proj = SphereProjector::with_sphere(sp);
    proj.project_circle(ci);
    as_line2d(&proj)
}

/// OCCT ProjLib::Project(const gp_Torus& To, const gp_Pnt& P) (L166-171):
/// ElSLib::Parameters(To, P, U, V).
pub fn project_torus_pnt(to: &ToroidalSurface, p: DVec3) -> DVec2 {
    to.world_to_uv(p)
}

/// OCCT ProjLib::Project(const gp_Torus& To, const gp_Circ& Ci) (L175-179) —
/// OCCT returns Proj.Line().
pub fn project_torus_circ(to: &ToroidalSurface, ci: &Circle3) -> Line2d {
    let mut proj = TorusProjector::with_torus(to);
    proj.project_circle(ci);
    as_line2d(&proj)
}

/// OCCT ProjLib::MakePCurveOfType(const ProjLib_ProjectedCurve& PC,
/// handle(Geom2d_Curve)& C2D) (L183-213): builds the 2D curve of the type
/// carried by the projection result.  Bezier / Other throw
/// Standard_NotImplemented (mirrored as a panic).
///
/// Architecture difference: OCCT consumes a ProjLib_ProjectedCurve; rcad's
/// projection result carrier is the [`Projector`] (the ProjLib_Projector
/// encoding holding the type + result).
pub fn make_pcurve_of_type(pc: &Projector) -> Option<Curve2d> {
    match pc.get_type() {
        CurveType::Line => Some(Curve2d::Line(as_line2d_of(pc))),
        CurveType::Circle => Some(Curve2d::Circle(as_circle2d_of(pc))),
        CurveType::Ellipse => Some(Curve2d::Ellipse(as_ellipse2d_of(pc))),
        CurveType::Parabola => {
            let parab = pc.parabola();
            Some(Curve2d::Parabola(Parabola2d {
                origin: DVec2::new(parab.vertex.x, parab.vertex.y),
                axis_dir: DVec2::new(parab.axis_dir.x, parab.axis_dir.y).normalize_or_zero(),
                focal_param: parab.focal_param,
            }))
        }
        CurveType::Hyperbola => {
            let hypr = pc.hyperbola();
            Some(Curve2d::Hyperbola(Hyperbola2d {
                center: DVec2::new(hypr.center.x, hypr.center.y),
                major_dir: DVec2::new(hypr.major_dir.x, hypr.major_dir.y).normalize_or_zero(),
                semi_major: hypr.semi_major,
                semi_minor: hypr.semi_minor,
            }))
        }
        CurveType::BSpline => pc.bspline().cloned(),
        // OCCT: case GeomAbs_BezierCurve / GeomAbs_OtherCurve / default:
        //   throw Standard_NotImplemented("ProjLib::MakePCurveOfType");
        CurveType::Bezier | CurveType::Other => {
            panic!("Standard_NotImplemented: ProjLib::MakePCurveOfType")
        }
    }
}

/// OCCT ProjLib::IsAnaSurf(const handle(Adaptor3d_Surface)& theAS)
/// (L217-233): true when the surface type is one of the analytic quadrics.
pub fn is_ana_surf(the_as: &dyn Adaptor3dSurface) -> bool {
    matches!(
        the_as.get_type(),
        GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus
    )
}

// ---------------------------------------------------------------------------
// Result unpack helpers — the projector stores 2D-embedded 3D forms; the
// ProjLib::Project overloads return the matching 2D types.
// ---------------------------------------------------------------------------

fn as_line2d_of(proj: &Projector) -> Line2d {
    let l = proj.line();
    Line2d {
        origin: DVec2::new(l.origin.x, l.origin.y),
        direction: DVec2::new(l.direction.x, l.direction.y).normalize_or_zero(),
    }
}

fn as_circle2d_of(proj: &Projector) -> Circle2d {
    let c = proj.circle();
    Circle2d {
        center: DVec2::new(c.center.x, c.center.y),
        x_dir: DVec2::new(c.x_dir.x, c.x_dir.y),
        y_dir: DVec2::new(c.y_dir.x, c.y_dir.y),
        radius: c.radius,
    }
}

fn as_ellipse2d_of(proj: &Projector) -> Ellipse2d {
    let e = proj.ellipse();
    Ellipse2d {
        center: DVec2::new(e.center.x, e.center.y),
        major_dir: DVec2::new(e.major_dir.x, e.major_dir.y).normalize_or_zero(),
        minor_dir: DVec2::new(-e.major_dir.y, e.major_dir.x),
        major_radius: e.major_radius,
        minor_radius: e.minor_radius,
    }
}

fn as_line2d<P>(proj: &P) -> Line2d
where
    P: ProjectorAccess,
{
    as_line2d_of(proj.projector())
}

fn as_circle2d<P>(proj: &P) -> Circle2d
where
    P: ProjectorAccess,
{
    as_circle2d_of(proj.projector())
}

fn as_ellipse2d<P>(proj: &P) -> Ellipse2d
where
    P: ProjectorAccess,
{
    as_ellipse2d_of(proj.projector())
}

/// Access to the projection result common to the ProjLib_Plane / Cylinder /
/// Cone / Sphere / Torus members (the OCCT members all expose the
/// ProjLib_Projector base result).
trait ProjectorAccess {
    fn projector(&self) -> &Projector;
}

impl ProjectorAccess for PlaneProjector {
    fn projector(&self) -> &Projector {
        self.projector()
    }
}

impl ProjectorAccess for CylinderProjector {
    fn projector(&self) -> &Projector {
        self.projector()
    }
}

impl ProjectorAccess for ConeProjector {
    fn projector(&self) -> &Projector {
        self.projector()
    }
}

impl ProjectorAccess for SphereProjector {
    fn projector(&self) -> &Projector {
        self.projector()
    }
}

impl ProjectorAccess for TorusProjector {
    fn projector(&self) -> &Projector {
        self.projector()
    }
}
