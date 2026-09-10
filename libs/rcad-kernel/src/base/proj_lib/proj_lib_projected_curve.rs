//! ProjLib_ProjectedCurve support translation (TKGeomBase/ProjLib) — part A:
//! the analytic ProjLib_Plane / Cylinder / Cone / Sphere / Torus members and
//! the ProjLib_ProjectedCurve.cxx static helpers.
//!
//! Sources (OCCT 8.0):
//!   - ProjLib_ProjectedCurve.cxx L52-159 (static helpers) + L242-270 (the
//!     Project dispatch)
//!   - ProjLib_Plane.cxx L26-169, ProjLib_Cylinder.cxx L29-175,
//!     ProjLib_Cone.cxx L32-180, ProjLib_Sphere.cxx L36-246,
//!     ProjLib_Torus.cxx L30-193, ProjLib_Projector.cxx L124-148 (the base
//!     Project bodies) + L166-251 (UFrame / VFrame)
//!
//! The consumed GeomAdaptor encodings (`GeomCurveAdaptor` /
//! `GeomSurfaceAdaptor` and their geometry-accessor traits) are the canonical
//! 1:1 translations of OCCT GeomAdaptor_Curve / GeomAdaptor_Surface that live
//! in the `geom_adaptor_curve` / `geom_adaptor_surface` modules; they stay
//! importable from this module path (the original consumer home).
//!
//! Architecture differences (rcad encodings):
//!   - OCCT ProjLib_* members derive from ProjLib_Projector and hold the
//!     gp_* 2D results (myLin / myCirc / ...).  rcad keeps the module's
//!     established encoding: the base result is the [`Projector`] struct of
//!     [`super`] with the 2D results 2D-embedded in its 3D forms (gp_Lin2d
//!     -> Line3 with the x/y payload; gp_Circ2d -> Circle3 with the 2D frame
//!     in x_dir/y_dir and the z normal) — see `project.rs` header.  The
//!     `myResult = P` copies of Perform translate as struct moves.
//!   - OCCT gp_Pln/gp_Cylinder/gp_Sphere/gp_Cone/gp_Torus payloads map onto
//!     the kernel analytic surface structs (Plane / CylindricalSurface /
//!     SphericalSurface / ConicalSurface / ToroidalSurface).  rcad
//!     `ConicalSurface::apex` is the gp_Cone Location() reference point
//!     (radius == RefRadius there); the OCCT Apex() is
//!     [`ConicalSurface::apex_point`].

use glam::{DVec2, DVec3};

use super::adaptor::Adaptor3dSurface;
use super::{CurveType, Projector};
use crate::core::precision;
use crate::geom::{
    Circle3, ConicalSurface, Curve2d, CylindricalSurface, Ellipse3, Hyperbola3, Line3, Parabola3,
    Plane, SphericalSurface, ToroidalSurface,
};

pub(crate) const TWO_PI: f64 = std::f64::consts::TAU;

// =========================================================================
// OCCT gp_Dir / gp_Lin2d small helpers (gp_Dir.hxx / gp_Trsf2d semantics)
// =========================================================================

/// OCCT gp_Dir::AngleWithRef(TheOther, VRef) — the angle in [0, 2*pi) of
/// `a` towards `b`, counted positive around the reference direction.
pub(crate) fn dir_angle_with_ref(a: DVec3, b: DVec3, v_ref: DVec3) -> f64 {
    let cross = a.cross(b);
    let ang = a.dot(b).clamp(-1.0, 1.0).acos();
    if v_ref.dot(cross) >= 0.0 {
        ang
    } else {
        TWO_PI - ang
    }
}

/// OCCT gp_Dir::IsNormal(TheOther, AngularTolerance) — |dot| <= tolerance.
fn dir_is_normal(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    a.dot(b).abs() <= angular_tolerance
}

/// OCCT gp_Dir::IsParallel(TheOther, AngularTolerance) — the cross modulus
/// of the unit directions is compared with the angular tolerance.
fn dir_is_parallel(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    a.cross(b).length() <= angular_tolerance
}

/// OCCT gp_Dir2d::IsEqual(Other, AngularTolerance).
fn dir2d_is_equal(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    a.dot(b) >= 1.0 - angular_tolerance
}

/// OCCT gp_Dir2d::IsOpposite(Other, AngularTolerance).
fn dir2d_is_opposite(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    a.dot(b) <= -1.0 + angular_tolerance
}

/// The 2D line payload of the [`Projector`] result (gp_Lin2d encoding:
/// origin/direction in the x/y slots).
#[inline]
pub(crate) fn proj_line_location(line: &Line3) -> DVec2 {
    DVec2::new(line.origin.x, line.origin.y)
}

#[inline]
pub(crate) fn proj_line_direction(line: &Line3) -> DVec2 {
    DVec2::new(line.direction.x, line.direction.y)
}

/// OCCT gp_Lin2d::Value(U) — Location + U * Direction.
pub(crate) fn proj_line_value(line: &Line3, u: f64) -> DVec2 {
    proj_line_location(line) + proj_line_direction(line) * u
}

/// OCCT ElCLib::Value(U, gp_Lin2d) — the same evaluation on the result
/// payload.
pub(crate) fn elclib_value_lin2d(line: &Line3, u: f64) -> DVec2 {
    proj_line_value(line, u)
}

/// OCCT gp_Lin2d::Translate(gp_Vec2d) on the result payload.
pub(crate) fn proj_line_translate(line: Line3, t: DVec2) -> Line3 {
    Line3::new(line.origin + DVec3::new(t.x, t.y, 0.0), line.direction)
}

/// OCCT gp_Trsf2d SetMirror(gp_Ax2d) applied to a gp_Lin2d — the reflection
/// of the location and the direction about the 2D axis.
fn mirror_lin2d_about_axis(line: Line3, axis_loc: DVec2, axis_dir: DVec2) -> Line3 {
    let loc = proj_line_location(&line);
    let dir = proj_line_direction(&line);
    let d = loc - axis_loc;
    let new_loc = axis_loc + axis_dir * (2.0 * d.dot(axis_dir)) - d;
    let new_dir = axis_dir * (2.0 * dir.dot(axis_dir)) - dir;
    Line3::new(
        DVec3::new(new_loc.x, new_loc.y, 0.0),
        DVec3::new(new_dir.x, new_dir.y, 0.0),
    )
}

// =========================================================================
// OCCT GeomAdaptor_Curve / GeomAdaptor_Surface (TKG3d) — the canonical 1:1
// translations live in the `geom_adaptor_curve` / `geom_adaptor_surface`
// modules; the consumed names stay importable from this module path (the
// original consumer home).
// =========================================================================

pub use super::geom_adaptor_curve::{Adaptor3dCurveGeom, GeomCurveAdaptor};
pub use super::geom_adaptor_surface::{Adaptor3dSurfaceGeom, GeomSurfaceAdaptor};

// =========================================================================
// OCCT GeomAdaptor_Surface — canonical 1:1 lives in `geom_adaptor_surface`.
// =========================================================================

/// Field copy of the [`Projector`] base result (the `myResult = P` move of
/// Perform; `Projector` itself carries no derive).
pub(crate) fn clone_projector(p: &Projector) -> Projector {
    Projector {
        proj_type: p.proj_type,
        lin: p.lin,
        circ: p.circ,
        elips: p.elips,
        hypr: p.hypr,
        parab: p.parab,
        bspline: p.bspline.clone(),
        is_periodic: p.is_periodic,
        is_done: p.is_done,
    }
}

macro_rules! impl_projlib_clone {
    ($t:ident, $base_field:ident) => {
        impl Clone for $t {
            fn clone(&self) -> Self {
                $t {
                    projector: clone_projector(&self.projector),
                    $base_field: self.$base_field,
                }
            }
        }
    };
}

impl_projlib_clone!(ProjLibPlane, my_plane);
impl_projlib_clone!(ProjLibCylinder, my_cylinder);
impl_projlib_clone!(ProjLibCone, my_cone);
impl_projlib_clone!(ProjLibSphere, my_sphere);
impl_projlib_clone!(ProjLibTorus, my_torus);

// =========================================================================
// OCCT ProjLib_Projector base Project bodies (ProjLib_Projector.cxx L124-148)
// =========================================================================

/// The OCCT ProjLib_Projector::Project overloads — the base bodies set the
/// result type to OtherCurve (the derived member overrides where an analytic
/// answer exists).  Shared by every ProjLib_* member through
/// [`ProjLibProjectorOverloads`].
pub trait ProjLibProjectorOverloads {
    /// The ProjLib_Projector base result.
    fn projector(&mut self) -> &mut Projector;

    /// OCCT ProjLib_Projector::Project(const gp_Lin&) (L124-128).
    fn project_lin(&mut self, _l: &Line3) {
        self.projector().set_type(CurveType::Other);
    }

    /// OCCT ProjLib_Projector::Project(const gp_Circ&) (L130-134).
    fn project_circ(&mut self, _c: &Circle3) {
        self.projector().set_type(CurveType::Other);
    }

    /// OCCT ProjLib_Projector::Project(const gp_Elips&) (L136-140).
    fn project_elips(&mut self, _e: &Ellipse3) {
        self.projector().set_type(CurveType::Other);
    }

    /// OCCT ProjLib_Projector::Project(const gp_Parab&) (L142-146).
    fn project_parab(&mut self, _p: &Parabola3) {
        self.projector().set_type(CurveType::Other);
    }

    /// OCCT ProjLib_Projector::Project(const gp_Hypr&) (L148).
    fn project_hypr(&mut self, _h: &Hyperbola3) {
        self.projector().set_type(CurveType::Other);
    }
}

// =========================================================================
// OCCT ProjLib_Plane (ProjLib_Plane.cxx L26-169)
// =========================================================================

/// OCCT ProjLib_Plane — the analytic projector onto a plane.
pub struct ProjLibPlane {
    /// OCCT ProjLib_Projector base.
    pub projector: Projector,
    /// OCCT: gp_Pln myPlane.
    pub my_plane: Plane,
}

/// OCCT ProjLib_Plane.cxx L86-92 — static EvalPnt2d(P, Pl): the plane UV of
/// the 3D point.
fn eval_pnt2d_plane(p: DVec3, pl: &Plane) -> DVec2 {
    let op = p - pl.origin;
    DVec2::new(op.dot(pl.u_dir), op.dot(pl.v_dir))
}

/// OCCT ProjLib_Plane.cxx L94-99 — static EvalDir2d(D, Pl): the plane UV of
/// the 3D direction.
fn eval_dir2d_plane(d: DVec3, pl: &Plane) -> DVec2 {
    DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir))
}

impl ProjLibPlane {
    /// OCCT ProjLib_Plane(const gp_Pln& Pl) -> Init(Pl) (L30-34, L77-84).
    pub fn new(pl: &Plane) -> Self {
        ProjLibPlane {
            projector: Projector::new(),
            my_plane: *pl,
        }
    }

    /// OCCT Init(Pl) (L77-84): myType = OtherCurve, isDone = false,
    /// myIsPeriodic = false, myPlane = Pl.
    pub fn init(&mut self, pl: &Plane) {
        self.projector = Projector::new();
        self.my_plane = *pl;
    }
}

/// The OCCT ProjLib_Plane::Project overloads (the virtual overrides of the
/// ProjLib_Projector base).
impl ProjLibProjectorOverloads for ProjLibPlane {
    fn projector(&mut self) -> &mut Projector {
        &mut self.projector
    }

    /// OCCT ProjLib_Plane::Project(const gp_Lin& L) (L101-107).
    fn project_lin(&mut self, l: &Line3) {
        let p2d = eval_pnt2d_plane(l.origin, &self.my_plane);
        let d2d = eval_dir2d_plane(l.direction, &self.my_plane).normalize_or_zero();
        self.projector.lin = Line3::new(
            DVec3::new(p2d.x, p2d.y, 0.0),
            DVec3::new(d2d.x, d2d.y, 0.0),
        );
        self.projector.proj_type = CurveType::Line;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Plane::Project(const gp_Circ& C) (L110-125).
    fn project_circ(&mut self, c: &Circle3) {
        let p2d = eval_pnt2d_plane(c.center, &self.my_plane);
        let x2d = eval_dir2d_plane(c.x_dir, &self.my_plane).normalize_or_zero();
        let y2d = eval_dir2d_plane(c.y_dir, &self.my_plane).normalize_or_zero();
        // OCCT: myCirc = gp_Circ2d(gp_Ax22d(P2d, X2d, Y2d), C.Radius()).
        self.projector.circ = Circle3 {
            center: DVec3::new(p2d.x, p2d.y, 0.0),
            normal: DVec3::Z,
            x_dir: DVec3::new(x2d.x, x2d.y, 0.0),
            y_dir: DVec3::new(y2d.x, y2d.y, 0.0),
            radius: c.radius,
        };
        self.projector.proj_type = CurveType::Circle;
        // OCCT: myIsPeriodic = true.
        self.projector.is_periodic = true;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Plane::Project(const gp_Elips& E) (L127-141).
    fn project_elips(&mut self, e: &Ellipse3) {
        let p2d = eval_pnt2d_plane(e.center, &self.my_plane);
        let x2d = eval_dir2d_plane(e.major_dir, &self.my_plane).normalize_or_zero();
        // OCCT: Y2d = EvalDir2d(E.Position().YDirection()) — the minor axis
        // (rcad Ellipse3 carries normal + major_dir; YDir = ZDir ^ XDir).
        let y_dir = e.normal.cross(e.major_dir);
        let y2d = eval_dir2d_plane(y_dir, &self.my_plane).normalize_or_zero();
        // OCCT: myElips = gp_Elips2d(gp_Ax22d(P2d, X2d, Y2d), Major, Minor).
        self.projector.elips = Ellipse3 {
            center: DVec3::new(p2d.x, p2d.y, 0.0),
            normal: DVec3::Z,
            major_dir: DVec3::new(x2d.x, x2d.y, 0.0),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        };
        let _ = y2d;
        self.projector.proj_type = CurveType::Ellipse;
        self.projector.is_periodic = true;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Plane::Project(const gp_Parab& P) (L143-156).
    fn project_parab(&mut self, p: &Parabola3) {
        let p2d = eval_pnt2d_plane(p.vertex, &self.my_plane);
        let x2d = eval_dir2d_plane(p.axis_dir, &self.my_plane).normalize_or_zero();
        // OCCT: myParab = gp_Parab2d(gp_Ax22d(P2d, X2d, Y2d), P.Focal()).
        self.projector.parab = Parabola3 {
            vertex: DVec3::new(p2d.x, p2d.y, 0.0),
            normal: DVec3::Z,
            axis_dir: DVec3::new(x2d.x, x2d.y, 0.0),
            focal_param: p.focal_param,
        };
        self.projector.proj_type = CurveType::Parabola;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Plane::Project(const gp_Hypr& H) (L158-169).
    fn project_hypr(&mut self, h: &Hyperbola3) {
        let p2d = eval_pnt2d_plane(h.center, &self.my_plane);
        let x2d = eval_dir2d_plane(h.major_dir, &self.my_plane).normalize_or_zero();
        // OCCT: myHypr = gp_Hypr2d(gp_Ax22d(P2d, X2d, Y2d), Major, Minor).
        self.projector.hypr = Hyperbola3 {
            center: DVec3::new(p2d.x, p2d.y, 0.0),
            normal: DVec3::Z,
            major_dir: DVec3::new(x2d.x, x2d.y, 0.0),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        };
        self.projector.proj_type = CurveType::Hyperbola;
        self.projector.is_done = true;
    }
}

// =========================================================================
// OCCT ProjLib_Cylinder (ProjLib_Cylinder.cxx L29-175)
// =========================================================================

/// OCCT ProjLib_Cylinder — the analytic projector onto a cylinder.
pub struct ProjLibCylinder {
    /// OCCT ProjLib_Projector base.
    pub projector: Projector,
    /// OCCT: gp_Cylinder myCylinder.
    pub my_cylinder: CylindricalSurface,
}

impl ProjLibCylinder {
    /// OCCT ProjLib_Cylinder(const gp_Cylinder& Cyl) -> Init(Cyl)
    /// (L33-37, L64-71).
    pub fn new(cyl: &CylindricalSurface) -> Self {
        ProjLibCylinder {
            projector: Projector::new(),
            my_cylinder: *cyl,
        }
    }

    /// OCCT Init(Cyl) (L64-71).
    pub fn init(&mut self, cyl: &CylindricalSurface) {
        self.projector = Projector::new();
        self.my_cylinder = *cyl;
    }

    /// OCCT ProjLib_Cylinder.cxx L78-85 — static EvalPnt2d(P, Cy).
    fn eval_pnt2d(&self, p: DVec3) -> DVec2 {
        let op = p - self.my_cylinder.origin;
        let x = op.dot(self.my_cylinder.ref_dir);
        let y = op.dot(self.my_cylinder.y_axis());
        let z = op.dot(self.my_cylinder.axis);
        let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
            y.atan2(x)
        } else {
            0.0
        };
        DVec2::new(u, z)
    }
}

/// The OCCT ProjLib_Cylinder::Project overloads.
impl ProjLibProjectorOverloads for ProjLibCylinder {
    fn projector(&mut self) -> &mut Projector {
        &mut self.projector
    }

    /// OCCT ProjLib_Cylinder::Project(const gp_Lin& L) (L95-120).
    fn project_lin(&mut self, l: &Line3) {
        // OCCT L97-102: the line must be parallel to the cylinder axis.
        if l.direction.cross(self.my_cylinder.axis).length_squared()
            > precision::ANGULAR * precision::ANGULAR
        {
            return;
        }
        let mut p2d = self.eval_pnt2d(l.origin);
        if p2d.x < 0.0 {
            p2d.x += TWO_PI;
        }
        // OCCT L110-111: Signe = L.Direction().Dot(axis) > 0 ? 1 : -1.
        let signe = if l.direction.dot(self.my_cylinder.axis) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let d2d = DVec2::new(0.0, signe);
        self.projector.lin = Line3::new(
            DVec3::new(p2d.x, p2d.y, 0.0),
            DVec3::new(d2d.x, d2d.y, 0.0),
        );
        self.projector.proj_type = CurveType::Line;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Cylinder::Project(const gp_Circ& C) (L122-159).
    fn project_circ(&mut self, c: &Circle3) {
        // OCCT L124-131: the circle's normal must be parallel to the axis.
        if self.my_cylinder.axis.cross(c.normal).length_squared()
            > precision::ANGULAR * precision::ANGULAR
        {
            return;
        }
        // OCCT L135: ZCyl = XDirection ^ YDirection.
        let z_cyl = self.my_cylinder.ref_dir.cross(self.my_cylinder.y_axis()).normalize();
        // OCCT L137: U = XDirection.AngleWithRef(aCircPos.XDirection(), ZCyl).
        let u = dir_angle_with_ref(self.my_cylinder.ref_dir, c.x_dir, z_cyl);
        // OCCT L139-140: V = (Location -> C.Location()) . Direction.
        let op = c.center - self.my_cylinder.origin;
        let v = op.dot(self.my_cylinder.axis);
        // OCCT L144-149.
        let d2d = if z_cyl.dot(c.normal) > 0.0 {
            DVec2::new(1.0, 0.0)
        } else {
            DVec2::new(-1.0, 0.0)
        };
        self.projector.lin = Line3::new(
            DVec3::new(u, v, 0.0),
            DVec3::new(d2d.x, d2d.y, 0.0),
        );
        self.projector.proj_type = CurveType::Line;
        self.projector.is_done = true;
    }

    /// OCCT ProjLib_Cylinder::Project(const gp_Elips&) (L161-165) — delegated
    /// to approximation: the projector stays NOT done (the base body is not
    /// called).
    fn project_elips(&mut self, _e: &Ellipse3) {}

    // OCCT L167-173: Project(gp_Parab) / Project(gp_Hypr) call the
    // ProjLib_Projector base bodies — the trait defaults.
}

// =========================================================================
// OCCT ProjLib_Cone (ProjLib_Cone.cxx L32-180)
// =========================================================================

/// OCCT ProjLib_Cone — the analytic projector onto a cone.
pub struct ProjLibCone {
    /// OCCT ProjLib_Projector base.
    pub projector: Projector,
    /// OCCT: gp_Cone myCone.
    pub my_cone: ConicalSurface,
}

/// The cone frame Y direction (gp_Ax3 YDirection of the cone position).
fn cone_y_dir(cone: &ConicalSurface) -> DVec3 {
    cone.ref_dir.cross(cone.axis).normalize()
}

impl ProjLibCone {
    /// OCCT ProjLib_Cone(const gp_Cone& Co) -> Init(Co) (L36-40, L59-66).
    pub fn new(co: &ConicalSurface) -> Self {
        ProjLibCone {
            projector: Projector::new(),
            my_cone: *co,
        }
    }

    /// OCCT Init(Co) (L59-66).
    pub fn init(&mut self, co: &ConicalSurface) {
        self.projector = Projector::new();
        self.my_cone = *co;
    }

    /// OCCT ElSLib::ConeD1 (ElSLib.cxx) for the cone frame — the point and
    /// the partials at (U, V).
    fn cone_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let pos = &self.my_cone;
        let x_dir = pos.ref_dir;
        let y_dir = cone_y_dir(pos);
        let z_dir = pos.axis;
        let cos_u = u.cos();
        let sin_u = u.sin();
        let cos_a = pos.half_angle_rad.cos();
        let sin_a = pos.half_angle_rad.sin();
        let r = pos.radius + v * sin_a;
        let a3 = v * cos_a;
        let a1 = r * cos_u;
        let a2 = r * sin_u;
        let r1 = sin_a * cos_u;
        let r2 = sin_a * sin_u;
        let p = pos.apex + x_dir * a1 + y_dir * a2 + z_dir * a3;
        let vu = x_dir * (-a2) + y_dir * a1;
        let vv = x_dir * r1 + y_dir * r2 + z_dir * cos_a;
        (p, vu, vv)
    }
}

/// The OCCT ProjLib_Cone::Project overloads.
impl ProjLibProjectorOverloads for ProjLibCone {
    fn projector(&mut self) -> &mut Projector {
        &mut self.projector
    }

    /// OCCT ProjLib_Cone::Project(const gp_Lin& L) (L69-106).
    fn project_lin(&mut self, l: &Line3) {
        let mut a_pnt = l.origin;
        let an_apex = self.my_cone.apex_point();
        let mut a_delta_v = 0.0_f64;

        if a_pnt.distance(an_apex) <= precision::CONFUSION {
            // OCCT L76-79: take another point of the line, off the apex.
            a_pnt += l.direction;
            a_delta_v = 1.0;
        }

        // OCCT L82: ElSLib::ConeParameters(Position, RefRadius, SemiAngle,
        // aPnt, U, V).
        let (u, v) = crate::math::el::elslib_cone_parameters(
            a_pnt,
            self.my_cone.apex,
            self.my_cone.ref_dir,
            cone_y_dir(&self.my_cone),
            self.my_cone.axis,
            self.my_cone.radius,
            self.my_cone.half_angle_rad,
        );

        // OCCT L85-89: ElSLib::ConeD1(U, V, ...).
        let (_p, _vu, vv) = self.cone_d1(u, v);

        // OCCT L91-104: the line is projected only when parallel to the
        // U-isoline of the cone.
        let dv = vv.normalize();
        if dir_is_parallel(dv, l.direction, precision::ANGULAR) {
            let a_sign = if l.direction.dot(dv) > 0.0 { 1.0 } else { -1.0 };
            let p2d = DVec2::new(u, v - a_delta_v * a_sign);
            let d2d = DVec2::new(0.0, a_sign);
            self.projector.lin = Line3::new(
                DVec3::new(p2d.x, p2d.y, 0.0),
                DVec3::new(d2d.x, d2d.y, 0.0),
            );
            self.projector.proj_type = CurveType::Line;
            self.projector.is_done = true;
        }
    }

    /// OCCT ProjLib_Cone::Project(const gp_Circ& C) (L108-165).
    fn project_circ(&mut self, c: &Circle3) {
        // OCCT L113-118: the circle's axis must be parallel to the cone axis.
        let cone_pos_dir = self.my_cone.axis;
        if !dir_is_parallel(cone_pos_dir, c.normal, precision::ANGULAR) {
            self.projector.is_done = false;
            return;
        }
        // OCCT L120-122: ZCone = XDirection ^ YDirection; ZCir = Xc ^ Yc.
        let z_cone = self.my_cone.ref_dir.cross(cone_y_dir(&self.my_cone)).normalize();
        let z_cir = c.x_dir.cross(c.y_dir).normalize();

        // OCCT L124-126: the local direction cosines of the circle location.
        let x = self.my_cone.ref_dir.dot(c.x_dir);
        let y = cone_y_dir(&self.my_cone).dot(c.x_dir);
        let z = (c.center - self.my_cone.apex).dot(self.my_cone.axis);

        // OCCT L128-147: the ElSLib-style U with the wrong-side guards.
        let mut u;
        if x.abs() <= precision::ANGULAR && y.abs() <= precision::ANGULAR {
            u = 0.0;
        } else if -self.my_cone.radius > z * self.my_cone.half_angle_rad.tan() {
            u = (-y).atan2(-x);
        } else {
            u = y.atan2(x);
        }
        if u < 0.0 {
            u += TWO_PI;
        }

        // OCCT L149: V = z / cos(SemiAngle).
        let v = z / self.my_cone.half_angle_rad.cos();

        // OCCT L151-160.
        let d2d = if z_cone.dot(z_cir) > 0.0 {
            DVec2::new(1.0, 0.0)
        } else {
            DVec2::new(-1.0, 0.0)
        };

        self.projector.lin = Line3::new(
            DVec3::new(u, v, 0.0),
            DVec3::new(d2d.x, d2d.y, 0.0),
        );
        self.projector.proj_type = CurveType::Line;
        self.projector.is_done = true;
    }

    // OCCT L167-177: Project(gp_Elips/gp_Parab/gp_Hypr) call the
    // ProjLib_Projector base bodies — the trait defaults.
}

// =========================================================================
// OCCT ProjLib_Sphere (ProjLib_Sphere.cxx L36-246)
// =========================================================================

/// OCCT ProjLib_Sphere — the analytic projector onto a sphere.
pub struct ProjLibSphere {
    /// OCCT ProjLib_Projector base.
    pub projector: Projector,
    /// OCCT: gp_Sphere mySphere.
    pub my_sphere: SphericalSurface,
}

impl ProjLibSphere {
    /// OCCT ProjLib_Sphere(const gp_Sphere& Sp) -> Init(Sp) (L40-44, L55-62).
    pub fn new(sp: &SphericalSurface) -> Self {
        ProjLibSphere {
            projector: Projector::new(),
            my_sphere: *sp,
        }
    }

    /// OCCT Init(Sp) (L55-62).
    pub fn init(&mut self, sp: &SphericalSurface) {
        self.projector = Projector::new();
        self.my_sphere = *sp;
    }

    /// OCCT ProjLib_Sphere.cxx L65-93 — static EvalPnt2d(P, Sp).
    fn eval_pnt2d(p: DVec3, sp: &SphericalSurface) -> DVec2 {
        let x = p.dot(sp.ref_dir);
        let y = p.dot(sp.ref_dir_perp());
        let z = p.dot(sp.axis);
        let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
            let uu = y.atan2(x);
            crate::math::el::in_period(uu, 0.0, TWO_PI)
        } else {
            0.0
        };
        let z = z.clamp(-1.0, 1.0);
        let v = z.asin();
        DVec2::new(u, v)
    }

    /// OCCT ProjLib_Sphere::SetInBounds(U) (L203-246) — place the projected
    /// line inside the sphere V bounds (mirroring about the pole axis when
    /// the line escapes [-pi/2, pi/2]).
    pub fn set_in_bounds(&mut self, u: f64) {
        // OCCT L204: StdFail_NotDone_Raise_if(!isDone, "ProjLib_Sphere:SetInBounds").
        assert!(self.projector.is_done, "ProjLib_Sphere:SetInBounds");

        // OCCT L207-210: first set the Y of the first point in [-pi, pi].
        let mut lin = self.projector.lin;
        let y = elclib_value_lin2d(&lin, u).y;
        let new_y = crate::math::el::in_period(y, -std::f64::consts::PI, std::f64::consts::PI);
        lin = proj_line_translate(lin, DVec2::new(0.0, new_y - y));

        // OCCT L211-215: re-evaluate the point and the direction.
        let p = proj_line_value(&lin, u);
        let tol = 1.0e-7;
        let d2d = proj_line_direction(&lin).normalize_or_zero();

        // OCCT L218-226: the north-pole mirror axis.
        let axis: Option<(DVec2, DVec2)> =
            if (p.y - std::f64::consts::FRAC_PI_2 > tol)
                || ((p.y - std::f64::consts::FRAC_PI_2).abs() < tol
                    && dir2d_is_equal(d2d, DVec2::Y, tol))
            {
                Some((DVec2::new(0.0, std::f64::consts::FRAC_PI_2), DVec2::X))
            } else if (p.y + std::f64::consts::FRAC_PI_2 < -tol)
                || ((p.y + std::f64::consts::FRAC_PI_2).abs() < tol
                    && dir2d_is_opposite(d2d, DVec2::Y, tol))
            {
                // OCCT L227-235: the south-pole mirror axis.
                Some((DVec2::new(0.0, -std::f64::consts::FRAC_PI_2), DVec2::X))
            } else {
                // OCCT L237-238: nothing to do.
                None
            };

        if let Some((axis_loc, axis_dir)) = axis {
            // OCCT L239-241: Trsf.SetMirror(Axis); myLin.Transform(Trsf).
            lin = mirror_lin2d_about_axis(lin, axis_loc, axis_dir);
            // OCCT L243: myLin.Translate(gp_Vec2d(M_PI, 0.)).
            lin = proj_line_translate(lin, DVec2::new(std::f64::consts::PI, 0.0));

            // OCCT L246-249: adjust the U parameter into [0, 2*pi).
            let x = proj_line_value(&lin, u).x;
            let new_x = crate::math::el::in_period(x, 0.0, TWO_PI);
            lin = proj_line_translate(lin, DVec2::new(new_x - x, 0.0));
        }
        self.projector.lin = lin;
    }
}

/// The OCCT ProjLib_Sphere::Project overloads.
impl ProjLibProjectorOverloads for ProjLibSphere {
    fn projector(&mut self) -> &mut Projector {
        &mut self.projector
    }

    /// OCCT ProjLib_Sphere::Project(const gp_Circ& C) (L97-178).
    #[allow(unused_assignments)] // OCCT L121 default-initializes P2d2 before the branches
    fn project_circ(&mut self, c: &Circle3) {
        // OCCT L99-113: the circle and sphere frames.
        let o = self.my_sphere.center;
        let xc = c.x_dir;
        let yc = c.y_dir;
        let zc = xc.cross(yc).normalize();
        let xs = self.my_sphere.ref_dir;
        let ys = self.my_sphere.ref_dir_perp();
        let zs = self.my_sphere.axis;

        // OCCT L115-119: isIsoU = Zc.IsNormal(Zs) && O == C.Location();
        // isIsoV = Xc.IsNormal(Zs) && Yc.IsNormal(Zs).
        let tol = precision::CONFUSION;
        let is_iso_u = dir_is_normal(zc, zs, tol) && o.distance(c.center) <= tol;
        let is_iso_v = dir_is_normal(xc, zs, tol) && dir_is_normal(yc, zs, tol);

        let mut p2d1 = DVec2::ZERO;
        let mut p2d2 = DVec2::ZERO;
        let mut d2d = DVec2::X;

        if is_iso_u {
            // OCCT L126: myType = GeomAbs_Line.
            self.projector.proj_type = CurveType::Line;

            // OCCT L128-129.
            p2d1 = Self::eval_pnt2d(xc, &self.my_sphere);
            p2d2 = Self::eval_pnt2d(yc, &self.my_sphere);

            // OCCT L131-138: P1 on the apex of the sphere takes its U from
            // P2d2 (U is undefined at the pole).
            if (p2d1.y - std::f64::consts::FRAC_PI_2).abs() < precision::PCONFUSION
                || (p2d1.y + std::f64::consts::FRAC_PI_2).abs() < precision::PCONFUSION
            {
                p2d1.x = p2d2.x;
            } else if ((p2d1.x - p2d2.x).abs() - std::f64::consts::PI).abs()
                < precision::PCONFUSION
            {
                // OCCT L139-157: U2 = U1 + PI; assume U1 = U2 with the
                // mirrored V2.
                p2d2.x = p2d1.x;
                if p2d2.y < 0.0 {
                    p2d2.y = -std::f64::consts::PI - p2d2.y;
                } else {
                    p2d2.y = std::f64::consts::PI - p2d2.y;
                }
            } else {
                // OCCT L158-160.
                p2d2.x = p2d1.x;
            }

            // OCCT L162: D2d = gp_Dir2d(gp_Vec2d(P2d1, P2d2)).
            let v2d = p2d2 - p2d1;
            d2d = v2d.normalize_or_zero();
            self.projector.is_done = true;
        } else if is_iso_v {
            // OCCT L165-176: the V-iso (parallel circle).
            self.projector.proj_type = CurveType::Line;

            let mut u = dir_angle_with_ref(xs, xc, xs.cross(ys));
            if u < 0.0 {
                u += TWO_PI;
            }
            let z = (c.center - o).dot(zs);
            let v = (z / self.my_sphere.radius).asin();
            p2d1 = DVec2::new(u, v);
            d2d = DVec2::new(xc.cross(yc).dot(xs.cross(ys)), 0.0);
            self.projector.is_done = true;
        }

        // OCCT L178: myLin = gp_Lin2d(P2d1, D2d) — stored in every case (the
        // OCCT member keeps the default line when the projection is not done).
        self.projector.lin = Line3::new(
            DVec3::new(p2d1.x, p2d1.y, 0.0),
            DVec3::new(d2d.x, d2d.y, 0.0),
        );
    }

    // OCCT L181-198: Project(gp_Lin/gp_Elips/gp_Parab/gp_Hypr) call the
    // ProjLib_Projector base bodies — the trait defaults.
}
// =========================================================================
// OCCT ProjLib_Torus (ProjLib_Torus.cxx L30-193)
// =========================================================================

/// OCCT ProjLib_Torus — the analytic projector onto a torus.
pub struct ProjLibTorus {
    /// OCCT ProjLib_Projector base.
    pub projector: Projector,
    /// OCCT: gp_Torus myTorus.
    pub my_torus: ToroidalSurface,
}

impl ProjLibTorus {
    /// OCCT ProjLib_Torus(const gp_Torus& To) -> Init(To) (L34-38, L49-56).
    pub fn new(to: &ToroidalSurface) -> Self {
        ProjLibTorus {
            projector: Projector::new(),
            my_torus: *to,
        }
    }

    /// OCCT Init(To) (L49-56).
    pub fn init(&mut self, to: &ToroidalSurface) {
        self.projector = Projector::new();
        self.my_torus = *to;
    }

    /// OCCT ProjLib_Torus.cxx L57-76 — static EvalPnt2d(Ve, To).
    fn eval_pnt2d(ve: DVec3, to: &ToroidalSurface) -> DVec2 {
        let x = ve.dot(to.ref_dir);
        let y = ve.dot(to.axis.cross(to.ref_dir).normalize());
        let u = if x.abs() > precision::PCONFUSION || y.abs() > precision::PCONFUSION {
            y.atan2(x)
        } else {
            0.0
        };
        DVec2::new(u, 0.0)
    }
}

/// The OCCT ProjLib_Torus::Project overloads.
impl ProjLibProjectorOverloads for ProjLibTorus {
    fn projector(&mut self) -> &mut Projector {
        &mut self.projector
    }

    /// OCCT ProjLib_Torus::Project(const gp_Circ& C) (L81-173).
    fn project_circ(&mut self, c: &Circle3) {
        // OCCT L83-89: the frames.
        let xc = c.x_dir;
        let yc = c.y_dir;
        let xt = self.my_torus.ref_dir;
        let yt = self.my_torus.axis.cross(self.my_torus.ref_dir).normalize();
        let zt = self.my_torus.axis;
        let oc = c.center - self.my_torus.center;

        let mut p1;
        let mut d2;

        if oc.length() < precision::CONFUSION
            || dir_is_parallel(c.normal, self.my_torus.axis, precision::ANGULAR)
        {
            // OCCT L91-141: Iso V.
            p1 = Self::eval_pnt2d(xc, &self.my_torus);
            let mut p2 = Self::eval_pnt2d(yc, &self.my_torus);
            let z = oc.dot(self.my_torus.axis) / self.my_torus.minor_radius;

            let v = if z > 1.0 {
                std::f64::consts::FRAC_PI_2 // simple protection
            } else if z < -1.0 {
                -std::f64::consts::FRAC_PI_2
            } else {
                z.asin()
            };

            let v = if c.radius < self.my_torus.major_radius {
                std::f64::consts::PI - v
            } else if v < 0.0 {
                v + TWO_PI
            } else {
                v
            };
            p1.y = v;
            p2.y = v;
            // OCCT L130-137: normally |P1.X() - P2.X()| = PI/2; a crossed
            // period reverses the direction.
            let v2d = p2 - p1;
            let v2d = if (p1.x - p2.x).abs() > std::f64::consts::PI {
                -v2d
            } else {
                v2d
            };
            d2 = v2d.normalize_or_zero();
            if p1.x < 0.0 {
                p1.x = TWO_PI + p1.x;
            }
        } else {
            // OCCT L142-166: Iso U.
            let mut u = dir_angle_with_ref(xt, oc, xt.cross(yt));
            if u < 0.0 {
                u += TWO_PI;
            }
            let mut v1 = dir_angle_with_ref(oc, xc, oc.cross(zt));
            if v1 < 0.0 {
                v1 += TWO_PI;
            }
            p1 = DVec2::new(u, v1);
            // OCCT L160-165: D2 = DY2d reversed when the circle sense is
            // negative w.r.t. (OC ^ Zt).
            d2 = DVec2::Y;
            if oc.cross(zt).dot(xc.cross(yc)) < 0.0 {
                d2 = -d2;
            }
        }

        self.projector.lin = Line3::new(
            DVec3::new(p1.x, p1.y, 0.0),
            DVec3::new(d2.x, d2.y, 0.0),
        );
        self.projector.proj_type = CurveType::Line;
        // OCCT L171: isDone = true.
        self.projector.is_done = true;
    }

    // OCCT L175-193: Project(gp_Lin/gp_Elips/gp_Parab/gp_Hypr) call the
    // ProjLib_Projector base bodies — the trait defaults.
}
// =========================================================================
// OCCT ProjLib_Projector::UFrame / VFrame (ProjLib_Projector.cxx L166-251)
// =========================================================================

impl Projector {
    /// OCCT ProjLib_Projector::UFrame(CFirst, CLast, UFirst, Period)
    /// (L166-192) — place the projected LINE inside the U period starting
    /// at UFirst.
    pub fn u_frame(&mut self, c_first: f64, _c_last: f64, u_first: f64, period: f64) {
        if self.proj_type == CurveType::Line {
            // OCCT L180-181: PFirst = ElCLib::Value(CFirst, myLin).
            let p_first = elclib_value_lin2d(&self.lin, c_first);
            let u = p_first.x;
            let new_u = crate::math::el::in_period(u, u_first, u_first + period);
            // OCCT L183: myLin.Translate(gp_Vec2d(NewU - U, 0.)).
            self.lin = proj_line_translate(self.lin, DVec2::new(new_u - u, 0.0));
        }
    }

    /// OCCT ProjLib_Projector::VFrame(CFirst, CLast, VFirst, Period)
    /// (L200-228).
    pub fn v_frame(&mut self, c_first: f64, _c_last: f64, v_first: f64, period: f64) {
        if self.proj_type == CurveType::Line {
            let p_first = elclib_value_lin2d(&self.lin, c_first);
            let v = p_first.y;
            let new_v = crate::math::el::in_period(v, v_first, v_first + period);
            self.lin = proj_line_translate(self.lin, DVec2::new(0.0, new_v - v));
        }
    }

    /// OCCT ProjLib_Projector::SetBezier(C) (L234-238) — the rcad encoding
    /// stores the bezier payload in the shared curve slot.
    pub fn set_bezier(&mut self, c: crate::geom::BezierCurve2) {
        self.bspline = Some(crate::geom::Curve2d::Bezier(c));
    }

    /// OCCT ProjLib_Projector::Bezier() (L244-250) — None models the null
    /// handle.
    pub fn bezier(&self) -> Option<crate::geom::BezierCurve2> {
        match self.bspline.as_ref() {
            Some(crate::geom::Curve2d::Bezier(b)) => Some(b.clone()),
            _ => None,
        }
    }

    /// OCCT ProjLib_Projector::BSpline() (L217-224) — the bspline payload.
    pub fn bspline_curve(&self) -> Option<crate::geom::BSplineCurve2> {
        match self.bspline.as_ref() {
            Some(crate::geom::Curve2d::BSpline(b)) => Some(b.clone()),
            _ => None,
        }
    }
}

// =========================================================================
// OCCT ProjLib_ProjectedCurve.cxx static helpers (L52-159)
// =========================================================================

/// OCCT GeomAbs_IsoType.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsoType {
    /// OCCT GeomAbs_IsoU.
    IsoU,
    /// OCCT GeomAbs_IsoV.
    IsoV,
}

/// OCCT ProjLib_ProjectedCurve.cxx L54-63 — static ComputeTolU.
pub fn compute_tol_u(surf: &dyn Adaptor3dSurface, tolerance: f64) -> f64 {
    let mut a_tol_u = surf.u_resolution(tolerance);
    if surf.is_u_periodic() {
        a_tol_u = a_tol_u.min(0.01 * surf.u_period());
    }
    a_tol_u
}

/// OCCT ProjLib_ProjectedCurve.cxx L67-76 — static ComputeTolV.
pub fn compute_tol_v(surf: &dyn Adaptor3dSurface, tolerance: f64) -> f64 {
    let mut a_tol_v = surf.v_resolution(tolerance);
    if surf.is_v_periodic() {
        a_tol_v = a_tol_v.min(0.01 * surf.v_period());
    }
    a_tol_v
}

/// OCCT ProjLib_ProjectedCurve.cxx L80-126 — static IsoIsDeg: whether the
/// iso of type `it` at `param` keeps a derivative magnitude inside
/// [TolMin, TolMax] along the sampled opposite parameter.
pub(crate) fn iso_is_deg(
    s: &dyn Adaptor3dSurface,
    param: f64,
    it: IsoType,
    tol_min: f64,
    tol_max: f64,
) -> bool {
    let u1 = s.first_u_parameter();
    let u2 = s.last_u_parameter();
    let v1 = s.first_v_parameter();
    let v2 = s.last_v_parameter();
    let mut along = true;
    let mut d1_norm_max = 0.0_f64;
    match it {
        IsoType::IsoV => {
            let step = (u2 - u1) / 10.0;
            let mut t = u1;
            while t <= u2 {
                let (_p, d1u, _d1v) = s.d1(t, param);
                d1_norm_max = d1_norm_max.max(d1u.length());
                t += step;
            }
            if d1_norm_max > tol_max || d1_norm_max < tol_min {
                along = false;
            }
        }
        IsoType::IsoU => {
            let step = (v2 - v1) / 10.0;
            let mut t = v1;
            while t <= v2 {
                let (_p, _d1u, d1v) = s.d1(param, t);
                d1_norm_max = d1_norm_max.max(d1v.length());
                t += step;
            }
            if d1_norm_max > tol_max || d1_norm_max < tol_min {
                along = false;
            }
        }
    }
    along
}

/// OCCT ProjLib_ProjectedCurve.cxx L130-159 — static TrimC3d: move the curve
/// ends off the degeneracy pole.  `is_trimmed` / `singular_case` are the
/// two-slot OCCT out arrays.
pub(crate) fn trim_c3d(
    my_curve: &mut std::sync::Arc<dyn Adaptor3dCurveGeom>,
    is_trimmed: &mut [bool; 2],
    dt: f64,
    pole: DVec3,
    singular_case: &mut [i32; 2],
    number_of_singular_case: i32,
    tol_conf: f64,
) {
    let mut f = my_curve.first_parameter();
    let mut l = my_curve.last_parameter();

    let mut p = my_curve.value(f);
    if p.distance(pole) <= tol_conf {
        is_trimmed[0] = true;
        f += dt;
        *my_curve = my_curve.trim_geom(f, l, precision::CONFUSION);
        singular_case[0] = number_of_singular_case;
    }

    p = my_curve.value(l);
    if p.distance(pole) <= tol_conf {
        is_trimmed[1] = true;
        l -= dt;
        *my_curve = my_curve.trim_geom(f, l, precision::CONFUSION);
        singular_case[1] = number_of_singular_case;
    }
}

// =========================================================================
// OCCT ProjLib_ProjectedCurve.cxx L163-238 — static ExtendC2d
// =========================================================================

/// OCCT static ExtendC2d(aRes, t, dt, u1, u2, v1, v2, FirstOrLast,
/// NumberOfSingularCase) (ProjLib_ProjectedCurve.cxx L163-238): appends the
/// straight segment running from the curve end to the degenerate-surface
/// boundary onto the 2d result (consumed by the BSpline-surface branch
/// L511-522, the default branch L678-687 and the ComputeApprox fallback
/// L742-753 of Perform).
///
/// GAP leaf: the body needs `Geom2dConvert_CompCurveToBSplineCurve`
/// (Geom2dConvert_CompCurveToBSplineCurve.cxx L20-250, the C1
/// concat-onto-BSpline machinery with `Add(aSegment, aTol, anAfter)` +
/// `BSplineCurve()`), not translated in rcad-kernel
/// (base::geom2d_convert carries only the simplified
/// `compose_curves_to_bspline` helper); the anchor is preserved with the
/// OCCT failure path (aRes unchanged, caller keeps the unextended result —
/// the same state OCCT produces when the concat is refused).  Remaining
/// consumed pieces once the dependency lands: Curve2d D1 at the end
/// parameter, the boundary direction switch (L188-210), the parallel /
/// intersection parameter (L213-229) and the trimmed-segment concat
/// (L231-237).
#[allow(clippy::too_many_arguments)]
pub(crate) fn extend_c2d(
    _a_res: &mut Curve2d,
    _t: f64,
    _dt: f64,
    _u1: f64,
    _u2: f64,
    _v1: f64,
    _v2: f64,
    _first_or_last: i32,
    _number_of_singular_case: i32,
) {
    unimplemented!(
        "ProjLib_ProjectedCurve::ExtendC2d (OCCT L163-238) needs \
         Geom2dConvert_CompCurveToBSplineCurve (untranslated)"
    );
}

// =========================================================================
// OCCT ProjLib_ProjectedCurve.cxx L242-270 — static Project dispatch
// =========================================================================

/// OCCT static Project(P, C): dispatch the curve kind onto the projector
/// member.  The BSpline / Bezier / Offset / Other kinds fall through (the
/// approximation route); OCCT throws Standard_NoSuchObject on the `default`
/// arm which the rcad closed [`CurveType`] enum cannot express.
pub(crate) fn project_dispatch<P: ProjLibProjectorOverloads + ?Sized>(
    p: &mut P,
    c: &dyn Adaptor3dCurveGeom,
) {
    let c_type = c.get_type();
    match c_type {
        CurveType::Line => {
            let l = c.line();
            p.project_lin(&l);
        }
        CurveType::Circle => {
            let ci = c.circle();
            p.project_circ(&ci);
        }
        CurveType::Ellipse => {
            let e = c.ellipse();
            p.project_elips(&e);
        }
        CurveType::Hyperbola => {
            let h = c.hyperbola();
            p.project_hypr(&h);
        }
        CurveType::Parabola => {
            let pa = c.parabola();
            p.project_parab(&pa);
        }
        // OCCT L262-265: BSplineCurve / BezierCurve / OffsetCurve /
        // OtherCurve — try the approximation (nothing to do here).
        CurveType::BSpline | CurveType::Bezier | CurveType::Other => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_project_meridian_circle_is_line() {
        // A great circle through the poles (a meridian, iso-U) projects to
        // the iso-U line U = const.
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        };
        let mut circle = Circle3::new(DVec3::ZERO, DVec3::Y, 2.0);
        circle.x_dir = DVec3::Z;
        circle.y_dir = DVec3::X;
        let mut proj = ProjLibSphere::new(&sphere);
        proj.project_circ(&circle);
        assert!(proj.projector.is_done);
        assert_eq!(proj.projector.get_type(), CurveType::Line);
    }

    #[test]
    fn sphere_project_equator_is_iso_v_line() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
        };
        let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        let mut proj = ProjLibSphere::new(&sphere);
        proj.project_circ(&circle);
        assert!(proj.projector.is_done);
        assert_eq!(proj.projector.get_type(), CurveType::Line);
        // The equator projects to V = 0.
        assert!(proj_line_location(&proj.projector.lin).y.abs() < 1e-12);
    }

    #[test]
    fn sphere_non_iso_circle_not_done() {
        // A small circle off-axis is neither iso-U nor iso-V: the OCCT
        // member keeps type OtherCurve and isDone == false (the projection
        // then goes through the approximation route).
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
        };
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 0.5), DVec3::new(1.0, 0.0, 1.0), 0.5);
        let mut proj = ProjLibSphere::new(&sphere);
        proj.project_circ(&circle);
        assert!(!proj.projector.is_done);
    }

    #[test]
    fn cylinder_project_axis_parallel_line() {
        let cyl = CylindricalSurface::new_with_ref_dir(DVec3::ZERO, DVec3::Z, 5.0, DVec3::X);
        let line = Line3::new(DVec3::new(5.0, 0.0, 1.0), DVec3::Z);
        let mut proj = ProjLibCylinder::new(&cyl);
        proj.project_lin(&line);
        assert!(proj.projector.is_done);
        assert_eq!(proj.projector.get_type(), CurveType::Line);
    }

    #[test]
    fn plane_project_line() {
        let pl = Plane::new(DVec3::ZERO, DVec3::Z);
        let line = Line3::new(DVec3::new(1.0, 2.0, 0.0), DVec3::X);
        let mut proj = ProjLibPlane::new(&pl);
        proj.project_lin(&line);
        assert!(proj.projector.is_done);
        assert_eq!(proj.projector.get_type(), CurveType::Line);
        let loc = proj_line_location(&proj.projector.lin);
        assert!((loc.x - 1.0).abs() < 1e-12 && (loc.y - 2.0).abs() < 1e-12);
    }
}
