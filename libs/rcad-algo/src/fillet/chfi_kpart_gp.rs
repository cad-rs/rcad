//! OCCT gp / ElCLib / ElSLib / gce analytic primitives used by the ChFiKPart
//! particular-case translations (fillet/chamfer SurfData computation).
//!
//! Each item carries its OCCT anchor (file + lines) from:
//!   - FoundationClasses/TKMath/gp (gp_Ax2 / gp_Ax3 / gp_Dir / gp_Mat / gp_Vec)
//!   - FoundationClasses/TKMath/ElCLib/ElCLib.cxx
//!   - FoundationClasses/TKMath/ElSLib/ElSLib.cxx
//!   - ModelingData/TKGeomBase/gce/gce_MakeCirc.cxx
//!
//! These are the support primitives the ChFiKPart case code manipulates as
//! gp_Ax3 / gp_Circ / Geom_* values; the rcad Surface3 / Curve2d / Curve3
//! payloads are produced from them at the registration points (D6-style
//! conversion at `IndexSurfaceInDS` / `IndexCurveInDS` call sites).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Circle2d, Circle3, ConicalSurface, CylindricalSurface, Line2d, Plane, SphericalSurface,
    Surface3, ToroidalSurface,
};
use rcad_kernel::topo::topods::Orientation;

/// OCCT Precision::Confusion() (rcad-kernel core::precision::CONFUSION).
pub const CONFUSION: f64 = 1e-7;
/// OCCT Precision::Angular().
pub const ANGULAR: f64 = 1e-12;
/// OCCT Precision::PConfusion().
pub const P_CONFUSION: f64 = 1e-9;
/// OCCT gp::Resolution() (gp.hxx) = 1e-12.
pub const GP_RESOLUTION: f64 = 1e-12;
/// OCCT RealEpsilon() = DBL_EPSILON.
pub const REAL_EPSILON: f64 = f64::EPSILON;

// =========================================================================
// OCCT gp_Ax3 (gp_Ax3.hxx L40-131).  A full 3D coordinate system which can
// be left-handed after YReverse / ZReverse (Direct() = (X ^ Y) . Z > 0).
// The gp_Ax2 constructors are shared: gp_Ax2 is a right-handed Ax3 without
// handedness state (gp_Ax3.hxx L477-483 models Ax3(Ax2) exactly this way).
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct GpAx3 {
    pub location: DVec3,
    pub vxdir: DVec3,
    pub vydir: DVec3,
    pub vzdir: DVec3,
}

impl GpAx3 {
    /// OCCT gp_Ax2/gp_Ax3 (P, N, Vx) constructor (gp_Ax2.hxx L72-77):
    /// vxdir = (N ^ Vx) ^ N; vydir = N ^ vxdir.
    pub fn new_pn_vx(p: DVec3, n: DVec3, vx: DVec3) -> Self {
        let n = n.normalize();
        let vxdir = cross_cross(n, vx, n).normalize();
        let vydir = n.cross(vxdir).normalize();
        GpAx3 {
            location: p,
            vxdir,
            vydir,
            vzdir: n,
        }
    }

    /// OCCT gp_Ax2/gp_Ax3 (P, V) constructor (gp_Ax2.cxx L31-80): the X
    /// direction is the unit vector perpendicular to V with the smallest
    /// component rule — identical to rcad Plane::new (gp_Ax3.cxx L29-80).
    pub fn new_pn(p: DVec3, n: DVec3) -> Self {
        let pl = Plane::new(p, n);
        GpAx3 {
            location: p,
            vxdir: pl.u_dir,
            vydir: pl.v_dir,
            vzdir: pl.normal,
        }
    }

    /// OCCT gp_Ax3.hxx L128 — YReverse: vydir.Reverse().
    pub fn y_reverse(&mut self) {
        self.vydir = -self.vydir;
    }

    /// OCCT gp_Ax3.hxx L131 — ZReverse: axis direction reversed.
    pub fn z_reverse(&mut self) {
        self.vzdir = -self.vzdir;
    }

    /// OCCT gp_Ax3.hxx L208 — Direct(): (X ^ Y) . Z > 0.
    pub fn direct(&self) -> bool {
        self.vxdir.cross(self.vydir).dot(self.vzdir) > 0.0
    }

    /// OCCT gp_Ax3.hxx L486-493 — Ax2(): Z2 = Direct ? Z : -Z, then the
    /// (P, N, Vx) Ax2 construction recomputes X/Y from (Z2, X).
    pub fn ax2(&self) -> GpAx3 {
        let mut zz = self.vzdir;
        if !self.direct() {
            zz = -zz;
        }
        GpAx3::new_pn_vx(self.location, zz, self.vxdir)
    }

    /// OCCT gp_Ax3.hxx — SetLocation.
    pub fn set_location(&mut self, p: DVec3) {
        self.location = p;
    }

    /// OCCT gp_Ax3.hxx L84 — Ax3(P, N, Vx) ("theVx gives the XDirection").
    pub fn new(p: DVec3, n: DVec3, vx: DVec3) -> Self {
        Self::new_pn_vx(p, n, vx)
    }

    // gp_Ax3 accessors.
    pub fn x_direction(&self) -> DVec3 {
        self.vxdir
    }

    pub fn y_direction(&self) -> DVec3 {
        self.vydir
    }

    pub fn direction(&self) -> DVec3 {
        self.vzdir
    }
}

/// OCCT gp_XYZ::CrossCrossed(A2, A3) — (this ^ A2) ^ A3.
pub fn cross_cross(v: DVec3, a2: DVec3, a3: DVec3) -> DVec3 {
    v.cross(a2).cross(a3)
}

/// OCCT gp_Dir::Angle (gp_Dir.cxx L20-42) — the angle in [0, PI].
pub fn dir_angle(d1: DVec3, d2: DVec3) -> f64 {
    let cosinus = d1.dot(d2);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        return cosinus.acos();
    }
    let sinus = d1.cross(d2).length();
    if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    }
}

/// OCCT gp_Dir::AngleWithRef (gp_Dir.cxx L55-84) — the signed angle.
pub fn dir_angle_with_ref(d1: DVec3, d2: DVec3, vref: DVec3) -> f64 {
    let xyz = d1.cross(d2);
    let cosinus = d1.dot(d2);
    let sinus = xyz.length();
    let ang;
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        ang = cosinus.acos();
    } else if cosinus < 0.0 {
        ang = std::f64::consts::PI - sinus.asin();
    } else {
        ang = sinus.asin();
    }
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT gp_Mat::SetRotation (gp_Mat.cxx L122-159) — the Rodrigues rotation
/// matrix; applied to a vector this is gp_Vec::Rotate (gp_Vec.hxx L521-526,
/// gp_Trsf::SetRotation + VectorialPart).  The translation part of the Ax1
/// does not act on (free) vectors.
pub fn rotate_vec_around_ax1(v: DVec3, axis_dir: DVec3, ang: f64) -> DVec3 {
    let av = axis_dir.normalize();
    let a = av.x;
    let b = av.y;
    let c = av.z;
    let a_cos = ang.cos();
    let a_sin = ang.sin();
    let a_om_cos = 1.0 - a_cos;
    let a2 = a * a;
    let b2 = b * b;
    let c2 = c * c;
    let ab = a * b;
    let ac = a * c;
    let bc = b * c;
    let m00 = 1.0 + a_om_cos * (-(b2 + c2));
    let m01 = a_om_cos * ab - a_sin * c;
    let m02 = a_om_cos * ac + a_sin * b;
    let m10 = a_om_cos * ab + a_sin * c;
    let m11 = 1.0 + a_om_cos * (-(a2 + c2));
    let m12 = a_om_cos * bc - a_sin * a;
    let m20 = a_om_cos * ac - a_sin * b;
    let m21 = a_om_cos * bc + a_sin * a;
    let m22 = 1.0 + a_om_cos * (-(a2 + b2));
    DVec3::new(
        v.x * m00 + v.y * m01 + v.z * m02,
        v.x * m10 + v.y * m11 + v.z * m12,
        v.x * m20 + v.y * m21 + v.z * m22,
    )
}

// =========================================================================
// OCCT ElCLib (ElCLib.cxx) — line / circle evaluation.
// =========================================================================

/// OCCT ElCLib::LineD1 (ElCLib.cxx L229-234).
pub fn elclib_line_d1(u: f64, pos: DVec3, dir: DVec3) -> (DVec3, DVec3) {
    let p = dir * u + pos;
    (p, dir)
}

/// OCCT ElCLib::CircleD1 (ElCLib.cxx L239-248).
pub fn elclib_circle_d1(u: f64, pos: &GpAx3, radius: f64) -> (DVec3, DVec3) {
    let xc = radius * u.cos();
    let yc = radius * u.sin();
    let coord1 = pos.vxdir;
    let coord2 = pos.vydir;
    let p = coord1 * xc + coord2 * yc + pos.location;
    let v1 = coord1 * (-yc) + coord2 * xc;
    (p, v1)
}

/// OCCT ElCLib::LineParameter (ElCLib.cxx L1191-1194).
pub fn elclib_line_parameter(pos: DVec3, dir: DVec3, p: DVec3) -> f64 {
    (p - pos).dot(dir)
}

/// OCCT ElCLib::CircleParameter (ElCLib.cxx L1199-1218).
pub fn elclib_circle_parameter(pos: &GpAx3, p: DVec3) -> f64 {
    let avec = p - pos.location;
    if avec.length_squared() < GP_RESOLUTION {
        // coinciding points -> infinite number of parameters
        return 0.0;
    }
    let dir = pos.vzdir;
    // Project vector on circle's plane
    let a_v_proj = cross_cross(dir, avec, dir);
    if a_v_proj.length_squared() < GP_RESOLUTION {
        return 0.0;
    }
    // Angle between X direction and projected vector
    let mut teta = dir_angle_with_ref(pos.vxdir, a_v_proj, dir);
    normalize_angle_elclib(&mut teta);
    teta
}

/// OCCT ElCLib.cxx normalizeAngle (ElCLib.cxx L28-55, shared namespace
/// helper) — maps to [0, 2PI].
pub fn normalize_angle_elclib(angle: &mut f64) {
    let pipi = std::f64::consts::PI + std::f64::consts::PI;
    let negative_resolution = -P_CONFUSION;
    while *angle < negative_resolution {
        *angle += pipi;
    }
    while *angle > pipi * (1.0 + GP_RESOLUTION) {
        *angle -= pipi;
    }
    if *angle < 0.0 {
        *angle = 0.0;
    }
}

// =========================================================================
// OCCT ElSLib (ElSLib.cxx) — analytic surface evaluation on a gp_Ax3 frame.
// Value/D0/D1 bodies: ElSLib.cxx L61-155 (Value), L562-863 (D0/D1).
// Parameters bodies: ElSLib.cxx L1547-1700.
// =========================================================================

/// OCCT ElSLib::PlaneParameters (ElSLib.cxx L1547-1554): local coordinates.
pub fn elslib_plane_parameters(pos: &GpAx3, p: DVec3) -> (f64, f64) {
    let d = p - pos.location;
    (d.dot(pos.vxdir), d.dot(pos.vydir))
}

/// OCCT ElSLib::CylinderParameters (ElSLib.cxx L1558-1570).
pub fn elslib_cylinder_parameters(pos: &GpAx3, _radius: f64, p: DVec3) -> (f64, f64) {
    let d = p - pos.location;
    let mut u = d.dot(pos.vydir).atan2(d.dot(pos.vxdir));
    normalize_angle_elclib(&mut u);
    let v = d.dot(pos.vzdir);
    (u, v)
}

/// OCCT ElSLib::ConeParameters (ElSLib.cxx L1574-1611).
pub fn elslib_cone_parameters(pos: &GpAx3, radius: f64, sangle: f64, p: DVec3) -> (f64, f64) {
    let d = p - pos.location;
    let x = d.dot(pos.vxdir);
    let y = d.dot(pos.vydir);
    let z = d.dot(pos.vzdir);

    // Check if point is at the apex
    let u;
    if x.abs() < GP_RESOLUTION && y.abs() < GP_RESOLUTION {
        u = 0.0;
    } else if -radius > z * sangle.tan() {
        // the point is at the wrong side of the apex
        u = (-y).atan2(-x);
    } else {
        u = y.atan2(x);
    }
    let mut u = u;
    normalize_angle_elclib(&mut u);

    // V = sin(Sang) * (x cosU + y SinU - R) + z * cos(Sang)
    let v = sangle.sin() * (x * u.cos() + y * u.sin() - radius) + sangle.cos() * z;
    (u, v)
}

/// OCCT ElSLib::SphereParameters (ElSLib.cxx L1615-1645).
pub fn elslib_sphere_parameters(pos: &GpAx3, _radius: f64, p: DVec3) -> (f64, f64) {
    let d = p - pos.location;
    let x = d.dot(pos.vxdir);
    let y = d.dot(pos.vydir);
    let z = d.dot(pos.vzdir);
    let l = (x * x + y * y).sqrt();
    if l < GP_RESOLUTION {
        // point on axis Z of the sphere
        if z > 0.0 {
            (0.0, std::f64::consts::FRAC_PI_2)
        } else {
            (0.0, -std::f64::consts::FRAC_PI_2)
        }
    } else {
        let v = (z / l).atan();
        let mut u = y.atan2(x);
        normalize_angle_elclib(&mut u);
        (u, v)
    }
}

/// OCCT ElSLib::TorusParameters (ElSLib.cxx L1649-1698).
pub fn elslib_torus_parameters(
    pos: &GpAx3,
    major_radius: f64,
    minor_radius: f64,
    p: DVec3,
) -> (f64, f64) {
    let d = p - pos.location;
    let x = d.dot(pos.vxdir);
    let y = d.dot(pos.vydir);
    let z = d.dot(pos.vzdir);

    // all that to process case of Major < Minor.
    let mut u = y.atan2(x);
    if major_radius < minor_radius {
        let cosu = u.cos();
        let sinu = u.sin();
        let z2 = z * z;
        let min_r2 = minor_radius * minor_radius;
        let rcosu = major_radius * cosu;
        let rsinu = major_radius * sinu;
        let xm = x - rcosu;
        let ym = y - rsinu;
        let xp = x + rcosu;
        let yp = y + rsinu;
        let d1 = xm * xm + ym * ym + z2 - min_r2;
        let d2 = xp * xp + yp * yp + z2 - min_r2;
        let ad1 = d1.abs();
        let ad2 = d2.abs();
        if ad2 < ad1 {
            u += std::f64::consts::PI;
        }
    }
    normalize_angle_elclib(&mut u);
    let cosu = u.cos();
    let sinu = u.sin();
    let dx = DVec3::new(cosu, sinu, 0.0);
    let dpv = DVec3::new(x - major_radius * cosu, y - major_radius * sinu, z);
    let a_mag = dpv.length();
    let mut v;
    if a_mag <= GP_RESOLUTION {
        v = 0.0;
    } else {
        let dp = dpv / a_mag;
        // OCCT: dx.AngleWithRef(dP, dx ^ gp::DZ()); the gp_Dir arguments are
        // the 2D vectors embedded in XY.
        v = dir_angle_with_ref(dx, dp, dx.cross(DVec3::Z));
    }
    normalize_angle_elclib(&mut v);
    (u, v)
}

/// OCCT ElSLib::PlaneD0 (ElSLib.cxx L562-570).
pub fn elslib_plane_d0(u: f64, v: f64, pos: &GpAx3) -> DVec3 {
    pos.vxdir * u + pos.vydir * v + pos.location
}

/// OCCT ElSLib::PlaneD1 (ElSLib.cxx L666-683).
pub fn elslib_plane_d1(u: f64, v: f64, pos: &GpAx3) -> (DVec3, DVec3, DVec3) {
    let p = pos.vxdir * u + pos.vydir * v + pos.location;
    (p, pos.vxdir, pos.vydir)
}

/// OCCT ElSLib::CylinderD1 (ElSLib.cxx L732-753).
pub fn elslib_cylinder_d1(u: f64, v: f64, pos: &GpAx3, radius: f64) -> (DVec3, DVec3, DVec3) {
    let a1 = radius * u.cos();
    let a2 = radius * u.sin();
    let p = pos.vxdir * a1 + pos.vydir * a2 + pos.vzdir * v + pos.location;
    let vu = pos.vxdir * (-a2) + pos.vydir * a1;
    let vv = pos.vzdir;
    (p, vu, vv)
}

/// OCCT ElSLib::ConeD1 (ElSLib.cxx L687-728).
pub fn elslib_cone_d1(
    u: f64,
    v: f64,
    pos: &GpAx3,
    radius: f64,
    sangle: f64,
) -> (DVec3, DVec3, DVec3) {
    let cos_u = u.cos();
    let sin_u = u.sin();
    let cos_a = sangle.cos();
    let sin_a = sangle.sin();
    let r = radius + v * sin_a;
    let a3 = v * cos_a;
    let a1 = r * cos_u;
    let a2 = r * sin_u;
    let r1 = sin_a * cos_u;
    let r2 = sin_a * sin_u;
    let p = pos.vxdir * a1 + pos.vydir * a2 + pos.vzdir * a3 + pos.location;
    let vu = pos.vxdir * (-a2) + pos.vydir * a1;
    let vv = pos.vxdir * r1 + pos.vydir * r2 + pos.vzdir * cos_a;
    (p, vu, vv)
}

/// OCCT ElSLib::SphereD1 (ElSLib.cxx L757-793).
pub fn elslib_sphere_d1(u: f64, v: f64, pos: &GpAx3, radius: f64) -> (DVec3, DVec3, DVec3) {
    let cos_u = u.cos();
    let sin_u = u.sin();
    let r1 = radius * v.cos();
    let r2 = radius * v.sin();
    let a1 = r1 * cos_u;
    let a2 = r1 * sin_u;
    let a3 = r2 * cos_u;
    let a4 = r2 * sin_u;
    let p = pos.vxdir * a1 + pos.vydir * a2 + pos.vzdir * r2 + pos.location;
    let vu = pos.vxdir * (-a2) + pos.vydir * a1;
    let vv = pos.vxdir * (-a3) + pos.vydir * (-a4) + pos.vzdir * r1;
    (p, vu, vv)
}

/// OCCT ElSLib::TorusD0 / TorusD1 (ElSLib.cxx L628-662 / L797-863) — with the
/// OCC620 eps zeroing of the small components.
pub fn elslib_torus_d1(
    u: f64,
    v: f64,
    pos: &GpAx3,
    major_radius: f64,
    minor_radius: f64,
) -> (DVec3, DVec3, DVec3) {
    let cos_u = u.cos();
    let sin_u = u.sin();
    let r1 = minor_radius * v.cos();
    let r2 = minor_radius * v.sin();
    let r = major_radius + r1;
    let mut a1 = r * cos_u;
    let mut a2 = r * sin_u;
    let mut a3 = r2 * cos_u;
    let mut a4 = r2 * sin_u;
    // Modified by skv - Tue Sep 9 15:10:34 2003 OCC620 Begin
    let eps = 10.0 * (minor_radius + major_radius) * REAL_EPSILON;
    if a1.abs() <= eps {
        a1 = 0.0;
    }
    if a2.abs() <= eps {
        a2 = 0.0;
    }
    if a3.abs() <= eps {
        a3 = 0.0;
    }
    if a4.abs() <= eps {
        a4 = 0.0;
    }
    // Modified by skv - Tue Sep 9 15:10:35 2003 OCC620 End
    let p = pos.vxdir * a1 + pos.vydir * a2 + pos.vzdir * r2 + pos.location;
    let vu = pos.vxdir * (-a2) + pos.vydir * a1;
    let vv = pos.vxdir * (-a3) + pos.vydir * (-a4) + pos.vzdir * r1;
    (p, vu, vv)
}

/// OCCT ElSLib::TorusD0 (ElSLib.cxx L628-662).
pub fn elslib_torus_d0(
    u: f64,
    v: f64,
    pos: &GpAx3,
    major_radius: f64,
    minor_radius: f64,
) -> DVec3 {
    let (p, _vu, _vv) = elslib_torus_d1(u, v, pos, major_radius, minor_radius);
    p
}

/// OCCT ElSLib::SphereD0 (ElSLib.cxx L609-626).
pub fn elslib_sphere_d0(u: f64, v: f64, pos: &GpAx3, radius: f64) -> DVec3 {
    let r = radius * v.cos();
    let a3 = radius * v.sin();
    let a1 = r * u.cos();
    let a2 = r * u.sin();
    pos.vxdir * a1 + pos.vydir * a2 + pos.vzdir * a3 + pos.location
}

// =========================================================================
// OCCT gp_Circ — circle positioned by an (Ax2) frame.
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct GpCirc {
    pub pos: GpAx3,
    pub radius: f64,
}

impl GpCirc {
    /// OCCT gp_Circ(A2, Radius).
    pub fn new(pos: GpAx3, radius: f64) -> Self {
        GpCirc { pos, radius }
    }

    pub fn location(&self) -> DVec3 {
        self.pos.location
    }

    pub fn axis_direction(&self) -> DVec3 {
        self.pos.vzdir
    }

    pub fn x_direction(&self) -> DVec3 {
        self.pos.vxdir
    }

    pub fn y_direction(&self) -> DVec3 {
        self.pos.vydir
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// OCCT gp_Ax2::SetLocation (position change in place).
    pub fn set_location(&mut self, p: DVec3) {
        self.pos.location = p;
    }

    /// rcad conversion — OCCT gp_Circ -> Geom_Circle (Curve3::Circle).
    pub fn to_circle3(&self) -> Circle3 {
        Circle3 {
            center: self.pos.location,
            normal: self.pos.vzdir,
            x_dir: self.pos.vxdir,
            y_dir: self.pos.vydir,
            radius: self.radius,
        }
    }

    /// rcad conversion — OCCT Geom_Circle / gp_Circ from the rcad payload.
    pub fn from_circle3(c: &Circle3) -> Self {
        GpCirc {
            pos: GpAx3 {
                location: c.center,
                vxdir: c.x_dir,
                vydir: c.y_dir,
                vzdir: c.normal,
            },
            radius: c.radius,
        }
    }
}

// =========================================================================
// OCCT gce_MakeCirc (gce_MakeCirc.cxx L32-156) — circle through 3 points.
// The perpendicular-bisector intersection (OCCT Extrema_ExtElC between the
// bisector lines L1/L2) is resolved directly: the closest points of the two
// bisector lines coincide at the circumcenter, and OCCT takes their midpoint
// (L117-131) with radius = (dist1 + dist2 + dist3) / 3 (L152).  The error
// branches keep the OCCT order (L48-59 coincident/confused, L71-83 colinear,
// L94-114 intersection errors).
// =========================================================================
pub fn gce_make_circ(p1: DVec3, p2: DVec3, p3: DVec3) -> Result<GpCirc, &'static str> {
    let mut dist1 = p1.distance(p2);
    let mut dist2 = p1.distance(p3);
    let mut dist3 = p2.distance(p3);
    let a_resolution = GP_RESOLUTION;

    if (dist1 < a_resolution) && (dist2 < a_resolution) && (dist3 < a_resolution) {
        // OCCT L50-53: degenerate zero-radius circle at P1 with (N=X, Vx=Z).
        return Ok(GpCirc::new(GpAx3::new_pn_vx(p1, DVec3::X, DVec3::Z), 0.0));
    }
    if dist1 < a_resolution || dist2 < a_resolution {
        // OCCT L55-58: gce_ConfusedPoints.
        return Err("gce_ConfusedPoints");
    }

    let dir1 = (p2 - p1).normalize();
    let vdir2 = p3 - p2;

    // OCCT L69-75: gp_Lin(P1, Dir1).Distance(P3) < Resolution -> colinear.
    if (p3 - p1).cross(dir1).length() < a_resolution {
        return Err("gce_ColinearPoints");
    }

    let vdir3 = dir1.cross(vdir2);
    if vdir3.length_squared() < a_resolution {
        // OCCT L79-83: gce_ColinearPoints.
        return Err("gce_ColinearPoints");
    }

    // Circumcenter of the triangle (p1, p2, p3) — the intersection point of
    // the two perpendicular bisectors (OCCT L87-131, Extrema between L1/L2
    // with pInt = midpoint of the closest points).
    let v21 = p2 - p1;
    let v31 = p3 - p1;
    let nn = 2.0 * v21.cross(v31).length_squared();
    if nn.abs() < a_resolution {
        return Err("gce_IntersectionError");
    }
    let p_int = p1 + (v21.cross(v31).cross(v21) * v31.length_squared()
        + v31.cross(v21).cross(v31) * v21.length_squared())
        / nn;
    // OCCT L136-138: recompute the three distances from pInt.
    dist1 = p1.distance(p_int);
    dist2 = p2.distance(p_int);
    dist3 = p3.distance(p_int);
    if dist1 < a_resolution {
        // OCCT L140-146: degenerate zero-radius circle.
        return Ok(GpCirc::new(
            GpAx3::new_pn_vx(p_int, DVec3::X, DVec3::Z),
            0.0,
        ));
    }
    let dir1bis = (p1 - p_int).normalize();
    // OCCT L152: gp_Ax2(pInt, gp_Dir(VDir3), Dir1) with the averaged radius.
    let ax = GpAx3::new_pn_vx(p_int, vdir3.normalize(), dir1bis);
    Ok(GpCirc::new(ax, (dist1 + dist2 + dist3) / 3.0))
}

// =========================================================================
// OCCT Geom_ConicalSurface (TKG3d/Geom/Geom_ConicalSurface.cxx) — the
// working carrier mutated by the chamfer case code (VReverse / SetPosition
// / Position / SemiAngle / Axis).  Registered as a rcad Surface3 payload at
// the ChangeSurf point.
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct GpConicalSurface {
    pub pos: GpAx3,
    /// OCCT: double semiAngle.
    pub semi_angle: f64,
    /// OCCT: double radius (RefRadius at the Ax3 location).
    pub ref_radius: f64,
}

impl GpConicalSurface {
    /// OCCT Geom_ConicalSurface(A3, Ang, Radius).
    pub fn new(pos: GpAx3, semi_angle: f64, ref_radius: f64) -> Self {
        GpConicalSurface {
            pos,
            semi_angle,
            ref_radius,
        }
    }

    /// OCCT Geom_ConicalSurface.cxx L104-108 — VReverse:
    /// semiAngle = -semiAngle; pos.ZReverse().
    pub fn v_reverse(&mut self) {
        self.semi_angle = -self.semi_angle;
        self.pos.z_reverse();
    }

    /// OCCT Geom_ConicalSurface.hxx L112-115 — SemiAngle().
    pub fn semi_angle(&self) -> f64 {
        self.semi_angle
    }

    /// OCCT Geom_ElementarySurface — Position().
    pub fn position(&self) -> GpAx3 {
        self.pos
    }

    /// OCCT Geom_ElementarySurface — SetPosition(A3).
    pub fn set_position(&mut self, pos: GpAx3) {
        self.pos = pos;
    }

    /// OCCT Geom_ElementarySurface — Axis().Direction() (the possibly
    /// reversed axis after VReverse).
    pub fn axis_direction(&self) -> DVec3 {
        self.pos.vzdir
    }

    /// rcad conversion — the Surface3 payload for IndexSurfaceInDS.
    pub fn to_surface3(&self) -> Surface3 {
        Surface3::Cone(ConicalSurface {
            apex: self.pos.location,
            axis: self.pos.vzdir,
            radius: self.ref_radius,
            half_angle_rad: self.semi_angle,
            ref_dir: self.pos.vxdir,
        })
    }
}

// =========================================================================
// OCCT Geom_CylindricalSurface — the working carrier mutated by the case
// code.  VReverse is Geom_ElementarySurface::VReverse (cxx L51-54).
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct GpCylindricalSurface {
    pub pos: GpAx3,
    pub radius: f64,
}

impl GpCylindricalSurface {
    /// OCCT Geom_CylindricalSurface(A3, Radius).
    pub fn new(pos: GpAx3, radius: f64) -> Self {
        GpCylindricalSurface { pos, radius }
    }

    /// OCCT Geom_ElementarySurface.cxx L51-54 — VReverse: pos.ZReverse().
    pub fn v_reverse(&mut self) {
        self.pos.z_reverse();
    }

    /// OCCT Geom_ElementarySurface — Position().
    pub fn position(&self) -> GpAx3 {
        self.pos
    }

    /// OCCT Geom_ElementarySurface — SetPosition(A3).
    pub fn set_position(&mut self, pos: GpAx3) {
        self.pos = pos;
    }

    /// rcad conversion — the Surface3 payload for IndexSurfaceInDS.
    pub fn to_surface3(&self) -> Surface3 {
        Surface3::Cylinder(CylindricalSurface {
            origin: self.pos.location,
            axis: self.pos.vzdir,
            radius: self.radius,
            ref_dir: self.pos.vxdir,
            y_dir: Some(self.pos.vydir),
        })
    }
}

// =========================================================================
// rcad conversions for the remaining OCCT Geom elementary surfaces used by
// ChFiKPart (registered through ChFiKPart_IndexSurfaceInDS).
// =========================================================================

/// OCCT Geom_Plane(Ax3) -> rcad Surface3 payload.
pub fn surface3_plane(pos: &GpAx3) -> Surface3 {
    Surface3::Plane(Plane::with_axes(pos.location, pos.vzdir, pos.vxdir))
}

/// OCCT Geom_SphericalSurface(Ax3, Radius) -> rcad Surface3 payload.
pub fn surface3_sphere(pos: &GpAx3, radius: f64) -> Surface3 {
    Surface3::Sphere(SphericalSurface {
        center: pos.location,
        axis: pos.vzdir,
        radius,
        ref_dir: pos.vxdir,
    })
}

/// OCCT Geom_ToroidalSurface(Ax3, MajorRadius, MinorRadius) -> rcad payload.
pub fn surface3_torus(pos: &GpAx3, major_radius: f64, minor_radius: f64) -> Surface3 {
    Surface3::Torus(ToroidalSurface {
        center: pos.location,
        axis: pos.vzdir,
        ref_dir: pos.vxdir,
        major_radius,
        minor_radius,
    })
}

// =========================================================================
// OCCT Adaptor3d_Surface evaluation stand-ins.  The ChFiKPart corner /
// sphere entries evaluate the support surfaces S1/S2 through
// Value / D0 / D1 / Parameters; rcad evaluates the analytic Surface3
// through the same ElSLib formulas on the surface frame.
// =========================================================================

/// OCCT gp_Ax3 of the support surface (S->Plane().Position() etc.).
/// The rcad analytic payloads carry the full frame (origin/axis/ref_dir and
/// an optional explicit Y).
pub fn surface3_ax3(s: &Surface3) -> GpAx3 {
    match s {
        Surface3::Plane(pl) => GpAx3 {
            location: pl.origin,
            vxdir: pl.u_dir,
            vydir: pl.v_dir,
            vzdir: pl.normal,
        },
        Surface3::Cylinder(cy) => GpAx3 {
            location: cy.origin,
            vxdir: cy.ref_dir,
            vydir: cy.y_axis(),
            vzdir: cy.axis,
        },
        Surface3::Sphere(sp) => GpAx3 {
            location: sp.center,
            vxdir: sp.ref_dir,
            vydir: sp.axis.cross(sp.ref_dir).normalize(),
            vzdir: sp.axis,
        },
        Surface3::Cone(co) => GpAx3 {
            location: co.apex,
            vxdir: co.ref_dir,
            vydir: co.axis.cross(co.ref_dir).normalize(),
            vzdir: co.axis,
        },
        Surface3::Torus(to) => GpAx3 {
            location: to.center,
            vxdir: to.ref_dir,
            vydir: to.axis.cross(to.ref_dir).normalize(),
            vzdir: to.axis,
        },
        _ => panic!("Standard_NotImplemented: surface kind not carried by ChFiKPart"),
    }
}

/// OCCT Adaptor3d_Surface::Value(U, V) / D0.
pub fn surface3_d0(s: &Surface3, u: f64, v: f64) -> DVec3 {
    let pos = surface3_ax3(s);
    match s {
        Surface3::Plane(_) => elslib_plane_d0(u, v, &pos),
        Surface3::Cylinder(cy) => {
            let (_p, _vu, _vv) = elslib_cylinder_d1(u, v, &pos, cy.radius);
            _p
        }
        Surface3::Sphere(sp) => elslib_sphere_d0(u, v, &pos, sp.radius),
        Surface3::Cone(co) => {
            let (p, _vu, _vv) = elslib_cone_d1(u, v, &pos, co.radius, co.half_angle_rad);
            p
        }
        Surface3::Torus(to) => elslib_torus_d0(u, v, &pos, to.major_radius, to.minor_radius),
        _ => panic!("Standard_NotImplemented: surface kind not carried by ChFiKPart"),
    }
}

/// OCCT Adaptor3d_Surface::D1(U, V, P, D1U, D1V).
pub fn surface3_d1(s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
    let pos = surface3_ax3(s);
    match s {
        Surface3::Plane(_) => elslib_plane_d1(u, v, &pos),
        Surface3::Cylinder(cy) => elslib_cylinder_d1(u, v, &pos, cy.radius),
        Surface3::Sphere(sp) => elslib_sphere_d1(u, v, &pos, sp.radius),
        Surface3::Cone(co) => elslib_cone_d1(u, v, &pos, co.radius, co.half_angle_rad),
        Surface3::Torus(to) => elslib_torus_d1(u, v, &pos, to.major_radius, to.minor_radius),
        _ => panic!("Standard_NotImplemented: surface kind not carried by ChFiKPart"),
    }
}

/// OCCT ElSLib::Parameters(S, P, U, V) — dispatch on the surface kind
/// (PlaneParameters / CylinderParameters / ConeParameters / SphereParameters
/// / TorusParameters).
pub fn elslib_surface_parameters(s: &Surface3, p: DVec3) -> (f64, f64) {
    let pos = surface3_ax3(s);
    match s {
        Surface3::Plane(_) => elslib_plane_parameters(&pos, p),
        Surface3::Cylinder(cy) => elslib_cylinder_parameters(&pos, cy.radius, p),
        Surface3::Sphere(sp) => elslib_sphere_parameters(&pos, sp.radius, p),
        Surface3::Cone(co) => {
            elslib_cone_parameters(&pos, co.radius, co.half_angle_rad, p)
        }
        Surface3::Torus(to) => {
            elslib_torus_parameters(&pos, to.major_radius, to.minor_radius, p)
        }
        _ => panic!("Standard_NotImplemented: surface kind not carried by ChFiKPart"),
    }
}

/// OCCT GeomAbs_SurfaceType comparison support — the kind names used by the
/// ChFiKPart dispatch.  OCCT dispatches on S->GetType().
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    Other,
}

/// OCCT Adaptor3d_Surface::GetType() stand-in.
pub fn surface3_kind(s: &Surface3) -> SurfaceKind {
    match s {
        Surface3::Plane(_) => SurfaceKind::Plane,
        Surface3::Cylinder(_) => SurfaceKind::Cylinder,
        Surface3::Cone(_) => SurfaceKind::Cone,
        Surface3::Sphere(_) => SurfaceKind::Sphere,
        Surface3::Torus(_) => SurfaceKind::Torus,
        _ => SurfaceKind::Other,
    }
}

// =========================================================================
// 2D primitives: gp_Pnt2d -> DVec2, gp_Lin2d -> geom::Line2d,
// gp_Circ2d (gp_Ax22d) -> geom::Circle2d.
// =========================================================================

/// OCCT gp_Lin2d(P, D) as a pcurve payload (Geom2d_Line).
pub fn curve2d_line(origin: DVec2, direction: DVec2) -> rcad_kernel::geom::Curve2d {
    rcad_kernel::geom::Curve2d::Line(Line2d {
        origin,
        direction: direction.normalize_or_zero(),
    })
}

/// OCCT gp_Ax22d — 2D frame of a gp_Circ2d.
#[derive(Debug, Clone, Copy)]
pub struct GpAx22d {
    pub location: DVec2,
    pub vx: DVec2,
    pub vy: DVec2,
}

impl GpAx22d {
    /// OCCT gp_Ax22d(P, X, Y).
    pub fn new(p: DVec2, x: DVec2, y: DVec2) -> Self {
        GpAx22d {
            location: p,
            vx: x,
            vy: y,
        }
    }

    /// OCCT Geom2d_Circle(Ax22d, Radius) -> rcad pcurve payload.
    pub fn to_circle2d(&self, radius: f64) -> rcad_kernel::geom::Curve2d {
        rcad_kernel::geom::Curve2d::Circle(Circle2d {
            center: self.location,
            x_dir: self.vx,
            y_dir: self.vy,
            radius,
        })
    }
}

/// OCCT TopAbs::Reverse(O) (TopAbs.hxx).
pub fn topabs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}
