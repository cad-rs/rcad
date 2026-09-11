//! OCCT ElSLib isoparametric isoline constructors (TKMath/ElSLib) and the gp
//! frame operations they need — the support layer consumed by
//! Adaptor3d_CurveOnSurface::EvalKPart (Adaptor3d_CurveOnSurface.cxx
//! L1552-1732).
//!
//! 1:1 translation of ElSLib.cxx L1716-1854 (CylinderUIso / ConeUIso /
//! SphereUIso / TorusUIso / CylinderVIso / ConeVIso / SphereVIso /
//! TorusVIso), of the referenced ElSLib::CylinderD1 / ConeD1 derivative
//! forms (ElSLib.cxx L687-755) and of the gp_Ax2 / gp_Dir2d operations they
//! perform (gp_Ax2.hxx, gp_Dir2d.cxx, gp_Trsf.cxx / gp_Mat.cxx rotation).
//!
//! Architecture mapping: the OCCT `gp_Ax3` surface frame (Location /
//! XDirection / YDirection / Direction / Directness) is carried by the rcad
//! payload structs as (location or apex, ref_dir, derived or explicit y,
//! axis); [`Ax3View`] reconstructs that frame from a payload so the iso
//! bodies keep the OCCT code shape.

use glam::{DVec2, DVec3};

use crate::geom::{Circle3, Line3};

// ---------------------------------------------------------------------------
// Ax3View — the OCCT gp_Ax3 frame view over a payload
// ---------------------------------------------------------------------------

/// The gp_Ax3 frame fields consumed by the ElSLib iso constructors.
///
/// OCCT keeps XDirection / YDirection / Direction and the Directness flag;
/// `direct` is `(XDirection ^ YDirection) . Direction > 0` (gp_Ax3.hxx
/// L486-493 reads it in Ax2()).
#[derive(Debug, Clone, Copy)]
pub struct Ax3View {
    /// OCCT: gp_Pnt Location().
    pub location: DVec3,
    /// OCCT: gp_Dir Direction() — the main direction.
    pub direction: DVec3,
    /// OCCT: gp_Dir XDirection().
    pub x_direction: DVec3,
    /// OCCT: gp_Dir YDirection().
    pub y_direction: DVec3,
    /// OCCT: Direct() — the frame handedness.
    pub direct: bool,
}

impl Ax3View {
    /// The right-handed frame of a payload with the OCCT default
    /// YDirection = Direction ^ XDirection.
    pub fn from_axes(location: DVec3, direction: DVec3, x_direction: DVec3) -> Self {
        let direction = direction.normalize_or_zero();
        let x_direction = x_direction.normalize_or_zero();
        let y_direction = direction.cross(x_direction).normalize_or_zero();
        Ax3View {
            location,
            direction,
            x_direction,
            y_direction,
            direct: true,
        }
    }

    /// The explicit-Y (possibly left-handed) frame of a payload that carries
    /// both directions (the swept-lateral cylinder frame encoding).
    pub fn with_y_dir(location: DVec3, direction: DVec3, x_direction: DVec3, y: DVec3) -> Self {
        let direction = direction.normalize_or_zero();
        let x_direction = x_direction.normalize_or_zero();
        let y_direction = y.normalize_or_zero();
        let direct = x_direction.cross(y_direction).dot(direction) > 0.0;
        Ax3View {
            location,
            direction,
            x_direction,
            y_direction,
            direct,
        }
    }

    /// OCCT gp_Ax3::Ax2() (gp_Ax3.hxx L486-493): the Ax2 form with main
    /// direction reversed for the non-direct (left-handed) frame.
    pub fn ax2(&self) -> Ax2View {
        let zz = if self.direct {
            self.direction
        } else {
            -self.direction
        };
        Ax2View {
            location: self.location,
            direction: zz,
            x_direction: self.x_direction,
        }
    }
}

// ---------------------------------------------------------------------------
// Ax2View — the OCCT gp_Ax2 frame (Location / Direction / XDirection; the
// YDirection is Direction ^ XDirection)
// ---------------------------------------------------------------------------

/// The gp_Ax2 frame fields consumed by the iso constructors.
#[derive(Debug, Clone, Copy)]
pub struct Ax2View {
    /// OCCT: gp_Pnt Location().
    pub location: DVec3,
    /// OCCT: gp_Dir Direction() — the main direction.
    pub direction: DVec3,
    /// OCCT: gp_Dir XDirection().
    pub x_direction: DVec3,
}

impl Ax2View {
    /// OCCT gp_Ax2::YDirection — Direction ^ XDirection
    /// (gp_Ax2.hxx L73-80: vydir.Cross(vxdir) on the orthonormalized frame).
    pub fn y_direction(&self) -> DVec3 {
        self.direction.cross(self.x_direction).normalize_or_zero()
    }

    /// OCCT gp_Ax2::Translate(gp_Vec) — Location += V.
    pub fn translate(&mut self, v: DVec3) {
        self.location += v;
    }

    /// OCCT gp_Ax2::Rotate(gp_Ax1, Ang) (gp_Ax2.hxx L301-309): rotates the
    /// location, the X and Y directions, then re-derives the main direction
    /// as XDirection ^ YDirection.
    pub fn rotate(&mut self, a1_loc: DVec3, a1_dir: DVec3, ang: f64) {
        let a_temp = xyz_rotate(a1_loc, self.location, a1_dir, ang);
        let vx = xyz_rotate_dir(self.x_direction, a1_dir, ang);
        let vy = xyz_rotate_dir(self.y_direction(), a1_dir, ang);
        self.location = a_temp;
        self.x_direction = vx;
        self.direction = vx.cross(vy).normalize_or_zero();
    }

    /// OCCT gp_Ax2::SetDirection(gp_Dir) (gp_Ax2.hxx L548-571).  The rcad
    /// view keeps the YDirection derived (Direction ^ XDirection), which
    /// reproduces the OCCT vydir assignments in every branch.
    pub fn set_direction(&mut self, the_v: DVec3) {
        let the_v = the_v.normalize_or_zero();
        let vxdir = self.x_direction;
        let vydir = self.y_direction();
        let old_main = self.direction;
        let a = the_v.dot(vxdir);
        if ((a.abs() - 1.0).abs()) <= crate::core::precision::ANGULAR {
            if a > 0.0 {
                // OCCT: vxdir = vydir; vydir = axis.Direction().
                self.x_direction = vydir;
            } else {
                // OCCT: vxdir = axis.Direction().
                self.x_direction = old_main;
            }
            self.direction = the_v;
        } else {
            // OCCT: axis.SetDirection(theV);
            //       vxdir = theV.CrossCrossed(vxdir, theV)  (theV ^ (vxdir ^ theV));
            //       vydir = theV.Crossed(vxdir).
            self.direction = the_v;
            self.x_direction = the_v.cross(vxdir.cross(the_v)).normalize_or_zero();
        }
    }
}

/// OCCT gp_Circ(gp_Ax2, Radius) — the circle payload from an Ax2 frame.
pub fn circ_from_ax2(axes: &Ax2View, radius: f64) -> Circle3 {
    Circle3 {
        center: axes.location,
        normal: axes.direction,
        x_dir: axes.x_direction,
        y_dir: axes.y_direction(),
        radius,
    }
}

/// OCCT gp_Circ::Rotate(A1, Ang) (gp_Circ.hxx Rotate -> gp_Ax2::Rotate,
/// gp_Ax2.hxx L301-309) — the position frame rotated, the radius kept.
pub fn circ_rotated(circ: &Circle3, a1_loc: DVec3, a1_dir: DVec3, ang: f64) -> Circle3 {
    let mut axes = Ax2View {
        location: circ.center,
        direction: circ.normal,
        x_direction: circ.x_dir,
    };
    axes.rotate(a1_loc, a1_dir, ang);
    circ_from_ax2(&axes, circ.radius)
}

/// OCCT EvalKPart opposite-direction arm: `gp_Ax2 Ax = myCirc.Position();
/// Ax.SetDirection(Ax.Direction().Reversed()); myCirc.SetPosition(Ax)`
/// (Adaptor3d_CurveOnSurface.cxx L1597-1602 and peers) — the SetDirection
/// step is gp_Ax2::SetDirection (gp_Ax2.hxx L548-571).
pub fn circ_with_direction_reversed(circ: &Circle3) -> Circle3 {
    let mut axes = Ax2View {
        location: circ.center,
        direction: circ.normal,
        x_direction: circ.x_dir,
    };
    axes.set_direction(-axes.direction);
    circ_from_ax2(&axes, circ.radius)
}

// ---------------------------------------------------------------------------
// gp_XYZ rotation (gp_Trsf::SetRotation + gp_Mat::SetRotation)
// ---------------------------------------------------------------------------

/// OCCT gp_Mat::SetRotation(Axis, Ang) (gp_Mat.cxx L122-159) applied as
/// gp_Trsf::Transforms (gp_Trsf.cxx L90-99): R*(P - A1Loc) + A1Loc.
pub fn xyz_rotate(a1_loc: DVec3, p: DVec3, a1_dir: DVec3, ang: f64) -> DVec3 {
    let a_v = a1_dir.normalize_or_zero();
    let (a, b, c) = (a_v.x, a_v.y, a_v.z);
    let a_cos = ang.cos();
    let a_sin = ang.sin();
    let a_om_cos = 1.0 - a_cos;
    let (a2, b2, c2) = (a * a, b * b, c * c);
    let (ab, ac, bc) = (a * b, a * c, b * c);
    // Row-major rotation matrix, the OCCT myMat layout.
    let m00 = 1.0 + a_om_cos * (-(b2 + c2));
    let m01 = a_om_cos * ab - a_sin * c;
    let m02 = a_om_cos * ac + a_sin * b;
    let m10 = a_om_cos * ab + a_sin * c;
    let m11 = 1.0 + a_om_cos * (-(a2 + c2));
    let m12 = a_om_cos * bc - a_sin * a;
    let m20 = a_om_cos * ac - a_sin * b;
    let m21 = a_om_cos * bc + a_sin * a;
    let m22 = 1.0 + a_om_cos * (-(a2 + b2));
    let d = p - a1_loc;
    DVec3::new(
        m00 * d.x + m01 * d.y + m02 * d.z,
        m10 * d.x + m11 * d.y + m12 * d.z,
        m20 * d.x + m21 * d.y + m22 * d.z,
    ) + a1_loc
}

/// OCCT gp_Vec::Rotate(A1, Ang) (gp_Vec.hxx L521-526): the vectorial part
/// only (no translation).
pub fn xyz_rotate_dir(v: DVec3, a1_dir: DVec3, ang: f64) -> DVec3 {
    let a_v = a1_dir.normalize_or_zero();
    let (a, b, c) = (a_v.x, a_v.y, a_v.z);
    let a_cos = ang.cos();
    let a_sin = ang.sin();
    let a_om_cos = 1.0 - a_cos;
    let (a2, b2, c2) = (a * a, b * b, c * c);
    let (ab, ac, bc) = (a * b, a * c, b * c);
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
        m00 * v.x + m01 * v.y + m02 * v.z,
        m10 * v.x + m11 * v.y + m12 * v.z,
        m20 * v.x + m21 * v.y + m22 * v.z,
    )
}

// ---------------------------------------------------------------------------
// ElSLib::CylinderD1 / ConeD1 (the derivative forms used by the U-iso)
// ---------------------------------------------------------------------------

/// OCCT ElSLib::CylinderD1(U, V, Pos, Radius, P, Vu, Vv)
/// (ElSLib.cxx L732-755) — the point and the V derivative consumed by
/// CylinderUIso.
fn cylinder_d1(u: f64, v: f64, pos: &Ax3View, radius: f64) -> (DVec3, DVec3) {
    let x_dir = pos.x_direction;
    let y_dir = pos.y_direction;
    let z_dir = pos.direction;
    let p_loc = pos.location;
    let a1 = radius * u.cos();
    let a2 = radius * u.sin();
    let p = DVec3::new(
        a1 * x_dir.x + a2 * y_dir.x + v * z_dir.x + p_loc.x,
        a1 * x_dir.y + a2 * y_dir.y + v * z_dir.y + p_loc.y,
        a1 * x_dir.z + a2 * y_dir.z + v * z_dir.z + p_loc.z,
    );
    // OCCT: Vv = ZDir (the V derivative; the Vu part feeds only P's frame).
    (p, z_dir)
}

/// OCCT ElSLib::ConeD1(U, V, Pos, Radius, SAngle, P, Vu, Vv)
/// (ElSLib.cxx L687-730) — the point and the V derivative consumed by
/// ConeUIso.
fn cone_d1(u: f64, v: f64, pos: &Ax3View, radius: f64, s_angle: f64) -> (DVec3, DVec3) {
    let x_dir = pos.x_direction;
    let y_dir = pos.y_direction;
    let z_dir = pos.direction;
    let p_loc = pos.location;
    let cos_u = u.cos();
    let sin_u = u.sin();
    let cos_a = s_angle.cos();
    let sin_a = s_angle.sin();
    let r = radius + v * sin_a;
    let a3 = v * cos_a;
    let a1 = r * cos_u;
    let a2 = r * sin_u;
    let r1 = sin_a * cos_u;
    let r2 = sin_a * sin_u;
    let p = DVec3::new(
        a1 * x_dir.x + a2 * y_dir.x + a3 * z_dir.x + p_loc.x,
        a1 * x_dir.y + a2 * y_dir.y + a3 * z_dir.y + p_loc.y,
        a1 * x_dir.z + a2 * y_dir.z + a3 * z_dir.z + p_loc.z,
    );
    let vv = DVec3::new(
        r1 * x_dir.x + r2 * y_dir.x + cos_a * z_dir.x,
        r1 * x_dir.y + r2 * y_dir.y + cos_a * z_dir.y,
        r1 * x_dir.z + r2 * y_dir.z + cos_a * z_dir.z,
    );
    (p, vv)
}

// ---------------------------------------------------------------------------
// The U-iso constructors (ElSLib.cxx L1716-1766)
// ---------------------------------------------------------------------------

/// OCCT ElSLib::PlaneUIso(Pos, U) (ElSLib.cxx L1705-1712): a line through the
/// plane location along YDirection, translated by U * XDirection.
pub fn elslib_plane_u_iso(pos: &Ax3View, u: f64) -> Line3 {
    let origin = pos.location + u * pos.x_direction;
    Line3::new(origin, pos.y_direction)
}

/// OCCT ElSLib::PlaneVIso(Pos, V) (ElSLib.cxx L1770-1777): a line through the
/// plane location along XDirection, translated by V * YDirection.
pub fn elslib_plane_v_iso(pos: &Ax3View, v: f64) -> Line3 {
    let origin = pos.location + v * pos.y_direction;
    Line3::new(origin, pos.x_direction)
}

/// OCCT ElSLib::CylinderUIso(Pos, Radius, U) (ElSLib.cxx L1716-1723).
pub fn elslib_cylinder_u_iso(pos: &Ax3View, radius: f64, u: f64) -> Line3 {
    let (p, dv) = cylinder_d1(u, 0.0, pos, radius);
    Line3::new(p, dv)
}

/// OCCT ElSLib::ConeUIso(Pos, Radius, SAngle, U) (ElSLib.cxx L1727-1734).
pub fn elslib_cone_u_iso(pos: &Ax3View, radius: f64, s_angle: f64, u: f64) -> Line3 {
    let (p, dv) = cone_d1(u, 0.0, pos, radius, s_angle);
    Line3::new(p, dv)
}

/// OCCT ElSLib::SphereUIso(Pos, Radius, U) (ElSLib.cxx L1738-1747).
pub fn elslib_sphere_u_iso(pos: &Ax3View, radius: f64, u: f64) -> Circle3 {
    let dx = pos.x_direction;
    let dy = pos.y_direction;
    let dz = pos.direction;
    let cx = (u.cos() * dx + u.sin() * dy).normalize_or_zero();
    let axes = Ax2View {
        location: pos.location,
        // OCCT: cx.Crossed(dz)
        direction: cx.cross(dz).normalize_or_zero(),
        x_direction: cx,
    };
    circ_from_ax2(&axes, radius)
}

/// OCCT ElSLib::TorusUIso(Pos, MajorRadius, MinorRadius, U)
/// (ElSLib.cxx L1751-1766).
pub fn elslib_torus_u_iso(pos: &Ax3View, major_radius: f64, minor_radius: f64, u: f64) -> Circle3 {
    let dx = pos.x_direction;
    let dy = pos.y_direction;
    let dz = pos.direction;
    let cx = (u.cos() * dx + u.sin() * dy).normalize_or_zero();
    let mut axes = Ax2View {
        location: pos.location,
        direction: cx.cross(dz).normalize_or_zero(),
        x_direction: cx,
    };
    // OCCT: Ve = cx * MajorRadius; axes.Translate(Ve).
    axes.translate(cx * major_radius);
    circ_from_ax2(&axes, minor_radius)
}

// ---------------------------------------------------------------------------
// The V-iso constructors (ElSLib.cxx L1781-1854)
// ---------------------------------------------------------------------------

/// OCCT ElSLib::CylinderVIso(Pos, Radius, V) (ElSLib.cxx L1781-1789).
pub fn elslib_cylinder_v_iso(pos: &Ax3View, radius: f64, v: f64) -> Circle3 {
    let mut axes = pos.ax2();
    // OCCT: Ve = Pos.Direction() * V — the Ax3 main direction, not the Ax2
    // one (they coincide for a direct frame).
    let ve = pos.direction * v;
    axes.translate(ve);
    circ_from_ax2(&axes, radius)
}

/// OCCT ElSLib::ConeVIso(Pos, Radius, SAngle, V) (ElSLib.cxx L1793-1811).
pub fn elslib_cone_v_iso(pos: &Ax3View, radius: f64, s_angle: f64, v: f64) -> Circle3 {
    // OCCT: gp_Ax3 axes(Pos) — a copy; then translate by V*cos(SA)*Dir.
    let mut axes = *pos;
    let ve = pos.direction * (v * s_angle.cos());
    axes.location += ve;
    let mut r = radius + v * s_angle.sin();
    if r < 0.0 {
        // OCCT: axes.XReverse(); axes.YReverse();
        axes.x_direction = -axes.x_direction;
        axes.y_direction = -axes.y_direction;
        r = -r;
    }
    circ_from_ax2(&axes.ax2(), r)
}

/// OCCT ElSLib::SphereVIso(Pos, Radius, V) (ElSLib.cxx L1815-1832).
pub fn elslib_sphere_v_iso(pos: &Ax3View, radius: f64, v: f64) -> Circle3 {
    let mut axes = pos.ax2();
    // OCCT: Ve = Pos.Direction() * (Radius * sin(V)).
    let ve = pos.direction * (radius * v.sin());
    axes.translate(ve);
    let mut radius_iso = radius * v.cos();
    // OCCT #23170: the analytical continuation for |V| > PI/2 — flip the main
    // direction and the radius sign.
    if radius_iso < 0.0 {
        axes.set_direction(-axes.direction);
        radius_iso = -radius_iso;
    }
    circ_from_ax2(&axes, radius_iso)
}

/// OCCT ElSLib::TorusVIso(Pos, MajorRadius, MinorRadius, V)
/// (ElSLib.cxx L1836-1854).
pub fn elslib_torus_v_iso(pos: &Ax3View, major_radius: f64, minor_radius: f64, v: f64) -> Circle3 {
    // OCCT: gp_Ax3 axes = Pos.Ax2() — the Ax2 rebuilt into an Ax3 (direct).
    let mut axes = pos.ax2();
    // OCCT: Ve = Pos.Direction() * (MinorRadius * sin(V)).
    let ve = pos.direction * (minor_radius * v.sin());
    axes.translate(ve);
    let mut r = major_radius + minor_radius * v.cos();
    if r < 0.0 {
        // OCCT: axes.XReverse(); axes.YReverse(); — on the direct Ax3 view
        // this flips the Ax2 XDirection (the YDirection is derived).
        axes.x_direction = -axes.x_direction;
        r = -r;
    }
    circ_from_ax2(&axes, r)
}

// ---------------------------------------------------------------------------
// gp_Dir2d angular predicates (gp_Dir2d.cxx / gp_Dir2d.hxx)
// ---------------------------------------------------------------------------

/// OCCT gp_Dir2d::Angle(Other) (gp_Dir2d.cxx L26-65) — the angular value in
/// [-PI, PI] computed through acos/asin by quadrant.
pub fn dir2d_angle(coord: DVec2, other: DVec2) -> f64 {
    let cosinus = coord.dot(other);
    let sinus = coord.x * other.y - coord.y * other.x;
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        -std::f64::consts::PI - sinus.asin()
    }
}

/// OCCT gp_Dir2d::IsParallel(Other, AngularTolerance)
/// (gp_Dir2d.hxx L422-431).
pub fn dir2d_is_parallel(coord: DVec2, other: DVec2, angular_tolerance: f64) -> bool {
    let mut an_ang = dir2d_angle(coord, other);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Dir2d::IsOpposite(Other, AngularTolerance)
/// (gp_Dir2d.hxx L410-418).
pub fn dir2d_is_opposite(coord: DVec2, other: DVec2, angular_tolerance: f64) -> bool {
    let mut an_ang = dir2d_angle(coord, other);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    std::f64::consts::PI - an_ang <= angular_tolerance
}
