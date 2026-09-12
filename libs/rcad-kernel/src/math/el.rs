//! OCCT ElCLib + ElSLib: elementary curve and surface evaluation.
//!
//! Analytic evaluation of elementary curves (line, circle, ellipse, hyperbola,
//! parabola) and surfaces (plane, cylinder, sphere, cone, torus).
//!
//! OCCT source: src/FoundationClasses/TKMath/ElCLib/ElCLib.cxx
//!             src/FoundationClasses/TKMath/ElSLib/ElSLib.cxx

use crate::core::precision::is_infinite_value;
use glam::DVec3;
use std::f64::consts::PI;

// ══════════════════════════════════════════════════════════════════════════
// ElCLib — elementary curve evaluation
// ══════════════════════════════════════════════════════════════════════════

/// Line: P(t) = origin + t · direction
pub fn elclib_line_value(t: f64, origin: DVec3, direction: DVec3) -> DVec3 {
    origin + t * direction
}

pub fn elclib_line_d1(t: f64, origin: DVec3, direction: DVec3) -> (DVec3, DVec3) {
    (origin + t * direction, direction)
}

/// Circle: P(u) = center + radius · (cos(u)·xDir + sin(u)·yDir)
pub fn elclib_circle_value(u: f64, center: DVec3, x_dir: DVec3, y_dir: DVec3, radius: f64) -> DVec3 {
    center + radius * (u.cos() * x_dir + u.sin() * y_dir)
}

pub fn elclib_circle_d1(u: f64, center: DVec3, x_dir: DVec3, y_dir: DVec3, radius: f64) -> (DVec3, DVec3) {
    let p = center + radius * (u.cos() * x_dir + u.sin() * y_dir);
    let d = radius * (-u.sin() * x_dir + u.cos() * y_dir);
    (p, d)
}

pub fn elclib_circle_d2(u: f64, center: DVec3, x_dir: DVec3, y_dir: DVec3, radius: f64) -> (DVec3, DVec3, DVec3) {
    let (su, cu) = u.sin_cos();
    let p = center + radius * (cu * x_dir + su * y_dir);
    let d1 = radius * (-su * x_dir + cu * y_dir);
    let d2 = -radius * (cu * x_dir + su * y_dir); // = -(p - center)
    (p, d1, d2)
}

/// Ellipse: P(u) = center + majorR·cos(u)·majorDir + minorR·sin(u)·(normal × majorDir)
pub fn elclib_ellipse_value(
    u: f64, center: DVec3, major_dir: DVec3, normal: DVec3,
    major_radius: f64, minor_radius: f64,
) -> DVec3 {
    let y_ax = normal.cross(major_dir).normalize();
    center + major_radius * u.cos() * major_dir + minor_radius * u.sin() * y_ax
}

pub fn elclib_ellipse_d1(
    u: f64, center: DVec3, major_dir: DVec3, normal: DVec3,
    major_radius: f64, minor_radius: f64,
) -> (DVec3, DVec3) {
    let y_ax = normal.cross(major_dir).normalize();
    let (su, cu) = u.sin_cos();
    let p = center + major_radius * cu * major_dir + minor_radius * su * y_ax;
    let d1 = -major_radius * su * major_dir + minor_radius * cu * y_ax;
    (p, d1)
}

/// Hyperbola: P(u) = center + majorR·cosh(u)·majorDir + minorR·sinh(u)·(normal × majorDir)
pub fn elclib_hyperbola_value(
    u: f64, center: DVec3, major_dir: DVec3, normal: DVec3,
    major_radius: f64, minor_radius: f64,
) -> DVec3 {
    let y_ax = normal.cross(major_dir).normalize();
    center + major_radius * u.cosh() * major_dir + minor_radius * u.sinh() * y_ax
}

/// Parabola: P(u) = vertex + u²/(4f)·axisDir + u·(normal × axisDir)
pub fn elclib_parabola_value(
    u: f64, vertex: DVec3, axis_dir: DVec3, normal: DVec3, focal: f64,
) -> DVec3 {
    let y_ax = normal.cross(axis_dir).normalize();
    vertex + (u * u) / (4.0 * focal) * axis_dir + u * y_ax
}

// ══════════════════════════════════════════════════════════════════════════
// ElSLib — elementary surface evaluation
// ══════════════════════════════════════════════════════════════════════════

/// Plane: P(u,v) = origin + u·uDir + v·vDir
pub fn elslib_plane_value(u: f64, v: f64, origin: DVec3, u_dir: DVec3, v_dir: DVec3) -> DVec3 {
    origin + u * u_dir + v * v_dir
}

pub fn elslib_plane_d1(u: f64, v: f64, origin: DVec3, u_dir: DVec3, v_dir: DVec3) -> (DVec3, DVec3, DVec3) {
    (origin + u * u_dir + v * v_dir, u_dir, v_dir)
}

/// Cylinder: P(u,v) = origin + R·(cos(u)·refDir + sin(u)·(axis×refDir)) + v·axis
/// u ∈ [0, 2π), v ∈ ℝ (height along axis)
pub fn elslib_cylinder_value(
    u: f64, v: f64, origin: DVec3, axis: DVec3, ref_dir: DVec3, radius: f64,
) -> DVec3 {
    let x_ax = ref_dir.normalize_or_zero();
    let y_ax = axis.cross(x_ax).normalize();
    origin + radius * (u.cos() * x_ax + u.sin() * y_ax) + v * axis
}

pub fn elslib_cylinder_d1(
    u: f64, v: f64, origin: DVec3, axis: DVec3, ref_dir: DVec3, radius: f64,
) -> (DVec3, DVec3, DVec3) {
    let x_ax = ref_dir.normalize_or_zero();
    let y_ax = axis.cross(x_ax).normalize();
    let (su, cu) = u.sin_cos();
    let p = origin + radius * (cu * x_ax + su * y_ax) + v * axis;
    let dpu = radius * (-su * x_ax + cu * y_ax);
    (p, dpu, axis)
}

/// Sphere: P(u,v) = center + R·(sin(v)·(cos(u)·refDir + sin(u)·(axis×refDir)) + cos(v)·axis)
/// u ∈ [0, 2π) longitude, v ∈ [0, π] colatitude
pub fn elslib_sphere_value(
    u: f64, v: f64, center: DVec3, axis: DVec3, ref_dir: DVec3, radius: f64,
) -> DVec3 {
    let x_ax = ref_dir.normalize();
    let y_ax = axis.cross(x_ax).normalize();
    center + radius * (v.sin() * (u.cos() * x_ax + u.sin() * y_ax) + v.cos() * axis)
}

pub fn elslib_sphere_d1(
    u: f64, v: f64, center: DVec3, axis: DVec3, ref_dir: DVec3, radius: f64,
) -> (DVec3, DVec3, DVec3) {
    let x_ax = ref_dir.normalize();
    let y_ax = axis.cross(x_ax).normalize();
    let (su, cu) = u.sin_cos();
    let (sv, cv) = v.sin_cos();
    let radial = cu * x_ax + su * y_ax;
    let p = center + radius * (sv * radial + cv * axis);
    let dpu = radius * sv * (-su * x_ax + cu * y_ax);
    let dpv = radius * (cv * radial - sv * axis);
    (p, dpu, dpv)
}

/// Cone: P(u,v) = apex + axial_from_slant(v)·axis + radius_at_slant(v)·(cos(u)·xAx + sin(u)·yAx)
pub fn elslib_cone_value(
    u: f64, v: f64, apex: DVec3, axis: DVec3, half_angle: f64, radius: f64,
) -> DVec3 {
    let x_ax = crate::geom::any_perpendicular(axis);
    let y_ax = axis.cross(x_ax).normalize();
    let rad_at_v = radius + v * half_angle.tan();
    let ax_at_v = v; // reference circle at v=0, axial offset = v·cos(α)
    let axial = ax_at_v * half_angle.cos();
    apex + axial * axis + rad_at_v * (u.cos() * x_ax + u.sin() * y_ax)
}

/// Torus: P(u,v) = center + (R + r·cos(v))·(cos(u)·xDir + sin(u)·yDir) + r·sin(v)·axis
pub fn elslib_torus_value(
    u: f64, v: f64, center: DVec3, axis: DVec3, major_radius: f64, minor_radius: f64,
) -> DVec3 {
    let ref_dir = if axis.x.abs() > 1.0 - 1e-12 { DVec3::Z } else { DVec3::X };
    let x_ax = (ref_dir - axis * ref_dir.dot(axis)).normalize_or_zero();
    let y_ax = axis.cross(x_ax).normalize();
    let radial = u.cos() * x_ax + u.sin() * y_ax;
    center + (major_radius + minor_radius * v.cos()) * radial + minor_radius * v.sin() * axis
}

/// OCCT ElCLib::Parameter(const gp_Lin2d&, const gp_Pnt2d&) — the parameter
/// of the 2D point on the 2D line: (P - Location)·Direction.
pub fn elclib_line_parameter_2d(p: glam::DVec2, origin: glam::DVec2, direction: glam::DVec2) -> f64 {
    (p - origin).dot(direction)
}

// ---------------------------------------------------------------------------
// OCCT ElCLib::To3d (ElCLib.cxx L1339-1440) — the 2D payloads lifted into the
// 3D frame of a gp_Ax2 (passed as the kernel `Ax2`; its y_direction is the
// `N ^ X` derived by the 3-argument gp_Ax2 constructor).
// ---------------------------------------------------------------------------

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Pnt2d&) (ElCLib.cxx L1339-1346):
/// `Vxy = XDirection*P.X + YDirection*P.Y + Location`.
pub fn elclib_to3d_pnt(pos: &crate::math::gp::Ax2, p: glam::DVec2) -> DVec3 {
    pos.x_direction * p.x + pos.y_direction * p.y + pos.location
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Vec2d&) (ElCLib.cxx L1360-1370).
pub fn elclib_to3d_vec(pos: &crate::math::gp::Ax2, v: glam::DVec2) -> DVec3 {
    let vx = pos.x_direction * v.x;
    let vy = pos.y_direction * v.y;
    vx + vy
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Ax22d&) (ElCLib.cxx L1381-1388):
/// `gp_Ax2(P, VX.Crossed(VY), VX)` — the 3-argument gp_Ax2 constructor
/// re-orthogonalizes VX against the new main direction and derives Y.
pub fn elclib_to3d_ax22d(
    pos: &crate::math::gp::Ax2,
    location: glam::DVec2,
    x_dir: glam::DVec2,
    y_dir: glam::DVec2,
) -> crate::math::gp::Ax2 {
    let p = elclib_to3d_pnt(pos, location);
    let vx = elclib_to3d_vec(pos, x_dir);
    let vy = elclib_to3d_vec(pos, y_dir);
    crate::math::gp::Ax2::new(p, vx.cross(vy), vx)
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Lin2d&) (ElCLib.cxx L1391-1394).
pub fn elclib_to3d_line(
    pos: &crate::math::gp::Ax2,
    origin: glam::DVec2,
    direction: glam::DVec2,
) -> crate::geom::Line3 {
    crate::geom::Line3::new(elclib_to3d_pnt(pos, origin), elclib_to3d_vec(pos, direction))
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Circ2d&) (ElCLib.cxx L1396-1399):
/// `gp_Circ(To3d(Pos, C.Axis()), C.Radius())`.
pub fn elclib_to3d_circle(
    pos: &crate::math::gp::Ax2,
    c: &crate::geom::Circle2d,
) -> crate::geom::Circle3 {
    let ax = elclib_to3d_ax22d(pos, c.center, c.x_dir, c.y_dir);
    crate::geom::Circle3 {
        center: ax.location,
        normal: ax.direction,
        x_dir: ax.x_direction,
        // gp_Circ keeps the Ax2 the 3-argument constructor built: Y = N ^ X.
        y_dir: ax.direction.cross(ax.x_direction),
        radius: c.radius,
    }
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Elips2d&) (ElCLib.cxx L1401-1404).
pub fn elclib_to3d_ellipse(
    pos: &crate::math::gp::Ax2,
    e: &crate::geom::Ellipse2d,
) -> crate::geom::Ellipse3 {
    let ax = elclib_to3d_ax22d(pos, e.center, e.major_dir, e.minor_dir);
    crate::geom::Ellipse3 {
        center: ax.location,
        normal: ax.direction,
        major_dir: ax.x_direction,
        major_radius: e.major_radius,
        minor_radius: e.minor_radius,
    }
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Hypr2d&) (ElCLib.cxx L1406-1410):
/// `gp_Hypr(To3d(Pos, H.Axis()), H.MajorRadius(), H.MinorRadius())`. The
/// rcad `Hyperbola2d` carries no minor direction, so the Ax22d Y direction
/// uses the codebase conic convention (-major.y, major.x).
pub fn elclib_to3d_hyperbola(
    pos: &crate::math::gp::Ax2,
    h: &crate::geom::Hyperbola2d,
) -> crate::geom::Hyperbola3 {
    let y_dir = glam::DVec2::new(-h.major_dir.y, h.major_dir.x);
    let ax = elclib_to3d_ax22d(pos, h.center, h.major_dir, y_dir);
    crate::geom::Hyperbola3 {
        center: ax.location,
        normal: ax.direction,
        major_dir: ax.x_direction,
        semi_major: h.semi_major,
        semi_minor: h.semi_minor,
    }
}

/// OCCT ElCLib::To3d(const gp_Ax2&, const gp_Parab2d&) (ElCLib.cxx L1412-1415):
/// `gp_Parab(To3d(Pos, Prb.Axis()), Prb.Focal())`. Same convention note as
/// the hyperbola overload.
pub fn elclib_to3d_parabola(
    pos: &crate::math::gp::Ax2,
    p: &crate::geom::Parabola2d,
) -> crate::geom::Parabola3 {
    let y_dir = glam::DVec2::new(-p.axis_dir.y, p.axis_dir.x);
    let ax = elclib_to3d_ax22d(pos, p.origin, p.axis_dir, y_dir);
    crate::geom::Parabola3 {
        vertex: ax.location,
        normal: ax.direction,
        axis_dir: ax.x_direction,
        focal_param: p.focal_param,
    }
}

/// OCCT ElCLib::InPeriod (ElCLib.cxx L95-111) — the value of U in the
/// periodic range [UFirst, ULast].
pub fn in_period(u: f64, ufirst: f64, ulast: f64) -> f64 {
    // OCCT L101-105: Precision::IsInfinite on all three arguments — "In
    // order to avoid FLT_Overflow exception".  IsInfinite is the ±1e100
    // threshold (Precision.hxx L350-355), not the IEEE finiteness test, so
    // the ±Precision::Infinite() domain sentinels returned by the unbounded
    // curve producers take the early return exactly like OCCT.
    if is_infinite_value(u) || is_infinite_value(ufirst) || is_infinite_value(ulast) {
        return u;
    }

    let period = ulast - ufirst;

    // OCCT: aPeriod < Epsilon(theULast), Epsilon(V) = relative machine eps.
    if period < f64::EPSILON * ulast.abs() {
        return u;
    }

    (ufirst).max(u + period * ((ufirst - u) / period).ceil())
}

// ============================================================================
// ElSLib::Parameters — the inverse parameterisation (3D point -> UV)
// (ElSLib.cxx L1547-1641).  The OCCT frames arrive via gp_Trsf
// SetTransformation(gp_Ax3); here the caller supplies the frame axes
// directly (local X = (P-O)·xDir, Y = (P-O)·yDir, Z = (P-O)·zDir).
// ============================================================================

const PIPI: f64 = std::f64::consts::PI + std::f64::consts::PI;
/// OCCT ElSLib.cxx NEGATIVE_RESOLUTION = -Precision::Computational().
const NEGATIVE_RESOLUTION: f64 = -f64::EPSILON;
/// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 1e-15;

/// OCCT ElSLib.cxx normalizeAngle (L42-56) — normalize to [0, 2·PI],
/// preserving the exact 2·PI seam value.
fn normalize_angle(angle: &mut f64) {
    while *angle < NEGATIVE_RESOLUTION {
        *angle += PIPI;
    }
    while *angle > PIPI * (1.0 + GP_RESOLUTION) {
        *angle -= PIPI;
    }
    if *angle < 0.0 {
        *angle = 0.0;
    }
}

/// OCCT ElSLib::PlaneParameters (ElSLib.cxx L1547-1556) — U and V are the
/// local coordinates of P in the plane frame.
pub fn elslib_plane_parameters(
    p: DVec3,
    origin: DVec3,
    u_dir: DVec3,
    v_dir: DVec3,
) -> (f64, f64) {
    let d = p - origin;
    (d.dot(u_dir), d.dot(v_dir))
}

/// OCCT ElSLib::CylinderParameters (ElSLib.cxx L1558-1570) — U = atan2 of the
/// local point, V = local axial coordinate.  The radius argument is unused in
/// OCCT (kept for signature parity).
pub fn elslib_cylinder_parameters(
    p: DVec3,
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    axis: DVec3,
    _radius: f64,
) -> (f64, f64) {
    let d = p - origin;
    let x = d.dot(x_dir);
    let y = d.dot(y_dir);
    let mut u = y.atan2(x);
    normalize_angle(&mut u);
    let v = d.dot(axis);
    (u, v)
}

/// OCCT ElSLib::ConeParameters (ElSLib.cxx L1574-1613) — U from the local
/// angle with the wrong-side-of-apex guards, V measured along the cone
/// generatrix direction: V = sin(SAngle)·(x·cosU + y·sinU - R) +
/// cos(SAngle)·z.
pub fn elslib_cone_parameters(
    p: DVec3,
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    axis: DVec3,
    radius: f64,
    semi_angle: f64,
) -> (f64, f64) {
    let d = p - origin;
    let x = d.dot(x_dir);
    let y = d.dot(y_dir);
    let z = d.dot(axis);

    let mut u;
    if x.abs() < GP_RESOLUTION && y.abs() < GP_RESOLUTION {
        // The point is on the cone axis (apex).
        u = 0.0;
    } else if -radius > z * semi_angle.tan() {
        // The point is at the wrong side of the apex.
        u = (-y).atan2(-x);
    } else {
        u = y.atan2(x);
    }
    normalize_angle(&mut u);

    let v = semi_angle.sin() * (x * u.cos() + y * u.sin() - radius) + semi_angle.cos() * z;
    (u, v)
}

/// OCCT ElSLib::SphereParameters (ElSLib.cxx L1615-1641) — V = latitude from
/// the local polar distance, U = longitude; degenerate on-axis points get
/// V = ±PI/2 and U = 0.
pub fn elslib_sphere_parameters(
    p: DVec3,
    center: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    axis: DVec3,
) -> (f64, f64) {
    let d = p - center;
    let x = d.dot(x_dir);
    let y = d.dot(y_dir);
    let z = d.dot(axis);
    let l = (x * x + y * y).sqrt();
    if l < GP_RESOLUTION {
        // Point on the Z axis of the sphere.
        let v = if z > 0.0 {
            std::f64::consts::FRAC_PI_2
        } else {
            -std::f64::consts::FRAC_PI_2
        };
        (0.0, v)
    } else {
        let v = (z / l).atan();
        let mut u = y.atan2(x);
        normalize_angle(&mut u);
        (u, v)
    }
}

/// OCCT ElSLib::TorusParameters (ElSLib.cxx L1646-1697) — U = the major angle
/// around the torus axis (with the Major < Minor branch that flips the
/// nearest side), V = the minor angle of the projected point.
#[allow(clippy::many_single_char_names)]
pub fn elslib_torus_parameters(
    p: DVec3,
    center: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    axis: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> (f64, f64) {
    let d = p - center;
    let x = d.dot(x_dir);
    let y = d.dot(y_dir);
    let z = d.dot(axis);

    // All that to process the case of Major < Minor.
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
    normalize_angle(&mut u);
    let cosu = u.cos();
    let sinu = u.sin();
    // dx = (cosU, sinU, 0); V = dx.AngleWithRef(dP, dx ^ DZ) =
    // atan2(z, cosU·(x - R·cosU) + sinU·y).
    let dpx = x - major_radius * cosu;
    let dpy = y - major_radius * sinu;
    let a_mag = (dpx * dpx + dpy * dpy + z * z).sqrt();
    let mut v = if a_mag <= GP_RESOLUTION {
        0.0
    } else {
        (z).atan2(cosu * dpx + sinu * dpy)
    };
    normalize_angle(&mut v);
    (u, v)
}
