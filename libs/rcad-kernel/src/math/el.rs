//! OCCT ElCLib + ElSLib: elementary curve and surface evaluation.
//!
//! Analytic evaluation of elementary curves (line, circle, ellipse, hyperbola,
//! parabola) and surfaces (plane, cylinder, sphere, cone, torus).
//!
//! OCCT source: src/FoundationClasses/TKMath/ElCLib/ElCLib.cxx
//!             src/FoundationClasses/TKMath/ElSLib/ElSLib.cxx

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

/// OCCT ElCLib::InPeriod (ElCLib.cxx L95-111) — the value of U in the
/// periodic range [UFirst, ULast].
pub fn in_period(u: f64, ufirst: f64, ulast: f64) -> f64 {
    // In order to avoid FLT_Overflow exception.
    if !u.is_finite() || !ufirst.is_finite() || !ulast.is_finite() {
        return u;
    }

    let period = ulast - ufirst;

    // OCCT: aPeriod < Epsilon(theULast), Epsilon(V) = relative machine eps.
    if period < f64::EPSILON * ulast.abs() {
        return u;
    }

    (ufirst).max(u + period * ((ufirst - u) / period).ceil())
}
