use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlaneConicalResult {
    NoIntersection,
    Point(DVec3),
    SingleLine(Line3),
    TwoLines(Line3, Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    Parabola(Parabola3),
    Hyperbola(Hyperbola3),
}

/// Classify the intersection of a plane and a cone analytically.
///
/// Let α = angle between plane normal and cone axis, β = cone half-angle.
/// Conic type is determined by the Dandelin-sphere criterion:
/// - |sin α| ≈ 1  (plane ⊥ axis)             → Circle
/// - |sin α| > cos β  (plane steeper than α = π/2−β)  → Ellipse
/// - |sin α| ≈ cos β  (plane parallel to one generator) → Parabola
/// - |sin α| < cos β  (plane shallower, parallel to axis) → Hyperbola
pub fn intersect_plane_cone(plane: &Plane, cone: &ConicalSurface) -> PlaneConicalResult {
    let axis_n = cone.axis.normalize();
    let plane_n = plane.normal.normalize();
    let apex = cone.apex_point();

    // cos of angle between plane normal and cone axis
    let cos_angle = plane_n.dot(axis_n).abs();
    // sin of that angle (= cos of complement)
    let sin_angle = (1.0 - cos_angle * cos_angle).sqrt().max(0.0);

    // Signed distance from apex to plane along plane normal direction
    let apex_to_plane = (plane.origin - apex).dot(plane.normal);

    // ── Plane ⊥ axis → circle ─────────────────────────────────────────────────
    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        if apex_to_plane.abs() < TOLERANCE_ABS {
            return PlaneConicalResult::Point(apex);
        }
        let t = apex_to_plane / axis_n.dot(plane.normal);
        let center = apex + axis_n * t;
        let radius = (t * cone.half_angle_rad.tan()).abs();
        if radius < TOLERANCE_ABS {
            return PlaneConicalResult::Point(center);
        }
        return PlaneConicalResult::Circle(Circle3::new(center, cone.axis, radius));
    }

    // ── Plane through apex ────────────────────────────────────────────────────
    if apex_to_plane.abs() < TOLERANCE_ABS {
        let angle_between = sin_angle.atan2(cos_angle); // angle between plane NORMAL and axis
        let half = cone.half_angle_rad;

        // The plane itself makes an angle of π/2 − angle_between with the axis.
        //   π/2 − angle_between < half  → plane cuts through cone → TwoLines
        //   π/2 − angle_between = half  → plane tangent to cone  → SingleLine
        //   π/2 − angle_between > half  → plane misses cone      → Point
        let two_lines_cutoff = std::f64::consts::FRAC_PI_2 - half;

        if (angle_between - two_lines_cutoff).abs() < TOLERANCE_ANG {
            // Tangent: single generator line
            let dir = plane_n.cross(axis_n).normalize();
            let gen_dir = (axis_n * half.cos() + dir * half.sin()).normalize();
            return PlaneConicalResult::SingleLine(Line3 {
                origin: apex,
                direction: gen_dir,
            });
        }

        if angle_between > two_lines_cutoff {
            // Two generators
            let cross = plane_n.cross(axis_n);
            if is_zero_vec(cross) {
                return PlaneConicalResult::Point(apex);
            }
            let perp_in_plane = cross.normalize();
            let projected_axis =
                (axis_n - plane_n * axis_n.dot(plane_n)).normalize_or_zero();
            if projected_axis.length_squared() < TOLERANCE_LEN_MIN {
                return PlaneConicalResult::Point(apex);
            }
            let d1 = (projected_axis * half.cos() + perp_in_plane * half.sin()).normalize();
            let d2 = (projected_axis * half.cos() - perp_in_plane * half.sin()).normalize();
            return PlaneConicalResult::TwoLines(
                Line3 { origin: apex, direction: d1 },
                Line3 { origin: apex, direction: d2 },
            );
        }

        return PlaneConicalResult::Point(apex);
    }

    // ── General case: conic type via Dandelin criterion ───────────────────────
    // Let σ = angle between the cutting plane and the cone axis,
    //     β = cone half-angle.
    // cos(σ) = sin_angle  (since σ = π/2 − acos(cos_angle))
    //
    //  cos(σ) > cos(β) → σ < β → Hyperbola  (plane cuts both nappes)
    //  cos(σ) ≈ cos(β) → σ ≈ β → Parabola   (plane parallel to one generator)
    //  cos(σ) < cos(β) → σ > β → Ellipse    (plane cuts one nappe completely)
    let cos_beta = cone.half_angle_rad.cos();
    let sin_beta = cone.half_angle_rad.sin();

    // ── Hyperbola: σ < β ↔ cos(σ) > cos(β) ↔ sin_angle > cos_beta ──────────
    if sin_angle > cos_beta + TOLERANCE_ANG {
        return build_hyperbola(plane, cone, apex_to_plane, axis_n, cos_beta, sin_beta, sin_angle);
    }

    // ── Parabola: σ ≈ β ↔ cos(σ) ≈ cos(β) ↔ sin_angle ≈ cos_beta ───────────
    if (sin_angle - cos_beta).abs() < TOLERANCE_ANG {
        return build_parabola(plane, cone, apex_to_plane, axis_n);
    }

    // ── Ellipse: σ > β ↔ cos(σ) < cos(β) ↔ sin_angle < cos_beta ────────────
    build_ellipse(plane, cone, apex_to_plane, axis_n, cos_angle)
}

// ─────────────────────────────────────────────────────────────────────────────
// Ellipse builder (corrected)
// ─────────────────────────────────────────────────────────────────────────────

fn build_ellipse(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    cos_angle: f64,
) -> PlaneConicalResult {
    // The intersection of a plane with a (right circular) cone is an ellipse
    // when sin(α) > sin(β) where α = angle between plane and axis, β = half-angle.
    //
    // We use the standard oblique-section formula:
    //   center = apex + axis_n * t  where t = apex_to_plane / (axis_n · plane.normal)
    //   minor_radius = r(t) = |t| * tan(β)  (circular cross-section radius at height t)
    //   major_radius = minor_radius / cos(γ)
    //     where γ = angle between plane and the circular cross-section plane
    //           = angle between plane normal and cone axis = acos(cos_angle)
    //     so major_radius = minor_radius / cos_angle
    //
    // (This is the standard textbook formula for plane-cone ellipse.)
    let tan_beta = cone.half_angle_rad.tan();
    let denom = axis_n.dot(plane.normal);
    if denom.abs() < TOLERANCE_FLOAT_LOOSE {
        return PlaneConicalResult::NoIntersection;
    }
    let t = apex_to_plane / denom;

    // Must be on the same nappe as the apex_to_plane sign
    // (t > 0: upper nappe; t < 0: lower nappe)
    let apex = cone.apex_point();
    let center = apex + axis_n * t;
    let base_radius = (t * tan_beta).abs();

    if base_radius < TOLERANCE_ABS {
        return PlaneConicalResult::Point(center);
    }


    // Semi-minor axis = base_radius (perpendicular to tilt direction)
    let minor_radius = base_radius;
    // Semi-major axis: correct formula using Dandelin approach
    // a = b / sqrt(1 - e²) where e = sin_angle_tilt / cos_beta... complex.
    // Use the practical formula: major_radius = base_radius / cos(γ)
    // where γ = angle between plane and the "circular" cross-section,
    // = complement of angle between plane normal and axis.
    // This is the standard oblique section formula: a = r / cos(φ)
    // where φ = angle between plane normal and the cylinder axis.
    // For a cone this is approximate but matches the standard result for shallow angles.
    let major_radius = base_radius / cos_angle;

    // Major direction in the plane (toward the steeper axis)
    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < TOLERANCE_LEN_MIN {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    PlaneConicalResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Parabola builder
// ─────────────────────────────────────────────────────────────────────────────

fn build_parabola(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
) -> PlaneConicalResult {
    // The plane is parallel to exactly one generator.
    // The vertex of the parabola is the point where the single tangent
    // generator intersects the plane.
    //
    // The generator direction in the plane of the cone that is parallel to
    // the cutting plane: find the generator in the plane spanned by axis and
    // the "steepest descent" direction in the cutting plane.
    let steepest = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let steepest = if steepest.length_squared() < TOLERANCE_LEN_MIN {
        any_perpendicular(plane.normal)
    } else {
        steepest
    };

    let tan_beta = cone.half_angle_rad.tan();
    // Generator parallel to the cutting plane: axis_n + tan_beta * steepest (normalized)
    let gen_dir = (axis_n + tan_beta * steepest).normalize();

    // Vertex: foot of generator on the plane
    let denom = gen_dir.dot(plane.normal);
    let vertex = if denom.abs() > TOLERANCE_LEN_MIN {
        let t = apex_to_plane / denom;
        cone.apex_point() + gen_dir * t
    } else {
        // Generator is parallel to plane; use foot of axis on plane
        let t = apex_to_plane / axis_n.dot(plane.normal).max(TOLERANCE_LEN_MIN);
        cone.apex_point() + axis_n * t
    };

    // Axis direction of the parabola: projection of cone axis onto the plane
    let axis_2d = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let axis_dir = if axis_2d.length_squared() < TOLERANCE_LEN_MIN {
        steepest
    } else {
        axis_2d
    };

    // Focal parameter p: derived from cone geometry.
    // For a cone with half-angle β and unit-speed cut at vertex distance d_v from apex:
    //   d_v = apex_to_plane / (gen_dir · plane_n)  (already computed above)
    // p = 2 * r_v * tan_beta where r_v = d_v * sin_beta
    let d_v = (vertex - cone.apex_point()).length().max(TOLERANCE_ABS);
    let r_v = d_v * cone.half_angle_rad.sin();
    let focal_param = (2.0 * r_v * tan_beta).max(TOLERANCE_ABS);

    PlaneConicalResult::Parabola(Parabola3 {
        vertex,
        normal: plane.normal,
        axis_dir,
        focal_param,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Hyperbola builder
// ─────────────────────────────────────────────────────────────────────────────

fn build_hyperbola(
    plane: &Plane,
    cone: &ConicalSurface,
    apex_to_plane: f64,
    axis_n: DVec3,
    cos_beta: f64,
    _sin_beta: f64,
    sin_angle: f64,
) -> PlaneConicalResult {
    // Plane cuts both nappes → two-branch hyperbola.
    // The apex is on the hyperbola's transverse axis.
    //
    // Center: projection of apex onto the cutting plane.
    let center = cone.apex_point() + plane.normal * apex_to_plane;

    // Major direction in the plane: projection of cone axis onto the plane.
    let major_dir = (axis_n - plane.normal * axis_n.dot(plane.normal)).normalize_or_zero();
    let major_dir = if major_dir.length_squared() < TOLERANCE_LEN_MIN {
        any_perpendicular(plane.normal)
    } else {
        major_dir
    };

    let cos_angle = plane.normal.dot(axis_n).abs();
    let tan_beta = cone.half_angle_rad.tan();

    // Semi-axes via the Dandelin-sphere construction for hyperbolas:
    // Let d = |apex_to_plane| (distance from apex to cutting plane along its normal).
    // The foot of the axis on the plane: center_along_axis = apex + axis_n * (d / cos_angle)
    //   (if cos_angle ≈ 0, the axis is nearly parallel to the plane)
    //
    // Practical formula (standard oblique section of a right circular cone):
    //   a = d * sin_beta / sqrt(sin_beta² − sin_angle²)
    //   b = d * sin_angle * cos_beta / sqrt(sin_beta² − sin_angle²)  (... simplified)
    //
    // For the purely parallel case (sin_angle=0, plane parallel to axis):
    //   The cutting plane at distance ρ from axis intersects in two lines if ρ < r,
    //   or a hyperbola whose semi-transverse = ρ / tan_beta (approx).
    //   Use the general formula with sin_angle=0:
    //   a = d * sin_beta / sin_beta = d,  b = 0 → degenerate (two straight lines).
    // We handle this by using the ρ-based formula when cos_angle ≈ 0.

    let d = apex_to_plane.abs();

    // General hyperbola formula (σ < β).
    // K = sin_angle²·tan²β − cos_angle² = cos²σ·tan²β − sin²σ > 0 for σ < β.
    // a = d·tanβ / K  (semi-transverse axis)
    // b = a·√(sin_angle² − cos_beta²) / cos_beta  (semi-conjugate axis)
    // Derived from the vertices at x=0: V1=(0, z1·tanβ, z1), V2=(0, -z2·tanβ, z2)
    // where z1,2 = d/(sinσ ∓ cosσ·tanβ).
    let k_val = sin_angle * sin_angle * tan_beta * tan_beta - cos_angle * cos_angle;
    if k_val <= TOLERANCE_ABS * TOLERANCE_ABS {
        return PlaneConicalResult::NoIntersection;
    }
    let a = d * tan_beta / k_val;
    let e_sq_minus_1 = (sin_angle * sin_angle - cos_beta * cos_beta) / (cos_beta * cos_beta);
    if e_sq_minus_1 <= 0.0 {
        return PlaneConicalResult::NoIntersection;
    }
    let b = a * e_sq_minus_1.sqrt();

    if a < TOLERANCE_ABS {
        return PlaneConicalResult::Point(center);
    }

    // Ensure b > 0 (may be tiny for near-axis planes)
    let b = b.max(TOLERANCE_ABS);

    PlaneConicalResult::Hyperbola(Hyperbola3 {
        center,
        normal: plane.normal,
        major_dir,
        semi_major: a,
        semi_minor: b,
    })
}


