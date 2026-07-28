use glam::DVec3;
use rcad_kernel::geom::*;



#[derive(Debug, Clone)]
pub enum PlaneSphereResult {
    NoIntersection,
    /// Legacy variant: [`intersect_plane_sphere`] now returns a tiny [`Circle3`] instead so Pave
    /// always records a trim curve for grazing planes (OCCT `bcommon_simple/A4`).
    TangentPoint(DVec3),
    Circle(Circle3),
}

pub fn intersect_plane_sphere(plane: &Plane, sphere: &SphericalSurface) -> PlaneSphereResult {
    let signed_dist = (sphere.center - plane.origin).dot(plane.normal);
    let abs_dist = signed_dist.abs();

    if abs_dist > sphere.radius + rcad_kernel::rcad_kernel::CONFUSION {
        return PlaneSphereResult::NoIntersection;
    }

    let r_sq = (sphere.radius * sphere.radius - signed_dist * signed_dist).max(0.0);
    let mut circle_radius = r_sq.sqrt();
    // Former `TangentPoint` branch: only inflate when the plane is numerically tangent
    // (|d|鈮圧), not for every legitimately small intersection circle 鈥?inflating all small radii
    // can destabilize other booleans (e.g. sphere鈥揷ylinder in `occt_alignment`).
    let tangent_band = rcad_kernel::rcad_kernel::CONFUSION.max(TOLERANCE_COORD_SUB * sphere.radius.abs());
    let tangent_like = (abs_dist - sphere.radius).abs() < tangent_band;
    const MIN_R_TANGENT: f64 = rcad_kernel::APPROXIMATION;
    if tangent_like && circle_radius < MIN_R_TANGENT {
        circle_radius = MIN_R_TANGENT;
    }

    let center = sphere.center - plane.normal * signed_dist;

    PlaneSphereResult::Circle(Circle3::new(center, plane.normal, circle_radius))
}
