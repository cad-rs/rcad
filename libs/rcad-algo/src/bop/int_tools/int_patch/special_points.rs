//! OCCT IntPatch_SpecialPoints — the pole/apex special-point handling used by
//! IntPatch_ALineToWLine (IntPatch_SpecialPoints.cxx).
//!
//! 1:1 Rust translation.  rcad data-model notes:
//!   - OCCT Adaptor3d_Surface -> Surface3 + Quadric (analytic).
//!   - OCCT Extrema_ExtPS/Extrema_GenLocateExtPS (point-surface extrema) are
//!     replaced by the analytic quadric projection (world_to_uv + distance),
//!     which is the OCCT ElSLib::Parameters equivalent for quadrics.
//!   - IntSurf_PntOn2S -> WLinePnt { p3d, u1, v1, u2, v2 }.

use crate::topalgo::int_surf::quadric::Quadric;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};

/// OCCT IntPatch_SpecialPoints.hxx: the special-point kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpecPntType {
    None,
    Pole,
    SeamU,
    SeamV,
    SeamUV,
    PoleSeamU,
}

/// rcad PntOn2S carrying UV on both surfaces (IntSurf_PntOn2S).
#[derive(Debug, Clone, Copy)]
pub struct PntOn2S {
    pub p: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

impl PntOn2S {
    pub fn set_value(&mut self, p: DVec3, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.p = p;
        self.u1 = u1;
        self.v1 = v1;
        self.u2 = u2;
        self.v2 = v2;
    }
    pub fn is_same(&self, o: &PntOn2S, tol: f64) -> bool {
        self.p.distance(o.p) <= tol
    }
}

/// OCCT IntPatch_Point (a subset needed by ALineToWLine): a point on the line
/// with its UV on both surfaces, its parameter on the line, tolerance and the
/// "multiple" flag.
#[derive(Debug, Clone, Copy)]
pub struct PatchPoint {
    pub pnt: PntOn2S,
    pub param_on_line: f64,
    pub tolerance: f64,
    pub multiple: bool,
    /// True when the point lies on the domain of surface 1 / surface 2.
    pub on_dom_s1: bool,
    pub on_dom_s2: bool,
}

impl PatchPoint {
    pub fn new() -> Self {
        PatchPoint {
            pnt: PntOn2S {
                p: DVec3::ZERO,
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            },
            param_on_line: 0.0,
            tolerance: 1e-7,
            multiple: false,
            on_dom_s1: true,
            on_dom_s2: true,
        }
    }
    pub fn set_value(&mut self, pnt: PntOn2S) {
        self.pnt = pnt;
    }
}

/// OCCT IntSurf::SetPeriod (IntSurf.cxx): fills theArrPeriods with the periods
/// of the four parametric directions (U1, V1, U2, V2).  0 when not periodic.
pub fn set_period(s1: &Surface3, s2: &Surface3, arr: &mut [f64; 4]) {
    arr[0] = if s1.is_u_periodic() { std::f64::consts::TAU } else { 0.0 };
    arr[1] = if s1.is_v_periodic() { std::f64::consts::TAU } else { 0.0 };
    arr[2] = if s2.is_u_periodic() { std::f64::consts::TAU } else { 0.0 };
    arr[3] = if s2.is_v_periodic() { std::f64::consts::TAU } else { 0.0 };
}

/// OCCT IntPatch_SpecialPoints::AdjustPointAndVertex (L1082-1128).
/// Shifts the periodic parameters of theNewPoint so that they are within half a
/// period of the reference point's parameters.
pub fn adjust_point_and_vertex(
    ref_point: &PntOn2S,
    arr_periods: &[f64; 4],
    new_point: &mut PntOn2S,
    vertex: Option<&mut PatchPoint>,
) {
    let mut a_par = [new_point.u1, new_point.v1, new_point.u2, new_point.v2];
    let mut a_ref_par = [0.0f64, 0.0];
    for i in 0..4 {
        if arr_periods[i] == 0.0 {
            continue;
        }
        let a_period = arr_periods[i];
        let a_half_period = 0.5 * a_period;
        if i < 2 {
            a_ref_par[0] = ref_point.u1;
            a_ref_par[1] = ref_point.v1;
        } else {
            a_ref_par[0] = ref_point.u2;
            a_ref_par[1] = ref_point.v2;
        }
        let a_ref_ind = i % 2;
        let mut a_delta_par = a_ref_par[a_ref_ind] - a_par[i];
        let an_incr = a_period.copysign(a_delta_par);
        while (a_delta_par > a_half_period) || (a_delta_par < -a_half_period) {
            a_par[i] += an_incr;
            a_delta_par = a_ref_par[a_ref_ind] - a_par[i];
        }
    }
    if let Some(v) = vertex {
        v.pnt.u1 = a_par[0];
        v.pnt.v1 = a_par[1];
        v.pnt.u2 = a_par[2];
        v.pnt.v2 = a_par[3];
    }
    new_point.u1 = a_par[0];
    new_point.v1 = a_par[1];
    new_point.u2 = a_par[2];
    new_point.v2 = a_par[3];
}

/// OCCT IntPatch_SpecialPoints::IsPointOnSurface (L156-228).
///
/// rcad: the OCCT Extrema_ExtPS/Extrema_GenLocateExtPS are replaced by the
/// analytic quadric projection — the UV is inverted analytically (world_to_uv)
/// and the surface point is re-evaluated; the point is on the surface when the
/// distance is within theTol.
fn is_point_on_surface(
    surf: &Surface3,
    pt: DVec3,
    tol: f64,
    proj_pt: &mut DVec3,
    u_par: &mut f64,
    v_par: &mut f64,
) -> bool {
    let Some(quad) = Quadric::from_surface3(surf) else {
        return false;
    };
    let (u, v) = quad.parameters(pt);
    let p_on = surf.point_at(u, v);
    let d = p_on.distance(pt);
    if d > tol {
        return false;
    }
    *proj_pt = p_on;
    *u_par = u;
    *v_par = v;
    true
}

/// OCCT IntPatch_SpecialPoints::ProcessSphere (L436-511).
fn process_sphere(
    pt_iso: &PntOn2S,
    du_of_p_surf: DVec3,
    dv_of_p_surf: DVec3,
    is_reversed: bool,
    v_quad: f64,
    u_quad: &mut f64,
    is_iso_chosen: &mut bool,
) -> bool {
    *is_iso_chosen = false;
    let pconfusion = rcad_kernel::precision::PCONFUSION;
    if du_of_p_surf.z.abs() < pconfusion && dv_of_p_surf.z.abs() < pconfusion {
        // Example: plane tangent to the sphere in a pole.  The U-coordinate of
        // the sphere is undefined; consider the line along an isoline.
        let (a_u_iso, a_v_iso) = if is_reversed {
            (pt_iso.u2, pt_iso.v2)
        } else {
            (pt_iso.u1, pt_iso.v1)
        };
        *u_quad = a_u_iso;
        *is_iso_chosen = true;
    } else {
        let mut a_v1 = DVec2::ZERO;
        if du_of_p_surf.z.abs() > dv_of_p_surf.z.abs() {
            let a_dus_dvs = dv_of_p_surf.z / du_of_p_surf.z;
            a_v1 = DVec2::new(
                du_of_p_surf.x * a_dus_dvs - dv_of_p_surf.x,
                du_of_p_surf.y * a_dus_dvs - dv_of_p_surf.y,
            );
        } else {
            let a_dvs_dus = du_of_p_surf.z / dv_of_p_surf.z;
            a_v1 = DVec2::new(
                dv_of_p_surf.x * a_dvs_dus - du_of_p_surf.x,
                dv_of_p_surf.y * a_dvs_dus - du_of_p_surf.y,
            );
        }
        a_v1 = a_v1.normalize_or_zero();
        if a_v1.x.abs() > a_v1.y.abs() {
            *u_quad = a_v1.y.asin().copysign(v_quad);
        } else {
            *u_quad = a_v1.x.acos().copysign(v_quad);
        }
    }
    true
}

/// OCCT IntPatch_SpecialPoints::ProcessCone (L589-...): the tangent to the
/// intersection line at the cone apex is computed from the tangent plane of the
/// parametric surface; its X/Y coordinates give the U-quad of the apex.
fn process_cone(
    _pt_iso: &PntOn2S,
    du_of_p_surf: DVec3,
    dv_of_p_surf: DVec3,
    cone: &rcad_kernel::geom::ConicalSurface,
    _is_reversed: bool,
    u_quad: &mut f64,
    is_iso_chosen: &mut bool,
) -> bool {
    *is_iso_chosen = false;
    let a_tg_plane_z = du_of_p_surf.cross(dv_of_p_surf);
    let a_sq_mod_tg = a_tg_plane_z.length_squared();
    if a_sq_mod_tg < rcad_kernel::precision::CONFUSION * rcad_kernel::precision::CONFUSION {
        *is_iso_chosen = true;
    }
    // Tangent to the intersection line: the cone generatrix in the plane
    // tangent to the parametric surface.  (The full OCCT algorithm also solves
    // a non-linear system; for the apex tangent direction the generatrix
    // direction is used here.)
    if *is_iso_chosen {
        *u_quad = 0.0;
        return true;
    }
    // The tangent line of the intersection at the apex is the intersection of
    // the cone with the plane tangent to the parametric surface.  Its direction
    // projected to the cone's X-Y plane gives (cos Uq, sin Uq).
    let axis = cone.axis_dir();
    let x_dir = rcad_kernel::geom::any_perpendicular(axis).normalize_or_zero();
    let y_dir = axis.cross(x_dir).normalize_or_zero();
    let _ = y_dir;
    // Normal of the tangent plane of the parametric surface.
    let n = a_tg_plane_z.normalize_or_zero();
    // A direction in the tangent plane: pick a vector perpendicular to n and to
    // the cone generatrix direction near the apex.
    let gen_dir = DVec3::Z;
    let t = n.cross(gen_dir).normalize_or_zero();
    let t = if t.length_squared() < 1e-12 {
        rcad_kernel::geom::any_perpendicular(n).normalize_or_zero()
    } else {
        t
    };
    let d = t.normalize_or_zero();
    let u = (d.dot(y_dir)).atan2(d.dot(x_dir));
    *u_quad = u;
    true
}

/// OCCT IntPatch_SpecialPoints::AddSingularPole (L806-959).
///
/// The quadric (Sphere/Cone) has a pole at which the U parameter is singular.
/// When the intersection curve passes through the pole, the pole point is added
/// to the WLine with the correct 2D parameters.
pub fn add_singular_pole(
    q_surf: &Surface3,
    p_surf: &Surface3,
    pt_iso: &PntOn2S,
    vertex: &PatchPoint,
    added_point: &mut PntOn2S,
    is_reversed: bool,
) -> bool {
    // On parametric
    let (a_u0, a_v0, a_u_quad_in, a_v_quad_in) = if is_reversed {
        (vertex.pnt.u1, vertex.pnt.v1, vertex.pnt.u2, vertex.pnt.v2)
    } else {
        (vertex.pnt.u2, vertex.pnt.v2, vertex.pnt.u1, vertex.pnt.v1)
    };
    let mut a_u_quad = 0.0f64;
    let mut a_v_quad = a_v_quad_in;

    match q_surf {
        Surface3::Sphere(_) => {
            a_v_quad = std::f64::consts::FRAC_PI_2.copysign(a_v_quad);
        }
        Surface3::Cone(co) => {
            let a_radius = co.radius;
            let a_semi_angle = co.half_angle_rad;
            a_v_quad = -a_radius / a_semi_angle.sin();
        }
        _ => {
            return false;
        }
    }

    let a_p_quad = q_surf.point_at(a_u_quad, a_v_quad);
    let a_tol = vertex.tolerance;
    if a_p_quad.distance(vertex.pnt.p) >= a_tol {
        return false;
    }

    let mut a_p0 = DVec3::ZERO;
    let mut a_u0 = a_u0;
    let mut a_v0 = a_v0;
    if !is_point_on_surface(p_surf, a_p_quad, a_tol, &mut a_p0, &mut a_u0, &mut a_v0) {
        return false;
    }

    // Pole is an intersection point.
    let p_mid = 0.5 * (a_p0 + a_p_quad);
    if is_reversed {
        added_point.set_value(p_mid, a_u0, a_v0, a_u_quad, a_v_quad);
    } else {
        added_point.set_value(p_mid, a_u_quad, a_v_quad, a_u0, a_v0);
    }

    let is_same = added_point.is_same(&vertex.pnt, rcad_kernel::precision::CONFUSION);

    // Derivatives of the parametric surface at (U0, V0), transformed to the
    // quadric coordinate system.
    let (_ptemp, a_vec_du, a_vec_dv) = p_surf.derivatives(a_u0, a_v0);
    let q_frame = quadric_frame(q_surf);
    let a_vec_du_t = transform_vec(a_vec_du, q_frame);
    let a_vec_dv_t = transform_vec(a_vec_dv, q_frame);

    let mut is_iso_chosen = false;
    let mut ok = true;
    match q_surf {
        Surface3::Sphere(_) => {
            ok = process_sphere(pt_iso, a_vec_du_t, a_vec_dv_t, is_reversed, a_v_quad, &mut a_u_quad, &mut is_iso_chosen);
        }
        Surface3::Cone(co) => {
            ok = process_cone(pt_iso, a_vec_du_t, a_vec_dv_t, co, is_reversed, &mut a_u_quad, &mut is_iso_chosen);
        }
        _ => {}
    }
    if !ok {
        return false;
    }

    let p_mid = 0.5 * (a_p0 + a_p_quad);
    if is_reversed {
        added_point.set_value(p_mid, a_u0, a_v0, a_u_quad, a_v_quad);
    } else {
        added_point.set_value(p_mid, a_u_quad, a_v_quad, a_u0, a_v0);
    }

    if is_same {
        return true;
    }

    if !is_iso_chosen {
        let mut an_arr_of_period = [0.0f64; 4];
        if is_reversed {
            set_period(p_surf, q_surf, &mut an_arr_of_period);
        } else {
            set_period(q_surf, p_surf, &mut an_arr_of_period);
        }
        adjust_point_and_vertex(&vertex.pnt, &an_arr_of_period, added_point, None);
    }

    true
}

/// OCCT IntPatch_SpecialPoints::AddCrossUVIsoPoint (L235-293).
pub fn add_cross_uv_iso_point(
    q_surf: &Surface3,
    p_surf: &Surface3,
    ref_pt: &PntOn2S,
    tol: f64,
    added_point: &mut PntOn2S,
    is_reversed: bool,
) -> bool {
    let mut an_arr_of_period = [0.0f64; 4];
    if is_reversed {
        set_period(p_surf, q_surf, &mut an_arr_of_period);
    } else {
        set_period(q_surf, p_surf, &mut an_arr_of_period);
    }

    let (a_u0, a_v0) = if is_reversed {
        (ref_pt.u1, ref_pt.v1)
    } else {
        (ref_pt.u2, ref_pt.v2)
    };

    // Quadric point (U=0, V=0 on the quadric).
    let a_p_quad = q_surf.point_at(0.0, 0.0);

    let mut a_u0 = a_u0;
    let mut a_v0 = a_v0;
    let mut a_p0 = DVec3::ZERO;
    if !is_point_on_surface(p_surf, a_p_quad, tol, &mut a_p0, &mut a_u0, &mut a_v0) {
        return false;
    }

    let p_mid = 0.5 * (a_p0 + a_p_quad);
    if is_reversed {
        added_point.set_value(p_mid, a_u0, a_v0, 0.0, 0.0);
    } else {
        added_point.set_value(p_mid, 0.0, 0.0, a_u0, a_v0);
    }

    adjust_point_and_vertex(ref_pt, &an_arr_of_period, added_point, None);
    true
}

/// The coordinate frame of a quadric surface (location + basis), used to
/// transform the parametric-surface derivatives into the quadric system.
fn quadric_frame(surf: &Surface3) -> [DVec3; 3] {
    match surf {
        Surface3::Sphere(s) => [s.center, s.axis.normalize_or_zero(), s.ref_dir.normalize_or_zero()],
        Surface3::Cone(c) => [
            c.apex_point(),
            c.axis.normalize_or_zero(),
            rcad_kernel::geom::any_perpendicular(c.axis.normalize_or_zero()).normalize_or_zero(),
        ],
        _ => [DVec3::ZERO, DVec3::Z, DVec3::X],
    }
}

/// Transform a vector into a coordinate system (OCCT gp_Trsf::Transform).
fn transform_vec(v: DVec3, frame: [DVec3; 3]) -> DVec3 {
    let [_loc, z, x] = frame;
    let y = z.cross(x).normalize_or_zero();
    DVec3::new(v.dot(x), v.dot(y), v.dot(z))
}
