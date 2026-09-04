//! IntCurveSurface_InterUtils.pxx (L1-1637) — the shared helpers of the
//! HInter assembly: UV interval decomposition, sampling (DoSurface /
//! DoNewBounds), parameter clamping, transition classification, the
//! interference start-point pipeline (Collect / Sort / Process) and the
//! infinite-surface bound estimators (EstLim*).
//!
//! SectionPointToParameters (InterUtils.pxx L740-848) already lives in
//! [`crate::geomalgo::intf::section_point_to_parameters`]; the polyhedron /
//! polygon sides implement [`crate::geomalgo::intf::PolyhedronLike`] /
//! [`crate::geomalgo::intf::PolygonLike`].

use glam::DVec3;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Line2d, Line3, Plane};
use rcad_kernel::math::direct_polynomial_roots::DirectPolynomialRoots;
use rcad_kernel::math::el::{
    elslib_cone_parameters, elslib_cylinder_parameters, elslib_plane_parameters,
    elslib_sphere_parameters, in_period,
};
use rcad_kernel::math::function_set_root::FunctionSetRoot;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::precision::PCONFUSION;

use crate::geomalgo::int_imp::int_cs::IntCS;
use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};
use crate::geomalgo::int_patch::int_conic_quad::IntConicQuad;
use crate::geomalgo::int_patch::int_lin_torus::IntLinTorus;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};
use crate::geomalgo::intf::{section_point_to_parameters, PolyhedronLike, PolygonLike};

use super::{Adaptor3dCurveBasis, Adaptor3dSurfaceBasis, HCurveTool, HSurfaceTool, IntersectionPoint, SurfaceType, TransitionOnCurve};

/// OCCT IntCurveSurface_InterUtils THE_TOLERANCE_ANGULAIRE (pxx L55).
pub(crate) const THE_TOLERANCE_ANGULAIRE: f64 = 1.0e-12;
/// OCCT IntCurveSurface_InterUtils THE_TOLTANGENCY (pxx L56).
pub(crate) const THE_TOLTANGENCY: f64 = 0.00000001;

/// OCCT IntCurveSurface_InterUtils::ProjectIntersectAndEstLim (pxx L61-134)
/// — project theLine and its X-axis symmetric line to thePln and intersect
/// the resulting curves with theBasCurvProj projection; accumulate the
/// extreme parameters on the basis curve projection.
pub(crate) fn project_intersect_and_est_lim(
    the_line: &Line3,
    the_pln: &Plane,
    the_bas_curv_proj: &mut rcad_kernel::base::proj_lib::PlaneProjector,
    the_vmin: &mut f64,
    the_vmax: &mut f64,
    the_no_intersection: &mut bool,
) {
    let mut a_line_proj = rcad_kernel::base::proj_lib::PlaneProjector::with_plane(the_pln);
    a_line_proj.project_line(the_line);
    if !a_line_proj.projector().is_done() {
        return;
    }
    let proj = a_line_proj.projector();
    let lin = proj.line();
    let a_lin2d = Line2d::new(
        glam::DVec2::new(lin.origin.x, lin.origin.y),
        glam::DVec2::new(lin.direction.x, lin.direction.y),
    );

    // Make a second line X-axe symmetric to the first one.
    let a_p1 = a_lin2d.origin;
    let a_p2 = a_lin2d.origin + a_lin2d.direction;
    let a_p1sym = glam::DVec2::new(a_p1.x, -a_p1.y);
    let a_p2sym = glam::DVec2::new(a_p2.x, -a_p2.y);
    let a_lin2dsym = Line2d::new(a_p1sym, a_p2sym - a_p1sym);

    // Intersect projections.
    use rcad_kernel::base::int_ana2d::Conic2d;
    let a_con = Conic2d::from_line(&a_lin2d);
    let a_con_sym = Conic2d::from_line(&a_lin2dsym);
    use rcad_kernel::base::int_ana2d::AnaIntersection2d;
    let mut an_intersect = AnaIntersection2d::new();
    let mut an_intersect_sym = AnaIntersection2d::new();

    let bas_proj = the_bas_curv_proj.projector();
    match bas_proj.get_type() {
        CurveType::Line => {
            let l = lin2d_of(bas_proj.line());
            an_intersect_sym.perform_lin_conic(&l, &a_con_sym);
            an_intersect.perform_lin_conic(&l, &a_con);
        }
        CurveType::Hyperbola => {
            let h = hypr2d_of(bas_proj.hyperbola());
            an_intersect_sym.perform_hyperbola_conic(&h, &a_con_sym);
            an_intersect.perform_hyperbola_conic(&h, &a_con);
        }
        CurveType::Parabola => {
            let p = parab2d_of(bas_proj.parabola());
            an_intersect_sym.perform_parabola_conic(&p, &a_con_sym);
            an_intersect.perform_parabola_conic(&p, &a_con);
        }
        _ => return, // not infinite curve
    }

    // Retrieve params of intersections.
    let a_nb_int_pnt = if an_intersect.is_done() { an_intersect.nb_points() } else { 0 };
    let a_nb_int_pnt_sym = if an_intersect_sym.is_done() { an_intersect_sym.nb_points() } else { 0 };
    let a_nb_pnt = a_nb_int_pnt.max(a_nb_int_pnt_sym);

    if a_nb_pnt == 0 {
        *the_no_intersection = true;
        return;
    }
    for i_pnt in 1..=a_nb_pnt {
        if i_pnt <= a_nb_int_pnt {
            let a_int_pnt = an_intersect.point(i_pnt);
            let a_param = a_int_pnt.param_on_first();
            *the_vmin = the_vmin.min(a_param);
            *the_vmax = the_vmax.max(a_param);
        }
        if i_pnt <= a_nb_int_pnt_sym {
            let a_int_pnt = an_intersect_sym.point(i_pnt);
            let a_param = a_int_pnt.param_on_first();
            *the_vmin = the_vmin.min(a_param);
            *the_vmax = the_vmax.max(a_param);
        }
    }
}

/// The projected Line3 (z = 0 in the plane frame) as a Line2d.
fn lin2d_of(l: &Line3) -> Line2d {
    Line2d::new(
        glam::DVec2::new(l.origin.x, l.origin.y),
        glam::DVec2::new(l.direction.x, l.direction.y),
    )
}

/// The projected Parabola3 (z = 0) as a Parabola2d.
fn parab2d_of(p: &rcad_kernel::geom::Parabola3) -> rcad_kernel::geom::Parabola2d {
    rcad_kernel::geom::Parabola2d {
        origin: glam::DVec2::new(p.vertex.x, p.vertex.y),
        axis_dir: glam::DVec2::new(p.axis_dir.x, p.axis_dir.y),
        focal_param: p.focal_param,
    }
}

/// The projected Hyperbola3 (z = 0) as a Hyperbola2d.
fn hypr2d_of(h: &rcad_kernel::geom::Hyperbola3) -> rcad_kernel::geom::Hyperbola2d {
    rcad_kernel::geom::Hyperbola2d {
        center: glam::DVec2::new(h.center.x, h.center.y),
        major_dir: glam::DVec2::new(h.major_dir.x, h.major_dir.y),
        semi_major: h.semi_major,
        semi_minor: h.semi_minor,
    }
}

/// OCCT gp_Dir::IsParallel(D, AngularTolerance) — the angle between the two
/// directions (or its complement to PI) is below the tolerance; near the
/// tolerance scale sin(angle) == angle.
pub(crate) fn dir_is_parallel(d1: DVec3, d2: DVec3, tolang: f64) -> bool {
    d1.normalize_or_zero().cross(d2.normalize_or_zero()).length() <= tolang
}

/// OCCT gp_Pln::Rotated(Ax1, angle) — rotate the plane position (origin and
/// normal) around the axis by the Rodrigues formula.
fn plane_rotated(pln: &Plane, axis_loc: DVec3, axis_dir: DVec3, angle: f64) -> Plane {
    let d = axis_dir.normalize_or_zero();
    let (s, c) = angle.sin_cos();
    let rot = |v: DVec3| -> DVec3 {
        v * c + d.cross(v) * s + d * (d.dot(v)) * (1.0 - c)
    };
    Plane {
        origin: axis_loc + rot(pln.origin - axis_loc),
        normal: rot(pln.normal).normalize_or_zero(),
        u_dir: rot(pln.u_dir).normalize_or_zero(),
        v_dir: rot(pln.v_dir).normalize_or_zero(),
    }
}

/// OCCT IntCurveSurface_InterUtils::EstLimForInfSurf (pxx L137-143) —
/// generic fallback bound estimation.
pub(crate) fn est_lim_for_inf_surf(u1new: &mut f64, u2new: &mut f64, v1new: &mut f64, v2new: &mut f64) {
    *u1new = u1new.max(-1.0e10);
    *u2new = u2new.min(1.0e10);
    *v1new = v1new.max(-1.0e10);
    *v2new = v2new.min(1.0e10);
}

/// OCCT IntCurveSurface_InterUtils::EstLimForInfExtr (pxx L148-371).
#[allow(clippy::too_many_arguments)]
pub(crate) fn est_lim_for_inf_extr<S, ST: HSurfaceTool<Surface = S>>(
    line: &Line3,
    surface: &S,
    is_off_surf: bool,
    nbsu: usize,
    u1inf: bool,
    u2inf: bool,
    v1inf: bool,
    v2inf: bool,
    u1new: &mut f64,
    u2new: &mut f64,
    v1new: &mut f64,
    v2new: &mut f64,
    no_intersection: &mut bool,
) {
    *no_intersection = false;

    let a_dir_of_ext = if is_off_surf {
        <ST as HSurfaceTool>::basis_surface(surface).direction()
    } else {
        <ST as HSurfaceTool>::direction(surface)
    };

    let tolang = THE_TOLERANCE_ANGULAIRE;

    if dir_is_parallel(a_dir_of_ext, line.direction, tolang) {
        *no_intersection = true;
        return;
    }

    if (v1inf || v2inf) && !(u1inf || u2inf) {
        let mut vmin = f64::MAX;
        let mut vmax = -f64::MAX;
        let step = (*u2new - *u1new) / nbsu as f64;
        let mut u = *u1new;

        use rcad_kernel::base::extrema::line_line_extrema;
        for _i in 0..=nbsu {
            let a_p = <ST as HSurfaceTool>::d0(surface, u, 0.0);
            let a_l = Line3 {
                origin: a_p,
                direction: a_dir_of_ext,
            };

            // OCCT Extrema_ExtElC aExtr(aL, Line, tolang): IsDone/IsParallel.
            let a_extr = line_line_extrema(&a_l, line);
            if a_extr.is_empty() {
                return;
            }

            // OCCT IsParallel: the two lines share the direction (checked
            // above for aDirOfExt; kept structural).
            if dir_is_parallel(a_dir_of_ext, line.direction, tolang) {
                *no_intersection = true;
                return;
            }

            // aExtr.Points(1, aP1, aP2); v = aP1.Parameter().
            let v = a_extr[0].1;
            vmin = vmin.min(v);
            vmax = vmax.max(v);

            u += step;
        }

        vmin = vmin - vmin.abs() - 10.0;
        vmax = vmax + vmax.abs() + 10.0;

        *v1new = v1new.max(vmin);
        *v2new = v2new.min(vmax);
    } else if u1inf || u2inf {
        let mut umin = f64::MAX;
        let mut umax = -f64::MAX;
        let u0 = 0.0f64.max(*u1new).min(*u2new);
        let v0 = 0.0f64.max(*v1new).min(*v2new);
        let a_p = <ST as HSurfaceTool>::d0(surface, u0, v0);
        let a_ref_pln = Plane::new(a_p, a_dir_of_ext);

        let a_bas_curv = if is_off_surf {
            <ST as HSurfaceTool>::basis_surface(surface).basis_curve()
        } else {
            <ST as HSurfaceTool>::basis_curve(surface)
        };

        let mut projector = rcad_kernel::base::proj_lib::PlaneProjector::with_plane(&a_ref_pln);
        projector.project_line(line);

        if !projector.projector().is_done() {
            return;
        }

        let proj = projector.projector();
        let line2d = lin2d_of(proj.line());

        let a_curv_typ = a_bas_curv.get_type();

        use rcad_kernel::base::int_ana2d::{AnaIntersection2d, Conic2d};
        if a_curv_typ == CurveType::Line {
            let mut projector2 = rcad_kernel::base::proj_lib::PlaneProjector::with_plane(&a_ref_pln);
            projector2.project_line(&a_bas_curv.line());

            if !projector2.projector().is_done() {
                return;
            }

            let a_l2d = lin2d_of(projector2.projector().line());

            let mut an_inter = AnaIntersection2d::new();
            an_inter.perform_lin_lin(&line2d, &a_l2d);

            if !an_inter.is_done() {
                return;
            }

            if an_inter.is_empty() || an_inter.identical_elements() || an_inter.parallel_elements() {
                *no_intersection = true;
                return;
            }

            let an_int_pnt = an_inter.point(1);
            umin = an_int_pnt.param_on_second();
            umax = umin;
        } else if a_curv_typ == CurveType::Parabola || a_curv_typ == CurveType::Hyperbola {
            let a_con = Conic2d::from_line(&line2d);
            let mut an_inter = AnaIntersection2d::new();

            if a_curv_typ == CurveType::Parabola {
                let mut projector2 =
                    rcad_kernel::base::proj_lib::PlaneProjector::with_plane(&a_ref_pln);
                projector2.project_parabola(&a_bas_curv.parabola());
                if !projector2.projector().is_done() {
                    return;
                }

                let a_p2d = parab2d_of(projector2.projector().parabola());

                an_inter.perform_parabola_conic(&a_p2d, &a_con);
            } else {
                let mut projector2 =
                    rcad_kernel::base::proj_lib::PlaneProjector::with_plane(&a_ref_pln);
                projector2.project_hyperbola(&a_bas_curv.hyperbola());
                if !projector2.projector().is_done() {
                    return;
                }

                let a_h2d = hypr2d_of(projector2.projector().hyperbola());
                an_inter.perform_hyperbola_conic(&a_h2d, &a_con);
            }

            if !an_inter.is_done() {
                return;
            }

            if an_inter.is_empty() {
                *no_intersection = true;
                return;
            }

            let nbint = an_inter.nb_points();
            for i in 1..=nbint {
                let an_int_pnt = an_inter.point(i);
                umin = umin.min(an_int_pnt.param_on_first());
                umax = umax.max(an_int_pnt.param_on_first());
            }
        } else {
            return;
        }

        umin = umin - umin.abs() - 10.0;
        umax = umax + umax.abs() + 10.0;

        *u1new = u1new.max(umin);
        *u2new = u2new.min(umax);

        if v1inf || v2inf {
            // OCCT recursion with U1inf = U2inf = false; U1new/U2new are
            // passed by reference and only read in that branch.
            let mut saved_u1 = *u1new;
            let mut saved_u2 = *u2new;
            est_lim_for_inf_extr::<S, ST>(
                line,
                surface,
                is_off_surf,
                nbsu,
                false,
                false,
                v1inf,
                v2inf,
                &mut saved_u1,
                &mut saved_u2,
                v1new,
                v2new,
                no_intersection,
            );
            *u1new = saved_u1;
            *u2new = saved_u2;
        }
    }
}

/// OCCT IntCurveSurface_InterUtils::EstLimForInfRevl (pxx L376-494).
#[allow(clippy::too_many_arguments)]
pub(crate) fn est_lim_for_inf_revl<S, ST: HSurfaceTool<Surface = S>>(
    line: &Line3,
    surface: &S,
    u1inf: bool,
    u2inf: bool,
    v1inf: bool,
    v2inf: bool,
    u1new: &mut f64,
    u2new: &mut f64,
    v1new: &mut f64,
    v2new: &mut f64,
    no_intersection: &mut bool,
) {
    *no_intersection = false;

    if u1inf || u2inf {
        if u1inf {
            *u1new = 0.0f64.max(*u1new);
        } else {
            *u2new = (2.0 * std::f64::consts::PI).min(*u2new);
        }
        if !v1inf && !v2inf {
            return;
        }
    }

    let a_basis_curve = <ST as HSurfaceTool>::basis_curve(surface);
    let (rev_ax_loc, rev_ax_dir) = <ST as HSurfaceTool>::axe_of_revolution(surface);
    let a_x_vec = rev_ax_dir;
    let a_tol_ang = rcad_kernel::precision::ANGULAR;

    // Make plane to project a basis curve.
    let o = rev_ax_loc;
    let mut a_u = 0.0f64;
    let mut p = a_basis_curve.value(a_u);
    while o.distance_squared(p) <= rcad_kernel::precision::PCONFUSION
        || dir_is_parallel(a_x_vec, p - o, a_tol_ang)
    {
        a_u += 1.0;
        p = a_basis_curve.value(a_u);
        if a_u > 3.0 {
            // Basis curve is a line coinciding with aXVec, P is any not on
            // aXVec.
            p = DVec3::new(a_u, a_u + 1.0, a_u + 2.0);
            break;
        }
    }
    let mut a_n_vec = a_x_vec.cross(p - o);
    let mut a_pln = Plane::new(o, a_n_vec);

    // Project basic curve.
    let mut a_bas_curv_proj = rcad_kernel::base::proj_lib::PlaneProjector::with_plane(&a_pln);
    match a_basis_curve.get_type() {
        CurveType::Line => a_bas_curv_proj.project_line(&a_basis_curve.line()),
        CurveType::Hyperbola => a_bas_curv_proj.project_hyperbola(&a_basis_curve.hyperbola()),
        CurveType::Parabola => a_bas_curv_proj.project_parabola(&a_basis_curve.parabola()),
        _ => return, // not infinite curve
    }
    if !a_bas_curv_proj.projector().is_done() {
        return;
    }
    // Make plane to project Line.
    if dir_is_parallel(a_x_vec, line.direction, a_tol_ang) {
        p = line.origin;
        while o.distance_squared(p) <= rcad_kernel::precision::PCONFUSION {
            a_u += 1.0;
            p = DVec3::new(a_u, a_u + 1.0, a_u + 2.0); // any not on aXVec
        }
        a_n_vec = a_x_vec.cross(p - o);
    } else {
        a_n_vec = a_x_vec.cross(line.direction);
    }

    a_pln = Plane::new(o, a_n_vec);

    // Make a second plane perpendicular to the first one, rotated around
    // aXVec.
    let a_pln_prp = plane_rotated(&a_pln, o, a_x_vec, std::f64::consts::PI / 2.0);

    // Project Line and its X-axe symmetric one to plane and intersect the
    // resulting curves with the projection of the Basic Curve.
    let mut a_vmin = f64::MAX;
    let mut a_vmax = -f64::MAX;
    let mut a_no_int1 = false;
    let mut a_no_int2 = false;
    project_intersect_and_est_lim(
        line,
        &a_pln,
        &mut a_bas_curv_proj,
        &mut a_vmin,
        &mut a_vmax,
        &mut a_no_int1,
    );
    project_intersect_and_est_lim(
        line,
        &a_pln_prp,
        &mut a_bas_curv_proj,
        &mut a_vmin,
        &mut a_vmax,
        &mut a_no_int2,
    );

    if a_no_int1 && a_no_int2 {
        *no_intersection = true;
        return;
    }

    a_vmin = a_vmin - a_vmin.abs() - 10.0;
    a_vmax = a_vmax + a_vmax.abs() + 10.0;

    if v1inf {
        *v1new = a_vmin;
    }
    if v2inf {
        *v2new = a_vmax;
    }
}

/// OCCT IntCurveSurface_InterUtils::EstLimForInfOffs (pxx L499-734).
#[allow(clippy::too_many_arguments)]
pub(crate) fn est_lim_for_inf_offs<S, ST: HSurfaceTool<Surface = S>>(
    line: &Line3,
    surface: &S,
    nbsu: usize,
    u1inf: bool,
    u2inf: bool,
    v1inf: bool,
    v2inf: bool,
    u1new: &mut f64,
    u2new: &mut f64,
    v1new: &mut f64,
    v2new: &mut f64,
    no_intersection: &mut bool,
) {
    use crate::geomalgo::int_surf::quadric::Quadric;
    *no_intersection = false;

    let a_bas_surf = <ST as HSurfaceTool>::basis_surface(surface);
    let an_off_val = <ST as HSurfaceTool>::offset_value(surface);

    let a_type_of_bas_surf = a_bas_surf.get_type();

    // Case for plane, cylinder and cone — make equivalent surface.
    if a_type_of_bas_surf == SurfaceType::Plane {
        let mut a_pln = a_bas_surf.plane();
        let a_t = a_pln.u_dir.cross(a_pln.v_dir) * an_off_val;
        a_pln.origin += a_t;
        let lin_plane = IntConicQuad::new_line_plane(line, &a_pln, THE_TOLERANCE_ANGULAIRE);

        if !lin_plane.is_done() {
            return;
        }

        if lin_plane.is_parallel() || lin_plane.is_in_quadric() {
            *no_intersection = true;
            return;
        }

        let (u, v) = elslib_plane_parameters(lin_plane.point(1), a_pln.origin, a_pln.u_dir, a_pln.v_dir);
        *u1new = u1new.max(u - 10.0);
        *u2new = u2new.min(u + 10.0);
        *v1new = v1new.max(v - 10.0);
        *v2new = v2new.min(v + 10.0);
    } else if a_type_of_bas_surf == SurfaceType::Cylinder {
        let mut a_cyl = a_bas_surf.cylinder();

        let mut a_r = a_cyl.radius;
        // OCCT anA.Direct() — the right-handed frame offsets outward.
        let direct = a_cyl.y_axis().dot(a_cyl.axis.cross(a_cyl.ref_dir).normalize_or_zero()) >= 0.0;
        if direct {
            a_r += an_off_val;
        } else {
            a_r -= an_off_val;
        }

        if a_r >= THE_TOLTANGENCY {
            a_cyl.radius = a_r;
        } else if a_r <= -THE_TOLTANGENCY {
            // OCCT anA.Rotate(gp_Ax1(anA.Location(), anA.Direction()), PI):
            // the frame X/Y directions turn by PI around the axis, the axis
            // itself is invariant.
            a_cyl.ref_dir = -a_cyl.ref_dir;
            a_cyl.y_dir = a_cyl.y_dir.map(|y| -y);
            a_cyl.radius = -a_r;
        } else {
            *no_intersection = true;
            return;
        }

        let lin_cylinder =
            IntConicQuad::new_line_quadric(line, &Quadric::from_cylinder(&a_cyl));

        if !lin_cylinder.is_done() {
            return;
        }

        if lin_cylinder.is_parallel() || lin_cylinder.is_in_quadric() {
            *no_intersection = true;
            return;
        }

        let nbp = lin_cylinder.nb_points();
        let mut vmin = f64::MAX;
        let mut vmax = -f64::MAX;

        for i in 1..=nbp {
            let y_ax = a_cyl.y_axis();
            let (u, v) = elslib_cylinder_parameters(
                lin_cylinder.point(i),
                a_cyl.origin,
                a_cyl.ref_dir,
                y_ax,
                a_cyl.axis,
                a_cyl.radius,
            );
            vmin = vmin.min(v);
            vmax = vmax.max(v);
        }

        *v1new = v1new.max(vmin - vmin.abs() - 10.0);
        *v2new = v2new.min(vmax + vmax.abs() + 10.0);
    } else if a_type_of_bas_surf == SurfaceType::Cone {
        let mut a_con = a_bas_surf.cone();
        let an_ang = a_con.half_angle_rad;
        let a_r = a_con.radius + an_off_val * an_ang.cos();
        if a_r >= 0.0 {
            // OCCT anA.Translate(aZ) with aZ = -anOffVal·sin(ang)·Z.
            let a_z = a_con.axis * (-an_off_val * an_ang.sin());
            a_con.apex += a_z;
            a_con.radius = a_r;
            a_con.half_angle_rad = an_ang;
        } else {
            return;
        }

        let lin_cone = IntConicQuad::new_line_quadric(line, &Quadric::from_cone(&a_con));

        if !lin_cone.is_done() {
            return;
        }

        if lin_cone.is_parallel() || lin_cone.is_in_quadric() {
            *no_intersection = true;
            return;
        }

        let nbp = lin_cone.nb_points();
        let mut vmin = f64::MAX;
        let mut vmax = -f64::MAX;

        for i in 1..=nbp {
            let y_ax = a_con.axis.cross(a_con.ref_dir).normalize_or_zero();
            let (u, v) = elslib_cone_parameters(
                lin_cone.point(i),
                a_con.apex,
                a_con.ref_dir,
                y_ax,
                a_con.axis,
                a_con.radius,
                a_con.half_angle_rad,
            );
            vmin = vmin.min(v);
            vmax = vmax.max(v);
        }

        *v1new = v1new.max(vmin - vmin.abs() - 10.0);
        *v2new = v2new.min(vmax + vmax.abs() + 10.0);
    } else if a_type_of_bas_surf == SurfaceType::SurfaceOfExtrusion {
        let an_u1 = *u1new;
        let an_u2 = *u2new;

        est_lim_for_inf_extr::<S, ST>(
            line,
            surface,
            true,
            nbsu,
            u1inf,
            u2inf,
            v1inf,
            v2inf,
            u1new,
            u2new,
            v1new,
            v2new,
            no_intersection,
        );

        if *no_intersection {
            return;
        }

        if u1inf || u2inf {
            let a_bas_curv_type = a_bas_surf.basis_curve().get_type();
            if a_bas_curv_type == CurveType::Line {
                *u1new = an_u1.max(-1.0e10);
                *u2new = an_u2.min(1.0e10);
            } else if a_bas_curv_type == CurveType::Parabola {
                let a_prb = a_bas_surf.basis_curve().parabola();
                let a_f = a_prb.focal_param * 0.5;
                let d_u = 2.0e5 * a_f.sqrt();
                *u1new = an_u1.max(-d_u);
                *u2new = an_u2.min(d_u);
            } else if a_bas_curv_type == CurveType::Hyperbola {
                *u1new = an_u1.max(-30.0);
                *u2new = an_u2.min(30.0);
            } else {
                *u1new = an_u1.max(-1.0e10);
                *u2new = an_u2.min(1.0e10);
            }
        }
    } else if a_type_of_bas_surf == SurfaceType::SurfaceOfRevolution {
        let a_bas_curv_type = a_bas_surf.basis_curve().get_type();
        if a_bas_curv_type == CurveType::Line {
            *v1new = v1new.max(-1.0e10);
            *v2new = v2new.min(1.0e10);
        } else if a_bas_curv_type == CurveType::Parabola {
            let a_prb = a_bas_surf.basis_curve().parabola();
            let a_f = a_prb.focal_param * 0.5;
            let d_v = 2.0e5 * a_f.sqrt();
            *v1new = v1new.max(-d_v);
            *v2new = v2new.min(d_v);
        } else if a_bas_curv_type == CurveType::Hyperbola {
            *v1new = v1new.max(-30.0);
            *v2new = v2new.min(30.0);
        } else {
            *v1new = v1new.max(-1.0e10);
            *v2new = v2new.min(1.0e10);
        }
    } else {
        *v1new = v1new.max(-1.0e10);
        *v2new = v2new.min(1.0e10);
    }
}

/// OCCT IntCurveSurface_InterUtils::ComputeTransitions (pxx L855-895).
pub(crate) fn compute_transitions<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>>(
    curve: &C,
    w: f64,
    trans_on_curve: &mut TransitionOnCurve,
    surface: &S,
    u: f64,
    v: f64,
) {
    let (_p_surf, d1u_s, d1v_s) = <ST as HSurfaceTool>::d1(surface, u, v);
    let n_surf = d1u_s.cross(d1v_s);
    // OCCT CurveTool::D1(curve, w, Psurf, D1U) reuses Psurf as the point
    // output (discarded) and D1U as the tangent output.
    let (_p_c, d1u_c) = <CT as HCurveTool>::d1(curve, w);
    let norm = n_surf.length();
    if norm > THE_TOLERANCE_ANGULAIRE && d1u_c.length_squared() > THE_TOLERANCE_ANGULAIRE {
        let d1u_norm = d1u_c.normalize_or_zero();
        let mut cos_dir = n_surf.dot(d1u_norm);
        cos_dir /= norm;
        if -cos_dir > THE_TOLERANCE_ANGULAIRE {
            //--  --Curve--->    <----Surface----
            *trans_on_curve = TransitionOnCurve::In;
        } else if cos_dir > THE_TOLERANCE_ANGULAIRE {
            //--  --Curve--->  ----Surface-->
            *trans_on_curve = TransitionOnCurve::Out;
        } else {
            *trans_on_curve = TransitionOnCurve::Tangent;
        }
    } else {
        *trans_on_curve = TransitionOnCurve::Tangent;
    }
}

/// OCCT IntCurveSurface_InterUtils::ComputeParamsOnQuadric (pxx L900-925).
pub(crate) fn compute_params_on_quadric<S, ST: HSurfaceTool<Surface = S>>(
    surface: &S,
    p: DVec3,
    u: &mut f64,
    v: &mut f64,
) {
    let surface_type = <ST as HSurfaceTool>::get_type(surface);
    match surface_type {
        SurfaceType::Plane => {
            let pln = <ST as HSurfaceTool>::plane(surface);
            let r = elslib_plane_parameters(p, pln.origin, pln.u_dir, pln.v_dir);
            *u = r.0;
            *v = r.1;
        }
        SurfaceType::Cylinder => {
            let cyl = <ST as HSurfaceTool>::cylinder(surface);
            let r = elslib_cylinder_parameters(p, cyl.origin, cyl.ref_dir, cyl.y_axis(), cyl.axis, cyl.radius);
            *u = r.0;
            *v = r.1;
        }
        SurfaceType::Cone => {
            let cone = <ST as HSurfaceTool>::cone(surface);
            let r = elslib_cone_parameters(
                p,
                cone.apex,
                cone.ref_dir,
                cone.axis.cross(cone.ref_dir).normalize_or_zero(),
                cone.axis,
                cone.radius,
                cone.half_angle_rad,
            );
            *u = r.0;
            *v = r.1;
        }
        SurfaceType::Sphere => {
            let sph = <ST as HSurfaceTool>::sphere(surface);
            let r = elslib_sphere_parameters(p, sph.center, sph.ref_dir, sph.ref_dir_perp(), sph.axis);
            *u = r.0;
            *v = r.1;
        }
        _ => {}
    }
}

/// OCCT IntCurveSurface_InterUtils::DoSurface (pxx L930-983) — sample the
/// 50x50 grid into `the_pnts_on_surface` (row-major, index (iU*50 + iV),
/// mirroring the 1-based SetValue(iU+1, iV+1)) and fill the bounding box.
pub(crate) fn do_surface<S, ST: HSurfaceTool<Surface = S>>(
    the_surface: &S,
    the_u0: f64,
    the_u1: f64,
    the_v0: f64,
    the_v1: f64,
    the_pnts_on_surface: &mut Vec<DVec3>,
    the_box_surface: &mut BndBox,
    the_gap: &mut f64,
) {
    let mut u;
    let mut v;
    let d_u = (the_u1 - the_u0) / 50.0;
    let d_v = (the_v1 - the_v0) / 50.0;

    the_pnts_on_surface.clear();
    for i_u in 0..50usize {
        if i_u == 0 {
            u = the_u0;
        } else if i_u == 49 {
            u = the_u1;
        } else {
            u = the_u0 + d_u * (i_u as f64);
        }

        for i_v in 0..50usize {
            if i_v == 0 {
                v = the_v0;
            } else if i_v == 49 {
                v = the_v1;
            } else {
                v = the_v0 + d_v * (i_v as f64);
            }

            let a_pnt = <ST as HSurfaceTool>::d0(the_surface, u, v);
            the_box_surface.add_point(a_pnt);
            the_pnts_on_surface.push(a_pnt);
        }
    }
    let u_res = <ST as HSurfaceTool>::u_resolution(the_surface, d_u);
    let v_res = <ST as HSurfaceTool>::v_resolution(the_surface, d_v);
    *the_gap = u_res.max(v_res);
}

/// OCCT IntCurveSurface_InterUtils::DoNewBounds (pxx L988-1107) — the grid
/// is row-major, index (iU - 1) * 50 + (iV - 1) for the 1-based
/// thePntsOnSurface.Value(iU, iV).
pub(crate) fn do_new_bounds<S, ST: HSurfaceTool<Surface = S>>(
    the_surface: &S,
    the_u0: f64,
    the_u1: f64,
    the_v0: f64,
    the_v1: f64,
    the_pnts_on_surface: &[DVec3],
    the_x: &[f64; 3],
    the_y: &[f64; 3],
    the_z: &[f64; 3],
    the_bounds: &mut [f64; 4],
) {
    the_bounds[0] = the_u0;
    the_bounds[1] = the_u1;
    the_bounds[2] = the_v0;
    the_bounds[3] = the_v1;

    let is_u_closed = <ST as HSurfaceTool>::is_u_closed(the_surface) || <ST as HSurfaceTool>::is_u_periodic(the_surface);
    let is_v_closed = <ST as HSurfaceTool>::is_v_closed(the_surface) || <ST as HSurfaceTool>::is_v_periodic(the_surface);
    let check_u = !is_u_closed;
    let check_v = !is_v_closed;

    let mut i_u_min = 50usize;
    let mut i_v_min = 50usize;
    let mut i_u_max = 1usize;
    let mut i_v_max = 1usize;

    for i in 0..2usize {
        for j in 0..2usize {
            for k in 0..2usize {
                let a_point = DVec3::new(the_x[i], the_y[j], the_z[k]);
                let mut dist_min = 1.0e100;
                let mut di_u = 0usize;
                let mut di_v = 0usize;
                for i_u in 1..=50usize {
                    for i_v in 1..=50usize {
                        let a_p = the_pnts_on_surface[(i_u - 1) * 50 + (i_v - 1)];
                        let dist = a_p.distance_squared(a_point);
                        if dist < dist_min {
                            dist_min = dist;
                            di_u = i_u;
                            di_v = i_v;
                        }
                    }
                }
                if di_u > 0 && di_u < i_u_min {
                    i_u_min = di_u;
                }
                if di_u > 0 && di_u > i_u_max {
                    i_u_max = di_u;
                }
                if di_v > 0 && di_v < i_v_min {
                    i_v_min = di_v;
                }
                if di_v > 0 && di_v > i_v_max {
                    i_v_max = di_v;
                }
            }
        }
    }

    let d_u = (the_u1 - the_u0) / 50.0;
    let d_v = (the_v1 - the_v0) / 50.0;

    let mut u_smin = the_u0 + d_u * ((i_u_min - 1) as f64);
    let mut u_smax = the_u0 + d_u * ((i_u_max - 1) as f64);
    let mut v_smin = the_v0 + d_v * ((i_v_min - 1) as f64);
    let mut v_smax = the_v0 + d_v * ((i_v_max - 1) as f64);

    if u_smin > u_smax {
        std::mem::swap(&mut u_smax, &mut u_smin);
    }
    if v_smin > v_smax {
        std::mem::swap(&mut v_smax, &mut v_smin);
    }

    u_smin -= 1.5 * d_u;
    if u_smin < the_u0 {
        u_smin = the_u0;
    }
    u_smax += 1.5 * d_u;
    if u_smax > the_u1 {
        u_smax = the_u1;
    }
    v_smin -= 1.5 * d_v;
    if v_smin < the_v0 {
        v_smin = the_v0;
    }
    v_smax += 1.5 * d_v;
    if v_smax > the_v1 {
        v_smax = the_v1;
    }

    if check_u {
        the_bounds[0] = u_smin;
        the_bounds[1] = u_smax;
    }
    if check_v {
        the_bounds[2] = v_smin;
        the_bounds[3] = v_smax;
    }
}

/// OCCT IntCurveSurface_InterUtils::ComputeAppendPoint (pxx L1116-1176) —
/// parameter validation + transition; returns the point when it should be
/// appended.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_append_point<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>>(
    the_curve: &C,
    the_lw: f64,
    the_surface: &S,
    the_su: f64,
    the_sv: f64,
) -> Option<IntersectionPoint> {
    let w0 = <CT as HCurveTool>::first_parameter(the_curve);
    let w1 = <CT as HCurveTool>::last_parameter(the_curve);
    let u0 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
    let u1 = <ST as HSurfaceTool>::last_u_parameter(the_surface);
    let v0 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
    let v1 = <ST as HSurfaceTool>::last_v_parameter(the_surface);

    let mut w = the_lw;
    let mut u = the_su;
    let mut v = the_sv;

    let a_c_type = <CT as HCurveTool>::get_type(the_curve);

    if <CT as HCurveTool>::is_periodic(the_curve)
        || a_c_type == CurveType::Circle
        || a_c_type == CurveType::Ellipse
    {
        w = in_period(w, w0, w0 + <CT as HCurveTool>::period(the_curve));
    }

    if (w0 - w) >= THE_TOLTANGENCY || (w - w1) >= THE_TOLTANGENCY {
        return None;
    }

    let a_s_type = <ST as HSurfaceTool>::get_type(the_surface);
    if <ST as HSurfaceTool>::is_u_periodic(the_surface)
        || a_s_type == SurfaceType::Cylinder
        || a_s_type == SurfaceType::Cone
        || a_s_type == SurfaceType::Sphere
    {
        u = in_period(u, u0, u0 + <ST as HSurfaceTool>::u_period(the_surface));
    }

    if <ST as HSurfaceTool>::is_v_periodic(the_surface) {
        v = in_period(v, v0, v0 + <ST as HSurfaceTool>::v_period(the_surface));
    }

    if (u0 - u) >= THE_TOLTANGENCY || (u - u1) >= THE_TOLTANGENCY {
        return None;
    }
    if (v0 - v) >= THE_TOLTANGENCY || (v - v1) >= THE_TOLTANGENCY {
        return None;
    }

    let mut trans_on_curve = TransitionOnCurve::Tangent;
    compute_transitions::<C, CT, S, ST>(the_curve, w, &mut trans_on_curve, the_surface, u, v);
    let p = <CT as HCurveTool>::value(the_curve, w);
    Some(IntersectionPoint::with_values(p, u, v, w, trans_on_curve))
}

/// OCCT IntCurveSurface_InterUtils::ProcessIntAna (pxx L1187-1228) —
/// returns done and fills theIsParallel/thePoints.
#[allow(clippy::too_many_arguments)]
pub(crate) fn process_int_ana<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>>(
    the_curve: &C,
    the_surface: &S,
    the_int_ana: &IntConicQuad,
    the_is_parallel: &mut bool,
    the_points: &mut Vec<IntersectionPoint>,
) -> bool {
    *the_is_parallel = false;
    the_points.clear();

    if !the_int_ana.is_done() {
        return false;
    }

    if the_int_ana.is_in_quadric() || the_int_ana.is_parallel() {
        *the_is_parallel = true;
        return true;
    }

    let nbp = the_int_ana.nb_points();
    let mut u = 0.0;
    let mut v = 0.0;
    for i in 1..=nbp {
        let p = the_int_ana.point(i);
        let w = the_int_ana.param_on_conic(i);
        compute_params_on_quadric::<S, ST>(the_surface, p, &mut u, &mut v);

        if let Some(a_point) =
            compute_append_point::<C, CT, S, ST>(the_curve, w, the_surface, u, v)
        {
            the_points.push(a_point);
        }
    }
    true
}

/// The QuadCurvExactType template parameter of PerformCurveQuadric — the
/// constructor + results contract (the OCCT
/// IntCurveSurface_TheQuadCurvExactHInter API).
pub trait QuadCurvExactLike<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>> {
    /// OCCT TheQuadCurvExactHInter(S, C).
    fn construct(s: &S, c: &C) -> Self;
    /// OCCT IsDone().
    fn is_done(&self) -> bool;
    /// OCCT NbRoots().
    fn nb_roots(&self) -> usize;
    /// OCCT Root(Index) — 1-based.
    fn root(&self, index: usize) -> f64;
}

/// OCCT IntCurveSurface_InterUtils::PerformCurveQuadric (pxx L1238-1274).
pub(crate) fn perform_curve_quadric<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>, Q>(
    the_curve: &C,
    the_surface: &S,
    the_points: &mut Vec<IntersectionPoint>,
) where
    Q: QuadCurvExactLike<C, CT, S, ST>,
{
    the_points.clear();

    let quad_curv = Q::construct(the_surface, the_curve);
    if quad_curv.is_done() {
        let nb_roots = quad_curv.nb_roots();
        let mut u = 0.0;
        let mut v = 0.0;
        for i in 1..=nb_roots {
            let w = quad_curv.root(i);
            compute_params_on_quadric::<S, ST>(the_surface, <CT as HCurveTool>::value(the_curve, w), &mut u, &mut v);

            if let Some(a_point) =
                compute_append_point::<C, CT, S, ST>(the_curve, w, the_surface, u, v)
            {
                the_points.push(a_point);
            }
        }
    }
}

/// OCCT IntCurveSurface_InterUtils::ProcessLinTorus (pxx L1282-1315) —
/// returns false when the fallback to the polyhedron path is needed.
pub(crate) fn process_lin_torus<C: ?Sized, CT: HCurveTool<Curve = C>, S, ST: HSurfaceTool<Surface = S>>(
    the_line: &Line3,
    the_curve: &C,
    the_surface: &S,
    the_points: &mut Vec<IntersectionPoint>,
) -> bool {
    the_points.clear();

    let intlintorus = IntLinTorus::new_line_torus(the_line, &<ST as HSurfaceTool>::torus(the_surface));
    if !intlintorus.is_done() {
        return false;
    }

    let nbp = intlintorus.nb_points();
    for i in 1..=nbp {
        let w = intlintorus.param_on_line(i);
        let (fi, theta) = intlintorus.param_on_torus(i);

        if let Some(a_point) =
            compute_append_point::<C, CT, S, ST>(the_curve, w, the_surface, fi, theta)
        {
            the_points.push(a_point);
        }
    }
    true
}

/// OCCT IntCurveSurface_InterUtils::SortedStartPoints (pxx L1318-1339).
#[derive(Debug, Clone, Default)]
pub struct SortedStartPoints {
    pub tab_u: Vec<f64>,
    pub tab_v: Vec<f64>,
    pub tab_w: Vec<f64>,
}

impl SortedStartPoints {
    pub fn clear(&mut self) {
        self.tab_u.clear();
        self.tab_v.clear();
        self.tab_w.clear();
    }

    pub fn size(&self) -> usize {
        self.tab_u.len()
    }

    pub fn append(&mut self, the_u: f64, the_v: f64, the_w: f64) {
        self.tab_u.push(the_u);
        self.tab_v.push(the_v);
        self.tab_w.push(the_w);
    }
}

/// OCCT IntCurveSurface_InterUtils::CollectInterferencePoints (pxx
/// L1345-1376).
pub(crate) fn collect_interference_points(
    the_interference: &crate::geomalgo::intf_interference_polygon_polyhedron::InterferencePolygonPolyhedron,
    the_polyhedron: &ThePolyhedronOfHInter,
    the_polygon: &ThePolygonOfHInter,
    the_points: &mut SortedStartPoints,
) {
    the_points.clear();

    let nb_section_points = the_interference.interf.nb_section_points();
    let nb_tangent_zones = the_interference.interf.nb_tangent_zones();

    for i in 1..=nb_section_points {
        let sp = the_interference.interf.pnt_value(i);
        let (u, v, w) = section_point_to_parameters(sp, the_polyhedron, the_polygon);
        the_points.append(u, v, w);
    }

    for i in 1..=nb_tangent_zones {
        let tz = the_interference.interf.zone_value(i);
        let nbpnts = tz.number_of_points();
        for j in 1..=nbpnts {
            let sp = tz.get_point(j);
            let (u, v, w) = section_point_to_parameters(&sp, the_polyhedron, the_polygon);
            the_points.append(u, v, w);
        }
    }
}

/// OCCT IntCurveSurface_InterUtils::SortStartPoints (pxx L1380-1447) — the
/// three bubble sorts (by W, then U for same W, then V for same W and U)
/// with the ptol = 10·PConfusion W/U collapsing.
pub(crate) fn sort_start_points(the_points: &mut SortedStartPoints) {
    let nb_start_points = the_points.size();
    if nb_start_points == 0 {
        return;
    }

    let ptol = 10.0 * PCONFUSION;

    // Sort by W.
    loop {
        let mut triok = true;
        for i in 1..nb_start_points {
            let im1 = i - 1;
            if the_points.tab_w[i] < the_points.tab_w[im1] {
                the_points.tab_w.swap(i, im1);
                the_points.tab_u.swap(i, im1);
                the_points.tab_v.swap(i, im1);
                triok = false;
            }
        }
        if triok {
            break;
        }
    }

    // Sort by U for same W.
    loop {
        let mut triok = true;
        for i in 1..nb_start_points {
            let im1 = i - 1;
            if the_points.tab_w[i] - the_points.tab_w[im1] < ptol {
                the_points.tab_w[i] = the_points.tab_w[im1];
                if the_points.tab_u[i] < the_points.tab_u[im1] {
                    the_points.tab_u.swap(i, im1);
                    the_points.tab_v.swap(i, im1);
                    triok = false;
                }
            }
        }
        if triok {
            break;
        }
    }

    // Sort by V for same W and U.
    loop {
        let mut triok = true;
        for i in 1..nb_start_points {
            let im1 = i - 1;
            if the_points.tab_w[i] - the_points.tab_w[im1] < ptol
                && the_points.tab_u[i] - the_points.tab_u[im1] < ptol
            {
                the_points.tab_u[i] = the_points.tab_u[im1];
                if the_points.tab_v[i] < the_points.tab_v[im1] {
                    the_points.tab_v.swap(i, im1);
                    triok = false;
                }
            }
        }
        if triok {
            break;
        }
    }
}

/// OCCT IntCurveSurface_InterUtils::ProcessSortedPoints (pxx L1455-1519) —
/// drive the exact intersection from the sorted start points.
#[allow(clippy::too_many_arguments)]
pub(crate) fn process_sorted_points<C: ?Sized, CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>, S, ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>, F>(
    the_exact_inter: &mut IntCS<S, C, ST, CT, F>,
    the_rsnld: &mut FunctionSetRoot,
    the_points: &SortedStartPoints,
    the_u0: f64,
    the_u1: f64,
    the_v0: f64,
    the_v1: f64,
    the_winf: f64,
    the_wsup: f64,
    the_curve: &C,
    the_surface: &S,
    the_result: &mut Vec<IntersectionPoint>,
) where
    F: rcad_kernel::math::function_set_root::FunctionSetWithDerivatives
        + crate::geomalgo::int_imp::int_cs::ZerCSAccessors<S, C>,
{
    the_result.clear();

    let nb_start_points = the_points.size();
    if nb_start_points == 0 {
        return;
    }

    let ptol = 10.0 * PCONFUSION;
    let mut su = 0.0;
    let mut sv = 0.0;
    let mut sw = 0.0;

    for i in 0..nb_start_points {
        let u = the_points.tab_u[i];
        let v = the_points.tab_v[i];
        let w = the_points.tab_w[i];

        if i == 0 {
            su = u - 1.0;
        }

        if (u - su).abs() > ptol || (v - sv).abs() > ptol || (w - sw).abs() > ptol {
            the_exact_inter.perform(u, v, w, the_rsnld, the_u0, the_u1, the_v0, the_v1, the_winf, the_wsup);
            if the_exact_inter.is_done() && !the_exact_inter.is_empty() {
                let w = the_exact_inter.parameter_on_curve();
                let (u, v) = the_exact_inter.parameter_on_surface();

                if let Some(a_point) =
                    compute_append_point::<C, CT, S, ST>(the_curve, w, the_surface, u, v)
                {
                    the_result.push(a_point);
                }
            }
        }
        su = the_points.tab_u[i];
        sv = the_points.tab_v[i];
        sw = the_points.tab_w[i];
    }
}

/// OCCT IntCurveSurface_InterUtils::UVBounds (pxx L1522-1544).
#[derive(Debug, Clone, Copy, Default)]
pub struct UVBounds {
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

impl UVBounds {
    pub fn new(the_u0: f64, the_u1: f64, the_v0: f64, the_v1: f64) -> Self {
        UVBounds {
            u0: the_u0,
            u1: the_u1,
            v0: the_v0,
            v1: the_v1,
        }
    }
}

/// OCCT IntCurveSurface_InterUtils::DecomposeSurfaceIntervals (pxx
/// L1549-1610).
pub(crate) fn decompose_surface_intervals<S, ST: HSurfaceTool<Surface = S>>(
    the_surface: &S,
    the_intervals: &mut Vec<UVBounds>,
) {
    the_intervals.clear();

    let nb_u_on_s = <ST as HSurfaceTool>::nb_u_intervals(the_surface, GeomAbsShape::C2);
    let nb_v_on_s = <ST as HSurfaceTool>::nb_v_intervals(the_surface, GeomAbsShape::C2);

    if nb_u_on_s > 1 {
        let mut tab_u = vec![0.0f64; nb_u_on_s + 1];
        <ST as HSurfaceTool>::u_intervals(the_surface, &mut tab_u, GeomAbsShape::C2);

        for iu in 1..=nb_u_on_s {
            let u0 = tab_u[iu - 1];
            let u1 = tab_u[iu];

            if nb_v_on_s > 1 {
                let mut tab_v = vec![0.0f64; nb_v_on_s + 1];
                <ST as HSurfaceTool>::v_intervals(the_surface, &mut tab_v, GeomAbsShape::C2);
                for iv in 1..=nb_v_on_s {
                    let v0 = tab_v[iv - 1];
                    let v1 = tab_v[iv];
                    the_intervals.push(UVBounds::new(u0, u1, v0, v1));
                }
            } else {
                let v0 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
                let v1 = <ST as HSurfaceTool>::last_v_parameter(the_surface);
                the_intervals.push(UVBounds::new(u0, u1, v0, v1));
            }
        }
    } else if nb_v_on_s > 1 {
        let u0 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
        let u1 = <ST as HSurfaceTool>::last_u_parameter(the_surface);

        let mut tab_v = vec![0.0f64; nb_v_on_s + 1];
        <ST as HSurfaceTool>::v_intervals(the_surface, &mut tab_v, GeomAbsShape::C2);

        for iv in 1..=nb_v_on_s {
            let v0 = tab_v[iv - 1];
            let v1 = tab_v[iv];
            the_intervals.push(UVBounds::new(u0, u1, v0, v1));
        }
    } else {
        let u0 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
        let u1 = <ST as HSurfaceTool>::last_u_parameter(the_surface);
        let v0 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
        let v1 = <ST as HSurfaceTool>::last_v_parameter(the_surface);
        the_intervals.push(UVBounds::new(u0, u1, v0, v1));
    }
}


/// OCCT IntCurveSurface_InterUtils::ClampUVParameters (pxx L1614-1633) —
/// protection from double overflow (bug26525).
pub(crate) fn clamp_uv_parameters(the_u1: &mut f64, the_u2: &mut f64, the_v1: &mut f64, the_v2: &mut f64) {
    const THE_PARAM_LIMIT: f64 = 1.0e50;
    if *the_u1 < -THE_PARAM_LIMIT {
        *the_u1 = -THE_PARAM_LIMIT;
    }
    if *the_u2 > THE_PARAM_LIMIT {
        *the_u2 = THE_PARAM_LIMIT;
    }
    if *the_v1 < -THE_PARAM_LIMIT {
        *the_v1 = -THE_PARAM_LIMIT;
    }
    if *the_v2 > THE_PARAM_LIMIT {
        *the_v2 = THE_PARAM_LIMIT;
    }
}

