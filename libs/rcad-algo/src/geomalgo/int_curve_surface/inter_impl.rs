//! IntCurveSurface_Inter.pxx (L1-1006) — the engine namespace
//! IntCurveSurface_InterImpl.  Every function is generic over the curve /
//! surface types and their tool traits exactly like the OCCT templates; the
//! OCCT callback lambdas (ResetFunc / PerformBoundsFunc / AppendFunc / ...)
//! are the [`HInterHost`] methods they wrap.

use rcad_kernel::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::math::bnd::BoundSortBox;
use rcad_kernel::math::function_set_root::FunctionSetRoot;
use rcad_kernel::precision::is_infinite_value;

use crate::geomalgo::int_imp::int_cs::IntCS;
use crate::geomalgo::int_imp::zer_cs_par_func::ZerCSParFunc;
use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};
use crate::geomalgo::int_patch::int_conic_quad::IntConicQuad;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};
use crate::geomalgo::intf::IntfTool;
use crate::geomalgo::int_curv_surf::{ThePolyhedronToolOfHInter, ThePolygonToolOfHInter};
use crate::geomalgo::intf_interference_polygon_polyhedron::InterferencePolygonPolyhedron;

use super::inter_utils::{
    clamp_uv_parameters, collect_interference_points, decompose_surface_intervals,
    est_lim_for_inf_extr, est_lim_for_inf_offs, est_lim_for_inf_revl, est_lim_for_inf_surf,
    perform_curve_quadric, process_int_ana, process_lin_torus, process_sorted_points,
    sort_start_points, QuadCurvExactLike, SortedStartPoints,
};
use super::{HCurveTool, HInterHost, HSurfaceTool, IntersectionPoint, SurfaceType, UVBounds};

/// OCCT IntCurveSurface_InterImpl THE_TOLTANGENCY (pxx L49).
pub(crate) const IMPL_THE_TOLTANGENCY: f64 = 0.00000001;
/// OCCT IntCurveSurface_InterImpl THE_TOLERANCE_ANGULAIRE (pxx L50).
pub(crate) const IMPL_THE_TOLERANCE_ANGULAIRE: f64 = 1.0e-12;
/// OCCT IntCurveSurface_InterImpl THE_TOLERANCE (pxx L51).
pub(crate) const IMPL_THE_TOLERANCE: f64 = 0.00000001;
/// OCCT IntCurveSurface_InterImpl THE_NBSAMPLESONCIRCLE (pxx L52).
pub(crate) const THE_NBSAMPLESONCIRCLE: usize = 32;
/// OCCT IntCurveSurface_InterImpl THE_NBSAMPLESONELLIPSE (pxx L53).
pub(crate) const THE_NBSAMPLESONELLIPSE: usize = 32;
/// OCCT IntCurveSurface_InterImpl THE_NBSAMPLESONPARAB (pxx L54).
pub(crate) const THE_NBSAMPLESONPARAB: usize = 16;
/// OCCT IntCurveSurface_InterImpl THE_NBSAMPLESONHYPR (pxx L55).
pub(crate) const THE_NBSAMPLESONHYPR: usize = 32;

/// OCCT IntCurveSurface_InterImpl::Perform (pxx L63-86) — perform the
/// intersection decomposing the surface by C2 intervals.
pub(crate) fn perform<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_surface: &S,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    host.reset_fields();
    *host.done_flag() = true;

    let mut a_intervals: Vec<UVBounds> = Vec::new();
    decompose_surface_intervals::<S, ST>(the_surface, &mut a_intervals);

    for a_bounds in &a_intervals {
        host.perform_bounds(the_curve, the_surface, a_bounds.u0, a_bounds.v0, a_bounds.u1, a_bounds.v1);
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformBounds (pxx L97-182).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_bounds<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    // Protection from double type overflow (bug26525).
    let mut u_u1 = the_u1;
    let mut u_u2 = the_u2;
    let mut v_v1 = the_v1;
    let mut v_v2 = the_v2;
    clamp_uv_parameters(&mut u_u1, &mut u_u2, &mut v_v1, &mut v_v2);

    let a_curve_type = <CT as HCurveTool>::get_type(the_curve);

    match a_curve_type {
        rcad_kernel::base::proj_lib::CurveType::Line => {
            let line = <CT as HCurveTool>::line(the_curve);
            host.perform_conic_line(&line, the_curve, the_surface, u_u1, v_v1, u_u2, v_v2);
        }
        rcad_kernel::base::proj_lib::CurveType::Circle => {
            let circle = <CT as HCurveTool>::circle(the_curve);
            host.perform_conic_circle(&circle, the_curve, the_surface, u_u1, v_v1, u_u2, v_v2);
        }
        rcad_kernel::base::proj_lib::CurveType::Ellipse => {
            let ellipse = <CT as HCurveTool>::ellipse(the_curve);
            host.perform_conic_ellipse(&ellipse, the_curve, the_surface, u_u1, v_v1, u_u2, v_v2);
        }
        rcad_kernel::base::proj_lib::CurveType::Parabola => {
            let parabola = <CT as HCurveTool>::parabola(the_curve);
            host.perform_conic_parabola(&parabola, the_curve, the_surface, u_u1, v_v1, u_u2, v_v2);
        }
        rcad_kernel::base::proj_lib::CurveType::Hyperbola => {
            let hyperbola = <CT as HCurveTool>::hyperbola(the_curve);
            host.perform_conic_hyperbola(&hyperbola, the_curve, the_surface, u_u1, v_v1, u_u2, v_v2);
        }
        _ => {
            let nb_intervals_on_curve = <CT as HCurveTool>::nb_intervals(the_curve, rcad_kernel::math::GeomAbsShape::C2);
            let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
            if a_surface_type != SurfaceType::Plane
                && a_surface_type != SurfaceType::Cylinder
                && a_surface_type != SurfaceType::Cone
                && a_surface_type != SurfaceType::Sphere
            {
                if nb_intervals_on_curve > 1 {
                    let mut tab_w = vec![0.0f64; nb_intervals_on_curve + 1];
                    <CT as HCurveTool>::intervals(the_curve, &mut tab_w, rcad_kernel::math::GeomAbsShape::C2);
                    for i in 1..=nb_intervals_on_curve {
                        let u1 = tab_w[i - 1];
                        let u2 = tab_w[i];

                        let defl = 0.1;
                        let nb_min = 10usize;
                        let a_pars = <CT as HCurveTool>::sample_pars(the_curve, u1, u2, defl, nb_min);

                        let polygon = ThePolygonOfHInter::new_tool_params::<C, CT>(the_curve, &a_pars);
                        host.internal_perform_bounds(the_curve, &polygon, the_surface, u_u1, v_v1, u_u2, v_v2);
                    }
                } else {
                    let u1 = <CT as HCurveTool>::first_parameter(the_curve);
                    let u2 = <CT as HCurveTool>::last_parameter(the_curve);

                    let defl = 0.1;
                    let nb_min = 10usize;
                    let a_pars = <CT as HCurveTool>::sample_pars(the_curve, u1, u2, defl, nb_min);

                    let polygon = ThePolygonOfHInter::new_tool_params::<C, CT>(the_curve, &a_pars);
                    host.internal_perform_bounds(the_curve, &polygon, the_surface, u_u1, v_v1, u_u2, v_v2);
                }
            } else {
                host.internal_perform_curve_quadric(the_curve, the_surface);
            }
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformPolygon (pxx L185-218).
pub(crate) fn perform_polygon<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    host.reset_fields();
    *host.done_flag() = true;
    let u1 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
    let v1 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
    let u2 = <ST as HSurfaceTool>::last_u_parameter(the_surface);
    let v2 = <ST as HSurfaceTool>::last_v_parameter(the_surface);
    let mut nbsu = <ST as HSurfaceTool>::nb_samples_u(the_surface, u1, u2);
    let mut nbsv = <ST as HSurfaceTool>::nb_samples_v(the_surface, v1, v2);
    if nbsu > 40 {
        nbsu = 40;
    }
    if nbsv > 40 {
        nbsv = 40;
    }
    let polyhedron = ThePolyhedronOfHInter::new_tool::<S, ST>(the_surface, nbsu, nbsv, u1, v1, u2, v2);
    host.perform_polygon_polyhedron(the_curve, the_polygon, the_surface, &polyhedron);
}

/// OCCT IntCurveSurface_InterImpl::PerformPolyhedron (pxx L220-242).
pub(crate) fn perform_polyhedron<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    host.reset_fields();
    *host.done_flag() = true;
    let u1 = <CT as HCurveTool>::first_parameter(the_curve);
    let u2 = <CT as HCurveTool>::last_parameter(the_curve);
    let polygon = ThePolygonOfHInter::new_tool::<C, CT>(the_curve, <CT as HCurveTool>::nb_samples(the_curve, u1, u2));
    host.perform_polygon_polyhedron(the_curve, &polygon, the_surface, the_polyhedron);
}

/// OCCT IntCurveSurface_InterImpl::PerformPolygonPolyhedron (pxx L245-268).
pub(crate) fn perform_polygon_polyhedron<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    host.reset_fields();
    *host.done_flag() = true;
    let u1 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
    let v1 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
    let u2 = <ST as HSurfaceTool>::last_u_parameter(the_surface);
    let v2 = <ST as HSurfaceTool>::last_v_parameter(the_surface);
    host.internal_perform(the_curve, the_polygon, the_surface, the_polyhedron, u1, v1, u2, v2);
}

/// OCCT IntCurveSurface_InterImpl::PerformPolygonPolyhedronBSB (pxx
/// L271-295).
pub(crate) fn perform_polygon_polyhedron_bsb<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    the_bsb: &mut BoundSortBox,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    host.reset_fields();
    *host.done_flag() = true;
    let u1 = <ST as HSurfaceTool>::first_u_parameter(the_surface);
    let v1 = <ST as HSurfaceTool>::first_v_parameter(the_surface);
    let u2 = <ST as HSurfaceTool>::last_u_parameter(the_surface);
    let v2 = <ST as HSurfaceTool>::last_v_parameter(the_surface);
    host.internal_perform_bsb(the_curve, the_polygon, the_surface, the_polyhedron, u1, v1, u2, v2, the_bsb);
}

/// OCCT IntCurveSurface_InterImpl::InternalPerformBSB (pxx L298-355).
#[allow(clippy::too_many_arguments)]
pub(crate) fn internal_perform_bsb<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    the_u0: f64,
    the_v0: f64,
    the_u1: f64,
    the_v1: f64,
    the_bsb: &mut BoundSortBox,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let interference = InterferencePolygonPolyhedron::new_polygon_polyhedron_grid::<
        ThePolygonOfHInter,
        ThePolyhedronOfHInter,
        ThePolygonToolOfHInter,
        ThePolyhedronToolOfHInter,
    >(the_polygon, the_polyhedron, the_bsb);

    internal_perform_body::<C, CT, S, ST, H>(
        interference, the_curve, the_polygon, the_surface, the_polyhedron, the_u0, the_v0, the_u1, the_v1, host,
    );
}

/// The shared tail of InternalPerform / InternalPerformBSB (OCCT has it
/// duplicated in the two templates L319-355 and L368-414).  The
/// CSFunctionType / ExactInterType template parameters are the concrete
/// instantiation TheCSFunctionOfHInter = ZerCSParFunc and TheExactHInter =
/// IntCS (IntCurveSurface_TheCSFunctionOfHInter_0.cxx /
/// TheExactHInter_0.cxx).
#[allow(clippy::too_many_arguments)]
fn internal_perform_body<C: ?Sized, CT, S: ?Sized, ST, H>(
    interference: InterferencePolygonPolyhedron,
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    the_u0: f64,
    the_v0: f64,
    the_u1: f64,
    the_v1: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let mut theicsfunction = ZerCSParFunc::<S, C, ST, CT>::new(the_surface, the_curve);
    let mut rsnld = FunctionSetRoot::new(&theicsfunction, &[0.0; 3], 100);
    let mut intersection_exacte = IntCS::new(&mut theicsfunction, IMPL_THE_TOLTANGENCY);

    let winf = the_polygon.inf_parameter();
    let wsup = the_polygon.sup_parameter();

    let mut a_start_points = super::inter_utils::SortedStartPoints::default();
    collect_interference_points(&interference, the_polyhedron, the_polygon, &mut a_start_points);
    sort_start_points(&mut a_start_points);

    let mut a_result_points: Vec<IntersectionPoint> = Vec::new();
    process_sorted_points::<C, CT, S, ST, _>(
        &mut intersection_exacte,
        &mut rsnld,
        &a_start_points,
        the_u0,
        the_u1,
        the_v0,
        the_v1,
        winf,
        wsup,
        the_curve,
        the_surface,
        &mut a_result_points,
    );

    for a_point in &a_result_points {
        host.append(a_point);
    }
}

/// OCCT IntCurveSurface_InterImpl::InternalPerform (pxx L357-414).
#[allow(clippy::too_many_arguments)]
pub(crate) fn internal_perform<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_polyhedron: &ThePolyhedronOfHInter,
    the_u0: f64,
    the_v0: f64,
    the_u1: f64,
    the_v1: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let interference = InterferencePolygonPolyhedron::new_polygon_polyhedron::<
        ThePolygonOfHInter,
        ThePolyhedronOfHInter,
        ThePolygonToolOfHInter,
        ThePolyhedronToolOfHInter,
    >(the_polygon, the_polyhedron);

    internal_perform_body::<C, CT, S, ST, H>(
        interference, the_curve, the_polygon, the_surface, the_polyhedron, the_u0, the_v0, the_u1, the_v1, host,
    );
}

/// OCCT IntCurveSurface_InterImpl::InternalPerformCurveQuadric (pxx
/// L416-439).
pub(crate) fn internal_perform_curve_quadric<C: ?Sized, CT, S: ?Sized, ST, Q, H>(
    the_curve: &C,
    the_surface: &S,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    Q: QuadCurvExactLike<C, CT, S, ST>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let mut a_points: Vec<IntersectionPoint> = Vec::new();

    perform_curve_quadric::<C, CT, S, ST, Q>(the_curve, the_surface, &mut a_points);

    for a_point in &a_points {
        host.append(a_point);
    }
}

/// OCCT IntCurveSurface_InterImpl::InternalPerformPolygonBounds (pxx
/// L441-520).
#[allow(clippy::too_many_arguments)]
pub(crate) fn internal_perform_polygon_bounds<C: ?Sized, CT, S: ?Sized, ST, Q, H>(
    the_curve: &C,
    the_polygon: &ThePolygonOfHInter,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    Q: QuadCurvExactLike<C, CT, S, ST>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    if a_surface_type != SurfaceType::Plane
        && a_surface_type != SurfaceType::Cylinder
        && a_surface_type != SurfaceType::Cone
        && a_surface_type != SurfaceType::Sphere
    {
        if a_surface_type != SurfaceType::BSplineSurface {
            let mut nbsu = <ST as HSurfaceTool>::nb_samples_u(the_surface, the_u1, the_u2);
            let mut nbsv = <ST as HSurfaceTool>::nb_samples_v(the_surface, the_v1, the_v2);
            if nbsu > 40 {
                nbsu = 40;
            }
            if nbsv > 40 {
                nbsv = 40;
            }
            let polyhedron =
                ThePolyhedronOfHInter::new_tool::<S, ST>(the_surface, nbsu, nbsv, the_u1, the_v1, the_u2, the_v2);
            host.internal_perform(the_curve, the_polygon, the_surface, &polyhedron, the_u1, the_v1, the_u2, the_v2);
        } else {
            // OCCT L489-501: the BSplineSurface branch samples through
            // Adaptor3d_TopolTool::SamplePnts — the TopolTool dependency is
            // runway 2a-4.
            let _ = <ST as HSurfaceTool>::u_trim(the_surface, the_u1, the_u2, 1.0e-9);
        }
    } else {
        internal_perform_curve_quadric::<C, CT, S, ST, Q, H>(the_curve, the_surface, host);
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformConicSurfLine (pxx L522-708).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_conic_surf_line<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_line: &Line3,
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    use crate::geomalgo::int_surf::quadric::Quadric;
    use rcad_kernel::base::proj_lib::CurveType;

    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    let mut is_ana_processed = true;

    match a_surface_type {
        SurfaceType::Plane => {
            let lin_plane = IntConicQuad::new_line_plane(the_line, &<ST as HSurfaceTool>::plane(the_surface), IMPL_THE_TOLERANCE_ANGULAIRE);
            host.append_int_ana(the_curve, the_surface, &lin_plane);
        }
        SurfaceType::Cylinder => {
            // OCCT IntAna_IntConicQuad(theLine, gp_Cylinder) — the cylinder
            // converts implicitly to IntAna_Quadric.
            let lin_cylinder =
                IntConicQuad::new_line_quadric(the_line, &Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)));
            host.append_int_ana(the_curve, the_surface, &lin_cylinder);
        }
        SurfaceType::Sphere => {
            let lin_sphere =
                IntConicQuad::new_line_quadric(the_line, &Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)));
            host.append_int_ana(the_curve, the_surface, &lin_sphere);
        }
        SurfaceType::Torus => {
            let mut a_points: Vec<IntersectionPoint> = Vec::new();
            if process_lin_torus::<C, CT, S, ST>(the_line, the_curve, the_surface, &mut a_points) {
                for a_point in &a_points {
                    host.append(a_point);
                }
            } else {
                is_ana_processed = false;
            }
        }
        SurfaceType::Cone => {
            let correction = 1.0e5 * rcad_kernel::precision::ANGULAR;
            let cn = <ST as HSurfaceTool>::cone(the_surface);
            if cn.half_angle_rad.abs() < std::f64::consts::PI / 2.0 - correction {
                let lin_cone = IntConicQuad::new_line_quadric(the_line, &Quadric::from_cone(&cn));
                host.append_int_ana(the_curve, the_surface, &lin_cone);
            } else {
                is_ana_processed = false;
            }
        }
        _ => is_ana_processed = false,
    }

    if !is_ana_processed {
        let mut nbsu = <ST as HSurfaceTool>::nb_samples_u(the_surface, the_u1, the_u2);
        let mut nbsv = <ST as HSurfaceTool>::nb_samples_v(the_surface, the_v1, the_v2);
        if nbsu < 1 {
            nbsu = 1;
        }
        if nbsv < 1 {
            nbsv = 1;
        }

        let u1inf = is_infinite_value(the_u1);
        let u2inf = is_infinite_value(the_u2);
        let v1inf = is_infinite_value(the_v1);
        let v2inf = is_infinite_value(the_v2);

        let mut u1new = the_u1;
        let mut u2new = the_u2;
        let mut v1new = the_v1;
        let mut v2new = the_v2;
        let mut no_intersection = false;

        if u1inf || u2inf || v1inf || v2inf {
            if a_surface_type == SurfaceType::SurfaceOfExtrusion {
                est_lim_for_inf_extr::<S, ST>(
                    the_line,
                    the_surface,
                    false,
                    nbsu,
                    u1inf,
                    u2inf,
                    v1inf,
                    v2inf,
                    &mut u1new,
                    &mut u2new,
                    &mut v1new,
                    &mut v2new,
                    &mut no_intersection,
                );
            } else if a_surface_type == SurfaceType::SurfaceOfRevolution {
                est_lim_for_inf_revl::<S, ST>(
                    the_line,
                    the_surface,
                    u1inf,
                    u2inf,
                    v1inf,
                    v2inf,
                    &mut u1new,
                    &mut u2new,
                    &mut v1new,
                    &mut v2new,
                    &mut no_intersection,
                );
            } else if a_surface_type == SurfaceType::OffsetSurface {
                est_lim_for_inf_offs::<S, ST>(
                    the_line,
                    the_surface,
                    nbsu,
                    u1inf,
                    u2inf,
                    v1inf,
                    v2inf,
                    &mut u1new,
                    &mut u2new,
                    &mut v1new,
                    &mut v2new,
                    &mut no_intersection,
                );
            } else {
                est_lim_for_inf_surf(&mut u1new, &mut u2new, &mut v1new, &mut v2new);
            }
        }

        if no_intersection {
            return;
        }

        if nbsu < 20 {
            nbsu = 20;
        }
        if nbsv < 20 {
            nbsv = 20;
        }

        let polyhedron = ThePolyhedronOfHInter::new_tool::<S, ST>(the_surface, nbsu, nbsv, u1new, v1new, u2new, v2new);
        let mut bnd_tool = crate::geomalgo::intf::IntfTool::new();
        let mut box_line = BndBox::new();
        bnd_tool.lin_box(the_line, polyhedron.bounding(), &mut box_line);

        for nbseg in 1..=bnd_tool.nb_segments() {
            let mut pinf = bnd_tool.begin_param(nbseg);
            let mut psup = bnd_tool.end_param(nbseg);
            if (psup - pinf) < 1e-10 {
                pinf -= 1e-10;
                psup += 1e-10;
            }
            let polygon = ThePolygonOfHInter::new_tool_range::<C, CT>(the_curve, pinf, psup, 2);
            host.internal_perform(the_curve, &polygon, the_surface, &polyhedron, u1new, v1new, u2new, v2new);
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformConicSurfCircle (pxx L711-759).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_conic_surf_circle<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_circle: &Circle3,
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    use crate::geomalgo::int_surf::quadric::Quadric;

    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    match a_surface_type {
        SurfaceType::Plane => {
            let circ_plane = IntConicQuad::new_circle_plane(
                the_circle,
                &<ST as HSurfaceTool>::plane(the_surface),
                IMPL_THE_TOLERANCE_ANGULAIRE,
                IMPL_THE_TOLERANCE,
            );
            host.append_int_ana(the_curve, the_surface, &circ_plane);
        }
        SurfaceType::Cylinder => {
            let circ_cylinder =
                IntConicQuad::new_circle_quadric(the_circle, &Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)));
            host.append_int_ana(the_curve, the_surface, &circ_cylinder);
        }
        SurfaceType::Cone => {
            let circ_cone = IntConicQuad::new_circle_quadric(the_circle, &Quadric::from_cone(&<ST as HSurfaceTool>::cone(the_surface)));
            host.append_int_ana(the_curve, the_surface, &circ_cone);
        }
        SurfaceType::Sphere => {
            let circ_sphere =
                IntConicQuad::new_circle_quadric(the_circle, &Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)));
            host.append_int_ana(the_curve, the_surface, &circ_sphere);
        }
        _ => {
            let polygon = ThePolygonOfHInter::new_tool::<C, CT>(the_curve, THE_NBSAMPLESONCIRCLE);
            host.internal_perform_bounds(the_curve, &polygon, the_surface, the_u1, the_v1, the_u2, the_v2);
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformConicSurfEllipse (pxx L762-810).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_conic_surf_ellipse<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_ellipse: &Ellipse3,
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    use crate::geomalgo::int_surf::quadric::Quadric;

    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    match a_surface_type {
        SurfaceType::Plane => {
            // OCCT Perform(const gp_Elips&, const gp_Pln&, ...) delegates to
            // the quadric path (Inter.pxx; IntConicQuad.cxx L562-565).
            let ellipse_plane = IntConicQuad::new_ellipse_plane(the_ellipse, &<ST as HSurfaceTool>::plane(the_surface));
            host.append_int_ana(the_curve, the_surface, &ellipse_plane);
        }
        SurfaceType::Cylinder => {
            let ellipse_cylinder =
                IntConicQuad::new_ellipse_quadric(the_ellipse, &Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)));
            host.append_int_ana(the_curve, the_surface, &ellipse_cylinder);
        }
        SurfaceType::Cone => {
            let ellipse_cone =
                IntConicQuad::new_ellipse_quadric(the_ellipse, &Quadric::from_cone(&<ST as HSurfaceTool>::cone(the_surface)));
            host.append_int_ana(the_curve, the_surface, &ellipse_cone);
        }
        SurfaceType::Sphere => {
            let ellipse_sphere =
                IntConicQuad::new_ellipse_quadric(the_ellipse, &Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)));
            host.append_int_ana(the_curve, the_surface, &ellipse_sphere);
        }
        _ => {
            let polygon = ThePolygonOfHInter::new_tool::<C, CT>(the_curve, THE_NBSAMPLESONELLIPSE);
            host.internal_perform_bounds(the_curve, &polygon, the_surface, the_u1, the_v1, the_u2, the_v2);
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformConicSurfParabola (pxx L813-888).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_conic_surf_parabola<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_parab: &Parabola3,
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    use crate::geomalgo::int_surf::quadric::Quadric;

    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    match a_surface_type {
        SurfaceType::Plane => {
            // OCCT Perform(const gp_Parab&, const gp_Pln&, Tolang) delegates
            // to the quadric path (IntConicQuad.cxx L567-570).
            let parab_plane = IntConicQuad::new_parabola_plane(the_parab, &<ST as HSurfaceTool>::plane(the_surface));
            host.append_int_ana(the_curve, the_surface, &parab_plane);
        }
        SurfaceType::Cylinder => {
            let parab_cylinder =
                IntConicQuad::new_parabola_quadric(the_parab, &Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)));
            host.append_int_ana(the_curve, the_surface, &parab_cylinder);
        }
        SurfaceType::Cone => {
            let parab_cone =
                IntConicQuad::new_parabola_quadric(the_parab, &Quadric::from_cone(&<ST as HSurfaceTool>::cone(the_surface)));
            host.append_int_ana(the_curve, the_surface, &parab_cone);
        }
        SurfaceType::Sphere => {
            let parab_sphere =
                IntConicQuad::new_parabola_quadric(the_parab, &Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)));
            host.append_int_ana(the_curve, the_surface, &parab_sphere);
        }
        _ => {
            let mut nbsu = <ST as HSurfaceTool>::nb_samples_u(the_surface, the_u1, the_u2);
            let mut nbsv = <ST as HSurfaceTool>::nb_samples_v(the_surface, the_v1, the_v2);
            if nbsu > 40 {
                nbsu = 40;
            }
            if nbsv > 40 {
                nbsv = 40;
            }
            let polyhedron =
                ThePolyhedronOfHInter::new_tool::<S, ST>(the_surface, nbsu, nbsv, the_u1, the_v1, the_u2, the_v2);
            let mut bnd_tool = crate::geomalgo::intf::IntfTool::new();
            let mut box_parab = BndBox::new();
            bnd_tool.parab_box(the_parab, polyhedron.bounding(), &mut box_parab);
            for nbseg in 1..=bnd_tool.nb_segments() {
                let polygon = ThePolygonOfHInter::new_tool_range::<C, CT>(
                    the_curve,
                    bnd_tool.begin_param(nbseg),
                    bnd_tool.end_param(nbseg),
                    THE_NBSAMPLESONPARAB,
                );
                host.internal_perform(the_curve, &polygon, the_surface, &polyhedron, the_u1, the_v1, the_u2, the_v2);
            }
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::PerformConicSurfHyperbola (pxx L891-966).
#[allow(clippy::too_many_arguments)]
pub(crate) fn perform_conic_surf_hyperbola<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_hypr: &Hyperbola3,
    the_curve: &C,
    the_surface: &S,
    the_u1: f64,
    the_v1: f64,
    the_u2: f64,
    the_v2: f64,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    use crate::geomalgo::int_surf::quadric::Quadric;

    let a_surface_type = <ST as HSurfaceTool>::get_type(the_surface);
    match a_surface_type {
        SurfaceType::Plane => {
            let hypr_plane = IntConicQuad::new_hyperbola_plane(the_hypr, &<ST as HSurfaceTool>::plane(the_surface));
            host.append_int_ana(the_curve, the_surface, &hypr_plane);
        }
        SurfaceType::Cylinder => {
            let hypr_cylinder =
                IntConicQuad::new_hyperbola_quadric(the_hypr, &Quadric::from_cylinder(&<ST as HSurfaceTool>::cylinder(the_surface)));
            host.append_int_ana(the_curve, the_surface, &hypr_cylinder);
        }
        SurfaceType::Cone => {
            let hypr_cone =
                IntConicQuad::new_hyperbola_quadric(the_hypr, &Quadric::from_cone(&<ST as HSurfaceTool>::cone(the_surface)));
            host.append_int_ana(the_curve, the_surface, &hypr_cone);
        }
        SurfaceType::Sphere => {
            let hypr_sphere =
                IntConicQuad::new_hyperbola_quadric(the_hypr, &Quadric::from_sphere(&<ST as HSurfaceTool>::sphere(the_surface)));
            host.append_int_ana(the_curve, the_surface, &hypr_sphere);
        }
        _ => {
            let mut nbsu = <ST as HSurfaceTool>::nb_samples_u(the_surface, the_u1, the_u2);
            let mut nbsv = <ST as HSurfaceTool>::nb_samples_v(the_surface, the_v1, the_v2);
            if nbsu > 40 {
                nbsu = 40;
            }
            if nbsv > 40 {
                nbsv = 40;
            }
            let polyhedron =
                ThePolyhedronOfHInter::new_tool::<S, ST>(the_surface, nbsu, nbsv, the_u1, the_v1, the_u2, the_v2);
            let mut bnd_tool = crate::geomalgo::intf::IntfTool::new();
            let mut box_hypr = BndBox::new();
            bnd_tool.hypr_box(the_hypr, polyhedron.bounding(), &mut box_hypr);
            for nbseg in 1..=bnd_tool.nb_segments() {
                let polygon = ThePolygonOfHInter::new_tool_range::<C, CT>(
                    the_curve,
                    bnd_tool.begin_param(nbseg),
                    bnd_tool.end_param(nbseg),
                    THE_NBSAMPLESONHYPR,
                );
                host.internal_perform(the_curve, &polygon, the_surface, &polyhedron, the_u1, the_v1, the_u2, the_v2);
            }
        }
    }
}

/// OCCT IntCurveSurface_InterImpl::AppendIntAna (pxx L969-1002).
pub(crate) fn append_int_ana<C: ?Sized, CT, S: ?Sized, ST, H>(
    the_curve: &C,
    the_surface: &S,
    the_int_ana: &IntConicQuad,
    host: &mut H,
) where
    CT: HCurveTool<Curve = C> + CurveTool3d<Curve = C>,
    ST: HSurfaceTool<Surface = S> + PSurfaceTool<Surface = S>,
    H: HInterHost<C, CT, S, ST> + ?Sized,
{
    let mut a_is_parallel = false;
    let mut a_points: Vec<IntersectionPoint> = Vec::new();

    if process_int_ana::<C, CT, S, ST>(the_curve, the_surface, the_int_ana, &mut a_is_parallel, &mut a_points) {
        if a_is_parallel {
            *host.is_parallel_flag() = true;
        } else {
            for a_point in &a_points {
                host.append(a_point);
            }
        }
    }
}
