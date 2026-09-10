// OCCT ProjLib_CompProjectedCurve.cxx continuation (the Init / Perform
// machinery, the file-local statics and the Extrema GAP carriers).
//
// Anchored to the same source file as the parent module
// (`proj_lib_h_comp_projected_curve.rs`): the statics d1 / d2 / d2CurvOnSurf
// (cxx L138-302), ExactBound (L306-420), DichExactBound (L424-471),
// InitialPoint (L475-543), Init (L658-1225), Perform (L1229-1408),
// UpdateTripleByTrapCriteria call sites, BuildCurveSplits (L2242-2273),
// SplitOnDirection (L2276-2308), FindSplitPoint (L2312-2390) and the
// SplitDS structure (L82-115).
//
// GAP carriers (staged): Extrema_ExtCS and Extrema_ExtCC are not translated
// — the carriers below model the not-done output, which preserves the OCCT
// control flow (the Init MaxDist early-exit check and the seam-split search
// are skipped exactly as OCCT skips them when the extrema find nothing).
// Extrema_ExtPS routes to the kernel Surface3 engine through the
// kernel_surface() bridge (the OCCT class is parameterized by the adaptor,
// the rcad engine by the kernel surface the adaptor wraps).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::extrema::{POnCurve, POnSurface};
use rcad_kernel::base::extrema::ExtPS;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_cc::ExtremaExtCC;
use rcad_kernel::base::extrema_ext_cs::ExtremaExtCS;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface, Curve2dHandle, CurveOnSurface,
    SurfaceHandle,
};
use rcad_kernel::base::proj_lib::prj_resolve::PrjResolve;
use rcad_kernel::core::precision::{p_confusion, CONFUSION};
use rcad_kernel::geom::{Curve2d, Line2d};
use rcad_kernel::math::GeomAbsShape;

use super::CompProjectedCurve;
use super::{FUNC_TOL, CurveHandleAlias};

// ---------------------------------------------------------------------------
// gp::Resolution() — the DBL_MIN-scale singularity bound the statics use.
// ---------------------------------------------------------------------------

/// OCCT gp::Resolution() (gp.cxx) — 2.2250738585072014e-308.
pub(crate) const GP_RESOLUTION: f64 = 2.2250738585072014e-308;

// ---------------------------------------------------------------------------
// The Extrema_ExtPS of BuildCurveSplits (cxx L2250-2259)
// ---------------------------------------------------------------------------

/// The reusable Extrema_ExtPS of BuildCurveSplits (cxx L2250-2259: Initialize
/// + SetFlag(MIN)); each Perform(point) routes to the kernel Surface3 engine
/// through the kernel_surface bridge.
pub(crate) struct ReusedExtPS {
    surface: SurfaceHandle,
    uinf: f64,
    usup: f64,
    vinf: f64,
    vsup: f64,
    tol_u: f64,
    tol_v: f64,
    /// The kernel engine result of the last Perform.
    last: Option<ExtPS>,
}

impl ReusedExtPS {
    /// OCCT anExtPS.Initialize(S, FirstU, LastU, FirstV, LastV, TolU, TolV)
    /// + SetFlag(Extrema_ExtFlag_MIN).
    pub(crate) fn initialize(
        surface: SurfaceHandle,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        ReusedExtPS {
            surface,
            uinf,
            usup,
            vinf,
            vsup,
            tol_u,
            tol_v,
            last: None,
        }
    }

    /// OCCT Perform(P).
    pub(crate) fn perform(&mut self, point: DVec3) {
        // The kernel engine is Surface3-parameterized; the OCCT engine is
        // adaptor-parameterized (bridge as in contap).
        if let Some(kernel_surf) = self.surface.kernel_surface() {
            self.last = Some(ExtPS::with_domain(
                point,
                kernel_surf,
                self.uinf,
                self.usup,
                self.vinf,
                self.vsup,
                self.tol_u,
                self.tol_v,
            ));
        } else {
            self.last = None;
        }
    }

    /// OCCT IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.last.as_ref().map(|e| e.is_done()).unwrap_or(false)
    }

    /// OCCT NbExt().
    pub(crate) fn nb_ext(&self) -> usize {
        self.last.as_ref().map(|e| e.nb_ext()).unwrap_or(0)
    }

    /// OCCT SquareDistance(N).
    pub(crate) fn square_distance(&self, n: usize) -> f64 {
        self.last.as_ref().map(|e| e.square_distance(n)).unwrap_or(0.0)
    }

    /// OCCT Point(N).
    pub(crate) fn point(&self, n: usize) -> &POnSurface {
        self.last.as_ref().map(|e| e.point(n)).unwrap()
    }
}

// ---------------------------------------------------------------------------
// SplitDS (cxx L82-115)
// ---------------------------------------------------------------------------

/// OCCT SplitDS (cxx L82-115) — the split-point computation state.
pub(crate) struct SplitDS<'a> {
    /// OCCT: const handle(Adaptor3d_Curve) myCurve.
    my_curve: CurveHandleAlias,
    /// OCCT: const handle(Adaptor3d_Surface) mySurface.
    my_surface: SurfaceHandle,
    /// OCCT: NCollection_DynamicArray(double)& mySplits.
    my_splits: &'a mut Vec<f64>,
    my_per_min_param: f64,
    my_per_max_param: f64,
    /// 0 for U periodicity and 1 for V periodicity.
    my_periodic_dir: i32,
    my_ext_cc_curve1: Option<Arc<CurveOnSurface>>,
    my_ext_cc_last_2d_param: f64,
    my_ext_ps: Option<ReusedExtPS>,
}

impl<'a> SplitDS<'a> {
    pub(crate) fn new(
        the_curve: CurveHandleAlias,
        the_surface: SurfaceHandle,
        the_splits: &'a mut Vec<f64>,
    ) -> Self {
        SplitDS {
            my_curve: the_curve,
            my_surface: the_surface,
            my_splits: the_splits,
            my_per_min_param: 0.0,
            my_per_max_param: 0.0,
            my_periodic_dir: 0,
            my_ext_cc_curve1: None,
            my_ext_cc_last_2d_param: 0.0,
            my_ext_ps: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Statics: d1 / d2 / d2CurvOnSurf (cxx L138-302)
// ---------------------------------------------------------------------------

/// OCCT static d1 (cxx L138-164) — the (u, v) velocity of the projection.
pub(crate) fn d1(
    t: f64,
    u: f64,
    v: f64,
    v_out: &mut DVec2,
    curve: &dyn Adaptor3dCurve,
    surface: &dyn Adaptor3dSurface,
) {
    let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv) = surface.d2(u, v);
    let (c, dc1_t) = curve.d1(t);
    let ort = s - c; // OCCT gp_Vec Ort(C, S).

    let d_e_dt = DVec2::new(-dc1_t.dot(ds1_u), -dc1_t.dot(ds1_v));
    let d_e_du = DVec2::new(ds1_u.dot(ds1_u) + ort.dot(ds2_u), ds1_u.dot(ds1_v) + ort.dot(ds2_uv));
    let d_e_dv = DVec2::new(ds1_v.dot(ds1_u) + ort.dot(ds2_uv), ds1_v.dot(ds1_v) + ort.dot(ds2_v));

    let det = d_e_du.x * d_e_dv.y - d_e_du.y * d_e_dv.x;
    if det.abs() < GP_RESOLUTION {
        panic!("Standard_ConstructionError");
    }

    // OCCT gp_Mat2d M(gp_XY(dE_dv.Y()/det, -dE_du.Y()/det),
    //                 gp_XY(-dE_dv.X()/det, dE_du.X()/det)); rows below.
    let m_row1 = DVec2::new(d_e_dv.y / det, -d_e_du.y / det);
    let m_row2 = DVec2::new(-d_e_dv.x / det, d_e_du.x / det);

    // OCCT: V = -gp_Vec2d(gp_Vec2d(M.Row(1)) * dE_dt, gp_Vec2d(M.Row(2)) * dE_dt).
    *v_out = -DVec2::new(m_row1.dot(d_e_dt), m_row2.dot(d_e_dt));
}

/// OCCT static d2 (cxx L168-229) — the (u, v) velocity and acceleration.
pub(crate) fn d2(
    t: f64,
    u: f64,
    v: f64,
    v1_out: &mut DVec2,
    v2_out: &mut DVec2,
    curve: &dyn Adaptor3dCurve,
    surface: &dyn Adaptor3dSurface,
) {
    // OCCT L178: Surface->D3(u, v, S, DS1_u, DS1_v, DS2_u, DS2_v, DS2_uv,
    // DS3_u, DS3_v, DS3_uuv, DS3_uvv) — the D3 call provides every order.
    let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv, ds3_u, ds3_v, ds3_uuv, ds3_uvv) =
        surface.d3(u, v);
    let (_c, dc1_t, dc2_t) = curve.d2(t);
    let ort = s - c;

    let d_e_dt = DVec2::new(-dc1_t.dot(ds1_u), -dc1_t.dot(ds1_v));
    let d_e_du = DVec2::new(ds1_u.dot(ds1_u) + ort.dot(ds2_u), ds1_u.dot(ds1_v) + ort.dot(ds2_uv));
    let d_e_dv = DVec2::new(ds1_v.dot(ds1_u) + ort.dot(ds2_uv), ds1_v.dot(ds1_v) + ort.dot(ds2_v));

    let det = d_e_du.x * d_e_dv.y - d_e_du.y * d_e_dv.x;
    if det.abs() < GP_RESOLUTION {
        panic!("Standard_ConstructionError");
    }

    let m_row1 = DVec2::new(d_e_dv.y / det, -d_e_du.y / det);
    let m_row2 = DVec2::new(-d_e_dv.x / det, d_e_du.x / det);

    // First derivative.
    *v1_out = -DVec2::new(m_row1.dot(d_e_dt), m_row2.dot(d_e_dt));
    let v1 = *v1_out;

    // Second derivative.

    // Computation of d2E_dt2 = S1.
    let d2_e_dt = DVec2::new(-dc2_t.dot(ds1_u), -dc2_t.dot(ds1_v));

    // Computation of 2*(d2E/dtdX)(dX/dt) = S2.
    let d2_e1_dtd_x = DVec2::new(-dc1_t.dot(ds2_u), -dc1_t.dot(ds2_uv));
    let d2_e2_dtd_x = DVec2::new(-dc1_t.dot(ds2_uv), -dc1_t.dot(ds2_v));
    let s2 = 2.0 * DVec2::new(d2_e1_dtd_x.dot(v1), d2_e2_dtd_x.dot(v1));

    // Computation of (d2E/dX2)*(dX/dt)2 = S3.
    let tmp = 2.0 * ds1_u.dot(ds2_uv) + ds1_v.dot(ds2_u) + ort.dot(ds3_uuv);
    let row11 = DVec2::new(3.0 * ds1_u.dot(ds2_u) + ort.dot(ds3_u), tmp);
    let row12 = DVec2::new(tmp, ds2_v.dot(ds1_u) + 2.0 * ds1_v.dot(ds2_uv) + ort.dot(ds3_uvv));
    let tmp2 = 2.0 * ds2_uv.dot(ds1_v) + ds1_u.dot(ds2_v) + ort.dot(ds3_uvv);
    let row21 = DVec2::new(ds2_u.dot(ds1_v) + 2.0 * ds1_u.dot(ds2_uv) + ort.dot(ds3_uuv), tmp2);
    let row22 = DVec2::new(tmp2, 3.0 * ds1_v.dot(ds2_v) + ort.dot(ds3_v));

    let s3 = DVec2::new(
        v1.x * row11.dot(v1) + v1.y * row12.dot(v1),
        v1.x * row21.dot(v1) + v1.y * row22.dot(v1),
    );

    let sum = d2_e_dt + s2 + s3;

    *v2_out = -DVec2::new(m_row1.dot(sum), m_row2.dot(sum));
}

/// OCCT static d2CurvOnSurf (cxx L235-302) — the 3D velocity and
/// acceleration of the projected curve.
pub(crate) fn d2_curv_on_surf(
    t: f64,
    u: f64,
    v: f64,
    v1_out: &mut DVec3,
    v2_out: &mut DVec3,
    curve: &dyn Adaptor3dCurve,
    surface: &dyn Adaptor3dSurface,
) {
    // OCCT L246: Surface->D3(u, v, S, DS1_u, DS1_v, DS2_u, DS2_v, DS2_uv,
    // DS3_u, DS3_v, DS3_uuv, DS3_uvv) — the D3 call provides every order.
    let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv, ds3_u, ds3_v, ds3_uuv, ds3_uvv) =
        surface.d3(u, v);
    let (_c, dc1_t, dc2_t) = curve.d2(t);
    let ort = s - c;

    let d_e_dt = DVec2::new(-dc1_t.dot(ds1_u), -dc1_t.dot(ds1_v));
    let d_e_du = DVec2::new(ds1_u.dot(ds1_u) + ort.dot(ds2_u), ds1_u.dot(ds1_v) + ort.dot(ds2_uv));
    let d_e_dv = DVec2::new(ds1_v.dot(ds1_u) + ort.dot(ds2_uv), ds1_v.dot(ds1_v) + ort.dot(ds2_v));

    let det = d_e_du.x * d_e_dv.y - d_e_du.y * d_e_dv.x;
    if det.abs() < GP_RESOLUTION {
        panic!("Standard_ConstructionError");
    }

    let m_row1 = DVec2::new(d_e_dv.y / det, -d_e_du.y / det);
    let m_row2 = DVec2::new(-d_e_dv.x / det, d_e_du.x / det);

    // First derivative.
    let v12d = -DVec2::new(m_row1.dot(d_e_dt), m_row2.dot(d_e_dt));

    // Second derivative.
    let d2_e_dt = DVec2::new(-dc2_t.dot(ds1_u), -dc2_t.dot(ds1_v));
    let d2_e1_dtd_x = DVec2::new(-dc1_t.dot(ds2_u), -dc1_t.dot(ds2_uv));
    let d2_e2_dtd_x = DVec2::new(-dc1_t.dot(ds2_uv), -dc1_t.dot(ds2_v));
    let s2 = 2.0 * DVec2::new(d2_e1_dtd_x.dot(v12d), d2_e2_dtd_x.dot(v12d));

    let tmp = 2.0 * ds1_u.dot(ds2_uv) + ds1_v.dot(ds2_u) + ort.dot(ds3_uuv);
    let row11 = DVec2::new(3.0 * ds1_u.dot(ds2_u) + ort.dot(ds3_u), tmp);
    let row12 = DVec2::new(tmp, ds2_v.dot(ds1_u) + 2.0 * ds1_v.dot(ds2_uv) + ort.dot(ds3_uvv));
    let tmp2 = 2.0 * ds2_uv.dot(ds1_v) + ds1_u.dot(ds2_v) + ort.dot(ds3_uvv);
    let row21 = DVec2::new(ds2_u.dot(ds1_v) + 2.0 * ds1_u.dot(ds2_uv) + ort.dot(ds3_uuv), tmp2);
    let row22 = DVec2::new(tmp2, 3.0 * ds1_v.dot(ds2_v) + ort.dot(ds3_v));

    let s3 = DVec2::new(
        v12d.x * row11.dot(v12d) + v12d.y * row12.dot(v12d),
        v12d.x * row21.dot(v12d) + v12d.y * row22.dot(v12d),
    );

    let sum = d2_e_dt + s2 + s3;

    let v22d = -DVec2::new(m_row1.dot(sum), m_row2.dot(sum));

    // OCCT L299-301: the 3D composition.
    *v1_out = ds1_u * v12d.x + ds1_v * v12d.y;
    *v2_out = ds2_u * v12d.x * v12d.x
        + ds1_u * v22d.x
        + 2.0 * ds2_uv * v12d.x * v12d.y
        + ds2_v * v12d.y * v12d.y
        + ds1_v * v22d.y;
}

// ---------------------------------------------------------------------------
// ExactBound (cxx L306-420)
// ---------------------------------------------------------------------------

/// OCCT static ExactBound (cxx L306-420) — the exact boundary search.
/// `sol` carries the (t, u, v) triple in its X/Y/Z slots.
pub(crate) fn exact_bound(
    sol: &mut DVec3,
    not_sol: f64,
    tol: f64,
    tol_u: f64,
    tol_v: f64,
    curve: &dyn Adaptor3dCurve,
    surface: &dyn Adaptor3dSurface,
) -> bool {
    let u0 = sol.y;
    let v0 = sol.z;
    let first_u = surface.first_u_parameter();
    let last_u = surface.last_u_parameter();
    let first_v = surface.first_v_parameter();
    let last_v = surface.last_v_parameter();

    // Here we have to compute the boundary that projection is going to
    // intersect (cxx L322-349).
    let d2d = {
        let mut d = DVec2::ZERO;
        d1(sol.x, u0, v0, &mut d, curve, surface);
        d
    };
    // Here we assume that D2d != (0, 0).
    let (ru1, ru2, rv1, rv2);
    if d2d.x.abs() < GP_RESOLUTION {
        ru1 = f64::INFINITY; // OCCT Precision::Infinite()
        ru2 = f64::INFINITY;
        rv1 = v0 - first_v;
        rv2 = last_v - v0;
    } else if d2d.y.abs() < GP_RESOLUTION {
        ru1 = u0 - first_u;
        ru2 = last_u - u0;
        rv1 = f64::INFINITY;
        rv2 = f64::INFINITY;
    } else {
        ru1 = DVec2::new(u0, v0).distance(DVec2::new(first_u, v0 + (first_u - u0) * d2d.y / d2d.x));
        ru2 = DVec2::new(u0, v0).distance(DVec2::new(last_u, v0 + (last_u - u0) * d2d.y / d2d.x));
        rv1 = DVec2::new(u0, v0).distance(DVec2::new(u0 + (first_v - v0) * d2d.x / d2d.y, first_v));
        rv2 = DVec2::new(u0, v0).distance(DVec2::new(u0 + (last_v - v0) * d2d.x / d2d.y, last_v));
    }

    // OCCT L350-368: the descending sort by the Y slot.
    let mut seq: Vec<DVec3> = Vec::new();
    seq.push(DVec3::new(first_u, ru1, 2.0));
    seq.push(DVec3::new(last_u, ru2, 2.0));
    seq.push(DVec3::new(first_v, rv1, 3.0));
    seq.push(DVec3::new(last_v, rv2, 3.0));
    for i in 1..=3usize {
        for j in 0..(4 - i) {
            if seq[j].y < seq[j + 1].y {
                seq.swap(j, j + 1);
            }
        }
    }

    let t = sol.x;
    let t1 = sol.x.min(not_sol);
    let t2 = sol.x.max(not_sol);

    let mut is_done = false;
    while !seq.is_empty() {
        let p = *seq.last().unwrap();
        seq.pop();
        let fix = p.z as i32;
        let mut a_prj_ps = PrjResolve::new(curve, surface, fix);
        if fix == 2 {
            a_prj_ps.perform(
                t,
                p.x,
                v0,
                DVec2::new(tol, tol_v),
                DVec2::new(t1, surface.first_v_parameter()),
                DVec2::new(t2, surface.last_v_parameter()),
                FUNC_TOL,
                false,
            );
            if !a_prj_ps.is_done() {
                continue;
            }
            let p_ons = a_prj_ps.solution();
            *sol = DVec3::new(p_ons.x, p.x, p_ons.y);
            is_done = true;
            break;
        } else {
            a_prj_ps.perform(
                t,
                u0,
                p.x,
                DVec2::new(tol, tol_u),
                DVec2::new(t1, surface.first_u_parameter()),
                DVec2::new(t2, surface.last_u_parameter()),
                FUNC_TOL,
                false,
            );
            if !a_prj_ps.is_done() {
                continue;
            }
            let p_ons = a_prj_ps.solution();
            *sol = DVec3::new(p_ons.x, p_ons.y, p.x);
            is_done = true;
            break;
        }
    }

    is_done
}

// ---------------------------------------------------------------------------
// DichExactBound (cxx L424-471)
// ---------------------------------------------------------------------------

/// OCCT static DichExactBound (cxx L424-471) — the dichotomy boundary search.
pub(crate) fn dich_exact_bound(
    sol: &mut DVec3,
    not_sol: f64,
    tol: f64,
    tol_u: f64,
    tol_v: f64,
    curve: &dyn Adaptor3dCurve,
    surface: &dyn Adaptor3dSurface,
) {
    let mut u0 = sol.y;
    let mut v0 = sol.z;
    let mut a_prj_ps = PrjResolve::new(curve, surface, 1);

    let mut a_not_sol = not_sol;
    while (sol.x - a_not_sol).abs() > tol {
        let t = (sol.x + a_not_sol) / 2.0;
        a_prj_ps.perform(
            t,
            u0,
            v0,
            DVec2::new(tol_u, tol_v),
            DVec2::new(surface.first_u_parameter(), surface.first_v_parameter()),
            DVec2::new(surface.last_u_parameter(), surface.last_v_parameter()),
            FUNC_TOL,
            true,
        );

        if a_prj_ps.is_done() {
            let p_ons = a_prj_ps.solution();
            *sol = DVec3::new(t, p_ons.x, p_ons.y);
            u0 = sol.y;
            v0 = sol.z;
        } else {
            a_not_sol = t;
        }
    }
}

// ---------------------------------------------------------------------------
// InitialPoint (cxx L475-543)
// ---------------------------------------------------------------------------

/// OCCT static InitialPoint (cxx L475-543) — the projection start-point
/// search.  The Extrema_ExtPS (Extrema_ExtFlag_MIN) routes to the kernel
/// engine through the kernel_surface bridge.
#[allow(clippy::too_many_arguments)]
pub(crate) fn initial_point(
    point: DVec3,
    t: f64,
    c: &dyn Adaptor3dCurve,
    s: &dyn Adaptor3dSurface,
    tol_u: f64,
    tol_v: f64,
    u: &mut f64,
    v: &mut f64,
    the_max_dist: f64,
) -> bool {
    let mut a_prj_ps = PrjResolve::new(c, s, 1);

    // OCCT L488-496: Extrema_ExtPS(Point, *S, FirstU, LastU, FirstV, LastV,
    // TolU, TolV, Extrema_ExtFlag_MIN).
    let a_ext_ps = s.kernel_surface().map(|kernel_surf| {
        ExtPS::with_domain(
            point,
            kernel_surf,
            s.first_u_parameter(),
            s.last_u_parameter(),
            s.first_v_parameter(),
            s.last_v_parameter(),
            tol_u,
            tol_v,
        )
    });

    let mut argmin = 0usize;
    let mut a_max_dist = the_max_dist;
    if a_max_dist > 0.0 {
        a_max_dist *= a_max_dist;
    }
    if let Some(a_ext_ps) = &a_ext_ps {
        if a_ext_ps.is_done() && a_ext_ps.nb_ext() > 0 {
            // Search for the nearest solution which is also a normal
            // projection.
            let nend = a_ext_ps.nb_ext();
            for i in 1..=nend {
                if a_max_dist > 0.0 && a_max_dist < a_ext_ps.square_distance(i) {
                    continue;
                }
                let p_ons = a_ext_ps.point(i);
                let (par_u, par_v) = (p_ons.u, p_ons.v);
                a_prj_ps.perform(
                    t,
                    par_u,
                    par_v,
                    DVec2::new(tol_u, tol_v),
                    DVec2::new(s.first_u_parameter(), s.first_v_parameter()),
                    DVec2::new(s.last_u_parameter(), s.last_v_parameter()),
                    FUNC_TOL,
                    true,
                );
                if a_prj_ps.is_done()
                    && (argmin == 0 || a_ext_ps.square_distance(i) < a_ext_ps.square_distance(argmin))
                {
                    argmin = i;
                }
            }
        }
    }
    if argmin == 0 {
        false
    } else {
        let p_ons = a_ext_ps.as_ref().unwrap().point(argmin);
        *u = p_ons.u;
        *v = p_ons.v;
        true
    }
}

// ---------------------------------------------------------------------------
// BuildCurveSplits (cxx L2242-2273)
// ---------------------------------------------------------------------------

/// OCCT BuildCurveSplits (cxx L2242-2273) — the seam split points.
pub(crate) fn build_curve_splits(
    the_curve: CurveHandleAlias,
    the_surface: SurfaceHandle,
    the_tol_u: f64,
    the_tol_v: f64,
    the_splits: &mut Vec<f64>,
) {
    let mut a_ds = SplitDS::new(the_curve, the_surface, the_splits);

    let an_ext_ps = ReusedExtPS::initialize(
        a_ds.my_surface.clone(),
        a_ds.my_surface.first_u_parameter(),
        a_ds.my_surface.last_u_parameter(),
        a_ds.my_surface.first_v_parameter(),
        a_ds.my_surface.last_v_parameter(),
        the_tol_u,
        the_tol_v,
    );
    a_ds.my_ext_ps = Some(an_ext_ps);

    if a_ds.my_surface.is_u_periodic() {
        a_ds.my_periodic_dir = 0;
        split_on_direction(&mut a_ds);
    }
    if a_ds.my_surface.is_v_periodic() {
        a_ds.my_periodic_dir = 1;
        split_on_direction(&mut a_ds);
    }

    // OCCT: std::sort(..., Comparator) — ascending.
    the_splits.sort_by(|a, b| a.partial_cmp(b).unwrap());
}

// ---------------------------------------------------------------------------
// SplitOnDirection (cxx L2276-2308)
// ---------------------------------------------------------------------------

/// OCCT SplitOnDirection (cxx L2276-2308).
fn split_on_direction(the_split_ds: &mut SplitDS) {
    let surface = the_split_ds.my_surface.clone();
    let a_start_pnt = DVec2::new(surface.first_u_parameter(), surface.first_v_parameter());
    let a_dir = DVec2::new(
        the_split_ds.my_periodic_dir as f64,
        (1 - the_split_ds.my_periodic_dir) as f64,
    );

    the_split_ds.my_per_min_param = if the_split_ds.my_periodic_dir == 0 {
        surface.first_u_parameter()
    } else {
        surface.first_v_parameter()
    };
    the_split_ds.my_per_max_param = if the_split_ds.my_periodic_dir == 0 {
        surface.last_u_parameter()
    } else {
        surface.last_v_parameter()
    };
    let a_last_2d_param = if the_split_ds.my_periodic_dir != 0 {
        surface.last_u_parameter() - surface.first_u_parameter()
    } else {
        surface.last_v_parameter() - surface.first_v_parameter()
    };

    // Create line which is represent periodic border (cxx L2298-2303).  The
    // OCCT restricted Geom2dAdaptor_Curve(aC2GC, 0, aLast2DParam) is modeled
    // by the unrestricted adaptor (the range is carried by
    // myExtCCLast2DParam and handed to the extrema through SetRange).
    let a_c2gc = Curve2d::Line(Line2d {
        origin: a_start_pnt,
        direction: a_dir,
    });
    let a_c_ons = Arc::new(CurveOnSurface::new(
        super::geom2d_adaptor_curve(a_c2gc),
        the_split_ds.my_surface.clone(),
    ));
    the_split_ds.my_ext_cc_curve1 = Some(a_c_ons);
    the_split_ds.my_ext_cc_last_2d_param = a_last_2d_param;

    let (curve, _surface) = (the_split_ds.my_curve.clone(), ());
    let first_param = curve.first_parameter();
    let last_param = curve.last_parameter();
    find_split_point(the_split_ds, first_param, last_param);
}

// ---------------------------------------------------------------------------
// FindSplitPoint (cxx L2312-2390)
// ---------------------------------------------------------------------------

/// OCCT FindSplitPoint (cxx L2312-2390) — the recursive split search.
fn find_split_point(the_split_ds: &mut SplitDS, the_min_param: f64, the_max_param: f64) {
    // Make extrema copy to avoid dependencies between different levels of
    // the recursion.
    let curve1 = the_split_ds.my_ext_cc_curve1.clone().unwrap();
    let curve2 = the_split_ds.my_curve.clone();
    let a_tool1 = CurveToolHandle::other(curve1.as_ref());
    let a_tool2 = CurveToolHandle::other(curve2.as_ref());
    let mut an_ext_cc = ExtremaExtCC::new(1.0e-10, 1.0e-10);
    an_ext_cc.set_curve(1, &a_tool1);
    an_ext_cc.set_curve(2, &a_tool2);
    // Search only one solution since multiple invocations are needed.
    an_ext_cc.set_single_solution_flag(true);
    an_ext_cc.set_range(1, 0.0, the_split_ds.my_ext_cc_last_2d_param);
    an_ext_cc.set_range(2, the_min_param, the_max_param);
    an_ext_cc.perform();

    if an_ext_cc.is_done() && !an_ext_cc.is_parallel() {
        let a_nb_ext = an_ext_cc.nb_ext();
        for an_idx in 1..=a_nb_ext {
            let mut p_on_c1 = POnCurve {
                param: 0.0,
                point: DVec3::ZERO,
            };
            let mut p_on_c2 = POnCurve {
                param: 0.0,
                point: DVec3::ZERO,
            };
            an_ext_cc.points(an_idx, &mut p_on_c1, &mut p_on_c2);

            let ext_ps = the_split_ds.my_ext_ps.as_mut().unwrap();
            ext_ps.perform(p_on_c2.point);
            if !ext_ps.is_done() {
                return;
            }

            // Find point with the minimal Euclidean distance to avoid false
            // positive points detection.
            let mut a_min_idx: isize = -1;
            let mut a_min_sq_dist = f64::MAX; // OCCT RealLast()
            let a_nb_pext = ext_ps.nb_ext();
            for a_p_idx in 1..=a_nb_pext {
                let a_curr_sq_dist = ext_ps.square_distance(a_p_idx);
                if a_curr_sq_dist < a_min_sq_dist {
                    a_min_sq_dist = a_curr_sq_dist;
                    a_min_idx = a_p_idx as isize;
                }
            }

            // Check that is point will be projected to the periodic border.
            let a_p_ons_param = {
                let a_p_ons = ext_ps.point(a_min_idx as usize);
                let (u, v) = (a_p_ons.u, a_p_ons.v);
                if the_split_ds.my_periodic_dir != 0 {
                    v
                } else {
                    u
                }
            };

            if (a_p_ons_param - the_split_ds.my_per_min_param).abs() < p_confusion()
                || (a_p_ons_param - the_split_ds.my_per_max_param).abs() < p_confusion()
            {
                let a_param = p_on_c2.param;
                let a_cf_param = the_split_ds.my_curve.first_parameter();
                let a_cl_param = the_split_ds.my_curve.last_parameter();

                if a_param > a_cf_param + p_confusion() && a_param < a_cl_param - p_confusion() {
                    // Add only inner points.
                    the_split_ds.my_splits.push(a_param);
                }

                let a_delta_coeff = 0.01;
                let a_delta = (the_max_param - the_min_param + a_cl_param - a_cf_param) * a_delta_coeff;

                if a_param - a_delta > the_min_param + p_confusion() {
                    find_split_point(the_split_ds, the_min_param, a_param - a_delta);
                }

                if a_param + a_delta < the_max_param - p_confusion() {
                    find_split_point(the_split_ds, a_param + a_delta, the_max_param);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Init (cxx L658-1225)
// ---------------------------------------------------------------------------

/// OCCT ProjLib_CompProjectedCurve::Init (cxx L658-1225) — computes the set
/// of projected points and the continuous parts of the projected curves.
pub(crate) fn init_body(this: &mut CompProjectedCurve) {
    this.my_tab_int = None;
    let mut a_splits: Vec<f64> = Vec::new();

    let mut tol = 0.0f64; // Tolerance for ExactBound
    let mut nend = 0usize;
    let mut a_split_idx = 0usize;
    let mut from_last_u = false;
    let mut is_splits_computed = false;

    let a_tol_ext = p_confusion();
    let curve = this.my_curve.clone().expect("Init");
    let surface = this.my_surface.clone().expect("Init");
    let a_curve_tool = CurveToolHandle::other(curve.as_ref());
    let mut cext = ExtremaExtCS::new_curve_surface(
        &a_curve_tool,
        surface.as_ref(),
        a_tol_ext,
        a_tol_ext,
    );
    if cext.is_done() && cext.nb_ext() > 0 {
        // Search for the minimum solution.
        // Avoid usage of extrema result that can be wrong for extrusion.
        if this.my_max_dist > 0.0 && surface.get_type() != GeomAbsSurfaceType::SurfaceOfExtrusion {
            let mut min_val2 = cext.square_distance(1);

            nend = cext.nb_ext();
            for i in 2..=nend {
                if cext.square_distance(i) < min_val2 {
                    min_val2 = cext.square_distance(i);
                }
            }
            if min_val2 > this.my_max_dist * this.my_max_dist {
                return; // No near solution -> exit.
            }
        }
    }

    let first_u = curve.first_parameter();
    let last_u = curve.last_parameter();
    const GLOBAL_MIN_STEP: f64 = 1.0e-4;
    // <GlobalMinStep> is sufficiently small to provide solving from initial
    // point and, on the other hand, it is sufficiently large to avoid too
    // close solutions.
    let min_step = 0.01 * (last_u - first_u);
    let max_step = 0.1 * (last_u - first_u);
    let search_step = 10.0 * min_step;
    let mut step = search_step;

    let a_low_border = DVec2::new(surface.first_u_parameter(), surface.first_v_parameter());
    let a_upp_border = DVec2::new(surface.last_u_parameter(), surface.last_v_parameter());
    let a_tol = DVec2::new(this.my_tol_u, this.my_tol_v);
    let mut a_prj_ps = PrjResolve::new(curve.as_ref(), surface.as_ref(), 1);

    let mut t = first_u;
    let mut new_part;
    let mut prev_deb = 0.0f64;
    let mut same_deb = false;

    let mut triple = DVec3::ZERO;

    // Basic loop (cxx L720).
    while t <= last_u {
        // Search for the beginning of a new continuous part
        // to avoid infinite computation in some difficult cases.
        new_part = false;
        if t > first_u && (t - prev_deb).abs() <= p_confusion() {
            same_deb = true;
        }
        while t <= last_u && !new_part && !from_last_u && !same_deb {
            prev_deb = t;
            if t == last_u {
                from_last_u = true;
            }
            let mut initpoint = false;
            let mut u = 0.0f64;
            let mut v = 0.0f64;
            let mut c_point = DVec3::ZERO;

            // Search an initial point in the list of Extrema Curve-Surface
            // (cxx L742-766).
            if nend != 0 && !cext.is_parallel() {
                for i in 1..=nend {
                    let mut p1 = POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    let mut p2 = POnSurface {
                        u: 0.0,
                        v: 0.0,
                        point: DVec3::ZERO,
                    };
                    cext.points(i, &mut p1, &mut p2);
                    let par_t = p1.param;
                    let (par_u, par_v) = (p2.u, p2.v);

                    a_prj_ps.perform(par_t, par_u, par_v, a_tol, a_low_border, a_upp_border, FUNC_TOL, true);

                    if a_prj_ps.is_done()
                        && p1.param > first_u.max(t - step + p_confusion())
                        && p1.param <= t
                    {
                        t = par_t;
                        u = par_u;
                        v = par_v;
                        c_point = p1.point;
                        initpoint = true;
                        break;
                    }
                }
            }
            if !initpoint {
                c_point = curve.value(t);
                // PConfusion - use geometric tolerances in extrema /
                // optimization (cxx L773-774).
                initpoint = initial_point(
                    c_point,
                    t,
                    curve.as_ref(),
                    surface.as_ref(),
                    this.my_tol_u,
                    this.my_tol_v,
                    &mut u,
                    &mut v,
                    this.my_max_dist,
                );
            }
            if initpoint {
                // When U or V lie on surface joint in some cases we cannot
                // use them as initial point for aPrjPS, so we switch them
                // (cxx L786-830).
                if (surface.is_u_periodic()
                    && (a_upp_border.x - a_low_border.x - surface.u_period()).abs() < CONFUSION)
                    || (surface.is_v_periodic()
                        && (a_upp_border.y - a_low_border.y - surface.v_period()).abs() < CONFUSION)
                {
                    if (u - a_low_border.x).abs() < surface.u_resolution(p_confusion())
                        && surface.is_u_periodic()
                    {
                        let mut d = DVec2::ZERO;
                        d1(t, u, v, &mut d, curve.as_ref(), surface.as_ref());
                        if d.x < 0.0 {
                            u = a_upp_border.x;
                        }
                    } else if (u - a_upp_border.x).abs() < surface.u_resolution(p_confusion())
                        && surface.is_u_periodic()
                    {
                        let mut d = DVec2::ZERO;
                        d1(t, u, v, &mut d, curve.as_ref(), surface.as_ref());
                        if d.x > 0.0 {
                            u = a_low_border.x;
                        }
                    }

                    if (v - a_low_border.y).abs() < surface.v_resolution(p_confusion())
                        && surface.is_v_periodic()
                    {
                        let mut d = DVec2::ZERO;
                        d1(t, u, v, &mut d, curve.as_ref(), surface.as_ref());
                        if d.y < 0.0 {
                            v = a_upp_border.y;
                        }
                    } else if (v - a_upp_border.y).abs() <= surface.v_resolution(p_confusion())
                        && surface.is_v_periodic()
                    {
                        let mut d = DVec2::ZERO;
                        d1(t, u, v, &mut d, curve.as_ref(), surface.as_ref());
                        if d.y > 0.0 {
                            v = a_low_border.y;
                        }
                    }
                }

                if this.my_max_dist > 0.0 {
                    // Here we are going to stop if the distance between
                    // projection and corresponding curve point is greater
                    // than myMaxDist (cxx L832-846).
                    let p_ons = surface.value(u, v);
                    let d = c_point.distance(p_ons);
                    if d > this.my_max_dist {
                        this.my_sequence.as_mut().unwrap().clear();
                        this.my_nb_curves = 0;
                        return;
                    }
                }
                triple = DVec3::new(t, u, v);
                if t != first_u {
                    // Search for exact boundary point
