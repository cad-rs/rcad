//! OCCT GeomFill_Darboux (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Darboux.hxx + GeomFill_Darboux.cxx (whole file L38-561,
//! including the file-local FDeriv / DDeriv / NormalD0 / NormalD1 /
//! NormalD2 statics).
//!
//! Architecture difference: OCCT stores the
//! `Handle(Adaptor3d_CurveOnSurface)` in the inherited myTrimmed member and
//! down-casts it in D0/D1/D2 (`GetCurve()` / `GetSurface()`).  Rust has no
//! handle RTTI, so the pair is carried explicitly in
//! [`Darboux::my_curve_on_surface`] and installed through
//! [`Darboux::set_curve_on_surface`] (the OCCT `SetCurve(COS)` call form).
//! OCCT `Adaptor3d_Surface::DN` is analytic per surface type while rcad
//! `SurfaceEval` exposes D1/D2 only: the higher orders use the
//! finite-difference fallback of [`surface_dn`] (same precedent as the
//! `CurveEval::derivative_at` / `derivatives2` defaults).

use glam::DVec3;

use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval};
use rcad_kernel::math::cs_lib::{
    dnnormal, dnnuv_array2, normal_from_derivatives_mag, normal_max_order, NormalStatus,
};
use rcad_kernel::math::GeomAbsShape;

use super::trihedron_law::{TrihedronLaw, TrihedronLawBase};

/// OCCT GeomAbs_C0 (default continuity when rcad carries no continuity data).
const CONT_C0: u8 = 0;
/// OCCT GeomAbs_C1.
const CONT_C1: u8 = 1;

/// OCCT `Handle(Adaptor3d_CurveOnSurface)` data pair — the downcast target
/// of the inherited myTrimmed member (architecture difference, see module
/// doc).
#[derive(Debug, Clone)]
pub struct CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface::GetCurve().
    pub curve2d: Curve2d,
    /// OCCT Adaptor3d_CurveOnSurface::GetSurface().
    pub surface: Surface3,
}

/// OCCT GeomFill_Darboux.
#[derive(Debug, Clone, Default)]
pub struct Darboux {
    /// OCCT GeomFill_TrihedronLaw protected base (myCurve / myTrimmed).
    pub(crate) base: TrihedronLawBase,
    /// OCCT myTrimmed down-cast to Adaptor3d_CurveOnSurface.
    my_curve_on_surface: Option<CurveOnSurface>,
}

/// OCCT static FDeriv (L38-47) — computes (F/|F|)'.
fn f_deriv(f: DVec3, df: DVec3) -> DVec3 {
    let norma = f.length();
    (df - f * (f.dot(df)) / (norma * norma)) / norma
}

/// OCCT static DDeriv (L49-60) — computes (F/|F|)''.
fn d_deriv(f: DVec3, df: DVec3, d2f: DVec3) -> DVec3 {
    let norma = f.length();
    (d2f - 2.0 * df * (f.dot(df)) / (norma * norma)) / norma
        - f * ((df.dot(df) + f.dot(d2f) - 3.0 * (f.dot(df)) * (f.dot(df)) / (norma * norma))
            / (norma * norma * norma))
}

/// OCCT Adaptor3d_Surface::DN(U, V, Nu, Nv) finite-difference fallback for
/// the orders rcad `SurfaceEval` does not host analytically (architecture
/// difference, see module doc).
fn surface_dn(surf: &Surface3, u: f64, v: f64, nu: usize, nv: usize) -> DVec3 {
    const H: f64 = 1.0e-5;
    if nu == 0 && nv == 0 {
        return surf.point_at(u, v);
    }
    if nv > 0 {
        return (surface_dn(surf, u, v + H, nu, nv - 1) - surface_dn(surf, u, v - H, nu, nv - 1))
            / (2.0 * H);
    }
    (surface_dn(surf, u + H, v, nu - 1, nv) - surface_dn(surf, u - H, v, nu - 1, nv)) / (2.0 * H)
}

/// OCCT static NormalD0 (L62-152) — computes the normal to the surface at
/// (U, V) (with the high-order fallback on singular points).
fn normal_d0(u: f64, v: f64, surf: &Surface3, normal: &mut DVec3, order_u: &mut i32, order_v: &mut i32) {
    //  gp_Vec D1U,D1V,D2U,D2V,DUV;
    let d1u: DVec3;
    let d1v: DVec3;
    // OCCT: Cont = min(UContinuity, VContinuity).  Architecture difference:
    // rcad Surface3 carries no continuity rank; C0 is the conservative
    // default (preserves the OCCT throw on the undefined-normal path).
    let cont = CONT_C0;
    *order_u = 0;
    *order_v = 0;
    // OCCT: gp_Pnt P; Surf->D1(U, V, P, D1U, D1V).
    let (_p, du, dv) = surf.derivatives(u, v);
    d1u = du;
    d1v = dv;
    let mag_tol = 0.000000001;
    let (a_normal, mut n_status) = normal_from_derivatives_mag(d1u, d1v, mag_tol);
    if let Some(n) = a_normal {
        *normal = n;
    }

    if n_status != NormalStatus::Defined {
        if cont == CONT_C0 || cont == CONT_C1 {
            panic!("Geom_UndefinedValue: GeomFill_Darboux::NormalD0");
        }
        let max_order = 3;
        let mut der_nuv = vec![vec![DVec3::ZERO; (max_order + 1) as usize]; (max_order + 1) as usize];
        let mut der_surf = vec![
            vec![DVec3::ZERO; (max_order + 2) as usize];
            (max_order + 2) as usize
        ];
        // OCCT: Umin = Surf->FirstUParameter(); ... Vmax = Surf->LastVParameter().
        let domain = surf.default_domain();
        let (umin, umax, vmin, vmax) = (domain[0], domain[1], domain[2], domain[3]);
        for i in 1..=(max_order + 1) {
            der_surf[i as usize][0] = surface_dn(surf, u, v, i as usize, 0);
        }

        for i in 0..=(max_order + 1) {
            for j in 1..=(max_order + 1) {
                der_surf[i as usize][j as usize] = surface_dn(surf, u, v, i as usize, j as usize);
            }
        }

        for i in 0..=max_order {
            for j in 0..=max_order {
                der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, &der_surf);
            }
        }

        let (status, a_normal, ou, ov) = normal_max_order(
            max_order,
            &der_nuv,
            mag_tol,
            u,
            v,
            umin,
            umax,
            vmin,
            vmax,
        );
        n_status = status;
        if let Some(n) = a_normal {
            *normal = n;
        }
        *order_u = ou;
        *order_v = ov;

        if n_status != NormalStatus::Defined {
            panic!("Geom_UndefinedValue: GeomFill_Darboux::NormalD0");
        }
    }
}

/// OCCT static NormalD1 (L154-247) — computes the normal to the surface and
/// its first derivative.
#[allow(clippy::too_many_arguments)]
fn normal_d1(
    u: f64,
    v: f64,
    surf: &Surface3,
    normal: &mut DVec3,
    d1u_normal: &mut DVec3,
    d1v_normal: &mut DVec3,
) {
    // OCCT: gp_Vec d2u, d2v, d2uv; gp_Pnt P;
    // Surf->D2(U, V, P, D1UNormal, D1VNormal, d2u, d2v, d2uv).
    let (_p, du, dv, d2u, d2uv, d2v) = surf.derivatives2(u, v);
    *d1u_normal = du;
    *d1v_normal = dv;
    let mag_tol = 0.000000001;
    let (a_normal, mut n_status) = normal_from_derivatives_mag(du, dv, mag_tol);
    if let Some(n) = a_normal {
        *normal = n;
    }
    let max_order = if n_status == NormalStatus::Defined { 0 } else { 3 };

    let mut der_nuv = vec![
        vec![DVec3::ZERO; (max_order + 2) as usize];
        (max_order + 2) as usize
    ];
    let mut der_surf = vec![
        vec![DVec3::ZERO; (max_order + 3) as usize];
        (max_order + 3) as usize
    ];
    // OCCT: Umin = Surf->FirstUParameter(); ... Vmax = Surf->LastVParameter().
    let domain = surf.default_domain();
    let (umin, umax, vmin, vmax) = (domain[0], domain[1], domain[2], domain[3]);

    der_surf[1][0] = *d1u_normal;
    der_surf[0][1] = *d1v_normal;
    der_surf[1][1] = d2uv;
    der_surf[2][0] = d2u;
    der_surf[0][2] = d2v;
    for i in 0..=(max_order + 1) {
        for j in i..=(max_order + 2) {
            if i + j > 2 {
                der_surf[i as usize][j as usize] = surface_dn(surf, u, v, i as usize, j as usize);
                if i != j {
                    der_surf[j as usize][i as usize] =
                        surface_dn(surf, u, v, j as usize, i as usize);
                }
            }
        }
    }

    for i in 0..=(max_order + 1) {
        for j in 0..=(max_order + 1) {
            der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, &der_surf);
        }
    }

    let (status, a_normal, order_u, order_v) =
        normal_max_order(max_order, &der_nuv, mag_tol, u, v, umin, umax, vmin, vmax);
    n_status = status;
    if let Some(n) = a_normal {
        *normal = n;
    }
    if n_status != NormalStatus::Defined {
        panic!("Geom_UndefinedValue: GeomFill_Darboux::NormalD1");
    }

    *d1u_normal = dnnormal(1, 0, &der_nuv, order_u, order_v);
    *d1v_normal = dnnormal(0, 1, &der_nuv, order_u, order_v);
}

/// OCCT static NormalD2 (L249-349) — computes the normal to the surface and
/// its first and second derivatives.
#[allow(clippy::too_many_arguments)]
fn normal_d2(
    u: f64,
    v: f64,
    surf: &Surface3,
    normal: &mut DVec3,
    d1u_normal: &mut DVec3,
    d1v_normal: &mut DVec3,
    d2u_normal: &mut DVec3,
    d2v_normal: &mut DVec3,
    d2uv_normal: &mut DVec3,
) {
    // OCCT: gp_Vec d3u, d3uuv, d3uvv, d3v; gp_Pnt P;
    // Surf->D3(U, V, P, D1UNormal, D1VNormal, D2UNormal, D2VNormal,
    //          D2UVNormal, d3u, d3v, d3uuv, d3uvv).
    let (_p, du, dv, d2u, d2v, d2uv) = surf.derivatives2(u, v);
    *d1u_normal = du;
    *d1v_normal = dv;
    *d2u_normal = d2u;
    *d2v_normal = d2v;
    *d2uv_normal = d2uv;
    let d3u = surface_dn(surf, u, v, 3, 0);
    let d3v = surface_dn(surf, u, v, 0, 3);
    let d3uuv = surface_dn(surf, u, v, 2, 1);
    let d3uvv = surface_dn(surf, u, v, 1, 2);
    let mag_tol = 0.000000001;
    let (a_normal, mut n_status) = normal_from_derivatives_mag(du, dv, mag_tol);
    if let Some(n) = a_normal {
        *normal = n;
    }
    let max_order = if n_status == NormalStatus::Defined { 0 } else { 3 };

    let mut der_nuv = vec![
        vec![DVec3::ZERO; (max_order + 3) as usize];
        (max_order + 3) as usize
    ];
    let mut der_surf = vec![
        vec![DVec3::ZERO; (max_order + 4) as usize];
        (max_order + 4) as usize
    ];

    // OCCT: Umin = Surf->FirstUParameter(); ... Vmax = Surf->LastVParameter().
    let domain = surf.default_domain();
    let (umin, umax, vmin, vmax) = (domain[0], domain[1], domain[2], domain[3]);

    der_surf[1][0] = *d1u_normal;
    der_surf[0][1] = *d1v_normal;
    der_surf[1][1] = *d2uv_normal;
    der_surf[2][0] = *d2u_normal;
    der_surf[0][2] = *d2v_normal;
    der_surf[3][0] = d3u;
    der_surf[2][1] = d3uuv;
    der_surf[1][2] = d3uvv;
    der_surf[0][3] = d3v;
    for i in 0..=(max_order + 2) {
        for j in i..=(max_order + 3) {
            if i + j > 3 {
                der_surf[i as usize][j as usize] = surface_dn(surf, u, v, i as usize, j as usize);
                if i != j {
                    der_surf[j as usize][i as usize] =
                        surface_dn(surf, u, v, j as usize, i as usize);
                }
            }
        }
    }

    for i in 0..=(max_order + 2) {
        for j in 0..=(max_order + 2) {
            der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, &der_surf);
        }
    }

    let (status, a_normal, order_u, order_v) =
        normal_max_order(max_order, &der_nuv, mag_tol, u, v, umin, umax, vmin, vmax);
    n_status = status;
    if let Some(n) = a_normal {
        *normal = n;
    }
    if n_status != NormalStatus::Defined {
        panic!("Geom_UndefinedValue: GeomFill_Darboux::NormalD2");
    }

    *d1u_normal = dnnormal(1, 0, &der_nuv, order_u, order_v);
    *d1v_normal = dnnormal(0, 1, &der_nuv, order_u, order_v);
    *d2u_normal = dnnormal(2, 0, &der_nuv, order_u, order_v);
    *d2v_normal = dnnormal(0, 2, &der_nuv, order_u, order_v);
    *d2uv_normal = dnnormal(1, 1, &der_nuv, order_u, order_v);
}

/// OCCT Adaptor2d_Curve2d::FirstParameter on the pcurve (architecture
/// helper for the CurveOnSurface pair — pure math/tool function).
pub(crate) fn curve2d_first_parameter(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Trimmed(tc) => tc.t_min,
        other => other.default_domain()[0],
    }
}

/// OCCT Adaptor2d_Curve2d::LastParameter on the pcurve (architecture
/// helper for the CurveOnSurface pair — pure math/tool function).
pub(crate) fn curve2d_last_parameter(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Trimmed(tc) => tc.t_max,
        other => other.default_domain()[1],
    }
}

impl Darboux {
    /// OCCT GeomFill_Darboux::GeomFill_Darboux() (L351).
    pub fn new() -> Self {
        Darboux {
            base: TrihedronLawBase::default(),
            my_curve_on_surface: None,
        }
    }

    /// OCCT `theLaw->SetCurve(Handle(Adaptor3d_CurveOnSurface))` — installs
    /// the curve-on-surface pair consumed by D0/D1/D2 (architecture
    /// difference, see module doc).  The base myCurve/myTrimmed members are
    /// left unset because rcad has no Curve3 view of a pcurve pair.
    pub fn set_curve_on_surface(&mut self, curve_on_surface: CurveOnSurface) {
        self.my_curve_on_surface = Some(curve_on_surface);
    }
}

impl TrihedronLaw for Darboux {
    fn my_curve(&self) -> &Option<rcad_kernel::geom::Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<rcad_kernel::geom::Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: rcad_kernel::geom::Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<rcad_kernel::geom::Curve3>) {
        self.base.my_trimmed = c;
    }

    /// OCCT Copy (L353-361).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let mut copy = Darboux::new();
        copy.my_curve_on_surface = self.my_curve_on_surface.clone();
        if let Some(curve) = &self.base.my_curve {
            TrihedronLaw::set_curve(&mut copy, curve.clone());
        }
        Box::new(copy)
    }

    /// OCCT D0 (L363-392).
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        let mut order_u = 0i32;
        let mut order_v = 0i32;
        // OCCT: myCurve2d->D1(Param, C2d, D2d).
        let c2d = cos.curve2d.point_at(param);
        let d2d = cos.curve2d.derivative_at(param);

        //  Normal = dS_du.Crossed(dS_dv).Normalized();
        let mut normal_dir = DVec3::ZERO;
        normal_d0(c2d.x, c2d.y, &cos.surface, &mut normal_dir, &mut order_u, &mut order_v);
        *binormal = normal_dir;

        // OCCT: mySupport->D1(C2d.X(), C2d.Y(), S, dS_du, dS_dv).
        let (_s, ds_du, ds_dv) = cos.surface.derivatives(c2d.x, c2d.y);

        *tangent = d2d.x * ds_du + d2d.y * ds_dv;
        *tangent = tangent.normalize_or_zero();

        *normal = *binormal;
        *normal = normal.cross(*tangent);

        true
    }

    /// OCCT D1 (L394-434).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
    ) -> bool {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        // OCCT: myCurve2d->D2(Param, C2d, D2d, D2_2d).
        let c2d = cos.curve2d.point_at(param);
        let d2d = cos.curve2d.derivative_at(param);
        let d2_2d = cos.curve2d.derivative2_at(param);
        // OCCT: mySupport->D2(C2d.X(), C2d.Y(), S, dS_du, dS_dv, d2S_du,
        // d2S_dv, d2S_duv).
        let (_s, ds_du, ds_dv, d2s_du, d2s_duv, d2s_dv) = cos.surface.derivatives2(c2d.x, c2d.y);

        let f = d2d.x * ds_du + d2d.y * ds_dv;
        *tangent = f.normalize_or_zero();
        let df = d2_2d.x * ds_du + d2_2d.y * ds_dv + d2s_du * d2d.x * d2d.x
            + 2.0 * d2s_duv * d2d.x * d2d.y
            + d2s_dv * d2d.y * d2d.y;

        *dtangent = f_deriv(f, df);

        let mut normal_dir = DVec3::ZERO;
        let mut d1u_normal = DVec3::ZERO;
        let mut d1v_normal = DVec3::ZERO;
        normal_d1(
            c2d.x,
            c2d.y,
            &cos.surface,
            &mut normal_dir,
            &mut d1u_normal,
            &mut d1v_normal,
        );
        *binormal = normal_dir;
        *dbinormal = d1u_normal * d2d.x + d1v_normal * d2d.y;

        *normal = *binormal;
        *normal = normal.cross(*tangent);
        *dnormal = binormal.cross(*dtangent) + dbinormal.cross(*tangent);

        true
    }

    /// OCCT D2 (L436-512).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        d2tangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        d2normal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        d2binormal: &mut DVec3,
    ) -> bool {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        // OCCT: myCurve2d->D3(Param, C2d, D2d, D2_2d, D3_2d).
        let c2d = cos.curve2d.point_at(param);
        let d2d = cos.curve2d.derivative_at(param);
        let d2_2d = cos.curve2d.derivative2_at(param);
        let d3_2d = cos.curve2d.derivative3_at(param);
        // OCCT: mySupport->D3(C2d.X(), C2d.Y(), S, dS_du, dS_dv, d2S_du,
        // d2S_dv, d2S_duv, d3S_du, d3S_dv, d3S_duuv, d3S_duvv).
        let (_s, ds_du, ds_dv, d2s_du, d2s_dv, d2s_duv) = cos.surface.derivatives2(c2d.x, c2d.y);
        let d3s_du = surface_dn(&cos.surface, c2d.x, c2d.y, 3, 0);
        let d3s_dv = surface_dn(&cos.surface, c2d.x, c2d.y, 0, 3);
        let d3s_duuv = surface_dn(&cos.surface, c2d.x, c2d.y, 2, 1);
        let d3s_duvv = surface_dn(&cos.surface, c2d.x, c2d.y, 1, 2);

        let f = d2d.x * ds_du + d2d.y * ds_dv;
        *tangent = f.normalize_or_zero();
        let df = d2_2d.x * ds_du + d2_2d.y * ds_dv + d2s_du * d2d.x * d2d.x
            + 2.0 * d2s_duv * d2d.x * d2d.y
            + d2s_dv * d2d.y * d2d.y;

        let d2f = d3_2d.x * ds_du
            + d3_2d.y * ds_dv
            + d2_2d.x * (d2d.x * d2s_du + d2d.y * d2s_duv)
            + d2_2d.y * (d2d.x * d2s_duv + d2d.y * d2s_dv)
            + d2d.x * d2d.x * (d2d.x * d3s_du + d2d.y * d3s_duuv)
            + d2d.y * d2d.y * (d2d.x * d3s_duvv + d2d.y * d3s_dv)
            + 2.0
                * (d2d.x * d2_2d.x * d2s_du
                    + d2d.y * d2_2d.y * d2s_dv
                    + d2d.x * d2d.y * (d2d.x * d3s_duuv + d2d.y * d3s_duvv)
                    + d2s_duv * (d2_2d.x * d2d.y + d2d.x * d2_2d.y));

        *dtangent = f_deriv(f, df);
        *d2tangent = d_deriv(f, df, d2f);

        let mut normal_dir = DVec3::ZERO;
        let mut d1u_normal = DVec3::ZERO;
        let mut d1v_normal = DVec3::ZERO;
        let mut d2u_normal = DVec3::ZERO;
        let mut d2v_normal = DVec3::ZERO;
        let mut d2uv_normal = DVec3::ZERO;
        normal_d2(
            c2d.x,
            c2d.y,
            &cos.surface,
            &mut normal_dir,
            &mut d1u_normal,
            &mut d1v_normal,
            &mut d2u_normal,
            &mut d2v_normal,
            &mut d2uv_normal,
        );
        *binormal = normal_dir;
        *dbinormal = d1u_normal * d2d.x + d1v_normal * d2d.y;
        *d2binormal = d1u_normal * d2_2d.x
            + d1v_normal * d2_2d.y
            + d2u_normal * d2d.x * d2d.x
            + 2.0 * d2uv_normal * d2d.x * d2d.y
            + d2v_normal * d2d.y * d2d.y;

        *normal = *binormal;
        *normal = normal.cross(*tangent);
        *dnormal = binormal.cross(*dtangent) + dbinormal.cross(*tangent);
        *d2normal = binormal.cross(*d2tangent) + 2.0 * dbinormal.cross(*dtangent)
            + d2binormal.cross(*tangent);

        true
    }

    /// OCCT NbIntervals (L514-517) — myCurve->NbIntervals(S) on the curve
    /// on surface (architecture difference: served from the pcurve domain).
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        let _first = curve2d_first_parameter(&cos.curve2d);
        let _last = curve2d_last_parameter(&cos.curve2d);
        1
    }

    /// OCCT Intervals (L519-522).
    fn intervals(&self, t: &mut Vec<f64>, _s: GeomAbsShape) {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        t[0] = curve2d_first_parameter(&cos.curve2d);
        let n = t.len();
        t[n - 1] = curve2d_last_parameter(&cos.curve2d);
    }

    /// OCCT GetAverageLaw (L524-551).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        let cos = self.my_curve_on_surface.as_ref().expect("null myTrimmed (CurveOnSurface)");
        let num = 20; // order of digitalization
        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        *atangent = DVec3::ZERO;
        *anormal = DVec3::ZERO;
        *abinormal = DVec3::ZERO;
        let first = curve2d_first_parameter(&cos.curve2d);
        let last = curve2d_last_parameter(&cos.curve2d);
        let step = (last - first) / num as f64;
        let mut param;
        for i in 0..=num {
            param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            TrihedronLaw::d0(self, param, &mut t, &mut n, &mut bn);
            *atangent += t;
            *anormal += n;
            *abinormal += bn;
        }
        *atangent /= (num + 1) as f64;
        *anormal /= (num + 1) as f64;

        *atangent = atangent.normalize_or_zero();
        *abinormal = atangent.cross(*anormal).normalize_or_zero();
        *anormal = abinormal.cross(*atangent);
    }

    /// OCCT IsConstant (L553-556) — (myCurve->GetType() == GeomAbs_Line),
    /// served from the pcurve type (architecture difference).
    fn is_constant(&self) -> bool {
        matches!(
            self.my_curve_on_surface.as_ref().map(|c| &c.curve2d),
            Some(Curve2d::Line(_))
        )
    }

    /// OCCT IsOnlyBy3dCurve (L558-561).
    fn is_only_by3d_curve(&self) -> bool {
        false
    }
}
