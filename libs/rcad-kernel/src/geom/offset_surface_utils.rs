//! OCCT Geom_OffsetSurfaceUtils (TKG3d/Geom/Geom_OffsetSurfaceUtils.pxx,
//! L41-1833) — 1:1 translation of the internal helper namespace of
//! Geom_OffsetSurface.  The rcad re-hosts of the `Geom_OffsetSurface` class
//! bodies that consume the namespace (the `SetBasisSurface` unwrap, `Surface()`
//! and `EvalD0`/`EvalD1`/`EvalD2`/`EvalDN`) live in
//! `offset_surface_utils_b.rs` and are re-exported from here.
//!
//! Architecture differences (rcad value model):
//!   - `Geom_Surface::ResD1/ResD2/ResD3` (Geom_Surface.hxx nested structs) map
//!     to the [`ResD1`] value struct here (the other rcad modules carry the
//!     same data as tuples);
//!   - `occ::handle(Geom_Surface)` maps to [`Surface3`]; the OCCT null handle
//!     of `aOscSurf` / `aDummy` maps to `None`;
//!   - `occ::handle(Geom_OsculatingSurface) myOscSurf` maps to
//!     `Option<OsculatingSurface>` — the rcad OffsetSurface payload does NOT
//!     carry the `myEvalRep` / `myOscSurf` members, so both are recomputed on
//!     demand ([`offset_equivalent_surface`] and
//!     [`offset_surface_osculating`]); OCCT's
//!     `ClearEvalRepresentation()` state (UReverse/VReverse/Transform) has no
//!     rcad counterpart, so the rcad offset surface behaves as the freshly
//!     constructed OCCT one (the eval representation present whenever
//!     `Surface()` yields a surface);
//!   - `NCollection_Array2<gp_Vec>` of the derivative tables maps to
//!     `Vec<Vec<DVec3>>` — the OCCT arrays are 0-based
//!     (`NCollection_Array2<gp_Vec>(buffer, 0, N, 0, M)`), so the rcad
//!     indices match the OCCT ones directly;
//!   - the OCCT `EvalD1/EvalD2/EvalD3/EvalDN` of the basis and of the
//!     osculating surface are the exact per-type kernel derivatives of
//!     [`eval_d1`] / [`eval_d2`] / [`eval_d3`] / [`eval_dn`]
//!     (`GeomAdaptor_Surface::DN` re-host: ElSLib for the quadrics,
//!     `Geom_BSplineSurface::DN` via BSplCLib::Eval for the polynomial kinds).

use glam::DVec3;

use super::osculating_surface::OsculatingSurface;
use crate::core::precision::CONFUSION;
use crate::geom::{BSplineSurface, Surface3, SurfaceEval};
use crate::math::cs_lib::{
    dnnormal, dnnuv_array2, normal_from_derivatives_mag, normal_max_order, NormalStatus,
};

/// OCCT `Geom_OffsetSurfaceUtils::THE_D1_MAGNITUDE_TOL` (pxx L45).
pub const THE_D1_MAGNITUDE_TOL: f64 = 1.0e-9;

/// OCCT `Geom_Surface::ResD1` (Geom_Surface.hxx) — point + first partials.
#[derive(Debug, Clone, Copy)]
pub struct ResD1 {
    pub point: DVec3,
    pub d1u: DVec3,
    pub d1v: DVec3,
}

/// OCCT `Geom_Surface::ResD2` (Geom_Surface.hxx) — point + partials up to 2nd
/// order.  The member order is the OCCT one (`Point, D1U, D1V, D2U, D2V,
/// D2UV`); the [`SurfaceEval::derivatives2`] tuple form used by the rcad
/// surfaces is the same data with `D2UV` and `D2V` swapped.
#[derive(Debug, Clone, Copy)]
pub struct ResD2 {
    pub point: DVec3,
    pub d1u: DVec3,
    pub d1v: DVec3,
    pub d2u: DVec3,
    pub d2v: DVec3,
    pub d2uv: DVec3,
}

/// OCCT `Geom_OffsetSurfaceUtils::OsculatingInfo` (pxx L49-60).
#[derive(Debug, Clone, Copy, Default)]
pub struct OsculatingInfo {
    /// True if osculating along U direction.
    pub along_u: bool,
    /// True if osculating along V direction.
    pub along_v: bool,
    /// True if normal direction should be reversed.
    pub is_opposite: bool,
}

impl OsculatingInfo {
    /// OCCT `Sign()` (pxx L56).
    pub fn sign(&self) -> f64 {
        if (self.along_u || self.along_v) && self.is_opposite {
            -1.0
        } else {
            1.0
        }
    }

    /// OCCT `HasOsculating()` (pxx L59).
    pub fn has_osculating(&self) -> bool {
        self.along_u || self.along_v
    }
}

// =========================================================================
// The Geom_Surface evaluation leaves (EvalD1 / EvalD2 / EvalD3 / EvalDN)
// =========================================================================

/// OCCT `Geom_Surface::EvalDN(U, V, Nu, Nv)` over the rcad surface value —
/// the exact per-type kernel derivative of the `GeomAdaptor_Surface::DN`
/// engine (ElSLib::DN for the quadrics, `BSplSLib::DN` for the polynomial
/// kinds, i.e. `Geom_BSplineSurface::EvalDN` and `Geom_BezierSurface::EvalDN`).
pub fn eval_dn(the_s: &Surface3, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    match the_s {
        Surface3::Plane(_)
        | Surface3::Cylinder(_)
        | Surface3::Cone(_)
        | Surface3::Sphere(_)
        | Surface3::Torus(_) => the_s.dn(u, v, nu, nv),
        Surface3::BSpline(bs) => bspl_slib_dn(bs, u, v, nu, nv),
        Surface3::Bezier(bez) => crate::geom::eval::bezier_surface_dn(bez, u, v, nu, nv),
        Surface3::Offset(of) => offset_payload_eval_dn(of, u, v, nu, nv),
        Surface3::LinearExtrusion(le) => {
            crate::geom::extrusion_utils::linear_extrusion_eval_dn(le, u, v, nu, nv)
        }
        Surface3::Revolution(rev) => {
            crate::geom::revolution_utils::revolution_eval_dn(rev, u, v, nu, nv)
        }
        Surface3::Ellipsoid(el) => crate::geom::eval_c::ellipsoid_eval_dn(el, u, v, nu, nv),
        Surface3::Helicoid(h) => crate::geom::eval_c::helicoid_eval_dn(h, u, v, nu, nv),
        _ => panic!(
            "GAP: Geom_Surface::EvalDN (TKG3d/Geom) is not translated for this surface type \
             (the rcad GeomAdaptor_Surface DN engine covers the ElSLib surfaces, \
             Geom_BSplineSurface, Geom_BezierSurface, Geom_OffsetSurface, \
             Geom_SurfaceOfLinearExtrusion, Geom_SurfaceOfRevolution, \
             GeomEval_EllipsoidSurface and GeomEval_CircularHelicoidSurface) — \
             Geom_OffsetSurfaceUtils::ComputeDerivatives"
        ),
    }
}

/// OCCT `Geom_BSplineSurface::EvalDN` (Geom_BSplineSurface_1.cxx L279-313) —
/// the (Nu, Nv) derivative of a BSpline surface at (U, V).
///
/// Delegates to the single faithful translation of its engine,
/// `BSplSLib::DN` (BSplSLib.cxx L1519-1605), in [`crate::geom::eval_b`], which
/// carries the whole OCCT chain (`PrepareEval` + the two `BSplCLib::Bohm`
/// passes + the `(n1 * (d2 + 1) + n2)` pole extraction, and
/// `BSplSLib::RationalDerivative` for the rational branch).  An earlier
/// revision of this file kept a second, rational-less copy of the same body
/// here.
fn bspl_slib_dn(bs: &BSplineSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    crate::geom::eval::bspline_surface_dn(bs, u, v, nu, nv)
}

/// OCCT `Geom_Surface::EvalD1(U, V)` over the rcad surface value.
pub fn eval_d1(the_s: &Surface3, u: f64, v: f64) -> ResD1 {
    let (point, d1u, d1v) = super::osculating_surface::surface_d1(the_s, u, v);
    ResD1 { point, d1u, d1v }
}

/// OCCT `Geom_Surface::EvalD2(U, V)` over the rcad surface value.
pub fn eval_d2(the_s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
    (
        SurfaceEval::point_at(the_s, u, v),
        eval_dn(the_s, u, v, 1, 0),
        eval_dn(the_s, u, v, 0, 1),
        eval_dn(the_s, u, v, 2, 0),
        eval_dn(the_s, u, v, 1, 1),
        eval_dn(the_s, u, v, 0, 2),
    )
}

/// OCCT `Geom_Surface::EvalD3(U, V)` over the rcad surface value.
#[allow(clippy::type_complexity)]
pub fn eval_d3(
    the_s: &Surface3,
    u: f64,
    v: f64,
) -> (
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
    DVec3,
) {
    let d2 = eval_d2(the_s, u, v);
    (
        d2.0,
        d2.1,
        d2.2,
        d2.3,
        d2.4,
        d2.5,
        eval_dn(the_s, u, v, 3, 0),
        eval_dn(the_s, u, v, 0, 3),
        eval_dn(the_s, u, v, 2, 1),
        eval_dn(the_s, u, v, 1, 2),
    )
}

// =========================================================================
// The vector helpers (pxx L62-210)
// =========================================================================

/// OCCT `IsInfiniteCoord` (pxx L65-69).
pub fn is_infinite_coord(the_vec: DVec3) -> bool {
    crate::core::precision::is_infinite_value(the_vec.x)
        || crate::core::precision::is_infinite_value(the_vec.y)
        || crate::core::precision::is_infinite_value(the_vec.z)
}

/// OCCT `IsSingular` (pxx L76-96).
pub fn is_singular(the_d1u: DVec3, the_d1v: DVec3, the_tol: f64) -> bool {
    let mut a_d1u = the_d1u;
    let mut a_d1v = the_d1v;
    let a_d1u_norm2 = a_d1u.length_squared();
    let a_d1v_norm2 = a_d1v.length_squared();
    if a_d1u_norm2 > 1.0 {
        a_d1u /= a_d1u_norm2.sqrt();
    }
    if a_d1v_norm2 > 1.0 {
        a_d1v /= a_d1v_norm2.sqrt();
    }
    let a_norm = a_d1u.cross(a_d1v);
    a_norm.length_squared() <= the_tol * the_tol
}

/// OCCT `ComputeNormal` (pxx L104-130) — `Some(normal)` is the OCCT `true`
/// return; `None` is the singular (`false`) return.
pub fn compute_normal(the_d1u: DVec3, the_d1v: DVec3, the_tol: f64) -> Option<DVec3> {
    let mut a_d1u = the_d1u;
    let mut a_d1v = the_d1v;
    let a_d1u_norm2 = a_d1u.length_squared();
    let a_d1v_norm2 = a_d1v.length_squared();
    if a_d1u_norm2 > 1.0 {
        a_d1u /= a_d1u_norm2.sqrt();
    }
    if a_d1v_norm2 > 1.0 {
        a_d1v /= a_d1v_norm2.sqrt();
    }
    let mut the_normal = a_d1u.cross(a_d1v);
    if the_normal.length_squared() <= the_tol * the_tol {
        return None;
    }
    the_normal = the_normal.normalize();
    Some(the_normal)
}

/// OCCT `ComputeDNormalU` (pxx L139-159).
pub fn compute_dnormal_u(
    the_d1u: DVec3,
    the_d1v: DVec3,
    the_d2u: DVec3,
    the_d2uv: DVec3,
    the_normal: DVec3,
) -> DVec3 {
    let a_scale = the_d1u.cross(the_d1v).dot(the_normal);
    let a_n1u = DVec3::new(
        the_d2u.y * the_d1v.z + the_d1u.y * the_d2uv.z - the_d2u.z * the_d1v.y
            - the_d1u.z * the_d2uv.y,
        -(the_d2u.x * the_d1v.z + the_d1u.x * the_d2uv.z - the_d2u.z * the_d1v.x
            - the_d1u.z * the_d2uv.x),
        the_d2u.x * the_d1v.y + the_d1u.x * the_d2uv.y - the_d2u.y * the_d1v.x
            - the_d1u.y * the_d2uv.x,
    );
    let a_scale_u = a_n1u.dot(the_normal);
    (a_n1u - a_scale_u * the_normal) / a_scale
}

/// OCCT `ComputeDNormalV` (pxx L168-188).
pub fn compute_dnormal_v(
    the_d1u: DVec3,
    the_d1v: DVec3,
    the_d2v: DVec3,
    the_d2uv: DVec3,
    the_normal: DVec3,
) -> DVec3 {
    let a_scale = the_d1u.cross(the_d1v).dot(the_normal);
    let a_n1v = DVec3::new(
        the_d2uv.y * the_d1v.z + the_d2v.z * the_d1u.y - the_d2uv.z * the_d1v.y
            - the_d2v.y * the_d1u.z,
        -(the_d2uv.x * the_d1v.z + the_d2v.z * the_d1u.x - the_d2uv.z * the_d1v.x
            - the_d2v.x * the_d1u.z),
        the_d2uv.x * the_d1v.y + the_d2v.y * the_d1u.x - the_d2uv.y * the_d1v.x
            - the_d2v.x * the_d1u.y,
    );
    let a_scale_v = a_n1v.dot(the_normal);
    (a_n1v - a_scale_v * the_normal) / a_scale
}

/// OCCT `CalculateD0` (pxx L197-210) — on `Some`, `the_value` holds the offset
/// point; `None` is the singular `false` return.
pub fn calculate_d0(the_value: DVec3, the_d1u: DVec3, the_d1v: DVec3, the_offset: f64) -> Option<DVec3> {
    calculate_d0_signed(the_value, the_d1u, the_d1v, the_offset, 1.0)
}

/// OCCT `CalculateD0` with the explicit `theSign` factor.
pub fn calculate_d0_signed(
    the_value: DVec3,
    the_d1u: DVec3,
    the_d1v: DVec3,
    the_offset: f64,
    the_sign: f64,
) -> Option<DVec3> {
    let a_norm = compute_normal(the_d1u, the_d1v, THE_D1_MAGNITUDE_TOL)?;
    Some(the_value + the_offset * the_sign * a_norm)
}

/// OCCT `CalculateD1` (pxx L222-248) — the free helper of the namespace (no
/// caller inside the pxx; the D1 body of `EvaluateD1` spells the same math
/// inline).  Returns `(point, d1u, d1v)` on success.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn calculate_d1(
    the_value: DVec3,
    the_d1u: DVec3,
    the_d1v: DVec3,
    the_d2u: DVec3,
    the_d2v: DVec3,
    the_d2uv: DVec3,
    the_offset: f64,
) -> Option<(DVec3, DVec3, DVec3)> {
    calculate_d1_signed(the_value, the_d1u, the_d1v, the_d2u, the_d2v, the_d2uv, the_offset, 1.0)
}

/// OCCT `CalculateD1` with the explicit `theSign` factor.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn calculate_d1_signed(
    the_value: DVec3,
    the_d1u: DVec3,
    the_d1v: DVec3,
    the_d2u: DVec3,
    the_d2v: DVec3,
    the_d2uv: DVec3,
    the_offset: f64,
    the_sign: f64,
) -> Option<(DVec3, DVec3, DVec3)> {
    let a_norm = compute_normal(the_d1u, the_d1v, THE_D1_MAGNITUDE_TOL)?;
    let value = the_value + the_offset * the_sign * a_norm;
    let a_n1u = compute_dnormal_u(the_d1u, the_d1v, the_d2u, the_d2uv, a_norm);
    let a_n1v = compute_dnormal_v(the_d1u, the_d1v, the_d2v, the_d2uv, a_norm);
    Some((
        value,
        the_d1u + the_offset * the_sign * a_n1u,
        the_d1v + the_offset * the_sign * a_n1v,
    ))
}

// =========================================================================
// ComputeDerivatives (pxx L266-414)
// =========================================================================

/// OCCT `CSLib::DNNUV(theNu, theNv, theDerSurf1, theDerSurf2)` (CSLib.cxx
/// L411-430) — the derivative of the non-normalized vector
/// `N = dS1/du ^ dS2/dv` built from TWO surfaces' derivative tables (the
/// osculating-surface variant).  rcad `math/cs_lib.rs` carries the single
/// array forms only, so this body is re-hosted here.
fn dnnuv_2arrays(nu: i32, nv: i32, der_surf1: &VecTable, der_surf2: &VecTable) -> DVec3 {
    let mut a_result = DVec3::ZERO;
    for i in 0..=nu {
        for j in 0..=nv {
            let a_vg = der_surf1[(i + 1) as usize][j as usize];
            let a_vd = der_surf2[(nu - i) as usize][(nv + 1 - j) as usize];
            let a_cross = a_vg.cross(a_vd);
            let a_bin_coef = crate::math::plib::binomial(nu as usize, i as usize)
                * crate::math::plib::binomial(nv as usize, j as usize);
            a_result += a_bin_coef * a_cross;
        }
    }
    a_result
}

/// The 0-based `NCollection_Array2<gp_Vec>` of the OCCT derivative tables —
/// `value(i, j)` is the OCCT `theArray(i, j)`.
type VecTable = Vec<Vec<DVec3>>;

fn make_table(upper_i: usize, upper_j: usize) -> VecTable {
    vec![vec![DVec3::ZERO; upper_j + 1]; upper_i + 1]
}

/// OCCT `ComputeDerivatives` (pxx L266-414).  `the_max_order` is passed by
/// value in OCCT and zeroed in the osculating branch; `the_der_surf` is the
/// caller's (already partially filled) derivative table, `the_der_nuv` its
/// normal-derivative counterpart.
#[allow(clippy::too_many_arguments)]
fn compute_derivatives(
    mut the_max_order: i32,
    the_min_order: i32,
    the_u: f64,
    the_v: f64,
    the_basis_surf: &Surface3,
    the_nu: i32,
    the_nv: i32,
    the_along_u: bool,
    the_along_v: bool,
    the_osc_surf: Option<&Surface3>,
    the_der_nuv: &mut VecTable,
    the_der_surf: &mut VecTable,
) -> bool {
    if the_along_u || the_along_v {
        the_max_order = 0;
        let osc = the_osc_surf.expect("ComputeDerivatives: null osculating surface");
        // OCCT: NCollection_Array2<gp_Vec> DerSurfL(0, MaxOrder+NU+1,
        // 0, MaxOrder+NV+1).
        let mut der_surf_l = make_table(
            (the_max_order + the_nu + 1) as usize,
            (the_max_order + the_nv + 1) as usize,
        );
        match the_min_order {
            1 => {
                let a_d1 = eval_d1(osc, the_u, the_v);
                der_surf_l[1][0] = a_d1.d1u;
                der_surf_l[0][1] = a_d1.d1v;
            }
            2 => {
                let a_d2 = eval_d2(osc, the_u, the_v);
                der_surf_l[1][0] = a_d2.1;
                der_surf_l[0][1] = a_d2.2;
                der_surf_l[1][1] = a_d2.4;
                der_surf_l[2][0] = a_d2.3;
                der_surf_l[0][2] = a_d2.5;
            }
            3 => {
                let a_d3 = eval_d3(osc, the_u, the_v);
                der_surf_l[1][0] = a_d3.1;
                der_surf_l[0][1] = a_d3.2;
                der_surf_l[1][1] = a_d3.4;
                der_surf_l[2][0] = a_d3.3;
                der_surf_l[0][2] = a_d3.5;
                der_surf_l[3][0] = a_d3.6;
                der_surf_l[2][1] = a_d3.7;
                der_surf_l[1][2] = a_d3.8;
                der_surf_l[0][3] = a_d3.9;
            }
            _ => {}
        }

        if the_nu <= the_nv {
            for i in 0..=(the_max_order + 1 + the_nu) {
                for j in i..=(the_max_order + the_nv + 1) {
                    if i + j > the_min_order {
                        let a_osc_dn = eval_dn(osc, the_u, the_v, i, j);
                        let a_basis_dn = eval_dn(the_basis_surf, the_u, the_v, i, j);
                        der_surf_l[i as usize][j as usize] = a_osc_dn;
                        the_der_surf[i as usize][j as usize] = a_basis_dn;
                        if i != j && j <= the_nu + 1 {
                            let a_basis_dnji = eval_dn(the_basis_surf, the_u, the_v, j, i);
                            let a_osc_dnji = eval_dn(osc, the_u, the_v, j, i);
                            the_der_surf[j as usize][i as usize] = a_basis_dnji;
                            der_surf_l[j as usize][i as usize] = a_osc_dnji;
                        }
                    }
                }
            }
        } else {
            for j in 0..=(the_max_order + 1 + the_nv) {
                for i in j..=(the_max_order + the_nu + 1) {
                    if i + j > the_min_order {
                        let a_osc_dn = eval_dn(osc, the_u, the_v, i, j);
                        let a_basis_dn = eval_dn(the_basis_surf, the_u, the_v, i, j);
                        der_surf_l[i as usize][j as usize] = a_osc_dn;
                        the_der_surf[i as usize][j as usize] = a_basis_dn;
                        if i != j && i <= the_nv + 1 {
                            let a_basis_dnji = eval_dn(the_basis_surf, the_u, the_v, j, i);
                            let a_osc_dnji = eval_dn(osc, the_u, the_v, j, i);
                            the_der_surf[j as usize][i as usize] = a_basis_dnji;
                            der_surf_l[j as usize][i as usize] = a_osc_dnji;
                        }
                    }
                }
            }
        }
        for i in 0..=(the_max_order + the_nu) {
            for j in 0..=(the_max_order + the_nv) {
                if the_along_u {
                    the_der_nuv[i as usize][j as usize] =
                        dnnuv_2arrays(i, j, &der_surf_l, the_der_surf);
                }
                if the_along_v {
                    the_der_nuv[i as usize][j as usize] =
                        dnnuv_2arrays(i, j, the_der_surf, &der_surf_l);
                }
            }
        }
    } else {
        for i in 0..=(the_max_order + the_nu + 1) {
            for j in i..=(the_max_order + the_nv + 1) {
                if i + j > the_min_order {
                    let a_dn = eval_dn(the_basis_surf, the_u, the_v, i, j);
                    the_der_surf[i as usize][j as usize] = a_dn;
                    let upper_row = the_der_surf.len() as i32 - 1;
                    let upper_col = the_der_surf[0].len() as i32 - 1;
                    if i != j && j <= upper_row && i <= upper_col {
                        let a_dnji = eval_dn(the_basis_surf, the_u, the_v, j, i);
                        the_der_surf[j as usize][i as usize] = a_dnji;
                    }
                }
            }
        }
        for i in 0..=(the_max_order + the_nu) {
            for j in 0..=(the_max_order + the_nv) {
                the_der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, the_der_surf);
            }
        }
    }
    true
}

/// OCCT `ReplaceDerivative` (pxx L432-508) — attempts to replace a zero
/// derivative by stepping away and recomputing.  `(d1u, d1v, replaced)`.
#[allow(clippy::too_many_arguments)]
fn replace_derivative(
    the_u: f64,
    the_v: f64,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
    mut the_du: DVec3,
    mut the_dv: DVec3,
    the_square_tol: f64,
    the_basis_surf: &Surface3,
) -> (DVec3, DVec3, bool) {
    let is_replace_du = the_du.length_squared() < the_square_tol;
    let is_replace_dv = the_dv.length_squared() < the_square_tol;
    let mut is_replaced = false;

    // Only handle the case where exactly one derivative is zero.
    if is_replace_du != is_replace_dv {
        // Calculate the step along the non-zero derivative.
        let mut a_step;
        if is_replace_dv {
            a_step = CONFUSION * the_du.length();
            if a_step > the_umax - the_umin {
                a_step = (the_umax - the_umin) / 100.;
            }
        } else {
            a_step = CONFUSION * the_dv.length();
            if a_step > the_vmax - the_vmin {
                a_step = (the_vmax - the_vmin) / 100.;
            }
        }

        // Step away from current parametric coordinates and calculate
        // derivatives once again.  Replace zero derivative by the obtained.
        let mut a_step_sign = -1.0f64;
        while a_step_sign <= 1.0 && !is_replaced {
            let mut a_u = the_u;
            let mut a_v = the_v;
            if is_replace_dv {
                a_u = the_u + a_step_sign * a_step;
                if a_u < the_umin || a_u > the_umax {
                    a_step_sign += 2.0;
                    continue;
                }
            } else {
                a_v = the_v + a_step_sign * a_step;
                if a_v < the_vmin || a_v > the_vmax {
                    a_step_sign += 2.0;
                    continue;
                }
            }

            let a_d1_result = eval_d1(the_basis_surf, a_u, a_v);
            if is_replace_du && a_d1_result.d1u.length_squared() > the_square_tol {
                the_du = a_d1_result.d1u;
                is_replaced = true;
            }
            if is_replace_dv && a_d1_result.d1v.length_squared() > the_square_tol {
                the_dv = a_d1_result.d1v;
                is_replaced = true;
            }
            a_step_sign += 2.0;
        }
    }
    (the_du, the_dv, is_replaced)
}

/// OCCT `ShiftPoint` (pxx L527-568) — shifts the evaluation point towards the
/// center of the parametric space; `false` when the center is overpassed.
#[allow(clippy::too_many_arguments)]
fn shift_point(
    the_u_start: f64,
    the_v_start: f64,
    the_u: &mut f64,
    the_v: &mut f64,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
    the_is_u_periodic: bool,
    the_is_v_periodic: bool,
    the_d1u: DVec3,
    the_d1v: DVec3,
) -> bool {
    let is_u_singular = the_d1u.length_squared() < THE_D1_MAGNITUDE_TOL * THE_D1_MAGNITUDE_TOL;
    let is_v_singular = the_d1v.length_squared() < THE_D1_MAGNITUDE_TOL * THE_D1_MAGNITUDE_TOL;

    let a_dir_u = if the_is_u_periodic || (is_u_singular && !is_v_singular) {
        0.
    } else {
        0.5 * (the_umin + the_umax) - the_u_start
    };
    let a_dir_v = if the_is_v_periodic || (is_v_singular && !is_u_singular) {
        0.
    } else {
        0.5 * (the_vmin + the_vmax) - the_v_start
    };
    let a_dist = (a_dir_u * a_dir_u + a_dir_v * a_dir_v).sqrt();

    let a_du = *the_u - the_u_start;
    let a_dv = *the_v - the_v_start;
    let mut a_step = (2. * (a_du * a_du + a_dv * a_dv).sqrt()).max(crate::core::precision::PCONFUSION);
    if a_step >= a_dist {
        return false;
    }

    a_step /= a_dist;
    *the_u += a_dir_u * a_step;
    *the_v += a_dir_v * a_step;
    true
}

// =========================================================================
// EvaluateD0 / EvaluateD1 (pxx L586-1136)
// =========================================================================

/// OCCT `EvaluateD0` (pxx L586-750) — the pre-computed-D1 overload.  `None`
/// is the OCCT `false` return (the caller raises Geom_UndefinedValue).
#[allow(clippy::too_many_arguments)]
fn evaluate_d0_precomputed(
    the_u0: f64,
    the_v0: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_value_in: DVec3,
    the_d1u_in: DVec3,
    the_d1v_in: DVec3,
) -> Option<DVec3> {
    let a_u_start = the_u0;
    let a_v_start = the_v0;
    let mut the_u = the_u0;
    let mut the_v = the_v0;
    let bounds = SurfaceEval::default_domain(the_basis_surf);
    let (a_umin, a_umax, a_vmin, a_vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    let is_u_per = SurfaceEval::is_u_periodic(the_basis_surf);
    let is_v_per = SurfaceEval::is_v_periodic(the_basis_surf);

    let mut a_d1u = the_d1u_in;
    let mut a_d1v = the_d1v_in;
    let mut the_value = the_value_in;
    let mut is_first_iteration = true;

    loop {
        // For subsequent iterations, recompute D1 at the shifted point.
        if !is_first_iteration {
            let a_d1_result = eval_d1(the_basis_surf, the_u, the_v);
            the_value = a_d1_result.point;
            a_d1u = a_d1_result.d1u;
            a_d1v = a_d1_result.d1v;
        }
        is_first_iteration = false;

        if is_infinite_coord(a_d1u) || is_infinite_coord(a_d1v) {
            return None;
        }

        // Try the non-singular case first.
        if let Some(value) = calculate_d0(the_value, a_d1u, a_d1v, the_offset) {
            return Some(value);
        }

        // Singular case - query the osculating surface and use higher order
        // derivatives.
        let a_max_order = 3;
        let mut a_osc_info = OsculatingInfo::default();
        let mut a_osc_surf: Option<Surface3> = None;
        if let Some(query) = the_osc_query {
            let (along_u, opposite_u, l_u) = query.u_osculating_surface(the_u, the_v);
            let (along_v, opposite_v, l_v) = query.v_osculating_surface(the_u, the_v);
            a_osc_info.along_u = along_u;
            a_osc_info.along_v = along_v;
            a_osc_info.is_opposite = opposite_u || opposite_v;
            a_osc_surf = l_u.or(l_v).map(Surface3::BSpline);
        }

        let mut a_der_nuv = make_table(a_max_order as usize, a_max_order as usize);
        let mut a_der_surf = make_table(
            (a_max_order + 1) as usize,
            (a_max_order + 1) as usize,
        );
        a_der_surf[1][0] = a_d1u;
        a_der_surf[0][1] = a_d1v;

        // Use ComputeDerivatives which handles the osculating surface.
        let ok = if a_osc_info.has_osculating() && a_osc_surf.is_some() {
            compute_derivatives(
                a_max_order,
                1,
                the_u,
                the_v,
                the_basis_surf,
                0,
                0,
                a_osc_info.along_u,
                a_osc_info.along_v,
                a_osc_surf.as_ref(),
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        } else {
            compute_derivatives(
                a_max_order,
                1,
                the_u,
                the_v,
                the_basis_surf,
                0,
                0,
                false,
                false,
                None,
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        };
        if !ok {
            return None;
        }

        let (mut a_n_status, mut a_normal, _order_u, _order_v) = normal_max_order(
            a_max_order,
            &a_der_nuv,
            THE_D1_MAGNITUDE_TOL,
            the_u,
            the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );

        // Handle CSLib_InfinityOfSolutions by replacing the zero derivative.
        if a_n_status == NormalStatus::InfinityOfSolutions {
            let (a_new_du, a_new_dv, replaced) = replace_derivative(
                the_u,
                the_v,
                a_umin,
                a_umax,
                a_vmin,
                a_vmax,
                a_d1u,
                a_d1v,
                THE_D1_MAGNITUDE_TOL * THE_D1_MAGNITUDE_TOL,
                the_basis_surf,
            );
            if replaced {
                // OCCT: CSLib::Normal(aNewDU, aNewDV, THE_D1_MAGNITUDE_TOL,
                // aNStatus, aNormal) — the MagTol overload.
                let (normal, status) =
                    normal_from_derivatives_mag(a_new_du, a_new_dv, THE_D1_MAGNITUDE_TOL);
                a_n_status = status;
                a_normal = normal;
            }
        }

        if a_n_status == NormalStatus::Defined {
            let n = a_normal.expect("CSLib::Normal(Defined) without a direction");
            return Some(the_value + the_offset * a_osc_info.sign() * n);
        }

        // Try shifting the point towards the center - false when overpassed.
        if !shift_point(
            a_u_start,
            a_v_start,
            &mut the_u,
            &mut the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
            is_u_per,
            is_v_per,
            a_d1u,
            a_d1v,
        ) {
            return None;
        }
    }
}

/// OCCT `EvaluateD0` (pxx L766-783) — the convenience overload that computes
/// the basis D1 first.
pub fn evaluate_d0(
    the_u: f64,
    the_v: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
) -> Option<DVec3> {
    let a_basis_d1 = eval_d1(the_basis_surf, the_u, the_v);
    evaluate_d0_precomputed(
        the_u,
        the_v,
        the_basis_surf,
        the_offset,
        the_osc_query,
        a_basis_d1.point,
        a_basis_d1.d1u,
        a_basis_d1.d1v,
    )
}

/// OCCT `EvaluateD1` (pxx L804-1095) — the pre-computed-D2 overload.  `None`
/// is the OCCT `false` return (the caller raises Geom_UndefinedDerivative).
#[allow(clippy::too_many_arguments)]
fn evaluate_d1_precomputed(
    the_u0: f64,
    the_v0: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_value_in: DVec3,
    the_d1u_in: DVec3,
    the_d1v_in: DVec3,
    the_d2u_in: DVec3,
    the_d2v_in: DVec3,
    the_d2uv_in: DVec3,
) -> Option<ResD1> {
    let a_u_start = the_u0;
    let a_v_start = the_v0;
    let mut the_u = the_u0;
    let mut the_v = the_v0;
    let bounds = SurfaceEval::default_domain(the_basis_surf);
    let (a_umin, a_umax, a_vmin, a_vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    let is_u_per = SurfaceEval::is_u_periodic(the_basis_surf);
    let is_v_per = SurfaceEval::is_v_periodic(the_basis_surf);

    let mut the_value = the_value_in;
    let mut the_d1u = the_d1u_in;
    let mut the_d1v = the_d1v_in;
    let mut a_d2u = the_d2u_in;
    let mut a_d2v = the_d2v_in;
    let mut a_d2uv = the_d2uv_in;
    let mut is_first_iteration = true;

    loop {
        // For subsequent iterations, recompute D2 at the shifted point.
        if !is_first_iteration {
            let a_d2_result = eval_d2(the_basis_surf, the_u, the_v);
            the_value = a_d2_result.0;
            the_d1u = a_d2_result.1;
            the_d1v = a_d2_result.2;
            a_d2u = a_d2_result.3;
            a_d2v = a_d2_result.5;
            a_d2uv = a_d2_result.4;
        }
        is_first_iteration = false;

        if is_infinite_coord(the_d1u) || is_infinite_coord(the_d1v) {
            return None;
        }

        // Check if singular by normalizing derivatives and computing the
        // cross product.
        let mut a_d1u = the_d1u;
        let mut a_d1v = the_d1v;
        let a_d1u_norm2 = a_d1u.length_squared();
        let a_d1v_norm2 = a_d1v.length_squared();
        if a_d1u_norm2 > 1.0 {
            a_d1u /= a_d1u_norm2.sqrt();
        }
        if a_d1v_norm2 > 1.0 {
            a_d1v /= a_d1v_norm2.sqrt();
        }

        let mut is_singular = false;
        let a_max_order = 3;
        let mut a_norm = a_d1u.cross(a_d1v);

        // Query the osculating surface only if singular.
        let mut a_osc_info = OsculatingInfo::default();
        let mut a_osc_surf: Option<Surface3> = None;
        if a_norm.length_squared() <= THE_D1_MAGNITUDE_TOL * THE_D1_MAGNITUDE_TOL {
            if let Some(query) = the_osc_query {
                let (along_u, opposite_u, l_u) = query.u_osculating_surface(the_u, the_v);
                let (along_v, opposite_v, l_v) = query.v_osculating_surface(the_u, the_v);
                a_osc_info.along_u = along_u;
                a_osc_info.along_v = along_v;
                a_osc_info.is_opposite = opposite_u || opposite_v;
                a_osc_surf = l_u.or(l_v).map(Surface3::BSpline);
            }
            is_singular = true;
        }

        // Compute the sign factor.
        let a_sign = a_osc_info.sign();

        // Non-singular case: the direct formulas.
        if !is_singular {
            a_norm = a_norm.normalize();
            let value = the_value + the_offset * a_sign * a_norm;

            // Compute normal derivatives using the inline formulas.
            let a_n0 = a_norm;
            let a_scale = the_d1u.cross(the_d1v).dot(a_n0);
            let mut a_n1u = DVec3::new(
                a_d2u.y * the_d1v.z + the_d1u.y * a_d2uv.z - a_d2u.z * the_d1v.y
                    - the_d1u.z * a_d2uv.y,
                (a_d2u.x * the_d1v.z + the_d1u.x * a_d2uv.z - a_d2u.z * the_d1v.x
                    - the_d1u.z * a_d2uv.x)
                    * -1.0,
                a_d2u.x * the_d1v.y + the_d1u.x * a_d2uv.y - a_d2u.y * the_d1v.x
                    - the_d1u.y * a_d2uv.x,
            );
            let a_scale_u = a_n1u.dot(a_n0);
            a_n1u = (a_n1u - a_scale_u * a_n0) / a_scale;

            let mut a_n1v = DVec3::new(
                a_d2uv.y * the_d1v.z + a_d2v.z * the_d1u.y - a_d2uv.z * the_d1v.y
                    - a_d2v.y * the_d1u.z,
                (a_d2uv.x * the_d1v.z + a_d2v.z * the_d1u.x - a_d2uv.z * the_d1v.x
                    - a_d2v.x * the_d1u.z)
                    * -1.0,
                a_d2uv.x * the_d1v.y + a_d2v.y * the_d1u.x - a_d2uv.y * the_d1v.x
                    - a_d2v.x * the_d1u.y,
            );
            let a_scale_v = a_n1v.dot(a_n0);
            a_n1v = (a_n1v - a_scale_v * a_n0) / a_scale;

            return Some(ResD1 {
                point: value,
                d1u: the_d1u + the_offset * a_sign * a_n1u,
                d1v: the_d1v + the_offset * a_sign * a_n1v,
            });
        }

        // Singular case - use higher order derivatives.
        let mut a_der_nuv = make_table(
            (a_max_order + 1) as usize,
            (a_max_order + 1) as usize,
        );
        let mut a_der_surf = make_table(
            (a_max_order + 2) as usize,
            (a_max_order + 2) as usize,
        );
        a_der_surf[1][0] = the_d1u;
        a_der_surf[0][1] = the_d1v;
        a_der_surf[1][1] = a_d2uv;
        a_der_surf[2][0] = a_d2u;
        a_der_surf[0][2] = a_d2v;

        let run_compute_derivatives = |der_nuv: &mut VecTable, der_surf: &mut VecTable| -> bool {
            if a_osc_info.has_osculating() && a_osc_surf.is_some() {
                compute_derivatives(
                    a_max_order,
                    2,
                    the_u,
                    the_v,
                    the_basis_surf,
                    1,
                    1,
                    a_osc_info.along_u,
                    a_osc_info.along_v,
                    a_osc_surf.as_ref(),
                    der_nuv,
                    der_surf,
                )
            } else {
                compute_derivatives(
                    a_max_order,
                    2,
                    the_u,
                    the_v,
                    the_basis_surf,
                    1,
                    1,
                    false,
                    false,
                    None,
                    der_nuv,
                    der_surf,
                )
            }
        };
        if !run_compute_derivatives(&mut a_der_nuv, &mut a_der_surf) {
            return None;
        }

        let (mut a_n_status, mut a_normal, mut a_order_u, mut a_order_v) = normal_max_order(
            a_max_order,
            &a_der_nuv,
            THE_D1_MAGNITUDE_TOL,
            the_u,
            the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );

        // Handle CSLib_InfinityOfSolutions by replacing the zero derivative.
        if a_n_status == NormalStatus::InfinityOfSolutions {
            let (a_new_du, a_new_dv, replaced) = replace_derivative(
                the_u,
                the_v,
                a_umin,
                a_umax,
                a_vmin,
                a_vmax,
                the_d1u,
                the_d1v,
                THE_D1_MAGNITUDE_TOL * THE_D1_MAGNITUDE_TOL,
                the_basis_surf,
            );
            if replaced {
                // Re-compute with the replaced derivatives.
                a_der_surf[1][0] = a_new_du;
                a_der_surf[0][1] = a_new_dv;
                if !run_compute_derivatives(&mut a_der_nuv, &mut a_der_surf) {
                    return None;
                }
                let (status, normal, order_u, order_v) = normal_max_order(
                    a_max_order,
                    &a_der_nuv,
                    THE_D1_MAGNITUDE_TOL,
                    the_u,
                    the_v,
                    a_umin,
                    a_umax,
                    a_vmin,
                    a_vmax,
                );
                a_n_status = status;
                a_normal = normal;
                a_order_u = order_u;
                a_order_v = order_v;
            }
        }

        if a_n_status == NormalStatus::Defined {
            let n = a_normal.expect("CSLib::Normal(Defined) without a direction");
            // Compute the offset point.
            let value = the_value + the_offset * a_sign * n;
            // Compute D1 using CSLib.
            let mut d1u = dnnormal(1, 0, &a_der_nuv, a_order_u, a_order_v);
            let mut d1v = dnnormal(0, 1, &a_der_nuv, a_order_u, a_order_v);
            d1u = d1u * (the_offset * a_sign) + a_der_surf[1][0];
            d1v = d1v * (the_offset * a_sign) + a_der_surf[0][1];
            return Some(ResD1 {
                point: value,
                d1u,
                d1v,
            });
        }

        // Try shifting the point towards the center - false when overpassed.
        if !shift_point(
            a_u_start,
            a_v_start,
            &mut the_u,
            &mut the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
            is_u_per,
            is_v_per,
            the_d1u,
            the_d1v,
        ) {
            return None;
        }
    }
}

/// OCCT `EvaluateD1` (pxx L1112-1136) — the convenience overload that computes
/// the basis D2 first.
pub fn evaluate_d1(
    the_u: f64,
    the_v: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
) -> Option<ResD1> {
    let a_basis_d2 = eval_d2(the_basis_surf, the_u, the_v);
    evaluate_d1_precomputed(
        the_u,
        the_v,
        the_basis_surf,
        the_offset,
        the_osc_query,
        a_basis_d2.0,
        a_basis_d2.1,
        a_basis_d2.2,
        a_basis_d2.3,
        a_basis_d2.5,
        a_basis_d2.4,
    )
}

// =========================================================================
// Geom_OffsetSurfaceUtils::EvaluateD2 / EvaluateDN (pxx L1161-1831)
// =========================================================================

/// OCCT `EvaluateD2` (pxx L1161-1350) — the pre-computed-D3 overload.  `None`
/// is the OCCT `false` return (the caller raises Geom_UndefinedDerivative).
#[allow(clippy::too_many_arguments)]
fn evaluate_d2_precomputed(
    the_u0: f64,
    the_v0: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_value_in: DVec3,
    the_d1u_in: DVec3,
    the_d1v_in: DVec3,
    the_d2u_in: DVec3,
    the_d2v_in: DVec3,
    the_d2uv_in: DVec3,
    the_d3u_in: DVec3,
    the_d3v_in: DVec3,
    the_d3uuv_in: DVec3,
    the_d3uvv_in: DVec3,
) -> Option<ResD2> {
    let a_u_start = the_u0;
    let a_v_start = the_v0;
    let mut the_u = the_u0;
    let mut the_v = the_v0;
    let bounds = SurfaceEval::default_domain(the_basis_surf);
    let (a_umin, a_umax, a_vmin, a_vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    let is_u_per = SurfaceEval::is_u_periodic(the_basis_surf);
    let is_v_per = SurfaceEval::is_v_periodic(the_basis_surf);

    let mut the_value = the_value_in;
    let mut the_d1u = the_d1u_in;
    let mut the_d1v = the_d1v_in;
    let mut the_d2u = the_d2u_in;
    let mut the_d2v = the_d2v_in;
    let mut the_d2uv = the_d2uv_in;
    // Use pre-computed D3 for the first iteration.
    let mut a_d3u = the_d3u_in;
    let mut a_d3v = the_d3v_in;
    let mut a_d3uuv = the_d3uuv_in;
    let mut a_d3uvv = the_d3uvv_in;
    let mut is_first_iteration = true;

    loop {
        // For subsequent iterations, recompute D3 at the shifted point.
        if !is_first_iteration {
            let a_d3_result = eval_d3(the_basis_surf, the_u, the_v);
            the_value = a_d3_result.0;
            the_d1u = a_d3_result.1;
            the_d1v = a_d3_result.2;
            the_d2u = a_d3_result.3;
            the_d2v = a_d3_result.5;
            the_d2uv = a_d3_result.4;
            a_d3u = a_d3_result.6;
            a_d3v = a_d3_result.7;
            a_d3uuv = a_d3_result.8;
            a_d3uvv = a_d3_result.9;
        }
        is_first_iteration = false;

        if is_infinite_coord(the_d1u) || is_infinite_coord(the_d1v) {
            return None;
        }

        // Check if singular using CSLib::Normal on the first-order derivatives.
        let (_a_normal0, a_n_status0) =
            normal_from_derivatives_mag(the_d1u, the_d1v, THE_D1_MAGNITUDE_TOL);

        // MaxOrder = 0 for non-singular, 3 for singular.
        let a_max_order = if a_n_status0 == NormalStatus::Defined { 0 } else { 3 };

        // Get the osculating surface info (singular case only).
        let mut a_osc_info = OsculatingInfo::default();
        let mut a_osc_surf: Option<Surface3> = None;
        if a_n_status0 != NormalStatus::Defined {
            if let Some(query) = the_osc_query {
                let (along_u, opposite_u, l_u) = query.u_osculating_surface(the_u, the_v);
                let (along_v, opposite_v, l_v) = query.v_osculating_surface(the_u, the_v);
                a_osc_info.along_u = along_u;
                a_osc_info.along_v = along_v;
                a_osc_info.is_opposite = opposite_u || opposite_v;
                a_osc_surf = l_u.or(l_v).map(Surface3::BSpline);
            }
        }

        // NCollection_Array2<gp_Vec> aDerNUV(buffer, 0, aMaxOrder + 2,
        // 0, aMaxOrder + 2) / aDerSurf(buffer, 0, aMaxOrder + 3,
        // 0, aMaxOrder + 3).
        let mut a_der_nuv = make_table(
            (a_max_order + 2) as usize,
            (a_max_order + 2) as usize,
        );
        let mut a_der_surf = make_table(
            (a_max_order + 3) as usize,
            (a_max_order + 3) as usize,
        );

        a_der_surf[1][0] = the_d1u;
        a_der_surf[0][1] = the_d1v;
        a_der_surf[1][1] = the_d2uv;
        a_der_surf[2][0] = the_d2u;
        a_der_surf[0][2] = the_d2v;
        a_der_surf[3][0] = a_d3u;
        a_der_surf[2][1] = a_d3uuv;
        a_der_surf[1][2] = a_d3uvv;
        a_der_surf[0][3] = a_d3v;

        // Use ComputeDerivatives to populate DerNUV.
        let ok = if a_osc_info.has_osculating() && a_osc_surf.is_some() {
            compute_derivatives(
                a_max_order,
                3,
                the_u,
                the_v,
                the_basis_surf,
                2,
                2,
                a_osc_info.along_u,
                a_osc_info.along_v,
                a_osc_surf.as_ref(),
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        } else {
            compute_derivatives(
                a_max_order,
                3,
                the_u,
                the_v,
                the_basis_surf,
                2,
                2,
                false,
                false,
                None,
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        };
        if !ok {
            return None;
        }

        // Compute the normal using CSLib with MaxOrder.
        let (a_n_status, a_normal, a_order_u, a_order_v) = normal_max_order(
            a_max_order,
            &a_der_nuv,
            THE_D1_MAGNITUDE_TOL,
            the_u,
            the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );

        if a_n_status == NormalStatus::Defined {
            let n = a_normal.expect("CSLib::Normal(Defined) without a direction");
            let a_sign = the_offset * a_osc_info.sign();

            // Compute the offset point.
            let value = the_value + a_sign * n;

            // Compute D1 using CSLib::DNNormal.
            let d1u = a_der_surf[1][0] + dnnormal(1, 0, &a_der_nuv, a_order_u, a_order_v) * a_sign;
            let d1v = a_der_surf[0][1] + dnnormal(0, 1, &a_der_nuv, a_order_u, a_order_v) * a_sign;

            // For D2, re-fetch from the basis surface.
            let a_dn20 = eval_dn(the_basis_surf, the_u, the_v, 2, 0);
            let a_dn02 = eval_dn(the_basis_surf, the_u, the_v, 0, 2);
            let a_dn11 = eval_dn(the_basis_surf, the_u, the_v, 1, 1);
            let d2u = a_dn20 + dnnormal(2, 0, &a_der_nuv, a_order_u, a_order_v) * a_sign;
            let d2v = a_dn02 + dnnormal(0, 2, &a_der_nuv, a_order_u, a_order_v) * a_sign;
            let d2uv = a_dn11 + dnnormal(1, 1, &a_der_nuv, a_order_u, a_order_v) * a_sign;
            return Some(ResD2 {
                point: value,
                d1u,
                d1v,
                d2u,
                d2v,
                d2uv,
            });
        }

        // Try shifting the point towards the center - false when overpassed.
        if !shift_point(
            a_u_start,
            a_v_start,
            &mut the_u,
            &mut the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
            is_u_per,
            is_v_per,
            the_d1u,
            the_d1v,
        ) {
            return None;
        }
    }
}

/// OCCT `EvaluateD2` (pxx L1370-1404) — the convenience overload that computes
/// the basis D3 first.
pub fn evaluate_d2(
    the_u: f64,
    the_v: f64,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
) -> Option<ResD2> {
    let a_basis_d3 = eval_d3(the_basis_surf, the_u, the_v);
    evaluate_d2_precomputed(
        the_u,
        the_v,
        the_basis_surf,
        the_offset,
        the_osc_query,
        a_basis_d3.0,
        a_basis_d3.1,
        a_basis_d3.2,
        a_basis_d3.3,
        a_basis_d3.5,
        a_basis_d3.4,
        a_basis_d3.6,
        a_basis_d3.7,
        a_basis_d3.8,
        a_basis_d3.9,
    )
}

/// OCCT `EvaluateDN` (pxx L1637-1794) — the pre-computed-D1 overload.  `None`
/// is the OCCT `false` return (the caller raises Geom_UndefinedDerivative).
#[allow(clippy::too_many_arguments)]
fn evaluate_dn_precomputed(
    the_u0: f64,
    the_v0: f64,
    the_nu: i32,
    the_nv: i32,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_d1u_in: DVec3,
    the_d1v_in: DVec3,
) -> Option<DVec3> {
    let a_u_start = the_u0;
    let a_v_start = the_v0;
    let mut the_u = the_u0;
    let mut the_v = the_v0;
    let bounds = SurfaceEval::default_domain(the_basis_surf);
    let (a_umin, a_umax, a_vmin, a_vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    let is_u_per = SurfaceEval::is_u_periodic(the_basis_surf);
    let is_v_per = SurfaceEval::is_v_periodic(the_basis_surf);

    // Use pre-computed D1 for the first iteration.
    let mut a_d1u = the_d1u_in;
    let mut a_d1v = the_d1v_in;
    let mut is_first_iteration = true;

    loop {
        // For subsequent iterations, recompute D1 at the shifted point.
        if !is_first_iteration {
            let a_d1_result = eval_d1(the_basis_surf, the_u, the_v);
            a_d1u = a_d1_result.d1u;
            a_d1v = a_d1_result.d1v;
        }
        is_first_iteration = false;

        if is_infinite_coord(a_d1u) || is_infinite_coord(a_d1v) {
            return None;
        }

        // Check if singular to determine MaxOrder.
        let (_a_normal0, a_n_status0) =
            normal_from_derivatives_mag(a_d1u, a_d1v, THE_D1_MAGNITUDE_TOL);
        let a_max_order = if a_n_status0 == NormalStatus::Defined { 0 } else { 3 };

        // NCollection_Array2<gp_Vec> aDerNUV(buffer, 0, aMaxOrder + theNu,
        // 0, aMaxOrder + theNv) / aDerSurf(buffer, 0, aMaxOrder + theNu + 1,
        // 0, aMaxOrder + theNv + 1).
        let mut a_der_nuv = make_table(
            (a_max_order + the_nu) as usize,
            (a_max_order + the_nv) as usize,
        );
        let mut a_der_surf = make_table(
            (a_max_order + the_nu + 1) as usize,
            (a_max_order + the_nv + 1) as usize,
        );

        a_der_surf[1][0] = a_d1u;
        a_der_surf[0][1] = a_d1v;

        // Check the osculating surface only in the singular case.
        let mut a_osc_info = OsculatingInfo::default();
        let mut a_osc_surf: Option<Surface3> = None;
        if a_n_status0 != NormalStatus::Defined {
            if let Some(query) = the_osc_query {
                let (along_u, opposite_u, l_u) = query.u_osculating_surface(the_u, the_v);
                let (along_v, opposite_v, l_v) = query.v_osculating_surface(the_u, the_v);
                a_osc_info.along_u = along_u;
                a_osc_info.along_v = along_v;
                a_osc_info.is_opposite = opposite_u || opposite_v;
                a_osc_surf = l_u.or(l_v).map(Surface3::BSpline);
            }
        }

        // Use ComputeDerivatives.
        let ok = if a_osc_info.has_osculating() && a_osc_surf.is_some() {
            compute_derivatives(
                a_max_order,
                1,
                the_u,
                the_v,
                the_basis_surf,
                the_nu,
                the_nv,
                a_osc_info.along_u,
                a_osc_info.along_v,
                a_osc_surf.as_ref(),
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        } else {
            compute_derivatives(
                a_max_order,
                1,
                the_u,
                the_v,
                the_basis_surf,
                the_nu,
                the_nv,
                false,
                false,
                None,
                &mut a_der_nuv,
                &mut a_der_surf,
            )
        };
        if !ok {
            return None;
        }

        // Compute the normal with CSLib.
        let (a_n_status, a_normal, a_order_u, a_order_v) = normal_max_order(
            a_max_order,
            &a_der_nuv,
            THE_D1_MAGNITUDE_TOL,
            the_u,
            the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
        );

        if a_n_status == NormalStatus::Defined {
            let _n = a_normal.expect("CSLib::Normal(Defined) without a direction");
            let a_sign = the_offset * a_osc_info.sign();

            // Compute the DN result: basis DN + offset * DNNormal.
            let a_basis_dn = eval_dn(the_basis_surf, the_u, the_v, the_nu, the_nv);
            return Some(
                a_basis_dn + dnnormal(the_nu, the_nv, &a_der_nuv, a_order_u, a_order_v) * a_sign,
            );
        }

        // Try shifting the point towards the center - false when overpassed.
        if !shift_point(
            a_u_start,
            a_v_start,
            &mut the_u,
            &mut the_v,
            a_umin,
            a_umax,
            a_vmin,
            a_vmax,
            is_u_per,
            is_v_per,
            a_d1u,
            a_d1v,
        ) {
            return None;
        }
    }
}

/// OCCT `EvaluateDN` (pxx L1811-1831) — the convenience overload that computes
/// the basis D1 first.
pub fn evaluate_dn(
    the_u: f64,
    the_v: f64,
    the_nu: i32,
    the_nv: i32,
    the_basis_surf: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
) -> Option<DVec3> {
    let a_basis_d1 = eval_d1(the_basis_surf, the_u, the_v);
    evaluate_dn_precomputed(
        the_u,
        the_v,
        the_nu,
        the_nv,
        the_basis_surf,
        the_offset,
        the_osc_query,
        a_basis_d1.d1u,
        a_basis_d1.d1v,
    )
}

// The `Geom_OffsetSurface` class bodies (the `SetBasisSurface` unwrap,
// `Surface()` and the `EvalD0`/`EvalD1`/`EvalD2`/`EvalDN` re-hosts) live in
// `offset_surface_utils_b.rs` and are re-exported here so that the existing
// `crate::geom::offset_surface_utils::*` import paths keep working.
pub use super::offset_surface_utils_b::*;

#[cfg(test)]
mod eval_tests {
    use super::*;
    use crate::geom::{BezierSurface, CylindricalSurface, OffsetSurface, SurfaceEval};

    /// The unit quarter circle as a degree-2 rational Bezier: poles
    /// `{(1,0), (1,1), (0,1)}`, weights `{1, sqrt(2)/2, 1}`, extruded along Z
    /// (degree 1, unit weights) so that `S(u, v) = (X(u), Y(u), v)`.
    ///
    /// The exact derivatives at `u = 0` follow from the quotient rule on
    /// `N(u)/W(u)` over the Bernstein polynomials:
    ///   `P(0) = (1, 0)`, `P'(0) = (0, sqrt(2))`,
    ///   `P''(0) = (-2, -2 + 2*sqrt(2))`.
    fn arc_patches() -> (BSplineSurface, BezierSurface) {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let control_points = vec![
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 1.0)],
            vec![DVec3::new(1.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0)],
            vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(0.0, 1.0, 1.0)],
        ];
        let weights = vec![vec![1.0, 1.0], vec![s, s], vec![1.0, 1.0]];
        let bs = BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: control_points.clone(),
            weights: weights.clone(),
            is_periodic_u: false,
            is_periodic_v: false,
        };
        (bs, BezierSurface { control_points, weights })
    }

    /// `Geom_Surface::EvalDN` over a rational BSpline basis — the branch that
    /// used to panic with the `BSplSLib::RationalDerivative` GAP.
    #[test]
    fn eval_dn_rational_bspline_basis() {
        let (bs, _bez) = arc_patches();
        let s = Surface3::BSpline(bs);
        let d1u = eval_dn(&s, 0.0, 0.5, 1, 0);
        assert!(
            (d1u - DVec3::new(0.0, std::f64::consts::SQRT_2, 0.0)).length() < 1e-12,
            "d1u={d1u:?}"
        );
        let d1v = eval_dn(&s, 0.0, 0.5, 0, 1);
        assert!((d1v - DVec3::Z).length() < 1e-12, "d1v={d1v:?}");
        let d2u = eval_dn(&s, 0.0, 0.5, 2, 0);
        let want = DVec3::new(-2.0, -2.0 + 2.0 * std::f64::consts::SQRT_2, 0.0);
        assert!((d2u - want).length() < 1e-12, "d2u={d2u:?} want={want:?}");
        let d2uv = eval_dn(&s, 0.0, 0.5, 1, 1);
        assert!(d2uv.length() < 1e-12, "d2uv={d2uv:?}");
        let d2v = eval_dn(&s, 0.0, 0.5, 0, 2);
        assert!(d2v.length() < 1e-12, "d2v={d2v:?}");
    }

    /// `Geom_Surface::EvalDN` over a rational Bezier basis
    /// (`Geom_BezierSurface::EvalDN`) and its agreement with the BSpline form
    /// of the same surface.
    #[test]
    fn eval_dn_rational_bezier_basis() {
        let (bs, bez) = arc_patches();
        let s_bs = Surface3::BSpline(bs);
        let s_bez = Surface3::Bezier(bez);
        for (nu, nv) in [(1, 0), (0, 1), (2, 0), (1, 1)] {
            let a = eval_dn(&s_bs, 0.25, 0.5, nu, nv);
            let b = eval_dn(&s_bez, 0.25, 0.5, nu, nv);
            assert!((a - b).length() < 1e-12, "({nu},{nv}) bs={a:?} bez={b:?}");
        }
    }

    /// `Geom_Surface::EvalD2` and `EvalD1` on the rational kinds — the
    /// derivative tables `ComputeDerivatives` fills.
    #[test]
    fn eval_d1_d2_rational_basis() {
        let (bs, bez) = arc_patches();
        for s in [Surface3::BSpline(bs), Surface3::Bezier(bez)] {
            let d1 = eval_d1(&s, 0.0, 0.5);
            assert!((d1.point - DVec3::new(1.0, 0.0, 0.5)).length() < 1e-12);
            assert!(
                (d1.d1u - DVec3::new(0.0, std::f64::consts::SQRT_2, 0.0)).length() < 1e-12,
                "d1u={:?}",
                d1.d1u
            );
            let d2 = eval_d2(&s, 0.0, 0.5);
            assert!((d2.3 - eval_dn(&s, 0.0, 0.5, 2, 0)).length() < 1e-12);
            assert!((d2.4 - eval_dn(&s, 0.0, 0.5, 1, 1)).length() < 1e-12);
            assert!((d2.5 - eval_dn(&s, 0.0, 0.5, 0, 2)).length() < 1e-12);
        }
    }

    // ---------------------------------------------------------------------
    // Geom_OffsetSurface::EvalD0/EvalD1/EvalD2/EvalDN (cxx L338-497)
    // ---------------------------------------------------------------------

    /// A cylindrical basis short-circuits through `Geom_OffsetSurface::Surface()`
    /// (cxx L908-925): the equivalent surface is the cylinder of radius
    /// `R + aSign * offset`, whose exact derivatives are known in closed form
    /// (`P = ((R+d) cos u, (R+d) sin u, v)`).
    #[test]
    fn offset_payload_d1_d2_dn_cylinder_basis() {
        use crate::geom::offset_surface_utils as utils;
        let radius = 2.0;
        let offset = 0.5;
        let of = OffsetSurface {
            basis: Box::new(Surface3::Cylinder(CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::Z,
                radius,
                ref_dir: DVec3::X,
                y_dir: None,
            })),
            offset_distance: offset,
        };
        let r = radius + offset;
        let (u, v) = (0.3f64, 1.5f64);
        let (su, cu) = u.sin_cos();
        let want_p = DVec3::new(r * cu, r * su, v);
        let want_d1u = DVec3::new(-r * su, r * cu, 0.0);
        let want_d1v = DVec3::Z;

        let p = utils::offset_payload_eval_d0(&of, u, v);
        assert!((p - want_p).length() < 1e-14, "p={p:?} want={want_p:?}");
        let d1 = utils::offset_payload_eval_d1(&of, u, v);
        assert!((d1.point - want_p).length() < 1e-14);
        assert!((d1.d1u - want_d1u).length() < 1e-14, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-14, "d1v={:?}", d1.d1v);
        let d2 = utils::offset_payload_eval_d2(&of, u, v);
        assert!((d2.point - want_p).length() < 1e-13);
        assert!((d2.d1u - want_d1u).length() < 1e-13);
        assert!((d2.d1v - want_d1v).length() < 1e-13);
        let want_d2u = DVec3::new(-r * cu, -r * su, 0.0);
        // The equivalent-surface `derivatives2` of a cylinder is the trait's
        // finite-difference default (only `derivatives` is analytic there), so
        // the D2 terms are only accurate to the default step.
        assert!((d2.d2u - want_d2u).length() < 1e-5, "d2u={:?}", d2.d2u);
        assert!(d2.d2v.length() < 1e-5, "d2v={:?}", d2.d2v);
        assert!(d2.d2uv.length() < 1e-5, "d2uv={:?}", d2.d2uv);
        // Geom_Surface::EvalDN goes through the equivalent cylinder's ElSLib
        // form, which is exact.
        assert!((utils::offset_payload_eval_dn(&of, u, v, 1, 0) - want_d1u).length() < 1e-13);
        assert!((utils::offset_payload_eval_dn(&of, u, v, 0, 1) - want_d1v).length() < 1e-13);
        assert!((utils::offset_payload_eval_dn(&of, u, v, 2, 0) - want_d2u).length() < 1e-13);
    }

    /// The parabolic cylinder patch `S(u, v) = (u, u^2, v)` as a degree-2
    /// Bezier in U (poles `(0,0)`, `(0.5,0)`, `(1,1)` in the XY plane) and
    /// degree 1 in V.  A BSpline basis has no equivalent surface, so this
    /// exercises the `Geom_OffsetSurfaceUtils::EvaluateD0/D1/D2/DN` bodies:
    /// with `s = sqrt(1 + 4u^2)` the surface normal
    /// `N = dS/du ^ dS/dv = (2u, -1, 0)/s` gives
    ///   `P   = (u + 2du/s, u^2 - d/s, v)`
    ///   `Pu  = (1 + 2d/s^3, 2u + 4du/s^3, 0)`
    ///   `Puu = (-24du/s^5, 2 + 4d/s^3 - 48du^2/s^5, 0)`
    /// with `Pv = (0, 0, 1)` and all V/V-mixed second derivatives zero.
    #[test]
    fn offset_payload_d1_d2_dn_parabolic_bspline_basis() {
        use crate::geom::offset_surface_utils as utils;
        let basis = Surface3::BSpline(BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 0.0, 1.0)],
                vec![DVec3::new(0.5, 0.0, 0.0), DVec3::new(0.5, 0.0, 1.0)],
                vec![DVec3::new(1.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0], vec![1.0, 1.0]],
            is_periodic_u: false,
            is_periodic_v: false,
        });
        let d = 0.7;
        let of = OffsetSurface {
            basis: Box::new(basis),
            offset_distance: d,
        };
        let (u, v) = (0.4f64, 0.3f64);
        let s = (1.0 + 4.0 * u * u).sqrt();
        let s3 = s * s * s;
        let s5 = s3 * s * s;
        let want_p = DVec3::new(u + 2.0 * d * u / s, u * u - d / s, v);
        let want_d1u = DVec3::new(1.0 + 2.0 * d / s3, 2.0 * u + 4.0 * d * u / s3, 0.0);
        let want_d1v = DVec3::Z;
        let want_d2u = DVec3::new(
            -24.0 * d * u / s5,
            2.0 + 4.0 * d / s3 - 48.0 * d * u * u / s5,
            0.0,
        );

        let p = utils::offset_payload_eval_d0(&of, u, v);
        assert!((p - want_p).length() < 1e-12, "p={p:?} want={want_p:?}");
        let d1 = utils::offset_payload_eval_d1(&of, u, v);
        assert!((d1.point - want_p).length() < 1e-12);
        assert!((d1.d1u - want_d1u).length() < 1e-12, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-12, "d1v={:?}", d1.d1v);
        let d2 = utils::offset_payload_eval_d2(&of, u, v);
        assert!((d2.point - want_p).length() < 1e-11);
        assert!((d2.d1u - want_d1u).length() < 1e-11);
        assert!((d2.d1v - want_d1v).length() < 1e-11);
        assert!((d2.d2u - want_d2u).length() < 1e-11, "d2u={:?}", d2.d2u);
        assert!(d2.d2v.length() < 1e-11, "d2v={:?}", d2.d2v);
        assert!(d2.d2uv.length() < 1e-11, "d2uv={:?}", d2.d2uv);
        assert!(
            (utils::offset_payload_eval_dn(&of, u, v, 0, 1) - want_d1v).length() < 1e-11,
            "dn(0,1)={:?}",
            utils::offset_payload_eval_dn(&of, u, v, 0, 1)
        );
        assert!((utils::offset_payload_eval_dn(&of, u, v, 1, 1) - d2.d2uv).length() < 1e-11);
        // Nu > Nv: OCCT's derivative table loop only visits the `j >= i`
        // triangle (Geom_OffsetSurfaceUtils.pxx L389-404), so for a non-square
        // table the transposed cell `DerSurf(2,0)` is never written and
        // `CSLib::DNNormal` yields zero — the result is the basis derivative.
        // The translated body reproduces that asymmetry verbatim.
        let basis_d1u = DVec3::new(1.0, 2.0 * u, 0.0);
        assert!(
            (utils::offset_payload_eval_dn(&of, u, v, 1, 0) - basis_d1u).length() < 1e-11,
            "dn(1,0)={:?}",
            utils::offset_payload_eval_dn(&of, u, v, 1, 0)
        );
    }
}
