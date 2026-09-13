//! OCCT Geom_OffsetSurfaceUtils (TKG3d/Geom/Geom_OffsetSurfaceUtils.pxx,
//! L41-1833) — 1:1 translation of the internal helper namespace of
//! Geom_OffsetSurface, plus the rcad re-hosts of the Geom_OffsetSurface
//! bodies that consume it:
//!   - `Geom_OffsetSurface::Surface()` (Geom_OffsetSurface.cxx L867-993) —
//!     [`offset_equivalent_surface`], the equivalent non-offset surface;
//!   - `Geom_OffsetSurface::SetBasisSurface` unwrap (L143-170) —
//!     [`offset_basis_and_value`] and the `myOscSurf` construction condition
//!     (L258-267) — [`offset_surface_osculating`];
//!   - `Geom_OffsetSurface::EvalD0` / `EvalD1` (L338-389) —
//!     [`offset_surface_eval_d0`] / [`offset_surface_eval_d1`].
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
use crate::geom::{
    BSplineSurface, ConicalSurface, CylindricalSurface, OffsetSurface, Plane, SphericalSurface,
    Surface3, SurfaceEval, ToroidalSurface,
};
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
/// engine (ElSLib::DN for the quadrics, [`bspl_slib_dn`] for the polynomial
/// kind).
pub fn eval_dn(the_s: &Surface3, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    match the_s {
        Surface3::Plane(_)
        | Surface3::Cylinder(_)
        | Surface3::Cone(_)
        | Surface3::Sphere(_)
        | Surface3::Torus(_) => the_s.dn(u, v, nu, nv),
        Surface3::BSpline(bs) => bspl_slib_dn(bs, u, v, nu, nv),
        _ => panic!(
            "GAP: Geom_Surface::EvalDN (TKG3d/Geom) is not translated for this surface type \
             (the rcad GeomAdaptor_Surface DN engine covers ElSLib surfaces and \
             Geom_BSplineSurface only) — Geom_OffsetSurfaceUtils::ComputeDerivatives"
        ),
    }
}

/// OCCT BSplSLib::DN (BSplSLib.cxx L1519-1605) — the (Nu, Nv) derivative of a
/// BSpline surface at (U, V), the engine of `Geom_BSplineSurface::EvalDN`.
///
/// Faithful for the non-rational branch: `PrepareEval` (the shared
/// `math/bspl_lib.rs`-family re-host of
/// [`super::osculating_surface`]) + the two `BSplCLib::Bohm` passes + the
/// `(n1 * (d2 + 1) + n2)` pole extraction.
///
/// The rational branch of OCCT DN calls
/// `BSplSLib::RationalDerivative(UDegree, VDegree, Nu, Nv, Poles, Results,
/// false)` (BSplSLib.cxx L87-300), which is NOT translated in rcad: the
/// rational branch keeps the GAP panic.
///
/// Note: `rcad-kernel/src/geom/eval.rs::bspline_surface_dn` is a re-host of
/// the same leaf whose derivative-result buffers are sized for a single pole
/// (`PolesResult`/`WeightsResult` of length `dim`), so `eval_homogeneous`
/// overruns them for Nu >= 1 ("index out of bounds"); the derivative does not
/// exist there and this translation is used instead.
fn bspl_slib_dn(bs: &BSplineSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    let u_deg = bs.degree_u as i32;
    let v_deg = bs.degree_v as i32;
    let u_rat = bs.is_rational_u();
    let v_rat = bs.is_rational_v();
    // OCCT: UIndex = VIndex = 0 (the LocateParameter arm of PrepareEval).
    let pe = super::osculating_surface::prepare_eval(
        u,
        v,
        0,
        0,
        u_deg,
        v_deg,
        u_rat,
        v_rat,
        bs.is_periodic_u,
        bs.is_periodic_v,
        &bs.control_points,
        Some(&bs.weights),
        &bs.knots_u,
        &bs.knots_v,
    );
    let dim: usize = if pe.rational { 4 } else { 3 };
    let u_first = pe.flag_u_or_v;
    let d1 = pe.d1;
    let d2 = pe.d2;
    let n1 = if u_first { nu } else { nv };
    let n2 = if u_first { nv } else { nu };

    if !pe.rational && (nu > u_deg || nv > v_deg) {
        return DVec3::ZERO;
    }

    let mut poles = pe.dc.poles;
    let row = dim * (d2 + 1) as usize;
    super::osculating_surface::bspl_clib_bohm(pe.u1, d1, n1, &pe.dc.knots1, row, &mut poles);
    for k in 0..=n1.min(d1) {
        let off = (k * dim as i32 * (d2 + 1)) as usize;
        let end = off + row;
        super::osculating_surface::bspl_clib_bohm(pe.u2, d2, n2, &pe.dc.knots2, dim, &mut poles[off..end]);
    }

    if pe.rational {
        panic!(
            "GAP: BSplSLib::RationalDerivative (TKM/BSplSLib) is not translated — the rational \
             branch of BSplSLib::DN / Geom_BSplineSurface::EvalDN"
        );
    }
    let idx = ((n1 * (d2 + 1) + n2) * dim as i32) as usize;
    DVec3::new(poles[idx], poles[idx + 1], poles[idx + 2])
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
// Geom_OffsetSurface::SetBasisSurface unwrap + Surface() + EvalD0/EvalD1
// =========================================================================

/// OCCT `Geom_OffsetSurface::SetBasisSurface` unwrap (cxx L150-170): the
/// nested `Geom_RectangularTrimmedSurface` / `Geom_OffsetSurface` wrappers are
/// unwrapped and the offset values are accumulated into a single
/// `offsetValue`.  Returns `(aCheckingSurf, accumulated_offset)`.
pub fn offset_basis_and_value(the_surf: &Surface3, the_offset: f64) -> (Surface3, f64) {
    let mut a_checking_surf = the_surf.clone();
    let mut offset_value = the_offset;
    while matches!(
        a_checking_surf,
        Surface3::Trimmed(_) | Surface3::Offset(_)
    ) {
        if let Surface3::Trimmed(t) = &a_checking_surf {
            let basis = t.basis.as_ref().clone();
            a_checking_surf = basis;
        }
        if let Surface3::Offset(o) = &a_checking_surf {
            let basis = o.basis.as_ref().clone();
            offset_value += o.offset_distance;
            a_checking_surf = basis;
        }
    }
    (a_checking_surf, offset_value)
}

/// OCCT `Geom_OffsetSurface::Surface()` (cxx L867-993) — the equivalent
/// non-offset surface of this offset surface, or `None` when none exists.
///
/// The rcad counterpart of the OCCT `myEvalRep` member
/// (`GeomEval_RepSurfaceDesc::Full`, built by `makeFullSurfaceRep(Surface())`
/// in `SetBasisSurface` L255-256) is this recomputation; `directRepSurface`
/// (cxx L83-95) is therefore `offset_equivalent_surface(...) != None` for a
/// freshly constructed surface.
pub fn offset_equivalent_surface(the_surf: &Surface3, the_offset: f64) -> Option<Surface3> {
    if the_offset == 0.0 {
        // Direct case - no offset.
        return Some(the_surf.clone());
    }

    let tol = CONFUSION;
    let mut result: Option<Surface3> = None;

    // Handle trimmed surfaces - extract the basis surface and bounds.
    let (base, is_trimmed, u1, u2, v1, v2) = match the_surf {
        Surface3::Trimmed(t) => {
            let b = t.basis.as_ref();
            (b.clone(), true, t.trim[0], t.trim[1], t.trim[2], t.trim[3])
        }
        other => (other.clone(), false, 0., 0., 0., 0.),
    };

    // Handle canonical surfaces - compute the equivalent offset surface.
    // For direct orientation, offset is along the outward normal; for
    // indirect, it is reversed.
    match &base {
        Surface3::Plane(p) => {
            // Plane normal is already available as Position().Direction().
            let t = p.normal * the_offset;
            // OCCT Geom_Plane::Translated(T) — the translated geometry.
            result = Some(Surface3::Plane(Plane {
                origin: p.origin + t,
                ..p.clone()
            }));
        }
        Surface3::Cylinder(c) => {
            // OCCT: gp_Ax3 Axis = C->Position(); aSign = Axis.Direct() ? 1 : -1.
            let a_sign = ax3_direct(c.y_dir, c.ref_dir, c.axis);
            let radius = c.radius + a_sign * the_offset;
            if radius >= tol {
                result = Some(Surface3::Cylinder(CylindricalSurface {
                    radius,
                    ..c.clone()
                }));
            } else if radius <= -tol {
                // Negative radius: flip the X-axis to reverse the normal
                // orientation.
                result = Some(Surface3::Cylinder(CylindricalSurface {
                    radius: -radius,
                    ref_dir: -c.ref_dir,
                    y_dir: c.y_dir.map(|y| -y),
                    ..c.clone()
                }));
            }
            // else: degenerate surface - radius is too small.
        }
        Surface3::Cone(c) => {
            // OCCT: gp_Ax3 anAxis = C->Position(); aSign = anAxis.Direct() ? 1 : -1.
            let a_sign = ax3_direct(None, c.ref_dir, c.axis);
            let an_alpha = c.half_angle_rad;
            let a_cos = an_alpha.cos();
            let a_sin = an_alpha.sin();
            let a_radius = c.radius + a_sign * the_offset * a_cos;
            if a_radius >= 0. {
                // Translate the apex along the axis by the offset component
                // (anAxis.Translate(aZ) moves the gp_Ax3 location).
                let a_z = c.axis * (-a_sign * the_offset * a_sin);
                result = Some(Surface3::Cone(ConicalSurface {
                    apex: c.apex + a_z,
                    radius: a_radius,
                    ..c.clone()
                }));
            }
            // else: degenerate surface - radius is negative.
        }
        Surface3::Sphere(s) => {
            let a_sign = ax3_direct(None, s.ref_dir, s.axis);
            let radius = s.radius + a_sign * the_offset;
            if radius >= tol {
                result = Some(Surface3::Sphere(SphericalSurface {
                    radius,
                    ..s.clone()
                }));
            } else if radius <= -tol {
                // Negative radius: flip both X and Z axes to reverse the
                // normal orientation.
                result = Some(Surface3::Sphere(SphericalSurface {
                    radius: -radius,
                    axis: -s.axis,
                    ref_dir: -s.ref_dir,
                    ..s.clone()
                }));
            }
            // else: degenerate surface - radius is too small.
        }
        Surface3::Torus(t) => {
            let major_radius = t.major_radius;
            let a_sign = ax3_direct(None, t.ref_dir, t.axis);
            let minor_radius = t.minor_radius + a_sign * the_offset;
            // Only handle the non-self-intersecting torus
            // (MinorRadius <= MajorRadius).
            if minor_radius >= tol && minor_radius <= major_radius {
                result = Some(Surface3::Torus(ToroidalSurface {
                    minor_radius,
                    ..t.clone()
                }));
            } else if minor_radius <= -tol && -minor_radius <= major_radius {
                // Negative minor radius: flip the X-axis to reverse the
                // normal orientation.
                result = Some(Surface3::Torus(ToroidalSurface {
                    minor_radius: -minor_radius,
                    ref_dir: -t.ref_dir,
                    ..t.clone()
                }));
            }
            // else: degenerate or self-intersecting torus - no equivalent
            // surface.
        }
        _ => {}
    }

    // Trim the result if the basis surface was trimmed.
    if is_trimmed {
        if let Some(b) = result.take() {
            result = Some(Surface3::Trimmed(super::TrimmedSurface::new(
                b, u1, u2, v1, v2,
            )));
        }
    }

    result
}

/// OCCT `gp_Ax3::Direct()` for the rcad surface frames: the frame
/// (X = `ref_dir`, Y, Z = `axis`) is right-handed (gp_Ax3.hxx L57-58
/// `Y = Z ^ X` for the direct form; the indirect form reverses Y,
/// gp_Ax3.cxx L44-46).
///
/// Architecture difference: only [`CylindricalSurface`] carries an explicit
/// frame Y direction in the rcad payload (`y_dir`), so a cylinder handedness
/// is decoded from `det(ref_dir, y_dir, axis)`; the
/// [`ConicalSurface`] / [`SphericalSurface`] / [`ToroidalSurface`] payloads
/// store `ref_dir` + `axis` only and can therefore represent the right-handed
/// frames exclusively — the OCCT `aSign = -1` arms of the cone / sphere /
/// torus are unreachable in rcad.
fn ax3_direct(y_dir: Option<DVec3>, ref_dir: DVec3, axis: DVec3) -> f64 {
    match y_dir {
        None => 1.0,
        Some(y) => {
            if ref_dir.dot(y.cross(axis)) > 0.0 {
                1.0
            } else {
                -1.0
            }
        }
    }
}

/// OCCT `Geom_OffsetSurface::SetBasisSurface` osculating-surface construction
/// (cxx L258-267): `myOscSurf` is built only for a BSpline / Bezier basis
/// (after the wrapper unwrap), with the hard-coded `Precision::Confusion()`
/// tolerance of `SetBasisSurface`.  `None` is the OCCT null `myOscSurf`.
pub fn offset_surface_osculating(the_basis: &Surface3) -> Option<OsculatingSurface> {
    // OCCT: constexpr double Tol = Precision::Confusion();
    let tol = CONFUSION;
    match the_basis {
        Surface3::BSpline(_) | Surface3::Bezier(_) => {
            Some(OsculatingSurface::with_surface(the_basis, tol))
        }
        _ => None,
    }
}

/// OCCT `Geom_OffsetSurface::EvalD0` (cxx L338-358) — the point of the offset
/// surface, without the `GeomEval_RepUtils::TryEvalSurfaceD0(myEvalRep, ...)`
/// short circuit (that member is recomputed by
/// [`offset_equivalent_surface`]; the caller decides whether to take it).
pub fn offset_surface_eval_d0(
    the_basis: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_u: f64,
    the_v: f64,
) -> Option<DVec3> {
    evaluate_d0(the_u, the_v, the_basis, the_offset, the_osc_query)
}

/// OCCT `Geom_OffsetSurface::EvalD1` (cxx L362-389) — the point + first
/// partials of the offset surface, without the
/// `GeomEval_RepUtils::TryEvalSurfaceD1(myEvalRep, ...)` short circuit.
pub fn offset_surface_eval_d1(
    the_basis: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_u: f64,
    the_v: f64,
) -> Option<ResD1> {
    evaluate_d1(the_u, the_v, the_basis, the_offset, the_osc_query)
}

/// OCCT `Geom_OffsetSurface::Offset()`-driven payload convenience: the
/// rcad `Geom_OffsetSurface`-equivalent evaluation of a rcad
/// [`OffsetSurface`] payload (the `offset_basis_and_value` unwrap + the
/// `offset_surface_osculating` member rebuild).
pub fn offset_payload_osculating(of: &OffsetSurface) -> Option<OsculatingSurface> {
    let (a_checking_surf, _offset_value) =
        offset_basis_and_value(of.basis.as_ref(), of.offset_distance);
    offset_surface_osculating(&a_checking_surf)
}
