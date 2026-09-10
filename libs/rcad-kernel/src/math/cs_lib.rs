//! CSLib — surface normal computation utilities.
//!
//! Analogous to OCCT `CSLib` package in TKMath.
//! Provides functions for computing surface normals and their derivatives,
//! including handling of singular/degenerate cases where the standard
//! cross product D1U × D1V vanishes.
//!
//! OCCT reference: CSLib.hxx (FoundationClasses/TKMath/CSLib/)

#![allow(unused_variables)]

use glam::DVec3;

/// Status returned by basic normal computation from first derivatives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivativeStatus {
    /// Normal was successfully computed (D1U × D1V is non-zero and within tolerance).
    Done,
    /// D1U has zero length.
    D1UIsNull,
    /// D1V has zero length.
    D1VIsNull,
    /// D1U and D1V are parallel (cross product below sine tolerance).
    D1UD1VAreParallel,
}

/// Status returned by singular normal computation using second derivatives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalStatus {
    /// Normal is defined.
    Defined,
    /// Normal is undefined (all derivatives vanish).
    Undefined,
    /// Normal is uncertain (singular point, but a direction could be estimated).
    Singular,
    /// OCCT CSLib_InfinityOfSolutions — the first non-zero derivative of the
    /// normal changes sign in the domain, so the normal direction is not
    /// uniquely determined (used by the max-order overload).
    InfinityOfSolutions,
}

// =============================================================================
// Primary: Normal from first derivatives
// =============================================================================

/// Compute the normal direction of a surface from first partial derivatives.
///
/// The normal is `D1U × D1V` (cross product), normalized.
///
/// Returns `None` in the status when:
/// - D1U has null length
/// - D1V has null length
/// - D1U and D1V are parallel (cross product magnitude < the_sin_tol)
///
/// OCCT-aligned: CSLib::Normal(D1U, D1V, SinTol, Status, Normal)
pub fn normal_from_derivatives(
    d1u: DVec3,
    d1v: DVec3,
    sin_tol: f64,
) -> (Option<DVec3>, DerivativeStatus) {
    let len_u_sq = d1u.length_squared();
    let len_v_sq = d1v.length_squared();

    if len_u_sq < 1e-30 {
        return (None, DerivativeStatus::D1UIsNull);
    }
    if len_v_sq < 1e-30 {
        return (None, DerivativeStatus::D1VIsNull);
    }

    let cross = d1u.cross(d1v);
    let cross_len_sq = cross.length_squared();
    let sin_angle_sq = cross_len_sq / (len_u_sq * len_v_sq);

    if sin_angle_sq < sin_tol * sin_tol {
        // D1U and D1V are parallel (or nearly so)
        return (None, DerivativeStatus::D1UD1VAreParallel);
    }

    let normal = cross / cross_len_sq.sqrt();
    (Some(normal), DerivativeStatus::Done)
}

// =============================================================================
// Singular: Normal from first + second derivatives
// =============================================================================

/// Compute an approximate normal direction at a singular point where
/// the first derivatives are parallel or zero, using second derivatives.
///
/// Uses a limited Taylor expansion of the non-normalized normal
/// `N = D1U × D1V`. When N(u0,v0) is zero (cross product vanishes),
/// the leading terms   dN/du * du + dN/dv * dv   are examined to find
/// a non-zero direction.
///
/// - `d1u`, `d1v` — First partial derivatives at the point.
/// - `d2u`, `d2v`, `d2uv` — Second partial derivatives (d²S/du², d²S/dv², d²S/dudv).
/// - `sin_tol` — Sine tolerance for parallelism checks.
///
/// Returns `(Some(normal), NormalStatus::Singular)` when a rescue direction
/// is found, `(None, NormalStatus::Undefined)` when no normal can be
/// determined, or `(None, NormalStatus::Singular)` when the solution is
/// ambiguous (OCCT InfinityOfSolutions).
///
/// OCCT-aligned: CSLib::Normal(D1U, D1V, D2U, D2V, D2UV, SinTol, Done, Status, Normal)
pub fn normal_from_derivatives_with_hessian(
    d1u: DVec3,
    d1v: DVec3,
    d2u: DVec3,
    d2v: DVec3,
    d2uv: DVec3,
    sin_tol: f64,
) -> (Option<DVec3>, NormalStatus) {
    // OCCT CSLib.cxx: Normal(D1U, D1V, D2U, D2V, D2UV, SinTol, Done, Status, Normal).
    // dN/du = D2U ^ D1V + D1U ^ D2UV,  dN/dv = D2UV ^ D1V + D1U ^ D2V.
    let d1nu = d2u.cross(d1v) + d1u.cross(d2uv);
    let d1nv = d2uv.cross(d1v) + d1u.cross(d2v);

    let l_d1nu = d1nu.length_squared();
    let l_d1nv = d1nv.length_squared();

    let eps = f64::EPSILON; // OCCT RealEpsilon()

    if l_d1nu <= eps && l_d1nv <= eps {
        return (None, NormalStatus::Undefined); // D1NIsNull, Done=false
    }
    if l_d1nu < eps {
        return (Some(d1nv.normalize_or_zero()), NormalStatus::Singular); // D1NuIsNull
    }
    if l_d1nv < eps {
        return (Some(d1nu.normalize_or_zero()), NormalStatus::Singular); // D1NvIsNull
    }
    if (l_d1nv / l_d1nu) <= eps {
        return (None, NormalStatus::Undefined); // D1NvNuRatioIsNull, Done=false
    }
    if (l_d1nu / l_d1nv) <= eps {
        return (None, NormalStatus::Undefined); // D1NuNvRatioIsNull, Done=false
    }

    let d1n_cross = d1nu.cross(d1nv);
    let sin2 = d1n_cross.length_squared() / (l_d1nu * l_d1nv);

    if sin2 < sin_tol * sin_tol {
        (Some(d1nu.normalize_or_zero()), NormalStatus::Singular) // D1NuIsParallelD1Nv
    } else {
        (None, NormalStatus::Singular) // InfinityOfSolutions, Done=false
    }
}

// =============================================================================
// Simplified: Normal using magnitude tolerance
// =============================================================================

/// Compute surface normal using a simpler magnitude-based tolerance check.
///
/// If `|D1U × D1V| >= mag_tol` and both `|D1U| >= mag_tol` and `|D1V| >= mag_tol`,
/// the normal is defined.
///
/// OCCT-aligned: CSLib::Normal(D1U, D1V, MagTol, Status, Normal)
pub fn normal_from_derivatives_mag(
    d1u: DVec3,
    d1v: DVec3,
    mag_tol: f64,
) -> (Option<DVec3>, NormalStatus) {
    let len_u = d1u.length();
    let len_v = d1v.length();

    if len_u < mag_tol || len_v < mag_tol {
        return (None, NormalStatus::Undefined);
    }

    let cross = d1u.cross(d1v);
    let cross_len = cross.length();

    if cross_len < mag_tol {
        return (None, NormalStatus::Singular);
    }

    (Some(cross / cross_len), NormalStatus::Defined)
}

// =============================================================================
// Derivative helpers
// =============================================================================

/// Compute the derivative of order (nu, nv) of the non-normalized normal vector
/// `N = dS/du × dS/dv`.
///
/// `der_surf` is a flat slice containing surface derivative vectors at the
/// required orders: `der_surf[i * (nv+2) + j]` = d^(i+j)S / (du^i · dv^j)
/// for `i = 0..=nu+1`, `j = 0..=nv+1`.
///
/// OCCT-aligned: CSLib::DNNUV(Nu, Nv, theDerSurf)
pub fn dnnuv(nu: usize, nv: usize, der_surf: &[DVec3], stride: usize) -> DVec3 {
    // N = Su × Sv
    // d^(nu+nv)N / (du^nu · dv^nv) = sum over k=0..nu, l=0..nv of
    //   C(nu,k) * C(nv,l) * (d^(k+l+1)S/(du^(k+1)·dv^l)) × (d^(nu+nv-k-l)S/(du^(nu-k)·dv^(nv-l)))
    // where (k,l) gives the order on Su and (nu-k, nv-l) gives the order on Sv.

    // Precompute binomial coefficients up to max(nu+1, nv+1)
    let max_n = (nu + 1).max(nv + 1);
    let mut binom = vec![vec![0i64; max_n + 1]; max_n + 1];
    for n in 0..=max_n {
        binom[n][0] = 1;
        binom[n][n] = 1;
        for k in 1..n {
            binom[n][k] = binom[n - 1][k - 1] + binom[n - 1][k];
        }
    }

    let mut result = DVec3::ZERO;
    for k in 0..=nu {
        for l in 0..=nv {
            let su_idx = (k + 1) * stride + l; // d^(k+l+1)S/(du^(k+1)·dv^l)
            let sv_idx = (nu - k) * stride + (nv - l); // d^(nu+nv-k-l)S/(du^(nu-k)·dv^(nv-l))
            let coeff = (binom[nu][k] * binom[nv][l]) as f64;
            result += coeff * der_surf[su_idx].cross(der_surf[sv_idx]);
        }
    }
    result
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_basic_orthogonal() {
        let d1u = DVec3::X;
        let d1v = DVec3::Y;
        let (n, status) = normal_from_derivatives(d1u, d1v, 1e-9);
        assert_eq!(status, DerivativeStatus::Done);
        assert!(n.is_some());
        assert!((n.unwrap() - DVec3::Z).length() < 1e-10);
    }

    #[test]
    fn normal_parallel_derivatives() {
        let d1u = DVec3::X;
        let d1v = DVec3::X; // parallel
        let (n, status) = normal_from_derivatives(d1u, d1v, 1e-9);
        assert_eq!(status, DerivativeStatus::D1UD1VAreParallel);
        assert!(n.is_none());
    }

    #[test]
    fn normal_null_u() {
        let (n, status) = normal_from_derivatives(DVec3::ZERO, DVec3::Y, 1e-9);
        assert_eq!(status, DerivativeStatus::D1UIsNull);
        assert!(n.is_none());
    }

    #[test]
    fn normal_null_v() {
        let (n, status) = normal_from_derivatives(DVec3::X, DVec3::ZERO, 1e-9);
        assert_eq!(status, DerivativeStatus::D1VIsNull);
        assert!(n.is_none());
    }

    #[test]
    fn normal_hessian_rescue_singular() {
        // Singular point: D1U and D1V are parallel (not independent) so the
        // first-derivative normal vanishes; D2UV = Y rescues it.
        // OCCT CSLib::Normal: dN/du = D2U ^ D1V + D1U ^ D2UV = X ^ Y = Z,
        // so the approximate normal should be near Z.
        let d1u = DVec3::X;
        let d1v = DVec3::X; // parallel
        let d2u = DVec3::ZERO;
        let d2v = DVec3::ZERO;
        let d2uv = DVec3::Y; // d²S/dudv = Y
        let (n, _status) = normal_from_derivatives_with_hessian(
            d1u, d1v, d2u, d2v, d2uv, 1e-9,
        );
        assert!(n.is_some(), "expected approximate normal at singular point");
        let dir = n.unwrap();
        assert!((dir - DVec3::Z).length() < 1e-6, "normal should be near Z, got {dir:?}");
    }

    #[test]
    fn normal_mag_tol() {
        let d1u = DVec3::new(1.0, 0.0, 0.0);
        let d1v = DVec3::new(0.0, 1e-10, 0.0); // very small
        let (n, status) = normal_from_derivatives_mag(d1u, d1v, 1e-7);
        // d1v is shorter than mag_tol → undefined
        assert_eq!(status, NormalStatus::Undefined);
        assert!(n.is_none());
    }
}

// =============================================================================
// Array2 forms used by BlendFunc::ComputeNormal / ComputeDNormal
// (OCCT CSLib.cxx L183-560)
// =============================================================================

/// OCCT CSLib_NormalPolyDef (CSLib_NormalPolyDef.hxx/cxx) — the polynomial
/// F(X) = Sum_{i=0}^{k0} C(k0,i) * cos^i(X) * sin^(k0-i)(X) * li(i) used by
/// the max-order Normal overload to detect sign changes.
struct NormalPolyDef {
    k0: i32,
    tab_li: Vec<f64>,
}

impl NormalPolyDef {
    /// OCCT CSLib_NormalPolyDef(theK0, theLi) (CSLib_NormalPolyDef.cxx L23-31).
    fn new(k0: i32, li: &[f64]) -> Self {
        NormalPolyDef {
            k0,
            tab_li: li.to_vec(),
        }
    }

    /// OCCT Value(theX, theF) (CSLib_NormalPolyDef.cxx L33-54).
    fn value(&mut self, the_x: f64, the_f: &mut f64) -> bool {
        *the_f = 0.0;

        let a_cos = the_x.cos();
        let a_sin = the_x.sin();

        // At singular points (cos or sin near zero), return zero to avoid
        // numerical instability.  OCCT: RealSmall().
        if a_cos.abs() <= f64::MIN_POSITIVE || a_sin.abs() <= f64::MIN_POSITIVE {
            return true;
        }

        // F(X) = Sum_{i=0}^{k0} C(k0,i) * cos^i(X) * sin^(k0-i)(X) * li(i)
        for i in 0..=self.k0 {
            *the_f += crate::math::plib::binomial(self.k0 as usize, i as usize)
                * a_cos.powi(i)
                * a_sin.powi(self.k0 - i)
                * self.tab_li[i as usize];
        }

        true
    }
}

impl crate::math::root::function_all_roots::FunctionValue for NormalPolyDef {
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut f = 0.0;
        if NormalPolyDef::value(self, x, &mut f) {
            Some(f)
        } else {
            None
        }
    }
}

impl crate::math::root::function_all_roots::FunctionWithDerivative for NormalPolyDef {
    /// OCCT Derivative(theX, theD) (CSLib_NormalPolyDef.cxx L56-76).
    fn derivative(&mut self, the_x: f64) -> Option<f64> {
        let mut the_d = 0.0;

        let a_cos = the_x.cos();
        let a_sin = the_x.sin();

        if a_cos.abs() <= f64::MIN_POSITIVE || a_sin.abs() <= f64::MIN_POSITIVE {
            return Some(the_d);
        }

        // dF/dX = Sum C(k0,i) * cos^(i-1)(X) * sin^(k0-i-1)(X)
        //                              * (k0*cos^2(X) - i) * li(i)
        for i in 0..=self.k0 {
            the_d += crate::math::plib::binomial(self.k0 as usize, i as usize)
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (self.k0 as f64 * a_cos * a_cos - i as f64)
                * self.tab_li[i as usize];
        }

        Some(the_d)
    }

    /// OCCT Values(theX, theF, theD) (CSLib_NormalPolyDef.cxx L78-105).
    fn values(&mut self, the_x: f64) -> Option<(f64, f64)> {
        let mut the_f = 0.0;
        let mut the_d = 0.0;

        let a_cos = the_x.cos();
        let a_sin = the_x.sin();

        if a_cos.abs() <= f64::MIN_POSITIVE || a_sin.abs() <= f64::MIN_POSITIVE {
            return Some((the_f, the_d));
        }

        for i in 0..=self.k0 {
            let a_bin_coeff = crate::math::plib::binomial(self.k0 as usize, i as usize);
            let a_li_coeff = self.tab_li[i as usize];

            the_f += a_bin_coeff * a_cos.powi(i) * a_sin.powi(self.k0 - i) * a_li_coeff;
            the_d += a_bin_coeff
                * a_cos.powi(i - 1)
                * a_sin.powi(self.k0 - i - 1)
                * (self.k0 as f64 * a_cos * a_cos - i as f64)
                * a_li_coeff;
        }

        Some((the_f, the_d))
    }
}

/// OCCT gp_Vec::IsParallel(Other, AngularTolerance) — cross-product test.
fn vectors_are_parallel(a: DVec3, b: DVec3, angular_tol: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la <= f64::EPSILON || lb <= f64::EPSILON {
        return true;
    }
    let sin_angle = a.cross(b).length() / (la * lb);
    sin_angle <= angular_tol
}

/// OCCT CSLib::Normal(theMaxOrder, theDerNUV, theSinTol, theU, theV, theUmin,
/// theUmax, theVmin, theVmax, theStatus, theNormal, theOrderU, theOrderV)
/// (CSLib.cxx L183-371).  Returns (status, normal, order_u, order_v).
#[allow(clippy::too_many_arguments)]
pub fn normal_max_order(
    max_order: i32,
    der_nuv: &[Vec<DVec3>],
    sin_tol: f64,
    the_u: f64,
    the_v: f64,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
) -> (NormalStatus, Option<DVec3>, i32, i32) {
    const THE_PARALLEL_ANGULAR_TOL: f64 = 1e-6;
    const THE_MAX_ROOT_ITERATIONS: i32 = 200;
    const THE_ROOT_FINDING_TOL: f64 = 1e-5;

    let mut an_order: i32 = -1;
    let mut a_found_u_idx: i32 = 0;
    let mut a_found = false;
    let mut a_norme: f64 = 0.0;
    let mut a_d = DVec3::ZERO;

    // Find k0 such that all derivatives N = dS/du ^ dS/dv are null till
    // order k0-1.
    while !a_found && an_order < max_order {
        an_order += 1;
        a_found_u_idx = an_order;
        while a_found_u_idx >= 0 && !a_found {
            let a_v_idx = an_order - a_found_u_idx;
            a_d = der_nuv[a_found_u_idx as usize][a_v_idx as usize];
            a_norme = a_d.length();
            a_found = a_norme >= sin_tol;
            a_found_u_idx -= 1;
        }
    }

    let order_u = a_found_u_idx + 1;
    let order_v = an_order - order_u;

    // Vk0 is the first non-null derivative of N: the reference vector.
    if !a_found {
        return (NormalStatus::Singular, None, order_u, order_v);
    }

    if an_order == 0 {
        return (
            NormalStatus::Defined,
            Some(a_d / a_d.length()),
            order_u,
            order_v,
        );
    }

    let a_vk0 = der_nuv[order_u as usize][order_v as usize];
    let mut a_ratio = vec![0.0f64; (an_order + 1) as usize];

    // Calculate lambda_i ratios for each derivative at this order.
    let mut a_ratio_idx: i32 = 0;
    let mut is_defined = false;
    while a_ratio_idx <= an_order && !is_defined {
        let a_der_vec = der_nuv[a_ratio_idx as usize][(an_order - a_ratio_idx) as usize];
        if a_der_vec.length() <= sin_tol {
            a_ratio[a_ratio_idx as usize] = 0.0;
        } else if vectors_are_parallel(a_der_vec, a_vk0, THE_PARALLEL_ANGULAR_TOL) {
            let mut a_magnitude_ratio = a_der_vec.length() / a_vk0.length();
            if a_der_vec.dot(a_vk0) < 0.0 {
                a_magnitude_ratio = -a_magnitude_ratio;
            }
            a_ratio[a_ratio_idx as usize] = a_magnitude_ratio;
        } else {
            is_defined = true;
        }
        a_ratio_idx += 1;
    }

    if is_defined {
        return (
            NormalStatus::Defined,
            Some(a_d / a_d.length()),
            order_u,
            order_v,
        );
    }

    // All lambda_i exist - analyze the polynomial sign.
    let mut a_inf = -std::f64::consts::PI;
    let mut a_sup = std::f64::consts::PI;

    // Determine domain based on position (interior, edge, corner).
    let p_confusion: f64 = 1e-9; // OCCT Precision::PConfusion() = Confusion() * 0.01.
    let is_fu = (the_u - the_umin).abs() < p_confusion;
    let is_lu = (the_u - the_umax).abs() < p_confusion;
    let is_fv = (the_v - the_vmin).abs() < p_confusion;
    let is_lv = (the_v - the_vmax).abs() < p_confusion;

    if is_lu {
        a_inf = std::f64::consts::FRAC_PI_2;
        a_sup = 3.0 * std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_inf = std::f64::consts::PI;
        }
        if is_fv {
            a_sup = std::f64::consts::PI;
        }
    } else if is_fu {
        a_sup = std::f64::consts::FRAC_PI_2;
        a_inf = -std::f64::consts::FRAC_PI_2;
        if is_lv {
            a_sup = 0.0;
        }
        if is_fv {
            a_inf = 0.0;
        }
    } else if is_lv {
        a_inf = -std::f64::consts::PI;
        a_sup = 0.0;
    } else if is_fv {
        a_inf = 0.0;
        a_sup = std::f64::consts::PI;
    }

    let mut a_changes_sign = false;
    let mut a_vprec: f64 = 0.0;
    let mut a_vsuiv: f64 = 0.0;

    // Create polynomial and find its roots (OCCT math_FunctionRoots; K = 0).
    let mut a_poly = NormalPolyDef::new(an_order, &a_ratio);
    let a_find_roots = crate::math::root::function_all_roots::FunctionRoots::new(
        &mut a_poly,
        a_inf,
        a_sup,
        THE_MAX_ROOT_ITERATIONS,
        THE_ROOT_FINDING_TOL,
        1e-7, // Precision::Confusion()
        1e-7, // Precision::Confusion()
        0.0,
    );

    if a_find_roots.is_done() && a_find_roots.nb_solutions() > 0 {
        // Sort roots in ascending order.
        let a_nb_sol = a_find_roots.nb_solutions();
        let mut a_sol = vec![0.0f64; a_nb_sol + 2];

        for (idx, a_root) in a_sol.iter_mut().enumerate().take(a_nb_sol + 1).skip(1) {
            *a_root = a_find_roots.value(idx);
        }
        let roots_slice = &mut a_sol[1..a_nb_sol + 1];
        roots_slice.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Add domain limits.
        a_sol[0] = a_inf;
        a_sol[a_nb_sol + 1] = a_sup;

        // Check for sign changes between consecutive roots.
        let mut a_first = 0usize;
        for a_interval_idx in 0..=a_nb_sol {
            if (a_sol[a_interval_idx + 1] - a_sol[a_interval_idx]).abs() > p_confusion {
                let _ = a_poly.value(
                    (a_sol[a_interval_idx] + a_sol[a_interval_idx + 1]) / 2.0,
                    &mut a_vsuiv,
                );
                if a_first == 0 {
                    a_first = a_interval_idx;
                    a_changes_sign = false;
                    a_vprec = a_vsuiv;
                } else {
                    a_changes_sign = a_changes_sign || (a_vprec * a_vsuiv) < 0.0;
                    a_vprec = a_vsuiv;
                }
            }
        }
    } else {
        // No roots found, polynomial doesn't change sign.
        a_changes_sign = false;
        let _ = a_poly.value(a_inf, &mut a_vsuiv);
    }

    // Determine status based on polynomial sign.
    if a_changes_sign {
        (NormalStatus::InfinityOfSolutions, None, order_u, order_v)
    } else {
        let a_sign = if a_vsuiv > 0.0 { 1.0 } else { -1.0 };
        let n = a_vk0 / a_vk0.length();
        (NormalStatus::Defined, Some(a_sign * n), order_u, order_v)
    }
}

/// OCCT CSLib::DNNUV(theNu, theNv, theDerSurf) (CSLib.cxx L391-408) — the
/// derivative (theNu, theNv) of the non-normalized normal N = dS/du ^ dS/dv
/// from the Array2 of surface derivatives `der_surf[i][j]` = d^(i+j)S/du^i dv^j.
pub fn dnnuv_array2(nu: i32, nv: i32, der_surf: &[Vec<DVec3>]) -> DVec3 {
    let mut a_result = DVec3::ZERO;

    for i in 0..=nu {
        for j in 0..=nv {
            let a_vg = der_surf[(i + 1) as usize][j as usize];
            let a_vd = der_surf[(nu - i) as usize][(nv + 1 - j) as usize];
            let a_cross = a_vg.cross(a_vd);
            let a_bin_coef = crate::math::plib::binomial(nu as usize, i as usize)
                * crate::math::plib::binomial(nv as usize, j as usize);
            a_result += a_bin_coef * a_cross;
        }
    }

    a_result
}

/// OCCT CSLib::DNNormal(theNu, theNv, theDerNUV, theIduref, theIdvref)
/// (CSLib.cxx L436-560) — the derivative (theNu, theNv) of the unit normal
/// vector of direction theDerNUV(theIduref, theIdvref).
pub fn dnnormal(nu: i32, nv: i32, der_nuv: &[Vec<DVec3>], iduref: i32, idvref: i32) -> DVec3 {
    let a_kderiv = nu + nv;

    let mut a_der_vec_nor =
        vec![vec![DVec3::ZERO; (a_kderiv + 1) as usize]; (a_kderiv + 1) as usize];
    let mut a_tab_scal = vec![vec![0.0f64; (a_kderiv + 1) as usize]; (a_kderiv + 1) as usize];
    let mut a_tab_norm = vec![vec![0.0f64; (a_kderiv + 1) as usize]; (a_kderiv + 1) as usize];

    let bin = crate::math::plib::binomial;

    let a_der_nor0 = der_nuv[iduref as usize][idvref as usize];
    let a_der_nor_unit = a_der_nor0 / a_der_nor0.length();
    a_der_vec_nor[0][0] = a_der_nor_unit;

    let a_dnorm0 = a_der_nor0.dot(a_der_vec_nor[0][0]);
    a_tab_norm[0][0] = a_dnorm0;
    a_tab_scal[0][0] = 0.0;

    for a_mderiv in 1..=a_kderiv {
        for a_pderiv in 0..=a_mderiv {
            let a_qderiv = a_mderiv - a_pderiv;
            if a_pderiv > nu || a_qderiv > nv {
                continue;
            }

            // Compute n . derivative(p,q) of n
            let mut a_scal = 0.0;
            if a_pderiv > a_qderiv {
                for a_jderiv in 1..=a_qderiv {
                    a_scal -= bin(a_qderiv as usize, a_jderiv as usize)
                        * a_der_vec_nor[0][a_jderiv as usize].dot(
                            a_der_vec_nor[a_pderiv as usize][(a_qderiv - a_jderiv) as usize],
                        );
                }

                for a_jderiv in 0..a_qderiv {
                    a_scal -= bin(a_qderiv as usize, a_jderiv as usize)
                        * a_der_vec_nor[a_pderiv as usize][a_jderiv as usize]
                            .dot(a_der_vec_nor[0][(a_qderiv - a_jderiv) as usize]);
                }

                for a_ideriv in 1..a_pderiv {
                    for a_jderiv in 0..=a_qderiv {
                        a_scal -= bin(a_pderiv as usize, a_ideriv as usize)
                            * bin(a_qderiv as usize, a_jderiv as usize)
                            * a_der_vec_nor[a_ideriv as usize][a_jderiv as usize].dot(
                                a_der_vec_nor[(a_pderiv - a_ideriv) as usize]
                                    [(a_qderiv - a_jderiv) as usize],
                            );
                    }
                }
            } else {
                for a_ideriv in 1..=a_pderiv {
                    a_scal -= bin(a_pderiv as usize, a_ideriv as usize)
                        * a_der_vec_nor[a_ideriv as usize][0]
                            .dot(a_der_vec_nor[(a_pderiv - a_ideriv) as usize][a_qderiv as usize]);
                }

                for a_ideriv in 0..a_pderiv {
                    a_scal -= bin(a_pderiv as usize, a_ideriv as usize)
                        * a_der_vec_nor[a_ideriv as usize][a_qderiv as usize]
                            .dot(a_der_vec_nor[(a_pderiv - a_ideriv) as usize][0]);
                }

                for a_ideriv in 0..=a_pderiv {
                    for a_jderiv in 1..a_qderiv {
                        a_scal -= bin(a_pderiv as usize, a_ideriv as usize)
                            * bin(a_qderiv as usize, a_jderiv as usize)
                            * a_der_vec_nor[a_ideriv as usize][a_jderiv as usize].dot(
                                a_der_vec_nor[(a_pderiv - a_ideriv) as usize]
                                    [(a_qderiv - a_jderiv) as usize],
                            );
                    }
                }
            }
            a_tab_scal[a_pderiv as usize][a_qderiv as usize] = a_scal / 2.0;

            // Compute the derivative (n,p) of NUV Length.
            let mut a_dnorm = der_nuv[(a_pderiv + iduref) as usize][(a_qderiv + idvref) as usize]
                .dot(a_der_vec_nor[0][0]);

            for a_jderiv in 0..a_qderiv {
                a_dnorm -= bin((a_qderiv + idvref) as usize, (a_jderiv + idvref) as usize)
                    * a_tab_norm[a_pderiv as usize][a_jderiv as usize]
                    * a_tab_scal[0][(a_qderiv - a_jderiv) as usize];
            }

            for a_ideriv in 0..a_pderiv {
                for a_jderiv in 0..=a_qderiv {
                    a_dnorm -= bin((a_pderiv + iduref) as usize, (a_ideriv + iduref) as usize)
                        * bin((a_qderiv + idvref) as usize, (a_jderiv + idvref) as usize)
                        * a_tab_norm[a_ideriv as usize][a_jderiv as usize]
                        * a_tab_scal[(a_pderiv - a_ideriv) as usize]
                            [(a_qderiv - a_jderiv) as usize];
                }
            }
            a_tab_norm[a_pderiv as usize][a_qderiv as usize] = a_dnorm;

            // Compute derivative (p,q) of n.
            let mut a_der_nor =
                der_nuv[(a_pderiv + iduref) as usize][(a_qderiv + idvref) as usize];

            for a_jderiv in 1..=a_qderiv {
                a_der_nor -= bin((a_pderiv + iduref) as usize, iduref as usize)
                    * bin((a_qderiv + idvref) as usize, (a_jderiv + idvref) as usize)
                    * a_tab_norm[0][a_jderiv as usize]
                    * a_der_vec_nor[a_pderiv as usize][(a_qderiv - a_jderiv) as usize];
            }

            for a_ideriv in 1..=a_pderiv {
                for a_jderiv in 0..=a_qderiv {
                    a_der_nor -= bin((a_pderiv + iduref) as usize, (a_ideriv + iduref) as usize)
                        * bin((a_qderiv + idvref) as usize, (a_jderiv + idvref) as usize)
                        * a_tab_norm[a_ideriv as usize][a_jderiv as usize]
                        * a_der_vec_nor[(a_pderiv - a_ideriv) as usize]
                            [(a_qderiv - a_jderiv) as usize];
                }
            }

            a_der_nor /= bin((a_pderiv + iduref) as usize, iduref as usize)
                * bin((a_qderiv + idvref) as usize, idvref as usize)
                * a_tab_norm[0][0];
            a_der_vec_nor[a_pderiv as usize][a_qderiv as usize] = a_der_nor;
        }
    }

    a_der_vec_nor[nu as usize][nv as usize]
}
