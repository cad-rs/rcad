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
