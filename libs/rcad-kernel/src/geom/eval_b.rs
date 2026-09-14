//! OCCT BSplSLib (TKMath/BSplSLib) — the surface evaluation leaves
//! `BSplSLib::D0` / `D1` / `D2` / `DN` (BSplSLib.cxx L737-1605).
//!
//! These four bodies are the per-type kernel derivative engine shared by
//! `Geom_BSplineSurface::EvalD0/D1/D2/D3/DN` (Geom_BSplineSurface_1.cxx
//! L111-313) and `Geom_BezierSurface::EvalD0/D1/D2/D3/DN` (Geom_BezierSurface.cxx
//! L1416-1724), and therefore also by the `GeomAdaptor_Surface::EvalDN` arms
//! (GeomAdaptor_Surface.cxx L1731-1752 for BSpline, L1807-1813 for Bezier).
//!
//! The chain is the OCCT one, leaf by leaf:
//!   - `BSplSLib::PrepareEval` (BSplSLib.cxx L313-728) — the local
//!     `(d1+1) x (d2+1)` pole block and the two `2*Degree` knot windows of the
//!     current span; rcad re-host in
//!     [`crate::geom::osculating_surface::prepare_eval`];
//!   - `BSplCLib::Bohm` (BSplCLib.cxx L1197-1541) — the in-place Taylor
//!     scheme producing the derivatives along one direction; rcad re-host in
//!     [`crate::geom::osculating_surface::bspl_clib_bohm`];
//!   - `BSplCLib::Eval` (BSplCLib.cxx L865-1006) — the in-place de Boor corner
//!     cutting; rcad
//!     [`crate::math::bspl_lib::bspl_clib_eval_inplace`];
//!   - `BSplSLib::RationalDerivative` (BSplSLib.cxx L87-300) — the
//!     tensor-product rational quotient rule; rcad
//!     [`crate::math::bspl_lib::bspl_slib_rational_derivative`].
//!
//! Architecture differences (rcad value model):
//!   - `NCollection_Array2<gp_Pnt> Poles` / `NCollection_Array2<double>
//!     Weights` map to the rcad `[Vec<DVec3>]` / `[Vec<f64>]` grids
//!     (`[u_index][v_index]`), `BSplSLib::NoWeights()` to `None`;
//!   - `BSplSLib_DataContainer::ders` (`double ders[48]`) is the local
//!     derivative table; here it is a stack-free `Vec<f64>` sized
//!     `(N+1) * (M+1) * 3` for the `All == true` calls;
//!   - `Geom_UndefinedDerivative` (`Nu + Nv < 1 || Nu < 0 || Nv < 0`) belongs
//!     to the `Geom_*` WRAPPERS, not to the `BSplSLib` leaves re-hosted here:
//!     `BSplSLib::DN` (BSplSLib.cxx L1519) carries no order check at all.
//!     The rcad wrappers DO reproduce the throw, as an assert carrying the same
//!     message (the `Trimmed` arm in `eval.rs`, plus `extrusion_utils.rs`,
//!     `revolution_utils.rs`, `eval_c.rs`, `curve_dn.rs`).  An earlier version
//!     of this note claimed the guard "has no rcad counterpart ... is not
//!     reproduced" — that was false.  OCCT also places the guard on OPPOSITE
//!     sides of the eval-rep short circuit: `Geom_BSplineSurface::EvalDN`
//!     (Geom_BSplineSurface_1.cxx L279-285) checks BEFORE it, while
//!     `Geom_BezierSurface::EvalDN` (Geom_BezierSurface.cxx L1674-1677) checks
//!     AFTER — so a shared guard must not be hoisted above both;
//!   - `Geom_BezierSurface` hands `BSplSLib` the COMPACT knot form
//!     (`UKnots()` = {0, 1}, `UMultiplicities()` = {Degree+1, Degree+1},
//!     Geom_BezierSurface.cxx L2098-2160).  The rcad
//!     [`crate::geom::osculating_surface::prepare_eval`] re-host carries only
//!     the OCCT `Mults == nullptr` arm (its own note: the `Mults != nullptr`
//!     arm is unreachable from Geom_OsculatingSurface), so the Bezier arms
//!     here expand that compact form to the equivalent clamped FLAT knot
//!     sequence ({0} x (Degree+1) ++ {1} x (Degree+1), i.e.
//!     `BSplCLib::KnotSequence` of that pair).  A clamped BSpline is exactly
//!     the Bezier, and for the single span the flat arm's local window and
//!     pole block coincide with the compact arm's, so the two paths feed
//!     `BSplCLib::Bohm` bit-identical data.

use glam::DVec3;

/// OCCT `Geom_BezierSurface::UKnots()` / `UMultiplicities()` (Geom_BezierSurface.cxx
/// L2098-2160) — the compact Bezier knot form expanded to the equivalent
/// clamped flat knot sequence `{0} x (degree+1) ++ {1} x (degree+1)`
/// (`BSplCLib::KnotSequence` of `({0,1}, {degree+1, degree+1})`).
///
/// The rcad `BezierSurface` value carries no knot vector, so the Bezier
/// entries of this module rebuild it; see the module header for why the flat
/// form is the faithful re-host.
pub(crate) fn bezier_flat_knots(degree: usize) -> Vec<f64> {
    let mut knots = vec![0.0f64; degree + 1];
    knots.resize(2 * (degree + 1), 1.0f64);
    knots
}

/// OCCT `BSplSLib::D0` (BSplSLib.cxx L737-778) over the
/// [`crate::geom::osculating_surface::prepare_eval`] re-host — the point of a
/// BSpline/Bezier surface at `(u, v)` with the rational reduction
/// `P = N / W` (the `HomogeneousD0` numerator, BSplSLib.cxx L782-851).
///
/// `weights` is the OCCT `Weights` pointer: `None` is `BSplSLib::NoWeights()`
/// (`Geom_BSplineSurface::Weights()`/`Geom_BezierSurface` pass it exactly when
/// the surface is non-rational).
#[allow(clippy::too_many_arguments)]
pub(crate) fn bspl_slib_d0(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
) -> DVec3 {
    let pe = crate::geom::osculating_surface::prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots,
    );
    let dim: usize = if pe.rational { 4 } else { 3 };
    let mut poles = pe.dc.poles;
    // OCCT L835-836 (rational) / L845-846 (non-rational): both arms are the
    // same two in-place `BSplCLib::Eval` passes, the rational one carrying
    // `dim = 4` (the weight in the last slot).
    crate::math::bspl_lib::bspl_clib_eval_inplace(
        pe.u1,
        pe.d1,
        &pe.dc.knots1,
        dim * (pe.d2 + 1) as usize,
        &mut poles,
    );
    crate::math::bspl_lib::bspl_clib_eval_inplace(
        pe.u2,
        pe.d2,
        &pe.dc.knots2,
        dim,
        &mut poles,
    );
    if pe.rational {
        let w = poles[3];
        DVec3::new(poles[0] / w, poles[1] / w, poles[2] / w)
    } else {
        DVec3::new(poles[0], poles[1], poles[2])
    }
}

/// OCCT `BSplSLib::D1` (BSplSLib.cxx L855-970) — point + first partials.
///
/// `All == true` is the OCCT default of `BSplSLib::RationalDerivative`
/// (BSplSLib.hxx L171-177), which the D1/D2/D3 calls rely on by omitting the
/// last argument (L913, L1129, L1325): it makes the rational branch produce
/// the WHOLE `(N+1) * (M+1) * 3` table, whose `(1,0)` / `(0,1)` triples are
/// then read at the arm-dependent offsets of L915-916 / L942-943.
#[allow(clippy::too_many_arguments)]
pub(crate) fn bspl_slib_d1(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
) -> (DVec3, DVec3, DVec3) {
    let pe = crate::geom::osculating_surface::prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots,
    );
    let u_first = pe.flag_u_or_v;
    let d1 = pe.d1;
    let d2 = pe.d2;
    let dim: usize = if pe.rational { 4 } else { 3 };
    // OCCT L909 (rational) / L922 (non-rational): dim2 = dim * (d2 + 1).
    let dim2 = dim * (d2 + 1) as usize;

    let mut poles = pe.dc.poles;
    crate::geom::osculating_surface::bspl_clib_bohm(
        pe.u1, d1, 1, &pe.dc.knots1, dim2, &mut poles,
    );
    crate::geom::osculating_surface::bspl_clib_bohm(
        pe.u2, d2, 1, &pe.dc.knots2, dim, &mut poles,
    );
    // OCCT L912 / L925: Eval(u2, d2, knots2, dim, poles + dim2).
    crate::math::bspl_lib::bspl_clib_eval_inplace(
        pe.u2,
        d2,
        &pe.dc.knots2,
        dim,
        &mut poles[dim2..(dim2 << 1)],
    );

    if pe.rational {
        // OCCT L913-916: RationalDerivative(d1, d2, 1, 1, poles, ders) — the
        // `M1 = 2` table; the index meaning follows the first direction.
        let mut ders = vec![0.0f64; 2 * 2 * 3];
        crate::math::bspl_lib::bspl_slib_rational_derivative(
            d1, d2, 1, 1, &poles, &mut ders, true,
        );
        if u_first {
            (triple(&ders, 0), triple(&ders, 6), triple(&ders, 3))
        } else {
            (triple(&ders, 0), triple(&ders, 3), triple(&ders, 6))
        }
    } else {
        // OCCT L926-928 / L953-955: result = poles, resVu = poles + dim2,
        // resVv = poles + 3 (or swapped when V is the first direction).
        if u_first {
            (triple(&poles, 0), triple(&poles, dim2), triple(&poles, 3))
        } else {
            (triple(&poles, 0), triple(&poles, 3), triple(&poles, dim2))
        }
    }
}

/// OCCT `BSplSLib::D2` (BSplSLib.cxx L1063-1247) — the six second-order
/// values, in the rcad `SurfaceEval::derivatives2` order
/// `(P, dP/du, dP/dv, d2P/du2, d2P/dudv, d2P/dv2)` (OCCT returns
/// `Point, D1U, D1V, D2U, D2V, D2UV`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn bspl_slib_d2(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
    let pe = crate::geom::osculating_surface::prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots,
    );
    let u_first = pe.flag_u_or_v;
    let d1 = pe.d1;
    let d2 = pe.d2;
    let dim: usize = if pe.rational { 4 } else { 3 };
    let dim2 = dim * (d2 + 1) as usize;

    let mut poles = pe.dc.poles;
    // OCCT L1122-1128 (rational) / L1142-1148 (non-rational).
    crate::geom::osculating_surface::bspl_clib_bohm(
        pe.u1, d1, 2, &pe.dc.knots1, dim2, &mut poles,
    );
    crate::geom::osculating_surface::bspl_clib_bohm(
        pe.u2, d2, 2, &pe.dc.knots2, dim, &mut poles,
    );
    crate::geom::osculating_surface::bspl_clib_bohm(
        pe.u2,
        d2,
        1,
        &pe.dc.knots2,
        dim,
        &mut poles[dim2..(dim2 << 1)],
    );
    if d1 > 1 {
        // OCCT L1127 / L1147: Eval(u2, d2, knots2, dim, poles + 2*dim2).
        crate::math::bspl_lib::bspl_clib_eval_inplace(
            pe.u2,
            d2,
            &pe.dc.knots2,
            dim,
            &mut poles[(dim2 << 1)..(dim2 * 3)],
        );
    }

    if pe.rational {
        // OCCT L1129-1135 / L1184-1190: RationalDerivative(d1, d2, 2, 2) —
        // the `M1 = 3` table, offsets per first direction.
        let mut ders = vec![0.0f64; 3 * 3 * 3];
        crate::math::bspl_lib::bspl_slib_rational_derivative(
            d1, d2, 2, 2, &poles, &mut ders, true,
        );
        if u_first {
            (
                triple(&ders, 0),
                triple(&ders, 9),
                triple(&ders, 3),
                triple(&ders, 18),
                triple(&ders, 12),
                triple(&ders, 6),
            )
        } else {
            (
                triple(&ders, 0),
                triple(&ders, 3),
                triple(&ders, 9),
                triple(&ders, 6),
                triple(&ders, 12),
                triple(&ders, 18),
            )
        }
    } else {
        // OCCT L1149-1168 / L1204-1223: the out-of-degree second derivatives
        // collapse to the static zero vector `BSplSLib_zero`.
        let zero = DVec3::ZERO;
        if u_first {
            let v_uu = if u_degree <= 1 { zero } else { triple(&poles, dim2 << 1) };
            let v_vv = if v_degree <= 1 { zero } else { triple(&poles, 6) };
            let v_uv = triple(&poles, ((d2 << 1) + d2 + 6) as usize);
            (
                triple(&poles, 0),
                triple(&poles, dim2),
                triple(&poles, 3),
                v_uu,
                v_uv,
                v_vv,
            )
        } else {
            let v_uu = if u_degree <= 1 { zero } else { triple(&poles, 6) };
            let v_vv = if v_degree <= 1 { zero } else { triple(&poles, dim2 << 1) };
            let v_uv = triple(&poles, ((d2 << 1) + d2 + 6) as usize);
            (
                triple(&poles, 0),
                triple(&poles, 3),
                triple(&poles, dim2),
                v_uu,
                v_uv,
                v_vv,
            )
        }
    }
}

/// OCCT `BSplSLib::DN` (BSplSLib.cxx L1519-1605) — the `(Nu, Nv)` derivative
/// of a BSpline/Bezier surface at `(U, V)`.
///
/// Unlike D1/D2/D3, the OCCT DN call passes `All = false` explicitly
/// (L1594), so the rational branch yields only the requested `(n1, n2)` triple
/// (indices swapped when V is the first direction).
#[allow(clippy::too_many_arguments)]
pub(crate) fn bspl_slib_dn(
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
) -> DVec3 {
    let pe = crate::geom::osculating_surface::prepare_eval(
        u, v, u_index, v_index, u_degree, v_degree, u_rat, v_rat, u_per, v_per, poles, weights,
        u_knots, v_knots,
    );
    let dim: usize = if pe.rational { 4 } else { 3 };
    let u_first = pe.flag_u_or_v;
    let d1 = pe.d1;
    let d2 = pe.d2;
    // OCCT L1581-1582: n1 is the request of the first direction.
    let n1 = if u_first { nu } else { nv };
    let n2 = if u_first { nv } else { nu };

    // OCCT L1570-1579: the non-rational out-of-degree request returns zero.
    if !pe.rational && (nu > u_degree || nv > v_degree) {
        return DVec3::ZERO;
    }

    let mut poles = pe.dc.poles;
    let row = dim * (d2 + 1) as usize;
    // OCCT L1584: Bohm(u1, d1, n1, knots1, dim * (d2 + 1), poles).
    crate::geom::osculating_surface::bspl_clib_bohm(pe.u1, d1, n1, &pe.dc.knots1, row, &mut poles);
    // OCCT L1586-1589: the V window of every U-derivative row.
    for k in 0..=n1.min(d1) {
        let off = (k * dim as i32 * (d2 + 1)) as usize;
        let end = off + row;
        crate::geom::osculating_surface::bspl_clib_bohm(
            pe.u2,
            d2,
            n2,
            &pe.dc.knots2,
            dim,
            &mut poles[off..end],
        );
    }

    if pe.rational {
        // OCCT L1592-1596: RationalDerivative(d1, d2, n1, n2, poles, ders,
        // false) -> result = ders.
        let mut ders = vec![0.0f64; 3];
        crate::math::bspl_lib::bspl_slib_rational_derivative(
            d1, d2, n1, n2, &poles, &mut ders, false,
        );
        DVec3::new(ders[0], ders[1], ders[2])
    } else {
        // OCCT L1599: result = poles + (n1 * (d2 + 1) + n2) * dim.
        triple(&poles, ((n1 * (d2 + 1) + n2) * dim as i32) as usize)
    }
}

/// The 3 consecutive doubles of an OCCT `gp_Vec`/`gp_Pnt` slot.
fn triple(values: &[f64], index: usize) -> DVec3 {
    DVec3::new(values[index], values[index + 1], values[index + 2])
}
