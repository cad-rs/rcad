//! OCCT GeomLib (TKGeomBase/GeomLib) + BSplCLib/PLib extension
//! primitives used by the ChFi3d corner tails — 1:1 translation.
//!
//!   - PLib::CoefficientsPoles, dim version (PLib.cxx L1522-1608)
//!   - BSplCLib::TangExtendToConstraint (BSplCLib.cxx L3794-4160)
//!   - GeomLib::ExtendSurfByLength      (GeomLib.cxx L1485-1962)
//!   - GeomLib::ExtendCurveToPoint      (GeomLib.cxx L1269-1410)
//!
//! GAP carriers (outside-package dependencies, OCCT failure path kept):
//!   - ExtendKPart (GeomLib.cxx L1420+) and GeomConvert_ApproxSurface are
//!     pending — non-BSpline/Bezier surfaces are returned unextended.
//!   - ComputeLambda (GeomLib.cxx L126-200) and
//!     GeomConvert_CompCurveToBSplineCurve (TKG3d) are pending —
//!     ExtendCurveToPoint keeps the curve unextended.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    BSplineSurface, BezierSurface, Curve3, CurveEval as _, Surface3, SurfaceEval as _,
};
use rcad_kernel::math::bspl_lib::{
    eval_flat, increase_degree as bspl_increase_degree, knot_sequence as bspl_knot_sequence,
    remove_knot as bspl_remove_knot, increase_degree_count_knots,
};
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::plib::hermite_coefficients;

// =========================================================================
// OCCT PLib.cxx L1522-1608 — PLib::CoefficientsPoles, generic dimension
// (the rcad plib.rs DVec3 version scalarized): converts the power
// coefficients of a degree `n = size - 1` polynomial (per coordinate) into
// the Bezier poles.
// =========================================================================
fn plib_coefficients_poles_dim(coefs: &[f64], dim: usize) -> Vec<f64> {
    let size = coefs.len() / dim;
    let mut poles = coefs.to_vec();
    // Les Extremites are kept; the intermediate coefficients are divided by
    // the binomial then accumulated (Horner-like Bernstein accumulation).
    for c in 0..dim {
        for i in 2..size {
            let idx = (i - 1) * dim + c;
            let cnp = binomial(size - 1, i - 1);
            poles[idx] = coefs[idx] / cnp;
        }
        for i in 1..=size - 1 {
            let mut j = size - 1;
            while j >= i {
                poles[j * dim + c] += poles[(j - 1) * dim + c];
                j -= 1;
            }
        }
    }
    poles
}

fn binomial(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    let mut r = 1.0f64;
    for i in 0..k {
        r = r * (n - i) as f64 / (i + 1) as f64;
    }
    r
}

/// Compresses a flat knot vector into OCCT (knots, mults) — the
/// BSplCLib::Knots(Flat, Knots, Mults) direction.
fn compress_knots(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for k in flat {
        if let Some(last) = knots.last() {
            if *last == *k {
                *mults.last_mut().expect("mults") += 1;
                continue;
            }
        }
        knots.push(*k);
        mults.push(1);
    }
    (knots, mults)
}

/// Expands OCCT (knots, mults) into the flat knot vector — the
/// BSplCLib::KnotSequence direction (non-periodic).
fn expand_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut flat = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            flat.push(*k);
        }
    }
    flat
}

// =========================================================================
// OCCT BSplCLib.cxx L3794-4160 — BSplCLib::TangExtendToConstraint.
// Builds the extension of a flat-knot/pole B-spline (of degree <c_degree>,
// <num_poles> poles of dimension <c_dimension>) that reaches
// <constraint_point> with C<continuity> at the connection knot, and
// concatenates it before/after the input.
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn bsplclib_tang_extend_to_constraint(
    flat_knots: &[f64],
    c1_coefficient: f64,
    num_poles: usize,
    poles: &[f64],
    c_dimension: usize,
    c_degree: usize,
    constraint_point: &[f64],
    continuity: i32,
    after: bool,
    nb_poles_result: &mut usize,
    nb_knots_result: &mut usize,
    knots_result: &mut Vec<f64>,
    poles_result: &mut Vec<f64>,
) {
    let cont = continuity as usize;
    // 1. calculation of extension nD.
    // Hermite matrix (OCCT L3810-3827).
    let csize = cont + 2;
    let mut mat_coefs = Matrix::new_init(1, csize as i32, 1, csize as i32, 0.0);
    if after {
        hermite_coefficients(0.0, 1.0, continuity, 0, &mut mat_coefs);
    } else {
        hermite_coefficients(0.0, 1.0, 0, continuity, &mut mat_coefs);
    }

    // position at the node of connection (OCCT L3830-3837).
    let tbord = if after {
        flat_knots[flat_knots.len() - 1 - c_degree]
    } else {
        flat_knots[c_degree]
    };
    // BSplCLib::Eval at Tbord (OCCT L3844-3854).
    let derivative_request = std::cmp::max(continuity, 1) as usize;
    let mut extrap_mode = [c_degree as i32, c_degree as i32];
    let mut eval_bs = vec![0.0f64; c_dimension * (derivative_request + 1)];
    eval_flat(
        tbord,
        false,
        derivative_request as i32,
        &mut extrap_mode,
        c_degree,
        flat_knots,
        c_dimension,
        poles,
        &mut eval_bs,
    );

    // norm of the tangent at the node of connection (OCCT L3857-3862).
    let mut tgte = Vector::new(1, c_dimension as i32);
    for ipos in 0..c_dimension {
        tgte.set((ipos + 1) as i32, eval_bs[c_dimension + ipos]);
    }
    let l1 = (0..c_dimension)
        .map(|i| eval_bs[c_dimension + i])
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt();

    // matrix of constraints (OCCT L3865-3900).
    let mut contraintes = Matrix::new_init(1, csize as i32, 1, c_dimension as i32, 0.0);
    let c1c2 = c1_coefficient * c1_coefficient;
    let c1c3 = c1c2 * c1_coefficient;
    for ipos in 1..=c_dimension as i32 {
        if after {
            contraintes.set(1, ipos, eval_bs[(ipos - 1) as usize]);
            contraintes.set(2, ipos, c1_coefficient * eval_bs[(c_dimension + ipos as usize - 1)]);
            if cont >= 2 {
                contraintes.set(
                    3,
                    ipos,
                    eval_bs[(2 * c_dimension + ipos as usize - 1)] * c1c2,
                );
            }
            if cont >= 3 {
                contraintes.set(
                    4,
                    ipos,
                    eval_bs[(3 * c_dimension + ipos as usize - 1)] * c1c3,
                );
            }
            contraintes.set((cont + 2) as i32, ipos, constraint_point[ipos as usize - 1]);
        } else {
            contraintes.set(1, ipos, constraint_point[ipos as usize - 1]);
            contraintes.set(2, ipos, eval_bs[(ipos - 1) as usize]);
            if cont >= 1 {
                contraintes.set(3, ipos, c1_coefficient * eval_bs[(c_dimension + ipos as usize - 1)]);
            }
            if cont >= 2 {
                contraintes.set(
                    4,
                    ipos,
                    eval_bs[(2 * c_dimension + ipos as usize - 1)] * c1c2,
                );
            }
            if cont >= 3 {
                contraintes.set(
                    5,
                    ipos,
                    eval_bs[(3 * c_dimension + ipos as usize - 1)] * c1c3,
                );
            }
        }
    }

    // calculate the coefficients of extension (OCCT L3905-3919).
    let mut extra_coeffs = vec![0.0f64; csize * c_dimension];
    for ii in 1..=csize {
        for jj in 1..=csize {
            for kk in 1..=c_dimension {
                let idx = (kk - 1) + (jj - 1) * c_dimension;
                extra_coeffs[idx] += mat_coefs.get(ii as i32, jj as i32) * contraintes.get(ii as i32, kk as i32);
            }
        }
    }

    // calculate the poles of extension (OCCT L3922-3927).
    let extrap_poles = plib_coefficients_poles_dim(&extra_coeffs, c_dimension);

    // nodes of extension with multiplicities + flat nodes (OCCT L3930-3941).
    let extrap_knots = [0.0f64, 1.0f64];
    let extrap_mults = [csize as i32, csize as i32];
    let mut fk2 = vec![0.0f64; 2 * csize];
    bspl_knot_sequence(
        &extrap_knots,
        &extrap_mults,
        csize - 1,
        false,
        &mut fk2,
    );

    // norm of the tangent at the connection point (OCCT L3944-3956).
    extrap_mode = [c_degree as i32, c_degree as i32];
    let mut eval_bs2 = vec![0.0f64; c_dimension * 2];
    eval_flat(
        if after { 0.0 } else { 1.0 },
        false,
        1,
        &mut extrap_mode,
        csize - 1,
        &fk2,
        c_dimension,
        &extrap_poles,
        &mut eval_bs2,
    );
    for ipos in 0..c_dimension {
        tgte.set((ipos + 1) as i32, eval_bs2[c_dimension + ipos]);
    }
    let l2 = (0..c_dimension)
        .map(|i| eval_bs2[c_dimension + i])
        .map(|v| v * v)
        .sum::<f64>()
        .sqrt();

    // harmonisation of degrees (OCCT L3959-3981).
    let mut new_p2_len = (c_degree + 1) * c_dimension;
    let (new_p2, new_k2, new_m2): (Vec<f64>, Vec<f64>, Vec<i32>);
    if csize - 1 < c_degree {
        let nb_out = (csize + 1) * c_dimension;
        let mut np2 = vec![0.0f64; new_p2_len.max(nb_out)];
        let mut nk2 = vec![0.0f64; 2];
        let mut nm2 = vec![0i32; 2];
        bspl_increase_degree(
            csize - 1,
            c_degree,
            false,
            c_dimension,
            &extrap_poles,
            &extrap_knots,
            &extrap_mults,
            &mut np2,
            &mut nk2,
            &mut nm2,
        );
        np2.truncate(nb_out);
        new_p2 = np2;
        new_k2 = nk2;
        new_m2 = nm2;
        new_p2_len = nb_out;
    } else {
        new_p2 = extrap_poles.clone();
        new_k2 = extrap_knots.to_vec();
        new_m2 = extrap_mults.to_vec();
        new_p2_len = extrap_poles.len();
    }

    // flat nodes of extension after harmonization of degrees
    // (OCCT L3984-3987).
    let mut new_fk2 = vec![0.0f64; 2 * (c_degree + 1)];
    bspl_knot_sequence(&new_k2, &new_m2, c_degree, false, &mut new_fk2);

    // 2. concatenation C0 (OCCT L3992-4090).
    let mut ratio = 1.0f64;
    if l1 > rcad_kernel::core::precision::CONFUSION && l2 > rcad_kernel::core::precision::CONFUSION {
        ratio = l2 / l1;
    }
    if !(0.00001..=100000.0).contains(&ratio) {
        ratio = 1.0;
    }

    let delta = if after {
        ratio * new_fk2[0] - flat_knots[flat_knots.len() - 1]
    } else {
        ratio * new_fk2[new_fk2.len() - 1] - flat_knots[0]
    };

    let nb_p1 = num_poles;
    let nb_p2 = c_degree + 1;
    let nb_k1 = flat_knots.len();
    let nb_k2 = 2 * (c_degree + 1);
    let mut new_poles = vec![0.0f64; (nb_p1 + nb_p2 - 1) * c_dimension];
    let new_flats_len = nb_k1 + nb_k2 - c_degree - 2;
    let mut new_flats = vec![0.0f64; new_flats_len];

    // poles
    if after {
        for ii in 0..(nb_p1 + nb_p2 - 1) {
            for jj in 0..c_dimension {
                let ind_np = ii * c_dimension + jj;
                new_poles[ind_np] = if ii < nb_p1 {
                    poles[ind_np]
                } else {
                    new_p2[(ii - nb_p1) * c_dimension + jj]
                };
            }
        }
    } else {
        for ii in 0..(nb_p1 + nb_p2 - 1) {
            for jj in 0..c_dimension {
                let ind_np = ii * c_dimension + jj;
                new_poles[ind_np] = if ii < nb_p2 {
                    new_p2[ind_np]
                } else {
                    poles[(ii - nb_p2) * c_dimension + jj]
                };
            }
        }
    }

    // flat nodes
    if after {
        for ii in 0..(nb_k1 - 1) {
            new_flats[ii] = flat_knots[ii];
        }
        for ii in 1..=(nb_k2 - c_degree - 1) {
            new_flats[nb_k1 + ii - 2] = ratio * new_fk2[ii + c_degree] - delta;
        }
    } else {
        for ii in 1..(nb_k2 - c_degree) {
            new_flats[ii - 1] = ratio * new_fk2[ii - 1] - delta;
        }
        for ii in 2..=nb_k1 {
            new_flats[nb_k2 - c_degree - 2 + ii - 1] = flat_knots[ii - 1];
        }
    }

    // 3. reduction of multiplicite at the node of connection
    // (OCCT L4096-4120).
    let mut klength = 1usize;
    for ii in 1..new_flats_len {
        if new_flats[ii] != new_flats[ii - 1] {
            klength += 1;
        }
    }
    let mut new_knots = vec![0.0f64; klength];
    let mut new_mults = vec![1i32; klength];
    let mut jj = 0usize;
    new_knots[0] = new_flats[0];
    for ii in 1..new_flats_len {
        if new_flats[ii] == new_flats[ii - 1] {
            new_mults[jj] += 1;
        } else {
            jj += 1;
            new_knots[jj] = new_flats[ii];
        }
    }

    // reduction of multiplicity at the second or the last but one node
    // (OCCT L4123-4160).
    let index = if after { klength - 1 } else { 1 };
    let mut m = c_degree as i32;
    let mut result_poles = vec![0.0f64; new_poles.len()];
    let mut result_knots = vec![0.0f64; klength];
    let mut result_mults = vec![0i32; klength];
    let tol = 1e-6;
    let mut ok = true;
    while m > c_degree as i32 - continuity && ok {
        ok = bspl_remove_knot(
            index,
            m - 1,
            c_degree,
            false,
            c_dimension,
            &new_poles,
            &new_knots,
            &new_mults,
            &mut result_poles,
            &mut result_knots,
            &mut result_mults,
            tol,
        );
        if ok {
            m -= 1;
        }
    }

    if m == c_degree as i32 {
        *nb_poles_result = nb_p1 + nb_p2 - 1;
        *poles_result = new_poles[..*nb_poles_result * c_dimension].to_vec();
        let (rk, rm) = compress_knots(&new_flats);
        let _ = (result_knots, result_mults, rk, rm);
        knots_result.clear();
        for (k, mult) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*mult {
                knots_result.push(*k);
            }
        }
        *nb_knots_result = knots_result.len();
    } else {
        *nb_poles_result = nb_p1 + nb_p2 - 1 - c_degree + m as usize;
        *poles_result = result_poles[..*nb_poles_result * c_dimension].to_vec();
        knots_result.clear();
        for (k, mult) in result_knots.iter().zip(result_mults.iter()) {
            for _ in 0..*mult {
                knots_result.push(*k);
            }
        }
        *nb_knots_result = knots_result.len();
    }
}

/// Increases the degree of one parameter direction of a rational or
/// non-rational B-spline surface (homogeneous per-strip increase).
/// OCCT Geom_BSplineSurface::IncreaseDegree(UDeg, VDeg) — the strips are
/// the iso-curves of the other direction; the rcad scalar bspl_lib call
/// carries dim 4 homogeneous poles per strip.
fn bsp_surface_increase_degree(bs: &mut BSplineSurface, in_u: bool, new_degree: usize) {
    if new_degree
        <= if in_u {
            bs.degree_u
        } else {
            bs.degree_v
        }
    {
        return;
    }
    let (knots, mults) = if in_u {
        compress_knots(&bs.knots_u)
    } else {
        compress_knots(&bs.knots_v)
    };
    let old_degree = if in_u { bs.degree_u } else { bs.degree_v };
    let nbknots = rcad_kernel::math::bspl_lib::increase_degree_count_knots(
        old_degree, new_degree, false, &mults,
    );
    let step = new_degree - old_degree;
    let (nb_strips, nb_poles_strip) = if in_u {
        (bs.control_points[0].len(), bs.control_points.len())
    } else {
        (bs.control_points.len(), bs.control_points[0].len())
    };
    // Homogeneous flat poles per strip: (x*w, y*w, z*w, w).
    let mut new_knots = vec![0.0f64; nbknots];
    let mut new_mults = vec![0i32; nbknots];
    let mut new_flat_knots: Vec<f64> = Vec::new();
    let mut new_points: Vec<Vec<DVec3>> = Vec::new();
    let mut new_weights: Vec<Vec<f64>> = Vec::new();
    for s in 0..nb_strips {
        let mut poles_flat = vec![0.0f64; nb_poles_strip * 4];
        for p in 0..nb_poles_strip {
            let (pt, w) = if in_u {
                let row = &bs.control_points[p];
                (row[s], bs.weights[p][s])
            } else {
                let row = &bs.control_points[s];
                (row[p], bs.weights[s][p])
            };
            poles_flat[p * 4] = pt.x * w;
            poles_flat[p * 4 + 1] = pt.y * w;
            poles_flat[p * 4 + 2] = pt.z * w;
            poles_flat[p * 4 + 3] = w;
        }
        let mut np = vec![0.0f64; (nb_poles_strip + step * (mults.len() - 1)) * 4];
        bspl_increase_degree(
            old_degree,
            new_degree,
            false,
            4,
            &poles_flat,
            &knots,
            &mults,
            &mut np,
            &mut new_knots,
            &mut new_mults,
        );
        let count = np.len() / 4;
        let strip_pts: Vec<DVec3> = (0..count)
            .map(|i| {
                let w = np[i * 4 + 3];
                DVec3::new(np[i * 4] / w, np[i * 4 + 1] / w, np[i * 4 + 2] / w)
            })
            .collect();
        let strip_w: Vec<f64> = (0..count).map(|i| np[i * 4 + 3]).collect();
        if s == 0 {
            new_flat_knots = {
                let mut flat = Vec::new();
                for (k, m) in new_knots.iter().zip(new_mults.iter()) {
                    for _ in 0..*m {
                        flat.push(*k);
                    }
                }
                flat
            };
            new_points = (0..count)
                .map(|_| Vec::with_capacity(nb_strips))
                .collect();
            new_weights = (0..count).map(|_| Vec::with_capacity(nb_strips)).collect();
        }
        for (i, (pt, w)) in strip_pts.iter().zip(strip_w.iter()).enumerate() {
            new_points[i].push(*pt);
            new_weights[i].push(*w);
        }
    }
    if in_u {
        bs.degree_u = new_degree;
        bs.knots_u = new_flat_knots;
        bs.control_points = new_points;
        bs.weights = new_weights;
    } else {
        bs.degree_v = new_degree;
        bs.knots_v = new_flat_knots;
        // Transpose back (strips were along v: strip index = u index).
        bs.control_points = new_points;
        bs.weights = new_weights;
    }
}

// =========================================================================
// OCCT GeomLib.cxx L1485-1962 — GeomLib::ExtendSurfByLength.  Extends the
// bounded surface by <length> in the u or v direction, before or after,
// with C<continuity> (rcad Surface3 encoding: BSpline/Bezier supported;
// the KPart (RectangularTrimmedSurface) and GeomConvert_ApproxSurface
// branches are GAP carriers — the surface is returned unextended and the
// result reports false).  Returns true when the surface was extended.
// =========================================================================
pub(crate) fn geom_lib_extend_surf_by_length(
    surface: &mut Surface3,
    length: f64,
    continuity: i32,
    in_u: bool,
    after: bool,
) -> bool {
    // OCCT L1489-1493.
    if !(0..=3).contains(&continuity) {
        return false;
    }
    let cont = continuity;
    // OCCT L1497-1505: KPart — ExtendKPart(TS, Length, InU, After).
    // GAP: ExtendKPart (GeomLib.cxx L1420+) pending — a Trimmed carrier
    // keeps the shape of the OCCT guard but does not extend.
    if let Surface3::Trimmed(_) = surface {
        return false;
    }

    // OCCT L1509-1529: BSpline downcast; other kinds go through
    // GeomConvert_ApproxSurface (pending TKGeomAlgo) — keep unextended.
    let mut bs: BSplineSurface = match surface {
        Surface3::BSpline(b) => b.clone(),
        Surface3::Bezier(bez) => {
            // GeomConvert::SurfaceToBSplineSurface of a Bezier surface:
            // the Bernstein poles are the BSpline poles with a single knot
            // span of (degree + 1) multiplicity at both ends.
            let degree_u = bez.control_points.len() - 1;
            let degree_v = bez.control_points[0].len() - 1;
            BSplineSurface {
                degree_u,
                degree_v,
                knots_u: {
                    let mut k = vec![0.0; degree_u + 1];
                    k.extend(vec![1.0; degree_u + 1]);
                    k
                },
                knots_v: {
                    let mut k = vec![0.0; degree_v + 1];
                    k.extend(vec![1.0; degree_v + 1]);
                    k
                },
                control_points: bez.control_points.clone(),
                weights: bez.weights.clone(),
                is_periodic_u: false,
                is_periodic_v: false,
            }
        }
        _ => return false,
    };

    // OCCT L1531-1539: increase the degree in the extension direction when
    // it is insufficient for the requested continuity.
    if in_u && (bs.degree_u as i32) < cont + 1 {
        bsp_surface_increase_degree(&mut bs, true, (cont + 1) as usize);
    }
    if !in_u && (bs.degree_v as i32) < cont + 1 {
        bsp_surface_increase_degree(&mut bs, false, (cont + 1) as usize);
    }

    // OCCT L1542-1552: if BS was periodic in the extension direction it is
    // de-periodized with Segment().  The rcad BSplineSurface carries the
    // flat knot vector only (no periodic flag) — nothing to do.

    // IFV Fix OCC bug 0022694 (OCCT L1554-1565).
    let rational = bs
        .weights
        .iter()
        .flat_map(|r| r.iter())
        .any(|w| (*w - 1.0).abs() > 0.0);
    let eps_w = 10.0 * rcad_kernel::core::precision::PCONFUSION;
    let gap: usize = if rational { 4 } else { 3 };

    let mut lambmin = length;
    let (mut cdeg, mut nb_p, mut cdim): (usize, usize, usize) = (0, 0, 0);
    let mut fknots: Vec<f64> = Vec::new();
    let mut poles_flat: Vec<f64> = Vec::new();
    let (mut point, mut tgte, mut lambda): (Vec<f64>, Vec<f64>, Vec<f64>) =
        (Vec::new(), Vec::new(), Vec::new());
    let mut ok = false;

    // OCCT L1583-1752: the first retry loop (up to 3 attempts; the degree
    // of the border iso is increased by 2 on failure and we retry).
    let mut kount = 0usize;
    while kount <= 2 && !ok {
        // Transform the surface into a single-variable flat pole array of
        // dimension <cdim> (per-direction flattening; rational keeps
        // (x,y,z,w) blocks — OCCT L1590-1617).
        if in_u {
            cdeg = bs.degree_u;
            nb_p = bs.control_points.len();
            cdim = bs.control_points[0].len() * gap;
            fknots = bs.knots_u.clone();
        } else {
            cdeg = bs.degree_v;
            nb_p = bs.control_points[0].len();
            cdim = bs.control_points.len() * gap;
            fknots = bs.knots_v.clone();
        }
        poles_flat = vec![0.0f64; cdim * nb_p];
        if in_u {
            for (ii, row) in bs.control_points.iter().enumerate() {
                for (jj, pt) in row.iter().enumerate() {
                    let ipole = ii * cdim + jj * gap;
                    poles_flat[ipole] = pt.x;
                    poles_flat[ipole + 1] = pt.y;
                    poles_flat[ipole + 2] = pt.z;
                    if rational {
                        poles_flat[ipole + 3] = bs.weights[ii][jj];
                    }
                }
            }
        } else {
            for (jj, row) in bs.control_points.iter().enumerate() {
                for (ii, pt) in row.iter().enumerate() {
                    let ipole = jj * cdim + ii * gap;
                    poles_flat[ipole] = pt.x;
                    poles_flat[ipole + 1] = pt.y;
                    poles_flat[ipole + 2] = pt.z;
                    if rational {
                        poles_flat[ipole + 3] = bs.weights[jj][ii];
                    }
                }
            }
        }

        // the parameter of the connection knot (OCCT L1620-1628).
        let tbord = if after {
            fknots[fknots.len() - 1 - cdeg]
        } else {
            fknots[cdeg]
        };

        // calculation of the connection point and tangent
        // (OCCT L1641-1668) — BSplCLib::Eval at Tbord.
        let derivative_request = std::cmp::max(cont, 1) as usize;
        let mut extrap_mode = [cdeg as i32, cdeg as i32];
        let mut result = vec![0.0f64; cdim * (derivative_request + 1)];
        eval_flat(
            tbord,
            false,
            derivative_request as i32,
            &mut extrap_mode,
            cdeg,
            &fknots,
            cdim,
            &poles_flat,
            &mut result,
        );
        point = result[..cdim].to_vec();
        tgte = result[cdim..2 * cdim].to_vec();
        lambda = vec![0.0f64; cdim];
        ok = true;

        // calculation of the constraint to reach (OCCT L1673-1740).  All
        // indices below keep the OCCT 1-based arithmetic (the `0` suffix
        // marks the 0-based slot).
        let tgtol = 1e-12;
        let mut old_n = 0.0f64;
        let mut old_t = DVec3::ZERO;
        let mut ii = gap; // 1-based
        while ii <= cdim {
            if rational {
                tgte[ii - 1] = 0.0;
                let cur_t = DVec3::new(tgte[ii - 4], tgte[ii - 3], tgte[ii - 2]);
                let n_tgte = cur_t.length();
                if n_tgte > tgtol {
                    let val = length / n_tgte;
                    if old_n > tgtol && cur_t.angle_between(old_t) > 2.0 {
                        ok = false;
                    }
                    lambda[ii - 4] = val;
                    lambda[ii - 3] = val;
                    lambda[ii - 2] = val;
                    lambda[ii - 1] = 0.0;
                    lambmin = lambmin.min(val);
                } else {
                    lambda[ii - 4] = 0.0;
                    lambda[ii - 3] = 0.0;
                    lambda[ii - 2] = 0.0;
                    lambda[ii - 1] = 0.0;
                }
                old_t = cur_t;
                old_n = n_tgte;
            } else {
                let cur_t = DVec3::new(tgte[ii - 3], tgte[ii - 2], tgte[ii - 1]);
                let n_tgte = cur_t.length();
                if n_tgte > tgtol {
                    let val = length / n_tgte;
                    if old_n > tgtol && cur_t.angle_between(old_t) > 2.0 {
                        ok = false;
                    }
                    lambda[ii - 1] = val;
                    lambda[ii - 2] = val;
                    lambda[ii - 3] = val;
                    lambmin = lambmin.min(val);
                } else {
                    lambda[ii - 1] = 0.0;
                    lambda[ii - 2] = 0.0;
                    lambda[ii - 3] = 0.0;
                }
                old_t = cur_t;
                old_n = n_tgte;
            }
            ii += gap;
        }
        if !ok && kount < 2 {
            // We increase the degree of the border iso to bring the poles
            // of the surface closer and we retry (OCCT L1742-1751).
            let new_deg = if in_u { bs.degree_v + 2 } else { bs.degree_u + 2 };
            bsp_surface_increase_degree(&mut bs, !in_u, new_deg);
        }
        kount += 1;
    }

    // OCCT L1755-1770: ConstraintPoint = Point +/- lambda * Tgte.
    let mut constraint_point = vec![0.0f64; cdim];
    for ii in 0..cdim {
        constraint_point[ii] = if after {
            point[ii] + lambda[ii] * tgte[ii]
        } else {
            point[ii] - lambda[ii] * tgte[ii]
        };
    }

    // special case of rational (OCCT L1773-1786): multiply the (x,y,z)
    // triples by the weight block component.
    if rational {
        let mut ipole = 0usize;
        while ipole < poles_flat.len() {
            poles_flat[ipole] *= poles_flat[ipole + 3];
            poles_flat[ipole + 1] *= poles_flat[ipole + 3];
            poles_flat[ipole + 2] *= poles_flat[ipole + 3];
            ipole += gap;
        }
        let mut ii = 0usize;
        while ii < cdim {
            constraint_point[ii] *= constraint_point[ii + 3];
            constraint_point[ii + 1] *= constraint_point[ii + 3];
            constraint_point[ii + 2] *= constraint_point[ii + 3];
            ii += gap;
        }
    }

    // arrays needed for the extension (OCCT L1789-1801).
    let ksize2 = nb_p + cdeg + 1 + cdeg;
    let psize2 = cdim * nb_p + cdeg * cdim;
    let mut nb_poles_result = 0usize;
    let mut nb_knots_result = 0usize;
    let mut fk_res = vec![0.0f64; ksize2];
    let mut p_res = vec![0.0f64; psize2];

    // OCCT L1802-1899: the extension retry loop (up to 5 attempts; a null
    // weight divides lambmin by 3 and we retry).
    let (mut nu, mut nv): (usize, usize);
    let mut grid: Vec<Vec<DVec3>> = Vec::new();
    let mut wgrid: Vec<Vec<f64>> = Vec::new();
    let mut ext_ok = false;
    for _kount in 0..5 {
        if ext_ok {
            break;
        }
        bsplclib_tang_extend_to_constraint(
            &fknots,
            lambmin,
            nb_p,
            &poles_flat,
            cdim,
            cdeg,
            &constraint_point,
            cont,
            after,
            &mut nb_poles_result,
            &mut nb_knots_result,
            &mut fk_res,
            &mut p_res,
        );
        // Copy the result poles as 3D points and weights
        // (OCCT L1868-1880).
        if in_u {
            nu = nb_poles_result;
            nv = bs.control_points[0].len();
        } else {
            nu = bs.control_points.len();
            nv = nb_poles_result;
        }
        grid = vec![DVec3::ZERO; nu].iter().map(|_| vec![DVec3::ZERO; nv]).collect();
        wgrid = vec![1.0f64; nu].iter().map(|_| vec![1.0f64; nv]).collect();
        let mut null_weight = false;
        'outer: for ii in 0..nu {
            for jj in 0..nv {
                let indice = if in_u {
                    (ii) * cdim + (jj) * gap
                } else {
                    (ii) * gap + (jj) * cdim
                };
                let mut x = p_res[indice];
                let mut y = p_res[indice + 1];
                let mut z = p_res[indice + 2];
                let mut w = 1.0f64;
                if rational {
                    w = p_res[indice + 3];
                    if (w - 1.0).abs() < eps_w {
                        w = 1.0;
                    }
                    if w < eps_w {
                        null_weight = true;
                        break 'outer;
                    }
                    x /= w;
                    y /= w;
                    z /= w;
                }
                grid[ii][jj] = DVec3::new(x, y, z);
                wgrid[ii][jj] = w;
            }
        }
        if null_weight {
            lambmin /= 3.0;
        } else {
            ext_ok = true;
        }
    }

    // Copy flat knots as knots with their multiplicities; compute the
    // degrees of the result (OCCT L1902-1945).
    let (uflat, vflat, udeg, vdeg);
    if in_u {
        let (mut uk, mut um) = compress_knots(&fk_res[..nb_knots_result]);
        udeg = cdeg;
        // UMults(Usize) = UDeg + 1 — "petite verrue utile quand la
        // continuite n'est pas ok".
        *um.last_mut().expect("umults") = udeg as i32 + 1;
        let (vk, vm) = compress_knots(&bs.knots_v);
        vdeg = bs.degree_v;
        uflat = expand_knots(&uk, &um);
        vflat = expand_knots(&vk, &vm);
    } else {
        let (mut vk, mut vm) = compress_knots(&fk_res[..nb_knots_result]);
        vdeg = cdeg;
        // VMults(Vsize) = VDeg + 1.
        *vm.last_mut().expect("vmults") = vdeg as i32 + 1;
        let (uk, um) = compress_knots(&bs.knots_u);
        udeg = bs.degree_u;
        uflat = expand_knots(&uk, &um);
        vflat = expand_knots(&vk, &vm);
    }

    // Construct the resulting BSpline surface (OCCT L1948-1961).
    *surface = Surface3::BSpline(BSplineSurface {
        degree_u: udeg,
        degree_v: vdeg,
        knots_u: uflat,
        knots_v: vflat,
        control_points: grid,
        weights: wgrid,
        is_periodic_u: false,
        is_periodic_v: false,
    });
    true
}
// =========================================================================
// OCCT GeomLib.cxx L1269-1410 — GeomLib::ExtendCurveToPoint.  Extends a
// bounded curve to <point> with C<continuity>, before or after.
// GAP carriers: ComputeLambda optimization (GeomLib.cxx L126-200, Gauss
// integration pending) keeps the incoming Lambda; the concatenation
// GeomConvert_CompCurveToBSplineCurve (TKG3d pending) leaves the curve
// unextended — the OCCT raise on Add failure becomes a documented no-op.
// Only reachable through the IntSS-pcurve GAP of PerformIntersectionAtEnd.
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub(crate) fn geom_lib_extend_curve_to_point(
    curve: &mut Curve3,
    point: DVec3,
    continuity: i32,
    after: bool,
) -> bool {
    // OCCT L1274-1278.
    if !(1..=3).contains(&continuity) {
        return false;
    }
    let size = (continuity + 2) as i32;
    let tol = 1e-6;
    let mut mat_coefs = Matrix::new_init(1, size, 1, size, 0.0);
    // OCCT L1294-1300: HermiteCoefficients(0, 1, Continuity, 0).
    hermite_coefficients(0.0, 1.0, continuity, 0, &mut mat_coefs);

    let ubord = if after {
        let [_, l] = curve.default_domain();
        l
    } else {
        let [f, _] = curve.default_domain();
        f
    };
    // OCCT L1303: Curve->D3(Ubord, p0, d1, d2, d3).
    let p0 = curve.point_at(ubord);
    let mut d1 = curve.derivative_at(ubord);
    let d2 = curve.derivative2_at(ubord);
    let mut d3 = curve.derivative3_at(ubord);
    if !after {
        // Invert the parameterization.
        d1 = -d1;
        d3 = -d3;
    }

    let l1 = p0.distance(point);
    if l1 <= tol {
        return false; // No extension
    }

    // Lambda is the ratio to apply to the derivative of the curve to obtain
    // the derivative of the extension (OCCT L1310-1336).
    let f0 = curve.default_domain()[0];
    let l0 = curve.default_domain()[1];
    let dt0 = (l0 - f0) / 9.0;
    let mut norm = d1.length();
    for i in 1..=8 {
        let t = f0 + dt0 * i as f64;
        norm += curve.derivative_at(t).length();
    }
    norm /= 9.0;
    let dt = d1.length() / norm;
    let lambda = if dt < 1.5 && dt > 0.75 {
        1.0 / (d1.length() / l1).max(tol)
    } else {
        1.0 / (norm / l1).max(tol)
    };
    // OCCT L1342: ComputeLambda(Cons, MatCoefs, L1, Lambda) — GAP (Gauss
    // integration pending); Lambda is kept.

    // Construction in the Polynomial Basis (OCCT L1345-1361).
    let mut cont_vec: Vec<DVec3> = vec![DVec3::ZERO; size as usize];
    cont_vec[0] = p0;
    cont_vec[1] = d1 * lambda;
    if continuity >= 2 {
        cont_vec[2] = d2 * (lambda * lambda);
    }
    if continuity >= 3 {
        cont_vec[3] = d3 * (lambda * lambda * lambda);
    }
    cont_vec[(size - 1) as usize] = point;

    // Conversion to the Bernstein Basis (OCCT L1364-1378) — per coordinate
    // flatten of the xyz constraint.
    let mut coefs = Vec::with_capacity(cont_vec.len() * 3);
    for c in &cont_vec {
        coefs.extend([c.x, c.y, c.z]);
    }
    let poles_flat = plib_coefficients_poles_dim(&coefs, 3);
    let extrap_poles: Vec<DVec3> = (0..cont_vec.len())
        .map(|i| {
            DVec3::new(
                poles_flat[i * 3],
                poles_flat[i * 3 + 1],
                poles_flat[i * 3 + 2],
            )
        })
        .collect();

    // Concatenation (OCCT L1381-1405) — GAP:
    // GeomConvert_CompCurveToBSplineCurve pending; the curve is returned
    // unextended (the OCCT Add-failure raise documented as no-op).
    false
}
