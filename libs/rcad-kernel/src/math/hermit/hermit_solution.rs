//! OCCT Hermit package (ModelingData/TKGeomBase/Hermit/Hermit.cxx L1-1207).
//!
//! Used to reparameterize rational BSpline curves so that they can be
//! concatenated into C1 curves: it builds the 1D reparameterizing function
//! `a(u)` from a Hermite interpolation, adding knots and modifying poles so
//! that `a(u)*D(u)` has value 1 at umin/umax and zero derivative there.
//!
//! Translated here: the whole Geom2d arm consumed by
//! `Geom2dConvert::ConcatC1` — `HermiteCoeff(Geom2d)` (L91-139),
//! `SignDenom` (L146-153), `Polemax` (L160-185), `PolyTest(Geom2d)`
//! (L520-842), `InsertKnots` (L849-860), `MovePoles` (L867-880) and
//! `Hermit::Solution(Geom2d)` (L994-1098).
//!
//! NOT translated (no carrier yet): the `Geom_BSplineCurve` (3D) overloads —
//! `HermiteCoeff(Geom)` (L36-83), `PolyTest(Geom)` (L192-513),
//! `Solution(Geom)` (L884-988) and `Solutionbis` (L1102-1207).  Their bodies
//! are identical to the Geom2d forms modulo the curve type; rcad-kernel has
//! no OCCT-member-data `Geom_BSplineCurve` struct (only the legacy flat
//! `BSplineCurve3`), which is a prerequisite for hosting them.  The 3D
//! SameParameter arm will need this follow-up.

use super::bspl_eval_kernels::bspl_clib_d1_scalar;
use crate::base::geom2d_convert::bspline_curve::Geom2dBSplineCurve;
use crate::base::geom2d_convert::c1_concat::{curve_knots, curve_mults, curve_weights};
use crate::core::precision::CONFUSION;
use crate::math::bspl_lib::{at, locate_parameter_main, reparametrize};
use glam::DVec2;

/// OCCT static HermiteCoeff (Hermit.cxx L91-139) — the Geom2d form: computes
/// the four Hermite coefficients of degree 3 from `bs` (the weight function
/// evaluated as a scalar BSpline at u=0 and u=1) and stores them in `herm`.
fn hermite_coeff(bs: &Geom2dBSplineCurve, herm: &mut [f64; 4]) {
    // NCollection_Array1<double> Knots(BS->Knots());
    let mut knots = curve_knots(bs);
    let weights = curve_weights(bs);
    let mults = curve_mults(bs);

    let degree = bs.degree();
    let periodic = bs.is_periodic();
    let index0 = bs.first_uknot_index();
    let index1 = bs.last_uknot_index() - 1;

    // BSplCLib::Reparametrize(0.0, 1.0, Knots); // affinity on the nodal vector
    reparametrize(0.0, 1.0, &mut knots);

    let mut denom0 = 0.0f64;
    let mut deriv0 = 0.0f64;
    let mut denom1 = 0.0f64;
    let mut deriv1 = 0.0f64;
    // Evaluate the weight function w(u) and its derivative at u=0 and u=1:
    // Weights passed as Poles, BSplCLib::NoWeights() for the unweighted arm.
    bspl_clib_d1_scalar(
        0.0,
        index0,
        degree,
        periodic,
        &weights,
        None,
        &knots,
        &mults,
        &mut denom0,
        &mut deriv0,
    );
    bspl_clib_d1_scalar(
        1.0,
        index1,
        degree,
        periodic,
        &weights,
        None,
        &knots,
        &mults,
        &mut denom1,
        &mut deriv1,
    );
    // Hermit coefficients.
    herm[0] = 1.0 / denom0;
    herm[1] = -deriv0 / (denom0 * denom0);
    herm[2] = -deriv1 / (denom1 * denom1);
    herm[3] = 1.0 / denom1;
}

/// OCCT static SignDenom (Hermit.cxx L146-153) — the sign of Herm(0),
/// true = positive.
fn sign_denom(poles: &[DVec2; 4]) -> bool {
    poles[0].y >= 0.0
}

/// OCCT static Polemax (Hermit.cxx L160-185) — the indices of the poles with
/// minimum and maximum ordinates.
fn polemax(poles: &[DVec2; 4], min: &mut i32, max: &mut i32) {
    *min = 0;
    *max = 0;
    let mut min_val = poles[0].y;
    let mut max_val = poles[0].y;
    for i in 1..=(poles.len() as i32 - 1) {
        if poles[i as usize].y < min_val {
            min_val = poles[i as usize].y;
            *min = i;
        }
        if poles[i as usize].y > max_val {
            max_val = poles[i as usize].y;
            *max = i;
        }
    }
}

/// OCCT BSplCLib::LocateParameter(Degree, Knots, U, IsPeriodic, FromK1, ToK2,
/// KnotIndex, NewU) — the flat (no Mults) overload (BSplCLib.cxx L189-214).
fn locate_parameter_flat(
    _degree: i32,
    knots: &[f64],
    u: f64,
    is_periodic: bool,
    from_k1: i32,
    to_k2: i32,
    knot_index: &mut i32,
    new_u: &mut f64,
) {
    let (uf, ul) = if is_periodic {
        let k_upper = knots.len() as i32;
        (
            at(knots, 1 + _degree),
            at(knots, k_upper - _degree),
        )
    } else {
        (0.0f64, 1.0f64)
    };
    locate_parameter_main(knots, u, is_periodic, from_k1, to_k2, knot_index, new_u, uf, ul);
}

/// OCCT static PolyTest (Hermit.cxx L520-842) — the Geom2d form: computes the
/// knots U4 and U5 to insert to `a(u)`.
fn poly_test(
    herm: &[f64; 4],
    bs: &Geom2dBSplineCurve,
    u4: &mut f64,
    u5: &mut f64,
    boucle: &mut i32,
    tol_poles: f64,
    _tol_knots: f64,
    ux: f64,
    uy: f64,
) {
    let mut i1 = 0i32;
    let mut i2 = 0i32;
    let mut i3 = 0i32;
    let mut i4 = 0i32;
    let mut polesinit = [DVec2::ZERO; 4];
    let mut cas = 0i32;
    let mut mark = 0i32;
    let mut dercas = 0i32;
    let mut min = 0i32;
    let mut max = 0i32;

    *u4 = 0.0;
    *u5 = 1.0; // default value
    if ux != 1.0 {
        bs.locate_u(ux, 0.0, &mut i1, &mut i2, false); // localization of the inserted knots
        if uy != 0.0 {
            bs.locate_u(uy, 0.0, &mut i3, &mut i4, false);
        }
    }

    // definition and filling of the array of knots (L552-593).
    let mut knots: Vec<f64>;
    if i1 == i2 {
        if (i3 == i4) || (i3 == 0) {
            knots = vec![0.0; bs.nb_knots() as usize];
            for i in 1..=bs.nb_knots() {
                knots[(i - 1) as usize] = bs.knot(i);
            }
        } else {
            knots = vec![0.0; bs.nb_knots() as usize + 1];
            for i in 1..=bs.nb_knots() {
                knots[(i - 1) as usize] = bs.knot(i);
            }
            knots[bs.nb_knots() as usize] = uy;
        }
    } else {
        if (i3 == i4) || (i3 == 0) {
            knots = vec![0.0; bs.nb_knots() as usize + 1];
            for i in 1..=bs.nb_knots() {
                knots[(i - 1) as usize] = bs.knot(i);
            }
            knots[bs.nb_knots() as usize] = ux;
        } else {
            knots = vec![0.0; bs.nb_knots() as usize + 2];
            for i in 1..=bs.nb_knots() {
                knots[(i - 1) as usize] = bs.knot(i);
            }
            knots[bs.nb_knots() as usize] = ux;
            knots[bs.nb_knots() as usize + 1] = uy;
        }
    }

    // sort of the array of knots.
    knots.sort_by(|a, b| a.partial_cmp(b).expect("NaN knot in PolyTest"));

    // poles of the Hermite polynome in the BSpline form.
    polesinit[0] = DVec2::new(0.0, herm[0]);
    polesinit[1] = DVec2::new(0.0, herm[0] + herm[1] / 3.0);
    polesinit[2] = DVec2::new(0.0, herm[3] - herm[2] / 3.0);
    polesinit[3] = DVec2::new(0.0, herm[3]);

    // loop to check the tolerances on poles.
    if tol_poles != 0.0 {
        polemax(&polesinit, &mut min, &mut max);
        let polemin = polesinit[min as usize].y;
        let polemax = polesinit[max as usize].y;
        if (polemax >= (1.0 / tol_poles) * polemin)
            || ((polemin == 0.0) && (polemax >= (1.0 / tol_poles)))
        {
            if polesinit[0].y >= (1.0 / tol_poles) * polesinit[3].y
                || polesinit[0].y <= tol_poles * polesinit[3].y
            {
                panic!("Standard_DimensionError: Hermit Impossible Tolerance");
            }
            if (max == 0) || (max == 3) {
                for pole in polesinit.iter_mut() {
                    pole.y = pole.y - tol_poles * polemax;
                }
            } else if (max == 1) || (max == 2) {
                if (min == 0) || (min == 3) {
                    for pole in polesinit.iter_mut() {
                        pole.y = pole.y - (1.0 / tol_poles) * polemin;
                    }
                } else {
                    if (tol_poles * polemax < polesinit[0].y)
                        && (tol_poles * polemax < polesinit[3].y)
                    {
                        for pole in polesinit.iter_mut() {
                            pole.y = pole.y - tol_poles * polemax;
                        }
                        mark = 1;
                    }

                    if ((1.0 / tol_poles) * polemin > polesinit[0].y)
                        && ((1.0 / tol_poles) * polemin > polesinit[3].y)
                        && (mark == 0)
                    {
                        for pole in polesinit.iter_mut() {
                            pole.y = pole.y - (1.0 / tol_poles) * polemin;
                        }
                        mark = 1;
                    }
                    if mark == 0 {
                        let pole0 = polesinit[0].y;
                        let pole3 = polesinit[3].y;
                        if pole0 < pole3 {
                            let a = (pole3 / pole0).log10();
                            if *boucle == 2 {
                                for pole in polesinit.iter_mut() {
                                    pole.y = pole.y
                                        - (pole3
                                            * (10f64.powf(-0.5 * tol_poles.log10() - a / 2.0)));
                                }
                            } else if *boucle == 1 {
                                for pole in polesinit.iter_mut() {
                                    pole.y = pole.y
                                        - (pole0 * (10f64.powf(a / 2.0 + 0.5 * tol_poles.log10())));
                                }
                                dercas = 1;
                            }
                        }
                        if pole0 > pole3 {
                            let a = (pole0 / pole3).log10();
                            if *boucle == 2 {
                                for pole in polesinit.iter_mut() {
                                    pole.y = pole.y
                                        - (pole0
                                            * (10f64.powf(-0.5 * tol_poles.log10() - a / 2.0)));
                                }
                            } else if *boucle == 1 {
                                for pole in polesinit.iter_mut() {
                                    pole.y = pole.y
                                        - (pole3 * (10f64.powf(a / 2.0 + 0.5 * tol_poles.log10())));
                                }
                                dercas = 1;
                            }
                        }
                    }
                }
            }
        }
    } // end of the loop

    if !sign_denom(&polesinit) {
        // inversion of the polynome sign
        for pole in polesinit.iter_mut() {
            pole.y = -pole.y;
        }
    }

    // positivity loop
    if (polesinit[1].y < 0.0) && (polesinit[2].y >= 0.0) {
        let mut us1 = polesinit[0].y / (polesinit[0].y - polesinit[1].y);
        if *boucle == 2 {
            us1 *= at(&knots, 2);
        }
        if *boucle == 1 {
            if ux != 0.0 {
                us1 *= ux;
            }
        }
        let knots_len = knots.len() as i32;
        locate_parameter_flat(3, &knots, us1, false, 1, knots_len, &mut i1, &mut us1);
        if i1 < 2 {
            *u4 = us1;
        } else {
            *u4 = at(&knots, i1);
        }
    }

    if (polesinit[1].y >= 0.0) && (polesinit[2].y < 0.0) {
        let mut us2 = polesinit[2].y / (polesinit[2].y - polesinit[3].y);
        if *boucle == 2 {
            let knots_len = knots.len() as i32;
            us2 = at(&knots, knots_len - 1) + us2 * (1.0 - at(&knots, knots_len - 1));
        }
        if *boucle == 1 {
            if ux != 0.0 {
                us2 = uy + us2 * (1.0 - uy);
            }
        }
        let knots_len = knots.len() as i32;
        locate_parameter_flat(3, &knots, us2, false, 1, knots_len, &mut i1, &mut us2);
        if i1 >= (knots.len() as i32 - 1) {
            *u5 = us2;
        } else {
            *u5 = at(&knots, i1 + 1);
        }
    }

    if dercas == 1 {
        *boucle += 1;
    }

    if (polesinit[1].y < 0.0) && (polesinit[2].y < 0.0) {
        let mut us1 = polesinit[0].y / (polesinit[0].y - polesinit[1].y);
        let mut us2 = polesinit[2].y / (polesinit[2].y - polesinit[3].y);
        if *boucle != 0 {
            if ux != 0.0 {
                us1 *= ux;
                us2 = uy + us2 * (1.0 - uy);
            }
        }
        let knots_len = knots.len() as i32;
        if us2 <= us1 {
            locate_parameter_flat(3, &knots, us1, false, 1, knots_len, &mut i1, &mut us1);
            if at(&knots, i1) >= us2 {
                // insertion of one knot for the two poles
                *u4 = at(&knots, i1);
            } else {
                if i1 >= 2 {
                    // insertion to the left and to the right without a new knot
                    *u4 = at(&knots, i1);
                    locate_parameter_flat(3, &knots, us2, false, 1, knots_len, &mut i3, &mut us2);
                    if i3 < (bs.nb_knots() - 1) {
                        *u5 = at(&knots, i3 + 1);
                        cas = 1;
                    }
                }
                if cas == 0 {
                    // insertion of only one new knot
                    *u4 = (us1 + us2) / 2.0;
                }
            }
        } else {
            // insertion of two knots
            locate_parameter_flat(3, &knots, us1, false, 1, knots_len, &mut i1, &mut us1);
            if i1 >= 2 {
                *u4 = at(&knots, i1);
            } else {
                *u4 = us1;
            }
            locate_parameter_flat(3, &knots, us2, false, 1, knots_len, &mut i3, &mut us2);
            if i3 < (bs.nb_knots() - 1) {
                *u5 = at(&knots, i3 + 1);
            } else {
                *u5 = us2;
            }
        }
    }
}

/// OCCT static InsertKnots (Hermit.cxx L849-860) — inserts the knots in the
/// BS knot sequence if they are not null.
fn insert_knots(bs: &mut Geom2dBSplineCurve, u4: f64, u5: f64) {
    if u4 != 0.0 {
        // insertion of :0 knot if U4=0, 1 knot if U4=U5, 2 knots otherwise.
        // OCCT BS->InsertKnot(U4): M = 1, ParametricTolerance = 0, Add = false.
        bs.insert_knots(&[u4], &[1], 0.0, false);
    }
    if (u5 != 1.0) && (u5 != u4) {
        bs.insert_knots(&[u5], &[1], 0.0, false);
    }
}

/// OCCT static MovePoles (Hermit.cxx L867-880) — moves the poles above the
/// x axis by raising the non-constrained poles to the first pole level.
/// Value-semantics architecture note: OCCT mutates `BS` in place through
/// SetPole; here the pole array is snapshotted, the same loop runs on it,
/// and the curve is rebuilt (pole 1 is untouched by the OCCT loop, so the
/// snapshot is stable).
fn move_poles(bs: &Geom2dBSplineCurve) -> Geom2dBSplineCurve {
    let nb = bs.nb_poles_curve();
    let y1 = bs.pole(1).y;
    let mut poles: Vec<DVec2> = (1..=nb).map(|i| bs.pole(i)).collect();
    let mut i = 3i32;
    while i <= (nb - 2) {
        // P.SetCoord(1, Pole(i).Coord(1)); P.SetCoord(2, Pole(1).Coord(2)).
        poles[(i - 1) as usize].y = y1;
        i += 1;
    }
    crate::base::geom2d_convert::c1_concat::set_pole_value(bs, poles)
}

/// OCCT Hermit::Solution — the Geom2d overload (Hermit.cxx L994-1098):
/// returns the correct spline a(u) which will be multiplied with BS later.
/// OCCT default arguments are TolPoles = TolKnots = 0.000001 (Hermit.hxx
/// L58-61); pass them explicitly at the call sites.
pub fn solution(bs: &Geom2dBSplineCurve, tol_poles: f64, tol_knots: f64) -> Geom2dBSplineCurve {
    let mut herm = [0.0f64; 4];
    let mut upos1 = 0.0f64;
    let mut upos2 = 1.0f64; // positivity knots
    let mut ux = 0.0f64;
    let mut uy = 1.0f64;
    let mut utol1 = 0.0f64;
    let mut utol2 = 1.0f64; // tolerance knots
    let mut uint1 = 0.0f64;
    let mut uint2 = 1.0f64; // tolerance knots for the first loop
    let mut boucle = 1i32; // loop mark
    let knots = vec![0.0f64, 1.0];
    let multiplicities = vec![4i32, 4];
    let mut zeroboucle = 0i32;

    // computation of the Hermite coefficient.
    hermite_coeff(bs, &mut herm);

    let poles: Vec<DVec2> = vec![
        DVec2::new(0.0, herm[0]),
        DVec2::new(0.0, herm[0] + herm[1] / 3.0),
        DVec2::new(0.0, herm[3] - herm[2] / 3.0),
        DVec2::new(0.0, herm[3]),
    ];
    // creation of the basic BSpline without modif.
    let mut bs1 = Geom2dBSplineCurve::new(poles.clone(), knots.clone(), multiplicities.clone(), 3, false);
    let mut bs2 = Geom2dBSplineCurve::new(poles, knots, multiplicities, 3, false);

    // computation of the positivity knots.
    poly_test(
        &herm,
        bs,
        &mut upos1,
        &mut upos2,
        &mut zeroboucle,
        CONFUSION,
        CONFUSION,
        1.0,
        0.0,
    );
    insert_knots(&mut bs2, upos1, upos2); // and insertion

    if upos1 != 0.0 {
        if upos2 != 1.0 {
            ux = upos1.min(upos2);
            uy = upos1.max(upos2);
        } else {
            ux = upos1;
            uy = upos1;
        }
    } else {
        ux = upos2;
        uy = upos2;
    }

    // computation of the Hermite coefficient on the positive BSpline.
    let nb = bs2.nb_poles_curve();
    herm[0] = bs2.pole(1).y;
    herm[1] = 3.0 * (bs2.pole(2).y - bs2.pole(1).y);
    herm[2] = 3.0 * (bs2.pole(nb).y - bs2.pole(nb - 1).y);
    herm[3] = bs2.pole(nb).y;

    // computation of the tolerance knots.
    poly_test(&herm, bs, &mut utol1, &mut utol2, &mut boucle, tol_poles, tol_knots, ux, uy);
    insert_knots(&mut bs2, utol1, utol2); // and insertion

    if boucle == 2 {
        // insertion of two knots
        let nb = bs2.nb_poles_curve();
        herm[0] = bs2.pole(1).y;
        herm[1] = 3.0 * (bs2.pole(2).y - bs2.pole(1).y);
        herm[2] = 3.0 * (bs2.pole(nb).y - bs2.pole(nb - 1).y);
        herm[3] = bs2.pole(nb).y;
        if utol1 == 0.0 {
            uint2 = utol2;
            poly_test(&herm, bs, &mut utol1, &mut utol2, &mut boucle, tol_poles, tol_knots, uint2, 0.0);
        } else {
            uint1 = utol1;
            poly_test(&herm, bs, &mut utol1, &mut utol2, &mut boucle, tol_poles, tol_knots, uint1, 0.0);
        }
        insert_knots(&mut bs2, utol1, utol2);
    }
    if (bs2.knot(2) < tol_knots) || (bs2.knot(bs2.nb_knots() - 1) > (1.0 - tol_knots)) {
        // checking of the knots tolerance
        panic!("Standard_DimensionError: Hermit Impossible Tolerance");
    } else {
        if (upos2 == 1.0) && (utol2 == 1.0) && (uint2 == 1.0) {
            // test on the final inserted knots
            insert_knots(&mut bs1, bs2.knot(2), 1.0);
        } else {
            if (upos1 == 0.0) && (utol1 == 0.0) && (uint1 == 0.0) {
                insert_knots(&mut bs1, bs2.knot(bs2.nb_knots() - 1), 1.0);
            } else {
                insert_knots(&mut bs1, bs2.knot(bs2.nb_knots() - 1), bs2.knot(2));
            }
        }
        bs1 = move_poles(&bs1); // relocation of the no-contrained knots
    }
    bs1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hand-derived Hermite coefficients for a cubic rational Bezier whose
    /// weight function is w(u) = (1-u)^3*1 + 3u(1-u)^2*2 + 3u^2(1-u)*3 +
    /// u^3*4 (Bernstein form with weights [1,2,3,4]).
    ///
    /// Independent derivation: w(0) = w0 = 1, w'(0) = 3(w1-w0) = 3,
    /// w(1) = w3 = 4, w'(1) = 3(w3-w2) = 3.  HermiteCoeff returns
    /// [1/D0, -D0'/D0^2, -D1'/D1^2, 1/D1] = [1, -3, -3/16, 1/4].
    #[test]
    fn hermite_coeff_bernstein_weights() {
        let poles = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 0.0),
        ];
        let weights = vec![1.0, 2.0, 3.0, 4.0];
        let knots = vec![0.0, 1.0];
        let mults = vec![4, 4];
        let bs = Geom2dBSplineCurve::new_rational(poles, weights, knots, mults, 3, false);
        let mut herm = [0.0f64; 4];
        hermite_coeff(&bs, &mut herm);
        assert!((herm[0] - 1.0).abs() < 1e-15, "herm0={}", herm[0]);
        assert!((herm[1] + 3.0).abs() < 1e-15, "herm1={}", herm[1]);
        assert!((herm[2] + 3.0 / 16.0).abs() < 1e-15, "herm2={}", herm[2]);
        assert!((herm[3] - 0.25).abs() < 1e-15, "herm3={}", herm[3]);
    }

    /// End-to-end Solution on the cubic rational Bezier with weights
    /// [1, 3, 0.5, 2], fully hand-derived from the OCCT control flow:
    ///
    /// Herm = [1, -6, -4.5/4, 1/2] (w(0)=1, w'(0)=6, w(1)=2, w'(1)=4.5), so
    /// the basic Hermite poles are [1, -1, 0.875, 0.5].  The positivity pass
    /// finds pole 1 negative: Us1 = 1/(1+1) = 0.5 -> Upos1 = 0.5 (exact).
    /// The tolerance pass shifts all Hermite poles by TolPoles*Polemax
    /// (eps = 1e-6) and returns Utol1 = (1-eps)/2 (the boucle=1 scaling by
    /// Ux = 0.5) — inserted into BS2.  The final insertion puts the
    /// *updated* BS2 Knot(2) = (1-eps)/2 into BS1.  Inserting
    /// t = (1-eps)/2 into the cubic [1, -1, 0.875, 0.5] gives
    /// Q1 = 1-2t = eps, Q2 = -1/16 - 15eps/16, Q3 = 11/16 + 3eps/16,
    /// Q4 = 1/2, and MovePoles raises pole 3 (i = 3..=nb-2) to the FIRST
    /// pole level (Q0 = 1): y = [1, eps, 1, Q3, 1/2].
    /// Knots/mults: [0, (1-eps)/2, 1] x [4, 1, 4].
    #[test]
    fn solution_inserts_positivity_knot() {
        let eps = 1.0e-6;
        let poles = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(2.0, 0.0),
            DVec2::new(3.0, 0.0),
        ];
        let weights = vec![1.0, 3.0, 0.5, 2.0];
        let bs = Geom2dBSplineCurve::new_rational(poles, weights, vec![0.0, 1.0], vec![4, 4], 3, false);
        let a = solution(&bs, eps, eps);
        assert_eq!(a.degree(), 3);
        assert!(!a.is_rational());
        assert_eq!(a.nb_knots(), 3);
        assert!(a.knot(1).abs() < 1e-15);
        assert!((a.knot(2) - (1.0 - eps) / 2.0).abs() < 1e-12, "knot2={}", a.knot(2));
        assert!((a.knot(3) - 1.0).abs() < 1e-15);
        assert_eq!(a.multiplicity(1), 4);
        assert_eq!(a.multiplicity(2), 1);
        assert_eq!(a.multiplicity(3), 4);
        let expected = [1.0, eps, 1.0, 11.0 / 16.0 + (3.0 * eps) / 16.0, 0.5];
        for (i, &y) in expected.iter().enumerate() {
            let p = a.pole(i as i32 + 1);
            assert!(p.x.abs() < 1e-15, "pole{} x={}", i + 1, p.x);
            assert!(
                (p.y - y).abs() < 1e-12,
                "pole{} y={} expected {}",
                i + 1,
                p.y,
                y
            );
        }
    }
}
