//! OCCT Law package services (TKGeomAlgo/Law) — 1:1 port of Law.hxx
//! (L42-111) and Law.cxx (whole file L29-395), including the file-local
//! static helpers `eval1` (L86-120) and `eval2` (L271-294).
//!
//! Architecture mappings: `NCollection_Array1<double>` / `NCollection_Array1<int>`
//! map to 0-based `Vec<f64>` / `Vec<i32>` slices (the bspl_lib helpers index
//! them with OCCT 1-based semantics); `Adaptor3d_Curve` maps to rcad
//! [`Curve3`]; `handle<Law_BSpline>` maps to `Rc<RefCell<LawBSpline>>`.

use std::cell::RefCell;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::math::bspl_lib::{build_schoenberg_points, interpolate as bspl_interpolate};

use super::law_bspline::LawBSpline;
use super::law_bsp_func::LawBSpFunc;
use super::law_function::LawFunction as _;
use super::law_interpolate::LawInterpolate;
use super::law_linear::LawLinear;

/// OCCT Law::MixBnd(Lin) (Law.cxx L29-46) — builds a 1d bspline that is near
/// from Lin with null derivatives at the extremities.
pub fn law_mix_bnd(lin: &mut LawLinear) -> LawBSpFunc {
    let (mut f, mut l) = (0.0f64, 0.0f64);
    lin.bounds(&mut f, &mut l);
    // OCCT 1-based Knots(1,4) / Mults(1,4) — stored 0-based.
    let knots = [f, 0.75 * f + 0.25 * l, 0.25 * f + 0.75 * l, l];
    let mults = [4i32, 1, 1, 4];
    let pol = law_mix_bnd_poles(3, &knots, &mults, lin);
    let bs = Rc::new(RefCell::new(LawBSpline::new(&pol, &knots, &mults, 3, false)));
    let mut bsf = LawBSpFunc::new();
    // OCCT L44: bsf->SetCurve(bs);
    super::law_s::set_curve(&mut bsf, bs);
    bsf
}

/// OCCT Law::MixBnd(Degree, Knots, Mults, Lin) (Law.cxx L48-84) — builds the
/// poles of the 1d bspline that is near from Lin with null derivatives at
/// the extremities.
pub fn law_mix_bnd_poles(
    degree: usize,
    knots: &[f64],
    mults: &[i32],
    lin: &mut LawLinear,
) -> Vec<f64> {
    let mut nbfk = 0usize;
    for i in 0..mults.len() {
        nbfk += mults[i] as usize;
    }
    // OCCT 1-based fk(1, nbfk): fk(++k) = Knots(i).
    let mut fk: Vec<f64> = Vec::with_capacity(nbfk);
    for i in 0..mults.len() {
        for _j in 1..=mults[i] {
            fk.push(knots[i]);
        }
    }
    let nbpol = nbfk - degree - 1;
    let mut par = vec![0.0f64; nbpol];
    build_schoenberg_points(degree, &fk, &mut par);
    // OCCT: res = new NCollection_HArray1<double>(1, nbpol); pol = res->Array1.
    let mut pol = vec![0.0f64; nbpol];
    for i in 1..=nbpol {
        pol[i - 1] = lin.value(par[i - 1]);
    }
    // OCCT 1-based ord(1, nbpol): ord.Init(0).
    let ord = vec![0i32; nbpol];
    let _inversion_problem = bspl_interpolate(degree, &fk, &par, &ord, 1, &mut pol);
    if nbpol >= 4 {
        pol[1] = pol[0];
        pol[nbpol - 2] = pol[nbpol - 1];
    }
    pol
}

/// OCCT static eval1 (Law.cxx L86-120).
fn eval1(p: f64, first: f64, last: f64, piv: f64, nulr: bool) -> f64 {
    if (nulr && p >= piv) || (!nulr && p <= piv) {
        0.0
    } else if nulr {
        let mut a = piv - first;
        a *= a;
        a = 1.0 / a;
        let mut b = p - first;
        a *= b;
        b = piv - p;
        a *= b;
        a *= b;
        a
    } else {
        let mut a = last - piv;
        a *= a;
        a = 1.0 / a;
        let mut b = last - p;
        a *= b;
        b = p - piv;
        a *= b;
        a *= b;
        a
    }
}

/// OCCT Law::MixTgt(Degree, Knots, Mults, NulOnTheRight, Index)
/// (Law.cxx L122-157) — builds the poles of the 1d bspline that is null on
/// the right side of Knots(Index) (on the left if NulOnTheRight is false).
pub fn law_mix_tgt(
    degree: usize,
    knots: &[f64],
    mults: &[i32],
    nul_on_the_right: bool,
    index: usize, // OCCT 1-based knot index
) -> Vec<f64> {
    let first = knots[0]; // OCCT: Knots(Knots.Lower())
    let last = knots[knots.len() - 1]; // OCCT: Knots(Knots.Upper())
    let piv = knots[index - 1]; // OCCT: Knots(Index)
    let mut nbfk = 0usize;
    for i in 0..mults.len() {
        nbfk += mults[i] as usize;
    }
    let mut fk: Vec<f64> = Vec::with_capacity(nbfk);
    for i in 0..mults.len() {
        for _j in 1..=mults[i] {
            fk.push(knots[i]);
        }
    }
    let nbpol = nbfk - degree - 1;
    let mut par = vec![0.0f64; nbpol];
    build_schoenberg_points(degree, &fk, &mut par);
    let mut pol = vec![0.0f64; nbpol];
    for i in 1..=nbpol {
        pol[i - 1] = eval1(par[i - 1], first, last, piv, nul_on_the_right);
    }
    let ord = vec![0i32; nbpol];
    let _inversion_problem = bspl_interpolate(degree, &fk, &par, &ord, 1, &mut pol);
    pol
}

/// OCCT Adaptor3d_Curve::FirstParameter / LastParameter for a rcad Curve3
/// (the trimmed view reports its own range; otherwise the natural domain).
fn curve_first_parameter(curve: &Curve3) -> f64 {
    match curve {
        Curve3::Trimmed(tc) => tc.first,
        other => other.default_domain()[0],
    }
}

fn curve_last_parameter(curve: &Curve3) -> f64 {
    match curve {
        Curve3::Trimmed(tc) => tc.last,
        other => other.default_domain()[1],
    }
}

/// OCCT Law::Reparametrize (Law.cxx L159-269) — computes a 1d curve to
/// reparametrize a curve: an interpolation of NbPoints points calculated at
/// quasi constant abscissa.
#[allow(clippy::too_many_arguments)]
pub fn law_reparametrize(
    curve: &Curve3,
    first: f64,
    last: f64,
    has_df: bool,
    has_dl: bool,
    d_first: f64,
    d_last: f64,
    rev: bool,
    nb_points: usize,
) -> Rc<RefCell<LawBSpline>> {
    // On evalue la longeur approximative de la courbe.
    let mut dd_first = d_first;
    let mut dd_last = d_last;
    if has_df && rev {
        dd_first = -d_first;
    }
    if has_dl && rev {
        dd_last = -d_last;
    }
    let n2 = 2 * nb_points;
    // OCCT 1-based cumdist(1, 2*NbPoints) / ucourbe(1, 2*NbPoints).
    let mut cumdist = vec![0.0f64; n2];
    let mut ucourbe = vec![0.0f64; n2];
    let u1 = curve_first_parameter(curve);
    let u2 = curve_last_parameter(curve);
    let mut p1: DVec3;
    let mut p2: DVec3;
    let mut u: f64;
    let du: f64;
    let mut length = 0.0f64;
    if !rev {
        p1 = CurveEval::point_at(curve, u1);
        u = u1;
        du = (u2 - u1) / (2 * nb_points - 1) as f64;
    } else {
        p1 = CurveEval::point_at(curve, u2);
        u = u2;
        du = (u1 - u2) / (2 * nb_points - 1) as f64;
    }
    for i in 1..=n2 {
        p2 = CurveEval::point_at(curve, u);
        length += p1.distance(p2);
        cumdist[i - 1] = length;
        ucourbe[i - 1] = u;
        u += du;
        p1 = p2;
    }
    if rev {
        ucourbe[n2 - 1] = u1;
    } else {
        ucourbe[n2 - 1] = u2;
    }

    // OCCT 1-based point(1, NbPoints) / param(1, NbPoints).
    let mut point = vec![0.0f64; nb_points];
    let mut param = vec![0.0f64; nb_points];

    let d_corde = length / (nb_points - 1) as f64;
    let mut corde = d_corde;
    let mut index = 1usize; // OCCT 1-based
    let fac = 1.0 / (nb_points - 1) as f64;

    point[0] = ucourbe[0];
    param[0] = first;
    point[nb_points - 1] = ucourbe[n2 - 1];
    param[nb_points - 1] = last;

    for i in 2..nb_points {
        while cumdist[index - 1] < corde {
            index += 1;
        }
        let alpha = (corde - cumdist[index - 2]) / (cumdist[index - 1] - cumdist[index - 2]);
        u = ucourbe[index - 2] + alpha * (ucourbe[index - 1] - ucourbe[index - 2]);
        point[i - 1] = u;
        param[i - 1] = ((nb_points - i) as f64 * first + (i - 1) as f64 * last) * fac;
        corde = i as f64 * d_corde;
    }
    let mut inter = LawInterpolate::with_parameters(point, param, false, 1.0e-9);
    if has_df || has_dl {
        let mut tgt = vec![0.0f64; nb_points];
        let mut flag = vec![false; nb_points];
        if has_df {
            flag[0] = true;
            tgt[0] = dd_first;
        }
        if has_dl {
            flag[nb_points - 1] = true;
            tgt[nb_points - 1] = dd_last;
        }
        inter.load(&tgt, &flag);
    }
    inter.perform();
    if !inter.is_done() {
        panic!("Law::Reparametrize echec interpolation");
    }
    inter.curve()
}

/// OCCT static eval2 (Law.cxx L271-294) — the `first` / `last` parameters are
/// unnamed in OCCT (commented out in the definition) and are kept as unused
/// placeholders here.
#[allow(clippy::too_many_arguments)]
fn eval2(
    p: f64,
    _first: f64,
    _last: f64,
    mil: f64,
    hasfirst: bool,
    haslast: bool,
    bs1: &Option<Rc<RefCell<LawBSpline>>>,
    bs2: &Option<Rc<RefCell<LawBSpline>>>,
) -> f64 {
    if hasfirst && p < mil {
        bs1.as_ref().unwrap().borrow().value(p)
    } else if haslast && p > mil {
        bs2.as_ref().unwrap().borrow().value(p)
    } else {
        1.0
    }
}

/// OCCT Law::Scale (Law.cxx L296-352) — computes a 1d curve to scale a field
/// of tangency.  Value is 1. for t = (First+Last)/2.
pub fn law_scale(
    first: f64,
    last: f64,
    has_f: bool,
    has_l: bool,
    v_first: f64,
    v_last: f64,
) -> Rc<RefCell<LawBSpline>> {
    let milieu = 0.5 * (first + last);
    // OCCT 1-based knot(1,3) / fknot(1,10) / mult(1,3).
    let knot = [first, milieu, last];
    let mut fknot = vec![0.0f64; 10];
    let mult = [4i32, 2, 4];
    fknot[0] = first;
    fknot[1] = first;
    fknot[2] = first;
    fknot[3] = first;
    fknot[9] = last;
    fknot[8] = last;
    fknot[7] = last;
    fknot[6] = last;
    fknot[4] = milieu;
    fknot[5] = milieu;

    // OCCT 1-based pbs(1,4) / kbs(1,2) / mbs(1,2).
    let mut pbs = [0.0f64; 4];
    let mut bs1: Option<Rc<RefCell<LawBSpline>>> = None;
    let mut bs2: Option<Rc<RefCell<LawBSpline>>> = None;
    if has_f {
        pbs[0] = v_first;
        pbs[1] = v_first;
        pbs[2] = 1.0;
        pbs[3] = 1.0;
        let kbs = [first, milieu];
        let mbs = [4i32, 4];
        bs1 = Some(Rc::new(RefCell::new(LawBSpline::new(
            &pbs, &kbs, &mbs, 3, false,
        ))));
    }
    if has_l {
        pbs[0] = 1.0;
        pbs[1] = 1.0;
        pbs[2] = v_last;
        pbs[3] = v_last;
        let kbs = [milieu, last];
        let mbs = [4i32, 4];
        bs2 = Some(Rc::new(RefCell::new(LawBSpline::new(
            &pbs, &kbs, &mbs, 3, false,
        ))));
    }

    // OCCT 1-based pol(1,6) / par(1,6).
    let mut pol = vec![0.0f64; 6];
    let mut par = vec![0.0f64; 6];
    build_schoenberg_points(3, &fknot, &mut par);
    for i in 1..=6 {
        pol[i - 1] = eval2(
            par[i - 1], first, last, milieu, has_f, has_l, &bs1, &bs2,
        );
    }
    let ord = vec![0i32; 6];
    let _inversion_problem = bspl_interpolate(3, &fknot, &par, &ord, 1, &mut pol);
    let bs1 = Rc::new(RefCell::new(LawBSpline::new(&pol, &knot, &mult, 3, false)));
    bs1
}

/// OCCT Law::ScaleCub (Law.cxx L354-395).
pub fn law_scale_cub(
    first: f64,
    last: f64,
    has_f: bool,
    has_l: bool,
    v_first: f64,
    v_last: f64,
) -> Rc<RefCell<LawBSpline>> {
    let milieu = 0.5 * (first + last);
    // OCCT 1-based pol(1,5) / knot(1,3) / mult(1,3).
    let mut pol = vec![0.0f64; 5];
    let knot = [first, milieu, last];
    let mult = [4i32, 1, 4];

    if has_f {
        pol[0] = v_first;
        pol[1] = v_first;
    } else {
        pol[0] = 1.0;
        pol[1] = 1.0;
    }
    if has_l {
        pol[3] = v_last;
        pol[4] = v_last;
    } else {
        pol[3] = 1.0;
        pol[4] = 1.0;
    }

    pol[2] = 1.0;

    let bs = Rc::new(RefCell::new(LawBSpline::new(&pol, &knot, &mult, 3, false)));
    bs
}
