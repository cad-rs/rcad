//! OCCT GeomFill_Frenet (TKGeomAlgo/GeomFill) — 1:1 port of GeomFill_Frenet.hxx
//! (L32-108) + GeomFill_Frenet.cxx (whole file L25-1108).
//!
//! Architecture mappings: `Adaptor3d_Curve` -> rcad `Curve3`; the trimmed
//! law curve (`myTrimmed`, a `GeomAdaptor` on a `Geom_TrimmedCurve` with
//! Adjust=False) is a `Curve3::Trimmed` evaluated at unchanged parameters.
//! `Extrema_ExtPC` on the SnglrFunc-as-curve -> rcad `ExtPC::new_fn` with
//! the SnglrFunc evaluation closure.

use std::f64::consts::PI;

use glam::DVec3;

use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::geom::{Curve3, CurveEval, TrimmedCurve3};
use rcad_kernel::math::gp::Ax2;
use rcad_kernel::math::GeomAbsShape;

use super::sngrl_func::SnglrFunc;
use super::trihedron_law::{curve_first_parameter, curve_last_parameter, TrihedronLaw, TrihedronLawBase};

// OCCT statics (GeomFill_Frenet.cxx L40-42).
const NULL_TOL: f64 = 1.0e-10;
const MAX_SINGULAR: f64 = 1.0e-5;
const MAX_DERIV_ORDER: i32 = 3;

// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 2.2250738585072014e-308;

/// OCCT static FDeriv (L48-56): computes (F/|F|)'.
fn f_deriv(f: DVec3, df: DVec3) -> DVec3 {
    let norma = f.length();
    (df - f * (f.dot(df)) / (norma * norma)) / norma
}

/// OCCT static DDeriv (L62-75): computes (F/|F|)''.
fn d_deriv(f: DVec3, df: DVec3, d2f: DVec3) -> DVec3 {
    let norma = f.length();
    (d2f - 2.0 * df * (f.dot(df)) / (norma * norma)) / norma
        - f * ((df.length_squared() + f.dot(d2f) - 3.0 * (f.dot(df)) * (f.dot(df)) / (norma * norma))
            / (norma * norma * norma))
}

/// OCCT static CosAngle (L81-102): cosine between two vectors.
fn cos_angle(v1: DVec3, v2: DVec3) -> f64 {
    let a_tol = GP_RESOLUTION;
    let m1 = v1.length();
    let m2 = v2.length();
    if m1 <= a_tol || m2 <= a_tol {
        // Vectors are codirectional
        return 1.0;
    }
    let mut a_cang = v1.dot(v2) / (m1 * m2);
    if a_cang > 1.0 {
        a_cang = 1.0;
    }
    if a_cang < -1.0 {
        a_cang = -1.0;
    }
    a_cang
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488 -> gp_Dir.cxx L27-53).
fn gp_vec_angle(a: DVec3, b: DVec3) -> f64 {
    assert!(
        a.length() > GP_RESOLUTION && b.length() > GP_RESOLUTION,
        "gp_VectorWithNullMagnitude"
    );
    let da = a.normalize_or_zero();
    let db = b.normalize_or_zero();
    let cosinus = da.dot(db);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = da.cross(db).length();
        if cosinus < 0.0 {
            PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT Adaptor3d_Curve::NbIntervals for the continuity mapping used by
/// Frenet (GeomAbs_C2/C3/CN): only the BSpline base has subdivisions; the
/// elementary curves have one interval.
fn curve_nb_intervals(c: &Curve3, s: GeomAbsShape) -> usize {
    match c {
        Curve3::BSpline(bs) => {
            let _ = s;
            bs.c2_intervals().len() - 1
        }
        Curve3::Trimmed(tc) => curve_nb_intervals(&tc.curve, s),
        _ => 1,
    }
}

/// OCCT Adaptor3d_Curve::Intervals.
fn curve_intervals(c: &Curve3, s: GeomAbsShape) -> Vec<f64> {
    match c {
        Curve3::BSpline(bs) => {
            let _ = s;
            bs.c2_intervals()
        }
        Curve3::Trimmed(tc) => {
            // The trimmed curve keeps the base parametrization; OCCT reports
            // the trim bounds as the outer interval.
            let inner = curve_intervals(&tc.curve, s);
            let mut out: Vec<f64> = inner
                .into_iter()
                .filter(|&t| t > tc.first && t < tc.last)
                .collect();
            out.insert(0, tc.first);
            out.push(tc.last);
            out
        }
        _ => vec![curve_first_parameter(c), curve_last_parameter(c)],
    }
}

/// OCCT myTrimmed->D2 / D3 / DN evaluations — the trimmed view evaluates the
/// base curve at unchanged parameters.
fn law_d2(c: &Curve3, u: f64) -> (DVec3, DVec3, DVec3) {
    match c {
        Curve3::BSpline(bs) => {
            let p = bs.point_at(u);
            let d1 = bs.dn(u, 1);
            let d2 = bs.dn(u, 2);
            (p, d1, d2)
        }
        Curve3::Trimmed(tc) => law_d2(&tc.curve, u),
        other => {
            let p = other.point_at(u);
            let d1 = other.derivative_at(u);
            let h = 1e-6;
            let d2 = (other.derivative_at(u + h) - other.derivative_at(u - h)) / (2.0 * h);
            (p, d1, d2)
        }
    }
}

fn law_d3(c: &Curve3, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
    match c {
        Curve3::BSpline(bs) => {
            let p = bs.point_at(u);
            let d1 = bs.dn(u, 1);
            let d2 = bs.dn(u, 2);
            let d3 = bs.dn(u, 3);
            (p, d1, d2, d3)
        }
        Curve3::Trimmed(tc) => law_d3(&tc.curve, u),
        _ => unimplemented!("GeomFill_Frenet D3 for non-BSpline bases is anchor-out-of-scope"),
    }
}

fn law_dn(c: &Curve3, u: f64, n: usize) -> DVec3 {
    match c {
        Curve3::BSpline(bs) => bs.dn(u, n),
        Curve3::Trimmed(tc) => law_dn(&tc.curve, u, n),
        _ => unimplemented!("GeomFill_Frenet DN for non-BSpline bases is anchor-out-of-scope"),
    }
}

fn law_d0(c: &Curve3, u: f64) -> DVec3 {
    c.point_at(u)
}

/// OCCT GeomFill_Frenet.
#[derive(Debug, Clone)]
pub struct Frenet {
    pub(crate) base: TrihedronLawBase,
    /// OCCT member `P` — written by every D0/D1/D2 in OCCT, never read on
    /// these paths (the field is kept for form).
    _p: DVec3,
    my_sngl: Option<Vec<f64>>,
    my_sngl_len: Option<Vec<f64>>,
    is_sngl: bool,
}

impl Frenet {
    /// OCCT GeomFill_Frenet() (L105-108).
    pub fn new() -> Self {
        Frenet {
            base: TrihedronLawBase::default(),
            _p: DVec3::ZERO,
            my_sngl: None,
            my_sngl_len: None,
            is_sngl: false,
        }
    }

    /// OCCT Init (L138-355) — searches the curve singularities.
    pub fn init(&mut self) {
        let curve = self.base.my_curve.as_ref().unwrap().clone();
        let mut func = SnglrFunc::new(curve.clone());
        const TOL_F: f64 = 1.0e-10;
        const TOL: f64 = 10.0 * TOL_F;
        const TOL2: f64 = TOL * TOL;
        const PTOL: f64 = 1e-12; // Precision::PConfusion()

        // We want to determine if the curve has linear segments
        let nb_int_c2 = curve_nb_intervals(&curve, GeomAbsShape::C2);
        let my_c2_disc = curve_intervals(&curve, GeomAbsShape::C2);
        let mut is_lin = vec![true; nb_int_c2];
        let mut is_const = vec![true; nb_int_c2];
        let mut ave_func = vec![0.0f64; nb_int_c2];
        let nb_control = 10usize;
        let mut average = 0.0f64;
        for i in 1..=nb_int_c2 {
            let step = (my_c2_disc[i] - my_c2_disc[i - 1]) / nb_control as f64;
            is_lin[i - 1] = true;
            is_const[i - 1] = true;
            let mut c1 = DVec3::ZERO;
            for j in 1..=nb_control {
                let c = func.eval_d0(my_c2_disc[i - 1] + (j - 1) as f64 * step);
                if j == 1 {
                    c1 = c;
                }
                let modulus = c.length();
                if modulus > TOL {
                    is_lin[i - 1] = false;
                }
                average += modulus;
                if is_const[i - 1]
                    && ((c.x - c1.x).abs() > TOL
                        || (c.y - c1.y).abs() > TOL
                        || (c.z - c1.z).abs() > TOL)
                {
                    is_const[i - 1] = false;
                }
            }
            ave_func[i - 1] = average / nb_control as f64;
        }

        // Here we are looking for singularities
        let mut seq_array: Vec<Vec<f64>> = vec![Vec::new(); nb_int_c2];
        let mut sngl_seq: Vec<f64> = Vec::new();
        let origin = DVec3::ZERO;

        for i in 1..=nb_int_c2 {
            if !is_lin[i - 1] && !is_const[i - 1] {
                func.set_ratio(1.0 / ave_func[i - 1]); // Normalization
                let value = |u: f64| func.eval_d0(u);
                let ext = ExtPC::new_fn(
                    origin,
                    TOL_F,
                    my_c2_disc[i - 1],
                    my_c2_disc[i],
                    &value,
                );
                if ext.is_done() && ext.nb_ext() != 0 {
                    for j in 1..=ext.nb_ext() {
                        let value2 = ext.square_distance(j);
                        if value2 < TOL2 {
                            let t = ext.point(j).param;
                            seq_array[i - 1].push(t);
                        }
                    }
                }
                // sorting
                if !seq_array[i - 1].is_empty() {
                    seq_array[i - 1].sort_by(|a, b| a.partial_cmp(b).unwrap());
                }
            }
        }

        // Filling SnglSeq by first sets of roots
        for seq in &seq_array {
            for &t in seq {
                sngl_seq.push(t);
            }
        }

        // Extrema works bad, need to pass second time
        for i in 0..nb_int_c2 {
            if !seq_array[i].is_empty() {
                let mut local = seq_array[i].clone();
                local.insert(0, my_c2_disc[i]);
                local.push(my_c2_disc[i + 1]);
                func.set_ratio(1.0 / ave_func[i]);
                for j in 0..local.len() - 1 {
                    if local[j + 1] - local[j] > PTOL {
                        let value = |u: f64| func.eval_d0(u);
                        let ext = ExtPC::new_fn(origin, TOL_F, local[j], local[j + 1], &value);
                        if ext.is_done() {
                            for k in 1..=ext.nb_ext() {
                                let value2 = ext.square_distance(k);
                                if value2 < TOL2 {
                                    let t = ext.point(k).param;
                                    if t - local[j] > PTOL && local[j + 1] - t > PTOL {
                                        sngl_seq.push(t);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if !sngl_seq.is_empty() {
            // sorting
            sngl_seq.sort_by(|a, b| a.partial_cmp(b).unwrap());

            // discard repeating elements
            let mut found = true;
            let mut j = 0usize;
            while found {
                found = false;
                let mut i = j;
                while i + 1 < sngl_seq.len() {
                    if sngl_seq[i + 1] - sngl_seq[i] <= PTOL {
                        sngl_seq.remove(i + 1);
                        j = i;
                        found = true;
                        break;
                    }
                    i += 1;
                }
            }

            let my_sngl = sngl_seq.clone();

            // computation of length of singular interval
            let mut my_sngl_len = vec![0.0f64; my_sngl.len()];
            for i in 1..=my_sngl.len() {
                let (_, sngl_der, sngl_der2) = law_d2(&curve, my_sngl[i - 1]);
                let norm = sngl_der.length();
                if norm > GP_RESOLUTION {
                    my_sngl_len[i - 1] = (NULL_TOL / norm).min(MAX_SINGULAR);
                } else {
                    let norm = sngl_der2.length();
                    if norm > GP_RESOLUTION {
                        my_sngl_len[i - 1] = ((2.0 * NULL_TOL / norm).sqrt()).min(MAX_SINGULAR);
                    } else {
                        my_sngl_len[i - 1] = MAX_SINGULAR;
                    }
                }
            }

            if my_sngl.len() > 1 {
                // we have to merge singular points that have common parts of
                // singular intervals
                let mut tmp_seq: Vec<(f64, f64)> = Vec::new();
                tmp_seq.push((my_sngl[0], my_sngl_len[0]));
                for i in 1..my_sngl.len() {
                    let u12 = tmp_seq.last().unwrap().0 + tmp_seq.last().unwrap().1;
                    let u21 = my_sngl[i] - my_sngl_len[i];
                    if u12 >= u21 {
                        let u11 = tmp_seq.last().unwrap().0 - tmp_seq.last().unwrap().1;
                        let u22 = my_sngl[i] + my_sngl_len[i];
                        let len = tmp_seq.len();
                        tmp_seq[len - 1] = ((u11 + u22) / 2.0, (u22 - u11) / 2.0);
                    } else {
                        tmp_seq.push((my_sngl[i], my_sngl_len[i]));
                    }
                }
                self.my_sngl = Some(tmp_seq.iter().map(|&(x, _)| x).collect());
                self.my_sngl_len = Some(tmp_seq.iter().map(|&(_, y)| y).collect());
            } else {
                self.my_sngl = Some(my_sngl);
                self.my_sngl_len = Some(my_sngl_len);
            }
            self.is_sngl = true;
        } else {
            self.is_sngl = false;
        }
    }

    /// OCCT RotateTrihedron (L361-470) — revolves the trihedron to coincide
    /// "Tangent" and "NewTangent" axes.
    fn rotate_trihedron(
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
        new_tangent: DVec3,
    ) -> bool {
        let an_inf_cos = 1.0e-12f64.cos(); // cos(Precision::Angular())
        let a_tol = GP_RESOLUTION;

        let mut an_axis = tangent.cross(new_tangent);
        let nt = an_axis.length();
        if nt <= a_tol {
            // No rotation required
            return true;
        }
        an_axis /= nt; // Normalization

        let a_px = an_axis.x;
        let a_py = an_axis.y;
        let a_pz = an_axis.z;
        let a_cang = cos_angle(*tangent, new_tangent); // cosine

        let an_add_cang = 1.0 - a_cang;
        let a_sang = (1.0 - a_cang * a_cang).sqrt(); // sine

        let a_v11 = DVec3::new(
            an_add_cang * a_px * a_px + a_cang,
            an_add_cang * a_px * a_py - a_pz * a_sang,
            an_add_cang * a_px * a_pz + a_py * a_sang,
        );
        let a_v12 = DVec3::new(
            an_add_cang * a_px * a_px + a_cang,
            an_add_cang * a_px * a_py + a_pz * a_sang,
            an_add_cang * a_px * a_pz - a_py * a_sang,
        );
        let a_v21 = DVec3::new(
            an_add_cang * a_px * a_py + a_pz * a_sang,
            an_add_cang * a_py * a_py + a_cang,
            an_add_cang * a_py * a_pz - a_px * a_sang,
        );
        let a_v22 = DVec3::new(
            an_add_cang * a_px * a_py - a_pz * a_sang,
            an_add_cang * a_py * a_py + a_cang,
            an_add_cang * a_py * a_pz + a_px * a_sang,
        );
        let a_v31 = DVec3::new(
            an_add_cang * a_px * a_pz - a_py * a_sang,
            an_add_cang * a_py * a_pz + a_px * a_sang,
            an_add_cang * a_pz * a_pz + a_cang,
        );
        let a_v32 = DVec3::new(
            an_add_cang * a_px * a_pz + a_py * a_sang,
            an_add_cang * a_py * a_pz - a_px * a_sang,
            an_add_cang * a_pz * a_pz + a_cang,
        );

        let a_t1 = DVec3::new(
            tangent.dot(a_v11),
            tangent.dot(a_v21),
            tangent.dot(a_v31),
        );
        let a_t2 = DVec3::new(
            tangent.dot(a_v12),
            tangent.dot(a_v22),
            tangent.dot(a_v32),
        );

        if cos_angle(a_t1, new_tangent) >= cos_angle(a_t2, new_tangent) {
            *tangent = a_t1;
            *normal = DVec3::new(normal.dot(a_v11), normal.dot(a_v21), normal.dot(a_v31));
            *binormal = DVec3::new(binormal.dot(a_v11), binormal.dot(a_v21), binormal.dot(a_v31));
        } else {
            *tangent = a_t2;
            *normal = DVec3::new(normal.dot(a_v12), normal.dot(a_v22), normal.dot(a_v32));
            *binormal = DVec3::new(binormal.dot(a_v12), binormal.dot(a_v22), binormal.dot(a_v32));
        }

        cos_angle(*tangent, new_tangent) >= an_inf_cos
    }

    /// OCCT IsSingular (L993-1007).
    fn is_singular(&self, u: f64, index: &mut usize) -> bool {
        if !self.is_sngl {
            return false;
        }
        let my_sngl = self.my_sngl.as_ref().unwrap();
        let my_sngl_len = self.my_sngl_len.as_ref().unwrap();
        for i in 1..=my_sngl.len() {
            if (u - my_sngl[i - 1]).abs() < my_sngl_len[i - 1] {
                *index = i;
                return true;
            }
        }
        false
    }

    /// OCCT DoSingular (L1009-1081).
    #[allow(clippy::too_many_arguments)]
    fn do_singular(
        &self,
        u: f64,
        index: usize,
        tangent: &mut DVec3,
        binormal: &mut DVec3,
        n: &mut usize,
        k: &mut usize,
        tflag: &mut i32,
        bnflag: &mut i32,
        delta: &mut f64,
    ) -> bool {
        let max_n = 20usize;
        *delta = 0.0;
        let my_sngl_len = self.my_sngl_len.as_ref().unwrap();
        let mut h = 2.0 * my_sngl_len[index - 1];

        let (mut a, mut b) = (0.0f64, 0.0f64);
        let mut t = DVec3::ZERO;
        let mut nn = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        *tflag = 1;
        *bnflag = 1;
        self.get_interval(&mut a, &mut b);
        if u >= (a + b) / 2.0 {
            h = -h;
        }
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let mut i = 1usize;
        while i <= max_n {
            *tangent = law_dn(trimmed, u, i);
            if tangent.length() > 1e-7 {
                break;
            }
            i += 1;
        }
        if i > max_n {
            return false;
        }
        *tangent = tangent.normalize();
        *n = i;

        let mut i = i + 1;
        while i <= max_n {
            *binormal = tangent.cross(law_dn(trimmed, u, i));
            let magn = binormal.length();
            if magn > 1e-7 {
                // modified by jgv, 12.08.03 for OCC605
                let next_binormal = tangent.cross(law_dn(trimmed, u, i + 1));
                if next_binormal.length() > magn {
                    i += 1;
                    *binormal = next_binormal;
                }
                break;
            }
            i += 1;
        }
        if i > max_n {
            *delta = h;
            return false;
        }

        *binormal = binormal.normalize();
        *k = i;

        self.d0(u + h, &mut t, &mut nn, &mut bn);

        if gp_vec_angle(*tangent, t) > PI / 2.0 {
            *tflag = -1;
        }
        if gp_vec_angle(*binormal, bn) > PI / 2.0 {
            *bnflag = -1;
        }

        true
    }

    /// OCCT SingularD0 (L1083-1098).
    fn singular_d0(
        &self,
        param: f64,
        index: usize,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
        delta: &mut f64,
    ) -> bool {
        let mut n = 0usize;
        let mut k = 0usize;
        let mut tflag = 0i32;
        let mut bnflag = 0i32;
        if !self.do_singular(
            param,
            index,
            tangent,
            binormal,
            &mut n,
            &mut k,
            &mut tflag,
            &mut bnflag,
            delta,
        ) {
            return false;
        }
        let _ = (n, k);

        *tangent *= tflag as f64;
        *binormal *= bnflag as f64;
        *normal = *binormal;
        *normal = normal.cross(*tangent);

        true
    }

    /// OCCT SingularD1 (L1100-1142).
    #[allow(clippy::too_many_arguments)]
    fn singular_d1(
        &self,
        param: f64,
        index: usize,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        delta: &mut f64,
    ) -> bool {
        let mut n = 0usize;
        let mut k = 0usize;
        let mut tflag = 0i32;
        let mut bnflag = 0i32;
        if !self.do_singular(
            param,
            index,
            tangent,
            binormal,
            &mut n,
            &mut k,
            &mut tflag,
            &mut bnflag,
            delta,
        ) {
            return false;
        }

        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let f = law_dn(trimmed, param, n);
        let df = law_dn(trimmed, param, n + 1);
        *dtangent = f_deriv(f, df);

        let dtmp = law_dn(trimmed, param, k);
        let f = tangent.cross(dtmp);
        let df = dtangent.cross(dtmp) + tangent.cross(law_dn(trimmed, param, k + 1));
        *dbinormal = f_deriv(f, df);

        if tflag < 0 {
            *tangent = -*tangent;
            *dtangent = -*dtangent;
        }
        if bnflag < 0 {
            *binormal = -*binormal;
            *dbinormal = -*dbinormal;
        }

        *normal = binormal.cross(*tangent);
        *dnormal = dbinormal.cross(*tangent) + binormal.cross(*dtangent);

        true
    }

    /// OCCT SingularD2 (L1144-1205).
    #[allow(clippy::too_many_arguments)]
    fn singular_d2(
        &self,
        param: f64,
        index: usize,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        d2tangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        d2normal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        d2binormal: &mut DVec3,
        delta: &mut f64,
    ) -> bool {
        let mut n = 0usize;
        let mut k = 0usize;
        let mut tflag = 0i32;
        let mut bnflag = 0i32;
        if !self.do_singular(
            param,
            index,
            tangent,
            binormal,
            &mut n,
            &mut k,
            &mut tflag,
            &mut bnflag,
            delta,
        ) {
            return false;
        }

        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let f = law_dn(trimmed, param, n);
        let df = law_dn(trimmed, param, n + 1);
        let d2f = law_dn(trimmed, param, n + 2);
        *dtangent = f_deriv(f, df);
        *d2tangent = d_deriv(f, df, d2f);

        let dtmp1 = law_dn(trimmed, param, k);
        let dtmp2 = law_dn(trimmed, param, k + 1);
        let f = tangent.cross(dtmp1);
        let df = dtangent.cross(dtmp1) + tangent.cross(dtmp2);
        let d2f = d2tangent.cross(dtmp1) + 2.0 * dtangent.cross(dtmp2)
            + tangent.cross(law_dn(trimmed, param, k + 2));
        *dbinormal = f_deriv(f, df);
        *d2binormal = d_deriv(f, df, d2f);

        if tflag < 0 {
            *tangent = -*tangent;
            *dtangent = -*dtangent;
            *d2tangent = -*d2tangent;
        }
        if bnflag < 0 {
            *binormal = -*binormal;
            *dbinormal = -*dbinormal;
            *d2binormal = -*d2binormal;
        }

        *normal = binormal.cross(*tangent);
        *dnormal = dbinormal.cross(*tangent) + binormal.cross(*dtangent);
        *d2normal =
            d2binormal.cross(*tangent) + 2.0 * dbinormal.cross(*dtangent) + binormal.cross(*d2tangent);

        true
    }
}

impl Default for Frenet {
    fn default() -> Self {
        Self::new()
    }
}

impl TrihedronLaw for Frenet {
    fn my_curve(&self) -> &Option<Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<Curve3>) {
        self.base.my_trimmed = c;
    }

    /// OCCT Copy (L110-125).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let mut copy = Frenet::new();
        if let Some(curve) = &self.base.my_curve {
            TrihedronLaw::set_curve(&mut copy, curve.clone());
        }
        Box::new(copy)
    }

    /// OCCT SetCurve (L127-154).
    fn set_curve(&mut self, c: Curve3) -> bool {
        TrihedronLaw::set_curve(self, c.clone());
        // GeomAbs_Circle/Ellipse/Hyperbola/Parabola/Line — no problem;
        // the other types need a singularity search.
        let analytic = matches!(
            c,
            Curve3::Line(_)
                | Curve3::Circle(_)
                | Curve3::Ellipse(_)
                | Curve3::Hyperbola(_)
                | Curve3::Parabola(_)
        );
        if analytic {
            self.is_sngl = false;
        } else {
            // We have to search singularities
            self.init();
        }
        true
    }

    /// OCCT SetInterval (base) + trimmed view (Adjust = False).
    fn set_interval(&mut self, first: f64, last: f64) {
        let curve = self
            .base
            .my_curve
            .as_ref()
            .expect("GeomFill_TrihedronLaw::SetInterval with null curve")
            .clone();
        self.base.my_trimmed =
            Some(Curve3::Trimmed(TrimmedCurve3::new(curve, first, last)));
    }

    /// OCCT D0 (L472-628).  The near-singular branch recursively calls D0
    /// at a shifted parameter (OCCT L596-610) — the recursion terminates
    /// there through the non-degenerate branch.
    #[allow(unused_mut)]
    fn d0(&self, the_param: f64, tangent: &mut DVec3, normal: &mut DVec3, binormal: &mut DVec3) -> bool {
        let a_tol = GP_RESOLUTION;
        let mut index = 0usize;
        let mut delta = 0.0f64;
        if self.is_singular(the_param, &mut index)
            && self.singular_d0(the_param, index, tangent, normal, binormal, &mut delta)
        {
            return true;
        }

        let a_param = the_param + delta;
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let (_p, mut tg, mut bn) = law_d2(trimmed, a_param);

        let division_factor = 1.0e-3;
        let an_uinfium = curve_first_parameter(trimmed);
        let an_usupremum = curve_last_parameter(trimmed);
        let a_delta = (an_usupremum - an_uinfium) * division_factor;
        let mut ndu = tg.length();

        if ndu <= a_tol {
            let mut a_tn;
            // Derivative is approximated by Taylor-series
            let mut an_index: i32 = 1; // Derivative order
            let mut is_derive_found = false;
            loop {
                an_index += 1;
                a_tn = law_dn(trimmed, the_param, an_index as usize);
                ndu = a_tn.length();
                is_derive_found = ndu > a_tol;
                if is_derive_found || an_index >= MAX_DERIV_ORDER {
                    break;
                }
            }

            if is_derive_found {
                let u;
                if the_param - an_uinfium < a_delta {
                    u = the_param + a_delta;
                } else {
                    u = the_param - a_delta;
                }
                let p1 = law_d0(trimmed, the_param.min(u));
                let p2 = law_d0(trimmed, the_param.max(u));
                let v1 = p2 - p1;
                let a_dir_factor = a_tn.dot(v1);
                if a_dir_factor < 0.0 {
                    a_tn = -a_tn;
                }
            } else {
                // Derivative is approximated by three points
                let p1: DVec3;
                let p2: DVec3;
                let p3: DVec3;
                let is_parameter_grown;
                if the_param - an_uinfium < 2.0 * a_delta {
                    p1 = law_d0(trimmed, the_param);
                    p2 = law_d0(trimmed, the_param + a_delta);
                    p3 = law_d0(trimmed, the_param + 2.0 * a_delta);
                    is_parameter_grown = true;
                } else {
                    p1 = law_d0(trimmed, the_param - 2.0 * a_delta);
                    p2 = law_d0(trimmed, the_param - a_delta);
                    p3 = law_d0(trimmed, the_param);
                    is_parameter_grown = false;
                }
                let ptemp = DVec3::ZERO;
                let v1 = p1 - ptemp;
                let v2 = p2 - ptemp;
                let v3 = p3 - ptemp;
                a_tn = if is_parameter_grown {
                    -3.0 * v1 + 4.0 * v2 - v3
                } else {
                    v1 - 4.0 * v2 + 3.0 * v3
                };
            }
            // Recursive calling is used to determine the trihedron for a
            // point which is near to the given one.
            let ok = if the_param - an_uinfium < 10.0 * a_delta {
                self.d0(a_param + 10.0 * a_delta, tangent, normal, binormal)
            } else {
                self.d0(a_param - 10.0 * a_delta, tangent, normal, binormal)
            };
            if !ok {
                return false;
            }

            if !Self::rotate_trihedron(tangent, normal, binormal, a_tn) {
                return false;
            }
        } else {
            *tangent = tg / tg.length();
            *binormal = tangent.cross(bn);
            let norm = binormal.length();
            if norm <= GP_RESOLUTION {
                let axe = Ax2::from_direction(DVec3::ZERO, *tangent);
                *binormal = axe.y_direction;
            } else {
                *binormal = binormal.normalize();
            }
            *normal = *binormal;
            *normal = normal.cross(*tangent);
        }

        true
    }

    /// OCCT D1 (L630-682).
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
        let mut index = 0usize;
        let mut delta = 0.0f64;
        if self.is_singular(param, &mut index)
            && self.singular_d1(
                param,
                index,
                tangent,
                dtangent,
                normal,
                dnormal,
                binormal,
                dbinormal,
                &mut delta,
            )
        {
            return true;
        }

        let the_param = param + delta;
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let (_p, dc1, dc2, dc3) = law_d3(trimmed, the_param);
        *tangent = dc1.normalize();

        if tangent.cross(dc2).length() <= GP_RESOLUTION {
            let axe = Ax2::from_direction(DVec3::ZERO, *tangent);
            *normal = axe.x_direction;
            *binormal = axe.y_direction;
            *dtangent = DVec3::ZERO;
            *dnormal = DVec3::ZERO;
            *dbinormal = DVec3::ZERO;
            return true;
        }
        *binormal = tangent.cross(dc2).normalize();

        *normal = binormal.cross(*tangent);

        *dtangent = f_deriv(dc1, dc2);

        let instead_dc1 = tangent.cross(dc2);
        let instead_dc2 = dtangent.cross(dc2) + tangent.cross(dc3);
        *dbinormal = f_deriv(instead_dc1, instead_dc2);

        *dnormal = dbinormal.cross(*tangent) + binormal.cross(*dtangent);
        true
    }

    /// OCCT D2 (L684-758).
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
        let mut index = 0usize;
        let mut delta = 0.0f64;
        if self.is_singular(param, &mut index)
            && self.singular_d2(
                param,
                index,
                tangent,
                dtangent,
                d2tangent,
                normal,
                dnormal,
                d2normal,
                binormal,
                dbinormal,
                d2binormal,
                &mut delta,
            )
        {
            return true;
        }

        let the_param = param + delta;
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let (_p, dc1, dc2, dc3) = law_d3(trimmed, the_param);
        let dc4 = law_dn(trimmed, the_param, 4);
        *tangent = dc1.normalize();

        if tangent.cross(dc2).length() <= GP_RESOLUTION {
            let axe = Ax2::from_direction(DVec3::ZERO, *tangent);
            *normal = axe.x_direction;
            *binormal = axe.y_direction;
            *dtangent = DVec3::ZERO;
            *dnormal = DVec3::ZERO;
            *dbinormal = DVec3::ZERO;
            *d2tangent = DVec3::ZERO;
            *d2normal = DVec3::ZERO;
            *d2binormal = DVec3::ZERO;
            return true;
        }
        *binormal = tangent.cross(dc2).normalize();

        *normal = binormal.cross(*tangent);

        *dtangent = f_deriv(dc1, dc2);
        *d2tangent = d_deriv(dc1, dc2, dc3);

        let instead_dc1 = tangent.cross(dc2);
        let instead_dc2 = dtangent.cross(dc2) + tangent.cross(dc3);
        let instead_dc3 =
            d2tangent.cross(dc2) + 2.0 * dtangent.cross(dc3) + tangent.cross(dc4);
        *dbinormal = f_deriv(instead_dc1, instead_dc2);
        *d2binormal = d_deriv(instead_dc1, instead_dc2, instead_dc3);

        *dnormal = dbinormal.cross(*tangent) + binormal.cross(*dtangent);

        *d2normal = d2binormal.cross(*tangent) + 2.0 * dbinormal.cross(*dtangent)
            + binormal.cross(*d2tangent);

        true
    }

    /// OCCT NbIntervals (L760-800).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let tmp_s = match s {
            GeomAbsShape::C0 => GeomAbsShape::C2,
            GeomAbsShape::C1 => GeomAbsShape::C3,
            GeomAbsShape::C2 | GeomAbsShape::C3 | GeomAbsShape::CN => GeomAbsShape::CN,
        };
        let curve = self.base.my_curve.as_ref().unwrap();
        let nb_trimmed = curve_nb_intervals(curve, tmp_s);

        if !self.is_sngl {
            return nb_trimmed;
        }

        let trim_int = curve_intervals(curve, tmp_s);
        let mut fusion: Vec<f64> = Vec::new();
        fuse_intervals(&trim_int, self.my_sngl.as_ref().unwrap(), &mut fusion, 1e-12, true);

        fusion.len() - 1
    }

    /// OCCT Intervals (L802-845).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        let tmp_s = match s {
            GeomAbsShape::C0 => GeomAbsShape::C2,
            GeomAbsShape::C1 => GeomAbsShape::C3,
            GeomAbsShape::C2 | GeomAbsShape::C3 | GeomAbsShape::CN => GeomAbsShape::CN,
        };
        let curve = self.base.my_curve.as_ref().unwrap();

        if !self.is_sngl {
            *t = curve_intervals(curve, tmp_s);
            return;
        }

        let nb_trimmed = curve_nb_intervals(curve, tmp_s);
        let mut trim_int = vec![0.0f64; nb_trimmed + 1];
        trim_int.copy_from_slice(&curve_intervals(curve, tmp_s)[..nb_trimmed + 1]);

        let mut fusion: Vec<f64> = Vec::new();
        fuse_intervals(&trim_int, self.my_sngl.as_ref().unwrap(), &mut fusion, 1e-12, true);

        for i in 1..=fusion.len() {
            t[i - 1] = fusion[i - 1];
        }
    }

    /// OCCT GetAverageLaw (L847-879).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        let num = 20usize; // order of digitalization
        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut bn = DVec3::ZERO;
        *atangent = DVec3::ZERO;
        *anormal = DVec3::ZERO;
        *abinormal = DVec3::ZERO;
        let trimmed = self.base.my_trimmed.as_ref().unwrap();
        let first = curve_first_parameter(trimmed);
        let last = curve_last_parameter(trimmed);
        let step = (last - first) / num as f64;
        for i in 0..=num {
            let mut param = first + i as f64 * step;
            if param > last {
                param = last;
            }
            self.d0(param, &mut t, &mut n, &mut bn);
            *atangent += t;
            *anormal += n;
            *abinormal += bn;
        }
        *atangent /= (num + 1) as f64;
        *anormal /= (num + 1) as f64;

        *atangent = atangent.normalize();
        *abinormal = atangent.cross(*anormal).normalize();
        *anormal = abinormal.cross(*atangent);
    }

    /// OCCT IsConstant (L881-884).
    fn is_constant(&self) -> bool {
        matches!(self.base.my_curve, Some(Curve3::Line(_)))
    }

    /// OCCT IsOnlyBy3dCurve (L886-889).
    fn is_only_by3d_curve(&self) -> bool {
        true
    }
}

