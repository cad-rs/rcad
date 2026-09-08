//! OCCT GeomFill_PolynomialConvertor (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_PolynomialConvertor.hxx (members) + GeomFill_PolynomialConvertor.cxx
//! (whole file L29-365).
//!
//! GAP carrier: `plib_hermite_coefficients` (PLib::HermiteCoefficients,
//! TKMath/PLib/PLib.cxx L1404-1470) is not yet translated in rcad-kernel; it
//! is carried here pending relocation to the kernel PLib port.

use glam::DVec3;
use rcad_kernel::math::convert_comp_polynomial_to_poles::ConvertCompPolynomialToPoles;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::VecD;

use super::gp_mat::GpMat;

/// OCCT GeomFill_PolynomialConvertor (GeomFill_PolynomialConvertor.hxx
/// L59-93): Ordre = 8, BH(1, Ordre, 1, Ordre), myinit flag.
pub struct PolynomialConvertor {
    /// OCCT GeomFill_PolynomialConvertor::Ordre = 8.
    ordre: i32,
    /// OCCT myinit.
    myinit: bool,
    /// OCCT BH (math_Matrix(1, Ordre, 1, Ordre)).
    bh: Matrix,
}

impl PolynomialConvertor {
    /// OCCT GeomFill_PolynomialConvertor::GeomFill_PolynomialConvertor
    /// (L29-34): Ordre(8), myinit(false), BH(1, Ordre, 1, Ordre).
    pub fn new() -> Self {
        let ordre = 8;
        PolynomialConvertor {
            ordre,
            myinit: false,
            bh: Matrix::new(1, ordre, 1, ordre),
        }
    }

    /// OCCT Initialized (L36-39).
    pub fn initialized(&self) -> bool {
        self.myinit
    }

    /// OCCT Init (L41-97): BH = B * H where B is the monomial-to-BSpline
    /// conversion matrix on [-1,1], degree 7, and H is the Hermite
    /// coefficients matrix. Both are mathematical constants computed once.
    pub fn init(&mut self) {
        if self.myinit {
            return;
        }

        // OCCT: static const math_Matrix THE_BH_MATRIX = []() { ... }();
        // (computed once; the static-lambda form maps to a plain block).
        let the_bh_matrix: Matrix = {
            let an_ordre = 8usize;
            // NCollection_Array1<double> aCoeffs(1, anOrdre * anOrdre);
            let mut a_coeffs = vec![0.0f64; an_ordre * an_ordre];
            let an_inter = [-1.0f64, 1.0f64];
            let a_true_inter = [-1.0f64, 1.0f64];
            for ii in 1..=an_ordre {
                a_coeffs[ii - 1 + (ii - 1) * an_ordre] = 1.0;
            }

            // OCCT: Convert_CompPolynomialToPoles aConverter(anOrdre,
            // anOrdre - 1, anOrdre - 1, aCoeffs, anInter, aTrueInter)
            // — the "only one span" constructor: Dimension = 8,
            // MaxDegree = 7, Degree = 7; it delegates to the multi-curve
            // constructor with NumCurves = 1 and NumCoeffPerCurve(1) =
            // Degree + 1 (Convert_CompPolynomialToPoles.cxx L140-165).
            let a_converter = ConvertCompPolynomialToPoles::from_arrays(
                1,
                an_ordre,
                an_ordre - 1,
                &[],
                &[(an_ordre as i32 - 1) + 1],
                &a_coeffs,
                &an_inter,
                &a_true_inter,
            );
            // OCCT: const NCollection_Array2<double>& aPoles =
            // aConverter.Poles(); — (1..NbPoles, 1..Dimension), row-major
            // flat storage in rcad.
            let a_poles = a_converter.poles();
            let mut a_b = Matrix::new(1, an_ordre as i32, 1, an_ordre as i32);
            for jj in 1..=an_ordre as i32 {
                for ii in 1..=an_ordre as i32 {
                    let mut a_term = a_poles[((ii - 1) as usize) * an_ordre + ((jj - 1) as usize)];
                    if (a_term - 1.0).abs() < 1.0e-9 {
                        a_term = 1.0;
                    }
                    if (a_term + 1.0).abs() < 1.0e-9 {
                        a_term = -1.0;
                    }
                    a_b.set(ii, jj, a_term);
                }
            }

            // OCCT: math_Matrix aH(1, anOrdre, 1, anOrdre);
            // PLib::HermiteCoefficients(-1, 1, anOrdre/2 - 1, anOrdre/2 - 1, aH);
            let mut a_h = Matrix::new(1, an_ordre as i32, 1, an_ordre as i32);
            plib_hermite_coefficients(
                -1.0,
                1.0,
                an_ordre as i32 / 2 - 1,
                an_ordre as i32 / 2 - 1,
                &mut a_h,
            );
            // OCCT: aH.Transpose(); return math_Matrix(aB * aH).
            let a_h = a_h.transposed();
            a_b.multiplied(&a_h)
        };

        self.bh = the_bh_matrix;
        self.myinit = true;
    }

    /// OCCT Section (L99-148) — poles of the degree-7 polynomial arc.
    pub fn section(
        &self,
        first_pnt: DVec3,
        center: DVec3,
        dir: DVec3,
        angle: f64,
        poles: &mut [DVec3],
    ) {
        let ordre = self.ordre;
        let mut vx = Vector::new(1, ordre);
        let mut vy = Vector::new(1, ordre);
        let cos_b = angle.cos();
        let sin_b = angle.sin();
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let v1 = first_pnt - center;
        // OCCT: V2 = Dir ^ V1.
        let v2 = dir.cross(v1);
        let beta = angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;

        // Calcul de la transformation — gp_Mat M(V1.X(), V2.X(), 0, ...):
        // columns are V1, V2, 0.
        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);

        // Calcul des contraintes -----------
        vx.set(1, 1.0);
        vy.set(1, 0.0);
        vx.set(2, 0.0);
        vy.set(2, beta);
        vx.set(3, -beta2);
        vy.set(3, 0.0);
        vx.set(4, 0.0);
        vy.set(4, -beta3);
        vx.set(5, cos_b);
        vy.set(5, sin_b);
        vx.set(6, -beta * sin_b);
        vy.set(6, beta * cos_b);
        vx.set(7, -beta2 * cos_b);
        vy.set(7, -beta2 * sin_b);
        vx.set(8, beta3 * sin_b);
        vy.set(8, -beta3 * cos_b);

        // Calcul des poles — Px = BH * Vx; Py = BH * Vy.
        let px = self.bh.multiplied_vec(&vx);
        let py = self.bh.multiplied_vec(&vy);
        for ii in 1..=ordre {
            let mut pnt = DVec3::new(px.get(ii), py.get(ii), 0.0);
            // OCCT: pnt *= M (gp_XYZ row-vector multiply).
            pnt = GpMat::multiply_xyz_row(pnt, &m);
            pnt += center;
            poles[(ii - 1) as usize] = pnt;
        }
    }

    /// OCCT Section (L150-232) — poles + first derivatives.
    #[allow(clippy::too_many_arguments)]
    pub fn section_d1(
        &self,
        first_pnt: DVec3,
        d_first_pnt: DVec3,
        center: DVec3,
        d_center: DVec3,
        dir: DVec3,
        d_dir: DVec3,
        angle: f64,
        d_angle: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
    ) {
        let ordre = self.ordre;
        let mut vx = Vector::new(1, ordre);
        let mut vy = Vector::new(1, ordre);
        let mut dvx = Vector::new(1, ordre);
        let mut dvy = Vector::new(1, ordre);
        let cos_b = angle.cos();
        let sin_b = angle.sin();
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let v1 = first_pnt - center;
        let mut v2 = dir.cross(v1);
        let beta = angle / 2.0;
        let bprim = d_angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;

        // Calcul des transformations
        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);
        let v1_prim = d_first_pnt - d_center;
        v2 = d_dir.cross(v1) + dir.cross(v1_prim);
        let m_prim = GpMat::from_rows(
            v1_prim.x, v2.x, 0.0, v1_prim.y, v2.y, 0.0, v1_prim.z, v2.z, 0.0,
        );

        // Calcul des contraintes -----------
        vx.set(1, 1.0);
        vy.set(1, 0.0);
        vx.set(2, 0.0);
        vy.set(2, beta);
        vx.set(3, -beta2);
        vy.set(3, 0.0);
        vx.set(4, 0.0);
        vy.set(4, -beta3);
        vx.set(5, cos_b);
        vy.set(5, sin_b);
        vx.set(6, -beta * sin_b);
        vy.set(6, beta * cos_b);
        vx.set(7, -beta2 * cos_b);
        vy.set(7, -beta2 * sin_b);
        vx.set(8, beta3 * sin_b);
        vy.set(8, -beta3 * cos_b);

        let b_bprim = bprim * beta;
        let b2_bprim = bprim * beta2;
        dvx.set(1, 0.0);
        dvy.set(1, 0.0);
        dvx.set(2, 0.0);
        dvy.set(2, bprim);
        dvx.set(3, -2.0 * b_bprim);
        dvy.set(3, 0.0);
        dvx.set(4, 0.0);
        dvy.set(4, -3.0 * b2_bprim);
        dvx.set(5, -2.0 * bprim * sin_b);
        dvy.set(5, 2.0 * bprim * cos_b);
        dvx.set(6, -bprim * sin_b - 2.0 * b_bprim * cos_b);
        dvy.set(6, bprim * cos_b - 2.0 * b_bprim * sin_b);
        dvx.set(7, 2.0 * b_bprim * (-cos_b + beta * sin_b));
        dvy.set(7, -2.0 * b_bprim * (sin_b + beta * cos_b));
        dvx.set(8, b2_bprim * (3.0 * sin_b + 2.0 * beta * cos_b));
        dvy.set(8, b2_bprim * (2.0 * beta * sin_b - 3.0 * cos_b));

        // Calcul des poles
        let px = self.bh.multiplied_vec(&vx);
        let py = self.bh.multiplied_vec(&vy);
        let dpx = self.bh.multiplied_vec(&dvx);
        let dpy = self.bh.multiplied_vec(&dvy);

        for ii in 1..=ordre {
            let mut p = DVec3::new(px.get(ii), py.get(ii), 0.0);
            poles[(ii - 1) as usize] = m.multiplied_xyz(p) + center;
            // OCCT: P *= MPrim (gp_XYZ row-vector multiply).
            p = GpMat::multiply_xyz_row(p, &m_prim);
            let mut dp = DVec3::new(dpx.get(ii), dpy.get(ii), 0.0);
            // OCCT: DP *= M.
            dp = GpMat::multiply_xyz_row(dp, &m);
            // OCCT: aux.SetLinearForm(1, P, 1, DP, DCenter.XYZ()).
            let aux = p + dp + d_center;
            dpoles[(ii - 1) as usize] = aux;
        }
    }

    /// OCCT Section (L234-365) — poles + first and second derivatives.
    #[allow(clippy::too_many_arguments)]
    pub fn section_d2(
        &self,
        first_pnt: DVec3,
        d_first_pnt: DVec3,
        d2_first_pnt: DVec3,
        center: DVec3,
        d_center: DVec3,
        d2_center: DVec3,
        dir: DVec3,
        d_dir: DVec3,
        d2_dir: DVec3,
        angle: f64,
        d_angle: f64,
        d2_angle: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
    ) {
        let ordre = self.ordre;
        let mut vx = Vector::new(1, ordre);
        let mut vy = Vector::new(1, ordre);
        let mut dvx = Vector::new(1, ordre);
        let mut dvy = Vector::new(1, ordre);
        let mut d2vx = Vector::new(1, ordre);
        let mut d2vy = Vector::new(1, ordre);

        let cos_b = angle.cos();
        let sin_b = angle.sin();
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let v1 = first_pnt - center;
        let mut v2 = dir.cross(v1);
        let beta = angle / 2.0;
        let bprim = d_angle / 2.0;
        let mut bsecn = d2_angle / 2.0;
        // OCCT L263 repeats the assignment: bsecn = D2Angle / 2;
        bsecn = d2_angle / 2.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let bprim2 = bprim * bprim;

        // Calcul des transformations
        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);
        let v1_prim = d_first_pnt - d_center;
        v2 = d_dir.cross(v1) + dir.cross(v1_prim);
        let m_prim = GpMat::from_rows(
            v1_prim.x, v2.x, 0.0, v1_prim.y, v2.y, 0.0, v1_prim.z, v2.z, 0.0,
        );
        let v1_secn = d2_first_pnt - d2_center;
        let mut v2 = d_dir.cross(v1_prim);
        v2 *= 2.0;
        v2 += d2_dir.cross(v1) + dir.cross(v1_secn);
        let m_secn = GpMat::from_rows(
            v1_secn.x, v2.x, 0.0, v1_secn.y, v2.y, 0.0, v1_secn.z, v2.z, 0.0,
        );

        // Calcul des contraintes -----------
        vx.set(1, 1.0);
        vy.set(1, 0.0);
        vx.set(2, 0.0);
        vy.set(2, beta);
        vx.set(3, -beta2);
        vy.set(3, 0.0);
        vx.set(4, 0.0);
        vy.set(4, -beta3);
        vx.set(5, cos_b);
        vy.set(5, sin_b);
        vx.set(6, -beta * sin_b);
        vy.set(6, beta * cos_b);
        vx.set(7, -beta2 * cos_b);
        vy.set(7, -beta2 * sin_b);
        vx.set(8, beta3 * sin_b);
        vy.set(8, -beta3 * cos_b);

        let b_bprim = bprim * beta;
        let b2_bprim = bprim * beta2;
        let b_bsecn = bsecn * beta;
        dvx.set(1, 0.0);
        dvy.set(1, 0.0);
        dvx.set(2, 0.0);
        dvy.set(2, bprim);
        dvx.set(3, -2.0 * b_bprim);
        dvy.set(3, 0.0);
        dvx.set(4, 0.0);
        dvy.set(4, -3.0 * b2_bprim);
        dvx.set(5, -2.0 * bprim * sin_b);
        dvy.set(5, 2.0 * bprim * cos_b);
        dvx.set(6, -bprim * sin_b - 2.0 * b_bprim * cos_b);
        dvy.set(6, bprim * cos_b - 2.0 * b_bprim * sin_b);
        dvx.set(7, 2.0 * b_bprim * (-cos_b + beta * sin_b));
        dvy.set(7, -2.0 * b_bprim * (sin_b + beta * cos_b));
        dvx.set(8, b2_bprim * (3.0 * sin_b + 2.0 * beta * cos_b));
        dvy.set(8, b2_bprim * (2.0 * beta * sin_b - 3.0 * cos_b));

        d2vx.set(1, 0.0);
        d2vy.set(1, 0.0);
        d2vx.set(2, 0.0);
        d2vy.set(2, bsecn);
        d2vx.set(3, -2.0 * (bprim2 + b_bsecn));
        d2vy.set(3, 0.0);
        d2vx.set(4, 0.0);
        d2vy.set(4, -3.0 * beta * (2.0 * bprim2 + b_bsecn));
        d2vx.set(5, -2.0 * (bsecn * sin_b + 2.0 * bprim2 * cos_b));
        d2vy.set(5, 2.0 * (bsecn * cos_b - 2.0 * bprim2 * sin_b));
        d2vx.set(
            6,
            (4.0 * beta * bprim2 - bsecn) * sin_b - 2.0 * (2.0 * bprim2 + b_bsecn) * cos_b,
        );
        d2vy.set(
            6,
            (bsecn - 4.0 * beta * bprim2) * cos_b - 2.0 * (b_bsecn + 2.0 * bprim2) * sin_b,
        );

        let mut aux = 2.0 * (bprim2 + b_bsecn);
        d2vx.set(
            7,
            aux * (-cos_b + beta * sin_b)
                + 2.0 * beta * bprim2 * (2.0 * beta * cos_b + 3.0 * sin_b),
        );
        d2vy.set(
            7,
            -aux * (sin_b + beta * cos_b) - 2.0 * beta * bprim2 * (3.0 * cos_b - 2.0 * beta * sin_b),
        );

        aux = beta * (2.0 * bprim2 + b_bsecn);
        d2vx.set(
            8,
            aux * (3.0 * sin_b + 2.0 * beta * cos_b)
                + 4.0 * beta2 * bprim2 * (2.0 * cos_b - beta * sin_b),
        );
        d2vy.set(
            8,
            aux * (2.0 * beta * sin_b - 3.0 * cos_b)
                + 4.0 * beta2 * bprim2 * (2.0 * sin_b + beta * cos_b),
        );

        // Calcul des poles
        let px = self.bh.multiplied_vec(&vx);
        let py = self.bh.multiplied_vec(&vy);
        let dpx = self.bh.multiplied_vec(&dvx);
        let dpy = self.bh.multiplied_vec(&dvy);
        let d2px = self.bh.multiplied_vec(&d2vx);
        let d2py = self.bh.multiplied_vec(&d2vy);

        for ii in 1..=ordre {
            let mut p = DVec3::new(px.get(ii), py.get(ii), 0.0);
            let mut dp = DVec3::new(dpx.get(ii), dpy.get(ii), 0.0);
            let mut d2p = DVec3::new(d2px.get(ii), d2py.get(ii), 0.0);

            poles[(ii - 1) as usize] = m.multiplied_xyz(p) + center;
            // OCCT: auxyz.SetLinearForm(1, MPrim * P, 1, M * DP, DCenter.XYZ()).
            let auxyz = m_prim.multiplied_xyz(p) + m.multiplied_xyz(dp) + d_center;
            dpoles[(ii - 1) as usize] = auxyz;
            // OCCT: P *= MSecn; DP *= MPrim; D2P *= M (row-vector form).
            p = GpMat::multiply_xyz_row(p, &m_secn);
            dp = GpMat::multiply_xyz_row(dp, &m_prim);
            d2p = GpMat::multiply_xyz_row(d2p, &m);
            // OCCT: auxyz.SetLinearForm(1, P, 2, DP, 1, D2P, D2Center.XYZ()).
            let auxyz = p + 2.0 * dp + d2p + d2_center;
            d2poles[(ii - 1) as usize] = auxyz;
        }
    }
}

impl Default for PolynomialConvertor {
    fn default() -> Self {
        Self::new()
    }
}

/// GAP carrier: OCCT PLib::HermiteCoefficients (TKMath/PLib/PLib.cxx
/// L1404-1470). Carried here (crate-private) because rcad-kernel's PLib port
/// does not yet host it; move to `rcad_kernel::math::plib` when available.
pub(crate) fn plib_hermite_coefficients(
    first_parameter: f64,
    last_parameter: f64,
    first_order: i32,
    last_order: i32,
    matrix_coefs: &mut Matrix,
) -> bool {
    let nb_coeff = (first_order + last_order + 2) as i32;
    let mut iof: i32 = 0;
    let mut prod: f64;
    let mut t_borne = first_parameter;
    let mut coeff = Vector::new(1, nb_coeff);
    let mut b = Vector::new_init(1, nb_coeff, 0.0);
    let mut mat = Matrix::new_init(1, nb_coeff, 1, nb_coeff, 0.0);

    // Test de validites
    if first_order < 0 || last_order < 0 {
        return false;
    }
    let d1 = first_parameter.abs();
    let d2 = last_parameter.abs();
    if d1 > 100.0 || d2 > 100.0 {
        return false;
    }
    let d2 = d2 + d1;
    if d2 < 0.01 {
        return false;
    }
    if (last_parameter - first_parameter).abs() / d2 < 0.01 {
        return false;
    }

    // Calcul de la matrice a inverser (MAT)
    let ordre = [first_order + 1, last_order + 1];

    for cote in 0..=1i32 {
        for r in 1..=nb_coeff {
            coeff.set(r, 1.0);
        }

        for pp in 1..=ordre[cote as usize] {
            let ii = pp + iof;
            prod = 1.0;

            for jj in pp..=nb_coeff {
                // tout se passe dans les 3 lignes suivantes
                let coeff_jj = coeff.get(jj);
                mat.set(ii, jj, coeff_jj * prod);
                coeff.set(jj, coeff_jj * (jj - pp) as f64);
                prod *= t_borne;
            }
        }
        t_borne = last_parameter;
        iof = ordre[0];
    }

    // resolution du systemes
    let resol_coeff = MathGauss::with_min_pivot(&mat.data, 1.0e-10);
    if !resol_coeff.is_done() {
        return false;
    }

    for ii in 1..=nb_coeff {
        // OCCT: B(ii) = 1; ResolCoeff.Solve(B, Coeff) — the rcad MathGauss
        // solves in place, so B is copied into a scratch vector.
        b.set(ii, 1.0);
        let mut x = VecD::new(nb_coeff as usize);
        for r in 1..=nb_coeff {
            x.set(r as usize, b.get(r));
        }
        let mut x = x;
        resol_coeff.solve(&mut x);
        for r in 1..=nb_coeff {
            coeff.set(r, x.get(r as usize));
        }
        // OCCT: MatrixCoefs.SetRow(ii, Coeff).
        for c in 1..=nb_coeff {
            matrix_coefs.set(ii, c, coeff.get(c));
        }
        b.set(ii, 0.0);
    }
    true
}
