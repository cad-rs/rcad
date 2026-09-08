//! OCCT GeomFill_QuasiAngularConvertor (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_QuasiAngularConvertor.hxx (members) + GeomFill_QuasiAngularConvertor.cxx
//! (whole file L54-591, plus the NullAngle constant at L27).

use glam::DVec3;
use rcad_kernel::math::math_matrix::{Matrix, Vector};

use super::gp_mat::GpMat;

/// OCCT L27: #define NullAngle 1.e-6.
const NULL_ANGLE: f64 = 1.0e-6;

/// OCCT GeomFill_QuasiAngularConvertor (GeomFill_QuasiAngularConvertor.hxx
/// L58-94): myinit, B(1,7,1,7), Px/Py/W/Vx/Vy/Vw(1,7).
///
/// QuasiAngular is rational definition of cos(theta(t)) and sin(theta(t)) on
/// [-alpha, +alpha] with U(t) = 1 + b t^2, V(t) = t + c t^3, c = 1/3 + b and
/// b = -1/gamma^2 + gamma / (3 (tang gamma - gamma)), gamma = alpha / 2.
pub struct QuasiAngularConvertor {
    /// OCCT myinit.
    myinit: bool,
    /// OCCT B (math_Matrix(1, 7, 1, 7)).
    b: Matrix,
    /// OCCT Px (math_Vector(1, 7)).
    px: Vector,
    /// OCCT Py (math_Vector(1, 7)).
    py: Vector,
    /// OCCT W (math_Vector(1, 7)).
    w: Vector,
    /// OCCT Vx (math_Vector(1, 7)).
    vx: Vector,
    /// OCCT Vy (math_Vector(1, 7)).
    vy: Vector,
    /// OCCT Vw (math_Vector(1, 7)).
    vw: Vector,
}

impl QuasiAngularConvertor {
    /// OCCT GeomFill_QuasiAngularConvertor::GeomFill_QuasiAngularConvertor
    /// (L54-64): myinit(false), B(1,7,1,7), Px(1,7), Py(1,7), W(1,7),
    /// Vx(1,7), Vy(1,7), Vw(1,7).
    pub fn new() -> Self {
        QuasiAngularConvertor {
            myinit: false,
            b: Matrix::new(1, 7, 1, 7),
            px: Vector::new(1, 7),
            py: Vector::new(1, 7),
            w: Vector::new(1, 7),
            vx: Vector::new(1, 7),
            vy: Vector::new(1, 7),
            vw: Vector::new(1, 7),
        }
    }

    /// OCCT Initialized (L66-69).
    pub fn initialized(&self) -> bool {
        self.myinit
    }

    /// OCCT Init (L71-130): B is the monomial-to-BSpline conversion matrix
    /// on [-1,1], degree 6. It is a mathematical constant computed once.
    pub fn init(&mut self) {
        if self.myinit {
            return;
        }

        // OCCT: static const math_Matrix THE_BASIS = []() { ... }();
        let the_basis: Matrix = {
            let an_ordre = 7usize;
            let mut a_coeffs = vec![0.0f64; an_ordre * an_ordre];
            let an_inter = [-1.0f64, 1.0f64];
            let a_true_inter = [-1.0f64, 1.0f64];
            for ii in 1..=an_ordre {
                a_coeffs[ii - 1 + (ii - 1) * an_ordre] = 1.0;
            }

            // OCCT: Convert_CompPolynomialToPoles aConverter(anOrdre,
            // anOrdre - 1, anOrdre - 1, aCoeffs, anInter, aTrueInter) — the
            // "only one span" constructor (Dimension = 7, MaxDegree = 6,
            // Degree = 6), see Convert_CompPolynomialToPoles.cxx L140-165.
            let a_converter =
                rcad_kernel::math::convert_comp_polynomial_to_poles::ConvertCompPolynomialToPoles::from_arrays(
                    1,
                    an_ordre,
                    an_ordre - 1,
                    &[],
                    &[(an_ordre as i32 - 1) + 1],
                    &a_coeffs,
                    &an_inter,
                    &a_true_inter,
                );
            let a_poles = a_converter.poles();
            let mut a_result = Matrix::new(1, an_ordre as i32, 1, an_ordre as i32);
            for jj in 1..=an_ordre as i32 {
                for ii in 1..=an_ordre as i32 {
                    let mut a_term = a_poles[((ii - 1) as usize) * an_ordre + ((jj - 1) as usize)];
                    if (a_term - 1.0).abs() < 1.0e-9 {
                        a_term = 1.0;
                    }
                    if (a_term + 1.0).abs() < 1.0e-9 {
                        a_term = -1.0;
                    }
                    a_result.set(ii, jj, a_term);
                }
            }
            a_result
        };

        self.b = the_basis;

        // OCCT: Vx.Init(0); Vx(1) = 1; Vy.Init(0); Vy(2) = 2; Vw.Init(0);
        // Vw(1) = 1; myinit = true;
        for r in 1..=7 {
            self.vx.set(r, 0.0);
            self.vy.set(r, 0.0);
            self.vw.set(r, 0.0);
        }
        self.vx.set(1, 1.0);
        self.vy.set(2, 2.0);
        self.vw.set(1, 1.0);
        self.myinit = true;
    }

    /// OCCT Section (L132-216) — poles + weights of the QuasiAngular arc.
    pub fn section(
        &mut self,
        first_pnt: DVec3,
        center: DVec3,
        dir: DVec3,
        angle: f64,
        poles: &mut [DVec3],
        weights: &mut [f64],
    ) {
        // Calcul de la transformation
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let mut v1 = first_pnt - center;
        let mut rot = GpMat::identity();
        rot.set_rotation(dir, angle / 2.0);
        // OCCT: aux = V1.XYZ(); aux *= Rot; V1.SetXYZ(aux).
        v1 = GpMat::multiply_xyz_row(v1, &rot);
        // OCCT: V2 = Dir ^ V1.
        let v2 = dir.cross(v1);

        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);

        // Calcul des coeffs -----------
        let beta = angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;

        // OCCT: double b, tan_b (declared at L139; tan_b is only assigned in
        // the general branch, so it stays branch-local in Rust).
        let b: f64;
        if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
            if beta.abs() < NULL_ANGLE {
                let cf = 2.0 / (3.0 * 5.0 * 7.0);
                b = -(0.2 + cf * beta2) / (1.0 + 0.2 * beta2);
                // b = beta5 / cf;
            } else {
                let tan_b = beta.tan();
                let mut b_val = -1.0 / beta2;
                b_val += beta / (3.0 * (tan_b - beta));
                b = b_val;
            }
        } else {
            b = -1.0 / beta2;
        }
        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;

        // X = U*U - V*V
        self.vx.set(3, beta2 * (2.0 * b - 1.0));
        self.vx.set(5, beta4 * (b2 - 2.0 * c));
        self.vx.set(7, -beta6 * c2);

        // Y = 2*U*V
        self.vy.set(2, 2.0 * beta);
        self.vy.set(4, beta3 * 2.0 * (c + b));
        self.vy.set(6, 2.0 * beta5 * b * c);

        // W = U*U + V*V
        self.vw.set(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set(5, beta4 * (2.0 * c + b2));
        self.vw.set(7, beta6 * c2);

        // Calculs des poles — Px.Multiply(B, Vx); Py.Multiply(B, Vy);
        // W.Multiply(B, Vw).
        self.px = self.b.multiplied_vec(&self.vx);
        self.py = self.b.multiplied_vec(&self.vy);
        self.w = self.b.multiplied_vec(&self.vw);

        // Transfo
        for ii in 1..=7 {
            let wi = self.w.get(ii);
            let mut pnt = DVec3::new(self.px.get(ii) / wi, self.py.get(ii) / wi, 0.0);
            // OCCT: pnt *= M (gp_XYZ row-vector multiply).
            pnt = GpMat::multiply_xyz_row(pnt, &m);
            pnt += center;
            poles[(ii - 1) as usize] = pnt;
            weights[(ii - 1) as usize] = wi;
        }
    }

    /// OCCT Section (L218-365) — poles + weights, first derivatives.
    #[allow(clippy::too_many_arguments)]
    pub fn section_d1(
        &mut self,
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
        weights: &mut [f64],
        dweights: &mut [f64],
    ) {
        let ordre = 7;
        // OCCT L232: math_Vector DVx(1, Ordre), DVy(1, Ordre), DVw(1, Ordre),
        // DPx(1, Ordre), DPy(1, Ordre), DW(1, Ordre) — the DV* vectors are
        // created at first use below (Rust has no uninitialized locals).
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let mut v1 = first_pnt - center;

        // Calcul des transformations
        // La rotation s'ecrit I + sin(Ang) * D + (1. - cos(Ang)) * D*D
        // ou D est l'application x -> Dir ^ x
        let mut rot = GpMat::identity();
        rot.set_rotation(dir, angle / 2.0);
        // La derive s'ecrit donc :
        // AngPrim * (sin(Ang)*D*D + cos(Ang)*D)
        // + sin(Ang)*DPrim  + (1. - cos(Ang)) *(DPrim*D + D*DPrim)
        let sina = (angle / 2.0).sin();
        let cosa = (angle / 2.0).cos();
        let mut d = GpMat::identity();
        d.set_cross(dir);
        let mut d_prim = GpMat::identity();
        d_prim.set_cross(d_dir);

        let mut rot_prim = d.powered(2).multiplied_scalar(sina);
        rot_prim = rot_prim.added(&d.multiplied_scalar(cosa));
        rot_prim = rot_prim.multiplied_scalar(d_angle / 2.0);
        rot_prim = rot_prim.added(&d_prim.multiplied_scalar(sina));
        rot_prim = rot_prim.added(
            &d_prim
                .multiplied_mat(&d)
                .added(&d.multiplied_mat(&d_prim))
                .multiplied_scalar(1.0 - cosa),
        );

        let mut v1_prim = GpMat::multiply_xyz_row(d_first_pnt - d_center, &rot);
        v1_prim += GpMat::multiply_xyz_row(v1, &rot_prim);
        v1 = GpMat::multiply_xyz_row(v1, &rot);
        let v2 = dir.cross(v1);
        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);
        let v2 = d_dir.cross(v1) + dir.cross(v1_prim);
        let m_prim = GpMat::from_rows(
            v1_prim.x, v2.x, 0.0, v1_prim.y, v2.y, 0.0, v1_prim.z, v2.z, 0.0,
        );

        // Calcul des constante -----------
        let beta = angle / 4.0;
        let betaprim = d_angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;

        // OCCT L234-236: double b, tan_b, bpr, dtan_b (branch-local assigns
        // stay branch-local in Rust).
        let b: f64;
        let bpr: f64;
        if beta.abs() < NULL_ANGLE {
            // On calcul b par D.L
            let cf = 2.0 / (3.0 * 5.0 * 7.0);
            let num = 0.2 + cf * beta2;
            let denom = 1.0 + 0.2 * beta2;
            b = -num / denom;
            bpr = -2.0 * beta * betaprim * (cf * denom - 0.2 * num) / (denom * denom);
        } else {
            let mut b_val = -1.0 / beta2;
            let mut bpr_val = (2.0 * betaprim) / beta3;
            if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
                let tan_b = beta.tan();
                let dtan_b = betaprim * (1.0 + tan_b * tan_b);
                let b2 = tan_b - beta;
                b_val += beta / (3.0 * b2);
                bpr_val += (betaprim * tan_b - beta * dtan_b) / (3.0 * b2 * b2);
            }
            b = b_val;
            bpr = bpr_val;
        }

        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;

        // X = U*U - V*V
        self.vx.set(3, beta2 * (2.0 * b - 1.0));
        self.vx.set(5, beta4 * (b2 - 2.0 * c));
        self.vx.set(7, -beta6 * c2);
        // OCCT: DVx.Init(0).
        let mut dvx = Self::zero_vector();
        dvx.set(3, 2.0 * (beta * betaprim * (2.0 * b - 1.0) + bpr * beta2));
        dvx.set(5, 4.0 * beta3 * betaprim * (b2 - 2.0 * c) + 2.0 * beta4 * bpr * (b - 1.0));
        dvx.set(7, -6.0 * beta5 * betaprim * c2 - 2.0 * beta6 * bpr * c);

        // Y = 2*U*V
        self.vy.set(2, 2.0 * beta);
        self.vy.set(4, beta3 * 2.0 * (c + b));
        self.vy.set(6, 2.0 * beta5 * b * c);
        // OCCT: DVy.Init(0).
        let mut dvy = Self::zero_vector();
        dvy.set(2, 2.0 * betaprim);
        dvy.set(4, 6.0 * beta2 * betaprim * (b + c) + 4.0 * beta3 * bpr);
        dvy.set(6, 10.0 * beta4 * betaprim * b * c + 2.0 * beta5 * bpr * (b + c));

        // W = U*U + V*V
        self.vw.set(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set(5, beta4 * (2.0 * c + b2));
        self.vw.set(7, beta6 * c2);
        // OCCT: DVw.Init(0) — with the commented-out variants at L332/L334/L336
        // replaced by the active expressions.
        let mut dvw = Self::zero_vector();
        dvw.set(3, 2.0 * beta * (betaprim * (1.0 + 2.0 * b) + beta * bpr));
        dvw.set(5, 2.0 * beta3 * (2.0 * betaprim * (2.0 * c + b2) + beta * bpr * (b + 1.0)));
        dvw.set(7, 2.0 * beta5 * c * (3.0 * betaprim * c + beta * bpr));

        // Calcul des poles
        self.px = self.b.multiplied_vec(&self.vx);
        self.py = self.b.multiplied_vec(&self.vy);
        self.w = self.b.multiplied_vec(&self.vw);
        let dpx = self.b.multiplied_vec(&dvx);
        let dpy = self.b.multiplied_vec(&dvy);
        let dw = self.b.multiplied_vec(&dvw);

        for ii in 1..=ordre {
            let wi = self.w.get(ii);
            let mut p = DVec3::new(self.px.get(ii) / wi, self.py.get(ii) / wi, 0.0);
            let mut dp = DVec3::new(dpx.get(ii) / wi, dpy.get(ii) / wi, 0.0);
            dp -= (dw.get(ii) / wi) * p;

            poles[(ii - 1) as usize] = m.multiplied_xyz(p) + center;
            // OCCT: P *= MPrim; DP *= M.
            p = GpMat::multiply_xyz_row(p, &m_prim);
            dp = GpMat::multiply_xyz_row(dp, &m);
            // OCCT: aux.SetLinearForm(1, P, 1, DP, DCenter.XYZ()).
            let aux = p + dp + d_center;
            dpoles[(ii - 1) as usize] = aux;
            weights[(ii - 1) as usize] = wi;
            dweights[(ii - 1) as usize] = dw.get(ii);
        }
    }

    /// OCCT Section (L367-591) — poles + weights, first and second derivatives.
    #[allow(clippy::too_many_arguments)]
    pub fn section_d2(
        &mut self,
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
        weights: &mut [f64],
        dweights: &mut [f64],
        d2weights: &mut [f64],
    ) {
        let ordre = 7;
        // OCCT gp_Vec V1(Center, FirstPnt) = FirstPnt - Center.
        let v1 = first_pnt - center;

        // Calcul des transformations
        // La rotation s'ecrit I + sin(Ang) * D + (1. - cos(Ang)) * D*D
        // ou D est l'application x -> Dir ^ x
        let mut rot = GpMat::identity();
        rot.set_rotation(dir, angle / 2.0);
        let sina = (angle / 2.0).sin();
        let cosa = (angle / 2.0).cos();
        let mut d = GpMat::identity();
        d.set_cross(dir);
        let mut d_prim = GpMat::identity();
        d_prim.set_cross(d_dir);
        let mut d_secn = GpMat::identity();
        d_secn.set_cross(d2_dir);

        let ddp = d_prim.multiplied_mat(&d).added(&d.multiplied_mat(&d_prim));
        let mut rot_prim = d.powered(2).multiplied_scalar(sina);
        rot_prim = rot_prim.added(&d.multiplied_scalar(cosa));
        rot_prim = rot_prim.multiplied_scalar(d_angle / 2.0);
        rot_prim = rot_prim.added(&d_prim.multiplied_scalar(sina));
        rot_prim = rot_prim.added(&ddp.multiplied_scalar(1.0 - cosa));

        let mut rot_secn = d.powered(2).multiplied_scalar(sina);
        rot_secn = rot_secn.added(&d.multiplied_scalar(cosa));
        rot_secn = rot_secn.multiplied_scalar(d2_angle / 2.0);
        let mut maux = d.powered(2).multiplied_scalar(cosa);
        maux = maux.added(&d.multiplied_scalar(-sina));
        maux = maux.multiplied_scalar(d_angle / 2.0);
        maux = maux.added(&ddp.multiplied_scalar(2.0 * sina));
        maux = maux.added(&d_prim.multiplied_scalar(2.0 * cosa));
        maux = maux.multiplied_scalar(d_angle / 2.0);
        rot_secn = rot_secn.added(&maux);
        maux = d_secn.multiplied_mat(&d).added(&d.multiplied_mat(&d_secn));
        maux = maux.added(&d_prim.powered(2).multiplied_scalar(2.0));
        maux = maux.multiplied_scalar(1.0 - cosa);
        maux = maux.added(&d_secn.multiplied_scalar(sina));
        rot_secn = rot_secn.added(&maux);

        let mut v1_prim = d_first_pnt - d_center;
        let mut auxyz = GpMat::multiply_xyz_row(d2_first_pnt - d2_center, &rot);
        auxyz += 2.0 * GpMat::multiply_xyz_row(v1_prim, &rot_prim);
        auxyz += GpMat::multiply_xyz_row(v1, &rot_secn);
        let v1_secn = auxyz;
        let auxyz = GpMat::multiply_xyz_row(v1_prim, &rot) + GpMat::multiply_xyz_row(v1, &rot_prim);
        v1_prim = auxyz;
        let v1 = GpMat::multiply_xyz_row(v1, &rot);
        let v2 = dir.cross(v1);

        let m = GpMat::from_rows(v1.x, v2.x, 0.0, v1.y, v2.y, 0.0, v1.z, v2.z, 0.0);
        let v2 = d_dir.cross(v1) + dir.cross(v1_prim);
        let m_prim = GpMat::from_rows(
            v1_prim.x, v2.x, 0.0, v1_prim.y, v2.y, 0.0, v1_prim.z, v2.z, 0.0,
        );

        let mut v2 = d_dir.cross(v1_prim);
        v2 *= 2.0;
        v2 += d2_dir.cross(v1) + dir.cross(v1_secn);
        let m_secn = GpMat::from_rows(
            v1_secn.x, v2.x, 0.0, v1_secn.y, v2.y, 0.0, v1_secn.z, v2.z, 0.0,
        );

        // Calcul des coeff -----------
        let beta = angle / 4.0;
        let betaprim = d_angle / 4.0;
        let betasecn = d2_angle / 4.0;
        let beta2 = beta * beta;
        let beta3 = beta * beta2;
        let beta4 = beta2 * beta2;
        let beta5 = beta3 * beta2;
        let beta6 = beta3 * beta3;
        let betaprim2 = betaprim * betaprim;

        // OCCT L458-459: double tan_b, dtan_b, d2tan_b, b, bpr, bsc (the tan
        // helpers are branch-local in Rust).
        let b: f64;
        let bpr: f64;
        let bsc: f64;
        if beta.abs() < NULL_ANGLE {
            // On calcul b par D.L
            let cf = -2.0 / 21.0;
            let num = 0.2 + cf * beta2;
            let denom = 1.0 + 0.2 * beta2;
            let aux = (cf * denom - 0.2 * num) / (denom * denom);
            b = -num / denom;
            bpr = -2.0 * beta * betaprim * aux;
            bsc = 2.0 * aux * (betaprim2 + beta * betasecn - 2.0 * beta * betaprim2);
        } else {
            let mut b_val = -1.0 / beta2;
            let mut bpr_val = (2.0 * betaprim) / beta3;
            let mut bsc_val = (2.0 * betasecn - 6.0 * betaprim * (betaprim / beta)) / beta3;
            if (std::f64::consts::PI / 2.0 - beta) > NULL_ANGLE {
                let tan_b = beta.tan();
                let dtan_b = betaprim * (1.0 + tan_b * tan_b);
                let d2tan_b = betasecn * (1.0 + tan_b * tan_b) + 2.0 * betaprim * tan_b * dtan_b;
                let b2 = tan_b - beta;
                b_val += beta / (3.0 * b2);
                let aux = betaprim * tan_b - beta * dtan_b;
                bpr_val += aux / (3.0 * b2 * b2);
                let daux = betasecn * tan_b - beta * d2tan_b;
                bsc_val += (daux - 2.0 * aux * betaprim * tan_b * tan_b / b2) / (3.0 * b2 * b2);
            }
            b = b_val;
            bpr = bpr_val;
            bsc = bsc_val;
        }

        let c = 1.0 / 3.0 + b;
        let b2 = b * b;
        let c2 = c * c;
        let bpr2 = bpr * bpr;

        // X = U*U - V*V
        self.vx.set(3, beta2 * (2.0 * b - 1.0));
        self.vx.set(5, beta4 * (b2 - 2.0 * c));
        self.vx.set(7, -beta6 * c2);
        // OCCT: DVx.Init(0).
        let mut dvx = Self::zero_vector();
        dvx.set(3, 2.0 * (beta * betaprim * (2.0 * b - 1.0) + bpr * beta2));
        dvx.set(5, 4.0 * beta3 * betaprim * (b2 - 2.0 * c) + 2.0 * beta4 * bpr * (b - 1.0));
        dvx.set(7, -6.0 * beta5 * betaprim * c2 - 2.0 * beta6 * bpr * c);
        // OCCT: D2Vx.Init(0).
        let mut d2vx = Self::zero_vector();
        d2vx.set(
            3,
            2.0 * ((betaprim2 + beta * betasecn) * (2.0 * b - 1.0)
                + 8.0 * beta * betaprim * bpr
                + bsc * beta2),
        );
        d2vx.set(
            5,
            4.0 * (b2 - 2.0 * c) * (3.0 * beta2 * betaprim2 + beta3 * betasecn)
                + 16.0 * beta3 * betaprim * bpr * (b - 1.0)
                + 2.0 * beta4 * (bsc * (b - 1.0) + bpr2),
        );
        d2vx.set(
            7,
            -6.0 * c2 * (5.0 * beta4 * betaprim2 + beta5 * betasecn)
                - 24.0 * beta5 * betaprim * bpr * c
                - 2.0 * beta6 * (bsc * c + bpr2),
        );

        // Y = 2*U*V
        self.vy.set(2, 2.0 * beta);
        self.vy.set(4, beta3 * 2.0 * (c + b));
        self.vy.set(6, 2.0 * beta5 * b * c);
        // OCCT: DVy.Init(0).
        let mut dvy = Self::zero_vector();
        dvy.set(2, 2.0 * betaprim);
        dvy.set(4, 6.0 * beta2 * betaprim * (b + c) + 4.0 * beta3 * bpr);
        dvy.set(6, 10.0 * beta4 * betaprim * b * c + 2.0 * beta5 * bpr * (b + c));
        // OCCT: D2Vy.Init(0).
        let mut d2vy = Self::zero_vector();
        d2vy.set(2, 2.0 * betasecn);
        d2vy.set(
            4,
            6.0 * (b + c) * (2.0 * beta * betaprim2 + beta2 * betasecn)
                + 24.0 * beta2 * betaprim * bpr * (b + c)
                + 4.0 * beta3 * bsc,
        );
        d2vy.set(
            6,
            10.0 * b * c * (4.0 * beta3 * betaprim2 + beta4 * betasecn)
                + 40.0 * beta4 * betaprim * bpr * (b + c)
                + 2.0 * beta5 * (bsc * (b + c) + 2.0 * bpr2),
        );

        // W = U*U + V*V
        self.vw.set(3, beta2 * (1.0 + 2.0 * b));
        self.vw.set(5, beta4 * (2.0 * c + b2));
        self.vw.set(7, beta6 * c2);
        // OCCT: DVw.Init(0).
        let mut dvw = Self::zero_vector();
        dvw.set(3, 2.0 * (beta * betaprim * (1.0 + 2.0 * b) + beta2 * bpr));
        dvw.set(5, 4.0 * beta3 * betaprim * (2.0 * c + b2) + 2.0 * beta4 * bpr * (b + 1.0));
        dvw.set(7, 6.0 * beta5 * betaprim * c2 + 2.0 * beta6 * bpr * c);
        // OCCT: D2Vw.Init(0).
        let mut d2vw = Self::zero_vector();
        d2vw.set(
            3,
            2.0 * ((betaprim2 + beta * betasecn) * (2.0 * b + 1.0)
                + 8.0 * beta * betaprim * bpr
                + bsc * beta2),
        );
        // OCCT L550: the literal OCCT expression keeps "b + 11" (reproduced
        // verbatim).
        d2vw.set(
            5,
            4.0 * (b2 + 2.0 * c) * (3.0 * beta2 * betaprim2 + beta3 * betasecn)
                + 16.0 * beta3 * betaprim * bpr * (b + 11.0)
                + 2.0 * beta4 * (bsc * (b + 1.0) + bpr2),
        );
        d2vw.set(
            7,
            6.0 * c2 * (5.0 * beta4 * betaprim2 + beta5 * betasecn)
                + 24.0 * beta5 * betaprim * bpr * c
                + 2.0 * beta6 * (bsc * c + bpr2),
        );

        // Calcul des poles — Px = B * Vx; ... D2W.Multiply(B, D2Vw).
        self.px = self.b.multiplied_vec(&self.vx);
        self.py = self.b.multiplied_vec(&self.vy);
        self.w = self.b.multiplied_vec(&self.vw);
        let dpx = self.b.multiplied_vec(&dvx);
        let dpy = self.b.multiplied_vec(&dvy);
        let dw = self.b.multiplied_vec(&dvw);
        let d2px = self.b.multiplied_vec(&d2vx);
        let d2py = self.b.multiplied_vec(&d2vy);
        let d2w = self.b.multiplied_vec(&d2vw);

        for ii in 1..=ordre {
            let wi = self.w.get(ii);
            let dwi = dw.get(ii);
            let mut p = DVec3::new(self.px.get(ii) / wi, self.py.get(ii) / wi, 0.0);
            let mut dp = DVec3::new(dpx.get(ii) / wi, dpy.get(ii) / wi, 0.0);

            let mut d2p = DVec3::new(d2px.get(ii) / wi, d2py.get(ii) / wi, 0.0);
            d2p -= 2.0 * (dwi / wi) * dp;
            d2p += (2.0 * (dwi / wi) * (dwi / wi) - d2w.get(ii) / wi) * p;
            dp -= (dw.get(ii) / wi) * p;

            poles[(ii - 1) as usize] = m.multiplied_xyz(p) + center;
            // OCCT: auxyz.SetLinearForm(1, MPrim * P, 1, M * DP, DCenter.XYZ()).
            let auxyz = m_prim.multiplied_xyz(p) + m.multiplied_xyz(dp) + d_center;
            dpoles[(ii - 1) as usize] = auxyz;
            // OCCT: P *= MSecn; DP *= MPrim; D2P *= M.
            p = GpMat::multiply_xyz_row(p, &m_secn);
            dp = GpMat::multiply_xyz_row(dp, &m_prim);
            d2p = GpMat::multiply_xyz_row(d2p, &m);
            // OCCT: auxyz.SetLinearForm(1, P, 2, DP, 1, D2P, D2Center.XYZ()).
            let auxyz = p + 2.0 * dp + d2p + d2_center;
            d2poles[(ii - 1) as usize] = auxyz;
            weights[(ii - 1) as usize] = wi;
            dweights[(ii - 1) as usize] = dwi;
            d2weights[(ii - 1) as usize] = d2w.get(ii);
        }
    }

    /// Helper for the OCCT `math_Vector(1, Ordre)` + `Init(0)` pair.
    fn zero_vector() -> Vector {
        Vector::new_init(1, 7, 0.0)
    }
}

impl Default for QuasiAngularConvertor {
    fn default() -> Self {
        Self::new()
    }
}
