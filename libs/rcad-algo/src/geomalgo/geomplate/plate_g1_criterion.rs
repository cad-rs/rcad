//! OCCT GeomPlate_PlateG1Criterion (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_PlateG1Criterion.hxx (L27-53) and GeomPlate_PlateG1Criterion.cxx
//! (whole file L31-133), served through the AdvApp2Var_Criterion base-class
//! interface (adv_app2_var::criterion::Criterion).

use glam::{DVec3, DVec2};
use rcad_kernel::math::plib::eval_poly2_var;

use crate::geomalgo::adv_app2_var::criterion::{Criterion, CriterionRepartition, CriterionType};
use crate::geomalgo::adv_app2_var::context::Context;
use crate::geomalgo::adv_app2_var::patch::Patch;

/// OCCT AdvApp2Var_CriterionType (AdvApp2Var_CriterionType.hxx L25-29) —
/// re-exported from the AdvApp2Var criterion module.
pub use crate::geomalgo::adv_app2_var::criterion::{CriterionRepartition as AdvApp2VarCriterionRepartition, CriterionType as AdvApp2VarCriterionType};

/// OCCT GeomPlate_PlateG1Criterion (hxx L28-46).
#[derive(Debug, Clone)]
pub struct PlateG1Criterion {
    /// hxx L47-48 (private): NCollection_Sequence<gp_XY> myData.
    my_data: Vec<DVec2>,
    /// hxx L49: NCollection_Sequence<gp_XYZ> myXYZ.
    my_xyz: Vec<DVec3>,
    /// OCCT AdvApp2Var_Criterion base members (myMaxValue / myType /
    /// myRepartition, hxx L50-52 of AdvApp2Var_Criterion.hxx).
    my_max_value: f64,
    my_type: AdvApp2VarCriterionType,
    my_repartition: AdvApp2VarCriterionRepartition,
}

impl PlateG1Criterion {
    /// OCCT ctor (GeomPlate_PlateG1Criterion.cxx L31-42) — default args
    /// Type = AdvApp2Var_Absolute, Repart = AdvApp2Var_Regular (hxx L38-39).
    pub fn new(
        data: &[DVec2],
        g1data: &[DVec3],
        maximum: f64,
        ctype: AdvApp2VarCriterionType,
        repart: AdvApp2VarCriterionRepartition,
    ) -> Self {
        // myData = Data; myXYZ = G1Data; myMaxValue = Maximum;
        // myType = Type; myRepartition = Repart;
        PlateG1Criterion {
            my_data: data.to_vec(),
            my_xyz: g1data.to_vec(),
            my_max_value: maximum,
            my_type: ctype,
            my_repartition: repart,
        }
    }
}

impl Criterion for PlateG1Criterion {
    /// OCCT Value(P, C) (GeomPlate_PlateG1Criterion.cxx L46-125).
    fn value(&self, p: &mut Patch, c: &Context) {
        // double UInt[2], VInt[2]; int MaxNbCoeff[2], NbCoeff[2];
        // adrCoeff = (double*)&P.Coefficients(1, C)->ChangeArray1()(
        //     P.Coefficients(1, C)->Lower());
        let adr_coeff = p.coefficients(1, c);
        // MaxNbCoeff[0] = C.ULimit(); MaxNbCoeff[1] = C.VLimit();
        let max_nb_coeff = [c.u_limit(), c.v_limit()];
        // NbCoeff[0] = P.NbCoeffInU(); NbCoeff[1] = P.NbCoeffInV();
        let nb_coeff = [p.nb_coeff_in_u(), p.nb_coeff_in_v()];
        // UInt[0] = P.U0(); UInt[1] = P.U1(); VInt[0] = P.V0(); VInt[1] = P.V1();
        let uint = [p.u0(), p.u1()];
        let vint = [p.v0(), p.v1()];

        // double up, vp, ang = 0.;
        let mut ang = 0.0f64;

        // int dimension = 3 * NbCoeff[1];
        let dimension = 3 * nb_coeff[1];
        // NCollection_Array1<double> Patch(1, NbCoeff[0] * dimension);
        let mut patch = vec![0.0f64; (nb_coeff[0] * dimension) as usize];
        // NCollection_Array1<double> Curve(1, 2 * dimension) — declared but
        // unused in the OCCT G1 body.
        let _curve = vec![0.0f64; (2 * dimension) as usize];
        // NCollection_Array1<double> Point(1, 3);
        let mut point = vec![0.0f64; 3usize];

        // int k1, k2, pos, ll = 1;
        let mut ll = 1i32;
        for k1 in 1..=nb_coeff[0] {
            // JAG 99.04.29    pos = 3*(MaxNbCoeff[0])*(k1-1);
            // pos = 3 * (MaxNbCoeff[1]) * (k1 - 1);
            let mut pos = 3 * max_nb_coeff[1] * (k1 - 1);
            for k2 in 1..=nb_coeff[1] {
                // Patch(ll)     = adrCoeff[pos];
                // Patch(ll + 1) = adrCoeff[pos + 1];
                // Patch(ll + 2) = adrCoeff[pos + 2];
                // ll += 3; pos += 3;
                patch[(ll - 1) as usize] = adr_coeff[pos as usize];
                patch[ll as usize] = adr_coeff[(pos + 1) as usize];
                patch[(ll + 1) as usize] = adr_coeff[(pos + 2) as usize];
                ll += 3;
                pos += 3;
            }
        }

        // int i, NbCtr = myData.Length();
        let nb_ctr = self.my_data.len();
        for i in 1..=nb_ctr {
            // gp_Vec v1s, v2s, v3s;
            // gp_Vec v3h(myXYZ.Value(i).X(), myXYZ.Value(i).Y(), myXYZ.Value(i).Z());
            let v3h = self.my_xyz[i - 1];
            // gp_XY P2d = myData.Value(i);
            let p2d = self.my_data[i - 1];
            if uint[0] < p2d.x && p2d.x < uint[1] && vint[0] < p2d.y && p2d.y < vint[1] {
                // u,v recadres sur (-1,1)
                // up = (2 * P2d.X() - UInt[0] - UInt[1]) / (UInt[1] - UInt[0]);
                // vp = (2 * P2d.Y() - VInt[0] - VInt[1]) / (VInt[1] - VInt[0]);
                let up = (2.0 * p2d.x - uint[0] - uint[1]) / (uint[1] - uint[0]);
                let vp = (2.0 * p2d.y - vint[0] - vint[1]) / (vint[1] - vint[0]);
                // PLib::EvalPoly2Var(up, vp, 1, 0, NbCoeff[0] - 1,
                //     NbCoeff[1] - 1, 3, Coeffs[0], Digit[0]);
                eval_poly2_var(
                    up,
                    vp,
                    1,
                    0,
                    (nb_coeff[0] - 1) as usize,
                    (nb_coeff[1] - 1) as usize,
                    3,
                    &patch,
                    &mut point,
                );

                // v1s.SetCoord(1, Digit[0]); SetCoord(2, Digit[1]);
                // SetCoord(3, Digit[2]);
                let v1s = DVec3::new(point[0], point[1], point[2]);

                // PLib::EvalPoly2Var(up, vp, 0, 1, NbCoeff[0] - 1,
                //     NbCoeff[1] - 1, 3, Coeffs[0], Digit[0]);
                eval_poly2_var(
                    up,
                    vp,
                    0,
                    1,
                    (nb_coeff[0] - 1) as usize,
                    (nb_coeff[1] - 1) as usize,
                    3,
                    &patch,
                    &mut point,
                );

                // v2s.SetCoord(1, Digit[0]); SetCoord(2, Digit[1]);
                // SetCoord(3, Digit[2]);
                let v2s = DVec3::new(point[0], point[1], point[2]);

                // v3s = v1s ^ v2s;
                let v3s = v1s.cross(v2s);
                // if (v3s.Angle(v3h) > (M_PI / 2))
                // {
                //   if ((M_PI - v3s.Angle(v3h)) > ang) ang = (M_PI - v3s.Angle(v3h));
                // }
                // else
                // {
                //   if (v3s.Angle(v3h) > ang) ang = v3s.Angle(v3h);
                // }
                let a = gp_vec_angle(v3s, v3h);
                if a > (std::f64::consts::PI / 2.0) {
                    if (std::f64::consts::PI - a) > ang {
                        ang = std::f64::consts::PI - a;
                    }
                } else if a > ang {
                    ang = a;
                }
            }
        }
        // P.SetCritValue(ang);
        p.set_crit_value(ang);
    }

    /// OCCT IsSatisfied(P) (GeomPlate_PlateG1Criterion.cxx L130-133).
    fn is_satisfied(&self, p: &Patch) -> bool {
        p.crit_value() < self.my_max_value
    }

    /// OCCT AdvApp2Var_Criterion::MaxValue (AdvApp2Var_Criterion.cxx L23-26).
    fn max_value(&self) -> f64 {
        self.my_max_value
    }

    /// OCCT AdvApp2Var_Criterion::Type (AdvApp2Var_Criterion.cxx L30-33).
    fn crit_type(&self) -> CriterionType {
        match self.my_type {
            AdvApp2VarCriterionType::Absolute => CriterionType::Absolute,
            AdvApp2VarCriterionType::Relative => CriterionType::Relative,
        }
    }

    /// OCCT AdvApp2Var_Criterion::Repartition
    /// (AdvApp2Var_Criterion.cxx L37-40).
    fn repartition(&self) -> CriterionRepartition {
        match self.my_repartition {
            AdvApp2VarCriterionRepartition::Regular => CriterionRepartition::Regular,
            AdvApp2VarCriterionRepartition::Incremental => CriterionRepartition::Incremental,
        }
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-494) — raises
/// VectorWithNullMagnitude when a magnitude is null, then delegates to
/// gp_Dir::Angle.  (gp::Resolution() — OCCT Standard_Real twin of
/// DBL_MIN; modeled by f64::MIN_POSITIVE.)
fn gp_vec_angle(v: DVec3, other: DVec3) -> f64 {
    if v.length() <= f64::MIN_POSITIVE || other.length() <= f64::MIN_POSITIVE {
        panic!("gp_VectorWithNullMagnitude");
    }
    // OCCT gp_Dir::Angle (gp_Dir.cxx L27-52): above 45 degrees the arccos
    // gives the best precision; otherwise the arcsin.  In 3d the angular
    // values are always positive and lie between 0 and PI.
    let cosinus = v.dot(other) / (v.length() * other.length());
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = v.cross(other).length() / (v.length() * other.length());
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}
