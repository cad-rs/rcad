//! OCCT math_GaussSetIntegration + math_FunctionSet (FoundationClasses/
//! TKMath/math), hosted under `fem_tool` because the `math` module
//! registration is closed and this class is only consumed by the FEmTool
//! criteria (FEmTool_LinearTension / LinearFlexion / LinearJerk).
//!
//! 1:1 translation of `math_GaussSetIntegration.hxx` (L17-77),
//! `math_GaussSetIntegration.lxx` (L17-37) and `math_GaussSetIntegration.cxx`
//! (L17-103), plus `math_FunctionSet.hxx` (L25-70) as a trait.

use super::super::gauss_points::{gauss_points, gauss_weights, gauss_points_max};
use super::super::math_matrix::{IntegerVector, Vector};

/// OCCT math_FunctionSet (math_FunctionSet.hxx L25-70) - a set of N functions
/// of M independent variables.
pub trait FunctionSet {
    /// OCCT math_FunctionSet::NbVariables - returns the number of variables
    /// of the function.
    fn nb_variables(&self) -> i32;

    /// OCCT math_FunctionSet::NbEquations - returns the number of equations
    /// of the function.
    fn nb_equations(&self) -> i32;

    /// OCCT math_FunctionSet::Value(X, F) - computes the values F of the
    /// functions for the variable X.
    fn value(&mut self, x: &Vector, f: &mut Vector) -> bool;

    /// OCCT math_FunctionSet::GetStateNumber - by default returns 0.
    fn get_state_number(&self) -> i32 {
        0
    }
}

/// OCCT math_GaussSetIntegration (math_GaussSetIntegration.hxx L29-77) -
/// Gauss-Legendre integration of a set of N functions of M variables between
/// the bounds Lower and Upper (the case M > 1 is not implemented).
#[derive(Debug, Clone)]
pub struct GaussSetIntegration {
    /// OCCT `Val`.
    val: Vector,
    /// OCCT `Done`.
    done: bool,
}

impl GaussSetIntegration {
    /// OCCT math_GaussSetIntegration::math_GaussSetIntegration
    /// (math_GaussSetIntegration.cxx L28-100).
    pub fn new(f: &mut dyn FunctionSet, lower: &Vector, upper: &Vector, order: &IntegerVector) -> Self {
        let nb_equa = f.nb_equations();
        let nb_var = f.nb_variables();

        let mut f_val1 = Vector::new(1, nb_equa);
        let mut f_val2 = Vector::new(1, nb_equa);
        let mut tval = Vector::new(1, nb_var);

        // Verification
        assert!(
            nb_var == 1 && order.get(order.lower()) <= gauss_points_max() as i32,
            "Standard_NotImplemented: GaussSetIntegration "
        );

        // Initialisations
        let done = false;

        let xdeb = lower.get(lower.lower());
        let xfin = upper.get(upper.lower());
        let ordre = order.get(order.lower());
        let mut val = Vector::new(1, nb_equa);
        let mut gauss_p = Vector::new(1, ordre);
        let mut gauss_w = Vector::new(1, ordre);

        // Recuperation des points de Gauss dans le fichier GaussPoints.
        gauss_points(ordre as usize, &mut gauss_p.data);
        gauss_weights(ordre as usize, &mut gauss_w.data);

        // Changement de variable pour la mise a l'echelle [Lower, Upper] :
        let xm = 0.5 * (xdeb + xfin);
        let xr = 0.5 * (xfin - xdeb);

        let ind = ordre / 2;
        let ind1 = (ordre + 1) / 2;
        let mut is_ok;
        if ind1 > ind {
            // odder case
            tval.set(1, xm); // +  Xr * GaussP(ind1);
            is_ok = f.value(&tval, &mut val);
            if !is_ok {
                return GaussSetIntegration { val, done };
            }
            // Val *= GaussW(ind1);
            let w = gauss_w.get(ind1);
            for r in 1..=val.length() {
                val.set(r, val.get(r) * w);
            }
        } else {
            // Val.Init(0);
            for r in 1..=val.length() {
                val.set(r, 0.0);
            }
        }

        for i in 1..=ind {
            tval.set(1, xm + xr * gauss_p.get(i));
            is_ok = f.value(&tval, &mut f_val1);
            if !is_ok {
                return GaussSetIntegration { val, done };
            }
            tval.set(1, xm - xr * gauss_p.get(i));
            is_ok = f.value(&tval, &mut f_val2);
            if !is_ok {
                return GaussSetIntegration { val, done };
            }
            // FVal1 += FVal2;
            for r in 1..=nb_equa {
                f_val1.set(r, f_val1.get(r) + f_val2.get(r));
            }
            // FVal1 *= GaussW(i);
            let w = gauss_w.get(i);
            for r in 1..=nb_equa {
                f_val1.set(r, f_val1.get(r) * w);
            }
            // Val += FVal1;
            for r in 1..=nb_equa {
                val.set(r, val.get(r) + f_val1.get(r));
            }
        }
        // Val *= Xr;
        for r in 1..=val.length() {
            val.set(r, val.get(r) * xr);
        }

        GaussSetIntegration { val, done: true }
    }

    /// OCCT math_GaussSetIntegration::IsDone (lxx L27-30).
    #[inline]
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT math_GaussSetIntegration::Value (lxx L32-36).
    #[inline]
    pub fn value(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: Integration ");
        &self.val
    }
}
