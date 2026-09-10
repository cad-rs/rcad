// OCCT math_Uzawa (TKMath/math/math_Uzawa.hxx/.cxx/.lxx) — 1:1 Rust
// translation: solves Cont*X = Secont (the Nce first equations are equal
// equations and the Nci last equations are inequalities <) with
// minimization of Norme(X-X0) via Uzawa's dual iteration; the Nci == 0 case
// resolves directly through math_Crout on Cont*Transposed(Cont).
//
// Note (kept from OCCT): ctor 1 forwards (Cont.RowNumber(), 0) into
// Perform's (Nci, Nce) parameter slots and ctor 2 forwards its (Nce, Nci)
// parameters in the same positional order — the C++ source's parameter
// naming asymmetry is preserved verbatim.

use super::math_crout::Crout;
use super::math_matrix::{Matrix, Vector};

/// OCCT math_Uzawa.
#[derive(Debug, Clone)]
pub struct Uzawa {
    resul: Vector,
    erruza: Vector,
    errinit: Vector,
    vardua: Vector,
    ctcinv: Matrix,
    nb_iter: i32,
    done: bool,
}

impl Uzawa {
    /// OCCT math_Uzawa(Cont, Secont, StartingPoint, EpsLix = 1.0e-06,
    /// EpsLic = 1.0e-06, NbIterations = 500) (math_Uzawa.cxx L16-31).
    pub fn new(
        cont: &Matrix,
        secont: &Vector,
        starting_point: &Vector,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut u = Uzawa::make_members(cont);
        u.perform(
            cont,
            secont,
            starting_point,
            cont.row_number(),
            0,
            eps_lix,
            eps_lic,
            nb_iterations,
        );
        u
    }

    /// OCCT math_Uzawa(Cont, Secont, StartingPoint, Nce, Nci, EpsLix,
    /// EpsLic, NbIterations) (math_Uzawa.cxx L33-48).
    #[allow(clippy::too_many_arguments)]
    pub fn with_constraints(
        cont: &Matrix,
        secont: &Vector,
        starting_point: &Vector,
        nce: i32,
        nci: i32,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut u = Uzawa::make_members(cont);
        u.perform(cont, secont, starting_point, nce, nci, eps_lix, eps_lic, nb_iterations);
        u
    }

    /// The shared member-initializer list (math_Uzawa.cxx L17-29 / L34-46).
    fn make_members(cont: &Matrix) -> Self {
        Uzawa {
            resul: Vector::new(1, cont.col_number()),
            erruza: Vector::new(1, cont.col_number()),
            errinit: Vector::new(1, cont.col_number()),
            vardua: Vector::new(1, cont.row_number()),
            ctcinv: Matrix::new(1, cont.row_number(), 1, cont.row_number()),
            nb_iter: 0,
            done: false,
        }
    }

    /// OCCT Perform(Cont, Secont, StartingPoint, Nci, Nce, EpsLix, EpsLic,
    /// NbIterations) (math_Uzawa.cxx L50-241).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        cont: &Matrix,
        secont: &Vector,
        starting_point: &Vector,
        nci: i32,
        nce: i32,
        eps_lix: f64,
        eps_lic: f64,
        nb_iterations: i32,
    ) {
        let coef = 1.0 / 2.0f64.sqrt();
        let nlig = cont.row_number();
        let ncol = cont.col_number();
        // Standard_DimensionError_Raise_if((Secont.Length() != Nlig)
        //     || ((Nce + Nci) != Nlig), " ").
        assert!(
            secont.length() == nlig && (nce + nci) == nlig,
            "Standard_DimensionError: math_Uzawa::Perform"
        );
        // Calcul du vecteur Cont*X0 - D:  (erreur initiale)
        //==================================================
        for i in 1..=nlig {
            let mut v = cont.get(i, 1) * starting_point.get(1) - secont.get(i);
            for j in 2..=ncol {
                v += cont.get(i, j) * starting_point.get(j);
            }
            self.errinit.set(i, v);
        }
        if nci == 0 {
            // cas de resolution directe
            self.nb_iter = 1; //==========================
            // Calcul de Cont*T(Cont)
            for i in 1..=nlig {
                for j in 1..=i {
                    // a utiliser avec Crout.
                    let mut v = cont.get(i, 1) * cont.get(j, 1);
                    for k in 2..=ncol {
                        v += cont.get(i, k) * cont.get(j, k);
                    }
                    self.ctcinv.set(i, j, v);
                }
            }
            // Calcul de l inverse de CTCinv :
            //================================
            let inv = Crout::new(&self.ctcinv, 1.0e-20); // utilisation de Crout.
            self.ctcinv = inv.inverse().clone();
            for i in 1..=nlig {
                let mut scale = self.ctcinv.get(i, 1) * self.errinit.get(1);
                for j in 2..=i {
                    scale += self.ctcinv.get(i, j) * self.errinit.get(j);
                }
                for j in (i + 1)..=nlig {
                    scale += self.ctcinv.get(j, i) * self.errinit.get(j);
                }
                self.vardua.set(i, scale);
            }
            for i in 1..=ncol {
                let mut v = -cont.get(1, i) * self.vardua.get(1);
                for j in 2..=nlig {
                    v -= cont.get(j, i) * self.vardua.get(j);
                }
                self.erruza.set(i, v);
            }
            // restitution des valeurs calculees:
            //===================================
            // Resul = StartingPoint + Erruza
            for i in 1..=ncol {
                let v = starting_point.get(i) + self.erruza.get(i);
                self.resul.set(i, v);
            }
            self.done = true;
        } else {
            // Initialisation des variables duales.
            //=====================================
            for i in 1..=nlig {
                if i <= nce {
                    self.vardua.set(i, 0.0);
                } else {
                    self.vardua.set(i, 1.0);
                }
            }
            // Calcul du coefficient Rho:
            //===========================
            let mut normat = 0.0;
            for i in 1..=nlig {
                let mut normli = cont.get(i, 1) * cont.get(i, 1);
                for j in 2..=ncol {
                    normli += cont.get(i, j) * cont.get(i, j);
                }
                normat += normli;
            }
            let rho = coef / normat;
            // Boucle des iterations de la methode d Uzawa.
            //=============================================
            let mut xmax = 0.0;
            let mut errmax;
            self.nb_iter = 1;
            while self.nb_iter <= nb_iterations {
                errmax = 0.0;
                for i in 1..=ncol {
                    let xian = self.erruza.get(i);
                    let mut v = -cont.get(1, i) * self.vardua.get(1);
                    for j in 2..=nlig {
                        v -= cont.get(j, i) * self.vardua.get(j);
                    }
                    self.erruza.set(i, v);
                    if self.nb_iter > 1 {
                        let diff = (self.erruza.get(i) - xian).abs();
                        if i == 1 {
                            xmax = diff;
                        }
                        xmax = xmax.max(diff);
                    }
                }
                // Calcul de Xmu a l iteration NbIter et evaluation de l erreur
                // sur la verification des contraintes.
                //=============================================================
                for i in 1..=nlig {
                    let mut err = cont.get(i, 1) * self.erruza.get(1) + self.errinit.get(i);
                    for j in 2..=ncol {
                        err += cont.get(i, j) * self.erruza.get(j);
                    }
                    let err1;
                    if i <= nce {
                        let v = self.vardua.get(i) + rho * err;
                        self.vardua.set(i, v);
                        err1 = (rho * err).abs();
                    } else {
                        let xmuian = self.vardua.get(i);
                        let v = (self.vardua.get(i) + rho * err).max(0.0);
                        self.vardua.set(i, v);
                        err1 = (v - xmuian).abs();
                    }
                    if i == 1 {
                        errmax = err1;
                    }
                    errmax = errmax.max(err1);
                }
                if self.nb_iter > 1 {
                    if xmax <= eps_lix {
                        if errmax <= eps_lic {
                            // Convergence atteinte dans Uzawa
                            self.done = true;
                        } else {
                            // convergence non atteinte pour le probleme dual
                            self.done = false;
                            return;
                        }
                        // Restitution des valeurs calculees
                        //==================================
                        // Resul = StartingPoint + Erruza
                        for i in 1..=ncol {
                            let v = starting_point.get(i) + self.erruza.get(i);
                            self.resul.set(i, v);
                        }
                        self.done = true;
                        return;
                    }
                }
                self.nb_iter += 1;
            } // fin de la boucle d iterations.
            self.done = false;
        }
    }

    /// OCCT Duale(V) (math_Uzawa.cxx L243-246) — the dual variables.
    pub fn duale(&self, v: &mut Vector) {
        for i in v.lower()..=v.upper() {
            let x = self.vardua.get(i);
            v.set(i, x);
        }
    }

    /// OCCT IsDone() (lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Value() (lxx) — the vector solution of the system.
    pub fn value(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: math_Uzawa::Value");
        &self.resul
    }

    /// OCCT InitialError() (lxx) — Cont*StartingPoint-Secont.
    pub fn initial_error(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: math_Uzawa::InitialError");
        &self.errinit
    }

    /// OCCT Error() (lxx) — the difference between X and the StartingPoint.
    pub fn error(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: math_Uzawa::Error");
        &self.erruza
    }

    /// OCCT NbIterations() (lxx) — the iterations really done.
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone: math_Uzawa::NbIterations");
        self.nb_iter
    }

    /// OCCT InverseCont() (lxx) — the inverse matrix of (C * Transposed(C)).
    pub fn inverse_cont(&self) -> &Matrix {
        assert!(self.done, "StdFail_NotDone: math_Uzawa::InverseCont");
        &self.ctcinv
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Equality-constrained projection: minimize |X - X0| subject to
    // Cont*X = Secont — ctor 1 iterative branch and ctor 2 direct branch.
    #[test]
    fn direct_resolution_projection() {
        // Cont = I (2x2): solution X must equal Secont exactly, with
        // StartingPoint = 0 the correction is Erruza = Secont.
        let mut cont = Matrix::new(1, 2, 1, 2);
        cont.set(1, 1, 1.0);
        cont.set(1, 2, 0.0);
        cont.set(2, 1, 0.0);
        cont.set(2, 2, 1.0);
        let mut secont = Vector::new(1, 2);
        secont.set(1, 3.0);
        secont.set(2, 4.0);
        let starting_point = Vector::new_init(1, 2, 0.0);

        // OCCT ctor 1 forwards Cont.RowNumber() into Perform's Nci slot, so
        // `new` ALWAYS takes the iterative branch (Nci != 0) with all rows
        // treated as inequalities (Nce = 0). For this input the starting
        // point already satisfies the inequalities (I*X0 <= Secont), so the
        // minimizer of |X - X0| is X0 itself and the duals clamp to 0 —
        // the faithful outcome of the positional semantics.
        let u = Uzawa::new(&cont, &secont, &starting_point, 1.0e-6, 1.0e-6, 500);
        assert!(u.is_done());
        assert!(u.value().get(1).abs() < 1.0e-5);
        assert!(u.value().get(2).abs() < 1.0e-5);

        // The direct (Crout) branch requires ctor 2 with its first
        // constraint count = 0 (it lands in Perform's Nci slot).
        let u2 = Uzawa::with_constraints(&cont, &secont, &starting_point, 0, 2, 1.0e-6, 1.0e-6, 500);
        assert!(u2.is_done());
        assert!((u2.value().get(1) - 3.0).abs() < 1.0e-10);
        assert!((u2.value().get(2) - 4.0).abs() < 1.0e-10);
        assert_eq!(u2.nb_iterations(), 1);

        // InitialError = Cont*X0 - Secont = -Secont.
        assert!((u.initial_error().get(1) + 3.0).abs() < 1.0e-10);
        assert!((u.initial_error().get(2) + 4.0).abs() < 1.0e-10);

        // InverseCont is only assembled by the direct branch; on u2 it must
        // be I (inverse of I*I = I).
        assert!((u2.inverse_cont().get(1, 1) - 1.0).abs() < 1.0e-10);
        assert!((u2.inverse_cont().get(2, 2) - 1.0).abs() < 1.0e-10);
        assert!((u2.inverse_cont().get(1, 2)).abs() < 1.0e-10);
    }

    // A 3x3 identity system through the same direct branch, exercising
    // CTCinv assembly and the triangular Vardua accumulation.
    #[test]
    fn direct_resolution_identity_3x3() {
        let mut cont = Matrix::new(1, 3, 1, 3);
        for i in 1..=3 {
            for j in 1..=3 {
                cont.set(i, j, if i == j { 2.0 } else { 0.0 });
            }
        }
        let mut secont = Vector::new(1, 3);
        secont.set(1, 2.0);
        secont.set(2, 4.0);
        secont.set(3, 6.0);
        let starting_point = Vector::new_init(1, 3, 0.0);

        // Direct branch through ctor 2 (first count = 0).
        let u = Uzawa::with_constraints(&cont, &secont, &starting_point, 0, 3, 1.0e-6, 1.0e-6, 500);
        assert!(u.is_done());
        assert!((u.value().get(1) - 1.0).abs() < 1.0e-10);
        assert!((u.value().get(2) - 2.0).abs() < 1.0e-10);
        assert!((u.value().get(3) - 3.0).abs() < 1.0e-10);
    }
}
