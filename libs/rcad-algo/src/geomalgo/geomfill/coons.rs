//! OCCT GeomFill_Coons (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Coons.cxx (whole file L25-201).
//!
//! Mapping: `NCollection_Array1<gp_Pnt>` -> `&[DVec3]`, `NCollection_Array2`
//! -> `Vec<Vec<DVec3>>` indexed `[i - 1][j - 1]` inside `1..=n` loops.
//! Dependencies: PLib::CoefficientsPoles (rcad plib::coefficients_poles) and
//! BSplCLib::IncreaseDegree for Bezier poles (rcad
//! bspl_lib::increase_degree with the Bezier knot vector).

use glam::DVec3;

use rcad_kernel::math::bspl_lib::increase_degree as bspl_increase_degree;
use rcad_kernel::math::plib::coefficients_poles;

use super::filling::FillingBase;

/// OCCT GeomFill_Coons.
#[derive(Debug, Clone, Default)]
pub struct Coons {
    pub(crate) base: FillingBase,
}

/// OCCT BSplCLib::IncreaseDegree(NewDegree, Poles, NoWeights, NewPoles,
/// NoWeights) for Bezier poles (BSplCLib.hxx L800-804, dispatched through
/// BSplCLib_BzSyntaxes.cxx L47-58 -> BSplCLib_CurveComputation.pxx
/// L2080-2099): the Bezier knot vector is {0, 1} with multiplicities
/// {Degree + 1, Degree + 1}.
fn increase_degree_bezier(new_degree: usize, poles: &[DVec3]) -> Vec<DVec3> {
    let degree = poles.len() - 1;
    let mut flat: Vec<f64> = Vec::with_capacity(poles.len() * 3);
    for p in poles {
        flat.extend([p.x, p.y, p.z]);
    }
    let knots = [0.0f64, 1.0];
    let mults = [degree as i32 + 1, degree as i32 + 1];
    // Bezier: one span; new pole count = Poles.Length() + (NewDegree - Degree).
    let nb_new_poles = poles.len() + (new_degree - degree);
    let mut new_flat = vec![0.0f64; nb_new_poles * 3];
    let mut new_knots = vec![0.0f64; 2];
    let mut new_mults = vec![0i32; 2];
    bspl_increase_degree(
        degree,
        new_degree,
        false,
        3,
        &flat,
        &knots,
        &mults,
        &mut new_flat,
        &mut new_knots,
        &mut new_mults,
    );
    new_flat
        .chunks(3)
        .map(|c| DVec3::new(c[0], c[1], c[2]))
        .collect()
}

impl Coons {
    /// OCCT GeomFill_Coons(P1, P2, P3, P4) (L31-37).
    pub fn new(p1: &[DVec3], p2: &[DVec3], p3: &[DVec3], p4: &[DVec3]) -> Self {
        let mut coons = Coons { base: FillingBase::new() };
        coons.init(p1, p2, p3, p4);
        coons
    }

    /// OCCT GeomFill_Coons(P1, P2, P3, P4, W1, W2, W3, W4) (L39-46).
    pub fn new_rational(
        p1: &[DVec3],
        p2: &[DVec3],
        p3: &[DVec3],
        p4: &[DVec3],
        w1: &[f64],
        w2: &[f64],
        w3: &[f64],
        w4: &[f64],
    ) -> Self {
        let mut coons = Coons { base: FillingBase::new() };
        coons.init_rational(p1, p2, p3, p4, w1, w2, w3, w4);
        coons
    }

    /// OCCT Init(P1, P2, P3, P4) (L48-139).
    pub fn init(&mut self, p1: &[DVec3], p2: &[DVec3], p3: &[DVec3], p4: &[DVec3]) {
        assert!(
            p1.len() == p3.len() && p2.len() == p4.len(),
            "Standard_DomainError"
        );
        let npolu = p1.len();
        let npolv = p2.len();
        self.base.is_rational = false;
        self.base.poles = vec![vec![DVec3::ZERO; npolv]; npolu];
        // The boundaries are not modified
        for i in 1..=npolu {
            self.base.poles[i - 1][0] = p1[i - 1];
            self.base.poles[i - 1][npolv - 1] = p3[i - 1];
        }
        for i in 1..=npolv {
            self.base.poles[0][i - 1] = p2[i - 1];
            self.base.poles[npolu - 1][i - 1] = p4[i - 1];
        }

        // Calcul des coefficients multiplicateurs
        // Coef(1) = (1, 0, 0); Coef(2) = (0, 0, 0); Coef(3) = (-3, 3, 0);
        // Coef(4) = (2, -2, 0)  (power coefficients of F(t) = 1 - 3t^2 + 2t^3).
        let coef = [
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(-3.0, 3.0, 0.0),
            DVec3::new(2.0, -2.0, 0.0),
        ];
        let pole = coefficients_poles(&coef);
        // OCCT `CoefU = Pole` / IncreaseDegree is reached only for NPolU >= 4
        // (GeomFill_BSplineCurves guards the CoonsStyle pole count).
        let coef_u = if npolu > 4 {
            increase_degree_bezier(npolu - 1, &pole)
        } else {
            pole.clone()
        };
        let coef_v = if npolv > 4 {
            increase_degree_bezier(npolv - 1, &pole)
        } else {
            pole.clone()
        };

        let mut fu = vec![0.0f64; npolu];
        let mut gu = vec![0.0f64; npolu];
        let mut fv = vec![0.0f64; npolv];
        let mut gv = vec![0.0f64; npolv];
        // OCCT: for (i = 2; i < NPolU; i++) — exclusive upper bound.
        for i in 2..npolu {
            fu[i - 1] = coef_u[i - 1].x;
            gu[i - 1] = coef_u[i - 1].y;
        }
        for i in 2..npolv {
            fv[i - 1] = coef_v[i - 1].x;
            gv[i - 1] = coef_v[i - 1].y;
        }

        // Calcul des poles interieurs
        for j in 2..npolv {
            for i in 2..npolu {
                let p = fv[j - 1] * self.base.poles[i - 1][0]
                    + gv[j - 1] * self.base.poles[i - 1][npolv - 1]
                    + fu[i - 1] * self.base.poles[0][j - 1]
                    + gu[i - 1] * self.base.poles[npolu - 1][j - 1]
                    - fu[i - 1] * fv[j - 1] * self.base.poles[0][0]
                    - fu[i - 1] * gv[j - 1] * self.base.poles[0][npolv - 1]
                    - gu[i - 1] * fv[j - 1] * self.base.poles[npolu - 1][0]
                    - gu[i - 1] * gv[j - 1] * self.base.poles[npolu - 1][npolv - 1];
                self.base.poles[i - 1][j - 1] = p;
            }
        }
    }

    /// OCCT Init(P1, P2, P3, P4, W1, W2, W3, W4) (L142-201).
    pub fn init_rational(
        &mut self,
        p1: &[DVec3],
        p2: &[DVec3],
        p3: &[DVec3],
        p4: &[DVec3],
        w1: &[f64],
        w2: &[f64],
        w3: &[f64],
        w4: &[f64],
    ) {
        assert!(
            w1.len() == w3.len() && w2.len() == w4.len(),
            "Standard_DomainError"
        );
        assert!(
            w1.len() == p1.len()
                && w2.len() == p2.len()
                && w3.len() == p3.len()
                && w4.len() == p4.len(),
            "Standard_DomainError"
        );
        self.init(p1, p2, p3, p4);
        self.base.is_rational = true;
        let npolu = w1.len();
        let npolv = w2.len();
        let nu = (npolu - 1) as f64;
        let nv = (npolv - 1) as f64;
        self.base.weights = vec![vec![0.0f64; npolv]; npolu];
        // The boundaries are not modified
        for i in 1..=npolu {
            self.base.weights[i - 1][0] = w1[i - 1];
            self.base.weights[i - 1][npolv - 1] = w3[i - 1];
        }
        for j in 2..=npolv - 1 {
            let pv = (j - 1) as f64 / nv;
            let pv1 = 1.0 - pv;
            self.base.weights[0][j - 1] = w4[j - 1];
            self.base.weights[npolu - 1][j - 1] = w2[j - 1];
            for i in 2..=npolu - 1 {
                let pu = (i - 1) as f64 / nu;
                let pu1 = 1.0 - pu;
                let w = pv1 * w1[i - 1]
                    + pv * w3[i - 1]
                    + pu * w2[j - 1]
                    + pu1 * w4[j - 1]
                    - (pu1 * pv1 * w1[0]
                        + pu * pv1 * w2[0]
                        + pu * pv * w3[npolu - 1]
                        + pu1 * pv * w4[npolv - 1]);
                self.base.weights[i - 1][j - 1] = w;
            }
        }
    }
}
