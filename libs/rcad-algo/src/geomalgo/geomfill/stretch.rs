//! OCCT GeomFill_Stretch (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Stretch.cxx (whole file L25-140).
//!
//! Mapping: `NCollection_Array1<gp_Pnt>` -> `&[DVec3]`, `NCollection_Array2`
//! -> `Vec<Vec<DVec3>>` indexed `[i - 1][j - 1]` inside `1..=n` loops.

use glam::DVec3;

use super::filling::FillingBase;

/// OCCT GeomFill_Stretch.
#[derive(Debug, Clone, Default)]
pub struct Stretch {
    pub(crate) base: FillingBase,
}

impl Stretch {
    /// OCCT GeomFill_Stretch(P1, P2, P3, P4) (L32-38).
    pub fn new(p1: &[DVec3], p2: &[DVec3], p3: &[DVec3], p4: &[DVec3]) -> Self {
        let mut stretch = Stretch { base: FillingBase::new() };
        stretch.init(p1, p2, p3, p4);
        stretch
    }

    /// OCCT GeomFill_Stretch(P1, P2, P3, P4, W1, W2, W3, W4) (L40-47).
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
        let mut stretch = Stretch { base: FillingBase::new() };
        stretch.init_rational(p1, p2, p3, p4, w1, w2, w3, w4);
        stretch
    }

    /// OCCT Init(P1, P2, P3, P4) (L49-89).
    pub fn init(&mut self, p1: &[DVec3], p2: &[DVec3], p3: &[DVec3], p4: &[DVec3]) {
        assert!(
            p1.len() == p3.len() && p2.len() == p4.len(),
            "Standard_DomainError"
        );
        let npolu = p1.len();
        let npolv = p2.len();
        self.base.is_rational = false;
        let nu = (npolu - 1) as f64;
        let nv = (npolv - 1) as f64;
        self.base.poles = vec![vec![DVec3::ZERO; npolv]; npolu];
        // The boundaries are not modified
        for i in 1..=npolu {
            self.base.poles[i - 1][0] = p1[i - 1];
            self.base.poles[i - 1][npolv - 1] = p3[i - 1];
        }
        for j in 2..=npolv - 1 {
            let pv = (j - 1) as f64 / nv;
            let pv1 = 1.0 - pv;
            self.base.poles[0][j - 1] = p4[j - 1];
            self.base.poles[npolu - 1][j - 1] = p2[j - 1];
            for i in 2..=npolu - 1 {
                let pu = (i - 1) as f64 / nu;
                let pu1 = 1.0 - pu;
                let p = pv1 * p1[i - 1]
                    + pv * p3[i - 1]
                    + pu * p2[j - 1]
                    + pu1 * p4[j - 1]
                    - (pu1 * pv1 * p1[0]
                        + pu * pv1 * p2[0]
                        + pu * pv * p3[npolu - 1]
                        + pu1 * pv * p4[npolv - 1]);
                self.base.poles[i - 1][j - 1] = p;
            }
        }
    }

    /// OCCT Init(P1, P2, P3, P4, W1, W2, W3, W4) (L91-140).
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
                let w = 0.5 * (pv1 * w1[i - 1] + pv * w3[i - 1] + pu * w2[j - 1] + pu1 * w4[j - 1]);
                self.base.weights[i - 1][j - 1] = w;
            }
        }
    }
}
