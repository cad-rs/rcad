//! OCCT GeomFill_Curved (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Curved.cxx (whole file L25-228).
//!
//! Mapping: `NCollection_Array1<gp_Pnt>` -> `&[DVec3]`, `NCollection_Array2`
//! -> `Vec<Vec<DVec3>>` indexed `[i - 1][j - 1]` inside `1..=n` loops.
//! gp_Pnt::Translated(gp_Vec(P2(1), P2(j))) -> `p2[j-1] - p2[0]` addition.

use glam::DVec3;

use super::filling::FillingBase;

/// OCCT GeomFill_Curved.
#[derive(Debug, Clone, Default)]
pub struct Curved {
    pub(crate) base: FillingBase,
}

impl Curved {
    /// OCCT GeomFill_Curved(P1, P2, P3, P4) (L29-35).
    pub fn new(p1: &[DVec3], p2: &[DVec3], p3: &[DVec3], p4: &[DVec3]) -> Self {
        let mut curved = Curved { base: FillingBase::new() };
        curved.init(p1, p2, p3, p4);
        curved
    }

    /// OCCT GeomFill_Curved(P1, P2, P3, P4, W1, W2, W3, W4) (L37-44).
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
        let mut curved = Curved { base: FillingBase::new() };
        curved.init_rational4(p1, p2, p3, p4, w1, w2, w3, w4);
        curved
    }

    /// OCCT GeomFill_Curved(P1, P2) (L46-52).
    pub fn new_two(p1: &[DVec3], p2: &[DVec3]) -> Self {
        let mut curved = Curved { base: FillingBase::new() };
        curved.init_two(p1, p2);
        curved
    }

    /// OCCT GeomFill_Curved(P1, P2, W1, W2) (L54-61).
    pub fn new_two_rational(p1: &[DVec3], p2: &[DVec3], w1: &[f64], w2: &[f64]) -> Self {
        let mut curved = Curved { base: FillingBase::new() };
        curved.init_two_rational(p1, p2, w1, w2);
        curved
    }

    /// OCCT Init(P1, P2, P3, P4) (L63-110).
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
            let mut pv = (j - 1) as f64 / nv;
            let mut pv1 = 1.0 - pv;
            pv /= 2.0;
            pv1 /= 2.0;
            self.base.poles[0][j - 1] = p4[j - 1];
            self.base.poles[npolu - 1][j - 1] = p2[j - 1];
            for i in 2..=npolu - 1 {
                let mut pu = (i - 1) as f64 / nu;
                let mut pu1 = 1.0 - pu;
                pu /= 2.0;
                pu1 /= 2.0;
                let p = pv1 * p1[i - 1] + pv * p3[i - 1] + pu * p2[j - 1] + pu1 * p4[j - 1];
                self.base.poles[i - 1][j - 1] = p;
            }
        }
    }

    /// OCCT Init(P1, P2, P3, P4, W1, W2, W3, W4) (L112-160).
    pub fn init_rational4(
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
            let mut pv = (j - 1) as f64 / nv;
            let mut pv1 = 1.0 - pv;
            pv /= 2.0;
            pv1 /= 2.0;
            self.base.weights[0][j - 1] = w4[j - 1];
            self.base.weights[npolu - 1][j - 1] = w2[j - 1];
            for i in 2..=npolu - 1 {
                let mut pu = (i - 1) as f64 / nu;
                let mut pu1 = 1.0 - pu;
                pu /= 2.0;
                pu1 /= 2.0;
                let w = pv1 * w1[i - 1] + pv * w3[i - 1] + pu * w2[j - 1] + pu1 * w4[j - 1];
                self.base.weights[i - 1][j - 1] = w;
            }
        }
    }

    /// OCCT Init(P1, P2) (L162-186).
    pub fn init_two(&mut self, p1: &[DVec3], p2: &[DVec3]) {
        let npolu = p1.len();
        let npolv = p2.len();
        self.base.is_rational = false;
        self.base.poles = vec![vec![DVec3::ZERO; npolv]; npolu];
        for j in 1..=npolv {
            // gp_Vec Tra(P2(1), P2(j))  ==  P2(j) - P2(1).
            let tra = p2[j - 1] - p2[0];
            for i in 1..=npolu {
                self.base.poles[i - 1][j - 1] = p1[i - 1] + tra;
            }
        }
    }

    /// OCCT Init(P1, P2, W1, W2) (L188-228).
    pub fn init_two_rational(&mut self, p1: &[DVec3], p2: &[DVec3], w1: &[f64], w2: &[f64]) {
        self.init_two(p1, p2);
        self.base.is_rational = true;
        // Initialisation des poids.
        let npolu = w1.len();
        let npolv = w2.len();
        self.base.weights = vec![vec![0.0f64; npolv]; npolu];
        for j in 1..=npolv {
            let factor = w2[j - 1] / w1[0];
            for i in 1..=npolu {
                self.base.weights[i - 1][j - 1] = w1[i - 1] * factor;
            }
        }
    }
}
