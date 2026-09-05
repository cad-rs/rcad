// OCCT HLRBRep_Data.cxx L77-492 (the file-static TableauRejection class).
//
// The OCCT member arrays are malloc/free row tables (UV[i] grown by SIZEUV
// slots, IndUV rows with the -1 sentinel, TabBit 32-bit words); rcad maps
// each malloc'd row to a Vec (the NCollection Array1 precedent) keeping the
// OCCT statement order, the -1 sentinel protocol and the Mask32 bit
// arithmetic bit-for-bit.  The OCCT_DEBUG statistics members
// (StNbLect/StNbEcr/StNbMax/StNbMoy/StNbMoyNonNul) are an #ifdef block and
// are not translated.
//
// The Data.my_reject field keeps the OCCT `new TableauRejection()` /
// `delete` shape through a unique-ownership raw pointer (the ctor and
// Destroy call sites live in [update]).

use crate::geomalgo::int_res2d::IntersectionPoint;

use super::MASK32;
use super::SIZEUV;

/// OCCT Standard_Real RealLast() (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;

/// OCCT class TableauRejection (cxx L77-492).
pub struct TableauRejection {
    /// OCCT `double** UV` — UV[i][j] holds the parameter (U on Ci) of the
    /// intersection of Ci with C(IndUV[j]).
    pub uv: Vec<Vec<f64>>,
    /// OCCT `int** IndUV` — IndUV[i][j] = J0 -> intersection between i and
    /// J0; the -1 entry terminates the filled prefix of each row.
    pub ind_uv: Vec<Vec<i32>>,
    /// OCCT `int* nbUV` — the allocated row width for line i.
    pub nb_uv: Vec<i32>,
    /// OCCT `int N`.
    pub n: i32,

    /// OCCT `long unsigned** TabBit`.
    pub tab_bit: Vec<Vec<u64>>,
    /// OCCT `int nTabBit`.
    pub n_tab_bit: i32,
}

impl Default for TableauRejection {
    fn default() -> Self {
        Self::new()
    }
}

impl TableauRejection {
    /// OCCT TableauRejection() (cxx L100-111).
    pub fn new() -> Self {
        TableauRejection {
            n: 0,
            n_tab_bit: 0,
            uv: Vec::new(),
            nb_uv: Vec::new(),
            ind_uv: Vec::new(),
            tab_bit: Vec::new(),
        }
    }

    /// OCCT SetDim(const int n) (cxx L114-146).
    pub fn set_dim(&mut self, n: i32) {
        if !self.uv.is_empty() {
            self.destroy();
        }
        self.n = n;
        // OCCT L127-129: UV/IndUV/nbUV = malloc(N * ...) — the rows are
        // allocated in the loops below.
        self.uv = (0..self.n).map(|_| Vec::new()).collect();
        self.ind_uv = (0..self.n).map(|_| Vec::new()).collect();
        self.nb_uv = vec![0; self.n as usize];
        //    for(int i=0;i<N;i++) {
        let mut i: i32;
        i = 0;
        while i < self.n {
            // OCCT L134: UV[i] = malloc(SIZEUV * sizeof(double)) — the OCCT
            // block is uninitialized; the neutral default is 0.
            self.uv[i as usize] = vec![0.0; SIZEUV];
            i += 1;
        }
        i = 0;
        while i < self.n {
            self.ind_uv[i as usize] = vec![0; SIZEUV];
            let mut k: i32 = 0;
            while k < SIZEUV as i32 {
                self.ind_uv[i as usize][k as usize] = -1;
                k += 1;
            }
            self.nb_uv[i as usize] = SIZEUV as i32;
            i += 1;
        }
        self.init_tab_bit(n);
    }

    /// OCCT ~TableauRejection() (cxx L149-153).
    pub fn drop_rejection(&mut self) {
        self.destroy();
    }

    /// OCCT Destroy() (cxx L156-241) — the #ifdef OCCT_DEBUG statistics
    /// block is not translated.
    pub fn destroy(&mut self) {
        if self.n != 0 {
            self.reset_tab_bit(self.n);
            //      for(int i=0;i<N;i++) {
            let mut i: i32;
            i = 0;
            while i < self.n {
                if !self.ind_uv[i as usize].is_empty() {
                    // OCCT L201-202: free(IndUV[i]); IndUV[i] = nullptr;
                    self.ind_uv[i as usize] = Vec::new();
                }
                i += 1;
            }
            i = 0;
            while i < self.n {
                if !self.uv[i as usize].is_empty() {
                    // OCCT L213-214: free(UV[i]); UV[i] = nullptr;
                    self.uv[i as usize] = Vec::new();
                }
                i += 1;
            }

            if !self.nb_uv.is_empty() {
                self.nb_uv = Vec::new();
            }
            if !self.ind_uv.is_empty() {
                self.ind_uv = Vec::new();
            }
            if !self.uv.is_empty() {
                self.uv = Vec::new();
            }
            self.n = 0;
        }
    }

    /// OCCT Set(int i0, int j0, const double u) (cxx L244-315).
    pub fn set(&mut self, i0: i32, j0: i32, u: f64) {
        let i0 = (i0 - 1) as usize; // OCCT: i0--
        let j0 = (j0 - 1) as usize; // OCCT: j0--
        let mut k: i32 = -1;
        //    for(int i=0; k==-1 && i<nbUV[i0]; i++) {
        let mut i: usize;
        i = 0;
        while k == -1 && i < self.nb_uv[i0] as usize {
            if self.ind_uv[i0][i] == -1 {
                k = i as i32;
            }
            i += 1;
        }
        if k == -1 {
            //-- on agrandit le tableau
            //--
            //-- declaration de la Nv ligne de taille : ancienne taille + SIZEUV
            //--
            let old_nb = self.nb_uv[i0] as usize;
            // OCCT L269-270: the new rows are malloc'd uninitialized; the
            // neutral default is 0 (every slot is overwritten below).
            let mut nv_ligne_uv: Vec<f64> = vec![0.0; old_nb + SIZEUV];
            let mut nv_ligne_ind: Vec<i32> = vec![0; old_nb + SIZEUV];
            //--
            //-- Recopie des anciennes valeurs ds la nouvelle ligne
            //--
            let mut ii: usize = 0;
            while ii < old_nb {
                nv_ligne_uv[ii] = self.uv[i0][ii];
                nv_ligne_ind[ii] = self.ind_uv[i0][ii];
                ii += 1;
            }

            //-- mise a jour de la nouvelle dimension ; free des anciennes
            //   lignes et affectation
            k = self.nb_uv[i0];
            self.nb_uv[i0] += SIZEUV as i32;
            self.uv[i0] = nv_ligne_uv; // OCCT L283-285: free(UV[i0]); UV[i0] = NvLigneUV;
            self.ind_uv[i0] = nv_ligne_ind; // OCCT L284/286: free(IndUV[i0]); IndUV[i0] = NvLigneInd;
            let mut kk = k;
            while kk < self.nb_uv[i0] {
                self.ind_uv[i0][kk as usize] = -1;
                kk += 1;
            }
        }
        self.ind_uv[i0][k as usize] = j0 as i32;
        self.uv[i0][k as usize] = u;

        //-- tri par ordre decroissant
        let mut tri_ok: bool;
        loop {
            tri_ok = true;
            let mut im1: usize = 0;
            i = 1;
            while self.ind_uv[i0][i] != -1 && i < self.nb_uv[i0] as usize {
                if self.ind_uv[i0][i] > self.ind_uv[i0][im1] {
                    tri_ok = false;
                    k = self.ind_uv[i0][i];
                    self.ind_uv[i0][i] = self.ind_uv[i0][im1];
                    self.ind_uv[i0][im1] = k;
                    let t = self.uv[i0][i];
                    self.uv[i0][i] = self.uv[i0][im1];
                    self.uv[i0][im1] = t;
                }
                i += 1;
                im1 += 1;
            }
            if tri_ok {
                break;
            }
        }
    }

    /// OCCT Get(int i0, int j0) (cxx L318-377) — the descending binary
    /// search; RealLast() when the pair is not found.
    pub fn get(&mut self, i0: i32, j0: i32) -> f64 {
        let i0 = (i0 - 1) as usize; // OCCT: i0--
        let j0 = (j0 - 1) as i32; // OCCT: j0--

        //-- ordre decroissant
        let mut a: i32 = 0;
        let mut b: i32 = self.nb_uv[i0] - 1;
        let mut ab: i32;
        if self.ind_uv[i0][a as usize] == -1 {
            return REAL_LAST;
        }
        if self.ind_uv[i0][a as usize] == j0 {
            return self.uv[i0][a as usize];
        }
        if self.ind_uv[i0][b as usize] == j0 {
            return self.uv[i0][b as usize];
        }
        while (self.ind_uv[i0][a as usize] > j0) && (self.ind_uv[i0][b as usize] < j0) {
            ab = (a + b) >> 1;
            if self.ind_uv[i0][ab as usize] < j0 {
                if b == ab {
                    return REAL_LAST;
                } else {
                    b = ab;
                }
            } else if self.ind_uv[i0][ab as usize] > j0 {
                if a == ab {
                    return REAL_LAST;
                } else {
                    a = ab;
                }
            } else {
                return self.uv[i0][ab as usize];
            }
        }

        REAL_LAST
    }

    /// OCCT ResetTabBit(const int nbedgs) (cxx L380-397).
    pub fn reset_tab_bit(&mut self, nbedgs: i32) {
        if !self.tab_bit.is_empty() {
            let mut i: i32 = 0;
            while i < nbedgs {
                if !self.tab_bit[i as usize].is_empty() {
                    // OCCT L389-390: free(TabBit[i]); TabBit[i] = nullptr;
                    self.tab_bit[i as usize] = Vec::new();
                }
                i += 1;
            }
            self.tab_bit = Vec::new();
            self.n_tab_bit = 0;
        }
    }

    /// OCCT InitTabBit(const int nbedgs) (cxx L400-419).
    pub fn init_tab_bit(&mut self, nbedgs: i32) {
        if !self.tab_bit.is_empty() && self.n_tab_bit != 0 {
            self.reset_tab_bit(self.n_tab_bit);
        }
        self.tab_bit = (0..nbedgs).map(|_| Vec::new()).collect();
        self.n_tab_bit = nbedgs;
        let n: usize = 1 + ((nbedgs >> 5) as usize);

        let mut i: i32 = 0;
        while i < nbedgs {
            // OCCT L413-417: TabBit[i] = malloc(n * sizeof(long unsigned))
            // zeroed by the loop; the calloc-equivalent Vec![0] reproduces it.
            self.tab_bit[i as usize] = vec![0u64; n];
            i += 1;
        }
    }

    /// OCCT SetNoIntersection(int i0, int i1) (cxx L422-436).
    pub fn set_no_intersection(&mut self, i0: i32, i1: i32) {
        let mut i0 = i0 - 1; // OCCT: i0--
        let mut i1 = i1 - 1; // OCCT: i1--
        if i0 > i1 {
            let t = i0;
            i0 = i1;
            i1 = t;
        }
        let c = (i1 >> 5) as usize;
        let o = (i1 & 31) as usize;
        self.tab_bit[i0 as usize][c] |= MASK32[o];
    }

    /// OCCT NoIntersection(int i0, int i1) (cxx L439-459).
    pub fn no_intersection(&mut self, i0: i32, i1: i32) -> bool {
        let mut i0 = i0 - 1; // OCCT: i0--
        let mut i1 = i1 - 1; // OCCT: i1--
        if i0 > i1 {
            let t = i0;
            i0 = i1;
            i1 = t;
        }
        let c = (i1 >> 5) as usize;
        let o = (i1 & 31) as usize;
        if self.tab_bit[i0 as usize][c] & MASK32[o] != 0 {
            return true;
        }
        false
    }

    /// OCCT SetIntersection(int i0, int i1, const IntRes2d_IntersectionPoint& IP)
    /// (cxx L462-477).
    pub fn set_intersection(&mut self, i0: i32, i1: i32, ip: &IntersectionPoint) {
        let t1 = ip.transition_of_first();
        let t2 = ip.transition_of_second();
        use super::super::super::super::geomalgo::int_res2d::{Position, TypeTrans};
        if t1.position_on_curve() == Position::Middle {
            if t2.position_on_curve() == Position::Middle {
                if t1.transition_type() == TypeTrans::In || t1.transition_type() == TypeTrans::Out {
                    self.set(i0, i1, ip.param_on_first());
                    self.set(i1, i0, ip.param_on_second());
                }
            }
        }
    }

    /// OCCT GetSingleIntersection(int i0, int i1, double& u, double& v)
    /// (cxx L480-491).
    pub fn get_single_intersection(&mut self, i0: i32, i1: i32, u: &mut f64, v: &mut f64) {
        *u = self.get(i0, i1);
        if *u != REAL_LAST {
            *v = self.get(i1, i0);
        } else {
            *v = REAL_LAST;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT anchor: SetDim allocates N rows of SIZEUV slots, all IndUV
    /// entries -1 (cxx L114-146); Get on an untouched row returns RealLast.
    #[test]
    fn tableau_rejection_set_dim_empty_get() {
        let mut tr = TableauRejection::new();
        tr.set_dim(3);
        assert_eq!(tr.n, 3);
        assert_eq!(tr.nb_uv.len(), 3);
        for i in 0..3usize {
            assert_eq!(tr.nb_uv[i], SIZEUV as i32);
            assert_eq!(tr.ind_uv[i][0], -1);
        }
        assert_eq!(tr.get(1, 2), f64::MAX, "RealLast on an empty row");
        assert_eq!(tr.tab_bit.len(), 3);
        // n = 1 + (nbedgs >> 5) = 1 for nbedgs <= 32.
        assert_eq!(tr.tab_bit[0].len(), 1);
        assert_eq!(tr.tab_bit[0][0], 0);
    }

    /// OCCT anchor: Set/Get round-trip (cxx L244-377) — the row is kept in
    /// descending IndUV order by the sort loop, and Get performs the binary
    /// search against it.
    #[test]
    fn tableau_rejection_set_get_round_trip() {
        let mut tr = TableauRejection::new();
        tr.set_dim(4);
        tr.set(1, 2, 0.5);
        assert_eq!(tr.get(1, 2), 0.5);
        assert_eq!(tr.get(1, 3), f64::MAX);
        assert_eq!(tr.get(2, 1), f64::MAX, "the table is not symmetric");

        tr.set(1, 3, 0.25);
        tr.set(1, 5, 0.75);
        // The descending order on the row: IndUV 5 (u=0.75), 3 (0.25), 2 (0.5).
        assert_eq!(tr.get(1, 5), 0.75);
        assert_eq!(tr.get(1, 3), 0.25);
        assert_eq!(tr.get(1, 2), 0.5);
        assert_eq!(tr.get(1, 4), f64::MAX);
        // The -1 sentinel is preserved after the filled prefix.
        assert!(tr.ind_uv[0].iter().any(|&v| v == -1));
    }

    /// OCCT anchor: SetNoIntersection / NoIntersection (cxx L422-459) — the
    /// Mask32 bit at (i0, i1) with the argument-order swap is bit-symmetric.
    #[test]
    fn tableau_rejection_no_intersection_bits() {
        let mut tr = TableauRejection::new();
        tr.set_dim(70);
        assert!(!tr.no_intersection(3, 67));
        tr.set_no_intersection(3, 67);
        assert!(tr.no_intersection(3, 67));
        // Symmetric query through the i0 > i1 swap branch.
        assert!(tr.no_intersection(67, 3));
        // Another pair in a second 32-bit word of the same row.
        tr.set_no_intersection(1, 2);
        assert!(tr.no_intersection(1, 2));
        assert!(tr.no_intersection(2, 1));
        assert!(!tr.no_intersection(2, 3));
        // init_tab_bit on a new SetDim resets the bits (cxx L119-122).
        tr.set_dim(70);
        assert!(!tr.no_intersection(3, 67));
    }

    /// OCCT anchor: SetIntersection (cxx L462-477) — only a Middle/Middle
    /// intersection with an IN or OUT transition on the first curve fills
    /// both directions; GetSingleIntersection (cxx L480-491) reads the
    /// mirrored parameters.
    #[test]
    fn tableau_rejection_set_intersection_and_single() {
        use crate::geomalgo::int_res2d::{Position, Transition, TypeTrans};

        let mut tr = TableauRejection::new();
        tr.set_dim(4);

        // A Middle/Middle IN point on the first curve.
        let mut ip = IntersectionPoint::empty();
        ip.set_values(
            glam::DVec2::new(0.25, 0.0),
            0.3,
            0.7,
            Transition::in_out(false, Position::Middle, TypeTrans::In),
            Transition::in_out(false, Position::Middle, TypeTrans::Out),
            false,
        );
        tr.set_intersection(1, 2, &ip);
        let mut u = 0.0;
        let mut v = 0.0;
        tr.get_single_intersection(1, 2, &mut u, &mut v);
        assert_eq!(u, 0.3);
        assert_eq!(v, 0.7);

        // A Head-position point is rejected by SetIntersection.
        let mut ip2 = IntersectionPoint::empty();
        ip2.set_values(
            glam::DVec2::ZERO,
            0.1,
            0.2,
            Transition::in_out(false, Position::Head, TypeTrans::In),
            Transition::in_out(false, Position::Middle, TypeTrans::Out),
            false,
        );
        tr.set_intersection(2, 3, &ip2);
        let mut u2 = 0.0;
        let mut v2 = 0.0;
        tr.get_single_intersection(2, 3, &mut u2, &mut v2);
        assert_eq!(u2, f64::MAX);
        assert_eq!(v2, f64::MAX);
    }

    /// OCCT anchor: Destroy (cxx L156-241) frees everything and Get after a
    /// re-SetDim behaves like a fresh table.
    #[test]
    fn tableau_rejection_destroy_resets() {
        let mut tr = TableauRejection::new();
        tr.set_dim(5);
        tr.set(1, 1, 1.0);
        tr.destroy();
        assert_eq!(tr.n, 0);
        assert!(tr.uv.is_empty() && tr.ind_uv.is_empty() && tr.nb_uv.is_empty());
        assert!(tr.tab_bit.is_empty() && tr.n_tab_bit == 0);
        // The destructor body runs once more without effect (OCCT dtor).
        tr.drop_rejection();
        assert_eq!(tr.n, 0);
    }
}
