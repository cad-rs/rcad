//! OCCT FEmTool_Assembly (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_Assembly.hxx` (L17-107) and
//! `FEmTool_Assembly.cxx` (L1-698).
//!
//! Container mapping: `NCollection_Sequence<NCollection_List<Handle(
//! NCollection_HArray1<double>)>>` (the constraint rows G) becomes
//! `Vec<Vec<Segment>>` with [`Segment`] mirroring the per-chunk
//! `HArray1<double>` (lower bound + storage); `NCollection_Sequence<double>`
//! C becomes `Vec<f64>`.  The `GHGt` nullable handle becomes `Option`.

use super::super::math_matrix::{Matrix, Vector};
use super::profile_matrix::ProfileMatrix;
use super::sparse_matrix::SparseMatrix;
use super::{AssemblingTable, IntArray2};

/// Mirror of one `NCollection_HArray1<double>` chunk of a constraint line
/// (G): values over the index interval [lower, upper].
#[derive(Debug, Clone)]
struct Segment {
    lower: i32,
    data: Vec<f64>,
}

impl Segment {
    /// OCCT new NCollection_HArray1<double>(Imin, Imax).
    fn new(lower: i32, upper: i32) -> Self {
        Segment {
            lower,
            data: vec![0.0; (upper - lower + 1).max(0) as usize],
        }
    }

    /// OCCT Init(0.).
    fn init(&mut self, value: f64) {
        self.data.fill(value);
    }

    #[inline]
    fn lower(&self) -> i32 {
        self.lower
    }

    #[inline]
    fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// OCCT Value(j).
    #[inline]
    fn value(&self, j: i32) -> f64 {
        self.data[(j - self.lower) as usize]
    }

    /// OCCT ChangeValue(j).
    #[inline]
    fn change_value(&mut self, j: i32) -> &mut f64 {
        &mut self.data[(j - self.lower) as usize]
    }
}

/// OCCT static MinIndex (FEmTool_Assembly.cxx L28-55).
fn min_index(table: &AssemblingTable) -> i32 {
    let mut imin;
    let diml = table.lower_row();
    let dimu = table.upper_row();
    let ell = table.lower_col();
    let elu = table.upper_col();

    let t = table.value(diml, ell);
    imin = t.value(t.lower());

    for dim in diml..=dimu {
        for el in ell..=elu {
            let t = table.value(dim, el);
            for nvar in t.lower()..=t.upper() {
                imin = imin.min(t.value(nvar));
            }
        }
    }
    imin
}

/// OCCT static MaxIndex (FEmTool_Assembly.cxx L59-86).
fn max_index(table: &AssemblingTable) -> i32 {
    let mut imax;
    let diml = table.lower_row();
    let dimu = table.upper_row();
    let ell = table.lower_col();
    let elu = table.upper_col();

    let t = table.value(diml, ell);
    imax = t.value(t.lower());

    for dim in diml..=dimu {
        for el in ell..=elu {
            let t = table.value(dim, el);
            for nvar in t.lower()..=t.upper() {
                imax = imax.max(t.value(nvar));
            }
        }
    }
    imax
}

/// OCCT FEmTool_Assembly - assemble and solve system from (one dimensional)
/// finite elements (FEmTool_Assembly.hxx L37-105).
#[derive(Debug)]
pub struct Assembly {
    /// OCCT myDepTable: NCollection_Array2<int>(1, Dependence.ColLength(),
    /// 1, Dependence.RowLength()).
    my_dep_table: IntArray2,
    /// OCCT myRefTable.
    my_ref_table: AssemblingTable,
    /// OCCT IsSolved.
    is_solved: bool,
    /// OCCT H: Handle(FEmTool_ProfileMatrix).
    h: ProfileMatrix,
    /// OCCT B: math_Vector(MinIndex(Table), MaxIndex(Table)).
    b: Vector,
    /// OCCT GHGt: nullable Handle(FEmTool_ProfileMatrix).
    gh_gt: Option<ProfileMatrix>,
    /// OCCT G: NCollection_Sequence<NCollection_List<Handle(HArray1<double>)>>
    /// (1-based sequence).
    g: Vec<Vec<Segment>>,
    /// OCCT C: NCollection_Sequence<double> (1-based sequence).
    c: Vec<f64>,
}

impl Assembly {
    /// OCCT FEmTool_Assembly::FEmTool_Assembly (cxx L90-137).
    pub fn new(dependence: &IntArray2, table: AssemblingTable) -> Self {
        // myDepTable(1, Dependence.ColLength(), 1, Dependence.RowLength());
        // myDepTable = Dependence;  (element-wise copy of the source values)
        let mut my_dep_table = IntArray2::new_init(
            1,
            dependence.col_length(),
            1,
            dependence.row_length(),
            0,
        );
        for i in 1..=dependence.col_length() {
            for j in 1..=dependence.row_length() {
                let v = dependence
                    .value(dependence.lower_row() + i - 1, dependence.lower_col() + j - 1);
                my_dep_table.set_value(i, j, v);
            }
        }
        let b_lower = min_index(&table);
        let b_upper = max_index(&table);
        let b = Vector::new(b_lower, b_upper);

        let mut first_indexes = vec![b.length(); b.length() as usize];
        // FirstIndexes.Init(B.Length());
        first_indexes.iter_mut().for_each(|v| *v = b.length());

        let i0 = 1 - b.lower();

        let diml = table.lower_row();
        let dimu = table.upper_row();
        let ell = table.lower_col();
        let elu = table.upper_col();

        for dim in diml..=dimu {
            for el in ell..=elu {
                let t = table.value(dim, el);
                let mut imin = t.value(t.lower()) + i0;

                for nvar in t.lower()..=t.upper() {
                    imin = imin.min(t.value(nvar) + i0);
                }

                for nvar in t.lower()..=t.upper() {
                    let i = t.value(nvar) + i0;
                    first_indexes[(i - 1) as usize] = first_indexes[(i - 1) as usize].min(imin);
                }
            }
        }

        let h = ProfileMatrix::new(&first_indexes);

        let mut a = Assembly {
            my_dep_table,
            my_ref_table: table,
            is_solved: false,
            h,
            b,
            gh_gt: None,
            g: Vec::new(),
            c: Vec::new(),
        };
        a.nullify_matrix();
        a.nullify_vector();
        a
    }

    /// OCCT FEmTool_Assembly::NullifyMatrix (cxx L139-143).
    pub fn nullify_matrix(&mut self) {
        self.h.init(0.0);
        self.is_solved = false;
    }

    /// OCCT FEmTool_Assembly::AddMatrix (cxx L147-179).
    pub fn add_matrix(&mut self, element: i32, dimension1: i32, dimension2: i32, mat: &Matrix) {
        assert!(
            self.my_dep_table.value(dimension1, dimension2) != 0,
            "Standard_DomainError: FEmTool_Assembly::AddMatrix"
        );

        let nvarl;
        let nvaru;
        {
            let t1 = self.my_ref_table.value(dimension1, element);
            let _t2 = self.my_ref_table.value(dimension2, element);

            nvarl = t1.lower();
            nvaru = t1.upper().min(nvarl + mat.row_number() - 1);
        }

        let i0 = 1 - self.b.lower();
        let i0m = mat.lower_row() - nvarl;
        let j0 = mat.lower_col() - nvarl;

        for i in nvarl..=nvaru {
            let big_i = self.my_ref_table.value(dimension1, element).value(i) + i0;
            let ii = i0m + i;
            for j in nvarl..=i {
                let big_j = self.my_ref_table.value(dimension2, element).value(j) + i0;
                let v = *self.h.change_value(big_i, big_j) + mat.get(ii, j0 + j);
                *self.h.change_value(big_i, big_j) = v;
            }
        }

        self.is_solved = false;
    }

    /// OCCT FEmTool_Assembly::NullifyVector (cxx L183-186).
    pub fn nullify_vector(&mut self) {
        for r in 1..=self.b.length() {
            self.b.set(self.b.lower() + r - 1, 0.0);
        }
    }

    /// OCCT FEmTool_Assembly::AddVector (cxx L190-203).
    pub fn add_vector(&mut self, element: i32, dimension: i32, vec: &Vector) {
        let t = self.my_ref_table.value(dimension, element);
        let nvarl = t.lower();
        let nvaru = t.upper().min(nvarl + vec.length() - 1);
        let i0 = vec.lower() - nvarl;

        for i in nvarl..=nvaru {
            let idx = t.value(i);
            let v = self.b.get(idx) + vec.get(i0 + i);
            self.b.set(idx, v);
        }
    }

    /// OCCT FEmTool_Assembly::Solve (cxx L207-395).
    pub fn solve(&mut self) -> bool {
        self.is_solved = self.h.decompose();

        if !self.g.is_empty() && self.is_solved {
            let nb_glob_var = self.nb_glob_var();
            let g_length = self.g.len() as i32;

            // calculating H-1 Gt
            let mut gh_gt_needs_build = true;
            if let Some(ghgt) = &self.gh_gt {
                if ghgt.row_number() == g_length {
                    gh_gt_needs_build = false;
                }
            }

            if gh_gt_needs_build {
                let mut first_indexes = vec![0i32; g_length as usize];
                //-----------------------------------------------------------------
                let mut h1 = super::IntArray2::new_init(1, nb_glob_var, 1, nb_glob_var, 1);
                let mut block_beg = 1;
                let mut i = 2;
                while i <= nb_glob_var {
                    let block_end = i - 1;
                    if !self.h.is_in_profile(i, block_end) {
                        // Maybe, begin of block
                        let mut block = true;
                        let mut j = i + 1;
                        while j <= nb_glob_var {
                            if self.h.is_in_profile(j, block_end) {
                                block = false;
                                break;
                            }
                            j += 1;
                        }
                        if block {
                            for j in i..=nb_glob_var {
                                for k in block_beg..=block_end {
                                    h1.set_value(j, k, 0);
                                    h1.set_value(k, j, 0);
                                }
                            }
                            block_beg = block_end + 1;
                        } else {
                            i = j;
                        }
                    }
                    i += 1;
                }

                for a_i in 1..=g_length {
                    let gi = &self.g[(a_i - 1) as usize];
                    'outer_j: for j in 1..=a_i {
                        let gj = &self.g[(j - 1) as usize];
                        let mut zero = true;
                        for a in gi.iter() {
                            for k in a.lower()..=a.upper() {
                                for b in gj.iter() {
                                    for l in b.lower()..=b.upper() {
                                        if h1.value(k, l) != 0 {
                                            zero = false;
                                            break;
                                        }
                                    }
                                    if !zero {
                                        break;
                                    }
                                }
                                if !zero {
                                    break;
                                }
                            }
                            if !zero {
                                break;
                            }
                        }
                        if !zero {
                            first_indexes[(a_i - 1) as usize] = j;
                            break 'outer_j;
                        }
                    }
                }
                //-----------------------------------------------------------------------
                self.gh_gt = Some(ProfileMatrix::new(&first_indexes));
            }

            let ghgt = self.gh_gt.as_mut().unwrap();
            ghgt.init(0.0);

            let mut gi = Vector::new(self.b.lower(), self.b.upper());
            let mut qi = Vector::new(self.b.lower(), self.b.upper());

            for a_i in 1..=g_length {
                let l = &self.g[(a_i - 1) as usize];
                // gi.Init(0.);
                for r in gi.lower()..=gi.upper() {
                    gi.set(r, 0.0);
                }
                // preparing i-th line of G (or column of Gt)
                for a in l.iter() {
                    for j in a.lower()..=a.upper() {
                        gi.set(j, a.value(j)); // gi - full line of G
                    }
                }
                // H*qi = gi, qi is column of H-1 Gt
                self.h.solve(&gi, &mut qi);

                // Calculation of product M = G H-1 Gt
                // for each i all elements of i-th column of M are calculated
                // for k >= i
                for k in a_i..=g_length {
                    if ghgt.is_in_profile(k, a_i) {
                        let mut m = 0.0; // m = M(k,i)

                        let al = &self.g[(k - 1) as usize];

                        for a in al.iter() {
                            for j in a.lower()..=a.upper() {
                                // scalar product of k-th line of G and i-th
                                // column of H-1 Gt
                                m += qi.get(j) * a.value(j);
                            }
                        }

                        *ghgt.change_value(k, a_i) = m;
                    }
                }
            }

            self.is_solved = ghgt.decompose();
        }

        self.is_solved
    }

    /// OCCT FEmTool_Assembly::Solution (cxx L399-464).
    pub fn solution(&self, solution: &mut Vector) {
        assert!(
            self.is_solved,
            "StdFail_NotDone: FEmTool_Assembly::Solution"
        );

        if self.g.is_empty() {
            self.h.solve(&self.b, solution);
        } else {
            let mut v1 = Vector::new(self.b.lower(), self.b.upper());
            self.h.solve(&self.b, &mut v1);

            let g_length = self.g.len() as i32;
            let mut l = Vector::new(1, g_length);
            let mut v2 = Vector::new(1, g_length);

            for a_i in 1..=g_length {
                let lg = &self.g[(a_i - 1) as usize];
                let mut m = 0.0;

                for a in lg.iter() {
                    for j in a.lower()..=a.upper() {
                        m += v1.get(j) * a.value(j); // scalar product, G v1
                    }
                }

                v2.set(a_i, m - self.c[(a_i - 1) as usize]);
            }

            let ghgt = self.gh_gt.as_ref().unwrap();
            ghgt.solve(&v2, &mut l); // Solving M*l = v2

            // v1 = B; Calculation v1 = B-Gt*l
            // v1(j) = B(j) - Gt(j,i)*l(i) = B(j) - G(i,j)*l(i)
            let mut v1 = self.b.clone();
            for a_i in 1..=g_length {
                let lg = &self.g[(a_i - 1) as usize];

                for a in lg.iter() {
                    for j in a.lower()..=a.upper() {
                        let v = v1.get(j) - l.get(a_i) * a.value(j);
                        v1.set(j, v);
                    }
                }
            }

            self.h.solve(&v1, solution);
        }
    }

    /// OCCT FEmTool_Assembly::NbGlobVar (cxx L466-470).
    pub fn nb_glob_var(&self) -> i32 {
        self.b.length()
    }

    /// OCCT FEmTool_Assembly::AssemblyTable (hxx L84-88).
    pub fn assembly_table(&self) -> &AssemblingTable {
        &self.my_ref_table
    }

    /// OCCT FEmTool_Assembly::ResetConstraint (cxx L478-482) - delete all
    /// constraints.
    pub fn reset_constraint(&mut self) {
        self.g.clear();
        self.c.clear();
    }

    /// OCCT FEmTool_Assembly::NullifyConstraint (cxx L484-499).
    pub fn nullify_constraint(&mut self) {
        for i in 1..=self.g.len() as i32 {
            self.c[(i - 1) as usize] = 0.0;

            for seg in self.g[(i - 1) as usize].iter_mut() {
                seg.init(0.0);
            }
        }
    }

    /// OCCT FEmTool_Assembly::AddConstraint (cxx L503-698).
    pub fn add_constraint(
        &mut self,
        indexof_constraint: i32,
        element: i32,
        dimension: i32,
        linear_form: &Vector,
        value: f64,
    ) {
        while (self.g.len() as i32) < indexof_constraint {
            // Add new lines in G
            self.g.push(Vec::new());
            self.c.push(0.0);
        }

        let indexes = self.my_ref_table.value(dimension, element).clone();
        let mut imax = 0;
        let mut imin = self.nb_glob_var();

        for i in indexes.lower()..=indexes.upper() {
            imin = imin.min(indexes.value(i));
            imax = imax.max(indexes.value(i));
        }

        let l = &mut self.g[(indexof_constraint - 1) as usize];

        // OCCT `Coeff` is a handle to the list element that receives the
        // additions; `coeff_idx` tracks the position of that element.
        let coeff_idx: usize;

        if l.is_empty() {
            let mut coeff = Segment::new(imin, imax);
            coeff.init(0.0);
            l.push(coeff);
            coeff_idx = l.len() - 1;
        } else {
            let mut s1 = 0;
            let mut s2 = 0;
            let mut aux1: Option<&Segment> = None;
            let mut aux2: Option<&Segment> = None;
            let mut i = 1i32;
            for seg in l.iter() {
                if imin >= seg.lower() {
                    s1 = i;
                    aux1 = Some(seg);
                    if imax <= seg.upper() {
                        s2 = s1;
                        break;
                    }
                }

                if imax <= seg.upper() {
                    s2 = i;
                    aux2 = Some(seg);
                }
                i += 1;
            }

            if s1 != s2 {
                let aux1 = aux1.expect("AddConstraint: Aux1");
                let aux2 = aux2.expect("AddConstraint: Aux2");
                if s1 == 0 {
                    if imax < aux2.lower() {
                        // inserting before first segment
                        let mut coeff = Segment::new(imin, imax);
                        coeff.init(0.0);
                        l.insert(0, coeff);
                        coeff_idx = 0;
                    } else {
                        // merge new and first segment
                        let mut coeff = Segment::new(imin, aux2.upper());
                        for i in imin..=aux2.lower() - 1 {
                            *coeff.change_value(i) = 0.0;
                        }
                        for i in aux2.lower()..=aux2.upper() {
                            *coeff.change_value(i) = aux2.value(i);
                        }
                        l[0] = coeff;
                        coeff_idx = 0;
                    }
                } else if s2 == 0 {
                    if imin > aux1.upper() {
                        // append new
                        let mut coeff = Segment::new(imin, imax);
                        coeff.init(0.0);
                        l.push(coeff);
                        coeff_idx = l.len() - 1;
                    } else {
                        // merge new and last segment
                        let mut coeff = Segment::new(aux1.lower(), imax);
                        for i in aux1.lower()..=aux1.upper() {
                            *coeff.change_value(i) = aux1.value(i);
                        }
                        for i in aux1.upper() + 1..=imax {
                            *coeff.change_value(i) = 0.0;
                        }
                        let last = l.len() - 1;
                        l[last] = coeff;
                        coeff_idx = last;
                    }
                } else if imin <= aux1.upper() && imax < aux2.lower() {
                    // merge s1 and new
                    let mut coeff = Segment::new(aux1.lower(), imax);
                    for i in aux1.lower()..=aux1.upper() {
                        *coeff.change_value(i) = aux1.value(i);
                    }
                    for i in aux1.upper() + 1..=imax {
                        *coeff.change_value(i) = 0.0;
                    }
                    l[(s1 - 1) as usize] = coeff;
                    coeff_idx = (s1 - 1) as usize;
                } else if imin > aux1.upper() && imax >= aux2.lower() {
                    // merge new and first segment
                    let mut coeff = Segment::new(imin, aux2.upper());
                    for i in imin..=aux2.lower() - 1 {
                        *coeff.change_value(i) = 0.0;
                    }
                    for i in aux2.lower()..=aux2.upper() {
                        *coeff.change_value(i) = aux2.value(i);
                    }
                    l[(s2 - 1) as usize] = coeff;
                    coeff_idx = (s2 - 1) as usize;
                } else if imin > aux1.upper() && imax < aux2.lower() {
                    // inserting new between s1 and s2
                    let mut coeff = Segment::new(imin, imax);
                    coeff.init(0.0);
                    l.insert(s1 as usize, coeff);
                    coeff_idx = s1 as usize;
                } else {
                    // merge s1, new, s2 and remove s2
                    let mut coeff = Segment::new(aux1.lower(), aux2.upper());
                    for i in aux1.lower()..=aux1.upper() {
                        *coeff.change_value(i) = aux1.value(i);
                    }
                    for i in aux1.upper() + 1..=aux2.lower() - 1 {
                        *coeff.change_value(i) = 0.0;
                    }
                    for i in aux2.lower()..=aux2.upper() {
                        *coeff.change_value(i) = aux2.value(i);
                    }
                    l[(s1 - 1) as usize] = coeff;
                    l.remove(s1 as usize);
                    coeff_idx = (s1 - 1) as usize;
                }
            } else {
                // The constraint falls inside the s1 segment itself: OCCT
                // keeps `Coeff` pointing at that (already existing) element.
                coeff_idx = (s1 - 1) as usize;
            }
        }

        // adding
        let mut j = linear_form.lower();
        for i in indexes.lower()..=indexes.upper() {
            *l[coeff_idx].change_value(indexes.value(i)) += linear_form.get(j);
            j += 1;
        }

        self.c[(indexof_constraint - 1) as usize] += value;
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::math_matrix::Matrix;
    use super::*;

    /// FEmTool_Assembly solve without constraints, hand-solved.
    /// Two elements, one dimension, global variables {1,2} and {2,3}:
    ///     H = [[2,1,0],[1,4,1],[0,1,2]]   (element blocks [[2,1],[1,2]])
    ///     B = [1, 5, 4]                    (AddVector [1,2] then [3,4])
    /// Gaussian elimination by hand gives x = [1/12, 5/6, 19/12]
    /// (det H = 12, all leading minors positive).
    #[test]
    fn assembly_solve_hand() {
        // Dependence: dimension 1 depends on itself.
        let dependence = IntArray2::new_init(1, 1, 1, 1, 1);

        // Assembly table: 1 dimension x 2 elements, local arrays (1,2).
        let mut table = AssemblingTable::new(1, 1, 1, 2, 1, 2);
        table.change_value(1, 1).set_value(1, 1);
        table.change_value(1, 1).set_value(2, 2);
        table.change_value(1, 2).set_value(1, 2);
        table.change_value(1, 2).set_value(2, 3);

        let mut assembly = Assembly::new(&dependence, table);
        assert_eq!(assembly.nb_glob_var(), 3);

        let mut mat = Matrix::new(1, 2, 1, 2);
        mat.set(1, 1, 2.0);
        mat.set(1, 2, 1.0);
        mat.set(2, 1, 1.0);
        mat.set(2, 2, 2.0);
        assembly.add_matrix(1, 1, 1, &mat);
        assembly.add_matrix(2, 1, 1, &mat);

        let mut vec = Vector::new(1, 2);
        vec.set(1, 1.0);
        vec.set(2, 2.0);
        assembly.add_vector(1, 1, &vec);
        vec.set(1, 3.0);
        vec.set(2, 4.0);
        assembly.add_vector(2, 1, &vec);

        assert!(assembly.solve());

        let mut solution = Vector::new(1, 3);
        assembly.solution(&mut solution);

        assert!((solution.get(1) - 1.0 / 12.0).abs() < 1e-13, "x1 = {}", solution.get(1));
        assert!((solution.get(2) - 5.0 / 6.0).abs() < 1e-13, "x2 = {}", solution.get(2));
        assert!((solution.get(3) - 19.0 / 12.0).abs() < 1e-13, "x3 = {}", solution.get(3));
    }
}
