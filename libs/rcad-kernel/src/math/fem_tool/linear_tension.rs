//! OCCT FEmTool_LinearTension (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_LinearTension.hxx` (L17-53) and
//! `FEmTool_LinearTension.cxx` (L1-241).  The C++ function-local statics
//! (`Order`, `MatrixElemts`) become a Mutex-guarded process-wide cache with
//! identical first-call/recompute semantics.

use std::sync::Mutex;

use super::super::gauss_points::gauss_points_max;
use super::super::math_matrix::{IntegerVector, Matrix, Vector};
use super::super::p_lib_jacobi::niv_constr;
use super::super::plib::hermit_jacobi::HermitJacobi;
use super::super::GeomAbsShape;
use super::elements_of_ref_matrix::ElementsOfRefMatrix;
use super::gauss_set_integration::GaussSetIntegration;
use super::elementary_criterion::{CriterionData, ElementaryCriterion};
use super::IntArray2;

/// OCCT WDeg (cxx L42): the reference matrices are computed once per
/// constraint order at working degree 14.
pub(crate) const W_DEG: i32 = 14;

/// OCCT static Order / MatrixElemts cache (cxx L41-42).
static MATRIX_ELEMTS: Mutex<(i32, Vec<f64>)> = Mutex::new((-333, Vec::new()));

/// OCCT FEmTool_LinearTension - criterium of LinearTension for
/// Hermit-Jacobi elements (hxx L32-51).
#[derive(Debug, Clone)]
pub struct LinearTension {
    /// OCCT protected base members (myCoeff / myFirst / myLast).
    data: CriterionData,
    /// OCCT RefMatrix: math_Matrix(0, WorkDegree, 0, WorkDegree).
    ref_matrix: Matrix,
    /// OCCT myOrder.
    my_order: i32,
}

impl LinearTension {
    /// OCCT FEmTool_LinearTension::FEmTool_LinearTension (cxx L36-76).
    pub fn new(work_degree: i32, constraint_order: GeomAbsShape) -> Self {
        let mut ref_matrix = Matrix::new(0, work_degree, 0, work_degree);

        let my_order = niv_constr(constraint_order) as i32;

        {
            let mut cache = MATRIX_ELEMTS.lock().unwrap();
            if my_order != cache.0 {
                // Calculating RefMatrix
                assert!(
                    work_degree <= W_DEG,
                    "Standard_ConstructionError: Degree too high"
                );
                let der_order = 1;
                let the_base = HermitJacobi::new(W_DEG, constraint_order);
                let mut elem = ElementsOfRefMatrix::new(&the_base, der_order);

                let max_degree = W_DEG + 1;
                let an_order = IntegerVector::new_init(
                    1,
                    1,
                    (4 * (max_degree / 2 + 1)).min(gauss_points_max() as i32),
                );
                let lower = Vector::new_init(1, 1, -1.0);
                let upper = Vector::new_init(1, 1, 1.0);

                let an_int = GaussSetIntegration::new(&mut elem, &lower, &upper, &an_order);
                cache.1 = an_int.value().data.v.clone();
                cache.0 = my_order;
            }

            // Fill RefMatrix from the flat 0-based MatrixElemts.
            let matrix_elemts = &cache.1;
            let mut ii = 0usize;
            for i in 0..=work_degree {
                ref_matrix.set(i, i, matrix_elemts[ii]);
                let mut jj = ii + 1;
                for j in (i + 1)..=work_degree {
                    let v = matrix_elemts[jj];
                    ref_matrix.set(j, i, v);
                    ref_matrix.set(i, j, v);
                    jj += 1;
                }
                ii += (W_DEG + 1 - i) as usize;
            }
        }

        LinearTension {
            data: CriterionData::default(),
            ref_matrix,
            my_order,
        }
    }

    fn coeff(&self) -> std::cell::Ref<'_, Matrix> {
        self.data
            .my_coeff
            .as_ref()
            .expect("FEmTool_LinearTension: myCoeff is null")
            .borrow()
    }
}

impl ElementaryCriterion for LinearTension {
    fn data(&self) -> &CriterionData {
        &self.data
    }

    fn data_mut(&mut self) -> &mut CriterionData {
        &mut self.data
    }

    /// OCCT FEmTool_LinearTension::DependenceTable (cxx L78-97).
    fn dependence_table(&self) -> IntArray2 {
        let coeff = self.coeff();
        let mut dep_tab = IntArray2::new_init(
            coeff.lower_col(),
            coeff.upper_col(),
            coeff.lower_col(),
            coeff.upper_col(),
            0,
        );
        for i in 1..=coeff.col_number() {
            dep_tab.set_value(i, i, 1);
        }

        dep_tab
    }

    /// OCCT FEmTool_LinearTension::Value (cxx L99-149).
    fn value(&mut self) -> f64 {
        let coeff = self.coeff();
        let my_first = self.data.my_first;
        let my_last = self.data.my_last;

        let deg = (coeff.row_number() - 1).min(self.ref_matrix.upper_row());
        let coeff_lower_row = coeff.lower_row();
        let deg_h = (2 * self.my_order + 1).min(deg);
        let nb_dim = coeff.col_number();

        let mut new_coeff = Matrix::new(1, nb_dim, 0, deg);

        let coeff_half = (my_last - my_first) / 2.0;
        let cteh3 = 2.0 / coeff_half;

        let mut j = 0.0;

        for i in 0..=deg_h {
            let k1 = if i <= self.my_order { i } else { i - self.my_order - 1 };
            let mfact = coeff_half.powi(k1);
            for dim in 1..=nb_dim {
                new_coeff.set(dim, i, coeff.get(coeff_lower_row + i, dim) * mfact);
            }
        }

        for i in (deg_h + 1)..=deg {
            for dim in 1..=nb_dim {
                new_coeff.set(dim, i, coeff.get(coeff_lower_row + i, dim));
            }
        }

        for dim in 1..=nb_dim {
            for i in 0..=deg {
                let mut jline = 0.5 * self.ref_matrix.get(i, i) * new_coeff.get(dim, i);

                for j2 in 0..i {
                    jline += self.ref_matrix.get(i, j2) * new_coeff.get(dim, j2);
                }

                j += jline * new_coeff.get(dim, i);
            }
        }

        cteh3 * j
    }

    /// OCCT FEmTool_LinearTension::Hessian (cxx L151-219).
    fn hessian(&mut self, dimension1: i32, dimension2: i32, h: &mut Matrix) {
        let dep_tab = self.dependence_table();

        assert!(
            dimension1 >= dep_tab.lower_row()
                && dimension1 <= dep_tab.upper_row()
                && dimension2 >= dep_tab.lower_col()
                && dimension2 <= dep_tab.upper_col(),
            "Standard_OutOfRange: FEmTool_LinearTension::Hessian"
        );

        assert!(
            dep_tab.value(dimension1, dimension2) != 0,
            "Standard_DomainError: FEmTool_LinearTension::Hessian"
        );

        let my_first = self.data.my_first;
        let my_last = self.data.my_last;

        let deg = self.ref_matrix.upper_row().min(h.row_number() - 1);
        let deg_h = (2 * self.my_order + 1).min(deg);

        let coeff_half = (my_last - my_first) / 2.0;
        let cteh3 = 2.0 / coeff_half;

        let i0 = h.lower_row();
        let j0 = h.lower_col();

        h.init(0.0);

        let mut i1 = i0;
        for i in 0..=deg_h {
            let k1 = if i <= self.my_order { i } else { i - self.my_order - 1 };
            let mfact = coeff_half.powi(k1) * cteh3;
            // Hermite*Hermite part of matrix
            let mut j1 = j0 + i;
            for j in i..=deg_h {
                let k2 = if j <= self.my_order { j } else { j - self.my_order - 1 };
                let v = mfact * coeff_half.powi(k2) * self.ref_matrix.get(i, j);
                h.set(i1, j1, v);
                if i != j {
                    h.set(j1, i1, v);
                }
                j1 += 1;
            }
            // Hermite*Jacobi part of matrix
            let mut j1 = j0 + deg_h + 1;
            for j in (deg_h + 1)..=deg {
                let v = mfact * self.ref_matrix.get(i, j);
                h.set(i1, j1, v);
                h.set(j1, i1, v);
                j1 += 1;
            }
            i1 += 1;
        }

        // Jacoby*Jacobi part of matrix
        let mut i1 = i0 + deg_h + 1;
        for i in (deg_h + 1)..=deg {
            let mut j1 = j0 + i;
            for j in i..=deg {
                let v = cteh3 * self.ref_matrix.get(i, j);
                h.set(i1, j1, v);
                if i != j {
                    h.set(j1, i1, v);
                }
                j1 += 1;
            }
            i1 += 1;
        }
    }

    /// OCCT FEmTool_LinearTension::Gradient (cxx L221-241).
    fn gradient(&mut self, dimension: i32, g: &mut Vector) {
        let (lower_col, upper_col, lower_row, nb_rows) = {
            let coeff = self.coeff();
            (coeff.lower_col(), coeff.upper_col(), coeff.lower_row(), coeff.row_number())
        };

        assert!(
            dimension >= lower_col && dimension <= upper_col,
            "Standard_OutOfRange: FEmTool_LinearTension::Gradient"
        );

        let deg = (g.length() - 1).min(nb_rows - 1);

        let mut x = Vector::new(0, deg);
        for i in 0..=deg {
            let v = self.coeff().get(lower_row + i, dimension);
            x.set(i, v);
        }

        let mut h = Matrix::new(0, deg, 0, deg);
        self.hessian(dimension, dimension, &mut h);

        // OCCT: G.Multiply(H, X)  (G = H*X)
        for gi in g.lower()..=g.upper() {
            let mut s = 0.0;
            for hj in h.lower_col()..=h.upper_col() {
                s += h.get(gi, hj) * x.get(hj);
            }
            g.set(gi, s);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::GeomAbsShape;
    use super::super::elementary_criterion::{CoeffHandle, ElementaryCriterion};
    use super::*;

    /// FEmTool_LinearTension::Value, hand-derived from the OCCT formula at a
    /// degenerate configuration.  Element on [0, 2] (coeff = (last-first)/2
    /// = 1, cteh3 = 2/coeff = 2), WorkDegree = 2, C0 (myOrder = 0), degree
    /// truncated to 1, so NewCoeff = Coeff (all mfact = coeff^0 = 1) and
    ///     Value = 2 * [c]^T [R] [c],
    /// with the Hermite block of RefMatrix hand-integrated from
    /// H00 = (1-t)/2, H01 = (1+t)/2:
    ///     R00 = int H00'^2 = 1/2,  R11 = int H01'^2 = 1/2,
    ///     R01 = int H00' H01' = -1/2.
    /// For Coeff = (0, 1): P(s) = H01(s) = (1+s)/2, so
    ///     Value = 2 * (1/2 * 1/2 * 1) = int (dP/ds)^2 ds = 1/2.
    /// For the constant curve Coeff = (1, 1): P' = 0 => Value = 0 exactly.
    #[test]
    fn linear_tension_value_hand() {
        let mut tension = LinearTension::new(2, GeomAbsShape::C0);
        tension.set_knots(0.0, 2.0);

        // Coefficient table: (degree+1) rows x 1 dimension.  The criterion
        // addresses rows from myCoeff->LowerRow() and dimension columns from
        // 1 (OCCT Value(j0 + i, dim), dim = 1..=RowLength).
        let coeff: CoeffHandle =
            std::rc::Rc::new(std::cell::RefCell::new(Matrix::new(0, 1, 1, 1)));
        {
            let mut c = coeff.borrow_mut();
            c.set(0, 1, 0.0);
            c.set(1, 1, 1.0);
        }
        tension.set_coeff(&coeff);

        let value = tension.value();
        assert!((value - 0.5).abs() < 1e-12, "value = {}", value);

        // Constant curve has exactly zero tension.
        {
            let mut c = coeff.borrow_mut();
            c.set(0, 1, 1.0);
            c.set(1, 1, 1.0);
        }
        let value = tension.value();
        assert!(value.abs() < 1e-12, "value = {}", value);
    }

    /// FEmTool_LinearTension::Gradient hand-check for the same degenerate
    /// configuration: G = H * X with the hand-derived Hessian
    /// H = 2 * cteh3/2-scaled R block = [[1, -1], [-1, 1]]
    /// (entry = 2 * mfact(i) * mfact(j) * R(i,j), all mfact = 1).
    /// For X = (0, 1): G = (-1, 1).
    #[test]
    fn linear_tension_gradient_hand() {
        let mut tension = LinearTension::new(2, GeomAbsShape::C0);
        tension.set_knots(0.0, 2.0);

        let coeff: CoeffHandle =
            std::rc::Rc::new(std::cell::RefCell::new(Matrix::new(0, 1, 1, 1)));
        {
            let mut c = coeff.borrow_mut();
            c.set(0, 1, 0.0);
            c.set(1, 1, 1.0);
        }
        tension.set_coeff(&coeff);

        let mut g = Vector::new(0, 1);
        tension.gradient(1, &mut g);

        assert!((g.get(0) - (-1.0)).abs() < 1e-12, "g0 = {}", g.get(0));
        assert!((g.get(1) - 1.0).abs() < 1e-12, "g1 = {}", g.get(1));
    }
}


