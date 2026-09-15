//! OCCT FEmTool_ElementsOfRefMatrix (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_ElementsOfRefMatrix.hxx` (L17-65) and
//! `FEmTool_ElementsOfRefMatrix.cxx` (L17-83).  Required by the
//! FEmTool_LinearTension / LinearFlexion / LinearJerk constructors, which
//! integrate it to build their reference matrices.

use super::super::math_matrix::Vector;
use super::gauss_set_integration::FunctionSet;
use super::super::plib::hermit_jacobi::HermitJacobi;

/// OCCT FEmTool_ElementsOfRefMatrix - functions needed for calculating
/// matrix elements of RefMatrix for linear criteriums (Tension, Flexion and
/// Jerk) by Gauss integration.  Each function from set gives value
/// Pi(u)'*Pj(u)' or Pi(u)''*Pj(u)'' or Pi(u)'''*Pj(u)''' for each i and j,
/// where Pi(u) is i-th basis function of expansion.
#[derive(Debug, Clone)]
pub struct ElementsOfRefMatrix {
    my_base: HermitJacobi,
    my_der_order: i32,
    my_nb_equations: i32,
}

impl ElementsOfRefMatrix {
    /// OCCT FEmTool_ElementsOfRefMatrix::FEmTool_ElementsOfRefMatrix
    /// (cxx L22-33).  OCCT copies the base (PLib_HermitJacobi myBase).
    pub fn new(the_base: &HermitJacobi, der_order: i32) -> Self {
        assert!(
            (0..=3).contains(&der_order),
            "Standard_ConstructionError: FEmTool_ElementsOfRefMatrix"
        );
        let my_nb_equations = (the_base.work_degree() + 2) * (the_base.work_degree() + 1) / 2;
        ElementsOfRefMatrix {
            my_base: the_base.clone(),
            my_der_order: der_order,
            my_nb_equations,
        }
    }

    /// OCCT FEmTool_ElementsOfRefMatrix::NbVariables (cxx L35-38).
    pub fn nb_variables(&self) -> i32 {
        1
    }

    /// OCCT FEmTool_ElementsOfRefMatrix::NbEquations (cxx L40-43).
    pub fn nb_equations(&self) -> i32 {
        self.my_nb_equations
    }

    /// OCCT FEmTool_ElementsOfRefMatrix::Value (cxx L45-82).
    pub fn value(&mut self, x: &Vector, f: &mut Vector) -> bool {
        assert!(
            f.length() >= self.my_nb_equations,
            "Standard_OutOfRange: FEmTool_ElementsOfRefMatrix::Value"
        );

        let u = x.get(x.lower());
        let work_degree = self.my_base.work_degree();
        let mut basis = vec![0.0f64; (work_degree + 1) as usize];
        // OCCT reuses one Aux array for all intermediate derivative slots
        // (D1(U, Aux, Basis), D2(U, Aux, Aux, Basis), ...); distinct scratch
        // buffers are equivalent here because the intermediate slots are
        // never read by this consumer.
        let mut aux = vec![0.0f64; (work_degree + 1) as usize];
        let mut aux2 = vec![0.0f64; (work_degree + 1) as usize];
        let mut aux3 = vec![0.0f64; (work_degree + 1) as usize];

        match self.my_der_order {
            0 => {
                self.my_base.d0(u, &mut basis);
            }
            1 => {
                self.my_base.d1(u, &mut aux, &mut basis);
            }
            2 => {
                self.my_base.d2(u, &mut aux, &mut aux2, &mut basis);
            }
            _ => {
                self.my_base
                    .d3(u, &mut aux, &mut aux2, &mut aux3, &mut basis);
            }
        }

        let mut ii = 0usize;
        for i in 0..=work_degree {
            for j in i..=work_degree {
                f.set(f.lower() + ii as i32, basis[i as usize] * basis[j as usize]);
                ii += 1;
            }
        }

        true
    }
}

impl FunctionSet for ElementsOfRefMatrix {
    fn nb_variables(&self) -> i32 {
        ElementsOfRefMatrix::nb_variables(self)
    }

    fn nb_equations(&self) -> i32 {
        ElementsOfRefMatrix::nb_equations(self)
    }

    fn value(&mut self, x: &Vector, f: &mut Vector) -> bool {
        ElementsOfRefMatrix::value(self, x, f)
    }
}
