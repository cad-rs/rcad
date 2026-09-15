//! OCCT FEmTool_ElementaryCriterion (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_ElementaryCriterion.hxx` (L17-68) and
//! `FEmTool_ElementaryCriterion.cxx` (L17-39).
//!
//! Architecture mapping: the OCCT abstract base class (Standard_Transient
//! with protected members `myCoeff` / `myFirst` / `myLast`) becomes the
//! [`ElementaryCriterion`] trait plus the shared [`CriterionData`] storage
//! that each concrete criterion embeds; the trait carries the two concrete
//! `Set` overloads as provided methods.  The OCCT shared handle
//! (`Handle(NCollection_HArray2<double>)`) becomes [`CoeffHandle`]
//! (`Rc<RefCell<Matrix>>`), preserving the OCCT aliasing semantics: the
//! caller keeps mutating the coefficients in place after `Set`.

use std::cell::RefCell;
use std::rc::Rc;

use super::super::math_matrix::{Matrix, Vector};
use super::IntArray2;

/// OCCT `Handle(NCollection_HArray2<double>) myCoeff` - shared, mutable
/// coefficient table (aliasing semantics preserved).
pub type CoeffHandle = Rc<RefCell<Matrix>>;

/// OCCT protected base-class data of FEmTool_ElementaryCriterion
/// (hxx L63-65): `myCoeff`, `myFirst`, `myLast`.
#[derive(Debug, Clone, Default)]
pub struct CriterionData {
    pub my_coeff: Option<CoeffHandle>,
    pub my_first: f64,
    pub my_last: f64,
}

/// OCCT FEmTool_ElementaryCriterion - defined J criteria used in
/// minimisation (hxx L31-60).
pub trait ElementaryCriterion {
    /// Access to the protected base-class members (trait substitute).
    fn data(&self) -> &CriterionData;

    /// Mutable access to the protected base-class members (trait substitute).
    fn data_mut(&mut self) -> &mut CriterionData;

    /// OCCT FEmTool_ElementaryCriterion::Set(Coeff) (cxx L29-32) - set the
    /// coefficient of the Element (the Curve).  The handle is shared, not
    /// copied.
    fn set_coeff(&mut self, coeff: &CoeffHandle) {
        self.data_mut().my_coeff = Some(coeff.clone());
    }

    /// OCCT FEmTool_ElementaryCriterion::Set(FirstKnot, LastKnot)
    /// (cxx L34-38) - set the definition interval of the Element.
    fn set_knots(&mut self, first_knot: f64, last_knot: f64) {
        self.data_mut().my_first = first_knot;
        self.data_mut().my_last = last_knot;
    }

    /// OCCT FEmTool_ElementaryCriterion::DependenceTable (hxx L45) - to know
    /// if two dimensions are independent.
    fn dependence_table(&self) -> IntArray2;

    /// OCCT FEmTool_ElementaryCriterion::Value (hxx L48) - to compute J(E)
    /// where E is the current Element.
    fn value(&mut self) -> f64;

    /// OCCT FEmTool_ElementaryCriterion::Hessian (hxx L54) - to compute the
    /// coefficients of the Hessian matrix of J(E) which are crossed
    /// derivatives in dimensions Dim1 and Dim2.
    fn hessian(&mut self, dimension1: i32, dimension2: i32, h: &mut Matrix);

    /// OCCT FEmTool_ElementaryCriterion::Gradient (hxx L58) - to compute the
    /// coefficients in the dimension Dim of the J(E)'s Gradient.
    fn gradient(&mut self, dimension: i32, g: &mut Vector);
}
