//! OCCT AppDef_SmoothCriterion (ModelingData/TKGeomBase/AppDef).
//!
//! 1:1 translation of `AppDef_SmoothCriterion.hxx` (L17-101) and
//! `AppDef_SmoothCriterion.cxx` (L17-22; the .cxx only carries the RTTI
//! implementation macro, no method bodies).
//!
//! Architecture mapping: the OCCT abstract base class (a Standard_Transient
//! with pure virtuals) becomes the [`SmoothCriterion`] trait, mirroring the
//! kernel `fem_tool::ElementaryCriterion` trait style.  The OCCT handles
//! become:
//! - `occ::handle<FEmTool_Curve>` -> [`CurveHandle`] (`Rc<RefCell<Curve>>`),
//!   preserving the OCCT handle aliasing (AppDef_Variational keeps mutating
//!   the curve through its own handle after `SetCurve`, and the criterion
//!   must see the updates);
//! - `occ::handle<NCollection_HArray1<double>>` (Parameters, and the
//!   `SetWeight(Array1)` argument) -> [`RealArray1`] by value.  OCCT only
//!   stores these (`myParameters = Parameters`), and the
//!   `NCollection_Array1` assignment `myPntWeight = Weight` copies bounds
//!   and data (NCollection_Array1.hxx L273, L325), so the Rust clone is the
//!   exact semantic counterpart;
//! - `handle<HArray2<handle<HArray1<int>>>>` (the assembly table) ->
//!   `fem_tool::AssemblingTable`, freshly built by `AssemblyTable()` like
//!   the OCCT `new NCollection_HArray2<...>`.

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::fem_tool::{AssemblingTable, Curve, IntArray2};
use rcad_kernel::math::math_matrix::{Matrix, Vector};

/// OCCT `occ::handle<FEmTool_Curve>` - shared, mutable curve (aliasing
/// semantics preserved).
pub type CurveHandle = Rc<RefCell<Curve>>;

/// Mirror of `NCollection_HArray1<double>`: a 1-based-storage double array
/// with an arbitrary lower bound (only the forms consumed by
/// AppDef_SmoothCriterion / AppDef_LinearCriteria are mirrored).
#[derive(Debug, Clone, Default)]
pub struct RealArray1 {
    lower: i32,
    data: Vec<f64>,
}

impl RealArray1 {
    /// OCCT new NCollection_HArray1<double>(Lower, Upper).
    pub fn new(lower: i32, upper: i32) -> Self {
        assert!(upper >= lower - 1, "NCollection_HArray1: bad range");
        RealArray1 {
            lower,
            data: vec![0.0; (upper - lower + 1).max(0) as usize],
        }
    }

    /// OCCT Init(Value).
    pub fn init(&mut self, value: f64) {
        self.data.fill(value);
    }

    /// OCCT Lower().
    #[inline]
    pub fn lower(&self) -> i32 {
        self.lower
    }

    /// OCCT Upper().
    #[inline]
    pub fn upper(&self) -> i32 {
        self.lower + self.data.len() as i32 - 1
    }

    /// OCCT Length().
    #[inline]
    pub fn length(&self) -> i32 {
        self.data.len() as i32
    }

    /// OCCT Value(Index).
    #[inline]
    pub fn value(&self, index: i32) -> f64 {
        self.data[(index - self.lower) as usize]
    }

    /// OCCT SetValue(Index, Value).
    #[inline]
    pub fn set_value(&mut self, index: i32, value: f64) {
        self.data[(index - self.lower) as usize] = value;
    }
}

/// OCCT AppDef_SmoothCriterion - defined criterion to smooth points in
/// curve (AppDef_SmoothCriterion.hxx L34-98).
pub trait SmoothCriterion {
    /// OCCT SetParameters (hxx L38-39).
    fn set_parameters(&mut self, parameters: &RealArray1);

    /// OCCT SetCurve (hxx L41).
    fn set_curve(&mut self, c: &CurveHandle);

    /// OCCT Curve() (hxx L45-50) - convenience accessor implemented on the
    /// base class by calling GetCurve; provided here as a default method.
    fn curve(&self) -> Option<CurveHandle> {
        let mut a_curve = None;
        self.get_curve(&mut a_curve);
        a_curve
    }

    /// OCCT GetCurve (hxx L52).
    fn get_curve(&self, c: &mut Option<CurveHandle>);

    /// OCCT SetEstimation (hxx L54).
    fn set_estimation(&mut self, e1: f64, e2: f64, e3: f64);

    /// OCCT EstLength (hxx L56) - the OCCT `double&` return becomes the
    /// mutable reference.
    fn est_length(&mut self) -> &mut f64;

    /// OCCT GetEstimation (hxx L58).
    fn get_estimation(&self, e1: &mut f64, e2: &mut f64, e3: &mut f64);

    /// OCCT AssemblyTable (hxx L60-61).
    fn assembly_table(&self) -> AssemblingTable;

    /// OCCT DependenceTable (hxx L63).
    fn dependence_table(&self) -> IntArray2;

    /// OCCT QualityValues (hxx L65-70).
    fn quality_values(
        &mut self,
        j1min: f64,
        j2min: f64,
        j3min: f64,
        j1: &mut f64,
        j2: &mut f64,
        j3: &mut f64,
    ) -> i32;

    /// OCCT ErrorValues (hxx L72-74).
    fn error_values(
        &mut self,
        max_error: &mut f64,
        quadratic_error: &mut f64,
        average_error: &mut f64,
    );

    /// OCCT Hessian (hxx L76-79).
    fn hessian(&mut self, element: i32, dimension1: i32, dimension2: i32, h: &mut Matrix);

    /// OCCT Gradient (hxx L81).
    fn gradient(&mut self, element: i32, dimension: i32, g: &mut Vector);

    /// OCCT InputVector (hxx L84-86) - convert the assembly Vector in an
    /// Curve.
    fn input_vector(&mut self, x: &Vector, ass_table: &AssemblingTable);

    /// OCCT SetWeight(QuadraticWeight, QualityWeight, percentJ1, percentJ2,
    /// percentJ3) (hxx L88-92).
    fn set_weight(
        &mut self,
        quadratic_weight: f64,
        quality_weight: f64,
        percent_j1: f64,
        percent_j2: f64,
        percent_j3: f64,
    );

    /// OCCT GetWeight (hxx L94).
    fn get_weight(&self, quadratic_weight: &mut f64, quality_weight: &mut f64);

    /// OCCT SetWeight(const NCollection_Array1<double>&) overload (hxx L96).
    /// Rust has no overloading; the array form is named after the OCCT
    /// argument.
    fn set_weight_array(&mut self, weight: &RealArray1);
}
