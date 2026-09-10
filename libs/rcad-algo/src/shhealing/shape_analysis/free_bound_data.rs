//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_FreeBoundData` (`ShapeAnalysis_FreeBoundData.hxx`
//! L17-145 + `.cxx` L1-98).
//!
//! Represents a free bound (a wire) and stores its properties: area of the
//! contour, perimeter, ratio of average length to average width, average
//! width, and the notches (narrow 'V'-like sub-contours) with their maximum
//! widths.
//!
//! Architecture bridges:
//! 1. `TopoDS_Wire` -> `rcad_kernel::topo_shape::Shape`.
//! 2. `occ::handle<NCollection_HSequence<TopoDS_Shape>>` -> `Vec<Shape>`.
//! 3. `NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher>
//!    myNotchesParams` -> a `HashMap<(u64, u32), f64>` keyed by the IsSame
//!    identity (TShape pointer + location) — `IsBound` / `Bind` / `Find`
//!    map to the entry API.

use std::collections::HashMap;

use rcad_kernel::topo_shape::Shape;

/// OCCT ShapeAnalysis_FreeBoundData (hxx L34-145).
#[derive(Clone)]
pub struct ShapeAnalysisFreeBoundData {
    /// OCCT `myBound`.
    my_bound: Shape,
    /// OCCT `myArea`.
    my_area: f64,
    /// OCCT `myPerimeter`.
    my_perimeter: f64,
    /// OCCT `myRatio`.
    my_ratio: f64,
    /// OCCT `myWidth`.
    my_width: f64,
    /// OCCT `myNotches`.
    my_notches: Vec<Shape>,
    /// OCCT `myNotchesParams`.
    my_notches_params: HashMap<(u64, u32), f64>,
}

impl Default for ShapeAnalysisFreeBoundData {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisFreeBoundData {
    /// OCCT ShapeAnalysis_FreeBoundData() (cxx L33-38): myNotches = new
    /// sequence; Clear().
    pub fn new() -> Self {
        let mut this = ShapeAnalysisFreeBoundData {
            my_bound: Shape::null(),
            my_area: 0.0,
            my_perimeter: 0.0,
            my_ratio: 0.0,
            my_width: 0.0,
            my_notches: Vec::new(),
            my_notches_params: HashMap::new(),
        };
        this.clear();
        this
    }

    /// OCCT ShapeAnalysis_FreeBoundData(freebound) (cxx L46-52).
    pub fn new_with_bound(freebound: &Shape) -> Self {
        let mut this = Self::new();
        this.set_free_bound(freebound);
        this
    }

    /// OCCT Clear() (cxx L59-69): clears all properties of the contour;
    /// the contour bound itself is not cleared.
    pub fn clear(&mut self) {
        self.my_area = -1.0;
        self.my_perimeter = -1.0;
        self.my_ratio = -1.0;
        self.my_width = -1.0;
        self.my_notches.clear();
        self.my_notches_params.clear();
    }

    /// OCCT SetFreeBound(freebound) (hxx L84-85 inline).
    pub fn set_free_bound(&mut self, freebound: &Shape) {
        self.my_bound = freebound.clone();
    }

    /// OCCT SetArea(area) (hxx L88-89 inline).
    pub fn set_area(&mut self, area: f64) {
        self.my_area = area;
    }

    /// OCCT SetPerimeter(perimeter) (hxx L92-93 inline).
    pub fn set_perimeter(&mut self, perimeter: f64) {
        self.my_perimeter = perimeter;
    }

    /// OCCT SetRatio(ratio) (hxx L96-97 inline).
    pub fn set_ratio(&mut self, ratio: f64) {
        self.my_ratio = ratio;
    }

    /// OCCT SetWidth(width) (hxx L100-101 inline).
    pub fn set_width(&mut self, width: f64) {
        self.my_width = width;
    }

    /// OCCT AddNotch(notch, width) (cxx L76-85).
    pub fn add_notch(&mut self, notch: &Shape, width: f64) {
        let key = (std::sync::Arc::as_ptr(&notch.data) as u64, notch.location);
        if self.my_notches_params.contains_key(&key) {
            return;
        }
        self.my_notches.push(notch.clone());
        self.my_notches_params.insert(key, width);
    }

    /// OCCT FreeBound() (hxx L104-105 inline).
    pub fn free_bound(&self) -> Shape {
        self.my_bound.clone()
    }

    /// OCCT Area() (hxx L108-109 inline).
    pub fn area(&self) -> f64 {
        self.my_area
    }

    /// OCCT Perimeter() (hxx L112-113 inline).
    pub fn perimeter(&self) -> f64 {
        self.my_perimeter
    }

    /// OCCT Ratio() (hxx L116-117 inline).
    pub fn ratio(&self) -> f64 {
        self.my_ratio
    }

    /// OCCT Width() (hxx L120-121 inline).
    pub fn width(&self) -> f64 {
        self.my_width
    }

    /// OCCT NbNotches() (hxx L124-125 inline).
    pub fn nb_notches(&self) -> i32 {
        self.my_notches.len() as i32
    }

    /// OCCT Notches() (hxx L128-129 inline).
    pub fn notches(&self) -> &[Shape] {
        &self.my_notches
    }

    /// OCCT Notch(index) (hxx L132-133 inline).
    pub fn notch(&self, index: i32) -> Shape {
        self.my_notches[(index - 1) as usize].clone()
    }

    /// OCCT NotchWidth(index) (cxx L90-95).
    pub fn notch_width(&self, index: i32) -> f64 {
        let wire = self.my_notches[(index - 1) as usize].clone();
        self.notch_width_of(&wire)
    }

    /// OCCT NotchWidth(notch) (cxx L100-104) — `myNotchesParams.Find(notch)`
    /// raises Standard_NoSuchObject on an unbound notch (the rcad panic).
    pub fn notch_width_of(&self, notch: &Shape) -> f64 {
        let key = (std::sync::Arc::as_ptr(&notch.data) as u64, notch.location);
        match self.my_notches_params.get(&key) {
            Some(&w) => w,
            None => panic!("Standard_NoSuchObject: ShapeAnalysis_FreeBoundData::NotchWidth"),
        }
    }
}
