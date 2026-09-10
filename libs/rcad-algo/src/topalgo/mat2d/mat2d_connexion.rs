//! OCCT MAT2d_Connexion — a connexion linking two lines of items in a set
//! of lines.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_Connexion.hxx L17-122, MAT2d_Connexion.cxx L26-274
//!
//! OCCT MAT2d_Connexion is a Standard_Transient manipulated through
//! occ::handle.  MAT2d_Circuit::DoubleLine mutates connexion objects that
//! are aliased inside MAT2d_MiniPath (theConnexions / theFather / thePath),
//! so the handle is Arc + RwLock (AGENTS.md Handle mapping).

use std::sync::{Arc, RwLock};

use glam::DVec2;

/// OCCT `occ::handle<MAT2d_Connexion>`.
pub type HandleMat2dConnexion = Arc<RwLock<Mat2dConnexion>>;

/// OCCT MAT2d_Connexion (MAT2d_Connexion.hxx L31-120).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat2dConnexion {
    line_a: i32,
    line_b: i32,
    item_a: i32,
    item_b: i32,
    distance: f64,
    parameter_on_a: f64,
    parameter_on_b: f64,
    point_a: DVec2,
    point_b: DVec2,
}

impl Mat2dConnexion {
    /// OCCT MAT2d_Connexion::MAT2d_Connexion() (cxx L26-35).
    pub fn new() -> Self {
        Mat2dConnexion {
            line_a: 0,
            line_b: 0,
            item_a: 0,
            item_b: 0,
            distance: 0.0,
            parameter_on_a: 0.0,
            parameter_on_b: 0.0,
            point_a: DVec2::ZERO,
            point_b: DVec2::ZERO,
        }
    }

    /// OCCT MAT2d_Connexion::MAT2d_Connexion(LineA, LineB, ItemA, ItemB,
    /// Distance, ParameterOnA, ParameterOnB, PointA, PointB) (cxx L39-58).
    #[allow(clippy::too_many_arguments)]
    pub fn new_full(
        line_a: i32,
        line_b: i32,
        item_a: i32,
        item_b: i32,
        distance: f64,
        parameter_on_a: f64,
        parameter_on_b: f64,
        point_a: DVec2,
        point_b: DVec2,
    ) -> Self {
        Mat2dConnexion {
            line_a,
            line_b,
            item_a,
            item_b,
            distance,
            parameter_on_a,
            parameter_on_b,
            point_a,
            point_b,
        }
    }

    /// OCCT MAT2d_Connexion::IndexFirstLine() const (cxx L62-65) — the index
    /// on the first line.
    pub fn index_first_line(&self) -> i32 {
        self.line_a
    }

    /// OCCT MAT2d_Connexion::IndexSecondLine() const (cxx L69-72).
    pub fn index_second_line(&self) -> i32 {
        self.line_b
    }

    /// OCCT MAT2d_Connexion::IndexItemOnFirst() const (cxx L76-79).
    pub fn index_item_on_first(&self) -> i32 {
        self.item_a
    }

    /// OCCT MAT2d_Connexion::IndexItemOnSecond() const (cxx L83-86).
    pub fn index_item_on_second(&self) -> i32 {
        self.item_b
    }

    /// OCCT MAT2d_Connexion::ParameterOnFirst() const (cxx L90-93).
    pub fn parameter_on_first(&self) -> f64 {
        self.parameter_on_a
    }

    /// OCCT MAT2d_Connexion::ParameterOnSecond() const (cxx L97-100).
    pub fn parameter_on_second(&self) -> f64 {
        self.parameter_on_b
    }

    /// OCCT MAT2d_Connexion::PointOnFirst() const (cxx L104-107).
    pub fn point_on_first(&self) -> DVec2 {
        self.point_a
    }

    /// OCCT MAT2d_Connexion::PointOnSecond() const (cxx L111-114).
    pub fn point_on_second(&self) -> DVec2 {
        self.point_b
    }

    /// OCCT MAT2d_Connexion::Distance() const (cxx L118-121).
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// OCCT MAT2d_Connexion::IndexFirstLine(anIndex) (cxx L125-128).
    pub fn set_index_first_line(&mut self, an_index: i32) {
        self.line_a = an_index;
    }

    /// OCCT MAT2d_Connexion::IndexSecondLine(anIndex) (cxx L132-135).
    pub fn set_index_second_line(&mut self, an_index: i32) {
        self.line_b = an_index;
    }

    /// OCCT MAT2d_Connexion::IndexItemOnFirst(anIndex) (cxx L139-143).
    pub fn set_index_item_on_first(&mut self, an_index: i32) {
        self.item_a = an_index;
    }

    /// OCCT MAT2d_Connexion::IndexItemOnSecond(anIndex) (cxx L146-150).
    pub fn set_index_item_on_second(&mut self, an_index: i32) {
        self.item_b = an_index;
    }

    /// OCCT MAT2d_Connexion::ParameterOnFirst(aParameter) (cxx L153-157).
    pub fn set_parameter_on_first(&mut self, a_parameter: f64) {
        self.parameter_on_a = a_parameter;
    }

    /// OCCT MAT2d_Connexion::ParameterOnSecond(aParameter) (cxx L160-164).
    pub fn set_parameter_on_second(&mut self, a_parameter: f64) {
        self.parameter_on_b = a_parameter;
    }

    /// OCCT MAT2d_Connexion::PointOnFirst(aPoint) (cxx L167-171).
    pub fn set_point_on_first(&mut self, a_point: DVec2) {
        self.point_a = a_point;
    }

    /// OCCT MAT2d_Connexion::PointOnSecond(aPoint) (cxx L174-178).
    pub fn set_point_on_second(&mut self, a_point: DVec2) {
        self.point_b = a_point;
    }

    /// OCCT MAT2d_Connexion::Distance(d) (cxx L181-184).
    pub fn set_distance(&mut self, d: f64) {
        self.distance = d;
    }

    /// OCCT MAT2d_Connexion::Reverse() const (cxx L188-199) — returns the
    /// reverse connexion: the firstpoint is the secondpoint, the secondpoint
    /// is the firstpoint.
    pub fn reverse(&self) -> HandleMat2dConnexion {
        Arc::new(RwLock::new(Mat2dConnexion::new_full(
            self.line_b,
            self.line_a,
            self.item_b,
            self.item_a,
            self.distance,
            self.parameter_on_b,
            self.parameter_on_a,
            self.point_b,
            self.point_a,
        )))
    }

    /// OCCT MAT2d_Connexion::IsAfter(C2, Sense) const (cxx L203-231).
    ///
    /// Returns true if my firstPoint is on the same line than the firstpoint
    /// of `a_connexion` and my firstpoint is after the firstpoint of
    /// `a_connexion` on the line.  `a_sense` = 1 if `a_connexion` is on the
    /// Left of its firstline, else `a_sense` = -1.
    pub fn is_after(&self, a_connexion: &HandleMat2dConnexion, a_sense: f64) -> bool {
        let c2 = a_connexion.read().unwrap();
        if self.line_a != c2.index_first_line() {
            return false;
        }

        if self.item_a > c2.index_item_on_first() {
            return true;
        } else if self.item_a == c2.index_item_on_first() {
            if self.parameter_on_a > c2.parameter_on_first() {
                return true;
            } else if self.parameter_on_a == c2.parameter_on_first() {
                // gp_Vec2d Vect1(C2->PointOnFirst(), C2->PointOnSecond());
                // gp_Vec2d Vect2(pointA, pointB);
                // (Vect1 ^ Vect2) — the 2d cross product (gp_Vec2d::Crossed).
                let vect1 = c2.point_on_second() - c2.point_on_first();
                let vect2 = self.point_b - self.point_a;
                let cross = vect1.x * vect2.y - vect1.y * vect2.x;
                if cross * a_sense > 0.0 {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT MAT2d_Connexion::Dump(Deep, Offset) const (cxx L246-274).
    pub fn dump(&self, _deep: i32, offset: i32) {
        // static void Indent(Offset) (cxx L233-242).
        fn indent(offset: i32) {
            for _ in 0..offset {
                print!(" ");
            }
        }
        let mut my_offset = offset;
        indent(offset);
        println!("MAT2d_Connexion :");
        my_offset += 1;
        indent(my_offset);
        println!("IndexFirstLine    :{}", self.line_a);
        indent(my_offset);
        println!("IndexSecondLine   :{}", self.line_b);
        indent(my_offset);
        println!("IndexItemOnFirst  :{}", self.item_a);
        indent(my_offset);
        println!("IndexItemOnSecond :{}", self.item_b);
        indent(my_offset);
        println!("ParameterOnFirst  :{}", self.parameter_on_a);
        indent(my_offset);
        println!("ParameterOnSecond :{}", self.parameter_on_b);
        indent(my_offset);
        println!("PointOnFirst      :");
        println!("  X = {}", self.point_a.x);
        println!("  Y = {}", self.point_a.y);
        indent(my_offset);
        println!("PointOnSecond     :");
        println!("  X = {}", self.point_b.x);
        println!("  Y = {}", self.point_b.y);
        indent(my_offset);
        println!("Distance          :{}", self.distance);
    }
}

impl Default for Mat2dConnexion {
    fn default() -> Self {
        Self::new()
    }
}
