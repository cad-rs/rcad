//! OCCT AdvApp2Var_Node (AdvApp2Var_Node.hxx + AdvApp2Var_Node.cxx) - used
//! to store constraints on a (Ui,Vj) point.
//!
//! Encoding notes (architecture):
//! - `Standard_Transient` base + RTTI: dropped (rcad has no handle runtime);
//!   the copy-assignment operator= and the deleted copy ctor become the
//!   derived `Clone` impl (the OCCT operator= copies the five fields, which
//!   is what `#[derive(Clone)]` does here).
//! - `gp_Pnt` -> `glam::DVec3`, `gp_XY` -> `glam::DVec2`.
//! - `NCollection_Array2<gp_Pnt>` / `NCollection_Array2<double>` ->
//!   [`Array2`] (bounds kept; AdvApp2Var_Node.cxx L37-38).

use glam::{DVec2, DVec3};

use super::nc_array::Array2;

/// OCCT AdvApp2Var_Node (AdvApp2Var_Node.hxx L27-88).
#[derive(Debug, Clone)]
pub struct Node {
    /// hxx L83: NCollection_Array2<gp_Pnt> myTruePoints.
    my_true_points: Array2<DVec3>,
    /// hxx L84: NCollection_Array2<double> myErrors.
    my_errors: Array2<f64>,
    /// hxx L85: gp_XY myCoord.
    my_coord: DVec2,
    /// hxx L86: int myOrdInU.
    my_ord_in_u: i32,
    /// hxx L87: int myOrdInV.
    my_ord_in_v: i32,
}

impl Node {
    /// OCCT initNodeFields (AdvApp2Var_Node.cxx L25-31) - shared
    /// initialization of myTruePoints/myErrors.
    fn init_node_fields(my_true_points: &mut Array2<DVec3>, my_errors: &mut Array2<f64>) {
        let a_zero_point = DVec3::new(0.0, 0.0, 0.0);
        my_true_points.init(a_zero_point);
        my_errors.init(0.0);
    }

    /// OCCT AdvApp2Var_Node() (AdvApp2Var_Node.cxx L36-43).
    pub fn new() -> Self {
        // myTruePoints(0, 2, 0, 2), myErrors(0, 2, 0, 2),
        // myOrdInU(2), myOrdInV(2); myCoord default (0,0).
        let mut my_true_points = Array2::new(0, 2, 0, 2);
        let mut my_errors = Array2::new(0, 2, 0, 2);
        Self::init_node_fields(&mut my_true_points, &mut my_errors);
        Node {
            my_true_points,
            my_errors,
            my_coord: DVec2::ZERO,
            my_ord_in_u: 2,
            my_ord_in_v: 2,
        }
    }

    /// OCCT AdvApp2Var_Node(const int iu, const int iv) (AdvApp2Var_Node.cxx
    /// L47-54).
    pub fn new_orders(iu: i32, iv: i32) -> Self {
        // myTruePoints(0, std::max(0, iu), 0, std::max(0, iv)), ...
        let mut my_true_points = Array2::new(0, iu.max(0), 0, iv.max(0));
        let mut my_errors = Array2::new(0, iu.max(0), 0, iv.max(0));
        Self::init_node_fields(&mut my_true_points, &mut my_errors);
        Node {
            my_true_points,
            my_errors,
            my_coord: DVec2::ZERO,
            my_ord_in_u: iu,
            my_ord_in_v: iv,
        }
    }

    /// OCCT AdvApp2Var_Node(const gp_XY& UV, const int iu, const int iv)
    /// (AdvApp2Var_Node.cxx L58-66).
    pub fn new_with_coord(uv: DVec2, iu: i32, iv: i32) -> Self {
        // myTruePoints(0, iu, 0, iv), ... myCoord(UV), ...
        // (no std::max in the OCCT ctor #3 - kept as-is).
        let mut my_true_points = Array2::new(0, iu, 0, iv);
        let mut my_errors = Array2::new(0, iu, 0, iv);
        Self::init_node_fields(&mut my_true_points, &mut my_errors);
        Node {
            my_true_points,
            my_errors,
            my_coord: uv,
            my_ord_in_u: iu,
            my_ord_in_v: iv,
        }
    }

    /// OCCT Coord() (hxx L38) - returns the coordinates (U,V) of the node.
    pub fn coord(&self) -> DVec2 {
        self.my_coord
    }

    /// OCCT SetCoord(x1, x2) (hxx L41-45) - changes the coordinates (U,V)
    /// to (x1,x2).
    pub fn set_coord(&mut self, x1: f64, x2: f64) {
        // myCoord.SetX(x1); myCoord.SetY(x2);
        self.my_coord.x = x1;
        self.my_coord.y = x2;
    }

    /// OCCT UOrder() (hxx L48) - returns the continuity order in U.
    pub fn u_order(&self) -> i32 {
        self.my_ord_in_u
    }

    /// OCCT VOrder() (hxx L51) - returns the continuity order in V.
    pub fn v_order(&self) -> i32 {
        self.my_ord_in_v
    }

    /// OCCT SetPoint(iu, iv, Pt) (hxx L54) - affects the value F(U,V) or its
    /// derivates on the node (U,V).
    pub fn set_point(&mut self, iu: i32, iv: i32, pt: DVec3) {
        self.my_true_points.set_value(iu, iv, pt);
    }

    /// OCCT Point(iu, iv) (hxx L57) - returns the value F(U,V) or its
    /// derivates on the node (U,V).
    pub fn point(&self, iu: i32, iv: i32) -> DVec3 {
        self.my_true_points.value(iu, iv)
    }

    /// OCCT SetError(iu, iv, error) (hxx L60-63) - affects the error between
    /// F(U,V) and its approximation.
    pub fn set_error(&mut self, iu: i32, iv: i32, error: f64) {
        self.my_errors.set_value(iu, iv, error);
    }

    /// OCCT Error(iu, iv) (hxx L66) - returns the error between F(U,V) and
    /// its approximation.
    pub fn error(&self, iu: i32, iv: i32) -> f64 {
        self.my_errors.value(iu, iv)
    }
}
