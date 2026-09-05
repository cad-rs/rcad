//! OCCT HLRBRep_BiPnt2D (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_BiPnt2D.hxx` (L27-108, fully inline).
//! `gp_Pnt2d` maps to `DVec2`, `TopoDS_Shape` to the kernel `Shape` (a
//! default-constructed null `TopoDS_Shape` maps to [`Shape::null`]).

use glam::DVec2;
use rcad_kernel::topods::Shape;

/// OCCT HLRBRep_BiPnt2D — the two end points of a projected segment with
/// its owning shape and the line-kind flags.
#[derive(Debug, Clone)]
pub struct BiPnt2D {
    my_p1: DVec2,
    my_p2: DVec2,
    my_shape: Shape,
    my_rg1_line: bool,
    my_rg_n_line: bool,
    my_out_line: bool,
    my_int_line: bool,
}

impl Default for BiPnt2D {
    fn default() -> Self {
        Self::new()
    }
}

impl BiPnt2D {
    /// OCCT HLRBRep_BiPnt2D() — hxx L32-38 (the shape stays a null
    /// TopoDS_Shape).
    pub fn new() -> Self {
        BiPnt2D {
            my_p1: DVec2::ZERO,
            my_p2: DVec2::ZERO,
            my_shape: Shape::null(),
            my_rg1_line: false,
            my_rg_n_line: false,
            my_out_line: false,
            my_int_line: false,
        }
    }

    /// OCCT HLRBRep_BiPnt2D(x1, y1, x2, y2, S, reg1, regn, outl, intl) —
    /// hxx L40-57.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        s: Shape,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        BiPnt2D {
            my_p1: DVec2::new(x1, y1),
            my_p2: DVec2::new(x2, y2),
            my_shape: s,
            my_rg1_line: reg1,
            my_rg_n_line: regn,
            my_out_line: outl,
            my_int_line: intl,
        }
    }

    /// OCCT HLRBRep_BiPnt2D(thePoint1, thePoint2, S, reg1, regn, outl,
    /// intl) — hxx L59-74 (the gp_XY overload).
    #[allow(clippy::too_many_arguments)]
    pub fn from_xy(
        the_point1: DVec2,
        the_point2: DVec2,
        s: Shape,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        BiPnt2D {
            my_p1: the_point1,
            my_p2: the_point2,
            my_shape: s,
            my_rg1_line: reg1,
            my_rg_n_line: regn,
            my_out_line: outl,
            my_int_line: intl,
        }
    }

    /// OCCT P1() — hxx L76.
    pub fn p1(&self) -> &DVec2 {
        &self.my_p1
    }

    /// OCCT P2() — hxx L78.
    pub fn p2(&self) -> &DVec2 {
        &self.my_p2
    }

    /// OCCT Shape() — hxx L80.
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT Shape(S) — hxx L82.
    pub fn set_shape(&mut self, s: Shape) {
        self.my_shape = s;
    }

    /// OCCT Rg1Line() — hxx L84.
    pub fn rg1_line(&self) -> bool {
        self.my_rg1_line
    }

    /// OCCT Rg1Line(B) — hxx L86.
    pub fn set_rg1_line(&mut self, b: bool) {
        self.my_rg1_line = b;
    }

    /// OCCT RgNLine() — hxx L88.
    pub fn rg_n_line(&self) -> bool {
        self.my_rg_n_line
    }

    /// OCCT RgNLine(B) — hxx L90.
    pub fn set_rg_n_line(&mut self, b: bool) {
        self.my_rg_n_line = b;
    }

    /// OCCT OutLine() — hxx L92.
    pub fn out_line(&self) -> bool {
        self.my_out_line
    }

    /// OCCT OutLine(B) — hxx L94.
    pub fn set_out_line(&mut self, b: bool) {
        self.my_out_line = b;
    }

    /// OCCT IntLine() — hxx L96.
    pub fn int_line(&self) -> bool {
        self.my_int_line
    }

    /// OCCT IntLine(B) — hxx L98.
    pub fn set_int_line(&mut self, b: bool) {
        self.my_int_line = b;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topods::Orientation;

    /// OCCT default ctor: the four flags are false and the shape is a null
    /// TopoDS_Shape (hxx L32-38).
    #[test]
    fn bi_pnt_2d_default_is_null_shape_and_false_flags() {
        let bp = BiPnt2D::new();
        assert_eq!(*bp.p1(), DVec2::ZERO);
        assert_eq!(*bp.p2(), DVec2::ZERO);
        assert!(bp.shape().is_null());
        assert!(!bp.rg1_line());
        assert!(!bp.rg_n_line());
        assert!(!bp.out_line());
        assert!(!bp.int_line());
    }

    /// OCCT scalar ctor + setter round-trip (hxx L40-57, L76-98).
    #[test]
    fn bi_pnt_2d_ctor_and_setter_roundtrip() {
        let shape = Shape::synthetic(3, Orientation::Reversed);
        let mut bp = BiPnt2D::from_parts(
            1.0, 2.0, 3.0, 4.0, shape.clone(), true, true, false, true,
        );
        assert_eq!(*bp.p1(), DVec2::new(1.0, 2.0));
        assert_eq!(*bp.p2(), DVec2::new(3.0, 4.0));
        assert_eq!(bp.shape().index, 3);
        assert!(bp.rg1_line());
        assert!(bp.rg_n_line());
        assert!(!bp.out_line());
        assert!(bp.int_line());
        bp.set_shape(Shape::null());
        assert!(bp.shape().is_null());
        bp.set_rg1_line(false);
        bp.set_rg_n_line(false);
        bp.set_out_line(true);
        bp.set_int_line(false);
        assert!(!bp.rg1_line());
        assert!(!bp.rg_n_line());
        assert!(bp.out_line());
        assert!(!bp.int_line());
    }

    /// OCCT gp_XY ctor (hxx L59-74) — identical field mapping.
    #[test]
    fn bi_pnt_2d_xy_ctor() {
        let shape = Shape::synthetic(4, Orientation::Forward);
        let bp = BiPnt2D::from_xy(
            DVec2::new(5.0, 6.0),
            DVec2::new(7.0, 8.0),
            shape,
            false,
            true,
            false,
            false,
        );
        assert_eq!(*bp.p1(), DVec2::new(5.0, 6.0));
        assert_eq!(*bp.p2(), DVec2::new(7.0, 8.0));
        assert_eq!(bp.shape().index, 4);
        assert!(!bp.rg1_line());
        assert!(bp.rg_n_line());
        assert!(!bp.out_line());
        assert!(!bp.int_line());
    }
}
