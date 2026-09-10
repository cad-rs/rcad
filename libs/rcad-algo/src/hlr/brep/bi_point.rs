//! OCCT HLRBRep_BiPoint (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_BiPnt.hxx` (L29-95, fully inline).  This is
//! the HLRBRep-owned point class and is distinct from the HLRAlgo_BiPoint
//! ([`crate::hlr::algo::bi_point::BiPoint`]).  `gp_Pnt` maps to `DVec3`,
//! `TopoDS_Shape` to the kernel `Shape` (a default-constructed null
//! `TopoDS_Shape` maps to [`Shape::null`]).

use glam::DVec3;
use rcad_kernel::topods::Shape;

/// OCCT HLRBRep_BiPoint — the two end points of a segment with its owning
/// shape and the line-kind flags.
#[derive(Debug, Clone)]
pub struct BiPoint {
    my_p1: DVec3,
    my_p2: DVec3,
    my_shape: Shape,
    my_rg1_line: bool,
    my_rg_n_line: bool,
    my_out_line: bool,
    my_int_line: bool,
}

impl Default for BiPoint {
    fn default() -> Self {
        Self::new()
    }
}

impl BiPoint {
    /// OCCT HLRBRep_BiPoint() — hxx L34-40 (the shape stays a null
    /// TopoDS_Shape).
    pub fn new() -> Self {
        BiPoint {
            my_p1: DVec3::ZERO,
            my_p2: DVec3::ZERO,
            my_shape: Shape::null(),
            my_rg1_line: false,
            my_rg_n_line: false,
            my_out_line: false,
            my_int_line: false,
        }
    }

    /// OCCT HLRBRep_BiPoint(x1, y1, z1, x2, y2, z2, S, reg1, regn, outl,
    /// intl) — hxx L42-61.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        x1: f64,
        y1: f64,
        z1: f64,
        x2: f64,
        y2: f64,
        z2: f64,
        s: Shape,
        reg1: bool,
        regn: bool,
        outl: bool,
        intl: bool,
    ) -> Self {
        BiPoint {
            my_p1: DVec3::new(x1, y1, z1),
            my_p2: DVec3::new(x2, y2, z2),
            my_shape: s,
            my_rg1_line: reg1,
            my_rg_n_line: regn,
            my_out_line: outl,
            my_int_line: intl,
        }
    }

    /// OCCT P1() — hxx L63.
    pub fn p1(&self) -> &DVec3 {
        &self.my_p1
    }

    /// OCCT P2() — hxx L65.
    pub fn p2(&self) -> &DVec3 {
        &self.my_p2
    }

    /// OCCT Shape() — hxx L67.
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT Shape(S) — hxx L69.
    pub fn set_shape(&mut self, s: Shape) {
        self.my_shape = s;
    }

    /// OCCT Rg1Line() — hxx L71.
    pub fn rg1_line(&self) -> bool {
        self.my_rg1_line
    }

    /// OCCT Rg1Line(B) — hxx L73.
    pub fn set_rg1_line(&mut self, b: bool) {
        self.my_rg1_line = b;
    }

    /// OCCT RgNLine() — hxx L75.
    pub fn rg_n_line(&self) -> bool {
        self.my_rg_n_line
    }

    /// OCCT RgNLine(B) — hxx L77.
    pub fn set_rg_n_line(&mut self, b: bool) {
        self.my_rg_n_line = b;
    }

    /// OCCT OutLine() — hxx L79.
    pub fn out_line(&self) -> bool {
        self.my_out_line
    }

    /// OCCT OutLine(B) — hxx L81.
    pub fn set_out_line(&mut self, b: bool) {
        self.my_out_line = b;
    }

    /// OCCT IntLine() — hxx L83.
    pub fn int_line(&self) -> bool {
        self.my_int_line
    }

    /// OCCT IntLine(B) — hxx L85.
    pub fn set_int_line(&mut self, b: bool) {
        self.my_int_line = b;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topods::Orientation;

    /// OCCT default ctor: the four flags are false and the shape is a null
    /// TopoDS_Shape (hxx L34-40).
    #[test]
    fn bi_point_default_is_null_shape_and_false_flags() {
        let bp = BiPoint::new();
        assert_eq!(*bp.p1(), DVec3::ZERO);
        assert_eq!(*bp.p2(), DVec3::ZERO);
        assert!(bp.shape().is_null());
        assert!(!bp.rg1_line());
        assert!(!bp.rg_n_line());
        assert!(!bp.out_line());
        assert!(!bp.int_line());
    }

    /// OCCT full ctor + setter round-trip (hxx L42-61, L63-85).
    #[test]
    fn bi_point_ctor_and_setter_roundtrip() {
        let shape = Shape::synthetic(7, Orientation::Forward);
        let mut bp = BiPoint::from_parts(
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, shape.clone(), true, false, true, false,
        );
        assert_eq!(*bp.p1(), DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(*bp.p2(), DVec3::new(4.0, 5.0, 6.0));
        assert_eq!(bp.shape().index, 7);
        assert!(bp.rg1_line());
        assert!(!bp.rg_n_line());
        assert!(bp.out_line());
        assert!(!bp.int_line());
        // Setters.
        bp.set_shape(Shape::null());
        assert!(bp.shape().is_null());
        bp.set_rg1_line(false);
        bp.set_rg_n_line(true);
        bp.set_out_line(false);
        bp.set_int_line(true);
        assert!(!bp.rg1_line());
        assert!(bp.rg_n_line());
        assert!(!bp.out_line());
        assert!(bp.int_line());
    }
}
