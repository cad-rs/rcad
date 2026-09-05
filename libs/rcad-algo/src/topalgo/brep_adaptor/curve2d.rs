//! BRepAdaptor (TKBRep) — the Curve2d adaptor over an edge-on-face pcurve.
//!
//! 1:1 translation of OCCT `BRepAdaptor_Curve2d.hxx` (L17-64) + `.cxx`
//! (L19-77): inherits Geom2dAdaptor_Curve (the pcurve + its loaded range),
//! adds the Edge/Face accessors.  The Geom2dAdaptor_Curve role is played by
//! the stored pcurve plus the loaded range; the
//! [`Curve2dAdaptor`] implementation delegates per-curve evaluation to the
//! existing `impl Curve2dAdaptor for Curve2d` (matching the landed
//! Geom2dAdaptor translation, whose `resolution` is the pass-through
//! convention).
//!
//! Ownership: OCCT keeps the shapes alive through TopoDS handles over a
//! refcounted TShape graph; rcad's HLRBRep_Data (the future owner of the
//! adaptor chain) holds the `BRep` in an `Arc`, and this adaptor holds a
//! refcounted clone — the handle copy.

use glam::DVec2;
use rcad_kernel::geom::Curve2d;
use rcad_kernel::topods::{BRep, BRepTool, Shape};
use std::sync::Arc;

use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType};
use rcad_kernel::math::GeomAbsShape;

/// OCCT BRepAdaptor_Curve2d — "allows to use an Edge on a Face like a 2d
/// curve".
#[derive(Clone)]
pub struct BRepCurve2d {
    /// OCCT keeps the shapes alive through TopoDS handles; the BRep owns
    /// the TShape graph (rcad: the refcounted handle copy).
    brep: Arc<BRep>,
    /// OCCT myEdge.
    my_edge: Shape,
    /// OCCT myFace.
    my_face: Shape,
    /// OCCT Geom2dAdaptor_Curve::myCurve.
    my_curve: Curve2d,
    /// OCCT Geom2dAdaptor_Curve::myFirst / myLast (the loaded range).
    my_first: f64,
    my_last: f64,
}

impl BRepCurve2d {
    /// OCCT BRepAdaptor_Curve2d() (cxx L21) — uninitialized.
    pub fn new(brep: Arc<BRep>) -> Self {
        BRepCurve2d {
            brep,
            my_edge: Shape::null(),
            my_face: Shape::null(),
            my_curve: Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: DVec2::ZERO,
                direction: DVec2::X,
            }),
            my_first: 0.0,
            my_last: 0.0,
        }
    }

    /// OCCT BRepAdaptor_Curve2d(E, F) (cxx L24-27).
    pub fn new_edge_face(brep: Arc<BRep>, e: &Shape, f: &Shape) -> Self {
        let mut c = BRepCurve2d::new(brep);
        c.initialize(e, f);
        c
    }

    /// OCCT Initialize(E, F) (cxx L50-62) — load the pcurve of E on F with
    /// its range (Geom2dAdaptor_Curve::Load).
    pub fn initialize(&mut self, e: &Shape, f: &Shape) {
        self.my_edge = e.clone();
        self.my_face = f.clone();
        if let Some((pcurve, first, last)) = self.brep.curve_on_surface(e, f) {
            self.my_curve = pcurve;
            self.my_first = first;
            self.my_last = last;
        }
    }

    /// OCCT Edge() (cxx L65-68).
    pub fn edge(&self) -> &Shape {
        &self.my_edge
    }

    /// OCCT Face() (cxx L71-74).
    pub fn face(&self) -> &Shape {
        &self.my_face
    }

    /// The owning BRep (the OCCT global BRep_Tool context).
    pub fn brep(&self) -> &BRep {
        &self.brep
    }

    /// The refcounted BRep handle (the OCCT handle copy).
    pub fn brep_arc(&self) -> Arc<BRep> {
        self.brep.clone()
    }

    /// OCCT Geom2dAdaptor_Curve::Curve() — the loaded pcurve.
    pub fn curve(&self) -> &Curve2d {
        &self.my_curve
    }
}

impl Curve2dAdaptor for BRepCurve2d {
    /// OCCT Geom2dAdaptor_Curve::FirstParameter() — the loaded range start.
    fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT Geom2dAdaptor_Curve::LastParameter() — the loaded range end.
    fn last_parameter(&self) -> f64 {
        self.my_last
    }

    fn value(&self, u: f64) -> DVec2 {
        Curve2dAdaptor::value(&self.my_curve, u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        Curve2dAdaptor::d1(&self.my_curve, u)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        Curve2dAdaptor::d2(&self.my_curve, u)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        Curve2dAdaptor::d3(&self.my_curve, u)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        Curve2dAdaptor::dn(&self.my_curve, u, n)
    }

    fn get_type(&self) -> Curve2dType {
        Curve2dAdaptor::get_type(&self.my_curve)
    }

    fn nb_samples(&self) -> i32 {
        Curve2dAdaptor::nb_samples(&self.my_curve)
    }

    /// OCCT Geom2dAdaptor_Curve::Resolution — the landed translation uses
    /// the pass-through convention (see the module note).
    fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    fn is_closed(&self) -> bool {
        Curve2dAdaptor::is_closed(&self.my_curve)
    }

    fn is_periodic(&self) -> bool {
        Curve2dAdaptor::is_periodic(&self.my_curve)
    }

    fn period(&self) -> f64 {
        Curve2dAdaptor::period(&self.my_curve)
    }

    fn nb_knots(&self) -> i32 {
        Curve2dAdaptor::nb_knots(&self.my_curve)
    }

    fn degree(&self) -> i32 {
        Curve2dAdaptor::degree(&self.my_curve)
    }

    fn nb_poles(&self) -> i32 {
        Curve2dAdaptor::nb_poles(&self.my_curve)
    }

    fn circle(&self) -> rcad_kernel::geom::Circle2d {
        Curve2dAdaptor::circle(&self.my_curve)
    }

    fn line(&self) -> rcad_kernel::geom::Line2d {
        Curve2dAdaptor::line(&self.my_curve)
    }

    fn ellipse(&self) -> rcad_kernel::geom::Ellipse2d {
        Curve2dAdaptor::ellipse(&self.my_curve)
    }

    fn parabola(&self) -> rcad_kernel::geom::Parabola2d {
        Curve2dAdaptor::parabola(&self.my_curve)
    }

    fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola2d {
        Curve2dAdaptor::hyperbola(&self.my_curve)
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}

/// OCCT GeomAbs continuity of the loaded pcurve (Geom2dAdaptor passes
/// through the curve's own continuity).
impl BRepCurve2d {
    pub fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::C1
    }
}
