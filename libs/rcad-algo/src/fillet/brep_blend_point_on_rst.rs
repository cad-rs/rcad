//! OCCT BRepBlend_PointOnRst (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_PointOnRst.hxx (L26-75) + BRepBlend_PointOnRst.lxx (L9-35).
//! The .cxx contains only the default constructor and SetArc
//! (BRepBlend_PointOnRst.cxx L19-38).
//!
//! Architecture mapping: `occ::handle<Adaptor2d_Curve2d>` maps to a cloned
//! rcad `Curve2d` value (rcad has no shared-handle adaptor layer for 2d
//! pcurves); `IntSurf_Transition` maps to
//! [`IntSurfTransition`](crate::geomalgo::int_patch::transitions::Transition).

use rcad_kernel::geom::Curve2d;

use crate::geomalgo::int_patch::transitions::Transition as IntSurfTransition;

/// OCCT BRepBlend_PointOnRst — definition of an imposed point on a
/// restriction (BRepBlend_PointOnRst.hxx L20).
#[derive(Debug, Clone)]
pub struct BRepBlendPointOnRst {
    arc: Option<Curve2d>,
    traline: IntSurfTransition,
    traarc: IntSurfTransition,
    prm: f64,
}

impl BRepBlendPointOnRst {
    /// OCCT BRepBlend_PointOnRst() (BRepBlend_PointOnRst.cxx L19-25) —
    /// empty constructor.
    pub fn new() -> Self {
        BRepBlendPointOnRst {
            arc: None,
            traline: IntSurfTransition::new(),
            traarc: IntSurfTransition::new(),
            prm: 0.0,
        }
    }

    /// OCCT BRepBlend_PointOnRst(A, Param, TLine, TArc) (L27-37).
    pub fn new_with_arc(
        a: &Curve2d,
        param: f64,
        t_line: IntSurfTransition,
        t_arc: IntSurfTransition,
    ) -> Self {
        BRepBlendPointOnRst {
            arc: Some(a.clone()),
            traline: t_line,
            traarc: t_arc,
            prm: param,
        }
    }

    /// OCCT SetArc(A, Param, TLine, TArc) (BRepBlend_PointOnRst.cxx L39-48).
    pub fn set_arc(
        &mut self,
        a: &Curve2d,
        param: f64,
        t_line: IntSurfTransition,
        t_arc: IntSurfTransition,
    ) {
        self.arc = Some(a.clone());
        self.traline = t_line;
        self.traarc = t_arc;
        self.prm = param;
    }

    /// OCCT Arc() (BRepBlend_PointOnRst.lxx L9-12).
    pub fn arc(&self) -> Option<&Curve2d> {
        self.arc.as_ref()
    }

    /// OCCT TransitionOnLine() (BRepBlend_PointOnRst.lxx L14-17).
    pub fn transition_on_line(&self) -> &IntSurfTransition {
        &self.traline
    }

    /// OCCT TransitionOnArc() (BRepBlend_PointOnRst.lxx L19-22).
    pub fn transition_on_arc(&self) -> &IntSurfTransition {
        &self.traarc
    }

    /// OCCT ParameterOnArc() (BRepBlend_PointOnRst.lxx L24-27).
    pub fn parameter_on_arc(&self) -> f64 {
        self.prm
    }
}
