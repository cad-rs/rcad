// OCCT Contap_Point (TKHLR) — a vertex on the contour line.
//
// Contap_Point.hxx L35-131 + Contap_Point.cxx L21-44 + Contap_Point.lxx
// L19-153.  The OCCT handle members `arc` (handle<Adaptor2d_Curve2d>) and
// `vtx` (handle<Adaptor3d_HVertex>) become Arc + the Adaptor3d HVertex
// value type; the Standard_DomainError guards map to panics with the same
// control flow.

use glam::DVec3;
use rcad_kernel::topo::topods::Orientation;

use crate::geomalgo::int_patch::transitions::Transition;
use crate::topalgo::adaptor3d::hvertex::HVertex;

/// OCCT `occ::handle<Adaptor2d_Curve2d>` as used by Contap data classes.
pub type Arc = std::sync::Arc<dyn crate::geomalgo::geom2d_int::Curve2dAdaptor>;

/// OCCT Contap_Point.
#[derive(Clone)]
pub struct Point {
    pt: DVec3,
    uparam: f64,
    vparam: f64,
    paraline: f64,
    onarc: bool,
    arc: Option<Arc>,
    traline: Transition,
    traarc: Transition,
    prmarc: f64,
    isvtx: bool,
    vtx: Option<HVertex>,
    ismult: bool,
    my_internal: bool,
}

impl Point {
    /// OCCT Contap_Point() (cxx L21-31) — empty constructor.
    pub fn new() -> Self {
        Point {
            pt: DVec3::ZERO,
            uparam: 0.0,
            vparam: 0.0,
            paraline: 0.0,
            onarc: false,
            arc: None,
            traline: Transition::new(),
            traarc: Transition::new(),
            prmarc: 0.0,
            isvtx: false,
            vtx: None,
            ismult: false,
            my_internal: false,
        }
    }

    /// OCCT Contap_Point(Pt, U, V) (cxx L33-44).
    pub fn with_uv(pt: DVec3, u: f64, v: f64) -> Self {
        Point {
            pt,
            uparam: u,
            vparam: v,
            paraline: 0.0,
            onarc: false,
            arc: None,
            traline: Transition::new(),
            traarc: Transition::new(),
            prmarc: 0.0,
            isvtx: false,
            vtx: None,
            ismult: false,
            my_internal: false,
        }
    }

    /// OCCT SetValue(Pt, U, V) (lxx L19-28).
    pub fn set_value(&mut self, pt: DVec3, u: f64, v: f64) {
        self.pt = pt;
        self.uparam = u;
        self.vparam = v;
        self.onarc = false;
        self.isvtx = false;
        self.ismult = false;
        self.my_internal = false;
    }

    /// OCCT SetParameter(Para) (lxx L30-34).
    pub fn set_parameter(&mut self, para: f64) {
        self.paraline = para;
    }

    /// OCCT SetVertex(V) (lxx L36-41).
    pub fn set_vertex(&mut self, v: HVertex) {
        self.isvtx = true;
        self.vtx = Some(v);
    }

    /// OCCT SetArc(A, Param, TLine, TArc) (lxx L43-54).
    pub fn set_arc(&mut self, a: Arc, param: f64, tline: Transition, tarc: Transition) {
        self.onarc = true;
        self.arc = Some(a);
        self.prmarc = param;
        self.traline = tline;
        self.traarc = tarc;
    }

    /// OCCT SetMultiple (lxx L56-59).
    pub fn set_multiple(&mut self) {
        self.ismult = true;
    }

    /// OCCT SetInternal (lxx L61-64).
    pub fn set_internal(&mut self) {
        self.my_internal = true;
    }

    /// OCCT IsMultiple (lxx L66-69).
    pub fn is_multiple(&self) -> bool {
        self.ismult
    }

    /// OCCT IsInternal (lxx L71-74).
    pub fn is_internal(&self) -> bool {
        self.my_internal
    }

    /// OCCT Value (lxx L76-80).
    pub fn value(&self) -> DVec3 {
        self.pt
    }

    /// OCCT ParameterOnLine (lxx L82-86).
    pub fn parameter_on_line(&self) -> f64 {
        self.paraline
    }

    /// OCCT Parameters(U1, V1) (lxx L88-93).
    pub fn parameters(&self) -> (f64, f64) {
        (self.uparam, self.vparam)
    }

    /// OCCT IsOnArc (lxx L95-98).
    pub fn is_on_arc(&self) -> bool {
        self.onarc
    }

    /// OCCT Arc (lxx L100-108) — raises Standard_DomainError when !onarc.
    pub fn arc(&self) -> Arc {
        if !self.onarc {
            panic!("Standard_DomainError: Contap_Point::Arc");
        }
        self.arc.clone().unwrap()
    }

    /// OCCT TransitionOnLine (lxx L110-118).
    pub fn transition_on_line(&self) -> Transition {
        if !self.onarc {
            panic!("Standard_DomainError: Contap_Point::TransitionOnLine");
        }
        self.traline
    }

    /// OCCT TransitionOnArc (lxx L120-128).
    pub fn transition_on_arc(&self) -> Transition {
        if !self.onarc {
            panic!("Standard_DomainError: Contap_Point::TransitionOnArc");
        }
        self.traarc
    }

    /// OCCT ParameterOnArc (lxx L130-138).
    pub fn parameter_on_arc(&self) -> f64 {
        if !self.onarc {
            panic!("Standard_DomainError: Contap_Point::ParameterOnArc");
        }
        self.prmarc
    }

    /// OCCT IsVertex (lxx L140-143).
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// OCCT Vertex (lxx L145-153).
    pub fn vertex(&self) -> &HVertex {
        if !self.isvtx {
            panic!("Standard_DomainError: Contap_Point::Vertex");
        }
        self.vtx.as_ref().unwrap()
    }
}

impl Default for Point {
    fn default() -> Self {
        Self::new()
    }
}

// The Orientation import is part of the HVertex value shape carried by
// Contap_Point::Vertex (Adaptor3d_HVertex::Orientation).
#[allow(unused)]
fn _orientation_shape(o: Orientation) -> Orientation {
    o
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec2;

    /// OCCT anchor: the default point carries no arc/vertex; SetValue
    /// resets the onarc/isvtx/ismult/internal flags; SetArc + accessors
    /// round-trip (lxx semantics).
    #[test]
    fn contap_point_lifecycle() {
        let mut p = Point::new();
        assert!(!p.is_on_arc());
        assert!(!p.is_vertex());
        assert!(!p.is_multiple());
        assert!(!p.is_internal());

        p.set_value(DVec3::new(1.0, 2.0, 3.0), 0.5, 0.25);
        assert_eq!(p.parameters(), (0.5, 0.25));
        assert_eq!(p.value(), DVec3::new(1.0, 2.0, 3.0));

        p.set_parameter(7.0);
        assert!((p.parameter_on_line() - 7.0).abs() < 1e-15);

        p.set_multiple();
        p.set_internal();
        assert!(p.is_multiple() && p.is_internal());

        // The DomainError path: Arc() without SetArc panics.
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let q = Point::new();
            let _ = q.arc();
        }))
        .is_err();
        assert!(panicked);
    }

    /// OCCT anchor: after SetValue, the previously set vertex flag is
    /// cleared (lxx L19-28 clears isvtx/ismult/myInternal).
    #[test]
    fn contap_point_set_value_clears_flags() {
        let mut p = Point::new();
        p.set_vertex(HVertex::new_with(DVec2::ZERO, Orientation::Forward, 1e-7));
        assert!(p.is_vertex());
        p.set_value(DVec3::ZERO, 0.0, 0.0);
        assert!(!p.is_vertex());
    }
}
