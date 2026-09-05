// OCCT Contap_ThePathPointOfTheSearch (TKHLR) — a starting point found on
// the surface domain by Contap_TheSearch.
//
// Contap_ThePathPointOfTheSearch.hxx L27-141 + _0.cxx L24-57.  The OCCT
// `occ::handle<Adaptor2d_Curve2d> arc` maps to [`Arc`], the vertex to the
// polymorphic HVertexHandle.

use rcad_kernel::geom::Point3;

use crate::topalgo::adaptor3d::hvertex::{HVertex, HVertexHandle};

use super::point::Arc;

/// OCCT Contap_ThePathPointOfTheSearch.
#[derive(Clone)]
pub struct ThePathPointOfTheSearch {
    point: Point3,
    tol: f64,
    isnew: bool,
    vtx: Option<HVertexHandle>,
    arc: Option<Arc>,
    param: f64,
}

impl ThePathPointOfTheSearch {
    /// OCCT Contap_ThePathPointOfTheSearch() (_0.cxx L24-29).
    pub fn new() -> Self {
        ThePathPointOfTheSearch {
            point: Point3::ZERO,
            tol: 0.0,
            isnew: true,
            vtx: None,
            arc: None,
            param: 0.0,
        }
    }

    /// OCCT Contap_ThePathPointOfTheSearch(P, Tol, V, A, Parameter)
    /// (_0.cxx L31-44).
    pub fn with_vertex(p: Point3, tol: f64, v: HVertexHandle, a: Arc, parameter: f64) -> Self {
        ThePathPointOfTheSearch {
            point: p,
            tol,
            isnew: false,
            vtx: Some(v),
            arc: Some(a),
            param: parameter,
        }
    }

    /// OCCT Contap_ThePathPointOfTheSearch(P, Tol, A, Parameter)
    /// (_0.cxx L46-57).
    pub fn without_vertex(p: Point3, tol: f64, a: Arc, parameter: f64) -> Self {
        ThePathPointOfTheSearch {
            point: p,
            tol,
            isnew: true,
            vtx: None,
            arc: Some(a),
            param: parameter,
        }
    }

    /// OCCT SetValue(P, Tol, V, A, Parameter) (hxx L81-93).
    pub fn set_value_vertex(
        &mut self,
        p: Point3,
        tol: f64,
        v: HVertexHandle,
        a: Arc,
        parameter: f64,
    ) {
        self.isnew = false;
        self.point = p;
        self.tol = tol;
        self.vtx = Some(v);
        self.arc = Some(a);
        self.param = parameter;
    }

    /// OCCT SetValue(P, Tol, A, Parameter) (hxx L95-105).
    pub fn set_value(&mut self, p: Point3, tol: f64, a: Arc, parameter: f64) {
        self.isnew = true;
        self.point = p;
        self.tol = tol;
        self.arc = Some(a);
        self.param = parameter;
    }

    /// OCCT Value (hxx L107-110).
    pub fn value(&self) -> Point3 {
        self.point
    }

    /// OCCT Tolerance (hxx L112-115).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT IsNew (hxx L117-120).
    pub fn is_new(&self) -> bool {
        self.isnew
    }

    /// OCCT Vertex (hxx L122-129) — raises Standard_DomainError when new.
    pub fn vertex(&self) -> &HVertexHandle {
        if self.isnew {
            panic!("Standard_DomainError: ThePathPointOfTheSearch::Vertex");
        }
        self.vtx.as_ref().unwrap()
    }

    /// OCCT Arc (hxx L131-134).
    pub fn arc(&self) -> Arc {
        self.arc.clone().expect("ThePathPointOfTheSearch::Arc")
    }

    /// OCCT Parameter (hxx L136-139).
    pub fn parameter(&self) -> f64 {
        self.param
    }
}

impl Default for ThePathPointOfTheSearch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::Curve2d;
    use rcad_kernel::topo::topods::Orientation;

    fn test_arc() -> Arc {
        std::sync::Arc::new(Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 0.0),
        }))
    }

    /// OCCT anchor: the vertex variant carries isnew=false, the plain
    /// variant isnew=true; Vertex() raises DomainError when new.
    #[test]
    fn path_point_vertex_semantics() {
        let vtx: HVertexHandle = std::sync::Arc::new(HVertex::new_with(
            DVec2::ZERO,
            Orientation::Forward,
            1e-7,
        ));
        let p = ThePathPointOfTheSearch::with_vertex(DVec3::ONE, 1e-5, vtx, test_arc(), 0.5);
        assert!(!p.is_new());
        assert!((p.parameter() - 0.5).abs() < 1e-15);
        assert!((p.tolerance() - 1e-5).abs() < 1e-15);

        let q = ThePathPointOfTheSearch::without_vertex(DVec3::ONE, 1e-5, test_arc(), 0.5);
        assert!(q.is_new());
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| q.vertex().clone())).is_err();
        assert!(panicked);
    }
}
