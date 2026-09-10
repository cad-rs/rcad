// OCCT `occ::handle<Adaptor3d_TopolTool>` — the domain of the contour
// engine.
//
// The Contap engines (IntStart_SearchOnBoundaries.gxx,
// IntStart_SearchInside.gxx, Contap_Contour.cxx) call the TopolTool
// through its abstract interface; this object-safe trait is that
// interface, and the generic Adaptor3d TopolTool
// (`topalgo::adaptor3d::topol_tool::TopolTool`) implements it, as does the
// BRepTopAdaptor_TopolTool (`topalgo::brep_top_adaptor::topol_tool_brep`).

use glam::DVec2;
use rcad_kernel::geom::Point3;
use rcad_kernel::topods::{Orientation, Shape, State};

use crate::hlr::contap::point::Arc;
use crate::topalgo::adaptor3d::hvertex::{HVertex, HVertexHandle};
use crate::topalgo::adaptor3d::topol_tool::TopolTool;
use crate::topalgo::adaptor2d::line2d::Line2dAdaptor;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;
use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;
use crate::geomalgo::int_curve_surface::HSurfaceTool;

/// The OCCT Adaptor3d_TopolTool interface as used by the Contap engines.
pub trait ContapDomain {
    /// OCCT Init().
    fn init(&mut self);
    /// OCCT More().
    fn more(&self) -> bool;
    /// OCCT Value() — the current restriction arc.
    fn arc(&self) -> Arc;
    /// OCCT Next().
    fn next(&mut self);
    /// OCCT Initialize(A) — attach the vertex iterator to the arc.
    fn initialize_arc(&mut self, a: &Arc);
    /// OCCT InitVertexIterator().
    fn init_vertex_iterator(&mut self);
    /// OCCT MoreVertex().
    fn more_vertex(&self) -> bool;
    /// OCCT Vertex() — the current arc-end vertex.
    fn vertex(&self) -> HVertexHandle;
    /// OCCT NextVertex().
    fn next_vertex(&mut self);
    /// OCCT Identical(V1, V2).
    fn identical(&self, v1: &HVertexHandle, v2: &HVertexHandle) -> bool;
    /// OCCT Classify(P, Tol, RecadreOnPeriodic).
    fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State;
    /// OCCT Orientation(A).
    fn orientation_arc(&self, a: &Arc) -> Orientation;
    /// OCCT Orientation(V).
    fn orientation_vertex(&self, v: &HVertexHandle) -> Orientation;
    /// OCCT Edge() — the BRep edge address when the domain carries real
    /// edges (BRepTopAdaptor_TopolTool); the base Adaptor3d version has
    /// no edge and returns None.
    fn edge(&self) -> Option<Shape>;
    /// OCCT NbSamplesU().
    fn nb_samples_u(&mut self) -> i32;
    /// OCCT NbSamplesV().
    fn nb_samples_v(&mut self) -> i32;
    /// OCCT NbSamples().
    fn nb_samples(&mut self) -> i32;
    /// OCCT SamplePoint(I, uv, uv3d) — 1-based.
    fn sample_point(&self, i: usize) -> (DVec2, Point3);
}

impl<'a, S, ST> ContapDomain for TopolTool<'a, S, ST>
where
    ST: HSurfaceTool<Surface = S>,
{
    fn init(&mut self) {
        TopolTool::init(self)
    }
    fn more(&self) -> bool {
        TopolTool::more(self)
    }
    fn arc(&self) -> Arc {
        // OCCT returns the handle to the restriction adaptor; the rcad
        // restrictions are Line2dAdaptor values — clone into the shared
        // handle (the OCCT handle copy).
        let a: Line2dAdaptor = *TopolTool::value(self);
        std::sync::Arc::new(a)
    }
    fn next(&mut self) {
        TopolTool::next(self)
    }
    fn initialize_arc(&mut self, a: &Arc) {
        TopolTool::initialize_curve(self, a.as_ref())
    }
    fn init_vertex_iterator(&mut self) {
        TopolTool::init_vertex_iterator(self)
    }
    fn more_vertex(&self) -> bool {
        TopolTool::more_vertex(self)
    }
    fn vertex(&self) -> HVertexHandle {
        std::sync::Arc::new(*TopolTool::vertex(self))
    }
    fn next_vertex(&mut self) {
        TopolTool::next_vertex(self)
    }
    fn identical(&self, v1: &HVertexHandle, v2: &HVertexHandle) -> bool {
        // OCCT Adaptor3d_TopolTool::Identical(V1, V2) = V1->IsSame(V2)
        // (Adaptor3d_TopolTool.cxx L616-619) — dispatched on the handles.
        v1.is_same(v2.as_ref())
    }
    fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        TopolTool::classify(self, p, tol, recadre_on_periodic)
    }
    fn orientation_arc(&self, a: &Arc) -> Orientation {
        TopolTool::orientation_curve(self, a.as_ref())
    }
    fn orientation_vertex(&self, v: &HVertexHandle) -> Orientation {
        // OCCT Adaptor3d_TopolTool::Orientation(V) = V->Orientation().
        v.orientation()
    }
    fn edge(&self) -> Option<Shape> {
        TopolTool::edge(self).map(|_| Shape::null())
    }
    fn nb_samples_u(&mut self) -> i32 {
        TopolTool::nb_samples_u(self)
    }
    fn nb_samples_v(&mut self) -> i32 {
        TopolTool::nb_samples_v(self)
    }
    fn nb_samples(&mut self) -> i32 {
        TopolTool::nb_samples(self)
    }
    fn sample_point(&self, i: usize) -> (DVec2, Point3) {
        TopolTool::sample_point(self, i)
    }
}

impl ContapDomain for BRepTopolTool {
    fn init(&mut self) {
        BRepTopolTool::init(self)
    }
    fn more(&self) -> bool {
        BRepTopolTool::more(self)
    }
    fn arc(&self) -> Arc {
        // OCCT Value() = down_cast<Adaptor2d_Curve2d>(myCIterator.Value());
        // the rcad arcs are BRepCurve2d clones in the shared handle.
        std::sync::Arc::new(BRepTopolTool::value(self).unwrap().clone())
    }
    fn next(&mut self) {
        BRepTopolTool::next(self)
    }
    fn initialize_arc(&mut self, a: &Arc) {
        // OCCT Initialize(C): myCurve = down_cast<BRepAdaptor_Curve2d>(C)
        // (BRepTopAdaptor_TopolTool.cxx L98-105).
        let c = a
            .as_any()
            .downcast_ref::<BRepCurve2d>()
            .expect("Standard_ConstructionError");
        BRepTopolTool::initialize_curve(self, c);
    }
    fn init_vertex_iterator(&mut self) {
        BRepTopolTool::init_vertex_iterator(self)
    }
    fn more_vertex(&self) -> bool {
        BRepTopolTool::more_vertex(self)
    }
    fn vertex(&self) -> HVertexHandle {
        // OCCT Vertex() = new BRepTopAdaptor_HVertex(V, myCurve)
        // (cxx L167-170); the rcad vertex owns its adaptor copy.
        std::sync::Arc::new(BRepTopolTool::vertex(self).unwrap())
    }
    fn next_vertex(&mut self) {
        BRepTopolTool::next_vertex(self)
    }
    fn identical(&self, v1: &HVertexHandle, v2: &HVertexHandle) -> bool {
        // OCCT does not override Identical — base V1->IsSame(V2)
        // dispatched to BRepTopAdaptor_HVertex::IsSame.
        v1.is_same(v2.as_ref())
    }
    fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        BRepTopolTool::classify(self, p, tol, recadre_on_periodic)
    }
    fn orientation_arc(&self, a: &Arc) -> Orientation {
        // OCCT Orientation(C): down_cast + brhc->Edge().Orientation()
        // (cxx L213-217).
        let c = a
            .as_any()
            .downcast_ref::<BRepCurve2d>()
            .expect("Standard_ConstructionError");
        BRepTopolTool::orientation_curve(self, c)
    }
    fn orientation_vertex(&self, v: &HVertexHandle) -> Orientation {
        // OCCT delegates to the base: V->Orientation() (cxx L221-224).
        v.orientation()
    }
    fn edge(&self) -> Option<Shape> {
        // OCCT Edge() = down_cast(myCIterator.Value())->Edge()
        // (cxx L138-143).
        BRepTopolTool::edge(self).map(|e| e.clone())
    }
    fn nb_samples_u(&mut self) -> i32 {
        BRepTopolTool::nb_samples_u(self)
    }
    fn nb_samples_v(&mut self) -> i32 {
        BRepTopolTool::nb_samples_v(self)
    }
    fn nb_samples(&mut self) -> i32 {
        BRepTopolTool::nb_samples(self)
    }
    fn sample_point(&self, i: usize) -> (DVec2, Point3) {
        BRepTopolTool::sample_point(self, i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::contap::h_cont_tool;
    use crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face;
    use glam::DVec3;

    /// OCCT anchor: the BRepTopAdaptor_TopolTool domain driven through the
    /// Contap engines' abstract interface — the arc iteration hands out
    /// down-castable pcurve handles, the vertex iterator yields
    /// BRepTopAdaptor_HVertex handles whose Parameter is
    /// BRep_Tool::Parameter (h_cont_tool::parameter dispatch), Identical is
    /// the TopoDS IsSame, Classify delegates to the FClass2d, and Edge()
    /// returns the current arc's edge (BRepTopAdaptor_TopolTool.cxx
    /// L98-224).
    #[test]
    fn contap_domain_brep_topol_tool_anchor() {
        let (brep, face) = square_face();
        let mut domain = BRepTopolTool::new_with_surface(std::sync::Arc::new(brep), &face);
        let domain: &mut dyn ContapDomain = &mut domain;

        domain.init();
        let mut arcs = Vec::new();
        let mut count = 0;
        while domain.more() {
            arcs.push(domain.arc());
            count += 1;
            domain.next();
        }
        assert_eq!(count, 4);

        // Engine order (SearchOnBoundaries): rewind to arc 0, then the
        // per-arc vertex iteration — Edge() reads myCIterator.Value() and
        // is valid inside the loop (cxx L138-143).
        domain.init();
        assert!(domain.more());
        let a0 = domain.arc();
        domain.initialize_arc(&a0);
        domain.init_vertex_iterator();
        assert!(domain.more_vertex());
        let v1 = domain.vertex();
        domain.next_vertex();
        assert!(domain.more_vertex());
        let v2 = domain.vertex();
        // BRepTopAdaptor_HVertex::Parameter = BRep_Tool::Parameter — the
        // first vertex of the bottom edge sits at the pcurve start.
        let p1 = h_cont_tool::parameter(v1.as_ref(), a0.as_ref());
        assert!((p1 - 0.0).abs() < 1e-12, "p1={}", p1);
        let p2 = h_cont_tool::parameter(v2.as_ref(), a0.as_ref());
        assert!((p2 - 1.0).abs() < 1e-12, "p2={}", p2);
        // Identical(V,V) topological IsSame; the two edge ends differ.
        assert!(domain.identical(&v1, &v1));
        assert!(!domain.identical(&v1, &v2));
        assert_eq!(v1.orientation(), Orientation::Forward);
        assert_eq!(domain.orientation_vertex(&v1), Orientation::Forward);
        assert_eq!(domain.orientation_vertex(&v2), Orientation::Reversed);
        domain.next_vertex();
        assert!(!domain.more_vertex());

        // Orientation(A) = the edge orientation (cxx L213-217).
        assert_eq!(domain.orientation_arc(&a0), Orientation::Forward);

        // Edge() = the current arc's edge (cxx L138-143).
        let e = domain.edge().expect("the BRep domain carries edges");
        assert!(e.is_same(a0
            .as_any()
            .downcast_ref::<BRepCurve2d>()
            .unwrap()
            .edge()));

        // Classify delegates to the FClass2d (cxx L174-187).
        assert_eq!(domain.classify(DVec2::new(0.5, 0.5), 1e-6, true), State::In);
        assert_eq!(domain.classify(DVec2::new(3.0, 3.0), 1e-6, true), State::Out);

        // Sample grid: plane clamps to 10x10 (cxx L396-472 / L558-566).
        assert_eq!(domain.nb_samples(), 100);
        let (uv, p3d) = domain.sample_point(0);
        assert!(uv.x > 0.0 && uv.y > 0.0);
        assert_eq!(p3d, DVec3::new(uv.x, uv.y, 0.0));
    }

    /// OCCT anchor: the base Adaptor3d domain vertex handle dispatches
    /// Value/Parameter on the shared trait (Adaptor3d_HVertex semantics
    /// preserved after the handle conversion).
    #[test]
    fn contap_domain_base_vertex_handle_anchor() {
        let v: HVertexHandle = std::sync::Arc::new(HVertex::new_with(
            DVec2::new(3.0, 2.0),
            Orientation::Reversed,
            1.0e-8,
        ));
        assert_eq!(v.value(), DVec2::new(3.0, 2.0));
        assert_eq!(v.orientation(), Orientation::Reversed);
        let other: HVertexHandle = std::sync::Arc::new(HVertex::new_with(
            DVec2::new(3.0, 2.0 + 1e-9),
            Orientation::Forward,
            1e-8,
        ));
        assert!(v.is_same(other.as_ref()));
        assert!(v.topo_vertex().is_none());
    }
}
