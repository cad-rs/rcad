// OCCT `occ::handle<Adaptor3d_TopolTool>` — the domain of the contour
// engine.
//
// The Contap engines (IntStart_SearchOnBoundaries.gxx,
// IntStart_SearchInside.gxx, Contap_Contour.cxx) call the TopolTool
// through its abstract interface; this object-safe trait is that
// interface, and the generic Adaptor3d TopolTool
// (`topalgo::adaptor3d::topol_tool::TopolTool`) implements it.
// The BRepTopAdaptor_TopolTool (real face wires) will implement it in
// Stage 3b/3g with the same contract.

use glam::DVec2;
use rcad_kernel::geom::Point3;
use rcad_kernel::topo::topods::{Orientation, State};

use crate::hlr::contap::point::Arc;
use crate::topalgo::adaptor3d::hvertex::HVertex;
use crate::topalgo::adaptor3d::topol_tool::TopolTool;
use crate::topalgo::adaptor2d::line2d::Line2dAdaptor;
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
    fn vertex(&self) -> HVertex;
    /// OCCT NextVertex().
    fn next_vertex(&mut self);
    /// OCCT Identical(V1, V2).
    fn identical(&self, v1: &HVertex, v2: &HVertex) -> bool;
    /// OCCT Classify(P, Tol, RecadreOnPeriodic).
    fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State;
    /// OCCT Orientation(A).
    fn orientation_arc(&self, a: &Arc) -> Orientation;
    /// OCCT Orientation(V).
    fn orientation_vertex(&self, v: &HVertex) -> Orientation;
    /// OCCT Edge() — the BRep edge address when the domain carries real
    /// edges (BRepTopAdaptor_TopolTool); the base Adaptor3d version has
    /// no edge and returns None.
    fn edge(&self) -> Option<()>;
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
    fn vertex(&self) -> HVertex {
        *TopolTool::vertex(self)
    }
    fn next_vertex(&mut self) {
        TopolTool::next_vertex(self)
    }
    fn identical(&self, v1: &HVertex, v2: &HVertex) -> bool {
        TopolTool::identical(self, v1, v2)
    }
    fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        TopolTool::classify(self, p, tol, recadre_on_periodic)
    }
    fn orientation_arc(&self, a: &Arc) -> Orientation {
        TopolTool::orientation_curve(self, a.as_ref())
    }
    fn orientation_vertex(&self, v: &HVertex) -> Orientation {
        TopolTool::orientation_vertex(self, v)
    }
    fn edge(&self) -> Option<()> {
        TopolTool::edge(self)
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
