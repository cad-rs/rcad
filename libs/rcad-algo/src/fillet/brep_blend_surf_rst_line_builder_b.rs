//! OCCT BRepBlend bridge layer for the two LineBuilders
//! (BRepBlend_SurfRstLineBuilder / BRepBlend_RstRstLineBuilder) — the
//! architecture-difference support that OCCT reaches through
//! `occ::handle` adaptors.  Split from brep_blend_surf_rst_line_builder.rs
//! per the 2000-line guideline; consumed by both LineBuilder modules.
//!
//! Contents (OCCT anchors inside):
//!   - [`RstArc`] — `occ::handle<Adaptor2d_Curve2d>` (the restriction arc,
//!     the BRepAdaptor_Curve2d view).
//!   - [`DomainVertex`] — `occ::handle<Adaptor3d_HVertex>` over the BRep
//!     domain (BRepTopAdaptor_HVertex).
//!   - [`DomainTool`] — the Adaptor3d_TopolTool domain protocol
//!     (Init/More/Next/Value/Initialize/vertex iteration) as an explicit
//!     iterator value.
//!   - [`blend_tool_*`] — BRepBlend_BlendTool
//!     (BlendTool.cxx L38-142 + .lxx L19-58).

use glam::DVec2;

use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::geom::{Curve2d, Curve2dEval as _};
use rcad_kernel::topo::topods::{self, BRepTool as _};
use rcad_kernel::topods::{Orientation, Shape};

use crate::geomalgo::geom2d_int::{Curve2dAdaptor, TheIntPCurvePCurveOfGInter};
use crate::geomalgo::int_res2d::Domain as Res2dDomain;
use crate::topalgo::adaptor2d::line2d::Line2dAdaptor;
use crate::topalgo::adaptor3d::hvertex::HVertex;

use super::chfi3d_builder_0::{brep_tool_parameter, topexp_face_edges};
use super::chfi3d_builder_2::BRepTopAdaptorTopolTool;


/// OCCT `occ::handle<Adaptor2d_Curve2d>` as consumed by the LineBuilder —
/// the restriction arc: the pcurve value plus the owning edge/face identity
/// (the BRepAdaptor_Curve2d view used by the vertex scans and by
/// BRepTopAdaptor_TopolTool::Orientation, BRepTopAdaptor_TopolTool.cxx
/// L213-217).
#[derive(Debug, Clone)]
pub(crate) struct RstArc {
    /// The 2d pcurve (the BRepAdaptor_Curve2d evaluation).
    pub curve: Curve2d,
    /// The owning BRep of the edge scan (empty for a bare pcurve).
    pub brep: topods::BRep,
    /// The owning edge (BRepAdaptor_Curve2d::Edge()); null for a bare pcurve.
    pub edge: Shape,
    /// The owning face (BRepAdaptor_Curve2d::Face()).
    pub face: Shape,
}

impl RstArc {
    /// A bare pcurve arc (no edge identity) — the consumer-side `rst`.
    pub fn bare(curve: &Curve2d) -> Self {
        RstArc {
            curve: curve.clone(),
            brep: topods::BRep::default(),
            edge: Shape::null(),
            face: Shape::null(),
        }
    }

    /// The null handle before Recadre assigns the arc (a degenerate line
    /// placeholder — OCCT holds a null handle here).
    pub fn empty() -> Self {
        RstArc::bare(&Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        }))
    }

    /// OCCT `operator !=` on the handles — the identity of the underlying
    /// edge (the handle identity in OCCT).
    pub fn is_same(&self, other: &Self) -> bool {
        !self.edge.is_null() && !other.edge.is_null() && self.edge.ptr_id() == other.edge.ptr_id()
    }

    /// OCCT FirstParameter() (the trimmed range; the rcad natural domain is
    /// the carrier of the range for trimmed curves).
    pub fn first_parameter(&self) -> f64 {
        self.curve.default_domain()[0]
    }

    /// OCCT LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.curve.default_domain()[1]
    }

    /// OCCT Value(W).
    pub fn value(&self, w: f64) -> DVec2 {
        self.curve.point_at(w)
    }

    /// OCCT D1(W, P, DP).
    pub fn d1(&self, w: f64) -> (DVec2, DVec2) {
        (self.curve.point_at(w), self.curve.derivative_at(w))
    }

    /// OCCT BRepTopAdaptor_TopolTool::Orientation(C)
    /// (BRepTopAdaptor_TopolTool.cxx L213-217) — the edge orientation;
    /// Forward for a bare pcurve.
    pub fn orientation(&self) -> Orientation {
        if self.edge.is_null() {
            Orientation::Forward
        } else {
            self.edge.orientation
        }
    }
}

/// OCCT `occ::handle<Adaptor3d_HVertex>` over the BRep domain — the
/// Adaptor3d_HVertex view plus the topological vertex identity
/// (BRepTopAdaptor_HVertex::myVtx).
#[derive(Debug, Clone)]
pub(crate) struct DomainVertex {
    /// The Adaptor3d_HVertex view.  Value() is the BRepTopAdaptor_HVertex
    /// dummy (RealFirst, RealFirst) — BRepTopAdaptor_HVertex.cxx L36-39.
    pub hv: HVertex,
    /// The topological vertex (BRepTopAdaptor_HVertex.cxx L188-192 — the
    /// IsSame identity); null for the empty vertex.
    pub vtx: Shape,
}

impl DomainVertex {
    /// The null handle before Recadre assigns the vertex.
    pub fn empty() -> Self {
        DomainVertex {
            hv: HVertex::new(),
            vtx: Shape::null(),
        }
    }

    /// OCCT BRepTopAdaptor_HVertex(V, C) — the bridge construction used by
    /// the domain vertex iterator.  Orientation is the vertex orientation in
    /// the edge; the Resolution seed is BRep_Tool::Tolerance(V) (the full
    /// BRepTopAdaptor_HVertex::Resolution UV-refinement chain is pending —
    /// see the file header).
    pub fn from_topo(v: &Shape) -> Self {
        let vdata = v.as_vertex().expect("domain vertex");
        let resolution = vdata.tolerance;
        DomainVertex {
            hv: HVertex::new_with(DVec2::new(-f64::MAX, -f64::MAX), v.orientation, resolution),
            vtx: v.clone(),
        }
    }

    /// OCCT Adaptor3d_HVertex::IsSame via BRepTopAdaptor_HVertex::IsSame
    /// (cxx L188-192) — the topological vertex identity.
    pub fn is_same(&self, other: &Self) -> bool {
        !self.vtx.is_null() && !other.vtx.is_null() && self.vtx.ptr_id() == other.vtx.ptr_id()
    }
}

/// OCCT Adaptor3d_TopolTool / BRepTopAdaptor_TopolTool — the domain protocol
/// (Init / More / Next / Value, Initialize(curve), InitVertexIterator /
/// MoreVertex / NextVertex / Vertex, Identical) as consumed by the
/// LineBuilder.  Architecture difference: the OCCT tool mutates itself
/// through the protocol; the rcad carrier (`BRepTopAdaptorTopolTool` stub)
/// is immutable, so the iteration state is this explicit value, rebuilt at
/// each `->Init()` call site.
pub(crate) struct DomainTool {
    /// The domain arcs (myShapes sequence / myCurve single-curve mode).
    arcs: Vec<RstArc>,
    /// The arc iterator position (1-based; 0 = not initialized).
    iarc: i32,
    /// The current-arc vertex list (BRepTopAdaptor_TopolTool::myVIterator
    /// over myCurve->Edge() vertices, cxx L148-152).
    vertices: Vec<DomainVertex>,
    /// The vertex iterator position (1-based; 0 = not initialized).
    ivtx: i32,
}

impl DomainTool {
    /// OCCT BRepTopAdaptor_TopolTool over a face — the domain arcs are the
    /// face edges with their pcurves (Adaptor3d_TopolTool::Initialize(S) +
    /// Init over the restrictions).
    pub fn face_domain(tool: &BRepTopAdaptorTopolTool) -> Self {
        let mut arcs = Vec::new();
        if let Some(hs) = &tool.surface {
            for e in topexp_face_edges(&hs.brep, &hs.face) {
                if let Some((pc, _, _)) = hs.brep.curve_on_surface(&e, &hs.face) {
                    arcs.push(RstArc {
                        curve: pc,
                        brep: hs.brep.clone(),
                        edge: e.clone(),
                        face: hs.face.clone(),
                    });
                }
            }
        }
        DomainTool {
            arcs,
            iarc: 0,
            vertices: Vec::new(),
            ivtx: 0,
        }
    }

    /// OCCT Adaptor3d_TopolTool::Initialize(C) — the single-restriction
    /// domain.
    pub fn single(arc: &RstArc) -> Self {
        DomainTool {
            arcs: vec![arc.clone()],
            iarc: 0,
            vertices: Vec::new(),
            ivtx: 0,
        }
    }

    /// OCCT Init() — restart the arc iterator.
    pub fn init(&mut self) {
        self.iarc = 1;
    }

    /// OCCT More().
    pub fn more(&self) -> bool {
        self.iarc >= 1 && (self.iarc as usize) <= self.arcs.len()
    }

    /// OCCT Next().
    pub fn next(&mut self) {
        self.iarc += 1;
    }

    /// OCCT Value() — the current arc.
    pub fn value(&self) -> RstArc {
        self.arcs[(self.iarc - 1) as usize].clone()
    }

    /// OCCT Initialize(C) — the single-arc domain (also the vertex-scan
    /// source for InitVertexIterator).
    pub fn initialize_arc(&mut self, arc: &RstArc) {
        self.arcs = vec![arc.clone()];
        self.iarc = 1;
    }

    /// OCCT BRepTopAdaptor_TopolTool::InitVertexIterator (cxx L148-152) —
    /// the topological vertices of the current curve's edge.  A bare pcurve
    /// (no edge identity) has no vertices — the pending consumer-side
    /// boundary (see the file header).
    pub fn init_vertex_iterator(&mut self) {
        self.vertices.clear();
        self.ivtx = 0;
        let idx = if self.iarc >= 1 && (self.iarc as usize) <= self.arcs.len() {
            self.iarc as usize - 1
        } else if !self.arcs.is_empty() {
            0
        } else {
            return;
        };
        let arc = self.arcs[idx].clone();
        if arc.edge.is_null() {
            return;
        }
        let ed = arc.edge.as_edge().expect("arc edge");
        // OCCT TopExp::Vertices of the edge — first then last.
        for v in [&ed.first, &ed.last] {
            if v.is_null() {
                continue;
            }
            self.vertices.push(DomainVertex::from_topo(v));
        }
    }

    /// OCCT MoreVertex().
    pub fn more_vertex(&self) -> bool {
        self.ivtx >= 1 && (self.ivtx as usize) <= self.vertices.len()
    }

    /// OCCT NextVertex().
    pub fn next_vertex(&mut self) {
        self.ivtx += 1;
    }

    /// OCCT Vertex() — the current vertex (the BRepTopAdaptor_HVertex
    /// construction, cxx L167-170).
    pub fn vertex(&self) -> DomainVertex {
        self.vertices[(self.ivtx - 1) as usize].clone()
    }
}

/// OCCT BRepBlend_BlendTool::Parameter(V, C) (BlendTool.lxx L46-51) —
/// BRep_Tool::Parameter(myVtx, Edge, Face).  With a bare arc/vertex (no edge
/// identity) the ElCLib stand-in of Adaptor3d_HVertex::Parameter is the
/// pending boundary.
pub(crate) fn blend_tool_parameter(v: &DomainVertex, arc: &RstArc) -> f64 {
    if !v.vtx.is_null() && !arc.edge.is_null() {
        brep_tool_parameter(&arc.brep, &v.vtx, &arc.edge)
    } else if !v.vtx.is_null() {
        v.hv.parameter(&arc.curve as &dyn Curve2dAdaptor)
    } else {
        0.0
    }
}

/// OCCT BRepBlend_BlendTool::Tolerance(V, A) (BlendTool.lxx L39-44) —
/// V->Resolution(A).  The stored seed (BRep_Tool::Tolerance(V)) is the
/// pending stand-in for the BRepTopAdaptor_HVertex::Resolution chain.
pub(crate) fn blend_tool_tolerance(v: &DomainVertex, arc: &RstArc) -> f64 {
    v.hv.resolution(&arc.curve as &dyn Curve2dAdaptor)
}

/// OCCT BRepBlend_BlendTool::Project(P, Surf, C, Paramproj, Dist)
/// (BlendTool.cxx L38-79) — the orthogonal projection of a point on a
/// curve; the endpoints seed then the extrema refinement (the Surf argument
/// is unused in OCCT; the Extrema_EPCOfExtPC2d sampling parameters are
/// carried by the rcad ExtPC2d).
pub(crate) fn blend_tool_project(p: DVec2, c: &RstArc, paramproj: &mut f64, dist: &mut f64) -> bool {
    *paramproj = c.first_parameter();
    let mut p2d = c.value(*paramproj);
    *dist = p2d.distance(p);

    let t = c.last_parameter();
    p2d = c.value(t);
    if p2d.distance(p) < *dist {
        *paramproj = t;
        *dist = p2d.distance(p);
    }

    let (first, last) = (c.first_parameter(), c.last_parameter());
    let extrema = ExtPC2d::new(p, &c.curve, 1.0e-5, first, last);
    if !extrema.is_done() {
        return true;
    }

    let nbext = extrema.nb_ext();
    let mut adist2 = *dist * *dist;
    for i in 1..=nbext {
        if extrema.square_distance(i) < adist2 {
            adist2 = extrema.square_distance(i);
            *paramproj = extrema.point(i).param;
        }
    }
    *dist = adist2.sqrt();

    true
}

/// OCCT BRepBlend_BlendTool::Inters(P1, P2, Surf, C, Param, Dist)
/// (BlendTool.cxx L85-120) — the intersection of a segment with a curve
/// (Geom2d_Line segment + Geom2dInt_GInter; the Surf argument is unused in
/// OCCT).
pub(crate) fn blend_tool_inters(
    p1: DVec2,
    p2: DVec2,
    c: &RstArc,
    param: &mut f64,
    dist: &mut f64,
) -> bool {
    let tol = 1.0e-8;
    let v = p2 - p1;
    let mag = v.length();
    if mag < tol {
        return false;
    }

    let d = v / mag;
    let seg = Line2dAdaptor::new_pnt_dir(p1, d, -0.01 * mag, 1.01 * mag);
    let d1 = Res2dDomain::bounded(
        seg.value(seg.first_parameter()),
        seg.first_parameter(),
        tol,
        seg.value(seg.last_parameter()),
        seg.last_parameter(),
        tol,
    );
    let cfirst = c.first_parameter();
    let clast = c.last_parameter();
    let d2 = Res2dDomain::bounded(c.value(cfirst), cfirst, tol, c.value(clast), clast, tol);

    let mut inter = TheIntPCurvePCurveOfGInter::new();
    inter.perform(&seg, &d1, &c.curve as &dyn Curve2dAdaptor, &d2, tol, tol);
    if !inter.base.is_done() {
        return false;
    }

    let nbint = inter.base.nb_points();
    if nbint == 0 {
        return false;
    }

    let ip = inter.base.point(1);
    *param = ip.param_on_second();
    *dist = p1.distance(ip.value());
    true
}

