//! OCCT BRepFill_OffsetWire — GAP placeholder types (TKMath MAT2d /
//! Bisector / MAT stack, BRepFill_TrimEdgeTool, BRepTools_Substitution) and
//! the free helpers of BRepFill_OffsetWire.cxx — split from offset_wire.rs
//! (file-size rule).  The class itself lives in offset_wire.rs.

use glam::{DVec2, DVec3};
use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use rcad_kernel::geom::{Plane, Surface3};
use crate::brep_fill::generator::{shape_key, shape_oriented, shape_reversed, ShapeKey};
use super::offset_wire::{GeomAbsJoinType, MatSide};

/// OCCT Precision::Confusion().
pub(super) const TOL_CONFUSION: f64 = CONFUSION;
/// OCCT Precision::PConfusion().
pub(super) const TOL_PCONFUSION: f64 = PCONFUSION;

// =============================================================================
// GAP placeholder types (see the file header)
// =============================================================================

/// GAP: MAT_Node (TKMath/MAT) — not translated (plan §0.6).
#[derive(Debug, Clone)]
pub struct MatNode;

impl MatNode {
    /// OCCT MAT_Node::Infinite().
    pub fn infinite(&self) -> bool {
        panic!("GAP: MAT_Node (TKMath/MAT) is not translated — see file header")
    }
    /// OCCT MAT_Node::Distance().
    pub fn distance(&self) -> f64 {
        panic!("GAP: MAT_Node (TKMath/MAT) is not translated — see file header")
    }
}

/// GAP: MAT_Arc (TKMath/MAT) — not translated (plan §0.6).
#[derive(Debug, Clone)]
pub struct MatArc;

impl MatArc {
    /// OCCT MAT_Arc::FirstElement().
    pub fn first_element(&self) -> Shape {
        panic!("GAP: MAT_Arc (TKMath/MAT) is not translated — see file header")
    }
    /// OCCT MAT_Arc::SecondElement().
    pub fn second_element(&self) -> Shape {
        panic!("GAP: MAT_Arc (TKMath/MAT) is not translated — see file header")
    }
    /// OCCT MAT_Arc::FirstNode().
    pub fn first_node(&self) -> MatNode {
        panic!("GAP: MAT_Arc (TKMath/MAT) is not translated — see file header")
    }
    /// OCCT MAT_Arc::SecondNode().
    pub fn second_node(&self) -> MatNode {
        panic!("GAP: MAT_Arc (TKMath/MAT) is not translated — see file header")
    }
}

/// GAP: MAT_Graph (TKMath/MAT) — not translated (plan §0.6).
#[derive(Debug, Clone)]
pub struct MatGraph;

impl MatGraph {
    /// OCCT MAT_Graph::NumberOfArcs().
    pub fn number_of_arcs(&self) -> usize {
        panic!("GAP: MAT_Graph (TKMath/MAT) is not translated — see file header")
    }
    /// OCCT MAT_Graph::Arc(i).
    pub fn arc(&self, _i: usize) -> MatArc {
        panic!("GAP: MAT_Graph (TKMath/MAT) is not translated — see file header")
    }
}

/// GAP: Bisector_Bisec (TKMath/Bisector) — not translated (plan §0.6).
#[derive(Debug, Clone)]
pub struct BisectorBisec;

impl BisectorBisec {
    /// OCCT Bisector_Bisec::Value() — the underlying 2D curve handle.
    pub fn value(&self) -> Curve2d {
        panic!("GAP: Bisector_Bisec (TKMath/Bisector) is not translated — see file header")
    }
}

/// GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) — not translated (plan §0.6).
#[derive(Debug, Clone, Default)]
pub struct BRepMAT2dBisectingLocus;

impl BRepMAT2dBisectingLocus {
    /// OCCT BRepMAT2d_BisectingLocus::NumberOfContours().
    pub fn number_of_contours(&self) -> usize {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::NumberOfElts(ic).
    pub fn number_of_elts(&self, _ic: usize) -> usize {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::BasicElt(ic, ie).
    pub fn basic_elt(&self, _ic: usize, _ie: usize) -> Shape {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::Graph().
    pub fn graph(&self) -> MatGraph {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::GeomBis(Arc, Reverse).
    pub fn geom_bis(&self, _a: &MatArc, _reverse: &mut bool) -> BisectorBisec {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::GeomElt(Node / BasicElt) — the 2D point.
    pub fn geom_elt_node(&self, _n: &MatNode) -> DVec2 {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    pub fn geom_elt(&self, _s: &Shape) -> Curve2d {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::Compute(...).
    #[allow(clippy::too_many_arguments)]
    pub fn compute(
        &mut self,
        _exp: &BRepMAT2dExplorer,
        _is: usize,
        _side: MatSide,
        _join: GeomAbsJoinType,
        _open: bool,
    ) {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_BisectingLocus::IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepMAT2d_BisectingLocus (TKMath/MAT2d) is not translated — see file header")
    }
}

/// GAP: BRepMAT2d_LinkTopoBilo (TKMath/MAT2d) — not translated (plan §0.6).
#[derive(Debug, Clone, Default)]
pub struct BRepMAT2dLinkTopoBilo;

impl BRepMAT2dLinkTopoBilo {
    /// OCCT BRepMAT2d_LinkTopoBilo::Perform(Exp, Locus).
    pub fn perform(&mut self, _exp: &BRepMAT2dExplorer, _locus: &BRepMAT2dBisectingLocus) {
        panic!("GAP: BRepMAT2d_LinkTopoBilo (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_LinkTopoBilo::GeneratingShape(BasicElt).
    pub fn generating_shape(&self, _e: &Shape) -> Shape {
        panic!("GAP: BRepMAT2d_LinkTopoBilo (TKMath/MAT2d) is not translated — see file header")
    }
}

/// GAP: BRepMAT2d_Explorer (TKMath/MAT2d) — not translated (plan §0.6).
#[derive(Debug, Clone, Default)]
pub struct BRepMAT2dExplorer;

impl BRepMAT2dExplorer {
    /// OCCT BRepMAT2d_Explorer::Perform(Face).
    pub fn perform(&mut self, _face: &Shape) {
        panic!("GAP: BRepMAT2d_Explorer (TKMath/MAT2d) is not translated — see file header")
    }
    /// OCCT BRepMAT2d_Explorer::ModifiedShape(S).
    pub fn modified_shape(&self, _s: &Shape) -> Shape {
        panic!("GAP: BRepMAT2d_Explorer (TKMath/MAT2d) is not translated — see file header")
    }
}

/// GAP: BRepFill_TrimEdgeTool (TKBool/BRepFill/BRepFill_TrimEdgeTool.cxx) —
/// its own file, not part of this stage (plan §0.6).
#[derive(Debug, Clone)]
pub struct BRepFillTrimEdgeTool;

impl BRepFillTrimEdgeTool {
    /// OCCT BRepFill_TrimEdgeTool::BRepFill_TrimEdgeTool(Bisec, E1, E2, Offset).
    pub fn new(_bisec: &BisectorBisec, _e1: &Curve2d, _e2: &Curve2d, _offset: f64) -> Self {
        panic!("GAP: BRepFill_TrimEdgeTool is not translated — see file header")
    }
    /// OCCT BRepFill_TrimEdgeTool::IntersectWith(E1, E2, S1, S2, VS, VE, Join,
    /// IsOpenResult, Params).
    #[allow(clippy::too_many_arguments)]
    pub fn intersect_with(
        &self,
        _e1: &Shape,
        _e2: &Shape,
        _s1: &Shape,
        _s2: &Shape,
        _vs: &Shape,
        _ve: &Shape,
        _join: GeomAbsJoinType,
        _open: bool,
        _params: &mut Vec<DVec3>,
    ) {
        panic!("GAP: BRepFill_TrimEdgeTool is not translated — see file header")
    }
    /// OCCT BRepFill_TrimEdgeTool::AddOrConfuse(Start, E1, E2, Params).
    pub fn add_or_confuse(&self, _start: bool, _e1: &Shape, _e2: &Shape, _params: &mut Vec<DVec3>) {
        panic!("GAP: BRepFill_TrimEdgeTool is not translated — see file header")
    }
    /// OCCT BRepFill_TrimEdgeTool::IsInside(P).
    pub fn is_inside(&self, _p: DVec2) -> bool {
        panic!("GAP: BRepFill_TrimEdgeTool is not translated — see file header")
    }
}

/// OCCT BRepTools_Substitution (TKTopAlgo/BRepTools) — GAP placeholder
/// (plan §0.6): each shape maps to its replacement list.
#[derive(Debug, Default)]
pub struct BRepToolsSubstitution;

impl BRepToolsSubstitution {
    /// OCCT BRepTools_Substitution::Substitute(Old, List).
    pub fn substitute(&mut self, _old: &Shape, _list: &[Shape]) {
        panic!("GAP: BRepTools_Substitution (TKTopAlgo/BRepTools) is not translated — see file header")
    }
    /// OCCT BRepTools_Substitution::Build(S).
    pub fn build(&mut self, _brep: &BRep, _s: &Shape) {
        panic!("GAP: BRepTools_Substitution (TKTopAlgo/BRepTools) is not translated — see file header")
    }
    /// OCCT BRepTools_Substitution::IsCopied(S).
    pub fn is_copied(&self, _s: &Shape) -> bool {
        panic!("GAP: BRepTools_Substitution (TKTopAlgo/BRepTools) is not translated — see file header")
    }
    /// OCCT BRepTools_Substitution::Copy(S).
    pub fn copy(&self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepTools_Substitution (TKTopAlgo/BRepTools) is not translated — see file header")
    }
    /// OCCT BRepTools_Substitution::Clear().
    pub fn clear(&mut self) {}
}

// =============================================================================
// Kernel-mapping helpers
// =============================================================================

pub(super) fn wire_edges(brep: &BRep, w: &Shape) -> Vec<Shape> {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

pub(super) fn face_wires(brep: &BRep, f: &Shape) -> Vec<Shape> {
    match f.data.as_ref() {
        TShape::Face(fd) => {
            let mut ws = vec![fd.outer_wire.clone()];
            ws.extend(fd.inner_wires.iter().cloned());
            ws
        }
        _ => Vec::new(),
    }
}

/// OCCT static EdgeVertices (L127-137): the vertices in traversal order.
pub(super) fn edge_vertices(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    if e.orientation == Orientation::Reversed {
        let ed = brep.edge(e.clone());
        (ed.last.clone(), ed.first.clone())
    } else {
        let ed = brep.edge(e.clone());
        (ed.first.clone(), ed.last.clone())
    }
}

/// OCCT BRep_Tool::Surface(F) — the face surface (rcad carries it in the
/// TFaceData; the TopLoc_Location is identity for these faces).
pub(super) fn brep_tool_surface(brep: &BRep, f: &Shape) -> Option<Surface3> {
    match f.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, f, l) — the pcurve bound to the face.
pub(super) fn brep_tool_curve_on_surface(brep: &BRep, e: &Shape, f: &Shape) -> Option<(Curve2d, f64, f64)> {
    let ed = brep.edge(e.clone());
    let fkey = (f.ptr_id(), f.location);
    // The seam representation carries two pcurves; CurveOnSurface returns the
    // first one.
    for r in &ed.representations {
        if let rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
            face,
            pcurve1,
            range,
            ..
        } = r
        {
            if *face == fkey {
                return Some((pcurve1.clone(), range[0], range[1]));
            }
        }
    }
    ed.pcurves.get(&fkey).map(|(c, a, b)| (c.clone(), *a, *b))
}

/// OCCT BRep_Tool::Tolerance (vertex / edge / face).
pub(super) fn brep_tool_tolerance(brep: &BRep, s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Pnt.
pub(super) fn vertex_point(brep: &BRep, v: &Shape) -> DVec3 {
    brep.vertex(v.clone()).point
}

/// OCCT BRep_Builder::UpdateVertex(V, P, Tol).
pub(super) fn update_vertex_point(brep: &mut BRep, v: &Shape, p: DVec3, tol: f64) {
    let vd = brep.vertex_mut(v.clone());
    vd.point = p;
    if tol > vd.tolerance {
        vd.tolerance = tol;
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, Tol).
pub(super) fn update_vertex_tolerance(brep: &mut BRep, v: &Shape, tol: f64) {
    let vd = brep.vertex_mut(v.clone());
    if tol > vd.tolerance {
        vd.tolerance = tol;
    }
}

/// OCCT BRepLib_MakeVertex(P) — always a NEW vertex TShape.
pub(super) fn make_vertex(brep: &mut BRep, p: DVec3) -> Shape {
    brep.add_tvertex_unique(p)
}

/// OCCT BRepLib_MakeEdge(V1, V2) — a straight edge between two vertices.
pub(super) fn make_edge_vertices(brep: &mut BRep, v1: &Shape, v2: &Shape) -> Shape {
    let p1 = vertex_point(brep, v1);
    let p2 = vertex_point(brep, v2);
    let dir = p2 - p1;
    let len = dir.length();
    let curve = if len < 1e-300 {
        None
    } else {
        Some(Curve3::Line(rcad_kernel::geom::Line3::new(p1, dir / len)))
    };
    brep.add_tedge(
        curve,
        shape_oriented(v1, Orientation::Forward),
        shape_oriented(v2, Orientation::Reversed),
        [0.0, len],
    )
}

/// OCCT BRep_Builder::Remove(A, B) for a wire child (rebuild the list).
pub(super) fn builder_remove_from_wire(brep: &mut BRep, w: &Shape, child: &Shape) {
    let idx = w.index;
    if let TShape::Wire(wd) = Arc::make_mut(&mut brep.tshapes[idx]) {
        wd.edges.retain(|e| !e.is_same(child));
        wd.my_shapes.retain(|e| !e.is_same(child));
    }
}

/// OCCT TopExp_Explorer(S, VERTEX/EDGE/WIRE) for the compound/wire levels
/// used here — flattened to the owned children.
pub(super) fn explored_children(brep: &BRep, s: &Shape, t: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    match s.data.as_ref() {
        TShape::Compound(children) | TShape::CompSolid(children) => {
            for c in children {
                if tshape_type(c.data.as_ref()) == t {
                    out.push(c.clone());
                }
                out.extend(explored_children(brep, c, t));
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                if t == ShapeType::Shell {
                    out.push(sh.clone());
                }
                out.extend(explored_children(brep, sh, t));
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                if t == ShapeType::Face {
                    out.push(f.clone());
                }
                out.extend(explored_children(brep, f, t));
            }
        }
        TShape::Face(_) => {
            if t == ShapeType::Wire {
                out.extend(face_wires(brep, s));
            }
        }
        TShape::Wire(_) => {
            if t == ShapeType::Edge {
                out.extend(wire_edges(brep, s));
            }
        }
        TShape::Edge(_) => {
            if t == ShapeType::Vertex {
                let ed = brep.edge(s.clone());
                out.push(ed.first.clone());
                out.push(ed.last.clone());
            }
        }
        _ => {}
    }
    out
}

/// The shape type of a TShape (TopoDS_Shape::ShapeType()).
pub(super) fn tshape_type(t: &TShape) -> ShapeType {
    match t {
        TShape::Vertex(_) => ShapeType::Vertex,
        TShape::Edge(_) => ShapeType::Edge,
        TShape::Wire(_) => ShapeType::Wire,
        TShape::Face(_) => ShapeType::Face,
        TShape::Shell(_) => ShapeType::Shell,
        TShape::Solid(_) => ShapeType::Solid,
        TShape::CompSolid(_) => ShapeType::CompSolid,
        TShape::Compound(_) => ShapeType::Compound,
    }
}

/// The wire Closed flag (TopoDS_Shape::Closed()).
pub(super) fn is_closed_wire(brep: &BRep, w: &Shape) -> bool {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.flags & tshape_flags::CLOSED != 0,
        _ => false,
    }
}

pub(super) fn set_wire_closed(brep: &mut BRep, w: &Shape, closed: bool) {
    let idx = w.index;
    if let TShape::Wire(wd) = Arc::make_mut(&mut brep.tshapes[idx]) {
        if closed {
            wd.flags |= tshape_flags::CLOSED;
        } else {
            wd.flags &= !tshape_flags::CLOSED;
        }
    }
}

// =============================================================================
// Static helpers (BRepFill_OffsetWire.cxx)
// =============================================================================

/// OCCT static PerformCurve (L101-108 decl, L2781-2823): discretize the
/// curve within the deflection (QuasiFleche).
pub(super) fn perform_curve(
    parameters: &mut Vec<f64>,
    points: &mut Vec<DVec3>,
    c: &Curve3,
    deflection: f64,
    u1: f64,
    u2: f64,
    epsilon: f64,
    nbmin: usize,
) -> bool {
    let uu1 = u1.min(u2);
    let uu2 = u1.max(u2);

    let (pdeb, ddeb) = (c.point_at(uu1), c.derivative_at(uu1));
    parameters.push(uu1);
    points.push(pdeb);

    let (pfin, dfin) = (c.point_at(uu2), c.derivative_at(uu2));

    let a_delta = uu2 - uu1;
    let a_dist = pdeb.distance(pfin);

    if (a_delta / a_dist) > 5.0e-14 {
        quasi_fleche(
            c,
            deflection * deflection,
            uu1,
            pdeb,
            ddeb,
            uu2,
            pfin,
            dfin,
            nbmin,
            epsilon * epsilon,
            parameters,
            points,
        );
    }

    true
}

/// OCCT static QuasiFleche (L88-99 decl, L2827-2922): recursive bounded-
/// deflection subdivision (IntWalk_IWalking::TestDeflection arrow estimate).
#[allow(clippy::too_many_arguments)]
pub(super) fn quasi_fleche(
    c: &Curve3,
    deflection2: f64,
    udeb: f64,
    pdeb: DVec3,
    vdeb: DVec3,
    ufin: f64,
    pfin: DVec3,
    vfin: DVec3,
    nbmin: usize,
    eps: f64,
    parameters: &mut Vec<f64>,
    points: &mut Vec<DVec3>,
) {
    let ptslength = points.len();
    let mut udelta = ufin - udeb;
    let (pdelta, vdelta);
    if nbmin > 2 {
        udelta /= (nbmin - 1) as f64;
        pdelta = c.point_at(udeb + udelta);
        vdelta = c.derivative_at(udeb + udelta);
    } else {
        pdelta = pfin;
        vdelta = vfin;
    }

    let norme = (pdelta - pdeb).length_squared();
    let mut the_fleche = 0.0f64;
    let mut flecheok = false;
    if norme > eps {
        // Evaluation of the arrow by interpolation. See IntWalk_IWalking::TestDeflection
        let n1 = vdeb.length_squared();
        let n2 = vdelta.length_squared();
        if n1 > eps && n2 > eps {
            let normediff = (vdeb / vdeb.length() - vdelta / vdelta.length()).length_squared();
            if normediff > eps {
                the_fleche = normediff * norme / 64.0;
                flecheok = true;
            }
        }
    }
    if !flecheok {
        let pmid = (pdeb + pdelta) / 2.0;
        let pverif = c.point_at(udeb + udelta / 2.0);
        the_fleche = pmid.distance_squared(pverif);
    }

    if the_fleche < deflection2 {
        parameters.push(udeb + udelta);
        points.push(pdelta);
    } else {
        quasi_fleche(
            c,
            deflection2,
            udeb,
            pdeb,
            vdeb,
            udeb + udelta,
            pdelta,
            vdelta,
            3,
            eps,
            parameters,
            points,
        );
    }

    if nbmin > 2 {
        quasi_fleche(
            c,
            deflection2,
            udeb + udelta,
            pdelta,
            vdelta,
            ufin,
            pfin,
            vfin,
            nbmin - (points.len() - ptslength),
            eps,
            parameters,
            points,
        );
    }
}

/// OCCT static CheckBadEdges (L110-114 decl, L2658-2777): collect the spine
/// edges whose curvature can exceed 1/Offset.  GAP: GeomLProp_CLProps2d
/// (TKMath/GeomLProp) is not translated (plan §0.6).
pub(super) fn check_bad_edges(
    _brep: &BRep,
    _spine: &Shape,
    _offset: f64,
    _locus: &BRepMAT2dBisectingLocus,
    _link: &BRepMAT2dLinkTopoBilo,
    _bad_edges: &mut Vec<Shape>,
) {
    panic!("GAP: CheckBadEdges requires GeomLProp_CLProps2d (TKMath, not translated) — see file header")
}

/// OCCT static CutEdge (L116-119 decl, L1930-2071): cut the edge at the
/// extrema of curvature and the inflexion points.  GAP: MAT2d_CutCurve
/// (TKMath/MAT2d) is not translated (plan §0.6).
pub(super) fn cut_edge(
    _brep: &mut BRep,
    _e: &Shape,
    _f: &Shape,
    _force_cut: i32,
    _cuts: &mut Vec<Shape>,
) -> usize {
    panic!("GAP: CutEdge requires MAT2d_CutCurve (TKMath/MAT2d, not translated) — see file header")
}

/// OCCT static CutCurve (L121-123 decl, L2076-2124): split a trimmed 2D
/// curve into nbParts parts.
pub(super) fn cut_curve(
    _c: &Curve2d,
    _nb_parts: usize,
    _the_curves: &mut Vec<Curve2d>,
) {
    panic!("GAP: CutCurve requires Geom2d_TrimmedCurve evaluation (see file header)")
}

/// OCCT static VertexFromNode (L139-143 decl, L2310-2344).
pub(super) fn vertex_from_node(
    brep: &mut BRep,
    a_node: &MatNode,
    offset: f64,
    pn: DVec2,
    map_node_vertex: &mut HashMap<ShapeKey, Shape>,
    vn: &mut Shape,
) -> bool {
    let status;
    let tol = TOL_CONFUSION;

    if !a_node.infinite() && (a_node.distance() - offset).abs() < tol {
        //------------------------------------------------
        // the Node gives a vertex on the offset
        //------------------------------------------------
        // OCCT L2325: MapNodeVertex.IsBound(aNode) — the map is keyed by the
        // MAT_Node handle.  rcad has no node identity yet (see the MatNode
        // gap note), so the bound test is approximated by the key set.
        if !map_node_vertex.is_empty() {
            // MapNodeVertex(aNode)
            *vn = map_node_vertex.values().next().cloned().unwrap();
        } else {
            // B.MakeVertex(VN); B.UpdateVertex(VN, P, Precision::Confusion());
            // MapNodeVertex.Bind(aNode, VN);
            let p = DVec3::new(pn.x, pn.y, 0.0);
            *vn = make_vertex(brep, p);
            update_vertex_point(brep, vn, p, TOL_CONFUSION);
            map_node_vertex.insert(ShapeKey(vn.ptr_id()), vn.clone());
        }
        status = true;
    } else {
        status = false;
    }

    status
}

/// OCCT static StoreInMap (L145-148 decl, L2348-2375).
pub(super) fn store_in_map(
    _brep: &BRep,
    v1: &Shape,
    v2: &Shape,
    map_vv: &mut Vec<(Shape, Shape)>,
) {
    let mut old_v = v1.clone();
    let mut new_v = v2.clone();

    // if (MapVV.Contains(V2)) NewV = MapVV.FindFromKey(V2);
    if let Some((_, val)) = map_vv.iter().find(|(k, _)| k.is_same(v2)) {
        new_v = val.clone();
    }

    // if (MapVV.Contains(V1)) MapVV.ChangeFromKey(V1) = NewV;
    for (k, val) in map_vv.iter_mut() {
        if k.is_same(v1) {
            *val = new_v.clone();
        }
    }

    for (_, val) in map_vv.iter_mut() {
        if val.is_same(&old_v) {
            *val = new_v.clone();
        }
    }
    let _ = &mut old_v;

    // MapVV.Add(V1, NewV)
    if !map_vv.iter().any(|(k, _)| k.is_same(v1)) {
        map_vv.push((v1.clone(), new_v));
    }
}

/// OCCT static IsInnerEdge (L161-164 decl, L2576-2602).
pub(super) fn is_inner_edge(brep: &BRep, pro_e: &Shape, all_spine: &Shape, tr_par1: &mut f64, tr_par2: &mut f64) -> bool {
    if tshape_type(pro_e.data.as_ref()) != ShapeType::Edge {
        return false;
    }

    let an_edge = pro_e.clone();

    // TopExp::MapShapesAndAncestors(AllSpine, VERTEX, EDGE, VEmap)
    let mut vemap: Vec<(Shape, Vec<Shape>)> = Vec::new();
    for w in face_wires(brep, all_spine) {
        for e in wire_edges(brep, &w) {
            let ed = brep.edge(e.clone());
            for v in [&ed.first, &ed.last] {
                match vemap.iter_mut().find(|(kv, _)| kv.is_same(v)) {
                    Some((_, l)) => l.push(e.clone()),
                    None => vemap.push((v.clone(), vec![e.clone()])),
                }
            }
        }
    }
    for (_, le) in &vemap {
        if le.len() == 1 && an_edge.is_same(&le[0]) {
            return false;
        }
    }

    let r = brep.edge(an_edge).range;
    *tr_par1 = r[0];
    *tr_par2 = r[1];
    true
}

/// OCCT static DoubleOrNotInside (L166 decl, L2609-2629): true when V appears
/// twice in LV or is not inside.
pub(super) fn double_or_not_inside(lv: &[Shape], v: &Shape) -> bool {
    let mut vu = false;
    for s in lv {
        if v.is_same(s) {
            if vu {
                return true;
            } else {
                vu = true;
            }
        }
    }
    !vu
}

/// OCCT static IsSmallClosedEdge (L168 decl, L2631-2656).
pub(super) fn is_small_closed_edge(brep: &BRep, an_edge: &Shape, a_vertex: &Shape) -> bool {
    let pv = vertex_point(brep, a_vertex);
    let pv2d = DVec2::new(pv.x, pv.y);

    // The pcurve of the edge (its first curve representation).
    let pcurve = {
        let ed = brep.edge(an_edge.clone());
        // The rcad pcurve map is keyed by face; the first bound pcurve is the
        // (unique) curve representation of these offset edges.
        match ed.pcurves.iter().next() {
            Some((_, (c, a, b))) => Some((c.clone(), *a, *b)),
            None => match ed.representations.first() {
                Some(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    pcurve1,
                    range,
                    ..
                }) => Some((pcurve1.clone(), range[0], range[1])),
                _ => None,
            },
        }
    };
    let Some((pcurve, fpar, lpar)) = pcurve else {
        return false;
    };
    let pfirst = pcurve.point_at(fpar);
    let plast = pcurve.point_at(lpar);
    let pmid = pcurve.point_at((fpar + lpar) * 0.5);

    let mut the_tol = brep_tool_tolerance(brep, a_vertex);
    the_tol *= 1.5;

    let dist1 = pfirst.distance(pv2d);
    let dist2 = plast.distance(pv2d);
    let dist3 = pmid.distance(pv2d);

    dist1 <= the_tol && dist2 <= the_tol && dist3 <= the_tol
}

/// OCCT CheckSmallParamOnEdge (L188 decl, L2924-2940).
pub fn check_small_param_on_edge(brep: &BRep, an_edge: &Shape) -> bool {
    let ed = brep.edge(an_edge.clone());
    // aList.IsEmpty() -> true
    let has_rep = !ed.pcurves.is_empty() || !ed.representations.is_empty();
    if has_rep {
        let (f, l) = if !ed.pcurves.is_empty() {
            let (_, (_, a, b)) = ed.pcurves.iter().next().unwrap();
            (*a, *b)
        } else {
            (ed.range[0], ed.range[1])
        };
        if (l - f).abs() < TOL_PCONFUSION {
            return false;
        }
    }
    true
}

/// OCCT static TrimEdge (L150-159 decl, L2379-2572): build the trimmed
/// offset edges from the vertex/parameter sequences.
#[allow(clippy::too_many_arguments)]
pub(super) fn trim_edge_offset(
    brep: &mut BRep,
    e: &Shape,
    pro_e: &Shape,
    all_spine: &Shape,
    detromp: &[Shape],
    the_ver: &mut Vec<Shape>,
    the_par: &mut Vec<f64>,
    s: &mut Vec<Shape>,
    map_vv: &mut Vec<(Shape, Shape)>,
    ind_of_e: i32,
) {
    let change = true;
    let _ = change;
    s.clear();

    //-----------------------------------------------------------
    // Parse two sequences depending on the parameter on the edge.
    //-----------------------------------------------------------
    // while (Change) bubble sort (L2397-2409)
    loop {
        let mut changed = false;
        for i in 0..the_par.len().saturating_sub(1) {
            if the_par[i] > the_par[i + 1] {
                the_par.swap(i, i + 1);
                the_ver.swap(i, i + 1);
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    //----------------------------------------------------------
    // If a vertex is not in the proofing, it is eliminated.
    //----------------------------------------------------------
    if !brep.edge(e.clone()).degenerated {
        let mut k = 0usize;
        while k < the_ver.len() {
            if double_or_not_inside(detromp, &the_ver[k]) {
                the_ver.remove(k);
                the_par.remove(k);
            } else {
                k += 1;
            }
        }
    }

    //----------------------------------------------------------
    // If a vertex_double appears twice in the proofing
    // the vertex is removed.
    // otherwise preserve only one of its representations.
    //----------------------------------------------------------
    if !brep.edge(e.clone()).degenerated {
        let a_par_tol = 2.0 * TOL_PCONFUSION;
        let mut k = 0usize;
        while k + 1 < the_ver.len() {
            if the_ver[k].is_same(&the_ver[k + 1]) || (the_par[k] - the_par[k + 1]).abs() <= a_par_tol {
                if k + 2 == the_ver.len() {
                    store_in_map(brep, &the_ver[k], &the_ver[k + 1], map_vv);
                    the_ver.remove(k);
                    the_par.remove(k);
                } else {
                    store_in_map(brep, &the_ver[k + 1], &the_ver[k], map_vv);
                    the_ver.remove(k + 1);
                    the_par.remove(k + 1);
                }
                // k--; — the OCCT loop re-scans the same index after k--
                continue;
            }
            k += 1;
        }
    }
    //-----------------------------------------------------------
    // Creation of edges.
    // the number of vertices should be even. The created edges
    // go from a vertex with uneven index i to vertex i+1;
    //-----------------------------------------------------------
    if ind_of_e == 1 || ind_of_e == -1 {
        // open result and extreme edges of result
        let mut new_edge = brep.empty_copied(e);
        let (v1, v2) = edge_vertices(brep, e);
        let (mut v1, mut v2) = (v1, v2);
        let pcv = brep_tool_curve_on_surface(brep, e, all_spine);
        let (mut fpar, mut lpar) = match &pcv {
            Some((_, a, b)) => (*a, *b),
            None => {
                let r = brep.edge(e.clone()).range;
                (r[0], r[1])
            }
        };
        let pcurve = pcv.map(|(c, _, _)| c);
        let mut tr_par1 = 0.0f64;
        let mut tr_par2 = 0f64;
        let to_trim_as_origin = is_inner_edge(brep, pro_e, all_spine, &mut tr_par1, &mut tr_par2);

        if ind_of_e == 1 {
            // first edge of open wire
            if new_edge.orientation == Orientation::Forward {
                if to_trim_as_origin {
                    fpar = tr_par1;
                    let tr_pnt2d = pcurve.as_ref().map(|c| c.point_at(fpar));
                    if let Some(tp) = tr_pnt2d {
                        v1 = make_vertex(brep, DVec3::new(tp.x, tp.y, 0.0));
                    }
                }
                // TheBuilder.Add(NewEdge, V1 FORWARD); Add(NewEdge, TheVer.First() REVERSED);
                // Range(NewEdge, fpar, ThePar.First())
                new_edge = rebuild_edge_with_vertices(
                    brep,
                    &new_edge,
                    shape_oriented(&v1, Orientation::Forward),
                    shape_oriented(&the_ver[0], Orientation::Reversed),
                    fpar,
                    the_par[0],
                );
            } else {
                if to_trim_as_origin {
                    lpar = tr_par2;
                    let tr_pnt2d = pcurve.as_ref().map(|c| c.point_at(lpar));
                    if let Some(tp) = tr_pnt2d {
                        v2 = make_vertex(brep, DVec3::new(tp.x, tp.y, 0.0));
                    }
                }
                new_edge = rebuild_edge_with_vertices(
                    brep,
                    &new_edge,
                    shape_oriented(&the_ver[0], Orientation::Reversed),
                    shape_oriented(&v2, Orientation::Forward),
                    the_par[0],
                    lpar,
                );
            }
        } else {
            // last edge of open wire
            if new_edge.orientation == Orientation::Forward {
                if to_trim_as_origin {
                    lpar = tr_par2;
                    let tr_pnt2d = pcurve.as_ref().map(|c| c.point_at(lpar));
                    if let Some(tp) = tr_pnt2d {
                        v2 = make_vertex(brep, DVec3::new(tp.x, tp.y, 0.0));
                    }
                }
                new_edge = rebuild_edge_with_vertices(
                    brep,
                    &new_edge,
                    shape_oriented(&the_ver[0], Orientation::Forward),
                    shape_oriented(&v2, Orientation::Reversed),
                    the_par[0],
                    lpar,
                );
            } else {
                if to_trim_as_origin {
                    fpar = tr_par1;
                    let tr_pnt2d = pcurve.as_ref().map(|c| c.point_at(fpar));
                    if let Some(tp) = tr_pnt2d {
                        v1 = make_vertex(brep, DVec3::new(tp.x, tp.y, 0.0));
                    }
                }
                new_edge = rebuild_edge_with_vertices(
                    brep,
                    &new_edge,
                    shape_oriented(&v1, Orientation::Reversed),
                    shape_oriented(&the_ver[0], Orientation::Forward),
                    fpar,
                    the_par[0],
                );
            }
        }
        s.push(new_edge);
    } else {
        let mut k = 0usize;
        while k + 1 < the_ver.len() {
            let new_edge0 = brep.empty_copied(e);
            let new_edge = if new_edge0.orientation == Orientation::Reversed {
                rebuild_edge_with_vertices(
                    brep,
                    &new_edge0,
                    shape_oriented(&the_ver[k], Orientation::Reversed),
                    shape_oriented(&the_ver[k + 1], Orientation::Forward),
                    the_par[k],
                    the_par[k + 1],
                )
            } else {
                rebuild_edge_with_vertices(
                    brep,
                    &new_edge0,
                    shape_oriented(&the_ver[k], Orientation::Forward),
                    shape_oriented(&the_ver[k + 1], Orientation::Reversed),
                    the_par[k],
                    the_par[k + 1],
                )
            };
            s.push(new_edge);
            k += 2;
        }
    }
}

/// OCCT TheBuilder.Add(NewEdge, V1); TheBuilder.Add(NewEdge, V2);
/// TheBuilder.Range(NewEdge, f, l) — rcad stores the endpoints and range in
/// the TEdgeData (B.Add/Range mapping).
pub(super) fn rebuild_edge_with_vertices(
    brep: &mut BRep,
    new_edge: &Shape,
    first: Shape,
    last: Shape,
    f: f64,
    l: f64,
) -> Shape {
    let ed = brep.edge_mut_inplace(new_edge.clone());
    ed.first = first;
    ed.last = last;
    ed.my_shapes = vec![ed.first.clone(), ed.last.clone()];
    ed.range = [f.min(l), f.max(l)];
    new_edge.clone()
}

/// OCCT static MakeCircle (L170-176 decl, L2130-2165): the offset circle at
/// a vertex.  GAP: BRepLib_MakeEdge(pcurve, plane) — pcurve-only edge
/// creation — is not available in the rcad kernel (plan §0.6).
pub(super) fn make_circle(
    _brep: &mut BRep,
    _e: &Shape,
    _v: &Shape,
    _f: &Shape,
    _offset: f64,
    _map: &mut Vec<(Shape, Vec<Shape>)>,
    _ref_plane: &Plane,
) {
    panic!("GAP: MakeCircle requires BRepLib_MakeEdge(Geom2d, Geom_Plane) pcurve-only edge creation — see file header")
}

/// OCCT static MakeOffset (L178-186 decl, L2169-2306): the parallel edge of
/// a spine edge.  GAP: Adaptor2d_OffsetCurve / Geom2d_OffsetCurve and the
/// pcurve-only BRepLib_MakeEdge are not translated (plan §0.6).
#[allow(clippy::too_many_arguments)]
pub(super) fn make_offset(
    _brep: &mut BRep,
    _e: &Shape,
    _f: &Shape,
    _offset: f64,
    _map: &mut Vec<(Shape, Vec<Shape>)>,
    _ref_plane: &Plane,
    _is_open_result: bool,
    _the_join_type: GeomAbsJoinType,
    _ends: &[Shape; 2],
) {
    panic!("GAP: MakeOffset requires Adaptor2d_OffsetCurve / Geom2d_OffsetCurve + BRepLib_MakeEdge(Geom2d, Geom_Plane) — see file header")
}

// OCCT static KPartCircle (L192-199 decl, L201-306).
#[allow(clippy::too_many_arguments)]
pub(super) fn k_part_circle(
    brep: &mut BRep,
    my_spine: &Shape,
    my_offset: f64,
    is_open_result: bool,
    alt: f64,
    my_shape: &mut Shape,
    my_map: &mut Vec<(Shape, Vec<Shape>)>,
    my_is_done: &mut bool,
) -> bool {
    // TopExp_Explorer(mySpine, TopAbs_EDGE): collect the edges; more than
    // one -> not a KPart.
    let mut e: Option<Shape> = None;
    for w in face_wires(brep, my_spine) {
        for ecur in wire_edges(brep, &w) {
            if e.is_some() {
                return false;
            }
            e = Some(ecur);
        }
    }
    let Some(e) = e else {
        return false;
    };

    let ed = brep.edge(e.clone());
    let (f, l) = (ed.range[0], ed.range[1]);
    let mut c = match &ed.curve {
        Some(c) => c.clone(),
        None => return false,
    };

    // OCCT L223-227: if (C->IsKind(Geom_TrimmedCurve)) C = Ct->BasisCurve();
    if let Curve3::Trimmed(t) = &c {
        c = (*t.curve).clone();
    }

    let is_circle = matches!(c, Curve3::Circle(_));
    let edge_closed = ed.flags & tshape_flags::CLOSED != 0;

    if (is_circle && edge_closed) || is_open_result {
        let mut an_offset = my_offset;

        // OCCT L234-275: the offset pcurve of the spine edge.  GAP: the
        // Geom2dAdaptor_Curve / Adaptor2d_OffsetCurve / Geom2d_OffsetCurve
        // family is not translated (plan §0.6).
        let _ = (&mut an_offset, f, l);
        panic!("GAP: KPartCircle requires Geom2dAdaptor_Curve / Adaptor2d_OffsetCurve pcurve offsets + BRepLib_MakeEdge(Geom2d, Geom_Plane) — see file header")
    }

    let _ = (alt, my_shape, my_map, my_is_done);
    false
}

// =============================================================================
