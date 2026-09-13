//! 1:1 translation of OCCT `ShapeFix_SplitTool`
//! (`TKShHealing/ShapeFix/ShapeFix_SplitTool.hxx` L17-91 +
//! `ShapeFix_SplitTool.cxx` L1-580, docket row `ShapeFix_SplitTool`).
//!
//! Function-count equation (OCCT SplitTool.cxx = rcad split_tool.rs):
//! `ShapeFix_SplitTool()` (L42) / `SplitEdge(edge,param,vert,face,...)`
//! (L46-145) / `SplitEdge(edge,param1,param2,vert,face,...)` (L149-202) /
//! `CutEdge` (L206-302) / `SplitEdge(edge,fp,V1,lp,V2,face,SeqE,...)`
//! (L306-579) — 5 OCCT functions = 5 rcad functions.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — OCCT `BRep_Tool`/`BRep_Builder` read and write
//!    TShapes through global accessors; the rcad equivalents take
//!    `brep: &mut BRep` (the ShapeFix_Edge precedent).
//! 2. OCCT overload disambiguation suffixes: `SplitEdge(param1,param2)` ->
//!    `split_edge_two_params`, `SplitEdge(fp,V1,lp,V2,SeqE)` ->
//!    `split_edge_sequence` (the ShapeFix_Wire precedent).
//! 3. `TopAbs_Orientation` -> `Orientation`; `TopoDS_Shape::Oriented(o)` ->
//!    the copy helper `shape_oriented` (the wire/mod.rs re-host).
//! 4. `BRep_Builder::Range(E,F,L,force)` — the rcad `set_edge_range` carries
//!    the (first,last) write; the OCCT `force` flag selects the tolerance
//!    update arm inside BRep_Builder, which the rcad kernel builder performs
//!    unconditionally (kernel bridge, annotated at the call sites).

use rcad_kernel::geom::{Curve2d, Curve3};
use rcad_kernel::precision::PCONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};
use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_build::brep_tool::topexp_explorer;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::edge::ShapeFixEdge;
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::{
    brep_tool_degenerated, brep_tool_pnt, brep_tool_same_parameter, shape_oriented,
};

// ---------------------------------------------------------------------------
// Shared BRep_Tool / BRep_Builder / BRepTools re-hosts (bridge #1; the
// wire_statics.rs precedents).  Used by split_tool / intersection_tool /
// face / shell.
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Surface(face, L) — the face surface and its location.
pub(crate) fn brep_tool_surface_loc(brep: &BRep, face: &Shape) -> (Option<rcad_kernel::geom::Surface3>, u32) {
    match brep.tshapes[face.index].as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fd.surface_location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Surface(face).
pub(crate) fn brep_tool_surface(brep: &BRep, face: &Shape) -> Option<rcad_kernel::geom::Surface3> {
    brep_tool_surface_loc(brep, face).0
}

/// OCCT BRep_Tool::Range(edge, first, last).
pub(crate) fn brep_tool_range(brep: &BRep, edge: &Shape) -> (f64, f64) {
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Curve(edge, L, f, l) — the 3d curve with its location.
pub(crate) fn brep_tool_curve_loc(brep: &BRep, edge: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), edge.location, ed.range[0], ed.range[1]),
        _ => (None, 0, 0.0, 0.0),
    }
}

/// OCCT BRepTools::Update(const TopoDS_Edge&) (BRepTools.cxx L375): empty.
pub(crate) fn brep_tools_update_edge(_brep: &mut BRep, _e: &Shape) {}

/// OCCT BRepTools::Update(const TopoDS_Face&) (BRepTools.cxx L383-390):
/// `if (!F.Checked()) { UpdateFaceUVPoints(F); F.TShape()->Checked(true); }`.
/// The UpdateFaceUVPoints pass (BRepTools.cxx L485-519) recomputes the stored
/// first/last UV points of every pcurve representation; the rcad
/// CurveRepresentation computes the UV points on read (kernel architecture),
/// so the pass reduces to the Checked flag.
pub(crate) fn brep_tools_update_face(brep: &mut BRep, f: &Shape) {
    let flags = match brep.tshapes.get(f.index).map(|ts| ts.as_ref()) {
        Some(TShape::Face(fd)) => fd.flags,
        _ => return,
    };
    if flags & rcad_kernel::topods::tshape_flags::CHECKED == 0 {
        if let TShape::Face(fd) = brep.tshapes[f.index].as_ref() {
            // SAFETY: the in-place TShape mutation through the shared pool
            // handle (the BRep_TShape pointer write of the OCCT builder).
            let ptr = std::sync::Arc::as_ptr(&brep.tshapes[f.index]) as *mut TShape;
            let _ = fd;
            let ts = unsafe { &mut *ptr };
            if let TShape::Face(fd2) = ts {
                fd2.flags |= rcad_kernel::topods::tshape_flags::CHECKED;
            }
        }
    }
}

/// OCCT BRepTools::Update(const TopoDS_Shape&) — the dispatch over the shape
/// type (BRepTools.cxx L371-429): Vertex/Edge/Wire are empty; Face is the
/// checked-flag pass; Shell/Solid/CompSolid/Compound explore down to faces.
pub(crate) fn brep_tools_update(brep: &mut BRep, s: &Shape) {
    match s.data.as_ref() {
        TShape::Vertex(_) | TShape::Edge(_) | TShape::Wire(_) => {}
        TShape::Face(_) => brep_tools_update_face(brep, s),
        TShape::Shell(_) | TShape::Solid(_) | TShape::CompSolid(_) | TShape::Compound(_) => {
            for f in topexp_explorer(brep, s, ShapeType::Face) {
                brep_tools_update_face(brep, &f);
            }
        }
    }
}

/// OCCT BRepTools::Compare(V1, V2) (BRepTools.cxx L522-538): the vertex
/// sameness by point distance against the tolerances.
pub(crate) fn brep_tools_compare(v1: &Shape, v2: &Shape) -> bool {
    if v1.is_same(&v2) {
        return true;
    }
    let p1 = brep_tool_pnt(v1);
    let p2 = brep_tool_pnt(v2);
    let l = p1.distance(p2);
    if l <= brep_tool_tolerance(v1) {
        return true;
    }
    l <= brep_tool_tolerance(v2)
}

/// OCCT TopExp::Vertices(edge, Vfirst, Vlast, CumOri) (TopExp.cxx L144-180):
/// the FORWARD-oriented vertex child becomes Vfirst, the REVERSED one
/// Vlast; `CumOri` composes the parent orientation (TopoDS_Iterator).
pub(crate) fn topexp_vertices_edge(e: &Shape) -> (Shape, Shape) {
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();
    if let TShape::Edge(ed) = e.data.as_ref() {
        for v in [&ed.first, &ed.last] {
            let mut vv = v.clone();
            // TopoDS_Iterator(E, CumOri): the child carries the composed
            // orientation.
            if e.orientation == Orientation::Reversed {
                vv.orientation = match vv.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }
            if vv.orientation == Orientation::Forward {
                vfirst = vv.clone();
            } else if vv.orientation == Orientation::Reversed {
                vlast = vv;
            }
        }
    }
    (vfirst, vlast)
}

/// OCCT TopExp::Vertices(wire, Vfirst, Vlast) (TopExp.cxx L184-238): the
/// vertex-parity walk over the wire edges — vertices seen once land in the
/// map, vertices seen twice (closed wire) drop out; an empty map means a
/// closed wire (both ends are the last V2), a two-entry map gives the
/// FORWARD and REVERSED ends.
pub(crate) fn topexp_vertices(brep: &mut BRep, w: &Shape) -> (Shape, Shape) {
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();
    // NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> vmap — the
    // IsSame identity keyed by (TShape pointer, location).
    let mut vmap: std::collections::HashMap<(u64, u32), Shape> = std::collections::HashMap::new();
    let mut v2 = Shape::null();
    for e in crate::shhealing::shape_build::brep_tool::iter_subshapes(brep, w, true, true) {
        let (mut v1, mut v2e) = if e.orientation == Orientation::Reversed {
            let (a, b) = topexp_vertices_edge(&e);
            (b, a)
        } else {
            topexp_vertices_edge(&e)
        };
        v2 = v2e.clone();
        v1.orientation = Orientation::Forward;
        v2e.orientation = Orientation::Reversed;
        // if (!vmap.Add(V1)) vmap.Remove(V1); — the parity map.
        let k1 = (v1.ptr_id(), v1.location);
        if vmap.remove(&k1).is_none() {
            vmap.insert(k1, v1);
        }
        let k2 = (v2e.ptr_id(), v2e.location);
        if vmap.remove(&k2).is_none() {
            vmap.insert(k2, v2e);
        }
    }
    if vmap.is_empty() {
        // closed
        vfirst = shape_oriented(&v2, Orientation::Forward);
        vlast = shape_oriented(&v2, Orientation::Reversed);
    } else if vmap.len() == 2 {
        // open
        for v in vmap.values() {
            if v.orientation == Orientation::Forward {
                vfirst = v.clone();
            }
        }
        for v in vmap.values() {
            if v.orientation == Orientation::Reversed {
                vlast = v.clone();
            }
        }
    }
    (vfirst, vlast)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_SplitTool.cxx L42 — the empty constructor.
// ---------------------------------------------------------------------------

/// OCCT ShapeFix_SplitTool (hxx L33-88): tool for splitting and cutting
/// edges; includes methods used in OverlappingTool and IntersectionTool.
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeFixSplitTool;

impl ShapeFixSplitTool {
    /// OCCT ShapeFix_SplitTool::ShapeFix_SplitTool() (cxx L42): empty.
    pub fn new() -> Self {
        ShapeFixSplitTool
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_SplitTool.cxx L46-145 — SplitEdge(edge,param,vert,face,
    // newE1,newE2,tol3d,tol2d).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_SplitTool::SplitEdge (cxx L46-145): splits the edge on
    /// two new edges using the new vertex `vert` and `param` — the parameter
    /// for splitting; `face` is necessary for pcurves and
    /// TransferParametersProj.
    #[allow(clippy::too_many_arguments)]
    pub fn split_edge(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        param: f64,
        vert: &Shape,
        face: &Shape,
        new_e1: &mut Shape,
        new_e2: &mut Shape,
        tol3d: f64,
        tol2d: f64,
    ) -> bool {
        // OCCT L55-58: sae.PCurve(edge, face, c2d, a, b, true).
        let sae = ShapeAnalysisEdge::new();
        let mut a = 0.0f64;
        let mut b = 0.0f64;
        let mut c2d: Option<Curve2d> = None;
        sae.pcurve_face(brep, edge, face, &mut c2d, &mut a, &mut b, true);
        // OCCT L59-62.
        if (a - param).abs() < tol2d || (b - param).abs() < tol2d {
            return false;
        }
        // OCCT L64: check distance between edge and new vertex.
        let p1;
        if brep_tool_same_parameter(edge) {
            // OCCT L68-74: c3d = BRep_Tool::Curve(edge, L, f, l).
            let (c3d, l, _f, _l) = brep_tool_curve_loc(brep, edge);
            let c3d = match c3d {
                Some(c) => c,
                // OCCT L70-73: return false.
                None => return false,
            };
            let mut val = rcad_kernel::geom::CurveEval::point_at(&c3d, param);
            // OCCT L75-78.
            if l != 0 {
                val = brep_loc_transform(brep, l, val);
            }
            p1 = val;
        } else {
            // OCCT L82-88: surf = BRep_Tool::Surface(face, L); sas->Value.
            let (surf, l) = brep_tool_surface_loc(brep, face);
            let sas = ShapeAnalysisSurface::new(surf.unwrap());
            let p2d = rcad_kernel::geom::Curve2dEval::point_at(c2d.as_ref().unwrap(), param);
            let mut val = sas.value(p2d);
            if l != 0 {
                val = brep_loc_transform(brep, l, val);
            }
            p1 = val;
        }
        // OCCT L90-96.
        let p2 = brep_tool_pnt(vert);
        if p1.distance(p2) > tol3d {
            // OCCT L94-95: BRep_Builder B; B.UpdateVertex(vert, P1.Distance(P2)).
            let mut b_builder = BRepBuilder::new();
            b_builder.update_vertex_tolerance(brep, vert.clone(), p1.distance(p2));
        }

        // OCCT L98-112.
        let mut transfer_parameters = ShapeAnalysisTransferParametersProj::new();
        transfer_parameters.base.set_max_tolerance(tol3d);
        transfer_parameters.init(brep, edge, face);
        let (first, last) = if a < b { (a, b) } else { (b, a) };

        // OCCT L114-121.
        let sbe = ShapeBuildEdge;
        let mut sfe = ShapeFixEdge::new();
        let orient = edge.orientation;
        let mut b_builder = BRepBuilder::new();
        let w_e = shape_oriented(edge, Orientation::Forward);
        let a_tmp_shape = shape_oriented(vert, Orientation::Reversed); // for porting
        *new_e1 = sbe.copy_replace_vertices(brep, &w_e, &sae.first_vertex(brep, &w_e), &a_tmp_shape);
        // OCCT L122.
        sbe.copy_pcurves(brep, new_e1, &w_e);
        // OCCT L123.
        transfer_parameters.transfer_range(brep, new_e1, first, param, true);
        // OCCT L124-125.
        b_builder.set_edge_same_range(brep, new_e1.clone(), false);
        sfe.fix_same_parameter(brep, new_e1, 0.0);
        // B.SameParameter(newE1,false);
        // OCCT L127-132.
        let a_tmp_shape = shape_oriented(vert, Orientation::Forward);
        *new_e2 = sbe.copy_replace_vertices(brep, &w_e, &a_tmp_shape, &sae.last_vertex(brep, &w_e));
        sbe.copy_pcurves(brep, new_e2, &w_e);
        transfer_parameters.transfer_range(brep, new_e2, param, last, true);
        b_builder.set_edge_same_range(brep, new_e2.clone(), false);
        sfe.fix_same_parameter(brep, new_e2, 0.0);
        // B.SameParameter(newE2,false);

        // OCCT L135-142.
        (*new_e1).orientation = orient;
        (*new_e2).orientation = orient;
        if orient == Orientation::Reversed {
            std::mem::swap(new_e1, new_e2);
        }

        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_SplitTool.cxx L149-202 — SplitEdge(edge,param1,param2,
    // vert,face,newE1,newE2,tol3d,tol2d).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_SplitTool::SplitEdge (cxx L149-202): splits the edge on
    /// two new edges using the new vertex `vert` and `param1`/`param2` — the
    /// parameters for splitting and cutting.
    #[allow(clippy::too_many_arguments)]
    pub fn split_edge_two_params(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        param1: f64,
        param2: f64,
        vert: &Shape,
        face: &Shape,
        new_e1: &mut Shape,
        new_e2: &mut Shape,
        tol3d: f64,
        tol2d: f64,
    ) -> bool {
        // OCCT L159.
        let param = (param1 + param2) / 2.0;
        // OCCT L160-200.
        if self.split_edge(brep, edge, param, vert, face, new_e1, new_e2, tol3d, tol2d) {
            // OCCT L162: cut new edges by param1 and param2.
            let mut is_cut_line = false;
            let mut fp1 = 0.0f64;
            let mut lp1 = 0.0f64;
            let mut fp2 = 0.0f64;
            let mut lp2 = 0.0f64;
            let sae = ShapeAnalysisEdge::new();
            let mut crv1: Option<Curve2d> = None;
            let mut crv2: Option<Curve2d> = None;
            if sae.pcurve_face(brep, new_e1, face, &mut crv1, &mut fp1, &mut lp1, false) {
                if sae.pcurve_face(brep, new_e2, face, &mut crv2, &mut fp2, &mut lp2, false) {
                    if lp1 == param {
                        if (lp1 - fp1) * (lp1 - param1) > 0.0 {
                            self.cut_edge(brep, new_e1, fp1, param1, face, &mut is_cut_line);
                            self.cut_edge(brep, new_e2, lp2, param2, face, &mut is_cut_line);
                        } else {
                            self.cut_edge(brep, new_e1, fp1, param2, face, &mut is_cut_line);
                            self.cut_edge(brep, new_e2, lp2, param1, face, &mut is_cut_line);
                        }
                    } else {
                        if (fp1 - lp1) * (fp1 - param1) > 0.0 {
                            self.cut_edge(brep, new_e1, lp1, param1, face, &mut is_cut_line);
                            self.cut_edge(brep, new_e2, fp2, param2, face, &mut is_cut_line);
                        } else {
                            self.cut_edge(brep, new_e1, lp1, param2, face, &mut is_cut_line);
                            self.cut_edge(brep, new_e2, fp2, param1, face, &mut is_cut_line);
                        }
                    }
                }
            }
            return true;
        }
        false
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_SplitTool.cxx L206-302 — CutEdge.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_SplitTool::CutEdge (cxx L206-302): cuts the edge by the
    /// parameters `pend` and `cut`.
    pub fn cut_edge(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        pend: f64,
        cut: f64,
        face: &Shape,
        iscutline: &mut bool,
    ) -> bool {
        // OCCT L212-215.
        if (cut - pend).abs() < 10.0 * PCONFUSION {
            return false;
        }
        // OCCT L216-218.
        let a_range = (cut - pend).abs();
        let (a0, b0) = brep_tool_range(brep, edge);
        let (mut a, mut b) = (a0, b0);
        // OCCT L219.
        *iscutline = false;
        // OCCT L220-223.
        if a_range < 10.0 * PCONFUSION {
            return false;
        }

        // OCCT L225: case pcurve is trimm of line.
        if !brep_tool_same_parameter(edge) {
            let sae = ShapeAnalysisEdge::new();
            let mut crv: Option<Curve2d> = None;
            let mut fp = 0.0f64;
            let mut lp = 0.0f64;
            if sae.pcurve_face(brep, edge, face, &mut crv, &mut fp, &mut lp, false) {
                // OCCT L233: Crv->IsKind(Geom2d_TrimmedCurve).
                if let Some(Curve2d::Trimmed(tc)) = crv.as_ref() {
                    // OCCT L236: tc->BasisCurve()->IsKind(Geom2d_Line).
                    if matches!(tc.curve.as_ref(), Curve2d::Line(_)) {
                        let mut b_builder = BRepBuilder::new();
                        // OCCT L239.
                        b_builder.set_edge_range(brep, edge.clone(), pend.min(cut), pend.max(cut));
                        if (pend - lp).abs() < PCONFUSION {
                            // cut from the beginning
                            // OCCT L242-249.
                            let cut3d = (cut - fp) * (b - a) / (lp - fp);
                            if cut3d <= PCONFUSION {
                                return false;
                            }
                            b_builder.set_edge_range(brep, edge.clone(), a + cut3d, b);
                            *iscutline = true;
                        } else if (pend - fp).abs() < PCONFUSION {
                            // cut from the end
                            // OCCT L252-258.
                            let cut3d = (lp - cut) * (b - a) / (lp - fp);
                            if cut3d <= PCONFUSION {
                                return false;
                            }
                            b_builder.set_edge_range(brep, edge.clone(), a, b - cut3d);
                            *iscutline = true;
                        }
                    }
                }
            }
            // OCCT L263.
            return true;
        }

        // OCCT L266-270: det-study on 03/12/01 checking the old and new ranges.
        if ((a - b).abs() - a_range).abs() < PCONFUSION {
            return false;
        }
        // OCCT L271-274.
        if a_range < 10.0 * PCONFUSION {
            return false;
        }

        // OCCT L276-280.
        let _c3d = brep_tool_curve_loc(brep, edge).0;
        let sac = ShapeAnalysisCurve;
        a = pend.min(cut);
        b = pend.max(cut);
        let (mut na, mut nb) = (a, b);

        // OCCT L282-295.
        let mut b_builder = BRepBuilder::new();
        if !brep_tool_degenerated(edge) && _c3d.is_some() && {
            let c = _c3d.as_ref().unwrap();
            sac.validate_range(c, &mut na, &mut nb, PCONFUSION) && (na != a || nb != b)
        } {
            b_builder.set_edge_range(brep, edge.clone(), na, nb);
            let sae = ShapeAnalysisEdge::new();
            if sae.has_pcurve_face(brep, edge, face) {
                b_builder.set_edge_same_range(brep, edge.clone(), false);
            }

            let mut sfe = ShapeFixEdge::new();
            sfe.fix_same_parameter(brep, edge, 0.0);
        } else {
            // OCCT L298.
            b_builder.set_edge_range(brep, edge.clone(), a, b);
        }

        // OCCT L301.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_SplitTool.cxx L306-579 — SplitEdge(edge,fp,V1,lp,V2,face,
    // SeqE,aNum,context,tol3d,tol2d).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_SplitTool::SplitEdge (cxx L306-579): splits the edge on
    /// two new edges using two new vertices V1 and V2 and two parameters fp
    /// and lp; `aNum` is the number of the edge in SeqE corresponding to
    /// [fp,lp].
    #[allow(clippy::too_many_arguments)]
    pub fn split_edge_sequence(
        &self,
        brep: &mut BRep,
        edge: &Shape,
        fp: f64,
        v1: &Shape,
        lp: f64,
        v2: &Shape,
        face: &Shape,
        seq_e: &mut Vec<Shape>,
        a_num: &mut i32,
        context: &ShapeBuildReShape,
        tol3d: f64,
        tol2d: f64,
    ) -> bool {
        // OCCT L318-321.
        if (lp - fp).abs() < tol2d {
            return false;
        }
        // OCCT L322-323.
        *a_num = 0;
        seq_e.clear();
        let mut b_builder = BRepBuilder::new();
        // OCCT L325-328.
        let sae = ShapeAnalysisEdge::new();
        let mut a = 0.0f64;
        let mut b = 0.0f64;
        let mut c2d: Option<Curve2d> = None;
        sae.pcurve_face(brep, edge, face, &mut c2d, &mut a, &mut b, true);
        // OCCT L329-338.
        let vf = sae.first_vertex(brep, edge);
        let vl = sae.last_vertex(brep, edge);
        let tol_vf = brep_tool_tolerance(&vf);
        let tol_vl = brep_tool_tolerance(&vl);
        let tol_v1 = brep_tool_tolerance(v1);
        let tol_v2 = brep_tool_tolerance(v2);
        let pvf = brep_tool_pnt(&vf);
        let pvl = brep_tool_pnt(&vl);
        let pv1 = brep_tool_pnt(v1);
        let pv2 = brep_tool_pnt(v2);

        // OCCT L340-352.
        let (par1, par2);
        let mut is_reverse = false;
        if (b - a) * (lp - fp) > 0.0 {
            par1 = fp;
            par2 = lp;
        } else {
            par1 = lp;
            par2 = fp;
            is_reverse = true;
        }

        let mut ctx = context.clone();

        // OCCT L354-424: fabs(a-par1)<=tol2d && fabs(b-par2)<=tol2d.
        if (a - par1).abs() <= tol2d && (b - par2).abs() <= tol2d {
            if is_reverse {
                // OCCT L356-388.
                let mut newtol = tol_vf + pvf.distance(pv2);
                if tol_v2 < newtol {
                    b_builder.update_vertex_tolerance(brep, v2.clone(), newtol);
                }
                if vf.orientation == v2.orientation {
                    ctx.replace(brep, &vf, v2);
                    // vf = v2 (the local copies are not reused after).
                } else {
                    let v2r = shape_oriented(v2, Orientation::Reversed);
                    ctx.replace(brep, &vf, &v2r);
                }
                newtol = tol_vl + pvl.distance(pv1);
                if tol_v1 < newtol {
                    b_builder.update_vertex_tolerance(brep, v1.clone(), newtol);
                }
                if vl.orientation == v1.orientation {
                    ctx.replace(brep, &vl, v1);
                } else {
                    let v1r = shape_oriented(v1, Orientation::Reversed);
                    ctx.replace(brep, &vl, &v1r);
                }
            } else {
                // OCCT L389-421.
                let mut newtol = tol_vf + pvf.distance(pv1);
                if tol_v1 < newtol {
                    b_builder.update_vertex_tolerance(brep, v1.clone(), newtol);
                }
                if vf.orientation == v1.orientation {
                    ctx.replace(brep, &vf, v1);
                } else {
                    let v1r = shape_oriented(v1, Orientation::Reversed);
                    ctx.replace(brep, &vf, &v1r);
                }
                newtol = tol_vl + pvl.distance(pv2);
                if tol_v2 < newtol {
                    b_builder.update_vertex_tolerance(brep, v2.clone(), newtol);
                }
                if vl.orientation == v2.orientation {
                    ctx.replace(brep, &vl, v2);
                } else {
                    let v2r = shape_oriented(v2, Orientation::Reversed);
                    ctx.replace(brep, &vl, &v2r);
                }
            }
            // OCCT L422-423.
            seq_e.push(edge.clone());
            *a_num = 1;
        }

        // OCCT L426-476: fabs(a-par1)<=tol2d && fabs(b-par2)>tol2d.
        if (a - par1).abs() <= tol2d && (b - par2).abs() > tol2d {
            let mut new_e1 = Shape::null();
            let mut new_e2 = Shape::null();
            if is_reverse {
                if !self.split_edge(brep, edge, par2, v1, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                // OCCT L435-449.
                let newtol = tol_vf + pvf.distance(pv2);
                if tol_v2 < newtol {
                    b_builder.update_vertex_tolerance(brep, v2.clone(), newtol);
                }
                if vf.orientation == v2.orientation {
                    ctx.replace(brep, &vf, v2);
                } else {
                    let v2r = shape_oriented(v2, Orientation::Reversed);
                    ctx.replace(brep, &vf, &v2r);
                }
            } else {
                if !self.split_edge(brep, edge, par2, v2, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                // OCCT L457-471.
                let newtol = tol_vf + pvf.distance(pv1);
                if tol_v1 < newtol {
                    b_builder.update_vertex_tolerance(brep, v1.clone(), newtol);
                }
                if vf.orientation == v1.orientation {
                    ctx.replace(brep, &vf, v1);
                } else {
                    let v1r = shape_oriented(v1, Orientation::Reversed);
                    ctx.replace(brep, &vf, &v1r);
                }
            }
            // OCCT L473-475.
            seq_e.push(new_e1);
            seq_e.push(new_e2);
            *a_num = 1;
        }

        // OCCT L478-528: fabs(a-par1)>tol2d && fabs(b-par2)<=tol2d.
        if (a - par1).abs() > tol2d && (b - par2).abs() <= tol2d {
            let mut new_e1 = Shape::null();
            let mut new_e2 = Shape::null();
            if is_reverse {
                if !self.split_edge(brep, edge, par1, v2, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                // OCCT L487-501.
                let newtol = tol_vl + pvl.distance(pv1);
                if tol_v1 < newtol {
                    b_builder.update_vertex_tolerance(brep, v1.clone(), newtol);
                }
                if vl.orientation == v1.orientation {
                    ctx.replace(brep, &vl, v1);
                } else {
                    let v1r = shape_oriented(v1, Orientation::Reversed);
                    ctx.replace(brep, &vl, &v1r);
                }
            } else {
                if !self.split_edge(brep, edge, par1, v1, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                // OCCT L509-523.
                let newtol = tol_vl + pvl.distance(pv2);
                if tol_v2 < newtol {
                    b_builder.update_vertex_tolerance(brep, v2.clone(), newtol);
                }
                if vl.orientation == v2.orientation {
                    ctx.replace(brep, &vl, v2);
                } else {
                    let v2r = shape_oriented(v2, Orientation::Reversed);
                    ctx.replace(brep, &vl, &v2r);
                }
            }
            // OCCT L525-527.
            seq_e.push(new_e1);
            seq_e.push(new_e2);
            *a_num = 2;
        }

        // OCCT L530-559: fabs(a-par1)>tol2d && fabs(b-par2)>tol2d.
        if (a - par1).abs() > tol2d && (b - par2).abs() > tol2d {
            let mut new_e1 = Shape::null();
            let mut new_e2 = Shape::null();
            let mut new_e3 = Shape::null();
            let mut new_e4 = Shape::null();
            if is_reverse {
                if !self.split_edge(brep, edge, par1, v2, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                if !self.split_edge(brep, &new_e2, par2, v1, face, &mut new_e3, &mut new_e4, tol3d, tol2d) {
                    return false;
                }
            } else {
                if !self.split_edge(brep, edge, par1, v1, face, &mut new_e1, &mut new_e2, tol3d, tol2d) {
                    return false;
                }
                if !self.split_edge(brep, &new_e2, par2, v2, face, &mut new_e3, &mut new_e4, tol3d, tol2d) {
                    return false;
                }
            }
            // OCCT L555-558.
            seq_e.push(new_e1);
            seq_e.push(new_e3);
            seq_e.push(new_e4);
            *a_num = 2;
        }

        // OCCT L561-564.
        if *a_num == 0 {
            return false;
        }

        // OCCT L566-571.
        let mut sewd = WireData::new();
        for s in seq_e.iter() {
            sewd.add_edge(s, 0);
        }
        let w = sewd.wire(brep);
        ctx.replace(brep, edge, &w);
        // OCCT L572-576.
        for e in topexp_explorer(brep, &w, ShapeType::Edge) {
            brep_tools_update_edge(brep, &e);
        }

        // OCCT L578.
        true
    }
}

/// OCCT `P1.Transformed(L.Transformation())` — the location transform.
pub(crate) fn brep_loc_transform(brep: &BRep, loc: u32, p: glam::DVec3) -> glam::DVec3 {
    brep.get_location(loc).transform_point3(p)
}
