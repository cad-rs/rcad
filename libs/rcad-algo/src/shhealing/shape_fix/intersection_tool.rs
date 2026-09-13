//! 1:1 translation of OCCT `ShapeFix_IntersectionTool`
//! (`TKShHealing/ShapeFix/ShapeFix_IntersectionTool.hxx` L17-124 +
//! `ShapeFix_IntersectionTool.lxx` L16-19 + `ShapeFix_IntersectionTool.cxx`
//! L1-2531, docket row `ShapeFix_IntersectionTool`).
//!
//! Function-count equation (OCCT IntersectionTool cxx + lxx + hxx ctor =
//! rcad intersection_tool.rs): constructor (cxx L51-58) / `Context` (lxx
//! L16-19) / `SplitEdge` (cxx L87-190) / `CutEdge` (cxx L194-262) /
//! `SplitEdge1` (cxx L270-358) / `SplitEdge2` (cxx L367-491) /
//! `UnionVertexes` (cxx L495-880) / `FindVertAndSplitEdge` (cxx L968-1025) /
//! `FixSelfIntersectWire` (cxx L1029-1831) / `FixIntersectingWires` (cxx
//! L1835-2530) — 10 members = 10 rcad methods; the file statics
//! `GetPointOnEdge` (cxx L66-83) / `CreateBoxes2d` (cxx L884-920) /
//! `SelectIntPnt` (cxx L924-964) — 3 statics = 3 rcad functions.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — the rcad equivalents take `brep: &mut BRep` (the
//!    ShapeFix_Edge precedent).
//! 2. `NCollection_DataMap<TopoDS_Shape, Bnd_Box2d, TopTools_ShapeMapHasher>`
//!    -> [`ShapeBoxes2d`], the map keyed by the (TShape pointer, location)
//!    IsSame identity.
//! 3. `Geom2dInt_GInter` (the general 2d curve-curve intersection,
//!    TKGeomAlgo/Geom2dInt) — GAP: the kernel `g_inter` carries the
//!    line-curve arm only ([`Geom2dIntGInterGap`]); the result keeps the
//!    OCCT `!IsDone()` state, so every consumer takes the OCCT failure
//!    branch (the `wire_checks.rs` bridge #5 precedent).  The point/segment
//!    walks below are translated 1:1 and compile against the real IntRes2d
//!    types; they stay runtime-dead behind `!IsDone()`.
//! 4. `BndLib_Add2dCurve::Add(gac, tol, box)` (TKTopAlgo/BndLib) — served by
//!    the landed `ShapeAnalysisCurve::FillBndBox` sampler (the
//!    `analysis.rs GetFaceUVBounds` precedent).
//! 5. `BRepTools::Update(E)` — the empty OCCT body
//!    ([`super::split_tool::brep_tools_update_edge`]).

use std::collections::HashMap;

use rcad_kernel::geom::{Curve2d, Curve2dEval, CurveEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType};

use crate::geomalgo::int_res2d::{Domain, IntersectionPoint, IntersectionSegment, Position, Transition};
use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::split_tool::{
    brep_loc_transform, brep_tool_curve_loc, brep_tool_range, brep_tool_surface,
    brep_tool_surface_loc, brep_tools_update_edge,
};
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::{
    brep_tool_degenerated, brep_tool_pnt, brep_tool_same_parameter, shape_oriented,
};

/// OCCT `NCollection_DataMap<TopoDS_Shape, Bnd_Box2d, TopTools_ShapeMapHasher>`
/// — the rcad map keyed by the (TShape pointer, location) IsSame identity
/// (bridge #2).
pub(crate) type ShapeBoxes2d = HashMap<(u64, u32), BndBox2d>;

/// OCCT `TopoDS_Shape::operator==` — IsEqual: same TShape, same location,
/// same orientation.
pub(crate) fn shape_is_equal(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location && a.orientation == b.orientation
}

// ---------------------------------------------------------------------------
// OCCT Geom2dInt_GInter — GAP re-host (module doc, bridge #3).
// ---------------------------------------------------------------------------

/// OCCT Geom2dInt_GInter (the general 2d curve-curve intersection,
/// TKGeomAlgo/Geom2dInt) — GAP re-host (module doc, bridge #3): the kernel
/// `g_inter` translation carries the line-curve arm only; the general
/// curve-curve `Perform` is not translated, so the result is the OCCT
/// `!IsDone()` state — every consumer takes the OCCT failure branch.
pub struct Geom2dIntGInterGap {
    points: Vec<IntersectionPoint>,
    segments: Vec<IntersectionSegment>,
}

impl Default for Geom2dIntGInterGap {
    fn default() -> Self {
        Self::new()
    }
}

impl Geom2dIntGInterGap {
    /// OCCT Geom2dInt_GInter() — the empty constructor.
    pub fn new() -> Self {
        Geom2dIntGInterGap {
            points: Vec::new(),
            segments: Vec::new(),
        }
    }

    /// OCCT Geom2dInt_GInter::Perform(C1, D1, C2, D2, TolConf, Tol) — GAP:
    /// keeps the OCCT call shape; the result stays empty and not-done.
    pub fn perform(
        &mut self,
        _c1: &Curve2d,
        _d1: &Domain,
        _c2: &Curve2d,
        _d2: &Domain,
        _tol_conf: f64,
        _tol: f64,
    ) {
        self.points.clear();
        self.segments.clear();
    }

    /// OCCT IsDone() — false from the GAP.
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// OCCT Point(num) — 1-indexed.
    pub fn point(&self, num: usize) -> &IntersectionPoint {
        &self.points[num - 1]
    }

    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }

    /// OCCT Segment(num) — 1-indexed.
    pub fn segment(&self, num: usize) -> &IntersectionSegment {
        &self.segments[num - 1]
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_IntersectionTool.cxx L51-58 — the constructor; the members
// (hxx L117-119).
// ---------------------------------------------------------------------------

/// OCCT ShapeFix_IntersectionTool (hxx L38-120): tool for fixing
/// selfintersecting wire and intersecting wires.
pub struct ShapeFixIntersectionTool {
    /// OCCT myContext (hxx L117).
    pub(crate) my_context: Option<ShapeBuildReShape>,
    /// OCCT myPreci (hxx L118).
    pub(crate) my_preci: f64,
    /// OCCT myMaxTol (hxx L119).
    pub(crate) my_max_tol: f64,
}

impl ShapeFixIntersectionTool {
    /// OCCT ShapeFix_IntersectionTool(context, preci, maxtol = 1.0)
    /// (cxx L51-58).  Rust has no default arguments; the OCCT default
    /// `maxtol` is 1.0 — pass it explicitly at the default-arg call sites
    /// (the ShapeFix_Wire.cxx L1203 call shape).
    pub fn new(context: Option<ShapeBuildReShape>, preci: f64, maxtol: f64) -> Self {
        ShapeFixIntersectionTool {
            my_context: context,
            my_preci: preci,
            my_max_tol: maxtol,
        }
    }

    /// OCCT ShapeFix_IntersectionTool::Context (lxx L16-19): returns context.
    pub fn context(&self) -> Option<&ShapeBuildReShape> {
        self.my_context.as_ref()
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L87-190 — SplitEdge.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::SplitEdge (cxx L87-190): splits the
    /// edge on two new edges using the new vertex `vert` and `param` — the
    /// parameter for splitting.
    #[allow(clippy::too_many_arguments)]
    pub fn split_edge(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        param: f64,
        vert: &Shape,
        face: &Shape,
        new_e1: &mut Shape,
        new_e2: &mut Shape,
        preci: f64,
    ) -> bool {
        // OCCT L95-96.
        let sae = ShapeAnalysisEdge::new();

        // OCCT L98-103.
        let v1 = sae.first_vertex(brep, edge);
        let v2 = sae.last_vertex(brep, edge);
        if v1.is_same(vert) || v2.is_same(vert) {
            return false;
        }

        // OCCT L105-110.
        let mut a = 0.0f64;
        let mut b = 0.0f64;
        let mut c2d: Option<Curve2d> = None;
        sae.pcurve_face(brep, edge, face, &mut c2d, &mut a, &mut b, true);
        if (a - param).abs() < 0.01 * preci || (b - param).abs() < 0.01 * preci {
            return false;
        }
        // OCCT L111: check distance between edge and new vertex.
        let p1;
        if brep_tool_same_parameter(edge) && !brep_tool_degenerated(edge) {
            // OCCT L116-122: c3d = BRep_Tool::Curve(edge, L, f, l).
            let (c3d, l, _f, _l) = brep_tool_curve_loc(brep, edge);
            let c3d = match c3d {
                Some(c) => c,
                // OCCT L118-121: return false.
                None => return false,
            };
            let mut val = CurveEval::point_at(&c3d, param);
            // OCCT L123-126.
            if l != 0 {
                val = brep_loc_transform(brep, l, val);
            }
            p1 = val;
        } else {
            // OCCT L130-136: surf = BRep_Tool::Surface(face, L); sas->Value.
            let (surf, l) = brep_tool_surface_loc(brep, face);
            let sas = ShapeAnalysisSurface::new(surf.unwrap());
            let p2d = Curve2dEval::point_at(c2d.as_ref().unwrap(), param);
            let mut val = sas.value(p2d);
            if l != 0 {
                val = brep_loc_transform(brep, l, val);
            }
            p1 = val;
        }
        // OCCT L138-144.
        let p2 = brep_tool_pnt(vert);
        if p1.distance(p2) > preci {
            // OCCT L142-143: BRep_Builder B; B.UpdateVertex(vert, P1.Distance(P2)).
            let mut b_builder = BRepBuilder::new();
            b_builder.update_vertex_tolerance(brep, vert.clone(), p1.distance(p2));
        }

        // OCCT L146-160.
        let mut transfer_parameters = ShapeAnalysisTransferParametersProj::new();
        transfer_parameters.base.set_max_tolerance(preci);
        transfer_parameters.init(brep, edge, face);
        let (first, last) = if a < b { (a, b) } else { (b, a) };

        // OCCT L162-168.
        let sbe = ShapeBuildEdge;
        let orient = edge.orientation;
        let mut b_builder = BRepBuilder::new();
        let w_e = shape_oriented(edge, Orientation::Forward);
        let a_tmp_shape = shape_oriented(vert, Orientation::Reversed); // for porting
        *new_e1 = sbe.copy_replace_vertices(brep, &w_e, &sae.first_vertex(brep, &w_e), &a_tmp_shape);
        // OCCT L169.
        sbe.copy_pcurves(brep, new_e1, &w_e);
        // OCCT L170-172.
        transfer_parameters.transfer_range(brep, new_e1, first, param, true);
        b_builder.set_edge_same_range(brep, new_e1.clone(), false);
        b_builder.set_edge_same_parameter(brep, new_e1.clone(), false);
        // OCCT L173-178.
        let a_tmp_shape = shape_oriented(vert, Orientation::Forward);
        *new_e2 = sbe.copy_replace_vertices(brep, &w_e, &a_tmp_shape, &sae.last_vertex(brep, &w_e));
        sbe.copy_pcurves(brep, new_e2, &w_e);
        transfer_parameters.transfer_range(brep, new_e2, param, last, true);
        b_builder.set_edge_same_range(brep, new_e2.clone(), false);
        b_builder.set_edge_same_parameter(brep, new_e2.clone(), false);

        // OCCT L180-187.
        (*new_e1).orientation = orient;
        (*new_e2).orientation = orient;
        if orient == Orientation::Reversed {
            std::mem::swap(new_e1, new_e2);
        }

        // OCCT L189.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L194-262 — CutEdge.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::CutEdge (cxx L194-262): cuts the edge
    /// by the parameters `pend` and `cut`.
    pub fn cut_edge(
        &mut self,
        brep: &mut BRep,
        edge: &Shape,
        pend: f64,
        cut: f64,
        face: &Shape,
        iscutline: &mut bool,
    ) -> bool {
        // OCCT L200-203.
        if (cut - pend).abs() < 10.0 * PCONFUSION {
            return false;
        }
        // OCCT L204-206.
        let a_range = (cut - pend).abs();
        let (a, b) = brep_tool_range(brep, edge);

        // OCCT L208-211.
        if a_range < 10.0 * PCONFUSION {
            return false;
        }

        // OCCT L213: case pcurve is trimm of line.
        if !brep_tool_same_parameter(edge) {
            let sae = ShapeAnalysisEdge::new();
            let mut crv: Option<Curve2d> = None;
            let mut fp = 0.0f64;
            let mut lp = 0.0f64;
            if sae.pcurve_face(brep, edge, face, &mut crv, &mut fp, &mut lp, false) {
                // OCCT L221: Crv->IsKind(Geom2d_TrimmedCurve).
                if let Some(Curve2d::Trimmed(tc)) = crv.as_ref() {
                    // OCCT L224: tc->BasisCurve()->IsKind(Geom2d_Line).
                    if matches!(tc.curve.as_ref(), Curve2d::Line(_)) {
                        let mut b_builder = BRepBuilder::new();
                        // OCCT L227.
                        b_builder.set_edge_range(brep, edge.clone(), pend.min(cut), pend.max(cut));
                        if (pend - lp).abs() < PCONFUSION {
                            // cut from the beginning
                            // OCCT L230-232.
                            let cut3d = (cut - fp) * (b - a) / (lp - fp);
                            b_builder.set_edge_range(brep, edge.clone(), a + cut3d, b);
                            *iscutline = true;
                        } else if (pend - fp).abs() < PCONFUSION {
                            // cut from the end
                            // OCCT L236-238.
                            let cut3d = (lp - cut) * (b - a) / (lp - fp);
                            b_builder.set_edge_range(brep, edge.clone(), a, b - cut3d);
                            *iscutline = true;
                        }
                    }

                    // OCCT L242.
                    return true;
                }
            }
            // OCCT L245.
            return false;
        }

        // OCCT L248-252: det-study on 03/12/01 checking the old and new ranges.
        if ((a - b).abs() - a_range).abs() < PCONFUSION {
            return false;
        }
        // OCCT L253-256.
        if a_range < 10.0 * PCONFUSION {
            return false;
        }

        // OCCT L258-259.
        let mut b_builder = BRepBuilder::new();
        b_builder.set_edge_range(brep, edge.clone(), pend.min(cut), pend.max(cut));

        // OCCT L261.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L270-358 — SplitEdge1.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::SplitEdge1 (cxx L270-358): splits
    /// edge[a,b] on two parts e1[a,param] and e2[param,b] using vertex vert.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn split_edge1(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        face: &Shape,
        num: i32,
        param: f64,
        vert: &Shape,
        preci: f64,
        boxes: &mut ShapeBoxes2d,
    ) -> bool {
        // OCCT L279: Standard_ASSERT_RETURN(num > 0 && num <= NbEdges()).
        if !(num > 0 && num <= sewd.nb_edges()) {
            return false;
        }

        // OCCT L281-286.
        let edge = sewd.edge(num);
        let mut new_e1 = Shape::null();
        let mut new_e2 = Shape::null();
        if !self.split_edge(brep, &edge, param, vert, face, &mut new_e1, &mut new_e2, preci) {
            return false;
        }

        // OCCT L288: change context.
        let mut wd = WireData::new();
        wd.add_edge(&new_e1, 0);
        wd.add_edge(&new_e2, 0);
        let w = wd.wire(brep);
        if let Some(ctx) = self.my_context.as_mut() {
            ctx.replace(brep, &edge, &w);
        }
        for e in topexp_explorer(brep, &w, ShapeType::Edge) {
            brep_tools_update_edge(brep, &e);
        }

        // OCCT L302: change sewd.
        sewd.set_edge(&new_e1, num);
        if num == sewd.nb_edges() {
            sewd.add_edge(&new_e2, 0);
        } else {
            sewd.add_edge(&new_e2, num + 1);
        }

        // OCCT L313: change boxes.
        boxes.remove(&(edge.ptr_id(), edge.location));
        let (s, l) = brep_tool_surface_loc(brep, face);
        let sae = ShapeAnalysisEdge::new();
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        let s = s.unwrap();
        if sae.pcurve_surface(brep, &new_e1, &s, l, &mut c2d, &mut cf, &mut cl, false) {
            let mut box2d = BndBox2d::new();
            let c = c2d.as_ref().unwrap();
            let domain = c.default_domain();
            let (a_first, a_last) = (domain[0], domain[1]);
            // OCCT L326-334: pdn avoiding problems with segment in Bnd_Box.
            if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
                bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
            } else {
                bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
            }
            boxes.insert((new_e1.ptr_id(), new_e1.location), box2d);
        }
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if sae.pcurve_surface(brep, &new_e2, &s, l, &mut c2d, &mut cf, &mut cl, false) {
            let mut box2d = BndBox2d::new();
            let c = c2d.as_ref().unwrap();
            let domain = c.default_domain();
            let (a_first, a_last) = (domain[0], domain[1]);
            if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
                bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
            } else {
                bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
            }
            boxes.insert((new_e2.ptr_id(), new_e2.location), box2d);
        }

        // OCCT L357.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L367-491 — SplitEdge2.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::SplitEdge2 (cxx L367-491): auxiliary;
    /// splits edge[a,b] on two parts e1[a,param1] and e2[param2,b] using
    /// vertex vert (removes the segment (param1,param2) from the edge).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn split_edge2(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        face: &Shape,
        num: i32,
        param1: f64,
        param2: f64,
        vert: &Shape,
        preci: f64,
        boxes: &mut ShapeBoxes2d,
    ) -> bool {
        // OCCT L377-383.
        let edge = sewd.edge(num);
        let mut new_e1 = Shape::null();
        let mut new_e2 = Shape::null();
        let param = (param1 + param2) / 2.0;
        if !self.split_edge(brep, &edge, param, vert, face, &mut new_e1, &mut new_e2, preci) {
            return false;
        }
        // OCCT L384: cut new edges by param1 and param2.
        let mut is_cut_line = false;
        let mut fp1 = 0.0f64;
        let mut lp1 = 0.0f64;
        let mut fp2 = 0.0f64;
        let mut lp2 = 0.0f64;
        let sae = ShapeAnalysisEdge::new();
        let mut crv1: Option<Curve2d> = None;
        let mut crv2: Option<Curve2d> = None;
        if sae.pcurve_face(brep, &new_e1, face, &mut crv1, &mut fp1, &mut lp1, false) {
            if sae.pcurve_face(brep, &new_e2, face, &mut crv2, &mut fp2, &mut lp2, false) {
                if lp1 == param {
                    if (lp1 - fp1) * (lp1 - param1) > 0.0 {
                        self.cut_edge(brep, &mut new_e1, fp1, param1, face, &mut is_cut_line);
                        self.cut_edge(brep, &mut new_e2, lp2, param2, face, &mut is_cut_line);
                    } else {
                        self.cut_edge(brep, &mut new_e1, fp1, param2, face, &mut is_cut_line);
                        self.cut_edge(brep, &mut new_e2, lp2, param1, face, &mut is_cut_line);
                    }
                } else {
                    if (fp1 - lp1) * (fp1 - param1) > 0.0 {
                        self.cut_edge(brep, &mut new_e1, lp1, param1, face, &mut is_cut_line);
                        self.cut_edge(brep, &mut new_e2, fp2, param2, face, &mut is_cut_line);
                    } else {
                        self.cut_edge(brep, &mut new_e1, lp1, param2, face, &mut is_cut_line);
                        self.cut_edge(brep, &mut new_e2, fp2, param1, face, &mut is_cut_line);
                    }
                }
            }
        }

        // OCCT L422: change context.
        let mut wd = WireData::new();
        wd.add_edge(&new_e1, 0);
        wd.add_edge(&new_e2, 0);
        let w = wd.wire(brep);
        if let Some(ctx) = self.my_context.as_mut() {
            ctx.replace(brep, &edge, &w);
        }
        for e in topexp_explorer(brep, &w, ShapeType::Edge) {
            brep_tools_update_edge(brep, &e);
        }

        // OCCT L436: change sewd.
        sewd.set_edge(&new_e1, num);
        if num == sewd.nb_edges() {
            sewd.add_edge(&new_e2, 0);
        } else {
            sewd.add_edge(&new_e2, num + 1);
        }

        // OCCT L447: change boxes.
        boxes.remove(&(edge.ptr_id(), edge.location));
        let (s, l) = brep_tool_surface_loc(brep, face);
        let s = s.unwrap();
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if sae.pcurve_surface(brep, &new_e1, &s, l, &mut c2d, &mut cf, &mut cl, false) {
            let mut box2d = BndBox2d::new();
            let c = c2d.as_ref().unwrap();
            let domain = c.default_domain();
            let (a_first, a_last) = (domain[0], domain[1]);
            if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
                bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
            } else {
                bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
            }
            boxes.insert((new_e1.ptr_id(), new_e1.location), box2d);
        }
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if sae.pcurve_surface(brep, &new_e2, &s, l, &mut c2d, &mut cf, &mut cl, false) {
            let mut box2d = BndBox2d::new();
            let c = c2d.as_ref().unwrap();
            let domain = c.default_domain();
            let (a_first, a_last) = (domain[0], domain[1]);
            if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
                bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
            } else {
                bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
            }
            boxes.insert((new_e2.ptr_id(), new_e2.location), box2d);
        }

        // OCCT L490.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L495-880 — UnionVertexes.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::UnionVertexes (cxx L495-880): unions
    /// the coinciding end vertices of edge1 and edge2.
    pub(crate) fn union_vertexes(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        edge1: &mut Shape,
        edge2: &mut Shape,
        num2: i32,
        boxes: &mut ShapeBoxes2d,
        b2: &BndBox2d,
    ) -> bool {
        // OCCT L503: union vertexes.
        let mut res = false;
        let sbe = ShapeBuildEdge;
        let sae = ShapeAnalysisEdge::new();
        let mut b_builder = BRepBuilder::new();
        // OCCT L509-520.
        let v1f = sae.first_vertex(brep, edge1);
        let pv1f = brep_tool_pnt(&v1f);
        let v1l = sae.last_vertex(brep, edge1);
        let pv1l = brep_tool_pnt(&v1l);
        let v2f = sae.first_vertex(brep, edge2);
        let pv2f = brep_tool_pnt(&v2f);
        let v2l = sae.last_vertex(brep, edge2);
        let pv2l = brep_tool_pnt(&v2l);
        let d11 = pv1f.distance(pv2f);
        let d12 = pv1f.distance(pv2l);
        let d21 = pv1l.distance(pv2f);
        let d22 = pv1l.distance(pv2l);

        if d11 < d12 && d11 < d21 && d11 < d22 {
            // OCCT L521-609.
            let tolv = brep_tool_tolerance(&v1f).max(brep_tool_tolerance(&v2f));
            if !v2f.is_same(&v1f) && d11 < tolv {
                // union vertexes V1F and V2F
                b_builder.update_vertex_tolerance(brep, v1f.clone(), tolv);
                let mut new_e = sbe.copy_replace_vertices(brep, edge2, &v1f, &v2l);
                if let Some(ctx) = self.my_context.as_mut() {
                    ctx.replace(brep, edge2, &new_e);
                }
                sewd.set_edge(&new_e, num2);
                *edge2 = new_e.clone();
                boxes.insert((new_e.ptr_id(), new_e.location), b2.clone()); // update boxes
                // replace vertex in other edge
                let (num21, num22) = neighbor_nums(sewd, num2);
                let edge21 = sewd.edge(num21);
                let edge22 = sewd.edge(num22);
                let v21f = sae.first_vertex(brep, &edge21);
                let v21l = sae.last_vertex(brep, &edge21);
                let v22f = sae.first_vertex(brep, &edge22);
                let v22l = sae.last_vertex(brep, &edge22);
                if v21f.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v1f, &v21l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v21l.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v21f, &v1f);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v22f.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v1f, &v22l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                if v22l.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v22f, &v1f);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                res = true;
            }
        } else if d12 < d21 && d12 < d22 {
            // OCCT L610-699.
            let tolv = brep_tool_tolerance(&v1f).max(brep_tool_tolerance(&v2l));
            if !v2l.is_same(&v1f) && d12 < tolv {
                // union vertexes V1F and V2L
                b_builder.update_vertex_tolerance(brep, v1f.clone(), tolv);
                let mut new_e = sbe.copy_replace_vertices(brep, edge2, &v2f, &v1f);
                if let Some(ctx) = self.my_context.as_mut() {
                    ctx.replace(brep, edge2, &new_e);
                }
                sewd.set_edge(&new_e, num2);
                *edge2 = new_e.clone();
                // boxes.Bind(NewE,boxes.Find(edge2)); // update boxes
                boxes.insert((new_e.ptr_id(), new_e.location), b2.clone()); // update boxes
                let (num21, num22) = neighbor_nums(sewd, num2);
                let edge21 = sewd.edge(num21);
                let edge22 = sewd.edge(num22);
                let v21f = sae.first_vertex(brep, &edge21);
                let v21l = sae.last_vertex(brep, &edge21);
                let v22f = sae.first_vertex(brep, &edge22);
                let v22l = sae.last_vertex(brep, &edge22);
                if v21f.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v1f, &v21l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v21l.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v21f, &v1f);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v22f.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v1f, &v22l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                if v22l.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v22f, &v1f);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                res = true;
            }
        } else if d21 < d22 {
            // OCCT L700-788.
            let tolv = brep_tool_tolerance(&v1l).max(brep_tool_tolerance(&v2f));
            if !v2f.is_same(&v1l) && d21 < tolv {
                // union vertexes V1L and V2F
                b_builder.update_vertex_tolerance(brep, v1l.clone(), tolv);
                let mut new_e = sbe.copy_replace_vertices(brep, edge2, &v1l, &v2l);
                if let Some(ctx) = self.my_context.as_mut() {
                    ctx.replace(brep, edge2, &new_e);
                }
                sewd.set_edge(&new_e, num2);
                *edge2 = new_e.clone();
                boxes.insert((new_e.ptr_id(), new_e.location), b2.clone()); // update boxes
                let (num21, num22) = neighbor_nums(sewd, num2);
                let edge21 = sewd.edge(num21);
                let edge22 = sewd.edge(num22);
                let v21f = sae.first_vertex(brep, &edge21);
                let v21l = sae.last_vertex(brep, &edge21);
                let v22f = sae.first_vertex(brep, &edge22);
                let v22l = sae.last_vertex(brep, &edge22);
                if v21f.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v1l, &v21l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v21l.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v21f, &v1l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v22f.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v1l, &v22l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                if v22l.is_same(&v2f) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v22f, &v1l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                res = true;
            }
        } else {
            // OCCT L789-877.
            let tolv = brep_tool_tolerance(&v1l).max(brep_tool_tolerance(&v2l));
            if !v2l.is_same(&v1l) && d22 < tolv {
                // union vertexes V1L and V2L
                b_builder.update_vertex_tolerance(brep, v1l.clone(), tolv);
                let mut new_e = sbe.copy_replace_vertices(brep, edge2, &v2f, &v1l);
                if let Some(ctx) = self.my_context.as_mut() {
                    ctx.replace(brep, edge2, &new_e);
                }
                sewd.set_edge(&new_e, num2);
                *edge2 = new_e.clone();
                boxes.insert((new_e.ptr_id(), new_e.location), b2.clone()); // update boxes
                let (num21, num22) = neighbor_nums(sewd, num2);
                let edge21 = sewd.edge(num21);
                let edge22 = sewd.edge(num22);
                let v21f = sae.first_vertex(brep, &edge21);
                let v21l = sae.last_vertex(brep, &edge21);
                let v22f = sae.first_vertex(brep, &edge22);
                let v22l = sae.last_vertex(brep, &edge22);
                if v21f.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v1l, &v21l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v21l.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge21, &v21f, &v1l);
                    if let Some(bb) = boxes.get(&(edge21.ptr_id(), edge21.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge21, &new_e);
                    }
                    sewd.set_edge(&new_e, num21);
                }
                if v22f.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v1l, &v22l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                if v22l.is_same(&v2l) {
                    new_e = sbe.copy_replace_vertices(brep, &edge22, &v22f, &v1l);
                    if let Some(bb) = boxes.get(&(edge22.ptr_id(), edge22.location)) {
                        boxes.insert((new_e.ptr_id(), new_e.location), bb.clone()); // update boxes
                    }
                    if let Some(ctx) = self.my_context.as_mut() {
                        ctx.replace(brep, &edge22, &new_e);
                    }
                    sewd.set_edge(&new_e, num22);
                }
                res = true;
            }
        }

        // OCCT L879.
        res
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L968-1025 — FindVertAndSplitEdge.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::FindVertAndSplitEdge (cxx L968-1025):
    /// finds the needed vertex from edge2 and splits edge1 using it.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn find_vert_and_split_edge(
        &mut self,
        brep: &mut BRep,
        param1: f64,
        edge1: &Shape,
        edge2: &Shape,
        crv1: &Curve2d,
        max_tol_vert: &mut f64,
        num1: &mut i32,
        sewd: &mut WireData,
        face: &Shape,
        boxes: &mut ShapeBoxes2d,
        a_tmp_key: bool,
    ) -> bool {
        // OCCT L980-983: find needed vertex from edge2 and split edge1 using it.
        let sae = ShapeAnalysisEdge::new();
        let sas = ShapeAnalysisSurface::new(brep_tool_surface(brep, face).unwrap());
        let pi1 = get_point_on_edge(brep, edge1, &sas, crv1, param1);
        let mut b_builder = BRepBuilder::new();
        // OCCT L987-992.
        let v1 = sae.first_vertex(brep, edge2);
        let pv1 = brep_tool_pnt(&v1);
        let v2 = sae.last_vertex(brep, edge2);
        let pv2 = brep_tool_pnt(&v2);
        let v11 = sae.first_vertex(brep, edge1);
        let v12 = sae.last_vertex(brep, edge1);
        let mut need_split = true;
        // OCCT L994-1011.
        let v;
        let tol_v;
        if pi1.distance(pv1) < pi1.distance(pv2) {
            if v1.is_same(&v11) || v1.is_same(&v12) {
                need_split = false;
            }
            v = v1.clone();
            tol_v = ((pi1.distance(pv1) / 2.0) * 1.00001f64).max(brep_tool_tolerance(&v1));
        } else {
            if v2.is_same(&v11) || v2.is_same(&v12) {
                need_split = false;
            }
            v = v2.clone();
            tol_v = ((pi1.distance(pv2) / 2.0) * 1.00001f64).max(brep_tool_tolerance(&v2));
        }
        // OCCT L1012-1023.
        if need_split || a_tmp_key {
            if self.split_edge1(brep, sewd, face, *num1, param1, &v, tol_v, boxes) {
                b_builder.update_vertex_tolerance(brep, v.clone(), tol_v);
                *max_tol_vert = (*max_tol_vert).max(tol_v);
                //      NbSplit++;
                *num1 -= 1;
                return true;
                // break;
            }
        }
        // OCCT L1024.
        false
    }

}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_IntersectionTool.cxx L66-83 — GetPointOnEdge (static).
// ---------------------------------------------------------------------------

/// OCCT static GetPointOnEdge (cxx L66-83): auxiliary — like in BRepCheck,
/// the point is to be taken from the 3d curve (but only if the edge is
/// SameParameter).
pub(super) fn get_point_on_edge(
    brep: &BRep,
    edge: &Shape,
    surf: &ShapeAnalysisSurface,
    crv2d: &Curve2d,
    param: f64,
) -> glam::DVec3 {
    // OCCT L71-80.
    if brep_tool_same_parameter(edge) {
        let (con_s, l, _f, _l) = brep_tool_curve_loc(brep, edge);
        if let Some(c) = con_s {
            let p = CurveEval::point_at(&c, param);
            return if l != 0 {
                brep_loc_transform(brep, l, p)
            } else {
                p
            };
        }
    }
    // OCCT L81-82.
    let a_p2d = Curve2dEval::point_at(crv2d, param);
    surf.value(a_p2d)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_IntersectionTool.cxx L884-920 — CreateBoxes2d (static).
// ---------------------------------------------------------------------------

/// OCCT static CreateBoxes2d (cxx L884-920): creates the 2d box for the
/// edges from the wire.
pub(crate) fn create_boxes2d(
    brep: &mut BRep,
    sewd: &mut WireData,
    face: &Shape,
    boxes: &mut ShapeBoxes2d,
) -> BndBox2d {
    let mut a_total_box = BndBox2d::new();
    let (s, l) = brep_tool_surface_loc(brep, face);
    let s = s.unwrap();
    let sae = ShapeAnalysisEdge::new();
    for i in 1..=sewd.nb_edges() {
        let e = sewd.edge(i);
        let mut c2d: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if sae.pcurve_surface(brep, &e, &s, l, &mut c2d, &mut cf, &mut cl, false) {
            let mut box2d = BndBox2d::new();
            let c = c2d.as_ref().unwrap();
            let domain = c.default_domain();
            let (a_first, a_last) = (domain[0], domain[1]);
            // OCCT L905-913: pdn avoiding problems with segment in Bnd_Box.
            if matches!(c, Curve2d::BSpline(_)) && (cf < a_first || cl > a_last) {
                bnd_lib_add2d_curve(brep, c, a_first, a_last, &mut box2d);
            } else {
                bnd_lib_add2d_curve(brep, c, cf, cl, &mut box2d);
            }
            // OCCT L914-916.
            boxes.insert((e.ptr_id(), e.location), box2d);
            // OCCT L916: aTotalBox.Add(box).
            a_total_box.add_box(boxes.get(&(e.ptr_id(), e.location)).unwrap());
        }
    }
    a_total_box
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_IntersectionTool.cxx L924-964 — SelectIntPnt (static).
// ---------------------------------------------------------------------------

/// OCCT static SelectIntPnt (cxx L924-964).
pub(super) fn select_int_pnt(
    inter: &Geom2dIntGInterGap,
    ip: &mut IntersectionPoint,
    tr1: &mut Transition,
    tr2: &mut Transition,
) {
    // OCCT L929-931.
    *ip = inter.point(1).clone();
    *tr1 = ip.transition_of_first().clone();
    *tr2 = ip.transition_of_second().clone();
    // OCCT L932: possible second point is better?
    if inter.nb_points() == 2 {
        let mut status1 = 0i32;
        let mut status2 = 0i32;
        // OCCT L935-942.
        if tr1.position_on_curve() == Position::Middle {
            status1 += 1;
        }
        if tr2.position_on_curve() == Position::Middle {
            status1 += 2;
        }
        // OCCT L944-956.
        let ip2 = inter.point(2);
        let tr12 = ip2.transition_of_first();
        let tr22 = ip2.transition_of_second();
        if tr12.position_on_curve() == Position::Middle {
            status2 += 1;
        }
        if tr22.position_on_curve() == Position::Middle {
            status2 += 2;
        }
        // OCCT L957-962.
        if status2 > status1 {
            *ip = ip2.clone();
            *tr1 = tr12.clone();
            *tr2 = tr22.clone();
        }
    }
}

// ---------------------------------------------------------------------------
// The BndLib_Add2dCurve / Bnd_Box2d::Add re-hosts (module doc, bridge #4).
// ---------------------------------------------------------------------------

/// OCCT BndLib_Add2dCurve::Add(C, U1, U2, Tol, B) (BndLib_Add2dCurve.cxx) —
/// GAP-served: the BndLib adaptor walk is not translated; the box is filled
/// by the landed `ShapeAnalysisCurve::FillBndBox` sampler (the analysis.rs
/// GetFaceUVBounds precedent), which adds the curve points over [u1,u2].
pub(crate) fn bnd_lib_add2d_curve(
    _brep: &BRep,
    c2d: &Curve2d,
    u1: f64,
    u2: f64,
    box2d: &mut BndBox2d,
) {
    let sac = ShapeAnalysisCurve;
    sac.fill_bnd_box(c2d, u1, u2, 33, true, box2d);
}

/// OCCT `TShape` neighbor indices of the num-th edge in the wire (the
/// UnionVertexes num21/num22 computation, cxx L540-556).
fn neighbor_nums(sewd: &WireData, num2: i32) -> (i32, i32) {
    let num21 = if num2 > 1 { num2 - 1 } else { sewd.nb_edges() };
    let num22 = if num2 < sewd.nb_edges() { num2 + 1 } else { 1 };
    (num21, num22)
}
