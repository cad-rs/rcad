//! 1:1 translation of OCCT `ShapeFix` package statics
//! (`TKShHealing/ShapeFix/ShapeFix.cxx` L1-818, docket row `ShapeFix` —
//! the `SameParameter` static is the W1-6 pull-forward, the remaining
//! statics complete the 818-LOC package file ahead of the W3 class batch).
//!
//! Function-count equation (OCCT ShapeFix.cxx = rcad shape_fix.rs):
//! `SameParameter` / `EncodeRegularity` / `RemoveSmallEdges` /
//! `FixVertexPosition` / `LeastEdgeSize` (public statics) +
//! `ReplaceVertex` / `getNearPoint` / `getNearestEdges` (file statics) —
//! 8 OCCT functions = 8 rcad functions.
//!
//! Architecture bridges (documented per call site):
//! - `Message_ProgressScope` — rcad carries no progress/abort channel; the
//!   scopes reduce to the [`MessageProgressScope`] no-abort bridge whose
//!   `More()` is always true, which keeps OCCT's non-aborted control flow
//!   exactly.
//! - `occ::handle<ShapeExtend_BasicMsgRegistrator>` — an
//!   `Option<&mut dyn BasicMsgRegistrator>` (the shape_extend trait).
//! - `Geom2dAdaptor_Curve` + `Adaptor3d_CurveOnSurface` +
//!   `GeomAdaptor_Surface(plane)` — the `Adaptor3d_CurveOnSurface::Value`
//!   chain reduces to `plane.point_at(pcurve.point_at(par))` over the
//!   location-applied face surface and the pcurve fetched by
//!   `BRep_Tool::CurveOnSurface`.
//! - `Bnd_Box` (LeastEdgeSize) — the three added points' component-wise
//!   min/max box, identical to OCCT's Bnd_Box::Get output.

use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType, TShape};

use crate::shhealing::shape_build::brep_tool::{iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::basic_msg_registrator::BasicMsgRegistrator;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::ShapeExtendStatus;
use crate::shhealing::shape_fix::edge::ShapeFixEdge;
use crate::shhealing::shape_fix::shape_fix_shape::ShapeFixShape;

/// OCCT Message_ProgressScope — the rcad architecture bridge: no abort
/// channel exists, so `More()` is always true and the while-loops of
/// SameParameter run to completion exactly as in OCCT's non-aborted path.
#[derive(Clone, Copy)]
pub struct MessageProgressScope;

impl MessageProgressScope {
    pub fn more(&self) -> bool {
        true
    }
    pub fn next(&self) {}
}

/// OCCT `Message_ProgressRange` — the no-abort bridge handle.
#[derive(Clone, Copy, Default)]
pub struct MessageProgressRange;

impl MessageProgressRange {
    /// OCCT Message_ProgressRange::UserBreak() — the no-abort bridge: the
    /// fix loops of `Perform` always run to completion.
    pub fn user_break(&self) -> bool {
        false
    }
}

/// OCCT `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` — the rcad indexed map keyed by (TShape
/// pointer, location index).  The TopTools_ShapeMapHasher IsSame identity
/// (TShape + Location) maps to the key pair; insertion order is kept by the
/// IndexMap.  Architecture note: OCCT compares the location value; rcad keys
/// by the location-table index, stable within one BRep.
pub(crate) type IndexedDataMapOfShapeListOfShape = IndexMap<(u64, u32), (Shape, Vec<Shape>)>;

/// Map key of a shape (TopTools_ShapeMapHasher: TShape pointer + location).
#[inline]
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsSame — same TShape with the same Location value
/// (rcad: the location-table index; see the map type note above).
#[inline]
pub(crate) fn occt_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-120) — re-host over
/// the shared explorer (architecture support; the same pattern as
/// feat/loc_ope_glued_shape.rs L67-81).
pub(crate) fn map_shapes_and_ancestors(
    brep: &mut BRep,
    the_s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexedDataMapOfShapeListOfShape,
) {
    for a_anc in topexp_explorer(brep, the_s, ta) {
        for a_exs in topexp_explorer(brep, &a_anc, ts) {
            let key = shape_key(&a_exs);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    for a_ex in topexp_explorer(brep, the_s, ts) {
        if topexp_explorer(brep, the_s, ta).is_empty() {
            let key = shape_key(&a_ex);
            m.entry(key).or_insert((a_ex, Vec::new()));
        }
    }
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) (TopExp.cxx L30-78) — the
/// no-CumOri form over the TopoDS_Iterator re-host (architecture re-host).
pub(crate) fn top_exp_vertices(brep: &mut BRep, e: &Shape) -> (Shape, Shape) {
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();
    for it in iter_subshapes(brep, e, false, true) {
        match it.orientation {
            Orientation::Forward => vfirst = it,
            Orientation::Reversed => vlast = it,
            _ => {}
        }
    }
    (vfirst, vlast)
}

// ===========================================================================
// OCCT ShapeFix.cxx L70-276: ShapeFix::SameParameter
// ===========================================================================

/// OCCT ShapeFix::SameParameter (ShapeFix.cxx L70-276).  Runs the per-edge
/// SameParameter fix (`ShapeFix_Edge::FixSameParameter`, the W1-6 GAP
/// carrier for the W3 class) over every edge of the shape, then the ":i2"
/// plane-edge tolerance update pass (L186-258).  `the_progress` is the
/// rcad no-abort bridge; `the_msg_reg` is None for the OCCT null handle.
pub fn same_parameter(
    brep: &mut BRep,
    shape: &Shape,
    enforce: bool,
    preci: f64,
    _the_progress: MessageProgressRange,
    mut the_msg_reg: Option<&mut dyn BasicMsgRegistrator>,
) -> bool {
    // L77-81: Calculate number of edges.
    let a_nb_edges = topexp_explorer(brep, shape, ShapeType::Edge).len();
    let _ = a_nb_edges;

    // L84-88: Calculate number of faces.
    let a_nb_faces = topexp_explorer(brep, shape, ShapeType::Face).len();
    let _ = a_nb_faces;

    // L90-92: aMapEF — edge -> faces ancestors map.
    let mut a_map_ef = IndexedDataMapOfShapeListOfShape::new();
    map_shapes_and_ancestors(brep, shape, ShapeType::Edge, ShapeType::Face, &mut a_map_ef);

    let mut builder = BRepBuilder::new();
    let mut status = true;
    let mut tol = preci;
    let iatol = tol > 0.0;
    // OCCT L101: Handle(ShapeFix_Edge) sfe = new ShapeFix_Edge — the W3
    // tranche 1 1:1 class (the W1-6 GAP carrier is retired).
    let mut sfe = ShapeFixEdge::new();

    // L102-103: the edge explorer and the done-message.
    let ex = topexp_explorer(brep, shape, ShapeType::Edge);
    let mut ex_pos = 0usize;
    let done_msg = MessageMsg::from_key("FixEdge.SameParameter.MSG0");

    // L106: the outer progress scope (2 steps).
    let a_ps_for_same_param = MessageProgressScope;

    {
        // L110: "Fixing edge" scope with aNbEdges steps.
        let a_ps = MessageProgressScope;

        while ex_pos < ex.len() {
            let mut ierr = 0;
            // L121-123: E = TopoDS::Edge(ex.Current()); ex.Next().
            let e = ex[ex_pos].clone();
            ex_pos += 1;

            // L125-128.
            if !iatol {
                tol = brep.tolerance(&e);
            }
            // L129-133.
            if enforce {
                builder.set_edge_same_range(brep, e.clone(), false);
                builder.set_edge_same_parameter(brep, e.clone(), false);
            }

            // L135-149: the faces sharing the edge; FixSameParameter per
            // face, or the face-less form (K2-SEP97).
            let a_list_of_faces = a_map_ef
                .get(&shape_key(&e))
                .map(|(_, l)| l.clone())
                .unwrap_or_default();
            if !a_list_of_faces.is_empty() {
                for a_f in &a_list_of_faces {
                    // L143: sfe->FixSameParameter(E, F) — the tolerance
                    // default is 0.0.
                    sfe.fix_same_parameter_face(brep, &e, a_f, 0.0);
                }
            } else {
                // L148: sfe->FixSameParameter(E).
                sfe.fix_same_parameter(brep, &e, 0.0);
            }

            // L151-157: BRep_Tool::SameParameter(E).
            let same_parameter_flag = match e.data.as_ref() {
                TShape::Edge(ed) => ed.same_parameter,
                _ => true,
            };
            if !same_parameter_flag {
                ierr = 1;
            }

            // L159-168.
            if ierr != 0 {
                status = false;
                builder.set_edge_same_range(brep, e.clone(), false);
                builder.set_edge_same_parameter(brep, e.clone(), false);
            } else if let Some(reg) = the_msg_reg.as_deref_mut() {
                if !sfe.status(ShapeExtendStatus::Ok) {
                    reg.send_shape(
                        &e,
                        &done_msg,
                        rcad_kernel::core::message::MessageGravity::Alert,
                    );
                }
            }

            // L171: aPS.Next().
            a_ps.next();
        } // -- end while

        // L174-178: the no-abort bridge never halts.
        let _ = a_ps_for_same_param;
    }

    {
        // L184: "Update tolerances" scope with aNbFaces steps.
        let a_ps = MessageProgressScope;

        // L186-188: :i2 abv 21 Aug 98: ProSTEP TR8 Motor.rle face 710 —
        // Update tolerance of edges on planes (no pcurves are stored).
        for exp in topexp_explorer(brep, shape, ShapeType::Face) {
            // L190-191: the face and its location-applied surface.
            let face = exp.clone();
            let surf = brep.face_surface_world(&face);

            // L193-206: the plane (one Geom_RectangularTrimmedSurface level
            // unwrapped to its basis).
            let plane = surf.and_then(|s| match s {
                rcad_kernel::geom::Surface3::Plane(p) => Some(p),
                rcad_kernel::geom::Surface3::Trimmed(ts) => match ts.basis.as_ref() {
                    rcad_kernel::geom::Surface3::Plane(p) => Some(p.clone()),
                    _ => None,
                },
                _ => None,
            });
            let Some(plane) = plane else {
                continue;
            };

            // L209-211: the edges of the face.
            for ed in topexp_explorer(brep, &face, ShapeType::Edge) {
                let edge = ed.clone();
                // L213: the 3D curve (BRep_Tool::Curve applies the location).
                let Some((crv, range)) = brep.edge_curve_world(&edge) else {
                    continue;
                };
                let (f, l) = (range[0], range[1]);

                // L219-223: the pcurve on the face.
                let Some((c2d, _f, _l)) = brep.curve_on_surface(&edge, &face) else {
                    continue;
                };

                // L224-225: Geom2dAdaptor_Curve + Adaptor3d_CurveOnSurface
                // over the plane (the architecture bridge — module doc).
                let tol0 = brep.tolerance(&edge);
                tol = tol0;
                let mut tol2 = tol * tol;
                let mut is_changed = false;
                const NCONTROL: i32 = 23;
                for i in 0..NCONTROL {
                    let par = (f * (NCONTROL as f64 - 1.0 - i as f64) + l * i as f64)
                        / (NCONTROL as f64 - 1.0);
                    // L235-236: pnt = crv->Value(par); prj = ACS.Value(par).
                    let pnt = crv.point_at(par);
                    let a_uv = c2d.point_at(par);
                    let prj = plane.point_at(a_uv.x, a_uv.y);
                    // L237-242.
                    let dist = pnt.distance_squared(prj);
                    if tol2 < dist {
                        tol2 = dist;
                        is_changed = true;
                    }
                }
                // L244-256.
                if is_changed {
                    // L246: coeff: see trj3_pm1-ct-203.stp #19681, edge 10.
                    tol = 1.00005 * tol2.sqrt();
                    if tol >= tol0 {
                        builder.update_edge_tolerance(brep, edge.clone(), tol);
                        for sv in iter_subshapes(brep, &edge, true, true) {
                            builder.update_vertex_tolerance(brep, sv, tol);
                        }
                    }
                }
            }
            // L188: exp.Next(), aPS.Next().
            a_ps.next();
        }
        // L259-263: the no-abort bridge never halts.
    }

    // L266-275: the OCCT_DEBUG report is compiled out; return the status.
    status
}

// ===========================================================================
// OCCT ShapeFix.cxx L280-283: ShapeFix::EncodeRegularity
// ===========================================================================

/// OCCT ShapeFix::EncodeRegularity (ShapeFix.cxx L280-283): runs
/// BRepLib::EncodeRegularity(shape, tolang).  GAP carrier: the kernel
/// BRepLib::EncodeRegularity port belongs to the docket section 4 gap 3
/// family (kernel completion item); the carrier keeps the call shape and
/// performs no update — replaced when the BRepLib batch lands.
pub fn encode_regularity(_brep: &mut BRep, shape: &Shape, tolang: f64) {
    let _ = (shape, tolang);
    // OCCT L282: BRepLib::EncodeRegularity(shape, tolang);
}

// ===========================================================================
// OCCT ShapeFix.cxx L287-309: ShapeFix::RemoveSmallEdges
// ===========================================================================

/// OCCT ShapeFix::RemoveSmallEdges (ShapeFix.cxx L287-309): removes edges
/// smaller than the tolerance through ShapeFix_Shape (the W3 tranche 4 1:1
/// class — the former W1-6 GAP carrier is retired, Rule 4).
pub fn remove_small_edges(
    brep: &mut BRep,
    shape: &mut Shape,
    tolerance: f64,
    context: &mut Option<ShapeBuildReShape>,
) -> Shape {
    // L291: occ::handle<ShapeFix_Shape> sfs = new ShapeFix_Shape.
    let mut sfs = ShapeFixShape::new();
    // L292.
    sfs.init(shape);
    // L293.
    sfs.set_precision(tolerance);
    // L294-296: the bool false assignments to the int& face modes.
    *sfs.fix_face_tool().fix_missing_seam_mode() = 0;
    *sfs.fix_face_tool().fix_orientation_mode() = 0;
    *sfs.fix_face_tool().fix_small_area_wire_mode() = 0;
    // L297 (L298 is the commented-out FixReorderMode line): ModifyTopologyMode
    // is a bool& mode.
    *sfs.fix_wire_tool().modify_topology_mode() = true;
    // L299-304: the bool assignments to the int& wire modes.
    *sfs.fix_wire_tool().fix_connected_mode() = 0;
    *sfs.fix_wire_tool().fix_edge_curves_mode() = 0;
    *sfs.fix_wire_tool().fix_degenerated_mode() = 0;
    *sfs.fix_wire_tool().fix_self_intersection_mode() = 0;
    *sfs.fix_wire_tool().fix_lacking_mode() = 0;
    *sfs.fix_wire_tool().fix_small_mode() = 1;
    // L305: sfs->Perform() — the OCCT no-arg default-progress form.
    sfs.perform(brep, MessageProgressRange::default());
    // L306: TopoDS_Shape result = sfs->Shape().
    let result = sfs.shape_result();
    // L307: context = sfs->Context().
    *context = sfs.base.my_context.clone();
    // L308: return result.
    result
}

// ===========================================================================
// OCCT ShapeFix.cxx L311-338: ReplaceVertex (file static)
// ===========================================================================

/// OCCT static ReplaceVertex (ShapeFix.cxx L315-338): auxiliary for
/// FixVertexPosition — rebuilds the edge on a fresh vertex at `the_p`.
fn replace_vertex(brep: &mut BRep, the_edge: &Shape, the_p: DVec3, the_fwd: bool) -> Shape {
    // L317-319: aNewVertex = BRep_Builder.MakeVertex(theP, Confusion()).
    let mut a_b = BRepBuilder::new();
    let a_new_vertex = a_b.add_vertex(brep, the_p, CONFUSION);
    // L320-330: aV1/aV2 by theFwd.
    let (a_v1, a_v2) = if the_fwd {
        let mut v = a_new_vertex.clone();
        v.orientation = Orientation::Forward;
        (Some(v), None)
    } else {
        let mut v = a_new_vertex;
        v.orientation = Orientation::Reversed;
        (None, Some(v))
    };
    // L331-337: ShapeBuild_Edge::CopyReplaceVertices over a FORWARD copy.
    // (rcad passes the null shape for the OCCT empty-handle vertex.)
    let a_sbe = ShapeBuildEdge;
    let mut e1 = the_edge.clone();
    let ori = e1.orientation;
    e1.orientation = Orientation::Forward;
    let a_new_edge = a_sbe.copy_replace_vertices(
        brep,
        &e1,
        a_v1.as_ref().unwrap_or(&Shape::null()),
        a_v2.as_ref().unwrap_or(&Shape::null()),
    );
    let mut a_new_edge = a_new_edge;
    a_new_edge.orientation = ori;
    a_new_edge
}

// ===========================================================================
// OCCT ShapeFix.cxx L340-376: getNearPoint (file static)
// ===========================================================================

/// OCCT static getNearPoint (ShapeFix.cxx L344-376): auxiliary for
/// FixVertexPosition — the closest pair between two point sequences.
fn get_near_point(a_seq1: &[DVec3], a_seq2: &[DVec3], acent: &mut DVec3) -> f64 {
    let mut ind1 = 0usize;
    let mut ind2 = 0usize;
    let mut mindist = f64::MAX;
    for i in 1..=a_seq1.len() {
        let p1 = a_seq1[i - 1];
        for j in 1..=a_seq2.len() {
            let p2 = a_seq2[j - 1];
            let d = p1.distance(p2);
            if (d - mindist).abs() <= CONFUSION {
                continue;
            }
            if d < mindist {
                mindist = d;
                ind1 = i;
                ind2 = j;
            }
        }
    }
    if ind1 != 0 && ind2 != 0 {
        *acent = (a_seq1[ind1 - 1] + a_seq2[ind2 - 1]) / 2.0;
    }
    mindist
}

// ===========================================================================
// OCCT ShapeFix.cxx L378-561: getNearestEdges (file static)
// ===========================================================================

/// OCCT static getNearestEdges (ShapeFix.cxx L382-561): auxiliary for
/// FixVertexPosition — splits the vertex edge list into suitable and
/// rejected edges around the vertex position.
#[allow(clippy::too_many_arguments)]
fn get_nearest_edges(
    brep: &mut BRep,
    the_ledges: &mut Vec<Shape>,
    the_vert: &Shape,
    the_suit_edges: &mut Vec<Shape>,
    the_reject_edges: &mut Vec<Shape>,
    the_tolerance: f64,
    thecentersuit: &mut DVec3,
    thecenterreject: &mut DVec3,
) -> bool {
    if the_ledges.is_empty() {
        return false;
    }
    // L394: aMapEdges — the processed-edge map.
    let mut a_map_edges: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();

    // L396-409: the first edge, its extremity vertices and curve.
    let mut atemp_list: Vec<Shape> = the_ledges.clone();
    let a_edge1 = atemp_list[0].clone();
    let (a_vert11, a_vert12) = top_exp_vertices(brep, &a_edge1);
    a_map_edges.insert(shape_key(&a_edge1));
    // L404-409: aCurve1 = BRep_Tool::Curve(aEdge1, aFirst1, aLast1); the
    // null-curve path returns false (L425-428).
    let Some((a_curve1, (a_first1, a_last1))) = brep
        .edge_curve_world(&a_edge1)
        .map(|(c, r)| (c, (r[0], r[1])))
    else {
        return false;
    };
    let mut p11 = DVec3::ZERO;
    let mut p12 = DVec3::ZERO;
    let is_first1 = occt_is_same(the_vert, &a_vert11);
    let is_same1 = occt_is_same(&a_vert11, &a_vert12);
    {
        // L410-424.
        if is_first1 {
            p11 = a_curve1.point_at(a_first1);
        } else if !is_same1 {
            p11 = a_curve1.point_at(a_last1);
        }
        if is_same1 {
            p12 = a_curve1.point_at(a_last1);
        }
    }

    // L429-433: the working lists.
    let mut aseqreject: Vec<Shape> = Vec::new();
    let mut aseqsuit: Vec<Shape> = Vec::new();

    let mut anum_loop = 0usize;
    let mut pos = 1usize;
    while pos < atemp_list.len() {
        let a_edge = atemp_list[pos].clone();
        if a_map_edges.contains(&shape_key(&a_edge)) {
            // L437-441: atempList.Remove(alIter); continue.
            atemp_list.remove(pos);
            continue;
        }

        // L443-446.
        let (a_vert1, a_vert2) = top_exp_vertices(brep, &a_edge);
        let is_first = occt_is_same(the_vert, &a_vert1);
        let is_same = occt_is_same(&a_vert1, &a_vert2);

        // L448-456: the closed-loop candidate handling.
        let is_loop = (occt_is_same(&a_vert1, &a_vert11) && occt_is_same(&a_vert2, &a_vert12))
            || (occt_is_same(&a_vert1, &a_vert12) && occt_is_same(&a_vert2, &a_vert11));
        if is_loop && atemp_list.len() > anum_loop {
            atemp_list.push(a_edge);
            atemp_list.remove(pos);
            anum_loop += 1;
            continue;
        }
        a_map_edges.insert(shape_key(&a_edge));

        // L458-460: aCurve = BRep_Tool::Curve(aEdge, aFirst, aLast).
        let (a_first, a_last, a_curve) = match brep.edge_curve_world(&a_edge) {
            Some((c, r)) => (r[0], r[1], Some(c)),
            None => (0.0, 0.0, None),
        };
        if let Some(a_curve) = a_curve {
            let mut p1;
            let mut p2 = DVec3::ZERO;
            if is_first {
                p1 = a_curve.point_at(a_first);
            } else {
                p1 = a_curve.point_at(a_last);
            }
            if is_same {
                p2 = a_curve.point_at(a_last);
            }
            let a_min_dist;
            let mut acent = DVec3::ZERO;
            if !is_same && !is_same1 {
                // L478-482.
                a_min_dist = p1.distance(p11);
                acent = (p1 + p11) / 2.0;
            } else {
                // L483-498.
                let mut a_seq1: Vec<DVec3> = Vec::new();
                let mut a_seq2: Vec<DVec3> = Vec::new();
                a_seq1.push(p11);
                if is_same1 {
                    a_seq1.push(p12);
                }
                a_seq2.push(p1);
                if is_same {
                    a_seq2.push(p2);
                }
                a_min_dist = get_near_point(&a_seq1, &a_seq2, &mut acent);
            }

            // L500-527.
            if a_min_dist > the_tolerance {
                if aseqreject.is_empty() {
                    *thecenterreject = acent;
                }
                aseqreject.push(a_edge);
            } else {
                if aseqsuit.is_empty() {
                    *thecentersuit = acent;
                    aseqsuit.push(a_edge);
                } else if !is_same1 {
                    aseqsuit.push(a_edge);
                } else if (*thecentersuit - acent).length() < the_tolerance {
                    aseqsuit.push(a_edge);
                } else {
                    aseqreject.push(a_edge);
                }
            }
        }
        // L529: atempList.Remove(alIter).
        atemp_list.remove(pos);
    }

    // L532-560: the recursive consume of the remaining list.
    let is_done = !aseqsuit.is_empty() || !aseqreject.is_empty();
    if is_done {
        if aseqsuit.is_empty() {
            the_reject_edges.push(a_edge1);
            the_ledges.remove(0);

            get_nearest_edges(
                brep,
                the_ledges,
                the_vert,
                the_suit_edges,
                the_reject_edges,
                the_tolerance,
                thecentersuit,
                thecenterreject,
            );
        } else {
            the_suit_edges.push(a_edge1);
            the_suit_edges.extend(aseqsuit);
            the_reject_edges.extend(aseqreject);
        }
    } else {
        the_reject_edges.push(a_edge1);
    }

    is_done
}

// ===========================================================================
// OCCT ShapeFix.cxx L563-789: ShapeFix::FixVertexPosition
// ===========================================================================

/// OCCT ShapeFix::FixVertexPosition (ShapeFix.cxx L565-789).
pub fn fix_vertex_position(
    brep: &mut BRep,
    theshape: &mut Shape,
    the_tolerance: f64,
    thecontext: &mut ShapeBuildReShape,
) -> bool {
    // L569-599: the vertex -> edges map.
    let mut a_map_vert_edges: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
    for a_exp1 in topexp_explorer(brep, theshape, ShapeType::Edge) {
        let mut n_v = 0;
        let a_exp3 = iter_subshapes(brep, &a_exp1, true, true);
        let mut a_vert1 = Shape::null();
        for a_vert in a_exp3 {
            n_v += 1;
            if n_v == 1 {
                a_vert1 = a_vert.clone();
            } else if occt_is_same(&a_vert1, &a_vert) {
                continue;
            }
            let key = shape_key(&a_vert);
            match a_map_vert_edges.get_mut(&key) {
                Some((_, list)) => list.push(a_exp1.clone()),
                None => {
                    a_map_vert_edges.insert(key, (a_vert, vec![a_exp1.clone()]));
                }
            }
        }
    }

    // L600-783: per vertex above the tolerance.
    let mut is_done = false;
    for i in 1..=a_map_vert_edges.len() {
        let (a_vert, a_ledges_all) = {
            let (_, (s, l)) = a_map_vert_edges.get_index(i - 1).unwrap();
            (s.clone(), l.clone())
        };
        let a_tol_vert = brep.tolerance(&a_vert);
        if a_tol_vert <= the_tolerance {
            continue;
        }

        // L611-614.
        let mut a_b1 = BRepBuilder::new();
        a_b1.update_vertex_tolerance(brep, a_vert.clone(), the_tolerance);
        let a_pvert = brep.vertex_position(&a_vert);
        let mut acenter = a_pvert;
        let mut acenterreject = a_pvert;

        // L616-623.
        let mut a_suit_edges: Vec<Shape> = Vec::new();
        let mut a_reject_edges: Vec<Shape> = Vec::new();
        let mut aledges = a_ledges_all.clone();
        if aledges.len() == 1 {
            continue;
        }
        // L625-636: the curve-distance check around the vertex.
        if !get_nearest_edges(
            brep,
            &mut aledges,
            &a_vert,
            &mut a_suit_edges,
            &mut a_reject_edges,
            the_tolerance,
            &mut acenter,
            &mut acenterreject,
        ) {
            continue;
        }

        // L639-719: update the vertex by the nearest point.
        let mut is_add = false;
        for k in 1..=a_suit_edges.len() {
            let a_edge_old = a_suit_edges[k - 1].clone();
            let (a_vert1, a_vert2) = top_exp_vertices(brep, &a_edge_old);

            let is_first = occt_is_same(&a_vert1, &a_vert);
            let is_last = occt_is_same(&a_vert2, &a_vert);
            if !is_first && !is_last {
                continue;
            }
            // L656: aEdge = TopoDS::Edge(thecontext->Apply(aEdgeOld)).
            let a_edge = thecontext.apply(brep, &a_edge_old, ShapeType::Shape);

            let (a_vert1n, a_vert2n) = top_exp_vertices(brep, &a_edge);
            let Some((a_curve, (a_first, a_last))) = brep
                .edge_curve_world(&a_edge)
                .map(|(c, r)| (c, (r[0], r[1])))
            else {
                continue;
            };
            {
                let p1 = a_curve.point_at(a_first);
                let p2 = a_curve.point_at(a_last);

                // L666-669: the same-vertex distance check.
                let is_replace =
                    occt_is_same(&a_vert1n, &a_vert2n) && p1.distance(p2) > the_tolerance;

                if is_first {
                    if k > 2 {
                        acenter += p1;
                        acenter /= 2.0;
                    }
                    if is_replace {
                        let enew = if p1.distance(acenter) < p2.distance(acenter) {
                            replace_vertex(brep, &a_edge, p2, false)
                        } else {
                            replace_vertex(brep, &a_edge, p1, true)
                        };
                        thecontext.replace(brep, &a_edge, &enew);
                        is_done = true;
                    }
                } else {
                    if k > 2 {
                        acenter += p2;
                        acenter /= 2.0;
                    }
                    if is_replace {
                        let enew = if p1.distance(acenter) < p2.distance(acenter) {
                            replace_vertex(brep, &a_edge, p2, false)
                        } else {
                            replace_vertex(brep, &a_edge, p1, true)
                        };
                        thecontext.replace(brep, &a_edge, &enew);
                        is_done = true;
                    }
                }

                is_add = true;
            }
        }

        // L722-734.
        if is_add && a_pvert.distance(acenter) > the_tolerance {
            let mut a_b = BRepBuilder::new();
            is_done = true;
            let mut a_new_vertex = a_b.add_vertex(brep, acenter, CONFUSION);
            a_new_vertex.orientation = a_vert.orientation;
            thecontext.replace(brep, &a_vert, &a_new_vertex);
        }

        // L736-782: the rejected edges — their extremities are replaced by
        // fresh vertices at the curve endpoints.
        for k in 1..=a_reject_edges.len() {
            let a_edge_old = a_reject_edges[k - 1].clone();
            let (a_vert1, a_vert2) = top_exp_vertices(brep, &a_edge_old);

            let is_first = occt_is_same(&a_vert1, &a_vert);
            let is_last = occt_is_same(&a_vert2, &a_vert);
            if !is_first && !is_last {
                continue;
            }
            let is_same = occt_is_same(&a_vert1, &a_vert2);
            let a_edge = thecontext.apply(brep, &a_edge_old, ShapeType::Shape);

            let Some((a_curve, (a_first, a_last))) = brep
                .edge_curve_world(&a_edge)
                .map(|(c, r)| (c, (r[0], r[1])))
            else {
                continue;
            };
            let p1 = a_curve.point_at(a_first);
            let p2 = a_curve.point_at(a_last);
            let mut enew = if is_first {
                let e = replace_vertex(brep, &a_edge, p1, true);
                if is_same {
                    replace_vertex(brep, &e, p2, false)
                } else {
                    e
                }
            } else {
                let e = replace_vertex(brep, &a_edge, p2, false);
                if is_same {
                    replace_vertex(brep, &e, p1, true)
                } else {
                    e
                }
            };
            enew.orientation = a_edge.orientation;
            thecontext.replace(brep, &a_edge, &enew);
            is_done = true;
        }
    }

    // L784-788.
    if is_done {
        *theshape = thecontext.apply(brep, theshape, ShapeType::Shape);
    }
    is_done
}

// ===========================================================================
// OCCT ShapeFix.cxx L791-818: ShapeFix::LeastEdgeSize
// ===========================================================================

/// OCCT ShapeFix::LeastEdgeSize (ShapeFix.cxx L793-818): the size of the
/// least edge (the squared component sum, square-rooted as in OCCT L816).
pub fn least_edge_size(brep: &mut BRep, the_shape: &Shape) -> f64 {
    let mut a_res = f64::MAX;
    for exp in topexp_explorer(brep, the_shape, ShapeType::Edge) {
        let edge = exp.clone();
        let Some((c3d, (first, last))) =
            brep.edge_curve_world(&edge).map(|(c, r)| (c, (r[0], r[1])))
        else {
            continue;
        };
        // L803-806: the three-point Bnd_Box.
        let bb = [
            c3d.point_at(first),
            c3d.point_at(last),
            c3d.point_at((last + first) / 2.0),
        ];
        let x1 = bb.iter().map(|p| p.x).fold(f64::MAX, f64::min);
        let x2 = bb.iter().map(|p| p.x).fold(f64::MIN, f64::max);
        let y1 = bb.iter().map(|p| p.y).fold(f64::MAX, f64::min);
        let y2 = bb.iter().map(|p| p.y).fold(f64::MIN, f64::max);
        let z1 = bb.iter().map(|p| p.z).fold(f64::MAX, f64::min);
        let z2 = bb.iter().map(|p| p.z).fold(f64::MIN, f64::max);
        // L807-813.
        let size = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1);
        if size < a_res {
            a_res = size;
        }
    }
    // L816: aRes = sqrt(aRes).
    a_res.sqrt()
}
