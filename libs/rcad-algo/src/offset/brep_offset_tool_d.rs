// OCCT BRepOffset_Tool.cxx L4422-L4659 — module d of the 1:1 translation
// (the ExtentFace tail: IsInOut / CorrectOrientation / CheckPlanesNormals /
// PerformPlanes / UpdateVertexTolerances; split from brep_offset_tool_c.rs
// to respect the 2000-line file limit; see the architecture-difference list
// in brep_offset_tool.rs — numbering continues there).
//
// This module also carries the BRepOffset_Analyse / BRepOffset_Interval GAP
// carriers shared with the Inter3d module (the OCCT dependency direction
// Inter3d -> Tool is kept).

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range, brep_tool_tolerance,
};

use super::brep_offset_tool::*;

/// OCCT BRepTopAdaptor_FClass2d(Face, Tol) — the rcad carrier over the
/// topalgo FClass2d + FaceShapeSource (architecture difference #29; the
/// loc_ope_wires_on_shape_b.rs #8 precedent).  The classifier is rebuilt
/// per Perform call (the OCCT class caches the UV polygons; the carrier
/// form keeps the OCCT call structure without the self-referential borrow).
pub(crate) struct FClass2dCarrier {
    face: Shape,
    surf: Surface3,
}

impl FClass2dCarrier {
    /// OCCT BRepTopAdaptor_FClass2d(F, Tol).
    pub fn new(face: &Shape) -> Self {
        let surf = face_surface_of(face)
            .expect("BRepTopAdaptor_FClass2d: the face carries no surface");
        FClass2dCarrier { face: face.clone(), surf }
    }

    /// OCCT BRepTopAdaptor_FClass2d::Perform(Puv).
    pub fn perform(&self, p_uv: DVec2) -> State {
        let locations = [glam::DAffine3::IDENTITY];
        let fss = crate::topalgo::shape_source::FaceShapeSource::new(
            &self.face,
            self.surf.clone(),
            &locations,
        );
        let fc = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
            &fss,
            0,
            rcad_kernel::precision::CONFUSION,
        );
        fc.perform(&fss, p_uv, false)
    }
}



// ---------------------------------------------------------------------------
// OCCT static IsInOut (cxx L4422-4438).
// ---------------------------------------------------------------------------

/// OCCT static IsInOut(FC, AC, S) (cxx L4422-4438) — the pcurve samples are
/// all classified S.  The GCPnts_QuasiUniformDeflection sampling is the
/// reduced uniform form (architecture difference #29); the classifier is
/// the topalgo FClass2d over FaceShapeSource (the loc_ope #8 precedent).
pub(crate) fn is_in_out(
    fc: &FClass2dCarrier,
    ac: &Curve2d,
    f: f64,
    l: f64,
    state: State,
) -> bool {
    // OCCT L4426-4427: Def = 100 * Confusion; QU(AC, Def).
    let n = 23usize;
    for i in 1..=n {
        let t = f + (l - f) * (i as f64) / (n as f64);
        let p = Curve2dEval::point_at(ac, t);
        // OCCT L4431-4434: if (FC.Perform(P) != S) return false.
        if fc.perform(p) != state {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::CorrectOrientation (hxx L195-200; cxx L4442-4505).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::CorrectOrientation(SI, NewEdges, AsDes,
/// InitOffset, Offset) (cxx L4442-4505).
pub fn correct_orientation(
    si: &Shape,
    new_edges: &OcctIndexedShapeMap,
    as_des: &mut crate::brep_algo::as_des::BRepAlgoAsDes,
    init_offset: &mut crate::brep_algo::image::BRepAlgoImage,
    offset: f64,
) {
    // OCCT L4450-4452: exp.Init(SI, FACE); f = 0., l = 0.
    let mut f_par = 0.0f64;
    let mut l_par = 0.0f64;

    for fi in explorer(si, ShapeType::Face, ShapeType::Shape) {
        // OCCT L4457-4458: LOF = InitOffset.Image(FI).
        let lof = init_offset.image(&fi);
        for of in &lof {
            // OCCT L4463: LOE = AsDes->ChangeDescendant(OF).
            let has_loe = as_des.has_descendant(of);
            if !has_loe {
                continue;
            }
            // OCCT L4465-4475: the YaInt probe.
            let loe_snapshot = as_des.descendant(of).to_vec();
            let mut ya_int = false;
            for oe in &loe_snapshot {
                if new_edges.contains(oe) {
                    ya_int = true;
                    break;
                }
            }
            if ya_int {
                // OCCT L4478-4479: FC = BRepTopAdaptor_FClass2d(the
                // FORWARD FI, Confusion).
                let mut fi_fwd = fi.clone();
                fi_fwd.orientation = Orientation::Forward;
                let fc = FClass2dCarrier::new(&fi_fwd);

                // OCCT L4482-4501: the orientation walk.
                let loe_snapshot = as_des.descendant(of).to_vec();
                for (pos, oe) in loe_snapshot.iter().enumerate() {
                    if new_edges.contains(oe) {
                        // OCCT L4487-4488: CO2d = CurveOnSurface(OE, OF);
                        // AC(CO2d, f, l).
                        let (co2d, cf, cl) = match brep_tool_curve_on_surface(oe, of) {
                            Some(t) => t,
                            None => continue,
                        };
                        f_par = cf;
                        l_par = cl;

                        if offset > 0.0 {
                            // OCCT L4492-4494: if (IsInOut(FC, AC,
                            // TopAbs_OUT)) OE.Reverse().
                            if is_in_out(&fc, &co2d, f_par, l_par, State::Out) {
                                // OCCT L4484: TopoDS_Shape& OE =
                                // itE.ChangeValue() — the in-place
                                // orientation flip of the descendant entry.
                                if let Some(loe) = as_des.change_descendant(of) {
                                    if let Some(oe_mut) = loe.get_mut(pos) {
                                        oe_mut.orientation = bat::top_abs_reverse(oe_mut.orientation);
                                    }
                                }
                            }
                        }
                        // OCCT L4497-4499: the commented-out else branch.
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::CheckPlanesNormals (hxx L206-208; cxx L4509-4533).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::CheckPlanesNormals(theFace1, theFace2, theTolAng)
/// (cxx L4509-4533) — compares the normal directions of the planar faces.
pub fn check_planes_normals(the_face1: &Shape, the_face2: &Shape, the_tol_ang: f64) -> bool {
    // OCCT L4513-4517: the planar probes.
    let s1 = face_surface_of(the_face1);
    let s2 = face_surface_of(the_face2);
    let (p1, p2) = match (&s1, &s2) {
        (Some(Surface3::Plane(a)), Some(Surface3::Plane(b))) => (a, b),
        _ => return false,
    };
    //
    // OCCT L4519-4523: aDN1 = the plane direction (reversed for REVERSED).
    let mut a_dn1 = p1.normal;
    if the_face1.orientation == Orientation::Reversed {
        a_dn1 = -a_dn1;
    }
    //
    // OCCT L4525-4529: aDN2 (idem).
    let mut a_dn2 = p2.normal;
    if the_face2.orientation == Orientation::Reversed {
        a_dn2 = -a_dn2;
    }
    //
    // OCCT L4531-4532: anAngle = aDN1.Angle(aDN2).
    let an_angle = a_dn1.angle_between(a_dn2);
    an_angle < the_tol_ang
}

// ---------------------------------------------------------------------------
// OCCT static PerformPlanes (cxx L4537-4602).
// ---------------------------------------------------------------------------

/// OCCT static PerformPlanes(theFace1, theFace2, theSide, theL1, theL2)
/// (cxx L4537-4602) — the IntTools_FaceFace plane/plane intersection used by
/// the planar fast path of BRepOffset_Tool::Inter3D (cxx L1453-1467).  The
/// bare-face IntTools_FaceFace::Perform form lives in the bop IntTools FaceFace
/// (crate::bop::int_tools::face_face::perform_face_face_planes); this leaf
/// builds the section edge from its first curve.
pub(crate) fn perform_planes(
    the_face1: &Shape,
    the_face2: &Shape,
    the_side: State,
    the_l1: &mut Vec<Shape>,
    the_l2: &mut Vec<Shape>,
) {
    // OCCT L4543-4544: theL1.Clear(); theL2.Clear().
    the_l1.clear();
    the_l2.clear();
    // Intersect the planes using IntTools_FaceFace directly
    // OCCT L4545-4548: IntTools_FaceFace aFF;
    //   aFF.SetParameters(true, true, true, Precision::Confusion());
    //   aFF.Perform(theFace1, theFace2).
    let (is_done, a_sc) = crate::bop::int_tools::face_face::perform_face_face_planes(
        the_face1,
        the_face2,
        rcad_kernel::precision::CONFUSION,
    );
    //
    // OCCT L4550-4553: if (!aFF.IsDone()) return.
    if !is_done {
        return;
    }
    // OCCT L4555-4559: const Sequence<IntTools_Curve>& aSC = aFF.Lines();
    //   if (aSC.IsEmpty()) return.
    if a_sc.is_empty() {
        return;
    }
    //
    // In Plane/Plane intersection only one curve is always produced.
    // Make the edge from this section curve.
    // OCCT L4563-4586.
    let a_e;
    let a_tf;
    let a_tl;
    {
        let mut a_bb_e = bat::builder_make_edge();
        // OCCT L4566-4568: const IntTools_Curve& aIC = aSC(1);
        //   const Handle(Geom_Curve)& aC3D = aIC.Curve();
        //   aBB.MakeEdge(aE, aC3D, aIC.Tolerance()).
        let a_ic = &a_sc[0];
        // OCCT BRep_Builder::MakeEdge(E, C, Tol) keeps the range of a bounded
        // curve (BRep_Builder.cxx UpdateEdge: the Curve3D range is the curve's
        // own first/last parameter).  The rcad IntersectionCurve carries the
        // trim separately (curve + t_range — architecture difference #30), so
        // the range is written here and re-affirmed at L4585.
        crate::offset::brep_offset_inter2d::builder_update_edge_curve(
            &mut a_bb_e,
            Some(a_ic.curve.clone()),
            a_ic.tolerance,
        );
        // OCCT L4569-4572: double aTF, aTL; gp_Pnt aPF, aPL;
        //   aIC.Bounds(aTF, aTL, aPF, aPL).
        a_tf = a_ic.t_range[0];
        a_tl = a_ic.t_range[1];
        let a_pf = CurveEval::point_at(&a_ic.curve, a_tf);
        let a_pl = CurveEval::point_at(&a_ic.curve, a_tl);
        bat::builder_range_edge(&mut a_bb_e, a_tf, a_tl);
        // OCCT L4574-4577: MakeVertex(aVF, aPF, aIC.Tolerance()); idem aVL;
        //   aVL.Orientation(TopAbs_REVERSED).
        let mut a_vf = bat::builder_make_vertex();
        bat::builder_update_vertex_point_tol(&mut a_vf, a_pf, a_ic.tolerance);
        let mut a_vl = bat::builder_make_vertex();
        bat::builder_update_vertex_point_tol(&mut a_vl, a_pl, a_ic.tolerance);
        a_vl.orientation = Orientation::Reversed;
        // OCCT L4579-4580: aBB.Add(aE, aVF); aBB.Add(aE, aVL).
        bat::builder_add_edge_vertex(&mut a_bb_e, &a_vf);
        bat::builder_add_edge_vertex(&mut a_bb_e, &a_vl);
        // OCCT L4582-4583: aBB.UpdateEdge(aE, aIC.FirstCurve2d(), theFace1,
        //   aIC.Tolerance()); aBB.UpdateEdge(aE, aIC.SecondCurve2d(),
        //   theFace2, aIC.Tolerance()).
        if let Some(c2d) = &a_ic.pcurve1 {
            bat::builder_update_edge_pcurve(&mut a_bb_e, c2d, the_face1, a_ic.tolerance);
        }
        if let Some(c2d) = &a_ic.pcurve2 {
            bat::builder_update_edge_pcurve(&mut a_bb_e, c2d, the_face2, a_ic.tolerance);
        }
        // OCCT L4585: aBB.Range(aE, aTF, aTL).
        bat::builder_range_edge(&mut a_bb_e, a_tf, a_tl);
        a_e = a_bb_e;
    }
    //
    // Orient section
    // OCCT L4589-4595: TopAbs_Orientation O1, O2;
    //   BRepOffset_Tool::OrientSection(aE, theFace1, theFace2, O1, O2);
    //   if (theSide == TopAbs_OUT) { O1 = TopAbs::Reverse(O1);
    //                                O2 = TopAbs::Reverse(O2); }
    let mut o1 = Orientation::Forward;
    let mut o2 = Orientation::Forward;
    orient_section(&a_e, the_face1, the_face2, &mut o1, &mut o2);
    if the_side == State::Out {
        o1 = bat::top_abs_reverse(o1);
        o2 = bat::top_abs_reverse(o2);
    }
    //
    // OCCT L4597: BRepLib::SameParameter(aE, Precision::Confusion(), true).
    crate::offset::brep_offset_inter2d::brep_lib_same_parameter(&a_e, rcad_kernel::precision::CONFUSION);
    //
    // Add edge to result
    // OCCT L4600-4601: theL1.Append(aE.Oriented(O1));
    //   theL2.Append(aE.Oriented(O2)).
    the_l1.push(bat::oriented(&a_e, o1));
    the_l2.push(bat::oriented(&a_e, o2));
}

// ---------------------------------------------------------------------------
// OCCT static UpdateVertexTolerances (cxx L4613-4659).
// ---------------------------------------------------------------------------

/// OCCT static UpdateVertexTolerances(theFace) (cxx L4613-4659).
pub(crate) fn update_vertex_tolerances(the_face: &Shape) {
    // OCCT L4615-4618: VEmap = MapShapesAndAncestors(VERTEX, EDGE).
    let mut vemak: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
        the_face,
        ShapeType::Vertex,
        ShapeType::Edge,
        &mut vemak,
    );

    for i in 1..=shape_indexed_data_map::extent(&vemak) {
        let a_vertex = shape_indexed_data_map::find_key_1(&vemak, i).clone();
        let elist = shape_indexed_data_map::value_1(&vemak, i).clone();
        let pnt_vtx = match bat::brep_tool_pnt(&a_vertex) {
            Some(t) => t,
            None => continue,
        };
        for it_value in &elist {
            let an_edge = it_value;
            let (v1, v2) = top_exp_vertices(an_edge);
            let (fpar, lpar) = brep_tool_range(an_edge);
            // OCCT L4633: aParam = (V1.IsSame(aVertex)) ? fpar : lpar.
            let a_param = if v1.is_same(&a_vertex) { fpar } else { lpar };
            if !brep_tool_degenerated(an_edge) {
                // OCCT L4636-4644: the 3d-curve distance probes.
                if let Some((c, _, _)) = brep_tool_curve(an_edge) {
                    let a_pnt = CurveEval::point_at(&c, a_param);
                    let a_dist = pnt_vtx.distance(a_pnt);
                    let mut vtx = a_vertex.clone();
                    bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    if v1.is_same(&v2) {
                        let a_pnt = CurveEval::point_at(&c, lpar);
                        let a_dist = pnt_vtx.distance(a_pnt);
                        bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    }
                }
            }
            // OCCT L4647-4656: the pcurve-on-face distance probes (the
            // BRepAdaptor_Curve(E, theFace) form) — the rcad form evaluates
            // the pcurve point through the face surface.
            if let Some((c2d, _, _)) = brep_tool_curve_on_surface(an_edge, the_face) {
                if let Some(surf) = face_surface_of(the_face) {
                    let p2d = Curve2dEval::point_at(&c2d, a_param);
                    let a_pnt = surf.point_at(p2d.x, p2d.y);
                    let a_dist = pnt_vtx.distance(a_pnt);
                    let mut vtx = a_vertex.clone();
                    bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    if v1.is_same(&v2) {
                        let p2d = Curve2dEval::point_at(&c2d, lpar);
                        let a_pnt = surf.point_at(p2d.x, p2d.y);
                        let a_dist = pnt_vtx.distance(a_pnt);
                        bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// BRepOffset_Analyse / BRepOffset_Interval GAP carriers (shared with the
// Inter3d module; the OCCT dependency direction Inter3d -> Tool keeps the
// carriers here).
// ---------------------------------------------------------------------------

// OCCT BRepOffset_Interval / BRepOffset_Analyse — the real bodies live in
// super::brep_offset_analyse (the E0 carrier-switch list; the local panic
// carriers are deleted).  Re-exported for the Inter3d / Tool consumers.
pub(crate) use super::brep_offset_analyse::{BRepOffsetAnalyse, BRepOffsetInterval};
