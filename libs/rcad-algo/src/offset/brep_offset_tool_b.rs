// OCCT BRepOffset_Tool.cxx L1441-L2576 — module b of the 1:1 translation
// (see the architecture-difference list in brep_offset_tool.rs; numbering
// continues there).
//
// Module b carries the class methods Inter3D / TryProject and the statics
// ExtentEdge / ProjectVertexOnEdge / Inter2d.
//
// D6 (the plan §0.6 key consumption): BRepOffset_Tool::Inter3D drives a
// full BOPAlgo_PaveFiller on the face pair (cxx L1475-1502) and reads the
// section edges back out of the BOPDS (InterfFF -> BOPDS_Curve ->
// PaveBlocks -> PB::Edge -> DS shape).  The rcad equivalent is the
// crate::bop::algo::pave_filler::PaveFiller pipeline: SetArguments +
// Perform, then the DS walk over `ds.interf_ff` (the BOPDS_InterfFF
// vector), `ds.intersection_curves[curves[i]]` (the BOPDS_Curve vector;
// its `pave_blocks` are the BOPDS_Curve::PaveBlocks and its `pcurve_on_a`
// / `pcurve_on_b` are the IntTools_Curve FirstCurve2d/SecondCurve2d), and
// the section edges through the DS ShapeSource (`p_ds.shape_at(pb.edge())`
// == `pDS->Shape(nSect)`).

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Plane, Surface3, SurfaceEval,
};
use rcad_kernel::topo::topods::{Orientation, State, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::tool as bat;
use crate::topalgo::shape_source::ShapeSource as _;
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_curve_on_surface;
use crate::brep_algo::tool::brep_tool_tolerance;

use super::brep_offset_tool::*;
use super::brep_offset_tool_d::{perform_planes, update_vertex_tolerances};

// ---------------------------------------------------------------------------
// Local builder re-hosts (thin BRep_Builder forms over the
// crate::brep_algo::tool layer).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::UpdateEdge(E, Tol) — the edge tolerance update
/// (BRep_TEdge::UpdateTolerance keeps the max).
fn b_update_edge_tolerance(e: &mut Shape, tol: f64) {
    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e.data) {
        ed.tolerance = ed.tolerance.max(tol);
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, P, E, Tol) — the parameter-on-edge
/// form (the rcad edge `vertex_params` map; the INTERNAL-vertex parameter
/// storage of the crate::brep_algo::tool layer).
pub(crate) fn b_update_vertex_on_edge(e_fwd: &mut Shape, v: &Shape, par: f64, tol: f64) {
    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e_fwd.data) {
        ed.vertex_params.insert(v.ptr_id(), par);
        ed.tolerance = ed.tolerance.max(tol);
    }
}

/// OCCT BRepLib::BuildCurves3d(S) (BRepLib.cxx L460-464) — no-op re-host
/// (the real body is translated in topalgo/brep_lib/build_curves3d.rs, but
/// the Inter3d pipeline carries no &mut BRep at these call sites; the
/// remaining GAP is the pool threading).
fn brep_lib_build_curves3d(_the_s: &Shape) {}

/// OCCT BRepLib::SameParameter(E, Tol, OnlySegments) — GAP no-op (the
/// brep_offset_offset_b.rs precedent).
fn brep_lib_same_parameter(_the_e: &Shape, _the_tol: f64) {}

/// OCCT BOPTools_AlgoTools2D::HasCurveOnSurface(E, F)
/// (BOPTools_AlgoTools2D.cxx L67-71) — the pcurve presence probe.
fn has_curve_on_surface(e: &Shape, f: &Shape) -> bool {
    brep_tool_curve_on_surface(e, f).is_some()
}

/// OCCT BOPTools_AlgoTools2D::AdjustPCurveOnFace(F, f, l, PC, PCNew,
/// Context) — GAP leaf: the rcad re-host (pave_filler_make_blocks.rs) is
/// DS-bound; the bare-face Tool form is staged.  The pass-through below
/// matches the OCCT behaviour on non-periodic domains (the adjustment is
/// the periodic-shift only; annotated per plan §0.6).
fn adjust_pcurve_on_face_gap(pcurve: &Curve2d) -> Curve2d {
    pcurve.clone()
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::Inter3D (hxx L84-91; cxx L1441-1941).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::Inter3D(F1, F2, L1, L2, Side, RefEdge, RefFace1,
/// RefFace2) (cxx L1441-1941) — the section between the two faces; the
/// edges of the solution are stored in L1 with the orientation on F1 and
/// in L2 with the orientation on F2.
#[allow(clippy::too_many_arguments)]
pub fn inter3d(
    f1: &Shape,
    f2: &Shape,
    l1: &mut Vec<Shape>,
    l2: &mut Vec<Shape>,
    side: State,
    ref_edge: &Shape,
    the_ref_face1: &Shape,
    the_ref_face2: &Shape,
) {
    // OCCT L1453-1467: the planar untrimmed fast path — PerformPlanes
    // without the PaveFiller.
    {
        let a_bas1 = face_surface_of(f1);
        let a_bas2 = face_surface_of(f2);
        let is_plane = |s: &Option<Surface3>| matches!(s, Some(Surface3::Plane(_)));
        if is_plane(&a_bas1) && is_plane(&a_bas2) {
            // OCCT L1456-1461: aBAS1.Initialize(F1, true) — the trimmed UV
            // domain; the IsInf probes use the trimmed (u-last, v-last)
            // parameters.  The rcad trimmed surface carries the trim box.
            let last_u_1 = trimmed_last_u(&a_bas1);
            let last_v_1 = trimmed_last_v(&a_bas1);
            if is_inf(last_u_1) && is_inf(last_v_1) {
                let last_u_2 = trimmed_last_u(&a_bas2);
                let last_v_2 = trimmed_last_v(&a_bas2);
                if is_inf(last_u_2) && is_inf(last_v_2) {
                    // OCCT L1463: PerformPlanes(F1, F2, Side, L1, L2).
                    perform_planes(f1, f2, side, l1, l2);
                    return;
                }
            }
        }
    }

    // create 3D curves on faces
    // OCCT L1470-1473: BRepLib::BuildCurves3d(F1/F2);
    // UpdateVertexTolerances(F1/F2).
    brep_lib_build_curves3d(f1);
    brep_lib_build_curves3d(f2);
    update_vertex_tolerances(f1);
    update_vertex_tolerances(f2);

    // OCCT L1475-1481 (D6): BOPAlgo_PaveFiller aPF;
    // NCollection_List<TopoDS_Shape> aLS; aLS.Append(F1); aLS.Append(F2);
    // aPF.SetArguments(aLS); aPF.Perform();
    //
    // rcad mapping (see the module header): the crate::bop PaveFiller
    // pipeline; the progress scope stands in for the default
    // Message_ProgressRange.
    let mut a_pf = PaveFiller::new();
    let mut a_ls: Vec<Shape> = Vec::new();
    a_ls.push(f1.clone());
    a_ls.push(f2.clone());
    a_pf.set_arguments(a_ls);
    let a_prog = rcad_kernel::message::NoopProgress;
    let a_ps = rcad_kernel::message::ProgressScope::new(&a_prog, "BRepOffset_Tool::Inter3D", 1);
    a_pf.perform(&a_ps);

    // OCCT L1483-1487: TrueEdges; if (!RefEdge.IsNull())
    // CheckIntersFF(aPF.PDS(), RefEdge, TrueEdges).
    let mut true_edges = OcctIndexedShapeMap::new();
    if !ref_edge.is_null() {
        check_inters_ff(a_pf.ds(), ref_edge, &mut true_edges);
    }

    // OCCT L1489-1490.
    let mut add_pcurve1 = true;
    let mut add_pcurve2 = true;

    // OCCT L1492-1495: pDS = aPF.PDS(); aFFs = pDS->InterfFF(); aNb =
    // aFFs.Length().
    //
    // DS-context annotation (the batch-2a form): the rcad DS deep-clones
    // the arguments; the F1/F2 handles the caller holds are NOT the DS
    // identities — the code below never matches faces by handle (the OCCT
    // loop walks every InterfFF entry), so only the section-edge reads go
    // through the DS (shape_at).
    let p_ds = a_pf.ds();
    let a_ffs = &p_ds.interf_ff;
    let a_nb = a_ffs.len();
    // Store Result
    // OCCT L1496-1500: L1.Clear(); L2.Clear(); O1, O2; BRep_Builder BB.
    l1.clear();
    l2.clear();
    let mut o1 = Orientation::Forward;
    let mut o2 = Orientation::Forward;
    //
    // OCCT L1502: aContext = aPF.Context() — the IntTools_Context consumed
    // only by the AdjustPCurveOnFace calls (GAP pass-through below).
    let _a_context = &a_pf.my_context;
    //
    // OCCT L1504-1591: the section-edge walk.
    for i in 0..a_nb {
        let a_ff = &a_ffs[i];
        let a_bcurves = &a_ff.curves;

        let a_nb_curves = a_bcurves.len();

        for j in 0..a_nb_curves {
            let a_bc = &p_ds.intersection_curves[a_bcurves[j]];
            let a_sect_edges = &a_bc.pave_blocks;

            for pb in a_sect_edges.iter() {
                // OCCT L1521-1523: nSect = aPB->Edge(); anEdge = pDS->Shape(nSect).
                let n_sect = pb.read().edge();
                let mut an_edge = p_ds.shape_at(n_sect);
                // OCCT L1524-1527: the TrueEdges filter.
                if !true_edges.is_empty() && !true_edges.contains(&an_edge) {
                    continue;
                }

                // OCCT L1529-1536: aC3DE = BRep_Tool::Curve(anEdge, f, l);
                // aC3DETrim = new Geom_TrimmedCurve(aC3DE, f, l).
                let (a_c3de, f_par, l_par) = match brep_tool_curve(&an_edge) {
                    Some(t) => t,
                    // OCCT: a null aC3DE keeps aC3DETrim null; no DS section
                    // edge reaches this form.
                    None => continue,
                };

                // OCCT L1538: aTolEdge = BRep_Tool::Tolerance(anEdge).
                let a_tol_edge = brep_tool_tolerance(&an_edge);

                // OCCT L1540-1558: the F1 pcurve attachment.
                if !has_curve_on_surface(&an_edge, f1) {
                    // OCCT L1542: aC2d = aBC.Curve().FirstCurve2d().
                    let mut a_c2d = match &a_bc.pcurve1 {
                        Some(c) => c.clone(),
                        None => {
                            // OCCT tolerates a null FirstCurve2d (the
                            // UpdateEdge with a null curve is skipped by
                            // the BB form); keep the OCCT skip.
                            continue;
                        }
                        #[allow(unreachable_patterns)]
                        _ => unreachable!(),
                    };
                    // OCCT L1543-1556: the AdjustPCurveOnFace forms — GAP
                    // pass-through (annotated above).
                    // if (aC3DE->IsPeriodic())
                    //   AdjustPCurveOnFace(F1, f, l, aC2d, aC2dNew, aContext);
                    // else
                    //   AdjustPCurveOnFace(F1, aC3DETrim, aC2d, aC2dNew, aContext);
                    if is_periodic_curve3(&a_c3de) {
                        a_c2d = adjust_pcurve_on_face_gap(&a_c2d);
                    } else {
                        a_c2d = adjust_pcurve_on_face_gap(&a_c2d);
                    }
                    // OCCT L1557: BB.UpdateEdge(anEdge, aC2d, F1, aTolEdge).
                    bat::builder_update_edge_pcurve(&mut an_edge, &a_c2d, f1, a_tol_edge);
                }

                // OCCT L1560-1578: the F2 pcurve attachment (idem).
                if !has_curve_on_surface(&an_edge, f2) {
                    // OCCT L1562: aC2d = aBC.Curve().SecondCurve2d().
                    let mut a_c2d = match &a_bc.pcurve2 {
                        Some(c) => c.clone(),
                        None => {
                            continue;
                        }
                    };
                    if is_periodic_curve3(&a_c3de) {
                        a_c2d = adjust_pcurve_on_face_gap(&a_c2d);
                    } else {
                        a_c2d = adjust_pcurve_on_face_gap(&a_c2d);
                    }
                    // OCCT L1577: BB.UpdateEdge(anEdge, aC2d, F2, aTolEdge).
                    bat::builder_update_edge_pcurve(&mut an_edge, &a_c2d, f2, a_tol_edge);
                }

                // OCCT L1580-1585: OrientSection + the Side reversal.
                orient_section(&an_edge, f1, f2, &mut o1, &mut o2);
                if side == State::Out {
                    o1 = bat::top_abs_reverse(o1);
                    o2 = bat::top_abs_reverse(o2);
                }

                // OCCT L1587-1588.
                l1.push(oriented(&an_edge, o1));
                l2.push(oriented(&an_edge, o2));
                let _ = (f_par, l_par);
            }
        }
    }

    // OCCT L1593-1622: aSameParTol; the isEl/addPCurve surface probes.
    let a_same_par_tol = rcad_kernel::precision::CONFUSION;
    let mut is_el1 = false;
    let mut is_el2 = false;

    let mut a_surf = face_surface_of(f1);
    if let Some(Surface3::Trimmed(ts)) = a_surf {
        a_surf = Some((*ts.basis).clone());
    }
    match &a_surf {
        Some(Surface3::Plane(_)) => {
            // OCCT L1601-1604: addPCurve1 = false.
            add_pcurve1 = false;
        }
        Some(s) if s.is_elementary() => {
            // OCCT L1605-1608: isEl1 = true.
            is_el1 = true;
        }
        _ => {}
    }

    let mut a_surf = face_surface_of(f2);
    if let Some(Surface3::Trimmed(ts)) = a_surf {
        a_surf = Some((*ts.basis).clone());
    }
    match &a_surf {
        Some(Surface3::Plane(_)) => {
            // OCCT L1615-1618: addPCurve2 = false.
            add_pcurve2 = false;
        }
        Some(s) if s.is_elementary() => {
            // OCCT L1619-1622: isEl2 = true.
            is_el2 = true;
        }
        _ => {}
    }

    // OCCT L1624-1702: remove the excess edges that are out of range.
    if l1.len() > 1 && (!is_el1 || !is_el2) && !the_ref_face1.is_null() {
        let (a_v1, a_v2) = top_exp_vertices(ref_edge);
        // only if RefEdge is open
        if !a_v1.is_same(&a_v2) {
            let a_ref_surf1 = face_surface_of(the_ref_face1);
            let a_ref_surf2 = face_surface_of(the_ref_face2);
            let any_closed = matches!(&a_ref_surf1, Some(s) if s.is_u_closed() || s.is_v_closed())
                || matches!(&a_ref_surf2, Some(s) if s.is_u_closed() || s.is_v_closed());
            if any_closed {
                // OCCT L1636-1640: MinAngleEdge; MinAngle; the RefEdge mid
                // point.
                let mut min_angle_edge = Shape::null();
                let mut min_angle = f64::INFINITY;
                let (ref_c, ref_f, ref_l) = match brep_tool_curve(ref_edge) {
                    Some(t) => t,
                    None => return,
                };
                let a_ref_pnt =
                    CurveEval::point_at(&ref_c, 0.5 * (ref_f + ref_l));

                // OCCT L1642-1678: the per-edge min-angle scan
                // (Extrema_ExtPC of aRefPnt onto the edge).
                let l1_snapshot = l1.clone();
                for it_value in &l1_snapshot {
                    let an_edge = it_value;
                    let (c, fe, le) = match brep_tool_curve(an_edge) {
                        Some(t) => t,
                        None => continue,
                    };
                    let a_mid_pnt_on_edge = CurveEval::point_at(&c, 0.5 * (fe + le));
                    let ref_to_mid = a_mid_pnt_on_edge - a_ref_pnt;

                    // OCCT L1647-1652: BRepAdaptor_Curve aBAcurve(anEdge);
                    // Extrema_ExtPC aProjector(aRefPnt, aBAcurve) — the real
                    // kernel body over the edge-range adaptor (TolF 1.0e-10).
                    let a_ba_adaptor = GeomCurveAdaptor::with_range(c.clone(), fe, le);
                    let a_ba_tool = CurveToolHandle::for_curve3(&c, &a_ba_adaptor, &a_ba_adaptor);
                    let a_projector = ExtremaExtPC::new_point_curve(a_ref_pnt, &a_ba_tool, 1.0e-10);
                    if a_projector.is_done() {
                        let mut imin = 0usize;
                        let mut min_sq_dist = f64::INFINITY;
                        for ind in 1..=a_projector.nb_ext() {
                            let a_sq_dist = a_projector.square_distance(ind);
                            if a_sq_dist < min_sq_dist {
                                min_sq_dist = a_sq_dist;
                                imin = ind;
                            }
                        }
                        if imin != 0 {
                            let a_projection_on_edge = a_projector.point(imin).point;
                            let ref_to_proj = a_projection_on_edge - a_ref_pnt;
                            let an_angle = ref_to_proj.angle_between(ref_to_mid);
                            if an_angle < min_angle {
                                min_angle = an_angle;
                                min_angle_edge = an_edge.clone();
                            }
                        }
                    }
                }

                // OCCT L1680-1699: keep only the MinAngleEdge entry.
                if !min_angle_edge.is_null() {
                    // OCCT L1685-1698: the while (itlist1.More()) walk with
                    // the paired Remove — the rcad index filter keeps the
                    // same pairing (L1 and L2 are element-wise parallel).
                    let mut keep: Vec<usize> = Vec::new();
                    for (idx, s) in l1.iter().enumerate() {
                        if s.is_same(&min_angle_edge) {
                            keep.push(idx);
                        }
                    }
                    let kept1: Vec<Shape> =
                        keep.iter().map(|&i| l1[i].clone()).collect();
                    let kept2: Vec<Shape> =
                        keep.iter().map(|&i| l2[i].clone()).collect();
                    *l1 = kept1;
                    *l2 = kept2;
                }
            } // OCCT L1700: } // if closed
        } // OCCT L1701: } // if (!aV1.IsSame(aV2))
    } // OCCT L1702: } // if (L1.Extent() > 1 && ...)

    // OCCT L1704-1930: the multi-edge assembly.
    if l1.len() > 1 && (!is_el1 || !is_el2) {
        let mut eseq: Vec<Shape> = Vec::new();
        let mut edges_for_concat: Vec<Shape> = Vec::new();

        if !true_edges.is_empty() {
            // OCCT L1709-1727: the TrueEdges assembly (1-based down loops).
            for i in (1..=true_edges.extent()).rev() {
                edges_for_concat.push(true_edges.at_1(i).clone());
            }
            let assembled_edge = assemble_edge(
                p_ds,
                f1,
                f2,
                add_pcurve1,
                add_pcurve2,
                &edges_for_concat,
            );
            if assembled_edge.is_null() {
                for i in (1..=true_edges.extent()).rev() {
                    eseq.push(true_edges.at_1(i).clone());
                }
            } else {
                eseq.push(assembled_edge);
            }
        } else {
            // OCCT L1728-1890: the wire-grouping assembly.
            let mut wseq: Vec<Shape> = Vec::new();
            let mut edges: Vec<Shape> = Vec::new();
            for it_value in l1.iter() {
                edges.push(it_value.clone());
            }
            // OCCT L1737-1790: while (!edges.IsEmpty()) — the AreConnex
            // wire grouping with the AngleWireEdge candidate selection.
            while !edges.is_empty() {
                let an_edge0 = edges[0].clone();
                let mut a_wire = bat::builder_make_wire();
                bat::builder_add_wire_edge(&mut a_wire, &an_edge0);
                let mut candidates: Vec<usize> = Vec::new();
                let mut res_wire = bat::builder_make_wire();
                let mut k = 0usize;
                for (wi, w) in wseq.iter().enumerate() {
                    // OCCT L1744-1752: resWire = TopoDS::Wire(wseq(k)).
                    res_wire = w.clone();
                    k = wi + 1; // OCCT k is 1-based
                    if are_connex_wire(&res_wire, &a_wire) {
                        candidates.push(1);
                        break;
                    }
                }
                if candidates.is_empty() {
                    wseq.push(a_wire.clone());
                    edges.remove(0);
                } else {
                    // OCCT L1758-1770: the candidate scan.
                    for j in 1..edges.len() {
                        // OCCT j is 2..=Length (1-based); rcad index j is
                        // the 0-based position (j == the 1-based index - 1
                        // offset handled below by the Candidates values).
                        let an_edge_j = edges[j].clone();
                        a_wire = bat::builder_make_wire();
                        bat::builder_add_wire_edge(&mut a_wire, &an_edge_j);
                        if are_connex_wire(&res_wire, &a_wire) {
                            // OCCT L1768: Candidates.Append(j) (1-based).
                            candidates.push(j + 1);
                        }
                    }
                    // OCCT L1771-1785: the min-angle candidate selection.
                    let mut minind = 1usize;
                    if candidates.len() > 1 {
                        let mut min_angle = f64::MAX;
                        for (cj, &cand) in candidates.iter().enumerate() {
                            let an_edge_c = edges[cand - 1].clone();
                            let an_angle = angle_wire_edge_wire(&res_wire, &an_edge_c);
                            if an_angle < min_angle {
                                min_angle = an_angle;
                                minind = cj + 1;
                            }
                        }
                    }
                    // OCCT L1786-1788: BB.Add(resWire, edges(Candidates(
                    // minind))); wseq(k) = resWire; edges.Remove(...).
                    let cand_edge = edges[candidates[minind - 1] - 1].clone();
                    bat::builder_add_wire_edge(&mut res_wire, &cand_edge);
                    wseq[k - 1] = res_wire;
                    edges.remove(candidates[minind - 1] - 1);
                }
            } // OCCT L1790: end of while

            // OCCT L1792-1889: for (i = 1; i <= wseq.Length(); i++).
            for i in 0..wseq.len() {
                let a_wire = wseq[i].clone();
                let mut a_local_edges_for_concat: Vec<Shape> = Vec::new();
                if shape_is_closed_flag(&a_wire) {
                    // OCCT L1796-1865: the closed-wire ordering walk.
                    let mut start_vertex = Shape::null();
                    let mut start_edge = Shape::null();
                    let mut start_found = false;
                    let mut elist: Vec<Shape> = Vec::new();

                    for e in bat::sub_shapes(&a_wire) {
                        let an_edge = e;
                        if start_found {
                            elist.push(an_edge);
                        } else {
                            let (v1, v2) = top_exp_vertices(&an_edge);
                            if !is_autonom_vertex(&v1, p_ds) {
                                start_vertex = v2;
                                start_edge = an_edge;
                                start_found = true;
                            } else if !is_autonom_vertex(&v2, p_ds) {
                                start_vertex = v1;
                                start_edge = an_edge;
                                start_found = true;
                            } else {
                                elist.push(an_edge);
                            }
                        }
                    }
                    if !start_found {
                        // OCCT L1833-1841.
                        start_edge = elist[0].clone();
                        elist.remove(0);
                        let (v1, _v2) = top_exp_vertices(&start_edge);
                        start_vertex = v1;
                    }
                    a_local_edges_for_concat.push(start_edge.clone());
                    // OCCT L1843-1865: while (!Elist.IsEmpty()).
                    while !elist.is_empty() {
                        let mut consumed = false;
                        for idx in 0..elist.len() {
                            let an_edge = elist[idx].clone();
                            let (v1, v2) = top_exp_vertices(&an_edge);
                            if v1.is_same(&start_vertex) {
                                start_vertex = v2;
                                a_local_edges_for_concat.push(an_edge);
                                elist.remove(idx);
                                consumed = true;
                                break;
                            } else if v2.is_same(&start_vertex) {
                                start_vertex = v1;
                                a_local_edges_for_concat.push(an_edge);
                                elist.remove(idx);
                                consumed = true;
                                break;
                            }
                        }
                        if !consumed {
                            // OCCT loops forever on a disconnected list;
                            // keep the OCCT termination (the sequences in
                            // the covered cases are connected).
                            break;
                        }
                    }
                } else {
                    // OCCT L1867-1874: the BRepTools_WireExplorer walk —
                    // the rcad sub-shape order (the WireExplorer order
                    // re-host is the same stored-edge order; annotated).
                    for e in bat::sub_shapes(&a_wire) {
                        a_local_edges_for_concat.push(e);
                    }
                }

                // OCCT L1876-1888: the per-wire assembly.
                let assembled_edge = assemble_edge(
                    p_ds,
                    f1,
                    f2,
                    add_pcurve1,
                    add_pcurve2,
                    &a_local_edges_for_concat,
                );
                if assembled_edge.is_null() {
                    for j in (0..a_local_edges_for_concat.len()).rev() {
                        eseq.push(a_local_edges_for_concat[j].clone());
                    }
                } else {
                    eseq.push(assembled_edge);
                }
            } // OCCT L1889: for (i = 1; i <= wseq.Length(); i++)
        } // OCCT L1890: end of else (when TrueEdges is empty)

        // OCCT L1892-1929: the tolerance-filtered rebuild of L1/L2.
        if eseq.len() < l1.len() {
            l1.clear();
            l2.clear();
            for i in 0..eseq.len() {
                let an_edge = eseq[i].clone();
                // OCCT L1900: BRepLib::SameParameter(anEdge, aSameParTol, true).
                brep_lib_same_parameter(&an_edge, a_same_par_tol);
                // OCCT L1901: EdgeTol = BRep_Tool::Tolerance(anEdge).
                let edge_tol = brep_tool_tolerance(&an_edge);
                if edge_tol > 1e-2 {
                    continue;
                }

                if edge_tol >= 1e-4 {
                    // OCCT L1912-1913: ReconstructPCurves(anEdge);
                    // BRepLib::SameParameter(anEdge, aSameParTol, true).
                    reconstruct_pcurves(&an_edge);
                    brep_lib_same_parameter(&an_edge, a_same_par_tol);
                }

                // OCCT L1919-1927.
                orient_section(&an_edge, f1, f2, &mut o1, &mut o2);
                if side == State::Out {
                    o1 = bat::top_abs_reverse(o1);
                    o2 = bat::top_abs_reverse(o2);
                }

                l1.push(oriented(&an_edge, o1));
                l2.push(oriented(&an_edge, o2));
            }
        }
    } // OCCT L1930: end of if (L1.Extent() > 1)
    else {
        // OCCT L1932-1940.
        let l1_snapshot = l1.clone();
        for it_value in &l1_snapshot {
            let an_edge = it_value;
            brep_lib_same_parameter(an_edge, a_same_par_tol);
        }
    }
}

/// OCCT IsInf(theVal) (cxx L4608-4611) — the value is close to TheInfini.
pub(crate) fn is_inf(the_val: f64) -> bool {
    the_val > THE_INFINI * 0.9
}

/// The trimmed-surface last-U probe of the Inter3D planar fast path (OCCT
/// BRepAdaptor_Surface(F, true)::LastUParameter) — the rcad trimmed box
/// form; a bare surface keeps its infinite default domain.
fn trimmed_last_u(s: &Option<Surface3>) -> f64 {
    match s {
        Some(Surface3::Trimmed(ts)) => ts.trim[1],
        Some(other) => other.default_domain()[1],
        None => 0.0,
    }
}

/// OCCT BRepAdaptor_Surface(F, true)::LastVParameter — idem.
fn trimmed_last_v(s: &Option<Surface3>) -> f64 {
    match s {
        Some(Surface3::Trimmed(ts)) => ts.trim[3],
        Some(other) => other.default_domain()[3],
        None => 0.0,
    }
}

/// OCCT aC3DE->IsPeriodic() — the rcad variant form (the circle / ellipse
/// BSpline-with-closed-knot cases collapse to the Circle/Ellipse variants;
/// the remaining rcad curves are non-periodic).
fn is_periodic_curve3(c: &Curve3) -> bool {
    matches!(c, Curve3::Circle(_) | Curve3::Ellipse(_))
}

/// The AreConnex static on bare wires (the module-a re-host exposed for the
/// Inter3D grouping loop).
fn are_connex_wire(w1: &Shape, w2: &Shape) -> bool {
    are_connex(w1, w2)
}

/// The AngleWireEdge static on bare wires (the module-a re-host exposed for
/// the Inter3D grouping loop).
fn angle_wire_edge_wire(a_wire: &Shape, an_edge: &Shape) -> f64 {
    angle_wire_edge(a_wire, an_edge)
}

/// OCCT aWire.Closed() — the CLOSED flag read.
fn shape_is_closed_flag(s: &Shape) -> bool {
    bat::shape_is_closed(s)
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::TryProject (hxx L97-103; cxx L1945-2006).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::TryProject(F1, F2, Edges, LInt1, LInt2, Side,
/// TolConf) (cxx L1945-2006) — finds whether the edges are laying on F1.
#[allow(clippy::too_many_arguments)]
pub fn try_project(
    f1: &Shape,
    f2: &Shape,
    edges: &[Shape],
    l_int1: &mut Vec<Shape>,
    l_int2: &mut Vec<Shape>,
    side: State,
    tol_conf: f64,
) -> bool {
    // try to find if the edges <Edges> are laying on the face F1.
    // OCCT L1955-1962.
    l_int1.clear();
    l_int2.clear();
    let mut is_ok = true;
    let mut ok = true;
    let mut o1 = Orientation::Forward;
    let mut o2 = Orientation::Forward;
    let bouchon = match face_surface_of(f1) {
        Some(s) => s,
        None => return ok,
    };

    // OCCT L1964-2004: the edge loop.
    for it_value in edges {
        // OCCT L1966-1978: the curve read with the BuildCurve3d fallback.
        let mut cur_e = it_value.clone();
        let (mut c, f, l) = match brep_tool_curve(&cur_e) {
            Some(t) => t,
            None => {
                // OCCT L1972: BRepLib::BuildCurve3d(CurE, Tolerance) — GAP
                // no-op (annotated above); the re-read stays null for a
                // degenerated edge and the OCCT continue applies.
                brep_lib_build_curves3d(&cur_e);
                match brep_tool_curve(&cur_e) {
                    Some(t) => t,
                    None => continue, // not 3d curve, can be degenerated, need to skip
                }
            }
        };
        // OCCT L1979-1983: C = new Geom_TrimmedCurve(C, f, l); the location
        // transform — the rcad surface/curve storage is location-baked
        // (identity in the covered cases; annotated).
        c = trimmed_curve3(&c, f, l);

        // OCCT L1984-1985: isOk = IsOnSurface(C, Bouchon, TolConf,
        // TolReached).
        let mut tol_reached = 0.0f64;
        is_ok = is_on_surface(&c, &bouchon, tol_conf, &mut tol_reached);

        if is_ok {
            // OCCT L1989: B.UpdateEdge(CurE, TolReached).
            b_update_edge_tolerance(&mut cur_e, tol_reached);
            // OCCT L1990: BuildPCurves(CurE, F1) — the module-a static.
            build_pcurves(&cur_e, f1);
            // OCCT L1991-1996.
            orient_section(&cur_e, f1, f2, &mut o1, &mut o2);
            if side == State::Out {
                o1 = bat::top_abs_reverse(o1);
                o2 = bat::top_abs_reverse(o2);
            }
            // OCCT L1997-1998.
            l_int1.push(oriented(&cur_e, o1));
            l_int2.push(oriented(&cur_e, o2));
        } else {
            // OCCT L2002: Ok = false.
            ok = false;
        }
    }
    // OCCT L2005: return Ok.
    ok
}

// ---------------------------------------------------------------------------
// OCCT static ExtentEdge (cxx L2072-2157) — the Tool-local (F, EF, E, NE)
// form (the Inter3d module carries its own two-argument static at cxx
// L63-87).
// ---------------------------------------------------------------------------

/// OCCT static ExtentEdge(F, EF, E, NE) (cxx L2072-2157) — the tangential
/// extension of E to the F surface bounds, then the new edge built on the
/// extended 2D curve.  The dependencies are the translated ones:
/// `GeomAPI::To3d` / `GeomAPI::To2d`
/// (`crate::geomalgo::geom_api`) and `GeomLib::ExtendCurveToPoint`
/// (`crate::geomalgo::geom_lib`).
pub(crate) fn extent_edge_tool(f: &Shape, ef: &Shape, e: &Shape, ne: &mut Shape) {
    // OCCT L2077-2081: CE = BRepAdaptor_Curve(E); Type; aLocalEdge =
    // E.EmptyCopied(); NE = TopoDS::Edge(aLocalEdge).
    let (ce_c, ce_f, ce_l) = match brep_tool_curve(e) {
        Some(t) => t,
        None => return,
    };
    let type_is_analytic = matches!(
        ce_c,
        Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_)
            | Curve3::Parabola(_)
    );
    *ne = empty_copied(e);

    // OCCT L2083-2087: the analytic early return.
    if type_is_analytic {
        return;
    }
    // Extension en tangence jusqu'au bord de la surface.
    // OCCT L2089-2098: PMax; S = the face surface; the bounds clamp.
    let p_max = 1.0e2f64;
    let s = match face_surface_of(f) {
        Some(v) => v,
        None => return,
    };
    let (mut umin, mut umax, mut vmin, mut vmax) = {
        let d = s.default_domain();
        (d[0], d[1], d[2], d[3])
    };
    umin = umin.max(-p_max);
    vmin = vmin.max(-p_max);
    umax = umax.min(p_max);
    vmax = vmax.min(p_max);

    // OCCT L2100-2101: C2d = BRep_Tool::CurveOnSurface(E, F, f, l).
    let c2d = match brep_tool_curve_on_surface(e, f) {
        Some((c, _, _)) => c,
        None => return,
    };

    // calcul point cible. ie point d'intersection du prolongement tangent
    // et des bords.
    // OCCT L2104-2119: the first-endpoint tangent target PF2d.
    let p = Curve2dEval::point_at(&c2d, ce_f);
    let mut tang = Curve2dEval::derivative_at(&c2d, ce_f);
    let mut tx = f64::INFINITY;
    let mut ty = f64::INFINITY;
    if tang.x.abs() > rcad_kernel::precision::CONFUSION {
        tx = ((umax - p.x) / tang.x).abs().min(((umin - p.x) / tang.x).abs());
    }
    if tang.y.abs() > rcad_kernel::precision::CONFUSION {
        ty = ((vmax - p.y) / tang.y).abs().min(((vmin - p.y) / tang.y).abs());
    }
    let tmin = tx.min(ty);
    tang = tang * tmin;
    let pf2d = DVec2::new(p.x - tang.x, p.y - tang.y);

    // OCCT L2121-2133: the last-endpoint tangent target PL2d.
    let p = Curve2dEval::point_at(&c2d, ce_l);
    let mut tang = Curve2dEval::derivative_at(&c2d, ce_l);
    let mut tx = f64::INFINITY;
    let mut ty = f64::INFINITY;
    if tang.x.abs() > rcad_kernel::precision::CONFUSION {
        tx = ((umax - p.x) / tang.x).abs().min(((umin - p.x) / tang.x).abs());
    }
    if tang.y.abs() > rcad_kernel::precision::CONFUSION {
        ty = ((vmax - p.y) / tang.y).abs().min(((vmin - p.y) / tang.y).abs());
    }
    let tmin = tx.min(ty);
    tang = tang * tmin;
    let pl2d = DVec2::new(p.x + tang.x, p.y + tang.y);

    // OCCT L2135: CC = GeomAPI::To3d(C2d, gp_Pln(gp::XOY())).
    let a_xoy = Plane::new(DVec3::ZERO, DVec3::Z);
    let cc = crate::geomalgo::geom_api::to3d(&c2d, &a_xoy);
    // OCCT L2137-2138: PF/PL — the 2D tangent targets lifted to z = 0.
    let pf = DVec3::new(pf2d.x, pf2d.y, 0.);
    let pl = DVec3::new(pl2d.x, pl2d.y, 0.);

    // OCCT L2140-2143: occ::down_cast<Geom_BoundedCurve>(CC); if
    // (ExtC.IsNull()) return; — the rcad Geom_BoundedCurve instantiations.
    if !matches!(cc, Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)) {
        return;
    }
    let mut ext_c = cc;

    // OCCT L2145-2146: GeomLib::ExtendCurveToPoint(ExtC, PF, 1, false);
    // GeomLib::ExtendCurveToPoint(ExtC, PL, 1, true);
    crate::geomalgo::geom_lib::extend_curve_to_point(&mut ext_c, pf, 1, false);
    crate::geomalgo::geom_lib::extend_curve_to_point(&mut ext_c, pl, 1, true);

    // OCCT L2148: CNE2d = GeomAPI::To2d(ExtC, gp_Pln(gp::XOY())).
    // The rcad None is the null result handle of GeomAPI::To2d (GeomAPI.cxx
    // L46-49); OCCT then dereferences it at L2152 (`B.Range(NE,
    // CNE2d->FirstParameter(), ...)`), which the expect below reports.
    let cne2d = crate::geomalgo::geom_api::to2d(&ext_c, &a_xoy)
        .expect("GeomAPI::To2d null result (OCCT L2152 null dereference)");

    // Construction de la nouvelle arrete (OCCT L2150-2156).
    // B.MakeEdge(NE);
    *ne = bat::builder_make_edge();
    // B.UpdateEdge(NE, CNE2d, EF, BRep_Tool::Tolerance(E));
    bat::builder_update_edge_pcurve(ne, &cne2d, ef, brep_tool_tolerance(e));
    // B.Range(NE, CNE2d->FirstParameter(), CNE2d->LastParameter());
    let a_cne_dom = Curve2dEval::default_domain(&cne2d);
    bat::builder_range_edge(ne, a_cne_dom[0], a_cne_dom[1]);
    // NE.Orientation(E.Orientation());
    ne.orientation = e.orientation;
}

// ---------------------------------------------------------------------------
// OCCT static ProjectVertexOnEdge (cxx L2161-2249).
// ---------------------------------------------------------------------------

/// OCCT static ProjectVertexOnEdge(V, E, TolConf) (cxx L2161-2249).
pub(crate) fn project_vertex_on_edge(v: &mut Shape, e: &Shape, tol_conf: f64) -> bool {
    // OCCT L2163-2167: f, l, U, L, found.
    let mut u = 0.0f64;
    let mut found = false;

    // OCCT L2169-2172: P = BRep_Tool::Pnt(V); C = BRepAdaptor_Curve(E).
    let p = match bat::brep_tool_pnt(v) {
        Some(t) => t,
        None => return false,
    };
    let (c, f, l) = match brep_tool_curve(e) {
        Some(t) => t,
        None => return false,
    };

    // OCCT L2174-2185: the FORWARD-endpoint equality.
    if v.orientation == Orientation::Forward {
        if f.abs() < rcad_kernel::precision::INFINITE_VALUE {
            let pf = CurveEval::point_at(&c, f);
            if pf.distance(p) <= tol_conf {
                u = f;
                found = true;
            }
        }
    }
    // OCCT L2186-2197: the REVERSED-endpoint equality.
    if v.orientation == Orientation::Reversed {
        if !found && l.abs() < rcad_kernel::precision::INFINITE_VALUE {
            let pl = CurveEval::point_at(&c, l);
            if pl.distance(p) <= tol_conf {
                u = l;
                found = true;
            }
        }
    }
    // OCCT L2198-2216: the Extrema_ExtPC projection — Extrema_ExtPC Proj(P, C)
    // over the edge-range adaptor (TolF 1.0e-10, the real kernel body).
    if !found {
        let a_proj_adaptor = GeomCurveAdaptor::with_range(c.clone(), f, l);
        let a_proj_tool = CurveToolHandle::for_curve3(&c, &a_proj_adaptor, &a_proj_adaptor);
        let proj = ExtremaExtPC::new_point_curve(p, &a_proj_tool, 1.0e-10);
        if proj.is_done() && proj.nb_ext() > 0 {
            let mut dist2_min = proj.square_distance(1);
            u = proj.point(1).param;
            for i in 2..=proj.nb_ext() {
                let dist2 = proj.square_distance(i);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    u = proj.point(i).param;
                }
            }
            found = true;
        }
    }

    // OCCT L2240-2247: the B.UpdateVertex form (the FORWARD-oriented edge
    // view; the INTERNAL-oriented vertex).
    if found {
        let mut ee = e.clone();
        ee.orientation = Orientation::Forward;
        let mut v_int = v.clone();
        v_int.orientation = Orientation::Internal;
        b_update_vertex_on_edge(&mut ee, &v_int, u, brep_tool_tolerance(e));
    }
    // OCCT L2248: return found.
    found
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::Inter2d (hxx L111-115; cxx L2253-2532).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::Inter2d(F, E1, E2, LV, TolConf) (cxx L2253-2532) —
/// the pcurve/pcurve intersection of E1 and E2 on F.
pub fn inter2d(
    f: &Shape,
    e1: &Shape,
    e2: &Shape,
    lv: &mut Vec<Shape>,
    tol_conf: f64,
) {
    // OCCT L2260-2261: fl1[2], fl2[2]; LV.Clear().
    let mut fl1 = [0.0f64; 2];
    let mut fl2 = [0.0f64; 2];
    lv.clear();

    // Si l edge a ete etendu les pcurves ne sont pas forcement
    // a jour.
    // OCCT L2265-2266: BuildPCurves(E1, F); BuildPCurves(E2, F).
    build_pcurves(e1, f);
    build_pcurves(e2, f);

    // OCCT L2282-2290: NbPC1/NbPC2 (the seam second pcurves).
    let mut nb_pc1 = 1usize;
    let mut nb_pc2 = 1usize;
    if bat::brep_tool_is_closed_on_surface(e1, f) {
        nb_pc1 += 1;
    }
    if bat::brep_tool_is_closed_on_surface(e2, f) {
        nb_pc2 += 1;
    }

    // OCCT L2292-2295: S = the face surface; C1, C2; YaSol; itry.
    let s = face_surface_of(f);
    let mut c1: Option<Curve2d> = None;
    let mut c2: Option<Curve2d> = None;
    let mut ya_sol = false;
    let mut itry = 0usize;

    // OCCT L2297-2490: while (!YaSol && itry < 2).
    while !ya_sol && itry < 2 {
        for i in 1..=nb_pc1 {
            // OCCT L2301-2309: C1 = CurveOnSurface(E1 [or Reversed], F,
            // fl1) — the rcad pcurve store is keyed by the face; the
            // reversed-edge seam lookup returns the same entry (annotated;
            // the closed-edge second representation is staged with the
            // pcurve store extension).
            let (cc1, rf1, rl1) = match brep_tool_curve_on_surface(e1, f) {
                Some(t) => t,
                None => return, // OCCT_DEBUG return form: no pcurve
            };
            c1 = Some(cc1);
            fl1[0] = rf1;
            fl1[1] = rl1;
            for j in 1..=nb_pc2 {
                // OCCT L2315-2323: C2 = CurveOnSurface(E2 [or Reversed], F, fl2).
                let (cc2, rf2, rl2) = match brep_tool_curve_on_surface(e2, f) {
                    Some(t) => t,
                    None => return,
                };
                c2 = Some(cc2);
                fl2[0] = rf2;
                fl2[1] = rl2;

                let mut u1 = 0.0f64;
                let mut u2 = 0.0f64;
                let mut p2d = DVec2::ZERO;
                let mut a_current_find = false;
                if itry == 1 {
                    // OCCT L2337-2343: the adaptor-range refresh.
                    if let Some(cc1) = &c1 {
                        let d = curve2d_domain(cc1, &fl1);
                        fl1[0] = d.0;
                        fl1[1] = d.1;
                    }
                    if let Some(cc2) = &c2 {
                        let d = curve2d_domain(cc2, &fl2);
                        fl2[0] = d.0;
                        fl2[1] = d.1;
                    }
                }

                // OCCT L2347-2433: the itry == 0 proximity branches.
                if itry == 0 {
                    let mut p1 = [DVec2::ZERO; 2];
                    let mut p2 = [DVec2::ZERO; 2];
                    if let Some(cc1) = &c1 {
                        p1[0] = Curve2dEval::point_at(cc1, fl1[0]);
                        p1[1] = Curve2dEval::point_at(cc1, fl1[1]);
                    }
                    if let Some(cc2) = &c2 {
                        p2[0] = Curve2dEval::point_at(cc2, fl2[0]);
                        p2[1] = Curve2dEval::point_at(cc2, fl2[1]);
                    }

                    // OCCT L2355-2373: the endpoint coincidence scan.
                    for i1 in 0..2usize {
                        for i2 in 0..2usize {
                            if fl1[i1].abs() < rcad_kernel::precision::INFINITE_VALUE
                                && fl2[i2].abs() < rcad_kernel::precision::INFINITE_VALUE
                            {
                                // gp_Pnt2d::IsEqual(P, TolConf) — the
                                // coordinate-wise tolerance form.
                                let d = p1[i1] - p2[i2];
                                if d.x.abs() <= tol_conf && d.y.abs() <= tol_conf {
                                    ya_sol = true;
                                    a_current_find = true;
                                    u1 = fl1[i1];
                                    u2 = fl2[i2];
                                    if let Some(cc1) = &c1 {
                                        p2d = Curve2dEval::point_at(cc1, u1);
                                    }
                                }
                            }
                        }
                    }
                    // OCCT L2374-2402: the Extrema_ExtPC2d of P1 onto AC2.
                    if !ya_sol {
                        for i1 in 0..2usize {
                            if let Some(cc2) = &c2 {
                                let (f2r, l2r) = (fl2[0], fl2[1]);
                                let mut extr = ExtPC2d::new(
                                    p1[i1],
                                    cc2,
                                    rcad_kernel::precision::CONFUSION,
                                    f2r,
                                    l2r,
                                );
                                extr.perform(p1[i1], cc2, f2r, l2r);
                                if extr.is_done() && extr.nb_ext() > 0 {
                                    let mut dist2_min = extr.square_distance(1);
                                    let mut index_min = 1usize;
                                    for ind in 2..=extr.nb_ext() {
                                        let dist2 = extr.square_distance(ind);
                                        if dist2 < dist2_min {
                                            dist2_min = dist2;
                                            index_min = ind;
                                        }
                                    }
                                    if dist2_min
                                        <= rcad_kernel::precision::CONFUSION
                                            * rcad_kernel::precision::CONFUSION
                                    {
                                        ya_sol = true;
                                        a_current_find = true;
                                        p2d = p1[i1];
                                        u1 = fl1[i1];
                                        u2 = extr.point(index_min).param;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    // OCCT L2404-2432: the Extrema_ExtPC2d of P2 onto AC1.
                    if !ya_sol {
                        for i2 in 0..2usize {
                            if let Some(cc1) = &c1 {
                                let (f1r, l1r) = (fl1[0], fl1[1]);
                                let mut extr = ExtPC2d::new(
                                    p2[i2],
                                    cc1,
                                    rcad_kernel::precision::CONFUSION,
                                    f1r,
                                    l1r,
                                );
                                extr.perform(p2[i2], cc1, f1r, l1r);
                                if extr.is_done() && extr.nb_ext() > 0 {
                                    let mut dist2_min = extr.square_distance(1);
                                    let mut index_min = 1usize;
                                    for ind in 2..=extr.nb_ext() {
                                        let dist2 = extr.square_distance(ind);
                                        if dist2 < dist2_min {
                                            dist2_min = dist2;
                                            index_min = ind;
                                        }
                                    }
                                    if dist2_min
                                        <= rcad_kernel::precision::CONFUSION
                                            * rcad_kernel::precision::CONFUSION
                                    {
                                        ya_sol = true;
                                        a_current_find = true;
                                        p2d = p2[i2];
                                        u2 = fl2[i2];
                                        u1 = extr.point(index_min).param;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L2436-2472: the Geom2dInt_GInter general branch —
                // GAP carrier (architecture difference #24): the general
                // 2d curve/curve intersection engine is not translated;
                // the OCCT branch keeps its structure behind the GAP (the
                // covered cases resolve through the proximity branches).
                if !ya_sol {
                    let _ = (u1, u2, p2d);
                    if general_2d_intersects_gap(&c1, &c2, tol_conf) {
                        // unreachable while the engine is a GAP.
                    }
                }
                // OCCT L2474-2486: the intersection-vertex build.
                if a_current_find {
                    let p = match &s {
                        Some(ss) => ss.point_at(p2d.x, p2d.y),
                        None => return,
                    };
                    // OCCT L2477: V = BRepLib_MakeVertex(P).
                    let mut v = bat::builder_make_vertex();
                    bat::builder_update_vertex_point_tol(&mut v, p, 0.0);
                    // OCCT L2478: V.Orientation(TopAbs_INTERNAL).
                    v.orientation = Orientation::Internal;
                    // OCCT L2479-2482: B.UpdateVertex(V, U1/U2, the
                    // FORWARD-oriented edges, TolConf).
                    let mut e1_fwd = e1.clone();
                    e1_fwd.orientation = Orientation::Forward;
                    b_update_vertex_on_edge(&mut e1_fwd, &v, u1, tol_conf);
                    let mut e2_fwd = e2.clone();
                    e2_fwd.orientation = Orientation::Forward;
                    b_update_vertex_on_edge(&mut e2_fwd, &v, u2, tol_conf);
                    // OCCT L2485: LV.Append(V).
                    lv.push(v);
                }
            }
        }
        // OCCT L2489: itry++.
        itry += 1;
    }

    // OCCT L2492-2524: keep only the extremity vertices.
    if lv.len() > 1 {
        //------------------------------------------------
        // garde seulement les vertex les plus proches du
        // debut et de la fin.
        //------------------------------------------------
        let mut vf = Shape::null();
        let mut vl = Shape::null();
        let mut u_min = f64::INFINITY;
        let mut u_max = f64::NEG_INFINITY;

        for it_value in lv.iter() {
            let cv = it_value.clone();
            // OCCT L2507-2508: U = BRep_Tool::Parameter(CV, the
            // FORWARD-oriented E1).
            let mut e1_fwd = e1.clone();
            e1_fwd.orientation = Orientation::Forward;
            let u = bat::brep_tool_parameter(&cv, &e1_fwd);
            if u < u_min {
                vf = cv.clone();
                u_min = u;
            }
            if u > u_max {
                vl = cv.clone();
                u_max = u;
            }
        }
        lv.clear();
        lv.push(vf);
        lv.push(vl);
    }
}

/// OCCT Geom2dAdaptor_Curve::FirstParameter/LastParameter (the bounded
/// range of the pcurve; the trimmed form carries it, the bare form keeps
/// the edge range).
fn curve2d_domain(c: &Curve2d, fallback: &[f64; 2]) -> (f64, f64) {
    match c {
        Curve2d::Trimmed(tc) => (tc.t_min, tc.t_max),
        _ => (fallback[0], fallback[1]),
    }
}

/// OCCT Geom2dInt_GInter(AC1, AC2, TolConf, TolConf) — GAP carrier
/// (architecture difference #24).  The general 2d curve/curve intersection
/// engine is not translated; the carrier keeps the OCCT call form and the
/// annotation fires for the cases the proximity branches do not resolve.
fn general_2d_intersects_gap(_c1: &Option<Curve2d>, _c2: &Option<Curve2d>, _tol: f64) -> bool {
    // OCCT: !Inter.IsEmpty() && (Inter.NbPoints() > 0 ||
    // Inter.NbSegments() > 0) — unreachable while the engine is a GAP.
    false
}

// ---------------------------------------------------------------------------
// GAP carriers shared with the sibling modules (architecture difference
// #24) — moved from module a to respect the 2000-line file limit; the OCCT
// call sites keep the single import path through the module-a glob.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// GAP carriers shared with the sibling modules (architecture difference
// #24) — declared here so the OCCT call sites keep a single import path.
// ---------------------------------------------------------------------------

/// OCCT Geom2dConvert_CompCurveToBSplineCurve (TKGeomBase/Geom2dConvert) —
/// GAP carrier (architecture difference #24).
pub(crate) struct Geom2dConvertCompCurveToBSplineCurve;

impl Geom2dConvertCompCurveToBSplineCurve {
    /// OCCT Geom2dConvert_CompCurveToBSplineCurve(BasisCurve).
    pub fn new(_basis_curve: &Curve2d) -> Self {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve (TKGeomBase not translated)");
    }
    /// OCCT Geom2dConvert_CompCurveToBSplineCurve::Add(Curve, Tol, After).
    pub fn add(&mut self, _curve: &Curve2d, _tol: f64, _after: bool) -> bool {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve::Add");
    }
    /// OCCT Geom2dConvert_CompCurveToBSplineCurve::BSplineCurve().
    pub fn bspline_curve(&self) -> Option<Curve2d> {
        panic!("GAP: Geom2dConvert_CompCurveToBSplineCurve::BSplineCurve");
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve (TKGeomBase/GeomConvert) — GAP
/// carrier (architecture difference #24).
pub(crate) struct GeomConvertCompCurveToBSplineCurve;

impl GeomConvertCompCurveToBSplineCurve {
    /// OCCT GeomConvert_CompCurveToBSplineCurve(BasisCurve).
    pub fn new(_basis_curve: &Curve3) -> Self {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve (TKGeomBase not translated)");
    }
    /// OCCT GeomConvert_CompCurveToBSplineCurve::Add(Curve, Tol, After).
    pub fn add(&mut self, _curve: &Curve3, _tol: f64, _after: bool) -> bool {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve::Add");
    }
    /// OCCT GeomConvert_CompCurveToBSplineCurve::BSplineCurve().
    pub fn bspline_curve(&self) -> Option<Curve3> {
        panic!("GAP: GeomConvert_CompCurveToBSplineCurve::BSplineCurve");
    }
}

/// OCCT Geom2dConvert_ApproxCurve (TKGeomBase/Geom2dConvert) — GAP carrier
/// (architecture difference #24).
pub(crate) struct Geom2dConvertApproxCurve;

impl Geom2dConvertApproxCurve {
    /// OCCT Geom2dConvert_ApproxCurve(Curve, Tol, Order, MaxSegments,
    /// MaxDegree).
    pub fn new(
        _curve: &Curve2d,
        _tol: f64,
        _order: GeomAbsShapeKind,
        _max_segments: i32,
        _max_degree: i32,
    ) -> Self {
        panic!("GAP: Geom2dConvert_ApproxCurve (TKGeomBase not translated)");
    }
    /// OCCT Geom2dConvert_ApproxCurve::HasResult().
    pub fn has_result(&self) -> bool {
        panic!("GAP: Geom2dConvert_ApproxCurve::HasResult");
    }
    /// OCCT Geom2dConvert_ApproxCurve::Curve().
    pub fn curve(&self) -> Curve2d {
        panic!("GAP: Geom2dConvert_ApproxCurve::Curve");
    }
}

/// OCCT GeomConvert_ApproxCurve (TKGeomBase/GeomConvert) — GAP carrier
/// (architecture difference #24).
pub(crate) struct GeomConvertApproxCurve;

impl GeomConvertApproxCurve {
    /// OCCT GeomConvert_ApproxCurve(Curve, Tol, Order, MaxSegments,
    /// MaxDegree).
    pub fn new(
        _curve: &Curve3,
        _tol: f64,
        _order: GeomAbsShapeKind,
        _max_segments: i32,
        _max_degree: i32,
    ) -> Self {
        panic!("GAP: GeomConvert_ApproxCurve (TKGeomBase not translated)");
    }
    /// OCCT GeomConvert_ApproxCurve::HasResult().
    pub fn has_result(&self) -> bool {
        panic!("GAP: GeomConvert_ApproxCurve::HasResult");
    }
    /// OCCT GeomConvert_ApproxCurve::Curve().
    pub fn curve(&self) -> Curve3 {
        panic!("GAP: GeomConvert_ApproxCurve::Curve");
    }
}

/// OCCT GeomAbs_Shape (the continuity kind) — the offset-module alias of the
/// GeomAbsShapeKind carrier from brep_offset_offset.rs.
pub(crate) use super::brep_offset_offset::GeomAbsShapeKind;

/// OCCT ProjLib_ProjectedCurve (TKTopAlgo/ProjLib) — GAP carrier
/// (architecture difference #24).
pub(crate) struct ProjLibProjectedCurve;

impl ProjLibProjectedCurve {
    /// OCCT ProjLib_ProjectedCurve(S, C, Tol).
    pub fn new(_s: &Surface3, _c: &Curve3, _tol: f64) -> Self {
        panic!("GAP: ProjLib_ProjectedCurve (TKTopAlgo not translated)");
    }
}

/// OCCT GCPnts_AbscissaPoint (TKGeomAlgo/GCPnts) — GAP carrier (architecture
/// difference #24).
pub(crate) struct GCPntsAbscissaPoint;

impl GCPntsAbscissaPoint {
    /// OCCT GCPnts_AbscissaPoint::Length(C).
    pub fn length(_c: &Curve3) -> f64 {
        panic!("GAP: GCPnts_AbscissaPoint::Length (TKGeomAlgo not translated)");
    }
}

/// OCCT ShapeCustom_Curve2d::ConvertToLine2d (TKShHealing/ShapeCustom) —
/// GAP carrier (architecture difference #24).
pub(crate) struct ShapeCustomCurve2d;

impl ShapeCustomCurve2d {
    /// OCCT ShapeCustom_Curve2d::ConvertToLine2d(Curve, f, l, TolConv,
    /// newFpar, newLpar, deviation).
    pub fn convert_to_line2d(
        _curve: &Curve2d,
        _f: f64,
        _l: f64,
        _tol_conv: f64,
        _new_fpar: &mut f64,
        _new_lpar: &mut f64,
        _deviation: &mut f64,
    ) -> Option<Line2d> {
        panic!("GAP: ShapeCustom_Curve2d::ConvertToLine2d (TKShHealing not translated)");
    }
}

/// OCCT GeomProjLib::Curve2d(C3d, f, l, Surface) (TKTopAlgo/GeomProjLib) —
/// GAP carrier (architecture difference #24; the
/// loc_ope_wires_on_shape_b.rs #10 precedent).
pub(crate) fn geom_proj_lib_curve2d(_c3d: &Curve3, _f: f64, _l: f64, _surf: &Surface3) -> Option<Curve2d> {
    panic!("GAP: GeomProjLib::Curve2d (TKTopAlgo not translated)");
}

// OCCT GeomAPI::To3d / GeomAPI::To2d (TKGeomAlgo/GeomAPI) live in their OCCT
// toolkit home now: `crate::geomalgo::geom_api::{to3d, to2d}`.  The former
// GAP carriers of this module are deleted.
