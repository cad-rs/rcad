// OCCT HLRTopoBRep_FaceIsoLiner (TKHLR/HLRTopoBRep/HLRTopoBRep_FaceIsoLiner.hxx
// L1-58 + .cxx L1-488) — builds the isoparametric curves of a face with a
// hatcher (Perform) and the (vertex, iso-edge) helpers MakeVertex /
// MakeIsoLine feeding the HLRTopoBRep_Data.
//
// Kernel-boundary note: the OCCT code reads BRep_Tool / BRepTools /
// BRepAdaptor_Surface / TopExp free functions over TopoDS handles; rcad keeps
// those reads on the landed TShape graph (`Shape::data`) — the established
// Shape-only convention of this crate (see the per-helper OCCT line
// references below).  The Perform signature is fixed by the DSFiller
// contract: no `BRep` handle crosses the boundary, so the pcurve key is
// `(face TShape pointer, face location)` and the location table is out of
// reach (identity-location HLR faces are unaffected).
//
// Standalone-shape note (OCCT BRep_Builder::MakeVertex / MakeEdge produce
// heap TShapes that are never "null"): rcad's `Shape::is_null` encoding marks
// index == usize::MAX, which a standalone `Shape::new` TShape would collide
// with; the fresh shapes here carry the `STANDALONE_INDEX` sentinel (largest
// non-null index) so the OCCT `IsNull` tests at cxx L292/L404 keep their
// meaning.

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Line2d, Plane, Point3, Surface3, SurfaceEval, TrimmedCurve2, Vec3,
};
use rcad_kernel::precision::{is_negative_infinite_value, is_positive_infinite_value, PCONFUSION};
use rcad_kernel::topods::{tshape_flags, Orientation, Shape, TEdgeData, TShape, TVertexData};

use crate::brep_algo::tool::brep_tool_tolerance;
use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::hatch::hatcher::Hatcher;
use crate::geomalgo::hatch::intersector::HatchIntersector;
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};

use super::data::Data;

// OCCT .cxx L43-47 — the module constants.
const INTERSECTOR_CONFUSION: f64 = 1.0e-10;
const INTERSECTOR_TANGENCY: f64 = 1.0e-10;
const HATCHER_CONFUSION2D: f64 = 1.0e-8;
const HATCHER_CONFUSION3D: f64 = 1.0e-8;
const INFINITE: f64 = 100.;

/// The index encoding of a standalone (heap TShape) shape built by the
/// BRep_Builder stand-ins below; never `Shape::is_null` (see the module
/// note).
const STANDALONE_INDEX: usize = usize::MAX - 1;

//=======================================================================
// Function : Perform
// Purpose  : Builds isoparametric curves with a hatcher.
//=======================================================================

/// OCCT HLRTopoBRep_FaceIsoLiner::Perform (cxx L54-415).
pub fn perform(fi: usize, f: &Shape, ds: &mut Data, nb_isos: usize) {
    // OCCT L59: (void)FI; // avoid compiler warning
    let _ = fi;

    // OCCT L61: double UMin, UMax, VMin, VMax, U1, U2; — the OCCT locals stay
    // uninitialized until assigned; the neutral 0.0 keeps Rust defined.
    let mut u_min;
    let mut u_max;
    let mut v_min;
    let mut v_max;
    let mut u1 = 0.0;
    let mut u2 = 0.0;
    // OCCT L62: int ne = 0;
    let mut ne: usize = 0;
    // OCCT L63: // BRep_Builder Builder; (commented out in OCCT).
    // OCCT L64: TopoDS_Edge Edge; — dead local kept for form.
    let _edge = Shape::null();
    // OCCT L65: TopExp_Explorer ExpEdges; — Init at L109 below.
    // OCCT L66-67: TopoDS_Face TF = F; TF.Orientation(TopAbs_FORWARD);
    let mut tf = f.clone();
    tf.orientation = Orientation::Forward;
    // OCCT L68-70: gp_Pnt2d P; gp_Pnt P1, P2; TopoDS_Vertex V1U, V2U, V1V, V2V;
    // (the OCCT default ctors give (0,0)/(0,0,0); the vertices start null and
    // may keep the value of a previous domain — the OCCT form).
    let mut p = DVec2::ZERO;
    let mut p1 = Point3::ZERO;
    let mut p2 = Point3::ZERO;
    let mut v1u = Shape::null();
    let mut v2u = Shape::null();
    let mut v1v = Shape::null();
    let mut v2v = Shape::null();

    // OCCT L72-73: the intersector and the hatcher (the OCCT ctor default
    // KeepSeg = false, hxx L48-49; rcad passes it explicitly).
    let intersector =
        HatchIntersector::with_tolerances(INTERSECTOR_CONFUSION, INTERSECTOR_TANGENCY);
    let mut hatcher = Hatcher::new(
        intersector,
        HATCHER_CONFUSION2D,
        HATCHER_CONFUSION3D,
        true,
        false,
    );

    // OCCT L75: BRepTools::UVBounds(TF, UMin, UMax, VMin, VMax);
    let (b_u_min, b_u_max, b_v_min, b_v_max) = uv_bounds(&tf);
    u_min = b_u_min;
    u_max = b_u_max;
    v_min = b_v_min;
    v_max = b_v_max;
    // OCCT L76-79.
    let infinite_u_min = is_negative_infinite_value(u_min);
    let infinite_u_max = is_positive_infinite_value(u_max);
    let infinite_v_min = is_negative_infinite_value(v_min);
    let infinite_v_max = is_positive_infinite_value(v_max);

    // OCCT L81-93.
    if infinite_u_min && infinite_u_max {
        u_min = -INFINITE;
        u_max = INFINITE;
    } else if infinite_u_min {
        u_min = u_max - INFINITE;
    } else if infinite_u_max {
        u_max = u_min + INFINITE;
    }

    // OCCT L95-107.
    if infinite_v_min && infinite_v_max {
        v_min = -INFINITE;
        v_max = INFINITE;
    } else if infinite_v_min {
        v_min = v_max - INFINITE;
    } else if infinite_v_max {
        v_max = v_min + INFINITE;
    }

    // OCCT L109-114: for (ExpEdges.Init(TF, TopAbs_EDGE); ExpEdges.More();
    //                   ExpEdges.Next()) { ne++; } — Edges of the face TF.
    let exp_edges = top_exp_edges(&tf);
    for _e in &exp_edges {
        ne += 1;
    }

    // OCCT L116-124: if (DS.FaceHasIntL(TF)) — OutLines built on face TF.
    if ds.face_has_int_l(&tf) {
        let int_l = ds.face_int_l(&tf);
        // NCollection_List<TopoDS_Shape>::Iterator itE.
        let mut it_e = 0;
        while it_e < int_l.len() {
            ne += 1;
            it_e += 1;
        }
    }

    // OCCT L126-127: NCollection_Array1<TopoDS_Shape> SH(1, ne);
    //                NCollection_Array1<bool> IL(1, ne); — internal OutLine.
    // The arrays are 1-based in OCCT; the [0] slots stay unused.  The C++
    // IL array is uninitialized, but every slot is assigned before its read;
    // the neutral false keeps Rust defined.
    let mut sh = vec![Shape::null(); ne + 1];
    let mut il = vec![false; ne + 1];

    // OCCT L129-154: the face edges as hatcher elements.
    let mut exp_i = 0;
    while exp_i < exp_edges.len() {
        // OCCT L132: const TopoDS_Edge& newE = TopoDS::Edge(ExpEdges.Current());
        let new_e = exp_edges[exp_i].clone();
        // OCCT L133: PC = BRep_Tool::CurveOnSurface(newE, TF, U1, U2);
        let (pc, c_u1, c_u2) = curve_on_surface(&new_e, &tf);
        u1 = c_u1;
        u2 = c_u2;
        let inde;
        // OCCT L134-138.
        if (Curve2dAdaptor::first_parameter(&pc) - u1).abs() <= PCONFUSION
            && (Curve2dAdaptor::last_parameter(&pc) - u2).abs() <= PCONFUSION
        {
            inde = hatcher.add_element_from_curve(&pc, new_e.orientation);
        } else {
            // OCCT L141-143: TPC = new Geom2d_TrimmedCurve(PC, U1, U2);
            //                Geom2dAdaptor_Curve aGAC(TPC);
            let tpc = Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(pc.clone()),
                t_min: u1,
                t_max: u2,
            });
            inde = hatcher.add_element(&tpc, new_e.orientation);
        }
        // OCCT L145: SH(IndE) = newE;
        sh[inde] = new_e.clone();
        // OCCT L146-153.
        if ds.is_out_l_face_edge(&tf, &new_e) {
            il[inde] = true;
        } else {
            il[inde] = false;
        }
        exp_i += 1;
    }

    // OCCT L156-178: if (DS.FaceHasIntL(TF)) — get the internal OutLines
    // built on face F.
    if ds.face_has_int_l(&tf) {
        let int_l = ds.face_int_l(&tf);
        // NCollection_List<TopoDS_Shape>::Iterator itE.
        let mut it_e = 0;
        while it_e < int_l.len() {
            // OCCT L162: const TopoDS_Edge& newE = TopoDS::Edge(itE.Value());
            let new_e = int_l[it_e].clone();
            // OCCT L163: PC = BRep_Tool::CurveOnSurface(newE, TF, U1, U2);
            let (pc, c_u1, c_u2) = curve_on_surface(&new_e, &tf);
            u1 = c_u1;
            u2 = c_u2;
            let inde;
            // OCCT L164-173 — the elements are added as TopAbs_INTERNAL.
            if (Curve2dAdaptor::first_parameter(&pc) - u1).abs() <= PCONFUSION
                && (Curve2dAdaptor::last_parameter(&pc) - u2).abs() <= PCONFUSION
            {
                inde = hatcher.add_element_from_curve(&pc, Orientation::Internal);
            } else {
                let tpc = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(pc.clone()),
                    t_min: u1,
                    t_max: u2,
                });
                inde = hatcher.add_element(&tpc, Orientation::Internal);
            }
            // OCCT L175-176.
            sh[inde] = new_e;
            il[inde] = true;
            it_e += 1;
        }
    }

    //-----------------------------------------------------------------------
    // Creation des hachures.
    //-----------------------------------------------------------------------

    // OCCT L184-185: BRepAdaptor_Surface Surface(TF);
    //                double Tolerance = BRep_Tool::Tolerance(TF);
    let surface = surface_adaptor(&tf);
    let tolerance = brep_tool_tolerance(&tf);

    // OCCT L187-191.
    let delta_u = (u_max - u_min).abs();
    let delta_v = (v_max - v_min).abs();
    let confusion = f64::min(delta_u, delta_v) * HATCHER_CONFUSION3D;
    hatcher.set_confusion3d(confusion);

    //-----------------------------------------------------------------------
    // Courbes Iso U.
    //-----------------------------------------------------------------------

    // OCCT L197: double StepU = DeltaU / (double)nbIsos;
    let step_u = delta_u / nb_isos as f64;
    if step_u > confusion {
        // OCCT L200-201: UPrm = UMin + StepU / 2.; Dir(gp_Dir2d::D::Y);
        let mut u_prm = u_min + step_u / 2.;
        let dir = DVec2::Y;
        // OCCT L203: for (IIso = 1; IIso <= nbIsos; IIso++)
        for _i_iso in 1..=nb_isos {
            // OCCT L205-206: Ori(UPrm, 0.); IsoLine = new Geom2d_Line(Ori, Dir);
            let ori = DVec2::new(u_prm, 0.);
            let iso_line = Curve2d::Line(Line2d::new(ori, dir));

            // OCCT L208-210: aGAC(IsoLine); IndH = AddHatching(aGAC); Trim(IndH);
            let ind_h = hatcher.add_hatching(&iso_line);
            hatcher.trim_hatching(ind_h);
            // OCCT L211-214.
            if hatcher.trim_done(ind_h) && !hatcher.trim_failed(ind_h) {
                // OCCT L213: Hatcher.ComputeDomains(IndH) — exactly one
                // hatching is bound at this point (RemHatching at the end of
                // every iteration), so the rcad no-arg ComputeDomains
                // computes the IndH-th one.
                hatcher.compute_domains();
            }
            // OCCT L215-241: the OCCT_DEBUG status dump is compiled out; the
            // failure path removes the hatching and continues.
            if !hatcher.is_done(ind_h) {
                hatcher.rem_hatching(ind_h);
                continue;
            }

            // OCCT L243: int NbDom = Hatcher.NbDomains(IndH);
            let nb_dom = hatcher.nb_domains(ind_h);
            if nb_dom > 0 {
                // OCCT L247: for (int IDom = 1; IDom <= NbDom; IDom++)
                for i_dom in 1..=nb_dom {
                    // OCCT L249: const HatchGen_Domain& Dom = Hatcher.Domain(IndH, IDom);
                    let dom = hatcher.domain(ind_h, i_dom);
                    // OCCT L250-251.
                    let u11 = if dom.has_first_point() {
                        dom.first_point().parameter()
                    } else {
                        v_min - INFINITE
                    };
                    let u21 = if dom.has_second_point() {
                        dom.second_point().parameter()
                    } else {
                        v_max + INFINITE
                    };
                    // OCCT L252-255: IsoLine->D0(U11, P); Surface.D0(P.X(), P.Y(), P1);
                    //                IsoLine->D0(U21, P); Surface.D0(P.X(), P.Y(), P2);
                    p = Curve2dEval::point_at(&iso_line, u11);
                    p1 = surface.value(p.x, p.y);
                    p = Curve2dEval::point_at(&iso_line, u21);
                    p2 = surface.value(p.x, p.y);
                    // OCCT L256-273: Iso U - Premier point.
                    if dom.has_first_point() {
                        let pnt_h = dom.first_point();
                        // OCCT L260: for (int IPntE = 1; IPntE <= PntH.NbPoints(); IPntE++)
                        for i_pnt_e in 1..=pnt_h.nb_points() {
                            // OCCT L262: const HatchGen_PointOnElement& PntE = PntH.Point(IPntE);
                            let pnt_e = pnt_h.point(i_pnt_e);
                            // OCCT L263-267.
                            v1u = make_vertex(
                                &sh[pnt_e.index() as usize],
                                p1,
                                pnt_e.parameter(),
                                tolerance,
                                ds,
                            );
                            // OCCT L268-271.
                            if il[pnt_e.index() as usize] {
                                ds.add_out_v(&v1u);
                            }
                        }
                    }
                    // OCCT L274-291: Iso U - Deuxieme point.
                    if dom.has_second_point() {
                        let pnt_h = dom.second_point();
                        for i_pnt_e in 1..=pnt_h.nb_points() {
                            let pnt_e = pnt_h.point(i_pnt_e);
                            // OCCT L281-285.
                            v2u = make_vertex(
                                &sh[pnt_e.index() as usize],
                                p2,
                                pnt_e.parameter(),
                                tolerance,
                                ds,
                            );
                            // OCCT L286-289.
                            if il[pnt_e.index() as usize] {
                                ds.add_out_v(&v2u);
                            }
                        }
                    }
                    // OCCT L292-295.
                    if !v1u.is_null() && !v2u.is_null() {
                        make_iso_line(f, &iso_line, &mut v1u, &mut v2u, u11, u21, tolerance, ds);
                    }
                }
            }

            // OCCT L299-300.
            hatcher.rem_hatching(ind_h);
            u_prm += step_u;
        }
    }

    //-----------------------------------------------------------------------
    // Courbes Iso V.
    //-----------------------------------------------------------------------

    // OCCT L308: double StepV = DeltaV / (double)nbIsos;
    let step_v = delta_v / nb_isos as f64;
    if step_v > confusion {
        // OCCT L311-312: VPrm = VMin + StepV / 2.; Dir(gp_Dir2d::D::X);
        let mut v_prm = v_min + step_v / 2.;
        let dir = DVec2::X;
        // OCCT L314: for (IIso = 1; IIso <= nbIsos; IIso++)
        for _i_iso in 1..=nb_isos {
            // OCCT L316-317: Ori(0., VPrm); IsoLine = new Geom2d_Line(Ori, Dir);
            let ori = DVec2::new(0., v_prm);
            let iso_line = Curve2d::Line(Line2d::new(ori, dir));

            // OCCT L319-321.
            let ind_h = hatcher.add_hatching(&iso_line);
            hatcher.trim_hatching(ind_h);
            // OCCT L322-325.
            if hatcher.trim_done(ind_h) && !hatcher.trim_failed(ind_h) {
                hatcher.compute_domains();
            }
            // OCCT L326-352: the OCCT_DEBUG status dump is compiled out.
            if !hatcher.is_done(ind_h) {
                hatcher.rem_hatching(ind_h);
                continue;
            }

            // OCCT L354.
            let nb_dom = hatcher.nb_domains(ind_h);
            if nb_dom > 0 {
                // OCCT L358.
                for i_dom in 1..=nb_dom {
                    // OCCT L360.
                    let dom = hatcher.domain(ind_h, i_dom);
                    // OCCT L361-362.
                    let u12 = if dom.has_first_point() {
                        dom.first_point().parameter()
                    } else {
                        v_min - INFINITE
                    };
                    let u22 = if dom.has_second_point() {
                        dom.second_point().parameter()
                    } else {
                        v_max + INFINITE
                    };
                    // OCCT L363-366.
                    p = Curve2dEval::point_at(&iso_line, u12);
                    p1 = surface.value(p.x, p.y);
                    p = Curve2dEval::point_at(&iso_line, u22);
                    p2 = surface.value(p.x, p.y);
                    // OCCT L367-385: Iso V - Premier point.
                    if dom.has_first_point() {
                        let pnt_h = dom.first_point();
                        for i_pnt_e in 1..=pnt_h.nb_points() {
                            let pnt_e = pnt_h.point(i_pnt_e);
                            // OCCT L374-378.
                            v1v = make_vertex(
                                &sh[pnt_e.index() as usize],
                                p1,
                                pnt_e.parameter(),
                                tolerance,
                                ds,
                            );
                            // OCCT L380-383.
                            if il[pnt_e.index() as usize] {
                                ds.add_out_v(&v1v);
                            }
                        }
                    }
                    // OCCT L386-403: Iso V - Deuxieme point.
                    if dom.has_second_point() {
                        let pnt_h = dom.second_point();
                        for i_pnt_e in 1..=pnt_h.nb_points() {
                            let pnt_e = pnt_h.point(i_pnt_e);
                            // OCCT L393-397.
                            v2v = make_vertex(
                                &sh[pnt_e.index() as usize],
                                p2,
                                pnt_e.parameter(),
                                tolerance,
                                ds,
                            );
                            // OCCT L398-401.
                            if il[pnt_e.index() as usize] {
                                ds.add_out_v(&v2v);
                            }
                        }
                    }
                    // OCCT L404-407.
                    if !v1v.is_null() && !v2v.is_null() {
                        make_iso_line(f, &iso_line, &mut v1v, &mut v2v, u12, u22, tolerance, ds);
                    }
                }
            }

            // OCCT L411-412.
            hatcher.rem_hatching(ind_h);
            v_prm += step_v;
        }
    }
}

//=================================================================================================

/// OCCT HLRTopoBRep_FaceIsoLiner::MakeVertex (cxx L419-463).
pub fn make_vertex(e: &Shape, p: Point3, par: f64, tol: f64, ds: &mut Data) -> Shape {
    // OCCT L425-427: TopoDS_Vertex V, VF, VL; BRep_Builder B;
    //                TopExp::Vertices(E, VF, VL);
    let mut v = Shape::null();
    let (v_f, v_l) = top_exp_vertices(e);
    // OCCT L428-431.
    if pnt_is_equal(p, brep_tool_pnt(&v_f), brep_tool_tolerance(&v_f)) {
        return v_f;
    }
    // OCCT L432-435.
    if pnt_is_equal(p, brep_tool_pnt(&v_l), brep_tool_tolerance(&v_l)) {
        return v_l;
    }

    // OCCT L437-453: for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex())
    ds.init_vertex(e);
    while ds.more_vertex() {
        // OCCT L439-440.
        let cur_v = ds.vertex();
        let cur_p = ds.parameter();
        // OCCT L441-445.
        if pnt_is_equal(p, brep_tool_pnt(&cur_v), brep_tool_tolerance(&cur_v)) {
            v = cur_v;
            break;
        } else if par < cur_p {
            // OCCT L448-451: B.MakeVertex(V, P, Tol);
            //                V.Orientation(TopAbs_INTERNAL);
            //                DS.InsertBefore(V, Par);
            v = brep_make_vertex(p, tol);
            v.orientation = Orientation::Internal;
            ds.insert_before(&v, par);
            break;
        }
        // the for-loop increment (skipped by the breaks above).
        ds.next_vertex();
    }

    // OCCT L455-460.
    if !ds.more_vertex() {
        v = brep_make_vertex(p, tol);
        v.orientation = Orientation::Internal;
        ds.append(&v, par);
    }

    // OCCT L462: return V;
    v
}

//=================================================================================================

/// OCCT HLRTopoBRep_FaceIsoLiner::MakeIsoLine (cxx L467-488).
pub fn make_iso_line(
    f: &Shape,
    iso: &Curve2d,
    v1: &mut Shape,
    v2: &mut Shape,
    u1: f64,
    u2: f64,
    tol: f64,
    ds: &mut Data,
) {
    // OCCT L476-478: BRep_Builder B; TopoDS_Edge E;
    //                E.Orientation(TopAbs_INTERNAL); — the orientation does
    //                not survive MakeEdge below (TopoDS_Builder::MakeShape,
    //                TopoDS_Builder.cxx L28-33, resets TopAbs_FORWARD); the
    //                literal sequence is kept and the built edge is FORWARD.
    let mut ed = brep_make_edge();
    // OCCT L479-480.
    v1.orientation = Orientation::Forward;
    v2.orientation = Orientation::Reversed;
    // OCCT L481: B.MakeEdge(E); — the fresh TEdgeData (above).
    // OCCT L482: B.UpdateEdge(E, Iso, F, Tol);
    brep_update_edge_pcurve(&mut ed, iso, f, tol);
    // OCCT L483: B.Add(E, V1);
    brep_add_edge_vertex(&mut ed, v1);
    // OCCT L484: B.UpdateVertex(V1, U1, E, Tol);
    brep_update_vertex_on_edge(&mut ed, v1, u1, tol);
    // OCCT L485: B.Add(E, V2);
    brep_add_edge_vertex(&mut ed, v2);
    // OCCT L486: B.UpdateVertex(V2, U2, E, Tol);
    brep_update_vertex_on_edge(&mut ed, v2, u2, tol);
    // The BRep_Builder conversation is complete: the TShape gets its Arc
    // handle (the heap allocation of OCCT).
    let e = Shape::from_parts(Arc::new(TShape::Edge(ed)), STANDALONE_INDEX, 0, Orientation::Forward);
    // OCCT L487: DS.AddIsoL(F).Append(E);
    ds.add_iso_l(f).push(e);
}

// ===========================================================================
//  Kernel-boundary stand-ins (the OCCT free functions over TopoDS handles).
//  Each helper cites its OCCT origin; the reads stay on the Shape's TShape
//  graph because the Perform contract carries no BRep handle.
// ===========================================================================

/// OCCT TopExp_Explorer(TF, TopAbs_EDGE) — the edges of every wire of the
/// face, each with the cumulative orientation
/// `TF.Orientation() x W.Orientation() x E.Orientation()`
/// (TopExp_Explorer.cxx, cumOri).  TF is FORWARD here (cxx L67).
fn top_exp_edges(tf: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    if let TShape::Face(fd) = &*tf.data {
        for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
            if let TShape::Wire(wd) = &*w.data {
                let w_or = tf.orientation.compose(w.orientation);
                for e in &wd.edges {
                    let mut ce = e.clone();
                    ce.orientation = w_or.compose(e.orientation);
                    out.push(ce);
                }
            }
        }
    }
    out
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) (TopExp.cxx L214-250, CumOri =
/// false) — the FORWARD child gives Vfirst, the REVERSED child gives Vlast.
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let mut v_first = Shape::null();
    let mut v_last = Shape::null();
    // TopoDS_Iterator ite(E) — the edge's direct children (myShapes).
    if let TShape::Edge(ed) = &*e.data {
        for sv in &ed.my_shapes {
            if let TShape::Vertex(_) = &*sv.data {
                match sv.orientation {
                    Orientation::Forward => v_first = sv.clone(),
                    Orientation::Reversed => v_last = sv.clone(),
                    _ => {}
                }
            }
        }
    }
    (v_first, v_last)
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, First, Last) (BRep_Tool.cxx L345-368)
/// — the pcurve of E on F with its range.  The rcad Shape-only read keys the
/// edge pcurve map by `(face TShape pointer, face location)` — the
/// established convention of the boundary collectors (fclass2d /
/// topol_tool_brep).  Unbound: the OCCT null handle would crash the caller —
/// the panic mirrors it.
fn curve_on_surface(e: &Shape, f: &Shape) -> (Curve2d, f64, f64) {
    if let TShape::Edge(ed) = &*e.data {
        let key = (f.ptr_id(), f.location);
        if let Some((pc, t1, t2)) = ed.pcurves.get(&key) {
            return (pc.clone(), *t1, *t2);
        }
    }
    panic!("Standard_NoSuchObject: BRep_Tool::CurveOnSurface");
}

/// OCCT BRep_Tool::Surface(F) — the face surface (no Location; the
/// identity-location convention).
fn face_surface(f: &Shape) -> Option<&Surface3> {
    match &*f.data {
        TShape::Face(fd) => fd.surface.as_ref(),
        _ => None,
    }
}

/// OCCT BRep_Tool::Pnt(V) (BRep_Tool.cxx) — the vertex 3D point (the
/// identity-location convention).
fn brep_tool_pnt(v: &Shape) -> Point3 {
    match &*v.data {
        TShape::Vertex(vd) => vd.point,
        _ => panic!("Standard_NoSuchObject: BRep_Tool::Pnt"),
    }
}

/// OCCT gp_Pnt::IsEqual(Other, Tolerance) (gp_Pnt.hxx L125-129) —
/// Distance(Other) <= Tolerance.
fn pnt_is_equal(p: Point3, other: Point3, tolerance: f64) -> bool {
    p.distance(other) <= tolerance
}

/// OCCT BRepTools::UVBounds(F, UMin, UMax, VMin, VMax) (BRepTools.cxx L64-80)
/// + AddUVBounds(F, B) (L126-157): the 2D box of every edge pcurve of the
/// FORWARD-oriented face (BndLib_Add2dCurve::Add, L185) clamped to the face
/// UV window for the non-periodic directions (L202-361; the B-spline
/// periodicity verification rides on the declared periodicity — the
/// established builder_face precedent), with the surface natural bounds as
/// the empty-box fallback (L139-153).  Returns (UMin, UMax, VMin, VMax).
fn uv_bounds(f: &Shape) -> (f64, f64, f64, f64) {
    // OCCT L128-129: TopoDS_Face F = FF; F.Orientation(TopAbs_FORWARD);
    let mut tf = f.clone();
    tf.orientation = Orientation::Forward;

    // OCCT L132-137: the per-edge boxes folded into aBox.
    let mut a_box: Option<[f64; 4]> = None;
    for e in top_exp_edges(&tf) {
        // OCCT L179-183: AddUVBounds(F, E, B) — the pcurve (a null pcurve
        // returns; rcad raises — see curve_on_surface).
        let (pc, t1, t2) = curve_on_surface(&e, &tf);
        // OCCT L185: BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC).
        let b = rcad_kernel::curve2d_bounding_box(&pc, t1, t2, 0.0);
        // OCCT L191-361: the non-periodic window clamp.
        let b = clamp_to_face_window(&tf, b);
        a_box = match a_box {
            None => Some(b),
            Some(prev) => Some([
                prev[0].min(b[0]),
                prev[1].max(b[1]),
                prev[2].min(b[2]),
                prev[3].max(b[3]),
            ]),
        };
    }

    // OCCT L139-153: the empty box takes the surface natural bounds; a null
    // surface leaves the box void (UVBounds then yields zeros, L76-79).
    let box_ = match a_box {
        Some(b) => b,
        None => match face_surface(&tf) {
            Some(surf) => SurfaceEval::default_domain(surf),
            None => return (0.0, 0.0, 0.0, 0.0),
        },
    };
    // OCCT L74: B.Get(UMin, VMin, UMax, VMax).
    (box_[0], box_[1], box_[2], box_[3])
}

/// OCCT BRepTools::AddUVBounds(F, E, B) L191-361 — clamp the pcurve box to
/// the face UV window for the non-periodic surface directions.
fn clamp_to_face_window(f: &Shape, mut b: [f64; 4]) -> [f64; 4] {
    if let TShape::Face(fd) = &*f.data {
        if let Some(uvd) = fd.uv_domain {
            if let Some(surf) = &fd.surface {
                if !SurfaceEval::is_u_periodic(surf) {
                    if b[0] < uvd[0] && uvd[0] < b[1] {
                        b[0] = uvd[0];
                    }
                    if b[0] < uvd[1] && uvd[1] < b[1] {
                        b[1] = uvd[1];
                    }
                }
                if !SurfaceEval::is_v_periodic(surf) {
                    if b[2] < uvd[2] && uvd[2] < b[3] {
                        b[2] = uvd[2];
                    }
                    if b[2] < uvd[3] && uvd[3] < b[3] {
                        b[3] = uvd[3];
                    }
                }
            }
        }
    }
    b
}

/// OCCT BRepAdaptor_Surface Surface(TF) (BRepAdaptor_Surface.cxx,
/// Initialize(F, R = true)) — the GeomAdaptor over the face surface with the
/// face UV window.  The rcad Shape-only stand-in reads the stored surface +
/// uv_domain (the BRep location table is unreachable from the Perform
/// signature; identity-location HLR faces are unaffected).
fn surface_adaptor(f: &Shape) -> GeomSurfaceAdapter {
    let world = face_surface(f).cloned().unwrap_or(Surface3::Plane(Plane {
        origin: Point3::ZERO,
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_dir: Vec3::new(1.0, 0.0, 0.0),
        v_dir: Vec3::new(0.0, 1.0, 0.0),
    }));
    let uv_domain = match &*f.data {
        TShape::Face(fd) => fd.uv_domain,
        _ => None,
    };
    match uv_domain {
        Some(d) => GeomSurfaceAdapter::with_domain(world, d),
        None => GeomSurfaceAdapter::new(world),
    }
}

// ===========================================================================
//  BRep_Builder stand-ins (the OCCT builder conversation of MakeVertex /
//  MakeIsoLine; the TShape data is assembled before the Arc handoff because
//  the landed TShape is immutable once shared).
// ===========================================================================

/// OCCT BRep_Builder::MakeVertex(V, P, Tol) (BRep_Builder.lxx L161-165 +
/// BRep_Builder.cxx L1203-1214) — a fresh heap TShape with the point and the
/// tolerance (UpdateVertex(V, P, Tol)).
fn brep_make_vertex(p: Point3, tol: f64) -> Shape {
    let tshape = Arc::new(TShape::Vertex(TVertexData {
        my_shapes: Vec::new(),
        flags: tshape_flags::FREE
            | tshape_flags::MODIFIED
            | tshape_flags::ORIENTABLE
            | tshape_flags::CLOSED
            | tshape_flags::CONVEX,
        point: p,
        tolerance: tol,
        points: Vec::new(),
    }));
    Shape::from_parts(tshape, STANDALONE_INDEX, 0, Orientation::Forward)
}

/// OCCT BRep_Builder::MakeEdge(E) (BRep_Builder.cxx L623-631 +
/// TopoDS_Builder::MakeShape L28-33) — the fresh BRep_TEdge defaults
/// (SameParameter/SameRange true, tolerance 0, range (0,0)).
fn brep_make_edge() -> TEdgeData {
    TEdgeData {
        my_shapes: Vec::new(),
        flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
        curve: None,
        first: Shape::null(),
        last: Shape::null(),
        range: [0.0, 0.0],
        degenerated: false,
        pcurves: indexmap::IndexMap::new(),
        representations: Vec::new(),
        vertex_params: std::collections::HashMap::new(),
        tolerance: 0.0,
        same_parameter: true,
        same_range: true,
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C, F, Tol) (BRep_Builder.cxx L679-701 +
/// the lxx (E, C, F, Tol) overload) — the pcurve is stored under
/// `L.Predivided(E.Location())` (a fresh edge has an identity location, so
/// the face location value); the fresh BRep_CurveOnSurface range is the
/// BRep_GCurve default (0,0) until UpdateVertex sets it.  The edge tolerance
/// is max-ed with Tol (BRep_TEdge::UpdateTolerance).
fn brep_update_edge_pcurve(ed: &mut TEdgeData, pc: &Curve2d, f: &Shape, tol: f64) {
    let key = (f.ptr_id(), f.location);
    ed.pcurves.insert(key, (pc.clone(), 0.0, 0.0));
    ed.tolerance = ed.tolerance.max(tol);
}

/// OCCT TopoDS_Builder::Add(E, V) (TopoDS_Builder.cxx L37-103) — the vertex
/// is appended to the edge TShape's myShapes with its own orientation (the
/// component is also frozen — Free(false) — out of reach through the shared
/// Arc<TShape>, and not consulted by the HLR consumers).  The rcad first /
/// last fields carry the FORWARD / REVERSED children (the kernel convention
/// of BRep_Tool::Vertices).
fn brep_add_edge_vertex(ed: &mut TEdgeData, v: &Shape) {
    ed.my_shapes.push(v.clone());
    if v.orientation == Orientation::Forward {
        ed.first = v.clone();
    } else {
        ed.last = v.clone();
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, Par, E, Tol) (BRep_Builder.cxx
/// L1319-1410) — the vertex is searched in the edge (IsSame) and its stored
/// orientation selects the parameter slot: FORWARD sets the curve
/// representation First, REVERSED sets Last.  The vertex tolerance is max-ed
/// with Tol.  The rcad vertex_params entry is the established encoding of
/// the vertex parameter read by BRepTool::parameter_on_edge.
fn brep_update_vertex_on_edge(ed: &mut TEdgeData, v: &Shape, par: f64, tol: f64) {
    // TopoDS_Iterator itv(E.Oriented(TopAbs_FORWARD)) — find V among the
    // edge's vertices.
    let mut ori = Orientation::Internal;
    for v_cur in &ed.my_shapes {
        if v.is_same(v_cur) {
            ori = v_cur.orientation;
            if ori == v.orientation {
                break;
            }
        }
    }
    match ori {
        Orientation::Forward => {
            // GC->First(Par) — the pcurve representation's first parameter.
            for (_, first, _) in ed.pcurves.values_mut() {
                *first = par;
            }
            ed.range[0] = par;
        }
        Orientation::Reversed => {
            // GC->Last(Par).
            for (_, _, last) in ed.pcurves.values_mut() {
                *last = par;
            }
            ed.range[1] = par;
        }
        _ => {}
    }
    ed.vertex_params.insert(v.ptr_id(), par);
    // TV->UpdateTolerance(Tol) — the max is out of reach through the shared
    // Arc<TShape>; the fresh MakeVertex shapes already carry Tol and the
    // shared ones are read-only here.
}

// ===========================================================================
//  Anchor tests (analytically assertable).
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::topo_brep::v_data::tests::face_outer_edge;
    use crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face;
    use rcad_kernel::topods::BRep;

    /// OCCT anchor: Perform on the unit square face with nbIsos = 1 — one U
    /// iso at x = 0.5 and one V iso at y = 0.5, both spanning the face
    /// (range [0, 1]), with the crossing vertices recorded on the edge lists.
    #[test]
    fn perform_iso_lines_square() {
        let (_brep, face) = square_face();
        let bottom = face_outer_edge(&face, 0);
        let right = face_outer_edge(&face, 1);
        let top = face_outer_edge(&face, 2);
        let left = face_outer_edge(&face, 3);
        let mut ds = Data::new();

        perform(0, &face, &mut ds, 1);

        // One U iso + one V iso.
        let iso_l = ds.face_iso_l(&face);
        assert_eq!(iso_l.len(), 2);
        let e_u = &iso_l[0];
        let e_v = &iso_l[1];

        // The U iso: pcurve x = 0.5 (origin (0.5, 0), direction Y), range
        // [0, 1] (the domain endpoints U11 = 0, U21 = 1).
        assert!(e_u.is_edge());
        // The MakeShape FORWARD reset (cxx L478 vs L481 — the OCCT quirk).
        assert_eq!(e_u.orientation, Orientation::Forward);
        let (pc_u, t1, t2) = curve_on_surface(e_u, &face);
        assert_eq!(t1, 0.0);
        assert_eq!(t2, 1.0);
        match &pc_u {
            Curve2d::Line(l) => {
                assert!((l.origin.x - 0.5).abs() < 1e-12);
                assert!(l.origin.y.abs() < 1e-12);
                assert!(l.direction.distance(DVec2::Y) < 1e-12);
            }
            _ => panic!("the U iso pcurve is a line"),
        }
        // The iso edge vertices: FORWARD at (0.5, 0, 0), REVERSED at
        // (0.5, 1, 0).
        let (v_f, v_l) = top_exp_vertices(e_u);
        assert_eq!(v_f.orientation, Orientation::Forward);
        assert!(pnt_is_equal(brep_tool_pnt(&v_f), Point3::new(0.5, 0.0, 0.0), 1e-12));
        assert_eq!(v_l.orientation, Orientation::Reversed);
        assert!(pnt_is_equal(brep_tool_pnt(&v_l), Point3::new(0.5, 1.0, 0.0), 1e-12));

        // The V iso: pcurve y = 0.5 (origin (0, 0.5), direction X), range
        // [0, 1], vertices (0, 0.5, 0) / (1, 0.5, 0).
        assert!(e_v.is_edge());
        assert_eq!(e_v.orientation, Orientation::Forward);
        let (pc_v, t1, t2) = curve_on_surface(e_v, &face);
        assert_eq!(t1, 0.0);
        assert_eq!(t2, 1.0);
        match &pc_v {
            Curve2d::Line(l) => {
                assert!(l.origin.x.abs() < 1e-12);
                assert!((l.origin.y - 0.5).abs() < 1e-12);
                assert!(l.direction.distance(DVec2::X) < 1e-12);
            }
            _ => panic!("the V iso pcurve is a line"),
        }
        let (v_f, v_l) = top_exp_vertices(e_v);
        assert_eq!(v_f.orientation, Orientation::Forward);
        assert!(pnt_is_equal(brep_tool_pnt(&v_f), Point3::new(0.0, 0.5, 0.0), 1e-12));
        assert_eq!(v_l.orientation, Orientation::Reversed);
        assert!(pnt_is_equal(brep_tool_pnt(&v_l), Point3::new(1.0, 0.5, 0.0), 1e-12));

        // The Data vertex lists: one fresh crossing vertex per boundary edge
        // (MakeVertex found no equal end vertex and the lists were empty).
        let check = |ds: &mut Data, e: &Shape, par: f64| {
            ds.init_vertex(e);
            assert!(ds.more_vertex());
            assert!((ds.parameter() - par).abs() < 1e-12);
            ds.next_vertex();
            assert!(!ds.more_vertex());
        };
        check(&mut ds, &bottom, 0.5);
        check(&mut ds, &top, 0.5);
        check(&mut ds, &left, 0.5);
        check(&mut ds, &right, 0.5);
    }

    /// OCCT anchor: Perform with a pcurve whose own parameter range already
    /// equals the edge range — the AddElement(handle) overload runs (cxx
    /// L134-138) instead of the Geom2d_TrimmedCurve re-trim (cxx L141-143);
    /// the resulting isolines are the same.
    #[test]
    fn perform_add_element_both_overloads() {
        let (mut brep, face) = square_face();
        let e0 = face_outer_edge(&face, 0);
        let key = (face.ptr_id(), face.location);
        {
            let ed = brep.edge_mut_inplace(e0.clone());
            let (pc, _, _) = ed.pcurves.get(&key).expect("the bottom pcurve");
            let wrapped = Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(pc.clone()),
                t_min: 0.0,
                t_max: 1.0,
            });
            ed.pcurves.insert(key, (wrapped, 0.0, 1.0));
        }
        let mut ds = Data::new();
        perform(0, &face, &mut ds, 1);
        let iso_l = ds.face_iso_l(&face);
        assert_eq!(iso_l.len(), 2);
        // The U iso geometry is unchanged (crossings at the parameter 0.5).
        let (v_f, _) = top_exp_vertices(&iso_l[0]);
        assert!(pnt_is_equal(brep_tool_pnt(&v_f), Point3::new(0.5, 0.0, 0.0), 1e-12));
    }

    /// OCCT anchor: MakeVertex (cxx L419-463) — the end-vertex equality
    /// returns, the fresh INTERNAL vertex is appended, the ordered
    /// InsertBefore keeps the list sorted and the round-trips (point /
    /// parameter / tolerance) hold.
    #[test]
    fn make_vertex_round_trip() {
        let (_brep, face) = square_face();
        let e0 = face_outer_edge(&face, 0);
        let (v_first, _) = top_exp_vertices(&e0);
        let mut ds = Data::new();

        // (a) P equals the FORWARD end vertex: returned as-is, before any
        // InitVertex (the edge list stays empty).
        let v = make_vertex(&e0, brep_tool_pnt(&v_first), 0.25, 1e-7, &mut ds);
        assert!(v.is_same(&v_first));
        ds.init_vertex(&e0);
        assert!(!ds.more_vertex());

        // (b) a fresh interior vertex: appended INTERNAL with the tolerance.
        let v1 = make_vertex(&e0, Point3::new(0.4, 0.0, 0.0), 0.4, 1e-7, &mut ds);
        assert!(!v1.is_null());
        assert_eq!(v1.orientation, Orientation::Internal);
        assert!((brep_tool_pnt(&v1).x - 0.4).abs() < 1e-12);
        assert!((brep_tool_tolerance(&v1) - 1e-7).abs() < 1e-15);
        ds.init_vertex(&e0);
        assert!(ds.more_vertex());
        assert!(ds.vertex().is_same(&v1));
        assert!((ds.parameter() - 0.4).abs() < 1e-12);
        ds.next_vertex();
        assert!(!ds.more_vertex());

        // (c) Par 0.25 < curP 0.4: InsertBefore keeps the list ordered.
        let v2 = make_vertex(&e0, Point3::new(0.25, 0.0, 0.0), 0.25, 1e-7, &mut ds);
        ds.init_vertex(&e0);
        assert!(ds.vertex().is_same(&v2));
        assert!((ds.parameter() - 0.25).abs() < 1e-12);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v1));
        assert!((ds.parameter() - 0.4).abs() < 1e-12);
        ds.next_vertex();
        assert!(!ds.more_vertex());

        // (d) P equals the stored v2: returned, no new entry.
        let v3 = make_vertex(&e0, Point3::new(0.25, 0.0, 0.0), 0.9, 1e-7, &mut ds);
        assert!(v3.is_same(&v2));
        ds.init_vertex(&e0);
        let mut count = 0;
        while ds.more_vertex() {
            count += 1;
            ds.next_vertex();
        }
        assert_eq!(count, 2);

        // (e) Par 0.9 beyond every stored parameter: appended at the end.
        let v4 = make_vertex(&e0, Point3::new(0.9, 0.0, 0.0), 0.9, 1e-7, &mut ds);
        assert_eq!(v4.orientation, Orientation::Internal);
        ds.init_vertex(&e0);
        assert!(ds.vertex().is_same(&v2));
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v1));
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v4));
        assert!((ds.parameter() - 0.9).abs() < 1e-12);
        ds.next_vertex();
        assert!(!ds.more_vertex());
    }

    /// OCCT anchor: MakeIsoLine (cxx L467-488) — the vertex orientations are
    /// forced FORWARD / REVERSED, the built edge carries the iso pcurve with
    /// the (U1, U2) range from the UpdateVertex calls and lands in the face
    /// IsoL list; the INTERNAL edge orientation set at L478 does not survive
    /// MakeEdge (the MakeShape FORWARD reset).
    #[test]
    fn make_iso_line_anchor() {
        let (_brep, face) = square_face();
        let mut ds = Data::new();

        let iso = Curve2d::Line(Line2d::new(DVec2::new(0.5, 0.0), DVec2::Y));
        let mut v1 = brep_make_vertex(Point3::new(0.5, 0.0, 0.0), 1e-7);
        let mut v2 = brep_make_vertex(Point3::new(0.5, 1.0, 0.0), 1e-7);
        // Mimic the INTERNAL orientation the MakeVertex callers see before
        // the MakeIsoLine resets (cxx L448/L457 + L479/L480).
        v1.orientation = Orientation::Internal;
        v2.orientation = Orientation::Internal;

        make_iso_line(&face, &iso, &mut v1, &mut v2, 0.25, 0.75, 1e-7, &mut ds);

        // V1/V2 are re-oriented in place (TopoDS_Vertex& V1, V2).
        assert_eq!(v1.orientation, Orientation::Forward);
        assert_eq!(v2.orientation, Orientation::Reversed);

        // The iso edge landed in DS.AddIsoL(F).
        let iso_l = ds.face_iso_l(&face);
        assert_eq!(iso_l.len(), 1);
        let e = &iso_l[0];
        assert!(e.is_edge());
        // The MakeShape FORWARD reset (the L478 INTERNAL is lost).
        assert_eq!(e.orientation, Orientation::Forward);
        assert!(!e.is_null());

        // The pcurve round-trip: the iso line and the (U1, U2) range.
        let (pc, t1, t2) = curve_on_surface(e, &face);
        assert!((t1 - 0.25).abs() < 1e-12);
        assert!((t2 - 0.75).abs() < 1e-12);
        match &pc {
            Curve2d::Line(l) => {
                assert!((l.origin.x - 0.5).abs() < 1e-12);
                assert!(l.origin.y.abs() < 1e-12);
                assert!(l.direction.distance(DVec2::Y) < 1e-12);
            }
            _ => panic!("the iso pcurve is a line"),
        }

        // The vertices and their orientations on the edge.
        let (v_f, v_l) = top_exp_vertices(e);
        assert!(v_f.is_same(&v1));
        assert_eq!(v_f.orientation, Orientation::Forward);
        assert!(v_l.is_same(&v2));
        assert_eq!(v_l.orientation, Orientation::Reversed);

        // UpdateEdge / UpdateVertex tolerances.
        assert!((brep_tool_tolerance(e) - 1e-7).abs() < 1e-15);
    }

    /// OCCT anchor: Perform bounds on the infinite-direction face — the UV
    /// bounds of a pcurve-less face fall back to the natural surface bounds;
    /// with nbIsos = 0 the loops stay empty and no iso appears.
    #[test]
    fn perform_nb_isos_zero_is_noop() {
        let (_brep, face) = square_face();
        let mut ds = Data::new();
        perform(0, &face, &mut ds, 0);
        assert!(!ds.face_has_iso_l(&face));
    }

    /// The brep handle is unused by the anchors (the shapes are extracted
    /// through the TShape graph); keep the import alive.
    #[allow(dead_code)]
    fn _brep_type_check(_b: &BRep) {}
}
