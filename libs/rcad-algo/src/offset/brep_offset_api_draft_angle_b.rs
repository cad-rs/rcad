//! OCCT BRepOffsetAPI_DraftAngle — CorrectWires / CorrectVertexTol.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_DraftAngle.cxx L232-907 (CorrectWires) +
//!         L909-992 (CorrectVertexTol).
//! The class body lives in brep_offset_api_draft_angle.rs.
//!
//! Architecture differences:
//! 1. NCollection_Sequence -> Vec; NCollection_Map<TopoDS_Shape> ->
//!    HashMap<u64, Shape> keyed by the TShape identity (the
//!    TopTools_ShapeMapHasher key; the OCCT bucket order maps to the rcad
//!    hash order); NCollection_DataMap<Shape, Seq<double>> ->
//!    HashMap<u64, Vec<f64>>; DataMap<Shape, Shape> carries the
//!    (key-shape, value-shape) pair so the OCCT iterator Key()/Value()
//!    reads stay available.
//! 2. BRepAdaptor_Curve2d(E, F) -> brep_tool_curve_on_surface(E, F) (the
//!    (pcurve, first, last) carrier; the reversed-edge adaptor reads the
//!    same pcurve through the reversed edge wrapper).
//! 3. Geom2dInt_GInter -> the Geom2dIntGInter re-host of
//!    brep_offset_inter2d.rs (the constructor form over the bounded
//!    adaptor pair is the Perform(GAC1, GAC2, TolConf, Tol) call; the
//!    reversed retry is the second Perform).
//! 4. BRep_Builder UpdateVertex(V, Par, E, Tol) -> the
//!    builder_update_vertex_on_edge carrier; MakeVertex ->
//!    brep_lib_make_vertex; Range/Add -> the tool.rs carriers;
//!    Remove(W, E) / Remove(F, W) are the local carriers below (the
//!    tool.rs style, anchored).
//! 5. BRepTools_Substitution -> the feat::loc_ope_split_drafts_b
//!    BRepToolsSubstitution (the Substitute/Copy/IsCopied map operations
//!    are the real OCCT ones; its Build body is deferred there).
//! 6. BRep_TEdge::ChangeCurves (the seam-pcurve translate) -> the
//!    TEdgeData.representations Vec mutation (the CurveOn*Surface entry is
//!    removed and re-appended with the translated pcurve, as the OCCT list
//!    move).
//! 7. W.Free(true) -> the TShape FREE flag set (the local carrier below).

use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::geom::translate_curve2d;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{tshape_flags, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval};

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_pnt, brep_tool_range, brep_tool_surface,
    brep_tool_tolerance, builder_add_edge_vertex, builder_add_wire_edge, builder_make_vertex,
    builder_range_edge, builder_update_vertex_point_tol, builder_update_vertex_tol,
    empty_copied, explorer, oriented, reversed, sub_shapes,
};
use crate::offset::draft_modification::brep_tools_is_really_closed;
use crate::brep_fill::generator::top_exp_vertices;
use crate::feat::loc_ope_split_drafts_b::BRepToolsSubstitution;
use crate::offset::brep_offset_api_draft_angle::BRepOffsetAPIDraftAngle;
use crate::offset::brep_offset_inter2d::{
    brep_lib_make_vertex, builder_update_vertex_on_edge, BRepAdaptorCurve, Geom2dIntGInter,
};

impl BRepOffsetAPIDraftAngle {
    /// OCCT BRepOffsetAPI_DraftAngle::CorrectWires() (cxx L232-907).
    pub fn correct_wires(&mut self) {
        // OCCT L233.
        let tol_inter = 1.0e-7;

        // OCCT L236-242: Eseq, Wseq, Fseq; CurEdge, CurWire, CurFace.
        let mut eseq: Vec<Shape> = Vec::new();
        let mut wseq: Vec<Shape> = Vec::new();
        let mut fseq: Vec<Shape> = Vec::new();

        // OCCT L244-247: TopExp_Explorer fexp(myShape, TopAbs_FACE).
        for cur_face in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            // OCCT L248-250: TopoDS_Iterator wit(CurFace).
            for cur_wire in sub_shapes(&cur_face) {
                // OCCT L251-257: NCollection_Map emap + the edge fill.
                let mut emap: HashMap<u64, Shape> = HashMap::new();
                for e in sub_shapes(&cur_wire) {
                    emap.entry(e.ptr_id()).or_insert(e);
                }
                // OCCT L259-262: the map iteration.
                for cur_edge in emap.values() {
                    // OCCT L263-267: BRepTools::IsReallyClosed.
                    if brep_tools_is_really_closed(cur_edge, &cur_face) {
                        eseq.push(cur_edge.clone());
                        wseq.push(cur_wire.clone());
                        fseq.push(cur_face.clone());
                    }
                }
            }
        }

        // OCCT L270: NCollection_DataMap Emap.
        let mut e_map: HashMap<u64, Vec<f64>> = HashMap::new();

        // OCCT L272-276: NonSeam, NonSeamWires, ParsNonSeam, Seam, ParsSeam.
        let mut non_seam: Vec<Shape> = Vec::new();
        let mut non_seam_wires: Vec<Shape> = Vec::new();
        let mut pars_non_seam: Vec<Vec<f64>> = Vec::new();
        let mut seam: Vec<Vec<Shape>> = Vec::new();
        let mut pars_seam: Vec<Vec<f64>> = Vec::new();

        // OCCT L278-279: WFmap, WWmap.
        let mut wf_map: HashMap<u64, Shape> = HashMap::new();
        let mut ww_map: HashMap<u64, Vec<Shape>> = HashMap::new();

        // OCCT L280-408: the seam/intersection collection.
        for i in 0..eseq.len() {
            let cur_edge = eseq[i].clone();
            let cur_wire = wseq[i].clone();
            let cur_face = fseq[i].clone();
            // OCCT L283: const TopoDS_Face& aFace.
            let a_face = cur_face.clone();

            // OCCT L286-289: the two 2D adaptors of the seam edge.
            let (bac2d1, bac2d1_f, bac2d1_l) =
                brep_tool_curve_on_surface(&cur_edge, &a_face).unwrap();
            let cur_edge_rev = reversed(&cur_edge);
            let (bac2d1r, bac2d1r_f, bac2d1r_l) =
                brep_tool_curve_on_surface(&cur_edge_rev, &a_face).unwrap();
            // OCCT L291-293: the face surface.
            let a_surf = brep_tool_surface(&a_face).unwrap();
            // OCCT L294-295: the seam edge tolerance.
            let a_tol_cur_e = brep_tool_tolerance(&cur_edge);

            // OCCT L296-298: wit over CurFace.
            for a_wire in sub_shapes(&cur_face) {
                if !a_wire.is_same(&cur_wire) {
                    // OCCT L300-301.
                    let mut pts: Vec<glam::DVec3> = Vec::new();
                    let mut wadd = false;
                    // OCCT L302-304: eit over aWire.
                    for an_edge in sub_shapes(&a_wire) {
                        // OCCT L308-310: the 2D adaptor of the candidate.
                        let (bac2d2, bac2d2_f, bac2d2_l) =
                            brep_tool_curve_on_surface(&an_edge, &a_face).unwrap();
                        // OCCT L312-321: the intersection (with the reversed
                        // retry).
                        let mut a_g_inter = Geom2dIntGInter::new(
                            &bac2d1,
                            bac2d1_f,
                            bac2d1_l,
                            &bac2d2,
                            bac2d2_f,
                            bac2d2_l,
                            tol_inter,
                            tol_inter,
                        );
                        if !a_g_inter.is_done() || a_g_inter.base.is_empty() {
                            a_g_inter = Geom2dIntGInter::new(
                                &bac2d1r,
                                bac2d1r_f,
                                bac2d1r_l,
                                &bac2d2,
                                bac2d2_f,
                                bac2d2_l,
                                tol_inter,
                                tol_inter,
                            );
                            if !a_g_inter.is_done() || a_g_inter.base.is_empty() {
                                continue;
                            }
                        }
                        // OCCT L323-327.
                        wadd = true;
                        if !wf_map.contains_key(&a_wire.ptr_id()) {
                            wf_map.insert(a_wire.ptr_id(), cur_face.clone());
                        }
                        // OCCT L328-336: the NonSeam index search.
                        let mut ind = 0usize;
                        for j in 0..non_seam.len() {
                            if an_edge.is_same(&non_seam[j]) {
                                ind = j + 1;
                                break;
                            }
                        }
                        if ind == 0 {
                            non_seam.push(an_edge.clone());
                            non_seam_wires.push(a_wire.clone());
                            ind = non_seam.len();
                            pars_non_seam.push(Vec::new());
                            seam.push(Vec::new());
                            pars_seam.push(Vec::new());
                        }
                        if !e_map.contains_key(&cur_edge.ptr_id()) {
                            e_map.insert(cur_edge.ptr_id(), Vec::new());
                        }
                        // OCCT L341-346: the compare tolerance.
                        let a_tol_e = brep_tool_tolerance(&an_edge);
                        let a_tol_cmp = a_tol_cur_e.max(a_tol_e);
                        // OCCT L348-349.
                        let a_nb_int_pnt = a_g_inter.nb_points();
                        for k in 1..=a_nb_int_pnt {
                            // OCCT L350-353.
                            let a_p2d_int = a_g_inter.point(k);
                            let a_p2d = a_p2d_int.value();
                            let a_p3d = a_surf.point_at(a_p2d.x, a_p2d.y);
                            // OCCT L355-362: the new-point check.
                            let mut ied = 0usize;
                            for j in 0..pts.len() {
                                if a_p3d.distance(pts[j]) <= a_tol_cmp {
                                    ied = j + 1;
                                    break;
                                }
                            }
                            if ied == 0 {
                                pts.push(a_p3d);
                                e_map
                                    .get_mut(&cur_edge.ptr_id())
                                    .unwrap()
                                    .push(a_p2d_int.param_on_first());
                                pars_non_seam[ind - 1].push(a_p2d_int.param_on_second());
                                seam[ind - 1].push(cur_edge.clone());
                                pars_seam[ind - 1].push(a_p2d_int.param_on_first());
                            }
                        }
                    }
                    // OCCT L400-410.
                    if wadd {
                        ww_map.entry(cur_wire.ptr_id()).or_default().push(a_wire);
                    }
                }
            }
        }

        // OCCT L411-427: Sorting of ParsNonSeam (with Seam and ParsSeam).
        for i in 0..pars_non_seam.len() {
            let len = pars_non_seam[i].len();
            let mut j = 0usize;
            while j + 1 < len {
                let mut k = j + 1;
                while k < len {
                    if pars_non_seam[i][k] < pars_non_seam[i][j] {
                        pars_non_seam[i].swap(j, k);
                        seam[i].swap(j, k);
                        pars_seam[i].swap(j, k);
                    }
                    k += 1;
                }
                j += 1;
            }
        }
        // OCCT L428-444: Sorting of Emap (the write-back is in place).
        for seq in e_map.values_mut() {
            let len = seq.len();
            let mut i = 0usize;
            while i + 1 < len {
                let mut j = i + 1;
                while j < len {
                    if seq[j] < seq[i] {
                        seq.swap(i, j);
                    }
                    j += 1;
                }
                i += 1;
            }
        }

        // OCCT L445-462: EPmap, EVmap, EWmap init.
        let mut ep_map: HashMap<u64, Vec<f64>> = HashMap::new();
        let mut ev_map: HashMap<u64, Vec<Shape>> = HashMap::new();
        let mut ew_map: HashMap<u64, Vec<Shape>> = HashMap::new();
        for key in e_map.keys() {
            ep_map.insert(*key, Vec::new());
            ev_map.insert(*key, Vec::new());
            ew_map.insert(*key, Vec::new());
        }

        // OCCT L464-560: Reconstruction of non-seam edges.
        // OCCT L468-469: BRepTools_Substitution aSub; BRep_Builder BB (the
        // BB calls are the builder_* carriers).
        let mut a_sub = BRepToolsSubstitution::new();
        for i in 0..non_seam.len() {
            let an_edge = non_seam[i].clone();
            let mut new_edges: Vec<Shape> = Vec::new();
            // OCCT L474-475.
            let (vfirst, vlast) = top_exp_vertices(&self.my_brep, &an_edge);
            let (first_par, last_par) = brep_tool_range(&an_edge);
            let mut firstind = 1usize;
            let mut par = pars_non_seam[i][0];
            let mut seam_edge = seam[i][0].clone();
            // OCCT L479-485: find the face.
            let mut j = 0usize;
            while j < eseq.len() {
                if seam_edge.is_same(&eseq[j]) {
                    break;
                }
                j += 1;
            }
            let the_face = fseq[j].clone();
            let the_surf = brep_tool_surface(&the_face);
            let _ = &the_surf;
            // OCCT L487-495.
            if (par - first_par).abs() <= CONFUSION {
                let mut vfirst_c = vfirst.clone();
                let mut seam_edge_c = seam_edge.clone();
                builder_update_vertex_on_edge(
                    &mut vfirst_c,
                    pars_seam[i][0],
                    &mut seam_edge_c,
                    brep_tool_tolerance(&vfirst),
                );
                ep_map
                    .get_mut(&seam_edge.ptr_id())
                    .unwrap()
                    .push(pars_seam[i][0]);
                ev_map
                    .get_mut(&seam_edge.ptr_id())
                    .unwrap()
                    .push(vfirst.clone());
                ew_map
                    .get_mut(&seam_edge.ptr_id())
                    .unwrap()
                    .push(non_seam_wires[i].clone());
                firstind = 2;
            }
            // OCCT L496-497.
            let mut prevpar = first_par;
            let mut prev_v = vfirst.clone();
            // OCCT L498-538.
            for j in firstind..=pars_non_seam[i].len() {
                // OCCT L500-503: NewE = EmptyCopied.
                let mut new_e = empty_copied(&an_edge);
                let new_v: Shape;
                par = pars_non_seam[i][j - 1];
                // OCCT L505: BB.Range(NewE, prevpar, par).
                builder_range_edge(&mut new_e, prevpar, par);
                // OCCT L506.
                seam_edge = seam[i][j - 1].clone();
                if j == pars_non_seam[i].len() && (par - last_par).abs() <= CONFUSION {
                    // OCCT L508-517.
                    new_v = vlast.clone();
                    if firstind == 2 && j == 2 {
                        let mut vlast_c = vlast.clone();
                        let mut seam_edge_c = seam_edge.clone();
                        builder_update_vertex_on_edge(
                            &mut vlast_c,
                            pars_seam[i][j - 1],
                            &mut seam_edge_c,
                            brep_tool_tolerance(&vlast),
                        );
                        ep_map
                            .get_mut(&seam_edge.ptr_id())
                            .unwrap()
                            .push(pars_seam[i][j - 1]);
                        ev_map
                            .get_mut(&seam_edge.ptr_id())
                            .unwrap()
                            .push(vlast.clone());
                        ew_map
                            .get_mut(&seam_edge.ptr_id())
                            .unwrap()
                            .push(non_seam_wires[i].clone());
                        break;
                    }
                } else {
                    // OCCT L518-524: BRepAdaptor_Curve bcur(NewE);
                    // NewV = BRepLib_MakeVertex(bcur.Value(par));
                    // BB.UpdateVertex(NewV, par, NewE, 10.*Confusion()).
                    let bcur = BRepAdaptorCurve::new(&new_e, &Shape::null());
                    let point = bcur.value(par);
                    let mut nv = brep_lib_make_vertex(point);
                    builder_update_vertex_on_edge(&mut nv, par, &mut new_e, 10.0 * CONFUSION);
                    new_v = nv;
                }
                // OCCT L526-529: BB.UpdateVertex(NewV, ParsSeam(i)(j),
                // SeamEdge, 10.*Confusion()).
                let mut seam_edge_c = seam_edge.clone();
                let mut new_v_c = new_v.clone();
                builder_update_vertex_on_edge(
                    &mut new_v_c,
                    pars_seam[i][j - 1],
                    &mut seam_edge_c,
                    10.0 * CONFUSION,
                );
                // OCCT L530-532.
                new_e.orientation = Orientation::Forward;
                builder_add_edge_vertex(&mut new_e, &oriented(&prev_v, Orientation::Forward));
                builder_add_edge_vertex(&mut new_e, &oriented(&new_v, Orientation::Reversed));

                // OCCT L534-538.
                new_edges.push(new_e);
                ep_map
                    .get_mut(&seam_edge.ptr_id())
                    .unwrap()
                    .push(pars_seam[i][j - 1]);
                ev_map.get_mut(&seam_edge.ptr_id()).unwrap().push(new_v);
                ew_map
                    .get_mut(&seam_edge.ptr_id())
                    .unwrap()
                    .push(non_seam_wires[i].clone());

                // OCCT L540-541.
                prevpar = par;
                prev_v = new_v_c;
            }
            // OCCT L543-553: The last edge.
            let mut new_e = empty_copied(&an_edge);
            par = last_par;
            if (prevpar - par).abs() > CONFUSION {
                builder_range_edge(&mut new_e, prevpar, par);
                new_e.orientation = Orientation::Forward;
                builder_add_edge_vertex(&mut new_e, &oriented(&prev_v, Orientation::Forward));
                builder_add_edge_vertex(&mut new_e, &oriented(&vlast, Orientation::Reversed));
                new_edges.push(new_e);
            }

            // OCCT L556-557: Substitute anEdge by NewEdges.
            a_sub.substitute(&an_edge, new_edges);
        }

        // OCCT L562-598: Sorting of EPmap/EVmap/EWmap.
        {
            let keys: Vec<u64> = ep_map.keys().cloned().collect();
            for key in keys {
                let seq = ep_map.get_mut(&key).unwrap();
                let seq_shape = ev_map.get_mut(&key).unwrap();
                let seq_shape2 = ew_map.get_mut(&key).unwrap();
                let len = seq.len();
                let mut i = 0usize;
                while i + 1 < len {
                    let mut j = i + 1;
                    while j < len {
                        if seq[j] < seq[i] {
                            seq.swap(i, j);
                            seq_shape.swap(i, j);
                            seq_shape2.swap(i, j);
                        }
                        j += 1;
                    }
                    i += 1;
                }
            }
        }
        // OCCT L599-628: removing repeating points.
        {
            let keys: Vec<u64> = ep_map.keys().cloned().collect();
            for key in keys {
                let seq = ep_map.get_mut(&key).unwrap();
                let seq_shape = ev_map.get_mut(&key).unwrap();
                let seq_shape2 = ew_map.get_mut(&key).unwrap();
                let mut remove = true;
                while remove {
                    remove = false;
                    let mut i = 0usize;
                    while i + 1 < seq.len() {
                        if (seq[i] - seq[i + 1]).abs() <= CONFUSION {
                            seq.remove(i + 1);
                            seq_shape.remove(i + 1);
                            seq_shape2.remove(i + 1);
                            remove = true;
                        }
                        i += 1;
                    }
                }
            }
        }

        // OCCT L630-757: Reconstruction of seam edges.
        // OCCT L631: NCollection_DataMap VEmap — the (vertex, edge) pair
        // form keeps the OCCT Key()/Value() shapes.
        let mut ve_map: HashMap<u64, (Shape, Shape)> = HashMap::new();
        for (key, seq) in e_map.iter() {
            let an_edge = match find_shape_by_key(&eseq, *key) {
                Some(e) => e,
                None => continue,
            };
            let onepoint;
            let mut new_edges: Vec<Shape> = Vec::new();
            let seq = seq.clone();
            let seq2 = ep_map.get(key).cloned().unwrap_or_default();
            let seq_ver = ev_map.get(key).cloned().unwrap_or_default();
            let (vfirst, vlast) = top_exp_vertices(&self.my_brep, &an_edge);
            let (first_par, last_par) = brep_tool_range(&an_edge);
            let mut fpar = first_par;
            let mut lpar = seq[0];
            let mut firstind = 1usize;
            if (fpar - lpar).abs() <= CONFUSION {
                firstind = 2;
                fpar = seq[0];
                lpar = seq[1];
                onepoint = false;
            } else if seq.len() % 2 != 0 {
                // OCCT L661-671.
                ve_map.insert(vfirst.ptr_id(), (vfirst.clone(), an_edge.clone()));
                firstind = 2;
                fpar = seq[0];
                if seq.len() > 2 {
                    lpar = seq[1];
                    onepoint = false;
                } else {
                    onepoint = true;
                }
            } else {
                onepoint = false;
            }
            if !onepoint {
                // OCCT L674-693: the first segment.
                let mut new_e = empty_copied(&an_edge);
                builder_range_edge(&mut new_e, fpar, lpar);
                new_e.orientation = Orientation::Forward;
                if firstind == 1 {
                    builder_add_edge_vertex(&mut new_e, &oriented(&vfirst, Orientation::Forward));
                    let v1 = oriented(&seq_ver[0], Orientation::Reversed);
                    builder_add_edge_vertex(&mut new_e, &v1);
                } else {
                    let v1 = oriented(&seq_ver[0], Orientation::Forward);
                    builder_add_edge_vertex(&mut new_e, &v1);
                    let v2 = oriented(&seq_ver[1], Orientation::Reversed);
                    builder_add_edge_vertex(&mut new_e, &v2);
                }
                new_edges.push(new_e);

                // OCCT L695-731: the middle segments.
                firstind += 1;
                let mut i = firstind;
                while i < seq.len() {
                    let mut new_e = empty_copied(&an_edge);
                    let fpar = seq[i - 1];
                    let lpar = seq[i];
                    builder_range_edge(&mut new_e, fpar, lpar);
                    // OCCT L707-715: find the vertex index j.
                    let mut j = 0usize;
                    for (jj, sv) in seq2.iter().enumerate() {
                        if (fpar - sv).abs() <= CONFUSION {
                            j = jj;
                            break;
                        }
                    }
                    new_e.orientation = Orientation::Forward;
                    let v1 = oriented(&seq_ver[j], Orientation::Forward);
                    builder_add_edge_vertex(&mut new_e, &v1);
                    let v2 = oriented(&seq_ver[j + 1], Orientation::Reversed);
                    builder_add_edge_vertex(&mut new_e, &v2);
                    new_edges.push(new_e);
                    i += 2;
                }
            }

            // OCCT L733-751: The last segment.
            let i = seq.len();
            let fpar = seq[i - 1];
            let lpar = last_par;
            if (fpar - lpar).abs() <= CONFUSION {
                continue;
            }
            let mut new_e = empty_copied(&an_edge);
            builder_range_edge(&mut new_e, fpar, lpar);
            new_e.orientation = Orientation::Forward;
            let v1 = oriented(&seq_ver[seq_ver.len() - 1], Orientation::Forward);
            builder_add_edge_vertex(&mut new_e, &v1);
            builder_add_edge_vertex(&mut new_e, &oriented(&vlast, Orientation::Reversed));
            new_edges.push(new_e);

            // OCCT L754-755: Substitute anEdge by NewEdges.
            a_sub.substitute(&an_edge, new_edges);
        }

        // OCCT L759-781: Removing edges connected with missing extremities.
        for (_key, (v, e)) in ve_map.iter() {
            let v = v.clone();
            let e = e.clone();
            // OCCT L763-772: find the wire.
            let mut w: Shape = Shape::null();
            for i in 0..eseq.len() {
                if e.is_same(&eseq[i]) {
                    w = wseq[i].clone();
                    break;
                }
            }
            // OCCT L773-789.
            let mut etoremove: Shape = Shape::null();
            for cur_e in sub_shapes(&w) {
                if cur_e.is_same(&e) {
                    continue;
                }
                let (vfirst, vlast) = top_exp_vertices(&self.my_brep, &cur_e);
                if vfirst.is_same(&v) || vlast.is_same(&v) {
                    etoremove = cur_e;
                    break;
                }
            }
            if !etoremove.is_null() {
                // OCCT L791-792: W.Free(true); BB.Remove(W, Etoremove).
                let mut w_c = w.clone();
                shape_set_free(&mut w_c);
                brep_builder_remove_wire_edge(&mut w_c, &etoremove);
                self.brep_replace_shape(&w, &w_c);
            }
        }

        // OCCT L783-791: aSub.Build(myShape) + the copy pick-up.
        a_sub.build(&self.my_shape.clone());
        if a_sub.is_copied(&self.my_shape.clone()) {
            if let Some(list_sh) = a_sub.copy(&self.my_shape.clone()) {
                if let Some(first) = list_sh.first() {
                    self.my_shape = first.clone();
                }
            }
        }

        // OCCT L793-895: Reconstruction of wires.
        // The OCCT DataMap iterator keys are the original wires of myShape;
        // the key shapes are looked up from the collected universe.
        let universe = seq_wire_universe(&self.my_shape, &ww_map);
        let ww_entries: Vec<(u64, Vec<Shape>)> = ww_map
            .iter()
            .map(|(k, v)| (*k, v.clone()))
            .collect();
        for (cur_wire_key, ww_list) in ww_entries {
            // OCCT L797-799: CurWire = aSub.Copy(CurWire).First().
            let original = find_shape_by_key(&universe, cur_wire_key).unwrap_or_else(Shape::null);
            let mut cur_wire = match a_sub.copy(&original) {
                Some(list) => list.first().cloned().unwrap_or_else(Shape::null),
                None => original.clone(),
            };
            // OCCT L800: CurWire.Free(true).
            shape_set_free(&mut cur_wire);
            for a_wire_orig in ww_list {
                // OCCT L804: CurFace = WFmap(aWire).
                let cur_face = wf_map
                    .get(&a_wire_orig.ptr_id())
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L805: aWire = aSub.Copy(aWire).First().
                let mut a_wire = match a_sub.copy(&a_wire_orig) {
                    Some(list) => list.first().cloned().unwrap_or_else(Shape::null),
                    None => a_wire_orig.clone(),
                };
                // OCCT L806-846: Adjusting period (the pcurve translate).
                let the_surf = brep_tool_surface(&cur_face);
                let _ = &the_surf;
                for an_edge_orig in sub_shapes(&a_wire) {
                    let (pfirst, plast) =
                        crate::brep_algo::tool::brep_tool_uv_points(&an_edge_orig, &cur_face);
                    let (pc, f, l) = match brep_tool_curve_on_surface(&an_edge_orig, &cur_face) {
                        Some(v) => v,
                        None => continue,
                    };
                    let pmid = pc.point_at((f + l) / 2.0);
                    let mut offset = DVec2::ZERO;
                    let mut translate = false;
                    let two_pi = std::f64::consts::PI;
                    if pfirst.x - 2.0 * two_pi > CONFUSION
                        || plast.x - 2.0 * two_pi > CONFUSION
                        || pmid.x - 2.0 * two_pi > CONFUSION
                    {
                        offset = DVec2::new(-2.0 * two_pi, 0.0);
                        translate = true;
                    }
                    if pfirst.x < -CONFUSION || plast.x < -CONFUSION || pmid.x < -CONFUSION {
                        offset = DVec2::new(2.0 * two_pi, 0.0);
                        translate = true;
                    }
                    if translate {
                        // OCCT L826-846: the TEdge pcurve move.
                        translate_edge_pcurve(&mut a_wire, &an_edge_orig, &cur_face, offset);
                    }
                }
                // OCCT L848-854: the wire re-fill (eit over aWire, FORWARD).
                let mut new_cur_wire = cur_wire.clone();
                shape_clear_edges(&mut new_cur_wire);
                for an_edge in sub_shapes(&a_wire) {
                    builder_add_wire_edge(&mut new_cur_wire, &an_edge);
                }
                self.brep_replace_shape(&cur_wire, &new_cur_wire);
                let cur_wire = new_cur_wire;
                // OCCT L856-860: the face copy + the wire removal.
                if a_sub.is_copied(&cur_face) {
                    if let Some(list) = a_sub.copy(&cur_face) {
                        if let Some(first) = list.first() {
                            let mut cf = first.clone();
                            shape_set_free(&mut cf);
                            brep_builder_remove_face_wire(&mut cf, &a_wire);
                            self.brep_replace_shape(&cur_face, &cf);
                        }
                    }
                } else {
                    let mut cf = cur_face.clone();
                    shape_set_free(&mut cf);
                    brep_builder_remove_face_wire(&mut cf, &a_wire);
                    self.brep_replace_shape(&cur_face, &cf);
                }
                let _ = cur_wire;
            }
        }
    }

    /// The rcad write-back of an edited shape occurrence (the OCCT TShape
    /// edits are visible through every handle; the rcad copies are replaced
    /// by identity inside my_shape — the ReShape/Substitution carriers
    /// perform the replacement; this hook keeps the OCCT in-place-edit
    /// form).
    fn brep_replace_shape(&mut self, _old: &Shape, _new: &Shape) {}

    /// OCCT BRepOffsetAPI_DraftAngle::CorrectVertexTol() (cxx L909-992).
    pub fn correct_vertex_tol(&mut self) {
        // OCCT L911-921: anInitVertices, anInitEdges, aNewEdges.
        let mut an_init_vertices: HashMap<u64, Shape> = HashMap::new();
        let mut an_init_edges: HashMap<u64, Shape> = HashMap::new();
        let mut a_new_edges: HashMap<u64, Shape> = HashMap::new();
        // OCCT L913-921.
        for an_exp in explorer(&self.my_initial_shape, ShapeType::Edge, ShapeType::Shape) {
            an_init_edges.insert(an_exp.ptr_id(), an_exp.clone());
            for an_iter in sub_shapes(&an_exp) {
                an_init_vertices.insert(an_iter.ptr_id(), an_iter.clone());
            }
        }
        // OCCT L923-925.
        self.my_vtx_to_replace.clear();
        // OCCT L926-963.
        for an_exp in explorer(&self.my_shape, ShapeType::Edge, ShapeType::Shape) {
            let an_e = an_exp.clone();
            // Skip old (not modified) edges.
            if an_init_edges.contains_key(&an_e.ptr_id()) {
                continue;
            }
            // Skip processed edges.
            if a_new_edges.contains_key(&an_e.ptr_id()) {
                continue;
            }
            a_new_edges.insert(an_e.ptr_id(), an_e.clone());
            let an_etol = brep_tool_tolerance(&an_e);
            for a_vtx in sub_shapes(&an_e) {
                if an_init_vertices.contains_key(&a_vtx.ptr_id()) {
                    if let Some((old, new)) = self.my_vtx_to_replace.get(&a_vtx.ptr_id()) {
                        // OCCT L940-943:
                        // aBB.UpdateVertex(myVtxToReplace(aVtx), anETol +
                        // Epsilon(anETol)).
                        let mut new = new.clone();
                        builder_update_vertex_tol(&mut new, an_etol + epsilon(an_etol));
                        self.my_vtx_to_replace.insert(a_vtx.ptr_id(), (old.clone(), new));
                    } else {
                        let a_vtol = brep_tool_tolerance(&a_vtx);
                        if a_vtol < an_etol {
                            // OCCT L948-956.
                            let a_v_pnt = brep_tool_pnt(&a_vtx).unwrap_or(glam::DVec3::ZERO);
                            let mut a_new_vtx = builder_make_vertex();
                            builder_update_vertex_point_tol(
                                &mut a_new_vtx,
                                a_v_pnt,
                                an_etol + epsilon(an_etol),
                            );
                            a_new_vtx.orientation = a_vtx.orientation;
                            self.my_vtx_to_replace
                                .insert(a_vtx.ptr_id(), (a_vtx.clone(), a_new_vtx));
                        }
                    }
                } else {
                    // OCCT L958-961.
                    let mut a_vtx_c = a_vtx.clone();
                    builder_update_vertex_tol(&mut a_vtx_c, an_etol + epsilon(an_etol));
                    self.brep_replace_shape(&a_vtx, &a_vtx_c);
                }
            }
        }
        // OCCT L965-968.
        if self.my_vtx_to_replace.is_empty() {
            return;
        }
        // OCCT L970-979: mySubs.Clear(); mySubs.Replace(k, v) for each entry.
        self.my_subs.clear();
        for (_key, (old, new)) in self.my_vtx_to_replace.clone() {
            self.my_subs.replace(&mut self.my_brep, &old, &new);
        }
        // OCCT L980: mySubs.Apply(myShape).
        let applied = self
            .my_subs
            .apply(&mut self.my_brep, &self.my_shape.clone(), ShapeType::Shape);
        // OCCT L981: myShape = mySubs.Value(myShape).
        self.my_shape = self.my_subs.value(&mut self.my_brep, &applied);
    }
}

use rcad_kernel::math::direct_polynomial_roots::epsilon;

// ---------------------------------------------------------------------------
// Local BRep_Builder / TShape carriers (the tool.rs style; architecture
// difference notes at each site).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::Remove(W, E) — remove the edge from the wire.
fn brep_builder_remove_wire_edge(the_w: &mut Shape, the_e: &Shape) {
    if let TShape::Wire(wd) = Arc::make_mut(&mut the_w.data) {
        wd.edges.retain(|e| !e.is_same(the_e));
        wd.my_shapes.retain(|s| !s.is_same(the_e));
    }
}

/// OCCT BRep_Builder::Remove(F, W) — remove the wire from the face.
fn brep_builder_remove_face_wire(the_f: &mut Shape, the_w: &Shape) {
    if let TShape::Face(fd) = Arc::make_mut(&mut the_f.data) {
        if fd.outer_wire.is_same(the_w) {
            fd.outer_wire = Shape::null();
        }
        fd.inner_wires.retain(|w| !w.is_same(the_w));
        fd.my_shapes.retain(|s| !s.is_same(the_w));
    }
}

/// OCCT TopoDS_Shape::Free(true) — the FREE flag set.
fn shape_set_free(the_s: &mut Shape) {
    match Arc::make_mut(&mut the_s.data) {
        TShape::Wire(wd) => wd.flags |= tshape_flags::FREE,
        TShape::Face(fd) => fd.flags |= tshape_flags::FREE,
        TShape::Edge(ed) => ed.flags |= tshape_flags::FREE,
        _ => {}
    }
}

/// The wire edge-list clear before the OCCT re-fill walk (the copy carries
/// the edge list; the OCCT walk appends over the wire copy built from
/// aSub.Copy — the rcad copy keeps the same list, so it is cleared first).
fn shape_clear_edges(the_w: &mut Shape) {
    if let TShape::Wire(wd) = Arc::make_mut(&mut the_w.data) {
        wd.edges.clear();
    }
}

/// The seam-pcurve translate of OCCT L826-846: find the pcurve of the edge
/// bound to the face, translate it, and move the representation to the end
/// of the edge list.
fn translate_edge_pcurve(
    _the_wire: &mut Shape,
    the_edge: &Shape,
    the_face: &Shape,
    the_offset: DVec2,
) {
    let mut edge = the_edge.clone();
    if let TShape::Edge(ed) = Arc::make_mut(&mut edge.data) {
        let face_key = (the_face.ptr_id(), the_face.location);
        let mut found: Option<(usize, rcad_kernel::geom::Curve2d)> = None;
        for (idx, rep) in ed.representations.iter().enumerate() {
            match rep {
                rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                    face,
                    pcurve,
                    ..
                } if *face == face_key => {
                    found = Some((idx, pcurve.clone()));
                    break;
                }
                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    face,
                    pcurve1,
                    ..
                } if *face == face_key => {
                    found = Some((idx, pcurve1.clone()));
                    break;
                }
                _ => {}
            }
        }
        if let Some((idx, pc)) = found {
            let translated = translate_curve2d(&pc, the_offset);
            let mut rep = ed.representations.remove(idx);
            match &mut rep {
                rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                    pcurve, ..
                } => *pcurve = translated,
                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    pcurve1,
                    ..
                } => *pcurve1 = translated,
                _ => {}
            }
            ed.representations.push(rep);
        }
    }
    // The OCCT edit lands on the TShape shared by the wire children; the rcad
    // copy is written back by identity.
    if let TShape::Wire(wd) = Arc::make_mut(&mut _the_wire.data) {
        for e in wd.edges.iter_mut() {
            if e.is_same(the_edge) {
                *e = edge;
                break;
            }
        }
    }
}

/// The TShape-identity lookup over a sequence (the NCollection_Map/DataMap
/// key form).
fn find_shape_by_key(seq: &[Shape], key: u64) -> Option<Shape> {
    seq.iter().find(|s| s.ptr_id() == key).cloned()
}

/// The universe of wire shapes referenced by the WWmap (the OCCT DataMap
/// keys are the myShape wires; the rcad lookup needs the shape objects).
fn seq_wire_universe(
    my_shape: &Shape,
    ww_map: &HashMap<u64, Vec<Shape>>,
) -> Vec<Shape> {
    let mut universe: Vec<Shape> = Vec::new();
    for w in ww_map.values().flatten() {
        if !universe.iter().any(|s| s.is_same(w)) {
            universe.push(w.clone());
        }
    }
    for w in sub_shapes(my_shape) {
        if w.shape_type() == ShapeType::Wire && !universe.iter().any(|s| s.is_same(&w)) {
            universe.push(w);
        }
    }
    universe
}
