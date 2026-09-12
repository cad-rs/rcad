// OCCT BRepOffset_MakeOffset.cxx — 1:1 translation, module d (see the
// module header of brep_offset_make_offset.rs for the split map and the
// architecture-difference list #38-#56).
//
// Module d carries the remaining file statics: UpdateInitOffset (cxx
// L3110-3146), ComputeMaxDist (L4167-4194), UpdateTolerance (L4196-4302),
// CorrectSolid (L4304-4382), checkSinglePoint (L4583-4626), RemoveShapes
// (L4628-4647), UpdateHistory (L4649-4672), TrimEdges (L4752-4924), TrimEdge
// (L4926-5034), GetEnlargedFaces (L5036-5062), BuildShellsCompleteInter
// (L5064-5229), GetSubShapes (L5386-5408),
// RemoveSeamAndDegeneratedEdges (L5550-5633), IsSolid (L5637-5645),
// AppendToList (L5647-5659).

use std::collections::HashMap;
use std::sync::Arc;

use glam::DVec3;
use rcad_kernel::geom::{CurveEval, Curve3, Surface3};
use rcad_kernel::topo::topods::{Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_inter2d_b::BRepOffsetInter2d;
use super::brep_offset_make_offset::{
    brep_check_edge_tolerance, brep_check_vertex_tolerance, brep_gprop_volume_properties,
    brep_lib_same_parameter_3, brep_tools_is_really_closed, bop_algo_tools_make_split_edge,
    find_parameter, top_exp_vertices_cum_ori, BRepOffset_Error,
    BOPAlgoMakerVolume, DataMapOfShapeListOfShape, DataMapOfShapeShape,
    IndexedDataMapOfShapeListOfShape, MapSF,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, top_exp_vertices,
    OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap,
};
use super::brep_offset_tool_d::BRepOffsetAnalyse;
use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool as bat;
use crate::brep_algo::tool::{brep_tool_pnt, shape_key};
use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve, brep_tool_degenerated, brep_tool_range, brep_tool_tolerance,
};
use rcad_kernel::core::precision::{
    CONFUSION as PRECISION_CONFUSION, SQUARE_CONFUSION as PRECISION_SQUARE_CONFUSION,
};

/// OCCT UpdateInitOffset (cxx L3110-3146) — Update and cleaning of
/// myInitOffset.
pub(crate) fn update_init_offset(
    my_init_offset: &mut BRepAlgoImage,
    my_image_offset: &BRepAlgoImage,
    my_offset_shape: &Shape,
    the_shape_type: ShapeType,
) {
    let mut niof = BRepAlgoImage::new();
    let roots = my_init_offset.roots().to_vec();
    for it in &roots {
        niof.set_root(it);
    }
    for it in &roots {
        let si = it;
        let mut li: Vec<Shape> = Vec::new();
        let mut l1: Vec<Shape> = Vec::new();
        my_init_offset.last_image(si, &mut l1);
        for it_l1 in &l1 {
            let o1 = it_l1;
            let mut l2: Vec<Shape> = Vec::new();
            my_image_offset.last_image(o1, &mut l2);
            li.extend(l2);
        }
        niof.bind_list(si, &li);
    }
    //  Modified by skv - Mon Apr  4 18:17:27 2005 Begin
    //  Supporting history.
    //   NIOF.Filter(myOffsetShape,TopAbs_FACE);
    niof.filter(my_offset_shape, the_shape_type);
    //  Modified by skv - Mon Apr  4 18:17:27 2005 End
    *my_init_offset = niof;
}

/// OCCT ComputeMaxDist (cxx L4167-4194).
pub(crate) fn compute_max_dist(
    the_plane: &rcad_kernel::geom::Plane,
    the_crv: &Curve3,
    the_first: f64,
    the_last: f64,
) -> f64 {
    let mut a_max_dist = 0.;
    let ncontrol: i32 = 23;
    for i in 0..ncontrol {
        let a_prm = ((ncontrol - 1 - i) as f64 * the_first + i as f64 * the_last)
            / (ncontrol - 1) as f64;
        let a_p = the_crv.point_at(a_prm);
        if a_p.x.is_infinite() || a_p.y.is_infinite() || a_p.z.is_infinite() {
            return f64::INFINITY;
        }
        // OCCT L4186: thePlane.SquareDistance(aP) — the gp_Pln signed-free
        // point-plane distance (the gp_Math form).
        let a_dist2 = gp_pln_square_distance(the_plane, a_p);
        if a_dist2 > a_max_dist {
            a_max_dist = a_dist2;
        }
    }
    a_max_dist.sqrt() * 1.05
}

/// OCCT gp_Pln::SquareDistance(P) (TKMath/gp) — the squared point-plane
/// distance (the plane normal is unit by the gp invariant).
fn gp_pln_square_distance(the_plane: &rcad_kernel::geom::Plane, p: DVec3) -> f64 {
    let d = (p - the_plane.origin).dot(the_plane.normal);
    d * d
}

/// OCCT UpdateTolerance (cxx L4196-4302).
pub(crate) fn update_tolerance(
    s: &mut Shape,
    faces: &OcctIndexedShapeMap,
    the_init_shape: &Shape,
) {
    let mut view: OcctShapeSet = HashMap::new();

    // The edges of caps are not modified.
    for j in 1..=faces.extent() {
        let f = faces.at_1(j).clone();
        for exp in bat::explorer(&f, ShapeType::Edge, ShapeType::Shape) {
            set_add(&mut view, &exp);
        }
    }

    // The edges of initial shape are  not modified
    let mut a_map_init_f: OcctShapeSet = HashMap::new();
    if !the_init_shape.is_null() {
        for an_exp_f in bat::explorer(the_init_shape, ShapeType::Face, ShapeType::Shape) {
            set_add(&mut a_map_init_f, &an_exp_f);
            for an_exp_e in bat::explorer(&an_exp_f, ShapeType::Edge, ShapeType::Shape) {
                set_add(&mut view, &an_exp_e);
                // OCCT L4232-4235: TopoDS_Iterator anItV(anExpE.Current()).
                for an_it_v in bat::sub_shapes(&an_exp_e) {
                    set_add(&mut view, &an_it_v);
                }
            }
        }
    }

    let mut tol;
    let an_exp_f = bat::explorer(s, ShapeType::Face, ShapeType::Shape);
    for f in &an_exp_f {
        if faces.contains(f) || set_contains(&a_map_init_f, f) {
            continue;
        }
        // OCCT L4252: BRepAdaptor_Surface aBAS(TopoDS::Face(F), false) — the
        // face surface carrier (architecture difference #54).
        let a_bas_surface = bat::brep_tool_surface(f);
        for exp in bat::explorer(f, ShapeType::Edge, ShapeType::Shape) {
            let mut e = exp;
            let mut is_updated = false;
            let a_curr_tol = brep_tool_tolerance(&e);
            if let Some(rcad_kernel::geom::Surface3::Plane(pln)) = &a_bas_surface {
                // Edge does not seem to have pcurve on plane,
                // so EdgeCorrector does not include it in tolerance calculation
                // OCCT L4263-4268.
                let (a_first, a_last) = brep_tool_range(&e);
                if let Some((a_crv, _, _)) = brep_tool_curve(&e) {
                    let a_max_dist = compute_max_dist(pln, &a_crv, a_first, a_last);
                    if a_max_dist > a_curr_tol {
                        update_edge_tolerance_host(&mut e, a_max_dist);
                        is_updated = true;
                    }
                }
            }
            if set_add(&mut view, &e) {
                // OCCT L4277: E.Locked(false).
                set_locked(&mut e, false);
                // OCCT L4278-4279: BRepCheck_Edge EdgeCorrector(E); Tol =
                // EdgeCorrector.Tolerance() (GAP leaf, arch. diff. #49).
                tol = brep_check_edge_tolerance(&e);
                if tol > a_curr_tol {
                    update_edge_tolerance_host(&mut e, tol);
                    is_updated = true;
                }
            }
            if is_updated {
                tol = brep_tool_tolerance(&e);
                // Update the vertices.
                // OCCT L4288: TopExp::Vertices(E, V[0], V[1]).
                let (v0, v1) = top_exp_vertices(&e);
                let vs: [Shape; 2] = [v0, v1];

                for v in &vs {
                    let mut v = v.clone();
                    // OCCT L4294: V[i].Locked(false).
                    set_locked(&mut v, false);
                    if set_add(&mut view, &v) {
                        // OCCT L4297-4302: TV->Tolerance(0.);
                        // BRepCheck_Vertex VertexCorrector(V[i]);
                        // B.UpdateVertex(V[i], VertexCorrector.Tolerance());
                        // (TV->ChangePoints()).Clear() — the rcad
                        // TShape::Vertex make_mut edits (arch. diff. #22).
                        set_vertex_tolerance(&mut v, 0.);
                        let vertex_tol = brep_check_vertex_tolerance(&v);
                        update_vertex_tolerance_host(&mut v, vertex_tol);
                        clear_vertex_points(&mut v);
                    }
                    update_vertex_tolerance_host(&mut v, tol);
                }
            }
        }
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol) — the rcad edge-tolerance edit
/// (the TShape::Edge make_mut form, arch. diff. #22).
fn update_edge_tolerance_host(the_e: &mut Shape, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, Tol) — the rcad vertex-tolerance edit
/// (the TShape::Vertex make_mut form, arch. diff. #22).
fn update_vertex_tolerance_host(the_v: &mut Shape, the_tol: f64) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
}

/// OCCT TopoDS_Shape::Locked(flag) — the rcad TShape flag edit.
fn set_locked(the_s: &mut Shape, the_flag: bool) {
    let d = Arc::make_mut(&mut the_s.data);
    match d {
        TShape::Vertex(vd) => set_flags_locked(&mut vd.flags, the_flag),
        TShape::Edge(ed) => set_flags_locked(&mut ed.flags, the_flag),
        TShape::Wire(wd) => set_flags_locked(&mut wd.flags, the_flag),
        TShape::Face(fd) => set_flags_locked(&mut fd.flags, the_flag),
        TShape::Shell(sd) => set_flags_locked(&mut sd.flags, the_flag),
        TShape::Solid(sd) => set_flags_locked(&mut sd.flags, the_flag),
        _ => {}
    }
}

fn set_flags_locked(flags: &mut u16, the_flag: bool) {
    if the_flag {
        *flags |= rcad_kernel::topo::topods::tshape_flags::LOCKED;
    } else {
        *flags &= !rcad_kernel::topo::topods::tshape_flags::LOCKED;
    }
}

/// OCCT Handle(BRep_TVertex)::Tolerance(0.) — the direct set form (the
/// UpdateTolerance reset before the VertexCorrector probe).
fn set_vertex_tolerance(the_v: &mut Shape, the_tol: f64) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.tolerance = the_tol;
    }
}

/// OCCT (TV->ChangePoints()).Clear().
fn clear_vertex_points(the_v: &mut Shape) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.points.clear();
    }
}

/// OCCT CorrectSolid (cxx L4304-4382).
pub(crate) fn correct_solid(
    brep: &mut rcad_kernel::topo::topods::BRep,
    the_sol: &mut Shape,
    the_sol_list: &mut Vec<Shape>,
) {
    let mut a_bb = rcad_kernel::topo::topods::BRepBuilder::new();
    let mut an_outer_shell = Shape::null();
    let mut a_vols: Vec<f64> = Vec::new();
    let mut a_vol_max = 0.;
    let mut an_outer_vol = 0.;

    for an_it in bat::sub_shapes(the_sol) {
        let a_sh = an_it;
        // OCCT L4322: BRepGProp::VolumeProperties(aSh, aVProps, true)
        // (GAP leaf, arch. diff. #51).
        let a_mass = brep_gprop_volume_properties(&a_sh);
        if a_mass.abs() > a_vol_max {
            an_outer_vol = a_mass;
            a_vol_max = an_outer_vol.abs();
            an_outer_shell = a_sh.clone();
        }
        a_vols.push(a_mass);
    }
    //
    if an_outer_vol.abs() < PRECISION_CONFUSION {
        return;
    }
    if an_outer_vol < 0. {
        // OCCT L4342: anOuterShell.Reverse().
        an_outer_shell = bat::reversed(&an_outer_shell);
    }
    // OCCT L4344-4346: TopoDS_Solid aNewSol; aBB.MakeSolid(aNewSol);
    // aNewSol.Closed(true); aBB.Add(aNewSol, anOuterShell) — the rcad
    // shells-Vec form (module-b arch. note).
    let mut a_new_sol = a_bb.make_solid(brep, vec![an_outer_shell.clone()]);
    bat::builder_set_closed(&mut a_new_sol, true);
    // OCCT L4348: BRepClass3d_SolidClassifier aSolClass(aNewSol) — the
    // crate::topalgo::brep_class3d translation (arch. diff. #50).
    let mut a_sol_class =
        crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(&a_new_sol);
    //
    for an_it in bat::sub_shapes(the_sol) {
        let a_v = a_vols.remove(0);
        let mut a_sh = an_it;
        if a_sh.is_same(&an_outer_shell) {
            continue;
        } else {
            // OCCT L4356-4358: TopExp_Explorer aVExp(aSh, TopAbs_VERTEX) —
            // the first vertex.
            let a_vtx = bat::explorer(&a_sh, ShapeType::Vertex, ShapeType::Shape)
                .into_iter()
                .next()
                .expect("CorrectSolid: shell without vertex");
            let a_p = brep_tool_pnt(&a_vtx).expect("BRep_Tool::Pnt null");
            a_sol_class.perform(a_p, brep_tool_tolerance(&a_vtx));
            // OCCT TopAbs_IN — the classifier state encoding (IN = 0).
            if a_sol_class.state() == 0 {
                if a_v > 0. {
                    a_sh = bat::reversed(&a_sh);
                }
                add_to_solid_host(brep, &mut a_new_sol, &a_sh);
            } else {
                if a_v < 0. {
                    a_sh = bat::reversed(&a_sh);
                }
                let mut a_sol = a_bb.make_solid(brep, vec![a_sh.clone()]);
                bat::builder_set_closed(&mut a_sol, true);
                the_sol_list.push(a_sol);
            }
        }
    }
    *the_sol = a_new_sol;
}

/// OCCT BRep_Builder::Add(Solid, Shell) — the rcad TShape::Solid make_mut
/// edit (arch. diff. #22; the incremental Add form of CorrectSolid).
fn add_to_solid_host(
    _brep: &mut rcad_kernel::topo::topods::BRep,
    the_sol: &mut Shape,
    the_sh: &Shape,
) {
    if let TShape::Solid(sd) = Arc::make_mut(&mut the_sol.data) {
        sd.shells.push(the_sh.clone());
    }
}

/// OCCT checkSinglePoint (cxx L4583-4626).
pub(crate) fn check_single_point(
    the_u_param: f64,
    the_v_param: f64,
    the_surf: &Surface3,
    the_bad_points: &[DVec3],
) -> BRepOffset_Error {
    use rcad_kernel::geom::SurfaceEval;
    let (a_pnt, a_d1u, a_d1v) = the_surf.derivatives(the_u_param, the_v_param);

    if a_d1u.length_squared() < PRECISION_SQUARE_CONFUSION
        || a_d1v.length_squared() < PRECISION_SQUARE_CONFUSION
    {
        let mut is_known_bad_pnt = false;
        for an_idx in 0..the_bad_points.len() {
            if a_pnt.distance_squared(the_bad_points[an_idx]) < PRECISION_SQUARE_CONFUSION {
                is_known_bad_pnt = true;
                break;
            }
        } // for(int anIdx  = theBadPoints.Lower();

        if !is_known_bad_pnt {
            return BRepOffset_Error::BadNormalsOnGeometry;
        } else {
            return BRepOffset_Error::NoError;
        }
    } //  if (aD1U.SquareMagnitude() < Precision::SquareConfusion() ||

    // OCCT L4618: aD1U.IsParallel(aD1V, Precision::Confusion()) — the
    // gp_Dir angular-parallel form.
    let an_angle = a_d1u.angle_between(a_d1v);
    if an_angle <= PRECISION_CONFUSION
        || (std::f64::consts::PI - an_angle) <= PRECISION_CONFUSION
    {
        // Isolines are collinear.
        return BRepOffset_Error::BadNormalsOnGeometry;
    }

    return BRepOffset_Error::NoError;
}

/// OCCT RemoveShapes (cxx L4628-4647).
pub(crate) fn remove_shapes(the_s: &mut Shape, the_ls: &[Shape]) {
    //
    let b_free = shape_is_free(the_s);
    bat::builder_set_free(the_s, true);
    //
    for a_it in the_ls {
        let a_si = a_it;
        remove_from_shape_host(the_s, a_si);
    }
    //
    bat::builder_set_free(the_s, b_free);
}

/// OCCT TopoDS_Shape::Free() — the rcad TShape flag read.
fn shape_is_free(the_s: &Shape) -> bool {
    match the_s.data.as_ref() {
        TShape::Vertex(vd) => vd.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        TShape::Edge(ed) => ed.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        TShape::Wire(wd) => wd.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        TShape::Face(fd) => fd.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        TShape::Shell(sd) => sd.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        TShape::Solid(sd) => sd.flags & rcad_kernel::topo::topods::tshape_flags::FREE != 0,
        _ => false,
    }
}

/// OCCT UpdateHistory (cxx L4649-4672).
pub(crate) fn update_history(
    the_lf: &[Shape],
    the_gf: &BOPAlgoMakerVolume,
    the_image: &mut BRepAlgoImage,
) {
    for a_it in the_lf {
        let a_f = a_it;
        let a_lf_im = the_gf.modified(a_f);
        if !a_lf_im.is_empty() {
            if the_image.has_image(a_f) {
                the_image.add_list(a_f, &a_lf_im);
            } else {
                the_image.bind_list(a_f, &a_lf_im);
            }
        }
    }
}

/// OCCT TrimEdges (cxx L4752-4924).
#[allow(clippy::too_many_arguments)]
pub(crate) fn trim_edges(
    the_brep: &rcad_kernel::topods::BRep,
    the_shape: &Shape,
    the_offset: f64,
    analyse: &BRepOffsetAnalyse,
    the_map_sf: &MapSF,
    the_mes: &mut DataMapOfShapeShape,
    the_build: &DataMapOfShapeShape,
    the_as_des: &mut BRepAlgoAsDes,
    the_as_des2d: &mut BRepAlgoAsDes,
    the_new_edges: &mut OcctIndexedShapeMap,
    the_e_trim_e_inf: &mut DataMapOfShapeShape,
    the_edges_origins: &mut DataMapOfShapeListOfShape,
) -> bool {
    let mut ne = Shape::null();

    let mut a_lfaces: Vec<Shape> = Vec::new();
    for exp in bat::explorer(the_shape, ShapeType::Face, ShapeType::Shape) {
        a_lfaces.push(exp);
    }

    let mut a_mf_generated: OcctShapeSet = HashMap::new();
    let mut a_dmef: IndexedDataMapOfShapeListOfShape = indexmap::IndexMap::new();
    for it in analyse.new_faces() {
        let a_fg = it;
        a_lfaces.push(a_fg.clone());
        set_add(&mut a_mf_generated, &a_fg);
        // OCCT L4789: TopExp::MapShapesAndUniqueAncestors(aFG, EDGE, FACE,
        // aDMEF).
        super::brep_offset_make_offset::top_exp_map_shapes_and_unique_ancestors(
            &a_fg,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_dmef,
        );
    }

    for it in a_lfaces.clone() {
        let fi = it;
        let mut nf = shape_data_map::value(the_map_sf, &fi).face();
        if shape_data_map::is_bound(the_mes, &nf) {
            nf = shape_data_map::find(the_mes, &nf);
        }
        //
        let mut view: OcctShapeSet = HashMap::new();
        let mut ve_map = OcctIndexedShapeMap::new();
        //
        // OCCT L4807-4808: TopExp::MapShapes(FI(FORWARD), EDGE/VERTEX,
        // VEmap).
        let fi_forward = bat::oriented(&fi, Orientation::Forward);
        for s in bat::explorer(&fi_forward, ShapeType::Edge, ShapeType::Shape) {
            ve_map.add(&s);
        }
        for s in bat::explorer(&fi_forward, ShapeType::Vertex, ShapeType::Shape) {
            ve_map.add(&s);
        }
        //
        let a_nb = ve_map.extent();
        for i in 1..=a_nb {
            let a_s = ve_map.at_1(i).clone();
            if !set_add(&mut view, &a_s) {
                continue;
            }
            //
            if shape_data_map::is_bound(the_build, &a_s) {
                ne = shape_data_map::find(the_build, &a_s);
                // keep connection to original edges
                for exp_c in bat::explorer(&ne, ShapeType::Edge, ShapeType::Shape) {
                    let nec = exp_c;
                    let p_le_or = shape_data_map_change_seek(the_edges_origins, &nec);
                    match p_le_or {
                        Some(le_or) => {
                            append_to_list_host(le_or, &a_s);
                        }
                        None => {
                            shape_data_map::bind(the_edges_origins, &nec, vec![a_s.clone()]);
                        }
                    }
                }
                // trim edges
                if ne.shape_type() == ShapeType::Edge {
                    if !the_new_edges.contains(&ne) {
                        the_new_edges.add(&ne);
                        let mut ne_trim = ne.clone();
                        if !trim_edge(&mut ne_trim, the_as_des2d, the_as_des, the_e_trim_e_inf) {
                            return false;
                        }
                    }
                } else {
                    //------------------------------------------------------------
                    // The Intersections are on several edges.
                    // The pieces without intersections with neighbors
                    // are removed from AsDes.
                    //------------------------------------------------------------
                    for exp_c in bat::explorer(&ne, ShapeType::Edge, ShapeType::Shape) {
                        let mut nec = exp_c;
                        if !the_new_edges.contains(&nec) {
                            the_new_edges.add(&nec);
                            if !the_as_des2d.descendant(&nec).is_empty() {
                                if !trim_edge(
                                    &mut nec,
                                    the_as_des2d,
                                    the_as_des,
                                    the_e_trim_e_inf,
                                ) {
                                    return false;
                                }
                            } else if the_as_des.has_ascendant(&nec) {
                                the_as_des.remove(&nec);
                            }
                        }
                    }
                }
            } else {
                if a_s.shape_type() != ShapeType::Edge {
                    continue;
                }
                if set_contains(&a_mf_generated, &fi)
                    && super::brep_offset_make_offset::shape_indexed_data_map_find(&a_dmef, &a_s)
                        .len()
                        == 1
                {
                    continue;
                }

                ne = shape_data_map::value(the_map_sf, &fi).generated(&a_s);
                //// modified by jgv, 19.12.03 for OCC4455 ////
                ne.orientation = a_s.orientation;
                //
                let p_le_or = shape_data_map_change_seek(the_edges_origins, &ne);
                match p_le_or {
                    Some(le_or) => {
                        append_to_list_host(le_or, &a_s);
                    }
                    None => {
                        shape_data_map::bind(the_edges_origins, &ne, vec![a_s.clone()]);
                    }
                }
                //
                if shape_data_map::is_bound(the_mes, &ne) {
                    ne = shape_data_map::find(the_mes, &ne);
                    ne.orientation = a_s.orientation;
                    if !the_new_edges.contains(&ne) {
                        the_new_edges.add(&ne);
                        let mut ne_trim = ne.clone();
                        if !trim_edge(&mut ne_trim, the_as_des2d, the_as_des, the_e_trim_e_inf) {
                            return false;
                        }
                    }
                } else {
                    // OCCT L4906-4913: BRepAdaptor_Curve aBAC(anEdge); type ==
                    // GeomAbs_Line → ExtentEdge + theETrimEInf.Bind.
                    let is_line = matches!(
                        bat::brep_tool_curve(&ne),
                        Some((Curve3::Line(_), _, _))
                    );
                    if is_line {
                        let mut a_new_edge = Shape::null();
                        let _ = BRepOffsetInter2d::extent_edge(the_brep, &ne, &mut a_new_edge, the_offset);
                        shape_data_map::bind(the_e_trim_e_inf, &ne, a_new_edge);
                    }
                }
                the_as_des.add(&nf, &ne);
            }
        }
    }
    true
}

/// OCCT TrimEdge (cxx L4926-5034).
pub(crate) fn trim_edge(
    ne: &mut Shape,
    as_des2d: &BRepAlgoAsDes,
    as_des: &mut BRepAlgoAsDes,
    the_e_trim_e_inf: &mut DataMapOfShapeShape,
) -> bool {
    let mut a_source_edge = Shape::null();
    let (mut v1, mut v2) = top_exp_vertices(ne);
    let (a_t1, a_t2) = brep_tool_range(ne);
    //
    bop_algo_tools_make_split_edge(ne, &v1, a_t1, &v2, a_t2, &mut a_source_edge);
    //
    let a_same_par_tol = PRECISION_CONFUSION;

    let mut u = 0.;
    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;

    let le = as_des2d.descendant(ne).to_vec();
    //
    let mut b_trim = false;
    //
    if le.len() > 1 {
        for it in &le {
            let mut v = it.clone();
            if ne.orientation == Orientation::Reversed {
                // OCCT L4960: V.Reverse().
                v.orientation = bat::top_abs_reverse(v.orientation);
            }
            // V.Orientation(TopAbs_INTERNAL);
            if !find_parameter(&v, ne, &mut u) {
                let (_the_curve, _f, _l) = match brep_tool_curve(ne) {
                    Some(x) => x,
                    None => {
                        // OCCT L4965: theCurve = BRep_Tool::Curve(NE, f, l) —
                        // a null curve reaches GeomAPI_ProjectPointOnCurve
                        // (null-handle dereference in OCCT).  rcad has no null
                        // Geom_Curve value, so the impossible-source state
                        // exits here.
                        return false;
                    }
                };
                let the_point = brep_tool_pnt(&v).unwrap_or(DVec3::ZERO);
                let projector =
                    GeomAPIProjectPointOnCurve::new_point_curve(the_point, &_the_curve);
                if projector.nb_points() == 0 {
                    return false;
                }
                u = projector.lower_distance_parameter();
            }
            if u < u_min {
                u_min = u;
                v1 = v.clone();
            }
            if u > u_max {
                u_max = u;
                v2 = v.clone();
            }
        }
        //
        if v1.is_null() || v2.is_null() {
            return false;
        }
        if !v1.is_same(&v2) {
            // OCCT L4982-4984: NE.Free(true); Or = NE.Orientation();
            // NE.Orientation(TopAbs_FORWARD).
            bat::builder_set_free(ne, true);
            let or = ne.orientation;
            ne.orientation = Orientation::Forward;
            let (vf, vl) = top_exp_vertices(ne);
            bat::builder_remove_edge_vertex(ne, &vf);
            bat::builder_remove_edge_vertex(ne, &vl);
            bat::builder_add_edge_vertex(ne, &bat::oriented(&v1, Orientation::Forward));
            bat::builder_add_edge_vertex(ne, &bat::oriented(&v2, Orientation::Reversed));
            bat::builder_range_edge(ne, u_min, u_max);
            ne.orientation = or;
            as_des.add(ne, &bat::oriented(&v1, Orientation::Forward));
            as_des.add(ne, &bat::oriented(&v2, Orientation::Reversed));
            // OCCT L4998: BRepLib::SameParameter(NE, aSameParTol, true).
            brep_lib_same_parameter_3(ne, a_same_par_tol);
            //
            b_trim = true;
        }
    }
    //
    if !b_trim {
        // OCCT L5007-5014: BRepAdaptor_Curve aBAC(NE); the GeomAbs_Line
        // probe (arch. diff. #54).
        let is_line = matches!(bat::brep_tool_curve(ne), Some((Curve3::Line(_), _, _)));
        if is_line && as_des.has_ascendant(ne) {
            as_des.remove(ne);
        }
    } else if !shape_data_map::is_bound(the_e_trim_e_inf, ne) {
        shape_data_map::bind(the_e_trim_e_inf, ne, a_source_edge);
    }
    true
}

/// OCCT GetEnlargedFaces (cxx L5036-5062).
pub(crate) fn get_enlarged_faces(
    the_faces: &[Shape],
    the_map_sf: &MapSF,
    the_mes: &DataMapOfShapeShape,
    the_faces_origins: &mut DataMapOfShapeShape,
    the_image: &mut BRepAlgoImage,
    the_lsf: &mut Vec<Shape>,
) {
    for it in the_faces {
        let fi = it;
        let ofi = shape_data_map::value(the_map_sf, fi).face();
        if shape_data_map::is_bound(the_mes, &ofi) {
            let a_local_face = shape_data_map::find(the_mes, &ofi);
            the_lsf.push(a_local_face.clone());
            the_image.set_root(&a_local_face);
            shape_data_map::bind(the_faces_origins, &a_local_face, fi.clone());
        }
    }
}

/// OCCT BuildShellsCompleteInter (cxx L5064-5229).
pub(crate) fn build_shells_complete_inter(
    the_lf: &[Shape],
    the_image: &mut BRepAlgoImage,
    the_shells: &mut Shape,
) -> bool {
    // make solids
    let mut a_mv1 = super::brep_offset_make_offset::BOPAlgoMakerVolume::new();
    a_mv1.set_arguments(the_lf);
    // we need to intersect the faces to process the tangential faces
    a_mv1.set_intersect(true);
    a_mv1.set_avoid_internal_shapes(true);
    a_mv1.perform();
    //
    let b_done = !a_mv1.has_errors();
    if !b_done {
        return b_done;
    }
    //
    update_history(the_lf, &a_mv1, the_image);
    //
    let a_result1 = a_mv1.shape();
    if a_result1.shape_type() == ShapeType::Solid {
        // result is the alone solid, nothing to do
        return get_sub_shapes(&a_result1, ShapeType::Shell, the_shells);
    }

    // Since the <theImage> object does not support multiple ancestors,
    // prepare local copy of the origins, which will be used to resolve
    // non-manifold solids produced by Maker Volume algorithm by comparison
    // of the normal directions of the split faces with their origins.
    let mut an_origins: ShapeDataMap<Vec<Shape>> = HashMap::new();
    for a_it_lr in the_image.roots().to_vec() {
        let a_fr = a_it_lr;

        // Find the last splits of the root face, including the ones
        // created during MakeVolume operation
        let mut a_lf_im: Vec<Shape> = Vec::new();
        the_image.last_image(&a_fr, &mut a_lf_im);

        for a_it_lf_im in &a_lf_im {
            let a_f_im = a_it_lf_im;
            let p_lf_or = shape_data_map_change_seek(&mut an_origins, a_f_im);
            match p_lf_or {
                Some(lf_or) => {
                    lf_or.push(a_fr.clone());
                }
                None => {
                    shape_data_map::bind(&mut an_origins, a_f_im, vec![a_fr.clone()]);
                }
            }
        }
    }

    // It is necessary to rebuild the solids, avoiding internal faces
    // Map faces to solids
    let mut a_dmfs: IndexedDataMapOfShapeListOfShape = indexmap::IndexMap::new();
    super::brep_offset_make_offset::top_exp_map_shapes_and_ancestors(
        &a_result1,
        ShapeType::Face,
        ShapeType::Solid,
        &mut a_dmfs,
    );
    //
    let a_nb = a_dmfs.len();
    let b_done = a_nb > 0;
    if !b_done {
        // unable to build any solid
        return b_done;
    }
    //
    // get faces attached to only one solid
    let mut a_lf: Vec<Shape> = Vec::new();
    for i in 0..a_nb {
        let (_k, (_key_shape, a_ls)) = a_dmfs.get_index(i).unwrap();
        if a_ls.len() == 1 {
            let a_f = _key_shape.clone();
            a_lf.push(a_f);
        }
    }
    //
    // make solids from the new list
    let mut a_mv2 = super::brep_offset_make_offset::BOPAlgoMakerVolume::new();
    a_mv2.set_arguments(&a_lf);
    // no need to intersect this time
    a_mv2.set_intersect(false);
    a_mv2.set_avoid_internal_shapes(true);
    a_mv2.perform();
    let b_done = !a_mv2.has_errors();
    if !b_done {
        return b_done;
    }
    //
    let a_result2 = a_mv2.shape();
    if a_result2.shape_type() == ShapeType::Solid {
        return get_sub_shapes(&a_result2, ShapeType::Shell, the_shells);
    }
    //
    let a_exp = bat::explorer(&a_result2, ShapeType::Face, ShapeType::Shape);
    let b_done = !a_exp.is_empty();
    if !b_done {
        return b_done;
    }
    //
    a_lf.clear();
    a_dmfs.clear();

    // the result is non-manifold - resolve it comparing normal
    // directions of the offset faces and original faces
    for a_f in &a_exp {
        let a_f = a_f.clone();
        let p_lf_or = shape_data_map::seek(&an_origins, &a_f);
        let lf_or = match p_lf_or {
            Some(x) => x.clone(),
            None => {
                // OCCT L5180-5183: Standard_ASSERT_INVOKE(...) — the OCCT
                // assertion-invoke path.
                panic!(
                    "BRepOffset_MakeOffset::BuildShellsCompleteInterSplit(): \
                     Origins map does not contain the split face"
                );
            }
        };
        // Check orientation
        for a_it_l_or in &lf_or {
            let a_f_or = a_it_l_or.clone();
            // OCCT L5191: BRepOffset_Tool::CheckPlanesNormals(aF, aFOr) —
            // the hxx default theTolAng = 1.e-8 (BRepOffset_Tool.hxx
            // L204-207).
            if super::brep_offset_tool_d::check_planes_normals(&a_f, &a_f_or, 1.0e-8) {
                a_lf.push(a_f.clone());
                break;
            }
        }
    }
    //
    // make solid from most outer faces with correct normal direction
    let mut a_mv3 = super::brep_offset_make_offset::BOPAlgoMakerVolume::new();
    a_mv3.set_arguments(&a_lf);
    a_mv3.set_intersect(false);
    a_mv3.set_avoid_internal_shapes(true);
    a_mv3.perform();
    let b_done = !a_mv3.has_errors();
    if !b_done {
        return b_done;
    }
    //
    let a_result3 = a_mv3.shape();
    get_sub_shapes(&a_result3, ShapeType::Shell, the_shells)
}

/// OCCT GetSubShapes (cxx L5386-5408).
pub(crate) fn get_sub_shapes(
    the_shape: &Shape,
    the_ss_type: ShapeType,
    the_result: &mut Shape,
) -> bool {
    let a_exp = bat::explorer(the_shape, the_ss_type, ShapeType::Shape);
    if a_exp.is_empty() {
        return false;
    }
    //
    let mut b = rcad_kernel::topo::topods::BRepBuilder::new();
    let mut brep = rcad_kernel::topo::topods::BRep::new();
    let mut a_result = b.make_compound(&mut brep, vec![]);
    //
    for a_ss in a_exp {
        b.add_to_compound(&mut brep, a_result.clone(), a_ss);
    }
    *the_result = a_result;
    true
}

/// OCCT RemoveSeamAndDegeneratedEdges (cxx L5550-5633).
pub(crate) fn remove_seam_and_degenerated_edges(the_face: &Shape, the_old_face: &Shape) {
    let mut a_face = the_face.clone();
    a_face.orientation = Orientation::Forward;

    let mut a_is_deg_or_seam_found = false;
    let mut a_esq: Vec<Shape> = Vec::new();
    for an_explo in bat::explorer(&a_face, ShapeType::Edge, ShapeType::Shape) {
        let an_edge = an_explo;
        if brep_tool_degenerated(&an_edge) || brep_tools_is_really_closed(&an_edge, the_old_face) {
            a_is_deg_or_seam_found = true;
        } else {
            a_esq.push(an_edge);
        }
    }

    if !a_is_deg_or_seam_found {
        return;
    }

    // Reconstruct wires
    let mut a_bb = rcad_kernel::topo::topods::BRepBuilder::new();
    let mut brep = rcad_kernel::topo::topods::BRep::new();
    let mut a_wlist: Vec<Shape> = Vec::new();
    for an_it_face in bat::sub_shapes(&a_face) {
        a_wlist.push(an_it_face);
    }

    bat::builder_set_free(&mut a_face, true);
    for an_itl in &a_wlist {
        // OCCT L5591: aBB.Remove(aFace, anItl.Value()) — the
        // BRep_Builder::Remove(Face, Wire) form (arch. diff. #22).
        remove_from_shape_host(&mut a_face, an_itl);
    }

    while !a_esq.is_empty() {
        let mut a_new_wire = a_bb.make_wire(&mut brep);
        let mut a_cur_edge = a_esq.remove(0);
        a_bb.add_to_wire(&mut brep, a_new_wire.clone(), a_cur_edge.clone());
        // OCCT L5600: TopExp::Vertices(aCurEdge, aFirstVertex, aCurVertex,
        // true) — with orientation.
        let (a_first_vertex, mut a_cur_vertex) = top_exp_vertices_cum_ori(&a_cur_edge);
        while !a_cur_vertex.is_same(&a_first_vertex) {
            let mut a_v1 = Shape::null();
            let mut a_v2 = Shape::null();
            let mut ind = 0usize;
            let mut found = false;
            while ind < a_esq.len() {
                a_cur_edge = a_esq[ind].clone();
                let (t_v1, t_v2) = top_exp_vertices_cum_ori(&a_cur_edge);
                a_v1 = t_v1;
                a_v2 = t_v2;
                ind += 1;
                if a_v1.is_same(&a_cur_vertex) {
                    found = true;
                    break;
                }
            }
            if !found {
                // error occurred: wire is not closed
                break;
            }
            let a_cur_edge = a_esq.remove(ind - 1);
            a_bb.add_to_wire(&mut brep, a_new_wire.clone(), a_cur_edge);
            a_cur_vertex = a_v2;
        }
        bat::builder_add_face_wire(&mut a_face, &a_new_wire);
    }
}

/// OCCT BRep_Builder::Remove(S, SS) — the generic container removal (the
/// rcad TShape make_mut edit, arch. diff. #22; the RemoveShapes /
/// RemoveSeamAndDegeneratedEdges forms).
fn remove_from_shape_host(the_s: &mut Shape, the_ss: &Shape) {
    let key = shape_key(the_ss);
    let d = Arc::make_mut(&mut the_s.data);
    match d {
        TShape::Wire(wd) => {
            wd.edges.retain(|c| shape_key(c) != key);
            wd.my_shapes.retain(|c| shape_key(c) != key);
        }
        TShape::Face(fd) => {
            if fd.outer_wire.is_same(the_ss) {
                fd.outer_wire = Shape::null();
            }
            fd.inner_wires.retain(|c| shape_key(c) != key);
            fd.my_shapes.retain(|c| shape_key(c) != key);
        }
        TShape::Shell(sd) => {
            sd.faces.retain(|c| shape_key(c) != key);
        }
        TShape::Compound(cd) => {
            cd.retain(|c| shape_key(c) != key);
        }
        TShape::Solid(sd) => {
            sd.shells.retain(|c| shape_key(c) != key);
        }
        _ => {}
    }
}

/// OCCT AppendToList (cxx L5647-5659).
pub(crate) fn append_to_list_host(the_list: &mut Vec<Shape>, the_shape: &Shape) {
    for a_it in the_list.iter() {
        let a_s = a_it;
        if a_s.is_same(the_shape) {
            return;
        }
    }
    the_list.push(the_shape.clone());
}


/// OCCT IsSolid (cxx L5637-5645).
pub(crate) fn is_solid(the_s: &Shape) -> bool {
    let a_exp = bat::explorer(the_s, ShapeType::Solid, ShapeType::Shape);
    !a_exp.is_empty()
}

/// OCCT DataMap::ChangeSeek(k) — the rcad form (the shape_data_map module
/// of brep_offset_tool.rs carries no ChangeSeek; the local form keeps the
/// OCCT call structure).
fn shape_data_map_change_seek<'a, V>(
    m: &'a mut ShapeDataMap<V>,
    k: &Shape,
) -> Option<&'a mut V> {
    m.get_mut(&shape_key(k)).map(|e| &mut e.1)
}
