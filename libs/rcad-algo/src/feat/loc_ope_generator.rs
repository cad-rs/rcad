// OCCT LocOpe_Generator.hxx L17-73 + LocOpe_Generator.cxx L17-1217 +
// LocOpe_Generator.lxx L17-65 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Generator.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Generator.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Generator.lxx
//
// OCCT inheritance chain (LocOpe_Generator.hxx L31): none —
// LocOpe_Generator is a standalone value class consuming a
// LocOpe_GeneratedShape handle (the abstract base maps to the
// LocOpeGeneratedShape trait, loc_ope_generated_shape.rs; implementors:
// LocOpe_GluedShape etc.).  (The stage plan's "Generator consumes
// BRepFill_Generator" note does not hold for this OCCT version — the .cxx
// consumes LocOpe_GeneratedShape, BRepAlgo_Loop, GeomProjLib and
// LocOpe_BuildShape; there is no BRepFill dependency — see the delivery
// report.)
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
//    (myModShapes / theEEMap / theFFMap) — HashMap<(TShape ptr, Location),
//    Vec<Shape>>; theEEMap / theFFMap carry the key shape alongside the
//    list (the OCCT iterators expose itf.Key()).  NCollection_Map
//    (theLeft / GEdg / GVtx / toRemove / mapTreated / EdgAdded) —
//    IndexMap / HashSet keyed by the same identity; the insertion order is
//    the deterministic stand-in for the OCCT bucket iteration order (same
//    reduction as loc_ope_build_shape.rs arch. diff. #3).
// 2. TopExp_Explorer is feat::brep_feat_builder::explorer;
//    TopoDS_Iterator is feat::brep_feat_builder::sub_shapes;
//    TopExp::MapShapesAndAncestors is
//    feat::loc_ope_glued_shape::map_shapes_and_ancestors;
//    TopExp::Vertices is top_exp_vertices (loc_ope_generator_b.rs).
// 3. BRep_Builder (TopoDS_Builder) — rcad's BRepBuilder over a local BRep
//    pool owned by Perform (loc_ope_build_shape.rs arch. diff. #1, the
//    AGENTS.md A1 vehicle).
// 4. BRep_Tool::Curve/Surface location branches: the feat pipeline shapes
//    carry identity locations (arch. diff. #1 of loc_ope_find_edges.rs),
//    so the "apply loc.Transformation()" steps are no-ops; the BRep_Tool
//    accessors are re-hosted in loc_ope_generator_b.rs.  Where the OCCT
//    source dereferences a null curve handle on a degenerated edge (cxx
//    L644 / L719 / L1027), the rcad translation skips the sub-shape (the
//    OCCT behavior is an unhandled null-handle dereference).
// 5. BRepTools::UVBounds (BRepTools.cxx L64-80) / IsReallyClosed
//    (BRepTools.cxx L1204-1220), GeomProjLib::Curve2d and the gp/ElCLib
//    statics are re-hosted in loc_ope_generator_b.rs (arch. diffs #4-#6
//    there).
// 6. BRepAlgo_Loop (TKBool, consumed at cxx L1074-1079) is NOT yet ported —
//    BRepAlgoLoop in loc_ope_generator_b.rs carries the consumed interface
//    with unimplemented bodies (the pending BRepAlgo port; the class is
//    1151 lines of its own package).
// 7. The OCCT dead locals `outw` / `newwire` (cxx L443; the MakeWire call
//    at L692 is commented out in the source) are not translated (same rule
//    as loc_ope_build_shape.rs arch. diff. #6).
// 8. Rust borrow discipline: Perform keeps the never-mutated myShape in a
//    local clone and the mutated myModShapes in a local written back before
//    returning; the OCCT "two identical pcurve blocks" (cxx L641-682 and
//    L1022-1065) share the make_pcurve_on_new_face helper below (same
//    statements; both OCCT occurrences are cited in its header).
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_Form::Perform uses
// LocOpe_Generator with a LocOpe_GeneratedShape operand.

use crate::feat::brep_feat_builder::{ explorer, sub_shapes };
use crate::feat::loc_ope_generated_shape::LocOpeGeneratedShape;
use crate::feat::loc_ope_generator_b::{
    basis_curve, basis_surface, brep_tool_curve, brep_tool_curve_on_surface,
    brep_tool_degenerated, brep_tool_parameter, brep_tool_range, brep_tool_surface,
    brep_tool_tolerance, brep_tool_uv_points, brep_tools_is_really_closed,
    brep_tools_uv_bounds, geomproj_lib_curve2d, new_parameter, shape_key,
    shape_reversed, standard_epsilon, tofuse_edge_edge, tofuse_edge_face_vertex,
    tofuse_face_face, top_abs_reverse, top_exp_vertices, with_orientation,
    BRepAlgoLoop,
};
use glam::{ DVec2, DVec3 };
use indexmap::IndexMap;
use rcad_kernel::geom::{
    translate_curve2d, Curve2d, Curve2dEval, CurveEval, Plane, Surface3, SurfaceEval,
};
use rcad_kernel::precision::ANGULAR;
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, BRepBuilder };
use std::collections::{ HashMap, HashSet };
use std::f64::consts::PI;

/// OCCT LocOpe_Generator (LocOpe_Generator.hxx L31-68).
pub struct LocOpeGenerator {
    my_shape: Shape, // OCCT: myShape
    // OCCT: myGen (hxx L63) — declared by OCCT and never assigned by the
    // .cxx (Perform works on its parameter); kept for the 1:1 member map.
    #[allow(dead_code)]
    my_gen: Option<std::sync::Arc<dyn LocOpeGeneratedShape>>,
    my_done: bool, // OCCT: myDone
    my_res: Option<Shape>, // OCCT: myRes (None = null shape)
    my_mod_shapes: HashMap<(u64, u32), Vec<Shape>>, // OCCT: myModShapes
}

impl Default for LocOpeGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeGenerator {
    /// OCCT LocOpe_Generator::LocOpe_Generator() (lxx L21-24).
    pub fn new() -> Self {
        LocOpeGenerator {
            my_shape: Shape::null(),
            my_gen: None,
            my_done: false,
            my_res: None,
            my_mod_shapes: HashMap::new(),
        }
    }

    /// OCCT LocOpe_Generator::LocOpe_Generator(S) (lxx L28-32).
    pub fn with_shape(s: &Shape) -> Self {
        LocOpeGenerator {
            my_shape: s.clone(),
            my_gen: None,
            my_done: false,
            my_res: None,
            my_mod_shapes: HashMap::new(),
        }
    }

    /// OCCT LocOpe_Generator::Init(S) (lxx L36-40).
    pub fn init(&mut self, s: &Shape) {
        self.my_shape = s.clone();
        self.my_done = false;
    }

    /// OCCT LocOpe_Generator::IsDone() (lxx L44-47).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_Generator::Shape() (lxx L51-54).
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT LocOpe_Generator::ResultingShape() (lxx L58-65).
    pub fn resulting_shape(&self) -> Option<&Shape> {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        self.my_res.as_ref()
    }

    /// OCCT LocOpe_Generator::Perform(G) (cxx L72-1217).
    #[allow(clippy::never_loop)]
    #[allow(unused_assignments)] // the c2d/c2d1/dir1/dir2 locals mirror the OCCT null handles / uninitialized gp_Vec (cxx L694, L805)
    pub fn perform(&mut self, g: &mut dyn LocOpeGeneratedShape) {
        // OCCT cxx L74-82.
        if self.my_shape.is_null() {
            panic!("Standard_NullObject");
        }
        self.my_done = false;
        self.my_res = None;
        self.my_mod_shapes.clear();

        // Rust borrow discipline (arch. diff. #8).
        let my_shape = self.my_shape.clone();
        let mut my_mod_shapes = std::mem::take(&mut self.my_mod_shapes);

        // OCCT cxx L84: ledges = G->GeneratingEdges().
        let ledges: Vec<Shape> = g.generating_edges().clone();

        // OCCT cxx L86-95.
        let mut the_left: IndexMap<(u64, u32), Shape> = IndexMap::new(); // "Faces a gauche"
        let mut g_edg: IndexMap<(u64, u32), Shape> = IndexMap::new(); // Edges generateurs
        let mut g_vtx: IndexMap<(u64, u32), Shape> = IndexMap::new(); // Vertex generateurs

        // OCCT cxx L97-131.
        for edg0 in &ledges {
            let edg = edg0.clone();
            g_edg.insert(shape_key(&edg), edg.clone());
            for vtx in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                if !g_vtx.contains_key(&shape_key(&vtx)) {
                    g_vtx.insert(shape_key(&vtx), vtx.clone());
                }
            }
            // OCCT cxx L110-130.
            for fac in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
                let mut hit = false;
                for e3 in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                    if e3.is_same(&edg) && e3.orientation == edg.orientation {
                        the_left.insert(shape_key(&fac), fac.clone());
                        if !my_mod_shapes.contains_key(&shape_key(&fac)) {
                            my_mod_shapes.insert(shape_key(&fac), Vec::new());
                        }
                        hit = true;
                        break;
                    }
                }
                // OCCT cxx L126-129: if (exp3.More()) break.
                if hit {
                    break;
                }
            }
        }

        // OCCT cxx L133-135: theEFMap.
        let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &my_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut the_ef_map,
        );

        // OCCT cxx L137-142: theEEMap / theFFMap / toRemove (the value
        // tuples carry the OCCT iterator key shape).
        let mut the_ee_map: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();
        let mut the_ff_map: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();
        let mut to_remove: HashSet<(u64, u32)> = HashSet::new();

        // OCCT cxx L144-224: "search for face fusions".
        let g_edg_keys: Vec<(u64, u32)> = g_edg.keys().copied().collect();
        for k in &g_edg_keys {
            let edg = g_edg[k].clone();
            // OCCT cxx L148-151.
            let Some((_, ef_faces)) = the_ef_map.get(&shape_key(&edg)) else {
                continue;
            };
            let ef_faces = ef_faces.clone();
            // OCCT cxx L152-158.
            let mut itl2_value: Option<Shape> = None;
            for fac in &ef_faces {
                if !the_left.contains_key(&shape_key(fac)) {
                    itl2_value = Some(fac.clone());
                    break;
                }
            }
            // OCCT cxx L159-223.
            if let Some(fac) = itl2_value {
                // OCCT cxx L164-165: facbis = G->Generated(edg).
                let facbis = g.generated_edge(&edg);
                if tofuse_face_face(&fac, &facbis) {
                    // OCCT cxx L168-195.
                    let mut facbisfound = false;
                    let ff_keys: Vec<(u64, u32)> = the_ff_map.keys().copied().collect();
                    for kf in &ff_keys {
                        let (kface, klist) = &the_ff_map[kf];
                        if kface.is_same(&fac) {
                            continue;
                        }
                        for s in klist.iter() {
                            if s.is_same(&facbis) {
                                facbisfound = true;
                                break;
                            }
                        }
                        if facbisfound {
                            the_ff_map.get_mut(kf).expect("theFFMap").1.push(fac.clone());
                            to_remove.insert(shape_key(&fac));
                            to_remove.insert(shape_key(&edg));
                            break;
                        }
                    }

                    // OCCT cxx L197-217.
                    if !facbisfound {
                        if !the_ff_map.contains_key(&shape_key(&fac)) {
                            the_ff_map.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                        }
                        let mut already = false;
                        for s in the_ff_map[&shape_key(&fac)].1.iter() {
                            if s.is_same(&facbis) {
                                already = true;
                                break;
                            }
                        }
                        if !already {
                            the_ff_map
                                .get_mut(&shape_key(&fac))
                                .expect("theFFMap")
                                .1
                                .push(facbis.clone());
                        }
                        to_remove.insert(shape_key(&edg));
                        to_remove.insert(shape_key(&facbis));
                    }
                } else {
                    // OCCT cxx L219-222: "face generee par edg : on la
                    // marque." (the body is commented out in the OCCT
                    // source).
                }
            } else {
                // OCCT cxx L160: edge "interne" au shell, ou bord libre.
            }
        }

        // OCCT cxx L226-265: "Il faut ici ajouter dans toRemove les edges de
        // connexites entre faces a fusionner avec une meme face de base".
        let ff_keys: Vec<(u64, u32)> = the_ff_map.keys().copied().collect();
        for kf in &ff_keys {
            let klist = the_ff_map[kf].1.clone();
            for itl_val in &klist {
                for ed0 in explorer(itl_val, ShapeType::Edge, ShapeType::Shape) {
                    let ed = ed0.clone();
                    if to_remove.contains(&shape_key(&ed)) {
                        continue;
                    }
                    for itl2 in &klist {
                        if !itl2.is_same(itl_val) {
                            let mut hit = false;
                            for e2 in explorer(itl2, ShapeType::Edge, ShapeType::Shape) {
                                if ed.is_same(&e2) {
                                    to_remove.insert(shape_key(&ed));
                                    hit = true;
                                    break;
                                }
                            }
                            // OCCT cxx L257-260: if (exp2.More()) break.
                            if hit {
                                break;
                            }
                        }
                    }
                }
            }
        }

        // OCCT cxx L267-270.
        let mut rebuild_face: Vec<Shape> = Vec::new();
        let mut map_treated: IndexMap<(u64, u32), Shape> = IndexMap::new();
        let mut dont_fuse: HashMap<(u64, u32), Shape> = HashMap::new();

        // OCCT cxx L272-433: the wire-modification pass over theFFMap.
        let ff_keys: Vec<(u64, u32)> = the_ff_map.keys().copied().collect();
        for kf in &ff_keys {
            let (fac, _gen_list) = {
                let e = &the_ff_map[kf];
                (e.0.clone(), e.1.clone())
            };

            for ed0 in explorer(&fac, ShapeType::Edge, ShapeType::Shape) {
                let edg = ed0.clone();
                if map_treated.contains_key(&shape_key(&edg)) {
                    continue; // OCCT cxx L281-283: "on saute l'edge".
                }

                map_treated.insert(shape_key(&edg), edg.clone());
                for vtx0 in explorer(&edg, ShapeType::Vertex, ShapeType::Shape) {
                    let vtx = vtx0.clone();
                    if !g_vtx.contains_key(&shape_key(&vtx)) {
                        continue;
                    }
                    // OCCT cxx L292: edgbis = G->Generated(vtx).
                    let edgbis = g.generated_vertex(&vtx);

                    // OCCT cxx L294-299.
                    if edgbis.is_null()
                        || brep_tool_degenerated(&edgbis)
                        || !tofuse_edge_edge(&edg, &edgbis)
                    {
                        continue;
                    }
                    // OCCT cxx L301-357: the really-closed branch.
                    if brep_tools_is_really_closed(&edg, &fac) {
                        if !the_ee_map.contains_key(&shape_key(&edg)) {
                            the_ee_map.insert(shape_key(&edg), (edg.clone(), Vec::new()));
                            the_ee_map
                                .get_mut(&shape_key(&edg))
                                .expect("theEEMap")
                                .1
                                .push(edgbis.clone());
                            // OCCT cxx L308: "toujours vrai pour edge double".
                            to_remove.insert(shape_key(&edgbis));
                            let mut fuse_edge = true;
                            let (_vf0, vl0) = top_exp_vertices(&edg);
                            let connect_last = vl0.is_same(&vtx);
                            for eee in explorer(
                                &with_orientation(&fac, Orientation::Forward),
                                ShapeType::Edge,
                                ShapeType::Shape,
                            ) {
                                let orient = eee.orientation;
                                if !eee.is_same(&edg) {
                                    let (vf, vl) = top_exp_vertices(&eee);
                                    if (vf.is_same(&vtx) || vl.is_same(&vtx))
                                        && !to_remove.contains(&shape_key(&eee))
                                    {
                                        fuse_edge = false;
                                        // OCCT cxx L323: "On recherche celui
                                        // qu'il ne faut pas fusionner".
                                        if (vf.is_same(&vtx) && orient == Orientation::Forward)
                                            || (vl.is_same(&vtx)
                                                && orient == Orientation::Reversed)
                                        {
                                            if connect_last {
                                                dont_fuse.insert(
                                                    shape_key(&edg),
                                                    with_orientation(&fac, Orientation::Forward),
                                                );
                                            } else {
                                                dont_fuse.insert(
                                                    shape_key(&edg),
                                                    with_orientation(&fac, Orientation::Reversed),
                                                );
                                            }
                                        } else if connect_last {
                                            dont_fuse.insert(
                                                shape_key(&edg),
                                                with_orientation(&fac, Orientation::Reversed),
                                            );
                                        } else {
                                            dont_fuse.insert(
                                                shape_key(&edg),
                                                with_orientation(&fac, Orientation::Forward),
                                            );
                                        }
                                        break;
                                    }
                                }
                            }
                            // OCCT cxx L352-355.
                            if fuse_edge {
                                to_remove.insert(shape_key(&vtx));
                            }
                        }
                    } else {
                        // OCCT cxx L359-428.
                        if !the_ee_map.contains_key(&shape_key(&edg)) {
                            the_ee_map.insert(shape_key(&edg), (edg.clone(), Vec::new()));
                        }
                        the_ee_map
                            .get_mut(&shape_key(&edg))
                            .expect("theEEMap")
                            .1
                            .push(edgbis.clone());
                        // OCCT cxx L370-384: the dedupe scan (the first
                        // append at L369 makes OK false — OCCT source as
                        // written).
                        let l_list = the_ee_map[&shape_key(&edg)].1.clone();
                        let mut ok = true;
                        for s in &l_list {
                            if s.is_same(&edgbis) {
                                ok = false;
                                break;
                            }
                        }
                        if ok {
                            the_ee_map
                                .get_mut(&shape_key(&edg))
                                .expect("theEEMap")
                                .1
                                .push(edgbis.clone());
                        }

                        // OCCT cxx L386.
                        let efac_list = the_ef_map
                            .get(&shape_key(&edg))
                            .expect("Standard_NoSuchObject")
                            .1
                            .clone();
                        // OCCT cxx L387.
                        let fuse_edge = tofuse_edge_face_vertex(&edg, &fac, &vtx, &to_remove);
                        if !fuse_edge {
                            dont_fuse.insert(shape_key(&edg), fac.clone());
                        } else {
                            for efac in &efac_list {
                                if !efac.is_same(&fac) {
                                    if the_ff_map.contains_key(&shape_key(efac)) {
                                        // OCCT cxx L400: "edge a fusionner".
                                        let fuse2 =
                                            tofuse_edge_face_vertex(&edg, efac, &vtx, &to_remove);
                                        if fuse2 {
                                            to_remove.insert(shape_key(&vtx));
                                        } else {
                                            if to_remove.contains(&shape_key(&vtx)) {
                                                to_remove.remove(&shape_key(&vtx));
                                            }
                                            dont_fuse.insert(shape_key(&edg), efac.clone());
                                        }
                                    } else {
                                        // OCCT cxx L416-417: "on marque comme
                                        // face a reconstruire".
                                        rebuild_face.push(efac.clone());
                                        if to_remove.contains(&shape_key(&vtx)) {
                                            to_remove.remove(&shape_key(&vtx));
                                        }
                                        dont_fuse.insert(shape_key(&edg), efac.clone());
                                    }

                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // OCCT cxx L435-439: RebuildFace -> theFFMap bindings.
        for s in &rebuild_face {
            the_ff_map.insert(shape_key(s), (s.clone(), Vec::new()));
        }

        // OCCT cxx L441-452: the shared locals.  The pool carries the
        // builder allocations (arch. diff. #3); the OCCT dead locals outw /
        // newwire are not translated (arch. diff. #7).
        let mut pool = BRep::new();
        let b = BRepBuilder::new();
        let mut b = b;
        let mut tol: f64 = 0.0;
        let mut prm: f64;
        let mut f: f64 = 0.0;
        let mut l: f64 = 0.0;
        let mut uminc: f64 = 0.0;
        let mut umaxc: f64 = 0.0;
        let mut pf = DVec2::ZERO;
        let mut pl = DVec2::ZERO;

        // OCCT cxx L454-522: "Fusion des edges".
        let ee_keys: Vec<(u64, u32)> = the_ee_map.keys().copied().collect();
        for k_ite in &ee_keys {
            let mut keep_new_edge = false;
            let (edg0, _) = the_ee_map[k_ite].clone();
            let edg = edg0;
            // OCCT cxx L461: BRep_Tool::Range(edg, f, l).
            let (f0, l0) = brep_tool_range(&edg);
            f = f0;
            l = l0;
            // OCCT cxx L462-465: newedg = Edge(edg.EmptyCopied()), FORWARD.
            let mut newedg = pool.empty_copied(&edg);
            newedg.orientation = Orientation::Forward;
            // OCCT cxx L466-511: the vertices of the FORWARD copy.
            for vtx0 in explorer(
                &with_orientation(&edg, Orientation::Forward),
                ShapeType::Vertex,
                ShapeType::Shape,
            ) {
                let vtx = vtx0.clone();
                // OCCT cxx L469.
                prm = brep_tool_parameter(&vtx, &edg);

                // OCCT cxx L471: newvtx.Nullify().
                let mut newvtx: Option<Shape> = None;
                // OCCT cxx L472-495: the edgbis containing vtx, and the
                // differing vertex on it.
                let ee_list = the_ee_map[k_ite].1.clone();
                let mut matched: Option<Shape> = None;
                for edgbis in &ee_list {
                    // OCCT cxx L475-481.
                    let mut hit = false;
                    for v2 in explorer(edgbis, ShapeType::Vertex, ShapeType::Shape) {
                        if v2.is_same(&vtx) {
                            hit = true;
                            break;
                        }
                    }
                    if hit {
                        // OCCT cxx L484-492.
                        for v2 in explorer(edgbis, ShapeType::Vertex, ShapeType::Shape) {
                            if !v2.is_same(&vtx) {
                                newvtx = Some(v2.clone());
                                prm = new_parameter(&edg, &vtx, &newedg, &v2);
                                break;
                            }
                        }
                        matched = Some(edgbis.clone());
                        break;
                    }
                }

                // OCCT cxx L497-510.
                if to_remove.contains(&shape_key(&vtx)) || (prm < l && prm > f) {
                    // OCCT cxx L499: B.Add(newedg, newvtx.Oriented(
                    // vtx.Orientation())) — a null newvtx raises
                    // Standard_NullObject in OCCT.
                    let nv = newvtx.expect("Standard_NullObject");
                    b.add_to_edge(
                        &mut pool,
                        newedg.clone(),
                        with_orientation(&nv, vtx.orientation),
                    );
                    let vtol = brep_tool_tolerance(&nv);
                    b.update_vertex_on_edge(&mut pool, nv.clone(), prm, newedg.clone(), vtol);
                    // OCCT cxx L502: toRemove.Add(itl.Value()) — i-e edgbis.
                    if let Some(m) = &matched {
                        to_remove.insert(shape_key(m));
                    }
                    keep_new_edge = true;
                } else {
                    // OCCT cxx L507-509.
                    b.add_to_edge(
                        &mut pool,
                        newedg.clone(),
                        with_orientation(&vtx, vtx.orientation),
                    );
                    let vtol = brep_tool_tolerance(&vtx);
                    b.update_vertex_on_edge(&mut pool, vtx.clone(), prm, newedg.clone(), vtol);
                }
            }
            // OCCT cxx L512-521.
            if keep_new_edge {
                if !my_mod_shapes.contains_key(&shape_key(&edg)) {
                    my_mod_shapes.insert(shape_key(&edg), Vec::new());
                }
                my_mod_shapes
                    .get_mut(&shape_key(&edg))
                    .expect("myModShapes")
                    .push(newedg.clone());
                to_remove.insert(shape_key(&edg));
            }
        }

        // OCCT cxx L524: EdgAdded.
        let mut edg_added: IndexMap<(u64, u32), Shape> = IndexMap::new();

        // OCCT cxx L526-1090: "Fusion des faces, ou reconstruction".
        let ff_keys: Vec<(u64, u32)> = the_ff_map.keys().copied().collect();
        for kf in &ff_keys {
            let (fac, gen_list) = {
                let e = &the_ff_map[kf];
                (e.0.clone(), e.1.clone())
            };
            // OCCT cxx L533-534.
            let mut mod_face = false;
            let mut listofedg: Vec<Shape> = Vec::new();

            // OCCT cxx L536.
            edg_added.clear();

            // OCCT cxx L538-545.
            let mut orface = Orientation::Forward;
            for fc in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
                if fc.is_same(&fac) {
                    orface = fc.orientation;
                    break;
                }
            }
            // OCCT cxx L546-549.
            let mut newface = pool.empty_copied(&fac);
            newface.orientation = Orientation::Forward;
            // OCCT cxx L550-555: S (trimmed -> basis) and the Geom_Plane
            // down-cast.
            let mut s_surf = brep_tool_surface(&fac);
            if let Some(sv) = s_surf.as_ref() {
                if matches!(sv, Surface3::Trimmed(_)) {
                    s_surf = Some(basis_surface(sv));
                }
            }
            let s_ref = s_surf.clone().expect("S");
            let p: Option<Plane> = match s_ref {
                Surface3::Plane(pln) => Some(pln),
                _ => None,
            };
            let s_for_pcurves = s_surf.clone().expect("S");

            // OCCT cxx L557-1071: the wire loop.
            for w0 in explorer(
                &with_orientation(&fac, Orientation::Forward),
                ShapeType::Wire,
                ShapeType::Shape,
            ) {
                let wir = w0.clone();
                // OCCT cxx L560-567.
                let mut wire_modified = false;
                for e2 in explorer(&wir, ShapeType::Edge, ShapeType::Shape) {
                    if to_remove.contains(&shape_key(&e2))
                        || my_mod_shapes.contains_key(&shape_key(&e2))
                    {
                        wire_modified = true;
                        break;
                    }
                }
                // OCCT cxx L568-575: "wire non modifie".
                if !wire_modified {
                    for e2 in explorer(&wir, ShapeType::Edge, ShapeType::Shape) {
                        listofedg.push(e2.clone());
                    }
                } else {
                    // OCCT cxx L578-688: if (!ModFace).
                    if !mod_face {
                        // OCCT cxx L581-586.
                        if p.is_none() {
                            let bounds = brep_tools_uv_bounds(&fac);
                            uminc = bounds[0];
                            umaxc = bounds[1];
                        }

                        // OCCT cxx L588-688: "premier passage".
                        for itl in &gen_list {
                            let mut facbis = itl.clone();
                            // OCCT cxx L593-616.
                            let of_list = g.oriented_faces().clone();
                            let mut found_of: Option<Shape> = None;
                            for s2 in &of_list {
                                if s2.is_same(&facbis) {
                                    found_of = Some(s2.clone());
                                    break;
                                }
                            }
                            if let Some(s2) = found_of {
                                facbis.orientation = s2.orientation;
                            } else {
                                // OCCT cxx L607.
                                for fc in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
                                    if fc.is_same(&facbis) {
                                        facbis.orientation = fc.orientation;
                                        break;
                                    }
                                }
                            }

                            // OCCT cxx L618-687.
                            for w3 in explorer(&facbis, ShapeType::Wire, ShapeType::Shape) {
                                // OCCT cxx L620-626.
                                let mut has_removed = false;
                                for e2 in explorer(&w3, ShapeType::Edge, ShapeType::Shape) {
                                    if to_remove.contains(&shape_key(&e2)) {
                                        has_removed = true;
                                        break;
                                    }
                                }
                                // OCCT cxx L627.
                                if !has_removed {
                                    // OCCT cxx L629-683 (theNew wire; its
                                    // B.Add to newface is commented out in
                                    // the OCCT source).
                                    for e2 in explorer(&w3, ShapeType::Edge, ShapeType::Shape) {
                                        let edg = e2.clone();
                                        let orient2 = orface.compose(edg.orientation);
                                        listofedg.push(with_orientation(&edg, orient2));
                                        edg_added.insert(shape_key(&edg), edg.clone());
                                        if p.is_none() {
                                            // OCCT cxx L641-682.
                                            make_pcurve_on_new_face(
                                                &edg, &facbis, &s_for_pcurves, &mut pool, &mut b,
                                                &newface, uminc, umaxc, &mut pf, &mut pl,
                                                &mut tol, &mut f, &mut l,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                        // OCCT cxx L688.
                        mod_face = true;
                    }

                    // OCCT cxx L691-923: "reconstruction du wire".
                    {
                        // OCCT cxx L694.
                        let mut c2d: Option<Curve2d> = None;
                        let mut c2d1: Option<Curve2d> = None;

                        // OCCT cxx L697-923.
                        for e2 in explorer(&wir, ShapeType::Edge, ShapeType::Shape) {
                            let edg = e2.clone();
                            let orient = edg.orientation;
                            // OCCT cxx L701-706.
                            if !to_remove.contains(&shape_key(&edg))
                                && !the_ee_map.contains_key(&shape_key(&edg))
                            {
                                listofedg.push(edg.clone());
                            } else if my_mod_shapes.contains_key(&shape_key(&edg))
                                || the_ee_map.contains_key(&shape_key(&edg))
                            {
                                // OCCT cxx L709-716.
                                let newedg = if my_mod_shapes.contains_key(&shape_key(&edg)) {
                                    my_mod_shapes[&shape_key(&edg)][0].clone()
                                } else {
                                    edg.clone()
                                };
                                listofedg.push(with_orientation(&newedg, orient));
                                // OCCT cxx L719-728.
                                let Some((c_raw, f0, l0)) = brep_tool_curve(&newedg) else {
                                    continue;
                                };
                                f = f0;
                                l = l0;
                                let mut c = basis_curve(&c_raw);
                                // OCCT cxx L729-764.
                                if p.is_none() {
                                    let oredonfafw = orient;
                                    let _ = oredonfafw;
                                    tol = brep_tool_tolerance(&newedg);
                                    // OCCT cxx L735-742.
                                    c2d =
                                        brep_tool_curve_on_surface(&edg, &fac).map(|(cv, _, _)| cv);
                                    if !brep_tools_is_really_closed(&edg, &fac) {
                                        if let Some(cv) = &c2d {
                                            b.update_edge_pcurve(
                                                &mut pool,
                                                newedg.clone(),
                                                cv.clone(),
                                                newface.clone(),
                                                tol,
                                            );
                                        }
                                    } else if c2d1.is_none() {
                                        c2d1 = c2d.clone();
                                    } else if orient == Orientation::Forward {
                                        b.update_edge_pcurve_closed(
                                            &mut pool,
                                            newedg.clone(),
                                            c2d.clone().expect("C2d"),
                                            c2d1.clone().expect("C2d1"),
                                            newface.clone(),
                                            f,
                                            l,
                                            tol,
                                        );
                                    } else {
                                        b.update_edge_pcurve_closed(
                                            &mut pool,
                                            newedg.clone(),
                                            c2d1.clone().expect("C2d1"),
                                            c2d.clone().expect("C2d"),
                                            newface.clone(),
                                            f,
                                            l,
                                            tol,
                                        );
                                    }
                                }
                                // OCCT cxx L765-801: AddPart.
                                let mut add_part = false;
                                if let Some(df) = dont_fuse.get(&shape_key(&edg)) {
                                    let df = df.clone();
                                    if !brep_tools_is_really_closed(&edg, &fac) {
                                        if df.is_same(&fac) {
                                            if my_mod_shapes.contains_key(&shape_key(&edg)) {
                                                add_part = true;
                                            }
                                        } else if !to_remove.contains(&shape_key(&edg)) {
                                            add_part = true;
                                        }
                                    } else if my_mod_shapes.contains_key(&shape_key(&edg)) {
                                        // OCCT cxx L785: "edg raccourci".
                                        if orient == df.orientation {
                                            add_part = true;
                                        }
                                    } else if orient == top_abs_reverse(df.orientation) {
                                        add_part = true;
                                    }
                                }
                                // OCCT cxx L802-920.
                                if add_part {
                                    let ee_list2 = the_ee_map[&shape_key(&edg)].1.clone();
                                    let mut dir1 = DVec3::ZERO;
                                    let mut dir2 = DVec3::ZERO;
                                    for e2v in &ee_list2 {
                                        if edg_added.contains_key(&shape_key(e2v)) {
                                            continue;
                                        }
                                        let edgbis = e2v.clone();
                                        // OCCT cxx L813-828: TopoDS_Iterator
                                        // it1(newedg), it2.
                                        let mut shared: Option<Shape> = None;
                                        'vtx_search: for v1 in sub_shapes(&newedg) {
                                            for v2 in sub_shapes(&edgbis) {
                                                if v1.is_same(&v2) {
                                                    shared = Some(v1.clone());
                                                    break 'vtx_search;
                                                }
                                            }
                                        }
                                        if let Some(sh_v) = &shared {
                                            // OCCT cxx L831-833: C->D1(prmvt,
                                            // ptbid, dir1) — the rcad
                                            // CurveEval::tangent_at carries the
                                            // unit tangent; the dot-sign test of
                                            // cxx L853 is scale-invariant.
                                            let prmvt = brep_tool_parameter(sh_v, &newedg);
                                            dir1 = c.tangent_at(prmvt);
                                            // OCCT cxx L835-846.
                                            let Some((cb_raw, f0, l0)) =
                                                brep_tool_curve(&edgbis)
                                            else {
                                                continue;
                                            };
                                            f = f0;
                                            l = l0;
                                            c = basis_curve(&cb_raw);
                                            let prmvt2 = brep_tool_parameter(sh_v, &edgbis);
                                            dir2 = c.tangent_at(prmvt2);
                                        } else {
                                            // OCCT cxx L850.
                                            dir1 = DVec3::new(1.0, 0.0, 0.0);
                                            dir2 = DVec3::new(1.0, 0.0, 0.0);
                                        }
                                        edg_added.insert(shape_key(&edgbis), edgbis.clone());
                                        // OCCT cxx L853-862.
                                        if dir1.dot(dir2) < 0.0 {
                                            listofedg.push(with_orientation(
                                                &edgbis,
                                                top_abs_reverse(orient),
                                            ));
                                        } else {
                                            listofedg.push(with_orientation(&edgbis, orient));
                                        }
                                        // OCCT cxx L864-919.
                                        if p.is_none() {
                                            // "C est la courbe de edgbis, f
                                            // et l s'y rapportent".
                                            let Some(mut ptc) =
                                                geomproj_lib_curve2d(&c, f, l, &s_for_pcurves, tol)
                                            else {
                                                continue;
                                            };
                                            if s_for_pcurves.is_u_periodic() {
                                                let uref;
                                                if dont_fuse.contains_key(&shape_key(&edg)) {
                                                    let oredge =
                                                        dont_fuse[&shape_key(&edg)].orientation;
                                                    if my_mod_shapes
                                                        .contains_key(&shape_key(&edg))
                                                    {
                                                        // OCCT cxx L877-886.
                                                        let (pfv, plv) = brep_tool_uv_points(
                                                            &with_orientation(&edg, oredge),
                                                            &fac,
                                                        );
                                                        pf = pfv;
                                                        pl = plv;
                                                    } else {
                                                        // OCCT cxx L890-899.
                                                        let (pfv, plv) = brep_tool_uv_points(
                                                            &with_orientation(
                                                                &edg,
                                                                top_abs_reverse(oredge),
                                                            ),
                                                            &fac,
                                                        );
                                                        pf = pfv;
                                                        pl = plv;
                                                    }
                                                } else {
                                                    let (pfv, plv) =
                                                        brep_tool_uv_points(&edg, &fac);
                                                    pf = pfv;
                                                    pl = plv;
                                                }
                                                uref = pf.x;
                                                let new_u = ptc.point_at(f).x;
                                                // OCCT cxx L912.
                                                if (new_u - uref).abs()
                                                    > standard_epsilon(surface_u_period(
                                                        &s_for_pcurves,
                                                    ))
                                                {
                                                    ptc = translate_curve2d(
                                                        &ptc,
                                                        DVec2::new(uref - new_u, 0.0),
                                                    );
                                                }
                                            }
                                            // OCCT cxx L918.
                                            b.update_edge_pcurve(
                                                &mut pool,
                                                edgbis.clone(),
                                                ptc.clone(),
                                                newface.clone(),
                                                tol,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // OCCT cxx L925-1069: "Recuperation des edges sur les
                    // faces a fusionner".
                    {
                        // OCCT cxx L931.
                        let mut includeinw = false;
                        for itl in &gen_list {
                            let mut facbis = itl.clone();
                            // OCCT cxx L935.
                            let mut genface = true;
                            // OCCT cxx L936-959.
                            let of_list = g.oriented_faces().clone();
                            let mut found_of: Option<Shape> = None;
                            for s2 in &of_list {
                                if s2.is_same(&facbis) {
                                    found_of = Some(s2.clone());
                                    break;
                                }
                            }
                            if let Some(s2) = found_of {
                                facbis.orientation = s2.orientation;
                            } else {
                                // OCCT cxx L949.
                                genface = false;
                                for fc in
                                    explorer(&my_shape, ShapeType::Face, ShapeType::Shape)
                                {
                                    if fc.is_same(&facbis) {
                                        facbis.orientation = fc.orientation;
                                        break;
                                    }
                                }
                            }

                            // OCCT cxx L961-1068.
                            for w3 in explorer(&facbis, ShapeType::Wire, ShapeType::Shape) {
                                // OCCT cxx L963-969.
                                let mut has_removed = false;
                                for e2 in explorer(&w3, ShapeType::Edge, ShapeType::Shape) {
                                    if to_remove.contains(&shape_key(&e2)) {
                                        has_removed = true;
                                        break;
                                    }
                                }
                                if genface {
                                    includeinw = true;
                                } else {
                                    // OCCT cxx L976-1004.
                                    if has_removed {
                                        for e2 in explorer(
                                            &w3,
                                            ShapeType::Edge,
                                            ShapeType::Shape,
                                        ) {
                                            if !to_remove.contains(&shape_key(&e2)) {
                                                continue;
                                            }
                                            let (vf, vl) = top_exp_vertices(&e2);
                                            let mut is_included_in_w = false;
                                            for v4 in explorer(
                                                &wir,
                                                ShapeType::Vertex,
                                                ShapeType::Shape,
                                            ) {
                                                if v4.is_same(&vf) || v4.is_same(&vl) {
                                                    is_included_in_w = true;
                                                    break;
                                                }
                                            }
                                            if is_included_in_w {
                                                break;
                                            }
                                        }
                                    }
                                }

                                // OCCT cxx L1007-1010.
                                if !includeinw {
                                    continue;
                                }

                                // OCCT cxx L1012-1067.
                                for e2 in explorer(&w3, ShapeType::Edge, ShapeType::Shape) {
                                    let edg = e2.clone();
                                    if !to_remove.contains(&shape_key(&edg))
                                        && !edg_added.contains_key(&shape_key(&edg))
                                    {
                                        let orient2 = orface.compose(edg.orientation);
                                        listofedg.push(with_orientation(&edg, orient2));
                                        edg_added.insert(shape_key(&edg), edg.clone());
                                        if p.is_none() {
                                            // OCCT cxx L1022-1065.
                                            make_pcurve_on_new_face(
                                                &edg, &facbis, &s_for_pcurves, &mut pool, &mut b,
                                                &newface, uminc, umaxc, &mut pf, &mut pl,
                                                &mut tol, &mut f, &mut l,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // OCCT cxx L1072-1089: BRepAlgo_Loop (architecture difference
            // #6 — the port is pending; the interface is consumed 1:1).
            if !listofedg.is_empty() {
                let mut loop_ = BRepAlgoLoop::new();
                loop_.init(&newface);
                loop_.add_const_edges(&listofedg);
                loop_.perform();
                loop_.wires_to_faces();
                let listoffaces = loop_.new_faces();
                to_remove.insert(shape_key(&fac));
                my_mod_shapes.insert(shape_key(&fac), listoffaces.clone());
                for s in &gen_list {
                    my_mod_shapes.insert(shape_key(s), listoffaces.clone());
                }
            }
        }

        // OCCT cxx L435-439 note: the RebuildFace bindings precede this
        // loop in the OCCT order (L435-439 before L527) — they are placed
        // there above.

        // OCCT cxx L1092-1129: the commented-out JAG 16.09.96 block (no
        // translation).

        // OCCT cxx L1131-1216: the final assembly.
        // OCCT cxx L1135: FaceRefOri.
        let mut face_ref_ori: Option<Shape> = None;

        // OCCT cxx L1137-1162.
        let mut lfres: Vec<Shape> = Vec::new();
        for fc in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
            let fac = fc.clone();
            if !the_left.contains_key(&shape_key(&fac)) {
                if !to_remove.contains(&shape_key(&fac)) {
                    lfres.push(fac.clone());
                    if !my_mod_shapes.contains_key(&shape_key(&fac)) {
                        my_mod_shapes.insert(shape_key(&fac), Vec::new());
                    }
                    my_mod_shapes
                        .get_mut(&shape_key(&fac))
                        .expect("myModShapes")
                        .push(fac.clone());
                    if face_ref_ori.is_none() {
                        face_ref_ori = Some(fac);
                    }
                } else if my_mod_shapes.contains_key(&shape_key(&fac)) {
                    // OCCT cxx L1159.
                    let first = my_mod_shapes[&shape_key(&fac)][0].clone();
                    lfres.push(with_orientation(&first, fac.orientation));
                }
            }
        }

        // OCCT cxx L1164-1187.
        let orsolid = my_shape.orientation;
        let of_list = g.oriented_faces().clone();
        for s in &of_list {
            let fac = s.clone();
            if to_remove.contains(&shape_key(&fac)) {
                continue;
            }
            if orsolid == Orientation::Forward {
                lfres.push(fac.clone());
            } else {
                lfres.push(shape_reversed(&fac));
            }
            if !my_mod_shapes.contains_key(&shape_key(&fac)) {
                my_mod_shapes.insert(shape_key(&fac), Vec::new());
            }
            my_mod_shapes
                .get_mut(&shape_key(&fac))
                .expect("myModShapes")
                .push(fac.clone());
        }

        // OCCT cxx L1189-1190.
        let bs = crate::feat::loc_ope_build_shape::LocOpeBuildShape::with_faces(&lfres);
        let mut my_res = bs.shape().cloned();

        // OCCT cxx L1192-1212: "Suite debug du 06.11.96".
        if let Some(res) = &my_res {
            if res.shape_type() == ShapeType::Solid {
                let mut do_flip = false;
                if let Some(r0) = &face_ref_ori {
                    for fc in explorer(res, ShapeType::Face, ShapeType::Shape) {
                        if fc.is_same(r0) {
                            if fc.orientation != r0.orientation {
                                do_flip = true;
                            }
                            break;
                        }
                    }
                }
                if do_flip {
                    // OCCT cxx L1205-1210.
                    let shells = explorer(res, ShapeType::Shell, ShapeType::Shape);
                    let shell_rev = shape_reversed(&shells[0]);
                    let new_sol = b.make_solid(&mut pool, vec![shell_rev]);
                    my_res = Some(new_sol);
                }
            }
        }

        // OCCT cxx L1214: "recodage des regularites..." (no body in the
        // OCCT source).

        self.my_mod_shapes = my_mod_shapes;
        self.my_res = my_res;
        self.my_done = true;
    }

    /// OCCT LocOpe_Generator::DescendantFace(F) (cxx L1221-1230).
    pub fn descendant_face(&self, f: &Shape) -> &Vec<Shape> {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        // OCCT cxx L1229: return myModShapes(F) — the DataMap operator()
        // asserts IsBound.
        self.my_mod_shapes.get(&shape_key(f)).expect("Standard_NoSuchObject")
    }
}

/// OCCT cxx L641-682 and L1022-1065 (two identical occurrences): the pcurve
/// of edg projected onto the new face — the "on met les courbes 2d si on
/// n'est pas sur un plan" block with the U-recadrage into the facette
/// (arch. diff. #8).
#[allow(clippy::too_many_arguments)]
fn make_pcurve_on_new_face(
    edg: &Shape,
    facbis: &Shape,
    s: &Surface3,
    pool: &mut BRep,
    b: &mut BRepBuilder,
    newface: &Shape,
    uminc: f64,
    umaxc: f64,
    pf: &mut DVec2,
    pl: &mut DVec2,
    tol: &mut f64,
    f: &mut f64,
    l: &mut f64,
) {
    // OCCT cxx L643-653: tol / C / loc / trim (identity locations, arch.
    // diff. #4).
    *tol = brep_tool_tolerance(edg);
    let Some((c_raw, f0, l0)) = brep_tool_curve(edg) else {
        // The OCCT source dereferences the null curve handle here; the rcad
        // edge without a 3D curve (degenerated) is skipped (arch. diff. #4).
        return;
    };
    *f = f0;
    *l = l0;
    let c = basis_curve(&c_raw);
    // OCCT cxx L655.
    let Some(mut c2d) = geomproj_lib_curve2d(&c, *f, *l, s, *tol) else {
        return;
    };
    // OCCT cxx L657-673: "Tentative de recalage dans la facette".
    *pf = c2d.point_at(*f);
    *pl = c2d.point_at(*l);
    let tttol = ANGULAR;
    while pf.x.min(pl.x) >= umaxc - tttol {
        c2d = translate_curve2d(&c2d, DVec2::new(-2.0 * PI, 0.0));
        *pf = c2d.point_at(*f);
        *pl = c2d.point_at(*l);
    }
    while pf.x.max(pl.x) <= uminc + tttol {
        c2d = translate_curve2d(&c2d, DVec2::new(2.0 * PI, 0.0));
        *pf = c2d.point_at(*f);
        *pl = c2d.point_at(*l);
    }
    // OCCT cxx L675-678.
    if !brep_tools_is_really_closed(edg, facbis) {
        b.update_edge_pcurve(pool, edg.clone(), c2d.clone(), newface.clone(), *tol);
    }
}

/// OCCT Geom_Surface::UPeriod() restricted to the rcad analytic carriers
/// (architecture difference: the periodicity table lives on the OCCT Geom
/// classes; rcad matches the U-periodic analytic types —
/// Geom_CylindricalSurface / Geom_ToroidalSurface: 2*PI).
fn surface_u_period(s: &Surface3) -> f64 {
    match s {
        Surface3::Cylinder(_) | Surface3::Torus(_) => 2.0 * PI,
        _ => 0.0,
    }
}
