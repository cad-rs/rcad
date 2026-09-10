//! OCCT BRepAlgo_Loop — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_Loop.cxx
//!         (L57-1151) + BRepAlgo_Loop.hxx (L33-113).
//!
//! Architecture differences:
//! 1. NCollection_IndexedDataMap<TopoDS_Shape, List, ShapeMapHasher> ->
//!    indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)> (the entry keeps the
//!    key shape; RemoveKey is the OCCT last-element swap, hence
//!    swap_remove).
//! 2. NCollection_DataMap/NCollection_Map keyed by TopTools_ShapeMapHasher
//!    -> HashMap/HashSet keyed by ShapeKey (TShape + Location).
//! 3. NCollection_Sequence<TopoDS_Shape> -> Vec<Shape> (1-based OCCT
//!    indexing adjusted).
//! 4. BRep_Builder edits are in-place Arc::make_mut mutations (tool.rs);
//!    the OCCT TShape sharing makes the mutation visible through every
//!    handle, rcad mutates the owning copy (UpdateVEmap writes the edited
//!    edge back into every theVEmap entry that carries it).
//! 5. TopoDS::Vertex/Edge/Face casts are no-ops (rcad Shape is untyped).
//! 6. Locations are identity in this pipeline (feat arch. diff. #1).
//! 7. The OCCT_DEBUG_ALGO debug blocks are not translated.

use crate::brep_algo::face_restrictor::BRepAlgoFaceRestrictor;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_is_closed_edge, brep_tool_is_closed_on_surface,
    brep_tool_parameter, brep_tool_pnt, brep_tool_range, brep_tool_surface,
    brep_tool_tolerance, brep_tool_uv_points, builder_add_edge_vertex, builder_add_wire_edge,
    builder_make_wire,
    builder_range_edge, builder_remove_edge_vertex, builder_set_closed, builder_set_free,
    builder_update_vertex_point_tol, builder_update_vertex_tol, empty_copied, explorer,
    is_same_opt, oriented, reversed, shape_is_equal_opt, shape_key, shape_same_opt,
    sub_shapes,
    top_exp_vertices_raw, ShapeKey,
};
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::base::geom_lib::axe_of_inertia;
use rcad_kernel::geom::{Curve2dEval, Surface3};
use rcad_kernel::math::gp::Ax2;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::{HashMap, HashSet};

/// The OCCT vertex -> list-of-edges map type (NCollection_IndexedDataMap).
type VEMap = IndexMap<ShapeKey, (Shape, Vec<Shape>)>;

/// OCCT BRepAlgo_Loop (BRepAlgo_Loop.hxx L33-113) — builds the loops from a
/// set of edges on a face.
pub struct BRepAlgoLoop {
    /// OCCT: myFace.
    my_face: Shape,
    /// OCCT: myConstEdges.
    my_const_edges: Vec<Shape>,
    /// OCCT: myEdges.
    my_edges: Vec<Shape>,
    /// OCCT: myVerOnEdges.
    my_ver_on_edges: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: myNewWires.
    my_new_wires: Vec<Shape>,
    /// OCCT: myNewFaces.
    my_new_faces: Vec<Shape>,
    /// OCCT: myCutEdges.
    my_cut_edges: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: myVerticesForSubstitute.
    my_vertices_for_substitute: HashMap<ShapeKey, (Shape, Shape)>,
    /// OCCT: myImageVV.
    my_image_vv: BRepAlgoImage,
    /// OCCT: myTolConf.
    my_tol_conf: f64,
}

impl Default for BRepAlgoLoop {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoLoop {
    /// OCCT BRepAlgo_Loop::BRepAlgo_Loop() (cxx L57-60) — myTolConf(0.001).
    pub fn new() -> Self {
        BRepAlgoLoop {
            my_face: Shape::null(),
            my_const_edges: Vec::new(),
            my_edges: Vec::new(),
            my_ver_on_edges: HashMap::new(),
            my_new_wires: Vec::new(),
            my_new_faces: Vec::new(),
            my_cut_edges: HashMap::new(),
            my_vertices_for_substitute: HashMap::new(),
            my_image_vv: BRepAlgoImage::new(),
            my_tol_conf: 0.001,
        }
    }

    /// OCCT BRepAlgo_Loop::Init(F) (cxx L64-73) — init with F; the set of
    /// edges must have pcurves on F.
    pub fn init(&mut self, f: &Shape) {
        self.my_const_edges.clear();
        self.my_edges.clear();
        self.my_ver_on_edges.clear();
        self.my_new_wires.clear();
        self.my_new_faces.clear();
        self.my_cut_edges.clear();
        self.my_face = f.clone();
    }

    /// OCCT BRepAlgo_Loop::AddEdge(E, LV) (cxx L125-129) — adds E with LV;
    /// E will be copied and trimmed by vertices in LV.
    pub fn add_edge(&mut self, e: &Shape, lv: &[Shape]) {
        self.my_edges.push(e.clone());
        self.my_ver_on_edges.insert(shape_key(e), lv.to_vec());
    }

    /// OCCT BRepAlgo_Loop::AddConstEdge(E) (cxx L133-136) — adds E as const
    /// edge; E can be in the result.
    pub fn add_const_edge(&mut self, e: &Shape) {
        self.my_const_edges.push(e.clone());
    }

    /// OCCT BRepAlgo_Loop::AddConstEdges(LE) (cxx L140-147) — adds LE as a
    /// set of const edges.
    pub fn add_const_edges(&mut self, le: &[Shape]) {
        for itl in le {
            self.my_const_edges.push(itl.clone());
        }
    }

    /// OCCT BRepAlgo_Loop::SetImageVV(theImageVV) (cxx L151-154) — sets the
    /// Image Vertex - Vertex (the OCCT member copy is a Rust move).
    pub fn set_image_vv(&mut self, the_image_vv: BRepAlgoImage) {
        self.my_image_vv = the_image_vv;
    }

    // NOTE: the "value assigned is never read" warnings inside perform are
    // intentional — they mirror the OCCT reference flow (cxx L657-661 the
    // declared VF/CV/CE/EF locals are written before every read).
    #[allow(unused_assignments)]

    /// OCCT BRepAlgo_Loop::Perform() (cxx L581-786).
    pub fn perform(&mut self) {
        let mut ya_couture = false;

        // (OCCT_DEBUG_ALGO block L587-601 not translated.)

        //------------------------------------------------
        // Cut edges (OCCT L602-615).
        //------------------------------------------------
        let my_edges_snapshot = self.my_edges.clone();
        for itl in &my_edges_snapshot {
            let an_edge = itl; // TopoDS::Edge(itl.Value())
            let mut lce: Vec<Shape> = Vec::new();
            let p_vertices = self.my_ver_on_edges.get(&shape_key(an_edge)).cloned();
            if let Some(p_vertices) = p_vertices {
                self.cut_edge(an_edge, &p_vertices, &mut lce);
                self.my_cut_edges.insert(shape_key(an_edge), lce);
            }
        }

        //-----------------------------------
        // Construction map vertex => edges (OCCT L616-639).
        //-----------------------------------
        let mut mve: VEMap = IndexMap::new();

        // add cut edges.
        let mut emap: HashSet<ShapeKey> = HashSet::new();
        for itl in &my_edges_snapshot {
            let list_len = self.my_cut_edges.get(&shape_key(itl)).map(|l| l.len());
            if let Some(list_len) = list_len {
                for j in 0..list_len {
                    let e = &mut self.my_cut_edges.get_mut(&shape_key(itl)).expect("pLCE")[j];
                    if !emap.insert(shape_key(e)) {
                        continue;
                    }
                    store_in_mve(
                        &self.my_face,
                        e,
                        &mut mve,
                        &mut ya_couture,
                        &mut self.my_vertices_for_substitute,
                        self.my_tol_conf,
                    );
                }
            }
        }

        // add const edges (OCCT L641-652).
        // Sewn edges can be doubled or not in myConstEdges
        // => call only once StoreInMVE which should double them.
        let mut deja_vu: HashSet<ShapeKey> = HashSet::new();
        for j in 0..self.my_const_edges.len() {
            let e = &mut self.my_const_edges[j];
            if !deja_vu.insert(shape_key(e)) {
                continue;
            }
            store_in_mve(
                &self.my_face,
                e,
                &mut mve,
                &mut ya_couture,
                &mut self.my_vertices_for_substitute,
                self.my_tol_conf,
            );
        }

        //-----------------------------------------------
        // Construction of wires and new faces (OCCT L654-663).
        //----------------------------------------------
        let mut vf: Option<Shape> = None;
        let mut cv: Option<Shape> = None;
        let mut ce: Shape = Shape::null();
        let mut ne_shape: Shape = Shape::null();
        let mut ef: Shape = Shape::null();
        let mut end: bool;

        self.update_ve_map(&mut mve);

        let mut used_edges: HashSet<ShapeKey> = HashSet::new();

        while mve.len() > 0 {
            let mut nw = builder_make_wire();
            //--------------------------------
            // Removal of hanging edges (OCCT L670-673).
            //--------------------------------
            remove_pending_edges(&mut mve);

            if mve.len() == 0 {
                break;
            }
            //--------------------------------
            // Start edge (OCCT L679-682).
            //--------------------------------
            let (_, (_, first_list)) = mve.get_index(0).expect("MVE(1)");
            ce = first_list.first().expect("MVE(1).First()").clone();
            ef = ce.clone();
            let (v1, v2) = top_exp_vertices_raw(&ce);
            //--------------------------------
            // VF vertex start of new wire (OCCT L684-694).
            //--------------------------------
            if ce.orientation == Orientation::Forward {
                cv = v1.clone();
                vf = v1;
            } else {
                cv = v2.clone();
                vf = v2;
            }
            let contains_cv = match &cv {
                Some(c) => mve.contains_key(&shape_key(c)),
                None => false,
            };
            if !contains_cv {
                continue;
            }
            // OCCT L699-707: remove CE from the CV entry (IsEqual match).
            if let Some(c) = &cv {
                let list = &mut mve.get_mut(&shape_key(c)).expect("ChangeFromKey(CV)").1;
                if let Some(pos) = list.iter().position(|x| x.is_equal(&ce)) {
                    list.remove(pos);
                }
            }
            end = false;

            while !end {
                //-------------------------------
                // Construction of a wire (OCCT L710-724).
                //-------------------------------
                let (v1, v2) = top_exp_vertices_raw(&ce);
                if !is_same_opt(&cv, &v1) {
                    cv = v1;
                } else {
                    cv = v2;
                }

                builder_add_wire_edge(&mut nw, &ce);
                used_edges.insert(shape_key(&ce));

                let (contains_cv, cv_empty) = match &cv {
                    Some(c) => {
                        let k = shape_key(c);
                        (
                            mve.contains_key(&k),
                            mve.get(&k).map(|(_, l)| l.is_empty()).unwrap_or(true),
                        )
                    }
                    None => (false, true),
                };
                if !contains_cv || cv_empty {
                    end = true;
                } else {
                    let c = cv.clone().expect("CV");
                    let list = &mut mve.get_mut(&shape_key(&c)).expect("ChangeFromKey(CV)").1;
                    end = !select_edge(&self.my_face, &ce, &c, &mut ne_shape, list);
                    if !end {
                        ce = ne_shape.clone();
                        // OCCT L737-741.
                        let k = shape_key(&c);
                        if mve.get(&k).map(|(_, l)| l.is_empty()).unwrap_or(false) {
                            mve.swap_remove(&k);
                        }
                    }
                }
            }
            //--------------------------------------------------
            // Add new wire to the set of wires (OCCT L745-782).
            //------------------------------------------------
            if is_same_opt(&vf, &cv) {
                let vf_ref = vf.clone().unwrap_or_else(Shape::null);
                if same_pnt2d(&vf_ref, &ef, &ce, &self.my_face) {
                    builder_set_closed(&mut nw, true);
                    self.my_new_wires.push(nw.clone());
                } else if brep_tool_tolerance(&vf_ref) < self.my_tol_conf {
                    // OCCT L756-773.
                    let mut vf_mut = vf_ref.clone();
                    builder_update_vertex_tol(&mut vf_mut, self.my_tol_conf);
                    if same_pnt2d(&vf_ref, &ef, &ce, &self.my_face) {
                        builder_set_closed(&mut nw, true);
                        self.my_new_wires.push(nw.clone());
                    }
                    // (OCCT_DEBUG_ALGO else block not translated.)
                }
            }
            // (OCCT_DEBUG_ALGO else block L775-782 not translated.)
        }

        purge_new_edges(&mut self.my_cut_edges, &used_edges);
    }

    // NOTE: the "value assigned is never read" warnings inside cut_edge are
    // intentional — they mirror the OCCT reference flow (cxx L823-826 the
    // declared VF/VL locals are written before every read).
    #[allow(unused_assignments)]

    /// OCCT BRepAlgo_Loop::CutEdge(E, VOnE, NE) const (cxx L790-948) — cuts
    /// the edge E in several edges NE on the vertices VonE.
    pub fn cut_edge(&self, e: &Shape, v_on_e: &[Shape], ne: &mut Vec<Shape>) {
        let we = oriented(e, Orientation::Forward); // TopoDS::Edge(E.Oriented(FORWARD))

        let mut sv: Vec<Shape> = Vec::new();

        for it in v_on_e {
            sv.push(it.clone());
        }
        //--------------------------------
        // Parse vertices on the edge (OCCT L807-810).
        //--------------------------------
        bubble(&we, &mut sv);

        let nb_ver = sv.len();
        //----------------------------------------------------------------
        // Construction of new edges (OCCT L812-822).
        // Note : vertices at the extremities of edges are not
        //        onligatorily in the list of vertices.
        //----------------------------------------------------------------
        if sv.is_empty() {
            ne.push(e.clone());
            return;
        }
        let mut vf: Option<Shape> = None;
        let mut vl: Option<Shape> = None;
        let (f, l) = brep_tool_range(&we);
        let (vf0, vl0) = top_exp_vertices_raw(&we);
        vf = vf0;
        vl = vl0;

        if nb_ver == 2 {
            let sv1 = sv.first().expect("SV(1)");
            let sv2 = sv.get(1).expect("SV(2)");
            // OCCT L830: SV(1).IsEqual(VF) && SV(2).IsEqual(VL).
            if shape_is_equal_opt(sv1, &vf) && shape_is_equal_opt(sv2, &vl) {
                ne.push(e.clone());
                return;
            }
        }
        //----------------------------------------------------
        // Processing of closed edges (OCCT L836-854).
        // If a vertex of intersection is on the common vertex
        // it should appear at the beginning and end of SV.
        //----------------------------------------------------
        if vf.is_some() && is_same_opt(&vf, &vl) {
            let vcei = update_closed_edge(&we, &mut sv);
            if let Some(vcei) = vcei {
                vf = Some(oriented(&vcei, Orientation::Forward));
                vl = Some(oriented(&vcei, Orientation::Reversed));
            }
            let vf_prepend = vf.clone().unwrap_or_else(Shape::null);
            sv.insert(0, vf_prepend);
            let vl_append = vl.clone().unwrap_or_else(Shape::null);
            sv.push(vl_append);
        } else {
            //-----------------------------------------
            // Eventually all extremities of the edge (OCCT L855-868).
            //-----------------------------------------
            if vf.is_some() {
                let first_shape = sv.first().cloned().unwrap_or_else(Shape::null);
                if !shape_same_opt(vf.as_ref().expect("VF"), &Some(first_shape.clone())) {
                    let vf_prepend = vf.clone().unwrap_or_else(Shape::null);
                    sv.insert(0, vf_prepend);
                }
            }
            if vl.is_some() {
                let last_shape = sv.last().cloned().unwrap_or_else(Shape::null);
                if !shape_same_opt(vl.as_ref().expect("VL"), &Some(last_shape)) {
                    let vl_append = vl.clone().unwrap_or_else(Shape::null);
                    sv.push(vl_append);
                }
            }
        }

        while !sv.is_empty() {
            while !sv.is_empty() && sv.first().expect("SV.First()").orientation != Orientation::Forward {
                sv.remove(0);
            }
            if sv.is_empty() {
                break;
            }
            let v1 = sv.remove(0);
            if sv.is_empty() {
                break;
            }
            if sv.first().expect("SV.First()").orientation == Orientation::Reversed {
                let v2 = sv.remove(0);
                //-------------------------------------------
                // Copy the edge and restriction by V1 V2 (OCCT L890-897).
                //-------------------------------------------
                let mut new_edge = empty_copied(&we);
                let a_local_edge = oriented(&v1, Orientation::Forward);
                builder_add_edge_vertex(&mut new_edge, &a_local_edge);
                let a_local_edge = oriented(&v2, Orientation::Reversed);
                builder_add_edge_vertex(&mut new_edge, &a_local_edge);
                // OCCT L898-915.
                let u1 = if shape_same_opt(&v1, &vf) {
                    f
                } else {
                    brep_tool_parameter(&oriented(&v1, Orientation::Internal), &we)
                };
                let u2 = if shape_same_opt(&v2, &vl) {
                    l
                } else {
                    brep_tool_parameter(&oriented(&v2, Orientation::Internal), &we)
                };
                builder_range_edge(&mut new_edge, u1, u2);
                ne.push(oriented(&new_edge, e.orientation));
            }
        }

        // Remove edges with size <= tolerance (OCCT L921-947).
        let tol = 0.001; // 5.e-05; //5.e-07;
        let mut i = 0usize;
        while i < ne.len() {
            // skl : I change "E" to "EE" (OCCT L927).
            let ee = ne[i].clone();
            let (fpar, lpar) = brep_tool_range(&ee);
            if lpar - fpar <= rcad_kernel::precision::CONFUSION {
                ne.remove(i);
            } else {
                let (pf, pl) = brep_tool_uv_points(&ee, &self.my_face);
                if pf.distance(pl) <= tol && !brep_tool_is_closed_edge(&ee) {
                    ne.remove(i);
                } else {
                    i += 1;
                }
            }
        }
    }

    /// OCCT BRepAlgo_Loop::NewWires() const (cxx L952-955) — the list of
    /// wires performed; can be an empty list.
    pub fn new_wires(&self) -> &[Shape] {
        &self.my_new_wires
    }

    /// OCCT BRepAlgo_Loop::NewFaces() const (cxx L959-962) — the list of
    /// faces; WiresToFaces has to be called before; can be empty.
    pub fn new_faces(&self) -> &[Shape] {
        &self.my_new_faces
    }

    /// OCCT BRepAlgo_Loop::WiresToFaces() (cxx L966-992) — builds faces from
    /// the wires result.
    pub fn wires_to_faces(&mut self) {
        if !self.my_new_wires.is_empty() {
            let mut fr = BRepAlgoFaceRestrictor::new();
            let a_local_s = oriented(&self.my_face, Orientation::Forward);
            fr.init(&a_local_s, false, false);
            for it in &self.my_new_wires {
                fr.add(it);
            }

            fr.perform();

            if fr.is_done() {
                let ori_f = self.my_face.orientation;
                while fr.more() {
                    let current = fr.current();
                    self.my_new_faces.push(oriented(&current, ori_f));
                    fr.next();
                }
            }
        }
    }

    /// OCCT BRepAlgo_Loop::NewEdges(E) const (cxx L996-999) — the list of
    /// new edges built from an edge E; it can be empty.
    pub fn new_edges(&self, e: &Shape) -> &[Shape] {
        match self.my_cut_edges.get(&shape_key(e)) {
            Some(l) => l,
            None => panic!("BRepAlgo_Loop::NewEdges"),
        }
    }

    /// OCCT BRepAlgo_Loop::GetVerticesForSubstitute(VerVerMap) const
    /// (cxx L1003-1007) — the datamap of vertices with their substitutes.
    pub fn get_vertices_for_substitute(
        &self,
        ver_ver_map: &mut HashMap<ShapeKey, (Shape, Shape)>,
    ) {
        *ver_ver_map = self.my_vertices_for_substitute.clone();
    }

    /// OCCT BRepAlgo_Loop::VerticesForSubstitute(VerVerMap) (cxx L1011-1015).
    pub fn vertices_for_substitute(
        &mut self,
        ver_ver_map: &HashMap<ShapeKey, (Shape, Shape)>,
    ) {
        self.my_vertices_for_substitute = ver_ver_map.clone();
    }

    /// OCCT BRepAlgo_Loop::SetTolConf(theTolConf) (hxx L95) — the maximal
    /// tolerance used for comparing distances between vertices.
    pub fn set_tol_conf(&mut self, the_tol_conf: f64) {
        self.my_tol_conf = the_tol_conf;
    }

    /// OCCT BRepAlgo_Loop::GetTolConf() const (hxx L98).
    pub fn get_tol_conf(&self) -> f64 {
        self.my_tol_conf
    }

    /// OCCT BRepAlgo_Loop::UpdateVEmap(theVEmap) (cxx L1019-1151) — updates
    /// the VE map according to Image Vertex - Vertex.
    pub fn update_ve_map(&mut self, the_vemap: &mut VEMap) {
        let mut ver_lver: VEMap = IndexMap::new();

        for ii in 0..the_vemap.len() {
            let (_, (key_shape, a_elist)) = the_vemap.get_index(ii).expect("theVEmap(ii)");
            let a_vertex = key_shape; // TopoDS::Vertex(theVEmap.FindKey(ii))
            if a_elist.len() == 1 && self.my_image_vv.is_image(a_vertex) {
                let a_pro_vertex = self.my_image_vv.image_from(a_vertex).clone();
                let kp = shape_key(&a_pro_vertex);
                if ver_lver.contains_key(&kp) {
                    let (_, a_vlist) = ver_lver.get_mut(&kp).expect("VerLver");
                    a_vlist.push(oriented(a_vertex, Orientation::Forward));
                } else {
                    let a_vlist = vec![oriented(a_vertex, Orientation::Forward)];
                    ver_lver.insert(kp, (a_pro_vertex, a_vlist));
                }
            }
        }

        if ver_lver.is_empty() {
            return;
        }

        for ii in 0..ver_lver.len() {
            let (_, (_, a_vlist)) = ver_lver.get_index(ii).expect("VerLver(ii)");
            if a_vlist.len() == 1 {
                continue;
            }

            let mut a_max_tol = 0.0f64;
            let mut points: Vec<DVec3> = Vec::with_capacity(a_vlist.len());

            for itl in a_vlist {
                let a_vertex = itl; // TopoDS::Vertex(itl.Value())
                let a_tol = brep_tool_tolerance(a_vertex);
                a_max_tol = a_max_tol.max(a_tol);
                let a_pnt = brep_tool_pnt(a_vertex).unwrap_or(DVec3::ZERO);
                points.push(a_pnt);
            }

            // OCCT L1075: gp_Ax2 anAxis (the default OXYZ frame).
            let mut an_axis = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
            let mut is_singular = false;
            // OCCT L1077: GeomLib::AxeOfInertia(Points, anAxis, IsSingular)
            // — the 3-argument overload carries the 1.0e-7 default tolerance
            // (GeomLib.hxx AxeOfInertia).
            axe_of_inertia(&points, &mut an_axis, &mut is_singular, 1.0e-7);
            let _ = is_singular;
            let a_centre = an_axis.location;
            let mut a_max_dist = 0.0f64;
            for p in &points {
                let a_sq_dist = (a_centre - *p).length_squared();
                a_max_dist = a_max_dist.max(a_sq_dist);
            }
            a_max_dist = a_max_dist.sqrt();
            a_max_tol = a_max_tol.max(a_max_dist);

            // Find constant vertex (OCCT L1088-1108).
            let mut a_const_vertex: Option<Shape> = None;
            for itl in a_vlist.clone() {
                let a_vertex = itl;
                let a_elist = the_vemap
                    .get(&shape_key(&a_vertex))
                    .map(|(_, l)| l.clone())
                    .unwrap_or_default();
                if a_elist.is_empty() {
                    continue;
                }
                let an_edge = &a_elist[0];
                for itcedges in &self.my_const_edges {
                    if an_edge.is_same(itcedges) {
                        a_const_vertex = Some(a_vertex.clone());
                        break;
                    }
                }
                if a_const_vertex.is_some() {
                    break;
                }
            }
            let a_const_vertex = match a_const_vertex {
                Some(v) => v,
                None => a_vlist.first().expect("aVlist.First()").clone(),
            };
            // OCCT L1113: aBB.UpdateVertex(aConstVertex, aCentre, aMaxTol).
            let mut cv = a_const_vertex.clone();
            builder_update_vertex_point_tol(&mut cv, a_centre, a_max_tol);
            // OCCT TShape sharing: the updated constant vertex is visible
            // through every handle; write the edited copy back into every
            // theVEmap entry carrying it (architecture difference #4).
            replace_shape_in_vemap(the_vemap, &a_const_vertex, &cv);

            for itl in a_vlist.clone() {
                let a_vertex = itl;
                // OCCT L1118: aVertex.IsSame(aConstVertex) — the identity of
                // the ORIGINAL constant vertex (cv is its edited copy).
                if a_vertex.is_same(&a_const_vertex) {
                    continue;
                }

                let a_elist = the_vemap
                    .get(&shape_key(&a_vertex))
                    .map(|(_, l)| l.clone())
                    .unwrap_or_default();
                if a_elist.is_empty() {
                    continue;
                }
                // OCCT L1124-1131: the edge is edited in place (shared
                // TShape); rcad edits a copy and writes it back into every
                // theVEmap entry carrying the edge (architecture diff. #4).
                let mut an_edge = a_elist.first().expect("aElist.First()").clone();
                an_edge.orientation = Orientation::Forward;
                let (a_v1, a_v2) = top_exp_vertices_raw(&an_edge);
                let a_vertex_to_remove = if is_same_opt(&a_v1, &Some(a_vertex.clone())) {
                    a_v1
                } else {
                    a_v2
                };
                builder_set_free(&mut an_edge, true);
                if let Some(to_remove) = a_vertex_to_remove.clone() {
                    builder_remove_edge_vertex(&mut an_edge, &to_remove);
                    let new_v = oriented(&cv, to_remove.orientation);
                    builder_add_edge_vertex(&mut an_edge, &new_v);
                }
                let old_edge = a_elist.first().expect("aElist.First()").clone();
                replace_shape_in_vemap(the_vemap, &old_edge, &an_edge);
            }
        }

        // OCCT L1135-1144: Emap — the indexed map of all theVEmap edges.
        let mut emap: IndexMap<ShapeKey, Shape> = IndexMap::new();
        for ii in 0..the_vemap.len() {
            let (_, (_, a_elist)) = the_vemap.get_index(ii).expect("theVEmap(ii)");
            for itl in a_elist.clone() {
                emap.insert(shape_key(&itl), itl.clone());
            }
        }

        // OCCT L1146-1150: rebuild theVEmap through
        // TopExp::MapShapesAndAncestors(Emap(ii), VERTEX, EDGE, theVEmap).
        the_vemap.clear();
        for (_, e) in emap.into_iter() {
            map_shapes_and_ancestors(&e, ShapeType::Vertex, ShapeType::Edge, the_vemap);
        }
    }
}

/// Write-back helper of the OCCT shared-TShape edit (architecture
/// difference #4): every theVEmap entry whose stored key shape or list
/// shape IsSame(the_old) receives the_new (orientation/location preserved).
/// The replacement is idempotent: after the first pass no entry carries the
/// old TShape pointer, so re-running it is a no-op.  The residual aliasing
/// limitation is that vertex copies held INSIDE unedited edge TShapes keep
/// the pre-edit content, while the OCCT shared handles see every content
/// edit.
fn replace_shape_in_vemap(the_vemap: &mut VEMap, the_old: &Shape, the_new: &Shape) {
    for (_, entry) in the_vemap.iter_mut() {
        if entry.0.is_same(the_old) {
            let mut new_key = the_new.clone();
            new_key.orientation = entry.0.orientation;
            new_key.location = entry.0.location;
            entry.0 = new_key;
        }
        for s in entry.1.iter_mut() {
            if s.is_same(the_old) {
                let mut new_s = the_new.clone();
                new_s.orientation = s.orientation;
                new_s.location = s.location;
                *s = new_s;
            }
        }
    }
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) (TopExp.cxx L87-119; the
/// bop/algo/section.rs model): M[TS shape] = list of TA ancestors, then the
/// TS shapes without ancestors get an empty entry.
fn map_shapes_and_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut VEMap,
) {
    for a_anc in explorer(s, ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, ts, ShapeType::Shape) {
            let key = shape_key(&a_exs);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    for a_ex in explorer(s, ts, ta) {
        let key = shape_key(&a_ex);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// OCCT static Bubble (cxx L80-121) — orders the sequence of vertices by
/// increasing parameter.
fn bubble(e: &Shape, seq: &mut Vec<Shape>) {
    // Remove duplicates (OCCT L83-93; Seq(i) == Seq(j) is the IsEqual test).
    let mut i = 1usize;
    while i < seq.len() {
        let mut j = i + 1usize;
        while j <= seq.len() {
            if seq[i - 1].is_equal(&seq[j - 1]) {
                seq.remove(j - 1);
                // OCCT: j-- then the for-increment keeps the position.
                continue;
            }
            j += 1;
        }
        i += 1;
    }

    let mut invert = true;
    let nb_points = seq.len();

    while invert {
        invert = false;
        for i in 1..nb_points {
            let v1 = oriented(&seq[i - 1], Orientation::Internal);
            let v2 = oriented(&seq[i], Orientation::Internal);

            let u1 = brep_tool_parameter(&v1, e);
            let u2 = brep_tool_parameter(&v2, e);
            if u2 < u1 {
                seq.swap(i - 1, i);
                invert = true;
            }
        }
    }
}

/// OCCT static UpdateClosedEdge (cxx L163-226) — if the first or the last
/// vertex of intersection coincides with the closing vertex, it is removed
/// from SV; it will be added at the beginning and the end of SV by the
/// caller.  Returns the closing vertex (None = the OCCT null vertex).
fn update_closed_edge(e: &Shape, sv: &mut Vec<Shape>) -> Option<Shape> {
    let mut v_res: Option<Shape> = None;
    let mut on_start = false;
    let mut on_end = false;
    // modified by jgv, 13.04.04 for OCC5634 (OCCT L168-171).
    let (v1, _v2) = top_exp_vertices_raw(e);
    let tol = v1.as_ref().map(brep_tool_tolerance).unwrap_or(0.0);

    if sv.is_empty() {
        return v_res;
    }

    let vb0 = sv.first().expect("SV.First()").clone();
    let vb1 = sv.last().expect("SV.Last()").clone();
    let pc = v1.as_ref().and_then(brep_tool_pnt).unwrap_or(DVec3::ZERO);

    for (i, vb) in [&vb0, &vb1].into_iter().enumerate() {
        let p = brep_tool_pnt(vb).unwrap_or(DVec3::ZERO);
        // OCCT L185: P.IsEqual(PC, Tol) — gp_Pnt::IsEqual is Distance <= Tol.
        if (p - pc).length() <= tol {
            v_res = Some(vb.clone());
            if i == 0 {
                on_start = true;
            } else {
                on_end = true;
            }
        }
    }
    if on_start && on_end {
        if !vb0.is_same(&vb1) {
            // (OCCT_DEBUG_ALGO block L202-206 not translated.)
        } else {
            sv.remove(0);
            if !sv.is_empty() {
                sv.remove(sv.len() - 1);
            }
        }
    } else if on_start {
        sv.remove(0);
    } else if on_end {
        sv.remove(sv.len() - 1);
    }

    v_res
}

/// OCCT static RemovePendingEdges (cxx L230-296) — removes hanging edges.
fn remove_pending_edges(mve: &mut VEMap) {
    let mut ya_supress = true;

    while ya_supress {
        ya_supress = false;
        let mut v_to_remove: Vec<Shape> = Vec::new();
        let mut e_to_remove: HashSet<ShapeKey> = HashSet::new();

        for iv in 0..mve.len() {
            let (_, (a_vertex, an_edges)) = mve.get_index(iv).expect("MVE.FindKey(iV)");
            let a_vertex = a_vertex.clone();
            if an_edges.is_empty() {
                v_to_remove.push(a_vertex.clone());
            }
            if an_edges.len() == 1 {
                let e = &an_edges[0]; // TopoDS::Edge(anEdges.First())
                let (v1, v2) = top_exp_vertices_raw(e);
                if !is_same_opt(&v1, &v2) {
                    v_to_remove.push(a_vertex.clone());
                    e_to_remove.insert(shape_key(&an_edges[0]));
                }
            }
        }

        if !v_to_remove.is_empty() {
            ya_supress = true;
            for itl in &v_to_remove {
                mve.swap_remove(&shape_key(itl));
            }
            if !e_to_remove.is_empty() {
                for iv in 0..mve.len() {
                    let (_, (_, le)) = mve.get_index_mut(iv).expect("MVE.ChangeFromIndex(iV)");
                    let mut j = 0usize;
                    while j < le.len() {
                        if e_to_remove.contains(&shape_key(&le[j])) {
                            le.remove(j);
                        } else {
                            j += 1;
                        }
                    }
                }
            }
        }
    }
}

/// OCCT static SamePnt2d (cxx L300-329).
fn same_pnt2d(v: &Shape, e1: &Shape, e2: &Shape, f: &Shape) -> bool {
    let ff = oriented(f, Orientation::Forward); // TopoDS::Face(F.Oriented(FORWARD))
    let Some((c1, f1, l1)) = brep_tool_curve_on_surface(e1, &ff) else {
        // OCCT dereferences the null C1 (cxx L311); rcad returns false
        // (marked).
        return false;
    };
    let Some((c2, f2, l2)) = brep_tool_curve_on_surface(e2, &ff) else {
        // OCCT dereferences the null C2 (cxx L308); rcad returns false
        // (marked).
        return false;
    };
    let p1 = if e1.orientation == Orientation::Forward {
        c1.point_at(f1)
    } else {
        c1.point_at(l1)
    };

    let p2 = if e2.orientation == Orientation::Forward {
        c2.point_at(l2)
    } else {
        c2.point_at(f2)
    };
    let tol = 100.0 * brep_tool_tolerance(v);
    let dist = p1.distance(p2);
    dist < tol
}

/// OCCT static SelectEdge (cxx L339-443) — finds the edge NE connected to CE
/// by vertex CV in the list LE; NE is removed from the list.  If CE is also
/// in the list LE with the same orientation, it is removed from the list.
fn select_edge(
    f: &Shape,
    ce: &Shape,
    cv: &Shape,
    ne: &mut Shape,
    le: &mut Vec<Shape>,
) -> bool {
    *ne = Shape::null(); // OCCT L346: NE.Nullify()
    // OCCT L356-363: remove the CE-equal entry.
    if let Some(pos) = le.iter().position(|x| x.is_equal(ce)) {
        le.remove(pos);
    }
    if le.len() > 1 {
        //--------------------------------------------------------------
        // Several edges possible (OCCT L364-373).
        // - Test edges different from CE, Selection of edge for which CV
        //   has U,V closer to the face than corresponding to CE.
        // - If several edges give representation less than the tolerance,
        //   discrimination on tangents.
        //--------------------------------------------------------------
        let f_forward = oriented(f, Orientation::Forward);

        let Some((c_ce, f0, l0)) = brep_tool_curve_on_surface(ce, &f_forward) else {
            // OCCT dereferences the null curve (cxx L379); rcad returns
            // false (marked).
            return false;
        };
        let mut k = 1usize;
        let mut kmin = 0usize;
        let mut distmin = 100.0 * brep_tool_tolerance(cv);
        let u_ce = if ce.orientation == Orientation::Forward { l0 } else { f0 };

        let pv = c_ce.point_at(u_ce);

        for itl in le.iter() {
            let e = itl; // TopoDS::Edge(itl.Value())
            if !e.is_same(ce) {
                let Some((c, f1, l1)) = brep_tool_curve_on_surface(e, &f_forward) else {
                    // OCCT dereferences the null curve (cxx L399); rcad skips
                    // the candidate (marked).
                    k += 1;
                    continue;
                };
                let u = if e.orientation == Orientation::Forward { f1 } else { l1 };
                let p2 = c.point_at(u);
                let dist = pv.distance(p2);
                if dist <= distmin {
                    kmin = k;
                    distmin = dist;
                }
            }
            k += 1;
        }
        if kmin == 0 {
            return false;
        }

        // OCCT L423-431: the kmin-th entry becomes NE and is removed.
        *ne = le[kmin - 1].clone();
        le.remove(kmin - 1);
    } else if le.len() == 1 {
        *ne = le.first().expect("LE.First()").clone();
        le.remove(0);
    } else {
        return false;
    }
    true
}

/// OCCT static PurgeNewEdges (cxx L447-471).
fn purge_new_edges(
    new_edges: &mut HashMap<ShapeKey, Vec<Shape>>,
    used_edges: &HashSet<ShapeKey>,
) {
    for (_, lne) in new_edges.iter_mut() {
        let mut j = 0usize;
        while j < lne.len() {
            let ne_shape = &lne[j];
            if !used_edges.contains(&shape_key(ne_shape)) {
                lne.remove(j);
            } else {
                j += 1;
            }
        }
    }
}

/// OCCT NCollection_DataMap::Bind semantics for the substitute map: an
/// existing key's value is replaced.
fn bind_or_replace(map: &mut HashMap<ShapeKey, (Shape, Shape)>, key: &Shape, value: &Shape) {
    map.insert(shape_key(key), (key.clone(), value.clone()));
}

/// OCCT static StoreInMVE (cxx L475-577).
fn store_in_mve(
    f: &Shape,
    e: &mut Shape,
    mve: &mut VEMap,
    ya_couture: &mut bool,
    vertices_for_substitute: &mut HashMap<ShapeKey, (Shape, Shape)>,
    the_tol_conf: f64,
) {
    for iv in 0..mve.len() {
        let (_, (key_shape, _)) = mve.get_index(iv).expect("MVE.FindKey(iV)").clone();
        let v = key_shape; // TopoDS::Vertex(MVE.FindKey(iV))
        let p = brep_tool_pnt(&v).unwrap_or(DVec3::ZERO);
        let v_list = sub_shapes(e); // TopoDS_Iterator VerExp(E)
        for itl in &v_list {
            let v1 = itl; // TopoDS::Vertex(itl.Value())
            let p1 = brep_tool_pnt(v1).unwrap_or(DVec3::ZERO);
            // OCCT L504: P.IsEqual(P1, theTolConf) — Distance <= Tol.
            if (p - p1).length() <= the_tol_conf && !v.is_same(v1) {
                // OCCT L506: V.Orientation(V1.Orientation()).
                let mut v_new = v.clone();
                v_new.orientation = v1.orientation;

                // OCCT L507-539: the substitute bookkeeping.
                let kv1 = shape_key(v1);
                if vertices_for_substitute.contains_key(&kv1) {
                    let old_new_v = vertices_for_substitute.get(&kv1).expect("VfS(V1)").1.clone();
                    if !old_new_v.is_same(&v_new) {
                        bind_or_replace(vertices_for_substitute, &old_new_v, &v_new);
                        vertices_for_substitute.get_mut(&kv1).expect("VfS(V1)").1 = v_new.clone();
                    }
                } else {
                    let kv = shape_key(&v_new);
                    if vertices_for_substitute.contains_key(&kv) {
                        let new_new_v = vertices_for_substitute.get(&kv).expect("VfS(V)").1.clone();
                        if !new_new_v.is_same(v1) {
                            bind_or_replace(vertices_for_substitute, v1, &new_new_v);
                        }
                    } else {
                        bind_or_replace(vertices_for_substitute, v1, &v_new);
                        let keys: Vec<ShapeKey> = vertices_for_substitute.keys().copied().collect();
                        for k in keys {
                            let entry = vertices_for_substitute.get(&k).expect("VfS");
                            if entry.1.is_same(v1) {
                                vertices_for_substitute.get_mut(&k).expect("VfS").1 = v_new.clone();
                            }
                        }
                    }
                }
                // OCCT L540-542: E.Free(true); BB.Remove(E, V1); BB.Add(E, V).
                builder_set_free(e, true);
                builder_remove_edge_vertex(e, v1);
                builder_add_edge_vertex(e, &v_new);
            }
        }
    }

    // OCCT L547-565.
    let (v1, v2) = top_exp_vertices_raw(e);
    if v1.is_none() && v2.is_none() {
        *ya_couture = false;
        return;
    }
    if let Some(vv) = &v1 {
        let k = shape_key(vv);
        if !mve.contains_key(&k) {
            mve.insert(k, (vv.clone(), Vec::new()));
        }
        mve.get_mut(&k).expect("ChangeFromKey(V1)").1.push(e.clone());
    }
    // OCCT L558: if (!V1.IsSame(V2)).
    if !is_same_opt(&v1, &v2) {
        let vv = v2.clone().unwrap_or_else(Shape::null);
        let k = shape_key(&vv);
        if !mve.contains_key(&k) {
            mve.insert(k, (vv.clone(), Vec::new()));
        }
        mve.get_mut(&k).expect("ChangeFromKey(V2)").1.push(e.clone());
    }
    // OCCT L566-576: the seam-edge doubling.
    let s: Option<Surface3> = brep_tool_surface(f);
    let is_closed = match &s {
        Some(_) => brep_tool_is_closed_on_surface(e, f),
        // OCCT would use the null surface handle; rcad treats the IsClosed
        // test as false (marked).
        None => false,
    };
    if is_closed {
        let vv = v2.clone().unwrap_or_else(Shape::null);
        let k = shape_key(&vv);
        mve.get_mut(&k).expect("ChangeFromKey(V2)").1.push(reversed(e));
        if !is_same_opt(&v1, &v2) {
            let vv1 = v1.clone().unwrap_or_else(Shape::null);
            let k1 = shape_key(&vv1);
            mve.get_mut(&k1).expect("ChangeFromKey(V1)").1.push(reversed(e));
        }
        *ya_couture = true;
    }
}
