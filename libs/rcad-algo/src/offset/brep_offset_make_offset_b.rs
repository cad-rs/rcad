// OCCT BRepOffset_MakeOffset.cxx — 1:1 translation, module b (see the
// module header of brep_offset_make_offset.rs for the split map and the
// architecture-difference list #38-#56).
//
// Module b carries cxx L837-L2394: MakeOffsetShape / MakeThickSolid /
// MakeOffsetFaces / BuildOffsetByInter / ReplaceRoots / BuildOffsetByArc /
// ToContext / UpdateFaceOffset.
//
// Extra architecture note for this module: the OCCT
// `BRep_Builder::Add(Solid, Shell)` step form has no incremental rcad
// builder (BRepBuilder carries make_solid(brep, shells) only), so the
// solid-shell accumulation keeps the OCCT branch structure and constructs
// the solid at the point of the OCCT finalization.

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_inter2d::IndexedShapeMap;
use super::brep_offset_inter2d_b::BRepOffsetInter2d;
use super::brep_offset_inter3d::BRepOffsetInter3d;
use super::brep_offset_make_offset::{
    analyse_add_faces, analyse_add_faces_rt, analyse_edges, analyse_has_generated,
    analyse_perform, analyse_set_face_offset_map, analyse_set_offset_value,
    brep_lib_sort_faces, brep_lib_update_tolerances, normalize_steps, remove_corks,
    BRepOffset_Error, BRepOffsetMakeOffset, BRepToolsQuilt,
    DataMapOfShapeListOfShape, DataMapOfShapeShape, IndexedDataMapOfShapeListOfShape,
    MapSF,
};
use super::brep_offset_offset::{BRepOffsetStatus, GeomAbsShapeKind};
use super::brep_offset_offset_b::BRepOffsetOffset;
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, top_exp_vertices, OcctIndexedShapeMap, OcctShapeSet,
};
use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool as bat;
use crate::brep_fill::offset_wire::GeomAbsJoinType;
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_tolerance;
use crate::brep_algo::tool::shape_key;
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;

impl BRepOffsetMakeOffset {
    /// OCCT BRepOffset_MakeOffset::MakeOffsetShape (cxx L837-1081).
    pub fn make_offset_shape(&mut self) {
        self.my_done = false;
        //

        // check if shape consists of only planar faces
        self.my_is_planar = self.is_planar();

        self.set_faces();
        self.set_faces_with_offset();

        self.build_face_comp();

        //------------------------------------------
        // Construction of myShape without caps.
        //------------------------------------------
        if !self.my_faces.is_empty() {
            let mut my_shape = self.my_shape.clone();
            remove_corks(&mut my_shape, &mut self.my_original_faces);
            self.my_shape = my_shape;
            let mut my_face_comp = self.my_face_comp.clone();
            remove_corks(&mut my_face_comp, &mut self.my_faces);
            self.my_face_comp = my_face_comp;
        }

        // OCCT L883-890: Message_ProgressScope aPS(theRange, "Making offset
        // shape", 100); NCollection_Array1<double> aSteps(0, Last-1);
        // analyzeProgress(100., aSteps) — the flattened rcad scope; the
        // aPS.Next(aSteps(...)) plumbing is result-inert (architecture
        // difference #40).
        let a_prog = NoopProgress;
        let mut a_ps = ProgressScope::new(&a_prog, "Making offset shape", 100);
        let mut a_steps = vec![0f64; PIO_OPERATION_LAST_ARRAY];
        self.analyze_progress(100., &mut a_steps);
        let _ = &a_steps;

        if !self.check_input_data() || self.my_error != BRepOffset_Error::NoError {
            // There is error in input data.
            // Check Error() method.
            return;
        }
        self.my_error = BRepOffset_Error::NoError;
        let mut side = State::In;
        if self.my_offset < 0. {
            side = State::Out;
        }

        // ------------
        // Preanalyse.
        // ------------
        let mut my_tol = self.my_tol;
        super::brep_offset_make_offset::eval_max(&self.my_shape, &mut my_tol);
        self.my_tol = my_tol;
        // There are possible second variant: analytical continuation of arcsin.
        let tol_angle_coeff =
            (self.my_tol / ((self.my_offset * 0.5).abs() + rcad_kernel::core::precision::CONFUSION))
                .min(1.0);
        let tol_angle = 4. * tol_angle_coeff.asin();
        if (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_inter
            && self.my_is_planar
        {
            analyse_set_offset_value(&mut self.my_analyse, self.my_offset);
            analyse_set_face_offset_map(&mut self.my_analyse, &self.my_face_offset);
        }
        analyse_perform(&mut self.my_analyse, &self.my_face_comp, tol_angle);
        let an_e_exp = bat::explorer(&self.my_face_comp, ShapeType::Edge, ShapeType::Shape);
        for an_e in &an_e_exp {
            let a_li = self.my_analyse.type_(an_e);
            if a_li.is_empty() {
                continue;
            }
            if a_li.last().unwrap().type_of() == ChFiDS_TypeOfConcavity::Mixed {
                self.my_error = BRepOffset_Error::MixedConnectivity;
                return;
            }
        }
        if a_ps.user_break() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }
        //---------------------------------------------------
        // Construction of Offset from preanalysis.
        //---------------------------------------------------
        //----------------------------
        // MaJ of SD Face - Offset
        //----------------------------
        self.update_face_offset();

        if self.my_join == GeomAbsJoinType::Arc {
            self.build_offset_by_arc();
        } else if self.my_join == GeomAbsJoinType::Intersection {
            self.build_offset_by_inter();
        }
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //-----------------
        // Auto unwinding.
        //-----------------
        // if (mySelfInter)  SelfInter(Modif);
        //-----------------
        // Intersection 3d .
        //-----------------

        // OCCT L984: BRepOffset_Inter3d Inter(myAsDes, Side, myTol) — the
        // OCCT Handle aliasing maps to the take + boundary-resync form
        // (architecture difference #39).
        let mut inter =
            BRepOffsetInter3d::new(std::mem::take(&mut self.my_as_des), side, self.my_tol);
        self.intersection_3d(&mut inter);
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //-----------------
        // Intersection2D
        //-----------------
        // OCCT L1000-1001: Modif = Inter.TouchedFaces(); NewEdges =
        // Inter.NewEdges(); — the references are read-only after
        // Intersection3d; the rcad forms are the boundary clones
        // (architecture difference #39).
        // OCCT aliasing resync — the TShape-free swap form (the
        // BRepAlgoAsDes carrier has no Clone).
        std::mem::swap(&mut self.my_as_des, inter.as_des_mut());
        let modif = indexed_shape_map_clone(inter.touched_faces());
        let new_edges = indexed_shape_map_view(inter.new_edges());

        if !modif.is_empty() {
            self.intersection_2d(&modif, &new_edges);
            if self.my_error != BRepOffset_Error::NoError {
                return;
            }
        }

        //-------------------------------------------------------
        // Unwinding 2D and reconstruction of modified faces
        //----------------------------------------------------
        self.make_loops(&modif);
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //-----------------------------------------------------
        // Reconstruction of non modified faces sharing
        // reconstructed edges
        //------------------------------------------------------
        if !modif.is_empty() {
            self.make_faces(&modif);
            if self.my_error != BRepOffset_Error::NoError {
                return;
            }
        }

        if self.my_thickening {
            self.make_missing_walls();
            if self.my_error != BRepOffset_Error::NoError {
                return;
            }
        }

        //-------------------------
        // Construction of shells.
        //-------------------------
        self.make_shells();
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        if self.my_offset_shape.is_null() {
            // not done
            self.my_done = false;
            return;
        }
        //--------------
        // Unwinding 3D.
        //--------------
        self.select_shells();
        //----------------------------------
        // Remove INTERNAL edges if necessary
        //----------------------------------
        if self.my_remove_int_edges {
            self.remove_internal_edges();
        }
        //----------------------------------
        // Coding of regularities.
        //----------------------------------
        self.encode_regularity();
        //----------------------------------
        // Replace roots in history maps
        //----------------------------------
        self.replace_roots();
        //----------------------
        // Creation of solids.
        //----------------------
        self.make_solid();
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //-----------------------------
        // MAJ Tolerance edge and Vertex
        // ----------------------------
        if !self.my_offset_shape.is_null() {
            if self.my_thickening {
                let mut my_offset_shape = self.my_offset_shape.clone();
                super::brep_offset_make_offset_d::update_tolerance(
                    &mut my_offset_shape,
                    &self.my_faces,
                    &self.my_shape,
                );
                self.my_offset_shape = my_offset_shape;
            } else {
                let a_dummy = Shape::null();
                let mut my_offset_shape = self.my_offset_shape.clone();
                super::brep_offset_make_offset_d::update_tolerance(
                    &mut my_offset_shape,
                    &self.my_faces,
                    &a_dummy,
                );
                self.my_offset_shape = my_offset_shape;
            }
            brep_lib_update_tolerances(&mut self.my_offset_shape);
        }

        self.correct_conical_faces();

        // Result solid should be computed in MakeOffset scope.
        if self.my_thickening && self.my_is_perform_sewing {
            let mut a_sew = super::bi_tgte_blended::BRepBuilderAPISewing::new(self.my_tol);
            a_sew.add(&self.my_offset_shape);
            a_sew.perform();
            if a_ps.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            self.my_offset_shape = a_sew.sewed_shape();

            // Rebuild solid.
            // Offset shape expected to be really closed after sewing.
            bat::builder_set_closed(&mut self.my_offset_shape, true);
            self.make_solid();
            if self.my_error != BRepOffset_Error::NoError {
                return;
            }
        }

        self.my_done = true;
    }

    /// OCCT BRepOffset_MakeOffset::MakeThickSolid (cxx L1083-1179).
    pub fn make_thick_solid(&mut self) {
        //--------------------------------------------------------------
        // Construction of shell parallel to shell (initial without cap).
        //--------------------------------------------------------------
        self.make_offset_shape();

        if !self.my_done {
            // Save return code and myDone state.
            return;
        }

        //--------------------------------------------------------------------
        // Construction of a solid with the initial shell, parallel shell
        // limited by caps.
        //--------------------------------------------------------------------
        if !self.my_faces.is_empty() {
            let mut b = rcad_kernel::topo::topods::BRepBuilder::new();
            let mut nb_f = self.my_faces.extent();

            // OCCT L1106: B.MakeSolid(Res) — the rcad shells-Vec form
            // (module-b architecture note).
            let mut res_shells: Vec<Shape> = Vec::new();

            let mut glue = BRepToolsQuilt::new();
            for exp in bat::explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
                nb_f += 1;
                glue.add(&exp);
            }
            let mut ya_result = false;
            if !self.my_offset_shape.is_null() {
                for exp in
                    bat::explorer(&self.my_offset_shape, ShapeType::Face, ShapeType::Shape)
                {
                    ya_result = true;
                    glue.add(&bat::reversed(&exp));
                }
            }

            if !ya_result {
                self.my_done = false;
                self.my_error = BRepOffset_Error::UnknownError;
                return;
            }

            self.my_offset_shape = glue.shells();
            for exp in
                bat::explorer(&self.my_offset_shape, ShapeType::Shell, ShapeType::Shape)
            {
                // OCCT L1140: B.Add(Res, exp.Current()).
                res_shells.push(exp);
            }
            let res = b.make_solid(&mut self.my_brep, res_shells);
            let mut res = res;
            bat::builder_set_closed(&mut res, true);
            self.my_offset_shape = res;

            // Test of Validity of the result of thick Solid
            // more face than the initial solid.
            let mut nb_of = 0;
            for _exp in bat::explorer(&self.my_offset_shape, ShapeType::Face, ShapeType::Shape) {
                nb_of += 1;
            }
            if nb_of < nb_f {
                self.my_done = false;
                self.my_error = BRepOffset_Error::UnknownError;
                return;
            }
            if nb_of == nb_f {
                self.my_offset = 0.;
            }
        }

        if self.my_offset > 0. {
            // OCCT L1188: myOffsetShape.Reverse().
            self.my_offset_shape = bat::reversed(&self.my_offset_shape);
        }

        self.my_done = true;
    }

    /// OCCT BRepOffset_MakeOffset::MakeOffsetFaces (cxx L1202-1275).
    pub(crate) fn make_offset_faces(&mut self, the_map_sf: &mut MapSF) {
        let mut a_cur_offset: f64;
        //
        let mut shape_tgt: HashMap<crate::brep_algo::tool::ShapeKey, Shape> = HashMap::new();
        //
        let offset_outside = self.my_offset > 0.;
        //
        let mut a_lf: Vec<Shape> = Vec::new();
        brep_lib_sort_faces(&self.my_face_comp, &mut a_lf);
        //
        for it_lf in a_lf.clone() {
            let a_f = it_lf;
            a_cur_offset = if shape_data_map::is_bound(&self.my_face_offset, &a_f) {
                shape_data_map::find(&self.my_face_offset, &a_f)
            } else {
                self.my_offset
            };
            // OCCT L1236: BRepOffset_Offset OF(aF, aCurOffset, ShapeTgt,
            // OffsetOutside, myJoin) — the Created map is read-only in the
            // OCCT ctor (cxx L403-411); the join type folds into the Arc
            // flag (architecture difference #38).
            let mut of = BRepOffsetOffset::with_face_created(
                &a_f,
                a_cur_offset,
                &shape_tgt,
                offset_outside,
                self.my_join == GeomAbsJoinType::Arc,
            );
            let mut let_: Vec<Shape> = Vec::new();
            analyse_edges(
                &self.my_analyse,
                &a_f,
                ChFiDS_TypeOfConcavity::Tangential,
                &mut let_,
            );
            for itl in let_.clone() {
                let cur = itl;
                if !shape_tgt.contains_key(&shape_key(&cur))
                    && !analyse_has_generated(&self.my_analyse, &cur)
                {
                    // OCCT L1243-1244: OTE = OF.Generated(Cur);
                    // ShapeTgt.Bind(Cur, OF.Generated(Cur)).
                    let ote = of.generated(&cur);
                    shape_tgt.insert(shape_key(&cur), of.generated(&cur));
                    let (v1, v2) = top_exp_vertices(&cur);
                    let (ov1, ov2) = top_exp_vertices(&ote);
                    if !shape_tgt.contains_key(&shape_key(&v1)) {
                        let mut le: Vec<Shape> = Vec::new();
                        analyse_edges(
                            &self.my_analyse,
                            &v1,
                            ChFiDS_TypeOfConcavity::Tangential,
                            &mut le,
                        );
                        let la = self.my_analyse.ancestors(&v1);
                        if le.len() == la.len() {
                            shape_tgt.insert(shape_key(&v1), ov1.clone());
                        }
                    }
                    if !shape_tgt.contains_key(&shape_key(&v2)) {
                        // OCCT L1256: LE.Clear().
                        let mut le: Vec<Shape> = Vec::new();
                        analyse_edges(
                            &self.my_analyse,
                            &v2,
                            ChFiDS_TypeOfConcavity::Tangential,
                            &mut le,
                        );
                        let la = self.my_analyse.ancestors(&v2);
                        if le.len() == la.len() {
                            shape_tgt.insert(shape_key(&v2), ov2.clone());
                        }
                    }
                }
            }
            shape_data_map::bind(the_map_sf, &a_f, of);
        }
        //
        let a_new_faces = self.my_analyse.new_faces();
        for it in &a_new_faces {
            let a_f = it;
            let of = BRepOffsetOffset::with_face_created(
                a_f,
                0.0,
                &shape_tgt,
                offset_outside,
                self.my_join == GeomAbsJoinType::Arc,
            );
            shape_data_map::bind(the_map_sf, a_f, of);
        }
    }

    /// OCCT BRepOffset_MakeOffset::BuildOffsetByInter (cxx L1277-1823).
    pub(crate) fn build_offset_by_inter(&mut self) {
        let a_prog = NoopProgress;
        let mut a_ps_outer =
            ProgressScope::new(&a_prog, "Connect offset faces by intersection", 100);

        // just for better management and visualization of the progress steps
        // define a nested enum listing all the steps of the current method.
        // OCCT L1301-1311: enum BuildOffsetByInter_PISteps.
        const MAKE_OFFSET_FACES: usize = 0;
        const CONNEX_INT_BY_INT: usize = 1;
        const CONTEXT_INT_BY_INT: usize = 2;
        const INTERSECT_EDGES: usize = 3;
        const COMPLETE_EDGES_INTERSECTION: usize = 4;
        const BUILD_FACES: usize = 5;
        const FILL_HISTORY_FOR_OFFSETS: usize = 6;
        const FILL_HISTORY_FOR_DEEPENINGS: usize = 7;
        const BUILD_OFFSET_BY_INTER_LAST: usize = 8;

        let a_nb_faces = bat::sub_shapes(&self.my_face_comp).len() as f64
            + self.my_analyse.new_faces().len() as f64
            + self.my_faces.extent() as f64;
        let an_offsets_part = (bat::sub_shapes(&self.my_face_comp).len() as f64
            + self.my_analyse.new_faces().len() as f64)
            / a_nb_faces;
        let a_deepenings_part = self.my_faces.extent() as f64 / a_nb_faces;

        let mut a_steps = vec![0f64; BUILD_OFFSET_BY_INTER_LAST];
        {
            // OCCT L1319: aSteps.Init(0).
            for v in a_steps.iter_mut() {
                *v = 0.;
            }

            let is_inter = self.my_join == GeomAbsJoinType::Intersection;
            let a_face_inter = if is_inter { 25. } else { 50. };
            let a_build_faces = if is_inter { 50. } else { 25. };
            a_steps[MAKE_OFFSET_FACES] = 5.;
            a_steps[CONNEX_INT_BY_INT] = a_face_inter * an_offsets_part;
            a_steps[CONTEXT_INT_BY_INT] = a_face_inter * a_deepenings_part;
            a_steps[INTERSECT_EDGES] = 10.;
            a_steps[COMPLETE_EDGES_INTERSECTION] = 5.;
            a_steps[BUILD_FACES] = a_build_faces;
            a_steps[FILL_HISTORY_FOR_OFFSETS] = 5. * an_offsets_part;
            a_steps[FILL_HISTORY_FOR_DEEPENINGS] = 5. * a_deepenings_part;
            normalize_steps(100., &mut a_steps);
        }
        let _ = &a_steps;

        //--------------------------------------------------------
        // Construction of faces parallel to initial faces
        //--------------------------------------------------------
        let mut map_sf = MapSF::new();
        self.make_offset_faces(&mut map_sf);
        if a_ps_outer.user_break() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }
        //--------------------------------------------------------------------
        // MES   : Map of OffsetShape -> Extended Shapes.
        // Build : Map of Initial SS  -> OffsetShape build by Inter.
        //                               can be an edge or a compound of edges
        //---------------------------------------------------------------------
        let mut mes = DataMapOfShapeShape::new();
        let mut build = DataMapOfShapeShape::new();
        let mut failed: Vec<Shape> = Vec::new();
        let side = State::In;
        // OCCT L1339: occ::handle<BRepAlgo_AsDes> AsDes = new BRepAlgo_AsDes().
        let mut as_des = BRepAlgoAsDes::new();

        //-------------------------------------------------------------------
        // Extension of faces and calculation of new edges of intersection.
        //-------------------------------------------------------------------
        let extent_context = self.my_offset > 0.;

        // OCCT L1348: BRepOffset_Inter3d Inter3(AsDes, Side, myTol) — the
        // Handle aliasing resync (architecture difference #39).
        let mut inter3 =
            BRepOffsetInter3d::new(std::mem::take(&mut as_des), side, self.my_tol);
        // Intersection between parallel faces
        inter3.connex_int_by_int(
            &self.my_face_comp,
            &map_sf,
            &self.my_analyse,
            &mut mes,
            &mut build,
            &mut failed,
            &a_ps_outer,
            self.my_is_planar,
        );
        std::mem::swap(&mut as_des, inter3.as_des_mut());
        if a_ps_outer.user_break() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }
        // Intersection with caps.
        inter3.context_int_by_int(
            &self.my_faces,
            extent_context,
            &map_sf,
            &self.my_analyse,
            &mut mes,
            &mut build,
            &mut failed,
            &a_ps_outer,
            self.my_is_planar,
        );
        std::mem::swap(&mut as_des, inter3.as_des_mut());
        if a_ps_outer.user_break() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }

        let mut a_lfaces: Vec<Shape> = Vec::new();
        for exp in bat::explorer(&self.my_face_comp, ShapeType::Face, ShapeType::Shape) {
            a_lfaces.push(exp);
        }
        for it in self.my_analyse.new_faces() {
            a_lfaces.push(it);
        }
        //---------------------------------------------------------------------------------
        // Extension of neighbor edges of new edges and intersection between neighbors.
        //--------------------------------------------------------------------------------
        let mut as_des2d = BRepAlgoAsDes::new();
        self.intersect_edges(
            &a_lfaces,
            &mut map_sf,
            &mut mes,
            &mut build,
            &mut as_des,
            &mut as_des2d,
        );
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //-----------------------------------------------------------
        // Great restriction of new edges and update of AsDes.
        //------------------------------------------ ----------------
        let mut an_edges_origins = DataMapOfShapeListOfShape::new(); // offset edge - initial edges
        let mut new_edges = OcctIndexedShapeMap::new();
        let mut a_e_trim_e_inf = DataMapOfShapeShape::new(); // trimmed - not trimmed edges
        //
        // Map of edges obtained after FACE-FACE (offsetted) intersection.
        // Key1 is edge trimmed by intersection points with other edges;
        // Item is not-trimmed edge.
        if !super::brep_offset_make_offset_d::trim_edges(
            &self.my_face_comp,
            self.my_offset,
            &self.my_analyse,
            &map_sf,
            &mut mes,
            &build,
            &mut as_des,
            &mut as_des2d,
            &mut new_edges,
            &mut a_e_trim_e_inf,
            &mut an_edges_origins,
        ) {
            self.my_error = BRepOffset_Error::CannotTrimEdges;
            return;
        }
        //
        //---------------------------------
        // Intersection 2D on //
        //---------------------------------
        let mut a_dmvv = IndexedDataMapOfShapeListOfShape::new();
        let mut a_faces_origins = DataMapOfShapeShape::new(); // offset face - initial face
        let mut lfe: Vec<Shape> = Vec::new();
        let mut imoe = BRepAlgoImage::new();
        super::brep_offset_make_offset_d::get_enlarged_faces(
            &a_lfaces,
            &map_sf,
            &mes,
            &mut a_faces_origins,
            &mut imoe,
            &mut lfe,
        );
        //
        // OCCT L1454-1467: aPS2d / aPS2dOffsets — the flattened scopes
        // (architecture difference #40).
        for it_lfe in lfe.clone() {
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let nef = it_lfe;
            let a_curr_face_tol = brep_tool_tolerance(&nef);
            BRepOffsetInter2d::compute(
                &mut as_des,
                &nef,
                &indexed_shape_map_view(&new_edges),
                a_curr_face_tol,
                &self.my_edge_int_edges,
                &mut a_dmvv,
                (),
            );
        }
        //----------------------------------------------
        // Intersections 2d on caps.
        //----------------------------------------------
        for i in 1..=self.my_faces.extent() {
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let cork = self.my_faces.at_1(i).clone();
            let a_curr_face_tol = brep_tool_tolerance(&cork);
            BRepOffsetInter2d::compute(
                &mut as_des,
                &cork,
                &indexed_shape_map_view(&new_edges),
                a_curr_face_tol,
                &self.my_edge_int_edges,
                &mut a_dmvv,
                (),
            );
        }
        //
        let _ = BRepOffsetInter2d::fuse_vertices(&a_dmvv, &mut as_des, &mut self.my_image_vv);
        //-------------------------------
        // Unwinding of extended Faces.
        //-------------------------------
        //
        let mut a_mf_done: OcctShapeSet = HashMap::new();
        //
        if (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_inter
            && self.my_is_planar
        {
            // OCCT MakeOffset_1.cxx L9511-9533: the BuildSplitsOfExtendedFaces
            // wrapper — the local BRepOffset_BuildOffsetFaces tool.
            // [INTERFACE NOTE — REPORTED] SetAsDesInfo consumes the OCCT
            // shared-handle carrier Rc<RefCell<BRepAlgoAsDes>>; the owned
            // BRepAlgoAsDes of this class has no Clone, so the handle is the
            // engine-local form (the OCCT handle aliasing — architecture
            // difference #39).
            let as_des_handle =
                std::rc::Rc::new(std::cell::RefCell::new(BRepAlgoAsDes::new()));
            let mut a_bf_tool =
                super::brep_offset_make_offset_1::BRepOffsetBuildOffsetFaces::new();
            a_bf_tool.set_faces(&lfe);
            a_bf_tool.set_as_des_info(&as_des_handle);
            a_bf_tool.set_analysis(&self.my_analyse);
            a_bf_tool.set_edges_origins(&an_edges_origins);
            a_bf_tool.set_faces_origins(&a_faces_origins);
            a_bf_tool.set_inf_edges(&a_e_trim_e_inf);
            let a_prog_1 = rcad_kernel::core::message::NoopProgress;
            let mut a_ps_1 = rcad_kernel::core::message::ProgressScope::new(
                &a_prog_1,
                "BuildSplitsOfExtendedFaces",
                1,
            );
            a_bf_tool.build_splits_of_extended_faces(&mut imoe, &a_ps_1);
            if self.my_error != BRepOffset_Error::NoError {
                return;
            }
            //
            for a_it_lf in lfe.clone() {
                let a_s = a_it_lf;
                set_add(&mut a_mf_done, &a_s);
            }
        } else {
            self.my_make_loops
                .build(&lfe, &mut as_des, &mut imoe, &mut self.my_image_vv);
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
        }
        //---------------------------
        // MAJ SD. for faces //
        //---------------------------
        for it in a_lfaces.clone() {
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let fi = it;
            self.my_init_offset_face.set_root(&fi);
            let mut of = shape_data_map::value(&map_sf, &fi).face();
            if shape_data_map::is_bound(&mes, &of) {
                of = shape_data_map::find(&mes, &of);
                if imoe.has_image(&of) {
                    let lofe = imoe.image(&of);
                    self.my_init_offset_face.bind_list(&fi, &lofe);
                    for it_lf in lofe.clone() {
                        let ofe = it_lf;
                        self.my_image_offset.set_root(&ofe);
                        let mut view: OcctShapeSet = HashMap::new();
                        for exp2 in bat::explorer(
                            &bat::oriented(&ofe, Orientation::Forward),
                            ShapeType::Edge,
                            ShapeType::Shape,
                        ) {
                            let coe = exp2;

                            self.my_as_des.add(&ofe, &coe);
                            if set_add(&mut view, &coe) {
                                if !self.my_as_des.has_descendant(&coe) {
                                    let (cv1, cv2) = top_exp_vertices(&coe);
                                    if !cv1.is_null() {
                                        self.my_as_des
                                            .add(&coe, &bat::oriented(&cv1, Orientation::Forward));
                                    }
                                    if !cv2.is_null() {
                                        self.my_as_des.add(
                                            &coe,
                                            &bat::oriented(&cv2, Orientation::Reversed),
                                        );
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if set_contains(&a_mf_done, &of) {
                        continue;
                    }
                    //
                    self.my_init_offset_face.bind(&fi, &of);
                    self.my_image_offset.set_root(&of);
                    let le = as_des.descendant(&of).to_vec();
                    for it_lf in le.clone() {
                        let oe = it_lf;
                        if imoe.has_image(&oe) {
                            let loe = imoe.image(&oe);
                            for it_loe in &loe {
                                let coe = bat::oriented(it_loe, oe.orientation);
                                self.my_as_des.add(&of, &coe);

                                if !self.my_as_des.has_descendant(&coe) {
                                    let (cv1, cv2) = top_exp_vertices(&coe);
                                    if !cv1.is_null() {
                                        self.my_as_des
                                            .add(&coe, &bat::oriented(&cv1, Orientation::Forward));
                                    }
                                    if !cv2.is_null() {
                                        self.my_as_des.add(
                                            &coe,
                                            &bat::oriented(&cv2, Orientation::Reversed),
                                        );
                                    }
                                }
                            }
                        } else {
                            self.my_as_des.add(&of, &oe);

                            let lv = as_des.descendant(&oe).to_vec();
                            self.my_as_des.add_list(&oe, &lv);
                        }
                    }
                }
            } else {
                self.my_init_offset_face.bind(&fi, &of);
                self.my_image_offset.set_root(&of);
                let mut view: OcctShapeSet = HashMap::new();
                for exp2 in bat::explorer(
                    &bat::oriented(&of, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                ) {
                    let coe = exp2;
                    self.my_as_des.add(&of, &coe);

                    if set_add(&mut view, &coe) {
                        if !self.my_as_des.has_descendant(&coe) {
                            let (cv1, cv2) = top_exp_vertices(&coe);
                            if !cv1.is_null() {
                                self.my_as_des
                                    .add(&coe, &bat::oriented(&cv1, Orientation::Forward));
                            }
                            if !cv2.is_null() {
                                self.my_as_des
                                    .add(&coe, &bat::oriented(&cv2, Orientation::Reversed));
                            }
                        }
                    }
                }
            }
        }
        //  Modified by skv - Tue Mar 15 16:20:43 2005
        // Add methods for supporting history.
        let mut a_map_edges: OcctShapeSet = HashMap::new();

        for it in a_lfaces.clone() {
            let a_face_ref = it;
            for an_edge_ref in bat::explorer(
                &bat::oriented(&a_face_ref, Orientation::Forward),
                ShapeType::Edge,
                ShapeType::Shape,
            ) {
                if set_add(&mut a_map_edges, &an_edge_ref) {
                    self.my_init_offset_edge.set_root(&an_edge_ref);
                    if shape_data_map::is_bound(&build, &an_edge_ref) {
                        let a_new_shape = shape_data_map::find(&build, &an_edge_ref);

                        if a_new_shape.shape_type() == ShapeType::Edge {
                            if imoe.has_image(&a_new_shape) {
                                let a_list_new_e = imoe.image(&a_new_shape);

                                self.my_init_offset_edge.bind_list(&an_edge_ref, &a_list_new_e);
                            } else {
                                self.my_init_offset_edge.bind(&an_edge_ref, &a_new_shape);
                            }
                        } else {
                            // aNewShape != TopAbs_EDGE
                            let mut a_list_new_edge: Vec<Shape> = Vec::new();

                            for exp_c in
                                bat::explorer(&a_new_shape, ShapeType::Edge, ShapeType::Shape)
                            {
                                let a_res_edge = exp_c;

                                if imoe.has_image(&a_res_edge) {
                                    let a_list_new_e = imoe.image(&a_res_edge);

                                    for a_new_e_iter in a_list_new_e {
                                        a_list_new_edge.push(a_new_e_iter);
                                    }
                                } else {
                                    a_list_new_edge.push(a_res_edge);
                                }
                            }

                            self.my_init_offset_edge.bind_list(&an_edge_ref, &a_list_new_edge);
                        }
                    } else {
                        // Free boundary.
                        let mut a_new_edge =
                            shape_data_map::value(&map_sf, &a_face_ref).generated(&an_edge_ref);

                        if shape_data_map::is_bound(&mes, &a_new_edge) {
                            a_new_edge = shape_data_map::find(&mes, &a_new_edge);
                        }

                        if imoe.has_image(&a_new_edge) {
                            let a_list_new_e = imoe.image(&a_new_edge);

                            self.my_init_offset_edge.bind_list(&an_edge_ref, &a_list_new_e);
                        } else {
                            self.my_init_offset_edge.bind(&an_edge_ref, &a_new_edge);
                        }
                    }
                }
            }
        }
        //  Modified by skv - Tue Mar 15 16:20:43 2005

        //---------------------------
        // MAJ SD. for caps
        //---------------------------
        for i in 1..=self.my_faces.extent() {
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let cork = self.my_faces.at_1(i).clone();
            let le = as_des.descendant(&cork).to_vec();
            for it_lf in le {
                let oe = it_lf;
                if imoe.has_image(&oe) {
                    let loe = imoe.image(&oe);
                    for it_loe in &loe {
                        let coe = bat::oriented(it_loe, oe.orientation);
                        self.my_as_des.add(&cork, &coe);

                        if !self.my_as_des.has_descendant(&coe) {
                            let (cv1, cv2) = top_exp_vertices(&coe);
                            if !cv1.is_null() {
                                self.my_as_des
                                    .add(&coe, &bat::oriented(&cv1, Orientation::Forward));
                            }
                            if !cv2.is_null() {
                                self.my_as_des
                                    .add(&coe, &bat::oriented(&cv2, Orientation::Reversed));
                            }
                        }
                    }
                } else {
                    self.my_as_des.add(&cork, &oe);
                    if as_des.has_descendant(&oe) {
                        let lv = as_des.descendant(&oe).to_vec();
                        self.my_as_des.add_list(&oe, &lv);
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::ReplaceRoots (cxx L1825-1871).
    pub(crate) fn replace_roots(&mut self) {
        // Replace the artificial faces and edges in InitOffset maps with the original ones.
        let mut view: OcctShapeSet = HashMap::new();
        for an_exp_f in bat::explorer(&self.my_face_comp, ShapeType::Edge, ShapeType::Shape) {
            let a_f = an_exp_f;
            for an_exp_e in bat::explorer(&a_f, ShapeType::Edge, ShapeType::Shape) {
                let a_e = an_exp_e;
                if !set_add(&mut view, &a_e) {
                    continue;
                }

                let a_f_gen = self.my_analyse.generated(&a_e);
                if a_f_gen.is_null() {
                    continue;
                }

                self.my_init_offset_face.replace_root(&a_f_gen, &a_e);

                // OCCT L1877: for (TopoDS_Iterator itV(aE); ...) — the direct
                // children walk (architecture difference #56).
                for it_v in bat::sub_shapes(&a_e) {
                    let a_v = it_v;
                    if !set_add(&mut view, &a_v) {
                        continue;
                    }

                    let a_e_gen = self.my_analyse.generated(&a_v);
                    if a_e_gen.is_null() {
                        continue;
                    }

                    self.my_init_offset_edge.replace_root(&a_e_gen, &a_v);
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::BuildOffsetByArc (cxx L1893-2152).
    pub(crate) fn build_offset_by_arc(&mut self) {
        let mut done: OcctShapeSet = HashMap::new();
        //--------------------------------------------------------
        // Construction of faces parallel to initial faces
        //--------------------------------------------------------
        let mut map_sf = MapSF::new();
        self.make_offset_faces(&mut map_sf);
        if self.my_error != BRepOffset_Error::NoError {
            return;
        }
        //--------------------------------------------------------
        // Construction of tubes on edge.
        //--------------------------------------------------------
        let ot = if self.my_offset < 0. {
            ChFiDS_TypeOfConcavity::Concave
        } else {
            ChFiDS_TypeOfConcavity::Convex
        };

        for exp in bat::explorer(&self.my_face_comp, ShapeType::Edge, ShapeType::Shape) {
            let e = exp;
            if set_add(&mut done, &e) {
                let anc = self.my_analyse.ancestors(&e);
                if anc.len() == 2 {
                    let l = self.my_analyse.type_(&e);
                    if !l.is_empty() && l[0].type_of() == ot {
                        let mut cur_offset = self.my_offset;
                        if shape_data_map::is_bound(&self.my_face_offset, &anc[0]) {
                            cur_offset = shape_data_map::find(&self.my_face_offset, &anc[0]);
                        }
                        let e_on1 = shape_data_map::value(&map_sf, &anc[0]).generated(&e);
                        let e_on2 = shape_data_map::value(&map_sf, &anc[1]).generated(&e);
                        // find if exits tangent edges in the original shape
                        let (v1f, v1l) = top_exp_vertices(&e);
                        let mut tang_e: Vec<Shape> = Vec::new();
                        self.my_analyse.tangent_edges(&e, &v1f, &mut tang_e);
                        // find if the pipe on the tangent edges are soon created.
                        let mut e1f = Shape::null();
                        let mut find = false;
                        for itl in &tang_e {
                            if find {
                                break;
                            }
                            if shape_data_map::is_bound(&map_sf, itl) {
                                e1f = shape_data_map::value(&map_sf, itl).generated(&v1f);
                                find = true;
                            }
                        }
                        // OCCT L2020: TangE.Clear().
                        let mut tang_e: Vec<Shape> = Vec::new();
                        self.my_analyse.tangent_edges(&e, &v1l, &mut tang_e);
                        // find if the pipe on the tangent edges are soon created.
                        let mut e1l = Shape::null();
                        let mut find = false;
                        for itl in &tang_e {
                            if find {
                                break;
                            }
                            if shape_data_map::is_bound(&map_sf, itl) {
                                e1l = shape_data_map::value(&map_sf, itl).generated(&v1l);
                                find = true;
                            }
                        }
                        // OCCT L2032: BRepOffset_Offset OF(E, EOn1, EOn2,
                        // CurOffset, E1f, E1l) — the (Path, Edge1, Edge2,
                        // Offset, FirstEdge, LastEdge) ctor with the OCCT
                        // defaults Polynomial = false, Tol = 1.0e-4, Conti =
                        // GeomAbs_C1.
                        let of = BRepOffsetOffset::with_path_first_last(
                            &e,
                            &e_on1,
                            &e_on2,
                            cur_offset,
                            &e1f,
                            &e1l,
                            false,
                            1.0e-4,
                            GeomAbsShapeKind::C1,
                        );
                        shape_data_map::bind(&mut map_sf, &e, of);
                    }
                } else {
                    // ----------------------
                    // free border.
                    // ----------------------
                    let e_on1 = shape_data_map::value(&map_sf, &anc[0]).generated(&e);
                    self.my_init_offset_edge.set_root(&e); // skv: supporting history.
                    self.my_init_offset_edge.bind(&e, &e_on1);
                }
            }
        }

        //--------------------------------------------------------
        // Construction of spheres on vertex.
        //--------------------------------------------------------
        done.clear();
        for exp in bat::explorer(&self.my_face_comp, ShapeType::Vertex, ShapeType::Shape) {
            let v = exp;
            if set_add(&mut done, &v) {
                let la = self.my_analyse.ancestors(&v);
                let mut le: Vec<Shape> = Vec::new();
                analyse_edges(&self.my_analyse, &v, ot, &mut le);

                if le.len() >= 3 && le.len() == la.len() {
                    let mut loe: Vec<Shape> = Vec::new();
                    //--------------------------------------------------------
                    // Return connected edges on tubes.
                    //--------------------------------------------------------
                    for it in &le {
                        loe.push(bat::reversed(
                            &shape_data_map::value(&map_sf, it).generated(&v),
                        ));
                    }
                    //----------------------
                    // construction sphere.
                    //-----------------------
                    let lla = self.my_analyse.ancestors(&la[0]);
                    let ff = lla[0].clone();
                    let mut cur_offset = self.my_offset;
                    if shape_data_map::is_bound(&self.my_face_offset, &ff) {
                        cur_offset = shape_data_map::find(&self.my_face_offset, &ff);
                    }

                    // OCCT L2066: BRepOffset_Offset OF(V, LOE, CurOffset) —
                    // the OCCT defaults Polynomial = false, Tol = 1.0e-4,
                    // Conti = GeomAbs_C1.
                    let of = BRepOffsetOffset::with_vertex(
                        &v,
                        &loe,
                        cur_offset,
                        false,
                        1.0e-4,
                        GeomAbsShapeKind::C1,
                    );
                    shape_data_map::bind(&mut map_sf, &v, of);
                }
                //--------------------------------------------------------------
                // Particular processing if V is at least a free border.
                //-------------------------------------------------------------
                let mut lbf: Vec<Shape> = Vec::new();
                analyse_edges(
                    &self.my_analyse,
                    &v,
                    ChFiDS_TypeOfConcavity::FreeBound,
                    &mut lbf,
                );
                if !lbf.is_empty() {
                    let mut first = true;
                    for it in &le {
                        if first {
                            self.my_init_offset_edge.set_root(&v); // skv: supporting history.
                            let a_gen = shape_data_map::value(&map_sf, it).generated(&v);
                            self.my_init_offset_edge.bind(&v, &a_gen);
                            first = false;
                        } else {
                            let a_gen = shape_data_map::value(&map_sf, it).generated(&v);
                            self.my_init_offset_edge.add(&v, &a_gen);
                        }
                    }
                }
            }
        }

        //------------------------------------------------------------
        // Extension of parallel faces to the context.
        // Extended faces are ordered in DS and removed from MapSF.
        //------------------------------------------------------------
        if !self.my_faces.is_empty() {
            self.to_context(&mut map_sf);
        }

        //------------------------------------------------------
        // MAJ SD.
        //------------------------------------------------------
        let rt = if self.my_offset < 0. {
            ChFiDS_TypeOfConcavity::Convex
        } else {
            ChFiDS_TypeOfConcavity::Concave
        };
        // OCCT L2131: NCollection_DataMap::Iterator It(MapSF) — the rcad
        // HashMap order (architecture difference #55); the rcad borrow-split
        // iteration (the BRepOffsetOffset values are borrowed, not cloned —
        // the brep_offset_offset_b.rs carrier has no Clone).
        for (_k, (si, sf)) in map_sf.iter() {
            if sf.status() == BRepOffsetStatus::Reversed
                || sf.status() == BRepOffsetStatus::Degenerated
            {
                //------------------------------------------------
                // Degenerated or returned faces are not stored.
                //------------------------------------------------
                continue;
            }

            let of = sf.face();
            self.my_init_offset_face.bind(si, &of);
            self.my_init_offset_face.set_root(si); // Initial<-> Offset
            self.my_image_offset.set_root(&of); // FaceOffset root of images

            if si.shape_type() == ShapeType::Face {
                for exp in bat::explorer(
                    &bat::oriented(si, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                ) {
                    //--------------------------------------------------------------------
                    // To each face are associatedthe edges that restrict that
                    // The edges that do not generate tubes or are not tangent
                    // to two faces are removed.
                    //--------------------------------------------------------------------
                    let e = exp;
                    let l = self.my_analyse.type_(&e);
                    if !l.is_empty() && l[0].type_of() != rt {
                        let oo = e.orientation;
                        let a_gen = sf.generated(&e);
                        let oe = bat::oriented(&a_gen, oo);
                        self.my_as_des.add(&of, &oe);
                    }
                }
            } else {
                for exp in bat::explorer(
                    &bat::oriented(&of, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                ) {
                    self.my_as_des.add(&of, &exp);
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::ToContext (cxx L2171-2340).
    pub(crate) fn to_context(&mut self, map_sf: &mut MapSF) {
        let mut created = DataMapOfShapeShape::new();
        let mut mef = DataMapOfShapeShape::new();
        let mut faces_to_build = OcctIndexedShapeMap::new();

        let side = State::Out;

        //--------------------------------------------------------
        // Determine the edges and faces reconstructed by
        // intersection.
        //---------------------------------------------------------
        for j in 1..=self.my_faces.extent() {
            let cf = self.my_faces.at_1(j).clone();
            for exp in bat::explorer(
                &bat::oriented(&cf, Orientation::Forward),
                ShapeType::Edge,
                ShapeType::Shape,
            ) {
                let e = exp;
                if self.my_analyse.has_ancestor(&e) {
                    let lea = self.my_analyse.ancestors(&e);
                    for itl in &lea {
                        let of = shape_data_map::value(map_sf, itl);
                        let a_gen = of.generated(&e);
                        faces_to_build.add(itl);
                        shape_data_map::bind(&mut mef, &a_gen, cf.clone());
                    }
                    // OCCT L2258-2268: TopoDS_Vertex V[2]; TopExp::Vertices;
                    // the ancestor walk of each extremity.
                    let (v0, v1) = top_exp_vertices(&e);
                    let vs = [v0, v1];
                    for v in &vs {
                        let lva = self.my_analyse.ancestors(v);
                        for itl in &lva {
                            if shape_data_map::is_bound(map_sf, itl) {
                                let of = shape_data_map::value(map_sf, itl);
                                let a_gen = of.generated(v);
                                faces_to_build.add(itl);
                                shape_data_map::bind(&mut mef, &a_gen, cf.clone());
                            }
                        }
                    }
                }
            }
        }
        //---------------------------
        // Reconstruction of faces.
        //---------------------------
        let rt = if self.my_offset < 0. {
            ChFiDS_TypeOfConcavity::Convex
        } else {
            ChFiDS_TypeOfConcavity::Concave
        };

        for j in 1..=faces_to_build.extent() {
            let s = faces_to_build.at_1(j).clone();
            // OCCT L2297: BOF = MapSF(S) — the rcad borrow-split form: the
            // (face, edge -> generated edge) pairs are captured before the
            // UnBind (the BRepOffsetOffset carrier has no Clone).
            let f = shape_data_map::value(map_sf, &s).face();
            let edge_gens: Vec<(Shape, Shape)> = bat::explorer(
                &bat::oriented(&s, Orientation::Forward),
                ShapeType::Edge,
                ShapeType::Shape,
            )
            .iter()
            .map(|e| (e.clone(), shape_data_map::value(map_sf, &s).generated(e)))
            .collect();
            let mut nf = Shape::null();
            super::brep_offset_tool_c::extent_face(
                &f,
                &mut created,
                &mut mef,
                side,
                self.my_tol,
                &mut nf,
            );
            shape_data_map::un_bind(map_sf, &s);
            //--------------
            // MAJ SD.
            //--------------
            self.my_init_offset_face.bind(&s, &nf);
            self.my_init_offset_face.set_root(&s); // Initial<-> Offset
            self.my_image_offset.set_root(&nf);

            if s.shape_type() == ShapeType::Face {
                for (e, oe) in &edge_gens {
                    let e = e.clone();
                    let mut oe = oe.clone();
                    let l = self.my_analyse.type_(&e);
                    let or = e.orientation;
                    oe.orientation = or;
                    if !l.is_empty() && l[0].type_of() != rt {
                        if shape_data_map::is_bound(&created, &oe) {
                            let mut ne = shape_data_map::find(&created, &oe);
                            if ne.orientation == Orientation::Reversed {
                                ne.orientation = bat::top_abs_reverse(or);
                            } else {
                                ne.orientation = or;
                            }
                            self.my_as_des.add(&nf, &ne);
                        } else {
                            self.my_as_des.add(&nf, &oe);
                        }
                    }
                }
            } else {
                //------------------
                // Tube
                //---------------------
                for exp in bat::explorer(
                    &bat::oriented(&nf, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                ) {
                    self.my_as_des.add(&nf, &exp);
                }
            }
            // OCCT L2339: MapSF.UnBind(S) — the OCCT double-unbind form.
            shape_data_map::un_bind(map_sf, &s);
        }

        //------------------
        // MAJ free borders
        //------------------
        let created_items: Vec<(Shape, Shape)> = created
            .values()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (oe, mut ne) in created_items {
            if self.my_init_offset_edge.is_image(&oe) {
                let e = self.my_init_offset_edge.image_from(&oe).clone();
                let or = self.my_init_offset_edge.image(&e)[0].orientation;
                if ne.orientation == Orientation::Reversed {
                    ne.orientation = bat::top_abs_reverse(or);
                } else {
                    ne.orientation = or;
                }
                self.my_init_offset_edge.remove(&oe);
                self.my_init_offset_edge.bind(&e, &ne);
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::UpdateFaceOffset (cxx L2342-2394).
    pub(crate) fn update_face_offset(&mut self) {
        let mut m: OcctShapeSet = HashMap::new();
        // OCCT L2350: CopiedMap.Assign(myFaceOffset).
        let copied_map_items: Vec<(Shape, f64)> = self
            .my_face_offset
            .values()
            .map(|(k, v)| (k.clone(), *v))
            .collect();

        let rt = if self.my_offset < 0. {
            ChFiDS_TypeOfConcavity::Concave
        } else {
            ChFiDS_TypeOfConcavity::Convex
        };

        for (f, cur_offset) in copied_map_items {
            if !set_add(&mut m, &f) {
                continue;
            }
            let mut build = rcad_kernel::topo::topods::BRepBuilder::new();
            let mut co = build.make_compound(&mut self.my_brep, vec![]);
            let mut dummy: OcctShapeSet = HashMap::new();
            build.add_to_compound(&mut self.my_brep, co.clone(), f.clone());
            if self.my_join == GeomAbsJoinType::Arc {
                analyse_add_faces_rt(
                    &self.my_analyse,
                    &f,
                    &mut co,
                    &mut dummy,
                    ChFiDS_TypeOfConcavity::Tangential,
                    rt,
                );
            } else {
                analyse_add_faces(
                    &self.my_analyse,
                    &f,
                    &mut co,
                    &mut dummy,
                    ChFiDS_TypeOfConcavity::Tangential,
                );
            }

            for exp in bat::explorer(&co, ShapeType::Face, ShapeType::Shape) {
                let ff = exp;
                if !set_add(&mut m, &ff) {
                    continue;
                }
                if shape_data_map::is_bound(&self.my_face_offset, &ff) {
                    shape_data_map::un_bind(&mut self.my_face_offset, &ff);
                }
                shape_data_map::bind(&mut self.my_face_offset, &ff, cur_offset);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Local helpers of module b.
// ---------------------------------------------------------------------------

/// OCCT array size (0, PIOperation_Last - 1) — 8 entries.
const PIO_OPERATION_LAST_ARRAY: usize = 8;

/// The OcctIndexedShapeMap -> IndexedShapeMap view (the two dependency
/// modules carry the two rcad forms of the OCCT NCollection_IndexedMap;
/// architecture difference #38).
pub(crate) fn indexed_shape_map_view(m: &OcctIndexedShapeMap) -> IndexedShapeMap {
    let mut out = indexmap::IndexMap::new();
    for s in m.iter() {
        out.insert(shape_key(s), s.clone());
    }
    out
}

/// The OcctIndexedShapeMap boundary clone (the type carries no Clone
/// derive; the rebuild keeps the OCCT insertion order — architecture
/// difference #39).
pub(crate) fn indexed_shape_map_clone(m: &OcctIndexedShapeMap) -> OcctIndexedShapeMap {
    let mut out = OcctIndexedShapeMap::new();
    for s in m.iter() {
        out.add(s);
    }
    out
}
