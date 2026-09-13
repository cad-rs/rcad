//! OCCT BRepFill_Evolved — 1:1 translation (part 3: the remaining private
//! methods).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_Evolved.cxx
//! - PlanarPerform           L1267-1414
//! - VerticalPerform         L1418-1520
//! - PrepareProfile          L1556-1712
//! - PrepareSpine            L1716-1779
//! - Add                     L1844-1955
//! - Transfert               L1966-2053
//! - AddTopAndBottom         L2078-2222
//! - MakeSolid               L2229-2262
//! - MakePipe                L2266-2325
//! - MakeRevol               L2329-2400
//! - FindLocation            L2404-2437
//! - TransformInitWork       L2441-2445
//! - ContinuityOnOffsetEdge  L2452-2544

use glam::{DAffine3, DVec2, DVec3};

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{CurveEval, Line2d, Plane, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType};

use crate::brep_algo::face_restrictor::BRepAlgoFaceRestrictor;
use crate::brep_algo::r#loop::BRepAlgoLoop;
use crate::brep_algo::tool::{
    builder_set_closed, explorer, shape_is_closed, top_abs_reverse,
};
// OCCT BRepLib_FindSurface (TKTopAlgo/BRepLib_FindSurface.hxx / .cxx) — the
// single 1:1 body lives in crate::topalgo::brep_lib_find_surface; the
// brep_fill_axe GAP carrier was retired.
use crate::topalgo::brep_lib_find_surface::BRepLibFindSurface;
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::brep_fill::brep_fill_evolved::{
    builder_add_compound_shape, builder_add_face_wire, builder_add_solid_shell,
    builder_add_wire_edge, builder_make_compound, builder_make_face_surface, builder_make_solid,
    builder_make_wire, brep_fill_confusion, brep_lib_build_curves3d, brep_lib_same_parameter,
    brep_tool_degenerated, brep_tool_pnt, brep_tool_surface, edge_vertices,
    location_shape_moved, location_shape_set, wire_edges, BRepFillEvolved,
    BRepMAT2dBisectingLocusCarrier, BRepMAT2dLinkTopoBiloCarrier, BRepSweepPrismCarrier,
    BRepSweepRevolCarrier, DataMapOfShapeItem, DataMapOfShapeListOfShape,
};
use crate::brep_fill::brep_fill_offset_ancestors::BRepFillOffsetAncestors;
use crate::brep_fill::brep_fill_pipe::{top_exp_vertices, BRepFillPipe, GeomFillTrihedron};
use crate::brep_fill::brep_fill_trim_edge_tool::GeomAbsJoinType;
use crate::brep_fill::generator::{shape_key, shape_oriented, shape_reversed, ShapeKey};
use crate::brep_fill::offset_wire::BRepFillOffsetWire;
use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;
use crate::topalgo::brep_tools_modification::BRepToolsTrsfModification;
use crate::topalgo::brep_tools_modifier::BRepToolsModifier;
// OCCT BRepTools_Quilt (TKBRep/BRepTools/BRepTools_Quilt.hxx / .cxx) — the
// single 1:1 body lives in crate::topalgo::brep_tools_quilt; the local
// BRepToolsQuiltCarrier GAP stand-in was retired.
use crate::topalgo::brep_tools_quilt::BRepToolsQuilt;

/// OCCT TopAbs_IN (the SolidClassifier state code).
const TOPABS_IN: u8 = 0;

impl BRepFillEvolved {
    // -------------------------------------------------------------------
    // PlanarPerform (cxx L1267-1414)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::PlanarPerform(Sp, Pr, Locus, Link, Join)
    /// (L1267-1414).
    pub(super) fn planar_perform(
        &mut self,
        sp: &Shape,
        pr: &Shape,
        locus: &BRepMAT2dBisectingLocusCarrier,
        link: &mut BRepMAT2dLinkTopoBiloCarrier,
        join: GeomAbsJoinType,
    ) {
        // OCCT L1273-1277.
        let a_local_shape_oriented = shape_oriented(sp, Orientation::Forward);
        self.my_spine = a_local_shape_oriented;
        self.my_profile = pr.clone();
        self.my_map.clear();

        // OCCT L1279-1280.
        self.my_shape = builder_make_compound();

        // OCCT L1284-1288.
        let mut map_vp: DataMapOfShapeItem<Shape> = DataMapOfShapeItem::new();

        // The BRepFill_OffsetWire translation is pool-bound (its
        // PerformWithBiLo / GeneratedShapes signatures take the BRep) — the
        // local pool carries the offset construction.
        let mut paral_brep = BRep::new();

        // OCCT L1290-1413: for (ProfExp.Init(myProfile); ProfExp.More(); ...).
        for e in wire_edges(&self.my_profile.clone()) {
            let mut fr = BRepAlgoFaceRestrictor::new();
            let mut off_anc = BRepFillOffsetAncestors::new();

            // OCCT L1296-1302.
            let (v0, v1) = edge_vertices(&e);
            let v = [v0, v1];
            let alt = crate::brep_fill::brep_fill_evolved_d::altitud(&v[0]);
            let offset = [
                crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v[0]),
                crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v[1]),
            ];
            let is_min_v1 = offset[0] < offset[1];

            // OCCT L1304-1373: for (i = 0; i <= 1; i++).
            for i in 0..=1usize {
                if !map_vp.is_bound(&v[i]) {
                    //------------------------------------------------
                    // Calculate parallel lines corresponding to vertices.
                    //------------------------------------------------
                    // OCCT L1311-1313.
                    let mut paral = BRepFillOffsetWire::new();
                    // The offset_wire.rs translation binds its
                    // PerformWithBiLo to its local placeholder types; the
                    // forwarded call carries the GAP panic (the Locus /
                    // Link of the OCCT call are the carriers of this
                    // translation).
                    let ph_locus = crate::brep_fill::offset_wire::BRepMAT2dBisectingLocus;
                    let mut ph_link = crate::brep_fill::offset_wire::BRepMAT2dLinkTopoBilo;
                    paral.perform_with_bilo(
                        &mut paral_brep,
                        &self.my_spine.clone(),
                        offset[i],
                        &ph_locus,
                        &mut ph_link,
                        join,
                        alt,
                    );
                    off_anc.perform(&mut paral, &paral_brep);
                    let paral_shape = paral.shape().clone();
                    map_vp.bind(&v[i], paral_shape.clone());

                    //-----------------------------
                    // Update myMap (.)(V[i])
                    //-----------------------------
                    // OCCT L1318-1330.
                    for exp_current in explorer(&paral_shape, ShapeType::Edge, ShapeType::Shape) {
                        let wc = exp_current;
                        let gs = off_anc.ancestor(&wc);
                        if !self.my_map.is_bound(&gs) {
                            self.my_map.bind(&gs, DataMapOfShapeListOfShape::new());
                        }
                        if !self.my_map.find(&gs).is_bound(&v[i]) {
                            let generated = paral.generated_shapes(&paral_brep, &gs);
                            self.my_map.find_mut(&gs).bind(&v[i], generated);
                        }
                    }
                }
                let rest = map_vp.find(&v[i]).clone();

                // OCCT L1334-1338.
                let to_reverse = (is_min_v1 && (i == 1)) || (!is_min_v1 && (i == 0));

                // OCCT L1340-1373.
                if !rest.is_null() {
                    if rest.shape_type() == ShapeType::Wire {
                        if to_reverse {
                            let a_local_shape = shape_reversed(&rest);
                            fr.add(&a_local_shape);
                        } else {
                            fr.add(&rest);
                        }
                    } else {
                        for exp_current in explorer(&rest, ShapeType::Wire, ShapeType::Shape) {
                            let wcop = exp_current;
                            if to_reverse {
                                let a_local_shape = shape_reversed(&wcop);
                                fr.add(&a_local_shape);
                            } else {
                                fr.add(&wcop);
                            }
                        }
                    }
                }
            }

            //----------------------------------------------------
            // Construction of faces limited by parallels.
            // - set to the height of the support face.
            //----------------------------------------------------
            // OCCT L1380-1386.
            let t = DAffine3::from_translation(DVec3::new(0.0, 0.0, alt));
            let lt = self.my_locations.intern(t);
            let a_local_shape = location_shape_moved(&mut self.my_locations, &self.my_spine, lt);
            fr.init(&a_local_shape, false, false);
            fr.perform();

            // OCCT L1388-1412: for (; FR.More(); FR.Next()).
            while fr.more() {
                let f = fr.current();
                builder_add_compound_shape(&mut self.my_shape, &f);
                //---------------------------------------
                // Update myMap(.)(E)
                //---------------------------------------
                for exp_current in explorer(&f, ShapeType::Edge, ShapeType::Shape) {
                    let ce = exp_current;
                    if off_anc.has_ancestor(&ce) {
                        let init_e = off_anc.ancestor(&ce);
                        if !self.my_map.is_bound(&init_e) {
                            self.my_map.bind(&init_e, DataMapOfShapeListOfShape::new());
                        }
                        if !self.my_map.find(&init_e).is_bound(&e) {
                            self.my_map.find_mut(&init_e).bind(&e, Vec::new());
                        }
                        self.my_map.find_mut(&init_e).find_mut(&e).push(f.clone());
                    }
                }
                fr.next();
            }
        } // End loop on profile.
    }

    // -------------------------------------------------------------------
    // VerticalPerform (cxx L1418-1520)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::VerticalPerform(Sp, Pr, Locus, Link, Join)
    /// (L1418-1520).
    pub(super) fn vertical_perform(
        &mut self,
        sp: &Shape,
        pr: &Shape,
        locus: &BRepMAT2dBisectingLocusCarrier,
        link: &mut BRepMAT2dLinkTopoBiloCarrier,
        join: GeomAbsJoinType,
    ) {
        // OCCT L1424-1428.
        let a_local_shape = shape_oriented(sp, Orientation::Forward);
        self.my_spine = a_local_shape;
        self.my_profile = pr.clone();
        self.my_map.clear();

        // OCCT L1430-1431.
        self.my_shape = builder_make_compound();

        // OCCT L1433-1443.
        let mut first = true;
        let mut base = Shape::null();

        // The BRepFill_OffsetWire pool (the planar_perform precedent).
        let mut paral_brep = BRep::new();

        // OCCT L1444-1519.
        for e in wire_edges(&self.my_profile.clone()) {
            let (v1, v2) = edge_vertices(&e);
            let alt1 = crate::brep_fill::brep_fill_evolved_d::altitud(&v1);
            let alt2 = crate::brep_fill::brep_fill_evolved_d::altitud(&v2);

            if first {
                // OCCT L1453-1457.
                let mut offset = crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v1);
                if offset.abs() < brep_fill_confusion() {
                    offset = 0.0;
                }
                // OCCT L1458-1460.
                let mut paral = BRepFillOffsetWire::new();
                // The perform_with_bilo placeholder-type forwarding (the
                // planar_perform precedent).
                let ph_locus = crate::brep_fill::offset_wire::BRepMAT2dBisectingLocus;
                let mut ph_link = crate::brep_fill::offset_wire::BRepMAT2dLinkTopoBilo;
                paral.perform_with_bilo(
                    &mut paral_brep,
                    &self.my_spine.clone(),
                    offset,
                    &ph_locus,
                    &mut ph_link,
                    join,
                    alt1,
                );
                let mut off_anc = BRepFillOffsetAncestors::new();
                off_anc.perform(&mut paral, &paral_brep);
                base = paral.shape().clone();

                // MAJ myMap
                // OCCT L1463-1477.
                for exp_current in explorer(&base.clone(), ShapeType::Edge, ShapeType::Shape) {
                    let an_edge = exp_current;
                    let ae = off_anc.ancestor(&an_edge);
                    if !self.my_map.is_bound(&ae) {
                        self.my_map.bind(&ae, DataMapOfShapeListOfShape::new());
                    }
                    if !self.my_map.find(&ae).is_bound(&v1) {
                        self.my_map.find_mut(&ae).bind(&v1, Vec::new());
                    }
                    self.my_map.find_mut(&ae).find_mut(&v1).push(an_edge.clone());
                }
                first = false;
            }

            // OCCT L1481.
            let ps = BRepSweepPrismCarrier::new(
                &base.clone(),
                DVec3::new(0.0, 0.0, alt2 - alt1),
                false,
            );

            // OCCT L1483.
            base = ps.last_shape();

            // OCCT L1485-1488.
            for exp_current in explorer(&ps.shape(), ShapeType::Face, ShapeType::Shape) {
                builder_add_compound_shape(&mut self.my_shape, &exp_current);
            }

            // MAJ myMap
            // OCCT L1491-1518.
            let map_keys: Vec<Shape> = self.my_map.iter().map(|(ks, _)| ks.clone()).collect();
            for key in map_keys {
                let lof = self.my_map.find(&key).find(&v1).clone();
                if !self.my_map.find(&key).is_bound(&v2) {
                    self.my_map.find_mut(&key).bind(&v2, Vec::new());
                }
                if !self.my_map.find(&key).is_bound(&e) {
                    self.my_map.find_mut(&key).bind(&e, Vec::new());
                }

                for os in &lof {
                    self.my_map
                        .find_mut(&key)
                        .find_mut(&v2)
                        .push(ps.last_shape_of(os));
                    self.my_map.find_mut(&key).find_mut(&e).push(ps.shape_of(os));
                }
            }
        }
    }

    // -------------------------------------------------------------------
    // PrepareProfile (cxx L1556-1712)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::PrepareProfile(WorkProf, MapProf) const
    /// (L1556-1712) — projection of the profile on the working plane, cut
    /// at the extrema of the distance to Oz, isolation of the vertical and
    /// horizontal parts, reconstruction of the wires starting from the cut
    /// edges.
    pub(super) fn prepare_profile(
        &self,
        work_prof: &mut Vec<Shape>,
        map_prof: &mut DataMapOfShapeItem<Shape>,
    ) {
        // Supposedly the profile is located so that the only transformation
        // to be carried out is a projection on plane yOz.

        // initialise the projection Plane and the Line to evaluate the extrema.
        // OCCT L1564-1565: Plane = Geom_Plane(gp_Ax3(gp::YOZ()));
        // Line = Geom2d_Line(gp::OY2d()).
        let plane = Plane::new(DVec3::ZERO, DVec3::X);
        let line = Line2d::new(DVec2::ZERO, DVec2::new(0.0, 1.0));

        // Map initial vertex -> projected vertex.
        // OCCT L1568.
        let mut map_ver_ref_moved: DataMapOfShapeItem<Shape> = DataMapOfShapeItem::new();

        // OCCT L1570-1575: V1, V2, VRef1, VRef2; W; B; WP; B.MakeWire(W);
        // WP.Append(W).  The current wire is the last WP entry (the shared
        // TShape handle — the Arc::make_mut sharing rule).
        let mut wp: Vec<Shape> = Vec::new();
        wp.push(builder_make_wire());

        // OCCT L1577-1640: BRepTools_WireExplorer Exp(myProfile);
        // while (Exp.More()).
        let profile_edges = wire_edges(&self.my_profile.clone());
        for idx in 0..profile_edges.len() {
            let mut cuts: Vec<Shape> = Vec::new();
            let mut new_wire = false;
            let e = profile_edges[idx].clone();

            // Cut of the edge.
            // OCCT L1586.
            crate::brep_fill::brep_fill_evolved_d::cut_edge_prof(
                &e,
                &plane,
                &line,
                &mut cuts,
                &mut map_ver_ref_moved,
            );

            // OCCT L1588: EdgeVertices(E, VRef1, VRef2) — the OCCT result
            // is not consumed afterwards (the source as written).
            let (_v_ref1, _v_ref2) = edge_vertices(&e);

            // OCCT L1590-1633.
            if cuts.is_empty() {
                // Neither extrema nor intersections nor vertices on the axis.
                builder_add_wire_edge(wp.last_mut().expect("WP"), &e);
                map_prof.bind(&e, e.clone());
            } else {
                while !cuts.is_empty() {
                    let ne = cuts.first().cloned().expect("Cuts.First()");
                    map_prof.bind(&ne, e.clone());
                    let (v1, v2) = edge_vertices(&ne);
                    if !map_prof.is_bound(&v1) {
                        map_prof.bind(&v1, e.clone());
                    }
                    if !map_prof.is_bound(&v2) {
                        map_prof.bind(&v2, e.clone());
                    }

                    builder_add_wire_edge(wp.last_mut().expect("WP"), &ne);
                    cuts.remove(0);

                    if crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v2) < brep_fill_confusion()
                        && crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v1) > brep_fill_confusion()
                    {
                        // NE ends on axis OZ => new wire
                        if cuts.is_empty() {
                            // last part of the current edge
                            // If it is not the last edge of myProfile
                            // create a new wire.
                            new_wire = true;
                        } else {
                            // New wire.
                            wp.push(builder_make_wire());
                        }
                    }
                }
            }
            let more = idx + 1 < profile_edges.len();
            if more && new_wire {
                // OCCT L1635-1639.
                wp.push(builder_make_wire());
            }
        }

        // In the list of Wires, find edges generating plane or vertical vevo.
        // OCCT L1643-1703.
        for cur_w in wp {
            let mut ya_modif = false;
            for ew_current in wire_edges(&cur_w) {
                let ee = ew_current;
                if crate::brep_fill::brep_fill_evolved_d::is_vertical(&ee)
                    || crate::brep_fill::brep_fill_evolved_d::is_planar(&ee)
                {
                    ya_modif = true;
                    break;
                }
            }

            if ya_modif {
                // Status = 0 for the beginning
                //          3 vertical
                //          2 horizontal
                //          1 other
                let mut status = 0;

                for ew_current in wire_edges(&cur_w) {
                    let ee = ew_current;
                    if crate::brep_fill::brep_fill_evolved_d::is_vertical(&ee) {
                        if status != 3 {
                            work_prof.push(builder_make_wire());
                            status = 3;
                        }
                    } else if crate::brep_fill::brep_fill_evolved_d::is_planar(&ee) {
                        if status != 2 {
                            work_prof.push(builder_make_wire());
                            status = 2;
                        }
                    } else if status != 1 {
                        work_prof.push(builder_make_wire());
                        status = 1;
                    }
                    builder_add_wire_edge(work_prof.last_mut().expect("WorkProf"), &ee);
                }
            } else {
                work_prof.push(cur_w);
            }
        }

        // connect vertices modified in MapProf;
        // OCCT L1706-1711.
        let pairs: Vec<(Shape, Shape)> = map_ver_ref_moved
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (k, v) in pairs {
            // OCCT L1710: MapProf.Bind(gilbert.Value(), gilbert.Key()).
            map_prof.bind(&v, k);
        }
    }

    // -------------------------------------------------------------------
    // PrepareSpine (cxx L1716-1779)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::PrepareSpine(WorkSpine, MapSpine) const
    /// (L1716-1779) — cut of the spine edges at the extrema of curvature
    /// and at the inflexion points.
    pub(super) fn prepare_spine(
        &self,
        work_spine: &mut Shape,
        map_spine: &mut DataMapOfShapeItem<Shape>,
    ) {
        // OCCT L1726-1728.
        let s = brep_tool_surface(&self.my_spine.clone()).expect("BRep_Tool::Surface");
        let tol_f = brep_tool_tolerance(&self.my_spine.clone());
        *work_spine = builder_make_face_surface(&s, tol_f);

        // OCCT L1730-1775: TopoDS_Iterator IteF(mySpine) — the face wires.
        for ite_f_value in face_wires(&self.my_spine.clone()) {
            let mut nw = builder_make_wire();
            let is_closed = shape_is_closed(&ite_f_value);

            for ite_w_value in wire_edges(&ite_f_value) {
                let e = ite_w_value;
                let (v1, v2) = edge_vertices(&e);
                map_spine.bind(&v1, v1.clone());
                map_spine.bind(&v2, v2.clone());
                let mut cuts: Vec<Shape> = Vec::new();

                // Cut
                // OCCT L1747.
                crate::brep_fill::brep_fill_evolved_d::cut_edge(&e, &self.my_spine.clone(), &mut cuts);

                // OCCT L1749-1771.
                if cuts.is_empty() {
                    builder_add_wire_edge(&mut nw, &e);
                    map_spine.bind(&e, e.clone());
                } else {
                    for ite_cuts_value in &cuts {
                        let ne = ite_cuts_value.clone();
                        builder_add_wire_edge(&mut nw, &ne);
                        map_spine.bind(&ne, e.clone());
                        let (v1, v2) = edge_vertices(&ne);
                        if !map_spine.is_bound(&v1) {
                            map_spine.bind(&v1, e.clone());
                        }
                        if !map_spine.is_bound(&v2) {
                            map_spine.bind(&v2, e.clone());
                        }
                    }
                }
            }
            // OCCT L1773: NW.Closed(IsClosed).
            builder_set_closed(&mut nw, is_closed);
            builder_add_face_wire(work_spine, &nw);
        }

        // Construct curves 3D of the spine
        // OCCT L1778.
        brep_lib_build_curves3d(&work_spine.clone());
    }

    // -------------------------------------------------------------------
    // Add (cxx L1844-1955)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::Add(Vevo, Prof, Glue) (L1844-1955).
    pub(super) fn add(
        &mut self,
        vevo: &mut BRepFillEvolved,
        prof: &Shape,
        glue: &mut BRepToolsQuilt,
    ) {
        // OCCT L1859-1862.
        if vevo.shape().is_null() {
            return;
        }

        //-------------------------------------------------
        // Find wires common to <me> and <Vevo>.
        //-------------------------------------------------
        // OCCT L1868-1911.
        for ex_prof_current in explorer(prof, ShapeType::Vertex, ShapeType::Shape) {
            let vv = ex_prof_current;
            //---------------------------------------------------------------
            // Parse edges generated by VV in myMap if they existent
            // and Bind in Glue
            //---------------------------------------------------------------

            //------------------------------------------------- -------------
            // Note: the curves of of reinforced edges are in the same direction
            //          if one remains on the same edge.
            //          if one passes from left to the right they are inverted.
            //------------------------------------------------- -------------
            let mut commun = false;
            crate::brep_fill::brep_fill_evolved_d::relative(
                &self.my_profile.clone(),
                prof,
                &vv,
                &mut commun,
            );

            if commun {
                let spine_keys: Vec<Shape> = self.my_map.iter().map(|(ks, _)| ks.clone()).collect();
                for sp in spine_keys {
                    if self.my_map.find(&sp).is_bound(&vv)
                        && vevo.generated().is_bound(&sp)
                        && vevo.generated().find(&sp).is_bound(&vv)
                    {
                        let my_list = self.my_map.find(&sp).find(&vv).clone();
                        let vevo_list = vevo.generated_shapes(&sp, &vv);
                        for (me_s, ve_s) in my_list.iter().zip(vevo_list.iter()) {
                            let me = me_s.clone();
                            let ve = ve_s.clone();
                            let og = crate::brep_fill::brep_fill_evolved_d::compare(&me, &ve);
                            let a_local_shape = shape_oriented(&ve, Orientation::Forward);
                            let a_local_shape2 = shape_oriented(&me, og);
                            glue.bind_edge(&a_local_shape, &a_local_shape2);
                        }
                    }
                }
            }
        }
        glue.add(vevo.shape());

        //----------------------------------------------------------
        // Add map of elements generate in Vevo in myMap.
        //----------------------------------------------------------
        // OCCT L1917-1954.
        let vevo_entries: Vec<(Shape, DataMapOfShapeListOfShape)> = vevo
            .generated()
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (current_spine, inner) in vevo_entries {
            let inner_entries: Vec<(Shape, Vec<Shape>)> =
                inner.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
            for (current_prof, gen_shapes) in inner_entries {
                if !self.my_map.is_bound(&current_spine) {
                    //------------------------------------------------
                    // The element of spine is not yet present .
                    // => previous profile not on the border.
                    //-------------------------------------------------
                    self.my_map
                        .bind(&current_spine, DataMapOfShapeListOfShape::new());
                }
                if !self.my_map.find(&current_spine).is_bound(&current_prof) {
                    self.my_map
                        .find_mut(&current_spine)
                        .bind(&current_prof, Vec::new());
                    for itl in &gen_shapes {
                        // during Glue.Add the shared shapes are recreated.
                        if glue.is_copied(itl) {
                            let copy = glue.copy(itl);
                            self.my_map
                                .find_mut(&current_spine)
                                .find_mut(&current_prof)
                                .push(copy);
                        } else {
                            self.my_map
                                .find_mut(&current_spine)
                                .find_mut(&current_prof)
                                .push(itl.clone());
                        }
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------
    // Transfert (cxx L1966-2053)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::Transfert(Vevo, MapProf, MapSpine, LS,
    /// InitLS, InitLP) (L1966-2053).
    pub(super) fn transfert(
        &mut self,
        vevo: &mut BRepFillEvolved,
        map_prof: &DataMapOfShapeItem<Shape>,
        map_spine: &DataMapOfShapeItem<Shape>,
        ls: u32,
        init_ls: u32,
        init_lp: u32,
    ) {
        //----------------------------------------------------------------
        // Transfer the shape from Vevo in myShape and Reposition shapes.
        //----------------------------------------------------------------
        // OCCT L1977-1980.
        self.my_shape = vevo.shape().clone();
        self.my_spine = location_shape_set(&mut self.my_locations, &self.my_spine.clone(), init_ls);
        self.my_profile =
            location_shape_set(&mut self.my_locations, &self.my_profile.clone(), init_lp);
        self.my_shape = location_shape_moved(&mut self.my_locations, &self.my_shape.clone(), ls);

        //
        // Expecting for better, the Same Parameter is forced here
        //  ( Pb Sameparameter between YaPlanar and Tuyaux
        //
        // OCCT L1986-1994.
        for ex_current in explorer(&self.my_shape.clone(), ShapeType::Edge, ShapeType::Shape) {
            crate::brep_fill::brep_fill_evolved_d::builder_same_range(&ex_current, false);
            crate::brep_fill::brep_fill_evolved_d::builder_same_parameter(&ex_current, false);
            brep_lib_same_parameter(&ex_current);
        }

        //--------------------------------------------------------------
        // Transfer of myMap of Vevo into myMap.
        //--------------------------------------------------------------
        // OCCT L1999-2045.
        let vevo_entries: Vec<(Shape, DataMapOfShapeListOfShape)> = vevo
            .generated()
            .iter()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        for (ite_s_key, _inner_snapshot) in vevo_entries {
            let spine_shape = map_spine.find(&ite_s_key).clone();
            let mut initial_spine =
                location_shape_moved(&mut self.my_locations, &spine_shape, ls);

            // The mutable walk of OCCT ChangeFind(iteS.Key()) — the
            // generated lists are moved and re-appended.
            let inner_keys: Vec<Shape> = vevo
                .generated()
                .find(&ite_s_key)
                .iter()
                .map(|(k, _)| k.clone())
                .collect();
            for ite_p_key in inner_keys {
                let prof_shape = map_prof.find(&ite_p_key).clone();
                let mut initial_prof =
                    location_shape_set(&mut self.my_locations, &prof_shape, init_lp);
                let _ = &mut initial_prof;

                // OCCT L2025-2032: GenShapes =
                // MapVevo.ChangeFind(iteS.Key()).ChangeFind(iteP.Key()).
                let gen_shapes: Vec<Shape> = {
                    let list = vevo.generated().find_mut(&ite_s_key).find_mut(&ite_p_key);
                    for itl in list.iter_mut() {
                        *itl = location_shape_moved(&mut self.my_locations, itl, ls);
                    }
                    list.clone()
                };

                if !self.my_map.is_bound(&initial_spine) {
                    self.my_map
                        .bind(&initial_spine, DataMapOfShapeListOfShape::new());
                }
                if !self.my_map.find(&initial_spine).is_bound(&initial_prof) {
                    self.my_map
                        .find_mut(&initial_spine)
                        .bind(&initial_prof, Vec::new());
                }
                self.my_map
                    .find_mut(&initial_spine)
                    .find_mut(&initial_prof)
                    .extend(gen_shapes);
            }
        }
        //--------------------------------------------------------------
        // Transfer of Top and Bottom of Vevo in myTop and myBottom.
        //--------------------------------------------------------------
        // OCCT L2049-2052.
        let top = vevo.top().clone();
        self.my_top = location_shape_moved(&mut self.my_locations, &top, ls);
        let bottom = vevo.bottom().clone();
        self.my_bottom = location_shape_moved(&mut self.my_locations, &bottom, ls);
    }

    // -------------------------------------------------------------------
    // AddTopAndBottom (cxx L2078-2222)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::AddTopAndBottom(Glue) (L2078-2222).
    pub(super) fn add_top_and_bottom(&mut self, glue: &mut BRepToolsQuilt) {
        //  return first and last vertex of the profile.
        // OCCT L2081-2086.
        let (v0, v1) = top_exp_vertices(&self.my_profile.clone());
        let v = [v0, v1];
        if v[0].is_same(&v[1]) {
            return;
        }

        let mut to_reverse = false;
        for i in 0..=1usize {
            // OCCT L2093-2097: Loop; S = gp_Pln(0,0,1,-Altitud(V[i]));
            // F = BRepLib_MakeFace(S); Loop.Init(F).
            let mut loop_ = BRepAlgoLoop::new();
            let s = Plane::new(
                DVec3::new(0.0, 0.0, -crate::brep_fill::brep_fill_evolved_d::altitud(&v[i])),
                DVec3::Z,
            );
            let f = builder_make_face_surface(&Surface3::Plane(s), CONFUSION);
            loop_.init(&f);

            // OCCT L2099-2100: ExpSpine; View (NCollection_Map).
            let mut view: std::collections::HashSet<ShapeKey> = std::collections::HashSet::new();

            // OCCT L2102-2137.
            for exp_spine_current in
                explorer(&self.my_spine.clone(), ShapeType::Edge, ShapeType::Shape)
            {
                let es = exp_spine_current;
                let l = self.generated_shapes(&es, &v[i]);
                let mut compute_orientation = false;

                for itl_value in &l {
                    let e = itl_value.clone();

                    if !compute_orientation {
                        // OCCT L2114-2126: the BRepAdaptor_Curve D1
                        // orientation probe.
                        let (c1, fs, ls) =
                            crate::brep_algo::tool::brep_tool_curve(&es).expect("ES curve");
                        let (c2, f2, l2) =
                            crate::brep_algo::tool::brep_tool_curve(&e).expect("E curve");
                        let u = 0.3 * f2 + 0.7 * l2;
                        let us = 0.3 * fs + 0.7 * ls;
                        let v1 = CurveEval::derivative_at(&c1, us);
                        let v2 = CurveEval::derivative_at(&c2, u);
                        to_reverse = v1.dot(v2) < 0.0;
                        compute_orientation = true;
                    }

                    let mut or = es.orientation;
                    if to_reverse {
                        or = top_abs_reverse(or);
                    }
                    let a_local_shape = shape_oriented(&e, or);
                    loop_.add_const_edge(&a_local_shape);
                }
            }

            // OCCT L2140-2141.
            let pv = brep_tool_pnt(&v[i]);
            let is_out = pv.y < 0.0;

            // OCCT L2143-2185.
            for exp_spine_current in
                explorer(&self.my_spine.clone(), ShapeType::Vertex, ShapeType::Shape)
            {
                let es = exp_spine_current;
                // OCCT L2146: if (View.Add(ES)).
                if view.insert(shape_key(&es)) {
                    let l = self.generated_shapes(&es, &v[i]);
                    for itl_value in &l {
                        let e = itl_value.clone();
                        if !brep_tool_degenerated(&e) {
                            // the center of circle (ie vertex) is IN the cap if vertex IsOut
                            //                                    OUT                   !IsOut
                            // OCCT L2156-2173.
                            let (c, f, l_par) =
                                crate::brep_algo::tool::brep_tool_curve(&e).expect("E curve");
                            let u = 0.3 * f + 0.7 * l_par;
                            let p = brep_tool_pnt(&es);
                            let pc = CurveEval::point_at(&c, u);
                            let vc = CurveEval::derivative_at(&c, u);
                            let a_ppc = pc - p;
                            let prod = a_ppc.cross(vc);
                            if is_out {
                                to_reverse = prod.z < 0.0;
                            } else {
                                to_reverse = prod.z > 0.0;
                            }
                            let or = if to_reverse {
                                Orientation::Reversed
                            } else {
                                Orientation::Forward
                            };
                            let a_local_shape = shape_oriented(&e, or);
                            loop_.add_const_edge(&a_local_shape);
                        }
                    }
                }
            }

            // OCCT L2187-2212.
            loop_.perform();
            loop_.wires_to_faces();
            let l: Vec<Shape> = loop_.new_faces().to_vec();

            // Maj of myTop and myBottom for the history
            // and addition of constructed faces.
            let mut bouchon = builder_make_compound();
            let mut j = 0;

            for an_iter_l_value in &l {
                j += 1;
                glue.add(an_iter_l_value);
                if j == 1 && i == 0 {
                    self.my_top = an_iter_l_value.clone();
                }
                if j == 1 && i == 1 {
                    self.my_bottom = an_iter_l_value.clone();
                }
                builder_add_compound_shape(&mut bouchon, an_iter_l_value);
            }
            if i == 0 && j > 1 {
                self.my_top = bouchon.clone();
            }
            if i == 1 && j > 1 {
                self.my_bottom = bouchon;
            }
        }
    }

    // -------------------------------------------------------------------
    // MakeSolid (cxx L2229-2262)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::MakeSolid() (L2229-2262).
    pub(super) fn make_solid(&mut self) {
        let exp = explorer(&self.my_shape.clone(), ShapeType::Shell, ShapeType::Shape);
        let mut ish = 0;
        let mut res = builder_make_compound();
        let mut sol = builder_make_solid();

        for sh in &exp {
            // OCCT L2242-2251.
            let mut sol_cur = builder_make_solid();
            builder_add_solid_shell(&mut sol_cur, sh);
            let mut sc = SolidClassifier::new();
            sc.load(&sol_cur.clone());
            sc.perform_infinite_point(brep_fill_confusion());
            if sc.state() == TOPABS_IN {
                sol_cur = builder_make_solid();
                builder_add_solid_shell(&mut sol_cur, &shape_reversed(sh));
            }
            builder_add_compound_shape(&mut res, &sol_cur);
            sol = sol_cur;
            ish += 1;
        }
        if ish == 1 {
            self.my_shape = sol;
        } else {
            self.my_shape = res;
        }
    }

    // -------------------------------------------------------------------
    // MakePipe (cxx L2266-2325)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::MakePipe(SE, AxeRef) (L2266-2325).
    pub(super) fn make_pipe(&mut self, se: &Shape, axe_ref: &rcad_kernel::math::gp::Ax3) {
        // OCCT L2271-2278.
        let mut trsf = DAffine3::IDENTITY;
        if crate::brep_fill::brep_fill_evolved_d::side(&self.my_profile.clone(), brep_fill_confusion()) > 3 {
            // side right
            trsf = crate::brep_fill::brep_fill_evolved_d::trsf_set_rotation(&DVec3::Z, std::f64::consts::PI);
        }
        let dum_loc = self.my_locations.intern(trsf);
        let a_local_shape =
            location_shape_moved(&mut self.my_locations, &self.my_profile.clone(), dum_loc);
        let dummy_prof = crate::brep_fill::brep_fill_evolved_d::put_profil_at(
            &a_local_shape,
            axe_ref,
            se,
            &self.my_spine.clone(),
            true,
            &mut self.my_locations,
        );

        // Copy of the profile to avoid the accumulation of
        // locations on the Edges of myProfile!
        // OCCT L2287-2288: TrsfMod = new
        // BRepTools_TrsfModification(gp_Trsf()); Modif =
        // BRepTools_Modifier(DummyProf, TrsfMod) — the OCCT two-argument
        // constructor (Init + Perform of BRepTools_Modifier.cxx L71-79).
        let mut trsf_mod = BRepToolsTrsfModification::new(rcad_kernel::math::gp::Trsf::identity());
        let modif = BRepToolsModifier::with_shape_and_modification(&dummy_prof, &mut trsf_mod);
        // OCCT L2290: GenProf = TopoDS::Wire(Modif.ModifiedShape(DummyProf)).
        let gen_prof = modif.modified_shape(&dummy_prof);

        // OCCT L2292: Pipe = BRepFill_Pipe(BRepLib_MakeWire(SE), GenProf).
        let spine_wire = crate::brep_fill::brep_fill_evolved::brep_lib_make_wire_from_edge(se);
        let mut pipe =
            BRepFillPipe::new(&spine_wire, &gen_prof, GeomFillTrihedron::IsFrenet, false, false);

        //---------------------------------------------
        // Arrangement of Tubes in myMap.
        //---------------------------------------------
        // OCCT L2299-2324.
        let mut l: Vec<Shape> = Vec::new();
        let mut first_vertex = true;

        // OCCT L2303-2305: DataMap P; myMap.Bind(SE, P) — an empty map.
        self.my_map
            .bind(se, DataMapOfShapeListOfShape::new());

        let prof_edges = wire_edges(&self.my_profile.clone());
        let gen_prof_edges = wire_edges(&gen_prof);
        let nb = prof_edges.len().min(gen_prof_edges.len());
        for idx in 0..nb {
            let prof_current = prof_edges[idx].clone();
            let gen_prof_current = gen_prof_edges[idx].clone();

            let (vf, vl) = edge_vertices(&prof_current);
            let (vfg, vlg) = edge_vertices(&gen_prof_current);

            if first_vertex {
                self.my_map.find_mut(se).bind(&vf, l.clone());
                let edge = pipe.edge(se, &vfg);
                self.my_map.find_mut(se).find_mut(&vf).push(edge);
                first_vertex = false;
            }
            self.my_map.find_mut(se).bind(&vl, l.clone());
            let edge = pipe.edge(se, &vlg);
            self.my_map.find_mut(se).find_mut(&vl).push(edge);
            self.my_map.find_mut(se).bind(&prof_current, l.clone());
            let face = pipe.face(se, &gen_prof_current);
            self.my_map.find_mut(se).find_mut(&prof_current).push(face);
        }
    }

    // -------------------------------------------------------------------
    // MakeRevol (cxx L2329-2400)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::MakeRevol(SE, VLast, AxeRef) (L2329-2400).
    pub(super) fn make_revol(
        &mut self,
        se: &Shape,
        v_last: &Shape,
        axe_ref: &rcad_kernel::math::gp::Ax3,
    ) {
        // OCCT L2336-2343.
        let mut trsf = DAffine3::IDENTITY;
        if crate::brep_fill::brep_fill_evolved_d::side(&self.my_profile.clone(), brep_fill_confusion()) > 3 {
            // side right
            trsf = crate::brep_fill::brep_fill_evolved_d::trsf_set_rotation(&DVec3::Z, std::f64::consts::PI);
        }
        let dum_loc = self.my_locations.intern(trsf);
        let a_local_shape =
            location_shape_moved(&mut self.my_locations, &self.my_profile.clone(), dum_loc);
        let mut gen_prof = crate::brep_fill::brep_fill_evolved_d::put_profil_at(
            &a_local_shape,
            axe_ref,
            se,
            &self.my_spine.clone(),
            false,
            &mut self.my_locations,
        );

        // OCCT L2349: AxeRev = gp_Ax1(BRep_Tool::Pnt(VLast), -gp::DZ()).
        let axe_rev = rcad_kernel::math::gp::Ax1::new(brep_tool_pnt(v_last), -DVec3::Z);

        // Position of the sewing on the edge of the spine
        // so that the bissectrices didn't cross the sewings.
        // OCCT L2353-2356.
        let dummy =
            crate::brep_fill::brep_fill_evolved_d::trsf_set_rotation(&axe_rev.direction, 1.5 * std::f64::consts::PI);
        let dummy_loc = self.my_locations.intern(dummy);
        gen_prof = location_shape_moved(&mut self.my_locations, &gen_prof, dummy_loc);

        // OCCT L2358.
        let rev = BRepSweepRevolCarrier::new(&gen_prof.clone(), axe_rev, true);

        //--------------------------------------------
        // Arrangement of revolutions in myMap.
        //---------------------------------------------
        // OCCT L2363-2399.
        let mut l: Vec<Shape> = Vec::new();
        let mut first_vertex = true;

        self.my_map.bind(v_last, DataMapOfShapeListOfShape::new());

        let prof_edges = wire_edges(&self.my_profile.clone());
        let gen_prof_edges = wire_edges(&gen_prof);
        let nb = prof_edges.len().min(gen_prof_edges.len());
        for idx in 0..nb {
            let prof_current = prof_edges[idx].clone();
            let gen_prof_current = gen_prof_edges[idx].clone();

            let (vf, vl) = edge_vertices(&prof_current);
            let (vfg, vlg) = edge_vertices(&gen_prof_current);

            let or = gen_prof_current.orientation;

            if first_vertex {
                self.my_map.find_mut(v_last).bind(&vf, l.clone());
                let rv = rev.shape_of(&vfg);
                //      TopAbs_Orientation OO = TopAbs::Compose(RV.Orientation(),Or);
                let oo = rv.orientation;
                self.my_map
                    .find_mut(v_last)
                    .find_mut(&vf)
                    .push(shape_oriented(&rv, oo));
                first_vertex = false;
            }
            self.my_map.find_mut(v_last).bind(&prof_current, l.clone());
            let rf = rev.shape_of(&gen_prof_current);
            let oo = top_abs_compose(rf.orientation, or);
            self.my_map
                .find_mut(v_last)
                .find_mut(&prof_current)
                .push(shape_oriented(&rf, oo));
            self.my_map.find_mut(v_last).bind(&vl, l.clone());
            let rv = rev.shape_of(&vlg);
            //    OO = TopAbs::Compose(RV.Orientation(),Or);
            let oo = rv.orientation;
            self.my_map
                .find_mut(v_last)
                .find_mut(&vl)
                .push(shape_oriented(&rv, oo));
        }
    }

    // -------------------------------------------------------------------
    // FindLocation (cxx L2404-2437)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::FindLocation(Face) const (L2404-2437) — the
    /// location transforming the planar shape in the plane xOy (if the
    /// shape is not planar).
    pub(super) fn find_location(&mut self, face: &Shape) -> u32 {
        // OCCT L2406-2408.
        let mut s = brep_tool_surface(face).expect("BRep_Tool::Surface");

        let is_plane = matches!(s, Surface3::Plane(_));
        if !is_plane {
            // OCCT L2410-2422: BRepLib_FindSurface FS(Face, -1, true);
            // if (FS.Found()) { S = FS.Surface(); L = FS.Location(); }.
            // The FindSurface pool lookup resolves the pcurve owner face
            // through the brep it is given (topalgo/brep_lib_find_surface
            // bridge #1); BRepFill_Evolved carries no pool of its own, so the
            // face is passed with an empty one: the existing-surface probe
            // finds nothing and the OCCT least-squares fallback (cxx
            // L422-579) fits the plane from the shape's own curves.
            let fs = BRepLibFindSurface::new(
                &mut rcad_kernel::topods::BRep::new(),
                face,
                -1.0,
                true,
            );
            if fs.found() {
                s = fs.surface().expect("BRepLib_FindSurface::Surface");
                // OCCT L2416: L = FS.Location() — the found surface carries
                // its location in the rcad encoding (architecture difference).
            } else {
                // OCCT L2420: throw Standard_NoSuchObject.
                panic!("Standard_NoSuchObject: BRepFill_Evolved : The Face is not planar");
            }
        }

        // OCCT L2424-2427: if (!L.IsIdentity()) S = S->Transformed(...) —
        // the face location is carried by the surface encoding (identity in
        // the pool-free world; the branch keeps the OCCT shape).

        // OCCT L2429-2430: P = down_cast<Geom_Plane>(S); Axis = P->Position().
        let axis = match &s {
            Surface3::Plane(p) => {
                rcad_kernel::math::gp::Ax3::from_pnt_n_vx(p.origin, p.normal, p.u_dir)
            }
            _ => {
                // OCCT: the down_cast is guaranteed by the Found branch.
                panic!("BRepFill_Evolved::FindLocation: not a plane");
            }
        };

        // OCCT L2432-2434.
        let axe_ref = rcad_kernel::math::gp::Ax3::from_pnt_n_vx(DVec3::ZERO, DVec3::Z, DVec3::X);
        let t = crate::brep_fill::brep_fill_evolved_d::trsf_set_transformation_between(&axe_ref, &axis);

        // OCCT L2436: return TopLoc_Location(T).
        self.my_locations.intern(t)
    }

    // -------------------------------------------------------------------
    // TransformInitWork (cxx L2441-2445)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::TransformInitWork(LS, LP) (L2441-2445) —
    /// applies LS to mySpine and LP to myProfil in order to set the Shapes
    /// in the work space.
    pub(super) fn transform_init_work(&mut self, ls: u32, lp: u32) {
        self.my_spine = location_shape_moved(&mut self.my_locations, &self.my_spine.clone(), ls);
        self.my_profile =
            location_shape_moved(&mut self.my_locations, &self.my_profile.clone(), lp);
    }

    // -------------------------------------------------------------------
    // ContinuityOnOffsetEdge (cxx L2452-2544)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::ContinuityOnOffsetEdge(WorkProf) (L2452-2544)
    /// — coding of regularities on edges parallel to CutVevo common to
    /// left and right parts of volevo.
    pub(super) fn continuity_on_offset_edge(&mut self, _work_prof: &Vec<Shape>) {
        // OCCT L2463-2470: WExp.Init(myProfile); FirstE = WExp.Current() —
        // the Current() read precedes the More() check (OCCT source as
        // written); CurrentVertex is the vertex connecting the current edge
        // to the previous one (the walk-first vertex — the
        // edge_vertices.0 reduction).
        let profile_edges = wire_edges(&self.my_profile.clone());
        let first_e = profile_edges
            .first()
            .cloned()
            .expect("BRepTools_WireExplorer::Current (the profile wire is empty)");
        let mut prec_e = first_e.clone();
        let (vf, mut v) = edge_vertices(&first_e);
        let mut cur_e = Shape::null();
        let mut vl = Shape::null();

        let rest_edges: Vec<Shape> = if profile_edges.len() > 1 {
            profile_edges[1..].to_vec()
        } else {
            Vec::new()
        };

        for w_exp_current in &rest_edges {
            cur_e = w_exp_current.clone();
            // OCCT L2475: V = WExp.CurrentVertex().
            v = edge_vertices(&cur_e).0;

            if crate::brep_fill::brep_fill_evolved_d::distance_to_oz(&v) <= brep_fill_confusion() {
                // the regularities are already coded on the edges of
                // elementary volevos
                // OCCT L2480-2484.
                let u1 = crate::brep_algo::tool::brep_tool_parameter(&v, &cur_e);
                let u2 = crate::brep_algo::tool::brep_tool_parameter(&v, &prec_e);
                let continuity =
                    crate::brep_fill::brep_fill_evolved_d::brep_lprop_continuity(&cur_e, &prec_e, u1, u2);

                if continuity >= 1 {
                    //-----------------------------------------------------
                    // Code continuity for all edges generated by V.
                    //-----------------------------------------------------
                    // OCCT L2491-2505.
                    let spine_keys: Vec<Shape> =
                        self.my_map.iter().map(|(ks, _)| ks.clone()).collect();
                    for sp in spine_keys {
                        if self.my_map.find(&sp).is_bound(&v)
                            && self.my_map.find(&sp).is_bound(&cur_e)
                            && self.my_map.find(&sp).is_bound(&prec_e)
                        {
                            if !self.my_map.find(&sp).find(&v).is_empty()
                                && !self.my_map.find(&sp).find(&cur_e).is_empty()
                                && !self.my_map.find(&sp).find(&prec_e).is_empty()
                            {
                                let e1 = self
                                    .my_map
                                    .find(&sp)
                                    .find(&v)
                                    .first()
                                    .cloned()
                                    .expect("myMap(SP)(V)");
                                let f1 = self
                                    .my_map
                                    .find(&sp)
                                    .find(&cur_e)
                                    .first()
                                    .cloned()
                                    .expect("myMap(SP)(CurE)");
                                let f2 = self
                                    .my_map
                                    .find(&sp)
                                    .find(&prec_e)
                                    .first()
                                    .cloned()
                                    .expect("myMap(SP)(PrecE)");
                                crate::brep_fill::brep_fill_evolved::builder_continuity(
                                    &e1, &f1, &f2, continuity,
                                );
                            }
                        }
                    }
                }
            }
            prec_e = cur_e.clone();
        }

        // OCCT L2511: EdgeVertices(PrecE, V, VL).
        let (v_out, vl_out) = edge_vertices(&prec_e);
        v = v_out;
        vl = vl_out;

        if vf.is_same(&vl) {
            // Closed profile.
            // OCCT L2516-2520.
            let u1 = crate::brep_algo::tool::brep_tool_parameter(&vf, &cur_e);
            let u2 = crate::brep_algo::tool::brep_tool_parameter(&vf, &first_e);
            let continuity =
                crate::brep_fill::brep_fill_evolved_d::brep_lprop_continuity(&cur_e, &first_e, u1, u2);

            if continuity >= 1 {
                //---------------------------------------------
                // Code continuity for all edges generated by V.
                //---------------------------------------------
                // OCCT L2527-2541.
                let spine_keys: Vec<Shape> =
                    self.my_map.iter().map(|(ks, _)| ks.clone()).collect();
                for sp in spine_keys {
                    if self.my_map.find(&sp).is_bound(&vf)
                        && self.my_map.find(&sp).is_bound(&cur_e)
                        && self.my_map.find(&sp).is_bound(&first_e)
                    {
                        if !self.my_map.find(&sp).find(&vf).is_empty()
                            && !self.my_map.find(&sp).find(&cur_e).is_empty()
                            && !self.my_map.find(&sp).find(&first_e).is_empty()
                        {
                            let e1 = self
                                .my_map
                                .find(&sp)
                                .find(&vf)
                                .first()
                                .cloned()
                                .expect("myMap(SP)(VF)");
                            let f1 = self
                                .my_map
                                .find(&sp)
                                .find(&cur_e)
                                .first()
                                .cloned()
                                .expect("myMap(SP)(CurE)");
                            let f2 = self
                                .my_map
                                .find(&sp)
                                .find(&first_e)
                                .first()
                                .cloned()
                                .expect("myMap(SP)(FirstE)");
                            crate::brep_fill::brep_fill_evolved::builder_continuity(
                                &e1, &f1, &f2, continuity,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// OCCT TopAbs::Compose(O1, O2) — the orientation composition
/// (Forward x Reversed = Reversed, Reversed x Reversed = Forward).
pub(super) fn top_abs_compose(o1: Orientation, o2: Orientation) -> Orientation {
    match (o1, o2) {
        (Orientation::Forward, _) => o2,
        (_, Orientation::Forward) => Orientation::Reversed,
        (Orientation::Reversed, Orientation::Reversed) => Orientation::Forward,
        // The Internal / External compositions never occur on the shapes
        // consumed here (TopAbs::Compose is total in OCCT; the rcad enum
        // carries the two extra TopAbs_State-like variants).
        _ => Orientation::Forward,
    }
}

/// OCCT TopoDS_Iterator(Face) — the face wires (the TFaceData reduction).
pub(super) fn face_wires(f: &Shape) -> Vec<Shape> {
    match f.data.as_ref() {
        rcad_kernel::topods::TShape::Face(fd) => {
            let mut wires = Vec::new();
            if !fd.outer_wire.is_null() {
                wires.push(fd.outer_wire.clone());
            }
            wires.extend(fd.inner_wires.iter().cloned());
            wires
        }
        _ => Vec::new(),
    }
}
