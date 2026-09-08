// OCCT BRepOffset_MakeOffset.cxx — 1:1 translation, module c (see the
// module header of brep_offset_make_offset.rs for the split map and the
// architecture-difference list #38-#56).
//
// Module c carries cxx L2396-L3108: CorrectConicalFaces / Intersection3D /
// Intersection2D / MakeLoops / MakeFaces.

use std::collections::HashMap;

use rcad_kernel::geom::Curve2dEval;
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_inter2d::DmvvMap;
use super::brep_offset_inter2d_b::BRepOffsetInter2d;
use super::brep_offset_inter3d::BRepOffsetInter3d;
use super::brep_offset_make_offset::{BRepOffset_Error, BRepOffsetMakeOffset};
use super::brep_offset_tool::OcctIndexedShapeMap;
use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::tool as bat;
use crate::brep_fill::offset_wire::GeomAbsJoinType;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range,
};

use super::brep_offset_make_offset::brep_lib_same_parameter_3;

use super::brep_offset_make_offset::{top_exp_first_vertex, top_exp_map_shapes};

// ===========================================================================
// Local re-hosts of the module (architecture differences #54, #57).
// ===========================================================================

/// OCCT BRepAdaptor_Curve(E, F) re-host — the edge 3D curve with its range
/// (the inter2d.rs #54 precedent; the GetType / Circle accessors are the
/// consumed surface of this module).
pub(crate) struct BRepAdaptorCurveC {
    my_curve: rcad_kernel::geom::Curve3,
    my_first: f64,
    my_last: f64,
}

impl BRepAdaptorCurveC {
    /// OCCT BRepAdaptor_Curve(E).
    pub fn new(the_e: &Shape) -> Self {
        match bat::brep_tool_curve(the_e) {
            Some((c, f, l)) => BRepAdaptorCurveC {
                my_curve: c,
                my_first: f,
                my_last: l,
            },
            None => panic!("GAP: BRepAdaptor_Curve(E) — the curve-on-surface 3D fallback"),
        }
    }

    /// OCCT BRepAdaptor_Curve::FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT BRepAdaptor_Curve::LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT BRepAdaptor_Curve::Value(U).
    pub fn value(&self, the_u: f64) -> glam::DVec3 {
        use rcad_kernel::geom::CurveEval;
        self.my_curve.point_at(the_u)
    }

    /// OCCT BRepAdaptor_Curve::GetType() — the GeomAbs_Liner /
    /// GeomAbs_Circle probes map onto the Curve3 variant match (arch.
    /// diff. #18/#54).
    pub fn get_type(&self) -> &rcad_kernel::geom::Curve3 {
        &self.my_curve
    }

    /// OCCT BRepAdaptor_Curve::Circle().
    pub fn circle(&self) -> rcad_kernel::geom::Circle3 {
        match &self.my_curve {
            rcad_kernel::geom::Curve3::Circle(c) => *c,
            _ => panic!("BRepAdaptor_Curve::Circle — the curve is not a circle"),
        }
    }
}

impl BRepOffsetMakeOffset {
    /// OCCT BRepOffset_MakeOffset::CorrectConicalFaces (cxx L2396-2853).
    pub(crate) fn correct_conical_faces(&mut self) {
        if self.my_offset_shape.is_null() {
            return;
        }
        //
        // OCCT L2410-2412: the Cones / Circs / Seams sequences — the rcad
        // Vec forms (appended only by the OCCT commented-out code).
        let _cones: Vec<Shape> = Vec::new();
        let _circs: Vec<Shape> = Vec::new();
        let _seams: Vec<Shape> = Vec::new();
        let tol_apex = 1.0e-5;

        let mut faces_of_cone:
            super::brep_offset_make_offset::DataMapOfShapeListOfShape = HashMap::new();
        if self.my_join == GeomAbsJoinType::Arc {
            for explo in
                bat::explorer(&self.my_offset_shape, ShapeType::Face, ShapeType::Shape)
            {
                let a_face = explo;
                // OCCT L2426: aSurf = BRep_Tool::Surface(aFace) — the face
                // surface carrier (the value is consumed by the degenerate
                // check only through the edge walk).
                let _a_surf = bat::brep_tool_surface(&a_face);

                let emap = top_exp_map_shapes(&a_face, ShapeType::Edge);
                for i in 1..=emap.extent() {
                    let an_edge = emap.at_1(i).clone();
                    if brep_tool_degenerated(&an_edge) {
                        // Check if anEdge is a really degenerated edge or not
                        // OCCT L2452-2458: BRepAdaptor_Curve BACurve(anEdge,
                        // aFace) — the 3D curve value path.
                        let ba_curve = BRepAdaptorCurveC::new(&an_edge);
                        let pfirst = ba_curve.value(ba_curve.first_parameter());
                        let plast = ba_curve.value(ba_curve.last_parameter());
                        let pmid = ba_curve
                            .value((ba_curve.first_parameter() + ba_curve.last_parameter()) / 2.);
                        if pfirst.distance(plast) <= tol_apex && pfirst.distance(pmid) <= tol_apex
                        {
                            continue;
                        }
                        let or_edge = self.my_init_offset_edge.root(&an_edge).clone();
                        let vf = top_exp_first_vertex(&or_edge);
                        if super::brep_offset_tool::shape_data_map::is_bound(
                            &faces_of_cone,
                            &vf,
                        ) {
                            // add a face to the existing list
                            let a_faces =
                                super::brep_offset_tool::shape_data_map::change_find(
                                    &mut faces_of_cone,
                                    &vf,
                                );
                            a_faces.push(a_face.clone());
                        } else {
                            // the vertex is not in the map => create a new key and items
                            let a_faces: Vec<Shape> = vec![a_face.clone()];
                            super::brep_offset_tool::shape_data_map::bind(
                                &mut faces_of_cone,
                                &vf,
                                a_faces,
                            );
                        }
                    }
                } // for (i = 1; i <= Emap.Extent(); i++)
            } // for (; fexp.More(); fexp.Next())
        } // if (myJoin == GeomAbs_Arc)

        // OCCT L2495: the FacesOfCone iterator — the rcad HashMap order
        // (architecture difference #55).
        let cone_items: Vec<(Shape, Vec<Shape>)> = faces_of_cone
            .values()
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();
        let mut bb = rcad_kernel::topo::topods::BRepBuilder::new();
        let mut is_modified = false;
        for (an_apex, faces) in &cone_items {
            // OCCT L2501-2503: gp_Sphere theSphere; aSphSurf; SphereWire.
            let mut the_sphere_center = glam::DVec3::ZERO;
            let mut the_sphere_axis = glam::DVec3::X;
            let mut the_sphere_ref_dir = glam::DVec3::Y;
            let mut the_sphere_has_position = false;
            let mut the_sphere_radius = 0f64;
            let mut sphere_wire = bb.make_wire(&mut self.my_brep);
            let an_apex = an_apex.clone();
            let faces = faces; // FacesOfCone(anApex);
            let mut is_first_face = true;
            let mut first_point = glam::DVec3::ZERO;
            let mut the_first_vertex = Shape::null();
            let mut cur_first_vertex = Shape::null();
            for it_faces in faces {
                let a_face = it_faces.clone(); // TopoDS::Face(Faces.First());
                let mut deg_edge = Shape::null(); // = TopoDS::Edge(DegEdges(aFace));
                for explo in bat::explorer(&a_face, ShapeType::Edge, ShapeType::Shape) {
                    deg_edge = explo;
                    if brep_tool_degenerated(&deg_edge) {
                        let or_edge = self.my_init_offset_edge.root(&deg_edge).clone();
                        let vf = top_exp_first_vertex(&or_edge);
                        if vf.is_same(&an_apex) {
                            break;
                        }
                    }
                }
                let cur_edge = bat::oriented(&deg_edge, Orientation::Forward);
                let mut cur_edge = cur_edge;
                bat::builder_set_degenerated(&mut cur_edge, false);
                // OCCT L2547-2548: BB.SameRange / BB.SameParameter(false).
                bb.set_edge_same_range(&mut self.my_brep, cur_edge.clone(), false);
                bb.set_edge_same_parameter(&mut self.my_brep, cur_edge.clone(), false);
                let mut fpnt = glam::DVec3::ZERO;
                let mut mpnt = glam::DVec3::ZERO;
                let mut lpnt = glam::DVec3::ZERO;
                super::brep_offset_make_offset::get_edge_points(
                    &cur_edge,
                    &a_face,
                    &mut fpnt,
                    &mut mpnt,
                    &mut lpnt,
                );
                let (f, l) = brep_tool_range(&cur_edge);
                let _ = (f, l);
                if is_first_face {
                    let a_vec1 = mpnt - fpnt;
                    let a_vec2 = lpnt - fpnt;
                    let a_norm = a_vec1.cross(a_vec2);
                    let the_apex =
                        bat::brep_tool_pnt(&an_apex).unwrap_or(glam::DVec3::ZERO);
                    let apex_to_fpnt = fpnt - the_apex;
                    // OCCT L2565-2566: Ydir = aNorm ^ ApexToFpnt; Xdir = Ydir
                    // ^ aNorm (the gp_Vec cross form).
                    let y_dir = a_norm.cross(apex_to_fpnt);
                    let x_dir = y_dir.cross(a_norm);
                    // OCCT L2568-2571: gp_Ax2 anAx2(theApex, Dir(aNorm),
                    // Dir(Xdir)); theSphere.SetRadius(myOffset);
                    // theSphere.SetPosition(gp_Ax3(anAx2)) — the Ax2/Ax3
                    // fold into the SphericalSurface fields (arch. diff.
                    // #53).
                    the_sphere_center = the_apex;
                    the_sphere_axis = a_norm;
                    the_sphere_ref_dir = x_dir;
                    the_sphere_radius = self.my_offset;
                    the_sphere_has_position = true;
                    first_point = fpnt;
                    // OCCT L2575-2576: BRepLib_MakeVertex(fPnt).
                    the_first_vertex = bb.add_vertex(
                        &mut self.my_brep,
                        fpnt,
                        rcad_kernel::core::precision::CONFUSION,
                    );
                    cur_first_vertex = the_first_vertex.clone();
                }

                let (v1, v2) = super::brep_offset_tool::top_exp_vertices(&cur_edge);
                let first_vert = cur_first_vertex.clone();
                let end_vert;
                if lpnt.distance(first_point) <= rcad_kernel::core::precision::CONFUSION {
                    end_vert = the_first_vertex.clone();
                } else {
                    end_vert = bb.add_vertex(
                        &mut self.my_brep,
                        lpnt,
                        rcad_kernel::core::precision::CONFUSION,
                    );
                }
                // OCCT L2591: CurEdge.Free(true).
                bat::builder_set_free(&mut cur_edge, true);
                bat::builder_remove_edge_vertex(&mut cur_edge, &v1);
                bat::builder_remove_edge_vertex(&mut cur_edge, &v2);
                bat::builder_add_edge_vertex(
                    &mut cur_edge,
                    &bat::oriented(&first_vert, Orientation::Forward),
                );
                bat::builder_add_edge_vertex(
                    &mut cur_edge,
                    &bat::oriented(&end_vert, Orientation::Reversed),
                );
                // take the curve from sphere an put it to the edge
                // OCCT L2600-2604: ElSLib::Parameters(theSphere, fPnt, Uf,
                // Vf); ElSLib::Parameters(theSphere, lPnt, Ul, Vl); the
                // |Ul| <= Confusion clamp to 2.PI.
                let (mut uf, vf) = elslib_sphere_parameters_host(
                    fpnt,
                    the_sphere_center,
                    the_sphere_ref_dir,
                    the_sphere_axis,
                );
                let (mut ul, _vl) = elslib_sphere_parameters_host(
                    lpnt,
                    the_sphere_center,
                    the_sphere_ref_dir,
                    the_sphere_axis,
                );
                if ul.abs() <= rcad_kernel::core::precision::CONFUSION {
                    ul = 2. * std::f64::consts::PI;
                }
                let _ = vf;
                // OCCT L2605: aCurv = aSphSurf->VIso(Vf) — the TKMath
                // ElSLib::SphereVIso iso constructor (GAP leaf, the bi_tgte
                // #28 precedent).
                let a_curv = elslib_sphere_viso(
                    the_sphere_center,
                    the_sphere_axis,
                    the_sphere_ref_dir,
                    the_sphere_radius,
                    vf,
                );
                // OCCT L2616: aTrimCurv = new Geom_TrimmedCurve(aCurv, Uf,
                // Ul).
                let a_trim_curv = rcad_kernel::geom::Curve3::Trimmed(
                    rcad_kernel::geom::TrimmedCurve3 {
                        curve: Box::new(a_curv.clone()),
                        first: uf,
                        last: ul,
                    },
                );
                // OCCT L2617: BB.UpdateEdge(CurEdge, aTrimCurv,
                // Precision::Confusion()).
                super::brep_offset_offset::update_edge_3d(
                    &cur_edge,
                    &a_trim_curv,
                    0,
                    rcad_kernel::core::precision::CONFUSION,
                );
                // OCCT L2618: BB.Range(CurEdge, Uf, Ul, true).
                bat::builder_range_edge(&mut cur_edge, uf, ul);
                // OCCT L2619-2621: the Vf iso pcurve
                // Geom2d_Line((0., Vf), gp::DX2d()) trimmed to [Uf, Ul].
                let the_trim_lin2d = rcad_kernel::geom::Curve2d::Trimmed(
                    rcad_kernel::geom::TrimmedCurve2 {
                        curve: Box::new(rcad_kernel::geom::Curve2d::Line(
                            rcad_kernel::geom::Line2d {
                                origin: glam::DVec2::new(0., vf),
                                direction: glam::DVec2::new(1., 0.),
                            },
                        )),
                        t_min: uf,
                        t_max: ul,
                    },
                );
                // OCCT L2622: BB.UpdateEdge(CurEdge, theTrimLin2d, aSphSurf,
                // L, Precision::Confusion()).
                bat::builder_update_edge_pcurve(
                    &mut cur_edge,
                    &the_trim_lin2d,
                    &a_face,
                    rcad_kernel::core::precision::CONFUSION,
                );
                // OCCT L2623: BB.Range(CurEdge, aSphSurf, L, Uf, Ul).
                bat::builder_range_edge_on_face(
                    &mut cur_edge,
                    &a_face,
                    uf,
                    ul,
                );
                // OCCT L2624: BRepLib::SameParameter(CurEdge).
                brep_lib_same_parameter_3(&cur_edge, rcad_kernel::core::precision::CONFUSION);
                bat::builder_add_wire_edge(&mut sphere_wire, &cur_edge);
                // Modifying correspondent edges in aFace: substitute vertices common with CurEdge
                // OCCT L2627-2629: BRepAdaptor_Curve2d BAc2d(CurEdge, aFace).
                let (ba_c2d, ba_f, ba_l) = brep_tool_curve_on_surface(&cur_edge, &a_face)
                    .expect("BRep_Tool::CurveOnSurface null");
                let fpnt2d = ba_c2d.point_at(ba_f);
                let lpnt2d = ba_c2d.point_at(ba_l);
                let emap = top_exp_map_shapes(&a_face, ShapeType::Edge);
                let mut ee: Vec<Shape> = Vec::new();
                for k in 1..=emap.extent() {
                    let an_edge = emap.at_1(k).clone();
                    if !brep_tool_degenerated(&an_edge) {
                        let (ev1, ev2) = super::brep_offset_tool::top_exp_vertices(&an_edge);
                        if ev1.is_same(&v1) || ev2.is_same(&v1) {
                            ee.push(an_edge);
                        }
                    }
                }
                for k in 0..ee.len() {
                    let eforward = bat::oriented(&ee[k], Orientation::Forward);
                    let mut eforward = eforward;
                    // OCCT L2662: Eforward.Free(true).
                    bat::builder_set_free(&mut eforward, true);
                    let (ev1, ev2) = super::brep_offset_tool::top_exp_vertices(&eforward);
                    let (ee_c, ee_f, ee_l) =
                        brep_tool_curve_on_surface(&eforward, &a_face)
                            .expect("BRep_Tool::CurveOnSurface null");
                    let p2d1 = ee_c.point_at(ee_f);
                    let p2d2 = ee_c.point_at(ee_l);
                    if ev1.is_same(&v1) {
                        let new_v = if p2d1.distance(fpnt2d)
                            <= rcad_kernel::core::precision::CONFUSION
                        {
                            first_vert.clone()
                        } else {
                            end_vert.clone()
                        };
                        bat::builder_remove_edge_vertex(&mut eforward, &ev1);
                        bat::builder_add_edge_vertex(
                            &mut eforward,
                            &bat::oriented(&new_v, Orientation::Forward),
                        );
                    } else {
                        let new_v = if p2d2.distance(fpnt2d)
                            <= rcad_kernel::core::precision::CONFUSION
                        {
                            first_vert.clone()
                        } else {
                            end_vert.clone()
                        };
                        bat::builder_remove_edge_vertex(&mut eforward, &ev2);
                        bat::builder_add_edge_vertex(
                            &mut eforward,
                            &bat::oriented(&new_v, Orientation::Reversed),
                        );
                    }
                }

                is_first_face = false;
                cur_first_vertex = end_vert;
            }
            // Building new spherical face
            // OCCT L2683: Ufirst = RealLast(), Ulast = RealFirst().
            let mut ufirst = f64::MAX;
            let mut ulast = f64::MIN;
            let mut p2d1 = glam::DVec2::ZERO;
            let mut p2d2 = glam::DVec2::ZERO;
            let mut edges_of_wire: Vec<Shape> = Vec::new();
            for itw in bat::sub_shapes(&sphere_wire) {
                let an_edge = itw;
                edges_of_wire.push(an_edge.clone());
                // OCCT L2691: aC2d = BRep_Tool::CurveOnSurface(anEdge,
                // aSphSurf, L, f, l) — the rcad curve-on-surface re-host
                // reads the pcurve stored on the edge (the face key is the
                // identity-location form).
                let a_c2d = match bat::brep_tool_first_curve_on_surface(&an_edge) {
                    Some((c, f, l)) => (c, f, l),
                    None => panic!("BRep_Tool::CurveOnSurface null"),
                };
                p2d1 = a_c2d.0.point_at(a_c2d.1);
                p2d2 = a_c2d.0.point_at(a_c2d.2);
                if p2d1.x < ufirst {
                    ufirst = p2d1.x;
                }
                if p2d1.x > ulast {
                    ulast = p2d1.x;
                }
                if p2d2.x < ufirst {
                    ufirst = p2d2.x;
                }
                if p2d2.x > ulast {
                    ulast = p2d2.x;
                }
            }
            let mut new_edges: Vec<Shape> = Vec::new();
            let mut first_edge = Shape::null();
            let mut remove_index: Option<usize> = None;
            for (idx, itl) in edges_of_wire.iter().enumerate() {
                first_edge = itl.clone();
                let a_c2d = match bat::brep_tool_first_curve_on_surface(&first_edge) {
                    Some((c, f, l)) => (c, f, l),
                    None => panic!("BRep_Tool::CurveOnSurface null"),
                };
                p2d1 = a_c2d.0.point_at(a_c2d.1);
                p2d2 = a_c2d.0.point_at(a_c2d.2);
                if (p2d1.x - ufirst).abs() <= rcad_kernel::core::precision::CONFUSION {
                    remove_index = Some(idx);
                    break;
                }
            }
            if let Some(idx) = remove_index {
                edges_of_wire.remove(idx);
            }
            new_edges.push(bat::oriented(&first_edge, Orientation::Forward));
            let (vf1, mut cur_vertex) =
                super::brep_offset_tool::top_exp_vertices(&first_edge);
            // OCCT L2726-2742: the wire chain walk with the list-iterator
            // removal form.
            let mut idx = 0usize;
            while idx < edges_of_wire.len() {
                let an_edge = edges_of_wire[idx].clone();
                let (v1, v2) = super::brep_offset_tool::top_exp_vertices(&an_edge);
                if v1.is_same(&cur_vertex) || v2.is_same(&cur_vertex) {
                    new_edges.push(bat::oriented(&an_edge, Orientation::Forward));
                    cur_vertex = if v1.is_same(&cur_vertex) { v2 } else { v1 };
                    edges_of_wire.remove(idx);
                } else {
                    idx += 1;
                }
            }

            let (vfirst, vlast);
            if p2d1.y > 0. {
                vfirst = p2d1.y;
                vlast = std::f64::consts::PI / 2.;
            } else {
                vfirst = -std::f64::consts::PI / 2.;
                vlast = p2d1.y;
            }
            // OCCT L2753: BRepLib_MakeFace(aSphSurf, Ufirst, Ulast, Vfirst,
            // Vlast, Precision::Confusion()) — the UV-bounds face maker
            // (GAP leaf; TKTopAlgo/BRepLib not translated).
            let mut new_spherical_face = brep_lib_make_face_uv(
                the_sphere_center,
                the_sphere_axis,
                the_sphere_ref_dir,
                the_sphere_radius,
                the_sphere_has_position,
                ufirst,
                ulast,
                vfirst,
                vlast,
                rcad_kernel::core::precision::CONFUSION,
            );
            let mut old_edge = Shape::null();
            let mut deg_edge = Shape::null();
            for explo in
                bat::explorer(&new_spherical_face, ShapeType::Edge, ShapeType::Shape)
            {
                deg_edge = explo;
                if brep_tool_degenerated(&deg_edge) {
                    break;
                }
            }
            let deg_vertex = top_exp_first_vertex(&deg_edge);
            for explo in
                bat::explorer(&new_spherical_face, ShapeType::Edge, ShapeType::Shape)
            {
                old_edge = explo;
                let (ov1, ov2) = super::brep_offset_tool::top_exp_vertices(&old_edge);
                if !ov1.is_same(&deg_vertex) && !ov2.is_same(&deg_vertex) {
                    break;
                }
            }
            let (v1, v2) = super::brep_offset_tool::top_exp_vertices(&old_edge);
            let lv1: Vec<Shape> = vec![bat::oriented(&vf1, Orientation::Forward)];
            let lv2: Vec<Shape> = vec![bat::oriented(&cur_vertex, Orientation::Forward)];
            // OCCT L2779-2790: BRepTools_Substitution theSubstitutor
            // (GAP carrier; TKTopAlgo/BRepTools not translated).
            let mut the_substitutor = BRepToolsSubstitution::new();
            the_substitutor.substitute(&bat::oriented(&v1, Orientation::Forward), &lv1);
            if !v1.is_same(&v2) {
                the_substitutor.substitute(&bat::oriented(&v2, Orientation::Forward), &lv2);
            }
            the_substitutor.substitute(&bat::oriented(&old_edge, Orientation::Forward), &new_edges);
            the_substitutor.build(&new_spherical_face);
            if the_substitutor.is_copied(&new_spherical_face) {
                let list_sh = the_substitutor.copy(&new_spherical_face);
                new_spherical_face = list_sh[0].clone();
            }

            // Adding NewSphericalFace to the shell
            let mut explo =
                bat::explorer(&self.my_offset_shape, ShapeType::Shell, ShapeType::Shape)
                    .into_iter();
            let the_shell = explo.next().expect("myOffsetShape has no shell");
            let mut the_shell = the_shell;
            bat::builder_set_free(&mut the_shell, true);
            let mut bb = rcad_kernel::topo::topods::BRepBuilder::new();
            bb.add_to_shell(&mut self.my_brep, the_shell.clone(), new_spherical_face);
            is_modified = true;
            if !bat::shape_is_closed(&the_shell) {
                if super::brep_offset_make_offset::brep_tool_is_closed(&the_shell) {
                    bat::builder_set_closed(&mut the_shell, true);
                }
            }
        }
        //
        if !is_modified {
            return;
        }
        //
        if self.my_shape.shape_type() == ShapeType::Solid || self.my_thickening {
            let mut nb_shell = 0;
            let mut bb = rcad_kernel::topo::topods::BRepBuilder::new();
            let mut nc = bb.make_compound(&mut self.my_brep, vec![]);
            let mut nc_count = 0usize;
            let mut s1 = Shape::null();

            // OCCT L2862-2864: TopoDS_Solid Sol; BB.MakeSolid(Sol);
            // Sol.Closed(true) — the shells-Vec form (module-b arch. note).
            let mut sol_shells: Vec<Shape> = Vec::new();
            for explo in
                bat::explorer(&self.my_offset_shape, ShapeType::Shell, ShapeType::Shape)
            {
                let mut sh = explo;
                nb_shell += 1;
                if bat::shape_is_closed(&sh) {
                    sol_shells.push(sh.clone());
                } else {
                    bb.add_to_compound(&mut self.my_brep, nc.clone(), sh.clone());
                    nc_count += 1;
                    if nb_shell == 1 {
                        s1 = sh.clone();
                    }
                }
            }
            let sol = bb.make_solid(&mut self.my_brep, sol_shells);
            let mut sol = sol;
            bat::builder_set_closed(&mut sol, true);
            let sol_children = bat::sub_shapes(&sol);
            let nbs = sol_children.len();
            let mut sol_is_null = nbs == 0;
            // Checking solid
            if nbs > 1 {
                // OCCT L2890-2891: BRepCheck_Analyzer aCheck(Sol, false).
                let a_check = super::brep_offset_make_offset::BRepCheckAnalyzer::new(&sol, false);
                if !a_check.is_valid() {
                    let mut a_sol_list: Vec<Shape> = Vec::new();
                    let mut sol_for_check = sol.clone();
                    super::brep_offset_make_offset_d::correct_solid(
                        &mut self.my_brep,
                        &mut sol_for_check,
                        &mut a_sol_list,
                    );
                    sol = sol_for_check;
                    if !a_sol_list.is_empty() {
                        bb.add_to_compound(&mut self.my_brep, nc.clone(), sol.clone());
                        nc_count += 1;
                        for a_sl_it in &a_sol_list {
                            bb.add_to_compound(&mut self.my_brep, nc.clone(), a_sl_it.clone());
                            nc_count += 1;
                        }
                        sol_is_null = true;
                    }
                }
            }
            //
            let nc_is_null = nc_count == 0;
            if !sol_is_null && !nc_is_null {
                bb.add_to_compound(&mut self.my_brep, nc.clone(), sol);
                self.my_offset_shape = nc;
            } else if sol_is_null && !nc_is_null {
                if nb_shell == 1 {
                    self.my_offset_shape = s1;
                } else {
                    self.my_offset_shape = nc;
                }
            } else if !sol_is_null && nc_is_null {
                self.my_offset_shape = sol;
            } else {
                self.my_offset_shape = nc;
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::Intersection3D (cxx L2855-2939).
    pub(crate) fn intersection_3d(&mut self, inter: &mut BRepOffsetInter3d) {
        // OCCT L2873-2875: Message_ProgressScope aPS(...) — the flattened
        // rcad scope (architecture difference #40).
        let a_prog = rcad_kernel::core::message::NoopProgress;
        let mut a_ps = rcad_kernel::core::message::ProgressScope::new(&a_prog, "Intersection 3D", 1);

        // In the Complete Intersection mode, implemented currently for planar
        // solids only, there is no need to intersect the faces here.
        // This intersection will be performed in the method BuildShellsCompleteInter
        // where the special treatment is applied to produced faces.
        //
        // Make sure to match the parameters in which the method
        // BuildShellsCompleteInter is called.
        if self.my_inter
            && (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_is_planar
            && !self.my_thickening
            && self.my_faces.is_empty()
            && super::brep_offset_make_offset_d::is_solid(&self.my_shape)
        {
            return;
        }

        let mut offset_faces: Vec<Shape> = Vec::new(); // list of faces // created.
        super::brep_offset_make_offset::make_list(
            &mut offset_faces,
            &self.my_init_offset_face,
            &self.my_faces,
        );

        if !self.my_faces.is_empty() {
            let in_side = self.my_offset < 0.; // Temporary
            // it is necessary to calculate Inside taking account of the concavity or convexity of edges
            // between the cap and the part.

            if self.my_join == GeomAbsJoinType::Arc {
                inter.context_int_by_arc(
                    &self.my_faces,
                    in_side,
                    &self.my_analyse,
                    &mut self.my_init_offset_face,
                    &mut self.my_init_offset_edge,
                    &a_ps,
                );
            }
        }
        if self.my_inter {
            //-------------
            // Complete.
            //-------------
            inter.complet_int(&offset_faces, &self.my_init_offset_face, &a_ps);
            if !a_ps_more() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            if self.my_join == GeomAbsJoinType::Intersection {
                super::brep_offset_tool_d::correct_orientation(
                    &self.my_face_comp,
                    inter.new_edges(),
                    &mut self.my_as_des,
                    &mut self.my_init_offset_face,
                    self.my_offset,
                );
            }
        } else {
            //--------------------------------
            // Only between neighbor faces.
            //--------------------------------
            inter.connex_int_by_arc(
                &offset_faces,
                &self.my_face_comp,
                &self.my_analyse,
                &self.my_init_offset_face,
                &a_ps,
            );
            if !a_ps_more() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::Intersection2D (cxx L2941-2989).
    pub(crate) fn intersection_2d(
        &mut self,
        modif: &OcctIndexedShapeMap,
        new_edges: &super::brep_offset_inter2d::IndexedShapeMap,
    ) {
        //--------------------------------------------------------
        // calculate intersections2d on faces concerned by
        // intersection3d
        //---------------------------------------------------------
        //-----------------------------------------------
        // Intersection of edges 2 by 2.
        //-----------------------------------------------
        let mut a_dmvv: DmvvMap = indexmap::IndexMap::new();
        for i in 1..=modif.extent() {
            if a_ps_more() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let f = modif.at_1(i).clone();
            BRepOffsetInter2d::compute(
                &mut self.my_as_des,
                &f,
                new_edges,
                self.my_tol,
                &self.my_edge_int_edges,
                &mut a_dmvv,
                (),
            );
        }
        //
        let _ = BRepOffsetInter2d::fuse_vertices(&a_dmvv, &mut self.my_as_des, &mut self.my_image_vv);
        //
    }

    /// OCCT BRepOffset_MakeOffset::MakeLoops (cxx L2991-3057).
    pub(crate) fn make_loops(&mut self, modif: &OcctIndexedShapeMap) {
        let mut lf: Vec<Shape> = Vec::new();
        let mut lc: Vec<Shape> = Vec::new();
        //-----------------------------------------
        // unwinding of faces // modified.
        //-----------------------------------------
        for i in 1..=modif.extent() {
            if !self.my_faces.contains(modif.at_1(i)) {
                lf.push(modif.at_1(i).clone());
            }
        }
        //
        if (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_inter
            && self.my_is_planar
        {
            // OCCT MakeOffset_1.cxx L9497-9508: the BuildSplitsOfTrimmedFaces
            // wrapper — the local BRepOffset_BuildOffsetFaces tool.
            // [INTERFACE NOTE — REPORTED] the SetAsDesInfo shared-handle form
            // (architecture difference #39; see brep_offset_make_offset_b.rs).
            let as_des_handle =
                std::rc::Rc::new(std::cell::RefCell::new(BRepAlgoAsDes::new()));
            let mut a_bf_tool =
                super::brep_offset_make_offset_1::BRepOffsetBuildOffsetFaces::new();
            a_bf_tool.set_faces(&lf);
            a_bf_tool.set_as_des_info(&as_des_handle);
            let a_prog_1 = rcad_kernel::core::message::NoopProgress;
            let mut a_ps_1 = rcad_kernel::core::message::ProgressScope::new(
                &a_prog_1,
                "BuildSplitsOfTrimmedFaces",
                lf.len(),
            );
            a_bf_tool.build_splits_of_trimmed_faces(&mut self.my_image_offset, &a_ps_1);
        } else {
            self.my_make_loops.build(
                &mut lf,
                &self.my_as_des,
                &self.my_image_offset,
                &mut self.my_image_vv,
            );
        }
        if a_ps_more() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }

        //-----------------------------------------
        // unwinding of caps.
        //-----------------------------------------
        for i in 1..=self.my_faces.extent() {
            lc.push(self.my_faces.at_1(i).clone());
        }

        let in_side = self.my_offset <= 0.;
        self.my_make_loops.build_on_context(
            &mut lc,
            &self.my_analyse,
            &self.my_as_des,
            &self.my_image_offset,
            in_side,
        );
    }

    /// OCCT BRepOffset_MakeOffset::MakeFaces (cxx L3059-3108).
    pub(crate) fn make_faces(&mut self, _modif: &OcctIndexedShapeMap) {
        let roots = self.my_init_offset_face.roots().to_vec();
        let mut lof: Vec<Shape> = Vec::new();
        //----------------------------------
        // Loop on all faces //.
        //----------------------------------
        for itr in &roots {
            let f = self.my_init_offset_face.image(itr)[0].clone();
            if !self.my_image_offset.has_image(&f) {
                lof.push(f);
            }
        }
        //
        if (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_inter
            && self.my_is_planar
        {
            // OCCT MakeOffset_1.cxx L9497-9508: the BuildSplitsOfTrimmedFaces
            // wrapper — the local BRepOffset_BuildOffsetFaces tool.
            // [INTERFACE NOTE — REPORTED] the SetAsDesInfo shared-handle form
            // (architecture difference #39; see brep_offset_make_offset_b.rs).
            let as_des_handle =
                std::rc::Rc::new(std::cell::RefCell::new(BRepAlgoAsDes::new()));
            let mut a_bf_tool =
                super::brep_offset_make_offset_1::BRepOffsetBuildOffsetFaces::new();
            a_bf_tool.set_faces(&lof);
            a_bf_tool.set_as_des_info(&as_des_handle);
            let a_prog_1 = rcad_kernel::core::message::NoopProgress;
            let mut a_ps_1 = rcad_kernel::core::message::ProgressScope::new(
                &a_prog_1,
                "BuildSplitsOfTrimmedFaces",
                lof.len(),
            );
            a_bf_tool.build_splits_of_trimmed_faces(&mut self.my_image_offset, &a_ps_1);
        } else {
            self.my_make_loops
                .build_faces(&mut lof, &self.my_as_des, &self.my_image_offset);
        }
        if a_ps_more() {
            self.my_error = BRepOffset_Error::UserBreak;
            return;
        }
    }
}

// ===========================================================================
// Local helpers of module c.
// ===========================================================================

/// OCCT ElSLib::Parameters(theSphere, P) — the rcad
/// elslib_sphere_parameters re-host keyed by the sphere frame
/// (center / ref_dir / axis).
pub(crate) fn elslib_sphere_parameters_host(
    p: glam::DVec3,
    center: glam::DVec3,
    ref_dir: glam::DVec3,
    axis: glam::DVec3,
) -> (f64, f64) {
    let y_dir = axis.cross(ref_dir);
    rcad_kernel::math::el::elslib_sphere_parameters(p, center, ref_dir, y_dir, axis)
}

/// OCCT Geom_SphericalSurface::VIso(V) (via ElSLib::SphereVIso) — GAP leaf
/// (architecture difference #52: the TKMath iso constructors are not
/// translated; the bi_tgte #28 precedent).
pub(crate) fn elslib_sphere_viso(
    _center: glam::DVec3,
    _axis: glam::DVec3,
    _ref_dir: glam::DVec3,
    _radius: f64,
    _v: f64,
) -> rcad_kernel::geom::Curve3 {
    panic!("GAP: ElSLib::SphereVIso / Geom_SphericalSurface::VIso (TKMath not translated)");
}

/// OCCT BRepLib_MakeFace(Surface, U1, U2, V1, V2, Tol) — GAP leaf
/// (architecture difference #52): the UV-bounds face maker of
/// CorrectConicalFaces.
#[allow(clippy::too_many_arguments)]
pub(crate) fn brep_lib_make_face_uv(
    _center: glam::DVec3,
    _axis: glam::DVec3,
    _ref_dir: glam::DVec3,
    _radius: f64,
    _has_position: bool,
    _u1: f64,
    _u2: f64,
    _v1: f64,
    _v2: f64,
    _tol: f64,
) -> Shape {
    panic!("GAP: BRepLib_MakeFace(S, U1, U2, V1, V2, Tol) (TKTopAlgo/BRepLib not translated)");
}

/// OCCT BRepTools_Substitution (TKTopAlgo/BRepTools) — GAP carrier
/// (architecture difference #52): the Substitutor of CorrectConicalFaces.
pub(crate) struct BRepToolsSubstitution;

impl BRepToolsSubstitution {
    /// OCCT BRepTools_Substitution::BRepTools_Substitution().
    pub fn new() -> Self {
        BRepToolsSubstitution
    }

    /// OCCT BRepTools_Substitution::Substitute(S, L).
    pub fn substitute(&mut self, _s: &Shape, _l: &[Shape]) {
        panic!("GAP: BRepTools_Substitution::Substitute (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Substitution::Build(S).
    pub fn build(&mut self, _s: &Shape) {
        panic!("GAP: BRepTools_Substitution::Build (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Substitution::IsCopied(S).
    pub fn is_copied(&self, _s: &Shape) -> bool {
        panic!("GAP: BRepTools_Substitution::IsCopied (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Substitution::Copy(S).
    pub fn copy(&self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepTools_Substitution::Copy (TKTopAlgo/BRepTools not translated)");
    }
}

/// The OCCT `!aPS.More()` user-break probe — the flattened NoopProgress
/// scope never breaks (architecture difference #40); the helper keeps the
/// OCCT branch structure at every former `aPS.More()` site without a scope.
fn a_ps_more() -> bool {
    false
}
