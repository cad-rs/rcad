// OCCT BRepOffset_MakeOffset.cxx — 1:1 translation, module e (see the
// module header of brep_offset_make_offset.rs for the split map and the
// architecture-difference list #38-#56).
//
// Module e carries cxx L3148-L4750 + L5231-L5659: MakeMissingWalls /
// MakeShells / MakeSolid / SelectShells / EncodeRegularity / CheckInputData /
// RemoveInternalEdges / IntersectEdges / Generated / Modified / IsDeleted /
// analyzeProgress / IsPlanar.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::core::precision::{
    CONFUSION as PRECISION_CONFUSION, PCONFUSION as PRECISION_PCONFUSION,
    SQUARE_CONFUSION as PRECISION_SQUARE_CONFUSION,
};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topo::topods::{
    BRep, BRepBuilder, GeomAbsShape, Orientation, ShapeType, State, TShape,
};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_inter2d::{DmvvMap, IndexedShapeMap};
use super::brep_offset_inter2d_b::BRepOffsetInter2d;
use super::brep_offset_make_offset::{
    brep_lib_build_curve3d_edge, brep_lib_build_curves3d_tol, brep_tools_uv_bounds,
    indexed_shape_map_remove_key, orientation_of_edge_in_face,
    top_exp_map_shapes_and_ancestors, BRepOffset_Error,
    BRepOffset_PIOperation, BRepOffsetMakeOffset, BRepToolsQuilt, DataMapOfShapeListOfShape,
    DataMapOfShapeShape, IndexedDataMapOfShapeListOfShape, MapSF,
};
use super::brep_offset_make_offset_c::BRepAdaptorCurveC;
use super::brep_offset_tool::{
    find_common_shapes, set_add, set_contains, shape_data_map, top_exp_vertices,
    OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap,
};
use super::brep_offset_tool_c::deboucle3d;
use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool as bat;
use crate::brep_fill::offset_wire::GeomAbsJoinType;
use crate::brep_algo::tool::{brep_tool_pnt, shape_key};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve, brep_tool_degenerated, brep_tool_range, brep_tool_tolerance,
};
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;

impl BRepOffsetMakeOffset {
    /// OCCT BRepOffset_MakeOffset::MakeMissingWalls (cxx L3148-3689).
    pub(crate) fn make_missing_walls(&mut self) {
        // OCCT L3151: Contours — Start vertex + list of connected edges
        // (free boundary).
        let mut contours = IndexedDataMapOfShapeListOfShape::new();
        // OCCT L3153-3154: MapEF — Edges of contours: edge + face.
        let mut map_ef = DataMapOfShapeShape::new();
        let offset_val = self.my_offset.abs();

        super::brep_offset_make_offset::fill_contours(
            &self.my_face_comp,
            &self.my_analyse,
            &mut contours,
            &mut map_ef,
        );

        // OCCT L3162: Message_ProgressScope aPS(theRange, "Making missing
        // walls", Contours.Extent()) — the flattened rcad scope (arch.
        // diff. #40).
        let a_prog = NoopProgress;
        let mut a_ps = ProgressScope::new(&a_prog, "Making missing walls", contours.len());
        for ic in 1..=contours.len() {
            a_ps.advance();
            if a_ps.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let start_vertex = super::brep_offset_make_offset::shape_indexed_data_map_find_key_1(
                &contours, ic,
            );
            let mut start_edge = Shape::null();
            let a_contour = super::brep_offset_make_offset::shape_indexed_data_map_value_1(
                &contours, ic,
            );
            let mut first_step = true;
            let mut prev_edge = Shape::null();
            let mut prev_vertex = start_vertex.clone();
            let mut is_build_from_scratch = false; // Problems with edges.
            for itl in a_contour.clone() {
                let mut an_edge = itl;
                let a_face_of_edge = shape_data_map::value(&map_ef, &an_edge).clone();

                // Check for offset existence.
                if !self.my_init_offset_edge.has_image(&an_edge) {
                    continue;
                }

                // Check for existence of two different vertices.
                let mut loe: Vec<Shape> = Vec::new();
                let mut loe2: Vec<Shape> = Vec::new();
                self.my_init_offset_edge.last_image(&an_edge, &mut loe);
                self.my_image_offset.last_image(loe.last().unwrap(), &mut loe2);
                let mut oe = loe2.last().unwrap().clone();
                let (v4_0, v3_0) = top_exp_vertices(&oe);
                let (mut v4, mut v3) = (v4_0, v3_0);
                let (v1_0, v2_0) = top_exp_vertices(&an_edge);
                let (mut v1, mut v2) = (v1_0, v2_0);
                let (a_f, a_l) = brep_tool_range(&an_edge);
                let a_c = brep_tool_curve(&an_edge);
                if v3.is_null() && v4.is_null() {
                    // Initially offset edge is created without vertices.
                    // Then edge is trimmed by intersection line between
                    // two adjacent extended offset faces and get vertices.
                    // When intersection lines are invalid for any reason,
                    // (one reason is mixed connectivity of faces)
                    // algorithm of cutting offset edge by intersection line
                    // can fail and offset edge cannot get vertices.
                    // Following workaround is only to avoid exception if V3 and V4 are Null
                    // Vertex points are invalid.
                    let an_oe_ori = oe.orientation;
                    oe.orientation = Orientation::Forward;
                    let an_oec = brep_tool_curve(&oe);
                    let mut a_bb = BRepBuilder::new();
                    let a_p1 = an_oec
                        .as_ref()
                        .map(|(c, _, _)| c.point_at(a_f))
                        .expect("BRep_Tool::Curve null");
                    let a_p2 = an_oec
                        .as_ref()
                        .map(|(c, _, _)| c.point_at(a_l))
                        .expect("BRep_Tool::Curve null");
                    let a_tol = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));
                    // OCCT L3244-3246: aBB.MakeVertex(anOEV1, aP1, aTol).
                    let mut an_oev1 = a_bb.add_vertex(&mut self.my_brep, a_p1, a_tol);
                    an_oev1.orientation = Orientation::Forward;
                    let mut an_oev2 = a_bb.add_vertex(&mut self.my_brep, a_p2, a_tol);
                    an_oev2.orientation = Orientation::Reversed;
                    bat::builder_add_edge_vertex(&mut oe, &an_oev1);
                    bat::builder_add_edge_vertex(&mut oe, &an_oev2);
                    bat::builder_range_edge(&mut oe, a_f, a_l);
                    oe.orientation = an_oe_ori;
                    let (v4_1, v3_1) = top_exp_vertices(&oe);
                    v4 = v4_1;
                    v3 = v3_1;
                }
                if let Some((a_c_v, _, _)) = &a_c {
                    // OCCT L3253: !aC->IsClosed() && !aC->IsPeriodic() — the
                    // rcad Curve3 carries no closed/periodic flags; the
                    // probe reduces to the parameter-range degeneracy of the
                    // extremities (arch. diff. #57).
                    let a_pnt_f = brep_tool_pnt(&v1).expect("BRep_Tool::Pnt null");
                    let a_pnt_l = brep_tool_pnt(&v2).expect("BRep_Tool::Pnt null");
                    let a_dist_e = a_pnt_f.distance_squared(a_pnt_l);
                    if a_dist_e < PRECISION_SQUARE_CONFUSION {
                        // Bad case: non closed, but vertexes mapped to same 3d point.
                        continue;
                    }

                    let an_edge_tol = brep_tool_tolerance(&an_edge);
                    if a_dist_e < an_edge_tol {
                        // Potential problems not detected via checkshape.
                        let a_pnt_of = brep_tool_pnt(&v4).expect("BRep_Tool::Pnt null");
                        let a_pnt_ol = brep_tool_pnt(&v3).expect("BRep_Tool::Pnt null");
                        if a_pnt_of.distance_squared(a_pnt_ol) > GP_RESOLUTION {
                            // To avoid computation of complex analytical continuation of Sin / ArcSin.
                            let a_sin_value =
                                (2. * an_edge_tol / a_pnt_of.distance(a_pnt_ol)).min(1.0);
                            let a_max_angle =
                                a_sin_value.asin().abs().min(std::f64::consts::PI / 4.); // Maximal angle.
                            let a_current_angle =
                                (a_pnt_l - a_pnt_f).angle_between(a_pnt_ol - a_pnt_of);
                            let is_line = matches!(a_c_v, Curve3::Line(_));
                            if is_line && a_current_angle.abs() > a_max_angle {
                                // anEdge not collinear to offset edge.
                                is_build_from_scratch = true;
                                self.my_is_perform_sewing = true;
                                continue;
                            }
                        }
                    }
                }

                let mut to_reverse = false;
                if !v1.is_same(&prev_vertex) {
                    std::mem::swap(&mut v1, &mut v2);
                    std::mem::swap(&mut v3, &mut v4);
                    to_reverse = true;
                }

                // OCCT L3303: OE.Orientation(TopAbs::Reverse(anEdge.Orientation())).
                oe.orientation = bat::top_abs_reverse(an_edge.orientation);
                let mut e3 = Shape::null();
                let mut e4 = Shape::null();
                let arc_on_v2 = (self.my_join == GeomAbsJoinType::Arc)
                    && self.my_init_offset_edge.has_image(&v2);
                if first_step || is_build_from_scratch {
                    // OCCT L3309: BRepLib_MakeEdge(V1, V4) — the vertex-pair
                    // edge form (no 3D curve, the [0, 1] parameter range of
                    // the OCCT BRepLib_MakeEdge(V, V) ctor).
                    e4 = BRepBuilder::new().add_edge(&mut self.my_brep, None, v1.clone(), v4.clone(), [0., 1.]);
                    if first_step {
                        start_edge = e4.clone();
                    }
                } else {
                    e4 = prev_edge.clone();
                }
                if v2.is_same(&start_vertex) && !arc_on_v2 {
                    e3 = start_edge.clone();
                } else {
                    e3 = BRepBuilder::new().add_edge(&mut self.my_brep, None, v2.clone(), v3.clone(), [0., 1.]);
                }
                // OCCT L3318: E4.Reverse().
                e4 = bat::reversed(&e4);

                if is_build_from_scratch {
                    e3 = bat::reversed(&e3);
                    e4 = bat::reversed(&e4);
                }

                let an_edge_fwd = bat::oriented(&an_edge, Orientation::Forward);
                let par_v1 = bat::brep_tool_parameter(&v1, &an_edge_fwd);
                let par_v2 = bat::brep_tool_parameter(&v2, &an_edge_fwd);
                let mut bb = BRepBuilder::new();
                let mut the_wire = bb.make_wire(&mut self.my_brep);
                if to_reverse {
                    bat::builder_add_wire_edge(&mut the_wire, &bat::reversed(&an_edge));
                    bat::builder_add_wire_edge(&mut the_wire, &bat::reversed(&e3));
                    bat::builder_add_wire_edge(&mut the_wire, &bat::reversed(&oe));
                    bat::builder_add_wire_edge(&mut the_wire, &bat::reversed(&e4));
                } else {
                    bat::builder_add_wire_edge(&mut the_wire, &an_edge);
                    bat::builder_add_wire_edge(&mut the_wire, &e3);
                    bat::builder_add_wire_edge(&mut the_wire, &oe);
                    bat::builder_add_wire_edge(&mut the_wire, &e4);
                }

                brep_lib_build_curves3d_tol(&the_wire, self.my_tol);
                bat::builder_set_closed(&mut the_wire, true);
                let mut new_face = Shape::null();
                let mut the_surf: Option<Surface3> = None;
                let ba_curve = BRepAdaptorCurveC::new(&an_edge);
                let ba_curve_oe = BRepAdaptorCurveC::new(&oe);
                let fpar = ba_curve.first_parameter();
                let lpar = ba_curve.last_parameter();
                let pon_e = ba_curve.value(fpar);
                let pon_oe = ba_curve_oe.value(fpar);
                let offset_dir =
                    super::brep_offset_make_offset::gce_make_dir(pon_e, pon_oe);
                // OCCT L3328-3329: EdgeLine2d / OELine2d / aLine2d /
                // aLine2d2 — the Option<Curve2d> carriers.
                let mut edge_line2d: Option<Curve2d> = None;
                let mut oe_line2d: Option<Curve2d> = None;
                let mut a_line2d: Option<Curve2d> = None;
                let mut a_line2d2: Option<Curve2d> = None;
                let mut is_planar = false;
                let a_circ = ba_curve.get_type();
                let a_circ_oe_t = ba_curve_oe.get_type();
                if matches!(a_circ, Curve3::Circle(_)) && matches!(a_circ_oe_t, Curve3::Circle(_))
                {
                    let a_circ = ba_curve.circle();
                    let a_circ_oe = ba_curve_oe.circle();
                    // OCCT L3345: gp_Lin anAxisLine(aCirc.Axis()) — the axis
                    // line (center, normal) of the circle.
                    let an_axis_line_origin = a_circ.center;
                    let an_axis_line_dir = a_circ.normal;
                    let mut circ_axis_dir = a_circ.normal;
                    let axes_parallel = {
                        let ang = an_axis_line_dir.angle_between(a_circ_oe.normal);
                        ang <= PRECISION_CONFUSION
                            || (std::f64::consts::PI - ang) <= PRECISION_CONFUSION
                    };
                    let axis_line_contains = {
                        // OCCT gp_Lin::Contains(P, Tol) — the point-line
                        // distance form.
                        let d = (a_circ_oe.center - an_axis_line_origin)
                            .cross(an_axis_line_dir)
                            .length();
                        d <= PRECISION_CONFUSION
                    };
                    if axes_parallel && axis_line_contains {
                        // cylinder, plane or cone
                        if (a_circ.radius - a_circ_oe.radius).abs() <= PRECISION_CONFUSION {
                            // case of cylinder
                            the_surf = Some(
                                super::brep_offset_make_offset::gc_make_cylindrical_surface(
                                    &a_circ,
                                ),
                            );
                        } else if a_circ.center.distance(a_circ_oe.center)
                            <= PRECISION_CONFUSION
                        {
                            // case of plane
                            is_planar = true;
                            //
                            let pon_el = ba_curve.value(lpar);
                            if pon_el.distance(pon_e) <= PRECISION_PCONFUSION {
                                let mut b_is_hole;
                                // OCCT L3378-3385: aE1/aE2, aW1/aW2, aPL,
                                // IntTools_FClass2d aClsf [D6 JUDGMENT
                                // REQUIRED: TKBool/IntTools].
                                let (a_e1, a_e2) = if a_circ.radius > a_circ_oe.radius {
                                    (an_edge.clone(), oe.clone())
                                } else {
                                    (oe.clone(), an_edge.clone())
                                };
                                //
                                let mut a_w1 = bb.make_wire(&mut self.my_brep);
                                bat::builder_add_wire_edge(&mut a_w1, &a_e1);
                                let mut a_w2 = bb.make_wire(&mut self.my_brep);
                                bat::builder_add_wire_edge(&mut a_w2, &a_e2);
                                //
                                // OCCT L3392: aPL = new Geom_Plane(aCirc.Location(), CircAxisDir).
                                let a_pl = Surface3::Plane(rcad_kernel::geom::Plane {
                                    origin: a_circ.center,
                                    normal: circ_axis_dir,
                                    u_dir: plane_u_dir_of(circ_axis_dir),
                                    v_dir: circ_axis_dir.cross(plane_u_dir_of(circ_axis_dir)),
                                });
                                for i in 0..2 {
                                    let (a_w, a_e) = if i == 0 {
                                        (&mut a_w1, &a_e1)
                                    } else {
                                        (&mut a_w2, &a_e2)
                                    };
                                    //
                                    // OCCT L3397: BB.MakeFace(aFace, aPL,
                                    // Precision::Confusion()).
                                    let mut a_face =
                                        bb.make_face(&mut self.my_brep, Some(a_pl.clone()), Shape::null());
                                    bat::builder_add_face_wire(&mut a_face, a_w);
                                    let mut a_clsf =
                                        super::brep_offset_make_offset::IntToolsFClass2d;
                                    a_clsf.init(&a_face, PRECISION_CONFUSION);
                                    b_is_hole = a_clsf.is_hole();
                                    if (b_is_hole && i == 0) || (!b_is_hole && i == 1) {
                                        // OCCT L3345: aW.Nullify(); BB.MakeWire(aW);
                                        // BB.Add(aW, aE.Reversed()).
                                        *a_w = bb.make_wire(&mut self.my_brep);
                                        bat::builder_add_wire_edge(a_w, &bat::reversed(a_e));
                                    }
                                    let _ = b_is_hole;
                                    b_is_hole = false;
                                }
                                //
                                let mut new_face_new = bb.make_face(
                                    &mut self.my_brep,
                                    Some(a_pl.clone()),
                                    Shape::null(),
                                );
                                bat::builder_add_face_wire(&mut new_face_new, &a_w1);
                                bat::builder_add_face_wire(&mut new_face_new, &a_w2);
                                new_face = new_face_new;
                            }
                        } else {
                            // case of cone
                            // OCCT L3427-3428: gce_MakeCone(L1, L2, R1, R2)
                            // (GAP leaf, arch. diff. #52).
                            let mut the_cone = super::brep_offset_make_offset::gce_make_cone(
                                a_circ.center,
                                a_circ_oe.center,
                                a_circ.radius,
                                a_circ_oe.radius,
                            );
                            // OCCT L3429: gp_Ax3 theAx3(aCirc.Position()).
                            let mut cone_axis = the_cone.axis;
                            if circ_axis_dir.dot(the_cone.axis) < 0. {
                                // OCCT L3431: theAx3.ZReverse().
                                cone_axis = -cone_axis;
                                the_cone.axis = cone_axis;
                                // OCCT L3433: CircAxisDir.Reverse().
                                circ_axis_dir = -circ_axis_dir;
                            }
                            the_surf = Some(Surface3::Cone(the_cone));
                        }
                        if !is_planar {
                            // OCCT L3445: EdgeLine2d = new Geom2d_Line((0.,
                            // 0.), DX2d).
                            edge_line2d = Some(line2d_of(0., 0., 1., 0.));
                            bb_update_edge_pcurve_on_surface(
                                &mut an_edge,
                                edge_line2d.as_ref().unwrap(),
                                the_surf.as_ref().unwrap(),
                                PRECISION_CONFUSION,
                            );
                            let coeff = if offset_dir.dot(circ_axis_dir) > 0. {
                                1.
                            } else {
                                -1.
                            };
                            oe_line2d = Some(line2d_of(0., offset_val * coeff, 1., 0.));
                            bb_update_edge_pcurve_on_surface(
                                &mut oe,
                                oe_line2d.as_ref().unwrap(),
                                the_surf.as_ref().unwrap(),
                                PRECISION_CONFUSION,
                            );
                            a_line2d = Some(line2d_of(par_v2, 0., 0., coeff));
                            a_line2d2 = Some(line2d_of(par_v1, 0., 0., coeff));
                            if e3.is_same(&e4) {
                                if coeff > 0. {
                                    bb_update_edge_pcurves_two_on_surface(
                                        &mut e3,
                                        a_line2d.as_ref().unwrap(),
                                        a_line2d2.as_ref().unwrap(),
                                        the_surf.as_ref().unwrap(),
                                        PRECISION_CONFUSION,
                                    );
                                } else {
                                    bb_update_edge_pcurves_two_on_surface(
                                        &mut e3,
                                        a_line2d2.as_ref().unwrap(),
                                        a_line2d.as_ref().unwrap(),
                                        the_surf.as_ref().unwrap(),
                                        PRECISION_CONFUSION,
                                    );
                                    // OCCT L3466-3472: theWire.Nullify();
                                    // BB.MakeWire(theWire); ... theWire.Closed(true).
                                    the_wire = bb.make_wire(&mut self.my_brep);
                                    bat::builder_add_wire_edge(
                                        &mut the_wire,
                                        &bat::oriented(&an_edge, Orientation::Reversed),
                                    );
                                    bat::builder_add_wire_edge(&mut the_wire, &e4);
                                    bat::builder_add_wire_edge(
                                        &mut the_wire,
                                        &bat::oriented(&oe, Orientation::Forward),
                                    );
                                    bat::builder_add_wire_edge(&mut the_wire, &e3);
                                    bat::builder_set_closed(&mut the_wire, true);
                                }
                            } else {
                                bb.set_edge_same_parameter(&mut self.my_brep, e3.clone(), false);
                                bb.set_edge_same_range(&mut self.my_brep, e3.clone(), false);
                                bb.set_edge_same_parameter(&mut self.my_brep, e4.clone(), false);
                                bb.set_edge_same_range(&mut self.my_brep, e4.clone(), false);
                                bb_update_edge_pcurve_on_surface(
                                    &mut e3,
                                    a_line2d.as_ref().unwrap(),
                                    the_surf.as_ref().unwrap(),
                                    PRECISION_CONFUSION,
                                );
                                bb_range_edge_on_surface(&mut e3, 0., offset_val);
                                bb_update_edge_pcurve_on_surface(
                                    &mut e4,
                                    a_line2d2.as_ref().unwrap(),
                                    the_surf.as_ref().unwrap(),
                                    PRECISION_CONFUSION,
                                );
                                bb_range_edge_on_surface(&mut e4, 0., offset_val);
                            }
                            // OCCT L3490: NewFace = BRepLib_MakeFace(theSurf,
                            // theWire).
                            new_face = bb.make_face(
                                &mut self.my_brep,
                                the_surf.clone(),
                                the_wire.clone(),
                            );
                        }
                    } // cylinder or cone
                } // if both edges are arcs of circles
                if new_face.is_null() {
                    let an_edge_tol = brep_tool_tolerance(&an_edge);
                    // Tolerances of input shape should not be increased by BRepLib_MakeFace
                    // OCCT L3498: BRepLib_FindSurface aFindPlane(theWire,
                    // anEdgeTol, true) (GAP carrier, arch. diff. #52).
                    let mut a_find_plane = super::brep_offset_make_offset::BRepLibFindSurface;
                    a_find_plane.init(&the_wire, an_edge_tol, true); // only plane
                    is_planar = false;
                    if a_find_plane.found() && a_find_plane.tolerance_reached() <= an_edge_tol {
                        let a_gc = brep_tool_curve(&an_edge);
                        let a_pln = match a_find_plane.surface() {
                            Surface3::Plane(p) => p,
                            _ => panic!("BRepLib_FindSurface::Surface is not a plane"),
                        };
                        let a_max_dist = match &a_gc {
                            Some((c, f, l)) => super::brep_offset_make_offset_d::compute_max_dist(
                                &a_pln, c, *f, *l,
                            ),
                            None => f64::INFINITY,
                        };
                        if a_max_dist <= an_edge_tol {
                            // OCCT L3509: BRepLib_MakeFace MF(aPln->Pln(),
                            // theWire).
                            let mf = bb.make_face(
                                &mut self.my_brep,
                                Some(Surface3::Plane(a_pln.clone())),
                                the_wire.clone(),
                            );
                            let mf_done = !mf.is_null();
                            if mf_done {
                                new_face = mf;
                                for an_it_e in bat::sub_shapes(&the_wire) {
                                    let mut an_e = an_it_e;
                                    if an_e.is_same(&an_edge) {
                                        continue;
                                    }
                                    if let Some((a_gc2, f2, l2)) = brep_tool_curve(&an_e) {
                                        let a_max_dist =
                                            super::brep_offset_make_offset_d::compute_max_dist(
                                                &a_pln, &a_gc2, f2, l2,
                                            );
                                        update_edge_tolerance_host_e(&mut an_e, a_max_dist);
                                    }
                                }
                                is_planar = true;
                            }
                        }
                    }
                    //
                    if !is_planar {
                        // Extrusion (by thrusections)
                        // OCCT L3532-3543: the trimmed curves + the
                        // GeomFill_Generator (GAP carrier, arch. diff. #52).
                        let mut thrusec_generator =
                            super::brep_offset_make_offset::GeomFillGenerator;
                        if let Some((c, _, _)) = brep_tool_curve(&an_edge) {
                            thrusec_generator.add_curve(&c);
                        }
                        if let Some((c, _, _)) = brep_tool_curve(&oe) {
                            thrusec_generator.add_curve(&c);
                        }
                        thrusec_generator.perform(PRECISION_PCONFUSION);
                        let the_surf_gen = thrusec_generator.surface();
                        // OCCT L3545: theSurf->Bounds(Uf, Ul, Vf, Vl).
                        let bounds = the_surf_gen.default_domain();
                        let (uf, ul, vf, vl) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                        the_surf = Some(the_surf_gen.clone());
                        edge_line2d = Some(line2d_of(0., vf, 1., 0.));
                        bb_update_edge_pcurve_on_surface(
                            &mut an_edge,
                            edge_line2d.as_ref().unwrap(),
                            &the_surf_gen,
                            PRECISION_CONFUSION,
                        );
                        oe_line2d = Some(line2d_of(0., vl, 1., 0.));
                        bb_update_edge_pcurve_on_surface(
                            &mut oe,
                            oe_line2d.as_ref().unwrap(),
                            &the_surf_gen,
                            PRECISION_CONFUSION,
                        );
                        let uon_v1 = if to_reverse { ul } else { uf };
                        let uon_v2 = if to_reverse { uf } else { ul };
                        a_line2d = Some(line2d_of(uon_v2, 0., 0., 1.));
                        a_line2d2 = Some(line2d_of(uon_v1, 0., 0., 1.));
                        if e3.is_same(&e4) {
                            bb_update_edge_pcurves_two_on_surface(
                                &mut e3,
                                a_line2d.as_ref().unwrap(),
                                a_line2d2.as_ref().unwrap(),
                                &the_surf_gen,
                                PRECISION_CONFUSION,
                            );
                            // OCCT L3570-3572: BSplC34 = theSurf->UIso(Uf);
                            // BB.UpdateEdge(E3, BSplC34, Tol);
                            // BB.Range(E3, Vf, Vl).
                            let bspl_c34 = surface_u_iso(&the_surf_gen, uf);
                            update_edge_3d_host(&mut e3, &bspl_c34, PRECISION_CONFUSION);
                            bat::builder_range_edge(&mut e3, vf, vl);
                        } else {
                            bb.set_edge_same_parameter(&mut self.my_brep, e3.clone(), false);
                            bb.set_edge_same_range(&mut self.my_brep, e3.clone(), false);
                            bb.set_edge_same_parameter(&mut self.my_brep, e4.clone(), false);
                            bb.set_edge_same_range(&mut self.my_brep, e4.clone(), false);
                            bb_update_edge_pcurve_on_surface(
                                &mut e3,
                                a_line2d.as_ref().unwrap(),
                                &the_surf_gen,
                                PRECISION_CONFUSION,
                            );
                            bb_range_edge_on_surface(&mut e3, vf, vl);
                            bb_update_edge_pcurve_on_surface(
                                &mut e4,
                                a_line2d2.as_ref().unwrap(),
                                &the_surf_gen,
                                PRECISION_CONFUSION,
                            );
                            bb_range_edge_on_surface(&mut e4, vf, vl);
                            let bspl_c3 = surface_u_iso(&the_surf_gen, uon_v2);
                            update_edge_3d_host(&mut e3, &bspl_c3, PRECISION_CONFUSION);
                            bat::builder_range_edge(&mut e3, vf, vl); // only for 3d curve
                            let bspl_c4 = surface_u_iso(&the_surf_gen, uon_v1);
                            update_edge_3d_host(&mut e4, &bspl_c4, PRECISION_CONFUSION);
                            bat::builder_range_edge(&mut e4, vf, vl); // only for 3d curve
                        }
                        // OCCT L3594: NewFace = BRepLib_MakeFace(theSurf,
                        // theWire).
                        new_face = bb.make_face(
                            &mut self.my_brep,
                            the_surf.clone(),
                            the_wire.clone(),
                        );
                    }
                }
                if !is_planar {
                    let fpar_oe = ba_curve_oe.first_parameter();
                    let lpar_oe = ba_curve_oe.last_parameter();
                    if (fpar - fpar_oe).abs() > PRECISION_CONFUSION {
                        let mut an_e4 = if to_reverse { e3.clone() } else { e4.clone() };
                        let fp2d = edge_line2d
                            .as_ref()
                            .expect("EdgeLine2d null")
                            .point_at(fpar);
                        let fp2d_oe = oe_line2d
                            .as_ref()
                            .expect("OELine2d null")
                            .point_at(fpar_oe);
                        // OCCT L3604: aLine2d2 = GC_MakeLine2d(fp2d,
                        // fp2dOE).Value() (GAP leaf, arch. diff. #52).
                        let line2d2 =
                            super::brep_offset_make_offset::gc_make_line2d(fp2d, fp2d_oe);
                        a_line2d2 = Some(line2d2.clone());
                        let first_par = 0.;
                        let last_par = fp2d.distance(fp2d_oe);
                        // OCCT L3607-3618: the
                        // Geom2dAdaptor_Curve / GeomAdaptor_Surface /
                        // Adaptor3d_CurveOnSurface + GeomLib::BuildCurve3d
                        // chain (the carriers of arch. diff. #52).
                        let con_s =
                            super::brep_offset_make_offset::Adaptor3dCurveOnSurface::new(
                                &line2d2,
                                the_surf.as_ref().unwrap(),
                            );
                        let mut a_curve: Option<Curve3> = None;
                        let mut max_deviation = 0.;
                        let mut average_deviation = 0.;
                        let _ = (&con_s, first_par, last_par);
                        super::brep_offset_make_offset::geom_lib_build_curve3d_cons(
                            &mut a_curve,
                            &mut max_deviation,
                            &mut average_deviation,
                        );
                        if let Some(a_curve) = a_curve {
                            update_edge_3d_host(&mut an_e4, &a_curve, max_deviation);
                        }
                        bb_update_edge_pcurve_on_surface(
                            &mut an_e4,
                            a_line2d2.as_ref().unwrap(),
                            the_surf.as_ref().unwrap(),
                            max_deviation,
                        );
                        bat::builder_range_edge(&mut an_e4, first_par, last_par);
                    }
                    if (lpar - lpar_oe).abs() > PRECISION_CONFUSION {
                        let mut an_e3 = if to_reverse { e4.clone() } else { e3.clone() };
                        let lp2d = edge_line2d
                            .as_ref()
                            .expect("EdgeLine2d null")
                            .point_at(lpar);
                        let lp2d_oe = oe_line2d
                            .as_ref()
                            .expect("OELine2d null")
                            .point_at(lpar_oe);
                        let line2d =
                            super::brep_offset_make_offset::gc_make_line2d(lp2d, lp2d_oe);
                        a_line2d = Some(line2d.clone());
                        let first_par = 0.;
                        let last_par = lp2d.distance(lp2d_oe);
                        let con_s =
                            super::brep_offset_make_offset::Adaptor3dCurveOnSurface::new(
                                &line2d,
                                the_surf.as_ref().unwrap(),
                            );
                        let mut a_curve: Option<Curve3> = None;
                        let mut max_deviation = 0.;
                        let mut average_deviation = 0.;
                        let _ = (&con_s, first_par, last_par);
                        super::brep_offset_make_offset::geom_lib_build_curve3d_cons(
                            &mut a_curve,
                            &mut max_deviation,
                            &mut average_deviation,
                        );
                        if let Some(a_curve) = a_curve {
                            update_edge_3d_host(&mut an_e3, &a_curve, max_deviation);
                        }
                        bb_update_edge_pcurve_on_surface(
                            &mut an_e3,
                            a_line2d.as_ref().unwrap(),
                            the_surf.as_ref().unwrap(),
                            max_deviation,
                        );
                        bat::builder_range_edge(&mut an_e3, first_par, last_par);
                    }
                }

                if !is_planar {
                    // For planar faces these operations are useless,
                    // because there are no curves on surface
                    super::brep_offset_make_offset::brep_lib_same_parameter_face(&mut new_face);
                    super::brep_offset_make_offset::brep_tools_update(&mut new_face);
                }
                // Check orientation
                let an_or = orientation_of_edge_in_face(&an_edge, &a_face_of_edge);
                let or_in_new_face = orientation_of_edge_in_face(&an_edge, &new_face);
                if or_in_new_face != bat::top_abs_reverse(an_or) {
                    new_face = bat::reversed(&new_face);
                }
                ///////////////////
                self.my_walls.push(new_face);
                if arc_on_v2 {
                    let mut an_arc = self.my_init_offset_edge.image(&v2)[0].clone();
                    let (mut arc_v1, mut arc_v2) = top_exp_vertices(&an_arc);
                    let mut arc_reverse = false;
                    if !arc_v1.is_same(&v3) {
                        std::mem::swap(&mut arc_v1, &mut arc_v2);
                        arc_reverse = true;
                    }
                    let mut ea1 = e3.clone();
                    // OCCT L3642: EA1.Reverse().
                    ea1 = bat::reversed(&ea1);
                    if to_reverse {
                        ea1 = bat::reversed(&ea1);
                    }
                    //////////////////////////////////////////////////////
                    let mut ea2;
                    if v2.is_same(&start_vertex) {
                        ea2 = start_edge.clone();
                    } else {
                        // OCCT L3650: BRepLib_MakeEdge(V2, arcV2).
                        ea2 = BRepBuilder::new().add_edge(
                            &mut self.my_brep,
                            None,
                            v2.clone(),
                            arc_v2.clone(),
                            [0., 1.],
                        );
                    }
                    an_arc.orientation = if arc_reverse {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    };
                    if ea1.orientation == Orientation::Reversed {
                        an_arc = bat::reversed(&an_arc);
                    }
                    ea2.orientation = bat::top_abs_reverse(ea1.orientation);
                    let mut arc_wire = bb.make_wire(&mut self.my_brep);
                    bat::builder_add_wire_edge(&mut arc_wire, &ea1);
                    bat::builder_add_wire_edge(&mut arc_wire, &an_arc);
                    bat::builder_add_wire_edge(&mut arc_wire, &ea2);
                    brep_lib_build_curves3d_tol(&arc_wire, self.my_tol);
                    bat::builder_set_closed(&mut arc_wire, true);
                    // OCCT L3659: BRepLib_MakeFace(arcWire, true) — the
                    // only-plane face maker (GAP leaf, arch. diff. #52).
                    let mut arc_face = super::brep_offset_make_offset::brep_lib_make_face_wire_only_plane(
                        &arc_wire,
                    );
                    super::brep_offset_make_offset::brep_tools_update(&mut arc_face);
                    self.my_walls.push(arc_face);
                    let cea2 = bat::oriented(&ea2, Orientation::Forward);
                    prev_edge = cea2;
                    prev_vertex = v2.clone();
                } else if is_build_from_scratch {
                    prev_edge = e4;
                    prev_vertex = v1;
                    is_build_from_scratch = false;
                } else {
                    prev_edge = e3;
                    prev_vertex = v2;
                }
                first_step = false;
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::MakeShells (cxx L3691-3792).
    pub(crate) fn make_shells(&mut self) {
        //
        // Prepare list of splits of the offset faces to make the shells
        let mut a_lsf: Vec<Shape> = Vec::new();
        let r = self.my_image_offset.roots().to_vec();
        //
        for it in &r {
            let mut a_f = it.clone();
            if self.my_thickening {
                // offsetted faces must change their orientations
                a_f = bat::reversed(&a_f);
            }
            //
            let mut image: Vec<Shape> = Vec::new();
            self.my_image_offset.last_image(&a_f, &mut image);
            for it2 in image {
                let a_f_im = it2;
                a_lsf.push(a_f_im);
            }
        }
        //
        if self.my_thickening {
            for explo in bat::explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
                let a_f = explo;
                a_lsf.push(a_f);
            }
            //
            for it in self.my_walls.clone() {
                let a_f = it;
                a_lsf.push(a_f);
            }
        }
        //
        if a_lsf.is_empty() {
            return;
        }
        //
        let mut b_done = false;
        if (self.my_join == GeomAbsJoinType::Intersection)
            && self.my_inter
            && !self.my_thickening
            && self.my_faces.is_empty()
            && super::brep_offset_make_offset_d::is_solid(&self.my_shape)
            && self.my_is_planar
        {
            //
            let mut a_shells = Shape::null();
            b_done = super::brep_offset_make_offset_d::build_shells_complete_inter(
                &a_lsf,
                &mut self.my_image_offset,
                &mut a_shells,
            );
            if b_done {
                self.my_offset_shape = a_shells;
            }
        }
        //
        if !b_done {
            let mut glue = BRepToolsQuilt::new();
            for a_it_ls in &a_lsf {
                glue.add(a_it_ls);
            }
            self.my_offset_shape = glue.shells();
        }
        //
        // Set correct value for closed flag
        for explo in bat::explorer(&self.my_offset_shape, ShapeType::Shell, ShapeType::Shape) {
            let mut a_s = explo;
            if !bat::shape_is_closed(&a_s) {
                if super::brep_offset_make_offset::brep_tool_is_closed(&a_s) {
                    bat::builder_set_closed(&mut a_s, true);
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::MakeSolid (cxx L3794-3894).
    pub(crate) fn make_solid(&mut self) {
        if self.my_offset_shape.is_null() {
            return;
        }
        //  Modified by skv - Mon Apr  4 18:17:27 2005 Begin
        //  Supporting history.
        let my_offset_shape = self.my_offset_shape.clone();
        super::brep_offset_make_offset_d::update_init_offset(
            &mut self.my_init_offset_face,
            &self.my_image_offset,
            &my_offset_shape,
            ShapeType::Face,
        );
        super::brep_offset_make_offset_d::update_init_offset(
            &mut self.my_init_offset_edge,
            &self.my_image_offset,
            &my_offset_shape,
            ShapeType::Edge,
        );
        //  Modified by skv - Mon Apr  4 18:17:27 2005 End
        let mut b = BRepBuilder::new();
        let mut nb_shell = 0;
        let mut nc = b.make_compound(&mut self.my_brep, vec![]);
        let mut nc_count = 0usize;
        let mut s1 = Shape::null();

        // OCCT L3824-3826: TopoDS_Solid Sol; B.MakeSolid(Sol);
        // Sol.Closed(true) — the rcad shells-Vec form (module-b arch.
        // note).
        let mut sol_shells: Vec<Shape> = Vec::new();
        let a_make_solid = self.my_shape.shape_type() == ShapeType::Solid || self.my_thickening;
        for exp in bat::explorer(&self.my_offset_shape, ShapeType::Shell, ShapeType::Shape) {
            if a_ps_user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
            let mut sh = exp;
            if self.my_thickening && self.my_offset > 0. {
                // OCCT L3839: Sh.Reverse().
                sh = bat::reversed(&sh);
            }
            nb_shell += 1;
            if bat::shape_is_closed(&sh) && a_make_solid {
                // OCCT L3843: B.Add(Sol, Sh).
                sol_shells.push(sh);
            } else {
                b.add_to_compound(&mut self.my_brep, nc.clone(), sh.clone());
                nc_count += 1;
                if nb_shell == 1 {
                    s1 = sh;
                }
            }
        }
        let sol = b.make_solid(&mut self.my_brep, sol_shells);
        let mut sol = sol;
        bat::builder_set_closed(&mut sol, true);
        let sol_children = bat::sub_shapes(&sol);
        let nbs = sol_children.len();
        let mut sol_is_null = nbs == 0;
        // Checking solid
        if nbs > 1 {
            // OCCT L3857: BRepCheck_Analyzer aCheck(Sol, false) (GAP
            // carrier, arch. diff. #49).
            let a_check = super::brep_offset_make_offset::BRepCheckAnalyzer::new(&sol, false);
            if !a_check.is_valid() {
                let mut a_sol_list: Vec<Shape> = Vec::new();
                super::brep_offset_make_offset_d::correct_solid(
                    &mut self.my_brep,
                    &mut sol,
                    &mut a_sol_list,
                );
                if !a_sol_list.is_empty() {
                    b.add_to_compound(&mut self.my_brep, nc.clone(), sol.clone());
                    nc_count += 1;
                    for a_sl_it in &a_sol_list {
                        b.add_to_compound(&mut self.my_brep, nc.clone(), a_sl_it.clone());
                        nc_count += 1;
                    }
                    sol_is_null = true;
                }
            }
        }
        let nc_is_null = nc_count == 0;
        if !sol_is_null && !nc_is_null {
            b.add_to_compound(&mut self.my_brep, nc.clone(), sol);
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

    /// OCCT BRepOffset_MakeOffset::SelectShells (cxx L3896-3936).
    pub(crate) fn select_shells(&mut self) {
        let mut free_edges: OcctShapeSet = HashMap::new();
        //-------------------------------------------------------------
        // FreeEdges all edges that can have free border in the
        // parallel shell
        // 1 - free borders of myShape .
        //-------------------------------------------------------------
        for exp in bat::explorer(&self.my_face_comp, ShapeType::Edge, ShapeType::Shape) {
            let e = exp;
            let la = self.my_analyse.ancestors(&e);
            if la.len() < 2 {
                // OCCT L3910: myAnalyse.Type(E).First().Type() — the OCCT
                // First() on the (Analyse-guaranteed non-empty) list.
                let t_e = self.my_analyse.type_(&e);
                let first_type = t_e.first().expect("Type(E).First() on empty list").my_type;
                if first_type == ChFiDS_TypeOfConcavity::FreeBound {
                    set_add(&mut free_edges, &e);
                }
            }
        }
        // myShape has free borders and there are no caps
        // no unwinding 3d.
        if !free_edges.is_empty() && self.my_faces.is_empty() {
            return;
        }

        self.my_offset_shape = deboucle3d(&self.my_offset_shape, &free_edges);
    }

    /// OCCT BRepOffset_MakeOffset::EncodeRegularity (cxx L3962-4165).
    pub(crate) fn encode_regularity(&mut self) {
        if self.my_offset_shape.is_null() {
            return;
        }
        // find edges G1 in the result
        let exp = bat::explorer(&self.my_offset_shape, ShapeType::Edge, ShapeType::Shape);

        let mut b = BRepBuilder::new();
        let mut ms: OcctShapeSet = HashMap::new();

        for oe in &exp {
            let mut oe = oe.clone();
            brep_lib_build_curve3d_edge(&oe, self.my_tol);
            let roe = oe.clone();

            if !set_add(&mut ms, &oe) {
                continue;
            }

            if self.my_image_offset.is_image(&oe) {
                oe = self.my_image_offset.root(&oe).clone();
            }

            let lof_of = self.my_as_des.ascendant(&roe).to_vec();

            if lof_of.len() != 2 {
                continue;
            }

            let f1 = lof_of[0].clone();
            let f2 = lof_of[1].clone();

            if f1.is_null() || f2.is_null() {
                continue;
            }

            let root1 = self.my_init_offset_face.root(&f1).clone();
            let root2 = self.my_init_offset_face.root(&f2).clone();

            let type1 = root1.shape_type();
            let type2 = root2.shape_type();

            if f1.is_same(&f2) {
                if bat::brep_tool_is_closed_on_surface(&oe, &f1) {
                    // Temporary Debug for the Bench.
                    // Check with YFR.
                    // In mode intersection, the edges are not coded in myInitOffsetEdge
                    // so, manage case by case
                    // Note DUB; for Hidden parts, it is NECESSARY to code CN
                    // Analytic Surfaces.
                    if self.my_join == GeomAbsJoinType::Intersection {
                        // OCCT L4002-4010: BRepAdaptor_Surface BS(F1, false);
                        // the analytic-surface probes (arch. diff. #54).
                        // OCCT L4003-4009: BRepAdaptor_Surface BS(F1, false); the
                        // analytic GetType probes (arch. diff. #54).
                        let s_type = bat::brep_tool_surface(&f1);
                        if matches!(
                            s_type.as_ref(),
                            Some(Surface3::Cylinder(_))
                                | Some(Surface3::Cone(_))
                                | Some(Surface3::Sphere(_))
                                | Some(Surface3::Torus(_))
                        ) {
                            b.continuity(&mut self.my_brep, &oe, &f1, &f1, GeomAbsShape::CN);
                        } else {
                            // See YFR : MaJ of myInitOffsetFace
                        }
                    } else if self.my_init_offset_edge.is_image(&roe) {
                        if type1 == ShapeType::Face && type2 == ShapeType::Face {
                            let ei = self.my_init_offset_edge.image_from(&roe).clone();
                            // OCCT L4021: GeomAbs_Shape Conti =
                            // BRep_Tool::Continuity(EI, FRoot, FRoot) — the
                            // draft_modification re-host.
                            let conti = crate::offset::draft_modification::brep_tool_continuity(
                                &ei, &root1, &root1,
                            );
                            if conti == GeomAbsShape::CN {
                                b.continuity(&mut self.my_brep, &oe, &f1, &f1, GeomAbsShape::CN);
                            } else if geom_abs_rank(&conti) > geom_abs_rank(&GeomAbsShape::C0) {
                                b.continuity(&mut self.my_brep, &oe, &f1, &f1, GeomAbsShape::G1);
                            }
                        }
                    }
                }
                continue;
            }

            //  code regularities G1 between :
            //    - sphere and tube : one root is a vertex, the other is an edge
            //                        and the vertex is included in the edge
            //    - face and tube   : one root is a face, the other an edge
            //                        and the edge is included in the face
            //    - face and face    : if two root faces are tangent in
            //                        the initial shape, they will be tangent in the offset shape
            //    - tube and tube  : if 2 edges generating tubes are
            //                        tangents, the 2 will be tangent either.
            if type1 == ShapeType::Edge && type2 == ShapeType::Vertex {
                let (v1, v2) = top_exp_vertices(&root1);
                if v1.is_same(&root2) || v2.is_same(&root2) {
                    b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                }
            } else if type1 == ShapeType::Vertex && type2 == ShapeType::Edge {
                let (v1, v2) = top_exp_vertices(&root2);
                if v1.is_same(&root1) || v2.is_same(&root1) {
                    b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                }
            } else if type1 == ShapeType::Face && type2 == ShapeType::Edge {
                for exp2 in bat::explorer(&root1, ShapeType::Edge, ShapeType::Shape) {
                    if exp2.is_same(&root2) {
                        b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                        break;
                    }
                }
            } else if type1 == ShapeType::Edge && type2 == ShapeType::Face {
                for exp2 in bat::explorer(&root2, ShapeType::Edge, ShapeType::Shape) {
                    if exp2.is_same(&root1) {
                        b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                        break;
                    }
                }
            } else if type1 == ShapeType::Face && type2 == ShapeType::Face {
                //  if two root faces are tangent in
                //  the initial shape, they will be tangent in the offset shape
                let mut le: Vec<Shape> = Vec::new();
                let mut lv: Vec<Shape> = Vec::new();
                find_common_shapes(&root1, &root2, &mut le, &mut lv);
                if le.len() == 1 {
                    let ed = le[0].clone();
                    if self.my_analyse.has_ancestor(&ed) {
                        let li = self.my_analyse.type_(&ed);
                        if li.len() == 1 && li[0].my_type == ChFiDS_TypeOfConcavity::Tangential {
                            b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                        }
                    }
                }
            } else if type1 == ShapeType::Edge && type2 == ShapeType::Edge {
                let mut lv: Vec<Shape> = Vec::new();
                let mut le: Vec<Shape> = Vec::new();
                find_common_shapes(&root1, &root2, &mut le, &mut lv);
                if lv.len() == 1 {
                    let mut led_tg: Vec<Shape> = Vec::new();
                    self.my_analyse
                        .tangent_edges(&root1, &lv[0], &mut led_tg);
                    for it in &led_tg {
                        if it.is_same(&root2) {
                            b.continuity(&mut self.my_brep, &oe, &f1, &f2, GeomAbsShape::G1);
                            break;
                        }
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::CheckInputData (cxx L4384-4510).
    pub fn check_input_data(&mut self) -> bool {
        // Set initial error state.
        self.my_error = BRepOffset_Error::NoError;
        self.my_bad_shape = Shape::null();
        // OCCT L4392: Message_ProgressScope aPS(theRange, nullptr, 1) — the
        // flattened rcad scope.
        let a_prog = NoopProgress;
        let a_ps = ProgressScope::new(&a_prog, "CheckInputData", 1);
        // Non-null offset.
        if self.my_offset.abs() <= self.my_tol {
            let mut is_found = false;
            for (_k, (_s, v)) in &self.my_face_offset {
                if v.abs() > self.my_tol {
                    is_found = true;
                    break;
                }
            }

            if !is_found {
                // No face with non-null offset found.
                self.my_error = BRepOffset_Error::NullOffset;
                return false;
            }
        }

        // Connectivity of input shape.
        if !super::brep_offset_make_offset::is_connected_shell(&self.my_face_comp) {
            self.my_error = BRepOffset_Error::NotConnectedShell;
            return false;
        }

        // Normals check and continuity check.
        let a_pnt_per_dim = 20; // 21 points on each dimension.
        let an_exp_sf = bat::explorer(&self.my_face_comp, ShapeType::Face, ShapeType::Shape);
        let mut a_presence_map: HashSet<u64> = HashSet::new();
        for a_f in &an_exp_sf {
            if a_ps.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return false;
            }

            if a_presence_map.contains(&a_f.ptr_id()) {
                // Not perform computations with partner shapes,
                // since they are contain same geometry.
                continue;
            }
            a_presence_map.insert(a_f.ptr_id());

            let a_surf = bat::brep_tool_surface(a_f);
            let (a_umin, a_umax, a_vmin, a_vmax) = brep_tools_uv_bounds(a_f);

            // Continuity check.
            // OCCT L4429-4433: aSurf->Continuity() == GeomAbs_C0 — the
            // Geom_Surface continuity probe (GAP leaf, arch. diff. #58).
            if geom_surface_continuity(a_surf.as_ref()) == GeomAbsShape::C0 {
                self.my_error = BRepOffset_Error::C0Geometry;
                return false;
            }

            // Get degenerated points, to avoid check them.
            let mut a_bad3d_pnts: Vec<glam::DVec3> = Vec::new();
            for an_exp_fe in bat::explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                let a_e = an_exp_fe;
                if brep_tool_degenerated(&a_e) {
                    let first_v = super::brep_offset_make_offset::top_exp_first_vertex(&a_e);
                    if let Some(p) = brep_tool_pnt(&first_v) {
                        a_bad3d_pnts.push(p);
                    }
                }
            }

            // Geometry grid check.
            for i in 0..=a_pnt_per_dim {
                let a_u_param = a_umin + (a_umax - a_umin) * i as f64 / a_pnt_per_dim as f64;
                for j in 0..=a_pnt_per_dim {
                    let a_v_param = a_vmin + (a_vmax - a_vmin) * j as f64 / a_pnt_per_dim as f64;

                    self.my_error = match &a_surf {
                        Some(s) => super::brep_offset_make_offset_d::check_single_point(
                            a_u_param, a_v_param, s, &a_bad3d_pnts,
                        ),
                        None => BRepOffset_Error::NoError,
                    };
                    if self.my_error != BRepOffset_Error::NoError {
                        return false;
                    }
                }
            }

            // Vertex list check.
            for an_exp_fv in bat::explorer(a_f, ShapeType::Vertex, ShapeType::Shape) {
                let a_v = an_exp_fv;
                // OCCT L4467: BRep_Tool::Parameters(aV, aF) — the vertex UV
                // parameters on the face (GAP leaf, arch. diff. #58).
                let a_pnt2d = brep_tool_parameters_vf(&a_v, a_f);

                self.my_error = match &a_surf {
                    Some(s) => super::brep_offset_make_offset_d::check_single_point(
                        a_pnt2d.x, a_pnt2d.y, s, &a_bad3d_pnts,
                    ),
                    None => BRepOffset_Error::NoError,
                };
                if self.my_error != BRepOffset_Error::NoError {
                    return false;
                }
            }
        }

        true
    }

    /// OCCT BRepOffset_MakeOffset::RemoveInternalEdges (cxx L4512-4560).
    pub(crate) fn remove_internal_edges(&mut self) {
        let mut a_dmelf = IndexedDataMapOfShapeListOfShape::new();
        //
        top_exp_map_shapes_and_ancestors(
            &self.my_offset_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_dmelf,
        );
        //
        let a_exp_f = bat::explorer(&self.my_offset_shape, ShapeType::Face, ShapeType::Shape);
        for a_f in &a_exp_f {
            let mut a_f = a_f.clone();
            //
            let mut a_liw: Vec<Shape> = Vec::new();
            //
            let a_exp_w = bat::explorer(&a_f, ShapeType::Wire, ShapeType::Shape);
            for a_w in &a_exp_w {
                let a_w = a_w.clone();
                //
                let mut b_remove_wire = true;
                let mut a_lie: Vec<Shape> = Vec::new();
                //
                let a_exp_e = bat::explorer(&a_w, ShapeType::Edge, ShapeType::Shape);
                for a_e in &a_exp_e {
                    let a_e = a_e.clone();
                    if a_e.orientation != Orientation::Internal {
                        b_remove_wire = false;
                        continue;
                    }
                    //
                    let a_lf = super::brep_offset_make_offset::shape_indexed_data_map_find(
                        &a_dmelf,
                        &a_e,
                    );
                    let b_remove_edge = a_lf.len() == 1;
                    if b_remove_edge {
                        a_lie.push(a_e);
                    } else {
                        b_remove_wire = false;
                    }
                }
                //
                if b_remove_wire {
                    a_liw.push(a_w);
                } else if !a_lie.is_empty() {
                    let mut a_w_mut = a_w.clone();
                    super::brep_offset_make_offset_d::remove_shapes(&mut a_w_mut, &a_lie);
                }
            }
            //
            if !a_liw.is_empty() {
                super::brep_offset_make_offset_d::remove_shapes(&mut a_f, &a_liw);
            }
        }
    }

    /// OCCT BRepOffset_MakeOffset::IntersectEdges (cxx L4674-4750).
    pub(crate) fn intersect_edges(
        &mut self,
        the_faces: &[Shape],
        the_map_sf: &mut MapSF,
        the_mes: &mut DataMapOfShapeShape,
        the_build: &DataMapOfShapeShape,
        the_as_des: &mut BRepAlgoAsDes,
        the_as_des2d: &mut BRepAlgoAsDes,
    ) {
        let mut a_dmvv: DmvvMap = indexmap::IndexMap::new();
        // intersect edges created from edges
        let mut a_mfv = IndexedShapeMap::new();
        // OCCT L4685: Message_ProgressScope aPSOuter(theRange, nullptr, 2) —
        // the flattened rcad scope.
        let a_prog = NoopProgress;
        let mut a_ps_outer = ProgressScope::new(&a_prog, "IntersectEdges", 2);
        for it in the_faces {
            let a_f = it.clone();
            let a_tol_f = brep_tool_tolerance(&a_f);
            // OCCT L4691-4703: BRepOffset_Inter2d::ConnexIntByInt(aF,
            // theMapSF(aF), ...) — the MES/Build DataMap-to-HashMap bridge
            // (architecture difference #59).
            let mut mes_view = datamap_view(the_mes);
            let build_view = datamap_view(the_build);
            // OCCT L4703: the Analyse argument — the Inter2d translation
            // carries its module-local Analyse GAP carrier (architecture
            // difference #60); theMapSF(aF) is the mutable entry form.
            let analyse_i2d = super::brep_offset_inter2d::BRepOffsetAnalyse;
            let ok = {
                let ofi = shape_data_map::change_find(the_map_sf, &a_f);
                BRepOffsetInter2d::connex_int_by_int(
                    &a_f,
                    ofi,
                &mut mes_view,
                &build_view,
                the_as_des,
                the_as_des2d,
                self.my_offset,
                a_tol_f,
                &analyse_i2d,
                &mut a_mfv,
                &mut self.my_image_vv,
                &mut self.my_edge_int_edges,
                &mut a_dmvv,
                (),
            )
            };
            datamap_merge_back(the_mes, mes_view);
            if !ok {
                self.my_error = BRepOffset_Error::CannotExtentEdge;
                return;
            }
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
        }
        // intersect edges created from vertices
        let a_nb_f = a_mfv.len();
        for i in 0..a_nb_f {
            // OCCT L4721: aF = TopoDS::Face(aMFV(i)).
            let a_f = a_mfv.get_index(i).unwrap().1.clone();
            let a_tol_f = brep_tool_tolerance(&a_f);
            let mut mes_view = datamap_view(the_mes);
            let build_view = datamap_view(the_build);
            // OCCT L4737: the Analyse argument — the Inter2d carrier (arch.
            // difference #60); theMapSF(aF) is the mutable entry form.
            let analyse_i2d = super::brep_offset_inter2d::BRepOffsetAnalyse;
            {
                let ofi = shape_data_map::change_find(the_map_sf, &a_f);
                BRepOffsetInter2d::connex_int_by_int_in_vert(
                    &a_f,
                    ofi,
                    &mut mes_view,
                    &build_view,
                    the_as_des,
                    the_as_des2d,
                    a_tol_f,
                    &analyse_i2d,
                    &mut a_dmvv,
                    (),
                );
            }
            datamap_merge_back(the_mes, mes_view);
            if a_ps_outer.user_break() {
                self.my_error = BRepOffset_Error::UserBreak;
                return;
            }
        }
        //
        // fuse vertices on edges
        if !BRepOffsetInter2d::fuse_vertices(&a_dmvv, the_as_des2d, &mut self.my_image_vv) {
            self.my_error = BRepOffset_Error::CannotFuseVertices;
            return;
        }
    }

    /// OCCT BRepOffset_MakeOffset::Generated (cxx L5231-5332).
    pub fn generated(&mut self, the_s: &Shape) -> &[Shape] {
        self.my_generated.clear();
        let a_type = the_s.shape_type();
        // OCCT L5240-5292: the switch with [[fallthrough]] — VERTEX falls
        // into EDGE falls into FACE; SOLID is a separate case.
        // case TopAbs_VERTEX:
        if a_type == ShapeType::Vertex {
            if self.my_analyse.has_ancestor(the_s) {
                let mut a_m_fence: OcctShapeSet = HashMap::new();
                let a_la = self.my_analyse.ancestors(the_s);
                let mut it_la = 0usize;
                while self.my_generated.is_empty() && it_la < a_la.len() {
                    let a_e = a_la[it_la].clone();
                    if !self.my_init_offset_edge.has_image(&a_e) {
                        it_la += 1;
                        continue;
                    }
                    let mut a_le_im: Vec<Shape> = Vec::new();
                    self.my_init_offset_edge.last_image(&a_e, &mut a_le_im);
                    let mut it_le_im = 0usize;
                    while self.my_generated.is_empty() && it_le_im < a_le_im.len() {
                        // OCCT L5264: TopoDS_Iterator itV(itLEIm.Value()).
                        for it_v in bat::sub_shapes(&a_le_im[it_le_im]) {
                            if !set_add(&mut a_m_fence, &it_v) {
                                self.my_generated.push(it_v.clone());
                                break;
                            }
                        }
                        it_le_im += 1;
                    }
                    it_la += 1;
                }
            }
        }
        // [[fallthrough]]
        // case TopAbs_EDGE:
        if matches!(a_type, ShapeType::Vertex | ShapeType::Edge) {
            if self.my_init_offset_edge.has_image(the_s) {
                let mut images: Vec<Shape> = Vec::new();
                self.my_init_offset_edge.last_image(the_s, &mut images);
                self.my_generated.extend(images);
            }
        }
        // [[fallthrough]]
        // case TopAbs_FACE:
        if matches!(
            a_type,
            ShapeType::Vertex | ShapeType::Edge | ShapeType::Face
        ) {
            let mut a_s = the_s.clone();
            if let Some(a_planface) = shape_data_map::seek(&self.my_face_planface_map, &a_s) {
                a_s = a_planface.clone();
            }

            if !self.my_faces.contains(&a_s) && self.my_init_offset_face.has_image(&a_s) {
                let mut images: Vec<Shape> = Vec::new();
                self.my_init_offset_face.last_image(&a_s, &mut images);
                self.my_generated.extend(images);

                if !self.my_faces.is_empty() {
                    // Reverse generated shapes in case of small solids.
                    // Useful only for faces without influence on others.
                    for it in self.my_generated.iter_mut() {
                        *it = bat::reversed(it);
                    }
                }
            }
        }
        // case TopAbs_SOLID:
        if a_type == ShapeType::Solid && the_s.is_same(&self.my_shape) {
            self.my_generated.push(self.my_offset_shape.clone());
        }

        if self.my_res_map.is_empty() {
            for s in bat::explorer(&self.my_offset_shape, ShapeType::Shape, ShapeType::Shape) {
                set_add(&mut self.my_res_map, &s);
            }
        }

        let res_map = &self.my_res_map;
        self.my_generated
            .retain(|s| res_map.contains_key(&shape_key(s)));

        &self.my_generated
    }

    /// OCCT BRepOffset_MakeOffset::Modified (cxx L5334-5367).
    pub fn modified(&mut self, the_shape: &Shape) -> &[Shape] {
        self.my_generated.clear();

        if the_shape.shape_type() == ShapeType::Face {
            let mut a_s = the_shape.clone();
            if let Some(a_planface) = shape_data_map::seek(&self.my_face_planface_map, &a_s) {
                a_s = a_planface.clone();
            }

            if self.my_faces.contains(&a_s) && self.my_init_offset_face.has_image(&a_s) {
                let mut images: Vec<Shape> = Vec::new();
                self.my_init_offset_face.last_image(&a_s, &mut images);
                self.my_generated.extend(images);

                if !self.my_faces.is_empty() {
                    // Reverse generated shapes in case of small solids.
                    // Useful only for faces without influence on others.
                    for it in self.my_generated.iter_mut() {
                        *it = bat::reversed(it);
                    }
                }
            }
        }

        &self.my_generated
    }

    /// OCCT BRepOffset_MakeOffset::IsDeleted (cxx L5369-5384).
    pub fn is_deleted(&mut self, the_s: &Shape) -> bool {
        if self.my_res_map.is_empty() {
            for s in bat::explorer(&self.my_offset_shape, ShapeType::Shape, ShapeType::Shape) {
                set_add(&mut self.my_res_map, &s);
            }
        }

        if self.my_res_map.contains_key(&shape_key(the_s)) {
            return false;
        }

        self.generated(the_s).is_empty() && self.modified(the_s).is_empty()
    }

    /// OCCT BRepOffset_MakeOffset::analyzeProgress (cxx L5410-5445).
    pub(crate) fn analyze_progress(&self, the_whole: f64, the_steps: &mut [f64]) {
        // OCCT L5414: theSteps.Init(0.0).
        for v in the_steps.iter_mut() {
            *v = 0.0;
        }

        // Set, approximately, the proportions for each operation.
        // It is not a problem that the sum of the set values will not
        // be equal to 100%, as the values will be normalized.
        // The main point is to make the proportions valid relatively each other.

        // Proportions will be different for different connection types
        let is_arc = self.my_join == GeomAbsJoinType::Arc;
        let is_planar_int_case = self.my_inter
            && !is_arc
            && self.my_is_planar
            && !self.my_thickening
            && self.my_faces.is_empty()
            && super::brep_offset_make_offset_d::is_solid(&self.my_shape);

        the_steps[BRepOffset_PIOperation::PIOperation_CheckInputData as usize] = 1.;
        the_steps[BRepOffset_PIOperation::PIOperation_Analyse as usize] = 2.;
        the_steps[BRepOffset_PIOperation::PIOperation_BuildOffsetBy as usize] =
            if is_planar_int_case {
                70.
            } else if is_arc {
                20.
            } else {
                50.
            };
        the_steps[BRepOffset_PIOperation::PIOperation_Intersection as usize] =
            if is_planar_int_case {
                0.
            } else if is_arc {
                50.
            } else {
                20.
            };
        if self.my_thickening {
            the_steps[BRepOffset_PIOperation::PIOperation_MakeMissingWalls as usize] = 5.;
        }
        the_steps[BRepOffset_PIOperation::PIOperation_MakeShells as usize] =
            if is_planar_int_case { 25. } else { 5. };
        the_steps[BRepOffset_PIOperation::PIOperation_MakeSolid as usize] = 5.;
        if self.my_is_perform_sewing && self.my_thickening {
            the_steps[BRepOffset_PIOperation::PIOperation_Sewing as usize] = 10.;
        }

        super::brep_offset_make_offset::normalize_steps(the_whole, the_steps);
    }

    /// OCCT BRepOffset_MakeOffset::IsPlanar (cxx L5447-5547).
    pub(crate) fn is_planar(&mut self) -> bool {
        let mut a_is_non_planar_found = false;
        let mut a_bb = BRepBuilder::new();

        let a_exp = bat::explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        for a_f in &a_exp {
            let a_f = a_f.clone();
            // OCCT L5458: BRepAdaptor_Surface aBAS(aF, false) — the plane
            // probe (arch. diff. #54: the rcad surface kind of the face
            // value).
            let a_bas_surface = bat::brep_tool_surface(&a_f);
            if matches!(a_bas_surface.as_ref(), Some(Surface3::Plane(_))) {
                continue;
            }

            if self.my_is_linearization_allowed {
                // define the toleance
                // (the OCCT face tolerance feeds BRepLib_MakeFace(F, P, Tol);
                // the rcad make_face form carries no tolerance argument).
                let a_tol_for_face = brep_tool_tolerance(&a_f);
                let _ = a_tol_for_face;

                // try to linearize
                // OCCT L5469-5470: GeomLib_IsPlanarSurface aPlanarityChecker
                // (GAP carrier, arch. diff. #52).
                let a_planarity_checker =
                    super::brep_offset_make_offset::GeomLibIsPlanarSurface::new(
                        a_bas_surface.as_ref().unwrap(),
                        PRECISION_CONFUSION,
                    );
                if a_planarity_checker.is_planar() {
                    let mut a_pln = a_planarity_checker.plan();
                    // OCCT L5474: aSurf->Bounds(u1, u2, v1, v2).
                    let bounds = a_bas_surface
                        .as_ref()
                        .unwrap()
                        .default_domain();
                    let (u1, u2, v1, v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                    let um;
                    let vm;
                    let is_inf1 = u1.is_infinite();
                    let is_inf2 = u2.is_infinite();
                    if !is_inf1 && !is_inf2 {
                        um = (u1 + u2) / 2.;
                    } else if is_inf1 && !is_inf2 {
                        um = u2 - 1.;
                    } else if !is_inf1 && is_inf2 {
                        um = u1 + 1.;
                    } else {
                        // isInf1 && isInf2
                        um = 0.;
                    }
                    let is_inf1 = v1.is_infinite();
                    let is_inf2 = v2.is_infinite();
                    if !is_inf1 && !is_inf2 {
                        vm = (v1 + v2) / 2.;
                    } else if is_inf1 && !is_inf2 {
                        vm = v2 - 1.;
                    } else if !is_inf1 && is_inf2 {
                        vm = v1 + 1.;
                    } else {
                        // isInf1 && isInf2
                        vm = 0.;
                    }
                    // OCCT L5495: aBAS.D1(um, vm, aP, aD1, aD2).
                    let (_a_p, a_d1, a_d2) = a_bas_surface
                        .as_ref()
                        .unwrap()
                        .derivatives(um, vm);
                    let a_norm = a_d1.cross(a_d2);
                    let mut a_pln_norm = a_pln.normal;
                    if a_norm.dot(a_pln_norm) < 0. {
                        // OCCT L5502-5504: aPlnNorm.Reverse(); gp_Ax1 anAx(...);
                        // aPln.SetAxis(anAx).
                        a_pln_norm = -a_pln_norm;
                        a_pln.normal = a_pln_norm;
                    }
                    // OCCT L5506-5511: aPlane = new Geom_Plane(aPln);
                    // aBB.MakeFace(aPlanarFace, aPlane, aTolForFace);
                    let a_plane = Surface3::Plane(a_pln);
                    let mut a_planar_face =
                        a_bb.make_face(&mut self.my_brep, Some(a_plane), Shape::null());
                    let a_face_forward = bat::oriented(&a_f, Orientation::Forward);
                    for an_it_face in bat::sub_shapes(&a_face_forward) {
                        let a_wire = an_it_face;
                        bat::builder_add_face_wire(&mut a_planar_face, &a_wire);
                    }
                    super::brep_offset_make_offset_d::remove_seam_and_degenerated_edges(
                        &a_planar_face,
                        &a_face_forward,
                    );
                    shape_data_map::bind(&mut self.my_face_planface_map, &a_f, a_planar_face.clone());
                    if self.my_faces.contains(&a_f) {
                        indexed_shape_map_remove_key(&mut self.my_faces, &a_f);
                        self.my_faces.add(&a_planar_face);
                    }
                } else {
                    a_is_non_planar_found = true;
                }
            } else {
                a_is_non_planar_found = true;
            }
        }

        !a_is_non_planar_found
    }
}

// ===========================================================================
// Local helpers of module e.
// ===========================================================================

/// OCCT gp::Resolution() (2.2250738585072014e-308) — the smallest positive
/// double (the offset_wire.rs precedent).
const GP_RESOLUTION: f64 = 2.2250738585072014e-308;

/// The Geom2d_Line constructor re-host (OCCT new Geom2d_Line(P, D)).
fn line2d_of(ox: f64, oy: f64, dx: f64, dy: f64) -> Curve2d {
    Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: DVec2::new(ox, oy),
        direction: DVec2::new(dx, dy),
    })
}

/// The rcad Plane U-direction stand-in for the Geom_Plane(P, N) Ax3 frame.
fn plane_u_dir_of(normal: glam::DVec3) -> glam::DVec3 {
    let a = if normal.x.abs() < 0.9 {
        glam::DVec3::X
    } else {
        glam::DVec3::Y
    };
    a.cross(normal).normalize_or_zero()
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) — the pcurve-on-surface
/// (face-less) form — GAP leaf (architecture difference #58: the rcad
/// pcurve storage is face-keyed; the surface-keyed storage of the
/// MakeMissingWalls forms is not translated).
fn bb_update_edge_pcurve_on_surface(the_e: &mut Shape, _the_c2d: &Curve2d, _the_s: &Surface3, _the_tol: f64) {
    let _ = the_e;
    panic!("GAP: BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) (surface-keyed pcurve storage)");
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) — the closed-surface
/// two-pcurve form — GAP leaf (architecture difference #58).
fn bb_update_edge_pcurves_two_on_surface(
    the_e: &mut Shape,
    _the_c1: &Curve2d,
    _the_c2: &Curve2d,
    _the_s: &Surface3,
    _the_tol: f64,
) {
    let _ = the_e;
    panic!("GAP: BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) (surface-keyed pcurve storage)");
}

/// OCCT BRep_Builder::Range(E, S, L, First, Last) — the pcurve range on the
/// surface — GAP leaf (architecture difference #58).
fn bb_range_edge_on_surface(the_e: &mut Shape, _the_first: f64, _the_last: f64) {
    let _ = the_e;
    panic!("GAP: BRep_Builder::Range(E, S, L, First, Last) (surface-keyed pcurve storage)");
}

/// OCCT BRep_Builder::UpdateEdge(E, C3d, Tol) — the 3D curve update of the
/// module-e forms (the TShape::Edge make_mut edit, arch. diff. #22).
fn update_edge_3d_host(the_e: &mut Shape, the_c: &Curve3, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.curve = Some(the_c.clone());
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol) — the tolerance-only form (the
/// wire-edge walk of MakeMissingWalls).
fn update_edge_tolerance_host_e(the_e: &mut Shape, the_tol: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT Geom_Surface::UIso(U) — GAP leaf (architecture difference #58: the
/// surface iso-construction is not translated).
fn surface_u_iso(_the_s: &Surface3, _the_u: f64) -> Curve3 {
    panic!("GAP: Geom_Surface::UIso (TKMath/Geom not translated)");
}

/// OCCT BRep_Tool::Parameters(V, F) — the vertex UV parameters on the face
/// — GAP leaf (architecture difference #58).
fn brep_tool_parameters_vf(_the_v: &Shape, _the_f: &Shape) -> DVec2 {
    panic!("GAP: BRep_Tool::Parameters(V, F) (TKTopAlgo/BRep not translated)");
}

/// OCCT Geom_Surface::Continuity() — GAP leaf (architecture difference #58:
/// the rcad Surface3 carries no continuity field; the OCCT default of the
/// elementary surfaces is GeomAbs_C2, so the C0 probe returns C2 here).
fn geom_surface_continuity(_the_s: Option<&Surface3>) -> GeomAbsShape {
    GeomAbsShape::C2
}

/// OCCT BRep_Tool::IsClosed(S) — the brep_algo re-host alias (the
/// shape_is_closed face-pair walk of the OCCT BRep_Tool::IsClosed(shell)).
pub(crate) fn brep_tool_is_closed_host(the_s: &Shape) -> bool {
    bat::shape_is_closed(the_s)
}

/// OCCT GeomAbs_Shape ranking (C0 < G1 < C1 < G2 < C2 < C3 < CN) — the
/// `Conti > GeomAbs_C0` comparison form.
fn geom_abs_rank(c: &GeomAbsShape) -> u8 {
    match c {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::G1 => 1,
        GeomAbsShape::C1 => 2,
        GeomAbsShape::G2 => 3,
        GeomAbsShape::C2 => 4,
        GeomAbsShape::C3 => 5,
        GeomAbsShape::CN => 6,
    }
}

/// OCCT TopExp::FirstVertex(E) — the module-e alias of the cum-ori form.
pub(crate) fn top_exp_first_vertex_host(edg: &Shape) -> Shape {
    super::brep_offset_make_offset::top_exp_first_vertex(edg)
}

/// OCCT TopExp::Vertices(E, V1, V2) — the module-e alias.
pub(crate) fn top_exp_vertices_host(edg: &Shape) -> (Shape, Shape) {
    top_exp_vertices(edg)
}

/// OCCT TopExp::Vertices(E, V1, V2) — the raw pair (module-e alias for the
/// tool re-host).
pub(crate) fn top_exp_vertices_of_host(edg: &Shape) -> (Shape, Shape) {
    top_exp_vertices(edg)
}

/// The OCCT `!aPS.More()` user-break probe — the flattened NoopProgress
/// scope never breaks (architecture difference #40).
fn a_ps_user_break() -> bool {
    false
}

/// OCCT DataMap<Shape, Shape> -> HashMap<ShapeKey, Shape> view (the
/// BRepOffset_Inter2d ConnexIntByInt forms; architecture difference #59).
fn datamap_view(m: &ShapeDataMap<Shape>) -> HashMap<crate::feat::loc_ope_wires_on_shape_b::ShapeKey, Shape> {
    m.iter().map(|(k, (_, v))| (*k, v.clone())).collect()
}

/// The write-back of the mutated HashMap view into the canonical DataMap
/// (architecture difference #59).
fn datamap_merge_back(
    m: &mut ShapeDataMap<Shape>,
    view: HashMap<crate::feat::loc_ope_wires_on_shape_b::ShapeKey, Shape>,
) {
    for (k, v) in view {
        let entry = m.entry(k).or_insert((Shape::null(), v.clone()));
        entry.1 = v;
    }
}
