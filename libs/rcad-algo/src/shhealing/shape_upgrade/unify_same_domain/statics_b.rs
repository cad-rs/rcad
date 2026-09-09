//! OCCT ShapeUpgrade_UnifySameDomain.cxx L536-1750 — the second file-statics
//! segment: `RelocatePCurvesToNewUorigin` (L536-697), `InsertWiresIntoFaces`
//! (L699-733), `FindCommonFace` (L735-776), `FindClosestPoints` (L778-825),
//! `ReconstructMissedSeam` (L829-918), `SameSurf` (L922-1053),
//! `TransformPCurves` (L1057-1266), `AddPCurves` (L1270-1294),
//! `AddOrdinaryEdges` (L1301-1352), `getCylinder` (L1356-1431),
//! `ClearRts` (L1435-1440), `GetNormalToSurface` (L1447-1514),
//! `IsSameDomain` (L1518-1642), `UpdateMapOfShapes` (L1646-1661),
//! `GlueEdgesWith3DCurves` (L1667-1750).

use std::collections::HashMap;
use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval, TrimmedCurve3};
use rcad_kernel::precision::{p_confusion, CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType, State, TShape};

use super::gap_deps::{
    dir_angle_3d, dir_angle_with_ref_3d, dir_is_parallel_3d, elclib_line_parameter_3d,
    geom2d_copy_untrim, geom2d_mirror_ox2d, geom2d_mirror_oy2d, geom2d_translate, geom2d_trimmed,
    geom_convert_c0_to_c1, geom_convert_comp_curve_add, geom_convert_concat_c1,
    geom_convert_curve_to_bspline, geom2d_convert_c0_to_c1, geom2d_convert_comp_curve_add,
    geom2d_convert_concat_c1, geom2d_convert_curve_to_bspline, precision_is_neg_inf,
    precision_is_pos_inf, pcurve_handle_same, surface_handle_same, BRepAdaptorCurve2d,
};
use super::statics_a::{
    get_curve_params, remove_edge_from_map, true_value_of_offset, REAL_FIRST, REAL_LAST,
};
use super::topexp::{
    brep_tool_curve, brep_tool_is_closed_edge_face, brep_tool_range, brep_tool_surface_loc,
    common_vertex, first_vertex, is_edge_degenerated, last_vertex, map_shapes_all,
    occt_is_same_shape, topexp_explorer, vertices,
};
use super::{map_add, shape_key, DataMapOfFacePlane, MapOfShape};
use crate::bop::int_tools::int_tools_fclass2d::IntToolsFClass2d;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes};
use crate::shhealing::shape_build::edge::builder_update_edge_pcurve_null;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

/// OCCT gp_Cylinder re-host (the `(location, direction, radius)` triple).
#[derive(Clone, Debug)]
pub struct GpCylinder {
    /// OCCT gp_Cylinder::Location().
    pub location: DVec3,
    /// OCCT gp_Cylinder::Position().Direction().
    pub direction: DVec3,
    /// OCCT gp_Cylinder::Radius().
    pub radius: f64,
}

/// OCCT `NCollection_DataMap<TopoDS_Shape, Handle(Geom2d_Curve), ...>` (the
/// EdgeNewPCurve local of RelocatePCurvesToNewUorigin).
pub type DataMapOfShapePCurve = HashMap<(u64, u32), (Shape, Curve2d)>;

/// The first/last parameter of the rcad 2D curve value over the stored knot
/// domain (the Geom2d_Curve FirstParameter/LastParameter reads).
pub fn first_param_of_2d(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Line(_) | Curve2d::Circle(_) => c.default_domain()[0],
        Curve2d::BSpline(b) => b.knots[b.degree],
        Curve2d::Bezier(_) => 0.0,
        Curve2d::Trimmed(t) => t.t_min,
        _ => c.default_domain()[0],
    }
}

/// See [`first_param_of_2d`].
pub fn last_param_of_2d(c: &Curve2d) -> f64 {
    match c {
        Curve2d::Line(_) | Curve2d::Circle(_) => c.default_domain()[1],
        Curve2d::BSpline(b) => b.knots[b.knots.len() - 1 - b.degree],
        Curve2d::Bezier(_) => 1.0,
        Curve2d::Trimmed(t) => t.t_max,
        _ => c.default_domain()[1],
    }
}

/// The coordinate selector of RelocatePCurvesToNewUorigin.
fn coord_of(ind: usize, p: DVec2) -> f64 {
    if ind == 1 {
        p.x
    } else {
        p.y
    }
}

fn coord_pair(ind: usize, a: DVec2, b: DVec2) -> (f64, f64) {
    (coord_of(ind, a), coord_of(ind, b))
}

/// OCCT static RelocatePCurvesToNewUorigin (cxx L536-697).
#[allow(clippy::too_many_arguments)]
pub fn relocate_pcurves_to_new_uorigin(
    brep: &mut BRep,
    the_edges: &[Shape],
    the_first_face: &Shape,
    the_ref_face: &Shape,
    the_coord_tol: f64,
    the_ind_coord: usize,
    the_period: f64,
    the_vemap: &mut super::IndexedDataMapOfShapeListOfShape,
    the_edge_new_pcurve: &mut DataMapOfShapePCurve,
    the_used_edges: &mut MapOfShape,
) {
    // OCCT L548-549.
    let mut an_edges_of_first_face = MapOfShape::new();
    map_shapes_all(brep, the_first_face, &mut an_edges_of_first_face);

    // OCCT L551: for (;;) // walk by contours.
    loop {
        // OCCT L553-603: the start edge with the minimum coordinate.
        let mut a_start_edge = Shape::null();
        let mut an_orientation = Orientation::Forward;
        let mut a_coord_min = REAL_LAST;
        for an_edge_index in 1..=the_edges.len() {
            let an_edge = the_edges[an_edge_index - 1].clone();
            if the_used_edges.contains_key(&shape_key(&an_edge)) {
                continue;
            }
            if !an_edges_of_first_face.contains_key(&shape_key(&an_edge)) {
                continue;
            }

            if a_start_edge.is_null() {
                a_start_edge = an_edge.clone();
                let (a_first_point, a_last_point) =
                    get_curve_params(brep, &a_start_edge, the_ref_face);
                let (fc, lc) = coord_pair(the_ind_coord, a_first_point, a_last_point);
                if fc < lc {
                    a_coord_min = fc;
                    an_orientation = Orientation::Forward;
                } else {
                    a_coord_min = lc;
                    an_orientation = Orientation::Reversed;
                }
            } else {
                let (a_first_point, a_last_point) = get_curve_params(brep, &an_edge, the_ref_face);
                let (fc, lc) = coord_pair(the_ind_coord, a_first_point, a_last_point);
                if fc < a_coord_min {
                    a_start_edge = an_edge.clone();
                    a_coord_min = fc;
                    an_orientation = Orientation::Forward;
                }
                if lc < a_coord_min {
                    a_start_edge = an_edge;
                    a_coord_min = lc;
                    an_orientation = Orientation::Reversed;
                }
            }
        }

        // OCCT L605-608: all contours are passed.
        if a_start_edge.is_null() {
            break;
        }

        // OCCT L610-621.
        let mut a_current_edge = a_start_edge.clone();
        let cur_pc = brep.curve_on_surface(&a_current_edge, the_ref_face);
        let cur_trimmed: Option<Curve2d> =
            cur_pc.as_ref().map(|(pc, f, l)| geom2d_trimmed(pc, *f, *l));
        if let Some(pc) = cur_trimmed.as_ref() {
            the_edge_new_pcurve
                .insert(shape_key(&a_current_edge), (a_current_edge.clone(), pc.clone()));
        }
        let (mut an_edge_start_param, mut an_edge_end_param) = match cur_pc.as_ref() {
            Some((_, f, l)) => (*f, *l),
            None => (0.0, 0.0),
        };
        if a_current_edge.orientation == Orientation::Reversed {
            std::mem::swap(&mut an_edge_start_param, &mut an_edge_end_param);
        }
        let mut cur_param = if an_orientation == Orientation::Forward {
            an_edge_end_param
        } else {
            an_edge_start_param
        };
        let mut cur_point = match cur_trimmed.as_ref() {
            Some(pc) => pc.point_at(cur_param),
            None => DVec2::ZERO,
        };

        // OCCT L623: for (;;) // collect pcurves of a contour.
        loop {
            // OCCT L625-628.
            if !remove_edge_from_map(brep, &a_current_edge, the_vemap) {
                break; // end of contour in 2d
            }
            map_add(the_used_edges, &a_current_edge);
            // OCCT L630-632.
            let cur_vertex = if an_orientation == Orientation::Forward {
                last_vertex(brep, &a_current_edge, true)
            } else {
                first_vertex(brep, &a_current_edge, true)
            };

            // OCCT L634-638.
            let elist = match the_vemap.get(&shape_key(&cur_vertex)) {
                Some((_, l)) if !l.is_empty() => l.clone(),
                _ => break, // end of contour in 3d
            };

            // OCCT L640-694.
            for an_edge in &elist {
                if occt_is_same_shape(an_edge, &a_current_edge) {
                    continue;
                }

                // OCCT L648-654.
                let a_first_vertex = if an_orientation == Orientation::Forward {
                    first_vertex(brep, an_edge, true)
                } else {
                    last_vertex(brep, an_edge, true)
                };
                if !occt_is_same_shape(&a_first_vertex, &cur_vertex) {
                    continue; // may be if CurVertex is deg.vertex
                }

                // OCCT L656-662.
                let Some((a_pc0, f, l)) = brep.curve_on_surface(an_edge, the_ref_face) else {
                    continue;
                };
                let mut a_pcurve = geom2d_trimmed(&a_pc0, f, l);
                let (mut sp, mut ep) = (f, l);
                if an_edge.orientation == Orientation::Reversed {
                    std::mem::swap(&mut sp, &mut ep);
                }
                // OCCT L663-670.
                let a_param = if an_orientation == Orientation::Forward {
                    sp
                } else {
                    ep
                };
                let a_point = a_pcurve.point_at(a_param);
                let an_offset =
                    coord_of(the_ind_coord, cur_point) - coord_of(the_ind_coord, a_point);
                if an_offset.abs() >= the_coord_tol
                    && (an_offset.abs() - the_period).abs() >= the_coord_tol
                {
                    continue; // may be if CurVertex is deg.vertex
                }

                // OCCT L672-687.
                if an_offset.abs() > the_period / 2.0 {
                    let an_offset = true_value_of_offset(an_offset, the_period);
                    let a_vec = if the_ind_coord == 1 {
                        DVec2::new(an_offset, 0.0)
                    } else {
                        DVec2::new(0.0, an_offset)
                    };
                    a_pcurve = geom2d_translate(&a_pcurve, a_vec);
                }
                // OCCT L688-693.
                the_edge_new_pcurve.insert(shape_key(an_edge), (an_edge.clone(), a_pcurve.clone()));
                a_current_edge = an_edge.clone();
                let cur_or = an_orientation.compose(a_current_edge.orientation);
                cur_param = if cur_or == Orientation::Forward {
                    last_param_of_2d(&a_pcurve)
                } else {
                    first_param_of_2d(&a_pcurve)
                };
                cur_point = a_pcurve.point_at(cur_param);
                break;
            }
        } // for (;;) (collect pcurves of a contour)
    } // for (;;) (walk by contours)
}

/// OCCT static InsertWiresIntoFaces (cxx L699-733).
pub fn insert_wires_into_faces(
    brep: &mut BRep,
    the_wires: &[Shape],
    the_faces: &[Shape],
    the_ref_face: &Shape,
) {
    for ii in 1..=the_wires.len() {
        let a_wire = the_wires[ii - 1].clone();
        // OCCT L707-711: the first edge of the wire and its pcurve midpoint.
        let it = iter_subshapes(brep, &a_wire, true, true);
        let Some(an_edge) = it.first() else { continue };
        let an_edge = an_edge.clone();
        let ba_curve2d = BRepAdaptorCurve2d::new(brep, &an_edge, the_ref_face);
        let a_pnt2d =
            ba_curve2d.value((ba_curve2d.first_parameter() + ba_curve2d.last_parameter()) / 2.0);
        // OCCT L712-723.
        let mut required_face = Shape::null();
        for jj in 1..=the_faces.len() {
            let a_face = the_faces[jj - 1].clone();
            let classifier =
                IntToolsFClass2d::new_face(Arc::new(brep.clone()), &a_face, CONFUSION);
            let a_status = classifier.perform(a_pnt2d, true);
            if a_status == State::In {
                required_face = a_face;
                required_face.orientation = Orientation::Forward;
                break;
            }
        }
        if !required_face.is_null() {
            // OCCT L726: BB.Add(RequiredFace, aWire).
            builder_add(brep, &required_face, &a_wire);
        } else {
            // OCCT L730: Standard_ASSERT_INVOKE — the assert failure raise.
            panic!("ShapeUpgrade_UnifySameDomain: wire remains unclassified");
        }
    }
}

/// OCCT static FindCommonFace (cxx L735-776).
pub fn find_common_face(
    brep: &mut BRep,
    the_edge1: &Shape,
    the_edge2: &Shape,
    the_vfmap: &super::IndexedDataMapOfShapeListOfShape,
    the_or_of_e1_on_face: &mut Orientation,
    the_or_of_e2_on_face: &mut Orientation,
) -> Shape {
    // OCCT L744-746.
    let a_vertex = common_vertex(brep, the_edge1, the_edge2).unwrap_or_else(Shape::null);
    let flist = the_vfmap
        .get(&shape_key(&a_vertex))
        .map(|(_, l)| l.clone())
        .unwrap_or_default();
    for itl in &flist {
        let mut a_face = itl.clone();
        a_face.orientation = Orientation::Forward;
        let mut e1found = false;
        let mut e2found = false;
        for an_edge in topexp_explorer(brep, &a_face, ShapeType::Edge) {
            if occt_is_same_shape(&an_edge, the_edge1) {
                e1found = true;
                *the_or_of_e1_on_face = an_edge.orientation;
            }
            if occt_is_same_shape(&an_edge, the_edge2) {
                e2found = true;
                *the_or_of_e2_on_face = an_edge.orientation;
            }
            if e1found && e2found {
                return a_face;
            }
        }
    }

    // OCCT L774-775: the null face.
    Shape::null()
}

/// OCCT static FindClosestPoints (cxx L778-825).
#[allow(clippy::too_many_arguments)]
pub fn find_closest_points(
    brep: &mut BRep,
    the_edge1: &Shape,
    the_edge2: &Shape,
    the_vfmap: &super::IndexedDataMapOfShapeListOfShape,
    the_common_face: &mut Shape,
    the_min_sq_dist: &mut f64,
    or_of_e1_on_face: &mut Orientation,
    or_of_e2_on_face: &mut Orientation,
    the_ind_on_e1: &mut i32,
    the_ind_on_e2: &mut i32,
    points_on_edge1: &mut [DVec2; 2],
    points_on_edge2: &mut [DVec2; 2],
) -> bool {
    // OCCT L793-797.
    *the_common_face = find_common_face(
        brep,
        the_edge1,
        the_edge2,
        the_vfmap,
        or_of_e1_on_face,
        or_of_e2_on_face,
    );
    if the_common_face.is_null() {
        return false;
    }

    // OCCT L799-807.
    let pcurve1 = brep.curve_on_surface(the_edge1, the_common_face);
    let pcurve2 = brep.curve_on_surface(the_edge2, the_common_face);
    let (Some((pc1, f1, l1)), Some((pc2, f2, l2))) = (pcurve1, pcurve2) else {
        return false;
    };
    points_on_edge1[0] = pc1.point_at(f1);
    points_on_edge1[1] = pc1.point_at(l1);
    points_on_edge2[0] = pc2.point_at(f2);
    points_on_edge2[1] = pc2.point_at(l2);
    // OCCT L808-823.
    *the_min_sq_dist = REAL_LAST;
    *the_ind_on_e1 = -1;
    *the_ind_on_e2 = -1;
    for ind1 in 0..2i32 {
        for ind2 in 0..2i32 {
            let a_sq_dist =
                points_on_edge1[ind1 as usize].distance_squared(points_on_edge2[ind2 as usize]);
            if a_sq_dist < *the_min_sq_dist {
                *the_min_sq_dist = a_sq_dist;
                *the_ind_on_e1 = ind1;
                *the_ind_on_e2 = ind2;
            }
        }
    }
    true
}

/// OCCT static ReconstructMissedSeam (cxx L829-918).
#[allow(clippy::too_many_arguments)]
pub fn reconstruct_missed_seam(
    brep: &mut BRep,
    the_removed_edges: &[Shape],
    the_fref_face: &Shape,
    the_cur_edge: &Shape,
    the_cur_vertex: &Shape,
    the_cur_point: DVec2,
    the_uperiod: f64,
    the_vperiod: f64,
    the_seam_edge: &mut Shape,
    the_next_point: &mut DVec2,
) {
    // OCCT L839.
    let _ref_surf = brep.face_surface_world(the_fref_face);

    // OCCT L841-842: find the seam edge among the removed edges.
    *the_seam_edge = Shape::null();
    for i in 1..=the_removed_edges.len() {
        let mut an_edge = the_removed_edges[i - 1].clone();
        if occt_is_same_shape(&an_edge, the_cur_edge) {
            continue;
        }
        let Some((mut a_pc, mut param1, mut param2)) =
            brep.curve_on_surface(&an_edge, the_fref_face)
        else {
            continue;
        };

        // OCCT L857-858.
        let (mut a_first_vertex, a_last_vertex) = vertices(brep, &an_edge, true);

        if (occt_is_same_shape(&a_first_vertex, the_cur_vertex)
            || occt_is_same_shape(&a_last_vertex, the_cur_vertex))
            && brep_tool_is_closed_edge_face(brep, &an_edge, the_fref_face)
        {
            // OCCT L863-866.
            let mut a_param = if an_edge.orientation == Orientation::Forward {
                param1
            } else {
                param2
            };
            let mut a_point = a_pc.point_at(a_param);
            let mut a_udiff = (a_point.x - the_cur_point.x).abs();
            let mut a_vdiff = (a_point.y - the_cur_point.y).abs();
            if (the_uperiod != 0.0 && a_udiff > the_uperiod / 2.0)
                || (the_vperiod != 0.0 && a_vdiff > the_vperiod / 2.0)
            {
                // OCCT L870-890.
                if occt_is_same_shape(&a_last_vertex, the_cur_vertex)
                    || (the_uperiod != 0.0 && the_vperiod != 0.0)
                {
                    an_edge.orientation = super::occt_reverse(an_edge.orientation);
                } else {
                    let an_ori = an_edge.orientation;
                    an_edge.orientation = Orientation::Forward;
                    let pc1 = brep
                        .curve_on_surface(&an_edge, the_fref_face)
                        .map(|(c, _, _)| c);
                    an_edge.orientation = super::occt_reverse(an_edge.orientation);
                    let pc2 = brep
                        .curve_on_surface(&an_edge, the_fref_face)
                        .map(|(c, _, _)| c);
                    an_edge.orientation = Orientation::Forward; // again FORWARD
                    let a_tol = brep.tolerance(&an_edge);
                    let (a_surf, a_loc) = brep_tool_surface_loc(brep, the_fref_face);
                    if let (Some(pc1), Some(pc2), Some(a_surf)) = (pc1, pc2, a_surf) {
                        // OCCT L888: aBB.UpdateEdge(anEdge, aPC2, aPC1, aSurf,
                        // aLoc, aTol) — the closed-surface pcurve pair (the
                        // kernel form carries the edge range).
                        let rng = brep_tool_range(brep, &an_edge);
                        let mut builder = BRepBuilder::new();
                        builder.update_edge_pcurve_closed(
                            brep,
                            an_edge.clone(),
                            pc2.clone(),
                            pc1.clone(),
                            the_fref_face.clone(),
                            rng[0],
                            rng[1],
                            a_tol,
                        );
                    }
                    an_edge.orientation = an_ori;
                    let _ = a_loc;
                }
                // OCCT L891-895.
                if let Some((pc, p1, p2)) = brep.curve_on_surface(&an_edge, the_fref_face) {
                    a_pc = pc;
                    param1 = p1;
                    param2 = p2;
                }
                a_param = if an_edge.orientation == Orientation::Forward {
                    param1
                } else {
                    param2
                };
                a_point = a_pc.point_at(a_param);
                a_udiff = (a_point.x - the_cur_point.x).abs();
                a_vdiff = (a_point.y - the_cur_point.y).abs();
            }
            // OCCT L897-906.
            if (the_uperiod == 0.0 || a_udiff < the_uperiod / 2.0)
                && (the_vperiod == 0.0 || a_vdiff < the_vperiod / 2.0)
            {
                a_first_vertex = first_vertex(brep, &an_edge, true);
                if occt_is_same_shape(&a_first_vertex, the_cur_vertex) {
                    *the_seam_edge = an_edge;
                    break;
                }
            }
        }
    }

    // OCCT L910-917.
    if !the_seam_edge.is_null() {
        if let Some((a_pc, param1, param2)) = brep.curve_on_surface(the_seam_edge, the_fref_face) {
            let a_param = if the_seam_edge.orientation == Orientation::Forward {
                param2
            } else {
                param1
            };
            *the_next_point = a_pc.point_at(a_param);
        }
    }
}

/// OCCT static SameSurf (cxx L922-1053): the sampled bounds/values surface
/// coincidence test (the handle-equality fast path lives in the caller).
pub fn same_surf(the_s1: &Surface3, the_s2: &Surface3) -> bool {
    // OCCT L924.
    const A_COEFS: [f64; 2] = [0.3399811, 0.7745966];

    // OCCT L926-928: the bounds.
    let (mut uf1, mut ul1, mut vf1, mut vl1) = {
        let d = the_s1.default_domain();
        (d[0], d[1], d[2], d[3])
    };
    let (uf2, ul2, vf2, vl2) = {
        let d = the_s2.default_domain();
        (d[0], d[1], d[2], d[3])
    };
    // OCCT L929.
    let a_p_tol = p_confusion();
    // OCCT L930-954 (u-first arm).
    if precision_is_neg_inf(uf1) {
        if !precision_is_neg_inf(uf2) {
            return false;
        } else {
            uf1 = (-1.0f64).min(ul1 - 1.0);
        }
    } else {
        if precision_is_neg_inf(uf2) {
            return false;
        }
        if (uf1 - uf2).abs() > a_p_tol {
            return false;
        }
    }
    // OCCT L956-980 (v-first arm).
    if precision_is_neg_inf(vf1) {
        if !precision_is_neg_inf(vf2) {
            return false;
        } else {
            vf1 = (-1.0f64).min(vl1 - 1.0);
        }
    } else {
        if precision_is_neg_inf(vf2) {
            return false;
        }
        if (vf1 - vf2).abs() > a_p_tol {
            return false;
        }
    }
    // OCCT L982-1006 (u-last arm).
    if precision_is_pos_inf(ul1) {
        if !precision_is_pos_inf(ul2) {
            return false;
        } else {
            ul1 = 1.0f64.max(uf1 + 1.0);
        }
    } else {
        if precision_is_pos_inf(ul2) {
            return false;
        }
        if (ul1 - ul2).abs() > a_p_tol {
            return false;
        }
    }
    // OCCT L1008-1032 (v-last arm).
    if precision_is_pos_inf(vl1) {
        if !precision_is_pos_inf(vl2) {
            return false;
        } else {
            vl1 = 1.0f64.max(vf1 + 1.0);
        }
    } else {
        if precision_is_pos_inf(vl2) {
            return false;
        }
        if (vl1 - vl2).abs() > a_p_tol {
            return false;
        }
    }

    // OCCT L1035-1050: the sampled surface values.
    let du = ul1 - uf1;
    let dv = vl1 - vf1;
    for i in 0..2 {
        let u = uf1 + A_COEFS[i] * du;
        for j in 0..2 {
            let v = vf1 + A_COEFS[j] * dv;
            let a_p1 = the_s1.point_at(u, v);
            let a_p2 = the_s2.point_at(u, v);
            if a_p1.distance(a_p2) > a_p_tol {
                return false;
            }
        }
    }

    true
}

/// OCCT static TransformPCurves (cxx L1057-1266).
pub fn transform_pcurves(
    brep: &mut BRep,
    the_ref_face: &Shape,
    the_face: &Shape,
    the_map_edges_with_temporary_pcurves: &mut MapOfShape,
) {
    // OCCT L1062-1066: the reference surface (one trimmed level unwrapped).
    let mut ref_surf = brep.face_surface_world(the_ref_face);
    if let Some(Surface3::Trimmed(ts)) = ref_surf.as_ref() {
        ref_surf = Some(ts.basis.as_ref().clone());
    }
    // OCCT L1068-1072.
    let mut surf_face = brep.face_surface_world(the_face);
    if let Some(Surface3::Trimmed(ts)) = surf_face.as_ref() {
        surf_face = Some(ts.basis.as_ref().clone());
    }
    let (Some(ref_surf), Some(surf_face)) = (ref_surf, surf_face) else {
        return;
    };

    // OCCT L1074-1077.
    let mut to_modify = false;
    let mut to_translate = false;
    let mut to_rotate = false;
    let mut x_reverse = false;
    let mut y_reverse = false;
    let mut to_project = false;
    let mut a_translation = 0.0f64;
    let mut an_angle = 0.0f64;

    // OCCT L1079-1084: the elementary-surface frames (the rcad value model
    // carries the frame through ref_dir/u_dir — see elementary_frame).
    let elem_surf_face = elementary_frame(&surf_face);
    let elem_ref_surf = elementary_frame(&ref_surf);

    if let (Some((sx_origin, sx_vdir, sx_xdir)), Some((rx_origin, rx_vdir, rx_xdir))) =
        (elem_surf_face, elem_ref_surf)
    {
        // OCCT L1089-1096.
        let a_param = elclib_line_parameter_3d(sx_vdir, sx_origin, rx_origin);
        if a_param.abs() > p_confusion() {
            a_translation = -a_param;
        }
        // OCCT L1098-1108: CrossProd = X ^ Y = the surface normal.
        if rx_vdir.dot(sx_vdir) < 0.0 {
            x_reverse = true;
        }
        // OCCT L1110-1114.
        let scal_prod = sx_vdir.dot(rx_vdir);
        if scal_prod < 0.0 {
            y_reverse = true;
        }
        // OCCT L1116-1128.
        if !x_reverse && !y_reverse {
            // OCCT L1119-1122: the Ax3 Direct() flip — the rcad frames are
            // direct by construction (architecture note).
            let dir_ref = rx_vdir;
            an_angle = dir_angle_with_ref_3d(rx_xdir, sx_xdir, dir_ref);
        } else {
            an_angle = dir_angle_3d(rx_xdir, sx_xdir);
        }

        // OCCT L1130-1134.
        to_rotate = an_angle.abs() > p_confusion();
        to_translate = a_translation.abs() > p_confusion();
        to_modify = to_translate || to_rotate || x_reverse || y_reverse;
    } else {
        // OCCT L1136-1142.
        if !same_surf(&ref_surf, &surf_face) {
            to_project = true;
        }
    }

    // OCCT L1144-1147.
    let mut a_emap = MapOfShape::new();
    for an_edge in topexp_explorer(brep, the_face, ShapeType::Edge) {
        let mut an_edge = an_edge.clone();
        // OCCT L1150-1153.
        if is_edge_degenerated(brep, &an_edge) && to_modify {
            continue;
        }
        // OCCT L1155-1158.
        if to_project && brep_tool_is_closed_edge_face(brep, &an_edge, the_face) {
            continue;
        }
        // OCCT L1160-1163.
        if !map_add(&mut a_emap, &an_edge) {
            continue;
        }

        // OCCT L1165.
        let an_or = an_edge.orientation;

        // OCCT L1167-1179: the pcurve pair on the reference face.
        let pc_on_ref = brep.curve_on_surface(&an_edge, the_ref_face);
        let mut pc2: Option<Curve2d> = None;
        if let Some((pcr, _, _)) = pc_on_ref.as_ref() {
            an_edge.orientation = super::occt_reverse(an_edge.orientation);
            pc2 = brep
                .curve_on_surface(&an_edge, the_ref_face)
                .map(|(c, _, _)| c);
            an_edge.orientation = super::occt_reverse(an_edge.orientation);
            if let Some(p2) = pc2.as_ref() {
                if pcurve_handle_same(pcr, p2) {
                    // OCCT L1176: one pcurve (no seam on the reference).
                    continue;
                }
            }
        }

        // OCCT L1181-1187: the pcurves on the face.
        an_edge.orientation = Orientation::Forward;
        let pc0 = brep.curve_on_surface(&an_edge, the_face).map(|(c, _, _)| c);
        an_edge.orientation = super::occt_reverse(an_edge.orientation);
        let pc1 = brep.curve_on_surface(&an_edge, the_face).map(|(c, _, _)| c);
        an_edge.orientation = Orientation::Forward;
        let nb_pcurves = match (pc0.as_ref(), pc1.as_ref()) {
            (Some(a), Some(b)) => {
                if pcurve_handle_same(a, b) {
                    1
                } else {
                    2
                }
            }
            _ => 1,
        };

        // OCCT L1189-1220: the transformed pcurves.
        let (fpar, lpar) = {
            let r = brep_tool_range(brep, &an_edge);
            (r[0], r[1])
        };
        let mut new_pcurves: Vec<Option<Curve2d>> = vec![None, None];
        for ii in 0..nb_pcurves {
            let src = if ii == 0 {
                pc0.clone()
            } else {
                pc1.clone()
            };
            if to_project {
                // OCCT L1192-1197: the projection onto the reference surface.
                if let Some((c3d, f, l)) = brep_tool_curve(brep, &an_edge) {
                    let a_c3d = rcad_kernel::geom::Curve3::Trimmed(TrimmedCurve3::new(c3d, f, l));
                    new_pcurves[ii] = rcad_kernel::base::geom_proj_lib::curve2d_simple(
                        &a_c3d,
                        fpar,
                        lpar,
                        &ref_surf,
                    );
                }
                let _ = an_or;
            } else if let Some(s) = src {
                // OCCT L1201: the copy.
                let mut npc = geom2d_copy_untrim(&s);
                // OCCT L1203-1219.
                if to_translate {
                    npc = geom2d_translate(&npc, DVec2::new(0.0, a_translation));
                }
                if y_reverse {
                    npc = geom2d_mirror_ox2d(&npc);
                }
                if x_reverse {
                    npc = geom2d_mirror_oy2d(&npc);
                    npc = geom2d_translate(&npc, DVec2::new(2.0 * std::f64::consts::PI, 0.0));
                }
                if to_rotate {
                    npc = geom2d_translate(&npc, DVec2::new(an_angle, 0.0));
                }
                new_pcurves[ii] = Some(npc);
            }
        }

        an_edge.orientation = Orientation::Forward;

        // OCCT L1224-1262: the UpdateEdge dispatch.
        let mut builder = BRepBuilder::new();
        if nb_pcurves == 1 {
            let is_u_closed = ref_surf.is_u_closed();
            let is_v_closed = ref_surf.is_v_closed();
            if pc2.is_none() || (!is_u_closed && !is_v_closed) {
                // OCCT L1228-1229.
                if let Some(npc) = new_pcurves[0].as_ref() {
                    builder.update_edge_pcurve(
                        brep,
                        an_edge.clone(),
                        npc.clone(),
                        the_ref_face.clone(),
                        0.0,
                    );
                    map_add(the_map_edges_with_temporary_pcurves, &an_edge);
                }
            } else if let (Some(pcr), Some(npc)) =
                (pc_on_ref.as_ref().map(|(c, _, _)| c), new_pcurves[0].as_ref())
            {
                // OCCT L1232-1256: the may-be-same-pcurve seam check.
                let (a_umin, a_umax, a_vmin, a_vmax) = {
                    let d = ref_surf.default_domain();
                    (d[0], d[1], d[2], d[3])
                };
                let a_uperiod = if is_u_closed {
                    a_umax - a_umin
                } else {
                    0.0
                };
                let a_vperiod = if is_v_closed {
                    a_vmax - a_vmin
                } else {
                    0.0
                };
                let a_p2d_on_pcurve1 = pcr.point_at(fpar);
                let a_p2d_on_pcurve2 = npc.point_at(fpar);
                if (a_uperiod != 0.0
                    && (a_p2d_on_pcurve1.x - a_p2d_on_pcurve2.x).abs() > a_uperiod / 2.0)
                    || (a_vperiod != 0.0
                        && (a_p2d_on_pcurve1.y - a_p2d_on_pcurve2.y).abs() > a_vperiod / 2.0)
                {
                    // OCCT L1243-1255.
                    builder_update_edge_pcurve_null(brep, &an_edge, the_ref_face);
                    if an_or == Orientation::Forward {
                        builder.update_edge_pcurve_closed(
                            brep,
                            an_edge.clone(),
                            npc.clone(),
                            pcr.clone(),
                            the_ref_face.clone(),
                            fpar,
                            lpar,
                            0.0,
                        );
                    } else {
                        builder.update_edge_pcurve_closed(
                            brep,
                            an_edge.clone(),
                            pcr.clone(),
                            npc.clone(),
                            the_ref_face.clone(),
                            fpar,
                            lpar,
                            0.0,
                        );
                    }
                    map_add(the_map_edges_with_temporary_pcurves, &an_edge);
                }
            }
        } else if let (Some(np0), Some(np1)) = (new_pcurves[0].as_ref(), new_pcurves[1].as_ref()) {
            // OCCT L1258-1262.
            builder.update_edge_pcurve_closed(
                brep,
                an_edge.clone(),
                np0.clone(),
                np1.clone(),
                the_ref_face.clone(),
                fpar,
                lpar,
                0.0,
            );
            map_add(the_map_edges_with_temporary_pcurves, &an_edge);
        }

        // OCCT L1264: BB.Range(anEdge, fpar, lpar).
        builder.set_edge_range(brep, an_edge.clone(), fpar, lpar);
    }
}

/// OCCT static AddPCurves (cxx L1270-1294).
pub fn add_pcurves(
    brep: &mut BRep,
    the_faces: &[Shape],
    the_ref_face: &Shape,
    the_map_edges_with_temporary_pcurves: &mut MapOfShape,
) {
    // OCCT L1275-1281: the planar reference short-circuit.
    let ref_surf = brep.face_surface_world(the_ref_face);
    let is_plane = match ref_surf.as_ref() {
        Some(Surface3::Plane(_)) => true,
        Some(Surface3::Trimmed(ts)) => matches!(ts.basis.as_ref(), Surface3::Plane(_)),
        _ => false,
    };
    if is_plane {
        return;
    }

    for i in 1..=the_faces.len() {
        let mut a_face = the_faces[i - 1].clone();
        a_face.orientation = Orientation::Forward;
        if occt_is_same_shape(&a_face, the_ref_face) {
            continue;
        }
        // OCCT L1292.
        transform_pcurves(
            brep,
            the_ref_face,
            &a_face,
            the_map_edges_with_temporary_pcurves,
        );
    }
}

/// OCCT static AddOrdinaryEdges (cxx L1301-1352): adds the edges of the
/// shape to the sequence; seams and equal edges are dropped.  Returns
/// whether one of the original edges was dropped.
pub fn add_ordinary_edges(
    brep: &mut BRep,
    edges: &mut Vec<Shape>,
    a_shape: &Shape,
    an_index: &mut i32,
    the_removed_edges: &mut Vec<Shape>,
) -> bool {
    // OCCT L1307: the edge map.
    let mut a_new_edges = MapOfShape::new();
    // OCCT L1309-1321: the edges without seams.
    for exp in topexp_explorer(brep, a_shape, ShapeType::Edge) {
        let edge = exp.clone();
        if a_new_edges.contains_key(&shape_key(&edge)) {
            a_new_edges.shift_remove(&shape_key(&edge));
            the_removed_edges.push(edge);
        } else {
            map_add(&mut a_new_edges, &edge);
        }
    }

    // OCCT L1323-1343: merge edges and drop seams.
    let mut is_dropped = false;
    let mut i = 1i32;
    while i <= edges.len() as i32 {
        let current = edges[(i - 1) as usize].clone();
        if a_new_edges.contains_key(&shape_key(&current)) {
            a_new_edges.shift_remove(&shape_key(&current));
            edges.remove((i - 1) as usize);
            the_removed_edges.push(current);
            i -= 1;

            if !is_dropped {
                is_dropped = true;
                *an_index = i;
            }
        }
        i += 1;
    }

    // OCCT L1345-1349: add the edges to the sequence.
    for (_, e) in a_new_edges.iter() {
        edges.push(e.clone());
    }

    is_dropped
}

/// OCCT static getCylinder (cxx L1356-1431).
pub fn get_cylinder(the_in_surface: &Surface3, the_out_cylinder: &mut GpCylinder) -> bool {
    let mut is_cylinder = false;

    match the_in_surface {
        // OCCT L1360-1367.
        Surface3::Cylinder(g) => {
            the_out_cylinder.location = g.origin;
            the_out_cylinder.direction = g.axis;
            the_out_cylinder.radius = g.radius;
            is_cylinder = true;
        }
        // OCCT L1368-1396.
        Surface3::Revolution(rs) => {
            let mut a_basis = rs.profile.as_ref().clone();
            while let rcad_kernel::geom::Curve3::Trimmed(tc) = a_basis {
                a_basis = (*tc.curve).clone();
            }
            if let rcad_kernel::geom::Curve3::Line(basis_line) = a_basis {
                let a_dir = rs.axis_dir;
                let a_basis_dir = basis_line.direction;
                if dir_is_parallel_3d(a_basis_dir, a_dir, CONFUSION.max(1e-12)) {
                    // OCCT L1387-1393: the cylinder.
                    let a_loc = rs.axis_origin;
                    // OCCT L1389: aR = aBasisLine->Lin().Distance(aLoc).
                    let v = a_loc - basis_line.origin;
                    let t = v.dot(basis_line.direction);
                    let a_r = (a_loc - (basis_line.origin + basis_line.direction * t)).length();
                    the_out_cylinder.location = a_loc;
                    the_out_cylinder.direction = a_dir;
                    the_out_cylinder.radius = a_r;
                    is_cylinder = true;
                }
            }
        }
        // OCCT L1397-1425.
        Surface3::LinearExtrusion(les) => {
            let mut a_basis = les.profile.as_ref().clone();
            while let rcad_kernel::geom::Curve3::Trimmed(tc) = a_basis {
                a_basis = (*tc.curve).clone();
            }
            if let rcad_kernel::geom::Curve3::Circle(basis_circle) = a_basis {
                let a_dir = les.direction;
                // OCCT L1413: the basis circle axis direction.
                let a_basis_dir = basis_circle.normal;
                if dir_is_parallel_3d(a_basis_dir, a_dir, CONFUSION.max(1e-12)) {
                    // OCCT L1416-1422: the cylinder.
                    let a_loc = basis_circle.center;
                    let a_r = basis_circle.radius;
                    the_out_cylinder.location = a_loc;
                    the_out_cylinder.direction = a_dir;
                    the_out_cylinder.radius = a_r;
                    is_cylinder = true;
                }
            }
        }
        // OCCT L1426-1428: the empty else.
        _ => {}
    }

    is_cylinder
}

/// OCCT static ClearRts (cxx L1435-1440).
pub fn clear_rts(a_surface: &Surface3) -> Surface3 {
    match a_surface {
        Surface3::Trimmed(a_rts) => a_rts.basis.as_ref().clone(),
        other => other.clone(),
    }
}

/// OCCT static GetNormalToSurface (cxx L1447-1514): the normal to the
/// surface by the given parameter on the edge.
pub fn get_normal_to_surface(
    brep: &mut BRep,
    the_face: &Shape,
    the_edge: &Shape,
    the_p: f64,
    the_normal: &mut DVec3,
) -> bool {
    // OCCT L1454-1481: the 2D curve (the in-face occurrence for seams).
    let a_c2d;
    if brep_tool_is_closed_edge_face(brep, the_edge, the_face) {
        // OCCT L1457-1474: find the edge in the FORWARD face.
        let mut a_face = the_face.clone();
        a_face.orientation = Orientation::Forward;
        let mut an_edge_in_face = Shape::null();
        for an_edge in topexp_explorer(brep, &a_face, ShapeType::Edge) {
            if occt_is_same_shape(&an_edge, the_edge) {
                an_edge_in_face = an_edge;
                break;
            }
        }
        if an_edge_in_face.is_null() {
            return false;
        }
        a_c2d = brep
            .curve_on_surface(&an_edge_in_face, &a_face)
            .map(|(c, _, _)| c);
    } else {
        a_c2d = brep
            .curve_on_surface(the_edge, the_face)
            .map(|(c, _, _)| c);
    }

    let Some(a_c2d) = a_c2d else {
        // OCCT L1483-1486.
        return false;
    };

    // OCCT L1489-1490: the 2D point.
    let a_p2d = a_c2d.point_at(the_p);

    // OCCT L1493-1497: the LOCAL surface D1.
    let (a_s, a_loc) = brep_tool_surface_loc(brep, the_face);
    let Some(a_s) = a_s else { return false };
    let (_a_p3d, a_du, a_dv) = a_s.derivatives(a_p2d.x, a_p2d.y);

    // OCCT L1499-1504.
    let mut a_v_normal = a_du.cross(a_dv);
    if a_v_normal.length() < CONFUSION {
        return false;
    }

    // OCCT L1506-1509.
    if the_face.orientation == Orientation::Reversed {
        a_v_normal = -a_v_normal;
    }

    // OCCT L1511-1513: the location transform (rotation part).
    let trsf = brep.get_location(a_loc);
    if trsf != glam::DAffine3::IDENTITY {
        a_v_normal = trsf.transform_vector3(a_v_normal);
    }
    *the_normal = a_v_normal.normalize_or_zero();
    true
}

/// OCCT static IsSameDomain (cxx L1518-1642).
pub fn is_same_domain(
    brep: &BRep,
    a_face: &Shape,
    a_checked_face: &Shape,
    the_lin_tol: f64,
    the_ang_tol: f64,
    the_face_plane_map: &mut DataMapOfFacePlane,
) -> bool {
    // OCCT L1524-1534: the same-handle fast path (the rcad identity proxy is
    // the surface value image + the location index).
    let (s1_loc, l1) = brep_tool_surface_loc(brep, a_face);
    let (s2_loc, l2) = brep_tool_surface_loc(brep, a_checked_face);
    if let (Some(sv1), Some(sv2)) = (&s1_loc, &s2_loc) {
        if surface_handle_same(sv1, sv2) && l1 == l2 {
            return true;
        }
    }

    // OCCT L1536-1540.
    let mut s1 = brep.face_surface_world(a_face);
    let mut s2 = brep.face_surface_world(a_checked_face);
    let (Some(sv1), Some(sv2)) = (s1.take(), s2.take()) else {
        return false;
    };
    let s1 = clear_rts(&sv1);
    let s2 = clear_rts(&sv2);

    // OCCT L1548-1582: the two-planar-surfaces case.
    let a_planarity_checker1 = rcad_kernel::base::geom_lib::IsPlanarSurface::new(&s1, the_lin_tol);
    if a_planarity_checker1.is_planar() {
        let a_planarity_checker2 =
            rcad_kernel::base::geom_lib::IsPlanarSurface::new(&s2, the_lin_tol);
        if a_planarity_checker2.is_planar() {
            let a_pln1 = a_planarity_checker1.plan().clone();
            let a_pln2 = a_planarity_checker2.plan().clone();

            // OCCT L1559-1560.
            let parallel = dir_is_parallel_3d(a_pln1.normal, a_pln2.normal, the_ang_tol);
            let distance = (a_pln1.origin - a_pln2.origin).dot(a_pln1.normal).abs();
            if parallel && distance < the_lin_tol {
                // OCCT L1562-1577: the shared Geom_Plane.
                let a_plane_of_faces = if let Some(p) = the_face_plane_map.get(&shape_key(a_face))
                {
                    p.clone()
                } else if let Some(p) = the_face_plane_map.get(&shape_key(a_checked_face)) {
                    p.clone()
                } else {
                    a_pln1.clone()
                };
                the_face_plane_map.insert(shape_key(a_face), a_plane_of_faces.clone());
                the_face_plane_map.insert(shape_key(a_checked_face), a_plane_of_faces);
                return true;
            }
        }
    }

    // OCCT L1584-1610: the two-elementary-surfaces case (the OCCT tool chain
    // GeomAdaptor_Surface + BRepTopAdaptor_TopolTool +
    // IntPatch_ImpImpIntersection — the rcad ImpImpIntersection port folds
    // the TopolTool, architecture note).
    if is_elementary(&s1) && is_elementary(&s2) {
        let uv1 = s1.default_domain();
        let uv2 = s2.default_domain();
        let mut an_ii_int =
            crate::geomalgo::int_patch::imp_imp_intersection::ImpImpIntersection::new();
        an_ii_int.perform(&s1, &s2, uv1, uv2, the_lin_tol, the_lin_tol);
        if !an_ii_int.is_done() || an_ii_int.is_empty() {
            return false;
        }
        return an_ii_int.tangent_faces();
    }

    // OCCT L1612-1639: the swept/cylinder case.
    let s1_cyl_or_swept = matches!(
        s1,
        Surface3::Cylinder(_) | Surface3::Revolution(_) | Surface3::LinearExtrusion(_)
    );
    let s2_cyl_or_swept = matches!(
        s2,
        Surface3::Cylinder(_) | Surface3::Revolution(_) | Surface3::LinearExtrusion(_)
    );
    if s1_cyl_or_swept && s2_cyl_or_swept {
        let mut a_cyl1 = GpCylinder {
            location: DVec3::ZERO,
            direction: DVec3::ONE,
            radius: 0.0,
        };
        let mut a_cyl2 = GpCylinder {
            location: DVec3::ZERO,
            direction: DVec3::ONE,
            radius: 0.0,
        };
        if get_cylinder(&s1, &mut a_cyl1) && get_cylinder(&s2, &mut a_cyl2) {
            // OCCT L1622.
            if (a_cyl1.radius - a_cyl2.radius).abs() < the_lin_tol {
                // OCCT L1624-1626.
                if dir_is_parallel_3d(
                    a_cyl1.direction,
                    a_cyl2.direction,
                    rcad_kernel::precision::ANGULAR,
                ) {
                    // OCCT L1628-1635.
                    let a_vec12 = a_cyl2.location - a_cyl1.location;
                    if a_vec12.length_squared() < the_lin_tol * the_lin_tol
                        || dir_is_parallel_3d(
                            a_vec12,
                            a_cyl1.direction,
                            rcad_kernel::precision::ANGULAR,
                        )
                    {
                        return true;
                    }
                }
            }
        }
    }

    // OCCT L1641.
    false
}

/// OCCT IsKind(STANDARD_TYPE(Geom_ElementarySurface)) over the rcad value
/// model (Plane/Cylinder/Sphere/Cone/Torus).
fn is_elementary(s: &Surface3) -> bool {
    matches!(
        s,
        Surface3::Plane(_)
            | Surface3::Cylinder(_)
            | Surface3::Sphere(_)
            | Surface3::Cone(_)
            | Surface3::Torus(_)
    )
}

/// The elementary-surface frame triple `(location, VDirection, XDirection)`
/// (the gp_Ax3 Position() read; the rcad value model carries the frame
/// through ref_dir/u_dir — architecture bridge).
pub fn elementary_frame(s: &Surface3) -> Option<(DVec3, DVec3, DVec3)> {
    use rcad_kernel::geom::any_perpendicular;
    match s {
        Surface3::Plane(p) => Some((p.origin, p.normal, p.u_dir)),
        Surface3::Cylinder(c) => Some((c.origin, c.axis, c.ref_dir)),
        Surface3::Sphere(sp) => Some((sp.center, sp.axis, sp.ref_dir)),
        Surface3::Cone(c) => Some((c.apex, c.axis, any_perpendicular(c.axis))),
        Surface3::Torus(t) => Some((t.center, t.axis, t.ref_dir)),
        _ => None,
    }
}

/// OCCT static UpdateMapOfShapes (cxx L1646-1661).
pub fn update_map_of_shapes(
    brep: &mut BRep,
    the_map_of_shapes: &mut MapOfShape,
    the_context: &mut ShapeBuildReShape,
) {
    // OCCT iterates the map adding context images; rcad walks the
    // insertion-ordered map by index (the additions append at the tail and
    // are visited, matching the OCCT tail-insertion walk).
    let mut i = 0usize;
    while i < the_map_of_shapes.len() {
        let (_, a_shape) = the_map_of_shapes.get_index(i).unwrap();
        let a_shape = a_shape.clone();
        let a_context_shape = the_context.apply(brep, &a_shape, ShapeType::Shape);
        if !occt_is_same_shape(&a_context_shape, &a_shape) {
            map_add(the_map_of_shapes, &a_context_shape);
        }
        i += 1;
    }
}

/// OCCT static GlueEdgesWith3DCurves (cxx L1667-1750): glues the 3D curves
/// of the edge chain.
pub fn glue_edges_with_3d_curves(
    brep: &mut BRep,
    a_chain: &[Shape],
    first_vertex: &Shape,
    last_vertex: &Shape,
) -> Shape {
    // OCCT L1672.
    let a_curve_count = a_chain.len();

    // OCCT L1674-1676.
    let mut a_prev_edge = a_chain[0].clone();
    let mut a_prev_vertex = first_vertex.clone();
    let mut a_max_tolerance = 0.0f64;
    // OCCT L1678-1680.
    let mut a_bsplines: Vec<rcad_kernel::geom::BSplineCurve3> = Vec::with_capacity(a_curve_count);
    let mut a_vertices_tolerances: Vec<f64> = vec![0.0; a_curve_count];

    for i in 1..=a_curve_count {
        let a_current_edge = a_chain[i - 1].clone();
        // OCCT L1684-1685.
        let (a_current_first_vertex, a_current_last_vertex) =
            vertices(brep, &a_current_edge, false);
        // OCCT L1686.
        let a_to_reverse = !occt_is_same_shape(&a_current_first_vertex, &a_prev_vertex);

        // OCCT L1688-1690.
        a_max_tolerance = a_max_tolerance
            .max(brep.tolerance(&a_current_first_vertex))
            .max(brep.tolerance(&a_current_last_vertex));
        if i > 1 {
            // OCCT L1692-1695.
            if let Some(a_common_vertex) = common_vertex(brep, &a_prev_edge, &a_current_edge) {
                a_vertices_tolerances[i - 2] = brep.tolerance(&a_common_vertex);
            }
        }

        // OCCT L1698-1703.
        let (a_current_curve, fp, lp) = match brep_tool_curve(brep, &a_current_edge) {
            Some((c, f, l)) => (c, f, l),
            None => continue,
        };
        let a_trimmed =
            rcad_kernel::geom::Curve3::Trimmed(TrimmedCurve3::new(a_current_curve, fp, lp));
        let mut bspl = geom_convert_curve_to_bspline(&a_trimmed, fp, lp);
        // OCCT L1704.
        geom_convert_c0_to_c1(&mut bspl);
        // OCCT L1705-1708.
        if a_to_reverse {
            bspl = bspl.reversed();
        }
        a_bsplines.push(bspl);

        // OCCT L1710-1711.
        a_prev_vertex = if a_to_reverse {
            a_current_first_vertex
        } else {
            a_current_last_vertex
        };
        a_prev_edge = a_current_edge;
    }

    // OCCT L1714-1725: the C1 concatenation.
    let mut a_concat_curves = geom_convert_concat_c1(&a_bsplines, CONFUSION);

    // OCCT L1727-1737.
    if a_concat_curves.len() > 1 {
        let mut head = a_concat_curves[0].clone();
        for c in a_concat_curves.iter().skip(1) {
            geom_convert_comp_curve_add(&mut head, c, a_max_tolerance);
        }
        a_concat_curves[0] = head;
    }

    // OCCT L1739-1749.
    let a_res_curve = a_concat_curves[0].clone();
    let res_first = a_res_curve.knots[a_res_curve.degree];
    let res_last = a_res_curve.knots[a_res_curve.knots.len() - 1 - a_res_curve.degree];
    let a_res_edge = super::gap_deps::brep_lib_make_edge(
        brep,
        &rcad_kernel::geom::Curve3::BSpline(a_res_curve),
        first_vertex,
        last_vertex,
        res_first,
        res_last,
    );
    let mut builder = BRepBuilder::new();
    // OCCT L1746-1747.
    builder.set_edge_same_range(brep, a_res_edge.clone(), false);
    builder.set_edge_same_parameter(brep, a_res_edge.clone(), false);
    // OCCT L1748: BRepLib::SameParameter(aResEdge, aMaxTolerance, true).
    crate::shhealing::shape_fix::shape_fix_gap_deps::brep_lib_same_parameter_edge(
        brep,
        &a_res_edge,
        a_max_tolerance,
        true,
    );
    a_res_edge
}
