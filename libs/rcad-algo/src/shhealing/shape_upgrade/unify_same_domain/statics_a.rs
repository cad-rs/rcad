//! OCCT ShapeUpgrade_UnifySameDomain.cxx L93-533 — the first file-statics
//! segment: `IsOnSingularity` (L93-106), `IsUiso` (L112-119), `IsLinear`
//! (L121-140), `SplitWire` forward declaration (L142-145; the body is
//! translated in `split_wire.rs` at cxx L4560-4687), `TrueValueOfOffset`
//! (L147-155), `GetFaceFromSeq` (L165-237), `UpdateBoundaries` (L241-266),
//! `TryMakeLine` (L268-304), `RemoveEdgeFromMap` (L310-342),
//! `ComputeMinEdgeSize` (L344-381), `FindCoordBounds` (L389-508),
//! `getCurveParams` (L518-532).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2dEval, Line2d};
use rcad_kernel::precision::{is_infinite_value, CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, ShapeType};

use super::gap_deps::{BRepAdaptorCurve, BRepAdaptorCurve2d, GeomAbsCurveType};
use super::topexp::{
    brep_tool_is_closed_edge_face, is_edge_degenerated, occt_is_same_shape, topexp_explorer,
    vertices,
};
use super::{map_add, shape_key, IndexedDataMapOfShapeListOfShape, MapOfShape};

/// OCCT RealLast() / RealFirst().
pub(crate) const REAL_LAST: f64 = f64::MAX;
pub(crate) const REAL_FIRST: f64 = -f64::MAX;

/// OCCT static IsOnSingularity (cxx L93-106): true when the list carries a
/// degenerated edge.
pub fn is_on_singularity(brep: &BRep, the_edge_list: &[Shape]) -> bool {
    for an_edge in the_edge_list {
        // OCCT L99: BRep_Tool::Degenerated(anEdge).
        if is_edge_degenerated(brep, an_edge) {
            return true;
        }
    }
    false
}

/// OCCT static IsUiso (cxx L112-119): only for seam edges assumed to be U-
/// or V- isolines.
pub fn is_uiso(brep: &BRep, the_edge: &Shape, the_face: &Shape) -> bool {
    let a_ba_curve2d = BRepAdaptorCurve2d::new(brep, the_edge, the_face);
    let (_a_p2d, a_vec) = a_ba_curve2d.d1(a_ba_curve2d.first_parameter());
    a_vec.y.abs() > a_vec.x.abs()
}

/// OCCT static IsLinear (cxx L121-140): the line detection for the two
/// mergeable analytic families; `the_dir` receives the direction.
pub fn is_linear(the_ba_curve: &BRepAdaptorCurve, the_dir: &mut DVec3) -> bool {
    // OCCT L123: GeomAbs_CurveType aType = theBAcurve.GetType().
    let a_type = the_ba_curve.get_type();

    if a_type == GeomAbsCurveType::Line {
        // OCCT L127: theDir = theBAcurve.Line().Position().Direction().
        if let Some(l) = the_ba_curve.line() {
            *the_dir = l.direction;
        }
        return true;
    }

    // OCCT L131-137: the two-pole BSpline/Bezier case.
    if (a_type == GeomAbsCurveType::BezierCurve || a_type == GeomAbsCurveType::BSplineCurve)
        && the_ba_curve.nb_poles() == 2
    {
        let a_first_pnt = the_ba_curve.value(the_ba_curve.first_parameter());
        let a_last_pnt = the_ba_curve.value(the_ba_curve.last_parameter());
        *the_dir = (a_last_pnt - a_first_pnt).normalize_or_zero();
        return true;
    }

    false
}

/// OCCT static TrueValueOfOffset (cxx L147-155).
pub fn true_value_of_offset(the_value: f64, the_period: f64) -> f64 {
    if the_value > 0.0 {
        return the_period;
    }
    -the_period
}

/// OCCT static GetFaceFromSeq (cxx L165-237): gets a face from the sequence
/// that has a common boundary with the already-used faces (preferring a
/// common 2D boundary, i.e. excluding seam edges).  Returns
/// `(face, theIndFace)`; a null face keeps OCCT's not-found path.
pub fn get_face_from_seq(
    brep: &mut BRep,
    the_faces: &[Shape],
    the_ref_face: &Shape,
    the_map_ef: &IndexedDataMapOfShapeListOfShape,
    the_used_faces: &mut MapOfShape,
    the_ind_face: &mut usize,
) -> Shape {
    // OCCT L174-179.
    if the_used_faces.is_empty() {
        *the_ind_face = 1;
        map_add(the_used_faces, &the_faces[0]);
        return the_faces[0].clone();
    }

    // OCCT L181-182.
    let mut a_seam_edge = Shape::null();
    let mut a_target_face = Shape::null();
    for ii in 1..=the_used_faces.len() {
        let a_used_face = the_used_faces.get_index(ii - 1).unwrap().1.clone();
        for an_edge in topexp_explorer(brep, &a_used_face, ShapeType::Edge) {
            // OCCT L190: theMapEF.FindFromKey(anEdge).
            let a_face_list = match the_map_ef.get(&shape_key(&an_edge)) {
                Some((_, l)) => l,
                None => continue,
            };
            if a_face_list.len() < 2 {
                continue;
            }

            // OCCT L196-198: the face on the other side of the edge.
            let a_face = if occt_is_same_shape(&a_face_list[0], &a_used_face) {
                a_face_list[a_face_list.len() - 1].clone()
            } else {
                a_face_list[0].clone()
            };

            if the_used_faces.contains_key(&shape_key(&a_face)) {
                continue;
            }

            // OCCT L205: BRep_Tool::IsClosed(anEdge, theRefFace).
            if brep_tool_is_closed_edge_face(brep, &an_edge, the_ref_face) {
                a_seam_edge = an_edge;
                continue;
            }

            a_target_face = a_face;
            break;
        }
        // OCCT L214-219: the seam-only fallback.
        if a_target_face.is_null() && !a_seam_edge.is_null() {
            if let Some((_, a_face_list)) = the_map_ef.get(&shape_key(&a_seam_edge)) {
                a_target_face = if occt_is_same_shape(&a_face_list[0], &a_used_face) {
                    a_face_list[a_face_list.len() - 1].clone()
                } else {
                    a_face_list[0].clone()
                };
            }
        }
        if !a_target_face.is_null() {
            break;
        }
    }

    // OCCT L226.
    map_add(the_used_faces, &a_target_face);
    // OCCT L227-234.
    for ii in 2..=the_faces.len() {
        if occt_is_same_shape(&the_faces[ii - 1], &a_target_face) {
            *the_ind_face = ii;
            break;
        }
    }

    a_target_face
}

/// OCCT static UpdateBoundaries (cxx L241-266).
pub fn update_boundaries(
    the_pcurve: &rcad_kernel::geom::Curve2d,
    the_first: f64,
    the_last: f64,
    the_ind_coord: usize,
    the_min_coord: &mut f64,
    the_max_coord: &mut f64,
) {
    // OCCT L248-249.
    const NB_SAMPLES: i32 = 4;
    let delta = (the_last - the_first) / NB_SAMPLES as f64;

    for i in 0..=NB_SAMPLES {
        // OCCT L253-254.
        let a_param = the_first + i as f64 * delta;
        let a_point = the_pcurve.point_at(a_param);
        let a_coord = if the_ind_coord == 1 {
            a_point.x
        } else {
            a_point.y
        };

        // OCCT L256-264.
        if a_coord < *the_min_coord {
            *the_min_coord = a_coord;
        }
        if a_coord > *the_max_coord {
            *the_max_coord = a_coord;
        }
    }
}

/// OCCT static TryMakeLine (cxx L268-304): builds the Geom2d_Line when the
/// pcurve is a straight segment (the out handle maps to the Option).
pub fn try_make_line(
    the_pcurve: &rcad_kernel::geom::Curve2d,
    the_first: f64,
    the_last: f64,
) -> Option<Line2d> {
    // OCCT L273-275.
    let a_first_pnt = the_pcurve.point_at(the_first);
    let a_last_pnt = the_pcurve.point_at(the_last);
    let a_vec = a_last_pnt - a_first_pnt;
    // OCCT L276-281.
    let a_sq_len = a_vec.length_squared();
    let a_sq_param_len = (the_last - the_first) * (the_last - the_first);
    if (a_sq_len - a_sq_param_len).abs() > CONFUSION {
        return None;
    }

    // OCCT L283-287.
    let a_dir = a_vec.normalize_or_zero();
    let an_origin = a_first_pnt - a_dir * the_first;
    let a_lin = Line2d::new(an_origin, a_dir);

    // OCCT L289-300.
    const NB_SAMPLES: i32 = 10;
    let a_delta = (the_last - the_first) / NB_SAMPLES as f64;
    for i in 1..NB_SAMPLES {
        let a_param = the_first + i as f64 * a_delta;
        let a_pnt = the_pcurve.point_at(a_param);
        // OCCT L295: aLin.Distance(aPnt).
        let v = a_pnt - a_lin.origin;
        let t = v.dot(a_lin.direction);
        let proj = a_lin.origin + a_lin.direction * t;
        let a_dist = a_pnt.distance(proj);
        if a_dist > CONFUSION {
            return None;
        }
    }

    // OCCT L302: theLine = new Geom2d_Line(aLin).
    Some(a_lin)
}

/// OCCT static RemoveEdgeFromMap (cxx L310-342): removes the specified edge
/// from the vertex-to-edge map.
pub fn remove_edge_from_map(
    brep: &mut BRep,
    the_edge: &Shape,
    the_vertex_to_edges: &mut IndexedDataMapOfShapeListOfShape,
) -> bool {
    // OCCT L315-318.
    let mut an_is_removed = false;
    let (a_first_vertex, a_last_vertex) = vertices(brep, the_edge, false);
    for a_vertex in [a_first_vertex, a_last_vertex] {
        // OCCT L321-324.
        let key = shape_key(&a_vertex);
        if !the_vertex_to_edges.contains_key(&key) {
            continue;
        }
        // OCCT L325-339: the remove-while-iterating walk.
        let (_, a_vertex_edges) = the_vertex_to_edges.get_mut(&key).unwrap();
        let mut it = 0usize;
        while it < a_vertex_edges.len() {
            let an_edge = a_vertex_edges[it].clone();
            if occt_is_same_shape(&an_edge, the_edge) {
                an_is_removed = true;
                a_vertex_edges.remove(it);
            } else {
                it += 1;
            }
        }
    }
    an_is_removed
}

/// OCCT static ComputeMinEdgeSize (cxx L344-381).
pub fn compute_min_edge_size(
    brep: &mut BRep,
    the_edges: &[Shape],
    the_ref_face: &Shape,
    the_edges_map: &mut MapOfShape,
) -> f64 {
    // OCCT L349.
    let mut min_size = REAL_LAST;

    for ind in 1..=the_edges.len() {
        let an_edge = the_edges[ind - 1].clone();
        // OCCT L354.
        map_add(the_edges_map, &an_edge);
        // OCCT L355-356.
        let (v1, v2) = vertices(brep, &an_edge, false);
        // OCCT L357-361.
        let ba_curve2d = BRepAdaptorCurve2d::new(brep, &an_edge, the_ref_face);
        if ba_curve2d.curve().is_none() {
            continue;
        }

        // OCCT L363-364.
        let first_p2d = ba_curve2d.value(ba_curve2d.first_parameter());
        let last_p2d = ba_curve2d.value(ba_curve2d.last_parameter());
        // OCCT L366-375.
        let a_sq_dist;
        if occt_is_same_shape(&v1, &v2) && !is_edge_degenerated(brep, &an_edge) {
            let mid_p2d = ba_curve2d
                .value((ba_curve2d.first_parameter() + ba_curve2d.last_parameter()) / 2.0);
            a_sq_dist = first_p2d.distance_squared(mid_p2d);
        } else {
            a_sq_dist = first_p2d.distance_squared(last_p2d);
        }

        // OCCT L377.
        min_size = min_size.min(a_sq_dist);
    }
    // OCCT L379-380.
    min_size.sqrt()
}

/// OCCT static FindCoordBounds (cxx L389-508): searching for the origin of
/// U (or V) in the 2D space.  Returns false when no curve on surface was
/// found (or no interval exists).
#[allow(clippy::too_many_arguments)]
pub fn find_coord_bounds(
    brep: &mut BRep,
    the_faces: &[Shape],
    the_ref_face: &Shape,
    the_map_ef: &IndexedDataMapOfShapeListOfShape,
    the_edges_map: &MapOfShape,
    the_ind_coord: usize,
    the_period: f64,
    the_min_coord: &mut f64,
    the_max_coord: &mut f64,
    the_number_of_intervals: &mut i32,
    the_ind_face_max: &mut i32,
) -> bool {
    // OCCT L403: aPairSeq.
    let mut a_pair_seq: Vec<(f64, f64)> = Vec::new();

    // OCCT L405-406.
    let mut a_simple_max = REAL_FIRST;
    *the_ind_face_max = 0;

    // OCCT L408-409.
    let mut a_used_faces = MapOfShape::new();
    let a_last_face = Shape::null(); // declared-unused in OCCT (cxx L409)
    let _ = a_last_face;

    for _ii in 1..=the_faces.len() {
        let mut an_ind_face: usize = 0;
        // OCCT L413-414: get a face bordering with the previous ones.
        let a_face = get_face_from_seq(
            brep,
            the_faces,
            the_ref_face,
            the_map_ef,
            &mut a_used_faces,
            &mut an_ind_face,
        );
        let mut a_min_coord = REAL_LAST;
        let mut a_max_coord = REAL_FIRST;
        for an_edge in topexp_explorer(brep, &a_face, ShapeType::Edge) {
            // OCCT L420-423.
            if !the_edges_map.contains_key(&shape_key(&an_edge)) {
                continue;
            }
            // OCCT L424-429.
            let Some((a_pcurve, fpar, lpar)) = brep.curve_on_surface(&an_edge, the_ref_face) else {
                return false;
            };
            update_boundaries(
                &a_pcurve,
                fpar,
                lpar,
                the_ind_coord,
                &mut a_min_coord,
                &mut a_max_coord,
            );
        }

        // OCCT L433-436.
        if is_infinite_value(a_min_coord) || is_infinite_value(a_max_coord) {
            continue;
        }

        // OCCT L438-442.
        if a_max_coord > a_simple_max {
            a_simple_max = a_max_coord;
            *the_ind_face_max = an_ind_face as i32;
        }

        // OCCT L444-476: insert the new interval into the sequence.
        let mut an_is_in_interval = false;
        for jj in 1..=a_pair_seq.len() {
            let a_local_min = a_pair_seq[jj - 1].0;
            let a_local_max = a_pair_seq[jj - 1].1;
            if a_min_coord >= a_local_min
                && a_min_coord <= a_local_max
                && a_max_coord >= a_local_min
                && a_max_coord <= a_local_max
            {
                an_is_in_interval = true;
                break;
            }

            if a_min_coord < a_local_min && a_max_coord >= a_local_min && a_max_coord <= a_local_max
            {
                a_pair_seq[jj - 1].0 = a_min_coord;
                an_is_in_interval = true;
                break;
            } else if a_min_coord < a_local_min && a_max_coord > a_local_max {
                a_pair_seq[jj - 1].0 = a_min_coord;
                a_pair_seq[jj - 1].1 = a_max_coord;
                an_is_in_interval = true;
                break;
            } else if a_min_coord >= a_local_min
                && a_min_coord <= a_local_max
                && a_max_coord > a_local_max
            {
                a_pair_seq[jj - 1].1 = a_max_coord;
                an_is_in_interval = true;
                break;
            }
        }
        // OCCT L477-488.
        if !an_is_in_interval {
            let an_interval = (a_min_coord, a_max_coord);
            if !a_pair_seq.is_empty() && a_max_coord < a_pair_seq[0].0 {
                a_pair_seq.insert(0, an_interval);
            } else {
                a_pair_seq.push(an_interval);
            }
        }
    }

    // OCCT L491.
    *the_number_of_intervals = a_pair_seq.len() as i32;

    // OCCT L493-504.
    if a_pair_seq.len() == 2 {
        *the_min_coord = a_pair_seq[1].0 - the_period;
    } else if !a_pair_seq.is_empty() {
        *the_min_coord = a_pair_seq[0].0;
    } else {
        return false;
    }

    // OCCT L506-507.
    *the_max_coord = a_pair_seq[0].1;
    true
}

/// OCCT static getCurveParams (cxx L518-532): the start and end points of
/// the edge in the parametric space of the face (orientation-aware).
pub fn get_curve_params(brep: &BRep, the_edge: &Shape, the_ref_face: &Shape) -> (DVec2, DVec2) {
    let a_curve_adaptor = BRepAdaptorCurve2d::new(brep, the_edge, the_ref_face);
    let mut a_first_param = a_curve_adaptor.first_parameter();
    let mut a_last_param = a_curve_adaptor.last_parameter();
    // OCCT L524-527.
    if the_edge.orientation != Orientation::Forward {
        std::mem::swap(&mut a_first_param, &mut a_last_param);
    }

    // OCCT L529-531.
    let a_first_point = a_curve_adaptor.value(a_first_param);
    let a_last_point = a_curve_adaptor.value(a_last_param);
    (a_first_point, a_last_point)
}
