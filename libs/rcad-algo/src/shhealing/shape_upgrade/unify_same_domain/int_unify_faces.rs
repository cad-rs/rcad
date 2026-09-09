//! OCCT ShapeUpgrade_UnifySameDomain.cxx L3185-4320 — `IntUnifyFaces` (the
//! per-shell / per-face-compound unification walk) and the static
//! helpers it inlines from the segment: the seam search, the periodic
//! relocation, the new-wire walk and the context merges.

use std::sync::Arc;

use crate::shhealing::shape_construct::gap_deps::{surface_u_period, surface_v_period};
use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Surface3, SurfaceEval, TrimmedSurface};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    tshape_flags, BRep, BRepBuilder, BRepTool, Orientation, ShapeType, TShape,
};

use super::gap_deps::{
    brep_lib_build_pcurve_for_edge_on_plane, brep_lib_continuity_of_faces, geom2d_trimmed,
    geom_bspline_surface_period_span_u, geom_bspline_surface_period_span_v,
    geom_convert_approx_surface, surface_handle_same,
};
use super::statics_a::{compute_min_edge_size, find_coord_bounds, is_on_singularity, is_uiso};
use super::statics_b::{
    add_ordinary_edges, add_pcurves, get_normal_to_surface, insert_wires_into_faces,
    is_same_domain, relocate_pcurves_to_new_uorigin, DataMapOfShapePCurve,
};
use super::topexp::{
    is_edge_degenerated, last_vertex, map_shapes_all, map_shapes_and_ancestors,
    map_shapes_and_unique_ancestors, occt_is_same_shape, topexp_explorer, vertices,
};
use super::unify_faces::is_same_sets;
use super::{
    map_add, occt_reverse, shape_key, DataMapOfShapeMapOfShape, IndexedDataMapOfShapeListOfShape,
    MapOfShape, ShapeUpgradeUnifySameDomain,
};
use crate::brep_sweep::tool_rehost::brep_tools_is_really_closed;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, set_flag_inplace};
use crate::shhealing::shape_build::edge::builder_update_edge_pcurve_null;

/// OCCT BRep_Builder::UpdateFace(F, S, L, Tol) (BRep_Builder.cxx L433-486) —
/// the surface/tolerance payload write through the shared handle.
fn brep_builder_update_face(brep: &mut BRep, face: &Shape, surf: Surface3, tol: f64) {
    // SAFETY: the in-place pattern of BRep::edge_mut_inplace (the healing
    // chain mutates the TShape through the shared handle).
    let ptr = Arc::as_ptr(&brep.tshapes[face.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    if let TShape::Face(fd) = ts {
        fd.surface = Some(surf);
        fd.tolerance = fd.tolerance.max(tol);
    }
}

/// OCCT BRepTools::OuterWire(F) — the rcad stored outer-wire read (the
/// brep_offset_api_middle_path.rs precedent; the TFaceData model carries
/// the outer wire explicitly).
fn outer_wire(brep: &BRep, the_face: &Shape) -> Shape {
    match &*brep.tshapes[the_face.index] {
        TShape::Face(fd) => fd.outer_wire.clone(),
        _ => Shape::null(),
    }
}

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L3185-4320: IntUnifyFaces.
    pub(crate) fn int_unify_faces(
        &mut self,
        brep: &mut BRep,
        the_inp_shape: &Shape,
        the_gmap_edge_faces: &IndexedDataMapOfShapeListOfShape,
        the_gmap_face_shells: &DataMapOfShapeMapOfShape,
        the_free_bound_map: &MapOfShape,
    ) {
        // OCCT L3193-3196: the local edge -> faces map.
        let mut a_map_edge_faces = IndexedDataMapOfShapeListOfShape::new();
        map_shapes_and_ancestors(
            brep,
            the_inp_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_map_edge_faces,
        );

        // OCCT L3198-3199: the processed map.
        let mut a_processed = MapOfShape::new();

        // OCCT L3201-3203: processing each face.
        for exp in topexp_explorer(brep, the_inp_shape, ShapeType::Face) {
            let a_face = exp.clone();
            if a_processed.contains_key(&shape_key(&a_face)) {
                continue;
            }

            // OCCT L3211-3216: the base surface (Bug 33894: no surface).
            let a_base_surface = match brep.face_surface_world(&a_face) {
                Some(s) => super::statics_b::clear_rts(&s),
                None => continue,
            };

            // OCCT L3218-3220: the boundary edges for the new face.
            let mut edges: Vec<Shape> = Vec::new();
            let mut removed_edges: Vec<Shape> = Vec::new();
            let mut dummy: i32 = 0;
            add_ordinary_edges(brep, &mut edges, &a_face, &mut dummy, &mut removed_edges);

            // OCCT L3225-3229: the faces to unify.
            let mut faces: Vec<Shape> = Vec::new();
            faces.push(a_face.clone());

            // OCCT L3231-3233.
            let ref_face_orientation = a_face.orientation;

            // OCCT L3235-3239: the original-surface reference face.
            let mut a_bb = BRepBuilder::new();
            let mut ref_face = a_bb.make_face(brep, Some(a_base_surface.clone()), Shape::null());
            ref_face.orientation = ref_face_orientation;
            // OCCT L3240-3241.
            let mut map_edges_with_temporary_pcurves = MapOfShape::new();

            // OCCT L3244-3245.
            let mut uperiod = if a_base_surface.is_u_periodic() {
                surface_u_period(&a_base_surface)
            } else {
                0.0
            };
            let mut vperiod = if a_base_surface.is_v_periodic() {
                surface_v_period(&a_base_surface)
            } else {
                0.0
            };

            // OCCT L3247-3249: the shells connected to the face.
            let pf_shells1 = the_gmap_face_shells.get(&shape_key(&a_face));

            // OCCT L3251-3350: the adjacent-face walk (the `i = dummy`
            // re-entry of OCCT L3342 becomes a while index).
            let mut i = 1usize;
            while i <= edges.len() {
                let edge = edges[i - 1].clone();
                if is_edge_degenerated(brep, &edge) {
                    i += 1;
                    continue;
                }

                // OCCT L3261-3269.
                let a_glist_len = the_gmap_edge_faces
                    .get(&shape_key(&edge))
                    .map(|(_, l)| l.len())
                    .unwrap_or(0);
                if !self.my_allow_internal
                    && (a_glist_len != 2
                        || self.my_keep_shapes.contains_key(&shape_key(&edge))
                        || the_free_bound_map.contains_key(&shape_key(&edge)))
                {
                    // OCCT L3267: non manifold case not processed.
                    i += 1;
                    continue;
                }
                // OCCT L3271-3276.
                let a_list = a_map_edge_faces
                    .get(&shape_key(&edge))
                    .map(|(_, l)| l.clone())
                    .unwrap_or_default();
                if a_list.len() < 2 {
                    i += 1;
                    continue;
                }

                // OCCT L3278-3283.
                if !self.my_safe_input_mode
                    && matches!(a_base_surface, Surface3::Plane(_) | Surface3::Trimmed(_))
                {
                    let is_plane = matches!(
                        super::statics_b::clear_rts(&a_base_surface),
                        Surface3::Plane(_)
                    );
                    if is_plane {
                        brep_lib_build_pcurve_for_edge_on_plane(brep, &edge, &a_face);
                    }
                }

                // OCCT L3285-3293: the face normal at the edge midpoint.
                let a_dn1: glam::DVec3;
                let range = super::topexp::brep_tool_range(brep, &edge);
                let a_t_mid = (range[0] + range[1]) * 0.5;
                let mut a_dn1_val = glam::DVec3::ONE;
                let b_check_normals =
                    get_normal_to_surface(brep, &a_face, &edge, a_t_mid, &mut a_dn1_val);
                a_dn1 = a_dn1_val;

                // OCCT L3295-3349.
                for a_checked_face in &a_list {
                    if occt_is_same_shape(a_checked_face, &a_face) {
                        continue;
                    }
                    if a_processed.contains_key(&shape_key(a_checked_face)) {
                        continue;
                    }

                    // OCCT L3311-3319: the shell-set guard.
                    let pf_shells2 = the_gmap_face_shells.get(&shape_key(a_checked_face));
                    if !is_same_sets(pf_shells1, pf_shells2) {
                        continue;
                    }

                    // OCCT L3321-3334: the normal agreement.
                    if b_check_normals {
                        let mut a_dn2 = glam::DVec3::ONE;
                        if get_normal_to_surface(brep, a_checked_face, &edge, a_t_mid, &mut a_dn2) {
                            let an_angle = super::gap_deps::dir_angle_3d(a_dn1, a_dn2);
                            if an_angle > self.my_ang_tol {
                                continue;
                            }
                        }
                    }
                    // OCCT L3336-3348.
                    if is_same_domain(
                        brep,
                        &a_face,
                        a_checked_face,
                        self.my_lin_tol,
                        self.my_ang_tol,
                        &mut self.my_face_plane_map,
                    ) {
                        if add_ordinary_edges(
                            brep,
                            &mut edges,
                            a_checked_face,
                            &mut dummy,
                            &mut removed_edges,
                        ) {
                            // OCCT L3341-3343: the sequence edges modified.
                            i = dummy as usize;
                        }
                        faces.push(a_checked_face.clone());
                        map_add(&mut a_processed, a_checked_face);
                        break;
                    }
                }
                i += 1;
            }

            // OCCT L3352-3479: the multi-face post-processing.
            if faces.len() > 1 {
                // OCCT L3354-3359: the planar reference surface update.
                if let Some(a_plane) = self.my_face_plane_map.get(&shape_key(&faces[0])) {
                    brep_builder_update_face(
                        brep,
                        &ref_face,
                        Surface3::Plane(a_plane.clone()),
                        CONFUSION,
                    );
                }
                // OCCT L3360-3363.
                let mut f_ref_face = ref_face.clone();
                f_ref_face.orientation = Orientation::Forward;
                add_pcurves(
                    brep,
                    &faces,
                    &f_ref_face,
                    &mut map_edges_with_temporary_pcurves,
                );

                // OCCT L3365-3373: the connectivity map for the faces.
                let mut a_map_ef = IndexedDataMapOfShapeListOfShape::new();
                for i in 1..=faces.len() {
                    map_shapes_and_ancestors(
                        brep,
                        &faces[i - 1],
                        ShapeType::Edge,
                        ShapeType::Face,
                        &mut a_map_ef,
                    );
                }
                // OCCT L3374-3389: the keep / multi-connected edges.
                let mut a_keep_edges: Vec<Shape> = Vec::new();
                for i in 1..=a_map_ef.len() {
                    let (a_e, a_lf) = {
                        let (_, (s, l)) = a_map_ef.get_index(i - 1).unwrap();
                        (s.clone(), l.clone())
                    };
                    if a_lf.len() == 2 {
                        let a_glf_len = the_gmap_edge_faces
                            .get(&shape_key(&a_e))
                            .map(|(_, l)| l.len())
                            .unwrap_or(0);
                        if a_glf_len > 2
                            || self.my_keep_shapes.contains_key(&shape_key(&a_e))
                            || the_free_bound_map.contains_key(&shape_key(&a_e))
                        {
                            a_keep_edges.push(a_e);
                        }
                    }
                }
                if !a_keep_edges.is_empty() {
                    if !self.my_allow_internal {
                        // OCCT L3394-3436: the avoid-faces pass.
                        let mut an_avoid_faces = MapOfShape::new();
                        for a_e in &a_keep_edges {
                            if let Some((_, a_lf)) = a_map_ef.get(&shape_key(a_e)) {
                                if let Some(f) = a_lf.first() {
                                    map_add(&mut an_avoid_faces, f);
                                }
                                if let Some(f) = a_lf.last() {
                                    map_add(&mut an_avoid_faces, f);
                                }
                            }
                        }
                        let mut i = 1usize;
                        while i <= faces.len() {
                            if an_avoid_faces.contains_key(&shape_key(&faces[i - 1])) {
                                // OCCT L3409-3434: the boundary update.
                                let mut has_connect_another_faces = false;
                                if let Some((_, a_lf)) = a_map_ef.get(&shape_key(&faces[i - 1])) {
                                    if a_lf.len() > 1 {
                                        for it in a_lf {
                                            if !an_avoid_faces.contains_key(&shape_key(it)) {
                                                has_connect_another_faces = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                                if !has_connect_another_faces {
                                    add_ordinary_edges(
                                        brep,
                                        &mut edges,
                                        &faces[i - 1],
                                        &mut dummy,
                                        &mut removed_edges,
                                    );
                                    faces.remove(i - 1);
                                    continue; // i unchanged (the OCCT i--)
                                }
                            }
                            i += 1;
                        }
                        // OCCT L3437-3467: the keep-edge containment pass.
                        if !faces.is_empty() {
                            let mut a_map_faces = MapOfShape::new();
                            for f in &faces {
                                map_add(&mut a_map_faces, f);
                            }
                            for it in &a_keep_edges {
                                let Some((_, a_lf)) = a_map_ef.get(&shape_key(it)) else {
                                    continue;
                                };
                                if a_lf.len() < 2 {
                                    continue;
                                }
                                let first_in = a_map_faces.contains_key(&shape_key(&a_lf[0]));
                                let last_in =
                                    a_map_faces.contains_key(&shape_key(&a_lf[a_lf.len() - 1]));
                                if first_in && last_in {
                                    let mut i = 0usize;
                                    while i < faces.len() {
                                        if occt_is_same_shape(&faces[i], &a_lf[0])
                                            || occt_is_same_shape(&faces[i], &a_lf[a_lf.len() - 1])
                                        {
                                            add_ordinary_edges(
                                                brep,
                                                &mut edges,
                                                &faces[i],
                                                &mut dummy,
                                                &mut removed_edges,
                                            );
                                            faces.remove(i);
                                            continue;
                                        }
                                        i += 1;
                                    }
                                }
                            }
                        }
                    } else {
                        // OCCT L3469-3478: the internal keep edges.
                        for it in &a_keep_edges {
                            let mut a_e = it.clone();
                            a_e.orientation = Orientation::Internal;
                            edges.push(a_e);
                        }
                    }
                }
            }

            // OCCT L3482-3489: the unique-ancestor map (fresh aMapEF).
            let mut a_map_ef = IndexedDataMapOfShapeListOfShape::new();
            for i in 1..=faces.len() {
                map_shapes_and_unique_ancestors(
                    brep,
                    &faces[i - 1],
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut a_map_ef,
                );
            }

            // OCCT L3491-3506: the edge orientation correction.
            let mut ii = 0usize;
            while ii < edges.len() {
                let an_edge = edges[ii].clone();
                let Some(ind_e) = a_map_ef.get_index_of(&shape_key(&an_edge)) else {
                    ii += 1;
                    continue;
                };
                let a_lf_len = a_map_ef.get_index(ind_e).unwrap().1 .1.len();
                if self.my_allow_internal
                    && self.my_keep_shapes.contains_key(&shape_key(&an_edge))
                    && a_lf_len == 2
                {
                    edges[ii].orientation = Orientation::Internal;
                }
                if edges[ii].orientation != Orientation::Internal {
                    // OCCT L3504: the map's representative occurrence.
                    edges[ii] = a_map_ef.get_index(ind_e).unwrap().1 .0.clone();
                }
                ii += 1;
            }

            // OCCT L3508-3523: exclude the internal edges.
            let mut internal_edges = MapOfShape::new();
            let mut ind_e = 0usize;
            while ind_e < edges.len() {
                if edges[ind_e].orientation == Orientation::Internal {
                    map_add(&mut internal_edges, &edges[ind_e]);
                    edges.remove(ind_e);
                } else {
                    ind_e += 1;
                }
            }

            // OCCT L3525-3533.
            if ref_face_orientation == Orientation::Reversed {
                for e in edges.iter_mut() {
                    e.orientation = occt_reverse(e.orientation);
                }
            }
            let mut f_ref_face = ref_face.clone();
            f_ref_face.orientation = Orientation::Forward;

            // OCCT L3535-4318: all faces collected — perform the union.
            if faces.len() > 1 {
                // OCCT L3538-3542.
                let mut edges_map = MapOfShape::new();
                let mut coord_tol =
                    compute_min_edge_size(brep, &edges, &f_ref_face, &mut edges_map);
                coord_tol /= 10.0;
                coord_tol = coord_tol.max(CONFUSION);

                // OCCT L3544-3551: the vertex -> edges map.
                let mut vemap = IndexedDataMapOfShapeListOfShape::new();
                for e in &edges {
                    map_shapes_and_unique_ancestors(
                        brep,
                        e,
                        ShapeType::Vertex,
                        ShapeType::Edge,
                        &mut vemap,
                    );
                }

                // OCCT L3553-3584: the seam search.
                let mut useam_found = false;
                let mut vseam_found = false;
                let mut edge_with_2pcurves = Shape::null();
                for f in &faces {
                    let face_ii = f.clone();
                    let an_outer_wire = outer_wire(brep, &face_ii);
                    for itw in iter_subshapes(brep, &an_outer_wire, true, true) {
                        let an_edge = itw;
                        if brep.is_edge_closed_on_face(&an_edge, &face_ii) {
                            if brep_tools_is_really_closed(&an_edge, &face_ii) {
                                if is_uiso(brep, &an_edge, &face_ii) {
                                    useam_found = true;
                                } else {
                                    vseam_found = true;
                                }
                            } else {
                                edge_with_2pcurves = an_edge;
                            }
                        }
                    }
                }
                let seam_found = useam_found || vseam_found;
                let _ = seam_found;

                // OCCT L3586-3596: the smoothness of the 2-pcurve edge.
                let mut a_is_edge_with_2pcurves_smooth = false;
                if self.my_concat_bsplines
                    && !edge_with_2pcurves.is_null()
                    && !(useam_found || vseam_found)
                {
                    let a_face_list = the_gmap_edge_faces
                        .get(&shape_key(&edge_with_2pcurves))
                        .map(|(_, l)| l.clone())
                        .unwrap_or_default();
                    if a_face_list.len() >= 2 {
                        let a_face1 = a_face_list[0].clone();
                        let a_face2 = a_face_list[a_face_list.len() - 1].clone();
                        let an_order_of_cont = brep_lib_continuity_of_faces(
                            brep,
                            &edge_with_2pcurves,
                            &a_face1,
                            &a_face2,
                            self.my_ang_tol,
                        );
                        use rcad_kernel::topods::GeomAbsShape as Gas;
                        a_is_edge_with_2pcurves_smooth = an_order_of_cont >= Gas::G1;
                    }
                }

                // OCCT L3598-3667: the periodic conversion of the reference.
                if a_is_edge_with_2pcurves_smooth {
                    // OCCT L3600-3616: the pcurve pair and the periodicity.
                    let mut a_p1: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut a_p2: Option<rcad_kernel::geom::Curve2d> = None;
                    if let Some((pc, a_first, a_last)) =
                        brep.curve_on_surface(&edge_with_2pcurves, &f_ref_face)
                    {
                        a_p1 = Some(pc);
                        let mut reversed = edge_with_2pcurves.clone();
                        reversed.orientation = occt_reverse(reversed.orientation);
                        a_p2 = brep
                            .curve_on_surface(&reversed, &f_ref_face)
                            .map(|(c, _, _)| c);
                        let _ = (a_first, a_last);
                    }
                    if let (Some(pc1), Some(pc2)) = (&a_p1, &a_p2) {
                        let a_pnt1 = pc1
                            .point_at(super::topexp::brep_tool_range(brep, &edge_with_2pcurves)[0]);
                        let a_pnt2 = pc2
                            .point_at(super::topexp::brep_tool_range(brep, &edge_with_2pcurves)[0]);
                        let an_is_uclosed =
                            (a_pnt1.x - a_pnt2.x).abs() > (a_pnt1.y - a_pnt2.y).abs();
                        let a_to_make_u_periodic = an_is_uclosed && uperiod == 0.0;
                        let a_to_make_v_periodic = !an_is_uclosed && vperiod == 0.0;

                        // OCCT L3618-3643.
                        if a_to_make_u_periodic || a_to_make_v_periodic {
                            let a_bspline_surface = match &a_base_surface {
                                Surface3::BSpline(b) => Some(b.clone()),
                                _ => {
                                    // OCCT L3624-3631: the approximation.
                                    Some(geom_convert_approx_surface(&a_base_surface))
                                }
                            };
                            if let Some(a_bs) = a_bspline_surface {
                                if a_to_make_u_periodic {
                                    // OCCT L3636-3637: SetUPeriodic.
                                    uperiod = geom_bspline_surface_period_span_u(&a_bs);
                                }
                                if a_to_make_v_periodic {
                                    // OCCT L3639-3640: SetVPeriodic.
                                    vperiod = geom_bspline_surface_period_span_v(&a_bs);
                                }

                                // OCCT L3646-3665: the reference rebuild.
                                if !surface_handle_same(
                                    &Surface3::BSpline(a_bs.clone()),
                                    &a_base_surface,
                                ) {
                                    let old_ref_face = ref_face.clone();
                                    ref_face = Shape::null();
                                    ref_face = a_bb.make_face(
                                        brep,
                                        Some(Surface3::BSpline(a_bs)),
                                        Shape::null(),
                                    );
                                    for e in &edges {
                                        let a_pcurve = brep
                                            .curve_on_surface(e, &old_ref_face)
                                            .map(|(c, _, _)| c);
                                        if map_edges_with_temporary_pcurves
                                            .contains_key(&shape_key(e))
                                        {
                                            builder_update_edge_pcurve_null(brep, e, &old_ref_face);
                                        }
                                        if let Some(pc) = a_pcurve {
                                            a_bb.update_edge_pcurve(
                                                brep,
                                                e.clone(),
                                                pc,
                                                ref_face.clone(),
                                                0.0,
                                            );
                                        }
                                    }
                                    f_ref_face = ref_face.clone();
                                    f_ref_face.orientation = Orientation::Forward;
                                }
                            }
                        }
                    }
                }

                // OCCT L3669-3803: the relocation to the new U-origin.
                let a_periods = [uperiod, vperiod];
                let mut an_is_seam_found = [useam_found, vseam_found];
                let a_surf_bounds = a_base_surface.default_domain();
                // OCCT L3673-3674: aSurfMin/aSurfMax from Bounds.
                let a_surf_min = [a_surf_bounds[0], a_surf_bounds[2]];
                let a_surf_max = [a_surf_bounds[1], a_surf_bounds[3]];
                let _ = a_surf_max;

                for ii in 0..2usize {
                    if a_periods[ii] != 0.0 {
                        // OCCT L3681.
                        if !an_is_seam_found[ii] {
                            // OCCT L3684-3699.
                            let mut a_min_coord = 0.0;
                            let mut a_max_coord = 0.0;
                            let mut a_number_of_intervals = 0i32;
                            let mut i_face_max = 0i32;
                            if !find_coord_bounds(
                                brep,
                                &faces,
                                &f_ref_face,
                                &a_map_ef,
                                &edges_map,
                                ii + 1,
                                a_periods[ii],
                                &mut a_min_coord,
                                &mut a_max_coord,
                                &mut a_number_of_intervals,
                                &mut i_face_max,
                            ) {
                                break;
                            }

                            // OCCT L3701-3704.
                            if a_max_coord - a_min_coord > a_periods[ii] - 1e-5 {
                                an_is_seam_found[ii] = true;
                            } else if a_number_of_intervals == 2 {
                                // OCCT L3705-3797.
                                let mut used_edges = MapOfShape::new();
                                let mut edge_new_pcurve = DataMapOfShapePCurve::new();

                                // OCCT L3710-3719.
                                if i_face_max >= 1 && (i_face_max as usize) <= faces.len() {
                                    relocate_pcurves_to_new_uorigin(
                                        brep,
                                        &edges,
                                        &faces[i_face_max as usize - 1],
                                        &f_ref_face,
                                        coord_tol,
                                        ii + 1,
                                        a_periods[ii],
                                        &mut vemap,
                                        &mut edge_new_pcurve,
                                        &mut used_edges,
                                    );
                                }

                                // OCCT L3721-3733: the unused edges.
                                for e in &edges {
                                    if !used_edges.contains_key(&shape_key(e)) {
                                        if let Some((pc, fpar, lpar)) =
                                            brep.curve_on_surface(e, &f_ref_face)
                                        {
                                            let pc = geom2d_trimmed(&pc, fpar, lpar);
                                            edge_new_pcurve.insert(shape_key(e), (e.clone(), pc));
                                        }
                                    }
                                }

                                // OCCT L3735-3740: restore the VEmap.
                                vemap.clear();
                                for e in &edges {
                                    map_shapes_and_unique_ancestors(
                                        brep,
                                        e,
                                        ShapeType::Vertex,
                                        ShapeType::Edge,
                                        &mut vemap,
                                    );
                                }

                                // OCCT L3742-3753: the new bounds.
                                let mut new_coord_min = super::statics_a::REAL_LAST;
                                let mut new_coord_max = super::statics_a::REAL_FIRST;
                                for e in &edges {
                                    if let Some((_, pc)) = edge_new_pcurve.get(&shape_key(e)) {
                                        super::statics_a::update_boundaries(
                                            pc,
                                            super::statics_b::first_param_of_2d(pc),
                                            super::statics_b::last_param_of_2d(pc),
                                            ii + 1,
                                            &mut new_coord_min,
                                            &mut new_coord_max,
                                        );
                                    }
                                }

                                // OCCT L3755-3797.
                                if new_coord_max - new_coord_min < a_periods[ii] - coord_tol
                                    && (-CONFUSION >= new_coord_min
                                        || new_coord_min >= a_periods[ii] + CONFUSION
                                        || -CONFUSION >= new_coord_max
                                        || new_coord_max >= a_periods[ii] + CONFUSION)
                                {
                                    // OCCT L3762-3782: the shifted surface.
                                    let rest_space_in_coord =
                                        a_periods[ii] - (new_coord_max - new_coord_min);
                                    let mut new_coord_origin =
                                        new_coord_min - rest_space_in_coord / 2.0;
                                    if new_coord_origin < a_surf_min[ii] {
                                        new_coord_origin = a_surf_min[ii];
                                    }
                                    let an_is_in_u = ii == 0;
                                    let new_surf = if new_coord_origin == a_surf_min[ii] {
                                        a_base_surface.clone()
                                    } else {
                                        let d = a_base_surface.default_domain();
                                        let (umin, umax, vmin, vmax) = (d[0], d[1], d[2], d[3]);
                                        let trim = if an_is_in_u {
                                            [
                                                new_coord_origin,
                                                new_coord_origin + a_periods[ii],
                                                vmin,
                                                vmax,
                                            ]
                                        } else {
                                            [
                                                umin,
                                                umax,
                                                new_coord_origin,
                                                new_coord_origin + a_periods[ii],
                                            ]
                                        };
                                        Surface3::Trimmed(TrimmedSurface {
                                            basis: Box::new(a_base_surface.clone()),
                                            trim,
                                        })
                                    };
                                    // OCCT L3783-3796: the reference rebuild.
                                    let old_ref_face = ref_face.clone();
                                    ref_face = Shape::null();
                                    ref_face = a_bb.make_face(brep, Some(new_surf), Shape::null());
                                    for e in &edges {
                                        if map_edges_with_temporary_pcurves
                                            .contains_key(&shape_key(e))
                                        {
                                            builder_update_edge_pcurve_null(brep, e, &old_ref_face);
                                        }
                                        if let Some((_, pc)) = edge_new_pcurve.get(&shape_key(e)) {
                                            a_bb.update_edge_pcurve(
                                                brep,
                                                e.clone(),
                                                pc.clone(),
                                                ref_face.clone(),
                                                0.0,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                useam_found = an_is_seam_found[0];
                vseam_found = an_is_seam_found[1];
                let _ = (useam_found, vseam_found);
                // OCCT L3804-3806.
                f_ref_face = ref_face.clone();
                f_ref_face.orientation = Orientation::Forward;

                let mut new_faces: Vec<Shape> = Vec::new();
                let mut new_wires: Vec<Shape> = Vec::new();

                // OCCT L3810-3829: the closed non-periodic periods.
                if uperiod == 0.0 || vperiod == 0.0 {
                    let (a_surf, _a_loc) = super::topexp::brep_tool_surface_loc(brep, &ref_face);
                    if let Some(a_surf) = a_surf {
                        let a_surf = super::statics_b::clear_rts(&a_surf);
                        let d = a_surf.default_domain();
                        let (ufirst, ulast, vfirst, vlast) = (d[0], d[1], d[2], d[3]);
                        if uperiod == 0.0 && a_surf.is_u_closed() {
                            uperiod = ulast - ufirst;
                        }
                        if vperiod == 0.0 && a_surf.is_v_closed() {
                            vperiod = vlast - vfirst;
                        }
                    }
                }

                let mut used_edges = MapOfShape::new();

                // OCCT L3833-3864: the minimal UV start.
                let mut face_umin = super::statics_a::REAL_LAST;
                let mut face_vmin = super::statics_a::REAL_LAST;
                for e in &edges {
                    let a_ba_curve = super::gap_deps::BRepAdaptorCurve2d::new(brep, e, &f_ref_face);
                    if a_ba_curve.curve().is_none() {
                        continue;
                    }
                    let a_first_point = a_ba_curve.value(a_ba_curve.first_parameter());
                    let a_last_point = a_ba_curve.value(a_ba_curve.last_parameter());
                    face_umin = face_umin.min(a_first_point.x).min(a_last_point.x);
                    face_vmin = face_vmin.min(a_first_point.y).min(a_last_point.y);
                }
                let _ = (face_umin, face_vmin);

                // OCCT L3866-4169: building the new wires and faces.
                while !edges.is_empty() {
                    // OCCT L3870-3877: the non-degenerated start edge.
                    let mut istart = 1usize;
                    let mut start_edge = edges[0].clone();
                    while is_edge_degenerated(brep, &start_edge) && istart < edges.len() {
                        istart += 1;
                        start_edge = edges[istart - 1].clone();
                    }

                    // OCCT L3879-3882.
                    let mut a_new_wire = a_bb.make_wire(brep);
                    a_bb.add_to_wire(brep, a_new_wire.clone(), start_edge.clone());
                    super::statics_a::remove_edge_from_map(brep, &start_edge, &mut vemap);
                    let mut splitting_vertices = MapOfShape::new();

                    // OCCT L3885-3892.
                    let start_pcurve = brep.curve_on_surface(&start_edge, &f_ref_face);
                    let Some((start_pc, fpar0, lpar0)) = start_pcurve.as_ref() else {
                        edges.remove(istart - 1);
                        continue;
                    };
                    let start_pc = start_pc.clone();
                    let (fpar0, lpar0) = (*fpar0, *lpar0);
                    let _ = (fpar0, lpar0);

                    // OCCT L3893-3905.
                    let (start_vertex, mut cur_vertex) = vertices(brep, &start_edge, true);
                    let (sp_first, sp_last) = (
                        start_pcurve.as_ref().unwrap().1,
                        start_pcurve.as_ref().unwrap().2,
                    );
                    let (start_param, cur_param0) =
                        if start_edge.orientation == Orientation::Forward {
                            (sp_first, sp_last)
                        } else {
                            (sp_last, sp_first)
                        };
                    let start_point = start_pc.point_at(start_param);
                    let mut cur_point = start_pc.point_at(cur_param0);
                    let _ = start_vertex;

                    // OCCT L3909-4116: the wire walk.
                    let mut cur_edge = start_edge.clone();
                    loop {
                        let mut next_edge = Shape::null();
                        let mut next_point = DVec2::ZERO;

                        // OCCT L3915-3921.
                        let elist: Vec<Shape> = vemap
                            .get(&shape_key(&cur_vertex))
                            .map(|(_, l)| l.clone())
                            .unwrap_or_default();
                        if elist.is_empty() && occt_is_same_shape(&cur_vertex, &start_vertex) {
                            // OCCT L3924-3947: the periodic seam branch.
                            if (uperiod != 0.0
                                && (start_point.x - cur_point.x).abs() > uperiod / 2.0)
                                || (vperiod != 0.0
                                    && (start_point.y - cur_point.y).abs() > vperiod / 2.0)
                            {
                                reconstruct_missed_seam_local(
                                    brep,
                                    &removed_edges,
                                    &f_ref_face,
                                    &cur_edge,
                                    &cur_vertex,
                                    cur_point,
                                    uperiod,
                                    vperiod,
                                    &mut next_edge,
                                    &mut next_point,
                                );
                            } else {
                                break; // end of wire
                            }
                        }

                        // OCCT L3950-4078.
                        if next_edge.is_null() {
                            let mut end_of_wire = false;

                            let an_is_on_singularity = is_on_singularity(brep, &elist);
                            if !an_is_on_singularity && elist.len() > 1 {
                                map_add(&mut splitting_vertices, &cur_vertex);
                            }

                            // OCCT L3960-3975.
                            let mut tmp_elist: Vec<Shape> = Vec::new();
                            for an_edge in &elist {
                                if used_edges.contains_key(&shape_key(an_edge)) {
                                    continue;
                                }
                                let a_first_vertex =
                                    super::topexp::first_vertex(brep, an_edge, true);
                                if !occt_is_same_shape(&a_first_vertex, &cur_vertex) {
                                    continue;
                                }
                                tmp_elist.push(an_edge.clone());
                            }
                            let mut true_elist: Vec<Shape> = Vec::new();
                            if tmp_elist.len() <= 1 || (uperiod != 0.0 || vperiod != 0.0) {
                                true_elist = tmp_elist;
                            } else {
                                // OCCT L3982-4017: the biggest-angle choice.
                                let mut max_angle = super::statics_a::REAL_FIRST;
                                let mut true_edge = Shape::null();
                                let cur_pc = brep.curve_on_surface(&cur_edge, &f_ref_face);
                                if let Some((cur_pc, fpar, lpar)) = cur_pc {
                                    let cur_param = if cur_edge.orientation == Orientation::Forward
                                    {
                                        lpar
                                    } else {
                                        fpar
                                    };
                                    let cur_point_v = cur_pc.point_at(cur_param);
                                    let cur_dir_v = cur_pc.derivative_at(cur_param);
                                    cur_point = cur_point_v;
                                    let mut cur_dir: DVec2 = cur_dir_v.normalize_or_zero();
                                    if cur_edge.orientation == Orientation::Reversed {
                                        cur_dir = -cur_dir;
                                    }
                                    for an_edge in &tmp_elist {
                                        let Some((a_pc, fpar, lpar)) =
                                            brep.curve_on_surface(an_edge, &f_ref_face)
                                        else {
                                            continue;
                                        };
                                        let a_param = if an_edge.orientation == Orientation::Forward
                                        {
                                            fpar
                                        } else {
                                            lpar
                                        };
                                        let a_dir_v = a_pc.derivative_at(a_param);
                                        let mut a_dir: DVec2 = a_dir_v.normalize_or_zero();
                                        if an_edge.orientation == Orientation::Reversed {
                                            a_dir = -a_dir;
                                        }
                                        let an_angle =
                                            super::gap_deps::dir_angle_2d(cur_dir, a_dir);
                                        if an_angle > max_angle {
                                            max_angle = an_angle;
                                            true_edge = an_edge.clone();
                                        }
                                    }
                                }
                                true_elist.push(true_edge);
                            }

                            // OCCT L4019-4072: the next-edge search.
                            for an_edge in &true_elist {
                                let Some((a_pc, fpar, lpar)) =
                                    brep.curve_on_surface(an_edge, &f_ref_face)
                                else {
                                    continue;
                                };
                                let a_param = if an_edge.orientation == Orientation::Forward {
                                    fpar
                                } else {
                                    lpar
                                };
                                let a_point = a_pc.point_at(a_param);
                                let diff_u = (a_point.x - cur_point.x).abs();
                                let diff_v = (a_point.y - cur_point.y).abs();
                                if uperiod != 0.0
                                    && diff_u > coord_tol
                                    && (diff_u - uperiod).abs() > coord_tol
                                {
                                    continue; // may be a deg.vertex
                                }
                                if vperiod != 0.0
                                    && diff_v > coord_tol
                                    && (diff_v - vperiod).abs() > coord_tol
                                {
                                    continue; // may be a deg.vertex
                                }

                                // OCCT L4039-4064: the periodic wrap check.
                                if (uperiod != 0.0 && diff_u > uperiod / 2.0)
                                    || (vperiod != 0.0 && diff_v > vperiod / 2.0)
                                {
                                    // OCCT literal quirk: StartOfNextEdge and
                                    // LastVertexOfSeam are declared but never
                                    // assigned in OCCT (cxx L4044-4045) — the
                                    // end check below reads the uninitialized
                                    // values; the rcad defaults keep that
                                    // verbatim.
                                    reconstruct_missed_seam_local(
                                        brep,
                                        &removed_edges,
                                        &f_ref_face,
                                        &cur_edge,
                                        &cur_vertex,
                                        cur_point,
                                        uperiod,
                                        vperiod,
                                        &mut next_edge,
                                        &mut next_point,
                                    );
                                    let start_of_next_edge = DVec2::ZERO;
                                    let last_vertex_of_seam = Shape::null();
                                    if occt_is_same_shape(&last_vertex_of_seam, &start_vertex)
                                        && (start_point.x - start_of_next_edge.x).abs()
                                            < uperiod / 2.0
                                    {
                                        end_of_wire = true;
                                    }
                                    break;
                                } else {
                                    // OCCT L4065-4071.
                                    next_edge = an_edge.clone();
                                    let last_param =
                                        if next_edge.orientation == Orientation::Forward {
                                            lpar
                                        } else {
                                            fpar
                                        };
                                    next_point = a_pc.point_at(last_param);
                                    break;
                                }
                            }

                            if end_of_wire {
                                break;
                            }
                        }

                        // OCCT L4080-4109.
                        if next_edge.is_null() {
                            if uperiod != 0.0 || vperiod != 0.0 {
                                if occt_is_same_shape(&cur_vertex, &start_vertex)
                                    && (uperiod == 0.0
                                        || (start_point.x - cur_point.x).abs() < uperiod / 2.0)
                                    && (vperiod == 0.0
                                        || (start_point.y - cur_point.y).abs() < vperiod / 2.0)
                                {
                                    break; // end of wire
                                }
                                reconstruct_missed_seam_local(
                                    brep,
                                    &removed_edges,
                                    &f_ref_face,
                                    &cur_edge,
                                    &cur_vertex,
                                    cur_point,
                                    uperiod,
                                    vperiod,
                                    &mut next_edge,
                                    &mut next_point,
                                );
                                if next_edge.is_null() {
                                    return;
                                }
                            } else {
                                return;
                            }
                        }

                        // OCCT L4110-4115.
                        cur_point = next_point;
                        cur_edge = next_edge;
                        cur_vertex = last_vertex(brep, &cur_edge, true);
                        a_bb.add_to_wire(brep, a_new_wire.clone(), cur_edge.clone());
                        map_add(&mut used_edges, &cur_edge);
                        super::statics_a::remove_edge_from_map(brep, &cur_edge, &mut vemap);
                    } // for (;;)

                    // OCCT L4118-4133.
                    set_flag_inplace(brep, &a_new_wire, tshape_flags::CLOSED, true);
                    map_add(&mut used_edges, &start_edge);
                    let mut ind = 0usize;
                    while ind < edges.len() {
                        if used_edges.contains_key(&shape_key(&edges[ind])) {
                            edges.remove(ind);
                        } else {
                            ind += 1;
                        }
                    }

                    // OCCT L4135-4168: face or hole classification.
                    let mut edge_on_bound_of_surf_found = false;
                    for itw in iter_subshapes(brep, &a_new_wire, true, true) {
                        if brep.is_edge_closed_on_face(&itw, &ref_face) {
                            edge_on_bound_of_surf_found = true;
                            break;
                        }
                    }
                    if edge_on_bound_of_surf_found {
                        // OCCT L4149-4156.
                        let (a_surf, a_loc) = super::topexp::brep_tool_surface_loc(brep, &ref_face);
                        if let Some(a_surf) = a_surf {
                            let mut a_result = a_bb.make_face(brep, Some(a_surf), Shape::null());
                            a_result.location = a_loc;
                            builder_add(brep, &a_result, &a_new_wire);
                            a_result.orientation = ref_face_orientation;
                            new_faces.push(a_result);
                        }
                    } else {
                        // OCCT L4157-4167: may be this wire is a hole.
                        if !splitting_vertices.is_empty() {
                            super::split_wire::split_wire(
                                brep,
                                &a_new_wire,
                                &f_ref_face,
                                &splitting_vertices,
                                &mut new_wires,
                            );
                        } else {
                            new_wires.push(a_new_wire);
                        }
                    }
                } // while (!edges.IsEmpty())

                // OCCT L4171-4228: the internal wires.
                let mut int_vemap = IndexedDataMapOfShapeListOfShape::new();
                for (_, ie) in internal_edges.iter() {
                    map_shapes_and_ancestors(
                        brep,
                        ie,
                        ShapeType::Vertex,
                        ShapeType::Edge,
                        &mut int_vemap,
                    );
                }
                let mut internal_wires: Vec<Shape> = Vec::new();
                while !internal_edges.is_empty() {
                    let first_key = *internal_edges.get_index(0).unwrap().0;
                    let a_first_edge = internal_edges.get_index(0).unwrap().1.clone();
                    internal_edges.shift_remove(&first_key);
                    super::statics_a::remove_edge_from_map(brep, &a_first_edge, &mut int_vemap);
                    let an_internal_wire = a_bb.make_wire(brep);
                    a_bb.add_to_wire(brep, an_internal_wire.clone(), a_first_edge.clone());
                    let mut end_edges = [a_first_edge.clone(), a_first_edge.clone()];
                    let (mut vv0, mut vv1) = vertices(brep, &a_first_edge, false);
                    loop {
                        if occt_is_same_shape(&vv0, &vv1) {
                            break; // closed wire
                        }
                        let mut found = false;
                        let vv = [vv0.clone(), vv1.clone()];
                        for ii in 0..2 {
                            let elist = int_vemap
                                .get(&shape_key(&vv[ii]))
                                .map(|(_, l)| l.clone())
                                .unwrap_or_default();
                            for an_edge in &elist {
                                if occt_is_same_shape(an_edge, &end_edges[ii]) {
                                    continue;
                                }
                                found = true;
                                let ek = shape_key(an_edge);
                                internal_edges.shift_remove(&ek);
                                super::statics_a::remove_edge_from_map(
                                    brep,
                                    an_edge,
                                    &mut int_vemap,
                                );
                                a_bb.add_to_wire(brep, an_internal_wire.clone(), an_edge.clone());
                                let (v1, v2) = vertices(brep, an_edge, false);
                                if occt_is_same_shape(&v1, &vv[ii]) {
                                    if ii == 0 {
                                        vv0 = v2;
                                    } else {
                                        vv1 = v2;
                                    }
                                } else if ii == 0 {
                                    vv0 = v1;
                                } else {
                                    vv1 = v1;
                                }
                                end_edges[ii] = an_edge.clone();
                                break;
                            }
                        }
                        if !found {
                            break; // end of open wire
                        }
                    }
                    internal_wires.push(an_internal_wire);
                }

                // OCCT L4230-4317: insert the new faces instead of the old.
                if new_faces.is_empty() {
                    // OCCT L4231-4252: one face without seam.
                    let (a_surf, a_loc) = super::topexp::brep_tool_surface_loc(brep, &ref_face);
                    if let Some(a_surf) = a_surf {
                        let mut a_result = a_bb.make_face(brep, Some(a_surf), Shape::null());
                        a_result.location = a_loc;
                        for w in &new_wires {
                            builder_add(brep, &a_result, w);
                        }
                        for w in &internal_wires {
                            builder_add(brep, &a_result, w);
                        }
                        a_result.orientation = ref_face_orientation;
                        self.my_context.merge(brep, &faces, &a_result);
                        for f in &faces {
                            self.my_face_new_face.insert(shape_key(f), a_result.clone());
                        }
                    }
                } else if new_faces.len() == 1 {
                    // OCCT L4254-4271.
                    let mut a_new_face = new_faces[0].clone();
                    a_new_face.orientation = Orientation::Forward;
                    for w in &new_wires {
                        builder_add(brep, &a_new_face, w);
                    }
                    for w in &internal_wires {
                        builder_add(brep, &a_new_face, w);
                    }
                    self.my_context.merge(brep, &faces, &new_faces[0]);
                    for f in &faces {
                        self.my_face_new_face
                            .insert(shape_key(f), new_faces[0].clone());
                    }
                } else {
                    // OCCT L4272-4317: the split result.
                    insert_wires_into_faces(brep, &new_wires, &new_faces, &ref_face);
                    insert_wires_into_faces(brep, &internal_wires, &new_faces, &ref_face);

                    let mut emaps: Vec<MapOfShape> = Vec::new();
                    for f in &faces {
                        let mut a_emap = MapOfShape::new();
                        map_shapes_all(brep, f, &mut a_emap);
                        emaps.push(a_emap);
                    }
                    for (ii, nf) in new_faces.iter().enumerate() {
                        let mut faces_for_this_face: Vec<Shape> = Vec::new();
                        let mut used_faces = MapOfShape::new();
                        for an_edge in topexp_explorer(brep, nf, ShapeType::Edge) {
                            if is_edge_degenerated(brep, &an_edge)
                                || brep.is_edge_closed_on_face(&an_edge, &ref_face)
                            {
                                continue;
                            }
                            let mut jj = 0usize;
                            while jj < emaps.len() && !emaps[jj].contains_key(&shape_key(&an_edge))
                            {
                                jj += 1;
                            }
                            if jj >= emaps.len() {
                                continue;
                            }
                            if map_add(&mut used_faces, &faces[jj]) {
                                faces_for_this_face.push(faces[jj].clone());
                            }
                        }
                        self.my_context.merge(brep, &faces_for_this_face, nf);
                        for f in &faces_for_this_face {
                            self.my_face_new_face
                                .insert(shape_key(f), new_faces[ii].clone());
                        }
                    }
                }
            } // if (faces.Length() > 1)
        } // end processing each face
    }
}

/// The module-local forward of the OCCT static ReconstructMissedSeam (the
/// full translation lives in statics_b.rs).
#[allow(clippy::too_many_arguments)]
fn reconstruct_missed_seam_local(
    brep: &mut BRep,
    removed_edges: &[Shape],
    f_ref_face: &Shape,
    cur_edge: &Shape,
    cur_vertex: &Shape,
    cur_point: DVec2,
    uperiod: f64,
    vperiod: f64,
    next_edge: &mut Shape,
    next_point: &mut DVec2,
) {
    super::statics_b::reconstruct_missed_seam(
        brep,
        removed_edges,
        f_ref_face,
        cur_edge,
        cur_vertex,
        cur_point,
        uperiod,
        vperiod,
        next_edge,
        next_point,
    );
}

