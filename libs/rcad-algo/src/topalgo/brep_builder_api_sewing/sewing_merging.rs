//! OCCT BRepBuilderAPI_Sewing.cxx — the merging group:
//! `replaceNMVertices` static (L3608-3705), `ReplaceEdge` static
//! (L3707-3787), `Merging` (L3789-4241), `GetSeqEdges` static (L4604-4652),
//! `GetFreeWires` (L4654-4707), `IsDegeneratedWire` static (L4709-4798),
//! `DegeneratedSection` static (L4800-4910), `EdgeProcessing`
//! (L4912-5001) and `EdgeRegularity` (L5003-5028).
//!
//! GAP note (BRepLib::EncodeRegularity(E, F1, F2, Tol)): the OCCT body
//! delegates to `ContinuityOfFaces` (BRepLib.cxx L2385-2565, the G1/G2
//! face-continuity classification) which is untranslated; the local re-host
//! keeps the call form and is a documented no-op, the same policy as the
//! crate `topalgo::brep_lib::BRepLib::same_parameter` stub.

use rcad_kernel::core::precision::PCONFUSION;
use rcad_kernel::geom::{Curve2d, CurveEval};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

use crate::brep_algo::tool as bat;

use super::{idx_list_get, iter_no_cumori, set_add, set_contains};
use super::BRepBuilderAPISewing;

/// OCCT BRepLib::EncodeRegularity(E, F1, F2, TolAng) (BRepLib.cxx L2570-2587)
/// — the local no-op re-host (see the module GAP note).
fn brep_lib_encode_regularity_edge(
    _brep: &mut BRep,
    _e: &Shape,
    _f1: &Shape,
    _f2: &Shape,
    _tol: f64,
) {
}

/// OCCT static replaceNMVertices(theEdge, theV1, theV2, theReShape)
/// (cxx L3608-3705).
fn replace_nm_vertices(
    brep: &mut BRep,
    the_edge: &Shape,
    the_v1: &Shape,
    the_v2: &Shape,
    the_re_shape: &mut crate::shhealing::shape_build::reshape::ShapeBuildReShape,
) {
    // To keep NM vertices on edge
    // OCCT L3612-3614.
    let mut a_seq_nm_vert: Vec<Shape> = Vec::new();
    let mut a_seq_nm_pars: Vec<f64> = Vec::new();
    let has_nm_vert =
        super::same_param::find_nm_vertices(brep, the_edge, &mut a_seq_nm_vert, &mut a_seq_nm_pars);
    // OCCT L3615-3618.
    if !has_nm_vert {
        return;
    }
    // OCCT L3619-3626.
    let (first, last) = bat::brep_tool_range(the_edge);
    let (_c3d, _a_loc, _f, _l) = match super::brep_tool_curve_loc(the_edge) {
        Some(v) => v,
        None => return,
    };
    let _ = (first, last);
    let mut a_ed_vert: Vec<Shape> = Vec::new();
    let mut a_ed_params: Vec<f64> = Vec::new();
    let nb = a_seq_nm_pars.len();

    // OCCT L3630-3670.
    for i in 1..=nb {
        let apar = a_seq_nm_pars[i - 1];
        if (apar - first).abs() <= PCONFUSION {
            // OCCT L3633-3636.
            let v = a_seq_nm_vert[i - 1].clone();
            the_re_shape.replace(brep, &v, the_v1);
            continue;
        }
        if (apar - last).abs() <= PCONFUSION {
            // OCCT L3637-3640.
            let v = a_seq_nm_vert[i - 1].clone();
            the_re_shape.replace(brep, &v, the_v2);
            continue;
        }
        let a_v = a_seq_nm_vert[i - 1].clone();
        // OCCT L3642-3660.
        let mut j = 1usize;
        let mut inserted = false;
        while j <= a_ed_params.len() {
            let apar2 = a_ed_params[j - 1];
            if (apar - apar2).abs() <= PCONFUSION {
                let v = a_ed_vert[j - 1].clone();
                the_re_shape.replace(brep, &a_v, &v);
                inserted = true;
                break;
            } else if apar < apar2 {
                let mut anew_v = brep.empty_copied(&a_v);
                a_ed_vert.insert(j - 1, anew_v.clone());
                a_ed_params.insert(j - 1, apar);
                // OCCT L3651-3656: append BRep_PointOnCurve(apar, c3d, aLoc)
                // to the new vertex's point representations.  Architecture
                // gap: the rcad PointRepresentation stores the curve by
                // pool index (no Curve3-valued form), so the parameter is
                // recorded with the sentinel curve index.
                if let TShape::Vertex(vd) = Arc::make_mut(&mut anew_v.data) {
                    vd.points.push(rcad_kernel::topo::topods::PointRepresentation::PointOnCurve {
                        curve: usize::MAX,
                        parameter: apar,
                        tolerance: 0.0,
                    });
                }
                the_re_shape.replace(brep, &a_v, &anew_v);
                inserted = true;
                break;
            } else {
                j += 1;
            }
        }
        // OCCT L3661-3675.
        if !inserted && j > a_ed_params.len() {
            let mut anew_v = brep.empty_copied(&a_v);
            a_ed_vert.push(anew_v.clone());
            a_ed_params.push(apar);
            if let TShape::Vertex(vd) = Arc::make_mut(&mut anew_v.data) {
                vd.points.push(rcad_kernel::topo::topods::PointRepresentation::PointOnCurve {
                    curve: usize::MAX,
                    parameter: apar,
                    tolerance: 0.0,
                });
            }
            the_re_shape.replace(brep, &a_v, &anew_v);
        }
    }

    // OCCT L3678-3704.
    let newnb = a_ed_params.len();
    if newnb < nb {
        let mut anew_edge = brep.empty_copied(&the_edge_alias(the_edge));
        let an_ori = the_edge.orientation;
        anew_edge.orientation = Orientation::Forward;
        let mut a_b = BRepBuilder::new();
        a_b.add_to_edge(brep, anew_edge.clone(), the_v1.clone());
        a_b.add_to_edge(brep, anew_edge.clone(), the_v2.clone());

        for i in 1..=a_ed_vert.len() {
            a_b.add_to_edge(brep, anew_edge.clone(), a_ed_vert[i - 1].clone());
        }
        anew_edge.orientation = an_ori;
        the_re_shape.replace(brep, the_edge, &anew_edge);
    }
}

/// TopoDS_Edge alias helper (the OCCT code passes the same handle; rcad
/// clones the handle).
fn the_edge_alias(e: &Shape) -> Shape {
    e.clone()
}

/// OCCT static ReplaceEdge(oldEdge, theNewShape, aReShape) (cxx L3707-3787).
fn replace_edge(
    brep: &mut BRep,
    old_edge: &Shape,
    the_new_shape: &Shape,
    a_re_shape: &mut crate::shhealing::shape_build::reshape::ShapeBuildReShape,
) {
    // OCCT L3711-3712.
    let old_shape = a_re_shape.apply(brep, old_edge, ShapeType::Shape);
    let new_shape = a_re_shape.apply(brep, the_new_shape, ShapeType::Shape);
    // OCCT L3713-3716.
    if old_shape.is_same(&new_shape) || a_re_shape.is_recorded(&new_shape) {
        return;
    }

    // OCCT L3718.
    a_re_shape.replace(brep, &old_shape, &new_shape);
    // OCCT L3719-3723.
    let (v1old, v2old) = bat::top_exp_vertices_raw(&old_shape);
    let v1old = v1old.unwrap_or_else(Shape::null);
    let v2old = v2old.unwrap_or_else(Shape::null);
    let or_old = old_shape.orientation;
    let mut or_new = or_old;
    let (mut v1new, mut v2new) = (Shape::null(), Shape::null());
    if new_shape.shape_type() == ShapeType::Edge {
        // OCCT L3725-3731.
        let a_en = new_shape.clone();
        let (a_v1, a_v2) = bat::top_exp_vertices_raw(&a_en);
        v1new = a_v1.unwrap_or_else(Shape::null);
        v2new = a_v2.unwrap_or_else(Shape::null);
        or_new = a_en.orientation;
        replace_nm_vertices(brep, &a_en, &v1new, &v2new, a_re_shape);
    } else if new_shape.shape_type() == ShapeType::Wire {
        // OCCT L3732-3747.
        for aex in bat::explorer(&new_shape, ShapeType::Edge, ShapeType::Shape) {
            let ed = aex;
            or_new = ed.orientation;
            let (a_v1, a_v2) = bat::top_exp_vertices_raw(&ed);
            let a_v1 = a_v1.unwrap_or_else(Shape::null);
            let a_v2 = a_v2.unwrap_or_else(Shape::null);
            replace_nm_vertices(brep, &ed, &a_v1, &a_v2, a_re_shape);
            if v1new.is_null() {
                v1new = a_v1;
            }
            v2new = a_v2;
        }
    }

    // OCCT L3750-3786.
    v1new.orientation = v1old.orientation;
    v2new.orientation = v2old.orientation;
    if v1old.is_same(&v2old) && !v1old.is_same(&v1new) && !a_re_shape.is_recorded(&v1new) {
        a_re_shape.replace(brep, &v1old, &v1new);
        return;
    }
    if or_old == or_new {
        v1new.orientation = v1old.orientation;
        v2new.orientation = v2old.orientation;
        if !v1old.is_same(&v1new) && !v1old.is_same(&v2new) && !a_re_shape.is_recorded(&v1new) {
            a_re_shape.replace(brep, &v1old, &v1new);
        }
        if !v2old.is_same(&v2new) && !v2old.is_same(&v1new) && !a_re_shape.is_recorded(&v2new) {
            a_re_shape.replace(brep, &v2old, &v2new);
        }
    } else {
        v1new.orientation = v2old.orientation;
        v2new.orientation = v1old.orientation;
        if !v1old.is_same(&v2new) && !v1old.is_same(&v1new) && !a_re_shape.is_recorded(&v2new) {
            a_re_shape.replace(brep, &v1old, &v2new);
        }
        if !v2old.is_same(&v2new) && !v2old.is_same(&v1new) && !a_re_shape.is_recorded(&v1new) {
            a_re_shape.replace(brep, &v2old, &v1new);
        }
    }
}

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::Merging(firstTime, theProgress) (cxx
    /// L3789-4241) — merges the free boundaries.
    pub(crate) fn merging(&mut self, brep: &mut BRep, _first_time: bool) {
        // OCCT L3792-3796: the myBoundFaces walk.
        let nb_bounds = self.my_bound_faces.len();
        for idx in 0..nb_bounds {
            let (bound, list_faces) = idx_list_get(&self.my_bound_faces, idx);

            // If bound was already merged - continue
            // OCCT L3803-3806.
            if set_contains(&self.my_merged_edges, &bound) {
                continue;
            }

            // OCCT L3808-3828: merge free edge - only vertices.
            if list_faces.is_empty() {
                let (no1, no2) = bat::top_exp_vertices_raw(&bound);
                let no1 = no1.unwrap_or_else(Shape::null);
                let no2 = no2.unwrap_or_else(Shape::null);
                let mut nno1 = no1.clone();
                let mut nno2 = no2.clone();
                if let Some((_, n)) = self.my_vertex_node_free.get(&bat::shape_key(&no1)) {
                    nno1 = n.clone();
                }
                if let Some((_, n)) = self.my_vertex_node_free.get(&bat::shape_key(&no2)) {
                    nno2 = n.clone();
                }
                if !no1.is_same(&nno1) {
                    let o = bat::oriented(&nno1, no1.orientation);
                    self.my_re_shape.replace(brep, &no1, &o);
                }
                if !no2.is_same(&nno2) {
                    let o = bat::oriented(&nno2, no2.orientation);
                    self.my_re_shape.replace(brep, &no2, &o);
                }
                set_add(&mut self.my_merged_edges, &bound);
                continue;
            }

            // Check for previous splitting, build replacing wire
            // OCCT L3831-3849.
            let mut bound_wire: Option<Shape> = None;
            let mut is_prev_split = false;
            let has_cutting_sections = self.my_bound_sections.contains_key(&bat::shape_key(&bound));
            if has_cutting_sections {
                let mut b = BRepBuilder::new();
                let mut wire = b.make_wire(brep);
                wire.orientation = bound.orientation;
                // Iterate on cutting sections
                let sections = self
                    .my_bound_sections
                    .get(&bat::shape_key(&bound))
                    .cloned()
                    .unwrap_or_default();
                for its in &sections {
                    let section = its.clone();
                    b.add_to_wire(brep, wire.clone(), section.clone());
                    if set_contains(&self.my_merged_edges, &section) {
                        is_prev_split = true;
                    }
                }
                bound_wire = Some(wire);
            }

            // Merge with bound
            // OCCT L3852-3853.
            let mut merged_with_bound: super::IdxShapeMap = super::IdxShapeMap::new();
            if !is_prev_split {
                // Obtain sequence of edges merged with bound
                // OCCT L3856-3858.
                let mut seq_merged_with_bound: Vec<Shape> = Vec::new();
                let mut seq_merged_with_bound_ori: Vec<bool> = Vec::new();
                if self.merged_nearest_edges(
                    brep,
                    &bound,
                    &mut seq_merged_with_bound,
                    &mut seq_merged_with_bound_ori,
                ) {
                    // Store bound in the map
                    // OCCT L3862.
                    super::idx_map_add(&mut merged_with_bound, &bound, bound.clone());
                    // Iterate on edges merged with bound
                    // OCCT L3865-3906.
                    let mut ii = 1usize;
                    while ii <= seq_merged_with_bound.len() {
                        let iedge = seq_merged_with_bound[ii - 1].clone();
                        // Remove edge if recorded as merged
                        let mut is_rejected = set_contains(&self.my_merged_edges, &iedge)
                            || merged_with_bound.contains_key(&bat::shape_key(&iedge));
                        if !is_rejected {
                            if let Some(sections) =
                                self.my_bound_sections.get(&bat::shape_key(&iedge)).cloned()
                            {
                                // Edge is split - check sections
                                for lit in &sections {
                                    if is_rejected {
                                        break;
                                    }
                                    let sec = lit.clone();
                                    is_rejected = set_contains(&self.my_merged_edges, &sec)
                                        || merged_with_bound.contains_key(&bat::shape_key(&sec));
                                }
                            }
                            if !is_rejected {
                                if let Some(bnd_v) =
                                    self.my_section_bound.get(&bat::shape_key(&iedge)).cloned()
                                {
                                    // Edge is a section - check bound
                                    let bnd = bnd_v;
                                    is_rejected = set_contains(&self.my_merged_edges, &bnd)
                                        || merged_with_bound.contains_key(&bat::shape_key(&bnd));
                                }
                            }
                        }
                        // To the next merged edge
                        if is_rejected {
                            // Remove rejected edge
                            seq_merged_with_bound.remove(ii - 1);
                            seq_merged_with_bound_ori.remove(ii - 1);
                        } else {
                            // Process accepted edge
                            super::idx_map_add(&mut merged_with_bound, &iedge, iedge.clone());
                            ii += 1;
                        }
                    }
                    // OCCT L3907-3991.
                    let mut nb_merged = seq_merged_with_bound.len();
                    if nb_merged != 0 {
                        // Create same parameter edge
                        let mut actually_merged: super::ShapeSet = super::ShapeSet::new();
                        let mut my_re_shape = std::mem::take(&mut self.my_re_shape);
                        let merged_edge = self.same_parameter_edge_seq(
                            brep,
                            &bound,
                            &seq_merged_with_bound,
                            &seq_merged_with_bound_ori,
                            &mut actually_merged,
                            &mut my_re_shape,
                        );
                        self.my_re_shape = my_re_shape;
                        let mut is_forward = false;
                        if !merged_edge.is_null() {
                            is_forward = merged_edge.orientation == Orientation::Forward;
                        }
                        // Process actually merged edges
                        // OCCT L3916-3941.
                        let mut nb_actually_merged = 0usize;
                        for ii in 1..=nb_merged {
                            let iedge = seq_merged_with_bound[ii - 1].clone();
                            if set_contains(&actually_merged, &iedge) {
                                nb_actually_merged += 1;
                                // Record merged edge in the map
                                let mut orient = iedge.orientation;
                                if !is_forward {
                                    orient = super::top_abs_reverse(orient);
                                }
                                if !seq_merged_with_bound_ori[ii - 1] {
                                    orient = super::top_abs_reverse(orient);
                                }
                                if let Some(entry) =
                                    merged_with_bound.get_mut(&bat::shape_key(&iedge))
                                {
                                    entry.1 = bat::oriented(&merged_edge, orient);
                                }
                            } else {
                                merged_with_bound.shift_remove(&bat::shape_key(&iedge));
                            }
                        }
                        // OCCT L3942-3956.
                        if nb_actually_merged != 0 {
                            // Record merged bound in the map
                            let mut orient = bound.orientation;
                            if !is_forward {
                                orient = super::top_abs_reverse(orient);
                            }
                            if let Some(entry) =
                                merged_with_bound.get_mut(&bat::shape_key(&bound))
                            {
                                entry.1 = bat::oriented(&merged_edge, orient);
                            }
                        }
                        nb_merged = nb_actually_merged;
                    }
                    // Remove bound from the map if not finally merged
                    // OCCT L3958-3962.
                    if nb_merged == 0 {
                        merged_with_bound.shift_remove(&bat::shape_key(&bound));
                    }
                }
            }
            // OCCT L3964: const bool isMerged = !MergedWithBound.IsEmpty();
            let is_merged = !merged_with_bound.is_empty();

            // Merge with cutting sections
            // OCCT L3967-4128.
            let mut sections_re_shape = crate::shhealing::shape_build::reshape::ShapeBuildReShape::new();
            let mut merged_with_sections: super::IdxShapeMap = super::IdxShapeMap::new();
            if has_cutting_sections {
                // Iterate on cutting sections
                let sections = self
                    .my_bound_sections
                    .get(&bat::shape_key(&bound))
                    .cloned()
                    .unwrap_or_default();
                for its in &sections {
                    // Retrieve cutting section
                    let section = its.clone();
                    // Skip section if already merged
                    if set_contains(&self.my_merged_edges, &section) {
                        continue;
                    }
                    // Merge cutting section
                    // OCCT L3991-4105.
                    let mut seq_merged_with_section: Vec<Shape> = Vec::new();
                    let mut seq_merged_with_section_ori: Vec<bool> = Vec::new();
                    if self.merged_nearest_edges(
                        brep,
                        &section,
                        &mut seq_merged_with_section,
                        &mut seq_merged_with_section_ori,
                    ) {
                        // Store section in the map
                        super::idx_map_add(&mut merged_with_sections, &section, section.clone());
                        // Iterate on edges merged with section
                        let mut ii = 1usize;
                        while ii <= seq_merged_with_section.len() {
                            let iedge = seq_merged_with_section[ii - 1].clone();
                            // Remove edge if recorded as merged
                            let mut is_rejected = set_contains(&self.my_merged_edges, &iedge)
                                || merged_with_sections.contains_key(&bat::shape_key(&iedge));
                            if !is_rejected {
                                if let Some(sections2) =
                                    self.my_bound_sections.get(&bat::shape_key(&iedge)).cloned()
                                {
                                    for lit in &sections2 {
                                        if is_rejected {
                                            break;
                                        }
                                        let sec = lit.clone();
                                        is_rejected = set_contains(&self.my_merged_edges, &sec)
                                            || merged_with_sections
                                                .contains_key(&bat::shape_key(&sec));
                                    }
                                }
                                if !is_rejected {
                                    if let Some(bnd_v) =
                                        self.my_section_bound.get(&bat::shape_key(&iedge)).cloned()
                                    {
                                        let bnd = bnd_v;
                                        is_rejected = set_contains(&self.my_merged_edges, &bnd)
                                            || merged_with_sections
                                                .contains_key(&bat::shape_key(&bnd));
                                    }
                                }
                            }
                            // To the next merged edge
                            if is_rejected {
                                seq_merged_with_section.remove(ii - 1);
                                seq_merged_with_section_ori.remove(ii - 1);
                            } else {
                                // Process accepted edge
                                super::idx_map_add(
                                    &mut merged_with_sections,
                                    &iedge,
                                    iedge.clone(),
                                );
                                ii += 1;
                            }
                        }
                        let mut nb_merged = seq_merged_with_section.len();
                        if nb_merged != 0 {
                            // Create same parameter edge
                            let mut actually_merged: super::ShapeSet = super::ShapeSet::new();
                            let merged_edge = self.same_parameter_edge_seq(
                                brep,
                                &section,
                                &seq_merged_with_section,
                                &seq_merged_with_section_ori,
                                &mut actually_merged,
                                &mut sections_re_shape,
                            );
                            let mut is_forward = false;
                            if !merged_edge.is_null() {
                                is_forward = merged_edge.orientation == Orientation::Forward;
                            }
                            // Process actually merged edges
                            let mut nb_actually_merged = 0usize;
                            for ii in 1..=nb_merged {
                                let iedge = seq_merged_with_section[ii - 1].clone();
                                if set_contains(&actually_merged, &iedge) {
                                    nb_actually_merged += 1;
                                    // Record merged edge in the map
                                    let mut orient = iedge.orientation;
                                    if !is_forward {
                                        orient = super::top_abs_reverse(orient);
                                    }
                                    if !seq_merged_with_section_ori[ii - 1] {
                                        orient = super::top_abs_reverse(orient);
                                    }
                                    let oedge = bat::oriented(&merged_edge, orient);
                                    if let Some(entry) =
                                        merged_with_sections.get_mut(&bat::shape_key(&iedge))
                                    {
                                        entry.1 = oedge.clone();
                                    }
                                    // OCCT L4053.
                                    let applied =
                                        self.my_re_shape.apply(brep, &iedge, ShapeType::Shape);
                                    replace_edge(brep, &applied, &oedge, &mut sections_re_shape);
                                } else {
                                    merged_with_sections.shift_remove(&bat::shape_key(&iedge));
                                }
                            }
                            // OCCT L4067-4081.
                            if nb_actually_merged != 0 {
                                // Record merged section in the map
                                let mut orient = section.orientation;
                                if !is_forward {
                                    orient = super::top_abs_reverse(orient);
                                }
                                let oedge = bat::oriented(&merged_edge, orient);
                                if let Some(entry) =
                                    merged_with_sections.get_mut(&bat::shape_key(&section))
                                {
                                    entry.1 = oedge.clone();
                                }
                                let applied =
                                    self.my_re_shape.apply(brep, &section, ShapeType::Shape);
                                replace_edge(brep, &applied, &oedge, &mut sections_re_shape);
                            }
                            nb_merged = nb_actually_merged;
                        }
                        // Remove section from the map if not finally merged
                        // OCCT L4084-4088.
                        if nb_merged == 0 {
                            merged_with_sections.shift_remove(&bat::shape_key(&section));
                        }
                    } else if is_merged {
                        // Reject merging of sections
                        // OCCT L4090-4094.
                        merged_with_sections.clear();
                        break;
                    }
                }
            }
            // OCCT L4129: const bool isMergedSplit = !MergedWithSections.IsEmpty();
            let is_merged_split = !merged_with_sections.is_empty();

            // OCCT L4131-4143.
            if !is_merged && !is_merged_split {
                // Nothing was merged in this iteration
                if is_prev_split {
                    // Replace previously split bound
                    let applied_bound =
                        self.my_re_shape.apply(brep, &bound, ShapeType::Shape);
                    let applied_wire =
                        self.my_re_shape.apply(brep, &bound_wire.clone().unwrap(), ShapeType::Shape);
                    self.my_re_shape.replace(brep, &applied_bound, &applied_wire);
                }
                continue;
            }

            // Set splitting flag
            // OCCT L4146.
            let mut is_splitted = (!is_merged && is_merged_split) || is_prev_split;

            // Choose between bound and sections merging
            // OCCT L4149-4193.
            if is_merged && is_merged_split && !is_prev_split {
                // Fill map of merged cutting sections
                let mut map_split_edges: super::ShapeSet = super::ShapeSet::new();
                for i in 0..merged_with_sections.len() {
                    let (edge, _) = super::idx_shape_get(&merged_with_sections, i);
                    set_add(&mut map_split_edges, &edge);
                }
                // Iterate on edges merged with bound
                for i in 0..merged_with_bound.len() {
                    // Retrieve edge merged with bound
                    let (edge, _) = super::idx_shape_get(&merged_with_bound, i);
                    // Remove edge from the map
                    if set_contains(&map_split_edges, &edge) {
                        map_split_edges.shift_remove(&bat::shape_key(&edge));
                    }
                    if let Some(sections) =
                        self.my_bound_sections.get(&bat::shape_key(&edge)).cloned()
                    {
                        // Edge has cutting sections
                        for its in &sections {
                            let sec = its.clone();
                            // Remove section from the map
                            if set_contains(&map_split_edges, &sec) {
                                map_split_edges.shift_remove(&bat::shape_key(&sec));
                            }
                        }
                    }
                }
                // Calculate section merging tolerance
                // OCCT L4174-4180.
                let mut min_split_tol = f64::MAX;
                for ii in 0..map_split_edges.len() {
                    let split_key = super::set_get(&map_split_edges, ii);
                    if let Some(entry) = merged_with_sections.get(&bat::shape_key(&split_key)) {
                        let edge = entry.1.clone();
                        min_split_tol = min_split_tol.min(bat::brep_tool_tolerance(&edge));
                    }
                }
                // Calculate bound merging tolerance
                // OCCT L4182-4188.
                let bound_edge = merged_with_bound
                    .get(&bat::shape_key(&bound))
                    .map(|(_, v)| v.clone())
                    .unwrap_or_else(Shape::null);
                let bound_edge_tol = bat::brep_tool_tolerance(&bound_edge);
                is_splitted =
                    (min_split_tol < bound_edge_tol + self.min_tolerance()) || self.my_nonmanifold;
                is_splitted = !map_split_edges.is_empty() && is_splitted;
            }

            // OCCT L4196-4236.
            if is_splitted {
                // Merging of cutting sections
                let applied_bound = self.my_re_shape.apply(brep, &bound, ShapeType::Shape);
                let applied_wire =
                    self.my_re_shape.apply(brep, &bound_wire.clone().unwrap(), ShapeType::Shape);
                self.my_re_shape.replace(brep, &applied_bound, &applied_wire);
                for i in 0..merged_with_sections.len() {
                    let (oldedge, newedge_v) = super::idx_shape_get(&merged_with_sections, i);
                    let newedge = sections_re_shape.apply(brep, &newedge_v, ShapeType::Shape);
                    let applied = self.my_re_shape.apply(brep, &oldedge, ShapeType::Shape);
                    replace_edge(brep, &applied, &newedge, &mut self.my_re_shape);
                    set_add(&mut self.my_merged_edges, &oldedge);
                    if self.my_bound_sections.contains_key(&bat::shape_key(&oldedge)) {
                        self.my_bound_sections.remove(&bat::shape_key(&oldedge));
                    }
                }
            } else {
                // Merging of initial bound
                for i in 0..merged_with_bound.len() {
                    let (oldedge, newedge) = super::idx_shape_get(&merged_with_bound, i);
                    let applied = self.my_re_shape.apply(brep, &oldedge, ShapeType::Shape);
                    replace_edge(brep, &applied, &newedge, &mut self.my_re_shape);
                    set_add(&mut self.my_merged_edges, &oldedge);
                    if self.my_bound_sections.contains_key(&bat::shape_key(&oldedge)) {
                        self.my_bound_sections.remove(&bat::shape_key(&oldedge));
                    }
                }
                if self.my_bound_sections.contains_key(&bat::shape_key(&bound)) {
                    self.my_bound_sections.remove(&bat::shape_key(&bound));
                }
                if !set_contains(&self.my_merged_edges, &bound) {
                    set_add(&mut self.my_merged_edges, &bound);
                }
            }
        }

        // OCCT L4238-4241.
        self.my_nb_vertices =
            (self.my_vertex_node.len() + self.my_vertex_node_free.len()) as i32;
        self.my_node_sections.clear();
        self.my_vertex_node.clear();
        self.my_vertex_node_free.clear();
        self.my_cutting_node.clear();
    }

    /// OCCT BRepBuilderAPI_Sewing::EdgeProcessing(theProgress) (cxx
    /// L4912-5001).
    pub(crate) fn edge_processing(&mut self, brep: &mut BRep) {
        // constructs sectionEdge
        // OCCT L4915-4916.
        let mut map_free_edges: super::ShapeSet = super::ShapeSet::new();
        let mut edge_face: super::DataShapeMap = super::DataShapeMap::new();
        let nb_bounds = self.my_bound_faces.len();
        for idx in 0..nb_bounds {
            let (bound, list_faces) = idx_list_get(&self.my_bound_faces, idx);
            if list_faces.len() == 1 {
                if let Some(sections) =
                    self.my_bound_sections.get(&bat::shape_key(&bound)).cloned()
                {
                    // OCCT L4928-4941.
                    for liter in &sections {
                        if !set_contains(&self.my_merged_edges, liter) {
                            let edge =
                                self.my_re_shape.apply(brep, liter, ShapeType::Shape);
                            if !set_contains(&map_free_edges, &edge) {
                                let face = list_faces[0].clone();
                                edge_face.insert(bat::shape_key(&edge), face);
                                set_add(&mut map_free_edges, &edge);
                            }
                        }
                    }
                } else if !set_contains(&self.my_merged_edges, &bound) {
                    // OCCT L4942-4954.
                    let edge = self.my_re_shape.apply(brep, &bound, ShapeType::Shape);
                    if !set_contains(&map_free_edges, &edge) {
                        let face = list_faces[0].clone();
                        edge_face.insert(bat::shape_key(&edge), face);
                        set_add(&mut map_free_edges, &edge);
                    }
                }
            }
        }

        // OCCT L4957-4998.
        if !map_free_edges.is_empty() {
            let mut seq_wires: Vec<Shape> = Vec::new();
            get_free_wires(brep, &mut map_free_edges, &mut seq_wires);
            for j in 1..=seq_wires.len() {
                let wire = seq_wires[j - 1].clone();
                if !is_degenerated_wire(brep, &wire) {
                    continue;
                }
                // OCCT L4972-4992.
                for i_e in iter_no_cumori(&wire) {
                    let applied = self.my_re_shape.apply(brep, &i_e, ShapeType::Shape);
                    let edge = applied;
                    let mut face = Shape::null();
                    if let Some(f) = edge_face.get(&bat::shape_key(&edge)) {
                        face = f.clone();
                    }
                    let degedge = degenerated_section(brep, &edge, &face);
                    if degedge.is_null() {
                        continue;
                    }
                    if !degedge.is_same(&edge) {
                        replace_edge(brep, &edge, &degedge, &mut self.my_re_shape);
                    }
                    if brep_tool_degenerated_checked(&degedge) {
                        set_add(&mut self.my_degenerated, &degedge);
                    }
                }
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::EdgeRegularity(theProgress) (cxx
    /// L5003-5028) — updates the Continuity flag on newly created edges.
    pub(crate) fn edge_regularity(&mut self, brep: &mut BRep) {
        // OCCT L5005-5007.
        let mut a_map_ef: super::IdxListMap = super::IdxListMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &self.my_sewed_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_map_ef,
        );

        // OCCT L5009-5021: the myMergedEdges walk.
        let merged: Vec<Shape> = self.my_merged_edges.values().cloned().collect();
        for a_me_it in &merged {
            let an_edge = self.my_re_shape.apply(brep, a_me_it, ShapeType::Shape);
            // encode regularity if and only if edges is shared by two faces
            if let Some((_, a_faces)) = a_map_ef.get(&bat::shape_key(&an_edge)) {
                if a_faces.len() == 2 {
                    brep_lib_encode_regularity_edge(
                        brep,
                        &an_edge,
                        &a_faces[0],
                        &a_faces[1],
                        1.0,
                    );
                }
            }
        }

        // OCCT L5023-5027.
        self.my_merged_edges.clear();
    }
}

/// OCCT `BRep_Tool::Degenerated(TopoDS::Edge(degedge))`.
fn brep_tool_degenerated_checked(e: &Shape) -> bool {
    super::brep_tool_degenerated(e)
}

/// OCCT static GetSeqEdges(edge, seqEdges, VertEdge) (cxx L4604-4652).
fn get_seq_edges(
    edge: &Shape,
    seq_edges: &mut Vec<Shape>,
    vert_edge: &super::DataListMap,
) {
    // OCCT L4606: int numV = 0.
    let mut num_v = 0usize;
    // OCCT L4607: TopoDS_Iterator Iv(edge, false).
    for iv in iter_no_cumori(edge) {
        let v1 = iv;
        num_v += 1;
        let list_edges = match vert_edge.get(&bat::shape_key(&v1)) {
            Some(l) => l.clone(),
            None => continue,
        };
        for l_it in &list_edges {
            let edge1 = l_it.clone();
            if edge1.is_same(edge) {
                continue;
            }
            // OCCT L4616-4626.
            let mut is_contained = false;
            let mut index = 1usize;
            for i in 1..=seq_edges.len() {
                if !is_contained {
                    is_contained = seq_edges[i - 1].is_same(&edge1);
                }
                if !is_contained && seq_edges[i - 1].is_same(edge) {
                    index = i;
                }
            }
            if !is_contained {
                if num_v == 1 {
                    seq_edges.insert(index - 1, edge1.clone());
                } else {
                    if index < seq_edges.len() {
                        seq_edges.insert(index, edge1.clone());
                    } else {
                        seq_edges.push(edge1.clone());
                    }
                }
                get_seq_edges(&edge1, seq_edges, vert_edge);
            }
        }
    }
}

/// OCCT BRepBuilderAPI_Sewing::GetFreeWires(MapFreeEdges, seqWires) (cxx
/// L4654-4707) — get wires from free edges.
pub(crate) fn get_free_wires(
    brep: &mut BRep,
    map_free_edges: &mut super::ShapeSet,
    seq_wires: &mut Vec<Shape>,
) {
    // OCCT L4655-4669.
    let mut vert_edge: super::DataListMap = super::DataListMap::new();
    let mut seq_free_edges: Vec<Shape> = Vec::new();
    for i in 0..map_free_edges.len() {
        let edge = super::set_get(map_free_edges, i);
        seq_free_edges.push(edge.clone());
        for iv in iter_no_cumori(&edge) {
            let v1 = iv;
            match vert_edge.get_mut(&bat::shape_key(&v1)) {
                Some(list) => list.push(edge.clone()),
                None => {
                    vert_edge.insert(bat::shape_key(&v1), vec![edge.clone()]);
                }
            }
        }
    }
    // OCCT L4672-4704.
    let mut b = BRepBuilder::new();
    for i in 1..=seq_free_edges.len() {
        let mut seq_edges: Vec<Shape> = Vec::new();
        let edge = seq_free_edges[i - 1].clone();
        if !set_contains(map_free_edges, &edge) {
            continue;
        }
        seq_edges.push(edge.clone());
        get_seq_edges(&edge, &mut seq_edges, &vert_edge);
        let wire = b.make_wire(brep);
        for j in 1..=seq_edges.len() {
            b.add_to_wire(brep, wire.clone(), seq_edges[j - 1].clone());
            map_free_edges.shift_remove(&bat::shape_key(&seq_edges[j - 1]));
        }
        seq_wires.push(wire.clone());
        if map_free_edges.is_empty() {
            break;
        }
    }
}

/// OCCT static IsDegeneratedWire(wire) (cxx L4709-4798).
fn is_degenerated_wire(brep: &BRep, wire: &Shape) -> bool {
    // OCCT L4710-4713.
    if wire.shape_type() != ShapeType::Wire {
        return false;
    }
    // Get maximal vertices tolerance
    // OCCT L4716-4722.
    let mut v1 = Shape::null();
    let mut v2 = Shape::null();
    let mut wire_length = 0.0f64;
    let mut nume = 0usize;
    let mut is_small = 0usize;
    // OCCT L4724: TopoDS_Iterator aIt(wire, false).
    for a_it in iter_no_cumori(wire) {
        nume += 1;
        let edge = a_it.clone();
        let (ve1, ve2) = bat::top_exp_vertices_raw(&edge);
        let ve1 = ve1.unwrap_or_else(Shape::null);
        let ve2 = ve2.unwrap_or_else(Shape::null);
        // OCCT L4729-4743: the endpoint chain walk.
        if nume == 1 {
            v1 = ve1.clone();
            v2 = ve2.clone();
        } else {
            if ve1.is_same(&v1) {
                v1 = ve2.clone();
            } else if ve1.is_same(&v2) {
                v2 = ve2.clone();
            }
            if ve2.is_same(&v1) {
                v1 = ve1.clone();
            } else if ve2.is_same(&v2) {
                v2 = ve1.clone();
            }
        }
        // OCCT L4745-4766.
        let c3d = super::brep_tool_curve_world(brep, &edge);
        if let Some((c3d, first, last)) = c3d {
            let pfirst = c3d.point_at(first);
            let plast = c3d.point_at(last);
            let pmid = c3d.point_at((first + last) * 0.5);
            let length = if pfirst.distance(plast) > pfirst.distance(pmid) {
                pfirst.distance(plast)
            } else {
                super::gcpnts_abscissa_length(&c3d, first, last)
            };
            let tole = bat::brep_tool_tolerance(&ve1) + bat::brep_tool_tolerance(&ve2);
            if length <= tole {
                is_small += 1;
            }
            wire_length += length;
        }
    }
    // OCCT L4779-4782.
    if is_small == nume {
        return true;
    }
    // OCCT L4784-4797.
    let tol = bat::brep_tool_tolerance(&v1) + bat::brep_tool_tolerance(&v2);
    wire_length <= tol
}

/// OCCT static DegeneratedSection(section, face) (cxx L4800-4910) — creates
/// a new degenerated edge if the section is degenerated.
fn degenerated_section(brep: &mut BRep, section: &Shape, face: &Shape) -> Shape {
    // Return if section is already degenerated
    // OCCT L4803-4806.
    if super::brep_tool_degenerated(section) {
        return section.clone();
    }

    // Retrieve edge curve
    // OCCT L4809-4822.
    let (c3d, first, last) = match super::brep_tool_curve_world(brep, section) {
        Some(v) => v,
        None => {
            // gka
            // OCCT L4812-4818.
            let edge1 = section.clone();
            let mut a_b = BRepBuilder::new();
            a_b.set_edge_degenerated(brep, edge1.clone(), true);
            return edge1;
        }
    };

    // Test if the new edge is degenerated
    // OCCT L4825-4830.
    let (v1, v2) = bat::top_exp_vertices_raw(section);
    let v1 = v1.unwrap_or_else(Shape::null);
    let v2 = v2.unwrap_or_else(Shape::null);

    // OCCT L4838-4840.
    let p1 = bat::brep_tool_pnt(&v1).unwrap_or(glam::DVec3::ZERO);
    let p3 = bat::brep_tool_pnt(&v2).unwrap_or(glam::DVec3::ZERO);

    // processing
    // OCCT L4853-4877.
    let mut a_builder = BRepBuilder::new();
    let edge = brep.empty_copied(section);
    if v1.is_same(&v2) {
        // OCCT L4857-4862.
        let an_edge = bat::oriented(&edge, Orientation::Forward);
        a_builder.add_to_edge(brep, an_edge.clone(), bat::oriented(&v1, Orientation::Forward));
        a_builder.add_to_edge(brep, an_edge, bat::oriented(&v2, Orientation::Reversed));
    } else {
        // OCCT L4863-4875.
        let p2 = c3d.point_at(0.5 * (first + last));
        let new_vertex = if p1.distance(p3) < bat::brep_tool_tolerance(&v1) {
            v1.clone()
        } else if p1.distance(p3) < bat::brep_tool_tolerance(&v2) {
            v2.clone()
        } else {
            let d1 = bat::brep_tool_tolerance(&v1) + p2.distance(p1);
            let d2 = bat::brep_tool_tolerance(&v2) + p2.distance(p3);
            let new_tolerance = d1.max(d2);
            super::make_vertex_at(brep, p2, new_tolerance)
        };
        let an_edge = bat::oriented(&edge, Orientation::Forward);
        a_builder.add_to_edge(brep, an_edge.clone(), bat::oriented(&new_vertex, Orientation::Forward));
        a_builder.add_to_edge(brep, an_edge, bat::oriented(&new_vertex, Orientation::Reversed));
    }

    // OCCT L4880-4882.
    let an_edge = bat::oriented(&edge, Orientation::Forward);
    a_builder.set_edge_range(brep, an_edge.clone(), first, last);
    a_builder.set_edge_degenerated(brep, edge.clone(), true);
    // OCCT L4883-4897.
    if !face.is_null() {
        let c2dt = bat::brep_tool_curve_on_surface(section, face);
        // OCCT L4886: aBuilder.UpdateEdge(edge, aC3dNew, 0) — the null-curve
        // update; the rcad carrier stores no null curve (no-op).
        let c2dn = bat::brep_tool_curve_on_surface(&edge, face);
        if c2dn.is_none() {
            if let Some((c2dt, _af, _al)) = c2dt {
                a_builder.update_edge_pcurve(brep, edge.clone(), c2dt, face.clone(), 0.0);
            }
        }
    }

    // OCCT L4909.
    edge
}

// The pcurve type anchor (Curve2d is consumed by the CreateSections group).
#[allow(unused_imports)]
use Curve2d as _Curve2dAnchor;
