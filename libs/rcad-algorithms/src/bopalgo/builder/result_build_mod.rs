use super::ResultBuilder;
use super::{BooleanBuilder, SourceSide};
use crate::bopalgo::builder::types::*;
use crate::bopalgo::builder::wire_splitter::{EdgeInfo, build_closed_wires};
use crate::bopalgo::{Alert, GlueEnum, Report};
use crate::bopds::ds::*;
use crate::boptools::bvh::{Aabb, BoxTree};
use crate::classify::{Classification, classify_point};
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use crate::inttools::context::Context;
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use glam::{DVec2, DVec3};
use indexmap::IndexMap;
use rcad_kernel::PCurve;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval, *};
use rcad_kernel::topods;
use rcad_kernel::topology::*;
use std::collections::{HashMap, HashSet, VecDeque};

impl<'a> BooleanBuilder<'a> {
    /// ✅ OCCT-aligned: BuildRC (BOPAlgo_BOP.cxx L597-881).
    pub(super) fn build_rc(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        let solids = std::mem::take(&mut result.tmp_solids);
        let sides: Vec<usize> = result.solid_side_origin.clone();

        // A. FUSE -- keep all split solids (fence-deduped)
        if self.my_operation == BooleanOpType::Union {
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut kept: Vec<Vec<usize>> = Vec::new();
            for s in &solids {
                if a_m_fence.insert(s.clone()) {
                    kept.push(s.clone());
                }
            }
            result.tmp_solids = kept;
            return;
        }

        // B. COMMON/CUT/CUT21: prepare building elements of arguments/tools
        let e_base = self.ds.vertex_count();
        let f_base = e_base + self.ds.edge_count() + self.ds.shape_info.iter().filter(|si| si.shape_type == rcad_kernel::topods::ShapeType::Wire && !si.is_new).count();

        // Map source shapes (V/E/F) per side
        let mut a_m_args: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_maps = [&mut a_m_args, &mut a_m_tools];

        for (side_idx, a_ms) in a_maps.iter_mut().enumerate() {
            let v_range = if side_idx == 0 {
                (0usize, self.ds.a_vertex_count)
            } else {
                (self.ds.a_vertex_count, self.ds.vertex_count())
            };
            let e_range = if side_idx == 0 {
                (0usize, self.ds.a_edge_count)
            } else {
                (self.ds.a_edge_count, self.ds.edge_count())
            };
            let f_range = if side_idx == 0 {
                (0usize, self.ds.a_face_count)
            } else {
                (self.ds.a_face_count, self.ds.face_count())
            };
            for vi in v_range.0..v_range.1 {
                a_ms.insert(vi);
            }
            for ei in e_range.0..e_range.1 {
                a_ms.insert(e_base + ei);
            }
            for fi in f_range.0..f_range.1 {
                a_ms.insert(f_base + fi);
            }
        }

        // Get splits of building elements (check myImages for split images)
        let mut a_m_args_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mset_args: Vec<std::collections::BTreeSet<usize>> = Vec::new();
        let mut a_mset_tools: Vec<std::collections::BTreeSet<usize>> = Vec::new();

        let mut im_maps = [&mut a_m_args_im, &mut a_m_tools_im];
        let mut set_maps = [&mut a_mset_args, &mut a_mset_tools];

        for (side_idx, (a_ms_im, a_mset)) in im_maps.iter_mut().zip(set_maps.iter_mut()).enumerate()
        {
            let a_ms = &a_maps[side_idx];
            let side_is_args = side_idx == 0;

            let mut sorted_elements: Vec<&usize> = a_ms.iter().collect();
            sorted_elements.sort();

            for &&flat_idx in &sorted_elements {
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                let local_idx = if is_edge {
                    flat_idx - e_base
                } else if is_face {
                    flat_idx - f_base
                } else {
                    flat_idx
                };

                if is_edge {
                    if self.ds.is_edge_degenerated(local_idx) {
                        continue;
                    }
                }

                let has_images = if is_edge {
                    self.my_images
                        .borrow()
                        .contains_key(&self.brep_sr(flat_idx))
                } else if is_face {
                    // Use my_images (OCCT myImages) instead of result.face_origins
                    // which may be consumed by fill_same_domain_faces.
                    let im = self
                        .my_images
                        .borrow()
                        .contains_key(&self.brep_sr(flat_idx));
                    im
                } else {
                    false
                };

                if has_images {
                    let (o_exp, sfi) = if side_is_args {
                        (ShapeOrigin::ShapeA, local_idx)
                    } else {
                        (ShapeOrigin::ShapeB, local_idx)
                    };

                    if is_face {
                        // Use my_images (OCCT myImages) instead of result.face_origins
                        // which may be consumed by fill_same_domain_faces.
                        if let Some(imgs) = self.my_images.borrow().get(&self.brep_sr(flat_idx)) {
                            for &sr in imgs {
                                a_ms_im.insert(f_base + sr.index);
                            }
                        }
                    } else if is_edge {
                        if let Some(imgs) = self.my_images.borrow().get(&self.brep_sr(flat_idx)) {
                            for &sr in imgs {
                                a_ms_im.insert(e_base + sr.index);
                            }
                        }
                    }
                } else {
                    a_ms_im.insert(flat_idx);

                    if is_face {
                        let mut a_st: std::collections::BTreeSet<usize> =
                            std::collections::BTreeSet::new();
                        let (o_exp2, sfi2) = if side_is_args {
                            (ShapeOrigin::ShapeA, local_idx)
                        } else {
                            (ShapeOrigin::ShapeB, local_idx)
                        };
                        a_st.insert(local_idx);
                        // Add sibling faces from the same shell
                        for (dfi2, df2) in self.ds.faces.iter().enumerate() {
                            if dfi2 != local_idx
                                && df2.origin == o_exp2
                                && df2.source_shell_idx == self.ds.source_shell_idx(local_idx)
                            {
                                a_st.insert(dfi2);
                            }
                        }
                        if !a_mset.contains(&a_st) {
                            a_mset.push(a_st);
                        }
                    }
                }
            }
        }

        let b_common = self.my_operation == BooleanOpType::Intersection;
        let b_cut21 = false; // --rcad: CUT21 not supported

        let a_m_it: &std::collections::HashSet<usize> =
            if b_cut21 { &a_m_tools_im } else { &a_m_args_im };
        let a_m_check: &std::collections::HashSet<usize> =
            if b_cut21 { &a_m_args_im } else { &a_m_tools_im };
        let a_mset_check: &Vec<std::collections::BTreeSet<usize>> =
            if b_cut21 { &a_mset_args } else { &a_mset_tools };

        // Expand sub-shapes for COMMON
        let a_m_it_exp: std::collections::HashSet<usize> = if b_common {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_it.iter().collect::<Vec<_>>() {
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.face_count() {
                        for &ei in self.ds.face_boundary_edges(local_fi) {
                            exp.insert(e_base + ei);
                        }
                        for &vi in self.ds.face_boundary_verts(local_fi) {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edge_count() {
                        exp.insert(self.ds.edge_start_vertex_ds(local_ei));
                        exp.insert(self.ds.edge_end_vertex_ds(local_ei));
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        } else {
            a_m_it.clone()
        };

        // Expand check side too
        let a_m_check_exp: std::collections::HashSet<usize> = {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_check.iter().collect::<Vec<_>>() {
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.face_count() {
                        for &ei in self.ds.face_boundary_edges(local_fi) {
                            exp.insert(e_base + ei);
                        }
                        for &vi in self.ds.face_boundary_verts(local_fi) {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edge_count() {
                        exp.insert(self.ds.edge_start_vertex_ds(local_ei));
                        exp.insert(self.ds.edge_end_vertex_ds(local_ei));
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        };

        // Compare building-element images and build keep set.
        let mut keep_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &&flat_idx in &a_m_it_exp.iter().collect::<Vec<_>>() {
            let mut b_contains = a_m_check_exp.contains(&flat_idx);
            let is_face = flat_idx >= f_base;
            if !b_contains && is_face {
                let local_fi = flat_idx - f_base;
                let mut a_st: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                if local_fi < self.ds.face_count() {
                    a_st.insert(local_fi);
                    for &vi in self.ds.face_boundary_verts(local_fi) {
                        a_st.insert(vi);
                    }
                    for &ei in self.ds.face_boundary_edges(local_fi) {
                        a_st.insert(e_base + ei);
                    }
                }
                b_contains = a_mset_check.iter().any(|s| s == &a_st);
            }
            let keep = if b_common { b_contains } else { !b_contains };
            if keep {
                keep_set.insert(flat_idx);
            }
        }

        // Filter result.tmp_solids: keep solids whose building elements pass.
        let mut kept_solids: Vec<Vec<usize>> = Vec::new();
        for (i, solid_shells) in solids.iter().enumerate() {
            let side = sides.get(i).copied().unwrap_or(0);
            // Check each solid: if ANY face's building element is in keep_set --keep.
            // A result solid is kept iff the source face(s) it was split from pass.
            let mut solid_keep = false;
            for &si in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(si) {
                    for &rfi in shell_faces {
                        // Check both index spaces that keep_set may use:
                        //   - f_base + rfi for split faces (images, line 174)
                        //   - f_base + dfi for unsplit faces (original, line 188)
                        let rfi_flat = f_base + rfi;
                        if keep_set.contains(&rfi_flat) {
                            solid_keep = true;
                            break;
                        }
                        let dfi_opt = match result.face_origins.get(rfi) {
                            Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f| {
                                f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi
                            }),
                            Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f| {
                                f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi
                            }),
                            _ => None,
                        };
                        if let Some(dfi) = dfi_opt {
                            let flat_fi = f_base + dfi;
                            if keep_set.contains(&flat_fi) {
                                solid_keep = true;
                                break;
                            }
                        }
                    }
                    if solid_keep {
                        break;
                    }
                }
            }
            if solid_keep {
                kept_solids.push(solid_shells.clone());
            }
        }

        //   rcad: OCCT re-iterates the compound by dimension (SOLID→SHELL→FACE) with
        //   a fence.  rcad solids are already at SOLID granularity; shell-count fence
        //   prevents duplicates (OCCT L799-804 fence at FACE+WIRE level).
        if b_common {
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut reordered: Vec<Vec<usize>> = Vec::new();
            for s in &kept_solids {
                if a_m_fence.insert(s.clone()) {
                    reordered.push(s.clone());
                }
            }
            kept_solids = reordered;
        }

        // Result edges are embedded in pre-assembled solids — no standalone DEs.

        result.tmp_solids = kept_solids;

        // Rebuild t_brep from DS faces in keep_set (OCCT L800-823: myRC).
        if !keep_set.is_empty() {
            let mut kept_refs: Vec<topods::ShapeRef> = Vec::new();
            let fr = self.my_face_refs.borrow();
            for (rfi, &sr) in fr.iter().enumerate() {
                if sr.is_null() {
                    continue;
                }
                let dsfi = if sr.index >= f_base && sr.index < f_base + self.ds.face_count() {
                    Some(sr.index - f_base)
                } else {
                    None
                };
                if let Some(dsfi) = dsfi {
                    if keep_set.contains(&(f_base + dsfi)) {
                        kept_refs.push(sr);
                    }
                }
            }
            if !kept_refs.is_empty() && kept_refs.len() < fr.len() {
                let mut nb = topods::BRep::new();
                nb.add_tshell(kept_refs);
                *t_brep = nb;
            }
        }
    }

    /// FillInternalShapes (Builder_3.cxx L622-887).
    /// Settles internal sub-shapes (vertices, edges) into result solids.
    ///
    /// rcad: internal V/E are marked via is_internal flag in DS.
    pub(super) fn detect_internal_voids(
        &self,
        result: &mut ResultBuilder,
        assignments: &[(usize, usize, &'static str)],
    ) {
        // Precompute DS face set, centroid, AABB per solid.
        //   OCCT operates on raw shells (myLoops); rcad operates on result.tmp_solids.
        let n_solids = result.tmp_solids.len();
        let mut ds_faces_of: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        let mut centroids: Vec<DVec3> = Vec::with_capacity(n_solids);
        let mut aabbs: Vec<Aabb> = Vec::with_capacity(n_solids);
        for si in 0..n_solids {
            let mut faces = Vec::new();
            let mut aabb = Aabb::empty();
            let mut centroid = DVec3::ZERO;
            for &sh in &result.tmp_solids[si] {
                if let Some(shell) = result.tmp_shells.get(sh) {
                    for &fi in shell {
                        if let Some(origin) = result.face_origins.get(fi) {
                            let ds_fi = match origin {
                                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f| {
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi
                                }),
                                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f| {
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi
                                }),
                                _ => None,
                            };
                            if let Some(dfi) = ds_fi {
                                faces.push(dfi);
                                for &vi in self.ds.face_boundary_verts(dfi) {
                                    if vi < self.ds.vertex_count() {
                                        aabb.expand_point(self.ds.vertex_point(vi));
                                    }
                                }
                                if fi < result.faces.len() {
                                    centroid = result.faces[fi].6;
                                }
                            }
                        }
                    }
                }
            }
            faces.sort_unstable();
            faces.dedup();
            ds_faces_of.push(faces);
            centroids.push(centroid);
            aabbs.push(aabb);
        }

        // Classify each shell as Growth or Hole
        let mut a_mhf: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut is_hole = vec![false; n_solids];

        for si in 0..n_solids {
            let is_growth = if !a_mhf.is_empty() {
                ds_faces_of[si].iter().any(|dfi| a_mhf.contains(dfi))
            } else {
                false
            };

            if !is_growth {
                let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
                let dead_faces: Vec<usize> = self
                    .ds
                    .faces
                    .iter()
                    .enumerate()
                    .filter(|(_, f)| match side {
                        0 => f.origin == ShapeOrigin::ShapeA,
                        _ => f.origin == ShapeOrigin::ShapeB,
                    })
                    .map(|(fi, _)| fi)
                    .collect();
                let class = classify_point(centroids[si], &dead_faces, self.ds);
                is_hole[si] = class == Classification::In;
            }

            if is_hole[si] {
                for &dfi in &ds_faces_of[si] {
                    a_mhf.insert(dfi);
                }
            }
        }

        let in_si: Vec<usize> = (0..n_solids).filter(|&i| is_hole[i]).collect();
        let out_si: Vec<usize> = (0..n_solids).filter(|&i| !is_hole[i]).collect();

        if in_si.is_empty() || out_si.is_empty() {
            return;
        }

        // Build BVH of hole shells
        let hole_key: Vec<usize> = in_si.clone();
        let hole_aabbs: Vec<Aabb> = in_si.iter().map(|&i| aabbs[i]).collect();
        let hole_bvh = crate::boptools::bvh::BoxTree::build(hole_key, hole_aabbs);

        // Classify holes against growth solids
        let mut a_hole_solid_map: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for &os in &out_si {
            let candidates = hole_bvh.query_aabb(&aabbs[os]);

            for &hole_idx in &candidates {
                let class = classify_point(centroids[hole_idx], &ds_faces_of[os], self.ds);
                if class != Classification::In && class != Classification::On {
                    continue;
                }

                // Select outermost containing solid
                use std::collections::hash_map::Entry;
                match a_hole_solid_map.entry(hole_idx) {
                    Entry::Occupied(mut e) => {
                        let prev_os = *e.get();
                        let prev_faces = &ds_faces_of[prev_os];
                        let os_inside = classify_point(centroids[os], prev_faces, self.ds);
                        if os_inside == Classification::In || os_inside == Classification::On {
                            e.insert(os);
                        }
                    }
                    Entry::Vacant(e) => {
                        e.insert(os);
                    }
                }
            }
        }

        // Build reverse map: solid → list of holes
        let mut solid_holes_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (&hole_idx, &os) in &a_hole_solid_map {
            solid_holes_map.entry(os).or_default().push(hole_idx);
        }

        // Add holes to solids
        let mut removed = vec![false; n_solids];
        for (&os, holes) in &solid_holes_map {
            for &hole_idx in holes {
                let void_shells = result.tmp_solids[hole_idx].clone();
                result.tmp_solids[os].extend(void_shells);
                removed[hole_idx] = true;
            }
        }
        let mut new_solids: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        for (si, solid) in result.tmp_solids.drain(..).enumerate() {
            if !removed[si] {
                new_solids.push(solid);
            }
        }
        result.tmp_solids = new_solids;
    }

    /// ✅ OCCT-aligned: FillInternalShapes (BOPAlgo_Builder_3.cxx L622-887).
    ///   Phase 1 (L648-709): Collect internal V/E/WIRE from arguments.
    ///   Phase 2 (L717-788): Internal V/E from source solids + build aMSx ancestry.
    ///   Phase 3 (L790-809): Filter shapes already attached via aMSx.
    ///   Phase 4 (L811-816): Early return if none.
    ///   Phase 5 (L820-877): Classify each internal shape against each split solid;
    ///     if IN --add as INTERNAL sub-shape (clone original if needed).
    pub(super) fn fill_internal_shapes(&self, result: &mut ResultBuilder) {
        //   rcad: adapted to Vec/HashSet equivalents.

        // === Phase 1: Shapes to process --collect from arguments (OCCT L648-709) ===
        //   OCCT L653-658: TreatCompound on each argument --flatten into aLSC.
        //   OCCT L660-681: filter VERTEX/EDGE/WIRE from aLSC --aLArgs.
        //   OCCT L684-709: for each aLArgs, check myImages.IsBound --aMSI (images or originals).
        //   rcad: DS vertices/edges with is_internal flag = sources.
        //   --TreatCompound: rcad treats DS V/E as already-flattened source shapes.
        //   --aMSI: maps shape-ref --true if it's an internal shape to process.
        let mut a_msi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // Collect internal vertices
        for (vi, v) in self.ds.vertices.iter().enumerate() {
            if v.is_internal {
                let v_ref = self.brep_sr(vi);
                if self.my_images.borrow().contains_key(&v_ref) {
                    for img in &self.my_images.borrow()[&v_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(vi);
                }
            }
        }
        // Collect internal edges
        for (ei, e) in self.ds.edges.iter().enumerate() {
            if e.is_internal {
                let e_ref = self.brep_sr(self.ds.vertex_count() + ei);
                if self.my_images.borrow().contains_key(&e_ref) {
                    for img in &self.my_images.borrow()[&e_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(ei);
                }
            }
        }

        // Build aMSx ancestry: internal shapes already on split-solid faces
        let mut a_msx: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        let mut a_lsd: Vec<usize> = Vec::new();

        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
            let mut edge_to_faces: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for &shi in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(shi) {
                    for &rfi in shell_faces {
                        if let Some(fe) = result.faces.get(rfi) {
                            for &(ei, _) in &fe.0 {
                                edge_to_faces.entry(ei).or_default().push(rfi);
                            }
                        }
                    }
                }
            }
            for (&ei, face_list) in &edge_to_faces {
                a_msx.entry(ei).or_default().extend(face_list);
                if ei < self.ds.edge_count() {
                    a_msx
                        .entry(self.ds.edge_start_vertex_ds(ei))
                        .or_default()
                        .push(ei);
                    a_msx
                        .entry(self.ds.edge_end_vertex_ds(ei))
                        .or_default()
                        .push(ei);
                }
            }
            a_lsd.push(si);
        }

        // Filter shapes already attached to split-solid faces
        let mut a_lsi: Vec<usize> = Vec::new();
        for &si in &a_msi {
            let is_attached = a_msx.get(&si).map_or(false, |anc| !anc.is_empty());
            if !is_attached {
                a_lsi.push(si);
            }
        }

        if a_lsi.is_empty() {
            return;
        }

        // Settle internal V/E into solids
        for &si in &a_lsd {
            let mut solid_ds_faces: Vec<usize> = Vec::new();
            if let Some(solid_shells) = result.tmp_solids.get(si) {
                for &shi in solid_shells {
                    if let Some(shell_faces) = result.tmp_shells.get(shi) {
                        for &rfi in shell_faces {
                            let dfi_opt = match result.face_origins.get(rfi) {
                                Some(FaceOrigin::FromA(sfi)) => {
                                    self.ds.faces.iter().position(|f| {
                                        f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi
                                    })
                                }
                                Some(FaceOrigin::FromB(sfi)) => {
                                    self.ds.faces.iter().position(|f| {
                                        f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi
                                    })
                                }
                                _ => None,
                            };
                            if let Some(dfi) = dfi_opt {
                                solid_ds_faces.push(dfi);
                            }
                        }
                    }
                }
            }
            solid_ds_faces.sort_unstable();
            solid_ds_faces.dedup();
            if solid_ds_faces.is_empty() {
                continue;
            }

            let mut i = 0usize;
            while i < a_lsi.len() {
                let si_idx = a_lsi[i];
                let pt = if si_idx < self.ds.vertex_count() {
                    self.ds.vertex_point(si_idx)
                } else {
                    let ei = si_idx;
                    if ei < self.ds.edge_count() {
                        (self.ds.vertex_point(self.ds.edge_start_vertex_ds(ei))
                            + self.ds.vertex_point(self.ds.edge_end_vertex_ds(ei)))
                            * 0.5
                    } else {
                        i += 1;
                        continue;
                    }
                };
                let a_state = classify_point(pt, &solid_ds_faces, self.ds);

                if a_state != Classification::In {
                    i += 1;
                    continue;
                }

                // Shape is IN — add as INTERNAL to the solid
                if let Some(&first_shi) = result.tmp_solids.get(si).and_then(|s| s.first()) {
                    if let Some(shell_faces) = result.tmp_shells.get(first_shi) {
                        if let Some(&first_rfi) = shell_faces.first() {
                            if first_rfi < result.face_internal_vtx.len() {
                                if si_idx < self.ds.vertex_count() {
                                    if !result.face_internal_vtx[first_rfi].contains(&si_idx) {
                                        result.face_internal_vtx[first_rfi].push(si_idx);
                                    }
                                }
                            }
                        }
                    }
                }

                a_lsi.swap_remove(i);
            }
        }
    }

    /// --FillImagesCompounds (Builder_1.cxx L197-342).
    ///   Phase 7: group result solids into COMPSOLID/COMPOUND hierarchy.
    ///
    /// OCCT flow:
    ///   L200-217 (FillImagesCompounds): Iterate source shapes for TopAbs_COMPOUND.
    ///     For each compound, call FillImagesCompound recursively.
    ///   L280-342 (FillImagesCompound): Recursively check each child for images.
    ///     If any child has images, build a new compound with image replacements.
    ///     Result stored in myImages[original_compound] = new_compound.
    ///
    /// rcad: records compound intent in ResultBuilder.  Actual compound
    ///   reconstruction happens after result.build() in build_with_history
    ///   (see the rebuild_compound_for_step post-step) because the result
    ///   BRep solids don't exist until build() is called.
    /// --FillImagesCompounds (Builder_1.cxx L197-217).
    ///
    /// OCCT FillImagesCompounds L197-217:
    ///   L200: aMFP fence map
    ///   L202-216: iterate source shapes, filter TopAbs_COMPOUND,
    ///             call FillImagesCompound(aC, aMFP)
    /// OCCT FillImagesCompound L280-342:
    ///   L290-293: fence --skip if processed
    ///   L296-308: check if any sub-shape has images
    ///   L309-312: if none modified --return
    ///   L314-341: build new compound from sub-shape (solid) images
    ///
    /// rcad: source compound solids are tracked in DS solid_images.
    ///   Compound reconstruction from result solids is deferred to
    ///   build_with_history's post-step (L6834-6840) because the
    ///   result BRep solids don't exist until ResultBuilder::build().
    /// --FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: aMFP fence map; NbSourceShapes --filter TopAbs_COMPOUND.
    ///   L280-293: FillImagesCompound --fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification --return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    /// FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: dispatcher with fence map; iterate source COMPOUND shapes.
    ///   L280-293: FillImagesCompound --fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification --return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    ///   --rcad: no compound nesting in DS.  Flat per-face source_compsolid_idx.
    ///     The recursive FillImagesCompound is collapsed to a single level.
    /// FillImagesCompounds (BOPAlgo_Builder_1.cxx L199-341).
    /// ✅ OCCT-aligned: FillImagesCompounds (BOPAlgo_Builder_1.cxx L197-216 + L280-342).
    /// Wraps sub-shape images into compound containers for non-destructive history.
    /// OCCT iterates source shapes, filters COMPOUND,
    ///   recursively processes each compound.  rcad: iterate compsolid groups from DS
    ///   faces (compounds not in shape_info for standard BRep inputs).
    ///
    ///   OCCT FillImagesCompound(theS, theMFP):
    ///     1. L290-293: Fence — theMFP.Add(theS) → skip if processed.
    ///     2. L295-312: Interference check — sub-shapes with myImages → bInterferred.
    ///        Sub-compounds → recurse.
    ///     3. L314-315: MakeContainer(TopAbs_COMPOUND, aCIm) for new compound.
    ///     4. L317-337: Add sub-shapes with orientation propagation (aSXIm.Orientation(aOrX)).
    ///     5. L339-341: myImages.Bind(theS, aLSIm).
    /// ✅ OCCT-aligned: FillImagesCompounds (BOPAlgo_Builder_1.cxx L197-216) + FillImagesCompound (L280-342).
    pub(super) fn fill_images_compounds(&self, result: &mut ResultBuilder) {
        let mut t = self.my_shape.borrow_mut();
        // OCCT: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFP — fence
        let mut a_mfp: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // Collect unique source_compsolid_idx from DS faces (equivalent to COMPOUND source shapes)
        let mut compound_groups: Vec<usize> = Vec::new();
        for df in &self.ds.faces {
            if let Some(csi) = df.source_compsolid_idx {
                // OCCT L290-293: theMFP.Add(theS) — fence
                if a_mfp.insert(csi) {
                    compound_groups.push(csi);
                }
            }
        }
        if compound_groups.is_empty() {
            return;
        }

        for &csi in &compound_groups {
            // OCCT L295-308: interference check — myImages.IsBound(aSx) for any sub-shape
            // rcad: check if any source solid under this compsolid has result solids
            let mut b_interferred = false;
            let sub_solid_indices: Vec<usize> = {
                let mut seen = std::collections::HashSet::new();
                self.ds
                    .faces
                    .iter()
                    .filter_map(|f| {
                        if f.source_compsolid_idx == Some(csi) {
                            f.source_solid_idx
                        } else {
                            None
                        }
                    })
                    .filter(|&ssi| seen.insert(ssi))
                    .collect()
            };
            for &ssi in &sub_solid_indices {
                // OCCT: myImages.IsBound(aSx) — check if sub-shape has split images
                let side = self
                    .ds
                    .faces
                    .iter()
                    .find(|f| f.source_solid_idx == Some(ssi))
                    .map(|f| match f.origin {
                        crate::bopds::ds::ShapeOrigin::ShapeA => 0,
                        crate::bopds::ds::ShapeOrigin::ShapeB => 1,
                    })
                    .unwrap_or(0);
                let has_images = result
                    .solid_side_origin
                    .iter()
                    .filter(|&&s| s == side)
                    .count()
                    > 0;
                if has_images {
                    b_interferred = true;
                    break;
                }
            }
            // OCCT L309-312: no interference → skip
            if !b_interferred {
                continue;
            }

            // OCCT L314-315: MakeContainer(TopAbs_COMPOUND, aCIm)
            // rcad: collect result solid ShapeRefs
            let mut a_c_im: Vec<usize> = Vec::new();
            for &ssi in &sub_solid_indices {
                let side = self
                    .ds
                    .faces
                    .iter()
                    .find(|f| f.source_solid_idx == Some(ssi))
                    .map(|f| match f.origin {
                        crate::bopds::ds::ShapeOrigin::ShapeA => 0,
                        crate::bopds::ds::ShapeOrigin::ShapeB => 1,
                    })
                    .unwrap_or(0);
                // OCCT L321: aOrX = aSX.Orientation() — get source orientation
                // rcad: orientation is per-face via FaceOrigin (no per-solid orientation)
                for (si, &s) in result.solid_side_origin.iter().enumerate() {
                    if s == side && !a_c_im.contains(&si) {
                        a_c_im.push(si);
                    }
                }
            }

            // OCCT L330: aBB.Add(aCIm, aSXIm) — add shapes to compound
            // OCCT L339-341: myImages.Bind(theS, aLSIm) — bind to myImages
            if !a_c_im.is_empty() {
                let solid_refs: Vec<topods::ShapeRef> = a_c_im
                    .iter()
                    .filter_map(|&si| self.my_solids.borrow().get(si).copied())
                    .collect();
                if !solid_refs.is_empty() {
                    let cmp_ref = t.add_tcompound(solid_refs);
                    // rcad: synthetic key (OCCT: original source shape via theS)
                    let ckey = topods::ShapeRef::synthetic(usize::MAX - a_c_im.len());
                    self.my_images
                        .borrow_mut()
                        .entry(ckey)
                        .or_default()
                        .push(cmp_ref);
                    self.my_compsolid_groups.borrow_mut().push(cmp_ref);
                }
            }
        }
    }

    /// Retrieve the EdgeInfo.is_inside status for the incoming edge at the given vertex.
    pub(super) fn incoming_edge_is_inside(
        &self,
        smart_map: &IndexMap<usize, Vec<EdgeInfo>>,
        vertex: usize,
        seg_idx: usize,
    ) -> bool {
        smart_map
            .get(&vertex)
            .and_then(|infos| infos.iter().find(|ei| ei.seg_idx == seg_idx && ei.in_flag))
            .map_or(false, |ei| ei.is_inside)
    }

    /// --face keep/discard policy (ComputeState --FillIn3DParts equivalent).
    ///   OCCT does NOT have a surface-type special case --ComputeState propagates
    ///   ON→IN/OUT based on face orientation + solid side, not surface type.
    /// --BOPAlgo_Builder::FillImagesFaces --face keep policy.
    ///   OCCT: after ComputeState returns IN/OUT/ON for a face against the other solid:
    ///     FUSE: keep OUT + ON
    ///     COMMON: keep IN + ON
    ///     CUT A-B:
    ///       face from A --keep if OUT or ON (A outside B)
    ///       face from B --keep if IN or ON (B inside A, the cut surface)
    pub(super) fn classification_keep_policy(
        &self,
        source: SourceSide,
        class: Classification,
        _fi: usize,
    ) -> bool {
        match self.my_operation {
            BooleanOpType::Intersection => {
                class == Classification::In || class == Classification::On
            }
            BooleanOpType::Difference => match source {
                SourceSide::A => class != Classification::In,
                SourceSide::B => class == Classification::In || class == Classification::On,
            },
            BooleanOpType::Union => class != Classification::In,
        }
    }

    /// --BuildResult --add split images to result (Builder_1.cxx L130-168).
    ///   OCCT: for each source shape of theType, if myImages bound --add images;
    ///   else add the original shape.  rcad: for Edge, creates topods edges in t_brep
    ///   (equivalent to OCCT's myShape) AND flat edge refs in result for face construction.
    ///   For Vertex/Wire/Shell/Solid, rcad handles these in other pipeline steps.
    /// BuildResult (Builder_1.cxx L130-168).
    ///   Add split images (or originals) of source shapes into the result.
    ///   OCCT L133: aMFence fence map.
    ///   L136-167: for each source argument of matching type --if myImages bound
    ///     --add all image shapes; else --add the original shape.
    /// --BuildResult (Builder_1.cxx L130-168).
    ///   Generic loop over myArguments matching OCCT form for ALL types.
    /// ✅ OCCT-aligned: BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    ///   For each argument (myArguments):
    ///     Path A — argument's ShapeType == theType: process directly.
    ///     Path B — TopExp_Explorer traverses sub-shapes of theType:
    ///       if myImages.IsBound(sub) → fence-add image shapes
    ///       else → fence-add the original sub-shape.
    ///   rcad: for Face with solid/compound arguments (Path B), DS face array
    ///   provides the equivalent of TopExp_Explorer, grouped by argument side.
    pub(super) fn build_result(
        &self,
        the_type: rcad_kernel::topods::ShapeType,
        result: &mut ResultBuilder,
    ) {
        let mut t = self.my_shape.borrow_mut();
        let mut a_m_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        if the_type == rcad_kernel::topods::ShapeType::Edge && self.my_edge_map.borrow().is_empty()
        {
            *self.my_edge_map.borrow_mut() =
                vec![rcad_kernel::topods::ShapeRef::NULL; self.ds.edge_count()];
        }

        let args = self.my_arguments.borrow();
        for (arg_idx, a_s) in args.iter().enumerate() {
            if a_s.shape_type(&*t) == the_type {
                // Path A — argument directly matches target type.
                let has_images = self.my_images.borrow().contains_key(a_s);
                if !has_images {
                    if a_m_fence.insert(a_s.index) {
                        self.add_to_result(*a_s, the_type, result, &mut *t);
                    }
                } else if let Some(imgs) = self.my_images.borrow().get(a_s) {
                    for &img_sr in imgs {
                        if a_m_fence.insert(img_sr.index) {
                            self.add_to_result(img_sr, the_type, result, &mut *t);
                        }
                    }
                }
            } else if the_type == rcad_kernel::topods::ShapeType::Face {
                // Path B — TopExp_Explorer on Face sub-shapes of this argument.
                // rcad: DS data replaces explorer; side maps argument index to origin.
                let side = if arg_idx == 0 {
                    ShapeOrigin::ShapeA
                } else {
                    ShapeOrigin::ShapeB
                };
                let f_base = self.ds.vertex_count() + self.ds.edge_count() + self.ds.shape_info.iter().filter(|si| si.shape_type == rcad_kernel::topods::ShapeType::Wire && !si.is_new).count();
                let side_offset = if side == ShapeOrigin::ShapeA {
                    0usize
                } else {
                    self.ds.a_face_count
                };
                for (fi, df) in self.ds.faces.iter().enumerate() {
                    if df.origin != side {
                        continue;
                    }
                    let src_sr = self.brep_sr(f_base + side_offset + df.source_face_idx);
                    let has_images = self.my_images.borrow().contains_key(&src_sr);
                    if !has_images {
                        // OCCT BuildResult: no images → add original shape as-is.
                        // BRep_Builder().Add(myShape, aS) adds the original TopoDS_Face.
                        // rcad: source face TShape at flat index src_sr.index.
                        if a_m_fence.insert(src_sr.index) {
                            self.add_to_result(src_sr, the_type, result, &mut *t);
                        }
                    } else if let Some(imgs) = self.my_images.borrow().get(&src_sr) {
                        for &img_sr in imgs {
                            if a_m_fence.insert(img_sr.index) {
                                self.add_to_result(img_sr, the_type, result, &mut *t);
                            }
                        }
                    }
                }
            } else if the_type == rcad_kernel::topods::ShapeType::Vertex {
                // Path B - all source DS vertices (OCCT: TopExp_Explorer VERTEX).
                for vi in 0..self.ds.vertex_count() {
                    let src_sr = self.brep_sr(vi);
                    let has_images = self.my_images.borrow().contains_key(&src_sr);
                    if !has_images {
                        if a_m_fence.insert(src_sr.index) {
                            self.add_to_result(src_sr, the_type, result, &mut *t);
                        }
                    } else if let Some(imgs) = self.my_images.borrow().get(&src_sr) {
                        for &img_sr in imgs {
                            if a_m_fence.insert(img_sr.index) {
                                self.add_to_result(img_sr, the_type, result, &mut *t);
                            }
                        }
                    }
                }
            } else if the_type == rcad_kernel::topods::ShapeType::Edge {
                // Path B - all source DS edges (OCCT: TopExp_Explorer EDGE).
                let e_base = self.ds.vertex_count();
                for ei in 0..self.ds.edge_count() {
                    let src_sr = self.brep_sr(e_base + ei);
                    let has_images = self.my_images.borrow().contains_key(&src_sr);
                    if !has_images {
                        if a_m_fence.insert(src_sr.index) {
                            self.add_to_result(src_sr, the_type, result, &mut *t);
                        }
                    } else if let Some(imgs) = self.my_images.borrow().get(&src_sr) {
                        for &img_sr in imgs {
                            if a_m_fence.insert(img_sr.index) {
                                self.add_to_result(img_sr, the_type, result, &mut *t);
                            }
                        }
                    }
                }
            } else if the_type == rcad_kernel::topods::ShapeType::Wire {
                // Path B - all source DS wires (OCCT: TopExp_Explorer WIRE).
                let w_base = self.ds.vertex_count() + self.ds.edge_count();
                for wi in 0..self.ds.shape_info.iter().filter(|si| si.shape_type == rcad_kernel::topods::ShapeType::Wire && !si.is_new).count() {
                    let src_sr = self.brep_sr(w_base + wi);
                    let has_images = self.my_images.borrow().contains_key(&src_sr);
                    if !has_images {
                        if a_m_fence.insert(src_sr.index) {
                            self.add_to_result(src_sr, the_type, result, &mut *t);
                        }
                    } else if let Some(imgs) = self.my_images.borrow().get(&src_sr) {
                        for &img_sr in imgs {
                            if a_m_fence.insert(img_sr.index) {
                                self.add_to_result(img_sr, the_type, result, &mut *t);
                            }
                        }
                    }
                }
            }
            // Path B for other types (e.g. Edge inside Solid) — no sub-shape
            // explorer implemented; the requirement is handled upstream.
        }

        // OCCT BuildResult(TopAbs_SOLID): builds result solids from split solids
        // created by BuildSplitSolids. rcad: process result.tmp_solids.
        if the_type == rcad_kernel::topods::ShapeType::Solid {
            let solids = std::mem::take(&mut result.tmp_solids);
            let mfr = self.my_face_refs.borrow();
            for shell_faces in &solids {
                let mut wire_refs: Vec<topods::ShapeRef> = Vec::new();
                for &fi in shell_faces {
                    if let Some(&face_sr) = mfr.get(fi) {
                        if face_sr.is_null() {
                            continue;
                        }
                        let face_tshape = &t.tshapes[face_sr.index];
                        if let topods::TShape::Face(fd) = &**face_tshape {
                            wire_refs.push(fd.outer_wire);
                        }
                    }
                }
                if !wire_refs.is_empty() {
                    let shell_sr = t.add_tshell(wire_refs);
                    t.add_tsolid(vec![shell_sr]);
                }
            }
            // OCCT BuildResult(SOLID): also add non-interfered solids from my_solids.
            for &solid_sr in self.my_solids.borrow().iter() {
                if a_m_fence.insert(solid_sr.index) {
                    self.add_to_result(solid_sr, the_type, result, &mut *t);
                }
            }
        }
    }

    /// add a ShapeRef to result (equivalent to BRep_Builder.Add(myShape, aS)).
    ///   Handles ALL shape types --populates ResultBuilder data structures.
    ///   Uses ensure_vertex_at/ensure_edge_at to place TShapes at flat indices
    ///   matching OCCT TopoDS_Shape identity (architecture alignment).
    pub(super) fn add_to_result(
        &self,
        a_s: rcad_kernel::topods::ShapeRef,
        topods_type: rcad_kernel::topods::ShapeType,
        result: &mut ResultBuilder,
        t: &mut rcad_kernel::topods::BRep,
    ) {
        match topods_type {
            rcad_kernel::topods::ShapeType::Vertex => {
                let vi = a_s.index;
                if vi < self.ds.vertex_count() {
                    // Create vertex TShape at flat index vi
                    t.ensure_vertex_at(vi, self.ds.vertex_point(vi));
                    result.add_ds_vertex(vi, self.ds.vertex_point(vi));
                }
            }
            rcad_kernel::topods::ShapeType::Edge => {
                let e_base = self.ds.vertex_count();
                let ei = a_s.index.saturating_sub(e_base);
                if ei < self.ds.edge_count() {
                    let edge = &self.ds.edges[ei];
                    // Ensure vertex TShapes exist at correct flat indices
                    let sv_sr = t.ensure_vertex_at(
                        edge.start_vertex,
                        self.ds.vertex_point(edge.start_vertex),
                    );
                    let ev_sr = t
                        .ensure_vertex_at(edge.end_vertex, self.ds.vertex_point(edge.end_vertex));
                    // Create edge TShape at flat index e_base + ei
                    let te = t.ensure_edge_at(
                        e_base + ei,
                        Some(edge.curve.clone()),
                        sv_sr,
                        ev_sr,
                        edge.t_range,
                    );
                    self.my_edge_map.borrow_mut()[ei] = te;
                }
            }
            rcad_kernel::topods::ShapeType::Wire => {
                let e_base = self.ds.vertex_count();
                let wi = a_s.index.saturating_sub(e_base + self.ds.edge_count());
                if wi < self.ds.shape_info.iter().filter(|si| si.shape_type == rcad_kernel::topods::ShapeType::Wire && !si.is_new).count() {
                    let w_ref = self.brep_sr(e_base + self.ds.edge_count() + wi);
                    let mut wire_edges = Vec::new();
                    let e_map = self.my_edge_map.borrow();
                    if let Some(imgs) = self.my_images.borrow().get(&w_ref) {
                        for &img_sr in imgs {
                            let nSpR = img_sr.index.saturating_sub(e_base);
                            if nSpR < e_map.len()
                                && e_map[nSpR] != rcad_kernel::topods::ShapeRef::NULL
                            {
                                wire_edges.push(e_map[nSpR]);
                            }
                        }
                    } else {
                        for &sub in &self.ds.shape_info[wi].sub_shapes {
                            // Convert flat shape index to per-type edge index
                            if sub >= self.ds.a_vertex_count {
                                let ei = sub - self.ds.a_vertex_count;
                                if ei < e_map.len() && e_map[ei] != rcad_kernel::topods::ShapeRef::NULL
                                {
                                    wire_edges.push(e_map[ei]);
                                }
                            }
                        }
                    }
                    drop(e_map);
                    if !wire_edges.is_empty() {
                        let sr = t.add_twire(wire_edges);
                        t.wire_mut(sr).flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                        self.my_wire_refs.borrow_mut().push(sr);
                    } else {
                        self.my_wire_refs
                            .borrow_mut()
                            .push(rcad_kernel::topods::ShapeRef::NULL);
                    }
                }
            }
            rcad_kernel::topods::ShapeType::Face => {
                // OCCT BuildResult: BRep_Builder().Add(myShape, aS) — adds original
                // face TShape. rcad: source face TShape created by
                // populate_source_shapes_in_t_brep at flat index a_s.index.
                // Ensure the face TShape exists in my_shape.
                if a_s.index >= t.tshapes.len() {
                    return;
                }
                if !self
                    .my_face_refs
                    .borrow()
                    .iter()
                    .any(|&r| !r.is_null() && r.index == a_s.index)
                {
                    self.my_face_refs.borrow_mut().push(a_s);
                }
            }
            rcad_kernel::topods::ShapeType::Shell => {
                if !self.my_shells.borrow().contains(&a_s) {
                    self.my_shells.borrow_mut().push(a_s);
                }
            }
            rcad_kernel::topods::ShapeType::Solid => {
                if !self.my_solids.borrow().contains(&a_s) {
                    self.my_solids.borrow_mut().push(a_s);
                }
            }
            rcad_kernel::topods::ShapeType::CompSolid
            | rcad_kernel::topods::ShapeType::Compound => {}
            rcad_kernel::topods::ShapeType::Shape => unreachable!(),
        }
    }

    /// ✅ OCCT-aligned: BuildBOP (BOPAlgo_BOP.cxx L885-926, state-based); ref: BOPAlgo_Builder.cxx L490-895.
    pub(super) fn build_bop(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        if self.has_errors() {
            return;
        }
        let my_face_refs = self.my_face_refs.borrow();
        let my_images = self.my_images.borrow();
        let my_in_parts = self.my_in_parts.borrow();

        let has_objects = self.ds.a_face_count > 0;
        let has_tools = self.ds.face_count() > self.ds.a_face_count;
        if !has_objects && !has_tools {
            self.my_report
                .borrow_mut()
                .add_alert(crate::bopalgo::Alert::TooFewArguments);
            return;
        }
        // State validation
        eprintln!(
            "[DBG_BOP] build_bop: n_faces_refs={} n_in_parts={}",
            my_face_refs.len(),
            my_in_parts.len()
        );
        let the_obj_state_in = matches!(self.my_operation, BooleanOpType::Intersection);
        let the_tools_state_in = matches!(
            self.my_operation,
            BooleanOpType::Intersection | BooleanOpType::Difference
        );
        // Build face maps per source solid.
        // rcad: iterate DS faces per side (A=0, B=1), collect face images with orientation.
        let f_base = self.ds.vertex_count() + self.ds.edge_count();
        let mut faces_with_ori: [Vec<topods::ShapeRef>; 2] = [Vec::new(), Vec::new()];
        let mut faces_fence: [std::collections::HashSet<u64>; 2] = [
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        ];
        let mut in_faces: [std::collections::HashSet<u64>; 2] = [
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        ];
        // Collect source solid indices per side for myInParts lookup
        let mut src_solids: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
        for (si, _sol_si) in self.ds.solid_shape_indices().iter().enumerate() {
            // Determine side from any face belonging to this solid
            let mut side = 2usize; // unknown
            for (fi, df) in self.ds.faces.iter().enumerate() {
                if df.source_solid_idx == Some(si) {
                    side = if df.origin == ShapeOrigin::ShapeA {
                        0
                    } else {
                        1
                    };
                    break;
                }
            }
            if side < 2 {
                src_solids[side].push(si);
            }
        }

        for side in 0..2 {
            // Collect all DS face indices for this side
            let side_origin = if side == 0 {
                ShapeOrigin::ShapeA
            } else {
                ShapeOrigin::ShapeB
            };
            for &src_si in &src_solids[side] {
                let mut solid_faces: Vec<usize> = Vec::new();
                for (fi, df) in self.ds.faces.iter().enumerate() {
                    if df.origin == side_origin && df.source_solid_idx == Some(src_si) {
                        solid_faces.push(fi);
                    }
                }
                for &fi in &solid_faces {
                    let face_sr = my_face_refs
                        .get(fi)
                        .copied()
                        .unwrap_or(topods::ShapeRef::synthetic(f_base + fi));
                    if face_sr.is_null() {
                        continue;
                    }

                    let has_images = my_images
                        .get(&face_sr)
                        .map_or(false, |imgs| !imgs.is_empty());
                    if has_images {
                        if let Some(imgs) = my_images.get(&face_sr) {
                            for &img_sr in imgs {
                                // IsSplitToReverse check
                                let need_reverse = if img_sr.index < t_brep.tshapes.len() {
                                    let orig_normal = self.face_approx_normal(fi);
                                    let split_normal = self.shape_ref_face_normal(img_sr, t_brep);
                                    split_normal.map_or(false, |sn| {
                                        crate::boptools::is_split_to_reverse(orig_normal, sn)
                                    })
                                } else {
                                    false
                                };
                                if need_reverse {
                                    let mut rev = img_sr;
                                    rev.orientation = topods::Orientation::Reversed;
                                    faces_with_ori[side].push(rev);
                                } else {
                                    faces_with_ori[side].push(img_sr);
                                }
                                faces_fence[side].insert(img_sr.ptr_id);
                            }
                        }
                    } else {
                        faces_with_ori[side].push(face_sr);
                        faces_fence[side].insert(face_sr.ptr_id);
                    }
                }
                // Collect IN faces from myInParts for this solid
                if let Some(in_face_indices) = my_in_parts.get(&src_si) {
                    for &in_fi in in_face_indices {
                        let in_sr = my_face_refs
                            .get(in_fi)
                            .copied()
                            .unwrap_or(topods::ShapeRef::synthetic(f_base + in_fi));
                        if !in_sr.is_null() {
                            in_faces[side].insert(in_sr.ptr_id);
                        }
                    }
                }
            }
        }

        // Face selection based on IN/OUT states
        let is_objects_in = the_obj_state_in;
        let is_tools_in = the_tools_state_in;
        let b_avoid_in = !is_objects_in && !is_tools_in;
        let b_avoid_in_for_both = is_objects_in != is_tools_in;
        let is_same_ori_needed = the_obj_state_in == the_tools_state_in;

        let mut a_m_res_faces_ori: Vec<topods::ShapeRef> = Vec::new();
        let mut a_m_res_faces_fence: std::collections::HashSet<u64> =
            std::collections::HashSet::new();
        let mut a_m_f_to_avoid: std::collections::HashSet<u64> = std::collections::HashSet::new();
        let mut a_m_fence: [std::collections::HashSet<u64>; 2] = [
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        ];
        let mut a_m_fence_ori: [std::collections::HashSet<u64>; 2] = [
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        ];

        for side in 0..2 {
            let b_take_in = if side == 0 {
                is_objects_in
            } else {
                is_tools_in
            };
            let opposite_side = 1 - side;
            for &f_sr in &faces_with_ori[side] {
                let f_id = f_sr.ptr_id;
                let is_in = in_faces[side].contains(&f_id);
                let is_in_opposite = in_faces[opposite_side].contains(&f_id);
                if b_avoid_in && (is_in || is_in_opposite) {
                    continue;
                }
                if b_avoid_in_for_both && is_in && is_in_opposite {
                    continue;
                }
                if !a_m_fence[side].insert(f_id) {
                    if !faces_fence[opposite_side].contains(&f_id) {
                        if b_take_in != is_same_ori_needed {
                            a_m_f_to_avoid.insert(f_id);
                        }
                    } else {
                        let is_same_ori = !a_m_fence_ori[side].insert(f_id);
                        if is_same_ori_needed == is_same_ori {
                            if a_m_res_faces_fence.insert(f_id) {
                                a_m_res_faces_ori.push(f_sr);
                            }
                        } else {
                            a_m_f_to_avoid.insert(f_id);
                        }
                        continue;
                    }
                }
                if !a_m_fence_ori[side].insert(f_id) {
                    continue;
                }
                if b_take_in == is_in_opposite {
                    if is_in {
                        a_m_res_faces_ori.push(f_sr);
                        let mut rev = f_sr;
                        rev.orientation = topods::Orientation::Reversed;
                        a_m_res_faces_ori.push(rev);
                    } else if b_take_in && !is_same_ori_needed {
                        let mut rev = f_sr;
                        rev.orientation = topods::Orientation::Reversed;
                        a_m_res_faces_ori.push(rev);
                    } else {
                        a_m_res_faces_ori.push(f_sr);
                    }
                    a_m_res_faces_fence.insert(f_id);
                }
            }
        }
        // Remove avoided faces
        let mut a_res_faces: Vec<topods::ShapeRef> = Vec::new();
        for &f_sr in &a_m_res_faces_ori {
            if !a_m_f_to_avoid.contains(&f_sr.ptr_id) {
                a_res_faces.push(f_sr);
            }
        }
        if a_res_faces.is_empty() {
            self.my_report
                .borrow_mut()
                .add_alert(crate::bopalgo::Alert::TooFewArguments);
            return;
        }
        // BuilderSolid from selected faces
        let mut bs_faces: Vec<usize> = Vec::new();
        for f_sr in &a_res_faces {
            for (dfi, &ref_sr) in my_face_refs.iter().enumerate() {
                if !ref_sr.is_null() && ref_sr.ptr_id == f_sr.ptr_id {
                    bs_faces.push(dfi);
                    break;
                }
            }
        }
        bs_faces.sort_unstable();
        bs_faces.dedup();

        let mut a_bs = crate::bopds::builder_solid::BuilderSolid::new();
        a_bs.set_shapes(&bs_faces);
        a_bs.perform(&self.ds);
        // Validate each resulting solid has ≥1 face from objects/tools
        let mut a_res_solids: Vec<topods::ShapeRef> = Vec::new();
        if !self.ds.solid_shape_indices().is_empty() || true {
            let a_areas: &[Vec<usize>] = a_bs.areas();
            for area_faces in a_areas {
                if area_faces.is_empty() {
                    continue;
                }
                let has_obj_or_tool_face = area_faces.iter().any(|&fi| {
                    let fi_sr = my_face_refs
                        .get(fi)
                        .copied()
                        .unwrap_or(topods::ShapeRef::synthetic(f_base + fi));
                    faces_fence[0].contains(&fi_sr.ptr_id) || faces_fence[1].contains(&fi_sr.ptr_id)
                });
                if !has_obj_or_tool_face {
                    continue;
                }
                let mut area_face_refs: Vec<topods::ShapeRef> = Vec::new();
                for &fi in area_faces {
                    let fi_sr = my_face_refs
                        .get(fi)
                        .copied()
                        .unwrap_or(topods::ShapeRef::synthetic(f_base + fi));
                    area_face_refs.push(fi_sr);
                }
                let shell_ref = t_brep.add_tshell(area_face_refs);
                let solid_ref = t_brep.add_tsolid(vec![shell_ref]);
                a_res_solids.push(solid_ref);
            }
        }
        // Collect unused faces → solids
        {
            let mut a_unused_fence: std::collections::HashSet<u64> =
                std::collections::HashSet::new();
            for solid_sr in &a_res_solids {
                if solid_sr.index < t_brep.tshapes.len() {
                    if let topods::TShape::Solid(ref sd) = *t_brep.tshapes[solid_sr.index] {
                        for &sh_sr in &sd.shells {
                            if sh_sr.index < t_brep.tshapes.len() {
                                if let topods::TShape::Shell(ref shd) = *t_brep.tshapes[sh_sr.index]
                                {
                                    for &face_sr in &shd.faces {
                                        a_unused_fence.insert(face_sr.ptr_id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            let mut a_unused_faces: Vec<topods::ShapeRef> = Vec::new();
            for &f_sr in &a_res_faces {
                if !a_unused_fence.contains(&f_sr.ptr_id) {
                    a_unused_faces.push(f_sr);
                }
            }
            if !a_unused_faces.is_empty() {
                let shell_ref = t_brep.add_tshell(a_unused_faces);
                let solid_ref = t_brep.add_tsolid(vec![shell_ref]);
                a_res_solids.push(solid_ref);
            }
        }
        // Combine solids into result compound
        for &sr in &a_res_solids {
            self.my_solids.borrow_mut().push(sr);
        }
    }

    /// Approximation of face normal for orientation comparison (OCCT IsSplitToReverse).
    /// Uses the first triangle or edges to estimate face normal direction.
    fn face_approx_normal(&self, fi: usize) -> glam::DVec3 {
        if let Some(face) = self.ds.faces.get(fi) {
            return face.surface.normal_at(0.5, 0.5);
        }
        glam::DVec3::Z
    }

    /// Lookup the normal of a ShapeRef that refers to a Face TShape in the BRep.
    fn shape_ref_face_normal(
        &self,
        sr: topods::ShapeRef,
        t_brep: &topods::BRep,
    ) -> Option<glam::DVec3> {
        if sr.index < t_brep.tshapes.len() {
            if let topods::TShape::Face(fd) = &*t_brep.tshapes[sr.index] {
                if let Some(ref surf) = fd.surface {
                    return Some(surf.normal_at(0.5, 0.5));
                }
            }
        }
        None
    }

    /// ✅ OCCT-aligned: CheckArgsForOpenSolid (BOPAlgo_BOP.cxx L1396-1470+).
    pub(super) fn check_args_for_open_solid(&self) -> bool {
        for (_sh_si, _shell_faces) in &self.ds.collect_source_shells() {
            // Check if the shell is closed (each edge appears twice).
            let mut ecount: std::collections::HashMap<usize, usize> =
                std::collections::HashMap::new();
            for &fi in _shell_faces {
                if let Some(face) = self.ds.faces.get(fi) {
                    for &ei in &face.boundary_edges {
                        *ecount.entry(ei).or_default() += 1;
                    }
                }
            }
            if ecount.values().any(|&c| c != 2) {
                return true; // non-closed shell found
            }
        }
        false // all shells closed
    }
}
