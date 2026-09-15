//! OCCT BRepBuilderAPI_Sewing.cxx — the cutting group:
//! `Cutting` (L4441-4602), `ProjectPointsOnCurve` (L5367-5478),
//! `CreateCuttingNodes` (L5480-5667) and `CreateSections` (L5669-5880).
//!
//! Re-hosts:
//! - `GeomAdaptor_Curve` + `Extrema_ExtPC` — the kernel `GeomCurveAdaptor`
//!   / `CurveToolHandle` / `ExtremaExtPC` bodies.
//! - `BndLib_Add3dCurve::Add(adptC, tol, box)` (Cutting L4501) — the 3d
//!   bounding box of the curve subrange; the rcad re-host samples the
//!   curve over the range (the pave_filler.rs `shrunk_range_bnd_box`
//!   precedent).

use glam::DVec3;
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::geom::{Curve2d, CurveEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation};

use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

use super::{idx_shape_get, make_vertex_at, set_add};
use super::BRepBuilderAPISewing;
use super::selectors::{BndBoxTree, BndBoxTreeFiller, BndBoxTreeSelector};

/// OCCT `BndLib_Add3dCurve::Add(adptC, tol, aGlobalBox)` (BndLib_Add3dCurve.cxx
/// L29-36) — the 3d bounding box of the curve subrange over [first, last];
/// the rcad re-host samples the curve (documented sampling reduction).
fn bnd_lib_add3d_curve(
    the_c: &rcad_kernel::geom::Curve3,
    first: f64,
    last: f64,
    tol: f64,
    b: &mut BndBox,
) {
    use rcad_kernel::geom::CurveEval;
    let n = 16usize;
    for i in 0..=n {
        let t = first + (last - first) * (i as f64) / (n as f64);
        let p = the_c.point_at(t);
        b.update(p.x, p.y, p.z, p.x, p.y, p.z);
    }
    b.enlarge(tol);
}

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::Cutting(theProgress) (cxx L4441-4602) —
    /// modifies myBoundSections, myNodeSections, myCuttingNode.
    pub(crate) fn cutting(&mut self, brep: &mut BRep) {
        // OCCT L4442-4447.
        let nb_vertices = self.my_vertex_node.len();
        if nb_vertices == 0 {
            return;
        }
        // Create a box tree with vertices
        // OCCT L4449-4458.
        let eps = self.my_tolerance * 0.5;
        let a_tree = BndBoxTree::new();
        let mut a_tree_filler = BndBoxTreeFiller::new(a_tree);
        let mut a_selector = BndBoxTreeSelector::new();
        for i in 1..=nb_vertices {
            let (node, _) = idx_shape_get(&self.my_vertex_node, i - 1);
            let pt = bat::brep_tool_pnt(&node).unwrap_or(DVec3::ZERO);
            let mut a_box = BndBox::from_point(pt);
            a_box.enlarge(eps);
            a_tree_filler.add(i as i32, a_box);
        }
        a_tree_filler.fill();

        // OCCT L4465-4468.
        // Iterate on all boundaries
        // OCCT L4470-4474.
        let nb_bounds = self.my_bound_faces.len();
        for idx in 0..nb_bounds {
            let (bound_v, list_faces) = super::idx_list_get(&self.my_bound_faces, idx);
            let bound = bound_v;
            // Do not cut floating edges
            // OCCT L4477-4480.
            if list_faces.is_empty() {
                continue;
            }
            // Obtain bound curve
            // OCCT L4482-4492.
            let (c3d, first, last) = match super::brep_tool_curve_world(brep, &bound) {
                Some(v) => v,
                None => continue,
            };
            // Create cutting sections
            // OCCT L4495-4586.
            let mut list_sections: Vec<Shape> = Vec::new();
            {
                // Obtain candidate vertices
                // OCCT L4498-4531.
                let mut candidate_vertices: super::ShapeSet = super::ShapeSet::new();
                {
                    // Create bounding box around curve
                    // OCCT L4500-4502.
                    let mut a_global_box = BndBox::new();
                    bnd_lib_add3d_curve(&c3d, first, last, self.my_tolerance, &mut a_global_box);
                    // Sort vertices to find candidates
                    // OCCT L4504-4506.
                    a_selector.set_current(a_global_box);
                    a_tree_filler.select(&mut a_selector);
                    // Skip bound if no node is in the boundind box
                    // OCCT L4507-4510.
                    if a_selector.res_ind().is_empty() {
                        continue;
                    }
                    // Retrieve bound nodes
                    // OCCT L4512-4515.
                    let (v1, v2) = bat::top_exp_vertices_raw(&bound);
                    let v1 = v1.unwrap_or_else(Shape::null);
                    let v2 = v2.unwrap_or_else(Shape::null);
                    let node1 = match self.my_vertex_node.get(&bat::shape_key(&v1)) {
                        Some((_, n)) => n.clone(),
                        None => continue,
                    };
                    let node2 = match self.my_vertex_node.get(&bat::shape_key(&v2)) {
                        Some((_, n)) => n.clone(),
                        None => continue,
                    };
                    // Fill map of candidate vertices
                    // OCCT L4517-4527.
                    for itl in a_selector.res_ind().clone() {
                        let index = itl;
                        let (vertex, node) = idx_shape_get(&self.my_vertex_node, (index - 1) as usize);
                        if !node.is_same(&node1) && !node.is_same(&node2) {
                            set_add(&mut candidate_vertices, &vertex);
                        }
                    }
                    a_selector.clear_res_list();
                    // OCCT: the V1/V2 locals live past the bracket scope —
                    // carried below through the explicit locals.
                    if !set_candidate_guard(&candidate_vertices) {
                        continue;
                    }
                    // Project vertices on curve
                    // OCCT L4533-4544.
                    let nb_candidates = candidate_vertices.len();
                    let mut arr_para = vec![0.0f64; nb_candidates];
                    let mut arr_dist = vec![0.0f64; nb_candidates];
                    let mut arr_pnt = vec![DVec3::ZERO; nb_candidates];
                    let mut arr_proj = vec![DVec3::ZERO; nb_candidates];
                    for j in 1..=nb_candidates {
                        let (cand, _) = super::set_entry(&candidate_vertices, j - 1);
                        arr_pnt[j - 1] = bat::brep_tool_pnt(&cand).unwrap_or(DVec3::ZERO);
                    }
                    self.project_points_on_curve(
                        &arr_pnt, &c3d, first, last, &mut arr_dist, &mut arr_para, &mut arr_proj,
                        true,
                    );
                    // Create cutting nodes
                    // OCCT L4547-4555.
                    let mut seq_node: Vec<Shape> = Vec::new();
                    let mut seq_para: Vec<f64> = Vec::new();
                    let (v1, v2) = {
                        let (a, b) = bat::top_exp_vertices_raw(&bound);
                        (
                            a.unwrap_or_else(Shape::null),
                            b.unwrap_or_else(Shape::null),
                        )
                    };
                    self.create_cutting_nodes(
                        brep,
                        &candidate_vertices,
                        &bound,
                        &v1,
                        &v2,
                        &arr_dist,
                        &arr_para,
                        // OCCT L4544-4552: the formal parameter is named arrPnt
                        // but the call passes arrProj — the PROJECTED points
                        // on the bound curve (they drive the closest-vertex
                        // search L5570-5577, the closeness test L5582 and the
                        // created cutting-vertex position L5603).
                        &arr_proj,
                        &mut seq_node,
                        &mut seq_para,
                    );
                    // OCCT L4556-4559.
                    if seq_para.is_empty() {
                        continue;
                    }
                    // Create cutting sections
                    self.create_sections(brep, &bound, &seq_node, &seq_para, &mut list_sections);
                }
            }
            // OCCT L4587-4610.
            if list_sections.len() > 1 {
                // modification of maps: myBoundSections
                for its in &list_sections {
                    let section = its.clone();
                    // Iterate on section vertices
                    for itv in bat::sub_shapes(&section) {
                        let mut vertex = itv;
                        // Convert vertex to node
                        if let Some((_, n)) = self.my_vertex_node.get(&bat::shape_key(&vertex)) {
                            vertex = n.clone();
                        }
                        // Update node sections
                        match self.my_node_sections.get_mut(&bat::shape_key(&vertex)) {
                            Some(list) => list.push(section.clone()),
                            None => {
                                self.my_node_sections
                                    .insert(bat::shape_key(&vertex), vec![section.clone()]);
                            }
                        }
                    }
                    // Store bound for section
                    self.my_section_bound
                        .insert(bat::shape_key(&section), bound.clone());
                }
                // Store split bound
                self.my_bound_sections
                    .insert(bat::shape_key(&bound), list_sections.clone());
            }
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::ProjectPointsOnCurve(arrPnt, c3d, first,
    /// last, arrDist, arrPara, arrProj, isConsiderEnds) (cxx L5367-5478) —
    /// projects points on curve.
    pub(crate) fn project_points_on_curve(
        &self,
        arr_pnt: &[DVec3],
        c3d: &rcad_kernel::geom::Curve3,
        first: f64,
        last: f64,
        arr_dist: &mut [f64],
        arr_para: &mut [f64],
        arr_proj: &mut [DVec3],
        is_consider_ends: bool,
    ) {
        // OCCT L5370: arrDist.Init(-1.0).
        for v in arr_dist.iter_mut() {
            *v = -1.0;
        }

        // OCCT L5373-5376.
        let pfirst = c3d.point_at(first);
        let plast = c3d.point_at(last);
        let find = 1usize; // (isConsiderEnds ? 1 : 2);
        let lind = arr_pnt.len(); // (isConsiderEnds ? ... : ... - 1);

        // OCCT L5379-5381.
        for i1 in find..=lind {
            let pt = arr_pnt[i1 - 1];
            // OCCT L5390: double worktol = myTolerance.  INFO (cxx
            // L5447-5456): OCCT wraps the projection in try/catch and on a
            // Standard_Failure from Perform or the result accessors
            // degrades worktol to MinTolerance() for the end-point
            // fallback below.  The rcad carrier (ExtremaExtPC) surfaces
            // those same conditions as panics, not catchable results
            // (rcad-kernel extrema_ext_pc.rs mirrors Standard_OutOfRange /
            // StdFail_NotDone as panic!), so the catch-and-degrade path is
            // unrepresentable here and worktol stays at myTolerance —
            // no invented recovery.
            let mut worktol = self.my_tolerance;
            let dist_f2 = pfirst.distance_squared(pt);
            let dist_l2 = plast.distance_squared(pt);
            let mut is_projected = false;

            // Project current point on curve
            // OCCT L5388-5443.
            {
                // OCCT L5390-5392: locProj.Initialize(GAC, first, last).
                let adaptor =
                    rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor::new(
                        c3d.clone(),
                    );
                let gac = rcad_kernel::base::extrema_curve_tool::CurveToolHandle::for_curve3(
                    c3d, &adaptor, &adaptor,
                );
                let mut loc_proj =
                    rcad_kernel::base::extrema_ext_pc::ExtremaExtPC::new_point_curve_ranged(
                        pt, &gac, first, last, 1.0e-10,
                    );
                if loc_proj.is_done() && loc_proj.nb_ext() > 0 {
                    let dist2_min = if is_consider_ends || i1 == find || i1 == lind {
                        dist_f2.min(dist_l2)
                    } else {
                        rcad_kernel::core::precision::INFINITE_VALUE
                    };
                    let mut ind_min: usize = 0;
                    let mut dist2_min = dist2_min;
                    for ind in 1..=loc_proj.nb_ext() {
                        let d_proj2 = loc_proj.square_distance(ind);
                        if d_proj2 < dist2_min {
                            ind_min = ind;
                            dist2_min = d_proj2;
                        }
                    }
                    if ind_min != 0 {
                        is_projected = true;
                        let p_on_c = loc_proj.point(ind_min);
                        let mut param_proj = p_on_c.param;
                        let mut pt_proj = c3d.point_at(param_proj);
                        let mut dist_proj2 = pt_proj.distance_squared(pt);
                        if !loc_proj.is_min(ind_min) {
                            if dist_f2.min(dist_l2) < dist2_min {
                                if dist_f2 < dist_l2 {
                                    param_proj = first;
                                    dist_proj2 = dist_f2;
                                    pt_proj = pfirst;
                                } else {
                                    param_proj = last;
                                    dist_proj2 = dist_l2;
                                    pt_proj = plast;
                                }
                            }
                        }
                        // OCCT L5425-5430.
                        if dist_proj2 < worktol * worktol || !is_consider_ends {
                            arr_dist[i1 - 1] = dist_proj2.sqrt();
                            arr_para[i1 - 1] = param_proj;
                            arr_proj[i1 - 1] = pt_proj;
                        }
                    }
                }
            }
            // OCCT L5445-5463.
            // INFO: the fallback threshold uses worktol, which OCCT may
            // have degraded to MinTolerance() by the L5447-5456 catch —
            // see the worktol declaration above for why rcad keeps it at
            // myTolerance (the carrier panics instead of throwing).
            if !is_projected && is_consider_ends {
                if dist_f2.min(dist_l2) < worktol * worktol {
                    if dist_f2 < dist_l2 {
                        arr_dist[i1 - 1] = dist_f2.sqrt();
                        arr_para[i1 - 1] = first;
                        arr_proj[i1 - 1] = pfirst;
                    } else {
                        arr_dist[i1 - 1] = dist_l2.sqrt();
                        arr_para[i1 - 1] = last;
                        arr_proj[i1 - 1] = plast;
                    }
                }
            }
            let _ = &mut worktol;
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::CreateCuttingNodes(MapVert, bound, vfirst,
    /// vlast, arrDist, arrPara, arrPnt, seqVert, seqPara) (cxx L5480-5667) —
    /// creates cutting vertices on projections.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn create_cutting_nodes(
        &mut self,
        brep: &mut BRep,
        map_vert: &super::ShapeSet,
        bound: &Shape,
        vfirst: &Shape,
        vlast: &Shape,
        arr_dist: &[f64],
        arr_para: &[f64],
        arr_pnt: &[DVec3],
        seq_vert: &mut Vec<Shape>,
        seq_para: &mut Vec<f64>,
    ) {
        // OCCT L5482-5483.
        let nb_proj = map_vert.len();

        // Reorder projections by distance
        // OCCT L5486-5504.
        let mut seq_ordered_index: Vec<usize> = Vec::new();
        {
            let mut seq_ordered_distance: Vec<f64> = Vec::new();
            for i in 1..=nb_proj {
                let dist_proj = arr_dist[i - 1];
                if dist_proj < 0.0 {
                    continue; // Skip vertex if not projected
                }
                let mut is_inserted = false;
                for j in 1..=seq_ordered_index.len() {
                    if is_inserted {
                        break;
                    }
                    is_inserted = dist_proj < seq_ordered_distance[j - 1];
                    if is_inserted {
                        seq_ordered_index.insert(j - 1, i - 1);
                        seq_ordered_distance.insert(j - 1, dist_proj);
                    }
                }
                if !is_inserted {
                    seq_ordered_index.push(i - 1);
                    seq_ordered_distance.push(dist_proj);
                }
            }
        }
        let nb_proj = seq_ordered_index.len();
        if nb_proj == 0 {
            return;
        }

        // OCCT L5511: BRep_Builder aBuilder.

        // Insert two initial vertices (to be removed later)
        // OCCT L5513-5531.
        let mut seq_dist: Vec<f64> = Vec::new();
        let mut seq_pnt: Vec<DVec3> = Vec::new();
        {
            // Retrieve bound curve
            let c3d = super::brep_tool_curve_world(brep, bound);
            if let Some((c3d, first, last)) = c3d {
                use rcad_kernel::geom::CurveEval;
                seq_vert.insert(0, vfirst.clone());
                seq_vert.push(vlast.clone());
                seq_para.insert(0, first);
                seq_para.push(last);
                seq_dist.insert(0, -1.0);
                seq_dist.push(-1.0);
                seq_pnt.insert(0, c3d.point_at(first));
                seq_pnt.push(c3d.point_at(last));
            }
        }

        // OCCT L5533-5625.
        let mut node_cutting_vertex: super::IdxShapeMap = super::IdxShapeMap::new();
        for i in 1..=nb_proj {
            let index = seq_ordered_index[i - 1];
            let dis_proj = arr_dist[index];
            let pnt_proj = arr_pnt[index];

            // Skip node if already bound to cutting vertex
            // OCCT L5537-5541.
            let (map_v, _) = super::set_entry(map_vert, index);
            let node = match self.my_vertex_node.get(&bat::shape_key(&map_v)) {
                Some((_, n)) => n.clone(),
                None => continue,
            };
            if node_cutting_vertex.contains_key(&bat::shape_key(&node)) {
                continue;
            }

            // Find the closest vertex
            // OCCT L5544-5552.
            let mut index_min = 1usize;
            let mut dist_min = if seq_pnt.is_empty() {
                f64::MAX
            } else {
                pnt_proj.distance(seq_pnt[0])
            };
            for j in 2..=seq_pnt.len() {
                let dist = pnt_proj.distance(seq_pnt[j - 1]);
                if dist < dist_min {
                    dist_min = dist;
                    index_min = j;
                }
            }

            // Check if current point is close to one of the existent
            // OCCT L5555-5572.
            if dist_min <= (dis_proj * 0.1).max(self.min_tolerance()) {
                // Check distance if close
                let jdist = seq_dist[index_min - 1];
                if jdist < 0.0 {
                    // Bind new cutting node (end vertex only)
                    seq_dist[index_min - 1] = dis_proj;
                    let cvertex = seq_vert[index_min - 1].clone();
                    super::idx_map_add(&mut node_cutting_vertex, &node, cvertex);
                } else {
                    // Bind secondary cutting nodes
                    super::idx_map_add(&mut node_cutting_vertex, &node, Shape::null());
                }
            } else {
                // Build new cutting vertex
                // OCCT L5573-5624.
                let cvertex = make_vertex_at(brep, pnt_proj, CONFUSION);
                // Bind new cutting vertex
                super::idx_map_add(&mut node_cutting_vertex, &node, cvertex.clone());
                // Insert cutting vertex in the sequences
                let par_proj = arr_para[index];
                for j in 2..=seq_para.len() {
                    if par_proj <= seq_para[j - 1] {
                        seq_vert.insert(j - 1, cvertex.clone());
                        seq_para.insert(j - 1, par_proj);
                        seq_dist.insert(j - 1, dis_proj);
                        seq_pnt.insert(j - 1, pnt_proj);
                        break;
                    }
                }
            }
        }

        // filling map for cutting nodes
        // OCCT L5628-5659.
        for i in 0..node_cutting_vertex.len() {
            let (node, cnode_v) = idx_shape_get(&node_cutting_vertex, i);
            let mut cnode = cnode_v;
            // Skip secondary nodes
            if cnode.is_null() {
                continue;
            }
            // Obtain vertex node
            if let Some((_, n)) = self.my_vertex_node.get(&bat::shape_key(&cnode)) {
                // This is an end vertex
                cnode = n.clone();
            } else {
                // Create link: cutting vertex -> node
                self.my_cutting_node
                    .insert(bat::shape_key(&cnode), vec![node.clone()]);
            }
            // Create link: node -> cutting vertex
            match self.my_cutting_node.get_mut(&bat::shape_key(&node)) {
                Some(list) => list.push(cnode.clone()),
                None => {
                    self.my_cutting_node
                        .insert(bat::shape_key(&node), vec![cnode.clone()]);
                }
            }
        }

        // Remove two initial vertices
        // OCCT L5662-5665.
        if !seq_vert.is_empty() {
            seq_vert.remove(0);
            seq_vert.pop();
        }
        if !seq_para.is_empty() {
            seq_para.remove(0);
            seq_para.pop();
        }
    }

    /// OCCT BRepBuilderAPI_Sewing::CreateSections(section, seqNode, seqPara,
    /// listEdge) (cxx L5669-5880) — performs cutting of bound.
    pub(crate) fn create_sections(
        &mut self,
        brep: &mut BRep,
        section: &Shape,
        seq_node: &[Shape],
        seq_para: &[f64],
        list_edge: &mut Vec<Shape>,
    ) {
        // OCCT L5670-5672.
        let sec = section;

        // To keep NM vertices on edge
        // OCCT L5675-5678.
        let mut a_seq_nm_vert: Vec<Shape> = Vec::new();
        let mut a_seq_nm_pars: Vec<f64> = Vec::new();
        super::same_param::find_nm_vertices(brep, sec, &mut a_seq_nm_vert, &mut a_seq_nm_pars);

        // OCCT L5681-5684.
        let (mut first, last) = bat::brep_tool_range(sec);

        // Create cutting sections
        // OCCT L5687-5726.
        let mut par1: f64;
        let mut par2: f64;
        let len = seq_para.len() + 1;
        for i in 1..=len {
            let mut edge = brep.empty_copied(sec);

            let (sec_first, sec_last) = bat::top_exp_vertices_raw(sec);
            let (v1, v2): (Shape, Shape);
            if i == 1 {
                par1 = first;
                par2 = seq_para[i - 1];
                v1 = sec_first.unwrap_or_else(Shape::null);
                v2 = seq_node[i - 1].clone();
            } else if i == len {
                par1 = seq_para[i - 2];
                par2 = last;
                v1 = seq_node[i - 2].clone();
                v2 = sec_last.unwrap_or_else(Shape::null);
            } else {
                par1 = seq_para[i - 2];
                par2 = seq_para[i - 1];
                v1 = seq_node[i - 2].clone();
                v2 = seq_node[i - 1].clone();
            }
            let _ = &mut first;

            // OCCT L5714-5718.
            let a_tmp_edge = bat::oriented(&edge, Orientation::Forward);
            let mut a_builder = BRepBuilder::new();
            a_builder.add_to_edge(brep, a_tmp_edge.clone(), bat::oriented(&v1, Orientation::Forward));
            a_builder.add_to_edge(brep, a_tmp_edge.clone(), bat::oriented(&v2, Orientation::Reversed));
            a_builder.set_edge_range(brep, a_tmp_edge.clone(), par1, par2);

            // OCCT L5721-5732.
            let mut k = 1usize;
            while k <= a_seq_nm_pars.len() {
                let apar = a_seq_nm_pars[k - 1];
                if apar >= par1 && apar <= par2 {
                    a_builder.add_to_edge(brep, a_tmp_edge.clone(), a_seq_nm_vert[k - 1].clone());
                    a_seq_nm_vert.remove(k - 1);
                    a_seq_nm_pars.remove(k - 1);
                } else {
                    k += 1;
                }
            }
            list_edge.push(edge.clone());
        }

        // OCCT L5735-5739.
        let list_faces = match self.my_bound_faces.get(&bat::shape_key(sec)) {
            Some((_, l)) => l.clone(),
            None => return,
        };
        if list_faces.is_empty() {
            return;
        }

        // OCCT L5741-5742.
        let tol_edge = bat::brep_tool_tolerance(sec);

        // Add cutting pcurves
        // OCCT L5745-5877.
        for itf in &list_faces {
            let fac = itf.clone();

            // Retrieve curve on surface
            // OCCT L5748-5752.
            let c2d = bat::brep_tool_curve_on_surface(sec, &fac);
            let (c2d, first2d, last2d) = match c2d {
                Some(v) => v,
                None => continue,
            };
            // OCCT L5755-5757.
            let mut c2d1: Option<Curve2d> = None;
            let is_seam = bat::brep_tool_is_closed_on_surface(sec, &fac);

            // gka fix for bug OCC12203 21.04.06 addition second curve for
            // seam edges
            // OCCT L5789-5800.
            if is_seam {
                let sec_rev = bat::reversed(sec);
                c2d1 = bat::brep_tool_curve_on_surface(&sec_rev, &fac)
                    .map(|(c, _f, _l)| c);
                if c2d1.is_none() {
                    continue;
                }
            }

            // Update cutting sections
            // OCCT L5813-5877.
            let mut a_builder = BRepBuilder::new();
            for ite in list_edge.clone() {
                // Retrieve cutting section
                // OCCT L5817-5819.
                let edge = ite.clone();
                let (par1e, par2e) = bat::brep_tool_range(&edge);
                let _ = (par1e, par2e);

                // Cut pcurve (the Segment() calls are commented in this
                // OCCT revision — the full-range copy is kept)
                // OCCT L5822: c2dNew = Copy().
                let c2d_new = c2d.clone();
                // OCCT L5824-5828.
                let mut c2d1_new: Option<Curve2d> = None;
                if let Some(c1) = &c2d1 {
                    c2d1_new = Some(c1.clone());
                }

                // OCCT L5844-5876.
                // INFO: in the else arm OCCT calls UpdateEdge with
                // (c2dNew, c2d1New) UNCONDITIONALLY (L5869/L5873); c2d1New
                // is non-null there because isSeam forces the L5797-5800
                // continue on a null c2d1, and the !isSeam state never
                // reaches the else arm (c2d1New is only set under isSeam) —
                // so the expect below never fires in reachable states.
                if !is_seam && c2d1_new.is_none() {
                    a_builder.update_edge_pcurve(
                        brep, edge.clone(), c2d_new.clone(), fac.clone(), tol_edge,
                    );
                } else {
                    let mut ori = edge.orientation;
                    if fac.orientation == Orientation::Reversed {
                        ori = super::top_abs_reverse(ori);
                    }

                    let c1n = c2d1_new.clone().expect(
                        "sewing CreateSections: null c2d1New in the seam UpdateEdge (unreachable; OCCT L5869/L5873 pass it unconditionally)",
                    );
                    if ori == Orientation::Forward {
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c2d_new.clone(), c1n, fac.clone(), tol_edge,
                        );
                    } else {
                        a_builder.update_edge_pcurve_closed(
                            brep, edge.clone(), c1n, c2d_new.clone(), fac.clone(), tol_edge,
                        );
                    }
                }
            }
            let _ = (first2d, last2d);
        }
        let _ = (PCONFUSION, &mut first);
    }
}

/// The candidate set guard (OCCT `continue` before the projection step when
/// no candidate survived the node filter).
fn set_candidate_guard(candidates: &super::ShapeSet) -> bool {
    !candidates.is_empty()
}
