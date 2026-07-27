use crate::bopds::ds::types::*;
use crate::bopds::face_info::FaceInfo;
use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::{Curve2dEval, Curve3, CurveEval, Line3, Plane, Surface3};
use rcad_kernel::topods;
use std::collections::HashMap;

/// Build DS from two topods::BRep shapes.
///
/// per-operand in dimension order: A_V, A_E, A_F, B_V, B_E, B_F.
pub fn new_from_topods(a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) -> DS {
    let tol = fuzzy_tol.max(TOLERANCE_ABS);
    let mut ds = DS {
        shapes: Vec::new(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
        vertex_shape_idx: Vec::new(),
        edge_shape_idx: Vec::new(),

        face_shape_idx: Vec::new(),
        interf_vv: Vec::new(),
        interf_ve: Vec::new(),
        interf_vf: Vec::new(),
        interf_ee: Vec::new(),
        interf_ef: Vec::new(),
        interf_ff: Vec::new(),
        interf_vz: Vec::new(),
        interf_ez: Vec::new(),
        interf_fz: Vec::new(),
        interf_zz: Vec::new(),
        intersection_curves: Vec::new(),
        shared_topology: SharedTopologyInfo::default(),
        shape_sd: ShapeSD::new(0, &SharedTopologyInfo::default()),
        common_blocks: Vec::new(),
        pave_blocks: Vec::new(),
        pave_blocks_pool: Vec::new(),
        face_info_pool: Vec::new(),
        locations: Vec::new(),
        increased_ss: std::collections::HashSet::new(),
        interf_tb: std::collections::HashSet::new(),
        map_ve: std::collections::HashMap::new(),
        shape_info: Vec::new(),
        nb_source_shapes: 0,
    };

    // per-operand loading A_V,A_E,A_F then B_V,B_E,B_F
    load_topods_brep(&mut ds, a, ShapeOrigin::ShapeA);
    load_topods_brep(&mut ds, b, ShapeOrigin::ShapeB);

    ds.compute_uv_boundaries();
    ds.build_face_reps();
    // ShapeInfo is built during init_shape_topo traversal (B2+B3 alignment).
    ds.nb_source_shapes = ds.shape_info.len();
    ds.build_map_ve();
    ds
}

/// Recursive DFS traversal of topods shape hierarchy.
/// Process Solid → Shell → Face → Wire → Edge → Vertex, pushing each shape
/// to ds.shapes[] and its type-specific flat array in depth-first order.
/// Skips already-visited shapes via visited set, maps topods ti → DS index.
#[allow(clippy::too_many_arguments)]
fn init_shape_topo(
    ds: &mut DS,
    brep: &topods::BRep,
    ti: usize,
    origin: ShapeOrigin,
    rank: usize,
    visited: &mut std::collections::HashSet<usize>,
    v_map: &mut HashMap<usize, usize>,
    e_map: &mut HashMap<usize, usize>,
    f_map: &mut HashMap<usize, usize>,
    w_map: &mut HashMap<usize, usize>,
    shell_map: &mut HashMap<usize, usize>,
    solid_map: &mut HashMap<usize, usize>,
    compsolid_map: &mut HashMap<usize, usize>,
    // Current shell counter for source_shell_idx assignment.
    shell_counter: &mut usize,
    // Current solid counter for source_solid_idx assignment.
    solid_counter: &mut usize,
) {
    if !visited.insert(ti) {
        return;
    }
    let ts = &brep.tshapes[ti];
    match &**ts {
        // OCCT BOPDS_DS.cxx L328-352: InitShape
        //   L342: myMapShapeIndex.Seek(aSubShape) — TopoDS identity check via
        //         map/hash, NOT spatial proximity. Two TopoDS_Vertex objects at
        //         the same position from different shapes are distinct and both
        //         are Appended (L344). rcad follows the same rule:
        //         visited.insert(ti) handles topological dedup (same OCCT shape
        //         visited twice); find_vertex_near (position-based dedup) has
        //         no OCCT equivalent and must NOT be called.
        //   L344: Append(aSubShape) — always creates new shape info entry.
        //         rcad: ds.push_vertex + v_map.insert + shape_info.push.
        topods::TShape::Vertex(vd) => {
            let vi = ds.vertex_count();
            ds.push_vertex(
                DSVertex { shape_idx: 0,
                    point: vd.point,
                    origin: Some(origin),
                    geom_tol: vd.tolerance,
                    is_internal: false,
                    location: 0,
                },
                Some(ts.clone()),
            );
            v_map.insert(ti, vi);
            // push_vertex keeps 1:1 between shapes[] and shape_info[]. However,
            // init_shape_topo processes sub-shapes recursively, and inner-wire
            // stub wires (push_wire with None tshape) create a TShape without a
            // ShapeInfo entry, breaking the 1:1 at this point.  The Wire/Face/
            // Shell/Solid handlers restore balance before returning.
            // push_vertex already created the ShapeInfo (keeps 1:1). Update it
            // with source-vertex data (not-new, proper rank, BRep box_gap).
            let last_si = ds.shape_info.len() - 1;
            ds.shape_info[last_si].has_brep = true;
            ds.shape_info[last_si].is_new = false;
            ds.shape_info[last_si].rank = rank;
            ds.shape_info[last_si].source_idx = vi;
            ds.shape_info[last_si].box_gap = vd.tolerance + TOLERANCE_ABS;
        }
        topods::TShape::Edge(ed) => {
            init_shape_topo(
                ds,
                brep,
                ed.first.index,
                origin,
                rank,
                visited,
                v_map,
                e_map,
                f_map,
                w_map,
                shell_map,
                solid_map,
                compsolid_map,
                shell_counter,
                solid_counter,
            );
            init_shape_topo(
                ds,
                brep,
                ed.last.index,
                origin,
                rank,
                visited,
                v_map,
                e_map,
                f_map,
                w_map,
                shell_map,
                solid_map,
                compsolid_map,
                shell_counter,
                solid_counter,
            );
            let start = v_map.get(&ed.first.index).copied().unwrap_or(0);
            let end = v_map.get(&ed.last.index).copied().unwrap_or(0);
            let curve = ed.curve.clone().unwrap_or_else(|| {
                let p0 = ds.vertex_point(start);
                let p1 = ds.vertex_point(end);
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: (p1 - p0).normalize_or_zero(),
                })
            });
            let t_range = ed.range;
            let mut vertex_params = HashMap::new();
            for (&vi, &param) in &ed.vertex_params {
                if let Some(&ds_vi) = v_map.get(&vi) {
                    vertex_params.insert(ds_vi, param);
                }
            }
            if !vertex_params.contains_key(&start) {
                vertex_params.insert(start, t_range[0]);
            }
            if !vertex_params.contains_key(&end) {
                vertex_params.insert(end, t_range[1]);
            }
            let is_geometric = ed.curve.is_some()
                && (matches!(
                    &curve,
                    Curve3::Line(_)
                        | Curve3::Circle(_)
                        | Curve3::Ellipse(_)
                        | Curve3::BSpline(_)
                        | Curve3::Bezier(_)
                ));
            let ds_ei = ds.edge_count();
            ds.push_edge(
                DSEdge { shape_idx: 0,
                    start_vertex: start,
                    end_vertex: end,
                    curve,
                    t_range,
                    origin,
                    geom_tol: ed.tolerance,
                    paves: Vec::new(),
                    pave_blocks: Vec::new(),
                    face_reps: Vec::new(),
                    is_internal: false,
                    vertex_params,
                    face_tolerances: Vec::new(),
                    is_geometric,
                    location: 0,
                },
                Some(ts.clone()),
            );
            e_map.insert(ti, ds_ei);
            // OCCT: PBs are lazily created via ChangePaveBlocks() on first access.
            // rcad: init is deferred to edge_pave_blocks_mut().

            // ShapeInfo is now created by push_edge() — no duplicate needed.
        }
        topods::TShape::Wire(wd) => {
            for esr in &wd.edges {
                init_shape_topo(
                    ds,
                    brep,
                    esr.index,
                    origin,
                    rank,
                    visited,
                    v_map,
                    e_map,
                    f_map,
                    w_map,
                    shell_map,
                    solid_map,
                    compsolid_map,
                    shell_counter,
                    solid_counter,
                );
            }
            let ds_edges: Vec<usize> = wd
                .edges
                .iter()
                .map(|esr| e_map.get(&esr.index).copied().unwrap_or(0))
                .collect();
            let sub: Vec<usize> = ds_edges
                .iter()
                .filter_map(|&ei| ds.edge_shape_idx.get(ei).copied())
                .collect();
            let wi = ds.push_wire(sub, Some(ts.clone()));
            w_map.insert(ti, wi);
        }
        topods::TShape::Face(fd) => {
            init_shape_topo(
                ds,
                brep,
                fd.outer_wire.index,
                origin,
                rank,
                visited,
                v_map,
                e_map,
                f_map,
                w_map,
                shell_map,
                solid_map,
                compsolid_map,
                shell_counter,
                solid_counter,
            );
            let outer_wire_data = brep.wire(fd.outer_wire);
            let boundary_edges_ordered = reorder_wire_topods(&outer_wire_data.edges, brep, e_map);
            let boundary_edges: Vec<usize> =
                boundary_edges_ordered.iter().map(|&(ei, _)| ei).collect();
            let boundary_edge_forwards: Vec<bool> =
                boundary_edges_ordered.iter().map(|&(_, fwd)| fwd).collect();
            let boundary_verts = {
                let edges = &outer_wire_data.edges;
                if edges.is_empty() {
                    Vec::new()
                } else if edges.len() == 1 {
                    let ed = brep.edge(edges[0]);
                    vec![
                        v_map.get(&ed.first.index).copied().unwrap_or(0),
                        v_map.get(&ed.last.index).copied().unwrap_or(0),
                    ]
                } else {
                    let mut verts = Vec::with_capacity(edges.len());
                    for i in 0..edges.len() {
                        let next_i = (i + 1) % edges.len();
                        let e = brep.edge(edges[i]);
                        let en = brep.edge(edges[next_i]);
                        let shared =
                            if e.first.index == en.first.index || e.first.index == en.last.index {
                                e.first.index
                            } else {
                                e.last.index
                            };
                        let non_shared = if shared == e.first.index {
                            e.last.index
                        } else {
                            e.first.index
                        };
                        verts.push(v_map.get(&non_shared).copied().unwrap_or(0));
                    }
                    verts
                }
            };
            // Update outer wire's ShapeInfo sub_shapes to traversal order.
            let ow_idx = w_map
                .get(&fd.outer_wire.index)
                .copied()
                .unwrap_or(usize::MAX);
            if ow_idx < ds.shape_info.len() {
                ds.shape_info[ow_idx].sub_shapes = boundary_edges
                    .iter()
                    .filter_map(|&ei| ds.edge_shape_idx.get(ei).copied())
                    .collect();
            }
            // Process inner wires
            for iw_sr in &fd.inner_wires {
                init_shape_topo(
                    ds,
                    brep,
                    iw_sr.index,
                    origin,
                    rank,
                    visited,
                    v_map,
                    e_map,
                    f_map,
                    w_map,
                    shell_map,
                    solid_map,
                    compsolid_map,
                    shell_counter,
                    solid_counter,
                );
            }
            let inner_boundary_edges: Vec<Vec<(usize, bool)>> = fd
                .inner_wires
                .iter()
                .map(|iw_sr| {
                    let iw_data = brep.wire(*iw_sr);
                    iw_data
                        .edges
                        .iter()
                        .map(|we_sr| {
                            (
                                e_map.get(&we_sr.index).copied().unwrap_or(0),
                                we_sr.orientation.is_forward(),
                            )
                        })
                        .collect()
                })
                .collect();
            let inner_wire_idxs: Vec<usize> = (0..fd.inner_wires.len())
                .map(|_| ds.push_wire(Vec::new(), None))
                .collect();
            for (ii, iw_sr) in fd.inner_wires.iter().enumerate() {
                let iw_data = brep.wire(*iw_sr);
                let iw_edges: Vec<usize> = iw_data
                    .edges
                    .iter()
                    .map(|we_sr| e_map.get(&we_sr.index).copied().unwrap_or(0))
                    .collect();
                let wi = inner_wire_idxs[ii];
                ds.shape_info[wi].sub_shapes = iw_edges
                    .iter()
                    .filter_map(|&ei| ds.edge_shape_idx.get(ei).copied())
                    .collect();
            }
            let surface = fd
                .surface
                .clone()
                .unwrap_or_else(|| Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)));
            let normal = match &surface {
                Surface3::Plane(p) => p.normal,
                _ => {
                    if boundary_verts.len() >= 3 {
                        let p0 = ds
                            .vertices
                            .get(boundary_verts[0])
                            .map(|v| v.point)
                            .unwrap_or(DVec3::ZERO);
                        let p1 = ds
                            .vertices
                            .get(boundary_verts[1])
                            .map(|v| v.point)
                            .unwrap_or(DVec3::ZERO);
                        let p2 = ds
                            .vertices
                            .get(boundary_verts[2])
                            .map(|v| v.point)
                            .unwrap_or(DVec3::ZERO);
                        (p1 - p0).cross(p2 - p0).normalize_or_zero()
                    } else {
                        DVec3::Z
                    }
                }
            };
            let ds_fi = ds.face_count();
            f_map.insert(ti, ds_fi);
            let bverts_copy = boundary_verts.clone();
            ds.push_face(
                DSFace { shape_idx: 0,
                    surface,
                    boundary_verts,
                    boundary_edges: boundary_edges.clone(),
                    boundary_edge_forwards,
                    inner_boundary_edges,
                    outer_wire_idx: Some(ow_idx),
                    inner_wire_idxs,
                    normal,
                    origin,
                    face_info: FaceInfo::default(),
                    source_face_idx: ds_fi,
                    geom_tol: fd.tolerance,
                    location: 0,
                    uv_boundary: None,
                    natural_restriction: fd.natural_restriction,
                    source_shell_idx: None,
                    source_solid_idx: None,
                    source_compsolid_idx: None,
                },
                Some(ts.clone()),
            );

            // ShapeInfo: sub_shapes = shapes[] indices of boundary edges
            let sub: Vec<usize> = boundary_edges
                .iter()
                .filter_map(|&ei| ds.edge_shape_idx.get(ei).copied())
                .collect();
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for &vi in &bverts_copy {
                if vi < ds.vertex_count() {
                    let p = ds.vertex_point(vi);
                    mn = mn.min(p);
                    mx = mx.max(p);
                }
            }
            // ShapeInfo must be 1:1 with shapes[].
            // debug_assert_eq is avoided: inner-wire stub wires create a TShape
            // (via push_wire) without a corresponding ShapeInfo entry.
            ds.shape_info.push(ShapeInfo {
                shape_type: rcad_kernel::topods::ShapeType::Face,
                sub_shapes: sub,
                flag: -1,
                reference: -1,
                has_brep: true,
                box_min: if mn.is_finite() { Some(mn) } else { None },
                box_max: if mx.is_finite() { Some(mx) } else { None },
                box_gap: fd.tolerance + TOLERANCE_ABS,
                is_new: false,
                rank,
                source_idx: ds_fi,
                is_internal: false,
                source_shell_idx: None,
                source_solid_idx: None,
                source_compsolid_idx: None,
            });
        }
        topods::TShape::Shell(shd) => {
            let prev_face_count = ds.face_count();
            for face_sr in &shd.faces {
                init_shape_topo(
                    ds,
                    brep,
                    face_sr.index,
                    origin,
                    rank,
                    visited,
                    v_map,
                    e_map,
                    f_map,
                    w_map,
                    shell_map,
                    solid_map,
                    compsolid_map,
                    shell_counter,
                    solid_counter,
                );
            }
            for fi in prev_face_count..ds.face_count() {
                ds.faces[fi].source_shell_idx = Some(*shell_counter);
            }
            let shell_face_idxs: Vec<usize> = (prev_face_count..ds.face_count()).collect();
            let shi = if shell_face_idxs.is_empty() {
                *shell_counter
            } else {
                let sub: Vec<usize> = shell_face_idxs
                    .iter()
                    .filter_map(|&fi| ds.face_shape_idx.get(fi).copied())
                    .collect();
                ds.push_shell(sub, Some(ts.clone()))
            };
            shell_map.insert(ti, shi);
            *shell_counter += 1;

        }
        topods::TShape::Solid(sd) => {
            let shell_start = *shell_counter;
            for shell_sr in &sd.shells {
                init_shape_topo(
                    ds,
                    brep,
                    shell_sr.index,
                    origin,
                    rank,
                    visited,
                    v_map,
                    e_map,
                    f_map,
                    w_map,
                    shell_map,
                    solid_map,
                    compsolid_map,
                    shell_counter,
                    solid_counter,
                );
            }
            for fi in 0..ds.face_count() {
                if ds.faces[fi]
                    .source_shell_idx
                    .map_or(false, |s| s >= shell_start)
                {
                    ds.faces[fi].source_solid_idx = Some(*solid_counter);
                }
            }
            let solid_shells: Vec<usize> = sd.shells
                .iter()
                .filter_map(|sr| shell_map.get(&sr.index).copied())
                .collect();
            let soi = if solid_shells.is_empty() {
                *solid_counter
            } else {
                ds.push_solid(solid_shells, Some(ts.clone()))
            };
            solid_map.insert(ti, soi);
            *solid_counter += 1;
        }
        topods::TShape::CompSolid(cs_solids) => {
            for sr in cs_solids {
                init_shape_topo(
                    ds,
                    brep,
                    sr.index,
                    origin,
                    rank,
                    visited,
                    v_map,
                    e_map,
                    f_map,
                    w_map,
                    shell_map,
                    solid_map,
                    compsolid_map,
                    shell_counter,
                    solid_counter,
                );
            }
            let ds_solid_indices: Vec<usize> = cs_solids
                .iter()
                .filter_map(|sr| solid_map.get(&sr.index).copied())
                .collect();
            let csi = if ds_solid_indices.is_empty() {
                *solid_counter
            } else {
                ds.push_compsolid(ds_solid_indices.clone(), Some(ts.clone()))
            };
            compsolid_map.insert(ti, csi);

            if !ds_solid_indices.is_empty() {
                // Update source_compsolid_idx
                for &si in &ds_solid_indices {
                    for fi in 0..ds.face_count() {
                        if ds.faces[fi].source_solid_idx == Some(si) {
                            ds.faces[fi].source_compsolid_idx = Some(csi);
                        }
                    }
                }
            }
        }
        topods::TShape::Compound(_) => {}
    }
}

/// Load one operand's shapes into DS using DFS hierarchy traversal.
/// Per-operand, calls init_shape_topo on top-level Solids/CompSolids, then
/// transfers pcurves as a post-processing step.
pub fn load_topods_brep(ds: &mut DS, brep: &topods::BRep, origin: ShapeOrigin) {
    let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut v_map: HashMap<usize, usize> = HashMap::new();
    let mut e_map: HashMap<usize, usize> = HashMap::new();
    let mut f_map: HashMap<usize, usize> = HashMap::new();
    let mut w_map: HashMap<usize, usize> = HashMap::new();
    let mut shell_map: HashMap<usize, usize> = HashMap::new();
    let mut solid_map: HashMap<usize, usize> = HashMap::new();
    let mut compsolid_map: HashMap<usize, usize> = HashMap::new();
    let mut shell_counter = 0usize;
    let mut solid_counter = 0usize;

    // Find top-level shapes: Solids and CompSolids
    // (these are the entry points for boolean operand shapes)
    for (ti, ts) in brep.tshapes.iter().enumerate() {
        match &**ts {
            topods::TShape::Solid(_) | topods::TShape::CompSolid(_) => {
                init_shape_topo(
                    ds,
                    brep,
                    ti,
                    origin,
                    if origin == ShapeOrigin::ShapeA { 0 } else { 1 },
                    &mut visited,
                    &mut v_map,
                    &mut e_map,
                    &mut f_map,
                    &mut w_map,
                    &mut shell_map,
                    &mut solid_map,
                    &mut compsolid_map,
                    &mut shell_counter,
                    &mut solid_counter,
                );
            }
            _ => {}
        }
    }

    // Post-process: transfer pcurves from topods edge data
    for (ti, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Edge(ed) = ts.as_ref() {
            let Some(&ds_ei) = e_map.get(&ti) else {
                continue;
            };
            for (&face_ti, (pcurve, param_start, param_end)) in &ed.pcurves {
                let Some(&ds_fi) = f_map.get(&face_ti) else {
                    continue;
                };
                if ds.edges[ds_ei]
                    .face_reps
                    .iter()
                    .any(|r| r.face_idx == ds_fi)
                {
                    continue;
                }
                let uv_start = pcurve.point_at(*param_start);
                let uv_end = pcurve.point_at(*param_end);
                let span = (uv_end - uv_start).length();
                if span < TOLERANCE_CLAMP_MIN || !span.is_finite() {
                    continue;
                }
                ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
                    face_idx: ds_fi,
                    pcurve: pcurve.clone(),
                    pcurve2: None,
                    pcurve_range: [*param_start, *param_end],
                    start_param: *param_start,
                    end_param: *param_end,
                });
            }
            for rep in &ed.representations {
                match rep {
                    topods::CurveRepresentation::CurveOnSurface {
                        face,
                        pcurve,
                        range,
                    } => {
                        let Some(&ds_fi) = f_map.get(face) else {
                            continue;
                        };
                        if ds.edges[ds_ei]
                            .face_reps
                            .iter()
                            .any(|r| r.face_idx == ds_fi)
                        {
                            continue;
                        }
                        let span = range[1] - range[0];
                        if span < TOLERANCE_CLAMP_MIN {
                            continue;
                        }
                        ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
                            face_idx: ds_fi,
                            pcurve: pcurve.clone(),
                            pcurve2: None,
                            pcurve_range: *range,
                            start_param: range[0],
                            end_param: range[1],
                        });
                    }
                    topods::CurveRepresentation::CurveOnClosedSurface {
                        face,
                        pcurve1,
                        pcurve2,
                        range,
                    } => {
                        let Some(&ds_fi) = f_map.get(face) else {
                            continue;
                        };
                        if ds.edges[ds_ei]
                            .face_reps
                            .iter()
                            .any(|r| r.face_idx == ds_fi)
                        {
                            continue;
                        }
                        let span = range[1] - range[0];
                        if span < TOLERANCE_CLAMP_MIN {
                            continue;
                        }
                        ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
                            face_idx: ds_fi,
                            pcurve: pcurve1.clone(),
                            pcurve2: Some(pcurve2.clone()),
                            pcurve_range: *range,
                            start_param: range[0],
                            end_param: range[1],
                        });
                    }
                    _ => {}
                }
            }
        }
    }
}

/// reorder wire edges by traversal order (TopExp_Explorer).
/// TopoDS version: reads edge vertex adjacency from topods::BRep.
/// Returns (DS edge index, forward_in_wire) pairs in traversal order.
pub fn reorder_wire_topods(
    wire_edges: &[topods::ShapeRef],
    brep: &topods::BRep,
    e_map: &HashMap<usize, usize>,
) -> Vec<(usize, bool)> {
    if wire_edges.len() <= 1 {
        return wire_edges
            .iter()
            .map(|we| {
                (
                    e_map.get(&we.index).copied().unwrap_or(0),
                    we.orientation.is_forward(),
                )
            })
            .collect();
    }
    let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, we_sr) in wire_edges.iter().enumerate() {
        let ed = brep.edge(*we_sr);
        adj.entry(ed.first.index).or_default().push(i);
        adj.entry(ed.last.index).or_default().push(i);
    }
    let first_ed = brep.edge(wire_edges[0]);
    let mut cur = first_ed.first.index;
    let mut used = vec![false; wire_edges.len()];
    let mut ordered = Vec::with_capacity(wire_edges.len());
    for _ in 0..wire_edges.len() {
        let next_i = adj
            .entry(cur)
            .or_default()
            .iter()
            .copied()
            .find(|&i| !used[i])
            .expect("wire is not closed -- broken topology");
        used[next_i] = true;
        let we_sr = &wire_edges[next_i];
        let ed = brep.edge(*we_sr);
        // orientation = Forward if traversal enters through first vertex,
        // Reversed if enters through last vertex (matches face wire winding).
        let is_fwd = ed.first.index == cur;
        ordered.push((e_map.get(&we_sr.index).copied().unwrap_or(0), is_fwd));
        cur = if is_fwd {
            ed.last.index
        } else {
            ed.first.index
        };
    }
    ordered
}

// ===== Tests: DS loading from topods::BRep + shapes array sync =====
#[cfg(test)]
mod tests {
    use super::*;
    use crate::bopds::ds::types::*;
    use crate::tolerance::TOLERANCE_ABS;
    use rcad_kernel::topods::{self, TShape};

    fn make_unit_box() -> topods::BRep {
        let mut br = topods::BRep::new();
        let v000 = br.add_tvertex(glam::DVec3::new(0.0, 0.0, 0.0));
        let v100 = br.add_tvertex(glam::DVec3::new(1.0, 0.0, 0.0));
        let v110 = br.add_tvertex(glam::DVec3::new(1.0, 1.0, 0.0));
        let v010 = br.add_tvertex(glam::DVec3::new(0.0, 1.0, 0.0));
        let v001 = br.add_tvertex(glam::DVec3::new(0.0, 0.0, 1.0));
        let v101 = br.add_tvertex(glam::DVec3::new(1.0, 0.0, 1.0));
        let v111 = br.add_tvertex(glam::DVec3::new(1.0, 1.0, 1.0));
        let v011 = br.add_tvertex(glam::DVec3::new(0.0, 1.0, 1.0));
        let e_bot = vec![
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::X,
                })),
                v000,
                v100,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Y,
                })),
                v100,
                v110,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: -glam::DVec3::X,
                })),
                v110,
                v010,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: -glam::DVec3::Y,
                })),
                v010,
                v000,
                [0.0, 1.0],
            ),
        ];
        let e_top = vec![
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::X,
                })),
                v001,
                v101,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Y,
                })),
                v101,
                v111,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: -glam::DVec3::X,
                })),
                v111,
                v011,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: -glam::DVec3::Y,
                })),
                v011,
                v001,
                [0.0, 1.0],
            ),
        ];
        let e_side = vec![
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Z,
                })),
                v000,
                v001,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Z,
                })),
                v100,
                v101,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Z,
                })),
                v110,
                v111,
                [0.0, 1.0],
            ),
            br.add_tedge(
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::Z,
                })),
                v010,
                v011,
                [0.0, 1.0],
            ),
        ];
        let _bottom_w = br.add_twire(e_bot.clone());
        let _top_w = br.add_twire(e_top.clone());
        let side_wires = [
            vec![e_bot[3], e_side[0], e_top[3], e_side[3]],
            vec![e_bot[0], e_side[1], e_top[0], e_side[0]],
            vec![e_bot[1], e_side[2], e_top[1], e_side[1]],
            vec![e_bot[2], e_side[3], e_top[2], e_side[2]],
        ];
        let plane = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
        ));
        let f_bottom = br.add_tface(
            Some(plane.clone()),
            _bottom_w,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        let f_top = br.add_tface(
            Some(plane.clone()),
            _top_w,
            vec![],
            None,
            None,
            vec![],
            false,
        );
        let mut all_faces = vec![f_bottom, f_top];
        for sw in &side_wires {
            let w = br.add_twire(sw.clone());
            let f = br.add_tface(Some(plane.clone()), w, vec![], None, None, vec![], false);
            all_faces.push(f);
        }
        let sh = br.add_tshell(all_faces);
        br.add_tsolid(vec![sh]);
        br
    }

    #[test]
    fn load_box_vertex_count() {
        let brep = make_unit_box();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        assert_eq!(ds.vertex_count(), 8, "box has 8 vertices");
    }

    #[test]
    fn load_box_edge_count() {
        let brep = make_unit_box();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        assert_eq!(ds.edge_count(), 12, "box has 12 edges");
    }

    #[test]
    fn load_box_face_count() {
        let brep = make_unit_box();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        assert_eq!(ds.face_count(), 6, "box has 6 faces");
    }

    #[test]
    fn load_box_shapes_sync() {
        let brep = make_unit_box();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        let expected_min = ds.vertex_count() + ds.edge_count() + ds.face_count();
        assert!(
            ds.shapes.len() >= expected_min,
            "shapes len {} >= verts+edges+faces {}",
            ds.shapes.len(),
            expected_min
        );
    }

    #[test]
    fn load_box_vertex_coords() {
        let brep = make_unit_box();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        assert!(
            ds.vertices
                .iter()
                .any(|v| v.point.distance_squared(glam::DVec3::ZERO) < 1e-10)
        );
        assert!(
            ds.vertices
                .iter()
                .any(|v| v.point.distance_squared(glam::DVec3::ONE) < 1e-10)
        );
    }

    #[test]
    fn push_vertex_maintains_sync() {
        let brep = make_unit_box();
        let mut ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        let nv_before = ds.vertex_count();
        let ns_before = ds.shapes.len();
        let vi = ds.push_vertex(
            DSVertex {
                shape_idx: 0,
                point: glam::DVec3::new(2.0, 2.0, 2.0),
                origin: None,
                geom_tol: TOLERANCE_ABS,
                is_internal: false,
                location: 0,
            },
            None,
        );
        assert_eq!(vi, nv_before, "push_vertex returns next vertex index");
        assert_eq!(ds.vertex_count(), nv_before + 1, "vertices array grew");
        assert_eq!(ds.shapes.len(), ns_before + 1, "shapes array grew");
        let new_sr_idx = ns_before;
        match &*ds.shapes[new_sr_idx] {
            TShape::Vertex(vd) => {
                assert!((vd.point - glam::DVec3::new(2.0, 2.0, 2.0)).length() < 1e-10);
            }
            _ => panic!("expected Vertex TShape"),
        }
    }

    #[test]
    fn push_edge_maintains_sync() {
        let brep = make_unit_box();
        let mut ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        let ne_before = ds.edge_count();
        let ns_before = ds.shapes.len();
        let ei = ds.push_edge(
            DSEdge {
                shape_idx: 0,
                start_vertex: 0,
                end_vertex: 1,
                curve: rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: glam::DVec3::ZERO,
                    direction: glam::DVec3::X,
                }),
                t_range: [0.0, 1.0],
                origin: ShapeOrigin::ShapeA,
                geom_tol: TOLERANCE_ABS,
                paves: Vec::new(),
                pave_blocks: Vec::new(),
                face_reps: Vec::new(),
                is_internal: false,
                vertex_params: HashMap::new(),
                face_tolerances: Vec::new(),
                is_geometric: true,
                location: 0,
            },
            None,
        );
        assert_eq!(ei, ne_before, "push_edge returns next edge index");
        assert_eq!(ds.edge_count(), ne_before + 1, "edges array grew");
        assert_eq!(ds.shapes.len(), ns_before + 1, "shapes array grew");
    }

    fn make_unit_sphere() -> topods::BRep {
        rcad_modeling::make_sphere_brep(glam::DVec3::ZERO, 1.0)
            .expect("Unit sphere creation failed")
    }

    #[test]
    fn load_sphere_faces() {
        let brep = make_unit_sphere();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        assert_eq!(ds.face_count(), 1, "sphere has 1 face");
    }

    #[test]
    fn load_sphere_surface_type() {
        let brep = make_unit_sphere();
        let ds = new_from_topods(&brep, &topods::BRep::new(), TOLERANCE_ABS);
        match ds.face_surface(0).unwrap() {
            rcad_kernel::geom::Surface3::Sphere(_) => {}
            other => panic!("expected Sphere surface, got {:?}", other),
        }
    }

    #[test]
    fn load_two_shapes_origin() {
        let sphere = make_unit_sphere();
        let box_brep = make_unit_box();
        let ds = new_from_topods(&sphere, &box_brep, TOLERANCE_ABS);
        let a_count = ds
            .faces
            .iter()
            .filter(|f| f.origin == ShapeOrigin::ShapeA)
            .count();
        let b_count = ds
            .faces
            .iter()
            .filter(|f| f.origin == ShapeOrigin::ShapeB)
            .count();
        assert_eq!(a_count, 1, "sphere (A) has 1 face");
        assert_eq!(b_count, 6, "box (B) has 6 faces");
    }

    #[test]
    fn load_two_shapes_vertex_edge_counts() {
        let sphere = make_unit_sphere();
        let box_brep = make_unit_box();
        let ds = new_from_topods(&sphere, &box_brep, TOLERANCE_ABS);
        assert!(
            ds.vertex_count() >= 8,
            "total vertices >= 8, got {}",
            ds.vertex_count()
        );
        assert!(
            ds.edge_count() >= 14,
            "total edges >= 14, got {}",
            ds.edge_count()
        );
        assert_eq!(ds.face_count(), 7, "total faces = 7 (6 box + 1 sphere)");
    }
}

// ───── OCCT-aligned per-type loading functions (prepareVertices/Edges/Faces/Solids) ─────
// Architecture diff: rcad loads all shape types in a single traversal (via load_topods_brep
// inside load_vertices_from_brep). The remaining prepare methods are structural no-ops since
// loading happens in the first pass. OCCT's BOPDS_DS Init uses four separate TopExp_Explorer
// passes (one per shape type), but rcad's flat topods::BRep storage allows a single scan.

/// OCCT BOPDS_DS::prepareVertices — load all shapes from a topods::BRep into DS.
/// Architecture diff: rcad loads all shape types in a single pass.
/// Captures per-operand A counts between the two operand loads.
pub fn load_vertices_from_brep(ds: &mut DS, brep: &topods::BRep, origin: ShapeOrigin) {
    // Guard: skip if shapes already loaded (re-entrant call from prepare_edges/faces/solids)
    if !ds.shape_info.is_empty()
        && ds
            .shape_info
            .iter()
            .any(|si| si.rank == (if origin == ShapeOrigin::ShapeA { 0 } else { 1 }))
    {
        return;
    }
    load_topods_brep(ds, brep, origin);
}

/// OCCT BOPDS_DS::prepareEdges — structural no-op (loading done in prepareVertices).
pub fn load_edges_from_brep(_ds: &mut DS, _brep: &topods::BRep, _origin: ShapeOrigin) {}

/// OCCT BOPDS_DS::prepareFaces — structural no-op (loading done in prepareVertices).
pub fn load_faces_from_brep(_ds: &mut DS, _brep: &topods::BRep, _origin: ShapeOrigin) {}

/// OCCT BOPDS_DS::prepareSolids — structural no-op (loading done in prepareVertices).
pub fn load_solids_from_brep(_ds: &mut DS, _brep: &topods::BRep, _origin: ShapeOrigin) {}
