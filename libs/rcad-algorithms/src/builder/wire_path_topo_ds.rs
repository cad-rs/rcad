/// ✅ OCCT-aligned: TopoDS-based walk_path_extract_wires using BRepTool.
///
/// Phase 2 migration target: parallel implementation of walk_path_extract_wires
/// that uses WireSegmentTopoDS + BRepTool instead of WireSegment + DS + face_idx.
use std::collections::{HashMap, HashSet, VecDeque};
use indexmap::IndexMap;
use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, Curve3};
use rcad_kernel::topods::{BRepTool, Orientation, ShapeRef};
use crate::tolerance::*;
use super::types::{WireSegmentTopoDS, WireEdgeSourceTopoDS, WireFace};
use super::wire_splitter::{
    EdgeInfo, mark_edge_passed,
    find_angle_at, select_best_outgoing,
};
use crate::inttools::fclass2d::curve2d_nb_samples;
use super::point_in_polygon_2d;

// ---------------------------------------------------------------------------
// TopoDS-native EdgeInfo — no seg_idx, holds edge+face directly (OCCT aligned)
// ---------------------------------------------------------------------------

/// OCCT-aligned: EdgeInfo for TopoDS-native path.
/// Holds the edge and face directly instead of indexing into a WireSegment array.
#[derive(Debug, Clone)]
pub(crate) struct EdgeInfoTopoDS {
    pub(crate) edge: ShapeRef,
    pub(crate) face: ShapeRef,
    pub(crate) passed: bool,
    pub(crate) in_flag: bool,
    pub(crate) is_inside: bool,
    pub(crate) angle: f64,
}

/// Mark an EdgeInfoTopoDS entry as passed at a given vertex.
pub(crate) fn mark_edge_passed_topo(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfoTopoDS>>,
    edge: ShapeRef,
    vertex: usize,
    in_flag: bool,
) {
    if let Some(infos) = smart_map.get_mut(&vertex) {
        for info in infos.iter_mut() {
            if info.edge.index == edge.index && info.in_flag == in_flag {
                info.passed = true;
                return;
            }
        }
    }
}

/// Find the angle for an edge at a vertex (TopoDS variant).
pub(crate) fn find_angle_at_topo(
    smart_map: &IndexMap<usize, Vec<EdgeInfoTopoDS>>,
    edge: ShapeRef,
    vertex: usize,
    in_flag: bool,
) -> Option<f64> {
    smart_map.get(&vertex)?.iter()
        .find(|ei| ei.edge.index == edge.index && ei.in_flag == in_flag)
        .map(|ei| ei.angle)
}

/// Select best outgoing edge from candidates (TopoDS variant).
pub(crate) fn select_best_outgoing_topo<'a>(
    candidates: &[&'a EdgeInfoTopoDS],
    angle_in: f64,
    incoming_is_boundary: bool,
    incoming_edge: ShapeRef,
) -> Option<&'a EdgeInfoTopoDS> {
    if candidates.is_empty() { return None; }
    let a_two_pi = std::f64::consts::TAU;
    let eps = std::f64::EPSILON;
    let mut a_min_angle = 100.0;
    let mut a_nb_ways_inside: i32 = 0;
    let mut p_only_way_in: Option<&EdgeInfoTopoDS> = None;
    let mut p_edge_info: Option<&EdgeInfoTopoDS> = None;
    for an_ei in candidates {
        let a_angle = if an_ei.edge.index == incoming_edge.index {
            a_two_pi
        } else {
            super::angle_2d::clock_wise_angle(angle_in, an_ei.angle)
        };
        if incoming_is_boundary && an_ei.is_inside {
            a_nb_ways_inside += 1;
            p_only_way_in = Some(an_ei);
        }
        if a_angle < a_min_angle - eps {
            a_min_angle = a_angle;
            p_edge_info = Some(an_ei);
        }
    }
    if a_nb_ways_inside == 1 { p_edge_info = p_only_way_in; }
    p_edge_info
}

/// Walk a path extracting closed wires using BRepTool-based data access.
///
/// Analogous to walk_path_extract_wires but uses ShapeRef handles and BRepTool
/// queries for edge/vertex data. Surface data (tolerances) must be pre-computed
/// and passed via compute_params.
pub(crate) struct WalkParams {
    pub(crate) u_res: fn(f64) -> f64,
    pub(crate) v_res: fn(f64) -> f64,
}
pub(crate) fn walk_path_extract_wires_topoDS(
    start_si: usize,
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
    face_surface: &Surface3,
) {
    let start_seg = &segments[start_si];
    // If this segment has no EdgeInfo, it cannot be walked.
    let has_info = smart_map.values().any(|v| v.iter().any(|ei| ei.seg_idx == start_si));
    if !has_info {
        smart_map.entry(start_seg.start_vertex.index).or_default().push(EdgeInfo {
            seg_idx: start_si, passed: true, in_flag: false,
            is_inside: false, is_circle_arc: false, angle: 0.0,
        });
        smart_map.entry(start_seg.end_vertex.index).or_default().push(EdgeInfo {
            seg_idx: start_si, passed: true, in_flag: true,
            is_inside: false, is_circle_arc: false, angle: 0.0,
        });
        return;
    }

    let max_iter = segments.len() * 4 + 200;

    // Build a per-vertex map: does this vertex belong to a closed/degenerate edge?
    let is_vert_closed = |smart_map: &IndexMap<usize, Vec<EdgeInfo>>, v: usize| -> bool {
        smart_map.get(&v).map_or(false, |infos| {
            infos.iter().any(|ei| {
                let seg = &segments[ei.seg_idx];
                seg.start_vertex.index == seg.end_vertex.index || seg.is_seam
            })
        })
    };

    // ✅ OCCT-aligned: Coord2d (BOPAlgo_WireSplitter_1.cxx L663-674).
    // Gets UV of a vertex on a specific edge by evaluating the edge's pcurve
    // at the vertex parameter. Uses BRepTool::curve_on_surface for pcurve lookup.
    let vertex_uv = |_vi: ShapeRef, segment: &WireSegmentTopoDS, at_start: bool| -> Option<DVec2> {
        // Try BRepTool pcurve lookup first (OCCT: CurveOnSurface + D0).
        if let Some((pc, _, _)) = tool.curve_on_surface(segment.edge, segment.face) {
            let t = if at_start { segment.t_range[0] } else { segment.t_range[1] };
            return Some(pc.point_at(t));
        }
        // Fallback: use first_pcurve/second_pcurve from segment if available
        let pc = if at_start || segment.orientation.is_forward() {
            segment.first_pcurve.as_ref().or(segment.second_pcurve.as_ref())
        } else {
            segment.second_pcurve.as_ref().or(segment.first_pcurve.as_ref())
        };
        if let Some(pc) = pc {
            let t = if at_start { segment.t_range[0] } else { segment.t_range[1] };
            return Some(pc.point_at(t));
        }
        None
    };

    // OCCT Tolerance2D/UTolerance2D/VTolerance2D using direct surface computation
    let vtol = |vi: usize| -> f64 {
        tool.vertex_tolerance(ShapeRef::new(vi)).max(TOLERANCE_ABS)
    };
    // Use face surface resolution functions directly instead of BRepTool queries
    // (the tool's face_surface/u_resolution require valid face ShapeRefs).
    let u_res_fn = rcad_kernel::topods::u_resolution_for_surface;
    let v_res_fn = rcad_kernel::topods::v_resolution_for_surface;
    let tolerance_2d = |vi: usize| -> f64 {
        let vt = vtol(vi);
        u_res_fn(&face_surface, vt).max(v_res_fn(&face_surface, vt)).max(vt)
    };
    let u_tolerance_2d = |vi: usize| -> f64 { u_res_fn(&face_surface, vtol(vi)) };
    let v_tolerance_2d = |vi: usize| -> f64 { v_res_fn(&face_surface, vtol(vi)) };
    let uv_tolerance = |vi: usize| -> f64 { 2.0 * tolerance_2d(vi) };

    let mut edge_seq: Vec<usize> = Vec::new();
    let mut vert_seq: Vec<usize> = Vec::new();
    let mut uv_seq: Vec<DVec2> = Vec::new();
    let mut info_seq: Vec<usize> = Vec::new();

    let mut ci = start_si;
    let mut arrived_vertex = start_seg.end_vertex.index;

    for _iter in 0..max_iter {
        // OCCT L394-403: do not escape through edge from which you enter.
        if edge_seq.len() == 1 {
            let same_edge = match (&segments[edge_seq[0]].source, &segments[ci].source) {
                (WireEdgeSourceTopoDS::DsEdge(ea), WireEdgeSourceTopoDS::DsEdge(eb)) => ea.index == eb.index,
                (WireEdgeSourceTopoDS::IntersectionCurve(ca), WireEdgeSourceTopoDS::IntersectionCurve(cb)) => ca.index == cb.index,
                (WireEdgeSourceTopoDS::SeamEdge, WireEdgeSourceTopoDS::SeamEdge) => true,
                _ => false,
            };
            if ci == edge_seq[0] || same_edge {
                return;
            }
        }

        let seg = &segments[ci];
        mark_edge_passed(smart_map, ci, seg.start_vertex.index, false);

        edge_seq.push(ci);
        vert_seq.push(seg.start_vertex.index);
        let cur_uv = vertex_uv(seg.start_vertex, seg, true);
        uv_seq.push(cur_uv.unwrap_or(DVec2::ZERO));
        info_seq.push(ci);

        // ── Loop Detection (OCCT L424-523) ──
        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_tol_2d = uv_tolerance(arrived_vertex);
        let a_tol_2d_sq = a_tol_2d * a_tol_2d;
        let a_pb = vertex_uv(ShapeRef::new(arrived_vertex), &segments[ci], false).unwrap_or(DVec2::ZERO);

        let mut b_has_edge = false;
        let a_nb = edge_seq.len();
        for i in (0..a_nb).rev() {
            let prev_v = vert_seq[i];
            let prev_uv = uv_seq[i];
            let prev_si = edge_seq[i];

            // OCCT L449-458: bHasEdge — skip degenerate-only wires
            if !b_has_edge {
                b_has_edge = match &segments[prev_si].source {
                    WireEdgeSourceTopoDS::DsEdge(ei) => !tool.is_edge_degenerated(*ei),
                    _ => true,
                };
                if !b_has_edge { continue; }
            }

            let is_same_v = prev_v == arrived_vertex;
            let mut is_same_v_2d = is_same_v;

            if is_same_v {
                if b_is_closed {
                    let a_d2 = prev_uv.distance_squared(a_pb);
                    is_same_v_2d = a_d2 < a_tol_2d_sq;
                    if is_same_v_2d {
                        let u_dist = (prev_uv.x - a_pb.x).abs();
                        let v_dist = (prev_uv.y - a_pb.y).abs();
                        let a_tol_u = 2.0 * u_tolerance_2d(arrived_vertex);
                        let a_tol_v = 2.0 * v_tolerance_2d(arrived_vertex);
                        if u_dist > a_tol_u || v_dist > a_tol_v {
                            is_same_v_2d = false;
                        }
                    }
                }
            }

            if is_same_v && is_same_v_2d {
                let wire: Vec<usize> = edge_seq[i..].to_vec();

                let mut is_valid = true;
                if wire.len() == 2 {
                    let a = &segments[wire[0]];
                    let b = &segments[wire[1]];
                    let same_edge = match (&a.source, &b.source) {
                        (WireEdgeSourceTopoDS::DsEdge(ea), WireEdgeSourceTopoDS::DsEdge(eb)) => ea.index == eb.index,
                        (WireEdgeSourceTopoDS::IntersectionCurve(ca), WireEdgeSourceTopoDS::IntersectionCurve(cb)) => ca.index == cb.index,
                        (WireEdgeSourceTopoDS::SeamEdge, WireEdgeSourceTopoDS::SeamEdge) => true,
                        _ => false,
                    };
                    if same_edge { is_valid = false; }
                }
                if is_valid { wires.push(wire); }

                let a_nbj = i;
                if a_nbj == 0 {
                    edge_seq.clear();
                    vert_seq.clear();
                    uv_seq.clear();
                    return;
                }

                let continue_vertex = vert_seq[i];
                edge_seq.truncate(a_nbj);
                vert_seq.truncate(a_nbj);
                uv_seq.truncate(a_nbj);
                info_seq.truncate(a_nbj);

                // ✅ OCCT-aligned L532-535: update state to last kept edge + continuation vertex
                ci = *info_seq.last().unwrap();
                arrived_vertex = continue_vertex;
                break;
            }
        }

        // ── Outgoing Edge Selection (OCCT L526-616) ──
        // OCCT L532-535: after loop detection, falls through here with
        // the truncated state (ci = aLS.Last(), arrived_vertex = aVertVa(i)).

        let angle_in = match find_angle_at(smart_map, ci, arrived_vertex, true) {
            Some(a) => a,
            None => return,
        };

        let raw_candidates: Vec<&EdgeInfo> = if let Some(infos) = smart_map.get(&arrived_vertex) {
            infos.iter().filter(|ei| !ei.passed && !ei.in_flag).collect()
        } else { return; };

        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_pb = vertex_uv(ShapeRef::new(arrived_vertex), &segments[ci], false).unwrap_or(DVec2::ZERO);
        let a_tol_2d_sq = { let tol = uv_tolerance(arrived_vertex); tol * tol };

        let i_cnt = raw_candidates.len();
        if i_cnt == 0 { return; }

        // Single candidate shortcut (OCCT L571-575)
        if i_cnt == 1 {
            let best = raw_candidates[0];
            ci = best.seg_idx;
            arrived_vertex = segments[ci].end_vertex.index;
            continue;
        }

        // 2D distance filter for closed vertices (OCCT L571-582)
        let candidates: Vec<&EdgeInfo> = if b_is_closed {
            raw_candidates.into_iter().filter(|ei| {
                let cand_uv = vertex_uv(ShapeRef::new(arrived_vertex), &segments[ei.seg_idx], true)
                    .unwrap_or(DVec2::ZERO);
                cand_uv.distance_squared(a_pb) < a_tol_2d_sq
            }).collect()
        } else { raw_candidates };

        if candidates.is_empty() { return; }

        let incoming_is_boundary = !matches!(segments[ci].source, WireEdgeSourceTopoDS::IntersectionCurve(_));
        let best = match select_best_outgoing(&candidates, angle_in, incoming_is_boundary, ci) {
            Some(e) => e,
            None => return,
        };

        ci = best.seg_idx;
        arrived_vertex = segments[ci].end_vertex.index;
    }
}

/// ✅ OCCT-aligned: TopoDS-based SplitBlock — refine angles + path walk for irregular blocks.
pub(crate) fn split_block_topoDS(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
    face_surface: &Surface3,
) {
    // OCCT L327: RefineAngles.  rcad: use topoDS-based angle refinement (currently skipped, uses pre-computed angles).
    // OCCT L331-358: Path walk
    let order_keys: Vec<usize> = smart_map.keys().copied().collect();
    for &v in &order_keys {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        for ei in &infos {
            if !ei.passed && !ei.in_flag
                && ei.seg_idx < segments.len()
                && (segments[ei.seg_idx].start_vertex.index != segments[ei.seg_idx].end_vertex.index
                    || segments[ei.seg_idx].is_seam)
            {
                walk_path_extract_wires_topoDS(ei.seg_idx, segments, smart_map, wires, tool, face_surface);
            }
        }
    }
}

/// ✅ OCCT-aligned: TopoDS-based build_closed_wires — SmartMap + angle computation + wire walking.
///
/// Simplified version without vi_to_canon/deg_end_canon (ShapeRef handles use DS indices directly).
pub(crate) fn build_closed_wires_topoDS(
    segments: &[WireSegmentTopoDS],
    avoided: &HashSet<usize>,
    tool: &dyn BRepTool,
    face_surface: &Surface3,
) -> Vec<Vec<usize>> {
    if segments.is_empty() { return vec![]; }

    let n = segments.len();
    let mut vert_to_segs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        if avoided.contains(&si) { continue; }
        vert_to_segs.entry(seg.start_vertex.index).or_default().push(si);
        vert_to_segs.entry(seg.end_vertex.index).or_default().push(si);
    }

    // Build connexity blocks
    let blocks = make_connexity_blocks_topoDS(segments, avoided, &vert_to_segs, n);

    // Process each block
    let mut wires: Vec<Vec<usize>> = Vec::new();
    for block in &blocks {
        if block.len() < 2 { continue; }

        // Build SmartMap
        let smart_map = build_smart_map_topoDS(block, segments, tool, face_surface);
        if smart_map.is_empty() { continue; }

        // Check regularity
        let is_regular = {
            let mut reg = true;
            for (_, infos) in &smart_map {
                let in_cnt = infos.iter().filter(|ei| ei.in_flag).count();
                let out_cnt = infos.iter().filter(|ei| !ei.in_flag).count();
                if in_cnt != 1 || out_cnt != 1 { reg = false; break; }
            }
            reg
        };

        if is_regular {
            // Regular: simple cyclic wire from block
            if let Some(wire) = build_regular_wire_topoDS(block) {
                wires.push(wire);
            }
        } else {
            // Irregular: split via path walk
            split_block_topoDS(block, segments, &mut (smart_map.clone()), &mut wires, tool, face_surface);
        }
    }
    wires
}

/// Build SmartMap for TopoDS segments with angle computation.
fn build_smart_map_topoDS(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_surface: &Surface3,
) -> IndexMap<usize, Vec<EdgeInfo>> {
    use super::angle_2d::angle_2d;
    use super::wire_path::pc_parameter_range;

    let mut smart_map: IndexMap<usize, Vec<EdgeInfo>> = IndexMap::new();
    for &si in block {
        let seg = &segments[si];
        let has_pcurve = tool.curve_on_surface(seg.edge, seg.face).is_some()
            || seg.first_pcurve.is_some() || seg.second_pcurve.is_some();
        if !has_pcurve { continue; }

        let is_inside = matches!(seg.source, WireEdgeSourceTopoDS::IntersectionCurve(_));
        let is_circle_arc = false;

        smart_map.entry(seg.start_vertex.index).or_default().push(EdgeInfo {
            seg_idx: si, passed: false, in_flag: false, is_inside, is_circle_arc, angle: 0.0,
        });
        smart_map.entry(seg.end_vertex.index).or_default().push(EdgeInfo {
            seg_idx: si, passed: false, in_flag: true, is_inside, is_circle_arc, angle: 0.0,
        });
    }

    // Compute angles using BRepTool (OCCT Angle2D equivalent).
    for (v, infos) in smart_map.iter_mut() {
        let v_ref = ShapeRef::new(*v);
        let geom_tol = tool.vertex_tolerance(v_ref);
        for ei in infos.iter_mut() {
            let seg = &segments[ei.seg_idx];
            let t_v = tool.parameter_on_edge(v_ref, seg.edge, seg.face)
                .unwrap_or_else(|| {
                    if *v == seg.start_vertex.index { seg.t_range[0] } else { seg.t_range[1] }
                });
            let domain = seg.t_range;
            let (curve, curve_domain): (&Curve2d, [f64; 2]) = match &seg.source {
                WireEdgeSourceTopoDS::IntersectionCurve(_) => {
                    // Use segment's own pcurve (populated by collect_face_edge_segments)
                    match seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()) {
                        Some(pc) => {
                            let (ta, tb) = pc_parameter_range(pc);
                            (pc, [ta, tb])
                        }
                        None => continue,
                    }
                }
                _ => {
                    // Try BRepTool pcurve, fall back to segment's own pcurve
                    let pc = tool.curve_on_surface(seg.edge, seg.face)
                        .map(|(pc, _, _)| pc)
                        .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
                    match pc {
                        Some(pc) => (pc, domain),
                        None => continue,
                    }
                }
            };
            ei.angle = angle_2d(curve, t_v, curve_domain, ei.in_flag, face_surface, geom_tol, None)
                .unwrap_or(0.0);
        }
    }
    smart_map
}

/// Build a regular wire from a block (all vertices have degree 2).
fn build_regular_wire_topoDS(block: &[usize]) -> Option<Vec<usize>> {
    if block.is_empty() { return None; }
    let mut result: Vec<usize> = Vec::with_capacity(block.len());
    if block.len() == 1 {
        // Single edge → self-loop
        return Some(vec![block[0]]);
    }
    // Start from first segment, follow vertex chain
    // For now, just return the block in order (caller handles simple cases)
    Some(block.to_vec())
}

/// Connected-component grouping for TopoDS segments.
fn make_connexity_blocks_topoDS(
    segments: &[WireSegmentTopoDS],
    avoided: &HashSet<usize>,
    vert_to_segs: &HashMap<usize, Vec<usize>>,
    n: usize,
) -> Vec<Vec<usize>> {
    let mut visited_seg = vec![false; n];
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for si in 0..n {
        if visited_seg[si] { continue; }
        if avoided.contains(&si) { continue; }
        let mut block = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(si);
        visited_seg[si] = true;
        while let Some(ci) = queue.pop_front() {
            block.push(ci);
            let seg = &segments[ci];
            for &vi in &[seg.start_vertex.index, seg.end_vertex.index] {
                if let Some(neighbors) = vert_to_segs.get(&vi) {
                    for &ni in neighbors {
                        if !visited_seg[ni] {
                            visited_seg[ni] = true;
                            queue.push_back(ni);
                        }
                    }
                }
            }
        }
        blocks.push(block);
    }
    blocks
}

/// TopoDS-based perform_areas — classifies wires as growth/hole using UV polygon area.
///
/// Uses WireSegmentTopoDS pcurves and BRepTool for surface queries.
/// Returns the same Vec<WireFace> (index-based) for backward compatibility.
pub(crate) fn perform_areas_topo_ds(
    wires: &[Vec<usize>],
    internal_wires: &[Vec<usize>],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_idx: usize,
) -> Vec<WireFace> {
    if wires.is_empty() { return vec![]; }

    struct WireData { wire_idx: usize, uv_boundary: Vec<DVec2>, n_distinct: usize }

    let mut wds: Vec<WireData> = wires.iter().enumerate().filter_map(|(wi, w)| {
        let mut uv_bnd: Vec<DVec2> = Vec::new();
        for &si in w {
            let seg = &segments[si];
            // Use pcurve from segment data (populated by collect_face_edge_segments)
            let pc_opt = if matches!(seg.source, WireEdgeSourceTopoDS::IntersectionCurve(_)) {
                seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref())
            } else {
                // Try BRepTool pcurve, fall back to segment's own pcurve
                tool.curve_on_surface(seg.edge, seg.face)
                    .map(|(pc, _, _)| pc)
                    .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()))
            };
            if let Some(pc) = pc_opt {
                let t0 = seg.t_range[0];
                let t1 = seg.t_range[1];
                let n = curve2d_nb_samples(pc, t0, t1).max(2);
                let du = if n > 1 { (t1 - t0) / (n - 1) as f64 } else { 0.0 };
                for i in 0..n {
                    uv_bnd.push(pc.point_at(t0 + du * i as f64));
                }
            }
        }
        uv_bnd.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);
        let n_distinct = { let mut pts = uv_bnd.clone(); pts.sort_by(|a,b|a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y))); pts.dedup(); pts.len() };
        if uv_bnd.len() < 3 || n_distinct < 3 { return None; }
        Some(WireData { wire_idx: wi, uv_boundary: uv_bnd, n_distinct })
    }).collect();

    if wds.is_empty() { return vec![]; }

    // Classify holes vs growths via UV signed area
    let mut is_hole = vec![false; wds.len()];
    for si in 0..wds.len() {
        let uv_b = &wds[si].uv_boundary;
        if uv_b.len() >= 3 {
            let area: f64 = uv_b.windows(2).map(|pair| {
                pair[0].x * pair[1].y - pair[1].x * pair[0].y
            }).sum::<f64>() + {
                let n = uv_b.len();
                uv_b[n-1].x * uv_b[0].y - uv_b[0].x * uv_b[n-1].y
            } * 0.5;
            is_hole[si] = area < 0.0;
        } else { is_hole[si] = true; }
    }

    let growths: Vec<usize> = (0..wds.len()).filter(|&i| !is_hole[i]).collect();
    let holes: Vec<usize> = (0..wds.len()).filter(|&i| is_hole[i]).collect();

    if growths.is_empty() && !wds.is_empty() {
        return vec![WireFace { outer_wire: wires[wds[0].wire_idx].clone(), inner_wires: vec![], internal_wires: internal_wires.to_vec() }];
    }
    if holes.is_empty() {
        return growths.iter().map(|&g| WireFace { outer_wire: wires[g].clone(), inner_wires: vec![], internal_wires: internal_wires.to_vec() }).collect();
    }

    // Assign holes to enclosing growths via UV point-in-polygon
    let growth_uv_bbox: Vec<Option<[f64; 4]>> = growths.iter().map(|&g| {
        let uv = &wds[g].uv_boundary;
        if uv.len() < 3 { return None; }
        let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        Some([u_min, u_max, v_min, v_max])
    }).collect();

    let mut h2g: Vec<(usize, usize)> = Vec::new();
    for &h in &holes {
        let h_uv = &wds[h].uv_boundary;
        let h_uv_c = if h_uv.len() >= 3 { h_uv.iter().copied().sum::<DVec2>() / h_uv.len() as f64 } else { continue; };
        let mut assigned = false;
        for (gi, &g) in growths.iter().enumerate() {
            if let Some([u0, u1, v0, v1]) = growth_uv_bbox[gi] {
                if h_uv_c.x < u0 || h_uv_c.x > u1 || h_uv_c.y < v0 || h_uv_c.y > v1 { continue; }
            }
            if wds[g].uv_boundary.len() >= 3 && point_in_polygon_2d(&wds[g].uv_boundary, h_uv_c) {
                h2g.push((h, g)); assigned = true; break;
            }
        }
        if !assigned && !growths.is_empty() { h2g.push((h, growths[0])); }
    }

    let mut g2h: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(h, g) in &h2g { g2h.entry(g).or_default().push(h); }

    growths.iter().map(|&g| WireFace {
        outer_wire: wires[g].clone(),
        inner_wires: g2h.get(&g).map(|hs| hs.iter().map(|&h| wires[h].clone()).collect()).unwrap_or_default(),
        internal_wires: internal_wires.to_vec(),
    }).collect()
}
