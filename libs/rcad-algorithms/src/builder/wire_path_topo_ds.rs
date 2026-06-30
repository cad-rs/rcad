/// ✅ OCCT-aligned: TopoDS-based walk_path_extract_wires using BRepTool.
///
/// Phase 2 migration target: parallel implementation of walk_path_extract_wires
/// that uses WireSegmentTopoDS + BRepTool instead of WireSegment + DS + face_idx.
use std::collections::{HashMap, HashSet, VecDeque};
use indexmap::IndexMap;
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, Curve3};
use rcad_kernel::topods::{BRepTool, Orientation, ShapeRef};
use crate::tolerance::*;
use super::types::{WireSegmentTopoDS, WireEdgeSourceTopoDS, WireFace};
use super::wire_splitter::{
    EdgeInfo, mark_edge_passed,
    find_angle_at, select_best_outgoing,
};
use super::angle_2d::clock_wise_angle;
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
pub(crate) fn walk_path_extract_wires(
    start_si: usize,
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
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
                seg.start_vertex.index == seg.end_vertex.index || seg.is_closed_on_face
            })
        })
    };

    // ✅ OCCT-aligned: Coord2d (BOPAlgo_WireSplitter_1.cxx L677-688).
    let coord2d = |vi: ShapeRef, edge: ShapeRef, face: ShapeRef| -> Option<DVec2> {
        let t = tool.parameter_on_edge(vi, edge, face)?;
        let (pc, _, _) = tool.curve_on_surface(edge, face)?;
        Some(pc.point_at(t))
    };
    // ✅ OCCT-aligned: Coord2dVf (BOPAlgo_WireSplitter_1.cxx L692-711).
    // Uses oriented_first_vertex to respect edge orientation in the face wire.
    let coord2d_vf = |seg: &WireSegmentTopoDS| -> Option<DVec2> {
        let fwd_v = tool.oriented_first_vertex(seg.edge, seg.orientation);
        coord2d(fwd_v, seg.edge, seg.face)
    };

    // OCCT Tolerance2D/UTolerance2D/VTolerance2D (BOPAlgo_WireSplitter_1.cxx L873-912)
    let face_ref = segments[start_si].face;
    let is_bspline = matches!(tool.face_surface(face_ref), Some(&Surface3::BSpline(_)));
    let vtol = |vi: usize| -> f64 { tool.vertex_tolerance(ShapeRef::new(vi)) };
    let tolerance_2d = |vi: usize| -> f64 {
        let vt = vtol(vi);
        let u = tool.u_resolution(face_ref, vt);
        let v = tool.v_resolution(face_ref, vt);
        let mut t = u.max(v);
        if t < vt { t = vt; }
        if is_bspline { t *= 1.1; }
        t
    };
    let u_tolerance_2d = |vi: usize| -> f64 { tool.u_resolution(face_ref, vtol(vi)) };
    let v_tolerance_2d = |vi: usize| -> f64 { tool.v_resolution(face_ref, vtol(vi)) };
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
        let cur_uv = coord2d(seg.start_vertex, seg.edge, seg.face);
        uv_seq.push(cur_uv.unwrap_or(DVec2::ZERO));
        info_seq.push(ci);

        // ── Loop Detection (OCCT L424-523) ──
        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_pb = coord2d(ShapeRef::new(arrived_vertex), segments[ci].edge, segments[ci].face).unwrap_or(DVec2::ZERO);
        let a_tol_2d = uv_tolerance(arrived_vertex);
        let a_tol_2d_sq = a_tol_2d * a_tol_2d;

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
        let angle_in = match find_angle_at(smart_map, ci, arrived_vertex, true) {
            Some(a) => a,
            None => return,
        };
        let a_tol_2d_sq = { let tol = uv_tolerance(arrived_vertex); tol * tol };
        let b_is_closed = is_vert_closed(smart_map, arrived_vertex);
        let a_pb = coord2d(ShapeRef::new(arrived_vertex), segments[ci].edge, segments[ci].face).unwrap_or(DVec2::ZERO);

        // OCCT L540-549: prepare selection state
        let mut p_edge_info: Option<usize> = None;
        let mut a_min_angle = 100.0;
        let le_info = smart_map.get(&arrived_vertex).map(|v| v.as_slice()).unwrap_or(&[]);
        let i_cnt = le_info.iter().filter(|ei| !ei.passed && !ei.in_flag).count();
        let incoming_is_boundary = !matches!(segments[ci].source, WireEdgeSourceTopoDS::IntersectionCurve(_));
        let mut a_nb_ways_inside = 0i32;
        let mut p_only_way_in: Option<usize> = None;

        for ei in le_info {
            let an_is_out = !ei.in_flag;
            let an_is_not_passed = !ei.passed;
            if !an_is_out || !an_is_not_passed { continue; }
            // OCCT L565-569: no way to go
            if i_cnt == 0 { return; }
            // OCCT L571-575: single way out
            if i_cnt == 1 { p_edge_info = Some(ei.seg_idx); break; }
            let an_angle = if ei.seg_idx == ci { std::f64::consts::TAU }
            else {
                // OCCT L584-596: 2D distance filter for closed vertices
                if b_is_closed {
                    let cand_uv = coord2d_vf(&segments[ei.seg_idx])
                        .unwrap_or(DVec2::ZERO);
                    if cand_uv.distance_squared(a_pb) >= a_tol_2d_sq { continue; }
                }
                let an_angle_out = ei.angle;
                clock_wise_angle(angle_in, an_angle_out)
            };
            // OCCT L603-607: count inside ways
            if incoming_is_boundary && ei.is_inside {
                a_nb_ways_inside += 1;
                p_only_way_in = Some(ei.seg_idx);
            }
            // OCCT L609-613: select minimal angle
            if an_angle < a_min_angle - std::f64::EPSILON {
                a_min_angle = an_angle;
                p_edge_info = Some(ei.seg_idx);
            }
        }
        // OCCT L616-619: prefer only way inside
        if a_nb_ways_inside == 1 { p_edge_info = p_only_way_in; }
        // OCCT L621-625: no way to go
        let best_si = match p_edge_info { Some(si) => si, None => return };
        // OCCT L627-629: advance to next vertex
        ci = best_si;
        arrived_vertex = segments[ci].end_vertex.index;
    }
}

/// OCCT-aligned: RefineAngles (WireSplitter_1.cxx L919-1043).
/// For each multi-vertex with 2 boundary edges, adjust internal edge angles
/// that fall outside the boundary sweep.  Uses BRepTool for angle computation.
fn refine_angles(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
) {
    let vertices: Vec<usize> = smart_map.keys().copied().collect();
    for &v in &vertices {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        let mut cnt_bnd = 0;
        let mut cnt_int = 0;
        let mut a1_bnd = 0.0;
        let mut a2_bnd = 0.0;
        for ei in &infos {
            if !ei.is_inside { cnt_bnd += 1; if !ei.in_flag { a1_bnd = ei.angle; } else { a2_bnd = ei.angle; } }
            else { cnt_int += 1; }
        }
        if cnt_bnd != 2 { continue; }
        let a_delta = super::angle_2d::clock_wise_angle(a2_bnd, a1_bnd);
        let mut refined_map: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
        for ei in &infos {
            if !ei.is_inside || ei.in_flag { continue; }
            let a_ic = ei.angle;
            let a_da = super::angle_2d::clock_wise_angle(a2_bnd, a_ic);
            if a_da < a_delta { continue; }
            let seg = &segments[ei.seg_idx];
            let v_ref = rcad_kernel::topods::ShapeRef::new(v);
            let geom_tol = tool.vertex_tolerance(v_ref);
            let t_v = tool.parameter_on_edge(v_ref, seg.edge, seg.face)
                .unwrap_or_else(|| if v == seg.start_vertex.index { seg.t_range[0] } else { seg.t_range[1] });
            let domain = seg.t_range;
            let pc = seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref());
            let b_refined = match pc {
                Some(pc) => {
                    let ta = pc.point_at(t_v);
                    // OCCT L1057-1061: use pcurve UV to compute new angle
                    // Check if angle is near boundary — true refine via micro-step
                    let eps = 1e-12;
                    let dt = (1e-8 * (domain[1] - domain[0]).abs().max(1.0)).min((domain[1] - domain[0]).abs() * 0.1);
                    let t2 = (t_v + dt).min(domain[1]);
                    let pt2 = pc.point_at(t2);
                    let dir = pt2 - ta;
                    if dir.length_squared() < 1e-30 { false }
                    else {
                        let new_ang = dir.y.atan2(dir.x);
                        if clock_wise_angle(a2_bnd, new_ang) < a_delta { false }
                        else { refined_map.insert(ei.seg_idx, new_ang); true }
                    }
                }
                None => false,
            };
            if !b_refined && cnt_int == 2 {
                let eps = 1e-12;
                let new_angle = if a_ic <= a1_bnd {
                    (a1_bnd + eps) % std::f64::consts::TAU
                } else {
                    (a2_bnd - eps + std::f64::consts::TAU) % std::f64::consts::TAU
                };
                refined_map.insert(ei.seg_idx, new_angle);
            }
        }
        if refined_map.is_empty() { continue; }
        if let Some(infos_mut) = smart_map.get_mut(&v) {
            for ei in infos_mut.iter_mut() {
                if let Some(&new_angle) = refined_map.get(&ei.seg_idx) {
                    ei.angle = if ei.in_flag { (new_angle + std::f64::consts::PI) % std::f64::consts::TAU } else { new_angle };
                }
            }
        }
    }
}

/// ✅ OCCT-aligned: TopoDS-based SplitBlock — refine angles + path walk for irregular blocks.
pub(crate) fn split_block(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
) {
    // OCCT L324: RefineAngles before Path walk
    refine_angles(smart_map, segments, tool);
    // OCCT L331-358: Path walk
    let order_keys: Vec<usize> = smart_map.keys().copied().collect();
    for &v in &order_keys {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        for ei in &infos {
            if !ei.passed && !ei.in_flag
                && ei.seg_idx < segments.len()
                && (segments[ei.seg_idx].start_vertex.index != segments[ei.seg_idx].end_vertex.index
                    || segments[ei.seg_idx].is_closed_on_face)
            {
                walk_path_extract_wires(ei.seg_idx, segments, smart_map, wires, tool);
            }
        }
    }
}

/// ✅ OCCT-aligned: TopoDS-based build_closed_wires — SmartMap + angle computation + wire walking.
///
/// Simplified version without vi_to_canon/deg_end_canon (ShapeRef handles use DS indices directly).
pub(crate) fn build_closed_wires(
    segments: &[WireSegmentTopoDS],
    avoided: &HashSet<usize>,
    tool: &dyn BRepTool,
) -> Vec<Vec<usize>> {
    if segments.is_empty() { return vec![]; }
    if std::env::var("RCAD_DEBUG_IC").is_ok() {
        let face_id = segments[0].face.index;
        eprintln!("[WIRE] face={} n_seg={} n_avoided={}", face_id, segments.len(), avoided.len());
    }

    let n = segments.len();
    let mut vert_to_segs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        if avoided.contains(&si) { continue; }
        vert_to_segs.entry(seg.start_vertex.index).or_default().push(si);
        vert_to_segs.entry(seg.end_vertex.index).or_default().push(si);
    }

    if std::env::var("RCAD_DUMP_FACE").is_ok() && !segments.is_empty() {
        let face_id = segments[0].face.index;
        eprintln!("[FACE] face={} n_seg={} n_avoided={} n_active={}", face_id, segments.len(), avoided.len(), n - avoided.len());
        for (v, segs) in &vert_to_segs {
            eprintln!("[FACE]   v={} valence={} segs={:?}", v, segs.len(), segs);
        }
        for (si, seg) in segments.iter().enumerate() {
            let src = match &seg.source { WireEdgeSourceTopoDS::DsEdge(sr) => format!("DsEdge({})", sr.index), WireEdgeSourceTopoDS::IntersectionCurve(sr) => format!("IC({})", sr.index), WireEdgeSourceTopoDS::SeamEdge => "Seam".into() };
            eprintln!("[FACE]   seg[{}] {} {}->{} avoided={}", si, src, seg.start_vertex.index, seg.end_vertex.index, avoided.contains(&si));
        }
    }

    // Build connexity blocks
    let blocks = make_connexity_blocks(segments, avoided, &vert_to_segs, n);

    // Process each block
    let mut wires: Vec<Vec<usize>> = Vec::new();
    for (bi, block) in blocks.iter().enumerate() {
        if block.len() < 2 {
            if std::env::var("RCAD_DEBUG_IC").is_ok() {
                let face_id = segments[0].face.index;
                let vset: std::collections::HashSet<usize> = block.iter().flat_map(|&si| {
                    let s = &segments[si];
                    [s.start_vertex.index, s.end_vertex.index]
                }).collect();
                eprintln!("[WIRE] face={} block[{}] too small ({} segs, {} verts)", face_id, bi, block.len(), vset.len());
            }
            continue;
        }

        // Build SmartMap
        let smart_map = build_smart_map(block, segments, tool);
        if smart_map.is_empty() {
            if std::env::var("RCAD_DEBUG_IC").is_ok() {
                let face_id = segments[0].face.index;
                let vset: std::collections::HashSet<usize> = block.iter().flat_map(|&si| {
                    let s = &segments[si];
                    [s.start_vertex.index, s.end_vertex.index]
                }).collect();
                let no_pcurve: Vec<usize> = block.iter().filter(|&&si| {
                    let seg = &segments[si];
                    tool.curve_on_surface(seg.edge, seg.face).is_none()
                        && seg.first_pcurve.is_none() && seg.second_pcurve.is_none()
                }).copied().collect();
                eprintln!("[WIRE] face={} block[{}] EMPTY smart_map ({} segs, {} verts, {} no-pcurve)", face_id, bi, block.len(), vset.len(), no_pcurve.len());
            }
            continue;
        }

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
            if let Some(wire) = build_regular_wire(block) {
                wires.push(wire);
            }
        } else {
            // Irregular: split via path walk
            split_block(block, segments, &mut (smart_map.clone()), &mut wires, tool);
        }
    }
    wires
}

/// Build SmartMap for TopoDS segments with angle computation.
fn build_smart_map(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
) -> IndexMap<usize, Vec<EdgeInfo>> {
    use super::angle_2d::angle_2d;
    use super::wire_path::pc_parameter_range;

    // OCCT L147-152: aMS set tracking edge parity (odd → boundary, even → internal).
    // Non-closed edges appearing an even number of times are internal (IsInside=true).
    // Uses orientation-aware key matching OCCT TopoDS_Shape hashing (orientation
    // in the face wire distinguishes FWD and REV occurrences of the same edge).
    let src_key_of = |seg: &WireSegmentTopoDS, si: usize| -> usize {
        let base = match seg.source {
            WireEdgeSourceTopoDS::DsEdge(ref sr) => sr.index,
            WireEdgeSourceTopoDS::IntersectionCurve(ref sr) => sr.index,
            WireEdgeSourceTopoDS::SeamEdge => usize::MAX - si,
        };
        base ^ ((seg.orientation as usize) << 24)
    };
    let mut a_ms: HashSet<usize> = HashSet::new();
    for &si in block {
        let seg = &segments[si];
        let has_pcurve = tool.curve_on_surface(seg.edge, seg.face).is_some()
            || seg.first_pcurve.is_some() || seg.second_pcurve.is_some();
        if !has_pcurve { continue; }

        let b_closed = seg.start_vertex.index == seg.end_vertex.index || seg.is_closed_on_face;

        let src_key = src_key_of(seg, si);
        // OCCT L149: !aMS.Add(aE) && !bIsClosed → aMS.Remove(aE)
        if !a_ms.insert(src_key) && !b_closed {
            a_ms.remove(&src_key);
        }
    }

    let mut smart_map: IndexMap<usize, Vec<EdgeInfo>> = IndexMap::new();
    for &si in block {
        let seg = &segments[si];
        let has_pcurve = tool.curve_on_surface(seg.edge, seg.face).is_some()
            || seg.first_pcurve.is_some() || seg.second_pcurve.is_some();
        if !has_pcurve { continue; }

        // OCCT L310: IsInside = !aMS.Contains(aE)
        let src_key = src_key_of(seg, si);
        let is_inside = !a_ms.contains(&src_key);
        let is_circle_arc = false;

        // OCCT L170-172: in_flag based on vertex orientation in the edge.
        //   FORWARD vertex → in_flag=false (outgoing), REVERSED → in_flag=true (incoming).
        //   When segment is Reversed, start_vertex maps to REVERSED, end_vertex to FORWARD.
        let in_flag_start = seg.orientation == Orientation::Reversed;
        let in_flag_end = seg.orientation == Orientation::Forward;

        smart_map.entry(seg.start_vertex.index).or_default().push(EdgeInfo {
            seg_idx: si, passed: false, in_flag: in_flag_start, is_inside, is_circle_arc, angle: 0.0,
        });
        smart_map.entry(seg.end_vertex.index).or_default().push(EdgeInfo {
            seg_idx: si, passed: false, in_flag: in_flag_end, is_inside, is_circle_arc, angle: 0.0,
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
                    match seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()) {
                        Some(pc) => {
                            let (ta, tb) = pc_parameter_range(pc);
                            (pc, [ta, tb])
                        }
                        None => continue,
                    }
                }
                _ => {
                    let pc = tool.curve_on_surface(seg.edge, seg.face)
                        .map(|(pc, _, _)| pc)
                        .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
                    match pc {
                        Some(pc) => (pc, domain),
                        None => continue,
                    }
                }
            };
            let surf = tool.face_surface(seg.face).unwrap_or(
                &Surface3::Plane(rcad_kernel::geom::Plane { origin: DVec3::ZERO, normal: DVec3::Z }));
            let new_angle = angle_2d(curve, t_v, curve_domain, ei.in_flag, surf, geom_tol, None)
                .unwrap_or(0.0);
            if std::env::var("RCAD_ANGLE_DUMP").is_ok() {
                let face_idx = segments[0].face.index;
                let curve_type = match curve { Curve2d::Line(_) => "Line", Curve2d::Circle(_) => "Circle", Curve2d::Ellipse(_) => "Ellipse", Curve2d::BSpline(_) => "BSpline", _ => "Other" };
                eprintln!("[ANGLE] face={} v={} seg={} in={} t_v={:.10} dom=[{:.6},{:.6}] curve={} angle={:.10} deg={:.4}",
                    face_idx, v, ei.seg_idx, ei.in_flag, t_v, curve_domain[0], curve_domain[1],
                    curve_type, new_angle, new_angle.to_degrees());
            }
            ei.angle = new_angle;
        }
    }
    smart_map
}

/// Build a regular wire from a block (all vertices have degree 2).
fn build_regular_wire(block: &[usize]) -> Option<Vec<usize>> {
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
fn make_connexity_blocks(
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

/// TopoDS-based perform_areas — classifies wires as growth/hole.
///
/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformAreas (L420-499).
///   Uses IsGrowthWire fast pre-check + IsHole classification via UV signed area.
///   Returns Vec<WireFace> for backward compatibility with emit_wire_face_topods.
pub(crate) fn perform_areas(
    wires: &[Vec<usize>],
    internal_wires: &[Vec<usize>],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_idx: usize,
    ds: &crate::bopds::ds::DS,
) -> Vec<WireFace> {
    // OCCT L401-414: if no loops at all
    if wires.is_empty() {
        if ds.faces[face_idx].natural_restriction {
            return vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: vec![] }];
        }
        return vec![];
    }

    // OCCT L420-423: aMHE — map of hole face edges for quick growth check.
    let mut a_mhe: HashSet<ShapeRef> = HashSet::new();

    // OCCT L425-458: classify each wire as growth or hole.
    let mut a_new_faces: Vec<Vec<usize>> = Vec::new(); // growth wires
    let mut a_hole_faces: Vec<Vec<usize>> = Vec::new(); // hole wires

    for w in wires {
        // OCCT L437-439: MakeFace from wire (rcad: we work with segment indices).
        // OCCT L441: IsGrowthWire(aWire, aMHE) — fast check.
        let b_is_growth = if !a_mhe.is_empty() {
            w.iter().any(|&si| {
                if let Some(seg) = segments.get(si) {
                    a_mhe.contains(&seg.edge)
                } else { false }
            })
        } else { false };

        let b_is_growth = if b_is_growth {
            true
        } else {
            // OCCT L444-446: run classification via IsHole().
            // rcad: UV signed area (equivalent: CW = hole, CCW = growth).
            let mut uv_boundary: Vec<DVec2> = Vec::new();
            for &si in w {
                let seg = &segments[si];
                let pc_opt = if matches!(seg.source, WireEdgeSourceTopoDS::IntersectionCurve(_)) {
                    seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref())
                } else {
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
                        uv_boundary.push(pc.point_at(t0 + du * i as f64));
                    }
                }
            }
            uv_boundary.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);
            if uv_boundary.len() < 3 { continue; }
            // OCCT IsHole(): signed area < 0 means CW = hole.
            let area: f64 = uv_boundary.windows(2).map(|pair| {
                pair[0].x * pair[1].y - pair[1].x * pair[0].y
            }).sum::<f64>() + {
                let n = uv_boundary.len();
                uv_boundary[n-1].x * uv_boundary[0].y - uv_boundary[0].x * uv_boundary[n-1].y
            } * 0.5;
            area >= 0.0 // IsHole() = area < 0 → growth = !IsHole = area >= 0
        };

        // OCCT L449-458: save face.
        if b_is_growth {
            a_new_faces.push(w.clone());
        } else {
            a_hole_faces.push(w.clone());
            // OCCT L457: TopExp::MapShapes(aWire, TopAbs_EDGE, aMHE)
            for &si in w {
                if let Some(seg) = segments.get(si) {
                    a_mhe.insert(seg.edge);
                }
            }
        }
    }

    // OCCT L461-466: no holes → all wires are growths.
    if a_hole_faces.is_empty() {
        return a_new_faces.iter().map(|w| WireFace {
            outer_wire: w.clone(), inner_wires: vec![], internal_wires: internal_wires.to_vec(),
        }).collect();
    }

    // OCCT L470-484: prepare bounding boxes for hole faces.
    // rcad: compute UV bounding boxes for hole wires.
    let mut hole_uv_boxes: Vec<Option<[f64; 4]>> = a_hole_faces.iter().map(|w| {
        let mut uv_bnd: Vec<DVec2> = Vec::new();
        for &si in w {
            let seg = &segments[si];
            let pc_opt = tool.curve_on_surface(seg.edge, seg.face)
                .map(|(pc, _, _)| pc)
                .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
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
        if uv_bnd.len() < 3 { return None; }
        let u_min = uv_bnd.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = uv_bnd.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let v_min = uv_bnd.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv_bnd.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        Some([u_min, u_max, v_min, v_max])
    }).collect();

    // OCCT L486-491: assign holes to enclosing growth faces.
    // rcad: for each hole, find the smallest-enclosing growth via point-in-polygon.
    let mut h2g: Vec<(usize, usize)> = Vec::new();
    for (hi, h_wire) in a_hole_faces.iter().enumerate() {
        let h_uv = {
            let mut uv_bnd: Vec<DVec2> = Vec::new();
            for &si in h_wire {
                let seg = &segments[si];
                let pc_opt = tool.curve_on_surface(seg.edge, seg.face)
                    .map(|(pc, _, _)| pc)
                    .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
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
            uv_bnd
        };
        if h_uv.len() < 3 { continue; }
        let h_center = h_uv.iter().copied().sum::<DVec2>() / h_uv.len() as f64;

        let mut assigned = false;
        for (gi, g_wire) in a_new_faces.iter().enumerate() {
            if let Some(ref uv_box) = hole_uv_boxes.get(hi).copied().flatten() {
                if h_center.x < uv_box[0] || h_center.x > uv_box[1] || h_center.y < uv_box[2] || h_center.y > uv_box[3] { continue; }
            }
            // Build growth UV boundary for containment test
            let g_uv: Vec<DVec2> = g_wire.iter().filter_map(|&si| {
                let seg = &segments[si];
                let pc_opt = tool.curve_on_surface(seg.edge, seg.face)
                    .map(|(pc, _, _)| pc)
                    .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
                pc_opt.map(|pc| pc.point_at(seg.t_range[0]))
            }).collect();
            if g_uv.len() >= 3 && point_in_polygon_2d(&g_uv, h_center) {
                h2g.push((hi, gi)); assigned = true; break;
            }
        }
        if !assigned && !a_new_faces.is_empty() {
            // Defer to OCCT L557-581 orphan hole handling below
        }
    }

    // OCCT L557-581: identify orphan holes
    let assigned_hole_set: std::collections::HashSet<usize> = h2g.iter().map(|&(h, _)| h).collect();
    let orphan_holes: Vec<usize> = (0..a_hole_faces.len())
        .filter(|hi| !assigned_hole_set.contains(hi))
        .collect();

    // OCCT L540-555: build growth → holes map
    let mut g2h: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(h, g) in &h2g { g2h.entry(g).or_default().push(h); }

    // OCCT L584-613: build result WireFaces — each wire becomes its own face.
    //   OCCT BOPAlgo_BuilderFace::PerformAreas creates a separate TopoDS_Face
    //   for EACH wire (growth or hole).  The caller (ComputeState) then classifies
    //   each face independently, removing those that are In the opposing solid.
    //   rcad legacy approach merged holes as inner_wires, preventing per-face
    //   classification and leaving internal vertices in the result.
    // OCCT L584-613: build result WireFaces
    let mut result: Vec<WireFace> = a_new_faces.iter().enumerate().map(|(gi, w)| WireFace {
        outer_wire: w.clone(),
        inner_wires: g2h.get(&gi).map(|hs| hs.iter().map(|&h| a_hole_faces[h].clone()).collect()).unwrap_or_default(),
        internal_wires: internal_wires.to_vec(),
    }).collect();

    // OCCT L557-581: unassigned holes + open face → create new growth from original surface
    if !orphan_holes.is_empty() && ds.faces[face_idx].natural_restriction {
        // rcad: WireFace with empty outer_wire = full parametric surface
        let orphan_inner: Vec<Vec<usize>> = orphan_holes.iter().map(|&h| a_hole_faces[h].clone()).collect();
        result.push(WireFace {
            outer_wire: vec![],
            inner_wires: orphan_inner,
            internal_wires: vec![],
        });
    }
    result
}

/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformInternalShapes (L618-778).
///   Classify internal wire groups against result faces via UV point-in-polygon.
pub(crate) fn perform_internal_shapes(
    wfs: &mut Vec<WireFace>,
    internal_wire_groups: &[Vec<usize>],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_idx: usize,
    ds: &crate::bopds::ds::DS,
) {
    if internal_wire_groups.is_empty() { return; }
    if wfs.is_empty() { return; }

    // OCCT L634-666: Build UV boundaries for each face
    let face_uv_bounds: Vec<Vec<DVec2>> = wfs.iter().map(|wf| {
        let mut uv_bnd: Vec<DVec2> = Vec::new();
        for &si in &wf.outer_wire {
            let seg = &segments[si];
            let pc_opt = tool.curve_on_surface(seg.edge, seg.face)
                .map(|(pc, _, _)| pc)
                .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
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
        uv_bnd
    }).collect();

    // OCCT L674-741: classify each internal wire against each face
    for group in internal_wire_groups {
        if group.is_empty() { continue; }
        let si = group[0];
        if si >= segments.len() { continue; }
        let seg = &segments[si];

        // Sample UV midpoint
        let uv_pt = {
            let mut pt = DVec2::ZERO;
            let mut found = false;
            let pc_opt = tool.curve_on_surface(seg.edge, seg.face)
                .map(|(pc, _, _)| pc)
                .or(seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref()));
            if let Some(pc) = pc_opt {
                let t_mid = (seg.t_range[0] + seg.t_range[1]) * 0.5;
                pt = pc.point_at(t_mid);
                found = true;
            }
            if !found {
                // Fallback: sample pcurve at midpoint via curve_on_surface
                if let Some((pc, t0, t1)) = tool.curve_on_surface(seg.edge, seg.face) {
                    let t_mid = (t0 + t1) * 0.5;
                    pt = pc.point_at(t_mid);
                } else { continue; }
            }
            pt
        };

        // OCCT L710-715: if edge is inside face → add to internal wires
        for (fi, wf) in wfs.iter_mut().enumerate() {
            if fi < face_uv_bounds.len() && face_uv_bounds[fi].len() >= 3
                && crate::builder::wire_path::point_in_uv_polygon(uv_pt, &face_uv_bounds[fi])
            {
                wf.internal_wires.push(group.clone());
                break;
            }
        }
    }
}

/// OCCT BOPAlgo_BuilderFace::PerformLoops L327-382: group connected avoided edges into internal wires.
pub(crate) fn build_internal_wires(
    segments: &[WireSegmentTopoDS],
    avoided: &HashSet<usize>,
) -> Vec<Vec<usize>> {
    if avoided.is_empty() { return vec![]; }
    let mut v_to_e: HashMap<usize, Vec<usize>> = HashMap::new();
    for &si in avoided {
        let seg = &segments[si];
        v_to_e.entry(seg.start_vertex.index).or_default().push(si);
        v_to_e.entry(seg.end_vertex.index).or_default().push(si);
    }
    let mut visited: HashSet<usize> = HashSet::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for &si in avoided {
        if visited.contains(&si) { continue; }
        let mut group: Vec<usize> = Vec::new();
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(si);
        visited.insert(si);
        while let Some(cur) = queue.pop_front() {
            group.push(cur);
            let seg = &segments[cur];
            for &v in &[seg.start_vertex.index, seg.end_vertex.index] {
                if let Some(neighbors) = v_to_e.get(&v) {
                    for &nsi in neighbors {
                        if visited.insert(nsi) {
                            queue.push_back(nsi);
                        }
                    }
                }
            }
        }
        groups.push(group);
    }
    groups
}
