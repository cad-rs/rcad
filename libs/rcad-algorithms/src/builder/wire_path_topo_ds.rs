/// TopoDS-based walk_path_extract_wires using BRepTool.
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
    find_angle_at,
};
use super::angle_2d::clock_wise_angle;
use crate::inttools::fclass2d::curve2d_nb_samples;
use super::point_in_polygon_2d;

// ---------------------------------------------------------------------------
// TopoDS-native EdgeInfo — no seg_idx, holds edge+face directly (OCCT aligned)
// ---------------------------------------------------------------------------

/// EdgeInfo for TopoDS-native path.
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
/// ✅ OCCT-aligned: BOPAlgo_WireSplitter::Path (WireSplitter_1.cxx L358-626).
/// Walks from start segment through smart_map forming closed wires by min-angle selection.
pub(crate) fn walk_path_extract_wires(
    start_si: usize,
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
) {
    let start_seg = &segments[start_si];

    // OCCT L382-384: init
    let mut a_va = start_seg.start_vertex.index;
    let mut a_e_outa = start_si;

    let mut edge_seq: Vec<usize> = Vec::new();
    let mut vert_seq: Vec<usize> = Vec::new();

    // OCCT Tolerance2D (WireSplitter_1.cxx L883-905)
    let face_ref = start_seg.face;
    let vtol = |vi: usize| -> f64 { tool.vertex_tolerance(ShapeRef::synthetic(vi)) };
    let u_tolerance = |vi: usize| -> f64 { tool.u_resolution(face_ref, vtol(vi)) };
    let v_tolerance = |vi: usize| -> f64 { tool.v_resolution(face_ref, vtol(vi)) };
    let tolerance_2d = |vi: usize| -> f64 {
        let vt = vtol(vi);
        let u = u_tolerance(vi);
        let v = v_tolerance(vi);
        let mut t = u.max(v);
        if t < vt { t = vt; }
        if matches!(tool.face_surface(face_ref), Some(&Surface3::BSpline(_))) { t *= 1.1; }
        t
    };

    // OCCT Coord2d (WireSplitter_1.cxx L677-683)
    let coord2d = |vi: usize, si: usize| -> DVec2 {
        let s = &segments[si];
        let t = tool.parameter_on_edge(ShapeRef::synthetic(vi), s.edge, s.face)
            .unwrap_or(if vi == s.start_vertex.index { s.t_range[0] } else { s.t_range[1] });
        let pc = tool.curve_on_surface(s.edge, s.face);
        pc.and_then(|(pc, _, _)| Some(pc.point_at(t))).unwrap_or(DVec2::ZERO)
    };

    // OCCT L392: for (;;)
    loop {
        // OCCT L394-403: Do not escape through edge from which you enter
        {
            let a_nb = edge_seq.len();
            if a_nb == 1 {
                let an_e_prev = edge_seq[a_nb - 1];
                if an_e_prev == a_e_outa {
                    return;
                }
            }
        }

        edge_seq.push(a_e_outa);

        let seg = &segments[a_e_outa];
        let a_vb = if seg.start_vertex.index == a_va { seg.end_vertex.index } else { seg.start_vertex.index };

        // OCCT L416-422: append to sequences (vert + UV coord)
        vert_seq.push(a_vb);
        let a_pb = coord2d(a_vb, a_e_outa);

        // OCCT L414: anEdgeInfo->SetPassed(true)
        mark_edge_passed(smart_map, a_e_outa, a_va, false);

        // OCCT L428-523: Loop detection
        let a_nb = edge_seq.len() - 1;  // number of entries BEFORE current append
        let a_tol_2d_2 = 2.0 * tolerance_2d(a_vb);
        let a_tol_2d_sq = a_tol_2d_2 * a_tol_2d_2;
        // Inline is_vert_closed to avoid borrow conflict
        let b_is_closed = smart_map.get(&a_vb).map_or(false, |infos| {
            infos.iter().any(|ei| {
                let s = &segments[ei.seg_idx];
                s.start_vertex.index == s.end_vertex.index || s.is_closed_on_face
            })
        });

        let mut loop_found = false;
        // OCCT L466: for (j = aNb; j >= 0; --j)
        // aNb is aLS.Length() before append (1-indexed). After append, aLS has aNb+1 entries.
        // j = aNb checks the last entry (the one just appended). In rcad 0-indexed,
        // the last entry is at index a_nb (edge_seq.len()-1). OCCT extracts aLS(j..aNb)
        // so j=aNb gives empty range → skip by using (0..a_nb).rev() (excludes a_nb).
        for j in (0..a_nb).rev() {
            // OCCT L468: aVertVa(j+1).IsSame(aVb)
            if vert_seq[j] != a_vb { continue; }

            // OCCT L472-474: if (aVertMap.IsBound(aVb) && aVertMap.Find(aVb))
            // → vertex belongs to closed/degenerated edge: check 2D proximity
            let is_same_v_2d = if b_is_closed {
                let a_d2 = coord2d(vert_seq[j], edge_seq[j]).distance_squared(a_pb);
                // OCCT L484-488: additional UV tolerance check
                if a_d2 < a_tol_2d_sq {
                    let prev_uv = coord2d(vert_seq[j], edge_seq[j]);
                    let u_dist = (prev_uv.x - a_pb.x).abs();
                    let v_dist = (prev_uv.y - a_pb.y).abs();
                    u_dist <= 2.0 * u_tolerance(a_vb) && v_dist <= 2.0 * v_tolerance(a_vb)
                } else { false }
            } else { true };

            // OCCT L476: if (anIsSameV && anIsSameV2d) → found loop
            if is_same_v_2d {
                // OCCT L478-496: MakeWire from aLS(j..aNb)
                let wire: Vec<usize> = edge_seq[j..].to_vec();
                let same_edge = wire.len() == 2
                    && match (&segments[wire[0]].source, &segments[wire[1]].source) {
                        (WireEdgeSourceTopoDS::DsEdge(a), WireEdgeSourceTopoDS::DsEdge(b)) => a.index == b.index,
                        (WireEdgeSourceTopoDS::IntersectionCurve(a), WireEdgeSourceTopoDS::IntersectionCurve(b)) => a.index == b.index,
                        (WireEdgeSourceTopoDS::SeamEdge, WireEdgeSourceTopoDS::SeamEdge) => true,
                        _ => false,
                    };
                if !same_edge { wires.push(wire); }

                // OCCT L505-523: backtrack
                if j == 0 { return; }
                a_va = vert_seq[j - 1];
                a_e_outa = edge_seq[j - 1];
                edge_seq.truncate(j);
                vert_seq.truncate(j);
                loop_found = true;
                break;
            }
        }
        if loop_found { continue; }

        // OCCT L526-616: Outgoing edge selection
        let i_cnt = smart_map.get(&a_vb).map_or(0, |infos| {
            infos.iter().filter(|ei| !ei.passed && !ei.in_flag).count()
        });
        if i_cnt == 0 { return; }

        let an_angle_in = find_angle_at(smart_map, a_e_outa, a_vb, true).unwrap_or(0.0);
        let incoming_is_boundary = !matches!(segments[a_e_outa].source, WireEdgeSourceTopoDS::IntersectionCurve(_));

        let le_info: Vec<&EdgeInfo> = smart_map.get(&a_vb)
            .map(|v| v.iter().filter(|ei| !ei.passed && !ei.in_flag).collect())
            .unwrap_or_default();

        if le_info.is_empty() { return; }

        if le_info.len() == 1 {
            a_va = a_vb;
            a_e_outa = le_info[0].seg_idx;
            continue;
        }

        // OCCT L526-624: Outgoing edge selection (inlined for 1:1 alignment)
        let a_two_pi = std::f64::consts::TAU;
        let eps = std::f64::EPSILON;
        let mut a_min_angle = 100.0;
        let mut a_nb_ways_inside: i32 = 0;
        let mut p_only_way_in: Option<&EdgeInfo> = None;
        let mut p_edge_info: Option<&EdgeInfo> = None;

        for an_ei in &le_info {
            // OCCT L541-542: anIsOut=!IsIn, anIsNotPassed=!Passed → already filtered above
            // OCCT L551-562: iCnt==0/1 handled above

            // OCCT L564-591: compute angle for this candidate
            let a_angle = if an_ei.seg_idx == a_e_outa {
                // OCCT L566: same edge → maximum angle (lowest priority)
                a_two_pi
            } else {
                // OCCT L570-582: bIsClosed 2D distance check (Coord2dVf)
                if b_is_closed {
                    let candidate_seg = &segments[an_ei.seg_idx];
                    let pc_opt = tool.curve_on_surface(candidate_seg.edge, candidate_seg.face)
                        .map(|(pc, _, _)| pc)
                        .or(candidate_seg.first_pcurve.as_ref()
                            .or(candidate_seg.second_pcurve.as_ref()));
                    if let Some(pc) = pc_opt {
                        let t_mid = (candidate_seg.t_range[0] + candidate_seg.t_range[1]) * 0.5;
                        let a_p2_dx = pc.point_at(t_mid);
                        let a_d2 = a_p2_dx.distance_squared(a_pb);
                        if a_d2 > a_tol_2d_sq {
                            continue; // OCCT L580-581
                        }
                    }
                }
                // OCCT L584-586: anAngleOut + ClockWiseAngle
                let an_angle_out = an_ei.angle;
                clock_wise_angle(an_angle_in, an_angle_out)
            };

            // OCCT L594-598: boundary/inside tracking
            if incoming_is_boundary && an_ei.is_inside {
                a_nb_ways_inside += 1;
                p_only_way_in = Some(an_ei);
            }

            // OCCT L600-604: minimum angle selection
            if a_angle < a_min_angle - eps {
                a_min_angle = a_angle;
                p_edge_info = Some(an_ei);
            }
        }

        // OCCT L607-610: single inside way → forced selection
        if a_nb_ways_inside == 1 {
            p_edge_info = p_only_way_in;
        }

        // OCCT L612-616: no way out → return
        if let Some(best_ei) = p_edge_info {
            a_va = a_vb;
            a_e_outa = best_ei.seg_idx;
        } else {
            return;
        }
    }
}

/// RefineAngles (WireSplitter_1.cxx L919-1043).
/// For each multi-vertex with 2 boundary edges, adjust internal edge angles
/// that fall outside the boundary sweep.  Uses BRepTool for angle computation.
fn refine_angles(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    a_ms: &HashSet<usize>,
    src_key_of: impl Fn(&WireSegmentTopoDS) -> usize,
) {
    let vertices: Vec<usize> = smart_map.keys().copied().collect();
    for &v in &vertices {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        // OCCT L308: aEI.SetIsInside(!aMS.Contains(aE))
        // IsInside: edge whose src_key is NOT in a_ms (appeared odd times) → interior edge
        {
            let infos_mut = smart_map.get_mut(&v).unwrap();
            for ei in infos_mut.iter_mut() {
                let seg = &segments[ei.seg_idx];
                let src_key = src_key_of(seg);
                ei.is_inside = !a_ms.contains(&src_key);
            }
        }
        // Refresh infos after mutation
        let infos = smart_map.get(&v).unwrap().clone();
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
            let v_ref = rcad_kernel::topods::ShapeRef::synthetic(v);
            let geom_tol = tool.vertex_tolerance(v_ref);
            let t_v = tool.parameter_on_edge(v_ref, seg.edge, seg.face)
                .unwrap_or_else(|| if v == seg.start_vertex.index { seg.t_range[0] } else { seg.t_range[1] });
            let domain = seg.t_range;
            let pc = seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref());
            let b_refined = match pc {
                Some(pc) => {
                    let ta = pc.point_at(t_v);
                    
                    // Check if angle is near boundary — true refine via micro-step
                    let eps = 1e-12;
                    let dt = (1e-8 * (domain[1] - domain[0]).abs().max(1.0)).min((domain[1] - domain[0]).abs() * 0.1);
                    let t2 = (t_v + dt).min(domain[1]);
                    let pt2 = pc.point_at(t2);
                    let dir = pt2 - ta;
                    if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE { false }
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

/// TopoDS-based SplitBlock — refine angles + path walk for irregular blocks.
pub(crate) fn split_block(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    wires: &mut Vec<Vec<usize>>,
    tool: &dyn BRepTool,
    a_ms: &HashSet<usize>,
) {
    // OCCT L197-227: vertex-level check + bNothingToDo
    // (rcad: skip bNothingToDo optimization — always run Path extraction)

    // OCCT L227-285: Path extraction
    let order_keys: Vec<usize> = smart_map.keys().copied().collect();

    // src_key_of: maps segment → physical edge key (same as in build_smart_map)
    let src_key_of = |seg: &WireSegmentTopoDS| -> usize {
        match seg.source {
            WireEdgeSourceTopoDS::DsEdge(ref sr) => sr.index,
            WireEdgeSourceTopoDS::IntersectionCurve(ref sr) => sr.index,
            WireEdgeSourceTopoDS::SeamEdge => usize::MAX,
        }
    };

    for &v in &order_keys {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        for ei in &infos {
            if !ei.passed && !ei.in_flag
                && ei.seg_idx < segments.len()
            {
                walk_path_extract_wires(ei.seg_idx, segments, smart_map, wires, tool);
            }
        }
    }

    // OCCT L298-318: Angles in mySmartMap (COMPUTE_ANGLES)
    // rcad: refine_angles matches OCCT's angle refinement after Path.
    // Must be called AFTER Path extraction to match OCCT order.
    refine_angles(smart_map, segments, tool, a_ms, src_key_of);
}

/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformLoops (BOPAlgo_BuilderFace.cxx L239-383).
/// TopoDS-based build_closed_wires — SmartMap + angle computation + wire walking.
///
/// Simplified version without vi_to_canon/deg_end_canon (ShapeRef handles use DS indices directly).
pub(crate) fn build_closed_wires(
    segments: &[WireSegmentTopoDS],
    avoided: &HashSet<usize>,
    tool: &dyn BRepTool,
) -> Vec<Vec<usize>> {
    if segments.is_empty() { return vec![]; }

    let n = segments.len();
    let mut vert_to_segs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        if avoided.contains(&si) { continue; }
        vert_to_segs.entry(seg.start_vertex.index).or_default().push(si);
        vert_to_segs.entry(seg.end_vertex.index).or_default().push(si);
    }

    if avoided.len() == n && !segments.is_empty() { return vec![]; }

    // Build connexity blocks
    let blocks = make_connexity_blocks(segments, avoided, &vert_to_segs, n);

    // Process each block
    let mut wires: Vec<Vec<usize>> = Vec::new();
    for (bi, block) in blocks.iter().enumerate() {
        if block.len() < 2 {
            continue;
        }

        // Build SmartMap
        let (smart_map, _a_vert_map, a_ms) = build_smart_map(block, segments, tool);
        if smart_map.is_empty() {
            continue;
        }

        // OCCT SplitBlock L227-285: Path extraction (always, no regular fast path)
        split_block(block, segments, &mut (smart_map.clone()), &mut wires, tool, &a_ms);
    }
    wires
}

/// ✅ OCCT-aligned: BOPAlgo_WireSplitter::SplitBlock (WireSplitter_1.cxx L136-195, L298-318).
/// Fills mySmartMap: for each edge in block, traverse vertices and build EdgeInfo list.
/// Also builds aVertMap (vertex closed/degenerated tracking).
fn build_smart_map(
    block: &[usize],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
) -> (IndexMap<usize, Vec<EdgeInfo>>, HashMap<usize, bool>, HashSet<usize>) {
    // OCCT L134: NCollection_Map<TopoDS_Shape> aMS — edge parity tracking
    // Orientation-independent key: same physical edge = same key (OCCT TopoDS_Shape ignores orientation)
    let src_key_of = |seg: &WireSegmentTopoDS| -> usize {
        match seg.source {
            WireEdgeSourceTopoDS::DsEdge(ref sr) => sr.index,
            WireEdgeSourceTopoDS::IntersectionCurve(ref sr) => sr.index,
            WireEdgeSourceTopoDS::SeamEdge => usize::MAX,
        }
    };
    let mut a_ms: HashSet<usize> = HashSet::new();

    let mut smart_map: IndexMap<usize, Vec<EdgeInfo>> = IndexMap::new();
    let mut a_vert_map: HashMap<usize, bool> = HashMap::new();

    // Track which physical edges have been added at each vertex, to deduplicate
    // FORWARD/REVERSED copies of the same physical edge matching OCCT TopoDS_Shape behavior.
    let mut vert_edge_added: HashMap<usize, HashSet<usize>> = HashMap::new();

    // OCCT L154: aV1 tracks first vertex for closed-edge detection
    let mut a_v1: Option<usize> = None;

    // OCCT L137-195: 1.Filling mySmartMap
    for &si in block {
        let seg = &segments[si];
        // OCCT L141-144: skip edges without pcurve on face
        let has_pcurve = tool.curve_on_surface(seg.edge, seg.face).is_some()
            || seg.first_pcurve.is_some() || seg.second_pcurve.is_some();
        if !has_pcurve { continue; }

        // OCCT L146: bIsClosed = Degenerated(aE) || IsClosed(aE, myFace)
        let b_is_closed = seg.start_vertex.index == seg.end_vertex.index
            || seg.is_closed_on_face;

        // OCCT L148-151: aMS parity (insert; if already exists and not closed, remove)
        let src_key = src_key_of(seg);
        if !a_ms.insert(src_key) && !b_is_closed {
            a_ms.remove(&src_key);
        }

        // OCCT L154-194: iterate vertices of the edge
        // aE has 2 vertices: i=0 (first), i=1 (last)
        // OCCT L168-171: aEI.SetEdge(aE); aOr = aV.Orientation(); bIsIN = (aOr == REVERSED)
        // In OCCT, TopExp_Explorer of an edge returns vertices in TShape order.
        // aV.Orientation() is the vertex's orientation WITHIN the edge:
        //   First vertex (i=0): FORWARD → bIsIN=false
        //   Last vertex (i=1): REVERSED → bIsIN=true
        // For REVERSED orientation of the edge in the face wire, the vertex
        // orientations are swapped by TopExp_Explorer:
        //   First vertex: REVERSED → bIsIN=true
        //   Last vertex: FORWARD → bIsIN=false

        // First vertex (i=0): orientation matches edge orientation in face
        {
            let a_v = seg.start_vertex.index;
            // OCCT: aV.Orientation() when edge iter starts.
            // For a segment: start_vertex orientation = FORWARD if seg is FORWARD, REVERSED if seg is REVERSED
            // OCCT TopExp: FORWARD edge → first vertex FORWARD(OUT), last REVERSED(IN)
            //              REVERSED edge → first vertex REVERSED(IN), last FORWARD(OUT)
            // rcad swaps start/end for REVERSED, so:
            //   REVERSED start = TShape last vertex → FORWARD → OUT
            //   REVERSED end = TShape first vertex → REVERSED → IN
            //   FORWARD start = TShape first vertex → FORWARD → OUT
            //   FORWARD end = TShape last vertex → REVERSED → IN
            // Result: b_is_in=false for start, true for end, regardless of orientation.
            let b_is_in = false;

            // Deduplicate: only add if this physical edge hasn't been added at this vertex
            let added = vert_edge_added.entry(a_v).or_default();
            if added.insert(src_key) {
                let entry = smart_map.entry(a_v).or_default();
                entry.push(EdgeInfo {
                    seg_idx: si, passed: false, in_flag: b_is_in,
                    is_inside: false, is_circle_arc: false, angle: 0.0,
                });
            }

            a_v1 = Some(a_v);

            // OCCT L183-193: aVertMap
            let closed = b_is_closed;
            a_vert_map.entry(a_v)
                .and_modify(|v| { if closed { *v = true; } })
                .or_insert(closed);
        }

        // Last vertex (i=1): opposite orientation
        {
            let a_v = seg.end_vertex.index;

            // Deduplicate: only add if this physical edge hasn't been added at this vertex
            let added = vert_edge_added.entry(a_v).or_default();
            if added.insert(src_key) {
                let entry = smart_map.entry(a_v).or_default();
                // For FORWARD: end = TShape last vertex -> REVERSED -> bIsIN=true
                // For REVERSED: end = TShape first vertex -> REVERSED -> bIsIN=true
                let b_is_in = true;
                entry.push(EdgeInfo {
                    seg_idx: si, passed: false, in_flag: b_is_in,
                    is_inside: false, is_circle_arc: false, angle: 0.0,
                });
            }

            // OCCT L178-181: bIsClosed = bIsClosed || aV1.IsSame(aV)
            let a_v1_v = a_v1.unwrap();
            let b_is_closed = b_is_closed || a_v1_v == a_v;

            // OCCT L183-193: aVertMap
            let closed = b_is_closed;
            a_vert_map.entry(a_v)
                .and_modify(|v| { if closed { *v = true; } })
                .or_insert(closed);
        }
    }

    // OCCT L298-318: 3.Angles in mySmartMap
    // Compute angles using BRepTool (equivalent to OCCT Angle2D)
    for (v, infos) in smart_map.iter_mut() {
        let v_ref = ShapeRef::synthetic(*v);
        let geom_tol = tool.vertex_tolerance(v_ref);
        for ei in infos.iter_mut() {
            // OCCT L311-316: compute Angle2D
            // IsInside flag (OCCT L308) moved to refine_angles (after Path walk, matching OCCT order).
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
                            let (ta, tb) = super::wire_path::pc_parameter_range(pc);
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
            let default_surf = Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z));
            let surf = tool.face_surface(seg.face).unwrap_or(&default_surf);
            let new_angle = super::angle_2d::angle_2d(curve, t_v, curve_domain, ei.in_flag, surf, geom_tol, None)
                .unwrap_or(0.0);
            ei.angle = new_angle;
        }
    }

    (smart_map, a_vert_map, a_ms)
}

/// Build a regular wire from a block (all vertices have degree 2).
/// OCCT BOPAlgo_WireSplitter::MakeWire (WireSplitter_1.cxx L735-800).
/// Follows vertex adjacency to build a properly ordered cyclic wire.
/// NOTE: Unlike OCCT's MakeWire which adds edges from the list as-is,
/// this function reorders edges into chain order. rcad's make_connexity_blocks
/// uses BFS which may not produce chain order, so explicit ordering is needed.
/// OCCT-aligned: BOPAlgo_WireSplitter.lxx MakeWire (BRep_Builder::Add loop)
/// + implicit edge ordering from BOPTools_AlgoTools::MakeConnexityBlocks.
fn build_regular_wire(block: &[usize], segments: &[WireSegmentTopoDS]) -> Option<Vec<usize>> {
    if block.is_empty() { return None; }
    if block.len() == 1 {
        // Single edge → self-loop
        return Some(vec![block[0]]);
    }

    // Build vertex → segment adjacency
    let mut vert_to_segs: HashMap<usize, Vec<usize>> = HashMap::new();
    for &si in block {
        let seg = &segments[si];
        vert_to_segs.entry(seg.start_vertex.index).or_default().push(si);
        vert_to_segs.entry(seg.end_vertex.index).or_default().push(si);
    }

    // For a regular block, each vertex has exactly 2 incident segments.
    // Walk from the first segment's start vertex, following end→start chain.
    let mut wire: Vec<usize> = Vec::with_capacity(block.len());

    let start_si = block[0];
    let start_seg = &segments[start_si];
    let start_v = start_seg.start_vertex.index;
    wire.push(start_si);

    // Track used segment indices (NOT block positions).
    // BUG HISTORY: used[] was previously indexed by block position, but
    // vert_to_segs uses segment indices. When block order != [0,1,2..],
    // the wrong used entry was checked, causing duplicate visits.
    let mut used_si: std::collections::HashSet<usize> = std::collections::HashSet::new();
    used_si.insert(start_si);

    let mut cur_v = start_seg.end_vertex.index;
    loop {
        // Find the next unused segment with cur_v as start vertex
        let next = vert_to_segs.get(&cur_v).and_then(|neighbors| {
            neighbors.iter().find(|&&ni| !used_si.contains(&ni) && segments[ni].start_vertex.index == cur_v)
                .or_else(|| {
                    // Try matching end vertex (reversed orientation)
                    neighbors.iter().find(|&&ni| !used_si.contains(&ni) && segments[ni].end_vertex.index == cur_v)
                })
                .copied()
        });

        match next {
            Some(si) => {
                used_si.insert(si);
                wire.push(si);
                let seg = &segments[si];
                if seg.start_vertex.index == cur_v {
                    cur_v = seg.end_vertex.index;
                } else {
                    cur_v = seg.start_vertex.index;
                }
                // Check if we returned to start
                if cur_v == start_v && wire.len() == block.len() {
                    break;
                }
            }
            None => break, // No more edges to follow
        }
    }

    if wire.len() != block.len() {
        // Fallback: return block as-is if we couldn't form a complete chain
        return Some(block.to_vec());
    }

    Some(wire)
}

/// Connected-component grouping for TopoDS segments.
/// ✅ OCCT-aligned: BOPTools_AlgoTools::MakeConnexityBlocks (BOPTools_AlgoTools.cxx L105-154).
/// Architecture note: OCCT uses TopExp_Explorer + MapShapesAndAncestors for BFS; rcad uses
/// pre-built vert_to_segs HashMap. Algorithm is equivalent (undirected edge adjacency BFS).
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
/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformAreas (BOPAlgo_BuilderFace.cxx L387-614).
/// Classifies wires as growth (outer boundary) or hole (inner boundary),
/// assigns holes to the nearest enclosing growth face, and handles orphan holes
/// on unbounded (open) faces.
pub(crate) fn perform_areas(
    wires: &[Vec<usize>],
    internal_wires: &[Vec<usize>],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_idx: usize,
    ds: &crate::bopds::ds::DS,
) -> Vec<WireFace> {
    
    if wires.is_empty() {
        if ds.face_natural_restriction(face_idx) {
            return vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: vec![] }];
        }
        return vec![];
    }

    
    let mut a_mhe: HashSet<ShapeRef> = HashSet::new();

    
    let mut a_new_faces: Vec<Vec<usize>> = Vec::new(); // growth wires
    let mut a_hole_faces: Vec<Vec<usize>> = Vec::new(); // hole wires

    for w in wires {
        
        
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

        
        if b_is_growth {
            a_new_faces.push(w.clone());
        } else {
            a_hole_faces.push(w.clone());
            
            for &si in w {
                if let Some(seg) = segments.get(si) {
                    a_mhe.insert(seg.edge);
                }
            }
        }
    }

    
    if a_hole_faces.is_empty() {
        return a_new_faces.iter().map(|w| WireFace {
            outer_wire: w.clone(), inner_wires: vec![], internal_wires: internal_wires.to_vec(),
        }).collect();
    }

    
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

    
    let assigned_hole_set: std::collections::HashSet<usize> = h2g.iter().map(|&(h, _)| h).collect();
    let orphan_holes: Vec<usize> = (0..a_hole_faces.len())
        .filter(|hi| !assigned_hole_set.contains(hi))
        .collect();

    
    let mut g2h: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(h, g) in &h2g { g2h.entry(g).or_default().push(h); }

    
    //   OCCT BOPAlgo_BuilderFace::PerformAreas creates a separate TopoDS_Face
    //   for EACH wire (growth or hole).  The caller (ComputeState) then classifies
    //   each face independently, removing those that are In the opposing solid.
    //   rcad legacy approach merged holes as inner_wires, preventing per-face
    //   classification and leaving internal vertices in the result.
    
    let mut result: Vec<WireFace> = a_new_faces.iter().enumerate().map(|(gi, w)| WireFace {
        outer_wire: w.clone(),
        inner_wires: g2h.get(&gi).map(|hs| hs.iter().map(|&h| a_hole_faces[h].clone()).collect()).unwrap_or_default(),
        internal_wires: internal_wires.to_vec(),
    }).collect();

    // OCCT L514-544: Add unused holes to the original face (if the original face is open/infinite)
    // OCCT: aHoleFaces.Extent() != aHoleFaceMap.Extent() → orphan holes exist
    // OCCT: aBoxF.IsOpen* → original face is unbounded (e.g. plane)
    // rcad: !natural_restriction → face is unbounded/open (e.g. plane, not sphere)
    if !orphan_holes.is_empty() && !ds.face_natural_restriction(face_idx) {
        // OCCT: Create new face from original surface (no outer wire = infinite surface)
        //        Add all orphan holes as inner wires
        let orphan_inner: Vec<Vec<usize>> = orphan_holes.iter()
            .map(|&hi| a_hole_faces[hi].clone()).collect();
        result.push(WireFace {
            outer_wire: vec![],
            inner_wires: orphan_inner,
            internal_wires: vec![],
        });
    }
    // OCCT: If the original face is bounded (closed), orphan holes are simply left
    //       unassigned — they were not enclosed by any growth.  Same as rcad.
    result
}

/// BOPAlgo_BuilderFace::PerformInternalShapes (L618-778).
///   Classify internal wire groups against result faces via UV point-in-polygon.
/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformInternalShapes (BOPAlgo_BuilderFace.cxx L618-778).
/// Architecture note: OCCT uses myLoopsInternal (TopoDS_Wire list) + myAreas (TopoDS_Face list)
/// with BVH tree + IsInside classification. rcad uses internal_wire_groups (precomputed segment index
/// groups) + UV polygon containment. The classification logic is functionally equivalent but the data
/// model differs (TopoDS_Shape vs DS segment index).
pub(crate) fn perform_internal_shapes(
    wfs: &mut Vec<WireFace>,
    internal_wire_groups: &[Vec<usize>],
    segments: &[WireSegmentTopoDS],
    tool: &dyn BRepTool,
    face_idx: usize,
    face_ref: rcad_kernel::topods::ShapeRef,
    ds: &crate::bopds::ds::DS,
) {
    if internal_wire_groups.is_empty() { return; }
    if wfs.is_empty() { return; }

    
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
        // Fallback: if pcurve sampling gave <3 points, build UV polygon by
        // projecting outer wire vertices to UV using the face surface.
        if uv_bnd.len() < 3 && !wf.outer_wire.is_empty() {
            uv_bnd.clear();
            if let Some(face_surf) = tool.face_surface(rcad_kernel::topods::ShapeRef::synthetic(face_idx)) {
                for &si in &wf.outer_wire {
                    let seg = &segments[si];
                    let pt = tool.vertex_position(seg.start_vertex);
                    if let Some(uv) = crate::builder::wire_splitter::world_to_uv(face_surf, pt) {
                        uv_bnd.push(uv);
                    }
                }
                uv_bnd.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);
            }
        }
        uv_bnd
    }).collect();


    
    // OCCT L670-671: aMEDone — fence to track classified internal edge groups.
    let mut a_me_done: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let n_total_groups = internal_wire_groups.len();

    for group in internal_wire_groups {
        if group.is_empty() { continue; }
        let si = group[0];
        if si >= segments.len() { continue; }
        let seg = &segments[si];

        // OCCT L692-716: check if group already classified
        // (OCCT uses an edge-level aMEDone map; rcad uses group-level tracking.)
        let group_id = si; // Use first segment index as group identifier (unique across groups).
        if !a_me_done.insert(group_id) {
            continue;
        }

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
                    found = true;
                }
            }
            if !found {
                // Second fallback: project internal wire vertex to UV via face surface
            if let Some(face_surf) = tool.face_surface(face_ref) {
                    let v_pt = tool.vertex_position(seg.start_vertex);
                    pt = crate::builder::wire_splitter::world_to_uv(face_surf, v_pt)
                        .unwrap_or(DVec2::ZERO);
                    found = true;
                }
            }
            if !found { continue; }
            pt
        };

        
        for (fi, wf) in wfs.iter_mut().enumerate() {
            if fi < face_uv_bounds.len() && face_uv_bounds[fi].len() >= 3
                && crate::builder::wire_path::point_in_uv_polygon(uv_pt, &face_uv_bounds[fi])
            {
                    wf.internal_wires.push(group.clone());

                // OCCT L736-740: early exit if all internal edge groups classified.
                if a_me_done.len() == n_total_groups {
                    return;
                }
                break;
            }
        }
    }

    // OCCT L743-777: handle unclassified internal edge groups.
    // Architecture note: OCCT uses MakeInternalWires to create TopoDS_Wire with INTERNAL edges
    // and adds them to the face via BRep_Builder. rcad appends unclassified groups to the
    // internal_wires of all wire faces. OCCT also emits an AlertFaceBuilderUnusedEdges warning.
    if a_me_done.len() < n_total_groups {
        // Add unclassified groups to the first wire face (OCCT adds to myFace directly).
        if let Some(wf) = wfs.first_mut() {
            for group in internal_wire_groups {
                if group.is_empty() { continue; }
                let group_id = group[0];
                if a_me_done.contains(&group_id) { continue; }
                wf.internal_wires.push(group.clone());
            }
        }
        // OCCT L777: AddWarning(new BOPAlgo_AlertFaceBuilderUnusedEdges(aWShape))
        // rcad: unclassified internal edges warning not yet integrated.
    }
}

/// ✅ OCCT-aligned: BOPAlgo_BuilderFace::PerformLoops L327-382: group connected avoided edges into internal wires.
/// Architecture note: OCCT uses TopExp::MapShapesAndAncestors + aMAdded for BFS; rcad uses
/// v_to_e HashMap + visited HashSet. Algorithm is equivalent (undirected edge adjacency BFS).
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
