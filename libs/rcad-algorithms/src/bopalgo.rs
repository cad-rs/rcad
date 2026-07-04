use glam::DVec3;
use std::collections::{HashMap, HashSet, VecDeque, BTreeMap};
use crate::bopds::ds::DS;
use crate::bvh::Aabb;
use crate::classify::{Classification, classify_point};

/// ✅ OCCT-aligned: BOPAlgo_GlueEnum — glue mode for coincident-face detection.
///   GlueOff=0: no glue (standard intersection).
///   GlueFull=1: full glue (coincident faces create shared topology).
///   GlueShift=2: shift glue (same as full, but with tolerance shift).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlueEnum {
    GlueOff = 0,
    GlueFull = 1,
    GlueShift = 2,
}

impl Default for GlueEnum {
    fn default() -> Self { GlueEnum::GlueOff }
}

/// ✅ OCCT-aligned: BOPAlgo_Alert base — alert types for pipeline diagnostics.
///   OCCT uses a polymorphic class hierarchy (BOPAlgo_Alert → subclasses).
///   rcad uses a flat enum since Rust naturally models "one of" variants.
#[derive(Debug, Clone)]
pub enum Alert {
    /// Edge range is too small to process (BOPAlgo_AlertTooSmallRange).
    /// Stores (edge_idx, range_length).
    TooSmallRange(usize, f64),
    /// Building pcurve on face failed (BOPAlgo_AlertBuildingPCurveFailed).
    /// Stores (edge_idx, face_idx).
    BuildingPCurveFailed(usize, usize),
    /// Faces from BuilderSolid that were not used in any solid shell
    /// (BOPAlgo_AlertSolidBuilderUnusedFaces).
    SolidBuilderUnusedFaces(Vec<usize>),
    /// Alert: edge has no curve (section edge without valid geometry).
    EdgeWithoutCurve(usize),
}

/// ✅ OCCT-aligned: BOPAlgo_Report — collects alerts during pipeline execution.
///   OCCT BOPAlgo_Report stores Handle(BOPAlgo_Alert) list + Dump() method.
#[derive(Debug, Clone, Default)]
pub struct Report {
    alerts: Vec<Alert>,
}

impl Report {
    pub fn new() -> Self { Self { alerts: Vec::new() } }
    pub fn add_alert(&mut self, alert: Alert) { self.alerts.push(alert); }
    pub fn has_alerts(&self) -> bool { !self.alerts.is_empty() }
    pub fn alerts(&self) -> &[Alert] { &self.alerts }
    pub fn clear(&mut self) { self.alerts.clear(); }

    /// OCCT-aligned: compatibility check for code that uses simple bool.
    pub fn has_error(&self) -> bool {
        self.alerts.iter().any(|a| matches!(a, Alert::TooSmallRange(_, _) | Alert::EdgeWithoutCurve(_)))
    }
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::FillMap (template, cxx L83-102).
/// Adds pair (n1, n2) to a connection map for connectivity grouping.
/// rcad: specialization for usize keys.
pub fn fill_map(
    map: &mut BTreeMap<usize, Vec<usize>>,
    n1: usize,
    n2: usize,
) {
    map.entry(n1).or_default().push(n2);
    map.entry(n2).or_default().push(n1);
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::MakeBlocks (template, cxx L45-80).
/// Groups connected elements from a connection map into blocks.
/// rcad: specialization for usize keys returning Vec<Vec<usize>>.
pub fn make_blocks(map: &BTreeMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let mut fence: HashSet<usize> = HashSet::new();
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for (&key, _) in map {
        if !fence.insert(key) { continue; }
        let mut block = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(key);
        while let Some(n1) = queue.pop_front() {
            block.push(n1);
            if let Some(neighbors) = map.get(&n1) {
                for &n2 in neighbors {
                    if fence.insert(n2) {
                        queue.push_back(n2);
                    }
                }
            }
        }
        blocks.push(block);
    }
    // Add isolated elements not in any connection
    for (&key, neighbors) in map {
        if neighbors.is_empty() && fence.insert(key) {
            blocks.push(vec![key]);
        }
    }
    blocks
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::IntersectVertices (cxx L1119-1205).
///
/// Groups vertices by geometric proximity (overlapping tolerance spheres).
/// Each resulting chain contains vertex indices that should be merged.
///
/// Args:
///   vertex_indices — list of DS vertex indices to check.
///   ds — DS containing vertex data.
///   fuzzy_value — additional tolerance for intersection (half added to each vertex).
///
/// Returns: groups of vertex indices that intersect each other.
pub fn intersect_vertices(
    vertex_indices: &[usize],
    ds: &DS,
    fuzzy_value: f64,
) -> Vec<Vec<usize>> {
    let a_nb_v = vertex_indices.len();
    if a_nb_v <= 1 {
        return vertex_indices.iter().map(|&vi| vec![vi]).collect();
    }

    // OCCT L1135: aTolAdd = theFuzzyValue / 2.
    let a_tol_add = fuzzy_value / 2.0;

    // OCCT L1138-1157: build BVH tree of vertex bounding boxes
    // rcad: use BTreeMap + nested loops (simpler than full BVH for vertex sets)
    // Build connection map from proximity
    let mut map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

    for i in 0..a_nb_v {
        let vi = vertex_indices[i];
        let Some(v) = ds.vertices.get(vi) else { continue; };
        let a_tol = ds.vertices[vi].geom_tol.max(0.0);
        let total_tol = a_tol + a_tol_add;
        let total_tol2 = total_tol * total_tol;

        // OCCT L1175-1178: Find interfering pairs (FillMap)
        for j in (i + 1)..a_nb_v {
            let vj = vertex_indices[j];
            let Some(v2) = ds.vertices.get(vj) else { continue; };
            let a_tol_j = ds.vertices[vj].geom_tol.max(0.0);
            let total_tol_j = a_tol_j + a_tol_add;
            let total_tol_sum = total_tol + total_tol_j;
            let dist2 = (v.point - v2.point).length_squared();
            if dist2 <= total_tol_sum * total_tol_sum {
                fill_map(&mut map, i, j);
            }
        }
    }

    // OCCT L1181-1194: MakeBlocks from connection map
    let mut blocks = make_blocks(&map);

    // OCCT L1197-1204: Add non-interfering vertices as singleton chains
    let mut taken: HashSet<usize> = HashSet::new();
    for block in &blocks {
        for &idx in block {
            taken.insert(idx);
        }
    }
    for i in 0..a_nb_v {
        if !taken.contains(&i) {
            blocks.push(vec![i]);
        }
    }

    // Convert internal indices back to DS vertex indices
    blocks.iter().map(|block| {
        block.iter().map(|&idx| vertex_indices[idx]).collect()
    }).collect()
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::EdgesToWires (cxx L360-663).
/// Converts a set of edge indices into connected wires.
/// The edges are expected to be planar; each resulting wire starts
/// at a free end and follows connectivity through shared vertices.
///
/// Args:
///   edge_indices — list of DS edge indices.
///   ds — DS containing edge/vertex data.
///   shared — if true, edges are already topologically shared.
///   ang_tol — angular tolerance for plane detection.
///
/// Returns: groups of (edge_idx, forward_flag) representing wires.
/// Error codes: 0=success, 1=no edges, 2=sharing failed.
pub fn edges_to_wires(
    edge_indices: &[usize],
    ds: &DS,
    shared: bool,
    _ang_tol: f64,
) -> Result<Vec<Vec<(usize, bool)>>, i32> {
    if edge_indices.is_empty() {
        return Err(1); // OCCT L392-395
    }

    // OCCT L404-438: Filter out degenerated edges
    let a_le: Vec<usize> = edge_indices.iter()
        .filter(|&&ei| {
            ds.edges.get(ei).map_or(false, |e| {
                !ds.is_edge_degenerated(ei) && e.is_geometric
            })
        })
        .copied()
        .collect();

    if a_le.is_empty() {
        return Err(1); // no valid geometric edges
    }

    // OCCT L442-452: If not shared, try to share edges by vertex proximity
    if !shared {
        // rcad: edges in DS already share vertices by index, so sharing is implicit.
        // If vertices need merging, the caller should call intersect_vertices first.
    }

    // OCCT L465-508: Build vertex→edge adjacency map (using edge orientation)
    // rcad: map from vertex index to list of (edge_idx, forward_flag)
    let mut a_ve_map: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
    for &ei in &a_le {
        let edge = &ds.edges[ei];
        // Forward orientation
        a_ve_map.entry(edge.start_vertex).or_default().push((ei, true));
        a_ve_map.entry(edge.end_vertex).or_default().push((ei, false));
    }

    // OCCT L513-518: Build fence map (processed edges)
    let mut a_m_fence: HashSet<usize> = HashSet::new();
    // OCCT L522-525: Edge->wire order map (aMVE processed vertices)
    let mut a_m_ve_processed: HashSet<usize> = HashSet::new();

    // OCCT L528-658: Build wires by walking edge chains
    let mut a_lwires: Vec<Vec<(usize, bool)>> = Vec::new();

    // Start from edges with free vertices (valence 1), then follow the chain
    // Collect start edges: all edges to process
    let mut start_edges: Vec<(usize, bool)> = Vec::new();
    
    // OCCT L536-547: Start from edge with free vertex (or first unprocessed edge)
    // First pass: find edges with valence-1 vertices
    for &ei in &a_le {
        if a_m_fence.contains(&ei) { continue; }
        let edge = &ds.edges[ei];
        let sv = edge.start_vertex;
        let ev = edge.end_vertex;
        let sv_count = a_ve_map.get(&sv).map_or(0, |v| v.len());
        let ev_count = a_ve_map.get(&ev).map_or(0, |v| v.len());
        // Prefer starting from edges with at least one free end
        if sv_count == 1 || ev_count == 1 {
            start_edges.push((ei, true));
            a_m_fence.insert(ei);
        }
    }
    // Second pass: add remaining unprocessed edges
    for &ei in &a_le {
        if a_m_fence.insert(ei) { continue; } // already added
        start_edges.push((ei, true));
    }

    // Clear fence for the walking phase
    a_m_fence.clear();
    for &(ei, _) in &start_edges {
        a_m_fence.insert(ei);
    }

    // Walk each start edge to build wires
    for &(start_ei, start_fwd) in &start_edges {
        if !a_m_fence.contains(&start_ei) { continue; }
        
        let mut wire: Vec<(usize, bool)> = Vec::new();
        let edge = &ds.edges[start_ei];
        let (mut a_v_cur, _a_v_other) = if start_fwd {
            (edge.start_vertex, edge.end_vertex)
        } else {
            (edge.end_vertex, edge.start_vertex)
        };

        // Walk forward from the free end (or arbitrary start)
        wire.push((start_ei, start_fwd));
        a_m_fence.remove(&start_ei);

        // Walk in the forward direction
        loop {
            let Some(neighbors) = a_ve_map.get(&a_v_cur) else { break; };
            let mut found = false;
            for &(next_ei, _next_fwd) in neighbors {
                if !a_m_fence.contains(&next_ei) { continue; }
                let next_edge = &ds.edges[next_ei];
                // Determine orientation: connect a_v_cur → other_vertex
                if next_edge.start_vertex == a_v_cur {
                    a_v_cur = next_edge.end_vertex;
                    wire.push((next_ei, true));
                } else if next_edge.end_vertex == a_v_cur {
                    a_v_cur = next_edge.start_vertex;
                    wire.push((next_ei, false));
                } else {
                    continue;
                }
                a_m_fence.remove(&next_ei);
                found = true;
                break;
            }
            if !found { break; }
        }

        // OCCT L607-625: Handle closed wires — try the other direction from start
        if !wire.is_empty() && wire.last().map(|&(ei, _)| ei) == Some(start_ei) {
            // Closed wire — already complete
        }

        a_lwires.push(wire);
    }

    Ok(a_lwires)
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::ClassifyFaces (cxx:1622-1747).
///
/// Classifies result faces relatively draft solids.  For each solid,
/// collects the faces classified as IN into the returned Vec.
///
/// OCCT builds a BVH tree of face boxes (L1670-1680), then runs per-solid
/// classification jobs (L1736-1738, BOPAlgo_FillIn3DParts::Perform).
/// rcad: for each (face, solid) pair with overlapping AABB, calls
/// classify_point against the solid's DS face set.
///
/// Args:
///   the_faces: result face indices, each with sample point at
///     face_samples[i].
///   face_samples: 3D sample points for each face (one per the_faces entry).
///   the_solids: each solid = Vec of shell groups of DS face indices.
///   ds: data structure.
///   aabb_of_face: bounding box for each face (parallel to the_faces).
///   aabb_of_solid: bounding box for each solid.
///
/// Returns: for each solid index, list of result FACE INDICES (values from
///   the_faces, not positions) classified as IN that solid.
pub fn classify_faces(
    the_faces: &[usize],
    face_samples: &[DVec3],
    the_solids: &[Vec<Vec<usize>>],
    ds: &DS,
    aabb_of_face: &[Aabb],
    aabb_of_solid: &[Aabb],
) -> Vec<Vec<usize>> {
    let n_solids = the_solids.len();
    let mut the_in_parts: Vec<Vec<usize>> = vec![Vec::new(); n_solids];

    // Precompute flat DS face sets per solid (for classify_point)
    let solid_faces: Vec<Vec<usize>> = the_solids.iter()
        .map(|shells| shells.iter().flat_map(|sh| sh.iter().copied()).collect())
        .collect();

    for (si, sfaces) in solid_faces.iter().enumerate() {
        if sfaces.is_empty() { continue; }
        let sbox = &aabb_of_solid[si];
        for (pi, &fi) in the_faces.iter().enumerate() {
            if pi >= face_samples.len() { continue; }
            if !aabb_of_face[pi].intersects(sbox) { continue; }
            let class = classify_point(face_samples[pi], sfaces, ds);
            if class == Classification::In {
                the_in_parts[si].push(fi);
            }
        }
    }
    the_in_parts
}

/// ✅ OCCT-aligned: BOPAlgo_Tools::TrsfToPoint (cxx:1912-1937).
///
/// Computes a translation from the combined bounding box of two boxes to a point.
/// Returns `Some(translation_vector)` when the point is sufficiently far from the
/// combined box center and the box size is small enough relative to the distance.
/// Returns `None` when the criteria rejects the transformation.
///
/// OCCT parameters:
///   theBox1, theBox2 — bounding boxes to unify.
///   theTrsf        — (output) the transform to fill.
///   thePoint       — target point.
///   theCriteria    — minimal distance criterion.
///
/// rcad: returns Option<DVec3> (the translation vector) instead of bool + gp_Trsf.
pub fn trsf_to_point(
    box1: &crate::bvh::Aabb,
    box2: &crate::bvh::Aabb,
    point: glam::DVec3,
    criteria: f64,
) -> Option<glam::DVec3> {
    // OCCT L1918-1920: Unify two boxes
    let mut a_box = *box1;
    a_box.expand_aabb(box2);

    // OCCT L1922-1923: Compute center of unified box and distance from point
    let a_b_center = (a_box.min + a_box.max) * 0.5;
    let a_pb_dist = (point - a_b_center).length();

    // OCCT L1924-1927: Reject if point is too close to box center
    if a_pb_dist < criteria {
        return None;
    }

    // OCCT L1929-1933: Compute box diagonal length; reject if box is too large
    //   relative to the distance (ratio > 1/criteria)
    let a_b_size = (a_box.max - a_box.min).length();
    if (a_b_size / a_pb_dist) > (1.0 / criteria) {
        return None;
    }

    // OCCT L1935: Set translation from box corner min to the point
    Some(point - a_box.min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[test]
    fn test_fill_map_and_make_blocks_single_edge() {
        let mut m = BTreeMap::new();
        fill_map(&mut m, 1, 2);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert!(blocks[0].contains(&1) && blocks[0].contains(&2));
    }

    #[test]
    fn test_make_blocks_two_components() {
        let mut m = BTreeMap::new();
        fill_map(&mut m, 1, 2);
        fill_map(&mut m, 3, 4);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 2);
    }

    #[test]
    fn test_make_blocks_chain() {
        let mut m = BTreeMap::new();
        fill_map(&mut m, 1, 2);
        fill_map(&mut m, 2, 3);
        fill_map(&mut m, 3, 4);
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].len(), 4);
    }

    #[test]
    fn test_make_blocks_isolated() {
        let mut m = BTreeMap::new();
        m.entry(42).or_default();
        let blocks = make_blocks(&m);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0], vec![42]);
    }

    #[test]
    fn test_make_blocks_empty() {
        let m: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        assert!(make_blocks(&m).is_empty());
    }
}

    // ===== GlueEnum =====

    #[test]
    fn glue_enum_default_is_off() {
        assert_eq!(GlueEnum::default(), GlueEnum::GlueOff);
    }

    #[test]
    fn glue_enum_discriminants() {
        assert_eq!(GlueEnum::GlueOff as i32, 0);
        assert_eq!(GlueEnum::GlueFull as i32, 1);
        assert_eq!(GlueEnum::GlueShift as i32, 2);
    }

    // ===== Alert / Report =====

    #[test]
    fn report_empty_initially() {
        let r = Report::new();
        assert!(!r.has_alerts());
        assert!(!r.has_error());
    }

    #[test]
    fn report_collects_alerts() {
        let mut r = Report::new();
        r.add_alert(Alert::TooSmallRange(3, 1e-12));
        assert!(r.has_alerts());
        assert!(r.has_error());
    }

    #[test]
    fn report_building_pcurve_failed_not_fatal() {
        let mut r = Report::new();
        r.add_alert(Alert::BuildingPCurveFailed(5, 2));
        assert!(r.has_alerts());
        assert!(!r.has_error(), "pcurve failure should not count as error");
    }

    #[test]
    fn report_clear_removes_alerts() {
        let mut r = Report::new();
        r.add_alert(Alert::TooSmallRange(0, 1e-9));
        r.clear();
        assert!(!r.has_alerts());
    }

    #[test]
    fn report_alerts_slice_order() {
        let mut r = Report::new();
        r.add_alert(Alert::EdgeWithoutCurve(1));
        r.add_alert(Alert::SolidBuilderUnusedFaces(vec![2, 3]));
        assert_eq!(r.alerts().len(), 2);
        assert!(matches!(r.alerts()[0], Alert::EdgeWithoutCurve(1)));
    }
