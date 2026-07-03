use std::collections::{HashMap, HashSet, VecDeque};
use indexmap::IndexMap;
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*;
use crate::bopds::ds::*; use crate::tolerance::*;
use crate::dbg_smartmap;
use super::types::{WireSegment, WireEdgeSource, WireFace, WireOrientation};
use super::wire_path::{pc_parameter_range, refine_angles, walk_path_extract_wires};
use crate::builder::point_in_polygon_2d;
use super::angle_2d::{dir_to_angle, angle_2d, clock_wise_angle};

/// 鉁?OCCT-aligned: Angle2D for seam edges (BOPAlgo_WireSplitter_1.cxx L768-840).
///
/// OCCT takes the edge's pcurve via BRep_Tool::CurveOnSurface, the vertex
/// parameter via BRep_Tool::Parameter, and calls Angle2D(aV, aE, aF, aGAS, bIsIN).
/// rcad equivalent: construct a Line pcurve along the surface isoline at the
/// parametric seam, then call angle_2d (which mirrors OCCT's dt/tol/step logic).
// 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
// 鈻?UV tangent angle computation
// 鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺愨晲鈺?
/// 鉁?OCCT-aligned Angle2D for DsEdge segments.
/// Evaluates the 3D curve at a micro-step near each vertex (OCCT
/// BOPAlgo_WireSplitter_1.cxx L768-840).  Maps both points to UV space
/// via world_to_uv, computes the direction.  Falls back to endpoint UV
/// difference when curve data is unavailable (plane is exact in both cases).
pub(crate) fn edge_uv_tangent(
    ds: &DS, sv: usize, ev: usize, surface: &Surface3,
    curve: Option<&Curve3>, t_range: Option<[f64; 2]>,
) -> (Option<f64>, Option<f64>) {
    // When curve data is available, use micro-step (OCCT Angle2D).
    // For plane surfaces, endpoint method is exact (linear pcurve).
    if let (Some(curve), Some(tr)) = (curve, t_range) {
        if !matches!(surface, Surface3::Plane(_)) {
            let fa = edge_angle_2d(curve, tr[0], tr, surface, false, ds.vertices[sv].geom_tol);
            let fb = edge_angle_2d(curve, tr[1], tr, surface, true, ds.vertices[ev].geom_tol);
            return (fa, fb);
        }
    }
    // Fallback: compute UV direction from endpoint UV difference.
    // Exact for plane surfaces; good approximation for small sub-edges.
    match surface {
        Surface3::Sphere(s) => {
            let uvs = s.world_to_uv(ds.vertices[sv].point);
            let uve = s.world_to_uv(ds.vertices[ev].point);
            let dir = uve - uvs;
            if dir.length_squared() < 1e-30 { return (None, None); }
            let a = dir_to_angle(dir);
            let na = a;
            (Some(na), Some((na + std::f64::consts::PI) % std::f64::consts::TAU))
        }
        Surface3::Plane(p) => {
            let x_axis = any_perpendicular(p.normal).normalize();
            let y_axis = p.normal.cross(x_axis).normalize();
            let local_s = ds.vertices[sv].point - p.origin;
            let local_e = ds.vertices[ev].point - p.origin;
            let uv_s = DVec2::new(local_s.dot(x_axis), local_s.dot(y_axis));
            let uv_e = DVec2::new(local_e.dot(x_axis), local_e.dot(y_axis));
            let dir = uv_e - uv_s;
            if dir.length_squared() < 1e-30 { return (None, None); }
            let a = dir_to_angle(dir);
            let na = a;
            (Some(na), Some((na + std::f64::consts::PI) % std::f64::consts::TAU))
        }
        _ => (None, None),
    }
}

/// 鉁?OCCT-aligned: micro-step Angle2D for a 3D curve mapped to face UV.
/// OCCT BOPAlgo_WireSplitter_1.cxx L768-840.  Evaluates the 3D curve at
/// t and t+dt, maps to UV via world_to_uv, returns UV direction angle.
pub(crate) fn edge_angle_2d(
    curve: &Curve3, t: f64, domain: [f64; 2],
    surface: &Surface3, b_is_in: bool, geom_tol: f64,
) -> Option<f64> {
    let range = (domain[1] - domain[0]).abs();
    if range < 1e-15 { return None; }
    // OCCT-aligned: Angle2D via 3D curve 鈫?UV mapping (WireSplitter_1.cxx L768-854).
    //   Unlike angle_2d.rs which takes a pcurve (Curve2d), this function evaluates
    //   the 3D curve and maps to UV via world_to_uv, then computes the UV tangent.
    let a_tol_2d = 2.0 * super::angle_2d::tolerance_2d(geom_tol, surface, None);
    let mut dt = a_tol_2d.max(1e-9); // OCCT L806: max with Precision::PConfusion (1e-9)
    // OCCT L808-821: curvature-aware adjustment for non-linear curves
    let eps = (1e-6 * range).max(1e-10);
    let tp = (t + eps).min(domain[1]);
    let tm = (t - eps).max(domain[0]);
    let p_p = curve.point_at(tp);
    let p_m = curve.point_at(tm);
    let d1 = p_p - p_m;
    let speed = d1.length();
    if speed > 1e-30 {
        let d1_n = d1 / speed;
        let d2 = p_p - 2.0 * curve.point_at(t) + p_m;
        let curvature = d1_n.cross(d2).length() / (speed * speed);
        if curvature > 1e-30 {
            let r_curv = 1.0 / curvature;
            let cos_phi = r_curv / (r_curv + a_tol_2d);
            if cos_phi < 1.0 {
                dt = dt.max(cos_phi.acos().max(1e-9));
            }
        }
    }
    // OCCT L824-834: clamp dt to 5% of range, floor at 5e-5
    let a_tx = 0.05 * range;
    let a_tx = if a_tx < 5e-5_f64 { (5e-5_f64).min(range * 0.5) } else { a_tx };
    if dt > a_tx { dt = a_tx; }
    // OCCT L822-829: step toward nearest curve end
    let t1 = if (t - domain[0]).abs() < (t - domain[1]).abs() {
        (t + dt).min(domain[1])
    } else {
        (t - dt).max(domain[0])
    };
    let p0 = curve.point_at(t);
    let p1 = curve.point_at(t1);
    let uv0 = world_to_uv(surface, p0)?;
    let uv1 = world_to_uv(surface, p1)?;
    let dir = if b_is_in { uv0 - uv1 } else { uv1 - uv0 };
    if dir.length_squared() < 1e-40 { return None; }
    Some(dir_to_angle(dir))
}

/// OCCT BOPTools_AlgoTools2D::IsEdgeIsoline (BOPTools_AlgoTools2D.cxx L669-700).
/// Checks if an edge's pcurve follows a U or V isoparametric line of the face surface.
/// Returns (is_uiso, is_v_iso) 鈥?true when the pcurve tangent at midpoint is
/// aligned with the U or V parametric direction respectively.
pub(crate) fn is_edge_isoline(pcurve: &Curve2d, range: [f64; 2]) -> (bool, bool) {
    let t_mid = (range[0] + range[1]) * 0.5;
    let eps = 1e-8 * (range[1] - range[0]).abs().max(1.0);
    let p_plus = pcurve.point_at((t_mid + eps).min(range[1]));
    let p_minus = pcurve.point_at((t_mid - eps).max(range[0]));
    let tangent = p_plus - p_minus;
    if tangent.length_squared() < 1e-30 {
        return (false, false);
    }
    let tangent_n = tangent.normalize();
    let tol = 1e-12;
    // OCCT L696-699: CrossMagnitude((0,1)) = |t.x|, CrossMagnitude((1,0)) = |t.y|
    let is_uiso = tangent_n.x.abs() <= tol;
    let is_v_iso = tangent_n.y.abs() <= tol;
    (is_uiso, is_v_iso)
}

/// Map a 3D point to UV space on a surface.  Returns None for unsupported
/// surface types (currently Sphere, Plane, Cylinder, Cone, Torus supported).
pub(crate) fn world_to_uv(surface: &Surface3, pt: DVec3) -> Option<DVec2> {
    match surface {
        Surface3::Sphere(s) => Some(s.world_to_uv(pt)),
        Surface3::Plane(p) => {
            let x_axis = any_perpendicular(p.normal).normalize();
            let y_axis = p.normal.cross(x_axis).normalize();
            let local = pt - p.origin;
            Some(DVec2::new(local.dot(x_axis), local.dot(y_axis)))
        }
        Surface3::Cylinder(c) => {
            let axis = c.axis.normalize_or_zero();
            if axis.length_squared() < 0.5 { return None; }
            let local = pt - c.origin;
            let v = local.dot(axis);
            let radial = local - axis * v;
            let u = radial.y.atan2(radial.x);
            let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
            Some(DVec2::new(u, v))
        }
        Surface3::Cone(c) => {
            let axis = c.axis_dir();
            let apex_to_pt = pt - c.apex;
            let v = apex_to_pt.dot(axis);
            let radial = apex_to_pt - axis * v;
            let u = radial.y.atan2(radial.x);
            let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
            Some(DVec2::new(u, v))
        }
        Surface3::Torus(t) => {
            let axis = t.axis.normalize_or_zero();
            if axis.length_squared() < 0.5 { return None; }
            let local = pt - t.center;
            let v = local.dot(axis);
            let radial = local - axis * v;
            let u = radial.y.atan2(radial.x);
            let tube_dir = radial.cross(axis).normalize_or_zero();
            let tube_local = local - radial;
            let w = tube_local.dot(tube_dir);
            // Simplified torus UV (OCCT uses analytic projection)
            let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
            Some(DVec2::new(u, w.atan2(t.minor_radius)))
        }
        // OCCT-aligned: numerical projection for BSpline/Bezier surfaces
        // (GeomAPI_ProjectPointOnSurf / Extrema_ExtPS).  Used by perform_areas
        // for hole-detection via UV boundary classification.
        Surface3::BSpline(_) | Surface3::Bezier(_) | Surface3::TriBezier(_) => {
            let proj = rcad_kernel::projection::closest_point_on_surface(surface, pt, 16);
            if proj.distance.is_finite() {
                Some(DVec2::new(proj.params.0, proj.params.1))
            } else {
                None
            }
        }
        _ => None,
    }
}

///  DS  ()
pub(crate) fn are_verts_coincident(ds: &DS, vi: usize, vj: usize) -> bool {
    if vi == vj { return true; }
    let d2 = ds.vertices[vi].point.distance_squared(ds.vertices[vj].point);
    d2 < TOLERANCE_ABS_SQ
}

// ================================================================
// OCCT-aligned: Angle2D  (BOPAlgo_WireSplitter_1.cxx L769-841)
// ================================================================

/// Convert a 2D direction vector to an angle in [0, 2).
///  OCCT  atan2(dir.y, dir.x)  [0, 2)
#[inline]
/// OCCT-aligned:  wire   BOPAlgo_WireSplitter
///    MakeConnexityBlocks + Path approach (PerformLoops L239-383)
///
/// Build closed wires:
///   1. MakeConnexityBlocks: BFS grouping by shared vertices
///   2. Regular block ( degree=2):
///   3. Irregular block ( degree>2 ): SmartMap + Path
// Returns (wires, internal_wires, vertex_positions) where vertex_positions
// maps canonical vertex indices (>= ds.vertices.len()) to their 3D position.
/// 鉁?OCCT-aligned: build canonical vertex map so different DS vertex indices
/// at the same 3D position map to one canonical index (OCCT BRep shares
/// TopoDS_Vertex).  Skips degenerate virtual-end vertices (>= ds.vertices.len()).
/// Extracted so the BuilderFace-level PerformShapesToAvoid and the WireSplitter
/// (build_closed_wires) agree on pole canonicalization.
pub(crate) fn build_vi_to_canon(segments: &[WireSegment], ds: &DS) -> Vec<usize> {
    let mut canon_vertices: Vec<DVec3> = Vec::new();
    let mut vi_to_canon: Vec<usize> = vec![usize::MAX; ds.vertices.len()];
    for seg in segments.iter() {
        if seg.end_vertex >= ds.vertices.len() { continue; } // skip deg (virtual end)
        for &vi in &[seg.start_vertex, seg.end_vertex] {
            if vi_to_canon[vi] != usize::MAX { continue; }
            let pt = ds.vertices[vi].point;
            let found = canon_vertices.iter().position(|c| c.distance_squared(pt) < TOLERANCE_ABS * TOLERANCE_ABS * 100_000_000.0);
            let canon = found.unwrap_or_else(|| { canon_vertices.push(pt); canon_vertices.len() - 1 });
            vi_to_canon[vi] = canon;
        }
    }
    vi_to_canon
}

/// 鉁?OCCT-aligned: physical-edge identity for a WireSegment.
/// Collapses FWD/REV of one physical edge (same source + same unordered
/// canonical endpoint pair) to ONE id, while keeping seam sub-edges that
/// share DsEdge(ei) but span different vertex pairs distinct.  This is the
/// rcad equivalent of TopoDS_Edge TShape identity used by BuilderFace's
/// aMVE (MapShapesAndAncestors VERTEX->EDGE).
pub(crate) fn physical_edge_id(seg: &WireSegment, vi_to_canon: &[usize], ds: &DS) -> (u8, usize, usize, usize) {
    let (tag, idx) = match &seg.source {
        WireEdgeSource::DsEdge(ei) => (0u8, *ei),
        WireEdgeSource::IntersectionCurve(ci) => (1u8, *ci),
        WireEdgeSource::SeamEdge => (2u8, 0),
    };
    let canon = |v: usize| vi_to_canon.get(v).copied().unwrap_or(v);
    let a = canon(seg.start_vertex);
    // Degenerate virtual end keeps its own id; never collapses with others.
    let b = if seg.end_vertex >= ds.vertices.len() { usize::MAX } else { canon(seg.end_vertex) };
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    (tag, idx, lo, hi)
}

/// 鉁?OCCT-aligned: WireSplitter / PerformLoops (BOPAlgo_WireSplitter).
///   OCCT BOPAlgo_WireSplitter organizes edges into ordered closed wires
///   by tracing 2D pcurves.  rcad: SmartMap-based edge-to-wire assembly
///   using canonical vertex indices and canonicalized edge connectivity.
pub(crate) fn build_closed_wires(segments: &mut Vec<WireSegment>, ds: &DS, face_idx: usize, avoided: &std::collections::HashSet<usize>) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, HashMap<usize, DVec3>) {
    if segments.is_empty() {
        return (vec![], vec![], HashMap::new());
    }

    let n = segments.len();

    // OCCT-aligned: canonicalize vertex indices so that different DS vertex
    // indices at the same 3D position map to a single canonical vertex.
    // OCCT BRep shares TopoDS_Vertex objects; rcad DS may assign different
    // indices to the same position (seam pole vs IC endpoint at pole).
    // 鉁?OCCT-aligned: vi_to_canon built skipping deg edges (build_vi_to_canon).
    let vi_to_canon: Vec<usize> = build_vi_to_canon(segments, ds);
    // Rebuild canon_vertices positions indexed by canonical id (deg_end_canon
    // below pushes new offset positions onto this vector).
    let mut canon_vertices: Vec<DVec3> = {
        let maxc = vi_to_canon.iter().filter(|&&c| c != usize::MAX).copied().max().map_or(0, |m| m + 1);
        let mut cv = vec![DVec3::ZERO; maxc];
        for vi in 0..ds.vertices.len() {
            let c = vi_to_canon[vi];
            if c != usize::MAX { cv[c] = ds.vertices[vi].point; }
        }
        cv
    };

    // Deg end canonical vertices with offset position, only for non-split seams.
    let seam_is_split = segments.iter().any(|s| {
        s.is_closed_on_face && matches!(&s.source, WireEdgeSource::DsEdge(ei)
            if ds.edges.get(*ei).map_or(0, |e| e.pave_blocks.len()) > 1)
    });
    let mut deg_end_canon: HashMap<usize, usize> = HashMap::new();
    if !seam_is_split {
        for (si, seg) in segments.iter().enumerate() {
            // OCCT-aligned: detect deg edges by virtual end vertex (>= ds.vertices.len())
            if seg.end_vertex >= ds.vertices.len() {
                let pt = ds.vertices[seg.start_vertex].point;
                canon_vertices.push(pt);
                deg_end_canon.insert(si, canon_vertices.len() - 1);
            }
        }
    }

    // OCCT-aligned: edge TShape dedup (BOPTools_AlgoTools.cxx L199-211).
    // IC section edges appear TWICE (FWD+REV) in the segment list.
    // The first appearance is the "primary" copy; the second is a duplicate.
    // Duplicate edges always make a block irregular.
    let mut seen_sources: HashSet<(u8, usize)> = HashSet::new();
    let mut duplicate_segs: HashSet<usize> = HashSet::new();
        for (si, seg) in segments.iter().enumerate() {
            if avoided.contains(&si) { continue; } // OCCT: avoided edges not in WireSplitter input
        // 鉁?OCCT-aligned: degenerate self-loop seam edges (sphere pole) appear twice in
        //   the WES (FORWARD+REVERSED) like any closed edge 鈥?not duplicates.  OCCT's
        //   bIsClosed guard (L148: !bIsClosed) preserves the second entry in aMS.
        if seg.is_closed_on_face && seg.start_vertex == seg.end_vertex { continue; }
        let variant = match &seg.source {
            WireEdgeSource::IntersectionCurve(ci) => (1u8, *ci),
            WireEdgeSource::DsEdge(ei) => (0u8, *ei),
            _ => continue,
        };
        if !seen_sources.insert(variant) {
            duplicate_segs.insert(si);
        }
    }

    // Build vertex鈫抯egments adjacency using CANONICAL vertex indices
    let mut vert_to_segs: HashMap<usize, Vec<usize>> = HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        if avoided.contains(&si) { continue; } // OCCT: avoided edges get no adjacency (not in aWES)
        let sv = vi_to_canon.get(seg.start_vertex).copied().unwrap_or(seg.start_vertex);
        let ev = deg_end_canon.get(&si).copied().unwrap_or_else(||
            vi_to_canon.get(seg.end_vertex).copied().unwrap_or(seg.end_vertex));
        vert_to_segs.entry(sv).or_default().push(si);
        vert_to_segs.entry(ev).or_default().push(si);
    }
    // 鉁?OCCT-aligned: DoSplitSEAMOnFace equivalent 鈥?reroute seam_rev
    //   through deg_end_canon vertices + remove redundant deg directions.
    //   This makes the seam+deg block regular (1 in + 1 out per vertex).
    let mut vertex_positions: HashMap<usize, DVec3> = HashMap::new();
    if deg_end_canon.len() == 2 {
        for &canon in deg_end_canon.values() {
            vertex_positions.insert(canon, canon_vertices[canon]);
        }
    }

    // MakeConnexityBlocks: BFS to find connected components
    let blocks = make_connexity_blocks(segments, &avoided, &vi_to_canon, &vert_to_segs, n);

    // Merge blocks that share canonical vertices (workaround for canonical
    // mapping precision issues that can split connected components).
    let mut merged_blocks: Vec<Vec<usize>> = Vec::new();
    {
        let n = blocks.len();
        let mut block_merged = vec![false; n];
        // Build vertex鈫抌lock index map using RAW vertex indices
        let mut v_to_b: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
        for (bi, b) in blocks.iter().enumerate() {
            for &si in b {
                let seg = &segments[si];
                v_to_b.entry(seg.start_vertex).or_default().push(bi);
                v_to_b.entry(seg.end_vertex).or_default().push(bi);
            }
        }
        for start_bi in 0..n {
            if block_merged[start_bi] { continue; }
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start_bi);
            block_merged[start_bi] = true;
            let mut merged: Vec<usize> = Vec::new();
            while let Some(bi) = queue.pop_front() {
                for &si in &blocks[bi] {
                    if !merged.contains(&si) { merged.push(si); }
                }
                // Find all blocks sharing ANY vertex with this block
                for &si in &blocks[bi] {
                    let seg = &segments[si];
                    for &vi in &[seg.start_vertex, seg.end_vertex] {
                        if let Some(neighbors) = v_to_b.get(&vi) {
                            for &nbi in neighbors {
                                if !block_merged[nbi] {
                                    block_merged[nbi] = true;
                                    queue.push_back(nbi);
                                }
                            }
                        }
                    }
                }
            }
            if !merged.is_empty() {
                merged_blocks.push(merged);
            }
        }
    }

    if std::env::var("RCAD_DEBUG_IC").is_ok() {
        eprintln!("[BLK_TRACE] fi={} n_merged_blocks={} n_total_segments={}", face_idx, merged_blocks.len(), segments.len());
        for (bi, b) in merged_blocks.iter().enumerate() {
            let seg_desc: Vec<String> = b.iter().map(|&si| {
                let seg = &segments[si];
                match &seg.source {
                    WireEdgeSource::DsEdge(ei) => format!("Ds({})", ei),
                    WireEdgeSource::IntersectionCurve(ci) => format!("IC({})", ci),
                    WireEdgeSource::SeamEdge => "Seam".to_string(),
                }
            }).collect();
            eprintln!("[BLK_TRACE]   block[{}] len={} segs=[{}]", bi, b.len(), seg_desc.join(","));
        }
    }

    // Process each block
    let mut wires: Vec<Vec<usize>> = Vec::new();
    let mut internal_wires: Vec<Vec<usize>> = Vec::new();

    for (bi, block) in merged_blocks.iter().enumerate() {
        if block.len() < 2 { continue; }
        if std::env::var("RCAD_DEBUG_IC").is_ok() {
            eprintln!("[BLK] fi={} bi={} n={}", face_idx, bi, block.len());
        }

        // 鉁?OCCT-aligned: aMS add/remove cycle (WireSplitter_1.cxx L149-152).
        //   OCCT: if (!aMS.Add(aE) && !bIsClosed) { aMS.Remove(aE); }
        //   Edges appearing exactly once 鈫?in aMS 鈫?IsInside=false (boundary).
        //   Edges appearing twice (FWD+REV) 鈫?removed from aMS 鈫?IsInside=true (interior).
        //   rcad: key by (source_type, canonical_idx) for identity.
        let mut a_ms: std::collections::HashSet<(u8, usize)> = std::collections::HashSet::new();
        for &si in block {
            let seg = &segments[si];
            let canon_key = match &seg.source {
                WireEdgeSource::DsEdge(ei) => (0u8, *ei),
                WireEdgeSource::IntersectionCurve(ci) => (1u8, *ci),
                WireEdgeSource::SeamEdge => continue,
            };
            let b_is_closed_check = seg.start_vertex == seg.end_vertex
                || (seg.is_closed_on_face && !matches!(&seg.source, WireEdgeSource::DsEdge(ei)
                    if ds.edges.get(*ei).map_or(false, |e| e.pave_blocks.len() > 1)));
            if !a_ms.insert(canon_key) && !b_is_closed_check {
                a_ms.remove(&canon_key);
            }
        }

        // 鉁?OCCT-aligned: Build SmartMap (WireSplitter_1.cxx L154-220).
        //    Always built first, used for BOTH regularity check and Path walk.
        let mut smart_map: IndexMap<usize, Vec<EdgeInfo>> = IndexMap::new();
        // OCCT L131: aVertMap 鈥?per-vertex closed edge flag (built inline)
        let mut vert_map: HashMap<usize, bool> = HashMap::new();
        for &si in block {
            let seg = &segments[si];
            // OCCT L142-144: skip edges without CurveOnSurface on the face
            let has_pcurve = match &seg.source {
                WireEdgeSource::IntersectionCurve(ci) => {
                    ds.intersection_curves.get(*ci).map_or(false, |ic| {
                        ic.pcurve_on_a.is_some() || ic.pcurve_on_b.is_some()
                    })
                }
                _ => seg.first_pcurve.is_some() || seg.second_pcurve.is_some(),
            };
            if !has_pcurve { continue; }

            // 鉁?OCCT-aligned: IsInside computed from aMS state (L310).
            //   Edges present in aMS (appeared once) 鈫?boundary (IsInside=false).
            //   Edges absent from aMS (appeared twice) 鈫?interior (IsInside=true).
            let canon_key = match &seg.source {
                WireEdgeSource::DsEdge(ei) => (0u8, *ei),
                WireEdgeSource::IntersectionCurve(ci) => (1u8, *ci),
                WireEdgeSource::SeamEdge => (2u8, 0),
            };
            let is_inside = !a_ms.contains(&canon_key);
            let is_circle_arc = is_inside && match &seg.source {
                WireEdgeSource::IntersectionCurve(ci) => {
                    ds.intersection_curves.get(*ci).map_or(false, |ic| {
                        matches!(&ic.curve, rcad_kernel::geom::Curve3::Circle(_))
                    })
                }
                _ => false,
            };

            // OCCT L147: bIsClosed = Degenerated(aE) || IsClosed(aE, myFace)
            let b_is_closed = seg.start_vertex == seg.end_vertex
                || (seg.is_closed_on_face && !matches!(&seg.source, WireEdgeSource::DsEdge(ei)
                    if ds.edges.get(*ei).map_or(false, |e| e.pave_blocks.len() > 1)));

            // OCCT L167-172: ONE EdgeInfo per vertex, in_flag from vertex orientation.
            // OCCT TopoDS_Edge: first vertex 鈫?FORWARD (out), second 鈫?REVERSED (in).
            // rcad: regardless of forward flag, start_vertex is geometrically where the
            // edge begins (OUT) and end_vertex is where it ends (IN).  Using forward flag
            // would create imbalanced IN/OUT for FWD+REV section edge pairs (v5 got 2 OUT,
            // v3 got 2 IN), breaking Path walking at irregular vertices.
            let sv_in = false;  // start_vertex 鈫?FORWARD/out
            smart_map.entry(seg.start_vertex).or_default().push(EdgeInfo {
                seg_idx: si, passed: false, in_flag: sv_in, is_inside, is_circle_arc, angle: 0.0,
            });
            let ev_in = true;   // end_vertex 鈫?REVERSED/in
            smart_map.entry(seg.end_vertex).or_default().push(EdgeInfo {
                seg_idx: si, passed: false, in_flag: ev_in, is_inside, is_circle_arc, angle: 0.0,
            });

            // OCCT L184-194: aVertMap 鈥?bind bIsClosed per vertex
            let e_sv = vert_map.entry(seg.start_vertex).or_default();
            *e_sv = *e_sv || b_is_closed;
            let e_ev = vert_map.entry(seg.end_vertex).or_default();
            *e_ev = *e_ev || b_is_closed;
        }

        // OCCT L298-319: compute angles for ALL EdgeInfo entries using Angle2D.
        // OCCT L316: aAngle = Angle2D(aVV, aE, myFace, aBAS, bIsIN, theContext)
        for (v, infos) in smart_map.iter_mut() {
            for ei in infos.iter_mut() {
                let seg = &segments[ei.seg_idx];
                // 鉁?OCCT-aligned: BRep_Tool::Parameter(aV, aE) 鈥?use vertex_params for DSEdges.
                let t_v = match &seg.source {
                    WireEdgeSource::DsEdge(ei) => {
                        ds.edges[*ei].vertex_param(*v).unwrap_or_else(|| {
                            if *v == seg.start_vertex { seg.t_range[0] } else { seg.t_range[1] }
                        })
                    }
                    _ => {
                        if *v == seg.start_vertex { seg.t_range[0] } else { seg.t_range[1] }
                    }
                };
                let domain = seg.t_range;
                let (curve, curve_domain) = match &seg.source {
                    WireEdgeSource::IntersectionCurve(ci) => {
                        let ic = &ds.intersection_curves[*ci];
                        if let Some(ref pc) = ic.pcurve_on_a {
                            let (ta, tb) = pc_parameter_range(pc);
                            (pc.clone(), [ta, tb])
                        } else if let Some(ref pc) = ic.pcurve_on_b {
                            let (ta, tb) = pc_parameter_range(pc);
                            (pc.clone(), [ta, tb])
                        } else { continue; }
                    }
                    _ => {
                        let pc = seg.first_pcurve.as_ref().or(seg.second_pcurve.as_ref());
                        match pc {
                            Some(pc) => (pc.clone(), domain),
                            None => continue,
                        }
                    }
                };
                ei.angle = angle_2d(&curve, t_v, curve_domain, ei.in_flag, &ds.faces[face_idx].surface, ds.vertices[*v].geom_tol, None)
                    .unwrap_or(0.0);
                if std::env::var("RCAD_DEBUG_SM").is_ok() && face_idx == 0 {
                    let seg = &segments[ei.seg_idx];
                    eprintln!("[SM] fi=0 v={} seg={} src={:?} pc={} ang={:.3} in={}", v, ei.seg_idx,
                        std::mem::discriminant(&seg.source),
                        seg.first_pcurve.is_some() || seg.second_pcurve.is_some(),
                        ei.angle, ei.in_flag);
                }
            }
        }

        // 鉁?OCCT-aligned: regularity check (L222-280) from SmartMap IN/OUT.
        //   Step 1 (L222-260): each vertex 1 IN + 1 OUT. Step 2 (L261-280):
        //   no duplicate edges.
        let mut is_regular = !block.iter().any(|&si| duplicate_segs.contains(&si));
        if is_regular {
            for (_, infos) in &smart_map {
                let in_cnt = infos.iter().filter(|ei| ei.in_flag).count();
                let out_cnt = infos.iter().filter(|ei| !ei.in_flag).count();
                if in_cnt != 1 || out_cnt != 1 {
                    is_regular = false;
                    break;
                }
            }
        }

        if is_regular {
            // OCCT L282-290: MakeWire 鈥?extract simple wire (no angles needed).
            if let Some(wire) = build_regular_wire(block, segments, &vert_to_segs, &vi_to_canon, &deg_end_canon) {
                wires.push(wire);
            }
        } else {
            // OCCT L292-358: SplitBlock 鈥?refine angles + path walk for irregular blocks.
            split_block(block, segments, &mut smart_map, ds, face_idx, &mut wires);
            if std::env::var("RCAD_DEBUG_IC").is_ok() {
                eprintln!("[BLK_WIRES] fi={} bi={} n_wires_in_block={}", face_idx, bi,
                    wires.iter().filter(|w| w.iter().any(|&si| block.contains(&si))).count());
            }
        }
    }

    (wires, internal_wires, vertex_positions)
}

/// 鉁?OCCT-aligned: BOPTools_AlgoTools::MakeConnexityBlocks(start_elements, VERTEX, EDGE).
///   Groups segments into connected components (blocks) by shared vertices.
///   Each block is a subset of segments that form a connected sub-graph.
pub(crate) fn make_connexity_blocks(
    segments: &[WireSegment],
    avoided: &std::collections::HashSet<usize>,
    vi_to_canon: &[usize],
    vert_to_segs: &HashMap<usize, Vec<usize>>,
    n: usize,
) -> Vec<Vec<usize>> {
    let mut visited_seg = vec![false; n];
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for si in 0..n {
        if visited_seg[si] { continue; }
        if avoided.contains(&si) { continue; }
        let mut block = Vec::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(si);
        visited_seg[si] = true;
        while let Some(ci) = queue.pop_front() {
            block.push(ci);
            let seg = &segments[ci];
            for &vi in &[seg.start_vertex, seg.end_vertex] {
                let cvi = vi_to_canon.get(vi).copied().unwrap_or(vi);
                if let Some(neighbors) = vert_to_segs.get(&cvi) {
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

/// 鉁?OCCT-aligned: SplitBlock (BOPAlgo_WireSplitter.cxx L292-358).
///   Handles irregular blocks (vertices with degree > 2) by refining edge
///   angles at multi-connected vertices, then walking paths through the
///   SmartMap to extract closed wires.  OCCT parallelizes per block;
///   rcad runs sequentially with the same angle + path logic.
pub(crate) fn split_block(
    block: &[usize],
    segments: &[WireSegment],
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    ds: &DS,
    face_idx: usize,
    wires: &mut Vec<Vec<usize>>,
) {
    // OCCT L327: RefineAngles 鈥?compute turning angles at each vertex.
    refine_angles(smart_map, segments, ds, face_idx);
    dbg_smartmap!("split_block", face_idx, smart_map);

    // OCCT L331-358: Path walk 鈥?iterate all vertices by insertion order,
    //   for each unpassed OUT entry start a new Path walk.
    let order_keys: Vec<usize> = smart_map.keys().copied().collect();
    for &v in &order_keys {
        let Some(infos) = smart_map.get(&v).cloned() else { continue; };
        for ei in &infos {
            if !ei.passed && !ei.in_flag
                && ei.seg_idx < segments.len()
                && (segments[ei.seg_idx].start_vertex != segments[ei.seg_idx].end_vertex
                    || segments[ei.seg_idx].is_closed_on_face)
            {
                walk_path_extract_wires(ei.seg_idx, segments, smart_map, wires, ds, face_idx);
            }
        }
    }
}

/// 鉁?OCCT-aligned: Regular block (degree=2) wire build.
pub(crate) fn build_regular_wire(
    block: &[usize],
    segments: &[WireSegment],
    vert_to_segs: &HashMap<usize, Vec<usize>>,
    vi_to_canon: &[usize],
    _deg_end_canon: &HashMap<usize, usize>,
) -> Option<Vec<usize>> {
    let cs = |seg: &WireSegment| vi_to_canon.get(seg.start_vertex).copied().unwrap_or(seg.start_vertex);
    let ce = |seg: &WireSegment| {
        // deg_end_canon is for specific seg indices; we don't have si here, use vi_to_canon
        vi_to_canon.get(seg.end_vertex).copied().unwrap_or(seg.end_vertex)
    };
    let block_set: std::collections::HashSet<usize> = block.iter().copied().collect();
    let mut visited = vec![false; segments.len()];
    let mut wire: Vec<usize> = Vec::new();

    let start_si = block[0];
    let start_seg = &segments[start_si];
    let start_vertex = cs(start_seg);
    let mut ci = start_si;
    let mut arrived_vertex = ce(start_seg);

    loop {
        visited[ci] = true;
        wire.push(ci);
        if arrived_vertex == start_vertex && wire.len() >= 2 { break; }

        let next = vert_to_segs.get(&arrived_vertex).and_then(|neighbors| {
            neighbors.iter().find(|&&ni| !visited[ni] && block_set.contains(&ni))
        }).copied();

        match next {
            Some(ni) => {
                let seg = &segments[ni];
                ci = ni;
                arrived_vertex = if cs(seg) == arrived_vertex { ce(seg) } else { cs(seg) };
            }
            None => break,
        }
    }

    if wire.len() >= 2 { Some(wire) } else { None }
}

/// OCCT-aligned: EdgeInfo  (BOPAlgo_WireSplitter.lxx L22-69)
#[derive(Debug, Clone)]
pub(crate) struct EdgeInfo {
    pub(crate) seg_idx: usize,
    pub(crate) passed: bool,
    pub(crate) in_flag: bool,
    pub(crate) is_inside: bool,
    pub(crate) is_circle_arc: bool,
    pub(crate) angle: f64,
}

// (SmartMap + Path moved into build_closed_wires 鈥?OCCT L154-358)

// ====================================================================
// 鉁?OCCT-aligned: PerformShapesToAvoid (BOPAlgo_BuilderFace.cxx L152-235)
//
// Face-level (BuilderFace) pass over the whole segment set (= myShapes).
// Builds the vertex->edge ancestor map (aMVE) using physical-edge identity
// so FWD/REV of one edge count as ONE edge, then repeatedly avoids:
//   - aNbE==1 dangling edges (non-degenerate)          (OCCT L198-210)
//   - aNbE==2 && aE2.IsSame(aE1) self-coincident edges  (OCCT L211-227)
// Returns the set of avoided SEGMENT indices (both FWD+REV of each avoided
// physical edge).  The caller excludes these from the WireSplitter input.
// ====================================================================
/// 鉁?OCCT-aligned: PerformShapesToAvoid (BOPAlgo_BuilderFace.cxx L152-235).
/// Physical edge identifier 鈥?groups WireSegments that belong to the same
/// source DS edge (or intersection curve).  OCCT-aligned: TopoDS_Edge identity.
pub(crate) type Pid = (u8, usize, usize, usize);

/// Build PID鈫抯egment and PID鈫抏ndpoint maps from WireSegments.
/// OCCT: segments grouped by their TopoDS_Edge identity.
pub(crate) fn build_pid_maps(
    segments: &[WireSegment],
    vi_to_canon: &[usize],
    ds: &DS,
) -> (
    std::collections::HashMap<Pid, Vec<usize>>,       // pid 鈫?segment indices
    std::collections::HashMap<Pid, (usize, usize)>,   // pid 鈫?(canon_start, canon_end)
) {
    let mut pid_segs: std::collections::HashMap<Pid, Vec<usize>> = std::collections::HashMap::new();
    let mut pid_endpoints: std::collections::HashMap<Pid, (usize, usize)> = std::collections::HashMap::new();
    let canon = |v: usize| -> usize {
        if v >= ds.vertices.len() { usize::MAX } else { vi_to_canon.get(v).copied().unwrap_or(v) }
    };
    for (si, seg) in segments.iter().enumerate() {
        let pid = physical_edge_id(seg, vi_to_canon, ds);
        pid_segs.entry(pid).or_default().push(si);
        pid_endpoints.entry(pid).or_insert_with(|| (canon(seg.start_vertex), canon(seg.end_vertex)));
    }
    (pid_segs, pid_endpoints)
}

/// Expand avoided PIDs to segment indices.  Called by split_face_and_emit_topo_ds
/// (OCCT: PerformLoops checks myShapesToAvoid.Contains, rcad needs segment indices).
pub(crate) fn expand_avoided_pids(
    avoided_pids: &std::collections::HashSet<Pid>,
    pid_segs: &std::collections::HashMap<Pid, Vec<usize>>,
) -> std::collections::HashSet<usize> {
    let mut segs = std::collections::HashSet::new();
    for pid in avoided_pids {
        if let Some(slist) = pid_segs.get(pid) {
            for &si in slist { segs.insert(si); }
        }
    }
    segs
}

/// 鉁?OCCT-aligned: PerformShapesToAvoid (BOPAlgo_BuilderFace.cxx L152-235).
/// Returns the set of physical edge PIDs to avoid (OCCT: myShapesToAvoid)
/// AND the pid鈫抯egment map for the caller to expand PIDs to segment indices.
pub(crate) fn perform_shapes_to_avoid(
    segments: &[WireSegment],
    vi_to_canon: &[usize],
    ds: &DS,
) -> (std::collections::HashSet<Pid>, std::collections::HashMap<Pid, Vec<usize>>) {
    let (pid_segs, pid_endpoints) = build_pid_maps(segments, vi_to_canon, ds);

    let is_degenerate = |pid: &Pid| -> bool {
        // Degenerate edge: virtual end vertex (b == usize::MAX in physical_edge_id).
        pid.3 == usize::MAX
    };

    let mut avoided_pids: std::collections::HashSet<Pid> = std::collections::HashSet::new();
    loop {
        let mut b_found = false;
        // Build ancestor map aMVE: vertex -> list of incident physical edge ids
        // (excluding already-avoided edges). A closed edge (a==b) is pushed twice.
        let mut anc: std::collections::HashMap<usize, Vec<Pid>> = std::collections::HashMap::new();
        for (&pid, &(a, b)) in &pid_endpoints {
            if avoided_pids.contains(&pid) { continue; }
            if a != usize::MAX { anc.entry(a).or_default().push(pid); }
            if b != usize::MAX { anc.entry(b).or_default().push(pid); }
        }
        for (&canon_vi, ids) in &anc {
            let a_nb_e = ids.len();
            // OCCT L204-207: INTERNAL vertices 鈫?skip (vertex inside face,
            //   not on boundary 鈥?dangling edges at internal vertices are valid).
            if canon_vi < ds.vertices.len() && ds.vertices[canon_vi].is_internal {
                continue;
            }
            if a_nb_e == 1 {
                // OCCT L198-210: dangling edge 鈫?avoid (skip degenerate).
                let pid = ids[0];
                if is_degenerate(&pid) { continue; }
                if avoided_pids.insert(pid) { b_found = true; }
            } else if a_nb_e == 2 && ids[0] == ids[1] {
                // OCCT L211-227: same edge twice at this vertex (self-coincident).
                let pid = ids[0];
                let (a, b) = pid_endpoints[&pid];
                if a == b { continue; } // OCCT L219-222: self-loop (closed) 鈫?keep
                if avoided_pids.insert(pid) { b_found = true; }
            }
        }
        if !b_found { break; }
    }

    (avoided_pids, pid_segs)
}

// ====================================================================
// OCCT-aligned: Assemble internal wires from avoided segments
// (BOPAlgo_BuilderFace.cxx L327-382)
// ====================================================================
/// 鉁?OCCT-aligned: PerformInternalShapes (BOPAlgo_BuilderFace.cxx L327-382).
/// 鉁?OCCT-aligned: BuilderFace::PerformInternalShapes (L618-735).
///   Classify avoided (internal) edges against each result WireFace,
///   assemble edges that fall INSIDE the face into per-face internal wires.
///
/// OCCT flow:
///   L642-663: Build BVH tree of 2D UV boxes for each edge
///   L674-716: For each result face, use BVH + IsInside 鈫?select internal edges
///   L718-735: MakeInternalWires (vertex-degree-based wire assembly) + add to face
///
/// rcad: for each WireFace, build 2D outer boundary polygon from segment pcurves,
///   classify each avoided segment's UV midpoint via 2D ray casting (point-in-polygon).
///   Segments inside the outer boundary (but not inside a hole) 鈫?assemble into
///   internal wires for that face.  Returns per-face internal wire segment groups:
///   `Vec<Vec<Vec<usize>>>` 鈥?outer index = WireFace index, inner = internal wires
///   for that face, each wire = Vec of segment indices.
pub(crate) fn assemble_internal_wires(
    avoided: &[usize],
    segments: &[WireSegment],
    wfs: &[WireFace],
) -> Vec<Vec<Vec<usize>>> {
    if avoided.is_empty() || wfs.is_empty() {
        return vec![vec![]; wfs.len()];
    }

    // OCCT L633-663: Build BVH tree 鈥?Bnd_Box2d per edge from pcurve UV bounds.
    //   OCCT: BRepTools::AddUVBounds(myFace, aE, aBoxE) samples the edge's
    //   pcurve on the face surface.  rcad: compute UV bounding box from segment's
    //   first_pcurve at sampled points (matching OCCT sampling density).
    let seg_uv_box: Vec<Option<[f64; 4]>> = avoided.iter().map(|&si| {
        let seg = &segments[si];
        seg.first_pcurve.as_ref().map(|pc| {
            let [t0, t1] = seg.t_range;
            let n_pts = 8usize; // OCCT IntTools_FClass2d uses NbSamples (鈮?)
            let mut u_min = f64::INFINITY; let mut u_max = f64::NEG_INFINITY;
            let mut v_min = f64::INFINITY; let mut v_max = f64::NEG_INFINITY;
            for k in 0..n_pts {
                let t = t0 + (t1 - t0) * k as f64 / (n_pts - 1) as f64;
                let uv = pc.point_at(t);
                u_min = u_min.min(uv.x); u_max = u_max.max(uv.x);
                v_min = v_min.min(uv.y); v_max = v_max.max(uv.y);
            }
            [u_min, u_max, v_min, v_max]
        })
    }).collect();

    // OCCT L674-716: For each WireFace, classify avoided segments via IsInside.
    let mut face_internal: Vec<Vec<usize>> = vec![Vec::new(); wfs.len()];

    for (fi, wf) in wfs.iter().enumerate() {
        // Build 2D outer boundary polygon from outer wire segments' pcurves.
        let outer_uv: Vec<DVec2> = wf.outer_wire.iter().filter_map(|&si| {
            if si >= segments.len() { return None; }
            let seg = &segments[si];
            seg.first_pcurve.as_ref().map(|pc| pc.point_at(seg.t_range[0]))
        }).collect();
        if outer_uv.len() < 3 { continue; }

        // Build 2D hole polygons for inner wires (to exclude segments inside holes).
        let hole_uvs: Vec<Vec<DVec2>> = wf.inner_wires.iter().map(|iw| {
            iw.iter().filter_map(|&si| {
                if si >= segments.len() { return None; }
                let seg = &segments[si];
                seg.first_pcurve.as_ref().map(|pc| pc.point_at(seg.t_range[0]))
            }).collect()
        }).filter(|poly: &Vec<DVec2>| poly.len() >= 3).collect();

        // OCCT L704-716: select edges inside this face via 2D ray casting.
        //   OCCT uses BVH box prefilter + IsInside.  rcad: box overlap prefilter
        //   + 2D point-in-polygon (equivalent to IsInside).
        for (ai, &si) in avoided.iter().enumerate() {
            let Some(box_e) = &seg_uv_box[ai] else { continue; };
            // OCCT L694-695: BVH box prefilter 鈥?skip if box doesn't overlap face.
            if outer_uv.iter().all(|p| p.x < box_e[0] || p.x > box_e[1] || p.y < box_e[2] || p.y > box_e[3]) {
                continue;
            }
            let uv_mid = DVec2::new(0.5 * (box_e[0] + box_e[1]), 0.5 * (box_e[2] + box_e[3]));
            if !point_in_polygon_2d(&outer_uv, uv_mid) { continue; }
            let in_hole = hole_uvs.iter().any(|hole| point_in_polygon_2d(hole, uv_mid));
            if in_hole { continue; }
            face_internal[fi].push(si);
        }
    }

    // OCCT L724-725: MakeInternalWires 鈥?per-face BFS assembly.
    let mut per_face_wires: Vec<Vec<Vec<usize>>> = vec![Vec::new(); wfs.len()];
    for (fi, assigned) in face_internal.iter().enumerate() {
        if assigned.is_empty() { continue; }
        let mut v_to_segs: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for &si in assigned {
            let seg = &segments[si];
            v_to_segs.entry(seg.start_vertex).or_default().push(si);
            v_to_segs.entry(seg.end_vertex).or_default().push(si);
        }
        let mut added = vec![false; segments.len()];
        for &start_si in assigned {
            if added[start_si] { continue; }
            let mut wire: Vec<usize> = Vec::new();
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start_si);
            added[start_si] = true;
            while let Some(si) = queue.pop_front() {
                wire.push(si);
                let seg = &segments[si];
                for &vtx in &[seg.start_vertex, seg.end_vertex] {
                    if let Some(neighbors) = v_to_segs.get(&vtx) {
                        for &ni in neighbors {
                            if !added[ni] {
                                added[ni] = true;
                                queue.push_back(ni);
                            }
                        }
                    }
                }
            }
            if !wire.is_empty() {
                per_face_wires[fi].push(wire);
            }
        }
    }
    per_face_wires
}

pub(crate) fn is_same_block_fwd_rev(a: &WireSegment, b: &WireSegment) -> bool {
    match (&a.source, &b.source) {
        (WireEdgeSource::DsEdge(ea), WireEdgeSource::DsEdge(eb)) => {
            ea == eb
            && a.start_vertex == b.end_vertex
            && a.end_vertex == b.start_vertex
        }
        // 鉁?OCCT-aligned: IntersectionCurve FWD+REV share curve index
        //    (TopoDS_Shape::IsSame check, WireSplitter_1.cxx L564-567).
        (WireEdgeSource::IntersectionCurve(ca), WireEdgeSource::IntersectionCurve(cb)) => {
            ca == cb
        }
        // 鉁?OCCT-aligned: SeamEdge FWD+REV (same seam, opposite directions).
        (WireEdgeSource::SeamEdge, WireEdgeSource::SeamEdge) => {
            a.is_closed_on_face && b.is_closed_on_face && a.orientation != b.orientation
        }
        _ => false,
    }
}

/// Check if a segment has been marked passed at a specific vertex with a specific in_flag.
pub(crate) fn is_seg_passed(smart_map: &IndexMap<usize, Vec<EdgeInfo>>, seg_idx: usize) -> bool {
    for infos in smart_map.values() {
        if infos.iter().any(|ei| ei.seg_idx == seg_idx && ei.passed) {
            return true;
        }
    }
    false
}

/// Mark the specific EdgeInfo AND its opposite-direction counterpart
/// (same physical edge, opposite in_flag) at the given vertex as passed.
/// OCCT has 1 entry per edge per vertex; rcad creates 2 (FWD+REV) that
/// must be treated as one physical edge.
pub(crate) fn mark_edge_passed_both_dirs(
    smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>,
    seg_idx: usize,
    vertex: usize,
    in_flag: bool,
    segments: &[WireSegment],
) {
    let Some(infos) = smart_map.get_mut(&vertex) else { return };
    let physical_key = match &segments[seg_idx].source {
        WireEdgeSource::DsEdge(ei) => (*ei, true),
        WireEdgeSource::IntersectionCurve(ci) => (*ci, false),
        WireEdgeSource::SeamEdge => return,
    };
    for info in infos.iter_mut() {
        let matches_physical = match (&segments[info.seg_idx].source, physical_key) {
            (WireEdgeSource::DsEdge(ei), (pe, true)) => *ei == pe,
            (WireEdgeSource::IntersectionCurve(ci), (pc, false)) => *ci == pc,
            _ => false,
        };
        if matches_physical {
            info.passed = true;
        }
    }
}

/// Mark only the specific EdgeInfo for a segment at a vertex+in_flag as passed.
pub(crate) fn mark_edge_passed(smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>, seg_idx: usize, vertex: usize, in_flag: bool) {
    if let Some(infos) = smart_map.get_mut(&vertex) {
        for info in infos.iter_mut() {
            if info.seg_idx == seg_idx && info.in_flag == in_flag {
                info.passed = true;
                return;
            }
        }
    }
}

/// Mark both orientations of a segment as passed (used for initial cleanup).
/// Not used during Path walking  use mark_edge_passed instead.
#[allow(dead_code)]
pub(crate) fn mark_seg_passed(smart_map: &mut IndexMap<usize, Vec<EdgeInfo>>, seg_idx: usize) {
    for infos in smart_map.values_mut() {
        for info in infos.iter_mut() {
            if info.seg_idx == seg_idx {
                info.passed = true;
            }
        }
    }
}

/// Find the EdgeInfo angle for a segment at a vertex with the given in_flag.
pub(crate) fn find_angle_at(smart_map: &IndexMap<usize, Vec<EdgeInfo>>, seg_idx: usize, vertex: usize, in_flag: bool) -> Option<f64> {
    smart_map.get(&vertex)?.iter()
        .find(|ei| ei.seg_idx == seg_idx && ei.in_flag == in_flag)
        .map(|ei| ei.angle)
}

/// Select the best outgoing edge at a vertex using ClockWiseAngle minimum selection.
/// (OCCT L622-660)
pub(crate) fn select_best_outgoing<'a>(
    candidates: &[&'a EdgeInfo],
    angle_in: f64,
    incoming_is_boundary: bool,
    incoming_ci: usize,
) -> Option<&'a EdgeInfo> {
    if candidates.is_empty() {
        return None;
    }
    let a_two_pi = std::f64::consts::TAU;
    let eps = std::f64::EPSILON; // OCCT: eps = Epsilon(1.)
    let mut a_min_angle = 100.0;
    let mut a_nb_ways_inside: i32 = 0;
    let mut p_only_way_in: Option<&EdgeInfo> = None;
    let mut p_edge_info: Option<&EdgeInfo> = None;
    for an_ei in candidates {
        let a_angle = if an_ei.seg_idx == incoming_ci {
            a_two_pi // OCCT L564-567: aE.IsSame(aEOuta) -> aTwoPI
        } else {
            clock_wise_angle(angle_in, an_ei.angle) // OCCT L585-586
        };
        if incoming_is_boundary && an_ei.is_inside {
            a_nb_ways_inside += 1; // OCCT L589-593
            p_only_way_in = Some(an_ei);
        }
        if a_angle < a_min_angle - eps {
            a_min_angle = a_angle; // OCCT L595-599
            p_edge_info = Some(an_ei);
        }
    }
    if a_nb_ways_inside == 1 {
        p_edge_info = p_only_way_in; // OCCT L602-604
    }
    p_edge_info
}

// ---------------------------------------------------------------------------
// TopoDS-based perform_shapes_to_avoid 鈥?BRepTool variant
// ---------------------------------------------------------------------------

/// Physical edge identity for WireSegmentTopoDS (BRepTool variant).
pub(crate) fn physical_edge_id_topo_ds(seg: &super::types::WireSegmentTopoDS) -> (u8, usize, usize, usize) {
    let (tag, idx) = match &seg.source {
        super::types::WireEdgeSourceTopoDS::DsEdge(e) => (0u8, e.index),
        super::types::WireEdgeSourceTopoDS::IntersectionCurve(c) => (1u8, c.index),
        super::types::WireEdgeSourceTopoDS::SeamEdge => (2u8, 0),
    };
    let a = seg.start_vertex.index;
    let b = seg.end_vertex.index;
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    (tag, idx, lo, hi)
}

/// Build PID maps for WireSegmentTopoDS (BRepTool variant).
pub(crate) fn build_pid_maps_topo_ds(
    segments: &[super::types::WireSegmentTopoDS],
) -> (
    std::collections::HashMap<Pid, Vec<usize>>,
    std::collections::HashMap<Pid, (usize, usize)>,
) {
    let mut pid_segs: std::collections::HashMap<Pid, Vec<usize>> = std::collections::HashMap::new();
    let mut pid_endpoints: std::collections::HashMap<Pid, (usize, usize)> = std::collections::HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        let pid = physical_edge_id_topo_ds(seg);
        pid_segs.entry(pid).or_default().push(si);
        pid_endpoints.entry(pid).or_insert_with(|| (seg.start_vertex.index, seg.end_vertex.index));
    }
    (pid_segs, pid_endpoints)
}

/// 鉁?OCCT-aligned: PerformShapesToAvoid for WireSegmentTopoDS.
///   Builds vertex鈫扨id adjacency map from segments (OCCT: TopExp::MapShapesAndAncestors),
///   then removes dangling edges (valence 1) and self-coincident edges in a fixed-point loop.
///   Returns avoided PIDs and pid鈫抯egment map.
pub(crate) fn perform_shapes_to_avoid_topo(
    segments: &[super::types::WireSegmentTopoDS],
    tool: &dyn rcad_kernel::topods::BRepTool,
) -> (std::collections::HashSet<Pid>, std::collections::HashMap<Pid, Vec<usize>>) {
    let (pid_segs, _) = build_pid_maps_topo_ds(segments);

    let is_degenerate = |pid: &Pid| -> bool {
        let (_tag, idx, a, b) = *pid;
        if a == b { return true; }
        if _tag == 0 {
            tool.is_edge_degenerated(rcad_kernel::topods::ShapeRef::new(idx))
        } else {
            false
        }
    };

    let is_internal_vertex = |vi: usize| -> bool {
        tool.vertex_orientation(rcad_kernel::topods::ShapeRef::new(vi))
            == rcad_kernel::topods::Orientation::Internal
    };

    // OCCT L152-156: aMVE = MapShapesAndAncestors(VERTEX, EDGE).
    //   OCCT uses vertex+orientation as key: an edge's start vertex has
    //   FORWARD orientation, end vertex has REVERSED.  FWD and REV of
    //   the same edge at the same vertex → different orientation →
    //   separate entries → not seen as self-coincident.
    //   rcad: track (Pid, is_outgoing) per vertex.  is_outgoing == true
    //   for start_vertex matches OCCT's FORWARD; false for end_vertex
    //   matches REVERSED.
    let mut a_mve: std::collections::HashMap<usize, Vec<(Pid, bool)>> = std::collections::HashMap::new();
    for (si, seg) in segments.iter().enumerate() {
        let pid = physical_edge_id_topo_ds(seg);
        a_mve.entry(seg.start_vertex.index).or_default().push((pid, true));
        a_mve.entry(seg.end_vertex.index).or_default().push((pid, false));
    }

    let mut avoided_pids: std::collections::HashSet<Pid> = std::collections::HashSet::new();

    // OCCT L182-228: fixed-point loop → avoid dangling edges (valence 1) and
    // self-coincident edges, updating aMVE after each avoidance.
    loop {
        let mut b_found = false;
        for (&v, entries) in &a_mve {
            // OCCT: MapShapesAndAncestors counts each edge ONCE per vertex
            // (with the vertex's orientation on the edge).  rcad may have
            // both FWD and REV of the same edge at the same vertex, which
            // share a single Pid.  Collapse by Pid for OCCT-equivalent
            // valence counting.  Skip degenerate Pids (self-loop edges) as
            // OCCT does (BRep_Tool::Degenerated -> continue).
            let mut unique: Vec<&Pid> = entries.iter().map(|(p, _)| p).filter(|p| !is_degenerate(p)).collect();
            unique.sort_unstable();
            unique.dedup();
            let val = unique.len();

            // OCCT L198-210: dangling edge (valence 1 at vertex) -> avoid
            if val == 1 {
                let pid = *unique[0];
                if is_internal_vertex(v) { continue; }
                if avoided_pids.insert(pid) { b_found = true; }
                continue;
            }
            // OCCT L211-227: self-coincident -- aE2.IsSame(aE1) matches
            // TShape regardless of orientation.  Two entries with the same
            // Pid at the same vertex (even with different directions) means
            // the same physical edge, caught here when entries.len() == 2
            // (only those 2 entries exist at this vertex).
            if val == 2 && entries.len() == 2
                && entries[0].0 == entries[1].0
            {
                let pid = entries[0].0;
                if avoided_pids.insert(pid) { b_found = true; }
            }
        }
        if !b_found { break; }
        // OCCT L230: rebuild aMVE without avoided edges
        a_mve.clear();
        for (si, seg) in segments.iter().enumerate() {
            let pid = physical_edge_id_topo_ds(seg);
            if avoided_pids.contains(&pid) { continue; }
            a_mve.entry(seg.start_vertex.index).or_default().push((pid, true));
            a_mve.entry(seg.end_vertex.index).or_default().push((pid, false));
        }
    }

    (avoided_pids, pid_segs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bopds::ds::DS;
    use rcad_kernel::geom::Plane;
    use rcad_modeling::make_box_brep;
    use std::collections::{HashMap, HashSet};

    #[test]
    fn are_verts_coincident_same_index() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        if ds.vertices.len() < 16 { return; }
        assert!(are_verts_coincident(&ds, 0, 0));
    }

    #[test]
    fn are_verts_coincident_distant() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(10.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        if ds.vertices.len() < 16 { return; }
        assert!(!are_verts_coincident(&ds, 0, ds.a_vertex_count));
    }

    #[test]
    fn world_to_uv_plane() {
        let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
        let surface = Surface3::Plane(plane);
        let uv = world_to_uv(&surface, DVec3::new(2.0, 3.0, 0.0));
        assert!(uv.is_some());
    }

    #[test]
    fn world_to_uv_plane_origin() {
        let plane = Plane { origin: DVec3::new(1.0, 2.0, 3.0), normal: DVec3::Z };
        let surface = Surface3::Plane(plane);
        let uv = world_to_uv(&surface, DVec3::new(1.0, 2.0, 3.0));
        assert!(uv.is_some());
        let uv = uv.unwrap();
        assert!(uv.length_squared() < 1e-20, "origin should map to (0,0)");
    }

    #[test]
    fn is_same_block_fwd_rev_true_for_same_edge() {
        let seg_a = WireSegment {
            start_vertex: 0, end_vertex: 1,
            source: WireEdgeSource::DsEdge(5),
            orientation: WireOrientation::Forward, is_closed_on_face: false,
            second_pcurve: None, first_pcurve: None,
            t_range: [0.0, 1.0],
        };
        let seg_b = WireSegment {
            start_vertex: 1, end_vertex: 0,
            source: WireEdgeSource::DsEdge(5),
            orientation: WireOrientation::Reversed, is_closed_on_face: false,
            second_pcurve: None, first_pcurve: None,
            t_range: [1.0, 0.0],
        };
        assert!(is_same_block_fwd_rev(&seg_a, &seg_b));
    }

    #[test]
    fn is_same_block_fwd_rev_false_for_different_edge() {
        let seg_a = WireSegment {
            start_vertex: 0, end_vertex: 1,
            source: WireEdgeSource::DsEdge(5),
            orientation: WireOrientation::Forward, is_closed_on_face: false,
            second_pcurve: None, first_pcurve: None,
            t_range: [0.0, 1.0],
        };
        let seg_b = WireSegment {
            start_vertex: 1, end_vertex: 0,
            source: WireEdgeSource::DsEdge(7),
            orientation: WireOrientation::Reversed, is_closed_on_face: false,
            second_pcurve: None, first_pcurve: None,
            t_range: [1.0, 0.0],
        };
        assert!(!is_same_block_fwd_rev(&seg_a, &seg_b));
    }

    #[test]
    fn build_vi_to_canon_empty() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        let canon = build_vi_to_canon(&[], &ds);
        // Returns vec of usize::MAX for all DS vertices, not empty
        assert!(canon.iter().all(|&c| c == usize::MAX),
            "all canonical indices should be unset for empty segments");
    }

    #[test]
    fn build_closed_wires_empty_input() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        let avoided = HashSet::new();
        let mut segments = Vec::new();
        let (outer, inner, vpos) = build_closed_wires(&mut segments, &ds, 0, &avoided);
        assert!(outer.is_empty());
        assert!(inner.is_empty());
        assert!(vpos.is_empty());
    }

    #[test]
    fn mark_seg_passed_then_is_passed() {
        let mut sm: IndexMap<usize, Vec<EdgeInfo>> = IndexMap::new();
        sm.entry(0).or_default().push(EdgeInfo {
            seg_idx: 5, passed: false, in_flag: true,
            is_inside: false, is_circle_arc: false, angle: 1.0,
        });
        assert!(!is_seg_passed(&sm, 5));
        mark_seg_passed(&mut sm, 5);
        assert!(is_seg_passed(&sm, 5));
    }

    #[test]
    fn edge_info_debug_display() {
        let ei = EdgeInfo {
            seg_idx: 5, passed: false, in_flag: true,
            is_inside: false, is_circle_arc: false, angle: 0.0,
        };
        let _ = format!("{:?}", ei);
    }

    #[test]
    fn edge_angle_2d_none_for_zero_range() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a = edge_angle_2d(&line, 0.5, [0.0, 0.0], &p, false, 1e-5);
        assert!(a.is_none());
    }

    #[test]
    fn edge_angle_2d_line_on_plane() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a = edge_angle_2d(&line, 0.5, [0.0, 1.0], &p, false, 1e-5);
        assert!(a.is_some());
    }

    #[test]
    fn edge_angle_2d_is_in_flips() {
        let line = Curve3::Line(Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let p = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        let a_out = edge_angle_2d(&line, 0.9, [0.0, 1.0], &p, false, 1e-5).unwrap();
        let a_in = edge_angle_2d(&line, 0.9, [0.0, 1.0], &p, true, 1e-5).unwrap();
        let diff = (a_out - a_in).abs();
        assert!((diff - std::f64::consts::PI).abs() < 0.01,
            "expected ~PI between in/out, got {}", diff);
    }
}
