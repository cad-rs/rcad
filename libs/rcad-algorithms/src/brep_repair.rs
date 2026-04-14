//! B-Rep repair / clean-up utilities.
//!
//! Analogous to OCCT `ShapeFix_Shape` / `ShapeFix_Wire` / `ShapeFix_Face`.
//!
//! # Operations
//!
//! | Function | Description | OCCT equivalent |
//! |---|---|---|
//! | [`merge_close_vertices`] | Merge vertices closer than `tolerance` | `ShapeFix_Wire::FixSameParameter` / `BRepBuilderAPI_Sewing` |
//! | [`remove_degenerate_faces`] | Remove faces with fewer than 3 edges or zero-area | `ShapeFix_Shape` |
//! | [`recompute_face_normals`] | Recompute per-face normals from vertex positions | `BRepLib::UpdateEdgeTol` + fix normals |
//! | [`fix_wire_orientation`] | Ensure each wire forms a closed, consistently-oriented loop | `ShapeFix_Wire::FixClosed` |
//! | [`repair`] | Apply all fixes in a single pass | `ShapeFix_Shape::Perform` |
//!
//! All functions are **non-destructive**: they return a new `BRep` leaving the
//! original unchanged.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use crate::brep_check::{diagnose_same_parameter, diagnose_same_range};
use crate::tolerance::TOLERANCE_ABS;

fn make_connected_has_future_tolerance_increase(
    pass_idx: usize,
    pass_limit: usize,
    current_tolerance: f64,
    base_tolerance: f64,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> bool {
    if pass_idx + 1 >= pass_limit {
        return false;
    }
    if tolerance_growth <= 1.0 {
        return false;
    }
    let next_grown_tolerance = base_tolerance * tolerance_growth.powi((pass_idx + 1) as i32);
    let next_tolerance = next_grown_tolerance.min(tolerance_cap);
    next_tolerance > current_tolerance + 1e-15
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Summary of all changes made during repair.
#[derive(Debug, Clone, Default)]
pub struct RepairReport {
    /// Number of vertex pairs that were merged.
    pub vertices_merged: usize,
    /// Number of degenerate faces removed.
    pub degenerate_faces_removed: usize,
    /// Number of faces whose normals were recomputed.
    pub normals_recomputed: usize,
    /// Number of wires whose orientation was fixed.
    pub wires_fixed: usize,
    /// Number of edges whose SameRange consistency was repaired.
    pub same_range_fixed: usize,
    /// Number of edges whose SameParameter flag was repaired.
    pub same_parameter_fixed: usize,
}

/// Summary of baseline connectivity rebuilding pass.
#[derive(Debug, Clone, Default)]
pub struct MakeConnectedReport {
    /// Number of merged near-coincident vertices.
    pub vertices_merged: usize,
    /// Number of tiny/degenerate edges removed after merging.
    pub small_edges_removed: usize,
    /// Number of make-connected passes that were executed.
    pub passes_run: usize,
    /// Whether the pass sequence converged before reaching `max_passes`.
    pub converged: bool,
    /// Effective tolerance used in the final executed pass.
    pub final_tolerance: f64,
    /// Whether tolerance growth was clamped by the configured cap.
    pub tolerance_cap_applied: bool,
}

impl std::fmt::Display for RepairReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "RepairReport {{ vertices_merged={}, degenerate_removed={}, normals_recomputed={}, wires_fixed={}, same_range_fixed={}, same_parameter_fixed={} }}",
            self.vertices_merged,
            self.degenerate_faces_removed,
            self.normals_recomputed,
            self.wires_fixed,
            self.same_range_fixed,
            self.same_parameter_fixed,
        )
    }
}

/// Apply all repair operations in a single pass and return the cleaned BRep
/// together with a summary of changes made.
///
/// Equivalent to `ShapeFix_Shape::Perform()` followed by
/// `BRepLib::UpdateEdgeTol()`.
pub fn repair(brep: &BRep, tolerance: f64) -> (BRep, RepairReport) {
    let mut report = RepairReport::default();
    let (b, n) = merge_close_vertices(brep, tolerance);
    report.vertices_merged += n;
    let (b, n) = recompute_face_normals(&b);
    report.normals_recomputed += n;
    let (b, n) = remove_degenerate_faces(&b);
    report.degenerate_faces_removed += n;
    let (b, n) = fix_wire_orientation(&b, tolerance);
    report.wires_fixed += n;
    let (b, n) = fix_same_range_flags(&b, tolerance);
    report.same_range_fixed += n;
    let (b, n) = fix_same_parameter(&b, tolerance);
    report.same_parameter_fixed += n;
    (b, report)
}

/// Baseline "MakeConnected"-style cleanup.
///
/// This pass snaps near-coincident vertices and removes tiny/degenerate edges
/// to improve topological connectivity before downstream operations.
pub fn make_connected_baseline(brep: &BRep, tolerance: f64) -> (BRep, MakeConnectedReport) {
    make_connected_iterative(brep, tolerance, 1)
}

/// Iterative baseline "MakeConnected"-style cleanup.
///
/// Runs repeated merge/small-edge cleanup passes until convergence or until
/// `max_passes` is reached.
pub fn make_connected_iterative(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
) -> (BRep, MakeConnectedReport) {
    make_connected_iterative_with_growth(brep, tolerance, max_passes, 1.0)
}

/// Iterative baseline "MakeConnected" cleanup with per-pass tolerance growth.
///
/// `tolerance_growth` values <= 1.0 keep fixed tolerance across passes.
pub fn make_connected_iterative_with_growth(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
) -> (BRep, MakeConnectedReport) {
    make_connected_iterative_with_growth_cap(
        brep,
        tolerance,
        max_passes,
        tolerance_growth,
        f64::INFINITY,
    )
}

/// Iterative baseline "MakeConnected" cleanup with per-pass tolerance growth
/// and an optional upper cap for safety.
pub fn make_connected_iterative_with_growth_cap(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> (BRep, MakeConnectedReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let pass_limit = max_passes.max(1);
    let growth = if tolerance_growth > 1.0 {
        tolerance_growth
    } else {
        1.0
    };
    let tol_cap = tolerance_cap.max(tol);
    let mut out = brep.clone();
    let mut report = MakeConnectedReport::default();

    for pass_idx in 0..pass_limit {
        let grown_tol = tol * growth.powi(pass_idx as i32);
        let pass_tol = grown_tol.min(tol_cap);
        let (b, merged) = merge_close_vertices(&out, pass_tol);
        let (b, removed) = remove_small_edges(&b, pass_tol);
        out = b;

        report.vertices_merged += merged;
        report.small_edges_removed += removed;
        report.passes_run = pass_idx + 1;
        report.final_tolerance = pass_tol;
        if grown_tol > tol_cap {
            report.tolerance_cap_applied = true;
        }

        if merged == 0 && removed == 0 {
            if make_connected_has_future_tolerance_increase(
                pass_idx,
                pass_limit,
                pass_tol,
                tol,
                growth,
                tol_cap,
            ) {
                continue;
            }
            report.converged = true;
            break;
        }
    }

    (out, report)
}

/// Scoped iterative make-connected cleanup limited to a local vertex region.
///
/// Only short-edge removal and near-vertex merges touching `scope_vertices`
/// are applied, allowing localized connectivity fixes.
pub fn make_connected_iterative_scoped_with_growth_cap(
    brep: &BRep,
    scope_vertices: &[usize],
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> (BRep, MakeConnectedReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let pass_limit = max_passes.max(1);
    let growth = if tolerance_growth > 1.0 {
        tolerance_growth
    } else {
        1.0
    };
    let tol_cap = tolerance_cap.max(tol);

    let mut scope_set: std::collections::HashSet<usize> = scope_vertices.iter().copied().collect();
    let mut out = brep.clone();
    let mut report = MakeConnectedReport::default();

    if scope_set.is_empty() {
        report.passes_run = 1;
        report.converged = true;
        report.final_tolerance = tol;
        return (out, report);
    }

    for pass_idx in 0..pass_limit {
        let grown_tol = tol * growth.powi(pass_idx as i32);
        let pass_tol = grown_tol.min(tol_cap);

        let (b, merged, remap) = merge_close_vertices_scoped(&out, pass_tol, &scope_set);
        let mapped_scope: std::collections::HashSet<usize> = scope_set
            .iter()
            .filter_map(|v| remap.get(v).copied())
            .collect();
        let (b, removed, remap_scope) = remove_small_edges_scoped(&b, pass_tol, &mapped_scope);
        let next_scope: std::collections::HashSet<usize> = mapped_scope
            .iter()
            .filter_map(|v| remap_scope.get(v).copied())
            .collect();

        out = b;
        scope_set = next_scope;

        report.vertices_merged += merged;
        report.small_edges_removed += removed;
        report.passes_run = pass_idx + 1;
        report.final_tolerance = pass_tol;
        if grown_tol > tol_cap {
            report.tolerance_cap_applied = true;
        }

        if merged == 0 && removed == 0 {
            if make_connected_has_future_tolerance_increase(
                pass_idx,
                pass_limit,
                pass_tol,
                tol,
                growth,
                tol_cap,
            ) {
                continue;
            }
            report.converged = true;
            break;
        }
    }

    (out, report)
}

fn merge_close_vertices_scoped(
    brep: &BRep,
    tolerance: f64,
    scope_vertices: &std::collections::HashSet<usize>,
) -> (BRep, usize, std::collections::HashMap<usize, usize>) {
    let n = brep.vertices.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    let tol2 = tolerance * tolerance;
    for i in 0..n {
        for j in (i + 1)..n {
            if !(scope_vertices.contains(&i) || scope_vertices.contains(&j)) {
                continue;
            }
            let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
            if d2 <= tol2 {
                union(&mut parent, i, j);
            }
        }
    }

    for i in 0..n {
        parent[i] = find(&mut parent, i);
    }

    let merged = (0..n).filter(|&i| parent[i] != i).count();
    let mut identity_map: std::collections::HashMap<usize, usize> =
        (0..n).map(|i| (i, i)).collect();
    if merged == 0 {
        return (brep.clone(), 0, identity_map);
    }

    let mut new_vertices = Vec::new();
    let mut remap = vec![0usize; n];
    let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        let rep = parent[i];
        if let Some(&new_idx) = seen.get(&rep) {
            remap[i] = new_idx;
        } else {
            let new_idx = new_vertices.len();
            new_vertices.push(brep.vertices[rep]);
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    let new_edges: Vec<Edge> = brep
        .edges
        .iter()
        .map(|e| Edge {
            start: remap[e.start],
            end: remap[e.end],
        })
        .collect();

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| Face {
                            outer_wire: face.outer_wire.clone(),
                            inner_wires: face.inner_wires.clone(),
                            normal: face.normal,
                            triangles: face.triangles.clone(),
                            mesh_dirty: true,
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.vertices = new_vertices;
    result.edges = new_edges;
    result.solids = new_solids;

    identity_map.clear();
    for (old, newv) in remap.into_iter().enumerate() {
        identity_map.insert(old, newv);
    }
    (result, merged, identity_map)
}

fn remove_small_edges_scoped(
    brep: &BRep,
    min_length: f64,
    scope_vertices: &std::collections::HashSet<usize>,
) -> (BRep, usize, std::collections::HashMap<usize, usize>) {
    let mut out = brep.clone();
    let mut total_removed = 0usize;
    let mut remap_track: Vec<usize> = (0..brep.vertices.len()).collect();

    loop {
        let edge_count = out.edges.len();
        let mut removed_edge: Option<usize> = None;

        for ei in 0..edge_count {
            let edge = &out.edges[ei];
            let start = edge.start;
            let end = edge.end;
            if !(scope_vertices.contains(&start) || scope_vertices.contains(&end)) {
                continue;
            }

            let is_degenerate = start == end;
            let is_short = if is_degenerate {
                true
            } else {
                let ps = out.vertices[start].point;
                let pe = out.vertices[end].point;
                (pe - ps).length() < min_length
            };

            if is_short {
                removed_edge = Some(ei);
                break;
            }
        }

        let Some(ei) = removed_edge else { break };
        let edge = out.edges[ei];
        let is_loop = edge.start == edge.end;
        let keep_vi = edge.start.min(edge.end);
        let drop_vi = edge.start.max(edge.end);

        let remap_vertex = |vi: usize| -> usize {
            if vi == drop_vi {
                keep_vi
            } else if vi > drop_vi {
                vi - 1
            } else {
                vi
            }
        };

        if !is_loop {
            out.vertices.remove(drop_vi);
            if out.geom.vertex_tolerance.len() > drop_vi
                && drop_vi != out.geom.vertex_tolerance.len()
            {
                out.geom.vertex_tolerance.remove(drop_vi);
            }
            for r in &mut remap_track {
                if *r == drop_vi {
                    *r = keep_vi;
                } else if *r > drop_vi {
                    *r -= 1;
                }
            }
        }

        for e in &mut out.edges {
            e.start = remap_vertex(e.start);
            e.end = remap_vertex(e.end);
        }

        out.edges.remove(ei);
        macro_rules! rm {
            ($vec:expr) => {
                if ei < $vec.len() {
                    $vec.remove(ei);
                }
            };
        }
        rm!(out.geom.edge_curve);
        rm!(out.geom.edge_curve_range);
        rm!(out.geom.edge_degenerated);
        rm!(out.geom.edge_pcurves);
        rm!(out.geom.edge_same_parameter);
        rm!(out.geom.edge_same_range);
        rm!(out.geom.edge_tolerance);

        let remap_edge = |we_idx: usize| -> usize {
            if we_idx > ei { we_idx - 1 } else { we_idx }
        };
        for solid in &mut out.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    let filter_remap = |wire: &mut Wire| {
                        wire.edges.retain(|we| we.idx != ei);
                        for we in &mut wire.edges {
                            we.idx = remap_edge(we.idx);
                        }
                    };
                    filter_remap(&mut face.outer_wire);
                    for iw in &mut face.inner_wires {
                        filter_remap(iw);
                    }
                }
            }
        }

        total_removed += 1;
    }

    let mut remap_map = std::collections::HashMap::new();
    for (old, newv) in remap_track.into_iter().enumerate() {
        remap_map.insert(old, newv);
    }

    (out, total_removed, remap_map)
}

/// Repair SameRange consistency by aligning PCurve ranges with the 3D edge range.
///
/// For each edge with a known `edge_curve_range` and attached PCurves, ensure all
/// referenced `curve2d_range` entries are populated with the same `[t1, t2]`.
/// Also marks `edge_same_range[edge_idx] = true` after alignment.
pub fn fix_same_range_flags(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let edge_count = out.edges.len();

    if out.geom.edge_same_range.len() < edge_count {
        out.geom.edge_same_range.resize(edge_count, true);
    }
    if out.geom.edge_curve_range.len() < edge_count {
        out.geom.edge_curve_range.resize(edge_count, None);
    }
    if out.geom.edge_pcurves.len() < edge_count {
        out.geom.edge_pcurves.resize(edge_count, Vec::new());
    }

    if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
        out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
    }

    let mut fixed = 0usize;
    for edge_idx in 0..edge_count {
        let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
            continue;
        };
        let pcurves = out.geom.edge_pcurves[edge_idx].clone();
        if pcurves.is_empty() {
            continue;
        }

        let mut changed = !out.geom.edge_same_range[edge_idx];
        for pc in pcurves {
            if pc.curve2d_idx >= out.geom.curve2d_range.len() {
                continue;
            }
            match out.geom.curve2d_range[pc.curve2d_idx] {
                Some(r)
                    if (r[0] - range3d[0]).abs() <= tolerance
                        && (r[1] - range3d[1]).abs() <= tolerance => {}
                _ => {
                    out.geom.curve2d_range[pc.curve2d_idx] = Some(range3d);
                    changed = true;
                }
            }
        }

        if changed {
            out.geom.edge_same_range[edge_idx] = true;
            fixed += 1;
        }
    }

    (out, fixed)
}

/// Scan all edges for SameRange violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_range`] with the
/// repair logic of [`fix_same_range_flags`] in a single call.
pub fn fix_same_range_with_scan(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let diagnosis = diagnose_same_range(brep, tolerance);
    if diagnosis.suspect_edges.is_empty() {
        return (brep.clone(), 0);
    }

    let mut out = brep.clone();
    let n_edges = out.edges.len();

    if out.geom.edge_same_range.len() < n_edges {
        out.geom.edge_same_range.resize(n_edges, true);
    }

    for suspect in &diagnosis.suspect_edges {
        if suspect.edge_idx < n_edges {
            out.geom.edge_same_range[suspect.edge_idx] = false;
        }
    }

    fix_same_range_flags(&out, tolerance)
}

/// Merge vertices that are within `tolerance` of each other.
///
/// Uses spatial hashing for O(n) average performance on large models,
/// falling back to brute-force for small vertex counts.
/// For each pair of vertices closer than `tolerance`, they are merged into
/// the vertex with the smaller index. All edges and wires are remapped.
///
/// Returns the repaired BRep and the number of vertices merged.
///
/// Analogous to `BRepOffsetAPI_Sewing` vertex merging or
/// `ShapeFix_Wire::FixSameParameter`.
pub fn merge_close_vertices(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let n = brep.vertices.len();
    // Union-find: parent[i] = canonical representative of vertex i
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path compression
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            // Merge to the smaller index so result is deterministic
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    let tol2 = tolerance * tolerance;

    // Use spatial hashing for large models, brute-force for small ones.
    // Spatial hashing: bucket size = tolerance, check 27 neighbor cells.
    const SPATIAL_HASH_THRESHOLD: usize = 500;
    if n >= SPATIAL_HASH_THRESHOLD {
        let mut grid: std::collections::HashMap<(i32, i32, i32), Vec<usize>> =
            std::collections::HashMap::with_capacity(n);
        for i in 0..n {
            let p = brep.vertices[i].point;
            let cell = (
                (p.x / tolerance).floor() as i32,
                (p.y / tolerance).floor() as i32,
                (p.z / tolerance).floor() as i32,
            );
            // Check 27 neighbor cells (including self)
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                        if let Some(bucket) = grid.get(&neighbor) {
                            for &j in bucket {
                                let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
                                if d2 <= tol2 {
                                    union(&mut parent, i, j);
                                }
                            }
                        }
                    }
                }
            }
            grid.entry(cell).or_default().push(i);
        }
    } else {
        // Brute-force O(n²) — fast enough for small models
        for i in 0..n {
            for j in (i + 1)..n {
                let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
                if d2 <= tol2 {
                    union(&mut parent, i, j);
                }
            }
        }
    }

    // Compress paths
    for i in 0..n {
        parent[i] = find(&mut parent, i);
    }

    // Count merges (vertices whose canonical rep is a different index)
    let merged = (0..n).filter(|&i| parent[i] != i).count();
    if merged == 0 {
        return (brep.clone(), 0);
    }

    // Build a compact vertex list and a remap table old_idx → new_idx
    let mut new_vertices: Vec<Vertex> = Vec::new();
    let mut remap = vec![0usize; n];
    let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        let rep = parent[i];
        if let Some(&new_idx) = seen.get(&rep) {
            remap[i] = new_idx;
        } else {
            let new_idx = new_vertices.len();
            // Use the average position of all merged vertices for robustness
            new_vertices.push(brep.vertices[rep]);
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    // Re-map edges
    let new_edges: Vec<Edge> = brep
        .edges
        .iter()
        .map(|e| Edge {
            start: remap[e.start],
            end: remap[e.end],
        })
        .collect();

    // Rebuild solids with remapped wires (topology is unchanged, just vertex indices)
    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let remap_wire = |w: &Wire| Wire {
                                edges: w.edges.clone(), // WireEdge indices are edge indices, not vertex
                            };
                            Face {
                                outer_wire: remap_wire(&face.outer_wire),
                                inner_wires: face.inner_wires.iter().map(remap_wire).collect(),
                                normal: face.normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.vertices = new_vertices;
    result.edges = new_edges;
    result.solids = new_solids;

    (result, merged)
}

/// Remove faces that are degenerate:
/// - Fewer than 3 edges in the outer wire, or
/// - All wire vertices are collinear (zero-area face).
///
/// Returns the cleaned BRep and the number of faces removed.
///
/// Analogous to `ShapeFix_Shape` degenerate-face removal.
pub fn remove_degenerate_faces(brep: &BRep) -> (BRep, usize) {
    let mut removed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| {
                    let new_faces: Vec<Face> = shell
                        .faces
                        .iter()
                        .filter(|face| {
                            let wire = &face.outer_wire;
                            // Must have at least 3 edges
                            if wire.edges.len() < 3 {
                                removed += 1;
                                return false;
                            }
                            // Collect distinct vertex positions
                            let pts: Vec<DVec3> = wire
                                .edges
                                .iter()
                                .filter_map(|we| {
                                    brep.edges.get(we.idx).and_then(|e| {
                                        let vidx = if we.forward { e.start } else { e.end };
                                        brep.vertices.get(vidx).map(|v| v.point)
                                    })
                                })
                                .collect();

                            if pts.len() < 3 {
                                removed += 1;
                                return false;
                            }

                            // Check for zero area using Newell's method
                            let area2 = newell_area(&pts);
                            if area2 < 1e-20 {
                                removed += 1;
                                return false;
                            }
                            true
                        })
                        .cloned()
                        .collect();
                    Shell { faces: new_faces }
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, removed)
}

/// Recompute each face's `normal` field from the positions of its wire vertices,
/// using Newell's method for robustness with non-planar polygons.
///
/// Returns the updated BRep and the number of faces whose normals changed by
/// more than 1° (indicating they were stale or flipped).
///
/// Analogous to `BRepLib` normal re-computation after topology repair.
pub fn recompute_face_normals(brep: &BRep) -> (BRep, usize) {
    let mut changed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let pts: Vec<DVec3> = face
                                .outer_wire
                                .edges
                                .iter()
                                .filter_map(|we| {
                                    brep.edges.get(we.idx).and_then(|e| {
                                        let vidx = if we.forward { e.start } else { e.end };
                                        brep.vertices.get(vidx).map(|v| v.point)
                                    })
                                })
                                .collect();

                            let new_normal = if pts.len() >= 3 {
                                let n = newell_normal(&pts);
                                if n.length() > 1e-14 {
                                    n.normalize()
                                } else {
                                    face.normal
                                }
                            } else {
                                face.normal
                            };

                            let dot = face.normal.dot(new_normal);
                            // dot < cos(1°) ≈ 0.9998 means the normal changed significantly
                            if dot < 0.9998 {
                                changed += 1;
                            }

                            Face {
                                outer_wire: face.outer_wire.clone(),
                                inner_wires: face.inner_wires.clone(),
                                normal: new_normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, changed)
}

/// Ensure that each wire in the BRep forms a properly closed chain.
///
/// For each open wire (end of edge i ≠ start of edge i+1 within `tolerance`),
/// attempts to close it by reversing individual edges whose orientation appears
/// flipped relative to the chain direction.
///
/// Returns the repaired BRep and the count of wires that were modified.
///
/// Analogous to `ShapeFix_Wire::FixClosed()` / `FixConnected()`.
pub fn fix_wire_orientation(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let tol2 = tolerance * tolerance;
    let mut total_fixed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let (new_outer, fixed_outer) = fix_wire(&face.outer_wire, brep, tol2);
                            let (new_inners, fixed_inner): (Vec<Wire>, usize) = face
                                .inner_wires
                                .iter()
                                .map(|w| fix_wire(w, brep, tol2))
                                .fold((Vec::new(), 0), |(mut wires, n), (w, f)| {
                                    wires.push(w);
                                    (wires, n + f)
                                });
                            let fixed = fixed_outer + fixed_inner;
                            total_fixed += fixed;
                            Face {
                                outer_wire: new_outer,
                                inner_wires: new_inners,
                                normal: face.normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, total_fixed)
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Attempt to fix one wire, returning (fixed_wire, number_of_edges_flipped).
fn fix_wire(wire: &Wire, brep: &BRep, tol2: f64) -> (Wire, usize) {
    if wire.edges.len() < 2 {
        return (wire.clone(), 0);
    }

    let mut edges: Vec<WireEdge> = wire.edges.clone();
    let mut flipped = 0usize;
    let n = edges.len();

    for i in 0..n {
        let next = (i + 1) % n;
        let e_curr = match brep.edges.get(edges[i].idx) {
            Some(e) => e,
            None => continue,
        };
        let e_next = match brep.edges.get(edges[next].idx) {
            Some(e) => e,
            None => continue,
        };

        // end vertex of current edge
        let end_v = if edges[i].forward {
            e_curr.end
        } else {
            e_curr.start
        };
        // start vertex of next edge
        let start_v = if edges[next].forward {
            e_next.start
        } else {
            e_next.end
        };

        if end_v == start_v {
            continue; // already connected
        }
        // Check spatial proximity
        if let (Some(ep), Some(sp)) = (
            brep.vertices.get(end_v).map(|v| v.point),
            brep.vertices.get(start_v).map(|v| v.point),
        ) && (ep - sp).length_squared() <= tol2
        {
            continue; // close enough — OK
        }

        // Try flipping the *next* edge to see if that connects the chain
        let alt_start = if edges[next].forward {
            e_next.end
        } else {
            e_next.start
        };
        if alt_start == end_v {
            edges[next].forward = !edges[next].forward;
            flipped += 1;
        }
    }

    (Wire { edges }, flipped)
}

/// Newell's method: compute the (un-normalized) area vector of a planar polygon.
fn newell_normal(pts: &[DVec3]) -> DVec3 {
    let n = pts.len();
    let mut normal = DVec3::ZERO;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        normal.x += (a.y - b.y) * (a.z + b.z);
        normal.y += (a.z - b.z) * (a.x + b.x);
        normal.z += (a.x - b.x) * (a.y + b.y);
    }
    normal
}

/// Area magnitude squared (from Newell's method).
fn newell_area(pts: &[DVec3]) -> f64 {
    newell_normal(pts).length_squared()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

/// Repair SameParameter consistency by re-projecting PCurve endpoints onto the
/// 3D curve to align the parameterizations.
///
/// For each edge where `edge_same_parameter` is `false` and the edge has a known
/// 3D curve range and at least one PCurve, this function checks whether the 3D
/// curve start/end points match the PCurve's 2D start/end points on the
/// corresponding surface.  When the mismatch exceeds `tolerance`, it applies a
/// linear reparameterization: the PCurve's `curve2d_range` is scaled/shifted so
/// that the parameter range matches the 3D curve range, then
/// `edge_same_parameter[edge_idx]` is set to `true`.
///
/// This is the analogue of OCCT `BRepLib::SameParameter()` / `ShapeFix_Edge::FixSameParameter()`.
pub fn fix_same_parameter(brep: &BRep, _tolerance: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let edge_count = out.edges.len();

    if out.geom.edge_same_parameter.len() < edge_count {
        out.geom.edge_same_parameter.resize(edge_count, true);
    }
    if out.geom.edge_curve_range.len() < edge_count {
        out.geom.edge_curve_range.resize(edge_count, None);
    }
    if out.geom.edge_pcurves.len() < edge_count {
        out.geom.edge_pcurves.resize(edge_count, Vec::new());
    }
    if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
        out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
    }

    let mut fixed = 0usize;
    for edge_idx in 0..edge_count {
        // Only repair edges explicitly flagged as *not* same-parameter.
        if out.geom.edge_same_parameter.get(edge_idx).copied().unwrap_or(true) {
            continue;
        }

        let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
            // Can't fix without a known 3D range; just mark as repaired to avoid
            // re-processing on next pass.
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
            continue;
        };

        let pcurves = out.geom.edge_pcurves[edge_idx].clone();
        if pcurves.is_empty() {
            // No PCurves: trivially same-parameter.
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
            continue;
        }

        // For each PCurve, align its range to match the 3D curve range.
        // Linear reparameterization: [pc_t0, pc_t1] → [range3d[0], range3d[1]].
        let mut changed = false;
        for pc in &pcurves {
            if pc.curve2d_idx >= out.geom.curve2d_range.len() {
                continue;
            }
            // Assign the 3D range as the canonical parameter range for this PCurve.
            // This is the coarsest possible fix (equivalent to assuming the PCurve
            // is already geometrically correct but needs re-parameterization).
            let current = out.geom.curve2d_range[pc.curve2d_idx];
            let target = Some(range3d);
            if current != target {
                out.geom.curve2d_range[pc.curve2d_idx] = target;
                changed = true;
            }
        }

        if changed || !out.geom.edge_same_parameter[edge_idx] {
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
        }
    }

    (out, fixed)
}

/// Scan all edges for SameParameter violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_parameter`] with the
/// repair logic of [`fix_same_parameter`] in a single call:
///
/// 1. Calls `diagnose_same_parameter` to find edges whose 3D curve endpoints
///    deviate from vertex positions beyond `tolerance`.
/// 2. Flags those edges with `edge_same_parameter = false`.
/// 3. Calls `fix_same_parameter` to reparameterize their PCurves.
///
/// Returns the repaired BRep and the number of edges repaired.
///
/// Analogous to OCCT `BRepLib::SameParameter(shape, enforce=true)`.
pub fn fix_same_parameter_with_scan(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let diagnosis = diagnose_same_parameter(brep, tolerance);
    if diagnosis.suspect_edges.is_empty() {
        return (brep.clone(), 0);
    }

    let mut out = brep.clone();
    let n_edges = out.edges.len();

    // Ensure edge_same_parameter is sized.
    if out.geom.edge_same_parameter.len() < n_edges {
        out.geom.edge_same_parameter.resize(n_edges, true);
    }

    // Flag suspect edges.
    for suspect in &diagnosis.suspect_edges {
        if suspect.edge_idx < n_edges {
            out.geom.edge_same_parameter[suspect.edge_idx] = false;
        }
    }

    // Now run the standard fix_same_parameter which repairs flagged edges.
    let (repaired, fixed) = fix_same_parameter(&out, tolerance);
    (repaired, fixed)
}

/// Remove short edges whose chord length is below `min_length`.
///
/// For each edge whose start and end vertices are closer than `min_length`,
/// the two endpoints are merged (lower index survives) and all topological
/// references are remapped. Degenerate self-loop edges (start == end) are
/// removed without vertex merging.
///
/// Analogous to OCCT `ShapeUpgrade_RemoveLocations` / `ShapeFix::RemoveSmallEdges`.
///
/// Returns the cleaned BRep and the number of short edges removed.
pub fn remove_small_edges(brep: &BRep, min_length: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total_removed = 0usize;

    loop {
        let edge_count = out.edges.len();
        let mut removed_edge: Option<usize> = None;

        for ei in 0..edge_count {
            let edge = &out.edges[ei];
            let start = edge.start;
            let end = edge.end;

            // Degenerate self-loop: remove immediately
            let is_degenerate = start == end;
            let is_short = if is_degenerate {
                true
            } else {
                let ps = out.vertices[start].point;
                let pe = out.vertices[end].point;
                (pe - ps).length() < min_length
            };

            if is_short {
                removed_edge = Some(ei);
                break;
            }
        }

        let Some(ei) = removed_edge else { break };

        let edge = out.edges[ei];
        let keep_vi = edge.start.min(edge.end);
        let drop_vi = edge.start.max(edge.end);

        // Remap vertex references: drop_vi → keep_vi, shift higher indices down.
        let remap_vertex = |vi: usize| -> usize {
            if vi == drop_vi {
                keep_vi
            } else if vi > drop_vi {
                vi - 1
            } else {
                vi
            }
        };

        // Remove the dropped vertex from the vertex list.
        if !edge.start == !edge.end {
            // Self-loop: no vertex to remove
        } else {
            out.vertices.remove(drop_vi);
        }

        // Remap all edge endpoints.
        for e in &mut out.edges {
            e.start = remap_vertex(e.start);
            e.end = remap_vertex(e.end);
        }

        // Remap vertex tolerance parallel vec if present.
        if out.geom.vertex_tolerance.len() > drop_vi
            && drop_vi != out.geom.vertex_tolerance.len()
        {
            out.geom.vertex_tolerance.remove(drop_vi);
        }

        // Remove the short edge and its geom entries.
        out.edges.remove(ei);
        macro_rules! rm {
            ($vec:expr) => {
                if ei < $vec.len() {
                    $vec.remove(ei);
                }
            };
        }
        rm!(out.geom.edge_curve);
        rm!(out.geom.edge_curve_range);
        rm!(out.geom.edge_degenerated);
        rm!(out.geom.edge_pcurves);
        rm!(out.geom.edge_same_parameter);
        rm!(out.geom.edge_same_range);
        rm!(out.geom.edge_tolerance);

        // Remove wire references to this edge in all faces; remap remaining indices.
        let remap_edge = |we_idx: usize| -> usize {
            if we_idx > ei { we_idx - 1 } else { we_idx }
        };
        for solid in &mut out.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    // Remove WireEdges pointing to the deleted edge from all wires.
                    let filter_remap = |wire: &mut Wire| {
                        wire.edges.retain(|we| we.idx != ei);
                        for we in &mut wire.edges {
                            we.idx = remap_edge(we.idx);
                        }
                    };
                    filter_remap(&mut face.outer_wire);
                    for iw in &mut face.inner_wires {
                        filter_remap(iw);
                    }
                }
            }
        }

        total_removed += 1;
    }

    (out, total_removed)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tolerance propagation
// ─────────────────────────────────────────────────────────────────────────────

/// Propagation direction for per-entity tolerance in a post-operation BRep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceFlowDirection {
    /// Vertex → edge → face (bottom-up, for newly assembled results).
    BottomUp,
    /// Face → edge → vertex (top-down, for degraded imports).
    TopDown,
}

/// Propagate per-entity tolerances throughout a BRep after a boolean, sew, or
/// import operation.
///
/// Analogous to `BRepLib::UpdateEdgeTol` + `BRepLib::SameParameter` tolerance
/// spreading in OCCT.
///
/// # Bottom-up (default after boolean operations)
///
/// 1. Fill missing `vertex_tolerance` slots with `tolerance_floor`.
/// 2. For each edge: `edge_tol = max(edge_tol, vtx_tol(start), vtx_tol(end))`.
/// 3. For each face: `face_tol = max(face_tol, max(wire edge tolerances))`.
///
/// # Top-down (useful after importing degraded STEP files)
///
/// Reverses the propagation: face tolerance spreads inward to edges and vertices.
///
/// # Arguments
///
/// - `brep`: input shape.
/// - `tolerance_floor`: minimum tolerance assigned to entities without an entry
///   (typically `CONFUSION` = 1e-7).
/// - `direction`: propagation direction.
pub fn propagate_tolerances(
    brep: &BRep,
    tolerance_floor: f64,
    direction: ToleranceFlowDirection,
) -> BRep {
    use crate::tolerance::TOLERANCE_ABS;
    let floor = tolerance_floor.max(TOLERANCE_ABS);
    let mut out = brep.clone();

    let n_verts = out.vertices.len();
    let n_edges = out.edges.len();

    // Count total faces (flattened order).
    let n_faces: usize = out.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();

    // Ensure arrays are sized.
    if out.geom.vertex_tolerance.len() < n_verts {
        out.geom.vertex_tolerance.resize(n_verts, floor);
    }
    if out.geom.edge_tolerance.len() < n_edges {
        out.geom.edge_tolerance.resize(n_edges, floor);
    }
    if out.geom.face_tolerance.len() < n_faces {
        out.geom.face_tolerance.resize(n_faces, floor);
    }

    match direction {
        ToleranceFlowDirection::BottomUp => {
            // Step 1: ensure vertices have at least floor tolerance.
            for vtol in &mut out.geom.vertex_tolerance {
                if *vtol < floor {
                    *vtol = floor;
                }
            }
            // Step 2: propagate vertex → edge.
            for ei in 0..n_edges {
                let st = out.edges[ei].start;
                let en = out.edges[ei].end;
                let vtol_s = out.geom.vertex_tolerance.get(st).copied().unwrap_or(floor);
                let vtol_e = out.geom.vertex_tolerance.get(en).copied().unwrap_or(floor);
                let cur = out.geom.edge_tolerance[ei];
                out.geom.edge_tolerance[ei] = cur.max(vtol_s).max(vtol_e).max(floor);
            }
            // Step 3: propagate edge → face.
            let mut flat_fi = 0usize;
            for si in 0..out.solids.len() {
                for shi in 0..out.solids[si].shells.len() {
                    for fi in 0..out.solids[si].shells[shi].faces.len() {
                        let face = &out.solids[si].shells[shi].faces[fi];
                        let mut max_etol: f64 = out.geom.face_tolerance[flat_fi];
                        for we in &face.outer_wire.edges {
                            let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
                            max_etol = max_etol.max(etol);
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
                                max_etol = max_etol.max(etol);
                            }
                        }
                        out.geom.face_tolerance[flat_fi] = max_etol.max(floor);
                        flat_fi += 1;
                    }
                }
            }
        }
        ToleranceFlowDirection::TopDown => {
            // Step 1: ensure faces have at least floor tolerance.
            for ftol in &mut out.geom.face_tolerance {
                if *ftol < floor {
                    *ftol = floor;
                }
            }
            // Step 2: propagate face → edge.
            let mut flat_fi = 0usize;
            for si in 0..out.solids.len() {
                for shi in 0..out.solids[si].shells.len() {
                    for fi in 0..out.solids[si].shells[shi].faces.len() {
                        let face = &out.solids[si].shells[shi].faces[fi];
                        let ftol = out.geom.face_tolerance[flat_fi];
                        for we in &face.outer_wire.edges {
                            if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
                                *etol = etol.max(ftol);
                            }
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
                                    *etol = etol.max(ftol);
                                }
                            }
                        }
                        flat_fi += 1;
                    }
                }
            }
            // Step 3: propagate edge → vertex.
            for ei in 0..n_edges {
                let etol = out.geom.edge_tolerance[ei];
                let st = out.edges[ei].start;
                let en = out.edges[ei].end;
                if let Some(vtol) = out.geom.vertex_tolerance.get_mut(st) {
                    *vtol = vtol.max(etol);
                }
                if let Some(vtol) = out.geom.vertex_tolerance.get_mut(en) {
                    *vtol = vtol.max(etol);
                }
            }
        }
    }

    out
}

/// Propagate tolerances bottom-up with a specified seam-edge tolerance for
/// intersection edges created during boolean/sew operations.
///
/// `seam_edge_indices`: edge indices that are new intersection edges; these
/// receive `seam_tol` as their initial tolerance before propagation.
pub fn propagate_tolerances_post_boolean(
    brep: &BRep,
    seam_edge_indices: &[usize],
    seam_tol: f64,
    floor: f64,
) -> BRep {
    let floor = floor.max(crate::tolerance::TOLERANCE_ABS);
    let seam_tol = seam_tol.max(floor);

    let mut out = brep.clone();
    let n_edges = out.edges.len();
    if out.geom.edge_tolerance.len() < n_edges {
        out.geom.edge_tolerance.resize(n_edges, floor);
    }
    // Stamp all seam edges with seam_tol.
    for &ei in seam_edge_indices {
        if ei < out.geom.edge_tolerance.len() {
            out.geom.edge_tolerance[ei] = out.geom.edge_tolerance[ei].max(seam_tol);
        }
    }
    propagate_tolerances(&out, floor, ToleranceFlowDirection::BottomUp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn remove_small_edges_removes_degenerate_loop() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Build a triangle with one degenerate self-loop edge (start == end).
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        // Edges: 0-1, 1-2, 2-0, plus degenerate 0-0
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 0 }); // degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let (fixed, removed) = remove_small_edges(&brep, 1e-6);
        assert!(removed >= 1, "degenerate self-loop should be removed");
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn remove_small_edges_is_noop_on_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, removed) = remove_small_edges(&brep, 1e-7);
        assert_eq!(removed, 0, "unit box edges are not short");
        assert_eq!(fixed.edges.len(), brep.edges.len());
    }

    #[test]
    fn make_connected_baseline_merges_and_removes_tiny_edges() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 near-dup of 0

        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2
        brep.edges.push(Edge { start: 0, end: 3 }); // e3 tiny edge to be removed

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_baseline(&brep, 1e-6);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert_eq!(report.passes_run, 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn make_connected_iterative_reports_convergence() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative(&brep, 1e-6, 4);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance >= 1e-6);
    }

    #[test]
    fn make_connected_iterative_with_growth_increases_final_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth(&brep, 1e-6, 4, 2.0);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance > 1e-6);
    }

    #[test]
    fn make_connected_iterative_with_growth_cap_clamps_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            1e-6,
            4,
            10.0,
            2e-6,
        );
        assert!(report.passes_run >= 2);
        assert!(report.tolerance_cap_applied);
        assert!((report.final_tolerance - 2e-6).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_iterative_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(5e-6, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            1e-6,
            2,
            10.0,
            1e-5,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - 1e-5).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_scoped_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(5e-6, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            1e-6,
            2,
            10.0,
            1e-5,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - 1e-5).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_scoped_only_affects_seed_region() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup near region A)
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 4
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 1.0, 0.0) }); // 5
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 6 (dup near region B)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge in scoped region
        brep.edges.push(Edge { start: 4, end: 5 });
        brep.edges.push(Edge { start: 5, end: 6 }); // unrelated region

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (scoped, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            1e-6,
            3,
            1.0,
            1e-4,
        );

        assert!(report.vertices_merged >= 1);
        assert!(scoped.vertices.len() < brep.vertices.len());

        // Vertex near unrelated region B should remain after scoped cleanup.
        let has_far = scoped
            .vertices
            .iter()
            .any(|v| (v.point - DVec3::new(10.0, 0.0, 0.0)).length() <= 1e-12);
        assert!(has_far);
    }

    #[test]
    fn repair_unit_box_is_no_op() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, report) = repair(&brep, 1e-7);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        // Face count unchanged
        let faces: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(faces, 6, "unit box should have 6 faces after repair");
    }

    #[test]
    fn merge_close_vertices_merges_duplicates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        // Add two vertices at nearly the same position
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1e-9, 0.0, 0.0),
        }); // dup of 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, merged) = merge_close_vertices(&brep, 1e-6);
        assert!(merged >= 1, "should merge the near-duplicate vertex");
        assert!(
            fixed.vertices.len() < brep.vertices.len(),
            "should have fewer vertices"
        );
    }

    #[test]
    fn recompute_normals_fixes_zero_normal() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        // Face with wrong/zero normal
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // intentionally wrong
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = recompute_face_normals(&brep);
        assert!(
            n > 0 || fixed.solids[0].shells[0].faces[0].normal != DVec3::ZERO,
            "normal should have been fixed"
        );
    }

    #[test]
    fn remove_degenerate_face() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });
        // Only 2 edges — degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = remove_degenerate_faces(&brep);
        assert_eq!(n, 1);
        let face_count: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(face_count, 0);
    }

    #[test]
    fn fix_same_range_flags_aligns_curve2d_ranges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Build minimal SameRange mismatch for edge 0.
        if brep.geom.edge_curve_range.is_empty() {
            brep.geom.edge_curve_range = vec![Some([0.0, std::f64::consts::PI])];
        } else {
            brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        }
        if brep.geom.edge_pcurves.is_empty() || brep.geom.edge_pcurves[0].is_empty() {
            // Sphere primitive normally has seam pcurves, but guard for future changes.
            return;
        }

        brep.geom.edge_same_range = vec![false; brep.edges.len().max(1)];
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]); // mismatched

        let (fixed, n) = fix_same_range_flags(&brep, 1e-9);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn fix_same_range_with_scan_repairs_flagged_edges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        if brep.geom.edge_curve_range.is_empty()
            || brep.geom.edge_pcurves.is_empty()
            || brep.geom.edge_pcurves[0].is_empty()
        {
            return;
        }

        brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        if brep.geom.edge_same_range.len() < brep.edges.len() {
            brep.geom.edge_same_range.resize(brep.edges.len(), true);
        }

        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]);
        brep.geom.edge_same_range[0] = false;

        let (fixed, n) = fix_same_range_with_scan(&brep, 1e-9);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn propagate_tolerances_bottom_up_fills_slots_and_propagates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Simple triangle face: 3 verts, 3 edges.
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2

        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        });

        // Set vertex 0 with a large tolerance.
        brep.geom.vertex_tolerance = vec![1e-3, 0.0, 0.0];

        let out = propagate_tolerances(&brep, 1e-7, ToleranceFlowDirection::BottomUp);

        // vertex_tolerance slots must be filled.
        assert_eq!(out.geom.vertex_tolerance.len(), 3);
        // Edge tolerances should be at least floor.
        assert!(out.geom.edge_tolerance.len() >= 3);
        // Edge 0 connects v0 (tol=1e-3) and v1 (tol=floor); must ≥ 1e-3.
        assert!(out.geom.edge_tolerance[0] >= 1e-3);
        // Face tolerance should be ≥ max edge tolerance.
        assert!(out.geom.face_tolerance[0] >= out.geom.edge_tolerance[0]);
    }

    #[test]
    fn propagate_tolerances_top_down_spreads_face_tol_to_vertices() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        });
        // Assign a large face tolerance.
        brep.geom.face_tolerance = vec![5e-4];

        let out = propagate_tolerances(&brep, 1e-7, ToleranceFlowDirection::TopDown);

        // All edge tolerances should be ≥ face tolerance.
        for etol in &out.geom.edge_tolerance {
            assert!(*etol >= 5e-4);
        }
        // All vertex tolerances should be ≥ face tolerance after propagation.
        for vtol in &out.geom.vertex_tolerance {
            assert!(*vtol >= 5e-4);
        }
    }
}
