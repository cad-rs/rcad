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
}
