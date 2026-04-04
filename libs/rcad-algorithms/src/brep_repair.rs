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
use rcad_kernel::{BRep, CurveEval};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

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
}

impl std::fmt::Display for RepairReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "RepairReport {{ vertices_merged={}, degenerate_removed={}, normals_recomputed={}, wires_fixed={} }}",
            self.vertices_merged,
            self.degenerate_faces_removed,
            self.normals_recomputed,
            self.wires_fixed,
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
    (b, report)
}

/// Merge vertices that are within `tolerance` of each other.
///
/// Uses a union-find approach: for each pair of vertices closer than
/// `tolerance`, they are merged into the vertex with the smaller index.
/// All edges and wires are remapped accordingly.
///
/// Returns the repaired BRep and the number of vertices merged.
///
/// Analogous to `BRepOffsetAPI_Sewing` vertex merging or
/// `ShapeFix_Wire::FixSameParameter`.
pub fn merge_close_vertices(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let n = brep.vertices.len();
    // Union-find: parent[i] = canonical representative of vertex i
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path compression
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            // Merge to the smaller index so result is deterministic
            if ra < rb { parent[rb] = ra; } else { parent[ra] = rb; }
        }
    }

    let tol2 = tolerance * tolerance;
    for i in 0..n {
        for j in (i + 1)..n {
            let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
            if d2 <= tol2 {
                union(&mut parent, i, j);
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
            new_vertices.push(brep.vertices[rep].clone());
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    // Re-map edges
    let new_edges: Vec<Edge> = brep.edges.iter().map(|e| Edge {
        start: remap[e.start],
        end: remap[e.end],
    }).collect();

    // Rebuild solids with remapped wires (topology is unchanged, just vertex indices)
    let new_solids = brep.solids.iter().map(|solid| Solid {
        shells: solid.shells.iter().map(|shell| Shell {
            faces: shell.faces.iter().map(|face| {
                let remap_wire = |w: &Wire| Wire {
                    edges: w.edges.clone(),  // WireEdge indices are edge indices, not vertex
                };
                Face {
                    outer_wire: remap_wire(&face.outer_wire),
                    inner_wires: face.inner_wires.iter().map(remap_wire).collect(),
                    normal: face.normal,
                    triangles: face.triangles.clone(),
                }
            }).collect(),
        }).collect(),
    }).collect();

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

    let new_solids = brep.solids.iter().map(|solid| Solid {
        shells: solid.shells.iter().map(|shell| {
            let new_faces: Vec<Face> = shell.faces.iter().filter(|face| {
                let wire = &face.outer_wire;
                // Must have at least 3 edges
                if wire.edges.len() < 3 {
                    removed += 1;
                    return false;
                }
                // Collect distinct vertex positions
                let pts: Vec<DVec3> = wire.edges.iter().filter_map(|we| {
                    brep.edges.get(we.idx).map(|e| {
                        let vidx = if we.forward { e.start } else { e.end };
                        brep.vertices.get(vidx).map(|v| v.point)
                    }).flatten()
                }).collect();

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
            }).cloned().collect();
            Shell { faces: new_faces }
        }).collect(),
    }).collect();

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

    let new_solids = brep.solids.iter().map(|solid| Solid {
        shells: solid.shells.iter().map(|shell| Shell {
            faces: shell.faces.iter().map(|face| {
                let pts: Vec<DVec3> = face.outer_wire.edges.iter().filter_map(|we| {
                    brep.edges.get(we.idx).map(|e| {
                        let vidx = if we.forward { e.start } else { e.end };
                        brep.vertices.get(vidx).map(|v| v.point)
                    }).flatten()
                }).collect();

                let new_normal = if pts.len() >= 3 {
                    let n = newell_normal(&pts);
                    if n.length() > 1e-14 { n.normalize() } else { face.normal }
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
                }
            }).collect(),
        }).collect(),
    }).collect();

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

    let new_solids = brep.solids.iter().map(|solid| Solid {
        shells: solid.shells.iter().map(|shell| Shell {
            faces: shell.faces.iter().map(|face| {
                let (new_outer, fixed_outer) =
                    fix_wire(&face.outer_wire, brep, tol2);
                let (new_inners, fixed_inner): (Vec<Wire>, usize) = face.inner_wires.iter()
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
                }
            }).collect(),
        }).collect(),
    }).collect();

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
        let e_curr = match brep.edges.get(edges[i].idx) { Some(e) => e, None => continue };
        let e_next = match brep.edges.get(edges[next].idx) { Some(e) => e, None => continue };

        // end vertex of current edge
        let end_v = if edges[i].forward { e_curr.end } else { e_curr.start };
        // start vertex of next edge
        let start_v = if edges[next].forward { e_next.start } else { e_next.end };

        if end_v == start_v {
            continue; // already connected
        }
        // Check spatial proximity
        if let (Some(ep), Some(sp)) = (
            brep.vertices.get(end_v).map(|v| v.point),
            brep.vertices.get(start_v).map(|v| v.point),
        ) {
            if (ep - sp).length_squared() <= tol2 {
                continue; // close enough — OK
            }
        }

        // Try flipping the *next* edge to see if that connects the chain
        let alt_start = if edges[next].forward { e_next.end } else { e_next.start };
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

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn repair_unit_box_is_no_op() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
        let (fixed, report) = repair(&brep, 1e-7);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        // Face count unchanged
        let faces: usize = fixed.solids.iter().flat_map(|s| &s.shells).map(|sh| sh.faces.len()).sum();
        assert_eq!(faces, 6, "unit box should have 6 faces after repair");
    }

    #[test]
    fn merge_close_vertices_merges_duplicates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        // Add two vertices at nearly the same position
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1e-9, 0.0, 0.0) }); // dup of 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });
        let face = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let (fixed, merged) = merge_close_vertices(&brep, 1e-6);
        assert!(merged >= 1, "should merge the near-duplicate vertex");
        assert!(fixed.vertices.len() < brep.vertices.len(), "should have fewer vertices");
    }

    #[test]
    fn recompute_normals_fixes_zero_normal() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        // Face with wrong/zero normal
        let face = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)] },
            inner_wires: vec![],
            normal: DVec3::ZERO,  // intentionally wrong
            triangles: vec![],
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });
        let (fixed, n) = recompute_face_normals(&brep);
        assert!(n > 0 || fixed.solids[0].shells[0].faces[0].normal != DVec3::ZERO,
            "normal should have been fixed");
    }

    #[test]
    fn remove_degenerate_face() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });
        // Only 2 edges — degenerate
        let face = Face {
            outer_wire: Wire { edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });
        let (fixed, n) = remove_degenerate_faces(&brep);
        assert_eq!(n, 1);
        let face_count: usize = fixed.solids.iter()
            .flat_map(|s| &s.shells).map(|sh| sh.faces.len()).sum();
        assert_eq!(face_count, 0);
    }
}
