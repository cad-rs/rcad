//! BRep validity checker.
//!
//! Analogous to OCCT `BRepCheck_Analyzer`. Checks structural and geometric
//! consistency of a BRep without modifying it.
//!
//! # Checks performed
//!
//! - **C1 Wire closure**: every wire must form a closed chain — the end vertex of
//!   each edge must equal the start vertex of the next edge.
//! - **C2 Face normal consistency**: each face's stored normal must not be a zero
//!   vector.
//! - **C3 Degenerate face**: faces with fewer than 3 wire edges are degenerate.
//! - **C4 Edge index validity**: WireEdge indices must be within bounds of
//!   `brep.edges`.
//! - **C5 Vertex index validity**: each edge's start/end indices must be within
//!   bounds of `brep.vertices`.
//! - **C6 Manifold topology**: each edge must be shared by exactly 2 faces
//!   (for closed manifold solids).
//! - **C7 Wire self-intersection**: a wire's edges must not share vertices
//!   except at consecutive junctions (no figure-8 or self-touching wires).

use glam::DVec3;
use rcad_kernel::BRep;

/// A single validity issue found during checking.
#[derive(Debug, Clone, PartialEq)]
pub enum CheckIssue {
    /// Wire is not closed: end vertex of edge `edge_idx` does not match start
    /// vertex of the next edge in the wire (solid `solid`, shell `shell`,
    /// face `face`, position `wire_pos`).
    OpenWire {
        solid: usize,
        shell: usize,
        face: usize,
        /// Index of the edge within the wire where the gap occurs.
        wire_pos: usize,
    },
    /// Face normal is a zero vector.
    ZeroNormal {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// Face outer wire has fewer than 3 edges.
    DegenerateFace {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// A WireEdge references an edge index that is out of bounds.
    InvalidEdgeIndex {
        solid: usize,
        shell: usize,
        face: usize,
        edge_idx: usize,
    },
    /// An edge references a vertex index that is out of bounds.
    InvalidVertexIndex { edge: usize, vertex_idx: usize },
    /// An edge is shared by more or fewer than 2 faces (non-manifold).
    NonManifoldEdge { edge_idx: usize, face_count: usize },
    /// A wire has self-intersecting topology: a vertex appears more than
    /// twice (once as start, once as end) in the same wire.
    SelfIntersectingWire {
        solid: usize,
        shell: usize,
        face: usize,
        wire_idx: usize,
        vertex: usize,
    },
    /// A wire's outer boundary edges intersect each other geometrically
    /// (non-adjacent edges in the wire cross in 3D space).
    ///
    /// This catches cases where a face wire forms a figure-eight or butterfly
    /// polygon rather than a simple closed loop.
    GeometricSelfIntersection {
        solid: usize,
        shell: usize,
        face: usize,
        /// Index of one of the crossing edges within the outer wire.
        edge_a: usize,
        /// Index of the other crossing edge within the outer wire.
        edge_b: usize,
    },
}

impl std::fmt::Display for CheckIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CheckIssue::OpenWire {
                solid,
                shell,
                face,
                wire_pos,
            } => write!(
                f,
                "OpenWire: solid={solid} shell={shell} face={face} at wire pos {wire_pos}"
            ),
            CheckIssue::ZeroNormal { solid, shell, face } => {
                write!(f, "ZeroNormal: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::DegenerateFace { solid, shell, face } => {
                write!(f, "DegenerateFace: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::InvalidEdgeIndex {
                solid,
                shell,
                face,
                edge_idx,
            } => write!(
                f,
                "InvalidEdgeIndex: solid={solid} shell={shell} face={face} edge={edge_idx}"
            ),
            CheckIssue::InvalidVertexIndex { edge, vertex_idx } => {
                write!(f, "InvalidVertexIndex: edge={edge} vertex={vertex_idx}")
            }
            CheckIssue::NonManifoldEdge { edge_idx, face_count } => {
                write!(f, "NonManifoldEdge: edge={edge_idx} shared by {face_count} faces (expected 2)")
            }
            CheckIssue::SelfIntersectingWire {
                solid,
                shell,
                face,
                wire_idx,
                vertex,
            } => write!(
                f,
                "SelfIntersectingWire: solid={solid} shell={shell} face={face} wire={wire_idx} vertex={vertex}"
            ),
            CheckIssue::GeometricSelfIntersection {
                solid,
                shell,
                face,
                edge_a,
                edge_b,
            } => write!(
                f,
                "GeometricSelfIntersection: solid={solid} shell={shell} face={face} edges {edge_a} and {edge_b} cross"
            ),
        }
    }
}

/// Result of a BRep validity check.
#[derive(Debug, Clone)]
pub struct CheckResult {
    pub issues: Vec<CheckIssue>,
}

impl CheckResult {
    /// Returns `true` if no issues were found.
    pub fn is_valid(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Check the validity of a BRep and return a `CheckResult` with any issues found.
///
/// Analogous to OCCT `BRepCheck_Analyzer::Perform()`.
pub fn check(brep: &BRep) -> CheckResult {
    let mut issues = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    // C5: edge vertex bounds
    for (eidx, edge) in brep.edges.iter().enumerate() {
        if edge.start >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.start,
            });
        }
        if edge.end >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.end,
            });
        }
    }

    // C6: manifold check — each edge must be shared by exactly 2 faces.
    // Count how many faces reference each edge across all solids/shells/faces.
    let mut edge_face_count: Vec<usize> = vec![0; n_edges];
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Count edges in outer wire
                for we in &face.outer_wire.edges {
                    if we.idx < n_edges {
                        edge_face_count[we.idx] += 1;
                    }
                }
                // Count edges in inner wires
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < n_edges {
                            edge_face_count[we.idx] += 1;
                        }
                    }
                }
            }
        }
    }
    for (eidx, &count) in edge_face_count.iter().enumerate() {
        if count != 2 {
            issues.push(CheckIssue::NonManifoldEdge {
                edge_idx: eidx,
                face_count: count,
            });
        }
    }

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                let wire = &face.outer_wire;

                // C2: zero normal
                if face.normal == DVec3::ZERO {
                    issues.push(CheckIssue::ZeroNormal {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                }

                // C3: degenerate face
                if wire.edges.len() < 3 {
                    issues.push(CheckIssue::DegenerateFace {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                    // Can't check wire closure for degenerate face
                    continue;
                }

                // C4: edge index bounds + collect start/end vertices for wire closure check
                let mut valid = true;
                let mut wire_verts: Vec<(usize, usize)> = Vec::new(); // (start_vidx, end_vidx)
                for we in &wire.edges {
                    if we.idx >= n_edges {
                        issues.push(CheckIssue::InvalidEdgeIndex {
                            solid: si,
                            shell: shi,
                            face: fi,
                            edge_idx: we.idx,
                        });
                        valid = false;
                    } else {
                        let edge = &brep.edges[we.idx];
                        let (sv, ev) = if we.forward {
                            (edge.start, edge.end)
                        } else {
                            (edge.end, edge.start)
                        };
                        wire_verts.push((sv, ev));
                    }
                }

                if !valid {
                    continue;
                }

                // C1: wire closure — end of edge[i] must match start of edge[i+1]
                let n = wire_verts.len();
                for i in 0..n {
                    let next = (i + 1) % n;
                    let end_v = wire_verts[i].1;
                    let start_v = wire_verts[next].0;
                    if end_v != start_v {
                        // Tolerance check: allow same position even if different vertex objects
                        let end_pt = brep.vertices[end_v].point;
                        let start_pt = brep.vertices[start_v].point;
                        if (end_pt - start_pt).length() > 1e-6 {
                            issues.push(CheckIssue::OpenWire {
                                solid: si,
                                shell: shi,
                                face: fi,
                                wire_pos: i,
                            });
                        }
                    }
                }

                // C7: wire self-intersection — each vertex should appear at most
                // twice in the wire (once as start of an edge, once as end of another).
                check_wire_self_intersection(
                    &wire_verts,
                    &brep.vertices,
                    si, shi, fi, 0, // outer wire index = 0
                    &mut issues,
                );

                // C8: geometric self-intersection — check if non-adjacent edges of
                // the outer wire cross each other in 3D space (projects to 2D via
                // the face plane for planar faces).
                check_geometric_self_intersection(
                    &wire_verts,
                    &brep.vertices,
                    si, shi, fi,
                    &mut issues,
                );

                // Check inner wires too
                for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
                    if inner_wire.edges.len() < 2 {
                        continue; // too few edges to self-intersect
                    }
                    let mut inner_verts: Vec<(usize, usize)> = Vec::new();
                    let mut inner_valid = true;
                    for we in &inner_wire.edges {
                        if we.idx >= n_edges {
                            issues.push(CheckIssue::InvalidEdgeIndex {
                                solid: si,
                                shell: shi,
                                face: fi,
                                edge_idx: we.idx,
                            });
                            inner_valid = false;
                        } else {
                            let edge = &brep.edges[we.idx];
                            let (sv, ev) = if we.forward {
                                (edge.start, edge.end)
                            } else {
                                (edge.end, edge.start)
                            };
                            inner_verts.push((sv, ev));
                        }
                    }
                    if !inner_valid {
                        continue;
                    }

                    // Inner wire closure check
                    let n_inner = inner_verts.len();
                    for i in 0..n_inner {
                        let next = (i + 1) % n_inner;
                        let end_v = inner_verts[i].1;
                        let start_v = inner_verts[next].0;
                        if end_v != start_v {
                            let end_pt = brep.vertices[end_v].point;
                            let start_pt = brep.vertices[start_v].point;
                            if (end_pt - start_pt).length() > 1e-6 {
                                issues.push(CheckIssue::OpenWire {
                                    solid: si,
                                    shell: shi,
                                    face: fi,
                                    wire_pos: i,
                                });
                            }
                        }
                    }

                    // Inner wire self-intersection
                    check_wire_self_intersection(
                        &inner_verts,
                        &brep.vertices,
                        si, shi, fi,
                        wi + 1, // inner wire indices start after outer
                        &mut issues,
                    );
                }
            }
        }
    }

    CheckResult { issues }
}

/// Check a single wire for self-intersecting topology.
///
/// A valid wire wire should have each vertex appear at most twice across
/// all edge endpoints: once as the start of some edge and once as the end
/// of another edge. If a vertex appears 3+ times, the wire self-intersects.
fn check_wire_self_intersection(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    wire_idx: usize,
    issues: &mut Vec<CheckIssue>,
) {
    use std::collections::HashMap;
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
    for &(sv, ev) in wire_verts {
        *vertex_count.entry(sv).or_insert(0) += 1;
        *vertex_count.entry(ev).or_insert(0) += 1;
    }
    // In a closed wire, each vertex should appear exactly twice (once as start,
    // once as end). Allow tolerance for vertices at the same position.
    for (&vidx, &count) in &vertex_count {
        if count > 2 {
            // Check if it's actually a geometric self-intersection (different
            // positions) or just the same vertex referenced multiple times
            // (which could be valid for some topologies).
            if vidx < vertices.len() {
                issues.push(CheckIssue::SelfIntersectingWire {
                    solid,
                    shell,
                    face,
                    wire_idx,
                    vertex: vidx,
                });
            }
        }
    }
}

/// Check whether non-adjacent edges of a wire intersect geometrically.
///
/// Projects the wire edge endpoints onto the face's 2D plane (using any two
/// non-collinear edges to form a local basis) and runs 2D segment intersection
/// tests on all non-adjacent edge pairs. Adjacent edges share an endpoint and
/// therefore trivially "intersect" at that endpoint — they are excluded.
///
/// This check only tests the start/end vertices; curved edges (circles, BSplines)
/// are approximated by their chord.
fn check_geometric_self_intersection(
    wire_verts: &[(usize, usize)],
    vertices: &[rcad_kernel::topology::Vertex],
    solid: usize,
    shell: usize,
    face: usize,
    issues: &mut Vec<CheckIssue>,
) {
    let n = wire_verts.len();
    if n < 4 {
        return; // Need at least 4 edges for non-adjacent crossing to exist
    }

    // Build the 3D segment list: each edge is (p0, p1).
    let segs: Vec<(DVec3, DVec3)> = wire_verts
        .iter()
        .map(|&(sv, ev)| {
            let p0 = vertices.get(sv).map(|v| v.point).unwrap_or(DVec3::ZERO);
            let p1 = vertices.get(ev).map(|v| v.point).unwrap_or(DVec3::ZERO);
            (p0, p1)
        })
        .collect();

    // Project to 2D using a local basis from the first non-degenerate edge.
    let (origin, axis_u, axis_v) = {
        let mut found = None;
        for i in 0..n {
            let d = segs[i].1 - segs[i].0;
            if d.length() > 1e-12 {
                let u = d.normalize();
                // Pick an arbitrary perpendicular
                let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                let v = u.cross(tmp).normalize();
                found = Some((segs[i].0, u, v));
                break;
            }
        }
        match found {
            Some(b) => b,
            None => return, // All edges are degenerate; skip
        }
    };

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(axis_u), d.dot(axis_v)]
    };

    let seg2d: Vec<([f64; 2], [f64; 2])> =
        segs.iter().map(|&(p0, p1)| (project(p0), project(p1))).collect();

    // Check all non-adjacent pairs: i and j are non-adjacent when |i - j| > 1
    // (mod n, so also check the wrap-around adjacency).
    for i in 0..n {
        for j in (i + 2)..n {
            // Skip adjacent pair at the wrap-around (edge n-1 and edge 0 share a vertex)
            if i == 0 && j == n - 1 {
                continue;
            }
            if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
                issues.push(CheckIssue::GeometricSelfIntersection {
                    solid,
                    shell,
                    face,
                    edge_a: i,
                    edge_b: j,
                });
                return; // Report only the first crossing per face
            }
        }
    }
}

/// Returns `true` if the open segment p1→p2 properly intersects segment p3→p4.
/// Returns `false` if they only share an endpoint (T-intersection) or don't cross.
fn segments_2d_properly_intersect(
    p1: [f64; 2],
    p2: [f64; 2],
    p3: [f64; 2],
    p4: [f64; 2],
) -> bool {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];

    let cross = d1[0] * d2[1] - d1[1] * d2[0];

    if cross.abs() < 1e-12 {
        return false; // Parallel or collinear
    }

    let dx = p3[0] - p1[0];
    let dy = p3[1] - p1[1];
    let t = (dx * d2[1] - dy * d2[0]) / cross;
    let s = (dx * d1[1] - dy * d1[0]) / cross;

    // Proper interior intersection: t and s must be strictly in (0, 1)
    let eps = 1e-9;
    t > eps && t < 1.0 - eps && s > eps && s < 1.0 - eps
}

// ── SameParameter diagnosis ───────────────────────────────────────────────────

/// A single edge whose 3D curve endpoints deviate from the vertex positions.
///
/// Analogous to the diagnostic output of `BRepCheck_Edge::SameParameter()` in OCCT.
#[derive(Debug, Clone)]
pub struct SuspectEdge {
    /// Index of the edge in `BRep.edges`.
    pub edge_idx: usize,
    /// Distance from `curve.point_at(t1)` to the start vertex position.
    pub start_gap: f64,
    /// Distance from `curve.point_at(t2)` to the end vertex position.
    pub end_gap: f64,
}

/// Report from [`diagnose_same_parameter`].
#[derive(Debug, Clone, Default)]
pub struct SameParameterDiagnosis {
    /// Edges whose curve endpoints deviate from vertex positions beyond `tolerance`.
    pub suspect_edges: Vec<SuspectEdge>,
}

impl SameParameterDiagnosis {
    /// Returns `true` if no SameParameter violations were found.
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// Scan all edges in `brep` for SameParameter violations.
///
/// For each edge that has a known 3D curve and curve range `[t1, t2]`, evaluates
/// the curve at `t1` and `t2` and compares the resulting points to the
/// corresponding vertex positions.  Edges whose endpoint gap exceeds `tolerance`
/// are reported as suspects and their `edge_same_parameter` flag is updated to
/// `false` so that a subsequent [`fix_same_parameter`] pass can repair them.
///
/// Analogous to `BRepLib::SameParameter()` diagnosis step in OCCT.
///
/// Returns a non-mutating diagnosis; call `fix_same_parameter_with_scan` to
/// also repair the flagged edges.
pub fn diagnose_same_parameter(brep: &BRep, tolerance: f64) -> SameParameterDiagnosis {
    use rcad_kernel::geom::CurveEval;

    let mut suspects = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    for edge_idx in 0..n_edges {
        let curve_idx = brep.geom.edge_curve.get(edge_idx).and_then(|c| *c);
        let Some(ci) = curve_idx else { continue };
        let Some(curve) = brep.geom.curves.get(ci) else { continue };
        let Some(range) = brep.geom.edge_curve_range.get(edge_idx).and_then(|r| *r) else { continue };

        let edge = &brep.edges[edge_idx];
        if edge.start >= n_verts || edge.end >= n_verts {
            continue;
        }

        let p_start = brep.vertices[edge.start].point;
        let p_end = brep.vertices[edge.end].point;

        let eval_start = curve.point_at(range[0]);
        let eval_end = curve.point_at(range[1]);

        let start_gap = (eval_start - p_start).length();
        let end_gap = (eval_end - p_end).length();

        if start_gap > tolerance || end_gap > tolerance {
            suspects.push(SuspectEdge {
                edge_idx,
                start_gap,
                end_gap,
            });
        }
    }

    SameParameterDiagnosis { suspect_edges: suspects }
}

// ── Shell topology analysis ───────────────────────────────────────────────────

/// Topology analysis report for a BRep's shell structure.
///
/// Analogous to `ShapeAnalysis_Shell` in OCCT's TKShHealing.
#[derive(Debug, Clone, Default)]
pub struct ShellTopologyReport {
    /// `true` if every edge is referenced by exactly 2 faces (no free edges).
    pub is_closed: bool,
    /// `true` if no edge is referenced by more than 2 faces.
    pub is_manifold: bool,
    /// Number of edges referenced by exactly 1 face (free / open edges).
    pub open_edge_count: usize,
    /// Number of edges referenced by 3 or more faces (non-manifold edges).
    pub non_manifold_edge_count: usize,
    /// Number of vertices not referenced by any edge in the BRep.
    pub isolated_vertex_count: usize,
    /// Total edge count in the BRep.
    pub total_edges: usize,
    /// Total face count across all solids/shells.
    pub total_faces: usize,
}

/// Analyze the shell topology of `brep`.
///
/// Counts free edges, non-manifold edges, and isolated vertices to determine
/// whether the BRep represents a closed manifold solid.
///
/// Analogous to `ShapeAnalysis_Shell::LoadUnorientedEdges()` + checks in OCCT.
pub fn analyze_shell_topology(brep: &BRep) -> ShellTopologyReport {
    let total_edges = brep.edges.len();

    // Count how many faces each edge is referenced by.
    let mut edge_face_count: Vec<usize> = vec![0; total_edges];
    let mut total_faces = 0usize;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                total_faces += 1;
                for we in &face.outer_wire.edges {
                    if we.idx < total_edges {
                        edge_face_count[we.idx] += 1;
                    }
                }
                for wire in &face.inner_wires {
                    for we in &wire.edges {
                        if we.idx < total_edges {
                            edge_face_count[we.idx] += 1;
                        }
                    }
                }
            }
        }
    }

    let open_edge_count = edge_face_count.iter().filter(|&&c| c == 1).count();
    let non_manifold_edge_count = edge_face_count.iter().filter(|&&c| c > 2).count();

    // Count isolated vertices: those not referenced by any edge as start or end.
    let mut vertex_used = vec![false; brep.vertices.len()];
    for edge in &brep.edges {
        if edge.start < vertex_used.len() {
            vertex_used[edge.start] = true;
        }
        if edge.end < vertex_used.len() {
            vertex_used[edge.end] = true;
        }
    }
    let isolated_vertex_count = vertex_used.iter().filter(|&&used| !used).count();

    ShellTopologyReport {
        is_closed: open_edge_count == 0,
        is_manifold: non_manifold_edge_count == 0,
        open_edge_count,
        non_manifold_edge_count,
        isolated_vertex_count,
        total_edges,
        total_faces,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::geom::{Curve3, Line3};

    #[test]
    fn analyze_shell_topology_unit_box_is_closed_manifold() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let report = analyze_shell_topology(&brep);
        assert!(report.is_closed, "unit box should be a closed shell");
        assert!(report.is_manifold, "unit box should be manifold");
        assert_eq!(report.open_edge_count, 0);
        assert_eq!(report.non_manifold_edge_count, 0);
        assert_eq!(report.total_faces, 6);
    }

    #[test]
    fn diagnose_same_parameter_clean_box_has_no_violations() {
        // A primitive box has no geom curves populated, so the diagnosis should
        // return empty (nothing to check = nothing flagged).
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let diagnosis = diagnose_same_parameter(&brep, 1e-7);
        assert!(
            diagnosis.is_clean(),
            "primitive box with no edge_curve entries should have no violations"
        );
    }

    #[test]
    fn diagnose_same_parameter_detects_mismatch() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Build a triangle with a Line3 curve whose range is deliberately mismatched.
        let mut brep = BRep::new();
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(rcad_kernel::topology::Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(rcad_kernel::topology::Edge { start: 0, end: 1 });
        brep.edges.push(rcad_kernel::topology::Edge { start: 1, end: 2 });
        brep.edges.push(rcad_kernel::topology::Edge { start: 2, end: 0 });

        // Edge 0: line from (0,0,0) toward (1,0,0), but with range [0, 999] — huge mismatch
        let ci = brep.geom.curves.len();
        brep.geom.curves.push(Curve3::Line(Line3 {
            origin: glam::DVec3::ZERO,
            direction: glam::DVec3::X,
        }));
        brep.geom.edge_curve.push(Some(ci));
        brep.geom.edge_curve_range.push(Some([0.0, 999.0])); // wrong range!
        brep.geom.edge_curve.push(None);
        brep.geom.edge_curve_range.push(None);
        brep.geom.edge_curve.push(None);
        brep.geom.edge_curve_range.push(None);

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: glam::DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let diagnosis = diagnose_same_parameter(&brep, 1e-6);
        assert!(!diagnosis.is_clean(), "mismatch should be detected");
        assert_eq!(diagnosis.suspect_edges[0].edge_idx, 0);
        assert!(diagnosis.suspect_edges[0].end_gap > 1.0, "end gap should be ~998");
    }

    #[test]
    fn unit_box_is_valid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check(&brep);
        assert!(
            result.is_valid(),
            "unit box should pass all checks; issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn open_wire_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with a deliberately open wire (gap between edge 1 end and edge 0 start)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        }); // 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        }); // 1
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 1.0, 0.0),
        }); // 2
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        }); // 3 (gap: wire goes 0→1→2 then 2→0 skips 3)

        // Edge 0: v0 → v1; Edge 1: v1 → v2; Edge 2: v2 → v0 (skips v3 — would close)
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // intentional gap: starts at v3 not v2

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2), // e2 starts at v3, but e1 ends at v2 → open
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(!result.is_valid(), "open wire should be detected");
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::OpenWire { .. }))
        );
    }

    #[test]
    fn degenerate_face_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });

        // Face with only 2 edges — degenerate
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

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::DegenerateFace { .. }))
        );
    }

    #[test]
    fn zero_normal_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        for p in [DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Z] {
            brep.vertices.push(Vertex { point: p });
        }
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // zero normal — invalid
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::ZeroNormal { .. })),
            "expected ZeroNormal issue"
        );
    }

    #[test]
    fn invalid_edge_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 }); // only edge 0 exists

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(99), // out-of-bounds
                    WireEdge::fwd(0),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidEdgeIndex { .. })),
            "expected InvalidEdgeIndex issue"
        );
    }

    #[test]
    fn invalid_vertex_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Vertex};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 99 }); // vertex 99 doesn't exist

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidVertexIndex { .. })),
            "expected InvalidVertexIndex issue"
        );
    }

    #[test]
    fn non_manifold_edge_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep where an edge is shared by 3 faces (non-manifold)
        let mut brep = BRep::new();
        // 4 vertices forming a square
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 1.0) });

        // 5 edges: 4 forming a square + 1 vertical
        brep.edges.push(Edge { start: 0, end: 1 }); // e0: bottom
        brep.edges.push(Edge { start: 1, end: 2 }); // e1: right
        brep.edges.push(Edge { start: 2, end: 3 }); // e2: top
        brep.edges.push(Edge { start: 3, end: 0 }); // e3: left
        brep.edges.push(Edge { start: 0, end: 4 }); // e4: vertical

        // 3 faces sharing edge e4 (vertical edge) — non-manifold
        // Face 1: uses e4
        let face1 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(4), WireEdge::fwd(0), WireEdge::rev(3)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        // Face 2: uses e4
        let face2 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::rev(4), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        // Face 3: uses e4 again — this makes e4 shared by 3 faces
        let face3 = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(4), WireEdge::fwd(3), WireEdge::rev(0)],
            },
            inner_wires: vec![],
            normal: DVec3::NEG_Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face1, face2, face3] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::NonManifoldEdge { edge_idx: 4, .. })),
            "expected NonManifoldEdge for edge 4, issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn self_intersecting_wire_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with a figure-8 wire: vertex 0 appears 3 times
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // v0 — center, appears 3x
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // v1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // v2
        brep.vertices.push(Vertex { point: DVec3::new(-1.0, 0.0, 0.0) }); // v3
        brep.vertices.push(Vertex { point: DVec3::new(0.0, -1.0, 0.0) }); // v4

        // Figure-8: v0→v1→v2→v0→v3→v4→v0 (v0 appears 3 times as start/end)
        brep.edges.push(Edge { start: 0, end: 1 }); // e0: v0→v1
        brep.edges.push(Edge { start: 1, end: 2 }); // e1: v1→v2
        brep.edges.push(Edge { start: 2, end: 0 }); // e2: v2→v0
        brep.edges.push(Edge { start: 0, end: 3 }); // e3: v0→v3
        brep.edges.push(Edge { start: 3, end: 4 }); // e4: v3→v4
        brep.edges.push(Edge { start: 4, end: 0 }); // e5: v4→v0

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                    WireEdge::fwd(3),
                    WireEdge::fwd(4),
                    WireEdge::fwd(5),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::SelfIntersectingWire { .. })),
            "expected SelfIntersectingWire issue, issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn inner_wire_open_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with an open inner wire (hole that doesn't close)
        let mut brep = BRep::new();
        // Outer wire: triangle
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(3.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.5, 3.0, 0.0) });
        // Inner wire vertices (don't close: v3→v4→v5, but v5≠v3)
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(2.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.5, 0.5, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2
        // Inner wire edges (open: e3: v3→v4, e4: v4→v5, e5: v5→v3 would close but we skip)
        brep.edges.push(Edge { start: 3, end: 4 }); // e3
        brep.edges.push(Edge { start: 4, end: 5 }); // e4
        // Intentionally missing: edge from v5 back to v3

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![Wire {
                edges: vec![WireEdge::fwd(3), WireEdge::fwd(4)], // open: v5≠v3
            }],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::OpenWire { .. })),
            "expected OpenWire for inner wire, issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn valid_box_passes_all_new_checks() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check(&brep);
        assert!(
            result.is_valid(),
            "unit box should pass all checks including manifold and self-intersection; issues: {:?}",
            result.issues
        );
    }
}
