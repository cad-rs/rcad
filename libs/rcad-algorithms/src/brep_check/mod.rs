//! BRep validity checker.
//!
//! Analogous to OCCT `BRepCheck_Analyzer`. Checks structural and geometric
//! consistency of a BRep without modifying it.
//!
//! # Checks performed
//!
//! - **C1 Wire closure**: every wire must form a closed chain — the end vertex of
//! each edge must equal the start vertex of the next edge.
//! - **C2 Face normal consistency**: each face's stored normal must not be a zero
//! vector.
//! - **C3 Degenerate face**: faces with fewer than 3 wire edges are degenerate.
//! - **C4 Edge index validity**: WireEdge indices must be within bounds of
//! `brep.edges`.
//! - **C5 Vertex index validity**: each edge's start/end indices must be within
//! bounds of `brep.vertices`.
//! - **C6 Manifold topology**: each edge must be shared by exactly 2 faces
//! (for closed manifold solids).
//! - **C7 Wire self-intersection**: a wire's edges must not share vertices
//! except at consecutive junctions (no figure-8 or self-touching wires).
//!
//! # Extended checks (OCCT BRepCheck_Analyzer equivalent)
//!
//! - **Surface continuity**: C0, C1, C2 continuity across adjacent faces
//! - **Curve-surface consistency**: 3D curve endpoints match surface evaluation
//! - **Edge-curve tolerance verification**: edge tolerance covers geometry deviation
//! - **Face-surface tolerance verification**: face tolerance covers surface deviation
//! - **Shell orientation consistency**: consistent normal orientation in shells
//! - **Solid closure verification**: all edges shared by exactly 2 faces
//! - **Wire orientation**: clockwise vs counter-clockwise validation
//! - **Nested wire validation**: inner loops properly contained within outer
//! - **Tolerance consistency**: adjacent faces have compatible tolerances
//! - **Vertex tolerance propagation**: vertices have appropriate tolerances
//! - **Aspect ratio checks**: face quality metrics
//! - **Degenerate geometry detection**: zero-length edges, collapsed faces
//! - **Sliver face detection**: very thin triangular faces
//! - **Small feature detection**: tiny faces, edges, vertices

use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::PCurve;
use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
use rcad_kernel::topods::{self, ShapeRef, TShape};

// ---- Backward-compat BRep navigation helpers ----
// These let analysis functions navigate topods::BRep using flat indices.

/// Find the ShapeRef for the n-th TShape::Solid (0-based).
fn ns_solid(brep: &rcad_kernel::BRep, n: usize) -> Option<ShapeRef> {
    let mut count = 0usize;
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if matches!(ts.as_ref(), TShape::Solid(_)) {
            if count == n {
                return Some(ShapeRef::synthetic(i));
            }
            count += 1;
        }
    }
    None
}

/// Get a reference to the n-th TShape::Solid's data.
fn ns_solid_data<'a>(brep: &'a rcad_kernel::BRep, n: usize) -> Option<&'a topods::TSolidData> {
    let sr = ns_solid(brep, n)?;
    match &*brep.tshapes[sr.index] {
        TShape::Solid(sd) => Some(sd),
        _ => None,
    }
}

/// Get shell data for a (solid_idx, shell_idx) pair.
fn ns_shell_data<'a>(
    brep: &'a rcad_kernel::BRep,
    solid_idx: usize,
    shell_idx: usize,
) -> Option<&'a topods::TShellData> {
    let sd = ns_solid_data(brep, solid_idx)?;
    let sh_sr = sd.shells.get(shell_idx)?;
    match &*brep.tshapes[sh_sr.index] {
        TShape::Shell(shd) => Some(shd),
        _ => None,
    }
}

/// Get face data for a (solid_idx, shell_idx, face_idx) triple.
fn ns_face_data<'a>(
    brep: &'a rcad_kernel::BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
) -> Option<&'a topods::TFaceData> {
    let shd = ns_shell_data(brep, solid_idx, shell_idx)?;
    let f_sr = shd.faces.get(face_idx)?;
    match &*brep.tshapes[f_sr.index] {
        TShape::Face(fd) => Some(fd),
        _ => None,
    }
}

/// Get wire data for a (solid_idx, shell_idx, face_idx, wire_idx) wire.
/// wire_idx = None -> outer wire, Some(i) -> inner wire.
fn ns_wire_data<'a>(
    brep: &'a rcad_kernel::BRep,
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    wire_idx: Option<usize>,
) -> Option<&'a topods::TWireData> {
    let fd = ns_face_data(brep, solid_idx, shell_idx, face_idx)?;
    let wr = match wire_idx {
        None => fd.outer_wire,
        Some(i) => *fd.inner_wires.get(i)?,
    };
    match &*brep.tshapes[wr.index] {
        TShape::Wire(wd) => Some(wd),
        _ => None,
    }
}

/// Get edge data by flat edge index (position in tshapes).
fn e_edge_data<'a>(brep: &'a rcad_kernel::BRep, edge_idx: usize) -> Option<&'a topods::TEdgeData> {
    let ts = brep.tshapes.get(edge_idx)?;
    match &**ts {
        TShape::Edge(ed) => Some(ed),
        _ => None,
    }
}

/// Look up the 3D curve for an edge.
fn edge_curve9(brep: &rcad_kernel::BRep, edge_idx: usize) -> Option<&rcad_kernel::geom::Curve3> {
    e_edge_data(brep, edge_idx)?.curve.as_ref()
}

/// Look up the parameter range for an edge's 3D curve.
fn edge_range9(brep: &rcad_kernel::BRep, edge_idx: usize) -> Option<[f64; 2]> {
    Some(e_edge_data(brep, edge_idx)?.range)
}

/// Check if an edge is degenerate.
fn edge_degenerated9(brep: &rcad_kernel::BRep, edge_idx: usize) -> bool {
    e_edge_data(brep, edge_idx)
        .map(|ed| ed.degenerated)
        .unwrap_or(false)
}

/// Find the tshape index of a face given its flat (solid, shell, face) coordinates.

// ---- Per-wire helpers for topology checking ----

/// Get vertex pairs for a wire's edges (oriented).
fn wire_vertex_pairs(brep: &rcad_kernel::BRep, wd: &topods::TWireData) -> Vec<(usize, usize)> {
    let mut verts = Vec::with_capacity(wd.edges.len());
    for wesr in &wd.edges {
        if let Some(ed) = e_edge_data(brep, wesr.index) {
            let (sv, ev) = if wesr.orientation.is_forward() {
                (ed.first.index, ed.last.index)
            } else {
                (ed.last.index, ed.first.index)
            };
            verts.push((sv, ev));
        }
    }
    verts
}

// ===========================================================?
// CheckIssue, CheckResult, and the main entry point
// ===========================================================?

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
    //    Geometry validation issues
    /// Surface continuity violation between adjacent faces.
    SurfaceContinuityViolation {
        solid: usize,
        face_a: usize,
        face_b: usize,
        shared_edge: usize,
        /// Expected continuity (0=C0, 1=C1, 2=C2)
        expected: u8,
        /// Actual continuity achieved
        actual: u8,
        /// Gap or angle deviation at the junction
        deviation: f64,
    },
    /// Curve-surface consistency violation: 3D curve doesn't match surface evaluation.
    CurveSurfaceMismatch {
        edge: usize,
        surface: usize,
        /// Maximum deviation between 3D curve and surface curve
        max_deviation: f64,
    },
    /// Edge tolerance insufficient to cover geometry deviation.
    EdgeToleranceViolation {
        edge: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    /// Face tolerance insufficient to cover surface deviation.
    FaceToleranceViolation {
        solid: usize,
        shell: usize,
        face: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    //    Topology validation issues
    /// Shell has inconsistent orientation (mixed inward/outward normals).
    ShellOrientationInconsistent {
        solid: usize,
        shell: usize,
        faces_with_inverted_normals: usize,
    },
    /// Solid is not closed (has boundary edges).
    SolidNotClosed {
        solid: usize,
        boundary_edge_count: usize,
    },
    /// Wire orientation is incorrect for its role (outer vs inner).
    WireOrientationIncorrect {
        solid: usize,
        shell: usize,
        face: usize,
        wire_idx: usize,
        /// true = should be CCW (outer), false = should be CW (inner)
        expected_ccw: bool,
        actual_ccw: bool,
    },
    /// Inner wire is not properly contained within outer wire.
    NestedWireViolation {
        solid: usize,
        shell: usize,
        face: usize,
        inner_wire_idx: usize,
        /// Number of inner wire vertices outside outer wire boundary
        vertices_outside: usize,
    },
    //    Tolerance issues
    /// Adjacent faces have inconsistent tolerances.
    ToleranceInconsistency {
        edge: usize,
        face_a: usize,
        face_b: usize,
        tolerance_a: f64,
        tolerance_b: f64,
        ratio: f64,
    },
    /// Vertex tolerance doesn't cover incident edge endpoints.
    VertexToleranceViolation {
        vertex: usize,
        stored_tolerance: f64,
        required_tolerance: f64,
    },
    //    Quality metric issues
    /// Face has poor aspect ratio.
    PoorAspectRatio {
        solid: usize,
        shell: usize,
        face: usize,
        aspect_ratio: f64,
    },
    /// Edge has near-zero length.
    DegenerateEdge { edge: usize, length: f64 },
    /// Face is a sliver (very thin).
    SliverFace {
        solid: usize,
        shell: usize,
        face: usize,
        area: f64,
        min_dimension: f64,
    },
    /// Small feature detected (tiny face or edge).
    SmallFeature {
        solid: usize,
        shell: usize,
        face: usize,
        feature_type: SmallFeatureType,
        size: f64,
    },
}

/// Type of small feature detected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SmallFeatureType {
    TinyFace,
    TinyEdge,
    TinyVertexGap,
}

impl std::fmt::Display for CheckIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CheckIssue::OpenWire {
                solid,
                shell,
                face,
                wire_pos,
            } => {
                write!(
                    f,
                    "OpenWire: solid={solid} shell={shell} face={face} at wire pos {wire_pos}"
                )
            }
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
            } => {
                write!(
                    f,
                    "InvalidEdgeIndex: solid={solid} shell={shell} face={face} edge={edge_idx}"
                )
            }
            CheckIssue::InvalidVertexIndex { edge, vertex_idx } => {
                write!(f, "InvalidVertexIndex: edge={edge} vertex={vertex_idx}")
            }
            CheckIssue::NonManifoldEdge {
                edge_idx,
                face_count,
            } => {
                write!(
                    f,
                    "NonManifoldEdge: edge={edge_idx} shared by {face_count} faces (expected 2)"
                )
            }
            CheckIssue::SelfIntersectingWire {
                solid,
                shell,
                face,
                wire_idx,
                vertex,
            } => {
                write!(
                    f,
                    "SelfIntersectingWire: solid={solid} shell={shell} face={face} wire={wire_idx} vertex={vertex}"
                )
            }
            CheckIssue::GeometricSelfIntersection {
                solid,
                shell,
                face,
                edge_a,
                edge_b,
            } => {
                write!(
                    f,
                    "GeometricSelfIntersection: solid={solid} shell={shell} face={face} edges {edge_a} and {edge_b} cross"
                )
            }
            CheckIssue::SurfaceContinuityViolation {
                solid,
                face_a,
                face_b,
                shared_edge,
                expected,
                actual,
                deviation,
            } => {
                write!(
                    f,
                    "SurfaceContinuityViolation: solid={solid} faces {face_a}/{face_b} edge={shared_edge} expected C{expected} got C{actual} deviation={deviation:.6e}"
                )
            }
            CheckIssue::CurveSurfaceMismatch {
                edge,
                surface,
                max_deviation,
            } => {
                write!(
                    f,
                    "CurveSurfaceMismatch: edge={edge} surface={surface} deviation={max_deviation:.6e}"
                )
            }
            CheckIssue::EdgeToleranceViolation {
                edge,
                stored_tolerance,
                required_tolerance,
            } => {
                write!(
                    f,
                    "EdgeToleranceViolation: edge={edge} stored={stored_tolerance:.6e} required={required_tolerance:.6e}"
                )
            }
            CheckIssue::FaceToleranceViolation {
                solid,
                shell,
                face,
                stored_tolerance,
                required_tolerance,
            } => {
                write!(
                    f,
                    "FaceToleranceViolation: solid={solid} shell={shell} face={face} stored={stored_tolerance:.6e} required={required_tolerance:.6e}"
                )
            }
            CheckIssue::ShellOrientationInconsistent {
                solid,
                shell,
                faces_with_inverted_normals,
            } => {
                write!(
                    f,
                    "ShellOrientationInconsistent: solid={solid} shell={shell} {faces_with_inverted_normals} inverted faces"
                )
            }
            CheckIssue::SolidNotClosed {
                solid,
                boundary_edge_count,
            } => {
                write!(
                    f,
                    "SolidNotClosed: solid={solid} {boundary_edge_count} boundary edges"
                )
            }
            CheckIssue::WireOrientationIncorrect {
                solid,
                shell,
                face,
                wire_idx,
                expected_ccw,
                actual_ccw,
            } => {
                let expected = if *expected_ccw { "CCW" } else { "CW" };
                let actual = if *actual_ccw { "CCW" } else { "CW" };
                write!(
                    f,
                    "WireOrientationIncorrect: solid={solid} shell={shell} face={face} wire={wire_idx} expected={expected} got={actual}"
                )
            }
            CheckIssue::NestedWireViolation {
                solid,
                shell,
                face,
                inner_wire_idx,
                vertices_outside,
            } => {
                write!(
                    f,
                    "NestedWireViolation: solid={solid} shell={shell} face={face} inner_wire={inner_wire_idx} {vertices_outside} vertices outside"
                )
            }
            CheckIssue::ToleranceInconsistency {
                edge,
                face_a,
                face_b,
                tolerance_a,
                tolerance_b,
                ratio,
            } => {
                write!(
                    f,
                    "ToleranceInconsistency: edge={edge} faces {face_a}/{face_b} tol_a={tolerance_a:.6e} tol_b={tolerance_b:.6e} ratio={ratio:.2}"
                )
            }
            CheckIssue::VertexToleranceViolation {
                vertex,
                stored_tolerance,
                required_tolerance,
            } => {
                write!(
                    f,
                    "VertexToleranceViolation: vertex={vertex} stored={stored_tolerance:.6e} required={required_tolerance:.6e}"
                )
            }
            CheckIssue::PoorAspectRatio {
                solid,
                shell,
                face,
                aspect_ratio,
            } => {
                write!(
                    f,
                    "PoorAspectRatio: solid={solid} shell={shell} face={face} ratio={aspect_ratio:.2}"
                )
            }
            CheckIssue::DegenerateEdge { edge, length } => {
                write!(f, "DegenerateEdge: edge={edge} length={length:.6e}")
            }
            CheckIssue::SliverFace {
                solid,
                shell,
                face,
                area,
                min_dimension,
            } => {
                write!(
                    f,
                    "SliverFace: solid={solid} shell={shell} face={face} area={area:.6e} min_dim={min_dimension:.6e}"
                )
            }
            CheckIssue::SmallFeature {
                solid,
                shell,
                face,
                feature_type,
                size,
            } => {
                let type_str = match feature_type {
                    SmallFeatureType::TinyFace => "TinyFace",
                    SmallFeatureType::TinyEdge => "TinyEdge",
                    SmallFeatureType::TinyVertexGap => "TinyVertexGap",
                };
                write!(
                    f,
                    "SmallFeature: solid={solid} shell={shell} face={face} type={type_str} size={size:.6e}"
                )
            }
        }
    }
}

/// Result of a BRep validity check.
#[derive(Debug, Clone, Default)]
pub struct CheckResult {
    pub issues: Vec<CheckIssue>,
}

impl CheckResult {
    pub fn is_valid(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Check the validity of a BRep and return a `CheckResult` with any issues found.
///
/// Analogous to OCCT `BRepCheck_Analyzer::Perform()`.
pub fn brep_check_analyze(brep: &rcad_kernel::BRep) -> CheckResult {
    let mut issues = Vec::new();
    let n_tshapes = brep.tshapes.len();

    // C5: edge vertex bounds
    for (eidx, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = &**ts else { continue };
        if ed.first.index >= n_tshapes {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: ed.first.index,
            });
        }
        if ed.last.index >= n_tshapes {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: ed.last.index,
            });
        }
    }

    // C6: manifold check — each edge must be shared by exactly 2 faces.
    let mut edge_face_count: Vec<usize> = vec![0; n_tshapes];
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    continue;
                };
                for wesr in &owd.edges {
                    if wesr.index < n_tshapes {
                        edge_face_count[wesr.index] += 1;
                    }
                }
                for iw_sr in &fd.inner_wires {
                    let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
                        continue;
                    };
                    for wesr in &iwd.edges {
                        if wesr.index < n_tshapes {
                            edge_face_count[wesr.index] += 1;
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

    // Per-face checks
    let mut si = 0usize;
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        let mut shi = 0usize;
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            let mut fi = 0usize;
            for face_sr in &shd.faces {
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else {
                    fi += 1;
                    continue;
                };

                // C2: zero normal — approximate from surface
                if fd
                    .surface
                    .as_ref()
                    .map(|s| s.normal_at(0.0, 0.0).length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE)
                    .unwrap_or(true)
                {
                    issues.push(CheckIssue::ZeroNormal {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                }

                // C3: degenerate face
                if wd.edges.len() < 3 {
                    issues.push(CheckIssue::DegenerateFace {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                    fi += 1;
                    continue;
                }

                // C4: edge index bounds + vertex pairs
                let wire_verts = wire_vertex_pairs(brep, wd);
                let mut has_invalid_edge = false;
                for wesr in &wd.edges {
                    if wesr.index >= n_tshapes {
                        issues.push(CheckIssue::InvalidEdgeIndex {
                            solid: si,
                            shell: shi,
                            face: fi,
                            edge_idx: wesr.index,
                        });
                        has_invalid_edge = true;
                    }
                }
                if has_invalid_edge || wire_verts.len() != wd.edges.len() {
                    fi += 1;
                    continue;
                }

                // C1: wire closure
                let n = wire_verts.len();
                for i in 0..n {
                    let next = (i + 1) % n;
                    let end_v = wire_verts[i].1;
                    let start_v = wire_verts[next].0;
                    if end_v != start_v {
                        let end_pt = brep.vertex_point(end_v).unwrap_or_default();
                        let start_pt = brep.vertex_point(start_v).unwrap_or_default();
                        if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
                            issues.push(CheckIssue::OpenWire {
                                solid: si,
                                shell: shi,
                                face: fi,
                                wire_pos: i,
                            });
                        }
                    }
                }

                // C7: wire self-intersection
                check_wire_self_intersection(&wire_verts, brep, si, shi, fi, 0, &mut issues);

                // C8: geometric self-intersection
                check_geometric_self_intersection(&wire_verts, brep, si, shi, fi, &mut issues);

                // Check inner wires
                for (wi, iw_sr) in fd.inner_wires.iter().enumerate() {
                    let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
                        continue;
                    };
                    let inner_verts = wire_vertex_pairs(brep, iwd);
                    if inner_verts.len() < 2 {
                        continue;
                    }

                    let n_inner = inner_verts.len();
                    for i in 0..n_inner {
                        let next = (i + 1) % n_inner;
                        let end_v = inner_verts[i].1;
                        let start_v = inner_verts[next].0;
                        if end_v != start_v {
                            let end_pt = brep.vertex_point(end_v).unwrap_or_default();
                            let start_pt = brep.vertex_point(start_v).unwrap_or_default();
                            if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
                                issues.push(CheckIssue::OpenWire {
                                    solid: si,
                                    shell: shi,
                                    face: fi,
                                    wire_pos: i,
                                });
                            }
                        }
                    }
                    check_wire_self_intersection(
                        &inner_verts,
                        brep,
                        si,
                        shi,
                        fi,
                        wi + 1,
                        &mut issues,
                    );
                }
                fi += 1;
            }
            shi += 1;
        }
        si += 1;
    }

    CheckResult { issues }
}

/// Convenience short alias for [`brep_check_analyze`].
pub fn check_brep(brep: &rcad_kernel::BRep) -> CheckResult {
    brep_check_analyze(brep)
}

/// Check a single wire for self-intersecting topology.
///
/// A valid wire wire should have each vertex appear at most twice across
/// all edge endpoints: once as the start of some edge and once as the end
/// of another edge. If a vertex appears 3+ times, the wire self-intersects.
///
/// Aligned with OCCT BRepCheck_Wire::SelfIntersection concept
/// (BRepCheck_Wire.cxx lines ~100-145).
fn check_wire_self_intersection(
    wire_verts: &[(usize, usize)],
    brep: &rcad_kernel::BRep,
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
    for (&vidx, &count) in &vertex_count {
        if count > 2 {
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

/// Check whether non-adjacent edges of a wire intersect geometrically.
///
/// Projects the wire edge endpoints onto the face's 2D plane (using any two
/// non-collinear edges to form a local basis) and runs 2D segment intersection
/// tests on all non-adjacent edge pairs.
fn check_geometric_self_intersection(
    wire_verts: &[(usize, usize)],
    brep: &rcad_kernel::BRep,
    solid: usize,
    shell: usize,
    face: usize,
    issues: &mut Vec<CheckIssue>,
) {
    let n = wire_verts.len();
    if n < 4 {
        return;
    }

    let segs: Vec<(DVec3, DVec3)> = wire_verts
        .iter()
        .map(|&(sv, ev)| {
            let p0 = brep.vertex_point(sv).unwrap_or(DVec3::ZERO);
            let p1 = brep.vertex_point(ev).unwrap_or(DVec3::ZERO);
            (p0, p1)
        })
        .collect();

    let (origin, axis_u, axis_v) = {
        let mut found = None;
        for i in 0..n {
            let d = segs[i].1 - segs[i].0;
            if d.length() > TOLERANCE_LEN_MIN {
                let u = d.normalize();
                let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                let v = u.cross(tmp).normalize();
                found = Some((segs[i].0, u, v));
                break;
            }
        }
        match found {
            Some(b) => b,
            None => return,
        }
    };

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(axis_u), d.dot(axis_v)]
    };

    let seg2d: Vec<([f64; 2], [f64; 2])> = segs
        .iter()
        .map(|&(p0, p1)| (project(p0), project(p1)))
        .collect();

    for i in 0..n {
        for j in (i + 2)..n {
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
                return;
            }
        }
    }
}

/// Returns `true` if the open segment p1-p2 properly intersects segment p3-p4.
fn segments_2d_properly_intersect(p1: [f64; 2], p2: [f64; 2], p3: [f64; 2], p4: [f64; 2]) -> bool {
    let d1 = [p2[0] - p1[0], p2[1] - p1[1]];
    let d2 = [p4[0] - p3[0], p4[1] - p3[1]];
    let cross = d1[0] * d2[1] - d1[1] * d2[0];
    if cross.abs() < TOLERANCE_LEN_MIN {
        return false;
    }
    let dx = p3[0] - p1[0];
    let dy = p3[1] - p1[1];
    let t = (dx * d2[1] - dy * d2[0]) / cross;
    let s = (dx * d1[1] - dy * d1[0]) / cross;
    let eps = TOLERANCE_COORD_SUB;
    t > eps && t < 1.0 - eps && s > eps && s < 1.0 - eps
}

fn count_geometric_self_intersections(wire_verts: &[(usize, usize)], v_points: &[DVec3]) -> usize {
    let n = wire_verts.len();
    if n < 4 {
        return 0;
    }
    let segs: Vec<(DVec3, DVec3)> = wire_verts
        .iter()
        .map(|&(sv, ev)| {
            let p0 = v_points.get(sv).copied().unwrap_or(DVec3::ZERO);
            let p1 = v_points.get(ev).copied().unwrap_or(DVec3::ZERO);
            (p0, p1)
        })
        .collect();
    let (origin, axis_u, axis_v) = {
        let mut found = None;
        for i in 0..n {
            let d = segs[i].1 - segs[i].0;
            if d.length() > TOLERANCE_LEN_MIN {
                let u = d.normalize();
                let tmp = if u.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                let v = u.cross(tmp).normalize();
                found = Some((segs[i].0, u, v));
                break;
            }
        }
        match found {
            Some(b) => b,
            None => return 0,
        }
    };
    let project = |p: DVec3| -> [f64; 2] {
        let d = p - origin;
        [d.dot(axis_u), d.dot(axis_v)]
    };
    let seg2d: Vec<([f64; 2], [f64; 2])> = segs
        .iter()
        .map(|&(p0, p1)| (project(p0), project(p1)))
        .collect();
    let mut count = 0usize;
    for i in 0..n {
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            if segments_2d_properly_intersect(seg2d[i].0, seg2d[i].1, seg2d[j].0, seg2d[j].1) {
                count += 1;
            }
        }
    }
    count
}

// ----- OCCT BRepCheck alignment: Shell/Wire/Face validation -----
//
// Functions in this section are aligned with OCCT's BRepCheck classes:
// BRepCheck_Shell.cxx  (Shell closure - each edge shared by exactly 2 faces)
// BRepCheck_Wire.cxx (Wire closure + self-intersection)
// BRepCheck_Face.cxx (Wire-on-surface check)
//
// OCCT source: $OCCT_SRC/src/BRepCheck/
// -----------------------------------------------------------------

/// Find a face by its flat (global) index across all solids/shells.
///
/// Returns `(solid_idx, shell_idx, face_idx, &TFaceData)` or `None` if the index
/// is out of range.
fn find_face_by_flat_idx<'a>(
    brep: &'a rcad_kernel::BRep,
    flat_idx: usize,
) -> Option<(usize, usize, usize, &'a topods::TFaceData)> {
    let mut idx = 0usize;
    let mut si = 0usize;
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        let mut shi = 0usize;
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            let mut fi = 0usize;
            for face_sr in &shd.faces {
                if idx == flat_idx {
                    let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                        return None;
                    };
                    return Some((si, shi, fi, fd));
                }
                idx += 1;
                fi += 1;
            }
            shi += 1;
        }
        si += 1;
    }
    None
}

/// BRepCheck_Wire::Closed equivalent.
///
/// Checks that every wire belonging to the face at `face_idx` (flat index)
/// forms a closed loop.
///
/// Aligned with OCCT BRepCheck_Wire::Closed (BRepCheck_Wire.cxx lines ~60-95)
pub fn check_wire_closed(brep: &rcad_kernel::BRep, face_idx: usize) -> bool {
    let (_, _, _, fd) = match find_face_by_flat_idx(brep, face_idx) {
        Some(f) => f,
        None => return false,
    };
    if !check_single_wire_closed(brep, fd.outer_wire) {
        return false;
    }
    for iw_sr in &fd.inner_wires {
        if !check_single_wire_closed(brep, *iw_sr) {
            return false;
        }
    }
    true
}

/// Internal helper: checks closure of a single wire (given by ShapeRef).
fn check_single_wire_closed(brep: &rcad_kernel::BRep, wire_sr: ShapeRef) -> bool {
    let TShape::Wire(wd) = &*brep.tshapes[wire_sr.index] else {
        return true;
    };
    let verts = wire_vertex_pairs(brep, wd);
    let n = verts.len();
    if n == 0 {
        return true;
    }
    if n == 1 {
        let (sv, ev) = verts[0];
        if sv == ev {
            return true;
        }
        let s_pt = brep.vertex_point(sv).unwrap_or_default();
        let e_pt = brep.vertex_point(ev).unwrap_or_default();
        return (s_pt - e_pt).length() <= TOLERANCE_MESH_LEGACY;
    }
    for i in 0..n {
        let next = (i + 1) % n;
        let (_, ev) = verts[i];
        let (sv, _) = verts[next];
        if ev != sv {
            let end_pt = brep.vertex_point(ev).unwrap_or_default();
            let start_pt = brep.vertex_point(sv).unwrap_or_default();
            if (end_pt - start_pt).length() > TOLERANCE_MESH_LEGACY {
                return false;
            }
        }
    }
    true
}

/// BRepCheck_Wire::SelfIntersection equivalent.
///
/// Returns a list of `(edge_idx_in_wire_a, edge_idx_in_wire_b)` pairs for
/// each topological self-intersection found in any wire of the face.
///
/// Aligned with OCCT BRepCheck_Wire::SelfIntersection (BRepCheck_Wire.cxx lines ~100-145)
pub fn check_wire_self_intersection_pairs(
    brep: &rcad_kernel::BRep,
    face_idx: usize,
) -> Vec<(usize, usize)> {
    let (_, _, _, fd) = match find_face_by_flat_idx(brep, face_idx) {
        Some(f) => f,
        None => return Vec::new(),
    };
    let mut result = Vec::new();
    let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
        return result;
    };
    result.extend(check_single_wire_self_intersection_pairs(brep, owd));
    for iw_sr in &fd.inner_wires {
        let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
            continue;
        };
        result.extend(check_single_wire_self_intersection_pairs(brep, iwd));
    }
    result
}

/// Internal helper: finds self-intersecting edge pairs in a single wire.
fn check_single_wire_self_intersection_pairs(
    brep: &rcad_kernel::BRep,
    wd: &topods::TWireData,
) -> Vec<(usize, usize)> {
    use std::collections::HashMap;
    let n = wd.edges.len();
    if n < 4 {
        return Vec::new();
    }
    let mut vertex_occurrences: HashMap<usize, Vec<(usize, bool)>> = HashMap::new();
    for (i, wesr) in wd.edges.iter().enumerate() {
        if let Some(ed) = e_edge_data(brep, wesr.index) {
            let (sv, ev) = if wesr.orientation.is_forward() {
                (ed.first.index, ed.last.index)
            } else {
                (ed.last.index, ed.first.index)
            };
            vertex_occurrences.entry(sv).or_default().push((i, true));
            vertex_occurrences.entry(ev).or_default().push((i, false));
        }
    }
    let mut pairs = Vec::new();
    for (&_vidx, occurrences) in &vertex_occurrences {
        if occurrences.len() <= 2 {
            continue;
        }
        let edge_positions: Vec<usize> = occurrences.iter().map(|(pos, _)| *pos).collect();
        for a in 0..edge_positions.len() {
            for b in (a + 1)..edge_positions.len() {
                let ea = edge_positions[a];
                let eb = edge_positions[b];
                let diff = if ea > eb { ea - eb } else { eb - ea };
                let is_adjacent = diff == 1 || (ea == 0 && eb == n - 1) || (eb == 0 && ea == n - 1);
                if is_adjacent {
                    continue;
                }
                let pair = if ea < eb { (ea, eb) } else { (eb, ea) };
                if !pairs.contains(&pair) {
                    pairs.push(pair);
                }
            }
        }
    }
    pairs
}

/// BRepCheck_Face::Intersection equivalent (wire-on-surface check).
///
/// Checks that every edge in the face's wires lies on the face surface
/// within the given tolerance.
///
/// Aligned with OCCT BRepCheck_Face::Intersection (BRepCheck_Face.cxx lines ~70-130)
pub fn check_face_wire_on_surface(
    brep: &rcad_kernel::BRep,
    face_idx: usize,
    tolerance: f64,
) -> bool {
    let (_, _, _, fd) = match find_face_by_flat_idx(brep, face_idx) {
        Some(f) => f,
        None => return false,
    };
    let Some(surface) = fd.surface.as_ref() else {
        return true;
    };
    for wesr in {
        let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
            return true;
        };
        &owd.edges
    } {
        if !check_edge_on_surface(brep, wesr.index, surface, tolerance) {
            return false;
        }
    }
    for iw_sr in &fd.inner_wires {
        let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
            continue;
        };
        for wesr in &iwd.edges {
            if !check_edge_on_surface(brep, wesr.index, surface, tolerance) {
                return false;
            }
        }
    }
    true
}

/// Check that a single edge's 3D curve lies on the given surface within tolerance.
fn check_edge_on_surface(
    brep: &rcad_kernel::BRep,
    edge_idx: usize,
    surface: &rcad_kernel::geom::Surface3,
    tolerance: f64,
) -> bool {
    let Some(ed) = e_edge_data(brep, edge_idx) else {
        return true;
    };
    let Some(curve) = ed.curve.as_ref() else {
        return true;
    };
    let range = ed.range;
    // Find pcurve for this surface by matching pcurves — we iterate and
    // check each pcurve's surface (the face that owns it).
    for (_pc_face_idx, (curve2d, _t1, _t2)) in &ed.pcurves {
        let sample_ts = [0.0, 0.5, 1.0];
        for &t_frac in &sample_ts {
            let t3 = range[0] + t_frac * (range[1] - range[0]);
            let p3d = curve.point_at(t3);
            let uv_range = ed.range;
            let t2 = uv_range[0] + t_frac * (uv_range[1] - uv_range[0]);
            let uv = curve2d.point_at(t2);
            let p_surf = surface.point_at(uv.x, uv.y);
            if (p3d - p_surf).length() > tolerance {
                return false;
            }
        }
    }
    true
}

//    SameParameter diagnosis

/// A single edge whose 3D curve endpoints deviate from the vertex positions.
#[derive(Debug, Clone)]
pub struct SuspectEdge {
    pub edge_idx: usize,
    pub start_gap: f64,
    pub end_gap: f64,
}

/// Report from [`diagnose_same_parameter`].
#[derive(Debug, Clone, Default)]
pub struct SameParameterDiagnosis {
    pub suspect_edges: Vec<SuspectEdge>,
}

impl SameParameterDiagnosis {
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// A single edge whose PCurve ranges deviate from the 3D edge range.
#[derive(Debug, Clone)]
pub struct SuspectSameRangeEdge {
    pub edge_idx: usize,
    pub mismatched_pcurves: usize,
    pub max_delta: f64,
}

/// Report from [`diagnose_same_range`].
#[derive(Debug, Clone, Default)]
pub struct SameRangeDiagnosis {
    pub suspect_edges: Vec<SuspectSameRangeEdge>,
}

impl SameRangeDiagnosis {
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// A single edge-PCurve pair whose UV-evaluated surface endpoints do not match
/// the edge's 3D curve endpoints.
#[derive(Debug, Clone)]
pub struct SuspectFaceSurfaceEdge {
    pub edge_idx: usize,
    pub pcurve_pos: usize,
    pub surface_idx: usize,
    pub start_gap: f64,
    pub end_gap: f64,
    pub max_gap: f64,
}

/// Report from [`diagnose_face_surface_consistency`].
#[derive(Debug, Clone, Default)]
pub struct FaceSurfaceConsistencyDiagnosis {
    pub suspect_edges: Vec<SuspectFaceSurfaceEdge>,
}

impl FaceSurfaceConsistencyDiagnosis {
    pub fn is_clean(&self) -> bool {
        self.suspect_edges.is_empty()
    }
}

/// Diagnose face-on-surface consistency via PCurves.
pub fn diagnose_face_surface_consistency(
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> FaceSurfaceConsistencyDiagnosis {
    use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
    let mut suspect_edges = Vec::new();

    for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = &**ts else { continue };
        let Some(curve3) = ed.curve.as_ref() else {
            continue;
        };
        let range3 = ed.range;

        for (pcurve_pos, (&pc_face_idx, (curve2d, t1, t2))) in ed.pcurves.iter().enumerate() {
            let TShape::Face(pc_fd) = &*brep.tshapes[pc_face_idx] else {
                continue;
            };
            let Some(surface) = pc_fd.surface.as_ref() else {
                continue;
            };

            let range2 = [*t1, *t2];

            let p3_start = curve3.point_at(range3[0]);
            let p3_end = curve3.point_at(range3[1]);
            let uv_start = curve2d.point_at(range2[0]);
            let uv_end = curve2d.point_at(range2[1]);
            let ps_start = surface.point_at(uv_start.x, uv_start.y);
            let ps_end = surface.point_at(uv_end.x, uv_end.y);

            let start_gap = (ps_start - p3_start).length();
            let end_gap = (ps_end - p3_end).length();
            let max_gap = start_gap.max(end_gap);

            if max_gap > tolerance {
                suspect_edges.push(SuspectFaceSurfaceEdge {
                    edge_idx,
                    pcurve_pos,
                    surface_idx: pc_face_idx,
                    start_gap,
                    end_gap,
                    max_gap,
                });
            }
        }
    }

    FaceSurfaceConsistencyDiagnosis { suspect_edges }
}

/// Per-wire diagnostics for gap and self-intersection analysis.
#[derive(Debug, Clone, Default)]
pub struct WireIssueReport {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    pub wire_idx: usize,
    pub edge_count: usize,
    pub open_gaps: usize,
    pub topological_self_intersections: usize,
    pub geometric_self_intersections: usize,
}

/// Aggregated report from [`analyze_wire_issues`].
#[derive(Debug, Clone, Default)]
pub struct WireAnalysisReport {
    pub wires: Vec<WireIssueReport>,
    pub total_open_gaps: usize,
    pub total_topological_self_intersections: usize,
    pub total_geometric_self_intersections: usize,
}

impl WireAnalysisReport {
    pub fn is_clean(&self) -> bool {
        self.total_open_gaps == 0
            && self.total_topological_self_intersections == 0
            && self.total_geometric_self_intersections == 0
    }
}

/// Analyze all face wires for gap and self-intersection issues.
///
/// This is a structured counterpart to checker issues C1/C7/C8 and is useful
/// for import diagnostics and healing reports.
/// Analyze wire issues in a BRep by examining topology directly.
/// Works on topods::BRep — no old BRep bridge needed.
pub fn analyze_wire_issues(brep: &topods::BRep, tolerance: f64) -> WireAnalysisReport {
    let tshapes = &brep.tshapes;
    let mut report = WireAnalysisReport::default();

    let v_points: Vec<DVec3> = tshapes
        .iter()
        .filter_map(|ts| {
            if let topods::TShape::Vertex(vd) = &**ts {
                Some(vd.point)
            } else {
                None
            }
        })
        .collect();

    let mut si = 0usize;
    for ts in tshapes {
        let topods::TShape::Solid(sd) = &**ts else {
            continue;
        };
        let mut shi = 0usize;
        for shell_sr in &sd.shells {
            let topods::TShape::Shell(shd) = &*tshapes[shell_sr.index] else {
                continue;
            };
            let mut fi = 0usize;
            for face_sr in &shd.faces {
                let topods::TShape::Face(fd) = &*tshapes[face_sr.index] else {
                    continue;
                };

                let mut all_wires: Vec<(usize, &topods::ShapeRef)> = Vec::new();
                all_wires.push((0, &fd.outer_wire));
                for (wi, inner) in fd.inner_wires.iter().enumerate() {
                    all_wires.push((wi + 1, inner));
                }

                for (wire_idx, wire_sr) in all_wires {
                    let topods::TShape::Wire(wd) = &*tshapes[wire_sr.index] else {
                        continue;
                    };
                    let mut wire_verts = Vec::with_capacity(wd.edges.len());
                    let mut valid = true;

                    for wesr in &wd.edges {
                        if wesr.index >= tshapes.len() {
                            valid = false;
                            break;
                        }
                        let topods::TShape::Edge(ed) = &*tshapes[wesr.index] else {
                            valid = false;
                            break;
                        };
                        let (sv, ev) = if wesr.orientation.is_forward() {
                            (ed.first.index, ed.last.index)
                        } else {
                            (ed.last.index, ed.first.index)
                        };
                        if sv >= v_points.len() || ev >= v_points.len() {
                            valid = false;
                            break;
                        }
                        wire_verts.push((sv, ev));
                    }
                    if !valid {
                        continue;
                    }

                    let mut open_gaps = 0usize;
                    let n = wire_verts.len();
                    if n > 1 {
                        for i in 0..n {
                            let next = (i + 1) % n;
                            let end_v = wire_verts[i].1;
                            let start_v = wire_verts[next].0;
                            if end_v != start_v {
                                let end_pt = v_points[end_v];
                                let start_pt = v_points[start_v];
                                if (end_pt - start_pt).length() > tolerance {
                                    open_gaps += 1;
                                }
                            }
                        }
                    }

                    let mut vertex_count: std::collections::HashMap<usize, usize> =
                        std::collections::HashMap::new();
                    for &(sv, ev) in &wire_verts {
                        *vertex_count.entry(sv).or_insert(0) += 1;
                        *vertex_count.entry(ev).or_insert(0) += 1;
                    }
                    let topological_self_intersections =
                        vertex_count.values().filter(|&&c| c > 2).count();
                    let geometric_self_intersections =
                        count_geometric_self_intersections(&wire_verts, &v_points);

                    if open_gaps > 0
                        || topological_self_intersections > 0
                        || geometric_self_intersections > 0
                    {
                        report.wires.push(WireIssueReport {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx,
                            edge_count: wire_verts.len(),
                            open_gaps,
                            topological_self_intersections,
                            geometric_self_intersections,
                        });
                    }
                    report.total_open_gaps += open_gaps;
                    report.total_topological_self_intersections += topological_self_intersections;
                    report.total_geometric_self_intersections += geometric_self_intersections;
                }
                fi += 1;
            }
            shi += 1;
        }
        si += 1;
    }
    report
}

/// Legacy internal helper: old flat BRep for tests that build flat structures.
#[cfg(test)]
pub(crate) fn analyze_wire_issues_flat(
    brep: &rcad_kernel::BRep,
    tolerance: f64,
) -> WireAnalysisReport {
    analyze_wire_issues(brep, tolerance)
}

/// Scan all edges in `brep` for SameParameter violations.
pub fn diagnose_same_parameter(brep: &rcad_kernel::BRep, tolerance: f64) -> SameParameterDiagnosis {
    use rcad_kernel::geom::CurveEval;
    let mut suspects = Vec::new();

    for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = &**ts else { continue };
        let Some(curve) = ed.curve.as_ref() else {
            continue;
        };
        let range = ed.range;

        let p_start = match brep.vertex_point(ed.first.index) {
            Some(p) => p,
            None => continue,
        };
        let p_end = match brep.vertex_point(ed.last.index) {
            Some(p) => p,
            None => continue,
        };

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

    SameParameterDiagnosis {
        suspect_edges: suspects,
    }
}

/// Scan all edges in `brep` for SameRange violations.
pub fn diagnose_same_range(brep: &rcad_kernel::BRep, tolerance: f64) -> SameRangeDiagnosis {
    let mut suspects = Vec::new();

    for (edge_idx, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = &**ts else { continue };
        let range3d = ed.range;
        if ed.pcurves.is_empty() {
            continue;
        }

        let mut mismatched_pcurves = 0usize;
        let mut max_delta = 0.0f64;

        for (_pc_face_idx, (_curve2d, t1, t2)) in &ed.pcurves {
            let d0 = (*t1 - range3d[0]).abs();
            let d1 = (*t2 - range3d[1]).abs();
            let d = d0.max(d1);
            if d > tolerance {
                mismatched_pcurves += 1;
                max_delta = max_delta.max(d);
            }
        }

        if mismatched_pcurves > 0 {
            suspects.push(SuspectSameRangeEdge {
                edge_idx,
                mismatched_pcurves,
                max_delta,
            });
        }
    }

    SameRangeDiagnosis {
        suspect_edges: suspects,
    }
}

//    Shell topology analysis

/// Topology analysis report for a BRep's shell structure.
#[derive(Debug, Clone, Default)]
pub struct ShellTopologyReport {
    pub is_closed: bool,
    pub is_manifold: bool,
    pub open_edge_count: usize,
    pub non_manifold_edge_count: usize,
    pub isolated_vertex_count: usize,
    pub total_edges: usize,
    pub total_faces: usize,
}

/// Analyze the shell topology of `brep`.
pub fn analyze_shell_topology(brep: &rcad_kernel::BRep) -> ShellTopologyReport {
    let total_edges = brep.edge_count();
    let mut edge_face_count: Vec<usize> = vec![0; brep.tshapes.len()];
    let mut total_faces = 0usize;

    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                total_faces += 1;
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    continue;
                };
                for wesr in &owd.edges {
                    if wesr.index < edge_face_count.len() {
                        edge_face_count[wesr.index] += 1;
                    }
                }
                for iw_sr in &fd.inner_wires {
                    let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
                        continue;
                    };
                    for wesr in &iwd.edges {
                        if wesr.index < edge_face_count.len() {
                            edge_face_count[wesr.index] += 1;
                        }
                    }
                }
            }
        }
    }

    let open_edge_count = edge_face_count.iter().filter(|&&c| c == 1).count();
    let non_manifold_edge_count = edge_face_count.iter().filter(|&&c| c > 2).count();

    let mut vertex_used = vec![false; brep.vertex_count()];
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Edge(ed) = &**ts else { continue };
        if ed.first.index < vertex_used.len() {
            vertex_used[ed.first.index] = true;
        }
        if ed.last.index < vertex_used.len() {
            vertex_used[ed.last.index] = true;
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

//    Euler characteristic analysis

/// Euler characteristic and topological genus for a single solid.
#[derive(Debug, Clone)]
pub struct EulerAnalysis {
    pub solid_idx: usize,
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub euler_number: i64,
    pub is_closed: bool,
    pub genus: Option<i64>,
}

/// Compute per-solid Euler analysis for every solid in `brep`.
pub fn euler_analysis(brep: &rcad_kernel::BRep) -> Vec<EulerAnalysis> {
    let mut results = Vec::new();
    let mut si = 0usize;

    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };

        let mut unique_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut face_count = 0usize;

        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                face_count += 1;
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    continue;
                };
                for wesr in &owd.edges {
                    unique_edges.insert(wesr.index);
                }
                for iw_sr in &fd.inner_wires {
                    let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
                        continue;
                    };
                    for wesr in &iwd.edges {
                        unique_edges.insert(wesr.index);
                    }
                }
            }
        }

        let mut unique_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &ei in &unique_edges {
            if let Some(ed) = e_edge_data(brep, ei) {
                unique_verts.insert(ed.first.index);
                unique_verts.insert(ed.last.index);
            }
        }

        let v = unique_verts.len();
        let e = unique_edges.len();
        let f = face_count;
        let euler_number = v as i64 - e as i64 + f as i64;

        let mut edge_face_count: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    continue;
                };
                for wesr in &owd.edges {
                    *edge_face_count.entry(wesr.index).or_insert(0) += 1;
                }
                for iw_sr in &fd.inner_wires {
                    let TShape::Wire(iwd) = &*brep.tshapes[iw_sr.index] else {
                        continue;
                    };
                    for wesr in &iwd.edges {
                        *edge_face_count.entry(wesr.index).or_insert(0) += 1;
                    }
                }
            }
        }
        let is_closed = edge_face_count.values().all(|&c| c == 2);

        let genus = if is_closed {
            let g = (2 - euler_number) / 2;
            if (2 - euler_number) % 2 == 0 && g >= 0 {
                Some(g)
            } else {
                None
            }
        } else {
            None
        };

        results.push(EulerAnalysis {
            solid_idx: si,
            vertices: v,
            edges: e,
            faces: f,
            euler_number,
            is_closed,
            genus,
        });
        si += 1;
    }

    results
}

//    Orientation consistency analysis

/// A face whose stored normal appears to point inward rather than outward.
#[derive(Debug, Clone)]
pub struct OrientationIssue {
    pub solid_idx: usize,
    pub face_idx: usize,
    pub dot_product: f64,
}

/// Report from [`check_orientation_consistency`].
#[derive(Debug, Clone, Default)]
pub struct OrientationReport {
    pub is_consistent: bool,
    pub issues: Vec<OrientationIssue>,
    pub consistent_face_count: usize,
    pub inconsistent_face_count: usize,
}

/// Check that every face's stored normal points outward from the solid interior.
pub fn check_orientation_consistency(brep: &rcad_kernel::BRep) -> OrientationReport {
    let mut issues = Vec::new();
    let mut consistent_face_count = 0usize;
    let mut inconsistent_face_count = 0usize;
    let mut flat_face_idx = 0usize;

    let mut si = 0usize;
    for ts in &brep.tshapes {
        let TShape::Solid(sd) = &**ts else { continue };

        // Compute solid centroid from all vertex points in solid's faces
        let mut solid_verts = std::collections::HashSet::new();
        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    continue;
                };
                for wesr in &owd.edges {
                    if let Some(ed) = e_edge_data(brep, wesr.index) {
                        solid_verts.insert(ed.first.index);
                        solid_verts.insert(ed.last.index);
                    }
                }
            }
        }

        if solid_verts.is_empty() {
            for shell_sr in &sd.shells {
                let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                    continue;
                };
                flat_face_idx += shd.faces.len();
            }
            si += 1;
            continue;
        }

        let solid_centroid: DVec3 = {
            let sum: DVec3 = solid_verts
                .iter()
                .filter_map(|&vi| brep.vertex_point(vi))
                .sum();
            sum / solid_verts.len() as f64
        };

        for shell_sr in &sd.shells {
            let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else {
                continue;
            };
            for face_sr in &shd.faces {
                let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else {
                    continue;
                };

                // Face centroid from first vertex of each outer-wire edge
                let TShape::Wire(owd) = &*brep.tshapes[fd.outer_wire.index] else {
                    flat_face_idx += 1;
                    continue;
                };
                let mut face_centroid = DVec3::ZERO;
                let mut n = 0usize;
                for wesr in &owd.edges {
                    if let Some(ed) = e_edge_data(brep, wesr.index) {
                        let vi = if wesr.orientation.is_forward() {
                            ed.first.index
                        } else {
                            ed.last.index
                        };
                        if let Some(pt) = brep.vertex_point(vi) {
                            face_centroid += pt;
                            n += 1;
                        }
                    }
                }
                if n == 0 {
                    flat_face_idx += 1;
                    continue;
                }
                face_centroid /= n as f64;

                let normal = fd
                    .surface
                    .as_ref()
                    .map(|s| s.normal_at(0.0, 0.0))
                    .unwrap_or(DVec3::Z);
                let outward = face_centroid - solid_centroid;
                let dot = normal.dot(outward);
                if dot >= 0.0 {
                    consistent_face_count += 1;
                } else {
                    inconsistent_face_count += 1;
                    issues.push(OrientationIssue {
                        solid_idx: si,
                        face_idx: flat_face_idx,
                        dot_product: dot,
                    });
                }
                flat_face_idx += 1;
            }
        }
        si += 1;
    }

    OrientationReport {
        is_consistent: issues.is_empty(),
        issues,
        consistent_face_count,
        inconsistent_face_count,
    }
}

//    Comprehensive richer validity analysis

/// Aggregated validity report combining all available checks.
///
/// This is the RCAD equivalent of OCCT's `BRepCheck_Analyzer` + `ShapeAnalysis`
/// combined output, giving a single entry-point for full BRep validation.
#[derive(Debug, Clone)]
pub struct RicherValidityReport {
    pub check_result: CheckResult,
    pub shell_topology: ShellTopologyReport,
    pub euler: Vec<EulerAnalysis>,
    pub orientation: OrientationReport,
    pub is_fully_valid: bool,
}

impl RicherValidityReport {
    pub fn summary(&self) -> String {
        let issues = self.check_result.issues.len();
        let euler_issues = self
            .euler
            .iter()
            .filter(|e| e.genus.is_none_or(|g| g < 0))
            .count();
        let orient_issues = self.orientation.inconsistent_face_count;
        if self.is_fully_valid {
            format!(
                "valid: {} solids, closed={}, manifold={}",
                self.euler.len(),
                self.shell_topology.is_closed,
                self.shell_topology.is_manifold,
            )
        } else {
            format!(
                "INVALID: {} structural issue(s), {} orientation inconsistency/ies, {} genus anomaly/ies",
                issues, orient_issues, euler_issues,
            )
        }
    }
}

/// Run all available validity checks on `brep` and return a consolidated report.
pub fn richer_validity_analysis(brep: &rcad_kernel::BRep) -> RicherValidityReport {
    let check_result = brep_check_analyze(brep);
    let shell_topology = analyze_shell_topology(brep);
    let euler = euler_analysis(brep);
    let orientation = check_orientation_consistency(brep);
    let is_fully_valid = check_result.is_valid() && orientation.is_consistent;

    RicherValidityReport {
        check_result,
        shell_topology,
        euler,
        orientation,
        is_fully_valid,
    }
}

//    Surface UV Analysis (ShapeAnalysis_Surface equivalent)

/// Report from surface UV domain analysis.
#[derive(Debug, Clone, Default)]
pub struct SurfaceAnalysisReport {
    pub faces_analyzed: usize,
    pub faces_with_uv_bounds_violation: Vec<UvBoundsViolation>,
    pub total_issues: usize,
}

impl SurfaceAnalysisReport {
    pub fn is_clean(&self) -> bool {
        self.total_issues == 0
    }
    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!("{} faces analyzed, no UV issues", self.faces_analyzed)
        } else {
            format!(
                "{} faces analyzed, {} UV bounds violations",
                self.faces_analyzed,
                self.faces_with_uv_bounds_violation.len()
            )
        }
    }
}

/// UV bounds violation for a single face.
#[derive(Debug, Clone)]
pub struct UvBoundsViolation {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    pub surface_type: String,
    pub expected_bounds: [f64; 4],
    pub actual_bounds: [f64; 4],
    pub violation: f64,
}

/// Quick structural validity check: at least one solid or face exists (topods::BRep).
pub fn is_valid_brep(brep: &rcad_kernel::topods::BRep) -> bool {
    brep.tshapes.iter().any(|ts| {
        matches!(
            ts.as_ref(),
            rcad_kernel::topods::TShape::Solid(_) | rcad_kernel::topods::TShape::Face(_)
        )
    })
}

include!("e1.rs");
