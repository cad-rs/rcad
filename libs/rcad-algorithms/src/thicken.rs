//! Shell thickening  ?analogous to OCCT `BRepOffsetAPI_MakeThickSolid`.
//!
//! # Algorithm
//!
//! 1. Identify boundary wires (edges appearing in exactly one face).
//! 2. Offset each face along its normal by the given thickness.
//! 3. For each boundary edge, create a lateral ruled face connecting
//! the original edge to the corresponding offset edge.
//! 4. Assemble offset faces + lateral faces into a closed solid.
//!
//! # Supported surfaces
//!
//! Plane, Sphere, Cylinder, Cone, Torus  ?each has a known parallel-surface
//! construction. B-spline and trimmed surfaces are skipped.
//!
//! # Features
//!
//! - **Face selection strategies**: Automatic selection for thin-wall features,
//! connectivity-based selection, area-based selection
//! - **Lateral face handling**: Configurable creation, tangency, and splitting
//! - **Thickness variation**: Variable thickness by face region with smooth transitions
//! - **Self-intersection handling**: Detection, automatic thickness reduction, warnings

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::topods::{Orientation, ShapeRef, TShape};
use rcad_kernel::topology::{Face, Shell, Wire, WireEdge};
use std::collections::{HashMap, HashSet};

use crate::tolerance::*;
use crate::triangulate::{TessellationParams, mesh_brep};

// ---------------------------------------------------------------------------
// Backward-compat helpers for accessing new topods::BRep data
// ---------------------------------------------------------------------------

/// Get vertex point by flat tshape index (replaces brep.vertices[vi].point).
fn brep_vertex_point(brep: &rcad_kernel::BRep, idx: usize) -> Option<DVec3> {
    brep.vertex_point(idx)
}

/// Get edge endpoint indices by flat tshape index (replaces brep.edges[ei]).
fn brep_edge_endpoints(brep: &rcad_kernel::BRep, idx: usize) -> Option<(usize, usize)> {
    brep.tshapes.get(idx).and_then(|ts| match ts.as_ref() {
        TShape::Edge(ed) => Some((ed.first.index, ed.last.index)),
        _ => None,
    })
}

/// Get face surface by flat face index in the first shell (replaces brep.geom.face_surface + brep.geom.surfaces).
fn brep_face_surface(brep: &rcad_kernel::BRep, flat_idx: usize) -> Option<Surface3> {
    let solid_idx = brep
        .tshapes
        .iter()
        .position(|ts| matches!(ts.as_ref(), TShape::Solid(_)))?;
    let sd = match &*brep.tshapes[solid_idx] {
        TShape::Solid(sd) => sd,
        _ => return None,
    };
    let shell_sr = sd.shells.first()?;
    let shelld = match &*brep.tshapes[shell_sr.index] {
        TShape::Shell(s) => s,
        _ => return None,
    };
    shelld
        .faces
        .get(flat_idx)
        .and_then(|face_sr| match &*brep.tshapes[face_sr.index] {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        })
}

/// Get first shell faces vec from a topods BRep (replaces brep.solids[0].shells[0].faces).
fn first_shell_face_count(brep: &rcad_kernel::BRep) -> Option<usize> {
    let solid_idx = brep
        .tshapes
        .iter()
        .position(|ts| matches!(ts.as_ref(), TShape::Solid(_)))?;
    let sd = match &*brep.tshapes[solid_idx] {
        TShape::Solid(sd) => sd,
        _ => return None,
    };
    let shell_sr = sd.shells.first()?;
    let shelld = match &*brep.tshapes[shell_sr.index] {
        TShape::Shell(s) => s,
        _ => return None,
    };
    Some(shelld.faces.len())
}

/// Build an old-style Shell from topods::BRep, extracting face/wire/edge data
/// with normals computed from face surfaces.
fn shell_from_brep(brep: &rcad_kernel::BRep) -> Option<Shell> {
    let solid_idx = brep
        .tshapes
        .iter()
        .position(|ts| matches!(ts.as_ref(), TShape::Solid(_)))?;
    let sd = match &*brep.tshapes[solid_idx] {
        TShape::Solid(sd) => sd,
        _ => return None,
    };
    let shell_sr = sd.shells.first()?;
    let shelld = match &*brep.tshapes[shell_sr.index] {
        TShape::Shell(s) => s,
        _ => return None,
    };
    let mut faces = Vec::new();
    for face_sr in &shelld.faces {
        let fd = match &*brep.tshapes[face_sr.index] {
            TShape::Face(fd) => fd,
            _ => continue,
        };
        let wd = match &*brep.tshapes[fd.outer_wire.index] {
            TShape::Wire(w) => w,
            _ => continue,
        };
        let outer_edges: Vec<WireEdge> = wd
            .edges
            .iter()
            .map(|e_sr| WireEdge::new(e_sr.index, e_sr.orientation == Orientation::Forward))
            .collect();
        let inner_wires: Vec<Wire> = fd
            .inner_wires
            .iter()
            .map(|iw_sr| {
                if iw_sr.index >= brep.tshapes.len() {
                    return Wire { edges: Vec::new() };
                }
                let iwd = match &*brep.tshapes[iw_sr.index] {
                    TShape::Wire(w) => w,
                    _ => return Wire { edges: Vec::new() },
                };
                Wire {
                    edges: iwd
                        .edges
                        .iter()
                        .map(|e_sr| {
                            WireEdge::new(e_sr.index, e_sr.orientation == Orientation::Forward)
                        })
                        .collect(),
                }
            })
            .collect();
        let normal = fd
            .surface
            .as_ref()
            .map(|s| s.normal_at(0.0, 0.0))
            .unwrap_or(DVec3::Z);
        faces.push(Face {
            outer_wire: Wire { edges: outer_edges },
            inner_wires,
            normal,
            triangles: Vec::new(),
            sample_point: fd.sample_point,
            mesh_dirty: true,
            surface_idx: Some(face_sr.index),
        });
    }
    Some(Shell { faces })
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Result Types
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Result of a thickening operation.
#[derive(Debug, Clone)]
pub struct ThickeningResult {
    /// The thickened solid as a new BRep.
    pub brep: rcad_kernel::BRep,
    /// Number of offset faces (one per input face).
    pub offset_faces: usize,
    /// Number of lateral faces connecting boundaries.
    pub lateral_faces: usize,
    /// Whether self-intersection was detected (thickness > half min face distance).
    pub self_intersection: bool,
    /// Warnings generated during the operation.
    pub warnings: Vec<ThickeningWarning>,
    /// If thickness was auto-reduced, this is the actual thickness used.
    pub actual_thickness: Option<f64>,
}

/// Warnings that can occur during thickening.
#[derive(Debug, Clone)]
pub enum ThickeningWarning {
    /// Self-intersection detected.
    SelfIntersection {
        /// Minimum distance between non-adjacent faces.
        min_distance: f64,
        /// Requested thickness.
        requested_thickness: f64,
    },
    /// Thickness was auto-reduced to avoid self-intersection.
    ThicknessAutoReduced {
        /// Original thickness.
        original: f64,
        /// Reduced thickness.
        reduced: f64,
    },
    /// Surface became degenerate during offset.
    DegenerateSurface {
        /// Face index.
        face_index: usize,
        /// Original surface type.
        surface_type: String,
    },
    /// Thin region detected where thickness may be problematic.
    ThinRegionDetected {
        /// Center of thin region.
        center: DVec3,
        /// Thickness at this region.
        thickness: f64,
    },
    /// Face was skipped due to unsupported geometry.
    SkippedFace {
        /// Face index.
        face_index: usize,
        /// Reason for skipping.
        reason: String,
    },
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Face Selection Strategies
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Strategy for selecting faces to remove during thickening.
#[derive(Debug, Clone)]
pub enum FaceSelectionStrategy {
    /// Use explicitly provided face indices.
    Explicit(Vec<usize>),
    /// Automatically select faces for thin-wall features.
    /// Selects faces that are roughly parallel pairs.
    AutoThinWall {
        /// Minimum area ratio for face pair consideration.
        area_ratio_threshold: f64,
        /// Maximum angle (radians) for parallel face consideration.
        parallel_angle_tolerance: f64,
    },
    /// Select faces by connectivity to a seed face.
    ByConnectivity {
        /// Seed face indices to start from.
        seed_faces: Vec<usize>,
        /// Maximum number of connected faces to select.
        max_faces: Option<usize>,
        /// Whether to stop at sharp edges (angle in radians).
        sharp_edge_angle: Option<f64>,
    },
    /// Select the N largest faces by area.
    ByArea {
        /// Number of faces to select.
        count: usize,
    },
    /// Select faces by normal direction.
    ByNormal {
        /// Target normal direction.
        direction: DVec3,
        /// Angle tolerance (radians).
        angle_tolerance: f64,
    },
    /// Select faces by planar property.
    PlanarOnly {
        /// Whether to include only planar faces.
        include_planar: bool,
    },
}

impl Default for FaceSelectionStrategy {
    fn default() -> Self {
        FaceSelectionStrategy::Explicit(Vec::new())
    }
}

impl FaceSelectionStrategy {
    /// Create an explicit selection strategy.
    pub fn explicit(indices: Vec<usize>) -> Self {
        FaceSelectionStrategy::Explicit(indices)
    }

    /// Create an auto thin-wall selection strategy.
    pub fn auto_thin_wall() -> Self {
        FaceSelectionStrategy::AutoThinWall {
            area_ratio_threshold: 0.5,
            parallel_angle_tolerance: 0.1, // ~5.7 degrees
        }
    }

    /// Create a connectivity-based selection strategy.
    pub fn by_connectivity(seed_faces: Vec<usize>) -> Self {
        FaceSelectionStrategy::ByConnectivity {
            seed_faces,
            max_faces: None,
            sharp_edge_angle: Some(std::f64::consts::PI / 4.0), // 45 degrees
        }
    }

    /// Create an area-based selection strategy.
    pub fn by_area(count: usize) -> Self {
        FaceSelectionStrategy::ByArea { count }
    }

    /// Create a normal-based selection strategy.
    pub fn by_normal(direction: DVec3, angle_tolerance: f64) -> Self {
        FaceSelectionStrategy::ByNormal {
            direction: direction.normalize(),
            angle_tolerance,
        }
    }
}

/// Select faces to remove based on the given strategy.
pub fn select_faces_for_removal(
    brep: &rcad_kernel::BRep,
    strategy: &FaceSelectionStrategy,
) -> Vec<usize> {
    let shell = match shell_from_brep(brep) {
        Some(s) => s,
        None => return Vec::new(),
    };

    match strategy {
        FaceSelectionStrategy::Explicit(indices) => indices
            .iter()
            .filter(|&&i| i < shell.faces.len())
            .copied()
            .collect(),

        FaceSelectionStrategy::AutoThinWall {
            area_ratio_threshold,
            parallel_angle_tolerance,
        } => select_faces_for_thin_wall(
            &shell,
            brep,
            *area_ratio_threshold,
            *parallel_angle_tolerance,
        ),

        FaceSelectionStrategy::ByConnectivity {
            seed_faces,
            max_faces,
            sharp_edge_angle,
        } => select_faces_by_connectivity(&shell, brep, seed_faces, *max_faces, *sharp_edge_angle),

        FaceSelectionStrategy::ByArea { count } => select_faces_by_area(&shell, brep, *count),

        FaceSelectionStrategy::ByNormal {
            direction,
            angle_tolerance,
        } => select_faces_by_normal(&shell, direction, *angle_tolerance),

        FaceSelectionStrategy::PlanarOnly { include_planar } => {
            if *include_planar {
                select_planar_faces(&shell, brep)
            } else {
                select_non_planar_faces(&shell, brep)
            }
        }
    }
}

/// Select faces that form thin-wall features (parallel pairs).
fn select_faces_for_thin_wall(
    shell: &Shell,
    brep: &rcad_kernel::BRep,
    area_ratio_threshold: f64,
    parallel_angle_tolerance: f64,
) -> Vec<usize> {
    let n = shell.faces.len();
    if n < 2 {
        return Vec::new();
    }

    // Compute face areas and normals
    let face_data: Vec<(usize, f64, DVec3)> = shell
        .faces
        .iter()
        .enumerate()
        .map(|(i, face)| {
            let area = compute_face_area(face, brep);
            (i, area, face.normal)
        })
        .collect();

    // Find parallel face pairs
    let mut selected: HashSet<usize> = HashSet::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let (_, area_i, normal_i) = &face_data[i];
            let (_, area_j, normal_j) = &face_data[j];

            // Check if faces are roughly parallel (normals are opposite)
            let dot = normal_i.dot(*normal_j);
            if dot > -1.0 + parallel_angle_tolerance {
                continue;
            }

            // Check area ratio
            let ratio = if *area_i > *area_j {
                area_j / area_i
            } else {
                area_i / area_j
            };

            if ratio > area_ratio_threshold {
                // This is a candidate thin-wall pair
                // Select the smaller face for removal
                if *area_i < *area_j {
                    selected.insert(i);
                } else {
                    selected.insert(j);
                }
            }
        }
    }

    selected.into_iter().collect()
}

/// Select faces by connectivity from seed faces.
fn select_faces_by_connectivity(
    shell: &Shell,
    _brep: &rcad_kernel::BRep,
    seed_faces: &[usize],
    max_faces: Option<usize>,
    sharp_edge_angle: Option<f64>,
) -> Vec<usize> {
    if seed_faces.is_empty() {
        return Vec::new();
    }

    let n = shell.faces.len();
    if n == 0 {
        return Vec::new();
    }

    // Build edge-to-face adjacency
    let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push(fi);
        }
    }

    // Build face adjacency through shared edges
    let mut face_adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for faces in edge_to_faces.values() {
        for &f1 in faces {
            for &f2 in faces {
                if f1 != f2 {
                    face_adjacency[f1].insert(f2);
                }
            }
        }
    }

    // BFS from seed faces
    let mut selected: HashSet<usize> = HashSet::new();
    let mut queue: Vec<usize> = seed_faces.to_vec();

    while let Some(fi) = queue.pop() {
        if fi >= n || selected.contains(&fi) {
            continue;
        }

        if let Some(max) = max_faces
            && selected.len() >= max
        {
            break;
        }

        selected.insert(fi);

        // Add adjacent faces
        if let Some(adjacent) = face_adjacency.get(fi) {
            for &adj in adjacent {
                if !selected.contains(&adj) {
                    // Check sharp edge condition
                    if let Some(angle_tol) = sharp_edge_angle {
                        let angle = compute_face_angle(fi, adj, shell);
                        if angle < angle_tol {
                            continue;
                        }
                    }
                    queue.push(adj);
                }
            }
        }
    }

    selected.into_iter().collect()
}

/// Select the N largest faces by area.
fn select_faces_by_area(shell: &Shell, brep: &rcad_kernel::BRep, count: usize) -> Vec<usize> {
    let mut face_areas: Vec<(usize, f64)> = shell
        .faces
        .iter()
        .enumerate()
        .map(|(i, face)| (i, compute_face_area(face, brep)))
        .collect();

    face_areas.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    face_areas.into_iter().take(count).map(|(i, _)| i).collect()
}

/// Select faces with normal matching the given direction.
fn select_faces_by_normal(shell: &Shell, direction: &DVec3, angle_tolerance: f64) -> Vec<usize> {
    shell
        .faces
        .iter()
        .enumerate()
        .filter(|(_, face)| {
            let angle = face.normal.dot(*direction).acos();
            angle < angle_tolerance
        })
        .map(|(i, _)| i)
        .collect()
}

/// Select planar faces.
fn select_planar_faces(shell: &Shell, brep: &rcad_kernel::BRep) -> Vec<usize> {
    shell
        .faces
        .iter()
        .enumerate()
        .filter(|(fi, _)| {
            brep_face_surface(brep, *fi).is_some_and(|s| matches!(s, Surface3::Plane(_)))
        })
        .map(|(i, _)| i)
        .collect()
}

/// Select non-planar faces.
fn select_non_planar_faces(shell: &Shell, brep: &rcad_kernel::BRep) -> Vec<usize> {
    shell
        .faces
        .iter()
        .enumerate()
        .filter(|(fi, _)| {
            !brep_face_surface(brep, *fi).is_some_and(|s| matches!(s, Surface3::Plane(_)))
        })
        .map(|(i, _)| i)
        .collect()
}

/// Compute the area of a face (approximate from vertex loop).
fn compute_face_area(face: &Face, brep: &rcad_kernel::BRep) -> f64 {
    let vertices: Vec<DVec3> = face
        .outer_wire
        .edges
        .iter()
        .filter_map(|we| brep_edge_endpoints(brep, we.idx))
        .filter_map(|(start, _end)| brep_vertex_point(brep, start))
        .collect();

    if vertices.len() < 3 {
        return 0.0;
    }

    // Compute area using shoelace formula in 3D
    // Project to the plane of the face
    let normal = face.normal;
    let mut area = 0.0;

    for i in 0..vertices.len() {
        let j = (i + 1) % vertices.len();
        let cross = (vertices[i] - vertices[0]).cross(vertices[j] - vertices[0]);
        area += cross.dot(normal);
    }

    (area * 0.5).abs()
}

/// Compute the angle between two adjacent faces.
fn compute_face_angle(fi: usize, fj: usize, shell: &Shell) -> f64 {
    let face_i = &shell.faces[fi];
    let face_j = &shell.faces[fj];

    // Angle between normals
    let dot = face_i.normal.dot(face_j.normal);
    dot.acos()
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Lateral Face Options
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Options for lateral face creation.
#[derive(Debug, Clone)]
pub struct LateralFaceOptions {
    /// Whether to create lateral faces.
    pub create: bool,
    /// Whether to ensure tangency with adjacent faces.
    pub ensure_tangency: bool,
    /// Whether to split lateral faces at sharp edges.
    pub split_at_sharp_edges: bool,
    /// Angle threshold for sharp edges (radians).
    pub sharp_edge_angle: f64,
    /// Maximum aspect ratio for lateral faces before splitting.
    pub max_aspect_ratio: Option<f64>,
    /// Whether to merge coplanar lateral faces.
    pub merge_coplanar: bool,
}

impl Default for LateralFaceOptions {
    fn default() -> Self {
        Self {
            create: true,
            ensure_tangency: false,
            split_at_sharp_edges: false,
            sharp_edge_angle: std::f64::consts::PI / 4.0, // 45 degrees
            max_aspect_ratio: None,
            merge_coplanar: false,
        }
    }
}

impl LateralFaceOptions {
    /// Create default lateral face options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Disable lateral face creation.
    pub fn none() -> Self {
        Self {
            create: false,
            ..Self::default()
        }
    }

    /// Enable tangency with adjacent faces.
    pub fn with_tangency(mut self) -> Self {
        self.ensure_tangency = true;
        self
    }

    /// Enable splitting at sharp edges.
    pub fn with_splitting(mut self, angle: f64) -> Self {
        self.split_at_sharp_edges = true;
        self.sharp_edge_angle = angle;
        self
    }

    /// Set maximum aspect ratio.
    pub fn with_max_aspect_ratio(mut self, ratio: f64) -> Self {
        self.max_aspect_ratio = Some(ratio);
        self
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Thickness Variation
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Thickness specification for a region or face.
#[derive(Debug, Clone)]
pub struct ThicknessSpec {
    /// Base thickness value.
    pub base: f64,
    /// Optional minimum thickness for auto-reduction.
    pub min: Option<f64>,
    /// Optional maximum thickness.
    pub max: Option<f64>,
    /// Transition mode to adjacent regions.
    pub transition: ThicknessTransition,
}

impl ThicknessSpec {
    /// Create a uniform thickness specification.
    pub fn uniform(value: f64) -> Self {
        Self {
            base: value,
            min: None,
            max: None,
            transition: ThicknessTransition::Sharp,
        }
    }

    /// Create a thickness specification with limits.
    pub fn with_limits(value: f64, min: f64, max: f64) -> Self {
        Self {
            base: value,
            min: Some(min),
            max: Some(max),
            transition: ThicknessTransition::Sharp,
        }
    }

    /// Set the transition mode.
    pub fn with_transition(mut self, transition: ThicknessTransition) -> Self {
        self.transition = transition;
        self
    }
}

/// Transition mode between thickness regions.
#[derive(Debug, Clone, Copy)]
pub enum ThicknessTransition {
    /// Sharp transition (default).
    Sharp,
    /// Linear interpolation over a distance.
    Linear {
        /// Distance over which to interpolate.
        distance: f64,
    },
    /// Smooth (cubic) interpolation.
    Smooth {
        /// Distance over which to interpolate.
        distance: f64,
    },
}

/// Thickness specification by face region.
#[derive(Debug, Clone)]
pub struct VariableThickness {
    /// Face-specific thickness values.
    pub face_thicknesses: Vec<(usize, f64)>,
    /// Default thickness for unspecified faces.
    pub default_thickness: f64,
    /// Transition mode between regions.
    pub transition: ThicknessTransition,
}

impl VariableThickness {
    /// Create a variable thickness specification.
    pub fn new(face_thicknesses: Vec<(usize, f64)>, default: f64) -> Self {
        Self {
            face_thicknesses,
            default_thickness: default,
            transition: ThicknessTransition::Sharp,
        }
    }

    /// Get thickness for a specific face.
    pub fn thickness_for_face(&self, face_index: usize) -> f64 {
        self.face_thicknesses
            .iter()
            .find(|(i, _)| *i == face_index)
            .map(|&(_, t)| t)
            .unwrap_or(self.default_thickness)
    }

    /// Set the transition mode.
    pub fn with_transition(mut self, transition: ThicknessTransition) -> Self {
        self.transition = transition;
        self
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Self-Intersection Handling
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Options for self-intersection handling.
#[derive(Debug, Clone)]
pub struct SelfIntersectionOptions {
    /// Whether to check for self-intersection.
    pub check: bool,
    /// Whether to automatically reduce thickness to avoid self-intersection.
    pub auto_reduce: bool,
    /// Minimum thickness after auto-reduction.
    pub min_thickness: f64,
    /// Warning threshold (fraction of min distance).
    pub warning_threshold: f64,
    /// Whether to abort on self-intersection.
    pub abort_on_detection: bool,
}

impl Default for SelfIntersectionOptions {
    fn default() -> Self {
        Self {
            check: true,
            auto_reduce: false,
            min_thickness: 0.01,
            warning_threshold: 0.8,
            abort_on_detection: false,
        }
    }
}

impl SelfIntersectionOptions {
    /// Create default options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Enable auto-reduction.
    pub fn with_auto_reduce(mut self, min_thickness: f64) -> Self {
        self.auto_reduce = true;
        self.min_thickness = min_thickness;
        self
    }

    /// Set warning threshold.
    pub fn with_warning_threshold(mut self, threshold: f64) -> Self {
        self.warning_threshold = threshold;
        self
    }

    /// Abort on self-intersection detection.
    pub fn abort_on_detection(mut self) -> Self {
        self.abort_on_detection = true;
        self
    }
}

/// Result of self-intersection analysis.
#[derive(Debug, Clone)]
pub struct SelfIntersectionAnalysis {
    /// Whether self-intersection would occur.
    pub would_intersect: bool,
    /// Minimum distance between non-adjacent faces.
    pub min_distance: f64,
    /// Safe thickness (half of min distance).
    pub safe_thickness: f64,
    /// Recommended thickness if auto-reduce is enabled.
    pub recommended_thickness: Option<f64>,
}

/// Analyze potential self-intersection for a given thickness.
pub fn analyze_self_intersection(
    brep: &rcad_kernel::BRep,
    thickness: f64,
) -> SelfIntersectionAnalysis {
    let shell = match shell_from_brep(brep) {
        Some(s) => s,
        None => {
            return SelfIntersectionAnalysis {
                would_intersect: false,
                min_distance: f64::MAX,
                safe_thickness: f64::MAX,
                recommended_thickness: Some(thickness),
            };
        }
    };

    // Compute face centroids
    let centroids: Vec<DVec3> = shell
        .faces
        .iter()
        .map(|face| compute_face_centroid(face, brep))
        .collect();

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            // Check if faces share an edge (adjacent)
            let share_edge = shell.faces[i].outer_wire.edges.iter().any(|we_i| {
                shell.faces[j]
                    .outer_wire
                    .edges
                    .iter()
                    .any(|we_j| we_i.idx == we_j.idx)
            });
            if share_edge {
                continue;
            }
            let dist = (centroids[i] - centroids[j]).length();
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    let safe_thickness = min_dist * 0.5;
    let would_intersect = thickness.abs() > safe_thickness;

    SelfIntersectionAnalysis {
        would_intersect,
        min_distance: min_dist,
        safe_thickness,
        recommended_thickness: if would_intersect {
            Some(safe_thickness * 0.95) // 5% safety margin
        } else {
            Some(thickness)
        },
    }
}

/// Compute the centroid of a face.
fn compute_face_centroid(face: &Face, brep: &rcad_kernel::BRep) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0;
    for we in &face.outer_wire.edges {
        if let Some((st, _en)) = brep_edge_endpoints(brep, we.idx) {
            if let Some(pt) = brep_vertex_point(brep, st) {
                sum += pt;
                count += 1;
            }
        }
    }
    if count > 0 {
        sum / count as f64
    } else {
        DVec3::ZERO
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Thickening Options
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Comprehensive options for thickening operations.
#[derive(Debug, Clone)]
pub struct ThickeningOptions {
    /// Thickness value (positive = outward, negative = inward).
    pub thickness: f64,
    /// Face selection strategy.
    pub face_selection: FaceSelectionStrategy,
    /// Lateral face options.
    pub lateral_faces: LateralFaceOptions,
    /// Variable thickness specification (optional).
    pub variable_thickness: Option<VariableThickness>,
    /// Self-intersection handling options.
    pub self_intersection: SelfIntersectionOptions,
    /// Geometric tolerance.
    pub tolerance: f64,
}

impl Default for ThickeningOptions {
    fn default() -> Self {
        Self {
            thickness: 0.1,
            face_selection: FaceSelectionStrategy::default(),
            lateral_faces: LateralFaceOptions::default(),
            variable_thickness: None,
            self_intersection: SelfIntersectionOptions::default(),
            tolerance: TOLERANCE_ABS,
        }
    }
}

impl ThickeningOptions {
    /// Create options with a given thickness.
    pub fn new(thickness: f64) -> Self {
        Self {
            thickness,
            ..Self::default()
        }
    }

    /// Set face selection strategy.
    pub fn with_face_selection(mut self, strategy: FaceSelectionStrategy) -> Self {
        self.face_selection = strategy;
        self
    }

    /// Set lateral face options.
    pub fn with_lateral_faces(mut self, options: LateralFaceOptions) -> Self {
        self.lateral_faces = options;
        self
    }

    /// Set variable thickness.
    pub fn with_variable_thickness(mut self, thickness: VariableThickness) -> Self {
        self.variable_thickness = Some(thickness);
        self
    }

    /// Set self-intersection options.
    pub fn with_self_intersection(mut self, options: SelfIntersectionOptions) -> Self {
        self.self_intersection = options;
        self
    }

    /// Set tolerance.
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Inline BRep builder helpers (avoids rcad_modeling dependency)
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

fn add_vertex(brep: &mut rcad_kernel::BRep, point: DVec3) -> usize {
    brep.add_tvertex(point).index
}

fn add_edge(
    brep: &mut rcad_kernel::BRep,
    curve: Curve3,
    t0: f64,
    t1: f64,
    v0: usize,
    v1: usize,
) -> usize {
    brep.add_edge_flat(v0, v1, Some(curve), [t0, t1])
}

fn add_face(
    brep: &mut rcad_kernel::BRep,
    surface: Surface3,
    outer: Wire,
    inner: Vec<Wire>,
    shell_sr: ShapeRef,
) -> usize {
    let outer_edge_refs: Vec<ShapeRef> = outer
        .edges
        .iter()
        .map(|we| {
            let orientation = if we.forward {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            ShapeRef::synthetic_with_orientation(we.idx, orientation)
        })
        .collect();
    let outer_wire = brep.add_twire(outer_edge_refs);
    let inner_wires: Vec<ShapeRef> = inner
        .iter()
        .map(|w| {
            let edge_refs: Vec<ShapeRef> = w
                .edges
                .iter()
                .map(|we| {
                    let orientation = if we.forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                    ShapeRef::synthetic_with_orientation(we.idx, orientation)
                })
                .collect();
            brep.add_twire(edge_refs)
        })
        .collect();
    let face_sr = brep.add_tface(
        Some(surface),
        outer_wire,
        inner_wires,
        None,
        None,
        vec![],
        true,
    );
    brep.shell_mut(shell_sr).faces.push(face_sr);
    brep.shell_mut(shell_sr).my_shapes.push(face_sr);
    face_sr.index
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Surface offset
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

fn offset_surface(surf: &Surface3, d: f64) -> Option<Surface3> {
    use rcad_kernel::geom::*;
    match surf {
        Surface3::Plane(p) => Some(Surface3::Plane(Plane::new(
            p.origin + p.normal * d,
            p.normal,
        ))),
        Surface3::Sphere(s) => {
            let r = s.radius + d;
            if r <= TOLERANCE_ABS {
                return None;
            }
            Some(Surface3::Sphere(SphericalSurface::new(s.center, s.axis, r)))
        }
        Surface3::Cylinder(c) => {
            let r = c.radius + d;
            if r <= TOLERANCE_ABS {
                return None;
            }
            Some(Surface3::Cylinder(CylindricalSurface {
                origin: c.origin,
                axis: c.axis,
                radius: r,
                ref_dir: c.ref_dir,
            }))
        }
        Surface3::Cone(c) => {
            let sin_a = c.half_angle_rad.sin();
            let shift = if sin_a.abs() > TOLERANCE_LINEAR_ULTRA_STRICT {
                d / sin_a
            } else {
                d
            };
            let new_r = c.radius + d;
            if new_r <= TOLERANCE_ABS {
                return None;
            }
            Some(Surface3::Cone(ConicalSurface {
                apex: c.apex - c.axis * shift,
                axis: c.axis,
                radius: new_r,
                half_angle_rad: c.half_angle_rad,
            }))
        }
        Surface3::Torus(t) => {
            let r = t.minor_radius + d;
            if r <= TOLERANCE_ABS {
                return None;
            }
            Some(Surface3::Torus(ToroidalSurface {
                center: t.center,
                axis: t.axis,
                major_radius: t.major_radius,
                minor_radius: r,
            }))
        }
        _ => None,
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Vertex normals
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

fn vertex_normal(shell: &Shell, brep: &rcad_kernel::BRep, vidx: usize) -> DVec3 {
    let mut n = DVec3::ZERO;
    let mut count = 0;
    for face in &shell.faces {
        let uses = face.outer_wire.edges.iter().any(|we| {
            let ed = match &*brep.tshapes[we.idx] {
                TShape::Edge(e) => e,
                _ => return false,
            };
            ed.first.index == vidx || ed.last.index == vidx
        });
        if uses {
            n += face.normal;
            count += 1;
        }
    }
    if count > 0 {
        (n / count as f64).normalize_or(DVec3::Z)
    } else {
        DVec3::Z
    }
}

fn vertex_normal_with_thickness(
    shell: &Shell,
    brep: &rcad_kernel::BRep,
    vidx: usize,
    thickness_map: &HashMap<usize, f64>,
    default_thickness: f64,
) -> DVec3 {
    let mut weighted_normal = DVec3::ZERO;
    let mut total_weight = 0.0;

    for (fi, face) in shell.faces.iter().enumerate() {
        let uses = face.outer_wire.edges.iter().any(|we| {
            let ed = match &*brep.tshapes[we.idx] {
                TShape::Edge(e) => e,
                _ => return false,
            };
            ed.first.index == vidx || ed.last.index == vidx
        });

        if uses {
            let thickness = thickness_map.get(&fi).copied().unwrap_or(default_thickness);
            let weight = thickness.abs();
            weighted_normal += face.normal * weight;
            total_weight += weight;
        }
    }

    if total_weight > 0.0 {
        (weighted_normal / total_weight).normalize_or(DVec3::Z)
    } else {
        DVec3::Z
    }
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Edge chaining
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

fn chain_edges(edge_indices: &[usize], edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    if edge_indices.is_empty() {
        return vec![];
    }
    let mut remaining: HashSet<usize> = edge_indices.iter().copied().collect();
    let mut loops = Vec::new();

    while let Some(&start_idx) = remaining.iter().next() {
        remaining.remove(&start_idx);
        let mut chain = vec![start_idx];
        let mut current_end = edges[start_idx].1;

        loop {
            let next = remaining
                .iter()
                .find(|&&ei| edges[ei].0 == current_end || edges[ei].1 == current_end)
                .copied();
            match next {
                Some(ei) => {
                    remaining.remove(&ei);
                    chain.push(ei);
                    let (e_st, e_en) = edges[ei];
                    current_end = if e_st == current_end { e_en } else { e_st };
                }
                None => break,
            }
        }
        if chain.len() >= 2 {
            loops.push(chain);
        }
    }
    loops
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Public API
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?

/// Thicken a solid using comprehensive options.
///
/// This is the main entry point for thickening operations, supporting:
/// - Multiple face selection strategies
/// - Lateral face configuration
/// - Variable thickness
/// - Self-intersection handling
pub fn thicken_solid(
    brep: &rcad_kernel::BRep,
    options: &ThickeningOptions,
) -> Option<ThickeningResult> {
    let shell = shell_from_brep(brep)?;
    if shell.faces.is_empty() {
        return None;
    }

    // Select faces to remove
    let removed_face_indices = select_faces_for_removal(brep, &options.face_selection);
    let removed_set: HashSet<usize> = removed_face_indices.iter().copied().collect();

    if removed_set.len() >= shell.faces.len() {
        return None; // can't remove all faces
    }

    // Determine thickness
    let mut thickness = options.thickness;
    let mut warnings = Vec::new();
    let mut actual_thickness: Option<f64> = None;

    // Handle variable thickness
    let face_thicknesses: HashMap<usize, f64> = if let Some(ref var) = options.variable_thickness {
        var.face_thicknesses.iter().map(|&(i, t)| (i, t)).collect()
    } else {
        HashMap::new()
    };

    // Self-intersection analysis
    if options.self_intersection.check && removed_face_indices.is_empty() {
        let analysis = analyze_self_intersection(brep, thickness);

        if analysis.would_intersect {
            if options.self_intersection.auto_reduce {
                if let Some(recommended) = analysis.recommended_thickness
                    && recommended >= options.self_intersection.min_thickness
                {
                    warnings.push(ThickeningWarning::ThicknessAutoReduced {
                        original: thickness,
                        reduced: recommended,
                    });
                    actual_thickness = Some(recommended);
                    thickness = recommended;
                }
            } else if options.self_intersection.abort_on_detection {
                return None;
            } else {
                warnings.push(ThickeningWarning::SelfIntersection {
                    min_distance: analysis.min_distance,
                    requested_thickness: thickness,
                });
            }
        }

        // Check warning threshold
        if thickness.abs() > analysis.safe_thickness * options.self_intersection.warning_threshold {
            warnings.push(ThickeningWarning::ThinRegionDetected {
                center: compute_face_centroid(&shell.faces[0], brep),
                thickness: analysis.safe_thickness * 2.0,
            });
        }
    }

    let d = thickness;

    // Build the "kept" shell (original faces minus removed)
    let kept_faces: Vec<(usize, &Face)> = shell
        .faces
        .iter()
        .enumerate()
        .filter(|(i, _)| !removed_set.contains(i))
        .collect();

    if kept_faces.is_empty() {
        return None;
    }

    // Find boundary edges of the kept shell
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for (_, face) in &kept_faces {
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }
    let boundary_edges: Vec<usize> = edge_use
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(idx, _)| idx)
        .collect();

    // Compute offset vertex positions
    let kept_shell = Shell {
        faces: kept_faces.iter().map(|(_, f)| (*f).clone()).collect(),
    };

    let new_pts: Vec<DVec3> = if let Some(ref var) = options.variable_thickness {
        // Variable thickness: use weighted normals
        (0..brep.vertex_count())
            .map(|i| {
                let pt = brep.vertex_point(i).unwrap();
                let n = vertex_normal_with_thickness(
                    &kept_shell,
                    brep,
                    i,
                    &face_thicknesses,
                    var.default_thickness,
                );
                let face_thickness = find_dominant_face_thickness(
                    i,
                    &kept_shell,
                    brep,
                    &face_thicknesses,
                    var.default_thickness,
                );
                pt + n * face_thickness
            })
            .collect()
    } else {
        (0..brep.vertex_count())
            .map(|i| {
                let pt = brep.vertex_point(i).unwrap();
                let n = vertex_normal(&kept_shell, brep, i);
                pt + n * d
            })
            .collect()
    };

    // Build result BRep with topods API
    let mut out = rcad_kernel::BRep::new();
    let out_shell = out.add_tshell(vec![]);
    let _out_solid = out.add_tsolid(vec![out_shell]);

    let mut orig_vidx: Vec<usize> = Vec::new();
    for i in 0..brep.vertex_count() {
        orig_vidx.push(add_vertex(&mut out, brep.vertex_point(i).unwrap()));
    }
    let mut off_vidx: Vec<usize> = Vec::new();
    for &p in &new_pts {
        off_vidx.push(add_vertex(&mut out, p));
    }

    // Offset kept faces
    let mut offset_face_count = 0;
    for &(fi, face) in &kept_faces {
        let off_surf = match brep_face_surface(brep, fi) {
            Some(s) => {
                let face_d = face_thicknesses.get(&fi).copied().unwrap_or(d);
                match offset_surface(&s, face_d) {
                    Some(os) => os,
                    None => {
                        warnings.push(ThickeningWarning::DegenerateSurface {
                            face_index: fi,
                            surface_type: format!("{:?}", s),
                        });
                        continue;
                    }
                }
            }
            None => continue,
        };

        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let ed = match &*brep.tshapes[we.idx] {
                TShape::Edge(e) => e,
                _ => continue,
            };
            let vs = off_vidx[ed.first.index];
            let ve = off_vidx[ed.last.index];
            let pvs = out.vertex_point(vs).unwrap();
            let pve = out.vertex_point(ve).unwrap();
            let dir = (pve - pvs).normalize_or(DVec3::X);
            let len = (pve - pvs).length();
            let curve = Curve3::Line(Line3 {
                origin: pvs,
                direction: dir,
            });
            let eidx = add_edge(&mut out, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(
            &mut out,
            off_surf,
            Wire { edges: wire_edges },
            Vec::new(),
            out_shell,
        );
        offset_face_count += 1;
    }

    if offset_face_count == 0 {
        return None;
    }

    // Create lateral faces
    let mut lateral_count = 0;
    if options.lateral_faces.create {
        // Build flat edge data from topods BRep
        let edges_flat: Vec<(usize, usize)> = (0..brep.edge_count())
            .map(|ei| match &*brep.tshapes[ei] {
                TShape::Edge(ed) => (ed.first.index, ed.last.index),
                _ => (usize::MAX, usize::MAX),
            })
            .collect();
        lateral_count = create_lateral_faces(
            &mut out,
            &boundary_edges,
            &edges_flat,
            &orig_vidx,
            &off_vidx,
            options,
            out_shell,
        );
    }

    // Triangulate
    mesh_brep(&mut out, &TessellationParams::default());

    // Final self-intersection check
    let self_intersection = if options.self_intersection.check
        && boundary_edges.is_empty()
        && removed_face_indices.is_empty()
    {
        detect_self_intersection(brep, thickness)
    } else {
        false
    };

    Some(ThickeningResult {
        brep: out,
        offset_faces: offset_face_count,
        lateral_faces: lateral_count,
        self_intersection,
        warnings,
        actual_thickness,
    })
}

/// Find the dominant face thickness for a vertex.
fn find_dominant_face_thickness(
    vidx: usize,
    shell: &Shell,
    brep: &rcad_kernel::BRep,
    thickness_map: &HashMap<usize, f64>,
    default: f64,
) -> f64 {
    let mut max_thickness = default;
    for (fi, face) in shell.faces.iter().enumerate() {
        let uses = face.outer_wire.edges.iter().any(|we| {
            let ed = match &*brep.tshapes[we.idx] {
                TShape::Edge(e) => e,
                _ => return false,
            };
            ed.first.index == vidx || ed.last.index == vidx
        });
        if uses {
            let t = thickness_map.get(&fi).copied().unwrap_or(default);
            if t.abs() > max_thickness.abs() {
                max_thickness = t;
            }
        }
    }
    max_thickness
}

/// Create lateral faces along boundary edges.
fn create_lateral_faces(
    out: &mut rcad_kernel::BRep,
    boundary_edges: &[usize],
    edges: &[(usize, usize)],
    orig_vidx: &[usize],
    off_vidx: &[usize],
    options: &ThickeningOptions,
    shell_sr: ShapeRef,
) -> usize {
    let loops = chain_edges(boundary_edges, edges);
    let mut lateral_count = 0;

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let (e_start, e_end) = edges[eidx];
            let o_vs = orig_vidx[e_start];
            let o_ve = orig_vidx[e_end];
            let f_vs = off_vidx[e_start];
            let f_ve = off_vidx[e_end];

            let p0 = out.vertex_point(o_vs).unwrap();
            let p1 = out.vertex_point(o_ve).unwrap();
            let p3 = out.vertex_point(f_vs).unwrap();

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
                continue;
            }

            // Check for splitting if enabled
            let should_split = if options.lateral_faces.split_at_sharp_edges {
                let edge_len = (p1 - p0).length();
                let thickness_len = (p3 - p0).length();
                let aspect_ratio = edge_len / thickness_len.max(TOLERANCE_LINEAR_ULTRA_STRICT);
                options
                    .lateral_faces
                    .max_aspect_ratio
                    .is_some_and(|max_ratio| aspect_ratio > max_ratio)
            } else {
                false
            };

            if should_split {
                // Split into two lateral faces
                let mid_orig = (p0 + p1) * 0.5;
                let mid_off = (p3 + out.vertex_point(f_ve).unwrap()) * 0.5;

                let mid_orig_vidx = add_vertex(out, mid_orig);
                let mid_off_vidx = add_vertex(out, mid_off);

                // First half
                lateral_count += create_single_lateral_face(
                    out,
                    o_vs,
                    mid_orig_vidx,
                    mid_off_vidx,
                    f_vs,
                    normal,
                    shell_sr,
                );
                // Second half
                lateral_count += create_single_lateral_face(
                    out,
                    mid_orig_vidx,
                    o_ve,
                    f_ve,
                    mid_off_vidx,
                    normal,
                    shell_sr,
                );
            } else {
                lateral_count +=
                    create_single_lateral_face(out, o_vs, o_ve, f_ve, f_vs, normal, shell_sr);
            }
        }
    }

    lateral_count
}

/// Create a single lateral face (quad).
fn create_single_lateral_face(
    out: &mut rcad_kernel::BRep,
    v0: usize,
    v1: usize,
    v2: usize,
    v3: usize,
    normal: DVec3,
    shell_sr: ShapeRef,
) -> usize {
    let p0 = out.vertex_point(v0).unwrap();

    let surf = Surface3::Plane(rcad_kernel::geom::Plane::new(
        p0,
        normal.normalize_or_zero(),
    ));

    let vseq = [v0, v1, v2, v3];
    let mut edges = Vec::new();
    for i in 0..4 {
        let s = vseq[i];
        let en = vseq[(i + 1) % 4];
        let ps = out.vertex_point(s).unwrap();
        let pe = out.vertex_point(en).unwrap();
        let dir = (pe - ps).normalize_or(DVec3::X);
        let len = (pe - ps).length();
        let curve = Curve3::Line(Line3 {
            origin: ps,
            direction: dir,
        });
        edges.push(WireEdge::fwd(add_edge(out, curve, 0.0, len, s, en)));
    }

    add_face(out, surf, Wire { edges }, Vec::new(), shell_sr);
    1
}

/// Thicken a solid by removing specified faces, offsetting the remaining
/// faces, and building lateral ruled faces at the removed-face boundaries.
///
/// This is analogous to OCCT `BRepOffsetAPI_MakeThickSolid`.
///
/// - `brep`: input solid (must have at least one shell with geometry).
/// - `removed_face_indices`: indices of faces to remove (relative to
/// `brep.solids[0].shells[0].faces`).
/// - `thickness`: positive = inward (material removed), negative = outward.
///
/// Returns `None` if all faces are removed, thickness is zero, or the offset
/// would create degenerate surfaces.
pub fn thick_solid_with_removed_faces(
    brep: &rcad_kernel::BRep,
    removed_face_indices: &[usize],
    thickness: f64,
) -> Option<ThickeningResult> {
    let options = ThickeningOptions::new(thickness).with_face_selection(
        FaceSelectionStrategy::explicit(removed_face_indices.to_vec()),
    );
    thicken_solid(brep, &options)
}

/// Detect self-intersection for closed-shell inward offsetting.
///
/// Computes the minimum distance between non-adjacent face centroids.
/// If `thickness > min_distance / 2`, the offset faces will self-intersect.
fn detect_self_intersection(brep: &rcad_kernel::BRep, thickness: f64) -> bool {
    let shell = match shell_from_brep(brep) {
        Some(s) => s,
        None => return false,
    };

    // Compute face centroids
    let centroids: Vec<DVec3> = shell
        .faces
        .iter()
        .map(|face| compute_face_centroid(face, brep))
        .collect();

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            // Check if faces share an edge (adjacent)
            let share_edge = shell.faces[i].outer_wire.edges.iter().any(|we_i| {
                shell.faces[j]
                    .outer_wire
                    .edges
                    .iter()
                    .any(|we_j| we_i.idx == we_j.idx)
            });
            if share_edge {
                continue;
            }
            let dist = (centroids[i] - centroids[j]).length();
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        return false; // no non-adjacent faces
    }

    thickness.abs() > min_dist * 0.5
}

/// Thicken an open shell by offsetting faces along their normals and
/// filling the gaps with lateral ruled faces.
///
/// The input BRep must have at least one face with populated surface data
/// (e.g. created via `make_box_brep` which populates analytic surfaces).
///
/// `thickness` > 0 offsets outward, < 0 offsets inward.
/// Returns `None` if the shell is closed, has no geometry, or the offset
/// would create degenerate surfaces.
pub fn thicken_shell(brep: &rcad_kernel::BRep, thickness: f64) -> Option<ThickeningResult> {
    if thickness.abs() < TOLERANCE_LEN_MIN {
        return None;
    }

    let shell = shell_from_brep(brep)?;
    if shell.faces.is_empty() {
        return None;
    }

    let d = thickness;

    // Find boundary edges
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }
    let boundary_edges: Vec<usize> = edge_use
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(idx, _)| idx)
        .collect();

    // Compute offset vertex positions
    let new_pts: Vec<DVec3> = (0..brep.vertex_count())
        .map(|i| {
            let pt = brep.vertex_point(i).unwrap();
            let n = vertex_normal(&shell, brep, i);
            pt + n * d
        })
        .collect();

    // Build result BRep with topods API
    let mut out = rcad_kernel::BRep::new();
    let out_shell = out.add_tshell(vec![]);
    let _out_solid = out.add_tsolid(vec![out_shell]);

    let mut orig_vidx: Vec<usize> = Vec::new();
    for i in 0..brep.vertex_count() {
        orig_vidx.push(add_vertex(&mut out, brep.vertex_point(i).unwrap()));
    }
    let mut off_vidx: Vec<usize> = Vec::new();
    for &p in &new_pts {
        off_vidx.push(add_vertex(&mut out, p));
    }

    // Offset faces
    let mut offset_face_count = 0;
    for (fi, _face) in shell.faces.iter().enumerate() {
        let off_surf = match brep_face_surface(brep, fi) {
            Some(s) => match offset_surface(&s, d) {
                Some(os) => os,
                None => continue,
            },
            None => continue,
        };
        let face = &shell.faces[fi];

        // Build wire from offset vertices
        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let ed = match &*brep.tshapes[we.idx] {
                TShape::Edge(e) => e,
                _ => continue,
            };
            let vs = off_vidx[ed.first.index];
            let ve = off_vidx[ed.last.index];
            let pvs = out.vertex_point(vs).unwrap();
            let pve = out.vertex_point(ve).unwrap();
            let dir = (pve - pvs).normalize_or(DVec3::X);
            let len = (pve - pvs).length();
            let curve = Curve3::Line(Line3 {
                origin: pvs,
                direction: dir,
            });
            let eidx = add_edge(&mut out, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(
            &mut out,
            off_surf,
            Wire { edges: wire_edges },
            Vec::new(),
            out_shell,
        );
        offset_face_count += 1;
    }

    if offset_face_count == 0 {
        return None;
    }

    // Lateral faces along boundary edges
    let mut lateral_count = 0;
    let edges_flat: Vec<(usize, usize)> = (0..brep.edge_count())
        .map(|ei| match &*brep.tshapes[ei] {
            TShape::Edge(ed) => (ed.first.index, ed.last.index),
            _ => (usize::MAX, usize::MAX),
        })
        .collect();
    let loops = chain_edges(&boundary_edges, &edges_flat);

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let (e_st, e_en) = edges_flat[eidx];
            let o_vs = orig_vidx[e_st];
            let o_ve = orig_vidx[e_en];
            let f_vs = off_vidx[e_st];
            let f_ve = off_vidx[e_en];

            let p0 = out.vertex_point(o_vs).unwrap();
            let p1 = out.vertex_point(o_ve).unwrap();
            let p3 = out.vertex_point(f_vs).unwrap();

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < TOLERANCE_LINEAR_ULTRA_STRICT {
                continue;
            }

            let surf = Surface3::Plane(rcad_kernel::geom::Plane::new(p0, p0.normalize_or_zero()));

            // Quad: orig_start -> orig_end -> off_end -> off_start
            let vseq = [o_vs, o_ve, f_ve, f_vs];
            let mut edges = Vec::new();
            for i in 0..4 {
                let s = vseq[i];
                let en = vseq[(i + 1) % 4];
                let ps = out.vertex_point(s).unwrap();
                let pe = out.vertex_point(en).unwrap();
                let dir = (pe - ps).normalize_or(DVec3::X);
                let len = (pe - ps).length();
                let curve = Curve3::Line(Line3 {
                    origin: ps,
                    direction: dir,
                });
                edges.push(WireEdge::fwd(add_edge(&mut out, curve, 0.0, len, s, en)));
            }

            add_face(&mut out, surf, Wire { edges }, Vec::new(), out_shell);
            lateral_count += 1;
        }
    }

    // Triangulate
    mesh_brep(&mut out, &TessellationParams::default());

    Some(ThickeningResult {
        brep: out,
        offset_faces: offset_face_count,
        lateral_faces: lateral_count,
        self_intersection: false,
        warnings: Vec::new(),
        actual_thickness: None,
    })
}

//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
// Tests
//  鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?鈧?
