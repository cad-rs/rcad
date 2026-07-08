// =============================================================================
// Shape Query Utilities
// =============================================================================

/// Determine the shape type of a BRep.
///
/// Returns the highest-level topological entity present:
/// - Compound if `brep.compound` is set
/// - CompSolid if `brep.compsolid` is set
/// - Solid if there are solids
/// - Shell if there are shells (no solids)
/// - etc.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::{get_shape_type, ShapeType};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// assert_eq!(get_shape_type(&brep), ShapeType::Solid);
///
/// let empty = BRep::new();
/// assert_eq!(get_shape_type(&empty), ShapeType::Empty);
/// ```
pub fn get_shape_type(brep: &rcad_kernel::BRep) -> ShapeType {
    if brep.compound.is_some() {
        return ShapeType::Compound;
    }
    if brep.compsolid.is_some() {
        return ShapeType::CompSolid;
    }
    if !brep.solids.is_empty() {
        // Check if solids have shells with faces
        let has_faces = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .any(|sh| !sh.faces.is_empty());
        if has_faces {
            return ShapeType::Solid;
        }
        // Check if there are empty shells
        let has_shells = brep.solids.iter().any(|s| !s.shells.is_empty());
        if has_shells {
            return ShapeType::Shell;
        }
        return ShapeType::Solid;
    }
    if !brep.edges.is_empty() {
        return ShapeType::Edge;
    }
    if !brep.vertices.is_empty() {
        return ShapeType::Vertex;
    }
    ShapeType::Empty
}

/// Get the outer wire of a face.
///
/// Returns a reference to the face's outer wire (boundary).
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Index of the face (flat index across all solids/shells)
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_outer_wire;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let outer_wire = get_outer_wire(&brep, 0).unwrap();
/// assert_eq!(outer_wire.edges.len(), 4); // Rectangle
/// ```
pub fn get_outer_wire(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<&Wire, BRepToolsError> {
    let (face, _) = get_face_by_flat_index(brep, face_idx)?;
    Ok(&face.outer_wire)
}

/// Get the inner wires (holes) of a face.
///
/// Returns references to the face's inner wires (holes/cutouts).
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Index of the face (flat index across all solids/shells)
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_inner_wires;
///
/// // A face with a hole would have inner_wires.len() > 0
/// let inner_wires = get_inner_wires(&brep, 0).unwrap();
/// for wire in inner_wires {
///     println!("Hole with {} edges", wire.edges.len());
/// }
/// ```
pub fn get_inner_wires(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<&[Wire], BRepToolsError> {
    let (face, _) = get_face_by_flat_index(brep, face_idx)?;
    Ok(&face.inner_wires)
}

/// Check if a shape is closed (forms a manifold solid).
///
/// A shape is closed if:
/// - It is a solid with a closed shell
/// - Each edge is shared by exactly two faces
/// - The shell encloses a finite volume
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::is_closed;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// assert!(is_closed(&brep));
///
/// let empty = BRep::new();
/// assert!(!is_closed(&empty));
/// ```
pub fn is_closed(brep: &rcad_kernel::BRep) -> bool {
    if brep.solids.is_empty() {
        return false;
    }

    for solid in &brep.solids {
        for shell in &solid.shells {
            if !is_shell_closed(brep, shell) {
                return false;
            }
        }
    }
    true
}

/// Check if a shell is closed by verifying edge manifoldness.
fn is_shell_closed(_brep: &rcad_kernel::BRep, shell: &Shell) -> bool {
    if shell.faces.is_empty() {
        return false;
    }

    // Count edge usage across all faces
    let mut edge_count = std::collections::HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_count.entry(we.idx).or_insert(0) += 1;
        }
        for inner in &face.inner_wires {
            for we in &inner.edges {
                *edge_count.entry(we.idx).or_insert(0) += 1;
            }
        }
    }

    // For a closed shell, each edge should appear exactly twice
    edge_count.values().all(|&count| count == 2)
}

// =============================================================================
// Geometry Utilities
// =============================================================================

/// Get the surface of a face.
///
/// Returns a reference to the surface geometry of the specified face.
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Flat index of the face across all solids/shells
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_surface;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let surface = get_surface(&brep, 0).unwrap();
/// // The surface of a box face is a plane
/// ```
pub fn get_surface(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<&Surface3, BRepToolsError> {
    let (_, flat_idx) = get_face_by_flat_index(brep, face_idx)?;

    match brep.geom.face_surface.get(flat_idx) {
        Some(Some(surf_idx)) => {
            brep.geom.surfaces.get(*surf_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "surface",
                    index: *surf_idx,
                })
        }
        Some(None) => Err(BRepToolsError::MissingGeometry {
            kind: "surface",
            index: face_idx,
        }),
        None => Err(BRepToolsError::InvalidIndex {
            kind: "face_surface",
            index: face_idx,
            max: brep.geom.face_surface.len(),
        }),
    }
}

/// Get the 3D curve of an edge.
///
/// Returns a reference to the 3D curve geometry of the specified edge.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge
/// * `edge_idx` - Index of the edge
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_curve;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let curve = get_curve(&brep, 0).unwrap();
/// ```
pub fn get_curve(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<&Curve3, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    match brep.geom.edge_curve.get(edge_idx) {
        Some(Some(curve_idx)) => {
            brep.geom.curves.get(*curve_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "curve",
                    index: *curve_idx,
                })
        }
        Some(None) => Err(BRepToolsError::MissingGeometry {
            kind: "curve",
            index: edge_idx,
        }),
        None => Err(BRepToolsError::InvalidIndex {
            kind: "edge_curve",
            index: edge_idx,
            max: brep.geom.edge_curve.len(),
        }),
    }
}

/// Get the parameter-space curve (pcurve) of an edge on a face.
///
/// Returns the 2D curve in the parameter space of the face's surface.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge and face
/// * `edge_idx` - Index of the edge
/// * `face_idx` - Flat index of the face
///
/// # Returns
///
/// A tuple containing:
/// - The 2D curve in the face's UV parameter space
/// - The surface index that the pcurve is defined on
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::get_pcurve;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// // Get pcurve of edge 0 on face 0
/// let (pcurve, surface_idx) = get_pcurve(&brep, 0, 0).unwrap();
/// ```
pub fn get_pcurve(brep: &rcad_kernel::BRep, edge_idx: usize, face_idx: usize) -> Result<(&Curve2d, usize), BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    let (_, _) = get_face_by_flat_index(brep, face_idx)?;

    // Get the surface index for this face
    let surf_idx = match brep.geom.face_surface.get(face_idx) {
        Some(Some(idx)) => *idx,
        _ => return Err(BRepToolsError::MissingGeometry {
            kind: "face_surface",
            index: face_idx,
        }),
    };

    // Find the pcurve for this edge on this surface
    let pcurves = brep.geom.edge_pcurves.get(edge_idx)
        .ok_or(BRepToolsError::MissingGeometry {
            kind: "edge_pcurves",
            index: edge_idx,
        })?;

    // Find the pcurve that matches this surface
    for pcurve in pcurves {
        if pcurve.surface_idx == surf_idx {
            let curve2d = brep.geom.curve2ds.get(pcurve.curve2d_idx)
                .ok_or(BRepToolsError::MissingGeometry {
                    kind: "curve2d",
                    index: pcurve.curve2d_idx,
                })?;
            return Ok((curve2d, surf_idx));
        }
    }

    // If no pcurve found, check if edge has single pcurve (common case)
    if pcurves.len() == 1 {
        let pcurve = &pcurves[0];
        let curve2d = brep.geom.curve2ds.get(pcurve.curve2d_idx)
            .ok_or(BRepToolsError::MissingGeometry {
                kind: "curve2d",
                index: pcurve.curve2d_idx,
            })?;
        return Ok((curve2d, pcurve.surface_idx));
    }

    Err(BRepToolsError::MissingGeometry {
        kind: "pcurve",
        index: edge_idx,
    })
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Get a face by its flat index (across all solids/shells).
///
/// Returns a tuple of (face reference, actual flat index used).
fn get_face_by_flat_index(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<(&Face, usize), BRepToolsError> {
    let mut current_idx = 0;

    for solid in &brep.solids {
        for shell in &solid.shells {
            if face_idx < current_idx + shell.faces.len() {
                let local_idx = face_idx - current_idx;
                return Ok((&shell.faces[local_idx], face_idx));
            }
            current_idx += shell.faces.len();
        }
    }

    Err(BRepToolsError::InvalidIndex {
        kind: "face",
        index: face_idx,
        max: current_idx,
    })
}

/// Get the parameter range of an edge's 3D curve.
///
/// Returns `[t_min, t_max]` if the edge has a curve with a defined range.
pub fn get_edge_range(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<Option<[f64; 2]>, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_curve_range.get(edge_idx).copied().flatten())
}

/// Check if an edge is degenerate (zero-length, like a pole).
pub fn is_edge_degenerate(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<bool, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false))
}

/// Get the tolerance of a vertex.
pub fn get_vertex_tolerance(brep: &rcad_kernel::BRep, vertex_idx: usize) -> Result<f64, BRepToolsError> {
    if vertex_idx >= brep.vertices.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "vertex",
            index: vertex_idx,
            max: brep.vertices.len(),
        });
    }

    Ok(brep.geom.vertex_tolerance.get(vertex_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

/// Get the tolerance of an edge.
pub fn get_edge_tolerance(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<f64, BRepToolsError> {
    if edge_idx >= brep.edges.len() {
        return Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.edges.len(),
        });
    }

    Ok(brep.geom.edge_tolerance.get(edge_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

/// Get the tolerance of a face.
pub fn get_face_tolerance(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<f64, BRepToolsError> {
    let (_, _) = get_face_by_flat_index(brep, face_idx)?;

    Ok(brep.geom.face_tolerance.get(face_idx).copied().unwrap_or(rcad_kernel::CONFUSION))
}

// =============================================================================
// Additional Shape Queries
// =============================================================================

/// Count the total number of faces in a BRep.
pub fn count_faces(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum()
}

/// Count the total number of edges in a BRep.
pub fn count_edges(brep: &rcad_kernel::BRep) -> usize {
    brep.edges.len()
}

/// Count the total number of vertices in a BRep.
pub fn count_vertices(brep: &rcad_kernel::BRep) -> usize {
    brep.vertices.len()
}

/// Count the total number of shells in a BRep.
pub fn count_shells(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.iter().map(|s| s.shells.len()).sum()
}

/// Count the total number of wires (outer + inner) across all faces in a BRep.
pub fn count_wires(brep: &rcad_kernel::BRep) -> usize {
    brep.solids.iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .map(|f| 1 + f.inner_wires.len())
        .sum()
}

/// Get the bounding box of a BRep.
///
/// Returns `[min_point, max_point]` or `None` if the BRep has no vertices.
pub fn bounding_box(brep: &rcad_kernel::BRep) -> Option<[DVec3; 2]> {
    if brep.vertices.is_empty() {
        return None;
    }

    let mut min_pt = brep.vertices[0].point;
    let mut max_pt = brep.vertices[0].point;

    for v in &brep.vertices[1..] {
        min_pt = min_pt.min(v.point);
        max_pt = max_pt.max(v.point);
    }

    Some([min_pt, max_pt])
}
