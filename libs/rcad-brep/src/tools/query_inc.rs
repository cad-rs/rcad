// =============================================================================
// Shape Query Utilities
// =============================================================================

/// Determine the shape type of a BRep.
///
/// Returns the highest-level topological entity present.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::{get_shape_type, ShapeType};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::new(); // topods
/// // empty => Empty
/// assert_eq!(get_shape_type(&brep), ShapeType::Empty);
/// ```
pub fn get_shape_type(brep: &rcad_kernel::BRep) -> ShapeType {
    for ts in &brep.tshapes {
        match &**ts {
            rcad_kernel::topods::TShape::Compound(_) => return ShapeType::Compound,
            rcad_kernel::topods::TShape::CompSolid(_) => return ShapeType::CompSolid,
            rcad_kernel::topods::TShape::Solid(sd) => {
                let has_faces = sd.shells.iter().any(|sr| {
                    if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                        !shd.faces.is_empty()
                    } else { false }
                });
                if has_faces { return ShapeType::Solid; }
                let has_shells = sd.shells.iter().any(|sr| {
                    matches!(&*brep.tshapes[sr.index], rcad_kernel::topods::TShape::Shell(_))
                });
                if has_shells { return ShapeType::Shell; }
                return ShapeType::Solid;
            }
            rcad_kernel::topods::TShape::Edge(_) => return ShapeType::Edge,
            rcad_kernel::topods::TShape::Vertex(_) => return ShapeType::Vertex,
            _ => {}
        }
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
pub fn get_outer_wire(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<&rcad_kernel::topods::TWireData, BRepToolsError> {
    let (_, _, fd) = get_face_by_flat_index(brep, face_idx)?;
    match &*brep.tshapes[fd.outer_wire.index] {
        rcad_kernel::topods::TShape::Wire(wd) => Ok(wd),
        _ => Err(BRepToolsError::InvalidIndex { kind: "wire", index: fd.outer_wire.index, max: brep.tshapes.len() }),
    }
}

/// Get the inner wires (holes) of a face.
///
/// Returns references to the face's inner wires (holes/cutouts).
///
/// # Arguments
///
/// * `brep` - The BRep containing the face
/// * `face_idx` - Index of the face (flat index across all solids/shells)
pub fn get_inner_wires(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<Vec<&rcad_kernel::topods::TWireData>, BRepToolsError> {
    let (_, _, fd) = get_face_by_flat_index(brep, face_idx)?;
    let mut wires = Vec::new();
    for iw_sr in &fd.inner_wires {
        match &*brep.tshapes[iw_sr.index] {
            rcad_kernel::topods::TShape::Wire(wd) => wires.push(wd),
            _ => return Err(BRepToolsError::InvalidIndex { kind: "wire", index: iw_sr.index, max: brep.tshapes.len() }),
        }
    }
    Ok(wires)
}

/// Check if a shape is closed (forms a manifold solid).
///
/// A shape is closed if:
/// - It is a solid with a closed shell
/// - Each edge is shared by exactly two faces
/// - The shell encloses a finite volume
pub fn is_closed(brep: &rcad_kernel::BRep) -> bool {
    let mut found_any = false;
    for ts in &brep.tshapes {
        if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    if !is_shell_closed(brep, shd) {
                        return false;
                    }
                    found_any = true;
                }
            }
        }
    }
    // Also check standalone shells
    for ts in &brep.tshapes {
        if let rcad_kernel::topods::TShape::Shell(shd) = &**ts {
            if !is_shell_closed(brep, shd) {
                return false;
            }
            found_any = true;
        }
    }
    found_any
}

/// Check if a shell is closed by verifying edge manifoldness.
fn is_shell_closed(brep: &rcad_kernel::BRep, shd: &rcad_kernel::topods::TShellData) -> bool {
    if shd.faces.is_empty() {
        return false;
    }

    // Count edge usage across all faces
    let mut edge_count = std::collections::HashMap::new();
    for face_sr in &shd.faces {
        if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
            if let rcad_kernel::topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for e_sr in &wd.edges {
                    *edge_count.entry(e_sr.index).or_insert(0) += 1;
                }
            }
            for iw_sr in &fd.inner_wires {
                if let rcad_kernel::topods::TShape::Wire(wd) = &*brep.tshapes[iw_sr.index] {
                    for e_sr in &wd.edges {
                        *edge_count.entry(e_sr.index).or_insert(0) += 1;
                    }
                }
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
pub fn get_surface(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<&Surface3, BRepToolsError> {
    let (_, _, fd) = get_face_by_flat_index(brep, face_idx)?;
    fd.surface.as_ref().ok_or(BRepToolsError::MissingGeometry {
        kind: "surface",
        index: face_idx,
    })
}

/// Get the 3D curve of an edge.
///
/// Returns a reference to the 3D curve geometry of the specified edge.
///
/// # Arguments
///
/// * `brep` - The BRep containing the edge
/// * `edge_idx` - Index of the edge (tshape index)
pub fn get_curve(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<&Curve3, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        rcad_kernel::topods::TShape::Edge(ed) => ed.curve.as_ref().ok_or(BRepToolsError::MissingGeometry {
            kind: "curve",
            index: edge_idx,
        }),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
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
/// * `edge_idx` - Index of the edge (tshape index)
/// * `face_idx` - Flat index of the face
///
/// # Returns
///
/// A tuple containing:
/// - The 2D curve in the face's UV parameter space
/// - The surface index that the pcurve is defined on
pub fn get_pcurve(brep: &rcad_kernel::BRep, edge_idx: usize, face_idx: usize) -> Result<(&Curve2d, usize), BRepToolsError> {
    let (face_ts_idx, _, fd) = get_face_by_flat_index(brep, face_idx)?;
    let _surf = fd.surface.as_ref().ok_or(BRepToolsError::MissingGeometry {
        kind: "surface",
        index: face_idx,
    })?;

    match &*brep.tshapes[edge_idx] {
        rcad_kernel::topods::TShape::Edge(ed) => {
            // Check pcurves map keyed by face tshape index
            if let Some((pc, _r1, _r2)) = ed.pcurves.get(&face_ts_idx) {
                return Ok((pc, face_ts_idx));
            }
            // Fallback: grab first pcurve if any
            if let Some((si, (pc, _, _))) = ed.pcurves.iter().next() {
                return Ok((pc, *si));
            }
            Err(BRepToolsError::MissingGeometry {
                kind: "pcurve",
                index: edge_idx,
            })
        }
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

// =============================================================================
// Helper Functions
// =============================================================================

/// Get a face by its flat index (across all solids/shells).
///
/// Returns a tuple of (face_tshape_index, TShellData reference, TFaceData reference).
fn get_face_by_flat_index<'a>(brep: &'a rcad_kernel::BRep, face_idx: usize) -> Result<(usize, &'a rcad_kernel::topods::TShellData, &'a rcad_kernel::topods::TFaceData), BRepToolsError> {
    let mut current_idx = 0;
    for ts in &brep.tshapes {
        if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let rcad_kernel::topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    if face_idx < current_idx + shd.faces.len() {
                        let local_idx = face_idx - current_idx;
                        let face_sr = &shd.faces[local_idx];
                        if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                            return Ok((face_sr.index, shd, fd));
                        }
                    }
                    current_idx += shd.faces.len();
                }
            }
        }
        if let rcad_kernel::topods::TShape::Shell(shd) = &**ts {
            if face_idx < current_idx + shd.faces.len() {
                let local_idx = face_idx - current_idx;
                let face_sr = &shd.faces[local_idx];
                if let rcad_kernel::topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                    return Ok((face_sr.index, shd, fd));
                }
            }
            current_idx += shd.faces.len();
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
    match &*brep.tshapes[edge_idx] {
        rcad_kernel::topods::TShape::Edge(ed) => Ok(Some(ed.range)),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Check if an edge is degenerate (zero-length, like a pole).
pub fn is_edge_degenerate(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<bool, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        rcad_kernel::topods::TShape::Edge(ed) => Ok(ed.degenerated),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the tolerance of a vertex.
pub fn get_vertex_tolerance(brep: &rcad_kernel::BRep, vertex_idx: usize) -> Result<f64, BRepToolsError> {
    match &*brep.tshapes[vertex_idx] {
        rcad_kernel::topods::TShape::Vertex(vd) => Ok(vd.tolerance),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "vertex",
            index: vertex_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the tolerance of an edge.
pub fn get_edge_tolerance(brep: &rcad_kernel::BRep, edge_idx: usize) -> Result<f64, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        rcad_kernel::topods::TShape::Edge(ed) => Ok(ed.tolerance),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the tolerance of a face.
pub fn get_face_tolerance(brep: &rcad_kernel::BRep, face_idx: usize) -> Result<f64, BRepToolsError> {
    let (_, _, fd) = get_face_by_flat_index(brep, face_idx)?;
    Ok(fd.tolerance)
}

// =============================================================================
// Additional Shape Queries
// =============================================================================

/// Count the total number of faces in a BRep.
pub fn count_faces(brep: &rcad_kernel::BRep) -> usize {
    brep.face_count()
}

/// Count the total number of edges in a BRep.
pub fn count_edges(brep: &rcad_kernel::BRep) -> usize {
    brep.edge_count()
}

/// Count the total number of vertices in a BRep.
pub fn count_vertices(brep: &rcad_kernel::BRep) -> usize {
    brep.vertex_count()
}

/// Count the total number of shells in a BRep.
pub fn count_shells(brep: &rcad_kernel::BRep) -> usize {
    count_shells_topods(brep)
}

/// Count the total number of wires (outer + inner) across all faces in a BRep.
pub fn count_wires(brep: &rcad_kernel::BRep) -> usize {
    count_wires_topods(brep)
}

/// Get the bounding box of a BRep.
///
/// Returns `[min_point, max_point]` or `None` if the BRep has no vertices.
pub fn bounding_box(brep: &rcad_kernel::BRep) -> Option<[DVec3; 2]> {
    let mut vertices: Vec<DVec3> = Vec::new();
    for ts in &brep.tshapes {
        if let rcad_kernel::topods::TShape::Vertex(vd) = &**ts {
            vertices.push(vd.point);
        }
    }
    if vertices.is_empty() {
        return None;
    }
    let mut min_pt = vertices[0];
    let mut max_pt = vertices[0];
    for p in &vertices[1..] {
        min_pt = min_pt.min(*p);
        max_pt = max_pt.max(*p);
    }
    Some([min_pt, max_pt])
}
