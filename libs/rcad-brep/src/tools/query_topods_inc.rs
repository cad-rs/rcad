// =============================================================================
// Topods-native BRep Tools (migration)
// =============================================================================
// These replace the old BRep query functions with topods::BRep equivalents.

/// Determine the shape type of a topods::BRep.
pub fn get_shape_type_topods(brep: &topods::BRep) -> ShapeType {
    for ts in &brep.tshapes {
        match &**ts {
            topods::TShape::Compound(_) => return ShapeType::Compound,
            topods::TShape::CompSolid(_) => return ShapeType::CompSolid,
            topods::TShape::Solid(sd) => {
                let has_faces = sd.shells.iter().any(|sr| {
                    if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
                        !shd.faces.is_empty()
                    } else { false }
                });
                if has_faces { return ShapeType::Solid; }
                let has_shells = sd.shells.iter().any(|sr| {
                    matches!(&*brep.tshapes[sr.index], topods::TShape::Shell(_))
                });
                if has_shells { return ShapeType::Shell; }
                return ShapeType::Solid;
            }
            topods::TShape::Edge(_) => return ShapeType::Edge,
            topods::TShape::Vertex(_) => return ShapeType::Vertex,
            _ => {}
        }
    }
    ShapeType::Empty
}

/// Count the total number of faces in a topods::BRep.
pub fn count_faces_topods(brep: &topods::BRep) -> usize {
    let mut count = 0;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    count += shd.faces.len();
                }
            }
        }
    }
    count
}

/// Count the total number of edges in a topods::BRep.
pub fn count_edges_topods(brep: &topods::BRep) -> usize {
    brep.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::Edge(_))).count()
}

/// Count the total number of vertices in a topods::BRep.
pub fn count_vertices_topods(brep: &topods::BRep) -> usize {
    brep.tshapes.iter().filter(|ts| matches!(&***ts, topods::TShape::Vertex(_))).count()
}

/// Count the total number of shells in a topods::BRep.
pub fn count_shells_topods(brep: &topods::BRep) -> usize {
    let mut count = 0;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            count += sd.shells.len();
        }
        // Also count standalone shells (not in a solid)
        if matches!(&**ts, topods::TShape::Shell(_)) {
            count += 1;
        }
    }
    count
}

/// Count the total number of wires across all faces in a topods::BRep.
pub fn count_wires_topods(brep: &topods::BRep) -> usize {
    let mut count = 0;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                            count += 1; // outer wire
                            count += fd.inner_wires.len(); // inner wires
                        }
                    }
                }
            }
        }
        if let topods::TShape::Shell(shd) = &**ts {
            for face_sr in &shd.faces {
                if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
                    count += 1;
                    count += fd.inner_wires.len();
                }
            }
        }
    }
    count
}

/// Get the bounding box of a topods::BRep.
pub fn bounding_box_topods(brep: &topods::BRep) -> Option<[DVec3; 2]> {
    let mut vertices: Vec<DVec3> = Vec::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Vertex(vd) = &**ts {
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

/// Check if a topods::BRep is a closed solid.
pub fn is_closed_topods(brep: &topods::BRep) -> bool {
    let mut found_any = false;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    if !is_shell_closed_topods(brep, shd) {
                        return false;
                    }
                    found_any = true;
                }
            }
        }
    }
    // Also check standalone shells (not inside a solid)
    for ts in &brep.tshapes {
        if let topods::TShape::Shell(shd) = &**ts {
            if !is_shell_closed_topods(brep, shd) {
                return false;
            }
            found_any = true;
        }
    }
    found_any
}

/// Check if a shell (topods) is closed.
fn is_shell_closed_topods(brep: &topods::BRep, shd: &topods::TShellData) -> bool {
    if shd.faces.is_empty() {
        return false;
    }
    let mut edge_count: HashMap<usize, usize> = HashMap::new();
    for face_sr in &shd.faces {
        if let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] {
            if let topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] {
                for sr in &wd.edges {
                    *edge_count.entry(sr.index).or_insert(0) += 1;
                }
            }
            for inner_sr in &fd.inner_wires {
                if let topods::TShape::Wire(wd) = &*brep.tshapes[inner_sr.index] {
                    for sr in &wd.edges {
                        *edge_count.entry(sr.index).or_insert(0) += 1;
                    }
                }
            }
        }
    }
    edge_count.values().all(|&count| count == 2)
}

/// Get the surface of a face (by face tshape index). Returns None if not a face or no surface.
pub fn get_face_surface_topods(brep: &topods::BRep, face_tshape_idx: usize) -> Result<&Surface3, BRepToolsError> {
    match &*brep.tshapes[face_tshape_idx] {
        topods::TShape::Face(fd) => fd.surface.as_ref().ok_or(BRepToolsError::MissingGeometry {
            kind: "surface",
            index: face_tshape_idx,
        }),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "face",
            index: face_tshape_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Find the first face tshape index matching a given face_index (sequential across solids/shells).
fn face_tshape_idx_by_flat_index(brep: &topods::BRep, face_idx: usize) -> Result<usize, BRepToolsError> {
    let mut current = 0usize;
    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    if face_idx < current + shd.faces.len() {
                        return Ok(shd.faces[face_idx - current].index);
                    }
                    current += shd.faces.len();
                }
            }
        }
        // Standalone shells
        if let topods::TShape::Shell(shd) = &**ts {
            if face_idx < current + shd.faces.len() {
                return Ok(shd.faces[face_idx - current].index);
            }
            current += shd.faces.len();
        }
    }
    Err(BRepToolsError::InvalidIndex { kind: "face", index: face_idx, max: current })
}

/// Get the 3D curve of an edge (by edge tshape index).
pub fn get_edge_curve_topods(brep: &topods::BRep, edge_idx: usize) -> Result<&Curve3, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        topods::TShape::Edge(ed) => ed.curve.as_ref().ok_or(BRepToolsError::MissingGeometry {
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

/// Get the tolerance of a vertex (by vertex tshape index).
pub fn get_vertex_tolerance_topods(brep: &topods::BRep, vertex_idx: usize) -> Result<f64, BRepToolsError> {
    match &*brep.tshapes[vertex_idx] {
        topods::TShape::Vertex(vd) => Ok(vd.tolerance),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "vertex",
            index: vertex_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the tolerance of an edge (by edge tshape index).
pub fn get_edge_tolerance_topods(brep: &topods::BRep, edge_idx: usize) -> Result<f64, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        topods::TShape::Edge(ed) => Ok(ed.tolerance),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the tolerance of a face (by face tshape index).
pub fn get_face_tolerance_topods(brep: &topods::BRep, face_idx: usize) -> Result<f64, BRepToolsError> {
    match &*brep.tshapes[face_idx] {
        topods::TShape::Face(fd) => Ok(fd.tolerance),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "face",
            index: face_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Get the parameter range of an edge's 3D curve (by edge tshape index).
pub fn get_edge_range_topods(brep: &topods::BRep, edge_idx: usize) -> Result<Option<[f64; 2]>, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        topods::TShape::Edge(ed) => Ok(Some(ed.range)),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}

/// Check if an edge is degenerate (by edge tshape index).
pub fn is_edge_degenerate_topods(brep: &topods::BRep, edge_idx: usize) -> Result<bool, BRepToolsError> {
    match &*brep.tshapes[edge_idx] {
        topods::TShape::Edge(ed) => Ok(ed.degenerated),
        _ => Err(BRepToolsError::InvalidIndex {
            kind: "edge",
            index: edge_idx,
            max: brep.tshapes.len(),
        }),
    }
}
