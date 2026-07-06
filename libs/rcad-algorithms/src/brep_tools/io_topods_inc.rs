// =============================================================================
// Topods-native I/O and remaining utilities (migration)
// =============================================================================
use std::io::Read;

/// Serialize a topods::BRep to a JSON string.
pub fn write_brep_to_string_topods(brep: &topods::BRep) -> Result<String, BRepToolsError> {
    serde_json::to_string_pretty(brep)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// Deserialize a topods::BRep from a JSON string.
pub fn read_brep_from_string_topods(s: &str) -> Result<topods::BRep, BRepToolsError> {
    serde_json::from_str(s)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))
}

/// Write a topods::BRep to a file as JSON.
pub fn write_brep_to_file_topods<P: AsRef<Path>>(brep: &topods::BRep, path: P) -> Result<(), BRepToolsError> {
    let json = write_brep_to_string_topods(brep)?;
    let mut file = File::create(path.as_ref())
        .map_err(|e| BRepToolsError::IoError(e.to_string()))?;
    file.write_all(json.as_bytes())
        .map_err(|e| BRepToolsError::IoError(e.to_string()))
}

/// Read a topods::BRep from a JSON file.
pub fn read_brep_from_file_topods<P: AsRef<Path>>(path: P) -> Result<topods::BRep, BRepToolsError> {
    let mut file = File::open(path.as_ref())
        .map_err(|e| BRepToolsError::IoError(e.to_string()))?;
    let mut s = String::new();
    file.read_to_string(&mut s)
        .map_err(|e| BRepToolsError::IoError(e.to_string()))?;
    read_brep_from_string_topods(&s)
}

/// Build an axis-aligned box as topods::BRep given origin and axes.
fn make_box_brep_topods_internal(
    origin: DVec3, x_dir: DVec3, y_dir: DVec3,
    w: f64, h: f64, d: f64,
) -> topods::BRep {
    let mut r = topods::BRep::new();

    let x = x_dir.normalize_or_zero();
    let y = y_dir.normalize_or_zero();
    let z = x.cross(y).normalize_or_zero();

    // 8 corners
    let corners = [
        origin,
        origin + x * w,
        origin + x * w + y * h,
        origin + y * h,
        origin + z * d,
        origin + x * w + z * d,
        origin + x * w + y * h + z * d,
        origin + y * h + z * d,
    ];

    // Add vertices
    let v: Vec<topods::ShapeRef> = corners.iter().map(|&p| r.add_tvertex(p)).collect();

    // 12 edges (6 faces x 4 edges, some shared)
    let face_loops = [
        [0, 1, 2, 3], // bottom
        [4, 5, 6, 7], // top
        [0, 1, 5, 4], // front
        [2, 3, 7, 6], // back
        [0, 3, 7, 4], // left
        [1, 2, 6, 5], // right
    ];

    let mut face_srs = Vec::new();
    for &fv in &face_loops {
        let mut edge_refs = Vec::new();
        for j in 0..4 {
            let a = fv[j];
            let b = fv[(j + 1) % 4];
            let edge_sr = r.add_tedge(None, v[a], v[b], [0.0, 1.0]);
            let orient = if a < b { topods::Orientation::Forward } else { topods::Orientation::Reversed };
            edge_refs.push(topods::ShapeRef::synthetic_with_orientation(edge_sr.index, orient));
        }
        let wire = r.add_twire(edge_refs);
        let face = r.add_tface(None, wire, vec![], None, None, vec![], true);
        face_srs.push(face);
    }

    let shell = r.add_tshell(face_srs);
    r.add_tsolid(vec![shell]);
    r
}

/// Create a half-space solid from a plane, extending in the normal direction.
pub fn make_face_half_space_topods(plane: &rcad_kernel::geom::Plane, bbox: &[DVec3; 2], normal_side: bool) -> topods::BRep {
    let [bmin, bmax] = *bbox;
    let diag = bmax - bmin;
    let margin = diag.length().max(1.0) * 2.0;

    let n = if normal_side { plane.normal } else { -plane.normal }.normalize();
    let abs = n.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z { DVec3::X }
    else if abs.y <= abs.z { DVec3::Y } else { DVec3::Z };
    let u = n.cross(candidate).normalize();
    let v = n.cross(u);

    let origin = plane.origin - u * (margin / 2.0) - v * (margin / 2.0)
        - if normal_side { DVec3::ZERO } else { n * margin };

    make_box_brep_topods_internal(origin, u, v, margin, margin, margin)
}

/// Get the surface of a face by its flat index.
pub fn get_surface_by_flat_index_topods(brep: &topods::BRep, face_idx: usize) -> Result<&Surface3, BRepToolsError> {
    let fi = face_tshape_idx_by_flat_index(brep, face_idx)?;
    get_face_surface_topods(brep, fi)
}

/// Get the 3D curve of an edge by its flat edge index.
pub fn get_curve_by_flat_index_topods(brep: &topods::BRep, edge_idx: usize) -> Result<&Curve3, BRepToolsError> {
    get_edge_curve_topods(brep, edge_idx)
}
