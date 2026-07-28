// =============================================================================
// BRep I/O Utilities
// =============================================================================

/// Serialize a BRep to a JSON string.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::write_brep_to_string;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let json = write_brep_to_string(&brep).unwrap();
/// assert!(json.contains("vertices"));
/// ```
pub fn write_brep_to_string(brep: &rcad_kernel::BRep) -> Result<String, BRepToolsError> {
    serde_json::to_string_pretty(brep)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// Deserialize a BRep from a JSON string.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_tools::{write_brep_to_string, read_brep_from_string};
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let json = write_brep_to_string(&brep).unwrap();
/// let restored = read_brep_from_string(&json).unwrap();
/// assert_eq!(brep.vertices.len(), restored.vertices.len());
/// ```
pub fn read_brep_from_string(s: &str) -> Result<rcad_kernel::BRep, BRepToolsError> {
    serde_json::from_str(s)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))
}

/// Write a BRep to a file as JSON.
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::write_brep_to_file;
/// use rcad_kernel::BRep;
///
/// let brep = BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// write_brep_to_file(&brep, "box.brep").unwrap();
/// ```
pub fn write_brep_to_file<P: AsRef<Path>>(brep: &rcad_kernel::BRep, path: P) -> Result<(), BRepToolsError> {
    let file = File::create(&path)
        .map_err(|e| BRepToolsError::IoError(format!("Failed to create file: {}", e)))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, brep)
        .map_err(|e| BRepToolsError::SerializationError(e.to_string()))
}

/// Read a BRep from a file.
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_tools::read_brep_from_file;
///
/// let brep = read_brep_from_file("box.brep").unwrap();
/// ```
pub fn read_brep_from_file<P: AsRef<Path>>(path: P) -> Result<rcad_kernel::BRep, BRepToolsError> {
    let file = File::open(&path)
        .map_err(|e| BRepToolsError::IoError(format!("Failed to open file: {}", e)))?;
    let reader = BufReader::new(file);
    serde_json::from_reader(reader)
        .map_err(|e| BRepToolsError::DeserializationError(e.to_string()))
}
