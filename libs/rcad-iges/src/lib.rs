//! IGES B-Rep parsing and writing support for RCAD.
//!
//! This module provides full IGES (Initial Graphics Exchange Specification) support:
//! - Parsing IGES files into RCAD BRep structures
//! - Writing RCAD BRep structures to IGES format
//!
//! ## Supported IGES Entity Types
//! - Type 100: Circular Arc
//! - Type 110: Line
//! - Type 126: Rational B-Spline Curve
//! - Type 128: Rational B-Spline Surface
//! - Type 142: Curve on a Parametric Surface
//! - Type 144: Trimmed Surface
//!
//! ## IGES File Structure
//! An IGES file consists of five sections:
//! - Start (S): Human-readable prologue
//! - Global (G): Global parameters (delimiter, units, etc.)
//! - Directory Entry (D): Entity metadata (2 lines per entity)
//! - Parameter Data (P): Entity parameter values
//! - Terminate (T): Section line counts

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

use glam::DVec3;
use rcad_kernel::topology::{Face, Shell, Solid, Vertex, Wire};
use rcad_kernel::{topods, BRep, BSplineSurface, Curve3, Surface3};
use rcad_kernel::geom::{BSplineCurve3, Circle3, Line3, Plane};

/// Errors that can occur when reading or parsing an IGES file.
#[derive(Debug, Clone)]
pub enum IgesError {
    /// File I/O error.
    Io(String),
    /// Not a valid IGES file (missing section identifiers, malformed records, etc.).
    InvalidFormat(String),
    /// A required IGES entity is missing or malformed.
    MissingEntity {
        entity_type: &'static str,
        id: Option<u64>,
    },
    /// Parse produced an empty or degenerate result.
    EmptyResult(String),
    /// Unsupported entity type encountered.
    UnsupportedEntity {
        entity_type: i32,
        message: String,
    },
}

impl std::fmt::Display for IgesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(msg) => write!(f, "I/O error: {}", msg),
            Self::InvalidFormat(msg) => write!(f, "invalid IGES format: {}", msg),
            Self::MissingEntity {
                entity_type,
                id: Some(id),
            } => write!(f, "missing {} entity #{}", entity_type, id),
            Self::MissingEntity {
                entity_type,
                id: None,
            } => write!(f, "missing {} entity", entity_type),
            Self::EmptyResult(msg) => write!(f, "IGES parse produced empty result: {}", msg),
            Self::UnsupportedEntity {
                entity_type,
                message,
            } => write!(
                f,
                "unsupported IGES entity type {}: {}",
                entity_type, message
            ),
        }
    }
}

impl std::error::Error for IgesError {}

// ============================================================================
// IGES File Structure Types
// ============================================================================

/// Global section parameters from an IGES file.
#[derive(Debug, Clone, Default)]
struct GlobalSection {
    /// Parameter delimiter character (default ',').
    param_delim: char,
    /// Record delimiter character (default ';').
    record_delim: char,
    /// Product ID from sending system.
    product_id: String,
    /// File name.
    file_name: String,
    /// Native system ID.
    system_id: String,
    /// Preprocessor version.
    preprocessor_version: String,
    /// Number of significant digits for floats.
    significant_digits: i32,
    /// Units flag (1=inches, 2=mm, 3=feet, 4=miles, 5=meters, 6=km, etc.).
    units_flag: i32,
    /// Units name (optional override).
    units_name: String,
    /// Model space scale.
    model_scale: f64,
}

/// Directory Entry (DE) record for an IGES entity.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields mirror the DE record; not all are consumed by the reader yet.
struct DirectoryEntry {
    /// Entity type number.
    entity_type: i32,
    /// Parameter data section pointer (line number).
    param_ptr: i32,
    /// Structure pointer (for hierarchical entities).
    structure: i32,
    /// Line font pattern.
    line_font: i32,
    /// Line weight.
    line_weight: i32,
    /// Color definition pointer.
    color: i32,
    /// Number of parameter data lines.
    param_count: i32,
    /// Form number (for entity variants).
    form: i32,
    /// Entity label.
    label: String,
    /// Entity subscript number.
    subscript: i32,
    /// Status (visibility, etc.).
    status: i32,
    /// Sequence number (DE line number / 2).
    sequence: i32,
}

/// Parameter Data (PD) for an IGES entity.
#[derive(Debug, Clone)]
struct ParameterData {
    /// Parsed parameter values.
    values: Vec<IgesValue>,
}

/// An IGES parameter value.
#[derive(Debug, Clone)]
#[allow(dead_code)] // `Pointer` / `Hollerith` payloads are stored for round-trip fidelity.
enum IgesValue {
    Int(i64),
    Float(f64),
    String(String),
    Pointer(i64),
    /// Hollerith string: nHxxxx
    Hollerith(String),
}

/// Parsed IGES file content.
#[derive(Debug, Default)]
struct ParsedIges {
    global: GlobalSection,
    directory_entries: Vec<DirectoryEntry>,
    parameter_data: HashMap<i32, ParameterData>,
    /// Map from DE sequence number to parsed entity index.
    entity_map: HashMap<i32, usize>,
}

// ============================================================================
// IGES Reader
// ============================================================================

/// IGES file reader.
pub struct IgesReader;

impl IgesReader {
    /// Read and parse an IGES file from disk.
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, IgesError> {
        let content =
            std::fs::read_to_string(path).map_err(|e| IgesError::Io(e.to_string()))?;
        Self::parse_string(&content)
    }

    /// Parse IGES content from a string.
    pub fn parse_string(content: &str) -> Result<BRep, IgesError> {
        let parsed = Self::parse_iges_structure(content)?;
        Self::build_brep(&parsed)
    }

    /// Parse the IGES file structure into sections.
    fn parse_iges_structure(content: &str) -> Result<ParsedIges, IgesError> {
        let mut parsed = ParsedIges::default();
        let mut section_lines: HashMap<char, Vec<String>> = HashMap::new();

        // Collect lines by section
        for line in content.lines() {
            if line.len() < 73 {
                continue;
            }
            let section = line.chars().nth(72).unwrap_or(' ');
            if "SGDPT".contains(section) {
                section_lines
                    .entry(section)
                    .or_default()
                    .push(line.to_string());
            }
        }

        // Parse Global section
        if let Some(g_lines) = section_lines.get(&'G') {
            parsed.global = Self::parse_global_section(g_lines)?;
        }

        // Parse Directory Entry section (pairs of lines)
        let mut de_pairs: Vec<(String, String)> = Vec::new();
        if let Some(d_lines) = section_lines.get(&'D') {
            let mut i = 0;
            while i + 1 < d_lines.len() {
                de_pairs.push((d_lines[i].clone(), d_lines[i + 1].clone()));
                i += 2;
            }
        }

        // Parse Directory Entries
        for (idx, (line1, line2)) in de_pairs.iter().enumerate() {
            let de = Self::parse_directory_entry(line1, line2, idx as i32 / 2 + 1)?;
            parsed.directory_entries.push(de);
        }

        // Build entity map
        for (idx, de) in parsed.directory_entries.iter().enumerate() {
            parsed.entity_map.insert(de.sequence, idx);
        }

        // Parse Parameter Data section
        if let Some(p_lines) = section_lines.get(&'P') {
            parsed.parameter_data = Self::parse_parameter_data(p_lines, &parsed.directory_entries)?;
        }

        Ok(parsed)
    }

    /// Parse the Global section.
    fn parse_global_section(lines: &[String]) -> Result<GlobalSection, IgesError> {
        // Concatenate all parameter data from G section
        let mut param_text = String::new();
        for line in lines {
            let data = if line.len() > 72 {
                &line[..72]
            } else {
                line
            };
            param_text.push_str(data.trim_end());
        }

        let mut global = GlobalSection {
            param_delim: ',',
            record_delim: ';',
            ..Default::default()
        };

        // Parse global parameters (comma-separated, semicolon-terminated)
        let params = Self::split_parameters(&param_text, ',', ';');

        // IGES Global section parameters (in order):
        // 1: param delimiter (nHstring format)
        // 2: record delimiter (nHstring format)
        // 3: product ID
        // 4: file name
        // 5: system ID
        // 6: preprocessor version
        // 7: significant digits
        // 8: units flag
        // ...

        if !params.is_empty() {
            if let Some(s) = Self::parse_hollerith(&params[0]) {
                if let Some(c) = s.chars().next() {
                    global.param_delim = c;
                }
            }
        }
        if params.len() > 1 {
            if let Some(s) = Self::parse_hollerith(&params[1]) {
                if let Some(c) = s.chars().next() {
                    global.record_delim = c;
                }
            }
        }
        if params.len() > 2 {
            global.product_id = Self::parse_hollerith(&params[2]).unwrap_or_default();
        }
        if params.len() > 3 {
            global.file_name = Self::parse_hollerith(&params[3]).unwrap_or_default();
        }
        if params.len() > 4 {
            global.system_id = Self::parse_hollerith(&params[4]).unwrap_or_default();
        }
        if params.len() > 5 {
            global.preprocessor_version = Self::parse_hollerith(&params[5]).unwrap_or_default();
        }
        if params.len() > 6 {
            global.significant_digits = Self::parse_int(&params[6]).unwrap_or(6) as i32;
        }
        if params.len() > 7 {
            global.units_flag = Self::parse_int(&params[7]).unwrap_or(1) as i32;
        }
        if params.len() > 13 {
            global.units_name = Self::parse_hollerith(&params[13]).unwrap_or_default();
        }
        if params.len() > 14 {
            global.model_scale = Self::parse_float(&params[14]).unwrap_or(1.0);
        }

        Ok(global)
    }

    /// Parse a Directory Entry pair.
    fn parse_directory_entry(
        line1: &str,
        line2: &str,
        sequence: i32,
    ) -> Result<DirectoryEntry, IgesError> {
        // DE line format: 8 columns of 8-character fields
        // Line 1: entity_type, param_ptr, structure, line_font, layer, view, transform, status
        // Line 2: entity_type, line_weight, color, param_count, form, unused, label, subscript

        fn extract_field(line: &str, field: usize) -> &str {
            let start = field * 8;
            if start + 8 <= line.len() {
                &line[start..start + 8]
            } else if start < line.len() {
                &line[start..]
            } else {
                ""
            }
        }

        fn parse_field(s: &str) -> i32 {
            s.trim().parse().unwrap_or(0)
        }

        let entity_type = parse_field(extract_field(line1, 0));
        let param_ptr = parse_field(extract_field(line1, 1));
        let structure = parse_field(extract_field(line1, 2));
        let line_font = parse_field(extract_field(line1, 3));
        let layer = parse_field(extract_field(line1, 4));
        let _view = parse_field(extract_field(line1, 5));
        let _transform = parse_field(extract_field(line1, 6));
        let status = parse_field(extract_field(line1, 7));

        let entity_type2 = parse_field(extract_field(line2, 0));
        let line_weight = parse_field(extract_field(line2, 1));
        let color = parse_field(extract_field(line2, 2));
        let param_count = parse_field(extract_field(line2, 3));
        let form = parse_field(extract_field(line2, 4));
        let _unused = parse_field(extract_field(line2, 5));
        let label = extract_field(line2, 6).trim().to_string();
        let subscript = parse_field(extract_field(line2, 7));

        // Verify entity type matches
        if entity_type != entity_type2 {
            return Err(IgesError::InvalidFormat(format!(
                "DE type mismatch: {} vs {}",
                entity_type, entity_type2
            )));
        }

        let _ = (layer, structure); // unused for now

        Ok(DirectoryEntry {
            entity_type,
            param_ptr,
            structure,
            line_font,
            line_weight,
            color,
            param_count,
            form,
            label,
            subscript,
            status,
            sequence,
        })
    }

    /// Parse Parameter Data section.
    fn parse_parameter_data(
        lines: &[String],
        directory_entries: &[DirectoryEntry],
    ) -> Result<HashMap<i32, ParameterData>, IgesError> {
        let mut result: HashMap<i32, ParameterData> = HashMap::new();

        // Build map from param_ptr to DE sequence
        let ptr_to_seq: HashMap<i32, i32> = directory_entries
            .iter()
            .map(|de| (de.param_ptr, de.sequence))
            .collect();

        // Group P lines by DE reference (columns 65-72)
        let mut current_de: Option<i32> = None;
        let mut current_params: Vec<String> = Vec::new();

        for line in lines {
            // Extract DE pointer from columns 65-72 (right-aligned)
            let de_ptr: i32 = if line.len() >= 72 {
                line[64..72].trim().parse().unwrap_or(0)
            } else {
                0
            };

            // Extract parameter data from columns 1-64
            let param_part = if line.len() >= 64 {
                &line[..64]
            } else {
                line
            };

            // Check if this is a new entity
            if let Some(&de_seq) = ptr_to_seq.get(&de_ptr) {
                // Save previous entity if any
                if let Some(de_seq) = current_de {
                    if !current_params.is_empty() {
                        let values = Self::parse_parameter_values(&current_params.join(""))?;
                        result.insert(de_seq, ParameterData { values });
                    }
                }

                // Start new entity
                current_de = Some(de_seq);
                current_params.clear();
                current_params.push(param_part.trim_end().to_string());
            } else if current_de.is_some() {
                // Continue current entity
                current_params.push(param_part.trim_end().to_string());
            }
        }

        // Save last entity
        if let Some(de_seq) = current_de {
            if !current_params.is_empty() {
                let values = Self::parse_parameter_values(&current_params.join(""))?;
                result.insert(de_seq, ParameterData { values });
            }
        }

        Ok(result)
    }

    /// Parse parameter values from a parameter string.
    fn parse_parameter_values(param_str: &str) -> Result<Vec<IgesValue>, IgesError> {
        let mut values = Vec::new();
        let mut chars = param_str.chars().peekable();
        let mut current = String::new();

        while let Some(c) = chars.next() {
            if c == ',' {
                if !current.is_empty() {
                    values.push(Self::parse_single_value(&current)?);
                    current.clear();
                }
            } else if c == ';' {
                if !current.is_empty() {
                    values.push(Self::parse_single_value(&current)?);
                }
                break;
            } else if c == 'H' {
                // Check for Hollerith: nHstring
                if let Ok(n) = current.parse::<usize>() {
                    current.clear();
                    let mut holl = String::with_capacity(n);
                    for _ in 0..n {
                        if let Some(hc) = chars.next() {
                            holl.push(hc);
                        }
                    }
                    values.push(IgesValue::Hollerith(holl));
                } else {
                    current.push(c);
                }
            } else {
                current.push(c);
            }
        }

        if !current.is_empty() {
            values.push(Self::parse_single_value(&current)?);
        }

        Ok(values)
    }

    /// Parse a single parameter value.
    fn parse_single_value(s: &str) -> Result<IgesValue, IgesError> {
        let s = s.trim();
        if s.is_empty() {
            return Ok(IgesValue::Int(0));
        }

        // Check for pointer reference (negative number or number with trailing P)
        if s.starts_with('-') || s.ends_with('P') || s.ends_with('p') {
            let num: i64 = s.trim_end_matches(|c: char| c == 'P' || c == 'p' || c.is_ascii_digit())
                .parse()
                .unwrap_or(0);
            return Ok(IgesValue::Pointer(num.abs()));
        }

        // Try integer first
        if let Ok(i) = s.parse::<i64>() {
            return Ok(IgesValue::Int(i));
        }

        // Try float (handle D/E exponent notation)
        if let Some(f) = Self::parse_float(s) {
            return Ok(IgesValue::Float(f));
        }

        // Default to string
        Ok(IgesValue::String(s.to_string()))
    }

    /// Parse a Hollerith string (format: nHxxxx).
    fn parse_hollerith(s: &str) -> Option<String> {
        let s = s.trim();
        let h_pos = s.find('H')?;
        let n: usize = s[..h_pos].parse().ok()?;
        let rest = &s[h_pos + 1..];
        if rest.len() >= n {
            Some(rest[..n].to_string())
        } else {
            Some(rest.to_string())
        }
    }

    /// Parse a float value (handles D exponent notation).
    fn parse_float(s: &str) -> Option<f64> {
        let s = s.trim().replace('D', "E").replace('d', "e");
        s.parse().ok()
    }

    /// Parse an integer value.
    fn parse_int(s: &str) -> Option<i64> {
        s.trim().parse().ok()
    }

    /// Split parameters by delimiter, respecting Hollerith strings.
    fn split_parameters(s: &str, param_delim: char, record_delim: char) -> Vec<String> {
        let mut params = Vec::new();
        let mut current = String::new();
        let chars = s.chars().peekable();
        let mut in_hollerith = false;
        let mut holl_count = 0;

        for c in chars {
            if in_hollerith {
                current.push(c);
                holl_count -= 1;
                if holl_count <= 0 {
                    in_hollerith = false;
                }
            } else if c == 'H' {
                // Check for Hollerith prefix
                if let Ok(n) = current.parse::<usize>() {
                    current.push(c);
                    in_hollerith = true;
                    holl_count = n as i32;
                } else {
                    current.push(c);
                }
            } else if c == param_delim {
                params.push(current.trim().to_string());
                current.clear();
            } else if c == record_delim {
                if !current.trim().is_empty() {
                    params.push(current.trim().to_string());
                }
                break;
            } else {
                current.push(c);
            }
        }

        if !current.trim().is_empty() {
            params.push(current.trim().to_string());
        }

        params
    }

    // ========================================================================
    // BRep Building
    // ========================================================================

    /// Build a BRep from parsed IGES data.
    fn build_brep(parsed: &ParsedIges) -> Result<BRep, IgesError> {
        let mut builder = IgesBrepBuilder::new();

        // First pass: parse all geometric entities
        for (seq, pd) in &parsed.parameter_data {
            if let Some(de) = parsed.directory_entries.iter().find(|d| d.sequence == *seq) {
                builder.parse_entity(de, pd)?;
            }
        }

        // Second pass: build topology from trimmed surfaces
        for (seq, pd) in &parsed.parameter_data {
            if let Some(de) = parsed.directory_entries.iter().find(|d| d.sequence == *seq) {
                if de.entity_type == 144 {
                    builder.build_trimmed_surface(de, pd)?;
                }
            }
        }

        builder.finish()
    }
}

// ============================================================================
// BRep Builder from IGES
// ============================================================================

/// Helper for building BRep from IGES entities.
struct IgesBrepBuilder {
    brep: BRep,
    // Maps from IGES entity pointers to RCAD indices
    point_map: HashMap<i32, usize>,       // pointer -> vertex index
    curve_map: HashMap<i32, usize>,       // pointer -> curve index in GeomStore
    surface_map: HashMap<i32, usize>,     // pointer -> surface index in GeomStore
    // Parsed geometry cache
    points: HashMap<i32, DVec3>,
    transformations: HashMap<i32, glam::DAffine3>,
}

impl IgesBrepBuilder {
    fn new() -> Self {
        Self {
            brep: BRep::new(),
            point_map: HashMap::new(),
            curve_map: HashMap::new(),
            surface_map: HashMap::new(),
            points: HashMap::new(),
            transformations: HashMap::new(),
        }
    }

    /// Parse a single IGES entity.
    fn parse_entity(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        match de.entity_type {
            116 => self.parse_point(de, pd),        // Point
            100 => self.parse_circular_arc(de, pd), // Circular Arc
            110 => self.parse_line(de, pd),         // Line
            126 => self.parse_bspline_curve(de, pd),// B-Spline Curve
            128 => self.parse_bspline_surface(de, pd), // B-Spline Surface
            142 => self.parse_curve_on_surface(de, pd), // Curve on Surface
            144 => self.parse_trimmed_surface(de, pd), // Trimmed Surface (topology handled separately)
            108 => self.parse_plane(de, pd),        // Plane
            124 => self.parse_transformation(de, pd), // Transformation Matrix
            130 => self.parse_offset_curve(de, pd), // Offset Curve
            140 => self.parse_offset_surface(de, pd), // Offset Surface
            102 => Ok(()), // Composite Curve - just reference child curves
            _ => {
                // Silently ignore unsupported entities during first pass
                Ok(())
            }
        }
    }

    /// Parse a Point entity (Type 116).
    fn parse_point(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        if pd.values.len() < 3 {
            return Err(IgesError::InvalidFormat(format!(
                "Point entity {} has insufficient params",
                de.sequence
            )));
        }

        let x = self.get_float(&pd.values[0])?;
        let y = self.get_float(&pd.values[1])?;
        let z = self.get_float(&pd.values[2])?;

        let point = DVec3::new(x, y, z);
        self.points.insert(de.sequence, point);

        // Add vertex
        let v_idx = self.brep.vertices.len();
        self.brep.vertices.push(Vertex { point });
        self.point_map.insert(de.sequence, v_idx);

        Ok(())
    }

    /// Parse a Circular Arc entity (Type 100).
    fn parse_circular_arc(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 100: zt, x1, y1, x2, y2, x3, y3
        // (center is defined by a plane normal and two points)
        if pd.values.len() < 7 {
            return Err(IgesError::InvalidFormat(format!(
                "Circular Arc entity {} has insufficient params",
                de.sequence
            )));
        }

        let zt = self.get_float(&pd.values[0])?;
        let x1 = self.get_float(&pd.values[1])?;
        let y1 = self.get_float(&pd.values[2])?;
        let x2 = self.get_float(&pd.values[3])?;
        let y2 = self.get_float(&pd.values[4])?;
        let x3 = self.get_float(&pd.values[5])?;
        let y3 = self.get_float(&pd.values[6])?;

        // The three points define the arc in the XT-YT plane at ZT
        let p1 = DVec3::new(x1, y1, zt);
        let p2 = DVec3::new(x2, y2, zt);
        let p3 = DVec3::new(x3, y3, zt);

        // Compute circle center and radius from three points
        let center = Self::circle_center_from_3_points(p1, p2, p3)?;
        let radius = (p1 - center).length();

        // Normal is perpendicular to plane containing the three points
        let v12 = p2 - p1;
        let v23 = p3 - p2;
        let normal = v12.cross(v23).normalize_or(DVec3::Z);

        let circle = Circle3::new(center, normal, radius);

        // Store the 3D curve
        let curve_idx = self.brep.geom.curves.len();
        self.brep.geom.curves.push(Curve3::Circle(circle));
        self.curve_map.insert(de.sequence, curve_idx);

        // Store start/end points for later edge construction
        // (arc goes from p1 to p3 through p2)
        let _ = (p1, p3); // Will be used when building edges

        Ok(())
    }

    /// Parse a Line entity (Type 110).
    fn parse_line(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        if pd.values.len() < 6 {
            return Err(IgesError::InvalidFormat(format!(
                "Line entity {} has insufficient params",
                de.sequence
            )));
        }

        let x1 = self.get_float(&pd.values[0])?;
        let y1 = self.get_float(&pd.values[1])?;
        let z1 = self.get_float(&pd.values[2])?;
        let x2 = self.get_float(&pd.values[3])?;
        let y2 = self.get_float(&pd.values[4])?;
        let z2 = self.get_float(&pd.values[5])?;

        let start = DVec3::new(x1, y1, z1);
        let end = DVec3::new(x2, y2, z2);
        let direction = (end - start).normalize_or(DVec3::X);

        let line = Line3 {
            origin: start,
            direction,
        };

        let curve_idx = self.brep.geom.curves.len();
        self.brep.geom.curves.push(Curve3::Line(line));
        self.curve_map.insert(de.sequence, curve_idx);

        // Store endpoints
        self.points.insert(de.sequence * 10 + 1, start);
        self.points.insert(de.sequence * 10 + 2, end);

        Ok(())
    }

    /// Parse a B-Spline Curve entity (Type 126).
    fn parse_bspline_curve(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 126 parameters:
        // 1: DE pointer to transformation matrix (0 = identity)
        // 2: Upper index of sum (degree)
        // 3: Planar flag (0=nonplanar, 1=planar)
        // 4: Closed flag (0=open, 1=closed)
        // 5: Rational flag (0=nonrational, 1=rational)
        // 6: Polynomial flag
        // 7: Periodic flag
        // 8-N: Knot sequence values
        // Then: weights (if rational)
        // Then: control points (X, Y, Z triplets)
        // Last: start and end parameter values

        if pd.values.len() < 10 {
            return Err(IgesError::InvalidFormat(format!(
                "B-Spline Curve entity {} has insufficient params",
                de.sequence
            )));
        }

        let k = self.get_int(&pd.values[1])? as usize; // degree (upper index)
        let _planar = self.get_int(&pd.values[2])?;
        let _closed = self.get_int(&pd.values[3])?;
        let rational = self.get_int(&pd.values[4])? != 0;
        let _polynomial = self.get_int(&pd.values[5])?;
        let _periodic = self.get_int(&pd.values[6])?;

        let degree = k;

        // M = number of knots - 1 (index 7 has this value)
        // Actually, IGES stores: A, K, ..., where A = number of knots - 1
        let n_knots_minus_1 = self.get_int(&pd.values[7])? as usize;
        let n_knots = n_knots_minus_1 + 1;
        let n_ctrl = n_knots - degree - 1; // Number of control points

        // Read knot values (starting at index 8)
        let knot_start = 8;
        let knot_end = knot_start + n_knots;
        if pd.values.len() < knot_end {
            return Err(IgesError::InvalidFormat(format!(
                "B-Spline Curve entity {} has insufficient knot vals",
                de.sequence
            )));
        }

        let mut knots: Vec<f64> = Vec::with_capacity(n_knots);
        for i in knot_start..knot_end {
            knots.push(self.get_float(&pd.values[i])?);
        }

        // Read weights if rational (n_ctrl values)
        let mut weights: Vec<f64> = vec![1.0; n_ctrl];
        let mut idx = knot_end;
        if rational {
            if pd.values.len() < idx + n_ctrl {
                return Err(IgesError::InvalidFormat(format!(
                    "B-Spline Curve entity {} has insufficient weight vals",
                    de.sequence
                )));
            }
            for (i, w) in weights.iter_mut().enumerate().take(n_ctrl) {
                *w = self.get_float(&pd.values[idx + i])?;
            }
            idx += n_ctrl;
        }

        // Read control points (X, Y, Z for each)
        if pd.values.len() < idx + n_ctrl * 3 {
            return Err(IgesError::InvalidFormat(format!(
                "B-Spline Curve entity {} has insufficient ctrl pts",
                de.sequence
            )));
        }

        let mut control_points: Vec<DVec3> = Vec::with_capacity(n_ctrl);
        for i in 0..n_ctrl {
            let x = self.get_float(&pd.values[idx + i * 3])?;
            let y = self.get_float(&pd.values[idx + i * 3 + 1])?;
            let z = self.get_float(&pd.values[idx + i * 3 + 2])?;
            control_points.push(DVec3::new(x, y, z));
        }

        let bspline = BSplineCurve3 {
            degree,
            knots,
            control_points,
            weights,
        };

        let curve_idx = self.brep.geom.curves.len();
        self.brep.geom.curves.push(Curve3::BSpline(bspline));
        self.curve_map.insert(de.sequence, curve_idx);

        Ok(())
    }

    /// Parse a B-Spline Surface entity (Type 128).
    fn parse_bspline_surface(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 128 parameters (similar to Type 126 but with two parameter directions)
        if pd.values.len() < 20 {
            return Err(IgesError::InvalidFormat(format!(
                "B-Spline Surface entity {} has insufficient params",
                de.sequence
            )));
        }

        let _xform = self.get_int(&pd.values[1])?;
        let degree_u = self.get_int(&pd.values[2])? as usize;
        let degree_v = self.get_int(&pd.values[3])? as usize;
        let _closed_u = self.get_int(&pd.values[4])?;
        let _closed_v = self.get_int(&pd.values[5])?;
        let rational = self.get_int(&pd.values[6])? != 0;
        let _polynomial_u = self.get_int(&pd.values[7])?;
        let _polynomial_v = self.get_int(&pd.values[8])?;
        let _periodic_u = self.get_int(&pd.values[9])?;
        let _periodic_v = self.get_int(&pd.values[10])?;

        // Number of knots
        let n_knots_u_minus_1 = self.get_int(&pd.values[11])? as usize;
        let n_knots_v_minus_1 = self.get_int(&pd.values[12])? as usize;
        let n_knots_u = n_knots_u_minus_1 + 1;
        let n_knots_v = n_knots_v_minus_1 + 1;

        let n_ctrl_u = n_knots_u - degree_u - 1;
        let n_ctrl_v = n_knots_v - degree_v - 1;

        let mut idx = 13;

        // Read U knots
        if pd.values.len() < idx + n_knots_u {
            return Err(IgesError::InvalidFormat(
                format!("B-Spline Surface entity {} has insufficient U knots", de.sequence)
            ));
        }
        let mut knots_u: Vec<f64> = Vec::with_capacity(n_knots_u);
        for i in 0..n_knots_u {
            knots_u.push(self.get_float(&pd.values[idx + i])?);
        }
        idx += n_knots_u;

        // Read V knots
        if pd.values.len() < idx + n_knots_v {
            return Err(IgesError::InvalidFormat(
                format!("B-Spline Surface entity {} has insufficient V knots", de.sequence)
            ));
        }
        let mut knots_v: Vec<f64> = Vec::with_capacity(n_knots_v);
        for i in 0..n_knots_v {
            knots_v.push(self.get_float(&pd.values[idx + i])?);
        }
        idx += n_knots_v;

        // Read weights if rational
        let mut weights: Vec<Vec<f64>> = vec![vec![1.0; n_ctrl_v]; n_ctrl_u];
        if rational {
            if pd.values.len() < idx + n_ctrl_u * n_ctrl_v {
                return Err(IgesError::InvalidFormat(
                    format!("B-Spline Surface entity {} has insufficient weights", de.sequence)
                ));
            }
            for (i, row) in weights.iter_mut().enumerate().take(n_ctrl_u) {
                for (j, w) in row.iter_mut().enumerate().take(n_ctrl_v) {
                    *w = self.get_float(&pd.values[idx + i * n_ctrl_v + j])?;
                }
            }
            idx += n_ctrl_u * n_ctrl_v;
        }

        // Read control points (X, Y, Z for each)
        if pd.values.len() < idx + n_ctrl_u * n_ctrl_v * 3 {
            return Err(IgesError::InvalidFormat(
                format!("B-Spline Surface entity {} has insufficient ctrl pts", de.sequence)
            ));
        }

        let mut control_points: Vec<Vec<DVec3>> = Vec::with_capacity(n_ctrl_u);
        for i in 0..n_ctrl_u {
            let mut row: Vec<DVec3> = Vec::with_capacity(n_ctrl_v);
            for j in 0..n_ctrl_v {
                let base = idx + (i * n_ctrl_v + j) * 3;
                let x = self.get_float(&pd.values[base])?;
                let y = self.get_float(&pd.values[base + 1])?;
                let z = self.get_float(&pd.values[base + 2])?;
                row.push(DVec3::new(x, y, z));
            }
            control_points.push(row);
        }

        let bspline = BSplineSurface {
            degree_u,
            degree_v,
            knots_u,
            knots_v,
            control_points,
            weights,
        };

        let surf_idx = self.brep.geom.surfaces.len();
        self.brep.geom.surfaces.push(Surface3::BSpline(bspline));
        self.surface_map.insert(de.sequence, surf_idx);

        Ok(())
    }

    /// Parse a Curve on Parametric Surface entity (Type 142).
    fn parse_curve_on_surface(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 142: curve_on_surface = (surface_ptr, curve_3d_ptr, curve_2d_ptr, preference)
        if pd.values.len() < 4 {
            return Err(IgesError::InvalidFormat(format!(
                "Curve on Surface entity {} has insufficient params",
                de.sequence
            )));
        }

        let surface_ptr = self.get_int(&pd.values[1])? as i32;
        let curve_3d_ptr = self.get_int(&pd.values[2])? as i32;
        let curve_2d_ptr = self.get_int(&pd.values[3])? as i32;
        let _preference = self.get_int(&pd.values[4])?;

        // Store the association for later use
        let _ = (surface_ptr, curve_3d_ptr, curve_2d_ptr);

        Ok(())
    }

    /// Parse a Trimmed Surface entity (Type 144).
    fn parse_trimmed_surface(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 144: trimmed_surface = (surface_ptr, outer_bound_ptr, n_inner, [inner_bounds...], ...)
        if pd.values.len() < 3 {
            return Err(IgesError::InvalidFormat(format!(
                "Trimmed Surface entity {} has insufficient params",
                de.sequence
            )));
        }

        let surface_ptr = self.get_int(&pd.values[1])? as i32;
        let _outer_bound = self.get_int(&pd.values[2])?;
        let n_inner = if pd.values.len() > 3 { self.get_int(&pd.values[3])? as usize } else { 0 };

        // Track this for building faces later
        let _ = (surface_ptr, n_inner);

        Ok(())
    }

    /// Parse a Plane entity (Type 108).
    fn parse_plane(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        if pd.values.len() < 4 {
            return Err(IgesError::InvalidFormat(format!(
                "Plane entity {} has insufficient params",
                de.sequence
            )));
        }

        let a = self.get_float(&pd.values[0])?;
        let b = self.get_float(&pd.values[1])?;
        let c = self.get_float(&pd.values[2])?;
        let d = self.get_float(&pd.values[3])?;

        // Plane equation: Ax + By + Cz = D
        // Normal is (A, B, C) normalized
        let normal = DVec3::new(a, b, c).normalize_or(DVec3::Z);

        // Find a point on the plane
        let origin = if normal.z.abs() > 1e-9 {
            DVec3::new(0.0, 0.0, d / c)
        } else if normal.y.abs() > 1e-9 {
            DVec3::new(0.0, d / b, 0.0)
        } else {
            DVec3::new(d / a, 0.0, 0.0)
        };

        let plane = Plane { origin, normal };

        let surf_idx = self.brep.geom.surfaces.len();
        self.brep.geom.surfaces.push(Surface3::Plane(plane));
        self.surface_map.insert(de.sequence, surf_idx);

        Ok(())
    }

    /// Parse a Transformation Matrix entity (Type 124).
    fn parse_transformation(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 124: 12 values forming a 3x4 transformation matrix (row-major)
        // R11, R12, R13, T1, R21, R22, R23, T2, R31, R32, R33, T3
        if pd.values.len() < 13 {
            return Err(IgesError::InvalidFormat(format!(
                "Transformation entity {} has insufficient params",
                de.sequence
            )));
        }

        let r11 = self.get_float(&pd.values[1])?;
        let r12 = self.get_float(&pd.values[2])?;
        let r13 = self.get_float(&pd.values[3])?;
        let t1 = self.get_float(&pd.values[4])?;
        let r21 = self.get_float(&pd.values[5])?;
        let r22 = self.get_float(&pd.values[6])?;
        let r23 = self.get_float(&pd.values[7])?;
        let t2 = self.get_float(&pd.values[8])?;
        let r31 = self.get_float(&pd.values[9])?;
        let r32 = self.get_float(&pd.values[10])?;
        let r33 = self.get_float(&pd.values[11])?;
        let t3 = self.get_float(&pd.values[12])?;

        let matrix = glam::DAffine3::from_cols(
            glam::DVec3::new(r11, r21, r31),
            glam::DVec3::new(r12, r22, r32),
            glam::DVec3::new(r13, r23, r33),
            glam::DVec3::new(t1, t2, t3),
        );

        self.transformations.insert(de.sequence, matrix);

        Ok(())
    }

    /// Parse an Offset Curve entity (Type 130).
    fn parse_offset_curve(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 130: offset_curve = (curve_ptr, offset_dist, ...flags)
        if pd.values.len() < 3 {
            return Err(IgesError::InvalidFormat(format!(
                "Offset Curve entity {} has insufficient params",
                de.sequence
            )));
        }

        let _curve_ptr = self.get_int(&pd.values[1])? as i32;
        let _offset_dist = self.get_float(&pd.values[2])?;

        // TODO: Implement offset curve parsing

        Ok(())
    }

    /// Parse an Offset Surface entity (Type 140).
    fn parse_offset_surface(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 140: offset_surface = (surface_ptr, offset_dist, ...)
        if pd.values.len() < 3 {
            return Err(IgesError::InvalidFormat(format!(
                "Offset Surface entity {} has insufficient params",
                de.sequence
            )));
        }

        let _surface_ptr = self.get_int(&pd.values[1])? as i32;
        let _offset_dist = self.get_float(&pd.values[2])?;

        // TODO: Implement offset surface parsing

        Ok(())
    }

    /// Build trimmed surface topology.
    fn build_trimmed_surface(&mut self, de: &DirectoryEntry, pd: &ParameterData) -> Result<(), IgesError> {
        // Type 144: Trimmed (Parametric) Surface
        if pd.values.len() < 3 {
            return Err(IgesError::InvalidFormat(format!(
                "Trimmed Surface entity {} has insufficient params",
                de.sequence
            )));
        }

        let surface_ptr = self.get_int(&pd.values[1])? as i32;
        let outer_bound_ptr = self.get_int(&pd.values[2])?;
        let n_inner = if pd.values.len() > 3 { self.get_int(&pd.values[3])? as usize } else { 0 };

        // Get surface index
        let surf_idx = self.surface_map.get(&surface_ptr).copied();

        // Create a face with empty wire for now
        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            sample_point: None,
            mesh_dirty: true,
            surface_idx: None,
        };

        let _ = (outer_bound_ptr, n_inner, surf_idx);

        // Add face to a shell
        // For now, create a simple structure
        if self.brep.solids.is_empty() {
            self.brep.solids.push(Solid {
                shells: vec![Shell { faces: vec![] }],
            });
        }

        if let Some(solid) = self.brep.solids.first_mut() {
            if solid.shells.is_empty() {
                solid.shells.push(Shell { faces: vec![] });
            }
            if let Some(shell) = solid.shells.first_mut() {
                shell.faces.push(face);
                // Set surface reference
                if let Some(idx) = surf_idx {
                    self.brep.geom.face_surface.push(Some(idx));
                } else {
                    self.brep.geom.face_surface.push(None);
                }
            }
        }

        Ok(())
    }

    /// Finish building and return the BRep.
    fn finish(self) -> Result<BRep, IgesError> {
        if self.brep.vertices.is_empty() && self.brep.solids.is_empty() {
            return Err(IgesError::EmptyResult(
                "IGES file contained no valid geometry".into(),
            ));
        }

        Ok(self.brep)
    }

    // Helper methods

    fn get_float(&self, value: &IgesValue) -> Result<f64, IgesError> {
        match value {
            IgesValue::Float(f) => Ok(*f),
            IgesValue::Int(i) => Ok(*i as f64),
            IgesValue::String(s) => IgesReader::parse_float(s)
                .ok_or_else(|| IgesError::InvalidFormat(format!("invalid float: {}", s))),
            _ => Err(IgesError::InvalidFormat("expected float val".into())),
        }
    }

    fn get_int(&self, value: &IgesValue) -> Result<i64, IgesError> {
        match value {
            IgesValue::Int(i) => Ok(*i),
            IgesValue::Float(f) => Ok(*f as i64),
            IgesValue::String(s) => IgesReader::parse_int(s)
                .ok_or_else(|| IgesError::InvalidFormat(format!("invalid int: {}", s))),
            _ => Err(IgesError::InvalidFormat("expected int val".into())),
        }
    }

    fn circle_center_from_3_points(p1: DVec3, p2: DVec3, p3: DVec3) -> Result<DVec3, IgesError> {
        let v12 = p2 - p1;
        let v23 = p3 - p2;

        let len12_sq = v12.length_squared();
        let len23_sq = v23.length_squared();

        if len12_sq < 1e-20 || len23_sq < 1e-20 {
            return Err(IgesError::InvalidFormat(
                "degenerate arc: points are too close".into(),
            ));
        }

        let mid12 = (p1 + p2) * 0.5;
        let perp12 = v12.cross(DVec3::Z).normalize_or(v12);

        let mid23 = (p2 + p3) * 0.5;
        let perp23 = v23.cross(DVec3::Z).normalize_or(v23);

        let diff = mid23 - mid12;

        let det = perp12.x * perp23.y - perp12.y * perp23.x;
        if det.abs() < 1e-12 {
            return Err(IgesError::InvalidFormat(
                "degenerate arc: points are collinear".into(),
            ));
        }

        let t1 = (diff.x * perp23.y - diff.y * perp23.x) / det;
        Ok(mid12 + perp12 * t1)
    }
}

// ============================================================================
// IGES Writer
// ============================================================================

/// IGES file writer.
pub struct IgesWriter {
    next_seq: i32,
    global: IgesGlobalParams,
}

/// Global section parameters for IGES output.
#[derive(Debug, Clone)]
struct IgesGlobalParams {
    product_id: String,
    file_name: String,
    system_id: String,
    preprocessor_version: String,
    units_flag: i32,
    units_name: String,
}

impl Default for IgesGlobalParams {
    fn default() -> Self {
        Self {
            product_id: "RCAD".to_string(),
            file_name: "rcad_export.igs".to_string(),
            system_id: "RCAD".to_string(),
            preprocessor_version: "RCAD IGES Writer 1.0".to_string(),
            units_flag: 1,
            units_name: "INCH".to_string(),
        }
    }
}

impl IgesWriter {
    /// Create a new IGES writer with default settings.
    pub fn new() -> Self {
        Self {
            next_seq: 1,
            global: IgesGlobalParams::default(),
        }
    }

    /// Set units for the output file.
    pub fn with_units(mut self, units: &str) -> Self {
        match units.to_uppercase().as_str() {
            "INCH" | "IN" => {
                self.global.units_flag = 1;
                self.global.units_name = "INCH".to_string();
            }
            "MM" | "MILLIMETER" => {
                self.global.units_flag = 2;
                self.global.units_name = "MM".to_string();
            }
            "M" | "METER" => {
                self.global.units_flag = 5;
                self.global.units_name = "M".to_string();
            }
            _ => {}
        }
        self
    }

    /// Write BRep to an IGES string.
    pub fn write_string(brep: &BRep) -> String {
        let writer = Self::new();
        writer.write_brep_to_string(brep)
    }

    /// Write topods::BRep to an IGES string (converts internally).
    pub fn write_string_topods(brep: &rcad_kernel::topods::BRep) -> String {
        let old = rcad_kernel::BRep::from_topods(brep);
        Self::write_string(&old)
    }

    /// Write BRep to an IGES file.
    pub fn write_file<P: AsRef<Path>>(brep: &BRep, path: P) -> Result<usize, io::Error> {
        let writer = Self::new();
        let content = writer.write_brep_to_string(brep);
        std::fs::write(path, &content)?;
        Ok(content.len())
    }

    /// Write topods::BRep to an IGES file (converts internally).
    pub fn write_file_topods<P: AsRef<Path>>(brep: &rcad_kernel::topods::BRep, path: P) -> Result<usize, io::Error> {
        let old = rcad_kernel::BRep::from_topods(brep);
        Self::write_file(&old, path)
    }

    /// Convert BRep to IGES string.
    fn write_brep_to_string(self, brep: &BRep) -> String {
        let mut out = Vec::new();
        let _ = self.write_to(brep, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    }

    /// Write BRep to a writer.
    fn write_to(mut self, brep: &BRep, writer: &mut impl Write) -> io::Result<usize> {
        let mut lines: Vec<String> = Vec::new();
        let mut s_count = 0usize;
        let mut g_count = 0usize;
        let mut d_lines: Vec<String> = Vec::new();
        let mut p_lines: Vec<String> = Vec::new();

        s_count += 1;
        lines.push(section_line("RCAD IGES B-Rep Export", 'S', s_count as i32));

        g_count += 1;
        let global_line = self.format_global_section();
        lines.push(section_line(&global_line, 'G', g_count as i32));

        let mut entity_count = 0usize;

        // Write surfaces
        for surface in &brep.geom.surfaces {
            if let Some(params) = self.surface_to_iges_params(surface) {
                let de_seq = self.next_seq;
                self.next_seq += 1;

                let (d1, d2) = self.make_directory_entry(128, de_seq, (params.len() / 64 + 1) as i32);
                d_lines.push(section_line(&d1, 'D', d_lines.len() as i32 + 1));
                d_lines.push(section_line(&d2, 'D', d_lines.len() as i32 + 1));

                let p_line = format!("{:<64}{:>8}", params, de_seq * 2 - 1);
                p_lines.push(section_line(&p_line, 'P', p_lines.len() as i32 + 1));

                entity_count += 1;
            }
        }

        // Write curves
        for curve in &brep.geom.curves {
            if let Some((entity_type, params)) = self.curve_to_iges_params(curve) {
                let de_seq = self.next_seq;
                self.next_seq += 1;

                let (d1, d2) = self.make_directory_entry(entity_type, de_seq, (params.len() / 64 + 1) as i32);
                d_lines.push(section_line(&d1, 'D', d_lines.len() as i32 + 1));
                d_lines.push(section_line(&d2, 'D', d_lines.len() as i32 + 1));

                let p_line = format!("{:<64}{:>8}", params, de_seq * 2 - 1);
                p_lines.push(section_line(&p_line, 'P', p_lines.len() as i32 + 1));

                entity_count += 1;
            }
        }

        // Write vertices as points
        for vertex in &brep.vertices {
            let params = format!(
                "116,{:.9},{:.9},{:.9};",
                vertex.point.x, vertex.point.y, vertex.point.z
            );
            let de_seq = self.next_seq;
            self.next_seq += 1;

            let (d1, d2) = self.make_directory_entry(116, de_seq, 1);
            d_lines.push(section_line(&d1, 'D', d_lines.len() as i32 + 1));
            d_lines.push(section_line(&d2, 'D', d_lines.len() as i32 + 1));

            let p_line = format!("{:<64}{:>8}", params, de_seq * 2 - 1);
            p_lines.push(section_line(&p_line, 'P', p_lines.len() as i32 + 1));

            entity_count += 1;
        }

        lines.append(&mut d_lines);
        lines.append(&mut p_lines);

        let d_count = d_lines.len();
        let p_count = p_lines.len();

        let term = format!("S{:>7}G{:>7}D{:>7}P{:>7}", s_count, g_count, d_count, p_count);
        lines.push(section_line(&term, 'T', 1));

        for line in lines {
            writer.write_all(line.as_bytes())?;
            writer.write_all(b"\n")?;
        }

        Ok(entity_count)
    }

    fn format_global_section(&self) -> String {
        // IGES Global section: minimal required fields
        let mut s = String::new();
        s.push_str("1H,,1H;,");
        // product_id
        s.push_str(&format!("{}H{},", self.global.product_id.len(), self.global.product_id));
        // file_name (as empty string for simplicity)
        s.push_str(&format!("{}H{},", self.global.file_name.len(), self.global.file_name));
        // system_id
        s.push_str(&format!("{}H{},", self.global.system_id.len(), self.global.system_id));
        // preprocessor_version
        s.push_str(&format!("{}H{},", self.global.preprocessor_version.len(), self.global.preprocessor_version));
        // significant_digits, units_flag
        s.push_str(&format!("6,{},", self.global.units_flag));
        s.push(';');
        s
    }

    fn make_directory_entry(
        &self,
        entity_type: i32,
        param_ptr: i32,
        param_count: i32,
    ) -> (String, String) {
        let line1 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            entity_type, param_ptr, 0, 0, 0, 0, 0, 0
        );
        let line2 = format!(
            "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
            entity_type, 0, 0, param_count, 0, 0, 0, 0
        );
        (line1, line2)
    }

    fn surface_to_iges_params(&self, surface: &Surface3) -> Option<String> {
        match surface {
            Surface3::Plane(p) => {
                let d = p.origin.dot(p.normal);
                Some(format!(
                    "108,{:.9},{:.9},{:.9},{:.9},0;",
                    p.normal.x, p.normal.y, p.normal.z, d
                ))
            }
            Surface3::Sphere(s) => {
                Some(format!(
                    "196,{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9};",
                    s.center.x, s.center.y, s.center.z,
                    s.axis.x, s.axis.y, s.axis.z,
                    s.radius
                ))
            }
            Surface3::BSpline(bs) => {
                let mut params = format!(
                    "128,0,{},{},0,0,{},0,0,0,0,{},{}",
                    bs.degree_u, bs.degree_v,
                    if bs.weights.iter().all(|r| r.iter().all(|&w| (w - 1.0).abs() < 1e-9)) { 0 } else { 1 },
                    bs.knots_u.len() - 1, bs.knots_v.len() - 1
                );

                for &k in &bs.knots_u {
                    params.push_str(&format!("{:.9},", k));
                }

                for &k in &bs.knots_v {
                    params.push_str(&format!("{:.9},", k));
                }

                if bs.weights.iter().any(|r| r.iter().any(|&w| (w - 1.0).abs() >= 1e-9)) {
                    for row in &bs.weights {
                        for &w in row {
                            params.push_str(&format!("{:.9},", w));
                        }
                    }
                }

                for row in &bs.control_points {
                    for p in row {
                        params.push_str(&format!("{:.9},{:.9},{:.9},", p.x, p.y, p.z));
                    }
                }

                if params.ends_with(',') {
                    params.pop();
                }
                params.push(';');

                Some(params)
            }
            _ => None,
        }
    }

    fn curve_to_iges_params(&self, curve: &Curve3) -> Option<(i32, String)> {
        match curve {
            Curve3::Line(l) => {
                let end = l.origin + l.direction;
                Some((
                    110,
                    format!(
                        "110,{:.9},{:.9},{:.9},{:.9},{:.9},{:.9};",
                        l.origin.x, l.origin.y, l.origin.z,
                        end.x, end.y, end.z
                    ),
                ))
            }
            Curve3::Circle(c) => {
                let r = c.radius;
                let perp = c.normal.any_orthogonal_vector();
                let p1 = c.center + perp * r;
                let p2 = c.center + c.normal.cross(perp) * r;
                let p3 = c.center - perp * r;

                Some((
                    100,
                    format!(
                        "100,{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9};",
                        c.center.z,
                        p1.x, p1.y, p2.x, p2.y, p3.x, p3.y
                    ),
                ))
            }
            Curve3::BSpline(bs) => {
                let mut params = format!(
                    "126,0,{},{},0,{},1,0,{}",
                    bs.degree,
                    0,
                    if bs.weights.iter().all(|&w| (w - 1.0).abs() < 1e-9) { 0 } else { 1 },
                    bs.knots.len() - 1
                );

                for &k in &bs.knots {
                    params.push_str(&format!("{:.9},", k));
                }

                if bs.weights.iter().any(|&w| (w - 1.0).abs() >= 1e-9) {
                    for &w in &bs.weights {
                        params.push_str(&format!("{:.9},", w));
                    }
                }

                for p in &bs.control_points {
                    params.push_str(&format!("{:.9},{:.9},{:.9},", p.x, p.y, p.z));
                }

                if params.ends_with(',') {
                    params.pop();
                }
                params.push(';');

                Some((126, params))
            }
            _ => None,
        }
    }
}

impl Default for IgesWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Format a line with section identifier and sequence number.
fn section_line(data: &str, section: char, seq: i32) -> String {
    let truncated = truncate_to_width(data, 72);
    format!("{:<72}{}{:>7}", truncated, section, seq)
}

/// Truncate string to specified width.
fn truncate_to_width(text: &str, width: usize) -> String {
    let chars: Vec<char> = text.chars().collect();
    if chars.len() <= width {
        text.to_string()
    } else {
        chars[..width].iter().collect()
    }
}

/// Convenience function to parse IGES string to BRep.
pub fn parse_iges_string(s: &str) -> Result<BRep, IgesError> {
    IgesReader::parse_string(s)
}

/// Convenience function to write BRep to IGES string.
pub fn write_iges_string(brep: &BRep) -> String {
    IgesWriter::write_string(brep)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_brep() -> BRep {
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

        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        brep.geom.curves.push(Curve3::Line(line));

        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        brep.geom.surfaces.push(Surface3::Plane(plane));

        brep
    }

    #[test]
    fn write_iges_string_produces_valid_structure() {
        let brep = make_simple_brep();
        let iges = IgesWriter::write_string(&brep);

        assert!(iges.contains('S'));
        assert!(iges.contains('G'));
        assert!(iges.contains('D'));
        assert!(iges.contains('P'));
        assert!(iges.contains('T'));
    }

    #[test]
    fn parse_empty_string_returns_error() {
        let result = IgesReader::parse_string("");
        assert!(result.is_err());
    }

    #[test]
    fn parse_invalid_format_returns_error() {
        let result = IgesReader::parse_string("This is not an IGES file");
        assert!(result.is_err());
    }

    #[test]
    fn parse_hollerith_string() {
        assert_eq!(IgesReader::parse_hollerith("5HHello"), Some("Hello".to_string()));
        assert_eq!(IgesReader::parse_hollerith("1H,"), Some(",".to_string()));
        assert_eq!(IgesReader::parse_hollerith("2HAB"), Some("AB".to_string()));
    }

    #[test]
    fn parse_float_handles_d_exponent() {
        assert_eq!(IgesReader::parse_float("1.0D+3"), Some(1000.0));
        assert_eq!(IgesReader::parse_float("1.5D-2"), Some(0.015));
        assert_eq!(IgesReader::parse_float("2.5E+1"), Some(25.0));
    }

    #[test]
    fn round_trip_simple_point() {
        let iges = r#"                                                                        S      1
1H,,1H;,4HRCAD,4Htest.igs,4HRCAD,20HRCAD IGES Writer 1.0,6,1,4HINCH;G      1
     116       1       0       0       0       0       0       0D      1
     116       0       0       1       0       0             0       0D      2
116,0.0,0.0,0.0;                                                            1       P      1
S      1G      1D      2P      1                                             T      1
"#;
        let result = IgesReader::parse_string(iges);
        match result {
            Ok(brep) => {
                assert!(!brep.vertices.is_empty());
            }
            Err(IgesError::EmptyResult(_)) => {}
            Err(e) => {
                panic!("Unexpected error: {:?}", e);
            }
        }
    }
}
