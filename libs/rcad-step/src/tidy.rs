//! StepTidy — STEP entity deduplication (OCCT TKDESTEP-aligned).
//!
//! Post-processes generated STEP text to find and merge duplicate entities:
//! - CARTESIAN_POINT:  same (x, y, z)
//! - DIRECTION:        same (dx, dy, dz)
//! - VECTOR:           same orientation + magnitude
//! - AXIS2_PLACEMENT_3D: same location + axis + ref_direction
//! - LINE:             same point + vector
//! - CIRCLE:           same placement + radius
//! - PLANE:            same placement
//!
//! OCCT source: src/DataExchange/TKDESTEP/GTests/StepTidy_*.cxx
//! These reducers follow the same logic but adapted for rcad-step's
//! direct text generation rather than OCCT's StepData_StepModel.

use std::collections::{BTreeMap, HashMap, HashSet};

/// A single STEP entity record parsed from the DATA section.
#[derive(Debug, Clone)]
struct EntityRecord {
    id: u64,
    type_name: String,
    body: String, // raw parameter string between '(' and ')' before ';'
}

/// Parsed representation of the STEP DATA section.
#[derive(Debug, Default)]
struct StepDataSection {
    records: Vec<EntityRecord>,
    id_to_idx: HashMap<u64, usize>,
}

/// Parsed entity bodies for each known type.
#[derive(Debug, Default)]
struct ResolvedEntities {
    cartesian_points: HashMap<u64, CartesianPointData>,
    directions: HashMap<u64, DirectionData>,
    vectors: HashMap<u64, VectorData>,
    axis2_placements: HashMap<u64, Axis2PlacementData>,
    lines: HashMap<u64, LineData>,
    circles: HashMap<u64, CircleData>,
    planes: HashMap<u64, PlaneData>,
}

#[derive(Debug, Clone, PartialEq)]
struct CartesianPointData {
    name: String,
    coords: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
struct DirectionData {
    name: String,
    ratios: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
struct VectorData {
    name: String,
    orientation_ref: u64,
    magnitude: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct Axis2PlacementData {
    name: String,
    location_ref: u64,
    axis_ref: u64,
    ref_direction_ref: u64,
}

#[derive(Debug, Clone, PartialEq)]
struct LineData {
    name: String,
    point_ref: u64,
    vector_ref: u64,
}

#[derive(Debug, Clone, PartialEq)]
struct CircleData {
    name: String,
    placement_ref: u64,
    radius: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct PlaneData {
    name: String,
    placement_ref: u64,
}

/// Parse the STEP DATA section into records.
fn parse_data_section(step_text: &str) -> Option<StepDataSection> {
    // Find DATA ... ENDSEC
    let data_start = step_text.find("DATA;")?;
    let data_start = data_start + 5; // skip "DATA;"
    let data_end = step_text[data_start..].find("ENDSEC").map(|p| data_start + p)?;
    let data_section = &step_text[data_start..data_end];

    let mut records = Vec::new();
    let mut pos = 0;
    let bytes = data_section.as_bytes();

    while pos < bytes.len() {
        // Skip whitespace/newlines
        while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos >= bytes.len() {
            break;
        }

        // Expect '#'
        if bytes[pos] != b'#' {
            pos += 1;
            continue;
        }
        pos += 1;

        // Read entity ID
        let id_start = pos;
        while pos < bytes.len() && bytes[pos].is_ascii_digit() {
            pos += 1;
        }
        if pos == id_start {
            continue;
        }
        let id: u64 = data_section[id_start..pos].parse().ok()?;

        // Expect '='
        while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos >= bytes.len() || bytes[pos] != b'=' {
            continue;
        }
        pos += 1;
        while pos < bytes.len() && bytes[pos].is_ascii_whitespace() {
            pos += 1;
        }

        // Read type name (until '(')
        let type_start = pos;
        while pos < bytes.len() && bytes[pos] != b'(' {
            pos += 1;
        }
        if pos >= bytes.len() {
            continue;
        }
        let type_name = data_section[type_start..pos].trim().to_string();
        pos += 1; // skip '('

        // Read body (balanced parentheses)
        let body_start = pos;
        let mut depth = 1;
        while pos < bytes.len() && depth > 0 {
            if bytes[pos] == b'(' {
                depth += 1;
            } else if bytes[pos] == b')' {
                depth -= 1;
            }
            pos += 1;
        }
        let body = data_section[body_start..pos - 1].trim().to_string();

        // Expect ';'
        while pos < bytes.len() && bytes[pos] != b';' {
            pos += 1;
        }
        if pos < bytes.len() {
            pos += 1; // skip ';'
        }

        records.push(EntityRecord {
            id,
            type_name,
            body,
        });
    }

    let id_to_idx: HashMap<u64, usize> = records
        .iter()
        .enumerate()
        .map(|(i, r)| (r.id, i))
        .collect();

    if records.is_empty() {
        None
    } else {
        Some(StepDataSection { records, id_to_idx })
    }
}

/// Parse a quoted string: `'name'`  ->  String.
fn parse_quoted_string(s: &str) -> Option<String> {
    let s = s.trim();
    if s.starts_with('\'') && (s.len() >= 2) {
        // Find matching end quote (handle escaped single quotes)
        let mut chars = s[1..].chars();
        let mut result = String::new();
        loop {
            match chars.next() {
                None => return None,
                Some('\'') => {
                    // Check for escaped quote ''
                    let next = chars.as_str().chars().next();
                    if next == Some('\'') {
                        result.push('\'');
                        chars.next(); // skip second '
                    } else {
                        break;
                    }
                }
                Some(c) => result.push(c),
            }
        }
        Some(result)
    } else {
        None
    }
}

/// Parse a reference: `#123`  ->  123.
fn parse_ref(s: &str) -> Option<u64> {
    let s = s.trim();
    if s.starts_with('#') {
        s[1..].trim().parse().ok()
    } else {
        None
    }
}

/// Parse a tuple of floats: `(1.0, 2.0, 3.0)`  ->  [f64; 3].
fn parse_f64_tuple_3(s: &str) -> Option<[f64; 3]> {
    let s = s.trim();
    if !s.starts_with('(') || !s.ends_with(')') {
        return None;
    }
    let inner = s[1..s.len() - 1].trim();
    let parts: Vec<&str> = inner.split(',').collect();
    if parts.len() != 3 {
        return None;
    }
    Some([
        parts[0].trim().parse().ok()?,
        parts[1].trim().parse().ok()?,
        parts[2].trim().parse().ok()?,
    ])
}

/// Parse a single float from a string.
fn parse_f64(s: &str) -> Option<f64> {
    s.trim().parse().ok()
}

/// Parse entity body for CARTESIAN_POINT: `'name',(x,y,z)`.
fn parse_cartesian_point(body: &str) -> Option<CartesianPointData> {
    let body = body.trim();
    // Find the comma separating name from coords: `',(` means `'` `,` `(`
    let name_end = body.find("',(")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    // Skip past `',` to reach `(x,y,z)`
    let coords_str = &body[name_end + 2..];
    let coords = parse_f64_tuple_3(coords_str)?;
    Some(CartesianPointData { name, coords })
}

/// Parse entity body for DIRECTION: `'name',(dx,dy,dz)`.
fn parse_direction(body: &str) -> Option<DirectionData> {
    let body = body.trim();
    let name_end = body.find("',(")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let ratios_str = &body[name_end + 2..];
    let ratios = parse_f64_tuple_3(ratios_str)?;
    Some(DirectionData { name, ratios })
}

/// Parse entity body for VECTOR: `'name',#ref,magnitude`.
fn parse_vector(body: &str) -> Option<VectorData> {
    let body = body.trim();
    // 'name',#ref,magnitude
    let name_end = body.find("',")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let rest = body[name_end + 2..].trim();
    let ref_end = rest.find(',')?;
    let orientation_ref = parse_ref(&rest[..ref_end])?;
    let magnitude = parse_f64(&rest[ref_end + 1..])?;
    Some(VectorData {
        name,
        orientation_ref,
        magnitude,
    })
}

/// Parse entity body for AXIS2_PLACEMENT_3D: `'name',#loc,#axis,#ref_dir`.
fn parse_axis2_placement(body: &str) -> Option<Axis2PlacementData> {
    let body = body.trim();
    let name_end = body.find("',")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let rest = body[name_end + 2..].trim();
    let parts: Vec<&str> = rest.split(',').collect();
    if parts.len() < 3 {
        return None;
    }
    let location_ref = parse_ref(parts[0])?;
    let axis_ref = parse_ref(parts[1])?;
    let ref_direction_ref = parse_ref(parts[2])?;
    Some(Axis2PlacementData {
        name,
        location_ref,
        axis_ref,
        ref_direction_ref,
    })
}

/// Parse entity body for LINE: `'name',#point,#vector`.
fn parse_line(body: &str) -> Option<LineData> {
    let body = body.trim();
    let name_end = body.find("',")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let rest = body[name_end + 2..].trim();
    let parts: Vec<&str> = rest.split(',').collect();
    if parts.len() < 2 {
        return None;
    }
    let point_ref = parse_ref(parts[0])?;
    let vector_ref = parse_ref(parts[1])?;
    Some(LineData {
        name,
        point_ref,
        vector_ref,
    })
}

/// Parse entity body for CIRCLE: `'name',#placement,radius`.
fn parse_circle(body: &str) -> Option<CircleData> {
    let body = body.trim();
    let name_end = body.find("',")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let rest = body[name_end + 2..].trim();
    let ref_end = rest.find(',')?;
    let placement_ref = parse_ref(&rest[..ref_end])?;
    let radius = parse_f64(&rest[ref_end + 1..])?;
    Some(CircleData {
        name,
        placement_ref,
        radius,
    })
}

/// Parse entity body for PLANE: `'name',#placement`.
fn parse_plane(body: &str) -> Option<PlaneData> {
    let body = body.trim();
    let name_end = body.find("',")?;
    let name = parse_quoted_string(&body[..name_end + 1])?;
    let rest = body[name_end + 2..].trim();
    let placement_ref = parse_ref(rest)?;
    Some(PlaneData {
        name,
        placement_ref,
    })
}

/// Resolve parsed entities from raw records.
fn resolve_entities(section: &StepDataSection) -> ResolvedEntities {
    let mut entities = ResolvedEntities::default();

    for record in &section.records {
        match record.type_name.as_str() {
            "CARTESIAN_POINT" => {
                if let Some(data) = parse_cartesian_point(&record.body) {
                    entities.cartesian_points.insert(record.id, data);
                }
            }
            "DIRECTION" => {
                if let Some(data) = parse_direction(&record.body) {
                    entities.directions.insert(record.id, data);
                }
            }
            "VECTOR" => {
                if let Some(data) = parse_vector(&record.body) {
                    entities.vectors.insert(record.id, data);
                }
            }
            "AXIS2_PLACEMENT_3D" => {
                if let Some(data) = parse_axis2_placement(&record.body) {
                    entities.axis2_placements.insert(record.id, data);
                }
            }
            "LINE" => {
                if let Some(data) = parse_line(&record.body) {
                    entities.lines.insert(record.id, data);
                }
            }
            "CIRCLE" => {
                if let Some(data) = parse_circle(&record.body) {
                    entities.circles.insert(record.id, data);
                }
            }
            "PLANE" => {
                if let Some(data) = parse_plane(&record.body) {
                    entities.planes.insert(record.id, data);
                }
            }
            _ => {}
        }
    }

    entities
}

/// Find the canonical ID for an entity, chasing replacement chain.
/// `replacements` maps old_id -> canonical_id.
fn canonical(replacements: &HashMap<u64, u64>, id: u64) -> u64 {
    let mut cur = id;
    while let Some(&next) = replacements.get(&cur) {
        if next == cur {
            break;
        }
        cur = next;
    }
    cur
}

/// Build a string key for a CartesianPoint that avoids f64 Hash issues.
fn cp_key(data: &CartesianPointData) -> String {
    format!("{}:{:.9}:{:.9}:{:.9}", data.name, data.coords[0], data.coords[1], data.coords[2])
}

/// Build a string key for a Direction.
fn dir_key(data: &DirectionData) -> String {
    format!("{}:{:.9}:{:.9}:{:.9}", data.name, data.ratios[0], data.ratios[1], data.ratios[2])
}

/// Build deduplication replacements for CARTESIAN_POINT entities.
fn dedup_cartesian_points(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.cartesian_points.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.cartesian_points.get(&id) {
            let key = cp_key(data);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for DIRECTION entities.
fn dedup_directions(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.directions.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.directions.get(&id) {
            let key = dir_key(data);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for VECTOR entities.
/// Two vectors are equal if they have the same name, reference the same (canonical)
/// direction, and have the same magnitude.
fn dedup_vectors(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
    dir_replacements: &HashMap<u64, u64>,
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.vectors.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.vectors.get(&id) {
            let canon_dir = canonical(dir_replacements, data.orientation_ref);
            let key = format!("{}:{}:{:.9}", data.name, canon_dir, data.magnitude);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for AXIS2_PLACEMENT_3D entities.
fn dedup_axis2_placements(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
    pt_replacements: &HashMap<u64, u64>,
    dir_replacements: &HashMap<u64, u64>,
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.axis2_placements.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.axis2_placements.get(&id) {
            let canon_loc = canonical(pt_replacements, data.location_ref);
            let canon_axis = canonical(dir_replacements, data.axis_ref);
            let canon_ref = canonical(dir_replacements, data.ref_direction_ref);
            let key = format!("{}:{}:{}:{}", data.name, canon_loc, canon_axis, canon_ref);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for LINE entities.
fn dedup_lines(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
    pt_replacements: &HashMap<u64, u64>,
    vec_replacements: &HashMap<u64, u64>,
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.lines.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.lines.get(&id) {
            let canon_pt = canonical(pt_replacements, data.point_ref);
            let canon_vec = canonical(vec_replacements, data.vector_ref);
            let key = format!("{}:{}:{}", data.name, canon_pt, canon_vec);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for CIRCLE entities.
fn dedup_circles(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
    a2p_replacements: &HashMap<u64, u64>,
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.circles.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.circles.get(&id) {
            let canon_a2p = canonical(a2p_replacements, data.placement_ref);
            let key = format!("{}:{}:{:.9}", data.name, canon_a2p, data.radius);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Build deduplication replacements for PLANE entities.
fn dedup_planes(
    entities: &ResolvedEntities,
    _records: &[EntityRecord],
    a2p_replacements: &HashMap<u64, u64>,
) -> HashMap<u64, u64> {
    let mut replacements = HashMap::new();
    let mut sig_to_id: HashMap<String, u64> = HashMap::new();

    let mut sorted: Vec<u64> = entities.planes.keys().copied().collect();
    sorted.sort();

    for &id in &sorted {
        if let Some(data) = entities.planes.get(&id) {
            let canon_a2p = canonical(a2p_replacements, data.placement_ref);
            let key = format!("{}:{}", data.name, canon_a2p);
            if let Some(&canonical_id) = sig_to_id.get(&key) {
                replacements.insert(id, canonical_id);
            } else {
                sig_to_id.insert(key, id);
            }
        }
    }

    replacements
}

/// Apply replacements to a STEP text.
/// Replaces all `#old_id` references with `#new_id` and removes
/// the definitions of replaced entities.
fn apply_replacements(step_text: &str, replacements: &HashMap<u64, u64>) -> String {
    if replacements.is_empty() {
        return step_text.to_string();
    }

    // Build a set of IDs that should be removed (all replaced entities)
    let removed: HashSet<u64> = replacements.keys().copied().collect();

    // Sort replacements by ID (descending) to avoid reference issues
    let mut sorted_repl: Vec<(u64, u64)> = replacements.iter().map(|(&k, &v)| (k, v)).collect();
    sorted_repl.sort_by(|a, b| b.0.cmp(&a.0));

    // Process line by line
    let mut result = String::with_capacity(step_text.len());
    // Remove definitions of replaced entities
    for line in step_text.lines() {
        let trimmed = line.trim();
        // Check if this line defines a removed entity: #ID = ...
        let is_removed_def = trimmed.starts_with('#')
            && trimmed.contains("=")
            && trimmed.ends_with(';')
            && {
                let id_str = trimmed[1..]
                    .split('=')
                    .next()
                    .unwrap_or("")
                    .trim();
                id_str.parse::<u64>().map(|id| removed.contains(&id)).unwrap_or(false)
            };

        if is_removed_def {
            continue; // skip this definition
        }

        result.push_str(line);
        result.push('\n');
    }

    // Now replace all references #old_id with #new_id in the remaining text
    // Sort by old_id descending to avoid nested replacement issues
    let mut text = result;
    for (old_id, new_id) in &sorted_repl {
        if old_id == new_id {
            continue;
        }
        let old_str = format!("#{}", old_id);
        let new_str = format!("#{}", new_id);
        text = text.replace(&old_str, &new_str);
    }

    text
}

/// Result of the deduplication process.
#[derive(Debug, Default)]
pub struct TidyReport {
    /// Number of entities removed (total across all types).
    pub removed_count: usize,
    /// Number of CARTESIAN_POINT entities removed.
    pub removed_cartesian_points: usize,
    /// Number of DIRECTION entities removed.
    pub removed_directions: usize,
    /// Number of VECTOR entities removed.
    pub removed_vectors: usize,
    /// Number of AXIS2_PLACEMENT_3D entities removed.
    pub removed_axis2_placements: usize,
    /// Number of LINE entities removed.
    pub removed_lines: usize,
    /// Number of CIRCLE entities removed.
    pub removed_circles: usize,
    /// Number of PLANE entities removed.
    pub removed_planes: usize,
}

/// Deduplicate entities in a STEP text.
///
/// This function parses the DATA section of a STEP file,
/// identifies duplicate entities (same type and same parameters),
/// and merges them by replacing all references to the duplicate
/// with references to the canonical entity.
pub fn deduplicate(step_text: &str) -> (String, TidyReport) {
    let Some(section) = parse_data_section(step_text) else {
        return (step_text.to_string(), TidyReport::default());
    };

    let entities = resolve_entities(&section);

    // Phase 1: dedup leaf types (no dependencies on other tidyable types)
    let cp_repl = dedup_cartesian_points(&entities, &section.records);
    let dir_repl = dedup_directions(&entities, &section.records);

    // Phase 2: dedup types that depend on leaves
    let vec_repl = dedup_vectors(&entities, &section.records, &dir_repl);
    let a2p_repl = dedup_axis2_placements(&entities, &section.records, &cp_repl, &dir_repl);

    // Phase 3: dedup types that depend on phase 2
    let line_repl = dedup_lines(&entities, &section.records, &cp_repl, &vec_repl);
    let circle_repl = dedup_circles(&entities, &section.records, &a2p_repl);
    let plane_repl = dedup_planes(&entities, &section.records, &a2p_repl);

    // Merge all replacements
    let mut all_replacements: HashMap<u64, u64> = HashMap::new();
    all_replacements.extend(cp_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(dir_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(vec_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(a2p_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(line_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(circle_repl.iter().map(|(&k, &v)| (k, v)));
    all_replacements.extend(plane_repl.iter().map(|(&k, &v)| (k, v)));

    let report = TidyReport {
        removed_count: all_replacements.len(),
        removed_cartesian_points: cp_repl.len(),
        removed_directions: dir_repl.len(),
        removed_vectors: vec_repl.len(),
        removed_axis2_placements: a2p_repl.len(),
        removed_lines: line_repl.len(),
        removed_circles: circle_repl.len(),
        removed_planes: plane_repl.len(),
    };

    let tidy_text = apply_replacements(step_text, &all_replacements);

    (tidy_text, report)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helper: build a minimal STEP file with the given DATA section ──

    fn make_step(data_entities: &[&str]) -> String {
        let header = "\
ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('RCAD test'),'2;1');
FILE_NAME('test.step','2026-04-02T00:00:00',(''),(''),'RCAD','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
";
        let footer = "ENDSEC;\nEND-ISO-10303-21;\n";
        let mut out = header.to_string();
        for e in data_entities {
            out.push_str(e);
            out.push('\n');
        }
        out.push_str(footer);
        out
    }

    // ── Helper to extract entity bodies ──

    fn entity_ids_in_text(text: &str) -> Vec<u64> {
        let mut ids = Vec::new();
        'lines: for line in text.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with('#') && trimmed.contains("=") && trimmed.ends_with(';') {
                let id_str = trimmed[1..].split('=').next().unwrap_or("").trim();
                if let Ok(id) = id_str.parse::<u64>() {
                    ids.push(id);
                }
            }
        }
        ids
    }

    fn entity_count_after(text: &str) -> usize {
        entity_ids_in_text(text).len()
    }

    fn count_type(text: &str, type_name: &str) -> usize {
        text.lines()
            .filter(|line| line.trim().contains(type_name))
            .count()
    }

    // ── CARTESIAN_POINT Reducer ──

    #[test]
    fn cp_different_names_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('FirstPt',(0.000000000,0.000000000,0.000000000));",
            "#2 = CARTESIAN_POINT('SecondPt',(0.000000000,0.000000000,0.000000000));",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_cartesian_points, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "CARTESIAN_POINT"), 2);
    }

    #[test]
    fn cp_identical_coords_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,1.000000000));",
            "#2 = CARTESIAN_POINT('',(0.000000000,0.000000000,1.000000000));",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_cartesian_points, 1, "identical points should be merged");
        assert_eq!(count_type(&tidy, "CARTESIAN_POINT"), 1);
    }

    #[test]
    fn cp_different_coords_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = CARTESIAN_POINT('',(1.000000000,2.000000000,3.000000000));",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_cartesian_points, 0);
        assert_eq!(count_type(&tidy, "CARTESIAN_POINT"), 2);
    }

    // ── DIRECTION Reducer ──

    #[test]
    fn dir_different_names_no_merge() {
        let step = make_step(&[
            "#1 = DIRECTION('dir1',(0.000000000,0.000000000,1.000000000));",
            "#2 = DIRECTION('dir2',(0.000000000,0.000000000,1.000000000));",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_directions, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "DIRECTION"), 2);
    }

    #[test]
    fn dir_identical_ratios_merge() {
        let step = make_step(&[
            "#1 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_directions, 1);
        assert_eq!(count_type(&tidy, "DIRECTION"), 1);
    }

    // ── VECTOR Reducer ──

    #[test]
    fn vec_different_names_no_merge() {
        let step = make_step(&[
            "#1 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#2 = VECTOR('vec1',#1,1.000000000);",
            "#3 = VECTOR('vec2',#1,1.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_vectors, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "VECTOR"), 2);
    }

    #[test]
    fn vec_identical_content_merge() {
        let step = make_step(&[
            "#1 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#2 = VECTOR('',#1,1.000000000);",
            "#3 = VECTOR('',#1,1.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_vectors, 1);
        assert_eq!(count_type(&tidy, "VECTOR"), 1);
    }

    #[test]
    fn vec_different_magnitude_no_merge() {
        let step = make_step(&[
            "#1 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#2 = VECTOR('',#1,1.000000000);",
            "#3 = VECTOR('',#1,2.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_vectors, 0);
        assert_eq!(count_type(&tidy, "VECTOR"), 2);
    }

    // ── AXIS2_PLACEMENT_3D Reducer ──

    #[test]
    fn a2p_different_names_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('Axis1',#1,#2,#3);",
            "#5 = AXIS2_PLACEMENT_3D('Axis2',#1,#2,#3);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_axis2_placements, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "AXIS2_PLACEMENT_3D"), 2);
    }

    #[test]
    fn a2p_identical_content_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_axis2_placements, 1);
        assert_eq!(count_type(&tidy, "AXIS2_PLACEMENT_3D"), 1);
    }

    #[test]
    fn a2p_different_location_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = CARTESIAN_POINT('',(1.000000000,0.000000000,0.000000000));",
            "#3 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#4 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#5 = AXIS2_PLACEMENT_3D('',#1,#3,#4);",
            "#6 = AXIS2_PLACEMENT_3D('',#2,#3,#4);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_axis2_placements, 0);
    }

    // ── LINE Reducer ──

    #[test]
    fn line_different_names_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = VECTOR('',#2,1.000000000);",
            "#4 = LINE('Line1',#1,#3);",
            "#5 = LINE('Line2',#1,#3);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_lines, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "LINE"), 2);
    }

    #[test]
    fn line_identical_content_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = VECTOR('',#2,1.000000000);",
            "#4 = LINE('',#1,#3);",
            "#5 = LINE('',#1,#3);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_lines, 1);
        assert_eq!(count_type(&tidy, "LINE"), 1);
    }

    // ── CIRCLE Reducer ──

    #[test]
    fn circle_different_names_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = CIRCLE('Circle1',#4,1.000000000);",
            "#6 = CIRCLE('Circle2',#4,1.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_circles, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "CIRCLE"), 2);
    }

    #[test]
    fn circle_identical_content_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = CIRCLE('',#4,1.000000000);",
            "#6 = CIRCLE('',#4,1.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_circles, 1);
        assert_eq!(count_type(&tidy, "CIRCLE"), 1);
    }

    #[test]
    fn circle_different_radius_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = CIRCLE('',#4,1.000000000);",
            "#6 = CIRCLE('',#4,2.000000000);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_circles, 0);
        assert_eq!(count_type(&tidy, "CIRCLE"), 2);
    }

    // ── PLANE Reducer ──

    #[test]
    fn plane_different_names_no_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = PLANE('Plane1',#4);",
            "#6 = PLANE('Plane2',#4);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_planes, 0, "different names should NOT be merged");
        assert_eq!(count_type(&tidy, "PLANE"), 2);
    }

    #[test]
    fn plane_identical_content_merge() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = PLANE('',#4);",
            "#6 = PLANE('',#4);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_planes, 1);
        assert_eq!(count_type(&tidy, "PLANE"), 1);
    }

    // ── Transitive Dedup: duplicate CartesianPoint → duplicate Axis2Placement → dedup ──

    #[test]
    fn transitive_cp_dedup_enables_a2p_dedup() {
        // Two Axis2Placement3d that use different but equal CartesianPoints
        // After CP dedup, the A2P should also be dedupable
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#3 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#4 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#5 = AXIS2_PLACEMENT_3D('',#1,#3,#4);",
            "#6 = AXIS2_PLACEMENT_3D('',#2,#3,#4);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_cartesian_points, 1);
        assert_eq!(report.removed_axis2_placements, 1);
        assert_eq!(count_type(&tidy, "CARTESIAN_POINT"), 1);
        assert_eq!(count_type(&tidy, "AXIS2_PLACEMENT_3D"), 1);
    }

    // ── Reference replacement integrity ──

    #[test]
    fn references_are_updated_correctly() {
        // Create two Axis2Placement3d (one is duplicate) and two PLANEs that reference them.
        // After dedup, all PLANEs should reference the same canonical A2P
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = DIRECTION('',(0.000000000,1.000000000,0.000000000));",
            "#4 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#5 = AXIS2_PLACEMENT_3D('',#1,#2,#3);",
            "#6 = PLANE('Plane1',#4);",
            "#7 = PLANE('Plane2',#5);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_axis2_placements, 1);
        assert_eq!(report.removed_planes, 0); // planes ref diff A2P names -> not merged

        // Both planes should still exist
        assert_eq!(count_type(&tidy, "PLANE"), 2);
        // Verify the text doesn't contain #5 anymore
        assert!(!tidy.contains("#5 ="), "removed entity #5 should not appear in output");
    }

    // ── Empty/no-op ──

    #[test]
    fn no_duplicates_does_nothing() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#3 = AXIS2_PLACEMENT_3D('',#1,#2,#1);",
        ]);
        let (tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_count, 0);
        assert_eq!(tidy, step);
    }

    #[test]
    fn empty_step_text_no_crash() {
        let text = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\nENDSEC;\nEND-ISO-10303-21;\n";
        let (tidy, report) = deduplicate(text);
        assert_eq!(report.removed_count, 0);
        assert_eq!(tidy, text);
    }

    // ── Full deduplication pipeline ──

    #[test]
    fn dedup_pipeline_removes_all_duplicates() {
        let step = make_step(&[
            "#1 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#2 = CARTESIAN_POINT('',(0.000000000,0.000000000,0.000000000));",
            "#3 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#4 = DIRECTION('',(0.000000000,0.000000000,1.000000000));",
            "#5 = VECTOR('',#3,1.000000000);",
            "#6 = VECTOR('',#3,1.000000000);",
            "#7 = AXIS2_PLACEMENT_3D('',#1,#3,#4);",
            "#8 = AXIS2_PLACEMENT_3D('',#2,#4,#3);",
            "#9 = LINE('',#1,#5);",
            "#10 = LINE('',#2,#6);",
            "#11 = CIRCLE('',#7,1.000000000);",
            "#12 = CIRCLE('',#8,1.000000000);",
            "#13 = PLANE('',#7);",
            "#14 = PLANE('',#8);",
        ]);
        let (_tidy, report) = deduplicate(&step);
        assert_eq!(report.removed_cartesian_points, 1);
        assert_eq!(report.removed_directions, 1);
        assert_eq!(report.removed_vectors, 1);
        assert_eq!(report.removed_axis2_placements, 1);
        assert_eq!(report.removed_lines, 1);
        assert_eq!(report.removed_circles, 1);
        assert_eq!(report.removed_planes, 1);
        assert_eq!(report.removed_count, 7);
    }
}
