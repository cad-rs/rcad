//! STEP assembly writer and reader.
//!
//! Writes a multi-component assembly where each component is a separate BRep
//! with an optional full affine transform and name.  Produces a STEP file with a
//! full NEXT_ASSEMBLY_USAGE_OCCURRENCE hierarchy so that importers (FreeCAD,
//! OCCT, etc.) can reconstruct the tree.
//!
//! # STEP assembly structure
//!
//! ```text
//! PRODUCT (assembly root)
//!   └─ PRODUCT_DEFINITION (assembly)
//!        └─ NEXT_ASSEMBLY_USAGE_OCCURRENCE → PRODUCT (component i)
//!                                               └─ PRODUCT_DEFINITION (component i)
//!                                                    └─ shape representation (geometry)
//! ```
//!
//! **Transform strategy**: transforms are baked into vertex/geometry coordinates
//! (`BRep::apply_transform`) rather than emitting `ITEM_DEFINED_TRANSFORMATION`
//! entities. This maximises compatibility with STEP readers that do not support
//! AP214 transformation entities.

use std::collections::{HashMap, HashSet};

use glam::{DAffine3, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::appearance::StepColor;

use crate::StepError;
use crate::writer::{ExportSelection, StepWriter};

// ─────────────────────────────────────────────────────────────────────────────
// AssemblyComponent
// ─────────────────────────────────────────────────────────────────────────────

/// A single component in an assembly.
///
/// The [`transform`][AssemblyComponent::transform] field is a full affine
/// transform (translation, rotation, uniform or non-uniform scale).  When
/// writing to STEP the transform is **baked** into vertex coordinates via
/// [`BRep::apply_transform`].
#[derive(Clone)]
pub struct AssemblyComponent {
    /// Human-readable part name.
    pub name: String,
    /// The geometry.
    pub brep: BRep,
    /// Full affine transform for this component (default: identity).
    pub transform: DAffine3,
    /// Optional RGB color for this component's faces.
    pub color: Option<rcad_kernel::appearance::Color>,
}

impl AssemblyComponent {
    /// Create a new component with an identity transform.
    pub fn new(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            brep,
            transform: DAffine3::IDENTITY,
            color: None,
        }
    }

    /// Set a full affine transform (replaces any previously set transform).
    pub fn with_transform(mut self, transform: DAffine3) -> Self {
        self.transform = transform;
        self
    }

    /// Convenience: set a pure translation transform.
    ///
    /// Equivalent to `with_transform(DAffine3::from_translation(t))`.
    pub fn with_translation(mut self, t: DVec3) -> Self {
        self.transform = DAffine3::from_translation(t);
        self
    }

    /// Set a color for this component.
    pub fn with_color(mut self, color: rcad_kernel::appearance::Color) -> Self {
        self.color = Some(color);
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Write
// ─────────────────────────────────────────────────────────────────────────────

/// Write a multi-component STEP assembly.
///
/// Each component is written as a separate `PRODUCT` + `PRODUCT_DEFINITION`
/// with its own geometry representation, linked into the root assembly via
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE`.
///
/// The [`AssemblyComponent::transform`] matrix is applied to vertex and geometry
/// coordinates before writing (baked into the STEP geometry, not stored as a
/// STEP transform entity).
pub fn write_assembly(assembly_name: &str, components: &[AssemblyComponent]) -> String {
    // Apply full DAffine3 transform into new BReps and collect colors.
    let prepared: Vec<(String, BRep, Option<rcad_kernel::appearance::Color>)> = components
        .iter()
        .map(|c| {
            let mut b = c.brep.clone();
            if c.transform != DAffine3::IDENTITY {
                b.apply_transform(c.transform);
            }
            (c.name.clone(), b, c.color)
        })
        .collect();

    write_assembly_step(assembly_name, &prepared)
}

fn push_record(records: &mut Vec<String>, next_id: &mut u64, body: String) -> u64 {
    let id = *next_id;
    *next_id += 1;
    records.push(format!("#{}={};", id, body));
    id
}

fn write_assembly_step(
    assembly_name: &str,
    components: &[(String, BRep, Option<rcad_kernel::appearance::Color>)],
) -> String {
    use std::fmt::Write as FmtWrite;

    let mut records: Vec<String> = Vec::new();
    let mut next_id: u64 = 1;

    macro_rules! push {
        ($body:expr) => {
            push_record(&mut records, &mut next_id, $body)
        };
    }

    // Shared context entities
    let app_ctx = push!("APPLICATION_CONTEXT('automotive_design')".to_string());
    push!(format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{})",
        app_ctx
    ));
    let prod_ctx = push!(format!(
        "PRODUCT_CONTEXT('part definition',#{},'mechanical')",
        app_ctx
    ));
    let def_ctx = push!(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{},'design')",
        app_ctx
    ));

    // Measurement context
    let len_unit = push!("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string());
    let rad_unit = push!("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
    let meas = push!(format!(
        "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})",
        rad_unit
    ));
    let dim_exp = push!("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)".to_string());
    let deg_unit = push!(format!(
        "( CONVERSION_BASED_UNIT('DEGREE',#{}) NAMED_UNIT(#{}) PLANE_ANGLE_UNIT() )",
        meas, dim_exp
    ));
    let sol_unit = push!("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string());
    let uncert = push!(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#{},'distance_accuracy_value','confusion accuracy')",
        len_unit
    ));
    let geom_ctx = push!(format!(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') )",
        uncert, len_unit, deg_unit, sol_unit
    ));

    // Assembly root product
    let asm_product = push!(format!(
        "PRODUCT('{}','{}','',( #{} ))",
        assembly_name, assembly_name, prod_ctx
    ));
    let asm_formation = push!(format!(
        "PRODUCT_DEFINITION_FORMATION('','',#{})",
        asm_product
    ));
    let asm_definition = push!(format!(
        "PRODUCT_DEFINITION('','',#{},#{})",
        asm_formation, def_ctx
    ));
    let asm_shape = push!(format!(
        "PRODUCT_DEFINITION_SHAPE('','',#{})",
        asm_definition
    ));
    let asm_rep = push!(format!(
        "SHAPE_REPRESENTATION('{}',(#{}),#{})",
        assembly_name, geom_ctx, geom_ctx
    ));
    push!(format!(
        "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
        asm_shape, asm_rep
    ));

    // Component products + NAUO links
    for (i, (comp_name, brep, comp_color)) in components.iter().enumerate() {
        // Write the component geometry as a standalone STEP string,
        // then inline its DATA section records with offset IDs.
        let colors = comp_color.map(|c| StepColor::new().with_solid_color(c));
        let comp_step = if let Some(sc) = &colors {
            StepWriter::write_string_colored(brep, sc)
        } else {
            StepWriter::write_string(
                brep,
                ExportSelection {
                    selected_faces: &[],
                    selected_edges: &[],
                },
            )
        };

        // Extract DATA records from component STEP string and re-number them
        let comp_records = extract_data_records(&comp_step);
        let id_offset = next_id - 1;

        // Re-number and collect component records
        let renumbered: Vec<String> = comp_records
            .iter()
            .map(|(orig_id, body)| {
                let new_id = orig_id + id_offset;
                let renumbered_body = renumber_refs(body, id_offset);
                format!("#{}={};", new_id, renumbered_body)
            })
            .collect();

        // Find the highest ID in the component (= shape_representation or similar)
        let comp_max_id = comp_records.iter().map(|(id, _)| *id).max().unwrap_or(1);
        let shape_rep_id = comp_max_id + id_offset; // last record is typically the shape repr

        // Advance next_id past all component records
        for record in &renumbered {
            records.push(record.clone());
        }
        next_id = shape_rep_id + 1;

        // Component product/definition
        let comp_product = push!(format!(
            "PRODUCT('{}','{}','',( #{} ))",
            comp_name, comp_name, prod_ctx
        ));
        let comp_formation = push!(format!(
            "PRODUCT_DEFINITION_FORMATION('','',#{})",
            comp_product
        ));
        let comp_definition = push!(format!(
            "PRODUCT_DEFINITION('','',#{},#{})",
            comp_formation, def_ctx
        ));
        let comp_pds = push!(format!(
            "PRODUCT_DEFINITION_SHAPE('','',#{})",
            comp_definition
        ));
        push!(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
            comp_pds, shape_rep_id
        ));

        // NEXT_ASSEMBLY_USAGE_OCCURRENCE: link component to assembly
        let nauo = push!(format!(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','',#{},#{},$)",
            i + 1,
            comp_name,
            asm_definition,
            comp_definition
        ));
        // PRODUCT_DEFINITION_SHAPE for the occurrence
        push!(format!(
            "PRODUCT_DEFINITION_SHAPE('Acme','occurrence shape',#{})",
            nauo
        ));
    }

    // Build output
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    let _ = writeln!(
        out,
        "FILE_DESCRIPTION(('RCAD assembly: {}'),'2;1');",
        assembly_name
    );
    out.push_str(
        "FILE_NAME('rcad_assembly.step','2026-04-11T00:00:00',(''),(''),'RCAD','RCAD','');\n",
    );
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");
    for record in &records {
        out.push_str(record);
        out.push('\n');
    }
    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Read
// ─────────────────────────────────────────────────────────────────────────────

/// Parse a STEP string that may contain a multi-component assembly and return
/// a flat list of [`AssemblyComponent`]s.
///
/// # Algorithm
///
/// 1. Scan the DATA section for `NEXT_ASSEMBLY_USAGE_OCCURRENCE` (NAUO) entities
///    to build a parent→child component map.
/// 2. For each NAUO child `PRODUCT_DEFINITION`, trace the chain:
///    `PRODUCT_DEFINITION` → `PRODUCT_DEFINITION_SHAPE` →
///    `SHAPE_DEFINITION_REPRESENTATION` → `SHAPE_REPRESENTATION`.
/// 3. BFS-collect all entity IDs reachable from that shape representation
///    (geometry, surfaces, curves, vertices, units, context, …).
/// 4. Build a self-contained sub-STEP string for each component and parse it
///    with [`crate::StepReader`] to obtain an isolated [`BRep`].
///
/// ## Fallback
///
/// If the file has no NAUO links (plain single-part STEP), the function returns
/// a single-element list containing the full parsed BRep with `transform =
/// IDENTITY`.
///
/// If a component's shape representation chain cannot be resolved (e.g. a
/// third-party file with non-standard structure), that component falls back to
/// the full merged BRep.
pub fn read_assembly(step: &str) -> Result<Vec<AssemblyComponent>, StepError> {
    let entity_map = parse_entity_map(step);
    let reverse_map = build_reverse_map(&entity_map);

    // Collect NAUO children: child PRODUCT_DEFINITION id + relation name
    let mut nauo_children: Vec<(u64, String)> = Vec::new();
    let mut has_nauo = false;

    for (_id, body) in &entity_map {
        if let Some(rest) = strip_entity_name(body, "NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            has_nauo = true;
            if let Some(args) = parse_args(rest) {
                let child_pd = parse_ref(args.get(4).copied().unwrap_or(""));
                let relation_name = unquote(args.get(1).copied().unwrap_or(""));
                if child_pd > 0 {
                    nauo_children.push((child_pd, relation_name));
                }
            }
        }
    }

    if !has_nauo || nauo_children.is_empty() {
        let brep = crate::StepReader::parse_string(step)?;
        let name = find_root_product_name(&entity_map).unwrap_or_else(|| "part".to_string());
        return Ok(vec![AssemblyComponent::new(name, brep)]);
    }

    // Parse the whole STEP once as a fallback for components whose geometry
    // cannot be isolated (e.g. third-party files with unusual structure).
    let merged_brep = crate::StepReader::parse_string(step)?;

    let mut components = Vec::new();
    for (pd_id, nauo_name) in &nauo_children {
        let name = resolve_product_name(*pd_id, &entity_map)
            .unwrap_or_else(|| nauo_name.clone());

        let brep = if let Some(sr_id) =
            find_shape_rep_for_pd(*pd_id, &entity_map, &reverse_map)
        {
            // BFS-collect all entities reachable from this component's shape
            // representation, then build a self-contained sub-STEP string.
            let reachable = collect_reachable(sr_id, &entity_map);
            let comp_step = build_component_step(&entity_map, &reachable);
            crate::StepReader::parse_string(&comp_step)
                .unwrap_or_else(|_| merged_brep.clone())
        } else {
            merged_brep.clone()
        };

        components.push(AssemblyComponent::new(name, brep));
    }

    Ok(components)
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal parsing helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Extract all `#N` reference IDs from a STEP entity body string.
fn extract_refs(body: &str) -> Vec<u64> {
    let mut refs = Vec::new();
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if i > start {
                if let Ok(n) = body[start..i].parse::<u64>() {
                    refs.push(n);
                }
            }
        } else {
            i += 1;
        }
    }
    refs
}

/// Build a reverse-reference map: referenced_id → list of entity IDs that
/// contain a `#referenced_id` in their body.
fn build_reverse_map(map: &HashMap<u64, String>) -> HashMap<u64, Vec<u64>> {
    let mut reverse: HashMap<u64, Vec<u64>> = HashMap::new();
    for (&id, body) in map {
        for ref_id in extract_refs(body) {
            reverse.entry(ref_id).or_default().push(id);
        }
    }
    reverse
}

/// BFS from `start_id`, following all `#N` references, and return the set of
/// all reachable entity IDs (including `start_id` itself).
fn collect_reachable(start_id: u64, map: &HashMap<u64, String>) -> HashSet<u64> {
    let mut visited: HashSet<u64> = HashSet::new();
    let mut queue = vec![start_id];
    while let Some(id) = queue.pop() {
        if !visited.insert(id) {
            continue;
        }
        if let Some(body) = map.get(&id) {
            for ref_id in extract_refs(body) {
                if !visited.contains(&ref_id) {
                    queue.push(ref_id);
                }
            }
        }
    }
    visited
}

/// Trace `PRODUCT_DEFINITION` → `PRODUCT_DEFINITION_SHAPE` →
/// `SHAPE_DEFINITION_REPRESENTATION` → shape representation ID.
///
/// Returns the shape representation entity ID, or `None` if the chain cannot
/// be resolved (e.g. the file uses a non-standard structure).
fn find_shape_rep_for_pd(
    pd_id: u64,
    map: &HashMap<u64, String>,
    reverse_map: &HashMap<u64, Vec<u64>>,
) -> Option<u64> {
    // Find PRODUCT_DEFINITION_SHAPE entities that reference pd_id
    let referencing = reverse_map.get(&pd_id)?;
    for &pds_id in referencing {
        let body = map.get(&pds_id)?;
        if !body.starts_with("PRODUCT_DEFINITION_SHAPE(") {
            continue;
        }
        // Find SHAPE_DEFINITION_REPRESENTATION entities that reference pds_id
        if let Some(referencing_pds) = reverse_map.get(&pds_id) {
            for &sdr_id in referencing_pds {
                if let Some(sdr_body) = map.get(&sdr_id) {
                    if let Some(args_str) =
                        sdr_body.strip_prefix("SHAPE_DEFINITION_REPRESENTATION(")
                    {
                        if let Some(args) = parse_args(args_str) {
                            let sr_id = parse_ref(args.get(1).copied().unwrap_or(""));
                            if sr_id > 0 {
                                return Some(sr_id);
                            }
                        }
                    }
                }
            }
        }
    }
    None
}

/// Build a self-contained STEP string containing only the entities in
/// `reachable`, using their original IDs from `map`.
fn build_component_step(map: &HashMap<u64, String>, reachable: &HashSet<u64>) -> String {
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    out.push_str("FILE_DESCRIPTION(('RCAD component'),'2;1');\n");
    out.push_str(
        "FILE_NAME('component.step','2026-04-11T00:00:00',(''),(''),'RCAD','RCAD','');\n",
    );
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");

    let mut ids: Vec<u64> = reachable.iter().copied().collect();
    ids.sort_unstable();
    for id in ids {
        if let Some(body) = map.get(&id) {
            out.push_str(&format!("#{}={};\n", id, body));
        }
    }

    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    out
}

/// Build a map from entity ID → entity body (the part after `#id=` and before `;`).
fn parse_entity_map(step: &str) -> std::collections::HashMap<u64, String> {
    let mut map = std::collections::HashMap::new();
    let mut in_data = false;
    for line in step.lines() {
        let line = line.trim();
        if line == "DATA;" {
            in_data = true;
            continue;
        }
        if line == "ENDSEC;" {
            in_data = false;
            continue;
        }
        if !in_data {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('#') {
            if let Some(eq) = stripped.find('=') {
                let id_str = &stripped[..eq];
                let body = stripped[eq + 1..].trim_end_matches(';');
                if let Ok(id) = id_str.parse::<u64>() {
                    map.insert(id, body.to_string());
                }
            }
        }
    }
    map
}

/// If `body` starts with `ENTITY_NAME(`, return the argument string (after the
/// opening paren, before the matching closing paren).
fn strip_entity_name<'a>(body: &'a str, entity: &str) -> Option<&'a str> {
    let prefix = entity;
    if body.starts_with(prefix) {
        body[prefix.len()..].strip_prefix('(')
    } else if let Some(inner) = body.strip_prefix('(') {
        // compound entity: ( ENTITY_NAME() ... )
        // Find the sub-entity by scanning for the name inside compound parens
        let inner = inner.trim();
        if inner.starts_with(entity) {
            inner[entity.len()..].strip_prefix('(')
        } else {
            None
        }
    } else {
        None
    }
}

/// Split a STEP argument list (without outer parens) on `,`, respecting nested
/// parens and string literals.
fn parse_args(args_str: &str) -> Option<Vec<&str>> {
    let mut result = Vec::new();
    let mut depth = 0i32;
    let mut in_str = false;
    let mut start = 0;
    let bytes = args_str.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'\'' if !in_str => in_str = true,
            b'\'' if in_str => in_str = false,
            b'(' if !in_str => depth += 1,
            b')' if !in_str => {
                if depth == 0 {
                    // closing paren of the whole arg list
                    result.push(args_str[start..i].trim());
                    return Some(result);
                }
                depth -= 1;
            }
            b',' if !in_str && depth == 0 => {
                result.push(args_str[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    if start < args_str.len() {
        result.push(args_str[start..].trim());
    }
    Some(result)
}

/// Parse `#N` → N, returning 0 on failure.
fn parse_ref(s: &str) -> u64 {
    s.trim().strip_prefix('#').and_then(|n| n.trim_end_matches(|c: char| !c.is_ascii_digit()).parse().ok()).unwrap_or(0)
}

/// Strip surrounding single-quotes from a STEP string literal.
fn unquote(s: &str) -> String {
    let s = s.trim();
    if s.starts_with('\'') && s.ends_with('\'') && s.len() >= 2 {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

/// Resolve PRODUCT_DEFINITION id → product name via PD → PD_FORMATION → PRODUCT.
fn resolve_product_name(
    pd_id: u64,
    map: &std::collections::HashMap<u64, String>,
) -> Option<String> {
    // PD body: PRODUCT_DEFINITION('','',#formation,#ctx)
    let pd_body = map.get(&pd_id)?;
    let pd_args = parse_args(pd_body.strip_prefix("PRODUCT_DEFINITION(")?.strip_suffix(')')?.as_ref())?;
    // Actually strip_entity_name handles compound; simpler approach:
    let formation_id = pd_args.get(2).map(|s| parse_ref(s))?;
    if formation_id == 0 {
        return None;
    }

    // Formation body: PRODUCT_DEFINITION_FORMATION('','',#product)
    let form_body = map.get(&formation_id)?;
    let form_args = parse_args(
        form_body
            .strip_prefix("PRODUCT_DEFINITION_FORMATION(")?.strip_suffix(')')?,
    )?;
    let prod_id = form_args.get(2).map(|s| parse_ref(s))?;
    if prod_id == 0 {
        return None;
    }

    // Product body: PRODUCT('id','name','desc',(#ctx))
    let prod_body = map.get(&prod_id)?;
    let prod_args = parse_args(prod_body.strip_prefix("PRODUCT(")?.strip_suffix(')')?.as_ref())?;
    // name is the second field (index 1)
    prod_args.get(1).map(|s| unquote(s))
}

/// Find the top-level PRODUCT name in a single-part STEP file (first PRODUCT entity).
fn find_root_product_name(map: &std::collections::HashMap<u64, String>) -> Option<String> {
    // Return the first PRODUCT entity's second field (name).
    let mut ids: Vec<u64> = map.keys().copied().collect();
    ids.sort();
    for id in ids {
        let body = &map[&id];
        if body.starts_with("PRODUCT(") {
            if let Some(args) = parse_args(body.strip_prefix("PRODUCT(")?.strip_suffix(')')?) {
                return args.get(1).map(|s| unquote(s));
            }
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Low-level helpers (retained from previous implementation)
// ─────────────────────────────────────────────────────────────────────────────

/// Extract (id, body) pairs from a STEP DATA section string.
fn extract_data_records(step: &str) -> Vec<(u64, String)> {
    let mut in_data = false;
    let mut result = Vec::new();
    for line in step.lines() {
        let line = line.trim();
        if line == "DATA;" {
            in_data = true;
            continue;
        }
        if line == "ENDSEC;" {
            in_data = false;
            continue;
        }
        if !in_data {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('#')
            && let Some(eq) = stripped.find('=')
        {
            let id_str = &stripped[..eq];
            let body = stripped[eq + 1..].trim_end_matches(';');
            if let Ok(id) = id_str.parse::<u64>() {
                result.push((id, body.to_string()));
            }
        }
    }
    result
}

/// Replace all `#N` references in a STEP entity body with `#(N + offset)`.
fn renumber_refs(body: &str, offset: u64) -> String {
    let mut out = String::with_capacity(body.len() + 16);
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if i > start {
                let num: u64 = body[start..i].parse().unwrap_or(0);
                out.push('#');
                out.push_str(&(num + offset).to_string());
            } else {
                out.push('#');
            }
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out
}
