//! STEP assembly writer.
//!
//! Writes a multi-component assembly where each component is a separate BRep
//! with an optional transform and name.  Produces a STEP file with a full
//! NEXT_ASSEMBLY_USAGE_OCCURRENCE hierarchy so that importers (FreeCAD, OCCT,
//! etc.) can reconstruct the tree.
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

use glam::DVec3;
use rcad_kernel::appearance::StepColor;
use rcad_kernel::BRep;

use crate::writer::{ExportSelection, StepWriter};

/// A single component in an assembly.
#[derive(Clone)]
pub struct AssemblyComponent {
    /// Human-readable part name.
    pub name: String,
    /// The geometry.
    pub brep: BRep,
    /// Optional translation offset applied to the component.
    pub translation: DVec3,
    /// Optional RGB color for this component's faces.
    pub color: Option<rcad_kernel::appearance::Color>,
}

impl AssemblyComponent {
    pub fn new(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            brep,
            translation: DVec3::ZERO,
            color: None,
        }
    }

    pub fn with_translation(mut self, t: DVec3) -> Self {
        self.translation = t;
        self
    }

    pub fn with_color(mut self, color: rcad_kernel::appearance::Color) -> Self {
        self.color = Some(color);
        self
    }
}

/// Write a multi-component STEP assembly.
///
/// Each component is written as a separate `PRODUCT` + `PRODUCT_DEFINITION`
/// with its own geometry representation, linked into the root assembly via
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE`.
///
/// Translations are applied to vertex coordinates directly (no STEP transform
/// entity), which maximises compatibility with readers that don't support
/// `ITEM_DEFINED_TRANSFORMATION`.
pub fn write_assembly(assembly_name: &str, components: &[AssemblyComponent]) -> String {
    // Apply translations into new BReps and collect colors
    let prepared: Vec<(String, BRep, Option<rcad_kernel::appearance::Color>)> = components
        .iter()
        .map(|c| {
            let mut b = c.brep.clone();
            if c.translation != DVec3::ZERO {
                for v in &mut b.vertices {
                    v.point += c.translation;
                }
            }
            (c.name.clone(), b, c.color)
        })
        .collect();

    // Write each component as a standalone STEP string, then merge into one file
    // with a shared header and a proper assembly hierarchy.
    // For simplicity, write each component separately and collect DATA sections,
    // then re-emit as a single assembly file with NAUO links.
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
        "PRODUCT_CONTEXT('part definition',#{},'mechanical')", app_ctx
    ));
    let def_ctx = push!(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{},'design')", app_ctx
    ));

    // Measurement context
    let len_unit = push!("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string());
    let rad_unit = push!("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
    let meas = push!(format!(
        "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})", rad_unit
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
        "PRODUCT('{}','{}','',( #{} ))", assembly_name, assembly_name, prod_ctx
    ));
    let asm_formation = push!(format!(
        "PRODUCT_DEFINITION_FORMATION('','',#{})", asm_product
    ));
    let asm_definition = push!(format!(
        "PRODUCT_DEFINITION('','',#{},#{})", asm_formation, def_ctx
    ));
    let asm_shape = push!(format!(
        "PRODUCT_DEFINITION_SHAPE('','',#{})", asm_definition
    ));
    let asm_rep = push!(format!(
        "SHAPE_REPRESENTATION('{}',(#{}),#{})", assembly_name, geom_ctx, geom_ctx
    ));
    push!(format!("SHAPE_DEFINITION_REPRESENTATION(#{},#{})", asm_shape, asm_rep));

    // Component products + NAUO links
    for (i, (comp_name, brep, comp_color)) in components.iter().enumerate() {
        // Write the component geometry as a standalone STEP string,
        // then inline its DATA section records with offset IDs.
        let colors = comp_color.map(|c| {
            StepColor::new().with_solid_color(c)
        });
        let comp_step = if let Some(sc) = &colors {
            StepWriter::write_string_colored(brep, sc)
        } else {
            StepWriter::write_string(brep, ExportSelection {
                selected_faces: &[], selected_edges: &[],
            })
        };

        // Extract DATA records from component STEP string and re-number them
        let comp_records = extract_data_records(&comp_step);
        let id_offset = next_id - 1;

        // Re-number and collect component records
        let renumbered: Vec<String> = comp_records.iter().map(|(orig_id, body)| {
            let new_id = orig_id + id_offset;
            let renumbered_body = renumber_refs(body, id_offset);
            format!("#{}={};", new_id, renumbered_body)
        }).collect();

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
            "PRODUCT('{}','{}','',( #{} ))", comp_name, comp_name, prod_ctx
        ));
        let comp_formation = push!(format!(
            "PRODUCT_DEFINITION_FORMATION('','',#{})", comp_product
        ));
        let comp_definition = push!(format!(
            "PRODUCT_DEFINITION('','',#{},#{})", comp_formation, def_ctx
        ));
        let comp_pds = push!(format!(
            "PRODUCT_DEFINITION_SHAPE('','',#{})", comp_definition
        ));
        push!(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})", comp_pds, shape_rep_id
        ));

        // NEXT_ASSEMBLY_USAGE_OCCURRENCE: link component to assembly
        let nauo = push!(format!(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','',#{},#{},$)",
            i + 1, comp_name, asm_definition, comp_definition
        ));
        // PRODUCT_DEFINITION_SHAPE for the occurrence
        push!(format!(
            "PRODUCT_DEFINITION_SHAPE('Acme','occurrence shape',#{})", nauo
        ));
    }

    // Build output
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    let _ = writeln!(out, "FILE_DESCRIPTION(('RCAD assembly: {}'),'2;1');", assembly_name);
    out.push_str("FILE_NAME('rcad_assembly.step','2026-04-03T00:00:00',(''),(''),'RCAD','RCAD','');\n");
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

/// Extract (id, body) pairs from a STEP DATA section string.
fn extract_data_records(step: &str) -> Vec<(u64, String)> {
    let mut in_data = false;
    let mut result = Vec::new();
    for line in step.lines() {
        let line = line.trim();
        if line == "DATA;" { in_data = true; continue; }
        if line == "ENDSEC;" { in_data = false; continue; }
        if !in_data { continue; }
        // Parse #id=body;
        if let Some(stripped) = line.strip_prefix('#') {
            if let Some(eq) = stripped.find('=') {
                let id_str = &stripped[..eq];
                let body = stripped[eq+1..].trim_end_matches(';');
                if let Ok(id) = id_str.parse::<u64>() {
                    result.push((id, body.to_string()));
                }
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
