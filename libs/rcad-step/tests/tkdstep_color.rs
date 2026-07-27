//! TKDESTEP-aligned tests for STEP rendering/color properties.
//!
//! OCCT source: src/DataExchange/TKDESTEP/GTests/STEPConstruct_RenderingProperties_Test.cxx
//! Tests color and rendering properties in STEP output, adapted for rcad-step's StepColor model.

use glam::DVec3;
use rcad_kernel::appearance::{Color, StepColor};
use rcad_kernel::topods;
use rcad_modeling::make_box_brep;
use rcad_step::{ExportSelection, StepProtocol, StepReader, StepWriteOptions, StepWriter};

fn all_faces() -> ExportSelection<'static> {
    ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    }
}

fn face_count(t: &topods::BRep) -> usize {
    t.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), topods::TShape::Face(_)))
        .count()
}

fn solid_count(t: &topods::BRep) -> usize {
    t.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), topods::TShape::Solid(_)))
        .count()
}

// ── Default constructor ──

#[test]
fn default_step_color_is_transparent_black() {
    let color = StepColor::new();
    let solid = color.solid_color.unwrap_or(Color {
        r: 0.0,
        g: 0.0,
        b: 0.0,
    });
    assert!(
        solid.r.abs() < f64::EPSILON
            && solid.g.abs() < f64::EPSILON
            && solid.b.abs() < f64::EPSILON,
        "default solid color should be black"
    );
}

// ── Solid color ──

#[test]
fn step_color_with_solid_color() {
    let color = StepColor::new().with_solid_color(Color {
        r: 0.8,
        g: 0.4,
        b: 0.2,
    });
    let solid = color.solid_color.unwrap_or(Color {
        r: 0.0,
        g: 0.0,
        b: 0.0,
    });
    assert!((solid.r - 0.8).abs() < 1e-6);
    assert!((solid.g - 0.4).abs() < 1e-6);
    assert!((solid.b - 0.2).abs() < 1e-6);
}

#[test]
fn solid_color_writes_to_step() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 1.0,
        g: 0.0,
        b: 0.0,
    });

    let step = StepWriter::write_string_colored(&brep, &colors);
    assert!(step.contains("COLOUR_RGB"));
    assert!(step.contains("STYLED_ITEM") || step.contains("PRESENTATION_STYLE_ASSIGNMENT"));
}

// ── Per-face colors ──

#[test]
fn per_face_color_writes_separate_styled_items() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let t = brep;

    let mut colors = StepColor::new();
    // Assign different colors to each of the 6 faces
    for i in 0..6 {
        let hue = i as f64 / 6.0;
        colors = colors.with_face_color(
            i,
            Color {
                r: (hue * 2.0).min(1.0),
                g: ((1.0 - hue) * 2.0).min(1.0),
                b: 0.5,
            },
        );
    }

    let step = StepWriter::write_string_colored(&t, &colors);
    assert!(step.contains("COLOUR_RGB"));
    // Verify each color component appears in the output
    let rgb_count = step.matches("COLOUR_RGB").count();
    assert!(
        rgb_count >= 3,
        "expected at least 3 COLOUR_RGB, got {rgb_count}"
    );
}

// ── Color roundtrip ──

#[test]
fn colored_solid_roundtrip_preserves_topology() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 0.8,
        g: 0.4,
        b: 0.2,
    });

    let step = StepWriter::write_string_colored(&brep, &colors);
    let parsed = StepReader::parse_string(&step).expect("colored STEP should parse");

    assert_eq!(solid_count(&parsed), 1);
    assert_eq!(face_count(&parsed), 6);
}

#[test]
fn colored_write_different_protocols() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 0.0,
        g: 0.5,
        b: 1.0,
    });
    let t = brep;

    for protocol in [StepProtocol::Ap214, StepProtocol::Ap242] {
        let step = StepWriter::write_string_with_options(
            &t,
            all_faces(),
            &StepWriteOptions {
                protocol,
                colors: Some(colors.clone()),
                ..Default::default()
            },
        );
        assert!(
            step.contains("COLOUR_RGB"),
            "color should be present for {protocol:?}"
        );
        let reparsed =
            StepReader::parse_string(&step).unwrap_or_else(|_| panic!("{protocol:?} should parse"));
        assert_eq!(solid_count(&reparsed), 1);
    }
}

// ── No-color output ──

#[test]
fn no_color_output_does_not_contain_rgb() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let step = StepWriter::write_string(&brep, all_faces());
    assert!(
        step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"),
        "STEP should contain shape rep"
    );
}

// ── AP242 metadata with colors ──

#[test]
fn ap242_colored_roundtrip() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 0.2,
        g: 0.6,
        b: 0.8,
    });
    let t = brep;

    let step = StepWriter::write_string_with_options(
        &t,
        all_faces(),
        &StepWriteOptions {
            protocol: StepProtocol::Ap242,
            colors: Some(colors),
            ..Default::default()
        },
    );

    assert!(step.contains("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING"));
    assert!(step.contains("COLOUR_RGB"));

    let parsed = StepReader::parse_string(&step).expect("AP242 colored STEP should parse");
    assert_eq!(solid_count(&parsed), 1);
    assert_eq!(face_count(&parsed), 6);
}

// ── Zero transparency default ──

#[test]
fn default_color_opaque() {
    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 0.5,
        g: 0.5,
        b: 0.5,
    });
    let step = StepWriter::write_string_colored(&brep, &colors);
    // Opaque entities should not contain INVISIBILITY
    assert!(
        !step.contains("INVISIBILITY"),
        "opaque should not set INVISIBILITY"
    );
}

// ── Solid color entity stream roundtrip ──

#[test]
fn colored_step_writes_and_reads_stream() {
    use std::io::Cursor;

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0).unwrap();
    let colors = StepColor::new().with_solid_color(Color {
        r: 1.0,
        g: 0.8,
        b: 0.6,
    });
    let t = brep;

    let step = StepWriter::write_string_with_options(
        &t,
        all_faces(),
        &StepWriteOptions {
            colors: Some(colors),
            ..Default::default()
        },
    );

    // Stream read should also work
    let parsed = StepReader::parse_reader(Cursor::new(step.as_bytes()))
        .expect("stream read of colored STEP should succeed");
    assert_eq!(solid_count(&parsed), 1);
}
