//! Integration tests for rcad-step assembly read/write.

use glam::{DAffine3, DVec3};
use rcad_modeling::make_box_brep;
use rcad_step::{AssemblyComponent, read_assembly, write_assembly};

fn make_box(origin: DVec3) -> rcad_kernel::BRep {
    let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    if origin != DVec3::ZERO {
        b.apply_transform(DAffine3::from_translation(origin));
    }
    b
}

/// Write an assembly with two named components, parse it back, and verify:
/// - component count == 2
/// - component names match
#[test]
fn write_read_assembly_component_count() {
    let comp_a = AssemblyComponent::new("box_a", make_box(DVec3::ZERO));
    let comp_b = AssemblyComponent::new("box_b", make_box(DVec3::new(5.0, 0.0, 0.0)));

    let step = write_assembly("test_asm", &[comp_a, comp_b]);

    // Basic structural check: NAUO should appear twice.
    let nauo_count = step
        .lines()
        .filter(|l| l.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
        .count();
    assert_eq!(nauo_count, 2, "expected 2 NAUO entries, got {}", nauo_count);

    let components = read_assembly(&step).expect("read_assembly failed");
    assert_eq!(
        components.len(),
        2,
        "expected 2 components, got {}",
        components.len()
    );

    let names: Vec<&str> = components.iter().map(|c| c.name.as_str()).collect();
    assert!(
        names.contains(&"box_a"),
        "expected 'box_a' in names: {:?}",
        names
    );
    assert!(
        names.contains(&"box_b"),
        "expected 'box_b' in names: {:?}",
        names
    );
}

/// Component with a translation transform: after write+read the geometry
/// (baked into vertices) should reflect the translated position.
#[test]
fn assembly_with_translation_baked_into_geometry() {
    let base_box = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let translation = DVec3::new(10.0, 0.0, 0.0);

    let comp = AssemblyComponent::new("shifted_box", base_box)
        .with_translation(translation);

    let step = write_assembly("shift_test", &[comp]);

    // After baking the transform in write_assembly, all vertices should be at x >= 10.
    let components = read_assembly(&step).expect("read_assembly");
    assert!(!components.is_empty());

    // The merged BRep returned by read_assembly contains baked geometry.
    let brep = &components[0].brep;
    for v in &brep.vertices {
        assert!(
            v.point.x >= 9.999,
            "vertex x should be >= 10 after baking translation, got {}",
            v.point.x
        );
    }
}

/// A plain single-part STEP file (no NAUO) parsed via read_assembly should
/// return exactly one component.
#[test]
fn single_part_step_returns_one_component() {
    use rcad_step::{ExportSelection, StepWriter};

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    let step = StepWriter::write_string(
        &brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );

    let components = read_assembly(&step).expect("read_assembly on single-part STEP");
    assert_eq!(
        components.len(),
        1,
        "single-part STEP should give 1 component, got {}",
        components.len()
    );
}

/// Full affine transform (rotation): write_assembly bakes DAffine3 rotation.
/// Verify the STEP string round-trip completes without panic.
#[test]
fn assembly_with_rotation_no_panic() {
    use std::f64::consts::FRAC_PI_4;

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let rotation = DAffine3::from_rotation_y(FRAC_PI_4);

    let comp = AssemblyComponent::new("rotated_box", brep).with_transform(rotation);
    let step = write_assembly("rotation_test", &[comp]);

    // Should not panic and should produce a valid STEP structure.
    assert!(step.contains("ISO-10303-21"));
    assert!(step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"));

    let components = read_assembly(&step).expect("read_assembly after rotation");
    assert_eq!(components.len(), 1);
}
