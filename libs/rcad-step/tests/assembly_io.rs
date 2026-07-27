//! Integration tests for rcad-step assembly read/write.

use glam::{DAffine3, DVec3};
use rcad_algorithms::{HealingMode, HealingOptions};
use rcad_kernel::topods;
use rcad_modeling::make_box_brep;
use rcad_step::{
    AssemblyComponent, AssemblyNode, read_assembly, read_assembly_tree,
    read_assembly_tree_with_healing, read_assembly_with_healing,
    read_assembly_with_healing_report_json, write_assembly, write_assembly_tree,
};

fn make_box(origin: DVec3) -> rcad_kernel::topods::BRep {
    let mut b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    if origin != DVec3::ZERO {
        b.apply_transform(DAffine3::from_translation(origin));
    }
    b.to_topods()
}

/// Write an assembly with two named components, parse it back, and verify:
/// - component count == 2
/// - component names match
#[test]
fn write_read_assembly_component_count() {
    let comp_a = AssemblyComponent::new("box_a", make_box(DVec3::ZERO));
    let comp_b = AssemblyComponent::new("box_b", make_box(DVec3::new(5.0, 0.0, 0.0)));

    let step = write_assembly("test_asm", &[comp_a, comp_b]);

    // Assembly NAUO (2) plus optional per-part placement NAUOs from `StepWriter`.
    let nauo_count = step
        .lines()
        .filter(|l| l.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
        .count();
    assert!(
        nauo_count >= 2,
        "expected at least 2 assembly NAUO lines, got {}",
        nauo_count
    );

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

    let comp =
        AssemblyComponent::new("shifted_box", base_box.to_topods()).with_translation(translation);

    let step = write_assembly("shift_test", &[comp]);

    // After baking the transform in write_assembly, all vertices should be at x >= 10.
    let components = read_assembly(&step).expect("read_assembly");
    assert!(!components.is_empty());

    // The merged BRep returned by read_assembly contains baked geometry.
    let brep = &components[0].brep;
    for ts in &brep.tshapes {
        if let topods::TShape::Vertex(vd) = &**ts {
            assert!(
                vd.point.x >= 9.999,
                "vertex x should be >= 10 after baking translation, got {}",
                vd.point.x
            );
        }
    }
}

/// A plain single-part STEP file (no NAUO) parsed via read_assembly should
/// return exactly one component.
#[test]
fn single_part_step_returns_one_component() {
    use rcad_step::{ExportSelection, StepWriter};

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    let step = StepWriter::write_string(
        &brep.to_topods(),
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

#[test]
fn assembly_with_rotation_no_panic() {
    use std::f64::consts::FRAC_PI_4;

    let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let rotation = DAffine3::from_rotation_y(FRAC_PI_4);

    let comp = AssemblyComponent::new("rotated_box", brep.to_topods()).with_transform(rotation);
    let step = write_assembly("rotation_test", &[comp]);

    // Should not panic and should produce a valid STEP structure.
    assert!(step.contains("ISO-10303-21"));
    assert!(step.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"));

    let components = read_assembly(&step).expect("read_assembly after rotation");
    assert_eq!(components.len(), 1);
}

#[test]
fn read_assembly_with_healing_reports_per_component() {
    let comp_a = AssemblyComponent::new("box_a", make_box(DVec3::ZERO));
    let comp_b = AssemblyComponent::new("box_b", make_box(DVec3::new(5.0, 0.0, 0.0)));

    let step = write_assembly("heal_asm", &[comp_a, comp_b]);
    let (components, reports) = read_assembly_with_healing(
        &step,
        HealingOptions {
            mode: HealingMode::AnalyzeOnly,
            ..HealingOptions::default()
        },
    )
    .expect("read_assembly_with_healing failed");

    assert_eq!(components.len(), 2);
    assert_eq!(reports.len(), 2);
    assert!(
        reports
            .iter()
            .all(|r| r.initial.is_valid() && r.final_result.is_valid())
    );
}

#[test]
fn read_assembly_with_healing_report_json_has_stable_schema() {
    let comp_a = AssemblyComponent::new("box_a", make_box(DVec3::ZERO));
    let comp_b = AssemblyComponent::new("box_b", make_box(DVec3::new(5.0, 0.0, 0.0)));

    let step = write_assembly("heal_json_asm", &[comp_a, comp_b]);
    let (components, reports, json) = read_assembly_with_healing_report_json(
        &step,
        HealingOptions {
            mode: HealingMode::AnalyzeOnly,
            ..HealingOptions::default()
        },
    )
    .expect("read_assembly_with_healing_report_json failed");

    assert_eq!(components.len(), 2);
    assert_eq!(reports.len(), 2);

    let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid report json");
    assert_eq!(parsed["schema"], "step.assembly.import.healing.v1");
    assert_eq!(parsed["component_count"], 2);
    assert_eq!(parsed["clean_components"], 2);
    assert_eq!(parsed["failed_components"], 0);
}

// ─── nested assembly tree tests ───────────────────────────────────────────────

/// Write a two-level nested assembly tree and verify the STEP structure.
#[test]
fn nested_assembly_tree_write_has_nauo() {
    let leaf_a = AssemblyNode::leaf("part_a", make_box(DVec3::ZERO));
    let leaf_b = AssemblyNode::leaf("part_b", make_box(DVec3::new(2.0, 0.0, 0.0)));
    let sub = AssemblyNode::branch("sub_asm", vec![leaf_a, leaf_b]);
    let leaf_c = AssemblyNode::leaf("part_c", make_box(DVec3::new(5.0, 0.0, 0.0)));
    let root = AssemblyNode::branch("root_asm", vec![sub, leaf_c]);

    let step = write_assembly_tree("root_asm", &root);

    // Should have 3 NAUO entries: sub_asm→part_a, sub_asm→part_b, root→sub_asm, root→part_c
    // Actually 4: root→sub_asm, root→part_c, sub_asm→part_a, sub_asm→part_b
    let nauo_count = step
        .lines()
        .filter(|l| l.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
        .count();
    assert!(
        nauo_count >= 3,
        "expected at least 3 NAUO entries, got {}",
        nauo_count
    );
    assert!(step.contains("part_a"), "should contain part_a");
    assert!(step.contains("part_b"), "should contain part_b");
    assert!(step.contains("part_c"), "should contain part_c");
    assert!(step.contains("sub_asm"), "should contain sub_asm");
}

/// Round-trip a nested assembly tree through STEP and verify the tree structure.
#[test]
fn nested_assembly_tree_round_trip() {
    let leaf_a = AssemblyNode::leaf("part_a", make_box(DVec3::ZERO));
    let leaf_b = AssemblyNode::leaf("part_b", make_box(DVec3::new(3.0, 0.0, 0.0)));
    let sub = AssemblyNode::branch("sub_asm", vec![leaf_a, leaf_b]);
    let root = AssemblyNode::branch("root_asm", vec![sub]);

    let step = write_assembly_tree("root_asm", &root);
    let tree = read_assembly_tree(&step).expect("read_assembly_tree failed");

    // Root should be a branch named "root_asm" with one child "sub_asm".
    assert_eq!(tree.name, "root_asm");
    assert!(tree.brep.is_none(), "root should be a branch (no geometry)");
    assert_eq!(tree.children.len(), 1, "root should have 1 child");

    let sub_node = &tree.children[0];
    assert_eq!(sub_node.name, "sub_asm");
    assert!(sub_node.brep.is_none(), "sub_asm should be a branch");
    assert_eq!(sub_node.children.len(), 2, "sub_asm should have 2 children");

    let names: Vec<&str> = sub_node.children.iter().map(|c| c.name.as_str()).collect();
    assert!(names.contains(&"part_a"), "sub_asm should contain part_a");
    assert!(names.contains(&"part_b"), "sub_asm should contain part_b");
}

#[test]
fn read_assembly_tree_with_healing_reports_leaf_nodes() {
    let leaf_a = AssemblyNode::leaf("part_a", make_box(DVec3::ZERO));
    let leaf_b = AssemblyNode::leaf("part_b", make_box(DVec3::new(3.0, 0.0, 0.0)));
    let sub = AssemblyNode::branch("sub_asm", vec![leaf_a, leaf_b]);
    let leaf_c = AssemblyNode::leaf("part_c", make_box(DVec3::new(7.0, 0.0, 0.0)));
    let root = AssemblyNode::branch("root_asm", vec![sub, leaf_c]);

    let step = write_assembly_tree("root_asm", &root);
    let (tree, reports) = read_assembly_tree_with_healing(
        &step,
        HealingOptions {
            mode: HealingMode::AnalyzeOnly,
            ..HealingOptions::default()
        },
    )
    .expect("read_assembly_tree_with_healing failed");

    assert_eq!(tree.name, "root_asm");
    assert_eq!(reports.len(), 3, "expected one report per leaf node");

    let report_names: Vec<&str> = reports.iter().map(|r| r.name.as_str()).collect();
    assert!(report_names.contains(&"part_a"));
    assert!(report_names.contains(&"part_b"));
    assert!(report_names.contains(&"part_c"));

    for report in &reports {
        assert!(
            report.report.initial.is_valid() && report.report.final_result.is_valid(),
            "healing report should remain valid for clean geometry"
        );
    }
}
