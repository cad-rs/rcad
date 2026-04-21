//! Example: Chamfer and fillet on BRep edges.
//!
//! Demonstrates:
//!   1. chamfer_edge: bevel a single edge of a box with a flat cut
//!   2. fillet_edge: round a single edge of a box with a cylindrical blend
//!   3. all_edges_chamfer: chamfer every edge of a box (applied sequentially)
//!
//! Run: cargo run --example phase_f_demo
//!
//! Output files (written to current directory):
//!   box_chamfer.step
//!   box_fillet.step
//!   box_all_chamfer.step

use glam::DVec3;
use rcad_modeling::{box_brep, chamfer_edge, fillet_edge};
use rcad_step::{ExportSelection, StepWriter};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn save(path: &str, content: &str) {
    std::fs::write(path, content).unwrap();
    println!("  → wrote {path}  ({} bytes)", content.len());
}

fn step_record_count(step: &str) -> usize {
    step.lines().filter(|l| l.starts_with('#')).count()
}

fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.solids
        .first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0)
}

// ── 1. Chamfer single edge ────────────────────────────────────────────────────

fn demo_chamfer() {
    println!("\n=== 1. Chamfer: edge 0 of 2×1.5×1 box (dist=0.2) ===");

    let base = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 1.5, 1.0).unwrap();
    println!("  Original faces: {}", face_count(&base));

    let result = chamfer_edge(&base, 0, 0.2).expect("chamfer_edge failed");
    println!(
        "  After chamfer: {} faces (expected 9 = 6 original + 1 chamfer + 2 closing)",
        face_count(&result)
    );

    let step = StepWriter::write_string(
        &result,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    save("box_chamfer.step", &step);
    println!("  STEP records: {}", step_record_count(&step));
}

// ── 2. Fillet single edge ─────────────────────────────────────────────────────

fn demo_fillet() {
    println!("\n=== 2. Fillet: edge 0 of 2×1.5×1 box (radius=0.2) ===");

    let base = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 1.5, 1.0).unwrap();
    println!("  Original faces: {}", face_count(&base));

    let result = fillet_edge(&base, 0, 0.2).expect("fillet_edge failed");
    println!(
        "  After fillet: {} faces (expected 9 = 6 original + 1 fillet + 2 closing)",
        face_count(&result)
    );

    let step = StepWriter::write_string(
        &result,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    save("box_fillet.step", &step);
    println!("  STEP records: {}", step_record_count(&step));
}

// ── 3. All edges chamfered ────────────────────────────────────────────────────

fn demo_all_edges_chamfer() {
    println!(
        "\n=== 3. All-edges chamfer: chamfer each original edge independently, export one ==="
    );

    let base = box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 1.5, 1.0).unwrap();
    let original_edge_count = base.edges.len();
    println!("  Original edge count: {original_edge_count}");

    // Chamfer each of the 12 box edges independently and collect results.
    let mut results = Vec::new();
    for ei in 0..original_edge_count {
        if let Ok(result) = chamfer_edge(&base, ei, 0.15) {
            results.push((ei, result));
        }
    }
    println!(
        "  Successfully chamfered {} of {original_edge_count} edges",
        results.len()
    );

    // Export the last successful result (all should be 9 faces each).
    if let Some((last_idx, last_result)) = results.last() {
        println!(
            "  Last chamfered edge {last_idx}: {} faces",
            face_count(last_result)
        );
        let step = StepWriter::write_string(
            last_result,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );
        save("box_all_chamfer.step", &step);
        println!("  STEP records: {}", step_record_count(&step));
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔═══════════════════════════════════════════════╗");
    println!("║            RCAD Chamfer / Fillet Demo        ║");
    println!("║      Chamfer · Fillet · All-edges Chamfer     ║");
    println!("╚═══════════════════════════════════════════════╝");

    demo_chamfer();
    demo_fillet();
    demo_all_edges_chamfer();

    println!("\n✓ Chamfer / fillet demo complete.");
    println!("  Output: box_chamfer.step  box_fillet.step  box_all_chamfer.step");
}
