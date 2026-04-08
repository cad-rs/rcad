//! Example: Section gallery — slice many shapes with cutting planes and export STEP.
//!
//! For each shape, takes 3 section cuts and reports the resulting polylines.
//! Also exports the original shape as a STEP file for visual reference.
//!
//! Run: cargo run --example section_gallery

use glam::DVec3;
use rcad_algorithms::section_polylines;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{BRep, PrimitiveSolid};
use rcad_modeling::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::{extrude, revolve};
use rcad_scene::append_brep;
use rcad_step::writer::{ExportSelection, StepWriter};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_3};

fn main() {
    // ── 1. Unit box — three axis-aligned sections ──────────────────────
    println!("=== 1. Unit Box Sections ===");
    let box1 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    write_step(&box1, "output_section_box.step");
    section_report(
        "x=0.5 cut",
        &box1,
        Plane {
            origin: DVec3::new(0.5, 0.0, 0.0),
            normal: DVec3::X,
        },
    );
    section_report(
        "y=0.5 cut",
        &box1,
        Plane {
            origin: DVec3::new(0.0, 0.5, 0.0),
            normal: DVec3::Y,
        },
    );
    section_report(
        "z=0.5 cut",
        &box1,
        Plane {
            origin: DVec3::new(0.0, 0.0, 0.5),
            normal: DVec3::Z,
        },
    );

    // ── 2. 2×3×4 box — diagonal cut ───────────────────────────────────
    println!("\n=== 2. 2×3×4 Box Diagonal Section ===");
    let box2 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 3.0,
        depth: 4.0,
    });
    write_step(&box2, "output_section_box234.step");
    // Diagonal plane cutting through x+z=3 (normal along XZ diagonal)
    let diag_normal = DVec3::new(1.0, 0.0, 1.0).normalize();
    section_report(
        "diagonal x+z=3",
        &box2,
        Plane {
            origin: DVec3::new(1.5, 0.0, 1.5),
            normal: diag_normal,
        },
    );
    // Horizontal mid-height
    section_report(
        "y=1.5 (mid height)",
        &box2,
        Plane {
            origin: DVec3::new(0.0, 1.5, 0.0),
            normal: DVec3::Y,
        },
    );

    // ── 3. Square pillar — section at mid-height ───────────────────────
    println!("\n=== 3. Square Pillar 1×1×5 Sections ===");
    let sq = polygon_profile(&unit_square());
    let pillar = extrude(&sq, 0, DVec3::Z, 5.0).unwrap();
    write_step(&pillar, "output_section_pillar.step");
    section_report(
        "z=2.5 (mid)",
        &pillar,
        Plane {
            origin: DVec3::new(0.0, 0.0, 2.5),
            normal: DVec3::Z,
        },
    );
    section_report(
        "x=0.5",
        &pillar,
        Plane {
            origin: DVec3::new(0.5, 0.0, 0.0),
            normal: DVec3::X,
        },
    );
    section_report(
        "45° through pillar",
        &pillar,
        Plane {
            origin: DVec3::new(0.5, 0.5, 2.5),
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        },
    );

    // ── 4. Equilateral triangular prism ───────────────────────────────
    println!("\n=== 4. Triangular Prism Sections ===");
    let tri = polygon_profile(&[
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(1.0, 1.732, 0.0),
    ]);
    let prism = extrude(&tri, 0, DVec3::Z, 4.0).unwrap();
    write_step(&prism, "output_section_prism.step");
    section_report(
        "z=2 (mid height)",
        &prism,
        Plane {
            origin: DVec3::new(0.0, 0.0, 2.0),
            normal: DVec3::Z,
        },
    );
    section_report(
        "x=1 (vertical)",
        &prism,
        Plane {
            origin: DVec3::new(1.0, 0.0, 0.0),
            normal: DVec3::X,
        },
    );
    section_report(
        "tilted through base edge",
        &prism,
        Plane {
            origin: DVec3::new(1.0, 0.5, 2.0),
            normal: DVec3::new(0.0, 1.0, 1.0).normalize(),
        },
    );

    // ── 5. Hexagonal bolt head ─────────────────────────────────────────
    println!("\n=== 5. Hexagonal Bolt Head Sections ===");
    let hex_pts: Vec<DVec3> = (0..6)
        .map(|i| {
            let a = FRAC_PI_3 * i as f64;
            DVec3::new(a.cos(), a.sin(), 0.0)
        })
        .collect();
    let hex = polygon_profile(&hex_pts);
    let bolt = extrude(&hex, 0, DVec3::Z, 0.6).unwrap();
    write_step(&bolt, "output_section_bolt.step");
    section_report(
        "z=0.3 (mid)",
        &bolt,
        Plane {
            origin: DVec3::new(0.0, 0.0, 0.3),
            normal: DVec3::Z,
        },
    );
    section_report(
        "x=0 (bisect)",
        &bolt,
        Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        },
    );

    // ── 6. L-beam cross section ────────────────────────────────────────
    println!("\n=== 6. L-Beam Sections ===");
    let l_pts = [
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(3.0, 0.0, 0.0),
        DVec3::new(3.0, 0.5, 0.0),
        DVec3::new(0.5, 0.5, 0.0),
        DVec3::new(0.5, 3.0, 0.0),
        DVec3::new(0.0, 3.0, 0.0),
    ];
    let l_prof = polygon_profile(&l_pts);
    let l_beam = extrude(&l_prof, 0, DVec3::Z, 6.0).unwrap();
    write_step(&l_beam, "output_section_l_beam.step");
    section_report(
        "z=3 (cross section)",
        &l_beam,
        Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
        },
    );
    section_report(
        "y=0.25 (through flange)",
        &l_beam,
        Plane {
            origin: DVec3::new(0.0, 0.25, 0.0),
            normal: DVec3::Y,
        },
    );
    section_report(
        "x=0.25 (through web)",
        &l_beam,
        Plane {
            origin: DVec3::new(0.25, 0.0, 0.0),
            normal: DVec3::X,
        },
    );

    // ── 7. Quarter-turn elbow (revolve 90°) ───────────────────────────
    println!("\n=== 7. Quarter Elbow Sections ===");
    let r_prof = polygon_profile(&[
        DVec3::new(3.0, -0.25, 0.0),
        DVec3::new(3.5, -0.25, 0.0),
        DVec3::new(3.5, 0.25, 0.0),
        DVec3::new(3.0, 0.25, 0.0),
    ]);
    let elbow = revolve(&r_prof, 0, DVec3::ZERO, DVec3::Z, FRAC_PI_2).unwrap();
    write_step(&elbow, "output_section_elbow.step");
    // Cut perpendicular to X at the elbow start (x=3.25)
    section_report(
        "x=3.25 (start face)",
        &elbow,
        Plane {
            origin: DVec3::new(3.25, 0.0, 0.0),
            normal: DVec3::X,
        },
    );
    // Cut perpendicular to Y at elbow end (y=3.25)
    section_report(
        "y=3.25 (end face)",
        &elbow,
        Plane {
            origin: DVec3::new(0.0, 3.25, 0.0),
            normal: DVec3::Y,
        },
    );
    // Oblique cut through the arc at 45°
    section_report(
        "45° bisect",
        &elbow,
        Plane {
            origin: DVec3::new(2.3, 2.3, 0.0),
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        },
    );

    // ── 8. Half-pipe (revolve 180°) — sections ───────────────────────
    println!("\n=== 8. Half-Pipe (Revolve 180°) Sections ===");
    let ring_prof = polygon_profile(&[
        DVec3::new(2.0, -0.3, 0.0),
        DVec3::new(2.6, -0.3, 0.0),
        DVec3::new(2.6, 0.3, 0.0),
        DVec3::new(2.0, 0.3, 0.0),
    ]);
    let half_pipe = revolve(&ring_prof, 0, DVec3::ZERO, DVec3::Y, std::f64::consts::PI).unwrap();
    write_step(&half_pipe, "output_section_halfpipe.step");
    // Planes that cut through the lateral faces
    section_report(
        "x=2.3 (outer wall)",
        &half_pipe,
        Plane {
            origin: DVec3::new(2.3, 0.0, 0.0),
            normal: DVec3::X,
        },
    );
    section_report(
        "z=-2.3 (at 90° rotation)",
        &half_pipe,
        Plane {
            origin: DVec3::new(0.0, 0.0, -2.3),
            normal: DVec3::Z,
        },
    );
    section_report(
        "y=0 (mid-height)",
        &half_pipe,
        Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        },
    );

    // ── 9. Stacked assembly — union of 3 extruded layers ──────────────
    println!("\n=== 9. Stacked Assembly Sections ===");
    let base_sq = polygon_profile(&[
        DVec3::new(-1.5, -1.5, 0.0),
        DVec3::new(1.5, -1.5, 0.0),
        DVec3::new(1.5, 1.5, 0.0),
        DVec3::new(-1.5, 1.5, 0.0),
    ]);
    let mut stack = extrude(&base_sq, 0, DVec3::Z, 0.5).unwrap();

    let mid_sq = polygon_profile(&[
        DVec3::new(-1.0, -1.0, 0.0),
        DVec3::new(1.0, -1.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(-1.0, 1.0, 0.0),
    ]);
    let mut mid = extrude(&mid_sq, 0, DVec3::Z, 0.5).unwrap();
    for v in &mut mid.vertices {
        v.point.z += 0.5;
    }

    let top_sq = polygon_profile(&[
        DVec3::new(-0.5, -0.5, 0.0),
        DVec3::new(0.5, -0.5, 0.0),
        DVec3::new(0.5, 0.5, 0.0),
        DVec3::new(-0.5, 0.5, 0.0),
    ]);
    let mut top = extrude(&top_sq, 0, DVec3::Z, 0.5).unwrap();
    for v in &mut top.vertices {
        v.point.z += 1.0;
    }

    append_brep(&mut stack, mid);
    append_brep(&mut stack, top);
    write_step(&stack, "output_section_stack.step");
    section_report(
        "z=0.25 (base level)",
        &stack,
        Plane {
            origin: DVec3::new(0.0, 0.0, 0.25),
            normal: DVec3::Z,
        },
    );
    section_report(
        "z=0.75 (mid level)",
        &stack,
        Plane {
            origin: DVec3::new(0.0, 0.0, 0.75),
            normal: DVec3::Z,
        },
    );
    section_report(
        "z=1.25 (top level)",
        &stack,
        Plane {
            origin: DVec3::new(0.0, 0.0, 1.25),
            normal: DVec3::Z,
        },
    );

    // ── 10. Plus / cross beam ─────────────────────────────────────────
    println!("\n=== 10. Plus Beam Sections ===");
    let th = 0.4_f64;
    let arm = 1.5_f64;
    let plus_pts = [
        DVec3::new(-th, -arm, 0.0),
        DVec3::new(th, -arm, 0.0),
        DVec3::new(th, -th, 0.0),
        DVec3::new(arm, -th, 0.0),
        DVec3::new(arm, th, 0.0),
        DVec3::new(th, th, 0.0),
        DVec3::new(th, arm, 0.0),
        DVec3::new(-th, arm, 0.0),
        DVec3::new(-th, th, 0.0),
        DVec3::new(-arm, th, 0.0),
        DVec3::new(-arm, -th, 0.0),
        DVec3::new(-th, -th, 0.0),
    ];
    let plus_prof = polygon_profile(&plus_pts);
    let plus_beam = extrude(&plus_prof, 0, DVec3::Z, 6.0).unwrap();
    write_step(&plus_beam, "output_section_plus.step");
    section_report(
        "z=3 (cross section)",
        &plus_beam,
        Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
        },
    );
    section_report(
        "y=0 (horizontal slice)",
        &plus_beam,
        Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        },
    );
    section_report(
        "45° oblique",
        &plus_beam,
        Plane {
            origin: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        },
    );

    println!("\nAll section STEP files exported.");
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn unit_square() -> Vec<DVec3> {
    vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ]
}

fn polygon_profile(pts: &[DVec3]) -> BRep {
    let n = pts.len();
    let mut brep = BRep::default();
    let vis: Vec<usize> = pts.iter().map(|&p| make_vertex(&mut brep, p)).collect();
    let mut wires = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = pts[i];
        let b = pts[j];
        let dir = (b - a).normalize_or_zero();
        let len = (b - a).length();
        let eidx = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            vis[i],
            vis[j],
        )
        .unwrap();
        wires.push(WireEdge::fwd(eidx));
    }
    let surface = Surface3::Plane(Plane {
        origin: pts[0],
        normal: DVec3::Z,
    });
    make_face(&mut brep, surface, make_wire(wires), vec![]).unwrap();
    brep
}

fn section_report(label: &str, brep: &BRep, plane: Plane) {
    let loops = section_polylines(brep, &plane);
    let total_pts: usize = loops.iter().map(|l| l.len()).sum();
    println!(
        "  {:<30}  {} loop(s), {} pts total",
        label,
        loops.len(),
        total_pts
    );
}

fn write_step(brep: &BRep, path: &str) {
    let step = StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    std::fs::write(path, step).expect("write STEP");
    println!("  -> {path}");
}
