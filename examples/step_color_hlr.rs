//! Example: Phase D — STEP color, STEP assembly, BSpline export, and HLR.
//!
//! Demonstrates:
//!   1. STEP colored export — solid-level and per-face colors
//!   2. STEP assembly — multi-BRep with translations and NAUO hierarchy
//!   3. BSpline curve export — quadratic B-spline arc in STEP
//!   4. HLR (Hidden-Line Removal) — isometric, front, and top SVG views
//!
//! Run: cargo run --example phase_d_demo
//!
//! Output files (written to current directory):
//!   colored_box.step
//!   multicolor_box.step
//!   assembly.step
//!   bspline_arc.step
//!   hlr_isometric.svg
//!   hlr_front.svg
//!   hlr_top.svg

use glam::DVec3;
use rcad_algorithms::{HlrCamera, hlr, hlr_to_svg};
use rcad_kernel::appearance::{Color, StepColor};
use rcad_kernel::{BRep, PrimitiveSolid};
use rcad_step::ExportSelection;
use rcad_step::StepWriter;
use rcad_step::{AssemblyComponent, write_assembly};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn save(path: &str, content: &str) {
    std::fs::write(path, content).unwrap();
    println!("  → wrote {path}  ({} bytes)", content.len());
}

// ── 1. Colored STEP export ────────────────────────────────────────────────────

fn demo_colored_step() {
    println!("\n=== 1. STEP Color Export ===");

    // 1a. Solid-level color (entire box is blue)
    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    let colors = StepColor::new().with_solid_color(Color::BLUE);
    let step = StepWriter::write_string_colored(&brep, &colors);
    save("colored_box.step", &step);
    println!(
        "  Solid color (BLUE): {} STEP records",
        step.lines().filter(|l| l.starts_with('#')).count()
    );

    // 1b. Per-face colors (6-face box with alternating colors)
    let face_colors = [
        Color::RED,
        Color::GREEN,
        Color::BLUE,
        Color::YELLOW,
        Color::CYAN,
        Color::MAGENTA,
    ];
    let mut step_color = StepColor::new().with_solid_color(Color::GRAY);
    for (i, &c) in face_colors.iter().enumerate() {
        step_color = step_color.with_face_color(i, c);
    }
    let step = StepWriter::write_string_colored(&brep, &step_color);
    save("multicolor_box.step", &step);
    println!(
        "  Per-face colors: {} STEP records",
        step.lines().filter(|l| l.starts_with('#')).count()
    );
}

// ── 2. STEP Assembly ─────────────────────────────────────────────────────────

fn demo_assembly() {
    println!("\n=== 2. STEP Assembly ===");

    // Three colored boxes in a row
    let box_brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });

    let components = vec![
        AssemblyComponent::new("RedPart", box_brep.clone()).with_color(Color::RED),
        AssemblyComponent::new("GreenPart", box_brep.clone())
            .with_translation(DVec3::new(1.5, 0.0, 0.0))
            .with_color(Color::GREEN),
        AssemblyComponent::new("BluePart", box_brep.clone())
            .with_translation(DVec3::new(3.0, 0.0, 0.0))
            .with_color(Color::BLUE),
    ];

    let step = write_assembly("ThreeBoxAssembly", &components);
    save("assembly.step", &step);
    println!(
        "  Assembly: {} components, {} STEP records",
        components.len(),
        step.lines().filter(|l| l.starts_with('#')).count()
    );
    let nauo_count = step
        .lines()
        .filter(|l| l.contains("NEXT_ASSEMBLY_USAGE_OCCURRENCE"))
        .count();
    println!("  NAUO links: {}", nauo_count);
}

// ── 3. BSpline STEP export ────────────────────────────────────────────────────

fn demo_bspline_step() {
    println!("\n=== 3. BSpline STEP Export ===");
    use rcad_kernel::geom::{BSplineCurve3, Curve3};
    use rcad_kernel::topology::{Wire, WireEdge};

    // Quadratic (degree 2) Bezier arc approximating a quarter circle in XY plane:
    // Control points for a quarter-circle arc via rational B-spline
    // Using standard conic section weights: w0=1, w1=cos(π/4)=√2/2, w2=1
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let bspline = BSplineCurve3 {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0],
    };

    // Build a minimal BRep containing a single edge with the B-spline curve
    let mut brep = BRep::new();
    brep.vertices.push(rcad_kernel::topology::Vertex {
        point: DVec3::new(1.0, 0.0, 0.0),
    });
    brep.vertices.push(rcad_kernel::topology::Vertex {
        point: DVec3::new(0.0, 1.0, 0.0),
    });
    let v0 = 0usize;
    let v1 = 1usize;
    let edge_idx = brep.edges.len();
    brep.edges
        .push(rcad_kernel::topology::Edge { start: v0, end: v1 });
    let curve_idx = brep.geom.curves.len();
    brep.geom.curves.push(Curve3::BSpline(bspline));
    brep.geom.edge_curve.push(Some(curve_idx));
    brep.geom.edge_curve_range.push(Some([0.0, 1.0]));
    brep.geom.edge_degenerated.push(false);

    let step = StepWriter::write_string(
        &brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[edge_idx],
        },
    );
    save("bspline_arc.step", &step);
    let has_bspline = step.contains("B_SPLINE_CURVE_WITH_KNOTS");
    println!("  B_SPLINE_CURVE_WITH_KNOTS present: {}", has_bspline);
    println!("  Degree-2 quadratic arc, control points: (1,0,0)→(1,1,0)→(0,1,0)");
}

// ── 4. HLR → SVG ─────────────────────────────────────────────────────────────

fn demo_hlr() {
    println!("\n=== 4. Hidden-Line Removal (HLR) ===");

    let brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 1.5,
        depth: 1.0,
    });

    let views = [
        ("hlr_isometric.svg", HlrCamera::isometric(8.0), "isometric"),
        ("hlr_front.svg", HlrCamera::front(8.0), "front"),
        ("hlr_top.svg", HlrCamera::top(8.0), "top"),
    ];

    for (filename, camera, label) in &views {
        let result = hlr(&brep, camera, 16);
        let svg = hlr_to_svg(&result, 120.0, 30.0);
        save(filename, &svg);
        println!(
            "  {label}: {} visible, {} hidden segments",
            result.visible().count(),
            result.hidden().count()
        );
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔═══════════════════════════════════════╗");
    println!("║         RCAD Phase D Demo             ║");
    println!("║  STEP Color · Assembly · BSpline · HLR║");
    println!("╚═══════════════════════════════════════╝");

    demo_colored_step();
    demo_assembly();
    demo_bspline_step();
    demo_hlr();

    println!("\n✓ Phase D demo complete. Check the generated .step and .svg files.");
}
