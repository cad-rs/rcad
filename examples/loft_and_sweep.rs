//! Example: Phase E — B-Spline Surface STEP Read, Loft, and Pipe Sweep.
//!
//! Demonstrates:
//!   1. B_SPLINE_SURFACE_WITH_KNOTS — synthetic STEP string parsed into BRep with analytic surface
//!   2. Loft — two square cross-sections connected into a ruled solid
//!   3. Pipe sweep — hexagonal profile swept along an S-curve spine
//!
//! Run: cargo run --example phase_e_demo
//!
//! Output files (written to current directory):
//!   loft_solid.step
//!   pipe_sweep.step

use glam::{DVec2, DVec3};
use rcad_modeling::{loft, sweep_pipe};
use rcad_step::{StepReader, StepWriter, ExportSelection};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn save(path: &str, content: &str) {
    std::fs::write(path, content).unwrap();
    println!("  → wrote {path}  ({} bytes)", content.len());
}

fn step_record_count(step: &str) -> usize {
    step.lines().filter(|l| l.starts_with('#')).count()
}

// ── 1. B-Spline Surface STEP read ─────────────────────────────────────────────

fn demo_bspline_surface_step_read() {
    println!("\n=== 1. B_SPLINE_SURFACE_WITH_KNOTS STEP Read ===");

    // Minimal valid STEP AP214 file containing a single bilinear (degree 1×1)
    // B-Spline surface patch over a 2×2 control grid.
    // Control points form a flat unit square in the XY plane.
    let step_str = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('B-Spline surface test'),'2;1');
FILE_NAME('bspline_surface_test.step','2026-04-03T00:00:00',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 3 1 1 }'));
ENDSEC;
DATA;
/* Control points: 2x2 grid (v-rows, u-cols) */
#1=CARTESIAN_POINT('',(0.0,0.0,0.0));
#2=CARTESIAN_POINT('',(1.0,0.0,0.0));
#3=CARTESIAN_POINT('',(0.0,1.0,0.0));
#4=CARTESIAN_POINT('',(1.0,1.0,0.0));
/* B-Spline surface: degree 1 in both u and v, 2x2 control grid */
#10=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.0,1.0),(0.0,1.0),.UNSPECIFIED.);
ENDSEC;
END-ISO-10303-21;
"#;

    match StepReader::parse_string(step_str) {
        Ok(brep) => {
            let n_surfaces = brep.geom.surfaces.len();
            println!("  Parsed OK — {} surface(s) in GeomStore", n_surfaces);
            if n_surfaces > 0 {
                let s = &brep.geom.surfaces[0];
                println!("  Surface type: {:?}", std::mem::discriminant(s));
            } else {
                println!("  (No full BRep topology — standalone surface entity only; see ADVANCED_FACE for full topology)");
            }
        }
        Err(e) => {
            println!("  Parse returned: {e}");
            println!("  (Expected: B-Spline surface parses without panicking)");
        }
    }
}

// ── 2. Loft ──────────────────────────────────────────────────────────────────

fn demo_loft() {
    println!("\n=== 2. Loft (square → square at height 2) ===");

    // Two square profiles at z=0 and z=2
    let profile_bot: Vec<DVec3> = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ];
    let profile_top: Vec<DVec3> = vec![
        DVec3::new(0.25, 0.25, 2.0),
        DVec3::new(0.75, 0.25, 2.0),
        DVec3::new(0.75, 0.75, 2.0),
        DVec3::new(0.25, 0.75, 2.0),
    ];

    let brep = loft(&[profile_bot, profile_top]).expect("loft failed");

    let n_faces = brep.solids.first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0);
    println!("  Faces: {} (expected 6 = 2 caps + 4 lateral)", n_faces);

    let step = StepWriter::write_string(&brep, ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    });
    save("loft_solid.step", &step);
    println!("  STEP records: {}", step_record_count(&step));
}

// ── 3. Pipe Sweep ─────────────────────────────────────────────────────────────

fn demo_pipe_sweep() {
    println!("\n=== 3. Pipe Sweep (hexagonal profile along S-curve spine) ===");

    // Hexagonal cross-section (radius 0.3)
    let r = 0.3_f64;
    let profile_2d: Vec<DVec2> = (0..6)
        .map(|i| {
            let a = i as f64 * std::f64::consts::TAU / 6.0;
            DVec2::new(r * a.cos(), r * a.sin())
        })
        .collect();

    // S-curve spine: goes forward in Z while snaking in X
    let n_spine = 12;
    let spine: Vec<DVec3> = (0..n_spine)
        .map(|i| {
            let t = i as f64 / (n_spine - 1) as f64;
            let z = t * 4.0;
            let x = (t * std::f64::consts::PI).sin() * 0.8;
            DVec3::new(x, 0.0, z)
        })
        .collect();

    let brep = sweep_pipe(&profile_2d, &spine).expect("sweep_pipe failed");

    let n_faces = brep.solids.first()
        .and_then(|s| s.shells.first())
        .map(|sh| sh.faces.len())
        .unwrap_or(0);
    let expected_lateral = (n_spine - 1) * 6; // (stations-1) × hex_sides
    println!("  Spine stations: {n_spine}, hex sides: 6");
    println!("  Faces: {n_faces} (expected {} lateral + 2 caps = {})", expected_lateral, expected_lateral + 2);

    let step = StepWriter::write_string(&brep, ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    });
    save("pipe_sweep.step", &step);
    println!("  STEP records: {}", step_record_count(&step));
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    println!("╔═══════════════════════════════════════════╗");
    println!("║            RCAD Phase E Demo              ║");
    println!("║   B-Spline Surface · Loft · Pipe Sweep   ║");
    println!("╚═══════════════════════════════════════════╝");

    demo_bspline_surface_step_read();
    demo_loft();
    demo_pipe_sweep();

    println!("\n✓ Phase E demo complete. Check loft_solid.step and pipe_sweep.step.");
}
