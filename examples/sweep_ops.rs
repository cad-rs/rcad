//! Example: Linear extrusion and revolution using the sweep API.
//!
//! Demonstrates `extrude()` and `revolve()` from rcad-modeling.
//!
//! Run: cargo run --example sweep_ops

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_modeling::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::{extrude, revolve};
use rcad_scene::append_brep;
use rcad_step::writer::{ExportSelection, StepWriter};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

fn main() {
    // ── 1. Unit square extruded 2 units along Z ────────────────────────
    println!("1. Square extrusion");
    let square = square_profile(0.0, 0.0, 1.0, 1.0);
    let pillar = extrude(&square, 0, DVec3::Z, 2.0).expect("extrude square");
    write_step(&pillar, "output_extrude_square.step");

    // ── 2. Rectangle extruded along Y (flat plate) ────────────────────
    println!("2. Rectangle extrusion (plate)");
    let rect = rect_profile(0.0, 0.0, 4.0, 0.5);
    let plate = extrude(&rect, 0, DVec3::Y, 0.2).expect("extrude rectangle");
    write_step(&plate, "output_extrude_plate.step");

    // ── 3. Triangle extruded along Z ──────────────────────────────────
    println!("3. Triangle extrusion (prism)");
    let tri = triangle_profile(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(1.0, 1.732, 0.0), // equilateral, side = 2
    );
    let prism = extrude(&tri, 0, DVec3::Z, 3.0).expect("extrude triangle");
    write_step(&prism, "output_extrude_prism.step");

    // ── 4. L-profile extrusion ────────────────────────────────────────
    println!("4. L-profile extrusion");
    let l_prof = l_profile(0.0, 0.0, 3.0, 3.0, 0.5);
    let l_beam = extrude(&l_prof, 0, DVec3::Z, 5.0).expect("extrude L");
    write_step(&l_beam, "output_extrude_l_beam.step");

    // ── 5. Hexagon extruded along Z (hex bolt head) ───────────────────
    println!("5. Hexagon extrusion (bolt head)");
    let hex = hexagon_profile(DVec3::ZERO, 1.0);
    let bolt_head = extrude(&hex, 0, DVec3::Z, 0.6).expect("extrude hex");
    write_step(&bolt_head, "output_extrude_hex.step");

    // ── 6. Square revolved 360° around Y axis (torus-like ring) ───────
    println!("6. Square revolve full 360°");
    let sq = square_profile(2.0, -0.5, 0.5, 0.5); // offset from axis
    let ring = revolve(&sq, 0, DVec3::ZERO, DVec3::Y, TAU).expect("revolve full");
    write_step(&ring, "output_revolve_ring.step");

    // ── 7. Triangle revolved 360° → cone-like shape ───────────────────
    println!("7. Triangle revolve 360° (cone-like)");
    let wedge = triangle_profile(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.5, 0.0, 0.0),
        DVec3::new(0.0, 3.0, 0.0),
    );
    let cone_like = revolve(&wedge, 0, DVec3::ZERO, DVec3::Y, TAU).expect("revolve triangle");
    write_step(&cone_like, "output_revolve_cone.step");

    // ── 8. Rectangle revolved 90° (quarter torus) ─────────────────────
    println!("8. Rectangle revolve 90° (elbow)");
    let r = rect_profile(3.0, -0.25, 0.5, 0.5);
    let elbow = revolve(&r, 0, DVec3::ZERO, DVec3::Z, FRAC_PI_2).expect("revolve 90");
    write_step(&elbow, "output_revolve_elbow.step");

    // ── 9. Rectangle revolved 180° (half pipe) ───────────────────────
    println!("9. Rectangle revolve 180° (half pipe)");
    let r2 = rect_profile(2.0, -0.25, 0.4, 0.5);
    let half_pipe = revolve(&r2, 0, DVec3::ZERO, DVec3::Z, PI).expect("revolve 180");
    write_step(&half_pipe, "output_revolve_half_pipe.step");

    // ── 10. Assembly: stacked extrusions ─────────────────────────────
    println!("10. Assembly: stacked extrusions");
    let base_sq = square_profile(-1.5, -1.5, 3.0, 3.0);
    let mut base = extrude(&base_sq, 0, DVec3::Z, 0.5).expect("base");

    let mid_sq = square_profile(-1.0, -1.0, 2.0, 2.0);
    let mut mid = extrude(&mid_sq, 0, DVec3::Z, 0.5).expect("mid");
    for v in &mut mid.vertices {
        v.point.z += 0.5;
    }

    let top_sq = square_profile(-0.5, -0.5, 1.0, 1.0);
    let mut top = extrude(&top_sq, 0, DVec3::Z, 0.5).expect("top");
    for v in &mut top.vertices {
        v.point.z += 1.0;
    }

    append_brep(&mut base, mid);
    append_brep(&mut base, top);
    write_step(&base, "output_stacked_extrusions.step");

    println!("\nExported 10 sweep operation STEP files.");
}

// ── Profile builders ──────────────────────────────────────────────────────────

/// Axis-aligned rectangle at (x, y, 0) with given width and height.
fn rect_profile(x: f64, y: f64, w: f64, h: f64) -> BRep {
    let pts = [
        DVec3::new(x, y, 0.0),
        DVec3::new(x + w, y, 0.0),
        DVec3::new(x + w, y + h, 0.0),
        DVec3::new(x, y + h, 0.0),
    ];
    polygon_profile(&pts)
}

/// Unit square profile.
fn square_profile(x: f64, y: f64, w: f64, h: f64) -> BRep {
    rect_profile(x, y, w, h)
}

/// Triangle profile from 3 points (in XY plane).
fn triangle_profile(a: DVec3, b: DVec3, c: DVec3) -> BRep {
    polygon_profile(&[a, b, c])
}

/// Regular hexagon centered at `center` in the XY plane with circumradius `r`.
fn hexagon_profile(center: DVec3, r: f64) -> BRep {
    let pts: Vec<DVec3> = (0..6)
        .map(|i| {
            let angle = std::f64::consts::FRAC_PI_3 * i as f64;
            center + DVec3::new(r * angle.cos(), r * angle.sin(), 0.0)
        })
        .collect();
    polygon_profile(&pts)
}

/// L-profile: outer rectangle (w x h) with a notch removed from top-right.
/// `t` is the thickness of both legs.
///
///   ┌──┐
///   │  │
///   │  └───┐
///   └──────┘
fn l_profile(x: f64, y: f64, w: f64, h: f64, t: f64) -> BRep {
    let pts = [
        DVec3::new(x, y, 0.0),
        DVec3::new(x + w, y, 0.0),
        DVec3::new(x + w, y + t, 0.0),
        DVec3::new(x + t, y + t, 0.0),
        DVec3::new(x + t, y + h, 0.0),
        DVec3::new(x, y + h, 0.0),
    ];
    polygon_profile(&pts)
}

/// Build a closed planar BRep face from an ordered polygon in the XY plane.
fn polygon_profile(pts: &[DVec3]) -> BRep {
    assert!(pts.len() >= 3, "polygon needs at least 3 points");
    let n = pts.len();
    let mut brep = BRep::default();

    let vis: Vec<usize> = pts.iter().map(|&p| make_vertex(&mut brep, p)).collect();

    let mut wire_edges = Vec::new();
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
        .expect("make_edge");
        wire_edges.push(WireEdge {
            idx: eidx,
            forward: true,
        });
    }

    let surface = Surface3::Plane(Plane {
        origin: pts[0],
        normal: DVec3::Z,
    });
    make_face(&mut brep, surface, make_wire(wire_edges), vec![]).expect("make_face");
    brep
}

fn write_step(brep: &BRep, path: &str) {
    let step = StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &[],
            selected_edges: &[],
        },
    );
    std::fs::write(path, step).expect("write STEP file");
    println!("  -> {path}");
}
