//! Example: Free-form BRep construction with BRepBuilder API.
//!
//! Shows how to use `make_vertex`, `make_edge`, `make_wire`, `make_face`,
//! `make_solid` to build shapes without relying on primitive constructors.
//!
//! Run: cargo run --example brep_builder_demo

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::BRep;
use rcad_modeling::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_scene::append_brep;
use rcad_step::writer::{ExportSelection, StepWriter};

fn main() {
    // ── 1. Manually built tetrahedron (4 triangular faces) ────────────
    println!("1. Tetrahedron (manually wired)");
    let tet = tetrahedron(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(2.0, 0.0, 0.0),
        DVec3::new(1.0, 1.732, 0.0),
        DVec3::new(1.0, 0.577, 1.633),
    );
    write_step(&tet, "output_tetrahedron.step");

    // ── 2. Flat diamond (rhombus in XY plane) ─────────────────────────
    println!("2. Diamond / rhombus face");
    let diamond = flat_diamond(DVec3::ZERO, 2.0, 1.0);
    write_step(&diamond, "output_diamond.step");

    // ── 3. Star polygon (6-pointed, two triangles) ────────────────────
    println!("3. Star / hexagram (two overlaid triangles)");
    let star = hexagram(DVec3::ZERO, 1.5);
    write_step(&star, "output_star.step");

    // ── 4. Stepped cross (plus sign in XY plane) ─────────────────────
    println!("4. Plus-sign / cross profile");
    let cross = plus_profile(DVec3::ZERO, 3.0, 1.0);
    write_step(&cross, "output_cross_profile.step");

    // ── 5. Pentagon face ──────────────────────────────────────────────
    println!("5. Regular pentagon");
    let penta = regular_polygon(DVec3::ZERO, 1.5, 5);
    write_step(&penta, "output_pentagon.step");

    // ── 6. Octagon face ───────────────────────────────────────────────
    println!("6. Regular octagon");
    let oct = regular_polygon(DVec3::ZERO, 1.5, 8);
    write_step(&oct, "output_octagon.step");

    // ── 7. Right-angle triangle with precise dimensions ────────────────
    println!("7. Right-angle triangle (3-4-5)");
    let right_tri = triangle_face(
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(3.0, 0.0, 0.0),
        DVec3::new(0.0, 4.0, 0.0),
    );
    write_step(&right_tri, "output_right_triangle.step");

    // ── 8. T-profile (structural beam cross-section) ──────────────────
    println!("8. T-profile face");
    let t_beam = t_profile(DVec3::ZERO, 4.0, 4.0, 0.5, 0.5);
    write_step(&t_beam, "output_t_profile.step");

    // ── 9. Assembly: three faces at different heights ─────────────────
    println!("9. Assembly: polygon stack at different heights");
    let mut stack = regular_polygon(DVec3::new(-3.0, 0.0, 0.0), 1.0, 3);

    let sq = flat_quad(
        DVec3::new(0.0, -0.75, 0.0),
        DVec3::new(1.5, -0.75, 0.0),
        DVec3::new(1.5, 0.75, 0.0),
        DVec3::new(0.0, 0.75, 0.0),
    );
    append_brep(&mut stack, sq);

    let hex = regular_polygon(DVec3::new(5.0, 0.0, 0.0), 1.0, 6);
    append_brep(&mut stack, hex);

    let oct = regular_polygon(DVec3::new(9.0, 0.0, 0.0), 1.0, 8);
    append_brep(&mut stack, oct);

    write_step(&stack, "output_polygon_gallery.step");

    println!("\nExported 9 BRepBuilder demo STEP files.");
}

// ── Shape builders ────────────────────────────────────────────────────────────

/// Build a regular N-gon face in the XY plane centered at `center`.
fn regular_polygon(center: DVec3, r: f64, n: usize) -> BRep {
    let pts: Vec<DVec3> = (0..n).map(|i| {
        let angle = std::f64::consts::TAU * i as f64 / n as f64;
        center + DVec3::new(r * angle.cos(), r * angle.sin(), 0.0)
    }).collect();
    planar_face_xy(&pts)
}

/// A flat quad from 4 explicit corner points.
fn flat_quad(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> BRep {
    planar_face_xy(&[a, b, c, d])
}

/// A flat triangle from 3 explicit points.
fn triangle_face(a: DVec3, b: DVec3, c: DVec3) -> BRep {
    planar_face_xy(&[a, b, c])
}

/// A diamond (rhombus): `w` = horizontal half-span, `h` = vertical half-span.
fn flat_diamond(center: DVec3, w: f64, h: f64) -> BRep {
    planar_face_xy(&[
        center + DVec3::new(w,  0.0, 0.0),
        center + DVec3::new(0.0, h,  0.0),
        center + DVec3::new(-w, 0.0, 0.0),
        center + DVec3::new(0.0, -h, 0.0),
    ])
}

/// A 6-pointed star (hexagram) built as two overlapping triangles (6 faces total).
fn hexagram(center: DVec3, r: f64) -> BRep {
    let inner = r * 0.5;
    // Star polygon — 12 alternating outer/inner vertices
    let pts: Vec<DVec3> = (0..12).map(|i| {
        let angle = std::f64::consts::TAU * i as f64 / 12.0;
        let radius = if i % 2 == 0 { r } else { inner };
        center + DVec3::new(radius * angle.cos(), radius * angle.sin(), 0.0)
    }).collect();
    planar_face_xy(&pts)
}

/// Plus/cross profile centered at `center`: arm_len × arm_len, thickness `t`.
fn plus_profile(center: DVec3, arm_len: f64, t: f64) -> BRep {
    let h = arm_len / 2.0;
    let th = t / 2.0;
    let c = center;
    // 12-vertex plus sign (CCW)
    let pts = [
        c + DVec3::new(-th, -h, 0.0),
        c + DVec3::new( th, -h, 0.0),
        c + DVec3::new( th, -th, 0.0),
        c + DVec3::new(  h, -th, 0.0),
        c + DVec3::new(  h,  th, 0.0),
        c + DVec3::new( th,  th, 0.0),
        c + DVec3::new( th,   h, 0.0),
        c + DVec3::new(-th,   h, 0.0),
        c + DVec3::new(-th,  th, 0.0),
        c + DVec3::new( -h,  th, 0.0),
        c + DVec3::new( -h, -th, 0.0),
        c + DVec3::new(-th, -th, 0.0),
    ];
    planar_face_xy(&pts)
}

/// T-profile: flange width `fw`, total height `fh`, web thickness `wt`, flange thickness `ft`.
///
///  ┌──────────────┐  ← flange (fw wide, ft thick)
///  └────┐    ┌────┘
///       │    │       ← web (wt wide, fh-ft tall)
///       └────┘
fn t_profile(center: DVec3, fw: f64, fh: f64, wt: f64, ft: f64) -> BRep {
    let hw = fw / 2.0;
    let hwt = wt / 2.0;
    let c = center;
    let pts = [
        c + DVec3::new(-hw,   0.0,     0.0),
        c + DVec3::new( hw,   0.0,     0.0),
        c + DVec3::new( hw,   ft,      0.0),
        c + DVec3::new( hwt,  ft,      0.0),
        c + DVec3::new( hwt,  fh,      0.0),
        c + DVec3::new(-hwt,  fh,      0.0),
        c + DVec3::new(-hwt,  ft,      0.0),
        c + DVec3::new(-hw,   ft,      0.0),
    ];
    planar_face_xy(&pts)
}

/// A tetrahedron: 4 triangular faces, each emitted as a separate face.
fn tetrahedron(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> BRep {
    let mut brep = BRep::default();

    // Emit the 4 triangular faces
    emit_tri_face(&mut brep, a, b, c);
    emit_tri_face(&mut brep, a, b, d);
    emit_tri_face(&mut brep, b, c, d);
    emit_tri_face(&mut brep, a, c, d);

    brep
}

/// Emit a single triangular face into an existing BRep.
fn emit_tri_face(brep: &mut BRep, a: DVec3, b: DVec3, c: DVec3) {
    let va = make_vertex(brep, a);
    let vb = make_vertex(brep, b);
    let vc = make_vertex(brep, c);

    let mk_edge = |brep: &mut BRep, p: DVec3, q: DVec3, vp: usize, vq: usize| {
        let dir = (q - p).normalize_or_zero();
        let len = (q - p).length();
        make_edge(brep, Curve3::Line(Line3 { origin: p, direction: dir }), 0.0, len, vp, vq)
            .expect("edge")
    };

    let e0 = mk_edge(brep, a, b, va, vb);
    let e1 = mk_edge(brep, b, c, vb, vc);
    let e2 = mk_edge(brep, c, a, vc, va);

    let normal = (b - a).cross(c - a).normalize_or_zero();
    let surface = Surface3::Plane(Plane { origin: a, normal });
    let wire = make_wire(vec![
        WireEdge { idx: e0, forward: true },
        WireEdge { idx: e1, forward: true },
        WireEdge { idx: e2, forward: true },
    ]);
    make_face(brep, surface, wire, vec![]).expect("face");
}

/// Build a planar face in the XY plane (normal = Z) from an ordered polygon.
fn planar_face_xy(pts: &[DVec3]) -> BRep {
    assert!(pts.len() >= 3);
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
            Curve3::Line(Line3 { origin: a, direction: dir }),
            0.0, len,
            vis[i], vis[j],
        ).expect("edge");
        wire_edges.push(WireEdge { idx: eidx, forward: true });
    }

    let surface = Surface3::Plane(Plane { origin: pts[0], normal: DVec3::Z });
    make_face(&mut brep, surface, make_wire(wire_edges), vec![]).expect("face");
    brep
}

fn write_step(brep: &BRep, path: &str) {
    let step = StepWriter::write_string(brep, ExportSelection {
        selected_faces: &[],
        selected_edges: &[],
    });
    std::fs::write(path, step).expect("write STEP file");
    println!("  -> {path}");
}
