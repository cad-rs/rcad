//! Example: Properties, validity checking, and section planes.
//!
//! Demonstrates:
//!   1. `surface_area`, `volume`, `centroid` from rcad-kernel
//!   2. `check` (BRepCheck) from rcad-algorithms
//!   3. `section_polylines` from rcad-algorithms
//!
//! Run: cargo run --example phase_c_demo

use glam::DVec3;
use rcad_algorithms::{check, section_polylines};
use rcad_kernel::geom::Plane;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{BRep, PrimitiveSolid, centroid, surface_area, volume};
use rcad_modeling::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::{extrude, revolve};
use std::f64::consts::TAU;

fn main() {
    // ── 1. Unit box properties ─────────────────────────────────────────
    println!("=== 1. Unit Box Properties ===");
    let unit_box = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    println!("  Surface area : {:.4}", surface_area(&unit_box)); // 6.0
    println!("  Volume       : {:.4}", volume(&unit_box)); // 1.0
    println!("  Centroid     : {:.4}", centroid(&unit_box)); // (0.5,0.5,0.5)

    // ── 2. 2×3×4 box properties ────────────────────────────────────────
    println!("\n=== 2. 2×3×4 Box Properties ===");
    let box234 = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 3.0,
        depth: 4.0,
    });
    println!("  Surface area : {:.4}", surface_area(&box234)); // 52.0
    println!("  Volume       : {:.4}", volume(&box234)); // 24.0
    println!("  Centroid     : {:.4}", centroid(&box234));

    // ── 3. Extruded square properties ─────────────────────────────────
    println!("\n=== 3. Extruded Square (1×1×2) Properties ===");
    let sq_profile = square_profile();
    let pillar = extrude(&sq_profile, 0, DVec3::Z, 2.0).expect("extrude");
    println!("  Surface area : {:.4}", surface_area(&pillar)); // 4*2 + 2*1 = 10
    println!("  Volume       : {:.4}", volume(&pillar)); // 2.0
    println!("  Centroid     : {:.4}", centroid(&pillar));

    // ── 4. BRepCheck on valid shape ────────────────────────────────────
    println!("\n=== 4. BRepCheck on Unit Box ===");
    let result = check(&unit_box);
    if result.is_valid() {
        println!("  ✓ Unit box passed all validity checks");
    } else {
        println!("  ✗ Issues found:");
        for issue in &result.issues {
            println!("    - {issue}");
        }
    }

    // ── 5. BRepCheck on extruded pillar ───────────────────────────────
    println!("\n=== 5. BRepCheck on Extruded Pillar ===");
    let result2 = check(&pillar);
    if result2.is_valid() {
        println!("  ✓ Extruded pillar passed all validity checks");
    } else {
        println!("  Issues ({} total):", result2.issues.len());
        for issue in &result2.issues {
            println!("    - {issue}");
        }
    }

    // ── 6. Section of unit box at z=0.5 ───────────────────────────────
    println!("\n=== 6. Section of Unit Box at z=0.5 ===");
    let plane_z = Plane {
        origin: DVec3::new(0.0, 0.0, 0.5),
        normal: DVec3::Z,
    };
    let loops = section_polylines(&unit_box, &plane_z);
    println!("  Loops found: {}", loops.len());
    for (i, lp) in loops.iter().enumerate() {
        println!("  Loop {i}: {} points", lp.len());
        for p in lp {
            println!("    ({:.3}, {:.3}, {:.3})", p.x, p.y, p.z);
        }
    }

    // ── 7. Section of 2×3×4 box at y=1.5 ─────────────────────────────
    println!("\n=== 7. Section of 2×3×4 Box at y=1.5 ===");
    let plane_y = Plane {
        origin: DVec3::new(0.0, 1.5, 0.0),
        normal: DVec3::Y,
    };
    let loops2 = section_polylines(&box234, &plane_y);
    println!("  Loops found: {}", loops2.len());
    for (i, lp) in loops2.iter().enumerate() {
        let xs: Vec<_> = lp.iter().map(|p| format!("{:.2}", p.x)).collect();
        let zs: Vec<_> = lp.iter().map(|p| format!("{:.2}", p.z)).collect();
        println!(
            "  Loop {i}: {} points  x=[{}..{}] z=[{}..{}]",
            lp.len(),
            xs.iter().min().unwrap_or(&"-".to_string()),
            xs.iter().max().unwrap_or(&"-".to_string()),
            zs.iter().min().unwrap_or(&"-".to_string()),
            zs.iter().max().unwrap_or(&"-".to_string()),
        );
    }

    // ── 8. Section of sphere at equator ───────────────────────────────
    println!("\n=== 8. Section of Revolve (ring) at z=0 ===");
    let sq = rect_profile(2.0, -0.25, 0.5, 0.5);
    let ring = revolve(&sq, 0, DVec3::ZERO, DVec3::Y, TAU).expect("revolve ring");
    let plane_eq = Plane {
        origin: DVec3::new(0.0, 0.0, 0.0),
        normal: DVec3::Z,
    };
    let loops3 = section_polylines(&ring, &plane_eq);
    println!("  Ring section loops: {}", loops3.len());

    // ── 9. Properties of sphere primitive ─────────────────────────────
    println!("\n=== 9. Sphere (r=1) Properties ===");
    let sphere = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
    // Sphere has no pre-triangulated data — properties use wire fan (minimal topology)
    println!(
        "  Faces        : {}",
        sphere
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    );
    let result3 = check(&sphere);
    // Sphere uses seam topology — degenerate by polygon checker standards
    println!("  BRepCheck issues: {}", result3.issues.len());
    for issue in &result3.issues {
        println!("    - {issue}");
    }

    println!("\nProperties demo complete.");
}

// ── Profile helpers ────────────────────────────────────────────────────────────

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
    let surface = Surface3::Plane(rcad_kernel::geom::Plane {
        origin: pts[0],
        normal: DVec3::Z,
    });
    make_face(&mut brep, surface, make_wire(wires), vec![]).unwrap();
    brep
}

fn square_profile() -> BRep {
    polygon_profile(&[
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(1.0, 0.0, 0.0),
        DVec3::new(1.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
    ])
}

fn rect_profile(x: f64, y: f64, w: f64, h: f64) -> BRep {
    polygon_profile(&[
        DVec3::new(x, y, 0.0),
        DVec3::new(x + w, y, 0.0),
        DVec3::new(x + w, y + h, 0.0),
        DVec3::new(x, y + h, 0.0),
    ])
}
