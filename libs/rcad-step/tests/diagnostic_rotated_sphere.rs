//! Diagnostic: rotated sphere boolean hang investigation.
//!
//! This test replicates OCCT bcut_simple A2 but with timeouts and
//! per-phase instrumentation to identify WHERE the boolean hangs.
//!
//! OCCT script:
//!   psphere s 1
//!   trotate s 0 0 0 0 0 1 -90    (Rz -90°)
//!   trotate s 0 0 0 0 1 0 -45    (Ry -45°)
//!   box b 1 1 1
//!   bcut result s b
//!   checkprops result -s 13.3517

use glam::{DAffine3, DVec3};
use rcad_algorithms::{boolean_op, BooleanOpType, total_surface_area};
use rcad_modeling::{make_box_brep, make_sphere_brep};
use std::time::Instant;

#[test]
fn diagnostic_rotated_sphere_a2() {
    // Phase 1: build sphere
    let s = make_sphere_brep(DVec3::ZERO, 1.0).expect("sphere");
    eprintln!("Sphere built: {} vertices", s.vertices.len());

    // Phase 2: rotate Rz(-90°) then Ry(-45°)
    let rz = DAffine3::from_rotation_z(-std::f64::consts::FRAC_PI_2);
    let ry = DAffine3::from_rotation_y(-std::f64::consts::FRAC_PI_4);
    let mut s = s;
    s.apply_transform(rz);
    s.apply_transform(ry);
    eprintln!("Sphere rotated");
    // Quick sanity: check center (should still be at origin since box is at origin)
    let nv = s.vertices.len() as f64;
    let center = s.vertices.iter().map(|v| v.point).reduce(|a, b| a + b).map(|s| s / nv);
    eprintln!("Sphere vertex count: {}, approx center: {:?}", s.vertices.len(), center);
    // Check surface geometry
    eprintln!("Sphere surfaces:");
    for (i, surf) in s.geom.surfaces.iter().enumerate() {
        eprintln!("  surface {i}: {surf:?}");
    }

    // Phase 3: build box
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box");
    eprintln!("Box built: {} vertices", b.vertices.len());

    // Phase 4: boolean difference
    eprintln!("Starting boolean_op(Difference)...");
    let start = Instant::now();
    let result = match boolean_op(BooleanOpType::Difference, &s, &b) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Boolean failed after {:.1}s: {e}", start.elapsed().as_secs_f64());
            panic!("boolean_op failed: {e}");
        }
    };
    let elapsed = start.elapsed();
    eprintln!("Boolean completed in {:.1}s", elapsed.as_secs_f64());

    // Phase 5: check result structure (before surface area)
    eprintln!("Result solids: {}", result.solids.len());
    for (si, solid) in result.solids.iter().enumerate() {
        eprintln!("  solid {si}: {} shells", solid.shells.len());
        for (shi, shell) in solid.shells.iter().enumerate() {
            eprintln!("    shell {shi}: {} faces", shell.faces.len());
            for (fi, face) in shell.faces.iter().enumerate() {
                eprintln!("      face {fi}: {} tris, {} outer_wire edges, {} inner_wires",
                    face.triangles.len(),
                    face.outer_wire.edges.len(),
                    face.inner_wires.len());
            }
        }
    }

    // Phase 6: per-face diagnostics
    eprintln!("\n--- Per-face diagnostics ---");
    use rcad_kernel::BRep;
    use rcad_kernel::geom::Surface3;
    use rcad_kernel::properties::face_surface_area;

    fn face_flat_iter(brep: &BRep) -> impl Iterator<Item = (usize, &rcad_kernel::topology::Face)> {
        let mut idx = 0;
        brep.solids.iter().flat_map(move |solid| {
            solid.shells.iter().flat_map(move |shell| {
                shell.faces.iter().map(move |f| {
                    let fi = idx;
                    idx += 1;
                    (fi, f)
                })
            })
        })
    }

    for (fi, face) in face_flat_iter(&result) {
        let surf_name = result.geom.face_surface.get(fi)
            .and_then(|o| *o)
            .and_then(|sidx| result.geom.surfaces.get(sidx))
            .map(|s| match s {
                Surface3::Plane(_) => "Plane",
                Surface3::Sphere(_) => "Sphere",
                Surface3::Cylinder(_) => "Cylinder",
                Surface3::Cone(_) => "Cone",
                Surface3::Torus(_) => "Torus",
                Surface3::BSpline(_) => "BSpline",
                _ => "Other",
            })
            .unwrap_or("None");
        eprintln!("  face {fi}: {surf_name}, {} outer edges, {} inner wires, {} tris",
            face.outer_wire.edges.len(),
            face.inner_wires.len(),
            face.triangles.len());
    }

    // Phase 7: total_surface_area
    eprintln!("Computing total_surface_area...");
    let sa_start = Instant::now();
    let sa = total_surface_area(&result);
    eprintln!("total_surface_area took {:.1}s, result: {sa:.4}", sa_start.elapsed().as_secs_f64());

    // Per-face area breakdown
    eprintln!("\n--- Per-face surface area ---");
    let mut per_face_areas: Vec<f64> = Vec::new();
    for (fi, face) in face_flat_iter(&result) {
        let a = rcad_kernel::properties::face_surface_area(&result, face, fi);
        per_face_areas.push(a);
        let surf_name = result.geom.face_surface.get(fi)
            .and_then(|o| *o)
            .and_then(|sidx| result.geom.surfaces.get(sidx))
            .map(|s| match s {
                Surface3::Plane(_) => "Plane",
                Surface3::Sphere(_) => "Sphere",
                _ => "Other",
            })
            .unwrap_or("None");
        eprintln!("  face {fi} ({surf_name}): area = {a:.4}, edges = {}", face.outer_wire.edges.len());
    }
    eprintln!("  Total from per-face: {:.4}", per_face_areas.iter().sum::<f64>());

    // The expected SA from OCCT is 13.3517
    let diff = (sa - 13.3517).abs();
    eprintln!("SA diff from expected: {diff:.4}");
    assert!(diff < 2.0, "SA mismatch: {sa:.4} vs expected 13.3517 (diff={diff:.4})");
}
