//! Example: Shell and solid offset operations.
//!
//! Demonstrates:
//!   1. Offset shell (outward/inward)
//!   2. Offset solid
//!   3. Hollow solid (thin-wall with openings)
//!   4. Thicken shell with lateral faces
//!   5. Self-intersection detection
//!
//! Run:
//!   cargo run -p rcad-examples --example offset_shell

use glam::DVec3;
use rcad_algorithms::{
    offset::{offset_shell, offset_solid, hollow_solid, offset_surface, OffsetOptions},
    thicken::{thicken_shell, thick_solid_with_removed_faces},
    geom_populate,
};
use rcad_kernel::{BRep, PrimitiveSolid};
use rcad_kernel::geom::{Surface3, Plane, SphericalSurface, CylindricalSurface};

fn separator(title: &str) {
    println!("\n──────────────────────────────────────────");
    println!("  {title}");
    println!("──────────────────────────────────────────");
}

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

fn make_box_with_geom() -> BRep {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    });
    geom_populate::populate_box_geom(&mut brep);
    brep
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. Surface offset
// ─────────────────────────────────────────────────────────────────────────────

fn demo_surface_offset() {
    separator("1. Surface Offset");

    // Plane offset
    {
        let plane = Surface3::Plane(Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let offset = offset_surface(&plane, 0.5).unwrap();
        if let Surface3::Plane(p) = offset {
            println!("  Plane offset by 0.5: origin.z = {:.4}", p.origin.z);
            assert!((p.origin.z - 0.5).abs() < 1e-9);
        }
    }

    // Sphere offset
    {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        let offset = offset_surface(&sphere, 0.5).unwrap();
        if let Surface3::Sphere(s) = offset {
            println!("  Sphere offset by 0.5: radius = {:.4} (was 2.0)", s.radius);
            assert!((s.radius - 2.5).abs() < 1e-9);
        }
    }

    // Cylinder offset
    {
        let cylinder = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let offset = offset_surface(&cylinder, 0.3).unwrap();
        if let Surface3::Cylinder(c) = offset {
            println!("  Cylinder offset by 0.3: radius = {:.4} (was 1.0)", c.radius);
            assert!((c.radius - 1.3).abs() < 1e-9);
        }
    }

    // Negative offset (shrinking)
    {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        });
        let offset = offset_surface(&sphere, -0.5).unwrap();
        if let Surface3::Sphere(s) = offset {
            println!("  Sphere offset by -0.5: radius = {:.4} (was 2.0)", s.radius);
            assert!((s.radius - 1.5).abs() < 1e-9);
        }
    }

    // Offset too large (degenerate)
    {
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let offset = offset_surface(&sphere, -2.0);
        println!("  Sphere radius 1.0 offset by -2.0: {:?}", if offset.is_none() { "None (degenerate)" } else { "Some" });
        assert!(offset.is_none(), "offset larger than radius should return None");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. Offset shell
// ─────────────────────────────────────────────────────────────────────────────

fn demo_offset_shell() {
    separator("2. Offset Shell");

    let brep = make_box_with_geom();
    let shell = &brep.solids[0].shells[0];
    let orig_faces = shell.faces.len();

    // Outward offset
    {
        let result = offset_shell(shell, &brep, 0.1);
        assert!(result.is_ok(), "offset_shell should succeed");
        let offset_brep = result.unwrap();
        let offset_faces = face_count(&offset_brep);
        println!("  Outward offset +0.1: {} faces (was {})", offset_faces, orig_faces);
        assert_eq!(offset_faces, orig_faces, "face count should be preserved");
    }

    // Inward offset
    {
        let result = offset_shell(shell, &brep, -0.1);
        assert!(result.is_ok(), "offset_shell with negative distance should succeed");
        let offset_brep = result.unwrap();
        let offset_faces = face_count(&offset_brep);
        println!("  Inward offset -0.1: {} faces (was {})", offset_faces, orig_faces);
    }

    // With options
    {
        let _opts = OffsetOptions::new(0.2)
            .with_tolerance(1e-6)
            .with_self_intersection_check(true);

        let result = offset_shell(shell, &brep, 0.2);
        assert!(result.is_ok());
        println!("  Offset with options: {} faces", face_count(&result.unwrap()));
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. Offset solid
// ─────────────────────────────────────────────────────────────────────────────

fn demo_offset_solid() {
    separator("3. Offset Solid");

    let brep = make_box_with_geom();
    let solid = &brep.solids[0];

    // Outward offset (thickening)
    {
        let result = offset_solid(solid, &brep, 0.2);
        assert!(result.is_ok(), "offset_solid should succeed");
        let offset_brep = result.unwrap();
        println!("  Solid outward offset +0.2:");
        println!("    Vertices: {} (was {})", offset_brep.vertices.len(), brep.vertices.len());
        println!("    Faces: {}", face_count(&offset_brep));
    }

    // Inward offset (shrink)
    {
        let result = offset_solid(solid, &brep, -0.2);
        assert!(result.is_ok(), "offset_solid with negative distance should succeed");
        let offset_brep = result.unwrap();
        println!("  Solid inward offset -0.2:");
        println!("    Vertices: {}", offset_brep.vertices.len());
        println!("    Faces: {}", face_count(&offset_brep));
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. Hollow solid
// ─────────────────────────────────────────────────────────────────────────────

fn demo_hollow_solid() {
    separator("4. Hollow Solid (Thin-Wall)");

    let brep = make_box_with_geom();
    let solid = &brep.solids[0];

    // Remove top face (index 5)
    {
        let result = hollow_solid(solid, &brep, 0.1, &[5]);
        assert!(result.is_ok(), "hollow_solid should succeed");
        let hollow = result.unwrap();
        let hollow_faces = face_count(&hollow);
        println!("  Hollow with top face removed:");
        println!("    Original faces: 6");
        println!("    Hollow faces: {} (kept: 5 offset + lateral)", hollow_faces);
        assert!(hollow_faces >= 5, "should have kept faces plus lateral faces");
    }

    // Remove top and bottom faces
    {
        let result = hollow_solid(solid, &brep, 0.1, &[0, 5]);
        assert!(result.is_ok(), "hollow_solid with multiple open faces should succeed");
        let hollow = result.unwrap();
        println!("  Hollow with top and bottom removed:");
        println!("    Faces: {}", face_count(&hollow));
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Thicken shell
// ─────────────────────────────────────────────────────────────────────────────

fn demo_thicken_shell() {
    separator("5. Thicken Shell");

    // Closed shell (no lateral faces)
    {
        let brep = make_box_with_geom();
        let result = thicken_shell(&brep, 0.1);
        assert!(result.is_some(), "thicken_shell should succeed for closed shell");
        let r = result.unwrap();
        println!("  Closed shell thickened:");
        println!("    Offset faces: {}", r.offset_faces);
        println!("    Lateral faces: {} (none for closed shell)", r.lateral_faces);
    }

    // Open shell (with lateral faces)
    {
        let mut brep = make_box_with_geom();
        // Remove top face to create open shell
        if let Some(s) = brep.solids.first_mut()
            && let Some(sh) = s.shells.first_mut()
                && sh.faces.len() > 1 {
                    sh.faces.pop();
                }

        let result = thicken_shell(&brep, 0.1);
        assert!(result.is_some(), "thicken_shell should succeed for open shell");
        let r = result.unwrap();
        println!("  Open shell thickened:");
        println!("    Offset faces: {}", r.offset_faces);
        println!("    Lateral faces: {} (created at boundary)", r.lateral_faces);
        assert!(r.lateral_faces > 0, "should create lateral faces for open shell");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Thick solid with removed faces
// ─────────────────────────────────────────────────────────────────────────────

fn demo_thick_solid_with_removed_faces() {
    separator("6. Thick Solid with Removed Faces");

    let brep = make_box_with_geom();

    // Remove one face
    {
        let result = thick_solid_with_removed_faces(&brep, &[5], 0.1);
        assert!(result.is_some(), "should succeed with one face removed");
        let r = result.unwrap();
        println!("  Remove top face, thickness 0.1:");
        println!("    Offset faces: {}", r.offset_faces);
        println!("    Lateral faces: {}", r.lateral_faces);
        println!("    Self-intersection: {}", r.self_intersection);
    }

    // Remove multiple faces
    {
        let result = thick_solid_with_removed_faces(&brep, &[0, 5], 0.1);
        assert!(result.is_some(), "should succeed with multiple faces removed");
        let r = result.unwrap();
        println!("  Remove top and bottom, thickness 0.1:");
        println!("    Offset faces: {}", r.offset_faces);
        println!("    Lateral faces: {}", r.lateral_faces);
    }

    // Self-intersection detection
    {
        let small_box = {
            let mut b = BRep::from_primitive(PrimitiveSolid::Box {
                width: 1.0,
                height: 1.0,
                depth: 1.0,
            });
            geom_populate::populate_box_geom(&mut b);
            b
        };

        // Thickness > half the minimum dimension should detect self-intersection
        let result = thick_solid_with_removed_faces(&small_box, &[], 0.6);
        assert!(result.is_some(), "should produce result even with self-intersection");
        let r = result.unwrap();
        println!("  Small box (1x1x1), thickness 0.6 (> 0.5):");
        println!("    Self-intersection detected: {}", r.self_intersection);
        assert!(r.self_intersection, "should detect self-intersection for large thickness");
    }

    println!("  PASS");
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

fn main() {
    println!("=================================================");
    println!("  Offset and Shell Operations Demo");
    println!("=================================================");

    demo_surface_offset();
    demo_offset_shell();
    demo_offset_solid();
    demo_hollow_solid();
    demo_thicken_shell();
    demo_thick_solid_with_removed_faces();

    println!("\n=================================================");
    println!("  Offset: All demos completed successfully");
    println!("=================================================");
}
