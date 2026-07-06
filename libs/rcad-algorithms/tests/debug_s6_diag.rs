use rcad_algorithms::{boolean_op, BooleanOpType};
use rcad_kernel::{topods, surface_area, BRep};
use rcad_kernel::face_flat_iter;
use rcad_modeling::{make_box_brep, make_cylinder_brep};
use glam::DVec3;

#[test]
fn debug_s6_diag() {
    let cyl = make_cylinder_brep(DVec3::new(0.0, 0.0, 1.0), DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
    let bx = make_box_brep(DVec3::new(-0.8, -0.8, 0.0), DVec3::X, DVec3::Y, 1.8, 1.8, 2.0).unwrap();

    println!("Cylinder SA: {:.6} (expected: 6π ≈ 18.85)", surface_area(&cyl));
    println!("Box SA: {:.6}", surface_area(&bx));
    println!("Expected Result: ~10.602");

    let result = boolean_op(BooleanOpType::Difference, &cyl, &bx).unwrap();
    println!("\nResult SA: {:.6}", surface_area(&result));

    // Per-face via face_flat_iter
    for (fi, (face_idx, face)) in face_flat_iter(&result).enumerate() {
        let n_edges = face.outer_wire.edges.len();
        let inner = face.inner_wires.len();
        // Compute area via try_analytic_face_surface_area or fallback
        let area = estimate_face_area(&result, face_idx);
        println!("  Face {fi} ({n_edges} edges, {inner} iw): area ≈ {area:.6}");
    }

    // Also count faces per solid/shell
    for (si, sol) in result.solids.iter().enumerate() {
        for (shi, shell) in sol.shells.iter().enumerate() {
            println!("  Solid {si} Shell {shi}: {} faces",
                shell.faces.len());
        }
    }
}

fn estimate_face_area(brep: &BRep, fi: usize) -> f64 {
    // Quick estimate from triangles if available
    use rcad_kernel::Face;
    for sol in &brep.solids {
        let mut idx = 0;
        for shell in &sol.shells {
            for face in &shell.faces {
                if idx == fi {
                    if !face.triangles.is_empty() {
                        let mut total = 0.0;
                        for tri in &face.triangles {
                            let v0 = brep.vertices[tri[0]].point;
                            let v1 = brep.vertices[tri[1]].point;
                            let v2 = brep.vertices[tri[2]].point;
                            total += (v1 - v0).cross(v2 - v0).length() * 0.5;
                        }
                        return total;
                    }
                    return 0.0;
                }
                idx += 1;
            }
        }
    }
    0.0
}
