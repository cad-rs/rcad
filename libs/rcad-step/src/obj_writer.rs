//! OBJ mesh exporter for BRep solids.
//!
//! Writes the triangulated faces of a `BRep` as a Wavefront OBJ file.
//! Each face must already be triangulated (`Face.triangles` non-empty); faces
//! without triangles are skipped.
//!
//! Analogous to the mesh-export path in OCCT `RWMesh_FaceMeshComp`.

use std::io::{self, Write};

use rcad_kernel::BRep;

/// Write `brep` as Wavefront OBJ text to `writer`.
///
/// Vertex indices in the OBJ file are 1-based.  All triangles from all faces
/// of all solids are emitted.  Faces without pre-triangulated data are skipped.
///
/// Returns the number of triangles written.
pub fn write_obj(brep: &BRep, writer: &mut impl Write) -> io::Result<usize> {
    // Emit all unique vertices first.
    for v in &brep.vertices {
        writeln!(writer, "v {:.9} {:.9} {:.9}", v.point.x, v.point.y, v.point.z)?;
    }

    let mut total_tris = 0usize;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for &[i, j, k] in &face.triangles {
                    // OBJ is 1-based
                    writeln!(writer, "f {} {} {}", i + 1, j + 1, k + 1)?;
                    total_tris += 1;
                }
            }
        }
    }

    Ok(total_tris)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topology::{Face, Shell, Solid, Vertex, Wire};
    use glam::DVec3;

    /// Build a minimal triangulated BRep: one solid with a single face that has
    /// two triangles (a square split diagonally).
    fn make_triangulated_brep() -> BRep {
        let mut brep = BRep::new();
        // 4 vertices of a 1×1 square in XY plane
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 3

        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![[0, 1, 2], [0, 2, 3]],
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        brep
    }

    #[test]
    fn write_obj_produces_correct_output() {
        let brep = make_triangulated_brep();
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).expect("write_obj should succeed");
        let text = String::from_utf8(buf).unwrap();

        let f_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('f')).collect();
        let v_lines: Vec<&str> = text.lines().filter(|l| l.starts_with('v')).collect();

        assert_eq!(n, 2, "should return 2 triangles");
        assert_eq!(f_lines.len(), 2, "should have 2 'f' lines");
        assert_eq!(v_lines.len(), 4, "should have 4 'v' lines");

        // OBJ indices are 1-based
        assert!(text.contains("f 1 2 3"), "first triangle should be f 1 2 3");
        assert!(text.contains("f 1 3 4"), "second triangle should be f 1 3 4");
    }

    #[test]
    fn write_obj_empty_brep() {
        let brep = BRep::new();
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).expect("write_obj should handle empty BRep");
        assert_eq!(n, 0);
    }

    #[test]
    fn write_obj_face_without_triangles_is_skipped() {
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        let face = Face {
            outer_wire: Wire { edges: vec![] },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![], // no triangles
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let mut buf = Vec::new();
        let n = write_obj(&brep, &mut buf).unwrap();
        assert_eq!(n, 0, "face without triangles should produce 0 triangles");
    }
}
