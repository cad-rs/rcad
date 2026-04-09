/// Smoke tests for rcad-render's CPU tessellation path.
/// These do NOT require a GPU — they only test the mesh-building logic.
use rcad_kernel::BRep;
use rcad_kernel::PrimitiveSolid;
use rcad_render::Tessellator;

fn make_box_brep() -> BRep {
    BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    })
}

/// Tessellating an empty BRep should return an empty mesh without panicking.
#[test]
fn tessellate_empty_brep_no_panic() {
    let empty = BRep::default();
    let mesh = Tessellator::tessellate(&empty);
    assert!(mesh.vertices.is_empty(), "empty BRep should yield no vertices");
    assert!(mesh.indices.is_empty(), "empty BRep should yield no indices");
}

/// Tessellating a box primitive (no geometry populated) should produce vertices.
#[test]
fn tessellate_box_has_vertices() {
    let brep = make_box_brep();
    let mesh = Tessellator::tessellate(&brep);
    // The box primitive has 8 corner vertices
    assert_eq!(mesh.vertices.len(), 8, "unit box should have 8 vertices");
}

/// Tessellating a box with triangles populated should produce triangle indices.
#[test]
fn tessellate_box_with_triangles_has_indices() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    let mesh = Tessellator::tessellate(&brep);
    assert!(!mesh.vertices.is_empty(), "should have vertices");
    // Triangle indices come in sets of 3
    assert!(mesh.indices.len() % 3 == 0, "index count must be divisible by 3");
}

/// All triangle indices must be within bounds of the vertex buffer.
#[test]
fn tessellate_box_indices_in_bounds() {
    use rcad_algorithms::geom_populate;

    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: 1.0,
        height: 1.0,
        depth: 1.0,
    });
    geom_populate::populate_box_geom(&mut brep);

    let mesh = Tessellator::tessellate(&brep);
    let nv = mesh.vertices.len() as u32;
    for &idx in &mesh.indices {
        assert!(idx < nv, "triangle index {idx} out of bounds (nv={nv})");
    }
}
