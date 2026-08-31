//! Unit tests for the TKHelix port (HelixGeom approximation + HelixBRep
//! wire building), asserted against OCCT DRAWEXE nbshapes references
//! (helix/standard grid, default nbshapes semantics).

use super::commands;
use rcad_kernel::topo::topods::TShape;

/// Count (vertices, edges) reachable from the wires of a built helix BRep —
/// OCCT nbshapes walks the shape tree counting shared TShapes once.
fn count_wire_ve(brep: &rcad_kernel::topo::topods::BRep) -> (usize, usize, usize) {
    let mut verts = std::collections::BTreeSet::new();
    let mut edges = std::collections::BTreeSet::new();
    let mut wires = std::collections::BTreeSet::new();
    for ts in &brep.tshapes {
        if let TShape::Wire(wd) = ts.as_ref() {
            wires.insert(std::sync::Arc::as_ptr(ts) as u64);
            for e in &wd.edges {
                edges.insert(std::sync::Arc::as_ptr(&e.data) as u64);
                if let TShape::Edge(ed) = e.data.as_ref() {
                    verts.insert(std::sync::Arc::as_ptr(&ed.first.data) as u64);
                    verts.insert(std::sync::Arc::as_ptr(&ed.last.data) as u64);
                }
            }
        }
    }
    (verts.len(), edges.len(), wires.len())
}

#[test]
fn helix_standard_a1_cylindrical_wire_topology() {
    // OCCT: helix result 1 100 100 5 0 -> V=6, E=5, W=1 (5 full turns).
    let brep = commands::helix(100.0, &[100.0], &[5.0], &[false]).expect("helix A1");
    let (v, e, w) = count_wire_ve(&brep);
    assert_eq!((v, e, w), (6, 5, 1));
}

#[test]
fn helix_standard_b1_composite_one_part_wire_topology() {
    // OCCT: comphelix result 1 100 100 100 20 1 -> V=6, E=5, W=1.
    let brep = commands::comphelix(&[100.0, 100.0], &[100.0], &[20.0], &[true])
        .expect("comphelix B1");
    let (v, e, w) = count_wire_ve(&brep);
    assert_eq!((v, e, w), (6, 5, 1));
}

#[test]
fn helix_standard_c1_spiral_three_parts_wire_topology() {
    // OCCT: spiral result 3 100 20 20 60 20 2 6 2 0 0 0 -> V=11, E=10, W=1
    // (three parts joined by coincident vertices).
    let brep = commands::spiral(
        100.0,
        20.0,
        &[20.0, 60.0, 20.0],
        &[2.0, 6.0, 2.0],
        &[false, false, false],
    )
    .expect("spiral C1");
    let (v, e, w) = count_wire_ve(&brep);
    assert_eq!((v, e, w), (11, 10, 1));
}

#[test]
fn helix_standard_g1_comphelix2_topology() {
    // OCCT: comphelix2 result 3 100 84 36 20 10 1 10 2 60 2
    // -> V=65, E=64, W=1.
    let brep = commands::comphelix2(
        &[100.0, 84.0, 36.0, 20.0],
        &[10.0, 1.0, 10.0],
        &[2.0, 60.0, 2.0],
    )
    .expect("comphelix2 G1");
    let (v, e, w) = count_wire_ve(&brep);
    assert_eq!((v, e, w), (65, 64, 1));
}

/// OCCT BSplCLib::Interpolate regression: a clamped cubic interpolates f(t)=t
/// at the Schoenberg points exactly.
#[test]
fn bspl_interpolate_linear() {
    use rcad_kernel::math::bspl_lib::{build_schoenberg_points, interpolate};
    let degree = 3usize;
    let flat = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let mut params = vec![0.0; 4];
    build_schoenberg_points(degree, &flat, &mut params);
    let contact = vec![0i32; 4];
    let mut poles = params.clone(); // RHS = f(params) for f(t) = t.
    let err = interpolate(degree, &flat, &params, &contact, 1, &mut poles);
    assert_eq!(err, 0, "poles = {poles:?}, params = {params:?}");
    for (i, p) in poles.iter().enumerate() {
        assert!((p - params[i]).abs() < 1e-12, "poles = {poles:?}");
    }
}
