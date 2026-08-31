//! Unit tests for the ShapeExtend_WireData port (against OCCT semantics).

use super::WireData;
use glam::DVec3;
use rcad_kernel::topo::topods::{BRep, Orientation, Shape};

/// Build a BRep pool with a triangle wire (-1,0)->(9,0)->(1,1)->closed,
/// mirroring the heal `wire_tails_composed` A1 construction:
///   vertex v1 -1 0 0 ; vertex v2 9 0 0 ; vertex v3 1 1 0
///   polyvertex w v1 v2 v3 v1 ; mkplane s w
fn triangle_wire() -> (BRep, Shape) {
    let mut brep = BRep::new();
    let pts = [DVec3::new(-1.0, 0.0, 0.0), DVec3::new(9.0, 0.0, 0.0), DVec3::new(1.0, 1.0, 0.0)];
    let v: Vec<Shape> = pts.iter().map(|p| brep.add_tvertex(*p)).collect();
    let mut edges: Vec<Shape> = Vec::new();
    for i in 0..3 {
        let a = pts[i];
        let b = pts[(i + 1) % 3];
        let mut v1 = v[i].clone();
        v1.orientation = Orientation::Forward;
        let mut v2 = v[(i + 1) % 3].clone();
        v2.orientation = Orientation::Reversed;
        let curve = rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(a, b - a));
        edges.push(brep.add_tedge(Some(curve), v1, v2, [0.0, (b - a).length()]));
    }
    let wire = brep.add_twire(edges);
    (brep, wire)
}

#[test]
fn init_from_wire_and_basic_access() {
    let (mut brep, wire) = triangle_wire();
    let mut wd = WireData::new();
    assert!(wd.init(&wire, true, true));
    assert_eq!(wd.nb_edges(), 3);
    // Edge 1 endpoint order.
    let e1 = wd.edge(1);
    let ed = WireData::edge_data(&e1).unwrap();
    assert_eq!(
        match ed.first.data.as_ref() {
            rcad_kernel::topo::topods::TShape::Vertex(v) => v.point,
            _ => panic!("vertex expected"),
        },
        DVec3::new(-1.0, 0.0, 0.0)
    );
    // Negative rank returns the reversed edge.
    let e1r = wd.edge(-1);
    assert_eq!(e1r.orientation, Orientation::Reversed);
    // Round-trip through wire().
    let w2 = wd.wire(&mut brep);
    assert!(matches!(w2.data.as_ref(), rcad_kernel::topo::topods::TShape::Wire(_)));
}

#[test]
fn reverse_and_remove_semantics() {
    let (mut brep, wire) = triangle_wire();
    let mut wd = WireData::new();
    wd.init(&wire, true, true);

    let e0 = wd.edge(1);
    wd.reverse();
    assert_eq!(wd.nb_edges(), 3);
    // After Reverse the last edge is the reversed former first edge.
    let elast = wd.edge(wd.nb_edges());
    assert_eq!(elast.ptr_id(), e0.ptr_id());
    assert_eq!(elast.orientation, Orientation::Reversed);

    // Remove(rank) drops the rank; Remove(0) drops the last.
    wd.remove(0);
    assert_eq!(wd.nb_edges(), 2);
    wd.remove(1);
    assert_eq!(wd.nb_edges(), 1);
    let _ = &mut brep;
}

#[test]
fn seam_detection() {
    let (mut brep, wire) = triangle_wire();
    let mut wd = WireData::new();
    wd.init(&wire, true, true);
    // Duplicate the first edge REVERSED: same TShape, opposite orientation
    // -> a seam pair per OCCT ComputeSeams.
    let mut e1 = wd.edge(1);
    e1.orientation = Orientation::Reversed;
    wd.add_edge(&e1, 0);
    assert!(wd.is_seam(1));
    assert!(wd.is_seam(4));
    assert!(!wd.is_seam(2));
    let _ = &mut brep;
}
