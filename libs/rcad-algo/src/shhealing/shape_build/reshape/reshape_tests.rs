//! Unit tests for the ShapeBuild_ReShape port. The ValueLeaf/cycle cases are
//! ported from the OCCT gtest `BRepTools_ReShape_Test.cxx`; the wire/edge
//! rebuild cases pin the applyImpl semantics used by the ShapeFix stack.

use super::*;
use rcad_kernel::geom::{Curve3, Line3 as Line};
use rcad_kernel::topo::topods::{BRep, Orientation, Shape, ShapeType};
use std::sync::Arc;

fn make_vertex(brep: &mut BRep, x: f64, y: f64, z: f64) -> Shape {
    brep.add_tvertex(glam::DVec3::new(x, y, z))
}

fn make_edge(brep: &mut BRep, v1: &Shape, v2: &Shape) -> Shape {
    let p1 = v1.as_vertex().unwrap().point;
    let p2 = v2.as_vertex().unwrap().point;
    let dir = (p2 - p1).normalize();
    let curve = Curve3::Line(Line {
        origin: p1,
        direction: dir,
    });
    let len = (p2 - p1).length();
    brep.add_tedge(
        Some(curve),
        Shape {
            orientation: Orientation::Forward,
            ..v1.clone()
        },
        Shape {
            orientation: Orientation::Reversed,
            ..v2.clone()
        },
        [0.0, len],
    )
}

// OCCT gtest BRepTools_ReShapeTest.ValueLeaf_FollowsChainToLeaf: Value()
// returns the direct replacement only (one hop); ValueLeaf() follows the
// chain through intermediate replacements to its terminal shape.
#[test]
fn value_leaf_follows_chain_to_leaf() {
    let mut brep = BRep::new();
    let a = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let b = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let c = make_vertex(&mut brep, 2.0, 0.0, 0.0);

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &a, &b);
    rs.replace(&mut brep, &b, &c);

    assert!(
        rs.value(&mut brep, &a).is_same(&b),
        "Value() must yield only the direct replacement"
    );
    assert!(
        rs.value_leaf(&mut brep, &a).is_same(&c),
        "ValueLeaf() must walk the full chain"
    );
    assert!(rs.value_leaf(&mut brep, &b).is_same(&c));
    assert!(
        rs.value_leaf(&mut brep, &c).is_same(&c),
        "Terminal element is a fixpoint"
    );
}

// OCCT gtest BRepTools_ReShapeTest.ValueLeaf_UnrecordedShapeIsIdentity.
#[test]
fn value_leaf_unrecorded_shape_is_identity() {
    let mut brep = BRep::new();
    let a = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let rs = ShapeBuildReShape::new();
    assert!(rs.value_leaf(&mut brep, &a).is_same(&a));
}

// OCCT gtest BRepTools_ReShapeTest.ValueLeaf_ChainEndingInRemoveReturnsNull:
// a chain ending in a removal surfaces as a null shape.
#[test]
fn value_leaf_chain_ending_in_remove_returns_null() {
    let mut brep = BRep::new();
    let a = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let b = make_vertex(&mut brep, 1.0, 0.0, 0.0);

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &a, &b);
    rs.remove(&mut brep, &b);

    assert!(rs.value_leaf(&mut brep, &a).index == usize::MAX);
}

// OCCT gtest BRepTools_ReShapeTest.Replace_DirectCycleHandledByApply: A->B->A
// bindings are both recorded; the DFS in-flight guard breaks the cycle.
#[test]
fn replace_direct_cycle_handled_by_apply() {
    let mut brep = BRep::new();
    let a = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let b = make_vertex(&mut brep, 1.0, 0.0, 0.0);

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &a, &b);
    rs.replace(&mut brep, &b, &a);

    assert!(rs.value(&mut brep, &a).is_same(&b));
    assert!(rs.value(&mut brep, &b).is_same(&a));

    let a_result = rs.apply(&mut brep, &a, ShapeType::Shape);
    assert!(
        a_result.index != usize::MAX,
        "Apply must terminate via the DFS guard"
    );
}

// OCCT gtest BRepTools_ReShapeTest.Replace_LongerCycleHandledByApply.
#[test]
fn replace_longer_cycle_handled_by_apply() {
    let mut brep = BRep::new();
    let a = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let b = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let c = make_vertex(&mut brep, 2.0, 0.0, 0.0);

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &a, &b);
    rs.replace(&mut brep, &b, &c);
    rs.replace(&mut brep, &c, &a);

    assert!(rs.value(&mut brep, &a).is_same(&b));
    assert!(rs.value(&mut brep, &b).is_same(&c));
    assert!(rs.value(&mut brep, &c).is_same(&a));

    let a_result = rs.apply(&mut brep, &a, ShapeType::Shape);
    assert!(
        a_result.index != usize::MAX,
        "Apply must terminate via the DFS guard"
    );
}

// applyImpl rebuilds a wire whose sub-edge was replaced: the result is a NEW
// wire TShape carrying the replacement edge, status DONE3 set.
#[test]
fn apply_rebuilds_wire_with_replaced_edge() {
    let mut brep = BRep::new();
    let v1 = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let v2 = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let v3 = make_vertex(&mut brep, 1.0, 1.0, 0.0);
    let e1 = make_edge(&mut brep, &v1, &v2);
    let e2 = make_edge(&mut brep, &v2, &v3);
    let wire = brep.add_twire(vec![e1.clone(), e2.clone()]);

    // Replacement edge: same geometry, distinct TShape.
    let e1_new = make_edge(&mut brep, &v1, &v2);
    assert!(!e1_new.is_same(&e1));

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &e1, &e1_new);

    let res = rs.apply(&mut brep, &wire, ShapeType::Shape);
    assert!(res.shape_type() == ShapeType::Wire);
    assert!(!res.is_same(&wire), "the wire must be rebuilt");
    let edges = res.as_wire().unwrap().edges.clone();
    assert_eq!(edges.len(), 2);
    assert!(edges[0].is_same(&e1_new), "edge 1 replaced");
    assert!(edges[1].is_same(&e2), "edge 2 untouched");
    assert!(
        rs.status_flag(ShapeExtendStatus::Done3),
        "DONE3: some subshapes replaced"
    );
    assert!(!rs.status_flag(ShapeExtendStatus::Done1));
}

// applyImpl on a wire with a removed edge yields an (empty) rebuilt wire and
// DONE4; the original edge keeps its replacement-recorded state.
#[test]
fn apply_wire_with_removed_edge() {
    let mut brep = BRep::new();
    let v1 = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let v2 = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let e1 = make_edge(&mut brep, &v1, &v2);
    let wire = brep.add_twire(vec![e1.clone()]);

    let mut rs = ShapeBuildReShape::new();
    rs.remove(&mut brep, &e1);

    let res = rs.apply(&mut brep, &wire, ShapeType::Shape);
    assert!(res.shape_type() == ShapeType::Wire);
    assert_eq!(
        res.as_wire().unwrap().edges.len(),
        0,
        "the removed edge must be gone"
    );
    assert!(rs.status_flag(ShapeExtendStatus::Done4));
}

// applyImpl on an edge whose vertex was replaced keeps the curve, range and
// pcurves (CopyRanges + EmptyCopy representation copy), and re-adds the new
// vertex into the right slot.
#[test]
fn apply_edge_rebuild_preserves_geometry_and_pcurves() {
    let mut brep = BRep::new();
    let v1 = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let v2 = make_vertex(&mut brep, 2.0, 0.0, 0.0);
    let e = make_edge(&mut brep, &v1, &v2);

    // A pcurve on a face-like key (face ptr placeholder).
    let face_key = (0xdead_beefu64, 0u32);
    {
        use rcad_kernel::geom::{Curve2d, Line2d};
        let pc = Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        });
        brep.edge_mut_inplace(e.clone())
            .pcurves
            .insert(face_key, (pc, 0.0, 2.0));
    }
    let old_curve = e.as_edge().unwrap().curve.clone();

    // Replacement vertex: distinct TShape at the same position. rcad's
    // add_tvertex shares one TShape per position (identity cache), so the
    // copy comes from EmptyCopied (OCCT BRep_TVertex::EmptyCopy: Pnt +
    // Tolerance).
    let v1_new = brep.empty_copied(&v1);
    assert!(!v1_new.is_same(&v1));

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &v1, &v1_new);

    let res = rs.apply(&mut brep, &e, ShapeType::Edge);
    // until = EDGE: OCCT rank of EDGE (6) >= rank of EDGE (6) -> stop at the
    // Value() result... the edge itself is not replaced, so the rebuild must
    // come from Apply with until = SHAPE.
    let _ = res;
    let res = rs.apply(&mut brep, &e, ShapeType::Shape);
    assert!(res.shape_type() == ShapeType::Edge);
    assert!(!res.is_same(&e), "the edge must be rebuilt");
    let ed = res.as_edge().unwrap();
    assert!(ed.first.is_same(&v1_new), "new vertex in the first slot");
    assert!(ed.last.is_same(&v2), "last vertex untouched");
    assert!(ed.curve.is_some());
    assert_eq!(ed.range, [0.0, 2.0], "CopyRanges restores the 3D range");
    assert!(
        ed.pcurves.contains_key(&face_key),
        "EmptyCopy keeps pcurve representations"
    );
    assert_eq!(
        ed.pcurves.get(&face_key).unwrap().1..ed.pcurves.get(&face_key).unwrap().2,
        0.0..2.0,
        "CopyRanges restores the pcurve range"
    );
    let _ = old_curve;
}

// BRep_Tool::IsClosed on wires: a closed square wire is closed, an open
// two-edge wire is not.
#[test]
fn brep_tool_is_closed_wire() {
    let mut brep = BRep::new();
    let v1 = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let v2 = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let v3 = make_vertex(&mut brep, 1.0, 1.0, 0.0);
    let e1 = make_edge(&mut brep, &v1, &v2);
    let e2 = make_edge(&mut brep, &v2, &v3);
    let open_wire = brep.add_twire(vec![e1.clone(), e2.clone()]);
    assert!(!brep_tool_is_closed(&brep, &open_wire));

    let v4 = make_vertex(&mut brep, 0.0, 1.0, 0.0);
    let e3 = make_edge(&mut brep, &v3, &v4);
    let e4 = make_edge(&mut brep, &v4, &v1);
    let closed_wire = brep.add_twire(vec![e1, e2.clone(), e3, e4]);
    assert!(brep_tool_is_closed(&brep, &closed_wire));
    let _ = e2;
}

// The replacement map ignores orientation (TopTools_ShapeMapHasher):
// recording e1 Forward and querying the REVERSED wrapper resolves to the
// replacement, flipped back.
#[test]
fn replace_is_orientation_insensitive() {
    let mut brep = BRep::new();
    let v1 = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    let v2 = make_vertex(&mut brep, 1.0, 0.0, 0.0);
    let e1 = make_edge(&mut brep, &v1, &v2);
    let e1_new = make_edge(&mut brep, &v1, &v2);

    let mut rs = ShapeBuildReShape::new();
    rs.replace(&mut brep, &e1, &e1_new);
    assert!(rs.is_recorded(&e1));

    let mut e1_rev = e1.clone();
    e1_rev.orientation = Orientation::Reversed;
    let v = rs.value(&mut brep, &e1_rev);
    assert!(v.is_same(&e1_new));
    assert!(
        v.orientation == Orientation::Reversed,
        "the result is flipped back to the queried orientation"
    );
}

// The Arc identity of pool slots survives in-place container mutation: the
// handle returned by apply() and the pool slot share one TShape allocation.
#[test]
fn builder_add_preserves_arc_identity() {
    let mut brep = BRep::new();
    let c = brep.add_tcompound(Vec::new());
    let v = make_vertex(&mut brep, 0.0, 0.0, 0.0);
    builder_add(&mut brep, &c, &v);
    // The pool slot and the handle must be the same allocation.
    assert_eq!(
        Arc::as_ptr(&brep.tshapes[c.index]) as u64,
        c.ptr_id(),
        "Arc::make_mut must not have split the identity"
    );
    let children = raw_children_count(&brep, &c);
    assert_eq!(children, 1);
}

fn raw_children_count(brep: &BRep, s: &Shape) -> usize {
    crate::shhealing::shape_build::brep_tool::raw_subshapes(brep, s).len()
}
