// Quick test: verify GWedge produces correct 8V/12E/6F box
use rcad_kernel::topods::TShape;
use rcad_kernel::BRep;
use rcad_modeling::prim::prim::builder::PrimBuilder;
use rcad_modeling::prim::prim::gwedge::GWedge;
use std::sync::Arc;

#[test]
fn gwedge_box_counts() {
    let mut t = BRep::new();
    let builder = PrimBuilder::new(&mut t);
    let mut wedge = GWedge::new_box(builder, 1.0, 1.0, 1.0);
    let _shell = wedge.build_shell();
    let is_v = |ts: &Arc<TShape>| matches!(&**ts, TShape::Vertex(_));
    let is_e = |ts: &Arc<TShape>| matches!(&**ts, TShape::Edge(_));
    let is_w = |ts: &Arc<TShape>| matches!(&**ts, TShape::Wire(_));
    let is_f = |ts: &Arc<TShape>| matches!(&**ts, TShape::Face(_));
    let nv = t.tshapes.iter().filter(|ts| is_v(ts)).count();
    let ne = t.tshapes.iter().filter(|ts| is_e(ts)).count();
    let nw = t.tshapes.iter().filter(|ts| is_w(ts)).count();
    let nf = t.tshapes.iter().filter(|ts| is_f(ts)).count();
    eprintln!("GWedge box: V={} E={} W={} F={}", nv, ne, nw, nf);
    assert_eq!(nv, 8, "expected 8 vertices, got {}", nv);
    assert_eq!(ne, 12, "expected 12 edges, got {}", ne);
    assert_eq!(nw, 6, "expected 6 wires, got {}", nw);
    assert_eq!(nf, 6, "expected 6 faces, got {}", nf);
}

