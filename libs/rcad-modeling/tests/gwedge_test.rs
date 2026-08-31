// Quick test: verify GWedge produces correct 8V/12E/6F box
use rcad_kernel::topods::{self, TShape};
use rcad_kernel::BRep;
use rcad_modeling::prim::prim::builder::PrimBuilder;
use rcad_modeling::prim::prim::gwedge::GWedge;
use std::sync::Arc;

/// Count shapes in BRep directly (scanning tshapes).
fn count_shapes(t: &BRep) -> (usize, usize, usize, usize, usize) {
    let mut nv=0; let mut ne=0; let mut nw=0; let mut nf=0; let mut ns=0;
    for ts in &t.tshapes {
        match &**ts {
            TShape::Vertex(_) => nv+=1,
            TShape::Edge(_) => ne+=1,
            TShape::Wire(_) => nw+=1,
            TShape::Face(_) => nf+=1,
            TShape::Shell(_) => ns+=1,
            _ => {}
        }
    }
    (nv, ne, nw, nf, ns)
}

/// Recursively collect shapes reachable from root via sub_shapes_of.
fn reachable(t: &BRep, shape: &topods::Shape) -> Vec<usize> {
    let mut seen = vec![shape.index];
    for child in sub_shapes_of(t, shape) {
        for s in reachable(t, &child) {
            if !seen.contains(&s) { seen.push(s); }
        }
    }
    seen
}

fn sub_shapes_of(t: &BRep, s: &topods::Shape) -> Vec<topods::Shape> {
    let idx = s.index;
    if idx >= t.tshapes.len() { return vec![]; }
    let ts = &t.tshapes[idx];
    let cp = |sr: &topods::Shape| topods::Shape {
        data: sr.data.clone(), index: sr.index,
        orientation: sr.orientation, location: sr.location,
    };
    match &**ts {
        TShape::Vertex(_) => vec![],
        TShape::Edge(ed) => vec![cp(&ed.first), cp(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(cp).collect(),
        TShape::Face(fd) => {
            let mut v = vec![cp(&fd.outer_wire)];
            v.extend(fd.inner_wires.iter().map(cp));
            v
        }
        TShape::Shell(sd) => sd.faces.iter().map(cp).collect(),
        TShape::Solid(sd) => sd.shells.iter().map(cp).collect(),
        _ => vec![],
    }
}

#[test]
fn gwedge_box_counts() {
    let mut t = BRep::new();
    let builder = PrimBuilder::new(&mut t);
    let mut wedge = GWedge::new_box(builder, 1.0, 1.0, 1.0);
    let shell = wedge.build_shell();
    t.add_tsolid(vec![shell]);

    let (nv, ne, nw, nf, ns) = count_shapes(&t);
    eprintln!("GWedge box: V={} E={} W={} F={} S={}", nv, ne, nw, nf, ns);
    assert_eq!(nv, 8, "expected 8 vertices, got {}", nv);
    assert_eq!(ne, 12, "expected 12 edges, got {}", ne);
    assert_eq!(nw, 6, "expected 6 wires, got {}", nw);
    assert_eq!(nf, 6, "expected 6 faces, got {}", nf);

    // Check shell face count
    let last_sh = t.tshapes.iter().enumerate().rev().find(|(_,ts)| matches!(ts.as_ref(), TShape::Shell(_)));
    if let Some((si, _)) = last_sh {
        if let TShape::Shell(sd) = &*t.tshapes[si] {
            eprintln!("  shell[{}] faces={}", si, sd.faces.len());
            for (fi, f) in sd.faces.iter().enumerate() {
                eprintln!("    face {}: outer_wire idx={}", fi, f.index);
            }
        }
    }

    // Check topology reachability from the solid
    let solid = t.tshapes.len() - 1;
    let root = topods::Shape { data: t.tshapes[solid].clone(), index: solid,
        orientation: topods::Orientation::Forward, location: 0 };
    let reachable_ids = reachable(&t, &root);
    let rv = reachable_ids.iter().filter(|&&i| matches!(&*t.tshapes[i], TShape::Vertex(_))).count();
    let re = reachable_ids.iter().filter(|&&i| matches!(&*t.tshapes[i], TShape::Edge(_))).count();
    eprintln!("  reachable: V={} E={} (total {} shapes)", rv, re, reachable_ids.len());
    assert_eq!(rv, 8, "reachable vertices");
    assert_eq!(re, 12, "reachable edges");
}

