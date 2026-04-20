//! Dump edge valence after boolean union (two overlapping boxes, same setup as `boolean_ops`).
//!
//! Run: `cargo run --example boolean_union_open_edges -p rcad-examples`
//!
//! Interpreting output:
//! - In a **closed** manifold shell, every edge index should appear exactly **twice** (once per
//!   adjacent face, possibly outer+inner wires).
//! - `ref_count == 1` → boundary / missing partner face (typical “open box” symptom).
//! - `ref_count > 2` → non-manifold or duplicated seam.
//!
//! With the seam pass in `ResultBuilder::build` (union only), edge valence should be **2** for a
//! closed shell on this scenario. If you still see `ref_count != 2`, check T-junction handling in
//! `builder.rs` (`subdivide_edges_at_interior_vertices`).

use glam::DVec3;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bvh::Bvh;
use rcad_algorithms::builder::BooleanBuilder;
use rcad_algorithms::geom_populate::{populate_box_geom, recompute_plane_surfaces};
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::{BooleanOpType, boolean_op};
use rcad_kernel::BRep;
use rcad_kernel::PrimitiveSolid;
use std::collections::HashMap;

fn make_box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
    let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
        width: w,
        height: h,
        depth: d,
    });
    for v in &mut brep.vertices {
        v.point += DVec3::new(x, y, z);
    }
    populate_box_geom(&mut brep);
    brep
}

/// Count how many times each edge index appears in all wires of all faces in all shells.
fn edge_reference_counts(brep: &BRep) -> HashMap<usize, usize> {
    let mut m: HashMap<usize, usize> = HashMap::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    *m.entry(we.idx).or_insert(0) += 1;
                }
                for w in &face.inner_wires {
                    for we in &w.edges {
                        *m.entry(we.idx).or_insert(0) += 1;
                    }
                }
            }
        }
    }
    m
}

/// List (shell_face_flat_idx, outer|inner_wire) for one edge.
fn edge_face_mentions(brep: &BRep, edge_idx: usize) -> Vec<(usize, &'static str)> {
    let mut out = Vec::new();
    let mut flat = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let mut hit = |label: &'static str| {
                    out.push((flat, label));
                };
                for we in &face.outer_wire.edges {
                    if we.idx == edge_idx {
                        hit("outer");
                    }
                }
                for w in &face.inner_wires {
                    for we in &w.edges {
                        if we.idx == edge_idx {
                            hit("inner");
                        }
                    }
                }
                flat += 1;
            }
        }
    }
    out
}

fn edge_endpoints(brep: &BRep, edge_idx: usize) -> Option<(DVec3, DVec3)> {
    let e = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(e.start)?.point;
    let p1 = brep.vertices.get(e.end)?.point;
    Some((p0, p1))
}

fn print_report(brep: &BRep, label: &str) {
    let counts = edge_reference_counts(brep);
    let mut bad: Vec<(usize, usize)> = counts
        .iter()
        .filter(|(_, c)| **c != 2)
        .map(|(&i, &c)| (i, c))
        .collect();
    bad.sort_by_key(|x| x.0);

    println!("=== {label} ===");
    println!("vertices: {}  edges: {}  faces: {}", brep.vertices.len(), brep.edges.len(), {
        let mut n = 0usize;
        for s in &brep.solids {
            for sh in &s.shells {
                n += sh.faces.len();
            }
        }
        n
    });
    println!(
        "edges with ref_count != 2: {} (closed shell expects 0)",
        bad.len()
    );

    for (ei, c) in &bad {
        let (mid, len) = edge_endpoints(brep, *ei)
            .map(|(a, b)| (0.5 * (a + b), (b - a).length()))
            .unwrap_or((DVec3::ZERO, f64::NAN));
        let mentions = edge_face_mentions(brep, *ei);
        println!(
            "  edge {:4}  ref_count={}  len={:.6}  mid=({:.4},{:.4},{:.4})  faces={:?}",
            ei, c, len, mid.x, mid.y, mid.z, mentions
        );
    }
    println!();
}

fn main() {
    let a = make_box_at(0.0, 0.0, 0.0, 3.0, 2.0, 2.0);
    let b = make_box_at(1.5, 0.5, 0.5, 3.0, 1.0, 1.0);

    // 1) Same pipeline as `boolean_op` up to `BooleanBuilder::build()` (no recompute_plane_surfaces).
    {
        let mut ds = DS::new(&a, &b);
        let ba = Bvh::build(&a);
        let bb = Bvh::build(&b);
        let mut filler = PaveFiller::with_bvh(&mut ds, &ba, &bb);
        filler.perform();
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let brep = builder.build().expect("builder.build");
        print_report(&brep, "After BooleanBuilder::build() only");
    }

    // 2) After `recompute_plane_surfaces` (still what `boolean_op` does before returning).
    {
        let mut ds = DS::new(&a, &b);
        let ba = Bvh::build(&a);
        let bb = Bvh::build(&b);
        let mut filler = PaveFiller::with_bvh(&mut ds, &ba, &bb);
        filler.perform();
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let mut brep = builder.build().expect("builder.build");
        recompute_plane_surfaces(&mut brep);
        print_report(&brep, "After build() + recompute_plane_surfaces (matches boolean_op output)");
    }

    // 3) Public API `boolean_op` (includes recompute_plane_surfaces).
    {
        let brep = boolean_op(BooleanOpType::Union, &a, &b).expect("boolean_op");
        print_report(&brep, "boolean_op Union (public API)");
    }
}
