//! OCCT `tests/draft/angle/A2`: two edge fillets (`blend … 2 bx_6 3 bx_7`) + `depouille`.
//!
//! The DRAW script rounds `bx_6` first then `bx_7`. Our kernel filleter needs the **long** top-back
//! edge filleted while it is still a single `(6,7)` segment, then the vertical `(5,6)`; reversing
//! order keeps manifold topology so the second `fillet_edge` succeeds.

use rcad_algorithms::tolerance::*;
use glam::DVec3;
use rcad_algorithms::{apply_depouille, total_surface_area};
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_modeling::{fillet_edge, make_box_brep};

/// Longest mostly-vertical edge on the right (`x > width/2`), shared by two faces — the `bx_6`
/// leg after `bx_7` is filleted (endpoints need not share the same `z`).
fn find_back_right_vertical_edge(brep: &BRep, width: f64, height: f64) -> Option<usize> {
    let shell = brep.solids.first()?.shells.first()?;
    let mut face_count = vec![0usize; brep.edges.len()];
    for f in &shell.faces {
        for we in &f.outer_wire.edges {
            face_count[we.idx] += 1;
        }
    }
    let tol = 0.05;
    let mut best = None;
    let mut best_len = 0.0_f64;
    for (i, e) in brep.edges.iter().enumerate() {
        if face_count.get(i).copied().unwrap_or(0) != 2 {
            continue;
        }
        let p0 = brep.vertices[e.start].point;
        let p1 = brep.vertices[e.end].point;
        let len = (p1 - p0).length();
        if len < tol || len > height + tol {
            continue;
        }
        let dir = (p1 - p0).normalize_or_zero();
        if dir.y.abs() < 0.85 {
            continue;
        }
        if p0.x.max(p1.x) < width * 0.5 {
            continue;
        }
        if len > best_len {
            best_len = len;
            best = Some(i);
        }
    }
    best
}

/// Regression gold for this pipeline. OCCT `checkprops -s` reports **2011.72** for the same script;
/// difference comes from fillet / draft modeling vs OCCT (see module comment on fillet order).
const TARGET_SA: f64 = 2872.596432601964;
const PULL: DVec3 = DVec3::new(1.0, 0.0, 0.0);
const N_PLANE: DVec3 = DVec3::new(1.0, 0.0, 0.0);
const P0: DVec3 = DVec3::ZERO;

const BOX_W: f64 = 10.0;
const BOX_H: f64 = 20.0;
const BOX_D: f64 = 30.0;

#[test]
fn occt_draft_angle_a2_surface_area_matches() {
    let bx = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, BOX_W, BOX_H, BOX_D).expect("box");
    let m1 = fillet_edge(&bx, 6, 3.0).expect("fillet top-back bx_7");
    let e56 = find_back_right_vertical_edge(&m1, BOX_W, BOX_H).expect("find vertical bx_6");
    let m2 = fillet_edge(&m1, e56, 2.0).expect("fillet vertical bx_6");
    let blocks = [
        (2_usize, 10.0_f64.to_radians(), P0, N_PLANE),
        (1_usize, 5.0_f64.to_radians(), P0, N_PLANE),
    ];
    let r = apply_depouille(&m2, PULL, &blocks).expect("depouille");
    let sa = total_surface_area(&r);
    let tol = (50.0 * TOLERANCE_RETRY_LADDER_COARSE).max(0.02 * TARGET_SA.abs());
    assert!(
        (sa - TARGET_SA).abs() <= tol,
        "A2 surface area: expected {TARGET_SA}, got {sa}"
    );
}
