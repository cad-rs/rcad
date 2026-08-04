//! Smoke test for the Builder stage-by-stage API.
//!
//! Verifies the snapshot contract (stage sequence, field population) without
//! asserting boolean correctness — the rcad-algo Builder is mid-alignment, so
//! the pipeline may legitimately stop early. This file is hand-written and NOT
//! overwritten by tools/gen_builder_tests.py.

use glam::DVec3;
use rcad_algo::bop::algo::builder::{Builder, BooleanOpType};
use rcad_algo::bop::algo::pave_filler::PaveFiller;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, TShape};
use rcad_modeling::prim::primapi::make_box_brep;

/// Expected stage sequence (OCCT BOPAlgo_BOP::PerformInternal1 DUMP_STAGE names).
const STAGE_ORDER: &[&str] = &[
    "after_FillImagesVertices",
    "after_FillImagesEdges",
    "after_BuildResultWire",
    "after_FillImagesFaces",
    "after_BuildResultShell",
    "after_FillImagesSolids",
    "after_BuildResultCompSolid",
    "after_FillImagesCompounds",
    "after_PrepareHistory",
    "after_PostTreat",
];

fn root_shape(brep: &topods::BRep, location: u32) -> Shape {
    for (i, ts) in brep.tshapes.iter().enumerate().rev() {
        match &**ts {
            TShape::Solid(_) | TShape::Shell(_) => {
                return Shape::from_parts(ts.clone(), i, location, Orientation::Forward);
            }
            _ => {}
        }
    }
    panic!("no root Solid/Shell in BRep");
}

#[test]
fn builder_stage_api_contract() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("box a");
    let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box b");
    let mut filler = PaveFiller::new();
    filler.set_arguments(vec![root_shape(&a, 0), root_shape(&b, 1)]);
    filler.set_fuzzy_value(0.0);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);

    let mut builder = Builder::new(filler.ds(), BooleanOpType::Union, 0.0);
    builder.set_arguments(filler.ds().arguments.clone());
    let (_brep, snaps) = builder
        .build_with_history_stage_by_stage()
        .expect("stage api should not Err on early stop");
    assert!(!snaps.is_empty(), "expected at least one snapshot");
    for (i, s) in snaps.iter().enumerate() {
        assert_eq!(s.stage_name, STAGE_ORDER[i], "stage {} name out of order", i);
    }
    // The result BRep accumulates across the FillImages*/BuildResult stages
    // (s01-s08). The last two snapshots bracket BuildShape (s09), which fuses
    // the split solids into the boolean result and legitimately reduces the
    // entity counts, so the monotonic check covers only the accumulation prefix.
    if snaps.len() == STAGE_ORDER.len() {
        for w in snaps[..snaps.len() - 2].windows(2) {
            assert!(w[1].n_brep_faces >= w[0].n_brep_faces, "result faces non-decreasing");
            assert!(w[1].n_brep_vertices >= w[0].n_brep_vertices, "result vertices non-decreasing");
        }
    }
}
