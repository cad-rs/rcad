//! Generate PaveFiller pipeline dumps for bcommon_simple A1.
//! Run with: RCAD_DUMP_PIPELINE=1 RCAD_DUMP_DIR=./target/pipeline_dumps cargo test --test pipeline_dump_gen -- --nocapture

use glam::DVec3;
use rcad_algorithms::pave_filler::PaveFiller;
use rcad_algorithms::bopds::ds::DS;
use rcad_algorithms::bop_occt_ops::pave_fill;
use rcad_algorithms::BooleanOpType;
use rcad_kernel::topods;
use rcad_modeling::{make_box_brep, make_sphere_brep};

#[test]
fn generate_bcommon_simple_a1_dump() {
    // bcommon_simple A1: psphere s 1; box b 1 1 1; bcommon result s b
    let sphere = make_sphere_brep(DVec3::ZERO, 1.0).expect("psphere s 1");
    let box_b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("box b 1 1 1");

    // Run the PaveFiller pipeline with dump enabled
    let mut brep = topods::BRep::new();
    let (mut ds, face_refs, ic_edge_map) = rcad_algorithms::bop_occt_ops::pave_fill(
        &sphere, &box_b, true, &mut brep, rcad_algorithms::tolerance::TOLERANCE_ABS);

    println!("DS: V={} E={} F={} IC={} PB={} CB={}",
        ds.vertices.len(), ds.edges.len(), ds.faces.len(),
        ds.intersection_curves.len(),
        ds.pave_blocks.len(), ds.common_blocks.len());

    // Print key counts at each stage (stages are dumped by PaveFiller::perform internally)
    let n_vv = ds.interf_vv.len();
    let n_ve = ds.interf_ve.len();
    let n_ee = ds.interf_ee.len();
    let n_vf = ds.interf_vf.len();
    let n_ef = ds.interf_ef.len();
    let n_ff = ds.interf_ff.len();
    println!("Interfs: VV={} VE={} EE={} VF={} EF={} FF={}",
        n_vv, n_ve, n_ee, n_vf, n_ef, n_ff);

    // Count edges by type
    let n_new_edges = ds.edges.iter().filter(|e| e.origin == rcad_algorithms::bopds::ds::ShapeOrigin::None).count();
    println!("Edges: total={} new={}",
        ds.edges.len(), n_new_edges);
}
