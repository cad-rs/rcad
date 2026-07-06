//! Minimal repro for **box–box union** when operands use different extents along local Y and Z
//! (not only along the overlap direction X). This is the case where Python `rcad-py-ocp-tests`
//! sees **more faces from RCAD than from OCCT `BRepAlgoAPI_Fuse`** (~14 on OCCT for the same metric).
//!
//! ## How to use during kernel work
//! - **Regression snapshot**: [`ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_RAW`] pins the current raw
//!   `boolean_op` face tally; lower it when merge / same-domain passes improve.
//! - **Target parity** (ignored until fixed): [`asymmetric_two_box_union_matches_occt_fourteen_faces`].
//! - **Pipeline diff**: run
//!   `cargo test -p rcad-algorithms --test asymmetric_box_union_occt_gap asymmetric_two_box_union_raw_vs_simplified -- --nocapture`
//!   to print face counts before vs after `boolean_op_simplified` (default simplify options).

use glam::DVec3;
use rcad_algorithms::geom_populate::populate_box_geom;
use rcad_algorithms::{boolean_op, boolean_op_simplified, BooleanOpType, SimplifyOptions};
use rcad_kernel::properties::volume;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_modeling::make_box_brep;

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

/// Same layout as the probe in `rcad-py-ocp-tests` (overlapping boxes, different Y/Z sizes).
fn two_asymmetric_boxes() -> (BRep, BRep) {
    let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 3.5, 2.5, 2.2).expect("box a");
    let mut b =
        make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.8, 2.0, 3.0).expect("box b");
    populate_box_geom(&mut a);
    populate_box_geom(&mut b);
    (a, b)
}

#[test]
fn asymmetric_two_box_union_succeeds_nonempty_volume() {
    let (a, b) = two_asymmetric_boxes();
    let out = boolean_op(BooleanOpType::Union, &a, &b).expect("union");
    assert!(!out.solids.is_empty());
    assert!(face_count(&out) >= 6);
    let v = volume(&out);
    assert!(v.is_finite() && v > 0.0, "unexpected volume {v}");
}

/// OCCT `fuse` face count for this configuration (see `occt-test-gen` / `rcad-py-ocp-tests`).
const OCCT_FACE_COUNT_FUSE_REFERENCE: usize = 14;

/// Update when orthogonal / same-domain merging converges toward OCCT.
const ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_RAW: usize = 14;

/// Snapshot after [`boolean_op_simplified`] with [`SimplifyOptions::default`] (merge activity).
const ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_AFTER_DEFAULT_SIMPLIFY: usize = 14;

#[test]
fn asymmetric_two_box_union_face_count_regression() {
    let (a, b) = two_asymmetric_boxes();
    let raw = boolean_op(BooleanOpType::Union, &a, &b).expect("union");
    assert_eq!(
        face_count(&raw),
        ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_RAW,
        "raw boolean_op face count drifted; if this decreased toward {OCCT_FACE_COUNT_FUSE_REFERENCE}, \
         adjust the constant and enable/trim the ignored OCCT parity test"
    );
}

#[test]
fn asymmetric_two_box_union_raw_vs_simplified() {
    let (a, b) = two_asymmetric_boxes();
    let raw = boolean_op(BooleanOpType::Union, &a, &b).expect("union");
    let n_raw = face_count(&raw);
    let (simp, rep) = boolean_op_simplified(
        BooleanOpType::Union,
        &a,
        &b,
        SimplifyOptions::default(),
    )
    .expect("union simplified");
    let n_s = face_count(&simp);

    eprintln!(
        "asymmetric box-box union: faces raw={n_raw} simplified={n_s} | \
         simplify: same_domain_face_merges={} orthogonal_coplanar_fusions={} internal_faces_removed={} \
         wires_fixed={} vertices_merged={}",
        rep.same_domain_face_merges,
        rep.orthogonal_coplanar_fusions,
        rep.internal_faces_removed,
        rep.wires_fixed,
        rep.vertices_merged
    );

    assert_eq!(n_raw, ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_RAW);
    assert_eq!(
        n_s,
        ASYMMETRIC_BOX_UNION_RCAD_FACE_COUNT_AFTER_DEFAULT_SIMPLIFY,
        "simplified face tally drifted; OCCT reference remains {OCCT_FACE_COUNT_FUSE_REFERENCE}"
    );
    assert!(
        n_s <= n_raw,
        "simplified face count {n_s} should not exceed raw {n_raw}"
    );
}

#[test]
fn asymmetric_two_box_union_matches_occt_fourteen_faces() {
    let (a, b) = two_asymmetric_boxes();
    let raw = boolean_op(BooleanOpType::Union, &a, &b).expect("union");
    assert_eq!(
        face_count(&raw),
        OCCT_FACE_COUNT_FUSE_REFERENCE,
        "align RCAD union merges with OCCT BRepAlgoAPI_Fuse for asymmetric box operands"
    );
}
