//! Cross-configuration and structural checks for boolean pipelines.
//!
//! - BVH vs non-BVH pave paths (`boolean_op` vs `boolean_op_with_options` / `use_bvh: false`)
//! - Serial vs parallel history builds (`boolean_op_with_history` vs `boolean_op_par`)
//! - `A ∩ B` vs `B ∩ A` volume agreement
//! - [`validate_solid_closure`] on successful outputs (edge–manifold closure)
//!
//! `use_bvh: false` uses the same pave/build stages as the BVH path; post-processing now matches
//! [`boolean_op_pave_fill_build`] so volumes stay aligned (see `boolean_postprocess_pave_result` in `lib.rs`).

use rcad_algorithms::tolerance::*;
use std::collections::HashSet;

use glam::DVec3;
use rcad_algorithms::{
    boolean_op, boolean_op_par, boolean_op_with_history, boolean_op_with_options, bvh::Bvh,
    total_volume, validate_solid_closure, BooleanOpType, BooleanOptions,
};
use rcad_kernel::BRep;
use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep};

fn face_count(brep: &BRep) -> usize {
    brep.solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

/// Options that only toggle the pave implementation; no healing, simplify, or make-connected.
fn bare_options(use_bvh: bool) -> BooleanOptions {
    BooleanOptions {
        use_bvh,
        run_healing: false,
        run_simplify: false,
        run_make_connected: false,
        include_history: false,
        fuzzy_tol: 0.0,
        use_glue: false,
        ..Default::default()
    }
}

fn assert_volume_close(a: f64, b: f64, label: &str) {
    let scale = a.abs().max(b.abs()).max(1.0);
    let tol = (TOLERANCE_ABS * scale).max(TOLERANCE_COORD_SUB);
    assert!(
        (a - b).abs() <= tol,
        "{label}: volume mismatch {a} vs {b} (tol {tol})"
    );
}

fn assert_bvh_non_bvh_agree(op: BooleanOpType, a: &BRep, b: &BRep, label: &str) {
    let with_bvh = boolean_op_with_options(op, a, b, bare_options(true))
        .unwrap_or_else(|e| panic!("{label} with BVH: {e:?}"));
    let no_bvh = boolean_op_with_options(op, a, b, bare_options(false))
        .unwrap_or_else(|e| panic!("{label} without BVH: {e:?}"));

    let v0 = total_volume(&with_bvh.0);
    let v1 = total_volume(&no_bvh.0);
    assert_volume_close(v0, v1, label);

    let f0 = face_count(&with_bvh.0);
    let f1 = face_count(&no_bvh.0);
    assert_eq!(
        f0, f1,
        "{label}: face count with_bvh={f0} no_bvh={f1}"
    );
}

#[test]
fn intersection_boxes_bvh_matches_non_bvh() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
    assert_bvh_non_bvh_agree(BooleanOpType::Intersection, &a, &b, "intersect boxes");
}

#[test]
fn difference_boxes_bvh_matches_non_bvh() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b");
    assert_bvh_non_bvh_agree(BooleanOpType::Difference, &a, &b, "cut boxes");
}

#[test]
fn intersection_sphere_box_bvh_matches_non_bvh() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    let b = make_box_brep(DVec3::new(-1.0, -1.0, -1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");
    assert_bvh_non_bvh_agree(BooleanOpType::Intersection, &s, &b, "sphere ∩ box");
}

/// Every face pair whose per-face AABBs intersect must appear in `Bvh::candidate_pairs`
/// (otherwise the pave face–face pass can drop real intersections).
#[test]
fn bvh_candidate_pairs_sound_sphere_cylinder_modeling() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    let c = make_cylinder_brep(DVec3::new(1.0, 0.0, -3.0), DVec3::Z, DVec3::X, 0.8, 6.0).expect("cyl");
    let bvh_a = Bvh::build(&s);
    let bvh_b = Bvh::build(&c);
    let cand: HashSet<(usize, usize)> = Bvh::candidate_pairs(&bvh_a, &bvh_b).into_iter().collect();
    let mut missing = Vec::new();
    for fa in 0..bvh_a.face_count() {
        for fb in 0..bvh_b.face_count() {
            let Some(aa) = bvh_a.face_aabb(fa) else { continue };
            let Some(ab) = bvh_b.face_aabb(fb) else { continue };
            if aa.intersects(ab) && !cand.contains(&(fa, fb)) {
                missing.push((fa, fb));
            }
        }
    }
    assert!(
        missing.is_empty(),
        "candidate_pairs missing {} aabb-intersecting face pair(s): {:?}",
        missing.len(),
        &missing[..missing.len().min(20)]
    );
}

#[test]
fn intersection_sphere_cylinder_bvh_matches_non_bvh() {
    let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
    let c = make_cylinder_brep(DVec3::new(1.0, 0.0, -3.0), DVec3::Z, DVec3::X, 0.8, 6.0).expect("cyl");
    assert_bvh_non_bvh_agree(BooleanOpType::Intersection, &s, &c, "sphere ∩ cylinder");
}

#[test]
fn intersection_commutes_volume_on_boxes() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
    let ab = boolean_op(BooleanOpType::Intersection, &a, &b).expect("a∩b");
    let ba = boolean_op(BooleanOpType::Intersection, &b, &a).expect("b∩a");
    assert_volume_close(
        total_volume(&ab),
        total_volume(&ba),
        "intersection symmetry volume",
    );
    assert_eq!(
        face_count(&ab),
        face_count(&ba),
        "intersection symmetry faces"
    );
}

#[test]
fn serial_and_parallel_nonunion_match_volume_and_faces() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");

    let (s_int, _) = boolean_op_with_history(BooleanOpType::Intersection, &a, &b).expect("serial ∩");
    let (p_int, _) = boolean_op_par(BooleanOpType::Intersection, &a, &b).expect("par ∩");
    assert_volume_close(
        total_volume(&s_int),
        total_volume(&p_int),
        "serial vs par intersection volume",
    );
    assert_eq!(face_count(&s_int), face_count(&p_int), "serial vs par ∩ faces");

    let (s_cut, _) = boolean_op_with_history(BooleanOpType::Difference, &a, &b).expect("serial cut");
    let (p_cut, _) = boolean_op_par(BooleanOpType::Difference, &a, &b).expect("par cut");
    assert_volume_close(
        total_volume(&s_cut),
        total_volume(&p_cut),
        "serial vs par difference volume",
    );
    assert_eq!(face_count(&s_cut), face_count(&p_cut), "serial vs par cut faces");

    let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).expect("torus");
    let cyl = make_cylinder_brep(DVec3::new(0.0, -3.0, 0.0), DVec3::Y, DVec3::X, 0.3, 6.0).expect("cyl");
    let (s_d, _) =
        boolean_op_with_history(BooleanOpType::Difference, &torus, &cyl).expect("serial torus-cyl");
    let (p_d, _) = boolean_op_par(BooleanOpType::Difference, &torus, &cyl).expect("par torus-cyl");
    assert_volume_close(
        total_volume(&s_d),
        total_volume(&p_d),
        "serial vs par torus−cylinder volume",
    );
    assert_eq!(
        face_count(&s_d),
        face_count(&p_d),
        "serial vs par torus−cylinder faces"
    );
}

#[test]
fn serial_and_parallel_union_match_volume_and_faces() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(1.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
    let (s_u, _) = boolean_op_with_history(BooleanOpType::Union, &a, &b).expect("serial ∪");
    let (p_u, _) = boolean_op_par(BooleanOpType::Union, &a, &b).expect("par ∪");
    assert_volume_close(
        total_volume(&s_u),
        total_volume(&p_u),
        "serial vs par union volume",
    );
    assert_eq!(face_count(&s_u), face_count(&p_u), "serial vs par ∪ faces");
}

/// Fuse path (`use_bvh: true`, default) vs full pave path (`use_bvh: false`) for union — both must agree.
#[test]
fn union_default_matches_non_bvh_pave_volume() {
    let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
    let b = make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
    let fused = boolean_op(BooleanOpType::Union, &a, &b).expect("fuse union");
    let paved = boolean_op_with_options(BooleanOpType::Union, &a, &b, bare_options(false))
        .expect("pave union")
        .0;
    assert_volume_close(
        total_volume(&fused),
        total_volume(&paved),
        "fuse vs pave union volume",
    );
    assert_eq!(
        face_count(&fused),
        face_count(&paved),
        "fuse vs pave union faces"
    );
}

/// Closed-shell edge-manifold check (same family as post-boolean validation in `boolean_op_with_options`).
/// Full `brep_check_analyze` wire-closure checks currently flag some valid boolean outputs; narrow scope.
#[test]
fn solid_closure_clean_after_representative_booleans() {
    let cases: &[(&str, fn() -> BRep)] = &[
        ("box ∩ box", || {
            let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
            let b = make_box_brep(DVec3::new(0.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
            boolean_op(BooleanOpType::Intersection, &a, &b).expect("∩")
        }),
        ("box − box", || {
            let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
            let b = make_box_brep(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).expect("b");
            boolean_op(BooleanOpType::Difference, &a, &b).expect("−")
        }),
        ("box ∪ box", || {
            let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("a");
            let b = make_box_brep(DVec3::new(1.5, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("b");
            boolean_op(BooleanOpType::Union, &a, &b).expect("∪")
        }),
        ("sphere ∩ box", || {
            let s = make_sphere_brep(DVec3::ZERO, 2.0).expect("sphere");
            let b = make_box_brep(DVec3::new(-1.0, -1.0, -1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).expect("box");
            boolean_op(BooleanOpType::Intersection, &s, &b).expect("sphere∩box")
        }),
    ];

    for (name, build) in cases.iter() {
        let brep = build();
        let r = validate_solid_closure(&brep);
        assert!(r.is_clean(), "{name}: validate_solid_closure failed: {:?}", r.issues);
    }
}
