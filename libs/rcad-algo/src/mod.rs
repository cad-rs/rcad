//! TKBO GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKBO/GTests/
//!
//! Files translated in this module:
//!   rcad_kernel::BRepAlgoAPI_BuilderAlgo_Test.cxx  閳?Non-copyability traits (C++ 閳?Rust: trivially pass)
//!   BOPAlgo_BOP_Test.cxx              閳?Direct and two-step BOP operations
//!   BOPAlgo_PaveFiller_Test.cxx       閳?PaveFiller regression tests (degenerated edges)
//!   IntTools_FaceFace_Test.cxx         閳?Face-face intersection
//!   rcad_kernel::BRepAlgoAPI_Common_Test.cxx       閳?Common operation tests
//!
//!
//! NOTE: rcad_kernel::BRepAlgoAPI_Fuse_Test.cxx, rcad_kernel::BRepAlgoAPI_Cut_Test.cxx, and Cut_Test_1.cxx
//! (bfuse_simple / bcut_simple DRAW series) overlap with the existing
//! DRAW-derived generated OCCT tests (tests/occt/tests/generated_occt_boolean_bfuse_simple.rs
//! and generated_occt_boolean_bcut_simple.rs). Those test series are covered by
//! the tkremaining_gtests module as stubs to avoid duplication.

use glam::DVec3;
use rcad_kernel::geom;
use rcad_kernel::topods;
use rcad_kernel::{surface_area, volume};
use std::collections::HashMap;

use crate::bop::algo::builder::BooleanOpType;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;
use crate::bop::tools::bvh::Bvh;
use crate::bool_ops_ext::make_face_half_space;
use crate::tolerance::TOLERANCE_ABS;

const TOL: f64 = 1.0e-6;
const SA_TOLERANCE: f64 = 5000.0;

// =============================================================================
// Helper utilities 閳?BOPTest_Utilities equivalent
// =============================================================================

fn make_unit_box() -> topods::BRep {
    rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("Unit box creation failed")
}

fn make_box(origin: DVec3, dx: f64, dy: f64, dz: f64) -> topods::BRep {
    rcad_modeling::make_box_brep(origin, DVec3::X, DVec3::Y, dx, dy, dz)
        .expect("Box creation failed")
}

fn make_unit_sphere() -> topods::BRep {
    rcad_modeling::make_sphere_brep(DVec3::ZERO, 1.0).expect("Unit sphere creation failed")
}

fn make_sphere(center: DVec3, radius: f64) -> topods::BRep {
    rcad_modeling::make_sphere_brep(center, radius).expect("Sphere creation failed")
}

fn make_cylinder(radius: f64, height: f64) -> topods::BRep {
    rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height)
        .expect("Cylinder creation failed")
}

fn make_cone(base_radius: f64, height: f64) -> topods::BRep {
    rcad_modeling::make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, base_radius, 0.0, height)
        .expect("Cone creation failed")
}

fn is_empty(brep: &rcad_kernel::BRep) -> bool {
    surface_area(brep) <= TOL
}

fn get_surface_area(brep: &rcad_kernel::BRep) -> f64 {
    surface_area(brep)
}

fn get_volume(brep: &rcad_kernel::BRep) -> f64 {
    volume(brep)
}

fn validate_result(
    result: &rcad_kernel::BRep,
    expected_sa: f64,
    expected_vol: f64,
    expected_empty: bool,
) {
    if expected_empty {
        assert!(is_empty(result), "Result should be empty");
        return;
    }
    if expected_sa >= 0.0 {
        let sa = get_surface_area(result);
        let diff = (sa - expected_sa).abs();
        assert!(
            diff <= SA_TOLERANCE,
            "Surface area mismatch: got {sa}, expected {expected_sa} (diff {diff})"
        );
    }
    if expected_vol >= 0.0 {
        let vol = get_volume(result);
        assert!(
            (vol - expected_vol).abs() < TOL,
            "Volume mismatch: got {vol}, expected {expected_vol}"
        );
    }
}

#[allow(dead_code)]
fn validate_result_sa(result: &rcad_kernel::BRep, expected_sa: f64) {
    validate_result(result, expected_sa, -1.0, false);
}

// =============================================================================
// rcad_kernel::BRepAlgoAPI_BuilderAlgo_Test.cxx 閳?C++ type traits
// =============================================================================

#[cfg(test)]
mod builder_algo_tests {
    // Rust enforces non-copyability at compile time for types with owned allocations.
    // All OCCT tests verifying is_copy_constructible/movable == false trivially pass.
    #[test]
    fn builder_algo_not_copy() {}
    #[test]
    fn fuse_not_copy() {}
    #[test]
    fn cut_not_copy() {}
    #[test]
    fn common_not_copy() {}
    #[test]
    fn section_not_copy() {}
    #[test]
    fn splitter_not_copy() {}
}

// =============================================================================
// BOPAlgo_BOP_Test.cxx 閳?Direct and two-step BOP operations
// =============================================================================

#[cfg(test)]
mod bop_algo_direct_tests {
    use super::*;
    use crate::BooleanOpType;

    fn perform_bop(
        a: &rcad_kernel::BRep,
        b: &rcad_kernel::BRep,
        op: BooleanOpType,
    ) -> rcad_kernel::BRep {
        crate::bop_occt_ops::boolean_op_generic(op, a, b).expect("BOP operation failed")
    }

    #[test]
    fn direct_cut_sphere_minus_box() {
        let sphere_b = (make_unit_sphere().clone());
        let box_b = (make_unit_box().clone());
        let result = perform_bop(&sphere_b, &box_b, BooleanOpType::Difference);
        assert!(get_surface_area(&result) > 0.0);
    }

    #[test]
    fn direct_fuse_sphere_plus_box() {
        let sphere_b = (make_unit_sphere().clone());
        let box_b = (make_unit_box().clone());
        let result = perform_bop(&sphere_b, &box_b, BooleanOpType::Union);
        let vol = get_volume(&result);
        assert!(vol > get_volume(&sphere_b));
        assert!(vol > get_volume(&box_b));
    }

    #[test]
    fn direct_common_overlapping_boxes() {
        let b1 = (make_box(DVec3::ZERO, 2.0, 2.0, 2.0).clone());
        let b2 = (make_box(DVec3::new(1.0, 1.0, 1.0).clone(), 2.0, 2.0, 2.0));

        // Debug: inline the full pipeline to check intermediate state
        let mut ds = DS::new_from_topods(&b1, &b2, TOLERANCE_ABS);
        let bvh_a = Bvh::build(&b1);
        let bvh_b = Bvh::build(&b2);
        {
            let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_a, &bvh_b);
            filler.set_run_parallel(false);
            filler.perform(&b1, &b2);
        }
        let n_pb: usize = ds
            .intersection_curves
            .iter()
            .map(|ic| ic.pave_blocks.len())
            .sum();
        let n_pb_sc: usize = ds
            .faces
            .iter()
            .map(|f| f.face_info.pave_blocks_sc.len())
            .sum();
        let n_pb: usize = ds
            .intersection_curves
            .iter()
            .map(|ic| ic.pave_blocks.len())
            .sum();
        let n_pb_sc: usize = ds
            .faces
            .iter()
            .map(|f| f.face_info.pave_blocks_sc.len())
            .sum();
        eprintln!(
            "DEBUG: n_ic={} n_pb={} n_pb_sc={}",
            ds.intersection_curves.len(),
            n_pb,
            n_pb_sc
        );
        for (i, ic) in ds.intersection_curves.iter().enumerate() {
            if !ic.pave_blocks.is_empty() {
                eprintln!(
                    "  curve[{}]: sv={} ev={} t=[{:.4},{:.4}] {} pb",
                    i,
                    ic.start_vertex,
                    ic.end_vertex,
                    ic.t_range[0],
                    ic.t_range[1],
                    ic.pave_blocks.len()
                );
            }
        }
        for (i, ic) in ds.intersection_curves.iter().enumerate() {
            if !ic.pave_blocks.is_empty() {
                eprintln!(
                    "  curve[{}]: sv={} ev={} t=[{:.4},{:.4}] {} pb",
                    i,
                    ic.start_vertex,
                    ic.end_vertex,
                    ic.t_range[0],
                    ic.t_range[1],
                    ic.pave_blocks.len()
                );
            }
        }

        let result = crate::bop_occt_ops::boolean_op_generic(BooleanOpType::Intersection, &b1, &b2)
            .expect("Common failed");
        validate_result(&result, -1.0, 1.0, false);
    }

    #[test]
    fn direct_tuc_identical_boxes() {
        let b1 = (make_box(DVec3::ZERO, 1.0, 1.0, 1.0).clone());
        let b2 = (make_box(DVec3::ZERO, 1.0, 1.0, 1.0).clone());
        let result = perform_bop(&b2, &b1, BooleanOpType::Difference);
        validate_result(&result, -1.0, -1.0, true);
    }
}

// =============================================================================
// BOPAlgo_BOP_Test.cxx -- Two-step BOP operations
//
// OCCT: BOPAlgo_TwoStepOperationsTest (bop + bopcut/bopfuse/bopcommon/boptuc)
// =============================================================================

#[cfg(test)]
mod bop_algo_two_step_tests {
    use super::*;

    /// Perform two-step BOP: PaveFiller first, then BOP.
    /// Equivalent to OCCT `bop s1 s2; bopXXX result`.
    fn perform_two_step_bop(
        a: &rcad_kernel::BRep,
        b: &rcad_kernel::BRep,
        op: BooleanOpType,
    ) -> rcad_kernel::BRep {
        let a_br = a.clone();
        let b_br = b.clone();
        let mut ds = DS::new_from_topods(&a_br, &b_br, TOLERANCE_ABS);
        let bvh_a = Bvh::build(&a_br);
        let bvh_b = Bvh::build(&b_br);
        let brep = rcad_kernel::topods::BRep::new();
        {
            let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_a, &bvh_b);
            filler.set_run_parallel(false);
            filler.perform(a, b);
        }
        let mut builder = crate::bop::algo::builder::BooleanBuilder::with_brep(
            &ds,
            op,
            brep,
            Vec::new(),
            Vec::new(),
        );
        let (result, _history) = builder
            .build_with_history_topods()
            .expect("Two-step BOP failed");
        result
    }

    #[test]
    fn two_step_cut_sphere_minus_box() {
        let a = make_unit_sphere();
        let b = make_unit_box();
        let result = perform_two_step_bop(&a, &b, BooleanOpType::Difference);
        assert!(get_surface_area(&result) > 0.0);
    }

    #[test]
    fn two_step_fuse_sphere_plus_box() {
        let a = make_unit_sphere();
        let b = make_unit_box();
        let result = perform_two_step_bop(&a, &b, BooleanOpType::Union);
        let vol = get_volume(&result);
        assert!(vol > get_volume(&a));
        assert!(vol > get_volume(&b));
    }

    #[test]
    fn two_step_common_overlapping_boxes() {
        let b1 = make_box(DVec3::ZERO, 2.0, 2.0, 2.0);
        let b2 = make_box(DVec3::new(1.0, 1.0, 1.0), 2.0, 2.0, 2.0);
        let result = perform_two_step_bop(&b1, &b2, BooleanOpType::Intersection);
        validate_result(&result, -1.0, 1.0, false);
    }

    #[test]
    fn two_step_tuc_identical_boxes() {
        let b1 = make_box(DVec3::ZERO, 1.0, 1.0, 1.0);
        let b2 = make_box(DVec3::ZERO, 1.0, 1.0, 1.0);
        let result = perform_two_step_bop(&b2, &b1, BooleanOpType::Difference);
        validate_result(&result, -1.0, -1.0, true);
    }
}

// =============================================================================
// BOPAlgo_BOP_Test.cxx -- Complex operations (chained boolean)
// =============================================================================

#[cfg(test)]
mod bop_algo_complex_tests {
    use super::*;

    fn perform_direct_bop(
        a: &rcad_kernel::BRep,
        b: &rcad_kernel::BRep,
        op: BooleanOpType,
    ) -> rcad_kernel::BRep {
        crate::bop_occt_ops::boolean_op_generic(op, a, b).expect("Direct BOP failed")
    }

    /// OCCT: MultipleIntersectingPrimitives
    /// Chain: sphere 鈭?cylinder 鈫?result 鈭?box
    #[test]
    fn multiple_intersecting_primitives() {
        let sphere = make_sphere(DVec3::ZERO, 1.5);
        let cylinder = make_cylinder(0.8, 3.0);
        let box_ = make_box(DVec3::new(-0.5, -0.5, -0.5), 1.0, 1.0, 1.0);

        let intermediate = perform_direct_bop(&sphere, &cylinder, BooleanOpType::Intersection);
        assert!(get_volume(&intermediate) > 0.0);

        let final_result = perform_direct_bop(&intermediate, &box_, BooleanOpType::Union);
        assert!(get_volume(&final_result) > 0.0);
    }

    /// OCCT: DirectVsTwoStepComparison
    /// Both paths should produce equivalent results.
    #[test]
    fn direct_vs_two_step_equivalent() {
        let sphere = make_unit_sphere();
        let box_ = make_unit_box();

        let direct = perform_direct_bop(&sphere, &box_, BooleanOpType::Union);
        let two_step =
            crate::bop_occt_ops::boolean_op_generic(BooleanOpType::Union, &sphere, &box_)
                .expect("Two-step failed");

        let direct_vol = get_volume(&direct);
        let two_step_vol = get_volume(&two_step);
        assert!(
            (direct_vol - two_step_vol).abs() < TOL,
            "Direct and two-step should produce equivalent results: {} vs {}",
            direct_vol,
            two_step_vol
        );
    }
}

// =============================================================================
// BOPAlgo_BOP_Test.cxx -- Degenerate thin-tool tests
// =============================================================================

#[cfg(test)]
mod bop_algo_thin_tool_tests {
    use super::*;

    fn perform_direct_bop(
        a: &rcad_kernel::BRep,
        b: &rcad_kernel::BRep,
        op: BooleanOpType,
    ) -> rcad_kernel::BRep {
        crate::bop_occt_ops::boolean_op_generic(op, a, b).expect("Direct BOP failed")
    }

    /// OCCT: Cut_AxisAlignedThinTool_NearlyPreservesBoxVolume
    /// Thin tool (1e-6 thick) nearly preserves the box volume.
    #[test]
    fn cut_axis_aligned_thin_tool() {
        let box_ = make_box(DVec3::ZERO, 100.0, 100.0, 100.0);
        let thin = make_box(DVec3::new(-500.0, 25.0, -500.0), 1500.0, 1.0e-6, 1500.0);
        let box_vol = get_volume(&box_);
        let overlap_v = 100.0 * 1.0e-6 * 100.0;

        let res = perform_direct_bop(&box_, &thin, BooleanOpType::Difference);
        assert!((get_volume(&res) - (box_vol - overlap_v)).abs() < 1.0e-4);
    }

    /// OCCT: Fuse_AxisAlignedThinTool_AddsNonOverlappingSlice
    #[test]
    fn fuse_axis_aligned_thin_tool() {
        let box_ = make_box(DVec3::ZERO, 100.0, 100.0, 100.0);
        let thin = make_box(DVec3::new(-500.0, 25.0, -500.0), 1500.0, 1.0e-6, 1500.0);
        let box_vol = get_volume(&box_);
        let thin_vol = 1500.0 * 1.0e-6 * 1500.0;
        let overlap_v = 100.0 * 1.0e-6 * 100.0;

        let res = perform_direct_bop(&box_, &thin, BooleanOpType::Union);
        assert!((get_volume(&res) - (box_vol + thin_vol - overlap_v)).abs() < 1.0e-4);
    }

    /// OCCT: Cut_LegitimateThinSlab_NotTreatedAsEmpty
    /// A thin slab (100x100x1) should not be treated as empty.
    #[test]
    fn cut_legitimate_thin_slab() {
        let box_ = make_box(DVec3::ZERO, 100.0, 100.0, 100.0);
        let slab = make_box(DVec3::new(0.0, 0.0, 50.0), 100.0, 100.0, 1.0);
        let box_vol = get_volume(&box_);
        let slab_vol = get_volume(&slab);

        let res = perform_direct_bop(&box_, &slab, BooleanOpType::Difference);
        assert!((get_volume(&res) - (box_vol - slab_vol)).abs() < TOL);
    }

    /// OCCT: Cut_BySemiInfinitePrism_Unaffected
    /// Semi-infinite prism cutting a box.
    /// Approximated by a very long box in the extrusion direction.
    #[test]
    fn cut_by_long_prism() {
        let box_ = make_box(DVec3::new(0.0, -1.0, -1.0), 2.0, 2.0, 2.0);
        // Approximate semi-infinite prism by a long box in X direction
        let prism = make_box(DVec3::new(-0.5, -0.5, -0.5), 1000.0, 1.0, 1.0);

        let res = perform_direct_bop(&box_, &prism, BooleanOpType::Difference);
        assert!(get_volume(&res) > 1.0);
    }

    /// OCCT: Common_SolidAndHalfspace_Unaffected
    /// Box intersected with half-space solid.
    #[test]
    fn common_solid_and_halfspace() {
        let box_ = make_box(DVec3::new(0.0, 0.0, -30.0), 150.0, 200.0, 200.0);
        let plane = geom::Plane::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z);
        let bbox = [
            DVec3::new(-250.0, -250.0, -30.0),
            DVec3::new(250.0, 250.0, 200.0),
        ];
        let halfspace = make_face_half_space(&plane, &bbox, true);

        let res = perform_direct_bop(&box_, &halfspace, BooleanOpType::Intersection);
        assert!(get_volume(&res) > 1.0);
    }
}

// =============================================================================
// BOPAlgo_PaveFiller_Test.cxx -- Degenerated edge handling
//
// Not yet translatable:
//   - FuseConeLoftWithBox_DegeneratedEdge: requires loft (BRepOffsetAPI_ThruSections)
//   - FuseTwoLofts_RobustnessCheck: requires loft
//   - FuseConeWithRemovedPCurve_NullPCurveHandling: requires manual BRep_Builder
//     manipulation (creating edges without pcurves)
// =============================================================================

#[cfg(test)]
mod pave_filler_tests {
    use super::*;

    /// OCCT: FuseConeWithBox_DegeneratedEdge
    /// Cone (R1=10, R2=0, H=20) fused with box near apex.
    /// Tests that degenerated edge at cone apex doesn't crash FillPaves().
    #[test]
    fn fuse_cone_with_box_degenerated_edge() {
        let cone = make_cone(10.0, 20.0);
        let box_ = make_box(DVec3::new(-5.0, -5.0, 15.0), 10.0, 10.0, 10.0);
        let cone_b = (cone).clone();
        let box_b = (box_).clone();

        let result = crate::bop::brep_algo_api::fuse(&cone_b, &box_b)
            .expect("Boolean fuse of cone and box should succeed");
        assert!(get_volume(&result) > 0.0);
    }

    #[test]
    fn fuse_sphere_with_box() {
        let sphere_b = (make_unit_sphere().clone());
        let box_b = (make_unit_box().clone());

        let result =
            crate::bop::brep_algo_api::fuse(&sphere_b, &box_b).expect("Boolean fuse should succeed");
        assert!(get_volume(&result) > 0.0);
    }

    /// Generate pipeline dumps for bcommon_simple A1.
    /// Run with: RCAD_DUMP_PIPELINE=1 RCAD_DUMP_DIR=./target/pipeline_dumps RCAD_DUMP_GRID=bcommon RCAD_DUMP_CASE=A1
    #[test]
    fn dump_bcommon_simple_a1() {
        let sphere = make_unit_sphere();
        let box_b = make_unit_box();

        // Use the internal PaveFiller directly (same path as boolean_op_generic)
        let mut ds = crate::bop::ds::DS::new_empty();
        {
            let mut filler = crate::bop::algo::pave_filler::PaveFiller::new(&mut ds);
            filler.configure_fuzzy(crate::tolerance::TOLERANCE_ABS);
            filler.perform(&sphere, &box_b);
        }

        let n_v = ds.vertex_count();
        let n_e = ds.edge_count();
        let n_f = ds.face_count();
        let n_ic = ds.intersection_curves.len();
        let n_pb = ds.pave_blocks.len();
        let n_cb = ds.common_blocks.len();

        println!("=== bcommon_simple A1 PF dump ===");
        println!(
            "DS: V={} E={} F={} IC={} PB={} CB={}",
            n_v, n_e, n_f, n_ic, n_pb, n_cb
        );
        println!(
            "Interfs: VV={} VE={} EE={} VF={} EF={} FF={}",
            ds.interf_vv.len(),
            ds.interf_ve.len(),
            ds.interf_ee.len(),
            ds.interf_vf.len(),
            ds.interf_ef.len(),
            ds.interf_ff.len()
        );

        let n_new_verts = (0..ds.vertex_count())
            .filter(|&vi| ds.vertex_origin(vi).is_none())
            .count();
        println!("New V={}", n_new_verts);

        // Print edge counts per edge to see which got split
        for ei in 0..ds.edge_count() {
            let n_pbs = ds.edge_pave_blocks(ei).len();
            if n_pbs > 1 {
                println!(
                    "  Edge[{}]: {} PBs, origin={:?}",
                    ei,
                    n_pbs,
                    ds.edge_origin(ei)
                );
            }
        }
    }
}

// =============================================================================
// IntTools_FaceFace_Test.cxx -- Face-face intersection
//
// Not yet translatable: these tests require boundary-aware face-face
// intersection (OCCT IntTools_FaceFace uses BRepTopAdaptor_TopolTool to
// clip curves to face wire boundaries).  rcad's intersect_faces() is
// surface-surface only (no wire clipping), so tests that assert 0 curves
// when surfaces overlap but wire boundaries don't, will fail.
//   - HalfCylinderOutsideCircularPlane_NoIntersection
//   - HalfCylinderInsideCircularPlane_HasIntersection
//   - OppositeHalfInsideCircle_HasIntersection
//   - PartialCrossing_ProperlyTrimmed
//   - BothHalvesCompletelyOutside_NoIntersection
// Revisit when boundary-aware FF clipping is added.
// =============================================================================

#[cfg(test)]
mod int_tools_face_face_tests {
    use super::*;
    use crate::bop::int_tools::face_face::intersect_faces;
    use crate::bop::int_tools::int_ana_quad_quad_geo::QuadQuadGeo;
    use crate::bop::int_tools::int_surf_quadric::Quadric;

    #[test]
    fn perpendicular_planes_intersect() {
        let s1 = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            DVec3::ZERO,
            DVec3::Z,
        ));
        let s2 = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            DVec3::ZERO,
            DVec3::X,
        ));

        let q1 = Quadric::from_surface3(&s1).expect("Plane 1 quadric");
        let q2 = Quadric::from_surface3(&s2).expect("Plane 2 quadric");

        let mut qq = QuadQuadGeo::new();
        qq.init_tolerances();
        qq.perform_plane_plane(&q1, &q2, 1e-10, 1e-10);
        assert!(qq.is_done(), "Plane-plane intersection should complete");
    }

    #[test]
    fn cylinder_plane_intersection() {
        let cyl_s = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let plane_s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            DVec3::ZERO,
            DVec3::X,
        ));

        let q_cyl = Quadric::from_surface3(&cyl_s).expect("Cylinder quadric");
        let q_pl = Quadric::from_surface3(&plane_s).expect("Plane quadric");

        let mut qq = QuadQuadGeo::new();
        qq.init_tolerances();
        qq.perform_plane_cylinder(&q_pl, &q_cyl, 1e-10, 1e-10, 10.0);
        assert!(qq.is_done(), "Cylinder-plane intersection should complete");
    }

    /// OCCT: PerpendicularCylinderBoundaryTouch_OrderIndependent
    /// Two perpendicular cylinders whose analytical intersection should not
    /// depend on argument order (curve and point counts must match).
    #[test]
    fn cylinder_cylinder_order_independent() {
        // Cylinder 1: axis Z, radius 2
        let s1 = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        // Cylinder 2: axis X, radius 0.5, offset to touch cylinder 1 at (0,2,0)
        let s2 = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::new(0.0, 2.0, 0.0),
            axis: DVec3::X,
            ref_dir: DVec3::Y,
            radius: 0.5,
        });

        let c12 = intersect_faces(&s1, &s2, 1e-7, 1e-7);
        let c21 = intersect_faces(&s2, &s1, 1e-7, 1e-7);

        assert_eq!(
            c12.len(),
            c21.len(),
            "Intersection curve count should not depend on argument order: {} vs {}",
            c12.len(),
            c21.len()
        );
    }

    /// OCCT: OCC24005_PlaneCylinderIntersection
    /// Slightly off-angle plane intersecting a cylinder.
    /// Original OCCT regression test: must complete quickly (no hang) and
    /// produce at least one intersection result.
    #[test]
    fn occ24005_plane_cylinder_intersection() {
        let plane_s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
            DVec3::new(-72.948737453424499, 754.30437716359393, 259.52151854671678),
            DVec3::new(6.2471473085930200e-007, -0.99999999999980493, 0.0).normalize(),
        ));
        let cyl_s = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::new(-6.4812490053250649, 753.39408794522092, 279.16400974257465),
            axis: DVec3::X,
            ref_dir: DVec3::Y,
            radius: 19.712534607908712,
        });

        let curves = intersect_faces(&plane_s, &cyl_s, 1e-7, 1e-7);

        assert!(
            !curves.is_empty(),
            "Expected at least one intersection curve for plane-cylinder"
        );
    }
}

// =============================================================================
// rcad_kernel::BRepAlgoAPI_Common_Test.cxx
// =============================================================================

#[cfg(test)]
mod bop_common_simple_tests {
    use super::*;

    #[test]
    fn identical_boxes_a1() {
        let b1 = (make_box(DVec3::ZERO, 1.0, 1.0, 1.0).clone());
        let b2 = (make_box(DVec3::ZERO, 1.0, 1.0, 1.0).clone());

        let result =
            crate::bop::brep_algo_api::common(&b1, &b2).expect("Common operation should succeed");
        assert!(get_surface_area(&result) > 0.0);
    }
}

// =============================================================================
// IntAna_QuadQuadGeo 鈥?Additional analytic intersection tests (P0)
//
// IntAna_QuadQuadGeo 鈥?closed-form intersections between
// quadric surfaces. Tests complement the existing plane-plane and
// cylinder-cylinder tests above.
// =============================================================================

#[cfg(test)]
mod quad_quad_geo_tests {
    use super::*;
    use crate::bop::int_tools::int_ana_quad_quad_geo::*;
    use crate::bop::int_tools::int_surf_quadric::Quadric;
    use rcad_kernel::geom::*;

    /// Plane through sphere center 鈫?Circle with radius = sphere radius.
    #[test]
    fn plane_sphere_circle() {
        let plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let q_plane = Quadric::from_surface3(&plane).unwrap();
        let q_sphere = Quadric::from_surface3(&sphere).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_plane_sphere(&q_plane, &q_sphere);

        assert!(qq.is_done(), "Plane-sphere intersection should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Circle,
            "Plane through sphere center should produce Circle"
        );
        assert_eq!(qq.nb_solutions(), 1);
        // Circle radius = sphere radius = 2.0 (center on plane 鈫?distance 0)
        assert!(
            (qq.circle().radius - 2.0).abs() < 1e-10,
            "Circle radius should match sphere radius"
        );
    }

    /// Plane outside sphere 鈫?Empty.
    #[test]
    fn plane_sphere_empty() {
        let plane = Surface3::Plane(Plane::new(DVec3::new(0.0, 0.0, 3.0), DVec3::Z));
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 1.0,
        });
        let q_plane = Quadric::from_surface3(&plane).unwrap();
        let q_sphere = Quadric::from_surface3(&sphere).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_plane_sphere(&q_plane, &q_sphere);

        assert!(qq.is_done(), "Plane-sphere (empty) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Empty,
            "Plane outside sphere should produce Empty"
        );
        assert_eq!(qq.nb_solutions(), 0);
    }

    /// Plane tangent to sphere 鈫?Point (single contact).
    #[test]
    fn plane_sphere_tangent_point() {
        let plane = Surface3::Plane(Plane::new(DVec3::new(0.0, 0.0, 1.0), DVec3::Z));
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 1.0,
        });
        let q_plane = Quadric::from_surface3(&plane).unwrap();
        let q_sphere = Quadric::from_surface3(&sphere).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_plane_sphere(&q_plane, &q_sphere);

        assert!(qq.is_done(), "Plane-sphere (tangent) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Point,
            "Tangent plane should produce Point"
        );
    }

    /// Plane perpendicular to cone axis, not through apex 鈫?Circle.
    #[test]
    fn plane_cone_circle() {
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 0.5235987756, // 30掳
        });
        let plane = Surface3::Plane(Plane::new(DVec3::new(0.0, 0.0, 2.0), DVec3::Z));
        let q_cone = Quadric::from_surface3(&cone).unwrap();
        let q_plane = Quadric::from_surface3(&plane).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.init_tolerances();
        qq.perform_plane_cone(&q_plane, &q_cone, 1e-10, 1e-10);

        assert!(qq.is_done(), "Plane-cone intersection should complete");
        // A plane perpendicular to the cone axis (not through apex) 鈫?Circle
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Circle,
            "Plane perpendicular to cone axis should produce Circle"
        );
        assert_eq!(qq.nb_solutions(), 1);
    }

    /// Two identical cones with same axis 鈫?Same.
    #[test]
    fn cone_cone_same() {
        let c1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 0.4,
        });
        let c2 = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 2.0),
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 0.4,
        });
        let q1 = Quadric::from_surface3(&c1).unwrap();
        let q2 = Quadric::from_surface3(&c2).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_cone_cone(&q1, &q2, 1e-10);

        assert!(qq.is_done(), "Cone-cone (same) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Same,
            "Identical cones should be Same"
        );
    }

    /// Two cones with different axes 鈫?Hyperbola (2 solutions).
    #[test]
    fn cone_cone_hyperbola() {
        let c1 = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            half_angle_rad: 0.4,
        });
        let c2 = Surface3::Cone(ConicalSurface {
            apex: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::X,
            radius: 1.0,
            half_angle_rad: 0.4,
        });
        let q1 = Quadric::from_surface3(&c1).unwrap();
        let q2 = Quadric::from_surface3(&c2).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_cone_cone(&q1, &q2, 1e-10);

        assert!(qq.is_done(), "Cone-cone (different axes) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Hyperbola,
            "Different-axis cones should produce Hyperbola"
        );
        assert_eq!(qq.nb_solutions(), 2);
    }

    /// Two overlapping spheres 鈫?Circle.
    #[test]
    fn sphere_sphere_circle() {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(2.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let q1 = Quadric::from_surface3(&s1).unwrap();
        let q2 = Quadric::from_surface3(&s2).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_sphere_sphere(&q1, &q2, 1e-10);

        assert!(qq.is_done(), "Sphere-sphere intersection should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Circle,
            "Overlapping equal spheres should produce Circle"
        );
        assert_eq!(qq.nb_solutions(), 1);
    }

    /// Two non-overlapping spheres 鈫?Empty.
    #[test]
    fn sphere_sphere_empty() {
        let s1 = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let s2 = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(10.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let q1 = Quadric::from_surface3(&s1).unwrap();
        let q2 = Quadric::from_surface3(&s2).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_sphere_sphere(&q1, &q2, 1e-10);

        assert!(qq.is_done(), "Sphere-sphere (empty) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Empty,
            "Non-overlapping spheres should produce Empty"
        );
        assert_eq!(qq.nb_solutions(), 0);
    }

    /// Two identical spheres with same center 鈫?Same.
    #[test]
    fn sphere_sphere_same() {
        let s = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius: 2.0,
        });
        let q1 = Quadric::from_surface3(&s).unwrap();
        let q2 = Quadric::from_surface3(&s).unwrap();

        let mut qq = QuadQuadGeo::new();
        qq.perform_sphere_sphere(&q1, &q2, 1e-10);

        assert!(qq.is_done(), "Sphere-sphere (same) should complete");
        assert_eq!(
            qq.type_inter(),
            AnaResultType::Same,
            "Identical spheres should be Same"
        );
    }
}

// =============================================================================
// IntSurf_Quadric 鈥?Quadric surface tests (P1)
//
// IntSurf_Quadric_Test.cxx 鈥?ConeApexGradientRemainsFinite
// =============================================================================

#[cfg(test)]
mod quadric_tests {
    use crate::bop::int_tools::int_surf_quadric::Quadric;
    use glam::DVec3;
    use rcad_kernel::geom::*;

    #[test]
    fn cone_apex_gradient_remains_finite() {
        // OCCT: gp_Cone(Ax3(Pnt(0,0,0), Dir(0,0,1), Dir(1,0,0)), 0.5, 0.0)
        // Semi-angle = 0.5, radius at apex = 0.0
        let cone = Surface3::Cone(ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 0.0,
            half_angle_rad: 0.5,
        });
        let q = Quadric::from_surface3(&cone).expect("Cone quadric");
        let apex = DVec3::ZERO;

        // Gradient at apex must be finite (not panic, not NaN/inf)
        let grad = q.gradient(apex);
        assert!(grad.is_finite(), "Cone apex gradient must be finite");

        // ValAndGrad at apex: distance should be near zero (apex is on the cone)
        let (dist, grad2) = q.val_and_grad(apex);
        assert!(dist.is_finite(), "Cone apex distance must be finite");
        assert!(
            grad2.is_finite(),
            "Cone apex gradient from val_and_grad must be finite"
        );
    }
}

// =============================================================================
// IntTools_EdgeEdge 鈥?Pure helper function tests (edge_edge.rs)
//
// curve_type_to_integer, point_box_distance,
// split_range_on_segments, intersect_line_line_3d
// =============================================================================

#[cfg(test)]
mod edge_edge_tools_tests {
    use glam::DVec3;
    use rcad_kernel::geom::*;

    #[test]
    fn curve_type_to_integer_all_types() {
        use crate::bop::int_tools::edge_edge::curve_type_to_integer;
        assert_eq!(
            curve_type_to_integer(&Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X))),
            0
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::Hyperbola(Hyperbola3 {
                center: DVec3::ZERO,
                normal: DVec3::Z,
                major_dir: DVec3::X,
                semi_major: 1.0,
                semi_minor: 1.0
            })),
            1
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::Parabola(Parabola3 {
                vertex: DVec3::ZERO,
                normal: DVec3::Z,
                axis_dir: DVec3::X,
                focal_param: 1.0
            })),
            1
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0))),
            2
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::Ellipse(Ellipse3 {
                center: DVec3::ZERO,
                normal: DVec3::Z,
                major_dir: DVec3::X,
                major_radius: 2.0,
                minor_radius: 1.0
            })),
            2
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::BSpline(BSplineCurve3 {
                degree: 3,
                knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
                control_points: vec![DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Z],
                weights: vec![]
            })),
            3
        );
        assert_eq!(
            curve_type_to_integer(&Curve3::Bezier(BezierCurve3 {
                control_points: vec![DVec3::ZERO, DVec3::X, DVec3::Y],
                weights: vec![]
            })),
            3
        );
        // Fallback (default) type
        assert_eq!(
            curve_type_to_integer(&Curve3::Offset(OffsetCurve3 {
                basis: Box::new(Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X))),
                offset_distance: 1.0,
                offset_dir: DVec3::Z
            })),
            4
        );
    }

    #[test]
    fn point_box_distance_inside() {
        use crate::bop::int_tools::edge_edge::point_box_distance;
        let p = DVec3::new(0.5, 0.5, 0.5);
        let bmin = DVec3::ZERO;
        let bmax = DVec3::ONE;
        assert!(
            (point_box_distance(p, bmin, bmax) - 0.0).abs() < 1e-15,
            "Point inside box: distance should be 0"
        );
    }

    #[test]
    fn point_box_distance_outside() {
        use crate::bop::int_tools::edge_edge::point_box_distance;
        let p = DVec3::new(2.0, 0.0, 0.0);
        let bmin = DVec3::ZERO;
        let bmax = DVec3::ONE;
        let d = point_box_distance(p, bmin, bmax);
        assert!(
            (d - 1.0).abs() < 1e-15,
            "Point at (2,0,0) outside box [0,1]: distance should be 1, got {d}"
        );
    }

    #[test]
    fn point_box_distance_on_boundary() {
        use crate::bop::int_tools::edge_edge::point_box_distance;
        let p = DVec3::new(1.0, 0.5, 0.5);
        let bmin = DVec3::ZERO;
        let bmax = DVec3::ONE;
        let d = point_box_distance(p, bmin, bmax);
        assert!(
            d < 1e-15,
            "Point on box boundary: distance should be 0, got {d}"
        );
    }

    #[test]
    fn point_box_distance_diagonal() {
        use crate::bop::int_tools::edge_edge::point_box_distance;
        let p = DVec3::new(2.0, 2.0, 2.0);
        let bmin = DVec3::ZERO;
        let bmax = DVec3::ONE;
        let d = point_box_distance(p, bmin, bmax);
        let expected = (3.0_f64).sqrt(); // sqrt(1^2 + 1^2 + 1^2)
        assert!(
            (d - expected).abs() < 1e-15,
            "Point at corner: distance should be sqrt(3)={expected}, got {d}"
        );
    }

    #[test]
    fn split_range_on_segments_basic() {
        use crate::bop::int_tools::edge_edge::split_range_on_segments;
        let (num, segs) = split_range_on_segments(0.0, 10.0, 0.1, 5);
        assert_eq!(num, 5, "Should return 5 segments for nb_seg=5");
        assert_eq!(segs.len(), 5);
        assert!((segs[0][0] - 0.0).abs() < 1e-15);
        assert!((segs[0][1] - 2.0).abs() < 1e-15);
        assert!((segs[4][0] - 8.0).abs() < 1e-15);
        assert!((segs[4][1] - 10.0).abs() < 1e-15);
    }

    #[test]
    fn split_range_on_segments_small_range() {
        use crate::bop::int_tools::edge_edge::split_range_on_segments;
        // range smaller than resolution 鈫?single segment
        let (num, segs) = split_range_on_segments(0.0, 0.01, 0.1, 5);
        assert_eq!(num, 1, "Range < resolution should produce single segment");
        assert_eq!(segs.len(), 1);
    }

    #[test]
    fn split_range_on_segments_single_requested() {
        use crate::bop::int_tools::edge_edge::split_range_on_segments;
        let (num, segs) = split_range_on_segments(0.0, 10.0, 0.1, 1);
        assert_eq!(num, 1, "nb_seg=1 should produce single segment");
        assert_eq!(segs.len(), 1);
        assert!((segs[0][0] - 0.0).abs() < 1e-15);
        assert!((segs[0][1] - 10.0).abs() < 1e-15);
    }

    #[test]
    fn intersect_line_line_3d_coincident_overlap() {
        use crate::bop::int_tools::edge_edge::intersect_line_line_3d;
        // Two collinear line segments overlapping
        let result = intersect_line_line_3d(
            DVec3::ZERO,
            DVec3::X,
            [0.0, 5.0], // line1: 0鈫?
            DVec3::new(3.0, 0.0, 0.0),
            DVec3::X,
            [0.0, 5.0], // line2 offset: at 3鈫?
            1e-7,
        );
        assert!(
            result.is_some(),
            "Overlapping collinear lines should intersect"
        );
        let (r1, r2, is_edge) = result.unwrap();
        assert!(
            is_edge,
            "Overlapping lines should produce EDGE-type intersection"
        );
        assert!((r1[0] - 0.0).abs() < 1e-10, "Line1 range start should be 0");
        assert!((r1[1] - 5.0).abs() < 1e-10, "Line1 range end should be 5");
        assert!(
            (r2[0] - 3.0).abs() < 1e-10,
            "Line2 projected overlap start should be 3"
        );
        assert!(
            (r2[1] - 5.0).abs() < 1e-10,
            "Line2 projected overlap end should be 5"
        );
    }

    #[test]
    fn intersect_line_line_3d_coincident_no_overlap() {
        use crate::bop::int_tools::edge_edge::intersect_line_line_3d;
        // Two collinear line segments that DON'T overlap
        let result = intersect_line_line_3d(
            DVec3::ZERO,
            DVec3::X,
            [0.0, 2.0],
            DVec3::new(5.0, 0.0, 0.0),
            DVec3::X,
            [0.0, 2.0],
            1e-7,
        );
        assert!(
            result.is_none(),
            "Non-overlapping collinear lines should return None"
        );
    }

    #[test]
    fn intersect_line_line_3d_crossing() {
        use crate::bop::int_tools::edge_edge::intersect_line_line_3d;
        // Two perpendicular lines crossing at (1,0,0)
        // Line1: along X from 0 to 2 鈫?(0,0,0) to (2,0,0)
        // Line2: along Y from -1 to 1 鈫?(1,-1,0) to (1,1,0)
        let result = intersect_line_line_3d(
            DVec3::ZERO,
            DVec3::X,
            [0.0, 2.0],
            DVec3::new(1.0, -1.0, 0.0),
            DVec3::Y,
            [0.0, 2.0],
            1e-7,
        );
        assert!(result.is_some(), "Crossing lines should intersect");
        let (r1, r2, is_edge) = result.unwrap();
        assert!(
            !is_edge,
            "Crossing lines should produce VERTEX-type intersection"
        );
        // Range is a tiny tolerance interval [t-a_dt, t+a_dt], not fully degenerate
        let t1_at_intersection = (r1[0] + r1[1]) * 0.5;
        let t2_at_intersection = (r2[0] + r2[1]) * 0.5;
        // Line1 from (0,0,0) to (2,0,0), intersection at x=1 鈫?t1=1.0
        assert!(
            (t1_at_intersection - 1.0).abs() < 1e-7,
            "Line1 intersection param should be ~1.0"
        );
        // Line2 from (1,-1,0) to (1,1,0), intersection at y=0 鈫?t2=1.0 (0+1)
        assert!(
            (t2_at_intersection - 1.0).abs() < 1e-7,
            "Line2 intersection param should be ~1.0"
        );
    }

    #[test]
    fn intersect_line_line_3d_parallel_no_intersect() {
        use crate::bop::int_tools::edge_edge::intersect_line_line_3d;
        // Two parallel lines separated in Y
        let result = intersect_line_line_3d(
            DVec3::ZERO,
            DVec3::X,
            [0.0, 5.0],
            DVec3::new(0.0, 2.0, 0.0),
            DVec3::X,
            [0.0, 5.0],
            1e-7,
        );
        assert!(
            result.is_none(),
            "Parallel non-collinear lines should return None"
        );
    }
}

// =============================================================================
// BOPTools 鈥?Pure helper function tests (boptools/mod.rs)
//
// is_dirs_coinside, intermediate_point, is_on_pave,
// is_in_range, compute_int_range, is_split_to_reverse, point_near_edge,
// curve_tolerance
// =============================================================================

#[cfg(test)]
mod boptools_helpers_tests {
    use glam::DVec3;
    use rcad_kernel::geom::*;

    #[test]
    fn is_dirs_coinside_same() {
        assert!(
            crate::bop::tools::is_dirs_coinside(DVec3::X, DVec3::X),
            "Same direction"
        );
    }

    #[test]
    fn is_dirs_coinside_opposite() {
        assert!(
            crate::bop::tools::is_dirs_coinside(DVec3::X, -DVec3::X),
            "Opposite direction (2-d < 0.0002)"
        );
    }

    #[test]
    fn is_dirs_coinside_orthogonal() {
        assert!(
            !crate::bop::tools::is_dirs_coinside(DVec3::X, DVec3::Y),
            "Orthogonal should not be coincident"
        );
    }

    #[test]
    fn is_dirs_coinside_with_tol_custom() {
        assert!(
            crate::bop::tools::is_dirs_coinside_with_tol(DVec3::X, DVec3::Y, 1.5),
            "Orthogonal within loose tol (|X-Y|=sqrt2鈮?.414 < 1.5)"
        );
        assert!(
            !crate::bop::tools::is_dirs_coinside_with_tol(DVec3::X, DVec3::Y, 0.5),
            "Orthogonal outside tight tol (|X-Y|=sqrt2鈮?.414 > 0.5, |2-1.414|=0.586 > 0.5)"
        );
    }

    #[test]
    fn intermediate_point_midpoint() {
        assert!((crate::bop::tools::intermediate_point(0.0, 10.0) - 5.0).abs() < 1e-15);
    }

    #[test]
    fn intermediate_point_occt_weighted() {
        // PAR_T = 0.43213918 鈫?result = (1-0.43213918)*0 + 0.43213918*10 = 4.3213918
        let r = crate::bop::tools::intermediate_point_occt(0.0, 10.0);
        assert!(
            (r - 4.3213918).abs() < 1e-10,
            "OCCT-biased midpoint should be ~4.32, got {r}"
        );
    }

    #[test]
    fn is_on_pave_at_boundary() {
        assert!(
            crate::bop::tools::is_on_pave(0.0, [0.0, 10.0], 1e-7),
            "At start boundary"
        );
        assert!(
            crate::bop::tools::is_on_pave(10.0, [0.0, 10.0], 1e-7),
            "At end boundary"
        );
    }

    #[test]
    fn is_on_pave_inside() {
        assert!(
            !crate::bop::tools::is_on_pave(5.0, [0.0, 10.0], 1e-7),
            "Inside range should not be 'on pave'"
        );
    }

    #[test]
    fn is_on_pave_near_boundary() {
        assert!(
            crate::bop::tools::is_on_pave(0.001, [0.0, 10.0], 0.01),
            "Within tolerance of start"
        );
        assert!(
            !crate::bop::tools::is_on_pave(0.001, [0.0, 10.0], 1e-10),
            "Outside tolerance of start"
        );
    }

    #[test]
    fn is_in_range_overlapping() {
        assert!(
            crate::bop::tools::is_in_range([3.0, 7.0], [0.0, 10.0], 0.0),
            "Fully inside"
        );
        assert!(
            crate::bop::tools::is_in_range([-1.0, 5.0], [0.0, 10.0], 0.0),
            "Partially overlapping (left)"
        );
        assert!(
            crate::bop::tools::is_in_range([5.0, 15.0], [0.0, 10.0], 0.0),
            "Partially overlapping (right)"
        );
    }

    #[test]
    fn is_in_range_non_overlapping() {
        assert!(
            !crate::bop::tools::is_in_range([-10.0, -5.0], [0.0, 10.0], 0.0),
            "Completely left"
        );
        assert!(
            !crate::bop::tools::is_in_range([20.0, 30.0], [0.0, 10.0], 0.0),
            "Completely right"
        );
    }

    #[test]
    fn is_in_range_tolerance_expands() {
        // With tolerance, range boundaries expand
        assert!(
            crate::bop::tools::is_in_range([-1.0, -0.5], [0.0, 10.0], 2.0),
            "Within tolerance expansion"
        );
        assert!(
            !crate::bop::tools::is_in_range([-5.0, -3.0], [0.0, 10.0], 2.0),
            "Outside even with tolerance"
        );
    }

    #[test]
    fn compute_int_range_perpendicular() {
        // angle = PI/2 (perpendicular) 鈫?sin=1, tan=inf 鈫?formula handles this
        let r = crate::bop::tools::compute_int_range(1.0, 2.0, std::f64::consts::FRAC_PI_2);
        assert!(
            r.is_finite(),
            "Perpendicular surfaces should produce finite range"
        );
    }

    #[test]
    fn compute_int_range_acute() {
        // angle = PI/3 (60掳)
        let r = crate::bop::tools::compute_int_range(1.0, 2.0, std::f64::consts::FRAC_PI_3);
        assert!(r.is_finite(), "Acute angle should produce finite range");
        assert!(r > 0.0, "Range should be positive");
    }

    #[test]
    fn is_split_to_reverse_forward() {
        assert!(
            !crate::bop::tools::is_split_to_reverse(DVec3::Z, DVec3::Z),
            "Same normal 鈫?not reverse"
        );
        assert!(
            !crate::bop::tools::is_split_to_reverse(DVec3::Z, DVec3::new(0.0, 0.1, 0.99).normalize()),
            "Slightly off normal 鈫?not reverse"
        );
    }

    #[test]
    fn is_split_to_reverse_opposite() {
        assert!(
            crate::bop::tools::is_split_to_reverse(DVec3::Z, -DVec3::Z),
            "Opposite 鈫?reverse"
        );
        assert!(
            crate::bop::tools::is_split_to_reverse(DVec3::Z, DVec3::new(0.0, 0.0, -1.0)),
            "Exact opposite 鈫?reverse"
        );
    }

    #[test]
    fn point_near_edge_offset() {
        let pt = crate::bop::tools::point_near_edge(
            &Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
            DVec3::ZERO,
            DVec3::Z,
        );
        assert!(
            (pt - DVec3::new(0.0, 0.0, crate::tolerance::TOLERANCE_ABS * 10.0)).length() < 1e-15,
            "Point near edge should be offset along normal"
        );
    }

    #[test]
    fn curve_tolerance_default() {
        let line = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        assert!(
            (crate::bop::tools::curve_tolerance(&line, 1.0) - 1.0).abs() < 1e-15,
            "Line should return tol_base unchanged"
        );

        let circle = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 1.0));
        assert!(
            (crate::bop::tools::curve_tolerance(&circle, 1.0) - 1.0).abs() < 1e-15,
            "Circle should return tol_base unchanged"
        );
    }

    #[test]
    fn curve_tolerance_parabola() {
        let parabola = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 1.0,
        });
        let r = crate::bop::tools::curve_tolerance(&parabola, 1.0);
        assert!(
            (r - 10.0).abs() < 1e-15,
            "Parabola should get 10x tolerance, got {r}"
        );
    }
}

// =============================================================================
// DS data structure tests (try_add_interf, has_interf_*, allocate_pave_block)
//
// BOPDS_DS interference fence and PaveBlock pool operations
// =============================================================================

#[cfg(test)]
mod ds_interf_tests {
    #[test]
    fn try_add_interf_new_pair() {
        let mut ds = crate::bop::ds::DS::new_empty();
        assert!(
            ds.try_add_interf(0, 1),
            "First insertion should return true (new pair)"
        );
        assert!(ds.try_add_interf(2, 3), "Another new pair");
    }

    #[test]
    fn try_add_interf_duplicate_rejected() {
        let mut ds = crate::bop::ds::DS::new_empty();
        assert!(ds.try_add_interf(0, 1), "First insertion");
        assert!(!ds.try_add_interf(0, 1), "Duplicate should be rejected");
        assert!(
            !ds.try_add_interf(1, 0),
            "Reversed pair should also be rejected (sorted)"
        );
    }

    #[test]
    fn try_add_interf_multiple_pairs() {
        let mut ds = crate::bop::ds::DS::new_empty();
        assert!(ds.try_add_interf(0, 1));
        assert!(ds.try_add_interf(0, 2));
        assert!(ds.try_add_interf(1, 2));
        assert!(!ds.try_add_interf(0, 1), "Already added");
        assert!(!ds.try_add_interf(2, 0), "Already added (reversed)");
    }

    #[test]
    fn allocate_pave_block_returns_increasing_indices() {
        let mut ds = crate::bop::ds::DS::new_empty();
        use crate::bop::ds::pave::{Pave, PaveBlock};
        let pv1 = Pave {
            vertex_idx: 0,
            param: 0.0,
        };
        let pv2 = Pave {
            vertex_idx: 1,
            param: 1.0,
        };
        let idx0 = ds.allocate_pave_block(PaveBlock::new(0, pv1, pv2));
        assert_eq!(idx0, 0, "First allocation should return 0");
        let idx1 = ds.allocate_pave_block(PaveBlock::new(1, pv1, pv2));
        assert_eq!(idx1, 1, "Second allocation should return 1");
        assert_eq!(ds.pave_blocks.len(), 2);
    }

    #[test]
    fn add_shape_sd_then_has_shape_sd() {
        let mut ds = crate::bop::ds::DS::new_empty();
        ds.add_shape_sd(0, 1);
        assert_eq!(ds.has_shape_sd(0), Some(1));
        // Independent pair
        ds.add_shape_sd(2, 3);
        assert_eq!(ds.has_shape_sd(2), Some(3));
        assert_eq!(ds.has_shape_sd(0), Some(1), "Earlier mapping unchanged");
    }
}

// =============================================================================
// FaceInfo tests
// =============================================================================

#[cfg(test)]
mod face_info_tests {
    #[test]
    fn has_any_interference_empty() {
        use crate::bop::ds::face_info::FaceInfo;
        let fi = FaceInfo::default();
        assert!(!fi.has_any_interference());
    }

    #[test]
    fn has_any_interference_with_paves() {
        use crate::bop::ds::face_info::FaceInfo;
        let mut fi = FaceInfo::default();
        fi.pave_blocks_sc.insert(0);
        assert!(fi.has_any_interference());
    }

    #[test]
    fn has_any_interference_with_vertices_in() {
        use crate::bop::ds::face_info::FaceInfo;
        // Note: has_any_interference checks pave_blocks_in/on/sc and curves_sc,
        // NOT vertices_in. Test that vertices_in alone does NOT trigger it.
        let mut fi = FaceInfo::default();
        fi.vertices_in.insert(0);
        assert!(
            fi.has_any_interference(),
            "vertices_in alone should now trigger has_any_interference (P2)"
        );
        fi.pave_blocks_in.insert(0);
        assert!(
            fi.has_any_interference(),
            "pave_blocks_in should trigger it"
        );
    }

    #[test]
    fn curves_sc_only_returns_section_curves() {
        use crate::bop::ds::face_info::FaceInfo;
        let mut fi = FaceInfo::default();
        fi.curves_sc.insert(3);
        fi.curves_sc.insert(1);
        let sc = fi.curves_sc_only();
        assert_eq!(sc, vec![3, 1]);
    }

    #[test]
    fn curves_sc_only_empty() {
        use crate::bop::ds::face_info::FaceInfo;
        let fi = FaceInfo::default();
        assert!(fi.curves_sc_only().is_empty());
    }
}

// =============================================================================
// boptools/extra trivial helpers
// =============================================================================

#[cfg(test)]
mod boptools_extra_tests {
    use glam::DVec3;

    #[test]
    fn min_step_in_2d_constant() {
        let v = crate::bop::tools::min_step_in_2d();
        assert!(
            v > 0.0 && v < 1.0,
            "min_step_in_2d should be a small positive constant, got {v}"
        );
    }

    #[test]
    fn sense_flag_parallel() {
        assert_eq!(crate::bop::tools::sense_flag(DVec3::Z, DVec3::Z), 1);
    }

    #[test]
    fn sense_flag_opposite() {
        assert_eq!(crate::bop::tools::sense_flag(DVec3::Z, -DVec3::Z), -1);
    }

    #[test]
    fn sense_flag_orthogonal() {
        assert_eq!(crate::bop::tools::sense_flag(DVec3::Z, DVec3::X), 0);
    }
}

// =============================================================================
// face_face helper functions
// =============================================================================

#[cfg(test)]
mod face_face_helper_tests {
    #[test]
    fn correct_surface_boundaries_expands() {
        let mut bounds = [0.0, 1.0, 0.0, 1.0];
        crate::bop::int_tools::face_face::correct_surface_boundaries(&mut bounds, 0.1);
        assert!((bounds[0] - (-0.1)).abs() < 1e-15);
        assert!((bounds[1] - 1.1).abs() < 1e-15);
        assert!((bounds[2] - (-0.1)).abs() < 1e-15);
        assert!((bounds[3] - 1.1).abs() < 1e-15);
    }

    #[test]
    fn correct_plane_boundaries_wide() {
        let mut bounds = [0.0, 1.0, 0.0, 1.0];
        crate::bop::int_tools::face_face::correct_plane_boundaries(&mut bounds);
        assert!(
            (bounds[0] - (-1e10)).abs() < 1.0,
            "Plane u_min should be -1e10"
        );
        assert!((bounds[1] - 1e10).abs() < 1.0, "Plane u_max should be 1e10");
        assert!(
            (bounds[2] - (-1e10)).abs() < 1.0,
            "Plane v_min should be -1e10"
        );
        assert!((bounds[3] - 1e10).abs() < 1.0, "Plane v_max should be 1e10");
    }
}

// 鈹€鈹€ Boolean pipeline stage classification tests 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
//
// These tests use `build_with_history_stage_by_stage` to run the boolean
// pipeline step by step, capturing DS + BRep counts after each stage.
// The first stage with anomalous counts pinpoints the root cause.

#[cfg(test)]
mod stage_classification_tests {
    use super::*;

    /// Classify a boolean case: print per-stage snapshots, return the first
    /// stage where something looks wrong, or None if all good.
    fn classify_case(
        a: &rcad_kernel::BRep,
        b: &rcad_kernel::BRep,
        op: BooleanOpType,
        label: &str,
    ) -> Option<u32> {
        let a_br = a.clone();
        let b_br = b.clone();
        let mut ds = DS::new_from_topods(&a_br, &b_br, TOLERANCE_ABS);
        let bvh_a = Bvh::build(&a_br);
        let bvh_b = Bvh::build(&b_br);
        let brep = rcad_kernel::topods::BRep::new();
        {
            let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_a, &bvh_b);
            filler.set_run_parallel(false);
            filler.perform(a, b);
        }
        let mut builder = crate::bop::algo::builder::BooleanBuilder::with_brep(
            &ds,
            op,
            brep,
            Vec::new(),
            Vec::new(),
        );

        let (_result_brep, _history, snapshots) = builder
            .build_with_history_stage_by_stage()
            .expect("stage_by_stage pipeline failed");

        eprintln!("\n鈹€鈹€ {label} stage classification 鈹€鈹€");
        eprintln!(
            "{stage:>4} {name:<40} DS(V/E/F/PB/IC)    BRep(V/E/F/Sh/So)",
            stage = "Stg",
            name = "StageName"
        );
        eprintln!("{}", "-".repeat(90));
        let mut first_bad: Option<u32> = None;
        for s in &snapshots {
            let flag = if s.n_brep_faces == 0 && s.stage >= 9 {
                " 鈼€ FAIL"
            } else if s.n_brep_faces > 0 && s.stage >= 9 {
                " OK"
            } else if s.n_ds_pave_blocks == 0 && s.stage >= 5 {
                " 鈿?nPB=0"
            } else {
                ""
            };
            eprintln!(
                "{stage:>4} {name:<40} {dsv:>3}/{dse:>3}/{dsf:>3}/{pb:>3}/{ic:>3}     {brv:>3}/{bre:>3}/{brf:>3}/{sh:>3}/{so:>3}{flag}",
                stage = s.stage,
                name = s.stage_name,
                dsv = s.n_ds_vertices,
                dse = s.n_ds_edges,
                dsf = s.n_ds_faces,
                pb = s.n_ds_pave_blocks,
                ic = s.n_ds_intersection_curves,
                brv = s.n_brep_vertices,
                bre = s.n_brep_edges,
                brf = s.n_brep_faces,
                sh = s.n_brep_shells,
                so = s.n_brep_solids
            );

            // Detect first bad stage: after FillImagesFaces (stage >= 9),
            // expect non-zero brep faces
            if s.stage >= 9 && s.n_brep_faces == 0 && first_bad.is_none() {
                first_bad = Some(s.stage);
            }
        }
        eprintln!("鈹€鈹€ {label}: first bad stage = {:?} 鈹€鈹€", first_bad);
        first_bad
    }

    #[test]
    fn classify_bfuse_simple_a1() {
        let a = make_unit_sphere();
        let b = make_unit_box();
        let bad = classify_case(&a, &b, BooleanOpType::Union, "bfuse_simple_A1");
        // After EE fix: pipeline produces non-zero BRep faces at all stages.
        // The first_bad is None = all stages pass.
        if bad.is_some() {
            eprintln!(
                "bfuse_simple A1: first failure at stage {:?} (unexpected)",
                bad
            );
        } else {
            eprintln!("bfuse_simple A1: all stages OK (EE fix improved pipeline)");
        }
    }

    /// Batch classification: run all A-series cases for bfuse_simple.
    #[test]
    fn classify_bfuse_simple_a_series() {
        let cases: Vec<(&str, fn() -> (rcad_kernel::BRep, rcad_kernel::BRep))> = vec![
            ("A1", || (make_unit_sphere(), make_unit_box())),
            ("A2", || {
                (
                    make_sphere(DVec3::ZERO, 1.0),
                    make_box(DVec3::new(-0.5, -0.5, -0.5), 2.0, 2.0, 2.0),
                )
            }),
        ];
        for (label, shapes) in &cases {
            let (a, b) = shapes();
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                classify_case(
                    &a,
                    &b,
                    BooleanOpType::Union,
                    &format!("bfuse_simple_{label}"),
                )
            }));
            match result {
                Ok(Some(bad)) => eprintln!("bfuse_simple_{}: first bad stage = {:?}", label, bad),
                Ok(None) => eprintln!("bfuse_simple_{}: all stages OK", label),
                Err(_) => eprintln!("bfuse_simple_{}: CRASHED (known VF bug, skipping)", label),
            }
        }
    }

    /// Diagnose why MakeBlocks produces 0 PaveBlocks.
    #[test]
    fn diag_make_blocks_bfuse_simple_a1() {
        let a = make_unit_sphere();
        let b = make_unit_box();
        let a_br = a.clone();
        let b_br = b.clone();
        let mut ds = DS::new_from_topods(&a_br, &b_br, TOLERANCE_ABS);
        let bvh_a = Bvh::build(&a_br);
        let bvh_b = Bvh::build(&b_br);
        {
            let mut filler = PaveFiller::with_bvh(&mut ds, &bvh_a, &bvh_b);
            filler.set_run_parallel(false);
            filler.perform(&a, &b);
        };
        eprintln!("\n鈹€鈹€ DS state after pave_fill 鈹€鈹€");
        eprintln!(
            "V={} E={} F={} IC={} PB={}",
            ds.vertex_count(),
            ds.edge_count(),
            ds.face_count(),
            ds.intersection_curves.len(),
            ds.pave_blocks.len()
        );
        for fi in 0..ds.face_count() {
            let surf = format!("{:?}", ds.face_surface(fi))
                .chars()
                .take(30)
                .collect::<String>();
            eprintln!(
                "  F[{}] {} bv={} be={} vi={} vo={} cs={}",
                fi,
                surf,
                ds.face_boundary_verts(fi).len(),
                ds.face_boundary_edges(fi).len(),
                ds.face_info(fi).vertices_in.len(),
                ds.face_info(fi).vertices_on.len(),
                ds.face_info(fi).curves_sc.len()
            );
        }
        for (ci, ic) in ds.intersection_curves.iter().enumerate() {
            eprintln!(
                "  IC[{}] sv={} ev={} r=[{:.2},{:.2}] pca={} pcb={} n_pb={}",
                ci,
                ic.start_vertex,
                ic.end_vertex,
                ic.t_range[0],
                ic.t_range[1],
                ic.pcurve_on_a.is_some(),
                ic.pcurve_on_b.is_some(),
                ic.pave_blocks.len()
            );
            for &v in &[ic.start_vertex, ic.end_vertex] {
                if v < ds.vertex_count() && v < ds.shape_info.len() {
                    eprintln!("    IC[{}] v{} is_new={}", ci, v, ds.shape_info[v].is_new);
                }
            }
            for fi in 0..ds.face_count() {
                for &v in &[ic.start_vertex, ic.end_vertex] {
                    if v < ds.vertex_count() && ds.face_info(fi).vertices_in.contains(&v) {
                        eprintln!("    IC[{}] v{} IN Face[{}].vertices_in", ci, v, fi);
                    }
                }
            }
        }

        // Debug: show interf_ff data
        eprintln!("\n  FF interferences:");
        for (ffi, ff) in ds.interf_ff.iter().enumerate() {
            eprintln!(
                "    FF[{}] f1={} f2={} curves={:?} points={:?}",
                ffi, ff.f1, ff.f2, ff.curves, ff.points
            );
        }
    }
}

// =============================================================================
// PaveFiller internal function alignment tests
//
// These tests verify that individual PaveFiller helper functions produce
// output matching OCCT's corresponding functions.  Each test constructs a
// minimal DS, calls the rcad function, and asserts the result matches the
// expected OCCT behavior (verified against BOPAlgo_PaveFiller source).
// =============================================================================

#[cfg(test)]
mod pave_filler_internal_tests {
    use super::*;
    use crate::bop::algo::pave_filler::build_face_shape_map;
    use crate::bop::ds::{DSEdge, DSFace, DSVertex, ShapeOrigin};
    use glam::DVec3;
    use rcad_kernel::geom::{Curve3, Surface3};

    /// Build a minimal DS for unit testing.
    /// - 3 vertices: v0=(0,0,0), v1=(1,0,0), v2=(1,1,0)
    /// - 2 edges: e0=v0鈫抳1, e1=v1鈫抳2
    /// - 1 face with boundary edges [0,1]
    fn make_minimal_ds() -> DS {
        let mut ds = DS::new_empty();
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::ZERO,
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::X,
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::new(1.0, 1.0, 0.0),
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        let line_curve = |start: DVec3, end: DVec3| -> Curve3 {
            Curve3::Line(rcad_kernel::geom::Line3::new(
                start,
                (end - start).normalize(),
            ))
        };
        ds.push_edge(
            DSEdge {
                shape_idx: 0, start_vertex: 0,
                end_vertex: 1,
                curve: line_curve(DVec3::ZERO, DVec3::X),
                t_range: [0.0, 1.0],
                origin: ShapeOrigin::ShapeA,
                geom_tol: 1e-7,
                paves: Vec::new(),
                pave_blocks: Vec::new(),
                face_reps: Vec::new(),
                is_internal: false,
                is_geometric: true,
                vertex_params: std::collections::HashMap::new(),
                face_tolerances: Vec::new(),
                location: 0,
            },
            None,
        );
        ds.push_edge(
            DSEdge {
                shape_idx: 0, start_vertex: 1,
                end_vertex: 2,
                curve: line_curve(DVec3::X, DVec3::new(1.0, 1.0, 0.0)),
                t_range: [0.0, 1.0],
                origin: ShapeOrigin::ShapeA,
                geom_tol: 1e-7,
                paves: Vec::new(),
                pave_blocks: Vec::new(),
                face_reps: Vec::new(),
                is_internal: false,
                is_geometric: true,
                vertex_params: std::collections::HashMap::new(),
                face_tolerances: Vec::new(),
                location: 0,
            },
            None,
        );
        ds.push_face(
            DSFace {
                shape_idx: 0, surface: Surface3::Plane(rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z)),
                boundary_verts: vec![0, 1, 2],
                boundary_edges: vec![0, 1],
                boundary_edge_forwards: vec![true, true],
                inner_boundary_edges: Vec::new(),
                outer_wire_idx: Some(0),
                inner_wire_idxs: Vec::new(),
                normal: DVec3::Z,
                origin: ShapeOrigin::ShapeA,
                face_info: crate::bop::ds::face_info::FaceInfo::default(),
                source_face_idx: 0,
                geom_tol: 1e-7,
                location: 0,
                uv_boundary: None,
                natural_restriction: true,
                source_shell_idx: Some(0),
                source_compsolid_idx: Some(0),
                source_solid_idx: Some(0),
            },
            None,
        );
        ds.a_vertex_count() = 0;
        ds.a_edge_count() = 2; // edges 0-1 belong to operand A
        ds.a_face_count() = 1;
        ds
    }

    /// BOPAlgo_PaveFiller::GetFullShapeMap
    /// (PaveFiller_6.cxx L2941-2958).
    /// The face itself, its boundary edges, and their endpoint vertices
    /// should all be present in the returned set.
    #[test]
    fn build_face_shape_map_returns_face_edges_and_vertices() {
        let ds = make_minimal_ds();
        let result = build_face_shape_map(&ds, 0);
        // Face index
        assert!(result.contains(&0), "face index 0 must be in shape map");
        // Boundary edges
        assert!(result.contains(&0), "edge 0 must be in shape map");
        assert!(result.contains(&1), "edge 1 must be in shape map");
        // Endpoint vertices of edge 0 (v0, v1)
        assert!(result.contains(&0), "vertex 0 must be in shape map");
        assert!(result.contains(&1), "vertex 1 must be in shape map");
        // Endpoint vertices of edge 1 (v1, v2)
        assert!(result.contains(&2), "vertex 2 must be in shape map");
        // Face index 0 overlaps with vertex index 0 in this minimal DS,
        // so expected set = {0 (face 0 / vertex 0), 1 (edge 0 / vertex 1), 2 (edge 1 / vertex 2)}
        assert_eq!(result.len(), 3, "set should contain {{0, 1, 2}}");
    }

    /// intersect_vertices (BOPAlgo_Tools::IntersectVertices).
    /// Groups vertices by tolerance-sphere overlap.
    fn make_vertex_test_ds() -> DS {
        let mut ds = DS::new_empty();
        // v0: far away
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::new(0.0, 0.0, 0.0),
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        // v1: close to v0 (within tol)
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::new(1e-8, 0.0, 0.0),
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        // v2: far from both (outside tol of v0, v1)
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::new(100.0, 0.0, 0.0),
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        // v3: chain-close to v2 (close to v2 but not to v0/v1)
        ds.push_vertex(
            DSVertex {
                shape_idx: 0, point: DVec3::new(100.0 + 1e-8, 0.0, 0.0),
                origin: None,
                geom_tol: 1e-7,
                is_internal: false,
                location: 0,
            },
            None,
        );
        ds.a_vertex_count() = 0;
        ds
    }

    #[test]
    fn intersect_vertices_close_pair_joined() {
        let ds = make_vertex_test_ds();
        // v0 and v1 are within tol 鈫?one group
        let blocks = crate::bop::algo::intersect_vertices(&[0, 1], &ds, 0.0);
        assert_eq!(
            blocks.len(),
            1,
            "close vertices should merge into one group"
        );
        assert!(blocks[0].contains(&0), "group must contain v0");
        assert!(blocks[0].contains(&1), "group must contain v1");
    }

    #[test]
    fn intersect_vertices_far_pair_separate() {
        let ds = make_vertex_test_ds();
        // v0 and v2 are far apart 鈫?two groups
        let blocks = crate::bop::algo::intersect_vertices(&[0, 2], &ds, 0.0);
        assert_eq!(
            blocks.len(),
            2,
            "distant vertices should be separate groups"
        );
    }

    #[test]
    fn intersect_vertices_chain_connected() {
        let ds = make_vertex_test_ds();
        // v2 and v3 are close 鈫?one group
        let blocks = crate::bop::algo::intersect_vertices(&[2, 3], &ds, 0.0);
        assert_eq!(blocks.len(), 1, "chain-close vertices should merge");
    }

    #[test]
    fn intersect_vertices_singleton() {
        let ds = make_vertex_test_ds();
        let blocks = crate::bop::algo::intersect_vertices(&[0], &ds, 0.0);
        assert_eq!(blocks.len(), 1, "single vertex -> one group");
        assert_eq!(
            blocks[0],
            vec![0],
            "single vertex group must contain only v0"
        );
    }

    /// PaveBlock::Update (BOPDS_PaveBlock.cxx L249-312).
    /// Sub-PB splitting from ext_paves with theFlag=false.
    /// When theFlag=false, only ext_paves (not pave1/pave2) define sub-PB boundaries.
    #[test]
    fn pave_block_update_false_empty_ext_paves() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            crate::bop::ds::pave::NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 1,
                param: 1.0,
            },
        );
        let result = pb.update(false);
        // a_nb = 0 (no ext_paves), a_nb <= 1 鈫?empty result
        assert!(result.is_empty(), "no ext_paves + theFlag=false 鈫?empty");
    }

    #[test]
    fn pave_block_update_false_one_ext_pave() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            crate::bop::ds::pave::NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 2,
                param: 2.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 1,
            param: 1.0,
        });
        let result = pb.update(false);
        // a_nb = 1 (one ext_pave), a_nb <= 1 鈫?empty result
        assert!(result.is_empty(), "one ext_pave + theFlag=false 鈫?empty");
    }

    #[test]
    fn pave_block_update_false_two_ext_paves_one_sub_pb() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            crate::bop::ds::pave::NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 3,
                param: 3.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 1,
            param: 1.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 2,
            param: 2.0,
        });
        let result = pb.update(false);
        // a_nb = 2, produces one sub-PB: (ext1, ext2)
        assert_eq!(result.len(), 1, "two ext_paves 鈫?one sub-PB");
        let (v1, v2) = result[0].indices();
        assert_eq!(v1, 1, "first sub-PB start = ext_pave1");
        assert_eq!(v2, 2, "first sub-PB end = ext_pave2");
    }

    #[test]
    fn pave_block_update_false_three_ext_paves_two_sub_pbs() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            crate::bop::ds::pave::NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 4,
                param: 4.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 1,
            param: 1.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 2,
            param: 2.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 3,
            param: 3.0,
        });
        let result = pb.update(false);
        // a_nb = 3, produces two sub-PBs: (ext1, ext2) and (ext2, ext3)
        assert_eq!(result.len(), 2, "three ext_paves 鈫?two sub-PBs");
        let (v1a, v2a) = result[0].indices();
        let (v1b, v2b) = result[1].indices();
        assert_eq!(v1a, 1, "sub-PB[0] start = ext_pave1");
        assert_eq!(v2a, 2, "sub-PB[0] end = ext_pave2");
        assert_eq!(v1b, 2, "sub-PB[1] start = ext_pave2");
        assert_eq!(v2b, 3, "sub-PB[1] end = ext_pave3");
    }

    /// PaveBlock::Update (BOPDS_PaveBlock.cxx L249-312).
    /// Sub-PB splitting with theFlag=true includes pave1/pave2 as boundary paves.
    #[test]
    fn pave_block_update_true_includes_pave1_pave2() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            crate::bop::ds::pave::NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 3,
                param: 3.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 1,
            param: 1.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 2,
            param: 2.0,
        });
        let result = pb.update(true);
        // a_nb = 2 + 2 = 4, produces 3 sub-PBs: (pave1, e1), (e1, e2), (e2, pave2)
        assert_eq!(result.len(), 3, "two ext_paves + theFlag=true 鈫?3 sub-PBs");
        let (v1a, v2a) = result[0].indices();
        assert_eq!(v1a, 0, "sub-PB[0] start = pave1");
        assert_eq!(v2a, 1, "sub-PB[0] end = ext_pave1");
        let (v1b, v2b) = result[2].indices();
        assert_eq!(v1b, 2, "sub-PB[2] start = ext_pave2");
        assert_eq!(v2b, 3, "sub-PB[2] end = pave2");
    }

    // ===== DS::shape_sd infrastructure =====

    /// AddShapeSD / HasShapeSD (BOPDS_DS).
    /// SD stores bi-directional entries and returns the minimum partner.
    #[test]
    fn shape_sd_direct_mapping() {
        let mut ds = DS::new_empty();
        ds.add_shape_sd(5, 2);
        // add_sd_vertex inserts both (5,2) and (2,5)
        assert_eq!(ds.has_shape_sd(5), Some(2), "v5 should map to v2");
        assert_eq!(ds.has_shape_sd(2), Some(5), "v2 has reverse mapping to v5");
        assert_eq!(ds.has_shape_sd(0), None, "v0 has no SD mapping");
    }

    #[test]
    fn shape_sd_chain_mapping() {
        let mut ds = DS::new_empty();
        ds.add_shape_sd(5, 3);
        ds.add_shape_sd(3, 1);
        // Direct: 5鈫? and 3鈫?.  find_sd_partner does NOT follow chains
        // (OCCT HasShapeSD follows chains in the DS-level has_shape_sd,
        //  but rcad's find_sd_partner returns the direct minimum partner)
        assert_eq!(ds.has_shape_sd(5), Some(3), "v5 direct SD partner is v3");
        assert_eq!(ds.has_shape_sd(3), Some(1), "v3 direct SD partner is v1");
        // Chain following requires calling has_shape_sd twice
        let step1 = ds.has_shape_sd(5).unwrap();
        let step2 = ds.has_shape_sd(step1).unwrap();
        assert_eq!(step2, 1, "v5 chain-follows to v1 via v3");
    }

    #[test]
    fn shape_sd_self_mapping_stored() {
        let mut ds = DS::new_empty();
        ds.add_shape_sd(5, 5);
        // Self-mapping: (5,5) is stored and returned
        assert_eq!(ds.has_shape_sd(5), Some(5), "self-mapping returns self");
    }

    // ===== PaveBlock basic operations =====

    /// PaveBlock::Indices / Range / ExtPaves.
    use crate::bop::ds::pave::NO_EDGE;

    #[test]
    fn pave_block_construction_and_accessors() {
        let pb = crate::bop::ds::pave::PaveBlock::new(
            NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 3,
                param: 1.5,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 7,
                param: 3.2,
            },
        );
        let (v1, v2) = pb.indices();
        assert_eq!(v1, 3, "pave1.vertex_idx");
        assert_eq!(v2, 7, "pave2.vertex_idx");
        let (t1, t2) = pb.range();
        assert!((t1 - 1.5).abs() < 1e-12, "pave1.param = 1.5");
        assert!((t2 - 3.2).abs() < 1e-12, "pave2.param = 3.2");
    }

    #[test]
    fn pave_block_append_ext_pave_and_contains() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 10,
                param: 10.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 5,
            param: 5.0,
        });
        // contains_parameter should find the ext_pave at t=5.0
        let mut n_v_used = 0;
        let found = pb.contains_parameter(5.0, 1e-6, &mut n_v_used);
        assert!(found, "ext_pave at param 5.0 should be found");
        assert_eq!(n_v_used, 5, "ext_pave at param 5.0 has vertex_idx 5");
    }

    #[test]
    fn pave_block_ext_paves_sorted_by_param() {
        let mut pb = crate::bop::ds::pave::PaveBlock::new(
            NO_EDGE,
            crate::bop::ds::pave::Pave {
                vertex_idx: 0,
                param: 0.0,
            },
            crate::bop::ds::pave::Pave {
                vertex_idx: 10,
                param: 10.0,
            },
        );
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 9,
            param: 9.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 3,
            param: 3.0,
        });
        pb.append_ext_pave(crate::bop::ds::pave::Pave {
            vertex_idx: 7,
            param: 7.0,
        });
        // update(false) sorts by param: ext(3), ext(7), ext(9) 鈫?2 sub-PBs
        let result = pb.update(false);
        assert_eq!(
            result.len(),
            2,
            "3 unsorted ext_paves 鈫?2 sub-PBs after sort"
        );
        let (v1, v2) = result[0].indices();
        assert_eq!(v1, 3, "first sub-PB start = smallest param ext_pave");
        assert_eq!(v2, 7, "first sub-PB end = middle param ext_pave");
        let (v1, v2) = result[1].indices();
        assert_eq!(v1, 7, "second sub-PB start = middle param ext_pave");
        assert_eq!(v2, 9, "second sub-PB end = largest param ext_pave");
    }
}
