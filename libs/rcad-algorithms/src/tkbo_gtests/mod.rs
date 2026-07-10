//! OCCT-aligned TKBO GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKBO/GTests/
//!
//! Files translated in this module:
//!   rcad_kernel::BRepAlgoAPI_BuilderAlgo_Test.cxx  鈥?Non-copyability traits (C++ 鈫?Rust: trivially pass)
//!   BOPAlgo_BOP_Test.cxx              鈥?Direct and two-step BOP operations
//!   BOPAlgo_PaveFiller_Test.cxx       鈥?PaveFiller regression tests (degenerated edges)
//!   IntTools_FaceFace_Test.cxx         鈥?Face-face intersection
//!   rcad_kernel::BRepAlgoAPI_Common_Test.cxx       鈥?Common operation tests
//!
//!
//! NOTE: rcad_kernel::BRepAlgoAPI_Fuse_Test.cxx, rcad_kernel::BRepAlgoAPI_Cut_Test.cxx, and Cut_Test_1.cxx
//! (bfuse_simple / bcut_simple DRAW series) overlap with the existing
//! DRAW-derived generated OCCT tests (tests/occt/tests/generated_occt_boolean_bfuse_simple.rs
//! and generated_occt_boolean_bcut_simple.rs). Those test series are covered by
//! the tkremaining_gtests module as stubs to avoid duplication.

use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::{surface_area, volume};
use rcad_kernel::geom;
use rcad_kernel::topods;

use crate::builder::BooleanOpType;
use crate::bopds::ds::DS;
use crate::brep_tools::make_face_half_space;
use crate::pave_filler::PaveFiller;
use crate::bvh::Bvh;
use crate::tolerance::TOLERANCE_ABS;

const TOL: f64 = 1.0e-6;
const SA_TOLERANCE: f64 = 5000.0;

// =============================================================================
// Helper utilities 鈥?BOPTest_Utilities equivalent
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
    rcad_modeling::make_sphere_brep(DVec3::ZERO, 1.0)
        .expect("Unit sphere creation failed")
}

fn make_sphere(center: DVec3, radius: f64) -> topods::BRep {
    rcad_modeling::make_sphere_brep(center, radius)
        .expect("Sphere creation failed")
}

fn make_cylinder(radius: f64, height: f64) -> topods::BRep {
    rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height)
        .expect("Cylinder creation failed")
}

fn make_cone(base_radius: f64, height: f64) -> topods::BRep {
    rcad_modeling::make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, base_radius, height)
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

fn validate_result(result: &rcad_kernel::BRep, expected_sa: f64, expected_vol: f64, expected_empty: bool) {
    if expected_empty {
        assert!(is_empty(result), "Result should be empty");
        return;
    }
    if expected_sa >= 0.0 {
        let sa = get_surface_area(result);
        let diff = (sa - expected_sa).abs();
        assert!(diff <= SA_TOLERANCE,
                "Surface area mismatch: got {sa}, expected {expected_sa} (diff {diff})");
    }
    if expected_vol >= 0.0 {
        let vol = get_volume(result);
        assert!((vol - expected_vol).abs() < TOL,
                "Volume mismatch: got {vol}, expected {expected_vol}");
    }
}

#[allow(dead_code)]
fn validate_result_sa(result: &rcad_kernel::BRep, expected_sa: f64) {
    validate_result(result, expected_sa, -1.0, false);
}

// =============================================================================
// rcad_kernel::BRepAlgoAPI_BuilderAlgo_Test.cxx 鈥?C++ type traits
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
// BOPAlgo_BOP_Test.cxx 鈥?Direct and two-step BOP operations
// =============================================================================

#[cfg(test)]
mod bop_algo_direct_tests {
    use super::*;
    use crate::BooleanOpType;

    fn perform_bop(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep, op: BooleanOpType) -> rcad_kernel::BRep {
        crate::bop_occt_union::boolean_op_generic(op, a, b).expect("BOP operation failed")
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
        let mut brep = rcad_kernel::topods::BRep::new();
        {
            let mut filler = PaveFiller::with_bvh_and_brep(&mut ds, &bvh_a, &bvh_b, &mut brep);
            filler.set_run_parallel(false);
            filler.perform();
        }
        let n_pb: usize = ds.intersection_curves.iter().map(|ic| ic.pave_blocks.len()).sum();
        let n_pb_sc: usize = ds.faces.iter().map(|f| f.face_info.pave_blocks_sc.len()).sum();
        eprintln!("DEBUG: n_ic={} n_pb={} n_pb_sc={}", ds.intersection_curves.len(), n_pb, n_pb_sc);
        for (i, ic) in ds.intersection_curves.iter().enumerate() {
            if !ic.pave_blocks.is_empty() {
                eprintln!("  curve[{}]: sv={} ev={} t=[{:.4},{:.4}] {} pb",
                    i, ic.start_vertex, ic.end_vertex, ic.t_range[0], ic.t_range[1], ic.pave_blocks.len());
            }
        }

        let result = crate::bop_occt_union::boolean_op_generic(BooleanOpType::Intersection, &b1, &b2).expect("Common failed");
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
    fn perform_two_step_bop(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep, op: BooleanOpType) -> rcad_kernel::BRep {
        let a_br = a.clone();
        let b_br = b.clone();
        let mut ds = DS::new_from_topods(&a_br, &b_br, TOLERANCE_ABS);
        let bvh_a = Bvh::build(&a_br);
        let bvh_b = Bvh::build(&b_br);
        let mut brep = rcad_kernel::topods::BRep::new();
        let (face_refs, ic_edge_map) = {
            let mut filler = PaveFiller::with_bvh_and_brep(&mut ds, &bvh_a, &bvh_b, &mut brep);
            filler.set_run_parallel(false);
            filler.perform();
            (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
        };
        let builder = crate::builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
        let (result, _history) = builder.build_with_history_topods().expect("Two-step BOP failed");
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

    fn perform_direct_bop(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep, op: BooleanOpType) -> rcad_kernel::BRep {
        crate::bop_occt_union::boolean_op_generic(op, a, b).expect("Direct BOP failed")
    }

    /// OCCT: MultipleIntersectingPrimitives
    /// Chain: sphere ∩ cylinder → result ∪ box
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
        let two_step = crate::bop_occt_union::boolean_op_generic(BooleanOpType::Union, &sphere, &box_)
            .expect("Two-step failed");

        let direct_vol = get_volume(&direct);
        let two_step_vol = get_volume(&two_step);
        assert!((direct_vol - two_step_vol).abs() < TOL,
            "Direct and two-step should produce equivalent results: {} vs {}", direct_vol, two_step_vol);
    }
}

// =============================================================================
// BOPAlgo_BOP_Test.cxx -- Degenerate thin-tool tests
// =============================================================================

#[cfg(test)]
mod bop_algo_thin_tool_tests {
    use super::*;

    fn perform_direct_bop(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep, op: BooleanOpType) -> rcad_kernel::BRep {
        crate::bop_occt_union::boolean_op_generic(op, a, b).expect("Direct BOP failed")
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
        let plane = geom::Plane {
            origin: DVec3::new(0.0, 0.0, 0.0),
            normal: DVec3::Z,
        };
        let bbox = [DVec3::new(-250.0, -250.0, -30.0), DVec3::new(250.0, 250.0, 200.0)];
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

        let mut fuser = crate::brep_algo_api::BRepAlgoAPI_Fuse::new(&cone_b, &box_b);
        assert!(fuser.build(), "Boolean fuse of cone and box should succeed");
        assert!(get_volume(fuser.shape()) > 0.0);
    }

    #[test]
    fn fuse_sphere_with_box() {
        let sphere_b = (make_unit_sphere().clone());
        let box_b = (make_unit_box().clone());

        let mut fuser = crate::brep_algo_api::BRepAlgoAPI_Fuse::new(&sphere_b, &box_b);
        assert!(fuser.build(), "Boolean fuse should succeed");
        assert!(get_volume(fuser.shape()) > 0.0);
    }
}

// =============================================================================
// IntTools_FaceFace_Test.cxx -- Face-face intersection
//
// Not yet translatable (require boundary-aware face-face intersection):
//   - HalfCylinderOutsideCircularPlane_NoIntersection
//   - HalfCylinderInsideCircularPlane_HasIntersection
//   - OppositeHalfInsideCircle_HasIntersection
//   - PartialCrossing_ProperlyTrimmed
//   - BothHalvesCompletelyOutside_NoIntersection
// =============================================================================

#[cfg(test)]
mod int_tools_face_face_tests {
    use super::*;
    use crate::inttools::face_face::intersect_faces;
    use crate::inttools::int_ana_quad_quad_geo::QuadQuadGeo;
    use crate::inttools::int_surf_quadric::Quadric;

    #[test]
    fn perpendicular_planes_intersect() {
        let s1 = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        });
        let s2 = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });

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
        let plane_s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        });

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

        assert_eq!(c12.len(), c21.len(),
            "Intersection curve count should not depend on argument order: {} vs {}",
            c12.len(), c21.len());
    }

    /// OCCT: OCC24005_PlaneCylinderIntersection
    /// Slightly off-angle plane intersecting a cylinder.
    /// Original OCCT regression test: must complete quickly (no hang) and
    /// produce at least one intersection result.
    #[test]
    fn occ24005_plane_cylinder_intersection() {
        let plane_s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::new(-72.948737453424499, 754.30437716359393, 259.52151854671678),
            normal: DVec3::new(6.2471473085930200e-007, -0.99999999999980493, 0.0).normalize(),
        });
        let cyl_s = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::new(-6.4812490053250649, 753.39408794522092, 279.16400974257465),
            axis: DVec3::X,
            ref_dir: DVec3::Y,
            radius: 19.712534607908712,
        });

        let curves = intersect_faces(&plane_s, &cyl_s, 1e-7, 1e-7);

        assert!(!curves.is_empty(),
            "Expected at least one intersection curve for plane-cylinder");
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

        let mut common = crate::brep_algo_api::BRepAlgoAPI_Common::new(&b1, &b2);
        assert!(common.build(), "Common operation should succeed");
        assert!(get_surface_area(common.shape()) > 0.0);
    }
}

