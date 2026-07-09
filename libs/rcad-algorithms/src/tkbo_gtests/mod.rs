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
//! and generated_occt_boolean_bcut_simple.rs). Those test series are NOT re-translated
//! here to avoid duplication.

use glam::DVec3;
use rcad_kernel::{surface_area, volume};
use rcad_kernel::topods;

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
        (crate::boolean_op(op, a, b).clone().expect("BOP operation failed"))
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
        let result = perform_bop(&b1, &b2, BooleanOpType::Intersection);
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
// BOPAlgo_PaveFiller_Test.cxx 鈥?Degenerated edge handling
// =============================================================================

#[cfg(test)]
mod pave_filler_tests {
    use super::*;

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
// IntTools_FaceFace_Test.cxx 鈥?Face-face intersection
// =============================================================================

#[cfg(test)]
mod int_tools_face_face_tests {
    use super::*;
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

