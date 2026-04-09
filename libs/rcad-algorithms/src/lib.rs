pub mod bopds;
pub mod brep_check;
pub mod brep_repair;
pub mod builder;
pub mod classify;
pub mod geom_populate;
pub mod history;
pub mod hlr;
pub mod imprint;
pub mod inttools;
pub mod pave_filler;
pub mod section;
pub mod tolerance;
pub mod triangulate;

use rcad_kernel::BRep;

pub use brep_check::{CheckIssue, CheckResult, check};
pub use brep_repair::{
    RepairReport, fix_wire_orientation, merge_close_vertices, recompute_face_normals,
    remove_degenerate_faces, repair,
};
pub use builder::{BooleanError, BooleanOpType};
pub use history::{BooleanHistory, FaceOrigin};
pub use hlr::{HlrCamera, HlrResult, HlrSegment, hlr, hlr_to_svg};
pub use imprint::{
    Gap, GapOverlapReport, ImprintResult, Overlap, detect_gaps_overlaps, imprint_brep,
};
pub use inttools::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces,
};
pub use section::{SectionCurve, section, section_curves, section_polylines};

/// Perform a boolean operation on two BReps.
///
/// Both BReps must have populated GeomStore (call
/// `geom_populate::populate_box_geom` first for box primitives).
pub fn boolean_op(op: BooleanOpType, a: &BRep, b: &BRep) -> Result<BRep, BooleanError> {
    // 1. Build the DS from both shapes
    let mut ds = bopds::ds::DS::new(a, b);

    // 2. Run PaveFiller — compute all interferences
    let mut filler = pave_filler::PaveFiller::new(&mut ds);
    filler.perform();

    // 3. Run Builder — classify and assemble result
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build()
}

/// Like [`boolean_op`] but also returns a [`BooleanHistory`] mapping each result
/// face back to its source in solid A or B.
pub fn boolean_op_with_history(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    let mut ds = bopds::ds::DS::new(a, b);
    let mut filler = pave_filler::PaveFiller::new(&mut ds);
    filler.perform();
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build_with_history()
}

/// Union two BReps and return both the result and face origin history.
pub fn union_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Union, a, b)
}

/// Intersect two BReps and return both the result and face origin history.
pub fn intersection_with_history(
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Intersection, a, b)
}

/// Subtract solid B from solid A and return both the result and face origin history.
pub fn difference_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Difference, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;
    use rcad_modeling::{make_box_brep, make_cylinder_brep, make_sphere_brep};

    fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: w,
            height: h,
            depth: d,
        });
        for v in &mut brep.vertices {
            v.point += DVec3::new(x, y, z);
        }
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    fn face_count(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }

    fn triangle_count(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .map(|f| f.triangles.len())
            .sum()
    }

    fn all_triangles_valid(brep: &BRep) -> bool {
        let nv = brep.vertices.len();
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .flat_map(|f| &f.triangles)
            .all(|tri| tri.iter().all(|&i| i < nv))
    }

    #[test]
    fn union_disjoint_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        // Disjoint: all 12 faces kept
        assert_eq!(face_count(&result), 12);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn intersection_disjoint_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // Disjoint: intersection is empty
        assert!(result.is_err());
    }

    #[test]
    fn union_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn intersection_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert!(face_count(&result) >= 6);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn difference_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn contained_box_difference() {
        // B completely inside A
        let a = box_at(0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn contained_box_intersection() {
        // B completely inside A → intersection is B
        let a = box_at(0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert_eq!(face_count(&result), 6); // B's 6 faces
        assert!(all_triangles_valid(&result));
    }

    // ─── Phase 4 edge case tests ───────────────────────────────────────

    #[test]
    fn identical_boxes_union() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn touching_face_union() {
        // Two boxes sharing a face (A right = B left)
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(1.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn touching_edge_union() {
        // Two boxes sharing an edge
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(1.0, 1.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert_eq!(face_count(&result), 12);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn non_unit_boxes_difference() {
        let a = box_at(0.0, 0.0, 0.0, 3.0, 2.0, 5.0);
        let b = box_at(1.0, 0.5, 1.0, 1.0, 1.0, 3.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 6);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn offset_3d_intersection() {
        let a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert!(face_count(&result) >= 6);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn difference_is_not_symmetric() {
        let a = box_at(0.0, 0.0, 0.0, 2.0, 1.0, 1.0);
        let b = box_at(1.0, 0.0, 0.0, 2.0, 1.0, 1.0);
        let a_minus_b = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        let b_minus_a = boolean_op(BooleanOpType::Difference, &b, &a).unwrap();
        assert!(face_count(&a_minus_b) > 0);
        assert!(face_count(&b_minus_a) > 0);
        assert!(all_triangles_valid(&a_minus_b));
        assert!(all_triangles_valid(&b_minus_a));
    }

    #[test]
    fn small_overlap_union() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.99, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn large_overlap_intersection() {
        let a = box_at(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        let b = box_at(0.1, 0.1, 0.1, 9.8, 9.8, 9.8);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert_eq!(face_count(&result), 6);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn classify_point_on_face() {
        use classify::Classification;
        let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        let ds = bopds::ds::DS::new(&brep, &rcad_kernel::BRep::new());
        let face_indices: Vec<usize> = (0..ds.faces.len())
            .filter(|&i| ds.faces[i].origin == bopds::ds::ShapeOrigin::ShapeA)
            .collect();
        let on_top = DVec3::new(1.0, 2.0, 1.0);
        assert_eq!(
            classify::classify_point(on_top, &face_indices, &ds),
            Classification::On
        );
    }

    #[test]
    fn triangulate_hexagon() {
        use triangulate::triangulate_polygon;
        let verts: Vec<DVec3> = (0..6)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 6.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 4);
        for tri in &tris {
            for &idx in tri {
                assert!(idx < 6);
            }
        }
    }

    // ─── Curved Boolean Tests ──────────────────────────────────────────────────

    #[test]
    fn boolean_box_sphere_intersection() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 1.5).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere intersection failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    fn boolean_box_sphere_difference() {
        // Small sphere inside a box — creates a hole
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b = make_sphere_brep(DVec3::new(2.0, 2.0, 2.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    fn boolean_box_sphere_union() {
        // Sphere protruding from box
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 2.5), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere union failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        let v_box = rcad_kernel::properties::volume(&a);
        let v_sphere = rcad_kernel::properties::volume(&b);
        assert!(v > v_box, "union should be larger than box");
        assert!(v > v_sphere, "union should be larger than sphere");
    }

    #[test]
    #[ignore = "sphere UV boundary degeneracy: the UV boundary [-π,0],[π,0],[π,π],[-π,π] \
                maps to only 2 distinct 3D points (poles), making sub-face boundary vertices \
                degenerate. Requires reparametrization of sphere domain or triangulation-based splitting."]
    fn boolean_sphere_sphere_intersection() {
        // Two overlapping unit spheres
        let a = make_sphere_brep(DVec3::new(-0.5, 0.0, 0.0), 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(0.5, 0.0, 0.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-sphere intersection failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        let v_sphere = rcad_kernel::properties::volume(&a);
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(
            v < v_sphere,
            "intersection should be smaller than one sphere"
        );
    }

    #[test]
    #[ignore = "sphere UV boundary degeneracy: same root cause as boolean_sphere_sphere_intersection; \
                requires reparametrization so the sphere's boundary polygon maps to a non-degenerate \
                3D polygon that can be split by the intersection curve."]
    fn boolean_sphere_sphere_difference() {
        // Large sphere minus small sphere
        let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-sphere difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        let v_large = rcad_kernel::properties::volume(&a);
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(v < v_large, "difference should be smaller than original");
    }

    #[test]
    fn boolean_box_cylinder_hole() {
        // Box minus a cylinder through it (classic hole)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        // Cylinder along Z axis through center of box
        let b =
            make_cylinder_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.5, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-cylinder difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    #[ignore = "cylinder-cylinder intersection: PaveFiller produces no intersection curves for \
                perpendicular cylinders (Steinmetz configuration); the FF pass marching fallback \
                needs Cylinder×Cylinder support. Returns DegenerateResult."]
    fn boolean_cylinder_cylinder_intersection() {
        // Two perpendicular cylinders (Steinmetz solid)
        let a =
            make_cylinder_brep(DVec3::new(0.0, -2.0, 0.0), DVec3::Y, DVec3::X, 1.0, 4.0).unwrap();
        let b =
            make_cylinder_brep(DVec3::new(-2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 4.0).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "cylinder-cylinder intersection failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        assert!(
            v > 0.0,
            "Steinmetz solid volume should be positive, got {v}"
        );
        // The Steinmetz solid of radius r has volume 16/3 * r^3, approx 5.33 for r=1
        assert!(
            (v - 16.0 / 3.0).abs() < 0.5,
            "Steinmetz solid volume should be ~5.33, got {v}"
        );
    }

    #[test]
    #[ignore = "volume conservation requires accurate volume computation for curved surface \
                result faces; the divergence-theorem volume currently returns 0 for sphere \
                faces due to missing triangulation in the primitive BRep"]
    fn volume_conservation_box_sphere() {
        // V(A∪B) ≈ V(A) + V(B) - V(A∩B), error < 5%
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.5), 1.0).unwrap();

        let union_result = boolean_op(BooleanOpType::Union, &a, &b);
        let inter_result = boolean_op(BooleanOpType::Intersection, &a, &b);

        assert!(
            union_result.is_ok(),
            "union failed: {:?}",
            union_result.err()
        );
        assert!(
            inter_result.is_ok(),
            "intersection failed: {:?}",
            inter_result.err()
        );

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_result.unwrap());
        let v_inter = rcad_kernel::properties::volume(&inter_result.unwrap());

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected;
        let error_pct = error * 100.0;
        assert!(
            error < 0.05,
            "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
        );
    }

    #[test]
    #[ignore = "sphere-sphere UV splitting and volume conservation both require more \
                complete implementation of the Boolean pipeline for curved surfaces"]
    fn volume_conservation_spheres() {
        // V(A∪B) ≈ V(A) + V(B) - V(A∩B), error < 5%
        let a = make_sphere_brep(DVec3::new(-0.5, 0.0, 0.0), 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(0.5, 0.0, 0.0), 1.0).unwrap();

        let union_result = boolean_op(BooleanOpType::Union, &a, &b);
        let inter_result = boolean_op(BooleanOpType::Intersection, &a, &b);

        assert!(
            union_result.is_ok(),
            "union failed: {:?}",
            union_result.err()
        );
        assert!(
            inter_result.is_ok(),
            "intersection failed: {:?}",
            inter_result.err()
        );

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_result.unwrap());
        let v_inter = rcad_kernel::properties::volume(&inter_result.unwrap());

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected;
        let error_pct = error * 100.0;
        assert!(
            error < 0.05,
            "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
        );
    }
}
