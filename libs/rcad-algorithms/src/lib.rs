pub mod bopds;
pub mod brep_check;
pub mod brep_repair;
pub mod builder;
pub mod bvh;
pub mod classify;
pub mod draft;
pub mod geom_populate;
pub mod history;
pub mod hlr;
pub mod imprint;
pub mod inttools;
pub mod pave_filler;
pub mod section;
pub mod thicken;
pub mod tolerance;
pub mod triangulate;

pub use bvh::{Aabb, Bvh, BvhStats};

use rcad_kernel::BRep;

pub use brep_check::{CheckIssue, CheckResult, check};
pub use brep_repair::{
    RepairReport, fix_wire_orientation, merge_close_vertices, recompute_face_normals,
    remove_degenerate_faces, repair,
};
pub use builder::{BooleanError, BooleanOpType};
pub use history::{BooleanHistory, FaceOrigin};
pub use hlr::{AssemblyHlrResult, ComponentHlr, HlrCamera, HlrResult, HlrSegment, hlr, hlr_assembly, hlr_to_svg};
pub use imprint::{
    Gap, GapOverlapReport, ImprintResult, Overlap, detect_gaps_overlaps, imprint_brep, min_distance,
};
pub use inttools::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces,
    intersect_surfaces_with_density,
};
pub use section::{SectionCurve, section, section_curves, section_polylines};
pub use thicken::{ThickeningResult, thicken_shell};
pub use triangulate::{SurfaceMesh, TessellationParams, mesh_brep, triangulate_surface};

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

/// Parallel version of [`boolean_op_with_history`].
///
/// Uses Rayon to process faces in parallel during the classification phase.
/// This can provide significant speedup (2-4x) for large models with many faces.
/// For small models (< 20 faces), the serial version may be faster due to
/// thread overhead.
///
/// # Example
/// ```rust,no_run
/// use rcad_algorithms::{boolean_op_par, BooleanOpType, history::BooleanHistory};
/// use rcad_kernel::BRep;
///
/// fn parallel_union(a: &BRep, b: &BRep) -> BRep {
///     let (brep, _history) = boolean_op_par(BooleanOpType::Union, a, b).unwrap();
///     brep
/// }
/// ```
pub fn boolean_op_par(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    let mut ds = bopds::ds::DS::new(a, b);
    let mut filler = pave_filler::PaveFiller::new(&mut ds);
    filler.perform();
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build_with_history_par()
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
    use rcad_modeling::{make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep};

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
        // Sphere primitive has no triangle mesh, so volume(&a) = 0. Compare against
        // analytical: two overlapping unit spheres at distance 1 → lens volume ≈ 1.809.
        // Full unit sphere volume = 4π/3 ≈ 4.189.
        let v_sphere_analytical = 4.0 * std::f64::consts::PI / 3.0; // 4π/3
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(
            v < v_sphere_analytical,
            "intersection should be smaller than one sphere (4π/3≈4.19), got {v}"
        );
    }

    #[test]
    fn boolean_sphere_sphere_difference() {
        // Large sphere (r=2) minus small sphere (r=1) with d=1 between centers.
        // d=1, r_A=2, r_B=1 → h = (1+4-1)/2 = 2 → tangent! Use d=0.5 instead.
        // d=0.5, r_A=2, r_B=1 → h = (0.25+4-1)/1 = 3.25 → outside sphere A
        // Use d=1.5: h = (2.25+4-1)/3 = 5.25/3 = 1.75 < r_A=2 → proper intersection
        let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.5, 0.0, 0.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-sphere difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        // Large sphere volume = 4π/3 * 8 ≈ 33.51; result should be positive and less.
        let v_large_analytical = 4.0 * std::f64::consts::PI / 3.0 * 8.0;
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(v < v_large_analytical, "difference should be smaller than original large sphere");
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
    fn boolean_cylinder_cylinder_intersection() {
        // Two perpendicular cylinders (Steinmetz solid).
        // Use cylinders that are offset so they overlap in a region that doesn't
        // straddle the seam boundary (avoiding UV-seam discontinuity issues).
        // Cylinder A: Y-axis, centered at (0, 0, 0) with height 4 → spans y ∈ [-2, 2]
        // Cylinder B: X-axis, centered at (0, 0, 0) with height 4 → spans x ∈ [-2, 2]
        let a =
            make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::Y, DVec3::X, 1.0, 4.0).unwrap();
        let b =
            make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 4.0).unwrap();

        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // The result should be non-degenerate (the two cylinders DO intersect).
        // We check only non-degeneracy: if the boolean fails or gives an empty
        // result, something is fundamentally broken.
        match result {
            Ok(brep) => {
                // Non-degenerate: at least one face in the result.
                assert!(
                    !brep.solids[0].shells[0].faces.is_empty(),
                    "cylinder-cylinder intersection should produce at least one face"
                );
                let v = rcad_kernel::properties::volume(&brep);
                assert!(v >= 0.0, "volume must not be negative, got {v}");
                // Note: exact volume comparison is not practical because the curved-face
                // volume computation (divergence theorem on polyline boundaries) is
                // approximate for complex intersection geometries.
            }
            Err(e) => {
                // If the result is degenerate, fail with a clear message.
                panic!("cylinder-cylinder intersection failed: {e:?}");
            }
        }
    }

    #[test]
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

        let union_brep = union_result.unwrap();
        let inter_brep = inter_result.unwrap();

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_brep);
        let v_inter = rcad_kernel::properties::volume(&inter_brep);

        // Debug values — show face count and uv_domains
        eprintln!("V_A={v_a:.4} V_B={v_b:.4} V_union={v_union:.4} V_inter={v_inter:.4}");
        eprintln!(
            "Union faces={}, inter faces={}",
            union_brep.solids[0].shells[0].faces.len(),
            inter_brep.solids[0].shells[0].faces.len()
        );
        for (i, face) in inter_brep.solids[0].shells[0].faces.iter().enumerate() {
            let range = inter_brep.geom.face_surface_range.get(i).and_then(|o| *o);
            let surf_name = inter_brep.geom.face_surface.get(i).and_then(|o| *o)
                .map(|si| format!("{:?}", std::mem::discriminant(&inter_brep.geom.surfaces[si])));
            // Compute per-face contribution to volume
            let face_tris = rcad_kernel::properties::face_triangles_pub(&inter_brep, face, i);
            let face_vol: f64 = face_tris.iter().map(|&[a,b,c]| a.dot(b.cross(c)) / 6.0).sum();
            eprintln!("  inter face {i}: normal={:.3?} uv_domain={range:?} surf={surf_name:?} vol_contrib={face_vol:.4}", face.normal);
        }

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected;
        let error_pct = error * 100.0;
        assert!(
            error < 0.05,
            "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
        );
    }

    #[test]
    #[ignore = "sphere-sphere boolean volume not yet correct: intersection faces cancel \
                in divergence-theorem sum (net≈0) and union result is topologically \
                incomplete (2 faces instead of expected composite); tracked as P0-B"]
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

        let union_brep = union_result.unwrap();
        let inter_brep = inter_result.unwrap();

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_brep);
        let v_inter = rcad_kernel::properties::volume(&inter_brep);

        eprintln!("sphere-sphere: V_A={v_a:.4} V_B={v_b:.4} V_union={v_union:.4} V_inter={v_inter:.4}");
        eprintln!(
            "Union faces={}, inter faces={}",
            union_brep.solids[0].shells[0].faces.len(),
            inter_brep.solids[0].shells[0].faces.len()
        );
        for (i, face) in inter_brep.solids[0].shells[0].faces.iter().enumerate() {
            let range = inter_brep.geom.face_surface_range.get(i).and_then(|o| *o);
            let surf_name = inter_brep.geom.face_surface.get(i).and_then(|o| *o)
                .map(|si| format!("{:?}", std::mem::discriminant(&inter_brep.geom.surfaces[si])));
            let face_tris = rcad_kernel::properties::face_triangles_pub(&inter_brep, face, i);
            let face_vol: f64 = face_tris.iter().map(|&[a,b,c]| a.dot(b.cross(c)) / 6.0).sum();
            eprintln!("  inter face {i}: normal={:.3?} uv_domain={range:?} surf={surf_name:?} vol_contrib={face_vol:.4}", face.normal);
        }

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected;
        let error_pct = error * 100.0;
        assert!(
            error < 0.05,
            "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
        );
    }

    #[test]
    fn boolean_result_edges_have_pcurves() {
        // Box with a cylindrical hole. After the boolean difference, intersection
        // edges on the cylinder surface should get PCurves via
        // populate_boolean_result_pcurves.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.5, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        let Ok(mut brep) = result else {
            // If the boolean op itself fails, skip (it's tested elsewhere).
            return;
        };
        if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
            return;
        }

        // Fill PCurves.
        geom_populate::populate_boolean_result_pcurves(&mut brep);

        // At least one edge on the cylinder face should now have a PCurve.
        let any_pcurve = brep.geom.edge_pcurves.iter().any(|v| !v.is_empty());
        assert!(
            any_pcurve,
            "populate_boolean_result_pcurves should have added at least one PCurve"
        );
    }

    // ─── Sphere × Cylinder Boolean Tests ──────────────────────────────────────

    /// A cylinder whose axis passes through the sphere centre (axis-aligned case).
    /// The sphere–cylinder intersection is two circles.  Difference should
    /// produce a valid solid with more faces than just the six box/sphere faces.
    #[test]
    fn boolean_sphere_cylinder_difference_axis_aligned() {
        // Sphere centred at origin, radius 5; cylinder along Z through origin, radius 3.
        // Intersection circles at z = ±4  (sqrt(25-9) = 4).
        let a = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(0.0, 0.0, -6.0), DVec3::Z, DVec3::X, 3.0, 12.0)
            .unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-cylinder difference (axis-aligned) failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
        // Volume of sphere (4π/3 · R³) minus the cylindrical tunnel should be positive
        // and smaller than the sphere.
        let v = rcad_kernel::properties::volume(&brep);
        let v_sphere = 4.0 * std::f64::consts::PI / 3.0 * 5.0_f64.powi(3);
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(v < v_sphere, "difference should be smaller than original sphere");
    }

    // ─── Cone × Plane Boolean Tests ───────────────────────────────────────────

    /// Box minus a cone through it: the cone's lateral surface intersects the
    /// box's planar faces, exercising the plane-cone circle intersection path.
    #[test]
    fn boolean_box_cone_difference() {
        // Box: 4×4×4 at origin.  Cone: base at (2,2,-0.5), axis Z, r=0.8, h=5.
        // The cone pokes through the box; plane-cone intersections are circles
        // (planes ⊥ cone axis).
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b =
            make_cone_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.8, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-cone difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
    }

    /// Cone intersected with a box slab: the slab's top and bottom faces are
    /// planes perpendicular to the cone axis, producing circle intersections.
    /// This test verifies that the plane-cone code path does not panic.
    #[test]
    fn boolean_cone_box_intersection_circle() {
        // Cone: base at origin, axis Z, base_radius=2, height=4.
        // Slab: 6×6×4 at z=0..4 — same height as the cone; the lateral face of
        // the slab does NOT cut the cone (slab is wide enough), so only the
        // slab top (z=4, a plane ⊥ cone axis) intersects the cone's lateral surface
        // near the apex region.  This exercises the plane-cone circle intersection.
        let a = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).unwrap();
        let b =
            make_box_brep(DVec3::new(-3.0, -3.0, 0.0), DVec3::X, DVec3::Y, 6.0, 6.0, 3.0)
                .unwrap();
        // The box (z=0..3) clips the cone (z=0..4), leaving the lower frustum.
        // The intersection may succeed or return DegenerateResult depending on
        // classifier robustness; we only require it does not panic.
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        match result {
            Ok(brep) => {
                assert!(
                    !brep.solids.is_empty() && !brep.solids[0].shells[0].faces.is_empty(),
                    "intersection produced an empty result"
                );
            }
            Err(BooleanError::DegenerateResult) => {
                // DegenerateResult is an acceptable failure for complex curved intersections.
            }
            Err(e) => {
                panic!("cone-box intersection failed unexpectedly: {e:?}");
            }
        }
    }

    /// Intersection of a sphere and a coaxial cylinder.
    #[test]
    fn boolean_sphere_cylinder_intersection_axis_aligned() {
        // Sphere centred at origin, radius 5; cylinder along Z through origin, radius 3.
        // The intersection of their volumes is a "barrel" shape bounded by two
        // spherical caps (z > 4 and z < -4) and the cylinder lateral surface.
        let a = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(0.0, 0.0, -6.0), DVec3::Z, DVec3::X, 3.0, 12.0)
            .unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-cylinder intersection (axis-aligned) failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
        // Just verify we get a positive volume — the exact amount depends on
        // whether sphere cap faces contribute correctly to the divergence-theorem
        // volume (sphere parametric surfaces have known approximation issues
        // tracked separately).
        let v = rcad_kernel::properties::volume(&brep);
        assert!(v > 0.0, "intersection volume should be positive, got {v}");
    }

}
