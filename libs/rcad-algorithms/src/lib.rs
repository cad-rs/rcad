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
pub use inttools::{SurfaceCurve, SurfaceSurfaceIntersection, intersect_surfaces};
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
}
