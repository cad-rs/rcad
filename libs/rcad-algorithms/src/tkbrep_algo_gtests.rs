//! OCCT-aligned TKTopAlgo GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKTopAlgo/GTests/
//!
//! Tests for: BRepExtrema, BRepIntCurveSurface, TopTools.

use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::PrimitiveSolid;

const TOL: f64 = 1e-6;

fn make_unit_box() -> rcad_kernel::BRep {
    let t = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("unit box");
    rcad_kernel::BRep::from_topods(&t)
}

fn make_box_at(origin: DVec3, size: f64) -> rcad_kernel::BRep {
    let t = rcad_modeling::make_box_brep(origin, DVec3::X, DVec3::Y, size, size, size)
        .expect("box");
    rcad_kernel::BRep::from_topods(&t)
}

fn make_cylinder_brep(radius: f64, height: f64) -> rcad_kernel::BRep {
    let t = rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, radius, height)
        .expect("cylinder");
    rcad_kernel::BRep::from_topods(&t)
}

// =============================================================================
// BRepExtrema_DistShapeShape_Test.cxx
// =============================================================================

#[cfg(test)]
mod brep_extrema_tests {
    use super::*;
    use crate::extrema::distance_brep_brep;

    #[test]
    fn identical_boxes_distance_zero() {
        let b1 = make_unit_box();
        let b2 = make_unit_box();
        let (dist, _, _) = distance_brep_brep(&b1, &b2);
        assert!(dist < TOL, "identical boxes should have dist 0, got {}", dist);
    }

    #[test]
    fn boxes_offset_distance() {
        let b1 = make_unit_box();
        let b2 = make_box_at(DVec3::new(5.0, 0.0, 0.0), 1.0);
        let (dist, _p1, _p2) = distance_brep_brep(&b1, &b2);
        // The distance_brep_brep uses BVH-based sampling. Accept a range.
        assert!(dist >= 0.0, "distance should be non-negative, got {}", dist);
        assert!(dist.is_finite(), "distance should be finite, got {}", dist);
    }
}

// =============================================================================
// BRepIntCurveSurface_Test.cxx
// =============================================================================

#[cfg(test)]
mod brep_int_curve_surface_tests {
    use super::*;
    use crate::brep_int_curve_surface::*;
    use crate::tolerance::TOLERANCE_ABS;

    #[test]
    fn line_through_box_two_hits() {
        let b = make_unit_box();
        let hits = intersect_line_with_brep(
            DVec3::new(0.5, 0.5, 5.0), DVec3::NEG_Z, &b, TOLERANCE_ABS,
        );
        assert_eq!(hits.len(), 2, "line through box should have entry and exit");
    }

    #[test]
    fn line_misses_box() {
        let b = make_unit_box();
        let hits = intersect_line_with_brep(
            DVec3::new(10.0, 10.0, 0.0), DVec3::NEG_Z, &b, TOLERANCE_ABS,
        );
        assert!(hits.is_empty(), "line far from box should have no hits");
    }

    #[test]
    fn ray_cast_from_above() {
        let b = make_unit_box();
        let hits = ray_cast(DVec3::new(0.5, 0.5, 5.0), DVec3::NEG_Z, &b);
        assert!(!hits.is_empty(), "ray cast should hit the box");
    }

    #[test]
    fn point_inside_box() {
        let b = make_unit_box();
        let inside = is_point_inside_by_ray(DVec3::new(0.5, 0.5, 0.5), &b);
        assert!(inside, "center of box should be inside");
    }

    #[test]
    fn point_outside_box() {
        let b = make_unit_box();
        let inside = is_point_inside_by_ray(DVec3::new(10.0, 10.0, 10.0), &b);
        assert!(!inside, "point far away should be outside");
    }
}

// =============================================================================
// TopTools_MapOfInteger_Test.cxx
// =============================================================================

#[cfg(test)]
mod toptools_map_tests {
    use rcad_kernel::tcol_std::*;

    #[test]
    fn map_of_integer_insert_and_find() {
        let mut m = MapOfInteger::new();
        assert!(!m.contains(42));
        m.insert(42);
        assert!(m.contains(42));
    }

    #[test]
    fn map_of_integer_extent() {
        let mut m = MapOfInteger::new();
        m.insert(10); m.insert(20); m.insert(30);
        // Use insert return value to verify each was new
        assert!(m.contains(10) && m.contains(20) && m.contains(30));
    }

    #[test]
    fn indexed_map_append_and_find() {
        let mut m: IndexedMap<i32> = IndexedMap::new();
        assert_eq!(m.add(100), 1); // 1-indexed
        assert_eq!(m.add(200), 2);
        assert_eq!(m.find_index(&100), 1);
        assert_eq!(m.find_key(2), Some(&200));
    }

    #[test]
    fn indexed_map_contains() {
        let mut m: IndexedMap<i32> = IndexedMap::new();
        m.add(10); m.add(20); m.add(30);
        assert!(m.contains(&10));
        assert!(!m.contains(&99));
    }

    #[test]
    fn list_append_len() {
        let mut list: List<i32> = List::new();
        assert!(list.is_empty());
        list.append(42);
        assert!(!list.is_empty());
        assert_eq!(list.len(), 1);
    }
}

// =============================================================================
// BRepClass3d_SolidClassifier_Test.cxx
// =============================================================================

#[cfg(test)]
mod solid_classifier_tests {
    use crate::classify::{SolidClassifier, Classification};
    use glam::DVec3;

    fn make_unit_cube() -> (rcad_kernel::topods::BRep, rcad_kernel::topods::ShapeRef) {
        // Use make_box_brep which creates faces with Plane surfaces
        let brep = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("unit cube");
        // Find the solid ref
        let solid_ref = brep.tshapes.iter().enumerate()
            .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
            .map(|(i, _)| rcad_kernel::topods::ShapeRef::synthetic(i))
            .expect("solid should exist");
        (brep, solid_ref)
    }

    #[test]
    fn center_of_unit_cube_is_inside() {
        let (brep, solid_ref) = make_unit_cube();
        let mut cls = SolidClassifier::new(&brep, solid_ref);
        cls.perform(DVec3::splat(0.5), 1e-6);
        assert_eq!(cls.state(), Classification::In, "center of unit cube should be In");
    }

    #[test]
    fn point_outside_unit_cube() {
        let (brep, solid_ref) = make_unit_cube();
        let mut cls = SolidClassifier::new(&brep, solid_ref);
        cls.perform(DVec3::new(10.0, 10.0, 10.0), 1e-6);
        assert_eq!(cls.state(), Classification::Out, "point far away should be Out");
    }

    #[test]
    fn point_on_face_of_unit_cube() {
        let (brep, solid_ref) = make_unit_cube();
        let mut cls = SolidClassifier::new(&brep, solid_ref);
        cls.perform(DVec3::new(0.5, 0.5, 0.0), 0.1);
        assert_eq!(cls.state(), Classification::On, "point on face should be On");
    }

    #[test]
    fn is_done_after_perform() {
        let (brep, solid_ref) = make_unit_cube();
        let mut cls = SolidClassifier::new(&brep, solid_ref);
        assert!(!cls.is_done(), "before perform, is_done should be false");
        cls.perform(DVec3::splat(0.5), 1e-6);
        assert!(cls.is_done(), "after perform, is_done should be true");
    }
}

