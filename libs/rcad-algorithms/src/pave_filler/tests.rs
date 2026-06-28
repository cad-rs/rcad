use super::*;
use rcad_kernel::{BRep, PrimitiveSolid};
use crate::bopds::ds::DS;

#[test]
fn sparse_test() {
    // Placeholder test to verify module compiles
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let mut ds = DS::new(&a, &b);
    let _filler = PaveFiller::new(&mut ds);
}

#[test]
fn test_section_attribute_defaults() {
    let attr = SectionAttribute::default();
    assert!(attr.approximation);
    assert!(attr.pcurve_on_s1);
    assert!(attr.pcurve_on_s2);
}

#[test]
fn test_section_attribute_custom() {
    let attr = SectionAttribute { approximation: false, pcurve_on_s1: true, pcurve_on_s2: false };
    assert!(!attr.approximation);
    assert!(attr.pcurve_on_s1);
    assert!(!attr.pcurve_on_s2);
}

#[test]
fn test_edge_range_distance() {
    let d = EdgeRangeDistance { first: 0.5, last: 1.5, distance: 0.01 };
    assert!((d.first - 0.5).abs() < 1e-15);
    assert!((d.last - 1.5).abs() < 1e-15);
    assert!((d.distance - 0.01).abs() < 1e-15);
}

#[test]
fn test_pavefiller_set_get_arguments() {
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let mut ds = DS::new(&a, &b);
    let mut pf = PaveFiller::new(&mut ds);
    assert!(pf.arguments().is_empty());
    pf.set_arguments(vec![a.clone(), b.clone()]);
    assert_eq!(pf.arguments().len(), 2);
}

#[test]
fn test_pavefiller_set_get_section_attribute() {
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let mut ds = DS::new(&a, &b);
    let mut pf = PaveFiller::new(&mut ds);
    assert!(pf.section_attribute.approximation);
    pf.set_section_attribute(SectionAttribute { approximation: false, pcurve_on_s1: false, pcurve_on_s2: false });
    assert!(!pf.section_attribute.approximation);
}

#[test]
fn test_pavefiller_set_get_avoid_build_pcurve() {
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let mut ds = DS::new(&a, &b);
    let mut pf = PaveFiller::new(&mut ds);
    assert!(!pf.is_avoid_build_pcurve());
    pf.set_avoid_build_pcurve(true);
    assert!(pf.is_avoid_build_pcurve());
    pf.set_avoid_build_pcurve(false);
    assert!(!pf.is_avoid_build_pcurve());
}

#[test]
fn test_pavefiller_new_initializes_all_fields() {
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let mut ds = DS::new(&a, &b);
    let pf = PaveFiller::new(&mut ds);
    assert!(pf.my_arguments.is_empty());
    assert!(pf.is_primary);
    assert!(!pf.avoid_build_pcurve);
    assert!(pf.fpbdone.is_empty());
    assert!(pf.verts_to_avoid_extension.is_empty());
    assert!(pf.distances.is_empty());
}

#[test]
fn test_shape_info_init() {
    let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let b = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
    let ds = crate::bopds::ds::DS::new(&a, &b);
    // shape_info should cover all shapes
    assert!(ds.nb_source_shapes > 0);
    assert_eq!(ds.shape_info.len(), ds.nb_source_shapes);
    // First entries are vertices
    assert_eq!(ds.shape_info_at(0).shape_type, rcad_kernel::topods::ShapeType::Vertex);
}

