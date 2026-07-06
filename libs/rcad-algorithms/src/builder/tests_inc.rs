#[cfg(test)]
mod tests {
 use crate::builder::*;
 use crate::bopds::ds::DS;
 use rcad_modeling::{make_box_brep};

 #[test]
 fn prepare_returns_empty_containers() {
 let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
 let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
 let ds = DS::new(&a, &b);
 let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
 let (t_brep, result) = builder.prepare();

 assert!(t_brep.tshapes.is_empty(), "t_brep should be empty after prepare");
 assert!(result.vertices.is_empty(), "result vertices should be empty");
 assert!(result.edges.is_empty(), "result edges should be empty");
 assert!(result.faces.is_empty(), "result faces should be empty");
 assert!(result.tmp_shells.is_empty());
 assert!(result.tmp_solids.is_empty());
 }

 #[test]
 fn minimal_box_union_pipeline_builds_result() {
 // Two tiny non-overlapping boxes  ?union should produce both boxes.
 let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
 let b = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

 let mut ds = DS::new(&a, &b);
 let mut t_brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let mut filler = crate::pave_filler::PaveFiller::new(&mut ds);
 filler.brep = Some(&mut t_brep);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };

 let builder = BooleanBuilder::with_brep(&ds, BooleanOpType::Union, t_brep, face_refs, ic_edge_map);
 let (brep, _history) = builder.build_with_history().expect("union should succeed");

 // Two disjoint boxes  ?12 faces total
 assert!(!brep.solids.is_empty(), "should produce at least one solid");
 let nf: usize = brep.solids.iter()
 .flat_map(|s| &s.shells)
 .map(|sh| sh.faces.len())
 .sum();
 assert!(nf >= 12, "expected >= 12 faces for two boxes, got {}", nf);
 }
}
