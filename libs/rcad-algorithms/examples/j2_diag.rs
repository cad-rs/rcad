use glam::{DAffine3, DVec3};
use rcad_modeling::{make_box_brep, make_cylinder_brep};
use rcad_algorithms::{boolean_op_with_retry, BooleanOpType, total_surface_area, brep_tools::extract_solids};
fn main() {
    let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 10.0, 10.0, 10.0).unwrap();
    let c = b.clone();
    let s = make_cylinder_brep(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, DVec3::X, 2.0, 4.0).unwrap();
    let mut s2 = s;
    s2.apply_transform(DAffine3::from_translation(DVec3::new(5.0, 5.0, -2.0)));
    let rr = boolean_op_with_retry(BooleanOpType::Difference, &c, &s2).unwrap();
    println!("bcut SA: {}", total_surface_area(&rr));
    let solids = extract_solids(&rr);
    println!("bcut solids: {}", solids.len());
    for (i, sol) in solids.iter().enumerate() {
        println!("solid[{}] SA: {}", i, total_surface_area(sol));
    }
    if !solids.is_empty() {
        let r = boolean_op_with_retry(BooleanOpType::Intersection, &solids[0], &c).unwrap();
        println!("bcommon SA: {}", total_surface_area(&r));
    }
}
