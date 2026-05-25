use rcad_modeling::make_cylinder_brep;
use rcad_algorithms::{boolean_op_with_retry, BooleanOpType, total_surface_area};
use glam::{DAffine3, DVec3};

fn main() {
    // Match the generated test: both at center (0,0,2)
    let b1 = make_cylinder_brep(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, DVec3::X, 1.0, 4.0).expect("b1");

    let pivot = DVec3::new(0.0, 0.0, 2.0);
    let mut b2 = make_cylinder_brep(DVec3::new(0.0, 0.0, 2.0), DVec3::Z, DVec3::X, 1.0, 4.0).expect("b2");

    let rot1 = DAffine3::from_axis_angle(DVec3::X, 90.0_f64.to_radians());
    let xf1 = DAffine3::from_translation(pivot) * rot1 * DAffine3::from_translation(-pivot);
    b2.apply_transform(xf1);

    let rot2 = DAffine3::from_axis_angle(DVec3::Y, 180.0_f64.to_radians());
    let xf2 = DAffine3::from_translation(pivot) * rot2 * DAffine3::from_translation(-pivot);
    b2.apply_transform(xf2);

    eprintln!("Starting bopfuse...");
    let result = boolean_op_with_retry(BooleanOpType::Union, &b1, &b2).expect("fuse");
    let sa = total_surface_area(&result);
    eprintln!("Result SA: {} (expected 46.8319)", sa);
    eprintln!("Difference: {} ({:.3}%)", (sa - 46.8319).abs(), (sa - 46.8319).abs() / 46.8319 * 100.0);
    eprintln!("Tolerance: {}", (0.15_f64).max(0.15 * 46.8319));
}
