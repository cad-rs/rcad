//! TKTopAlgo GTest translations — ported from rcad-algorithms.
//! OCCT source: src/ModelingAlgorithms/TKTopAlgo/GTests/
//! Adapted for rcad_kernel/rcad_algo APIs.

use glam::DVec3;
use rcad_kernel::topods;
use rcad_algo::topalgo::brep_class3d::solid_classifier::SolidClassifier;

fn make_unit_box_brep() -> (topods::BRep, topods::Shape) {
    let brep = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
        .expect("unit box");
    let solid_ref = brep.tshapes.iter().enumerate()
        .find(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Solid(_)))
        .map(|(i, _)| topods::Shape::new(brep.tshapes[i].clone(), 0, topods::Orientation::Forward))
        .expect("solid");
    (brep, solid_ref)
}

// BRepClass3d_SolidClassifier_Test.cxx
mod solid_classifier_tests {
    use super::*;

    fn point_state(x: f64, y: f64, z: f64, tol: f64) -> u8 {
        let (brep, solid_ref) = make_unit_box_brep();
        let mut cls = SolidClassifier::new();
        cls.load(&solid_ref);
        cls.perform(DVec3::new(x, y, z), tol);
        cls.state()
    }

    const INSIDE: u8 = 1;  // Classification::In
    const OUTSIDE: u8 = 2; // Classification::Out
    const ON: u8 = 3;      // Classification::On

    #[test] fn center_inside() { assert_eq!(point_state(0.5, 0.5, 0.5, 1e-6), INSIDE); }
    #[test] fn point_outside() { assert_eq!(point_state(10.0, 10.0, 10.0, 1e-6), OUTSIDE); }
    #[test] fn point_on_face() { assert_eq!(point_state(0.5, 0.5, 0.0, 0.1), ON); }
}
