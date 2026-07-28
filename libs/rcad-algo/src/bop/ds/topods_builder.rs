//! Build DS from topods::BRep shapes (OCCT: BOPDS_DS::Init equivalent).

use crate::bop::ds::DS;
use rcad_kernel::topods::{self, TShape};
use rcad_kernel::topo_shape::Shape;
use std::sync::Arc;

/// Build a new DS from two input BReps.
pub fn new_from_topods(a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) -> DS {
    let mut ds = DS::new();
    let shape_a = Shape::new(
        Arc::new(TShape::Compound(
            a.tshapes.iter().map(|ts| {
                Shape::new(ts.clone(), 0, topods::Orientation::Forward)
            }).collect()
        )),
        0, topods::Orientation::Forward,
    );
    let shape_b = Shape::new(
        Arc::new(TShape::Compound(
            b.tshapes.iter().map(|ts| {
                Shape::new(ts.clone(), 0, topods::Orientation::Forward)
            }).collect()
        )),
        1, topods::Orientation::Forward,
    );
    ds.set_arguments(vec![shape_a, shape_b]);
    ds.init(fuzzy_tol);
    ds
}
