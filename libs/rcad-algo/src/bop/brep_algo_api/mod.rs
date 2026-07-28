//! OCCT BRepAlgoAPI — high-level boolean operations.
//!
//! | Rust        | OCCT                       |
//! |-------------|----------------------------|
//! | fuse        | BRepAlgoAPI_Fuse           |
//! | common      | BRepAlgoAPI_Common         |
//! | cut         | BRepAlgoAPI_Cut            |

pub mod builder_operation;
pub mod section;
pub mod argument_analyzer;

use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::ds::DS;

/// Fuse (Union): combine two shapes.
pub fn fuse(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut ds = DS::new();
    ds.set_arguments(vec![
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                a.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 0, rcad_kernel::topods::Orientation::Forward,
        ),
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                b.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 1, rcad_kernel::topods::Orientation::Forward,
        ),
    ]);
    ds.init(1e-7);

    let mut builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
    builder.build().map_err(|_| BooleanError::InvalidOperation)
}

/// Common (Intersection): shared volume.
pub fn common(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut ds = DS::new();
    ds.set_arguments(vec![
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                a.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 0, rcad_kernel::topods::Orientation::Forward,
        ),
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                b.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 1, rcad_kernel::topods::Orientation::Forward,
        ),
    ]);
    ds.init(1e-7);
    let mut builder = BooleanBuilder::new(&ds, BooleanOpType::Intersection);
    builder.build().map_err(|_| BooleanError::InvalidOperation)
}

/// Cut (Difference): a minus b.
pub fn cut(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut ds = DS::new();
    ds.set_arguments(vec![
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                a.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 0, rcad_kernel::topods::Orientation::Forward,
        ),
        rcad_kernel::topo_shape::Shape::new(
            std::sync::Arc::new(rcad_kernel::topods::TShape::Compound(
                b.tshapes.iter().map(|ts| rcad_kernel::topo_shape::Shape::new(
                    ts.clone(), 0, rcad_kernel::topods::Orientation::Forward
                )).collect()
            )), 1, rcad_kernel::topods::Orientation::Forward,
        ),
    ]);
    ds.init(1e-7);
    let mut builder = BooleanBuilder::new(&ds, BooleanOpType::Difference);
    builder.build().map_err(|_| BooleanError::InvalidOperation)
}
