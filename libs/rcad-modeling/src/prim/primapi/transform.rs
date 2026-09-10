//! OCCT BRepBuilderAPI_Transform (TKTopAlgo) — apply a gp_Trsf to a shape.
//!
//! OCCT BRepBuilderAPI_Transform builds a copy of the shape and applies the
//! transformation (BRepBuilderAPI_Transform.cxx, via BRepTools_Modifier).
//! The rcad port applies the transform in place on the flat BRep pool
//! (topods::BRep::apply_transform transforms every vertex, edge curve, face
//! surface and stored location — equivalent to the OCCT modifier).

use rcad_kernel::math::gp::Trsf;
use rcad_kernel::topo::topods::BRep;

/// OCCT BRepBuilderAPI_Transform(shape, trsf) — apply `trsf` to the shape.
pub fn transform_brep(brep: &mut BRep, trsf: &Trsf) {
    brep.apply_transform(trsf.to_daffine3());
}
