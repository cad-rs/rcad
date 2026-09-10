// OCCT LocOpe_Operation.hxx L20-25 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Operation.hxx
//
// LocOpe_Operation is a plain C++ enum (no class, no .cxx). rcad maps the
// C++ enumerator list to a Rust enum with CamelCase variants (same mapping
// as TopAbs_EDGE -> ShapeType::Edge).
//
// first consumer: BRepFeat_Form family (3b) through LocOpe_Gluer::OpeType;
// also consumed by LocOpe_Gluer::myOpe.

/// OCCT LocOpe_Operation (LocOpe_Operation.hxx L20-25).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocOpeOperation {
    /// LocOpe_FUSE
    Fuse,
    /// LocOpe_CUT
    Cut,
    /// LocOpe_INVALID
    Invalid,
}
