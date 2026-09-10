// OCCT BRepFeat_Status.hxx L20-25 + BRepFeat_StatusError.hxx L20-51 +
// BRepFeat_PerfSelection.hxx L20-37 + BRepFeat.cxx L719-809 (Print) — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Status.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_StatusError.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_PerfSelection.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat.cxx (Print)
//
// The three OCCT enums are plain C++ enums (no class, no .cxx for Status /
// StatusError; PerfSelection has no .cxx either). rcad maps each C++
// enumerator list to a Rust enum with CamelCase variants (same mapping as
// TopAbs_EDGE -> ShapeType::Edge).
//
// The BRepFeat::Print static (BRepFeat.cxx L719-809) prints a
// BRepFeat_StatusError description on a stream; it is re-hosted below because
// the status enums live here (same re-hosting model as BopToolsSet in
// brep_feat_builder.rs). The Standard_OStream out-parameter maps to a
// returned &'static str (rcad has no stream machinery on this facade).
//
// first consumer: BRepFeat_MakeCylindricalHole (3a; its local copies were
// retired and it imports BRepFeatStatus from here), BRepFeat_Form family (3b).

/// OCCT BRepFeat_Status (BRepFeat_Status.hxx L20-25).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFeatStatus {
    /// BRepFeat_NoError
    NoError,
    /// BRepFeat_InvalidPlacement
    InvalidPlacement,
    /// BRepFeat_HoleTooLong
    HoleTooLong,
}

/// OCCT BRepFeat_StatusError (BRepFeat_StatusError.hxx L21-51) — describes
/// the error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum BRepFeatStatusError {
    OK,
    BadDirect,
    BadIntersect,
    EmptyBaryCurve,
    EmptyCutResult,
    FalseSide,
    IncDirection,
    IncSlidFace,
    IncParameter,
    IncTypes,
    IntervalOverlap,
    InvFirstShape,
    InvOption,
    InvShape,
    LocOpeNotDone,
    LocOpeInvNotDone,
    NoExtFace,
    NoFaceProf,
    NoGluer,
    NoIntersectF,
    NoIntersectU,
    NoParts,
    NoProjPt,
    NotInitialized,
    NotYetImplemented,
    NullRealTool,
    NullToolF,
    NullToolU,
}

/// OCCT BRepFeat_PerfSelection (BRepFeat_PerfSelection.hxx L30-37) — to
/// declare the type of selection semantics for local operation Perform
/// methods.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum BRepFeatPerfSelection {
    /// BRepFeat_NoSelection
    NoSelection,
    /// BRepFeat_SelectionFU — selection of a face up to which a local
    /// operation will be performed
    SelectionFU,
    /// BRepFeat_SelectionU — selection of a point up to which a local
    /// operation will be performed
    SelectionU,
    /// BRepFeat_SelectionSh — selection of a shape on which a local
    /// operation will be performed
    SelectionSh,
    /// BRepFeat_SelectionShU — selection of a shape up to which a local
    /// operation will be performed
    SelectionShU,
}

/// OCCT BRepFeat::Print(se, s) (BRepFeat.cxx L719-809) — prints the error
/// description of the status <se>. The OCCT stream appends are reduced to
/// the returned string slice (the caller appends it to its sink).
pub fn print(se: BRepFeatStatusError) -> &'static str {
    match se {
        BRepFeatStatusError::OK => "No error",
        BRepFeatStatusError::BadDirect => "Directions must be opposite",
        BRepFeatStatusError::BadIntersect => "Intersection failure",
        BRepFeatStatusError::EmptyBaryCurve => "Empty BaryCurve",
        BRepFeatStatusError::EmptyCutResult => "Failure in Cut : Empty resulting shape",
        BRepFeatStatusError::FalseSide => "Verify plane and wire orientation",
        BRepFeatStatusError::IncDirection => "Incoherent Direction for shapes From and Until",
        BRepFeatStatusError::IncSlidFace => "Sliding face not in Base shape",
        BRepFeatStatusError::IncParameter => "Incoherent Parameter : shape Until before shape From",
        BRepFeatStatusError::IncTypes => {
            "Invalid option for faces From and Until : 1 Support and 1 not"
        }
        BRepFeatStatusError::IntervalOverlap => "Shapes From and Until overlap",
        BRepFeatStatusError::InvFirstShape => "Invalid First shape : more than 1 face",
        BRepFeatStatusError::InvOption => "Invalid option",
        BRepFeatStatusError::InvShape => "Invalid shape",
        BRepFeatStatusError::LocOpeNotDone => "Local Operation not done",
        BRepFeatStatusError::LocOpeInvNotDone => "Local Operation : intersection line conflict",
        BRepFeatStatusError::NoExtFace => "No Extreme faces",
        BRepFeatStatusError::NoFaceProf => "No Face Profile",
        BRepFeatStatusError::NoGluer => "Gluer Failure",
        BRepFeatStatusError::NoIntersectF => "No intersection between Feature and shape From",
        BRepFeatStatusError::NoIntersectU => "No intersection between Feature and shape Until",
        BRepFeatStatusError::NoParts => "No parts of tool kept",
        BRepFeatStatusError::NoProjPt => "No projection points",
        BRepFeatStatusError::NotInitialized => "Fields not initialized",
        BRepFeatStatusError::NotYetImplemented => "Not yet implemented",
        BRepFeatStatusError::NullRealTool => "Real Tool : Null DPrism",
        BRepFeatStatusError::NullToolF => "Null Tool : Invalid type for shape Form",
        BRepFeatStatusError::NullToolU => "Null Tool : Invalid type for shape Until",
    }
}
