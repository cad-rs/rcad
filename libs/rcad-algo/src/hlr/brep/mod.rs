// OCCT HLRBRep (TKHLR) — the exact-HLR internal package. Stage 3a: the
// adaptor/tool layer consumed by the HLR Data (3f) and the contour engine
// (Contap).
//
//   - [`surface`]       HLRBRep_Surface (the face adaptor wrapper with the
//     projector-based IsSide/IsAbove queries).
//   - [`curve`]         HLRBRep_Curve (the 3D edge curve seen as a 2D
//     projection with optional perspective).
//   - [`b_surface_tool`] HLRBRep_BSurfaceTool (statics over the face
//     adaptor).
//   - [`b_curve_tool`]  HLRBRep_BCurveTool (statics over the edge adaptor).

pub mod b_curve_tool;
pub mod cl_props;
pub mod b_surface_tool;
pub mod curve;
pub mod line_tool;
pub mod surface;

pub use b_curve_tool::{b_curve_value, CurveView};
pub use b_surface_tool::{
    nb_samples_u_range, nb_samples_u_total, nb_samples_v_range, nb_samples_v_total,
};
pub use cl_props::CLProps;
pub use line_tool as line_tool_mod;
pub use surface::Surface;
