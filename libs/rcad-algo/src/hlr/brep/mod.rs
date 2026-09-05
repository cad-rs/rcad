// OCCT HLRBRep (TKHLR) — the exact-HLR internal package. Stage 3a: the
// adaptor/tool layer consumed by the HLR Data (3f) and the contour engine
// (Contap).  Stage 3b: the intersection layer (CurveTool / CInter /
// Intersector) over the projected curves.
//
//   - [`surface`]       HLRBRep_Surface (the face adaptor wrapper with the
//     projector-based IsSide/IsAbove queries).
//   - [`curve`]         HLRBRep_Curve (the 3D edge curve seen as a 2D
//     projection with optional perspective).
//   - [`b_surface_tool`] HLRBRep_BSurfaceTool (statics over the face
//     adaptor).
//   - [`b_curve_tool`]  HLRBRep_BCurveTool (statics over the edge adaptor).
//   - [`curve_tool`]    HLRBRep_CurveTool (the static tool over
//     HLRBRep_Curve for the CInter instantiation).
//   - [`c_inter`]       HLRBRep_CInter + the The*OfCInter instantiation
//     family over the projected curves.
//   - [`edge_data`]     HLRBRep_EdgeData (the per-edge record).
//   - [`intersector`]   HLRBRep_Intersector (the 2D/CS intersection layer).

pub mod b_curve_tool;
pub mod c_inter;
pub mod cl_props;
pub mod b_surface_tool;
pub mod curve;
pub mod curve_tool;
pub mod edge_data;
pub mod intersector;
pub mod line_tool;
pub mod sl_props;
pub mod inter_csurf;
pub mod surface_tool;
pub mod surface;

pub use b_curve_tool::{b_curve_value, CurveView};
pub use b_surface_tool::{
    nb_samples_u_range, nb_samples_u_total, nb_samples_v_range, nb_samples_v_total,
};
pub use cl_props::CLProps;
pub use c_inter::CInter;
pub use sl_props::SLProps;
pub use line_tool as line_tool_mod;
pub use surface::Surface;
