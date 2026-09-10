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
//   Stage 3e — the interference data structures consumed by the Hider:
//   - [`bi_point`]      HLRBRep_BiPoint (the 3D segment point record).
//   - [`bi_pnt_2d`]     HLRBRep_BiPnt2D (the projected segment point
//     record).
//   - [`area_limit`]    HLRBRep_AreaLimit (a vertex on the edge with the
//     state on the left and the right).
//   - [`edge_interference_tool`] HLRBRep_EdgeInterferenceTool (the
//     HLRBRep_Data view used to instantiate the interference lists).
//   - [`edge_ilist`]    HLRBRep_EdgeIList (the sorted interference list
//     tools).
//   - [`vertex_list`]   HLRBRep_VertexList (the merged boundary/interference
//     vertex iterator).
//   - [`face_iterator`] HLRBRep_FaceIterator (the face wires/edges explorer).
//   - [`edge_builder`]  HLRBRep_EdgeBuilder (the area splitter over the
//     vertex list).

pub mod area_limit;
pub mod b_curve_tool;
pub mod bi_pnt_2d;
pub mod bi_point;
pub mod c_inter;
pub mod cl_props;
pub mod b_surface_tool;
pub mod curve;
pub mod curve_tool;
pub mod edge_builder;
pub mod data;           // HLRBRep_Data (the HLR data structure, Stage 3f)
pub mod edge_face_tool; // HLRBRep_EdgeFaceTool
pub mod face_data;      // HLRBRep_FaceData
pub mod algo;           // HLRBRep_Algo (the user entry)
pub mod hider;          // HLRBRep_Hider
pub mod hlr_to_shape;   // HLRBRep_HLRToShape
pub mod internal_algo;  // HLRBRep_InternalAlgo
pub mod shape_bounds;   // HLRBRep_ShapeBounds
pub mod shape_to_hlr;   // HLRBRep_ShapeToHLR
pub mod edge_data;
pub mod edge_ilist;
pub mod edge_interference_tool;
pub mod face_iterator;
pub mod intersector;
pub mod line_tool;
pub mod sl_props;
pub mod inter_csurf;
pub mod surface_tool;
pub mod surface;
pub mod vertex_list;

pub use b_curve_tool::{b_curve_value, CurveView};
pub use b_surface_tool::{
    nb_samples_u_range, nb_samples_u_total, nb_samples_v_range, nb_samples_v_total,
};
pub use cl_props::CLProps;
pub use c_inter::CInter;
pub use sl_props::SLProps;
pub use line_tool as line_tool_mod;
pub use surface::Surface;
