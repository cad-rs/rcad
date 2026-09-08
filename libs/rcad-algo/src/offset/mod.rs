//! OCCT TKOffset — modules named 1:1 after the OCCT TKOffset packages:
//! `brep_offset_*` (BRepOffset), `brep_offset_api_*` (BRepOffsetAPI),
//! `bi_tgte_*` (BiTgte), `draft_*` (Draft).  Stage 2 per the port plan.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/

pub mod brep_offset_make_simple_offset;
pub mod bi_tgte_blended;
pub mod bi_tgte_contact;
pub mod bi_tgte_curve_on_edge;
pub mod bi_tgte_curve_on_vertex;
pub mod brep_offset_inter2d;
pub mod brep_offset_inter2d_b;
pub mod brep_offset_inter3d;
pub mod brep_offset_offset;
pub mod brep_offset_surface;
pub mod brep_offset_tool;
pub mod brep_offset_tool_b;
pub mod brep_offset_tool_c;
pub mod brep_offset_tool_d;
pub mod draft_modification;
pub mod draft_modification_1;
pub mod draft_modification_1_b;
pub mod draft_modification_1_c;
pub mod draft;
pub mod draft_edge_info;
pub mod draft_error_status;
pub mod draft_face_info;
pub mod draft_vertex_info;
pub mod brep_offset_offset_b;
