//! OCCT TKFeat — modules named 1:1 after the OCCT TKFeat packages:
//! `brep_feat_*` (BRepFeat), `loc_ope_*` (LocOpe).
//!
//! Stage 3 per the port plan (docs/tkfeat-fillet-offset-port-plan.md);
//! started in parallel with Stage 1 per ruling D7 (2026-09-07).

pub mod brep_feat_builder;
pub mod brep_feat_form;
pub mod brep_feat_form_2;
pub mod brep_feat_gluer;
pub mod brep_feat_split_shape;
pub mod brep_feat_status;
pub mod brep_feat_make_cylindrical_hole;
pub mod loc_ope_operation;
pub mod loc_ope_pnt_face;
pub mod loc_ope_find_edges;
pub mod loc_ope_find_edges_in_face;
pub mod loc_ope_cs_intersector;
pub mod loc_ope_curve_shape_intersector;
pub mod loc_ope_generated_shape;
pub mod loc_ope_glued_shape;
pub mod loc_ope_build_shape;
pub mod loc_ope_gluer;
pub mod loc_ope_build_wires;
pub mod loc_ope_spliter;
