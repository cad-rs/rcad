// OCCT BRepTopAdaptor — Topology adaptors for classification.
//
// OCCT ref: TKTopAlgo/BRepTopAdaptor/
// Contains FClass2d (2D face classifier), TopolTool, HVertex, Tool.

pub mod class2d;
pub mod fclass2d;
pub mod fclass2d_topol; // BRepTopAdaptor_FClass2d (the TopolTool classifier)
pub mod hvertex_brep;   // BRepTopAdaptor_HVertex
pub mod topol_tool_brep; // BRepTopAdaptor_TopolTool
pub mod tool;            // BRepTopAdaptor_Tool (TopolTool+HVertex cache for the MST map)
pub mod brep_adaptor_bridge {}
