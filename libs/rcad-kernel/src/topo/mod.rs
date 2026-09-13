pub mod topology;
pub mod topo_shape;
pub mod topods;
// OCCT BRep_Tool::Tolerance canonical body (TKBRep/BRep/BRep_Tool.cxx).
pub mod brep_tool;
pub mod topo_query;
pub mod topo_simplify;
pub mod brep_graph;
pub mod brep_lib;
// OCCT TopoDS_Builder / TopoDS_Iterator / BRep_Builder facility (TKBRep):
// the shape -> BRep pool adoption glue.
pub mod topo_builder;
