pub mod topology;
pub mod topo_shape;
pub mod topods;
// OCCT TKMath/Poly package — triangulation data layer (Poly core subset).
// Lives in rcad-kernel because BRep_TFace/BRep_PolygonOnTriangulation carry
// Poly values (OCCT: TKBRep depends on TKMath, never the reverse).
pub mod poly;
// OCCT BRep_Tool::Tolerance canonical body (TKBRep/BRep/BRep_Tool.cxx).
pub mod brep_tool;
pub mod topo_query;
pub mod topo_simplify;
pub mod brep_graph;
pub mod brep_lib;
// OCCT TopoDS_Builder / TopoDS_Iterator / BRep_Builder facility (TKBRep):
// the shape -> BRep pool adoption glue.
pub mod topo_builder;
