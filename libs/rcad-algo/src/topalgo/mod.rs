// OCCT TKTopAlgo — Topological Algorithms.
//
// OCCT: TKTopAlgo toolkit — topological algorithms belonging
// to TKTopAlgo layer: BRepClass3d, BRepClass, BRepExtrema, BRepBndLib, etc.

pub mod brep_bnd_lib;       // BRepBndLib
pub mod brep_class;         // BRepClass (FaceExplorer, Edge, FClassifier)
pub mod brep_class3d;       // BRepClass3d
pub mod shape_source;       // ShapeSource — DS shape access for the classifiers
pub mod brep_extrema;       // BRepExtrema
pub mod brep_int_curve_surface; // BRepIntCurveSurface
pub mod brep_lib;           // BRepLib
pub mod brep_top_adaptor;   // BRepTopAdaptor (FClass2d, TopolTool)
pub mod brep_check;         // BRepCheck (migrated from algo_ext)
pub mod gcpnts;             // GCPnts (QuasiUniformDeflection)
pub mod brep_copy;          // BRepBuilderAPI_Copy
pub mod thru_sections;      // BRepOffsetAPI_ThruSections (loft — BRepFill port pending)
