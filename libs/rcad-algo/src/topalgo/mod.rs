// OCCT TKTopAlgo — Topological Algorithms.
//
// OCCT: TKTopAlgo toolkit — topological algorithms belonging
// to TKTopAlgo layer: BRepClass3d, BRepClass, BRepExtrema, BRepBndLib, etc.

pub mod adaptor2d;          // Adaptor2d (Line2d — the TopolTool restriction adaptor)
pub mod adaptor3d;          // Adaptor3d (HVertex + TopolTool)
pub mod brep_bnd_lib;       // BRepBndLib
pub mod brep_class;         // BRepClass (FaceExplorer, Edge, FClassifier)
pub mod brep_class3d;       // BRepClass3d
pub mod shape_source;       // ShapeSource — DS shape access for the classifiers
pub mod brep_extrema;       // BRepExtrema
pub mod brep_int_curve_surface; // BRepIntCurveSurface
pub mod brep_lib;           // BRepLib
pub mod brep_adaptor;       // BRepAdaptor (Curve2d)
pub mod brep_top_adaptor;   // BRepTopAdaptor (FClass2d, TopolTool)
pub mod brep_check;         // BRepCheck (migrated from algo_ext)
pub mod mat2d;            // MAT2d (bisecting locus on 2d lines — OffsetWire engine)
pub mod bisector;          // Bisector (2d bisector curves — MAT2d feed)
pub mod mat;               // MAT (MAT graph — consumed by MAT2d)
pub mod gcpnts;             // GCPnts (QuasiUniformDeflection)
pub mod brep_copy;          // BRepBuilderAPI_Copy
pub mod brep_builderapi_transform; // BRepBuilderAPI_Transform + BRepTools_TrsfModification
pub mod brep_tools_modification;   // BRepTools_Modification + BRepTools_TrsfModification
pub mod brep_tools_modifier;       // BRepTools_Modifier
pub mod thru_sections;      // BRepOffsetAPI_ThruSections (loft — BRepFill port pending)
pub mod brep_tools_substitution;
pub mod brep_tools_quilt;
pub mod brep_lib_find_surface;
pub mod brep_lib_validate_edge;
pub mod brep_lib_encode_regularity;
pub mod int_curves_face_intersector; // IntCurvesFace (Intersector/ShapeIntersector)
