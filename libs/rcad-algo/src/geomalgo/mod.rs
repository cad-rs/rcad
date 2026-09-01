// OCCT TKGeomAlgo — Geometric Algorithms.
//
// OCCT: TKGeomAlgo toolkit — lower-level geometric algorithms: IntAna
// (analytic surface intersection), IntCurveSurface, IntPatch (patch-patch
// intersection), IntRes2d (2D intersection data types), IntSurf
// (surface-intersection support types), TopTrans (CurveTransition).
//
// These are consumed by the boolean pipeline (TKBO -> BOPAlgo_PaveFiller ->
// IntTools_FaceFace -> this module) and by topological classifiers
// (TKTopAlgo -> BRepClass -> IntRes2d).

pub mod int_patch;   // IntPatch + IntAna + IntCurveSurface intersection chain
pub mod int_res2d;   // IntRes2d (2D intersection data types)
pub mod int_surf;    // IntSurf (Quadric, LineOn2S, PntOn2S)
pub mod top_trans;   // TopTrans (CurveTransition, SurfaceTransition)
pub mod int_polyh;   // IntPolyh (triangle-triangle intersection support types)
pub mod intf;        // Intf (Intf_PIType, Intf_SectionPoint, Intf_Tool, InterUtils helpers)
pub mod int_curv_surf; // IntCurveSurface polygon/polyhedron sampling + IntPatch_Polyhedron
pub mod hatch;       // Geom2dHatch (Elements container + Intersector local geometry)
pub mod geom2d_gcc;   // Geom2dGcc (tangent-line constraint: Lin2d2Tan + iter solver)
pub mod approx_int;  // ApproxInt_KnotTools + ApproxInt_Approx (WLApprox) chain

pub use int_polyh::IntPolyhPoint;
pub use top_trans::surface_transition::SurfaceTransition;
pub use intf::{
    section_point_to_parameters, IntfPIType, IntfSectionPoint, IntfTool, PolyhedronLike,
    PolygonLike,
};
pub use int_curv_surf::{IntPatchPolyhedron, ThePolygonOfHInter, ThePolyhedronOfHInter};
