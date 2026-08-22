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
pub mod top_trans;   // TopTrans (CurveTransition)
pub mod approx_int;  // ApproxInt_KnotTools + ApproxInt_Approx (WLApprox) chain
