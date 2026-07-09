//! OCCT-aligned TKGeomAlgo GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKGeomAlgo/GTests/ (32 files, ~130 tests)
//!
//! This module covers geometric algorithm tests. Most are not yet translatable
//! because rcad lacks the required OCCT math libraries (Gordon surfaces,
//! GeomFill, IntPolyh, Plate, Hatching, etc.).
//!
//! Translatable subset (~12 tests) deferred until rcad has equivalent APIs:
//!   Geom2dAPI_InterCurveCurve_Test.cxx    — 2D curve intersection
//!   Geom2dAPI_Interpolate_Test.cxx        — 2D interpolation
//!   Geom2dAPI_PointsToBSpline_Test.cxx    — BSpline fitting 2D
//!   GeomAPI_PointsToBSpline_Test.cxx      — BSpline fitting 3D
//!   GeomAPI_PointsToBSplineSurface_Test.cxx — BSpline surface fitting
//!   GeomAPI_ProjectPointOnSurf_Test.cxx   — Point on surface projection
//!   GeomAPI_IntSS_Test.cxx                — Surface-surface intersection
//!   IntCurveSurface_IntersectionPoint_Test.cxx — Curve-surface intersection
//!   Geom2dConvert_BSplineCurveToBezierCurve_Test.cxx — BSpline->Bezier conversion
//!   Geom2dGcc_Circ2d2TanRad_Test.cxx      — Circle tangent construction
//!   Geom2dGcc_Circ2d3Tan_Test.cxx         — Circle tangent construction
//!   Geom2dGcc_Lin2d2Tan_Test.cxx          — Line tangent construction

// TODO: implement as stubs when needed
// For now, all TKGeomAlgo tests are skipped with a compile-time marker
#[cfg(test)]
#[allow(dead_code)]
mod tkgeom_algo_stub {
    #[test]
    fn all_tkgeom_algo_tests_deferred() {
        // 32 test files, ~130 tests — deferred until rcad has OCCT math equivalents
        assert!(true, "TKGeomAlgo tests deferred");
    }
}
