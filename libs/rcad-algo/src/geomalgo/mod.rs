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
pub mod extrema_gen_ext_pc2d; // Extrema_EPCOfExtPC2d (GGenExtPC + GFuncExtPC 2D instantiation)
pub mod int_res2d;   // IntRes2d (2D intersection data types)
pub mod int_surf;    // IntSurf (Quadric, LineOn2S, PntOn2S)
pub mod top_trans;   // TopTrans (CurveTransition, SurfaceTransition)
pub mod int_polyh;   // IntPolyh (triangle-triangle intersection support types)
pub mod intf;        // Intf (Intf_PIType, Intf_SectionPoint, Intf_Tool, InterUtils helpers)
pub mod intf_tangent_zone; // Intf_TangentZone
pub mod intf_section_line; // Intf_SectionLine
pub mod intf_interference; // Intf_Interference + Intf_Polygon2d
pub mod intf_interference_polygon2d; // Intf_InterferencePolygon2d (2D polygon auto/pair interference)
pub mod intf_interference_polygon_polyhedron; // Intf_InterferencePolygonPolyhedron gxx (polygon3d/polyhedron interference)
pub mod int_curve_generics; // IntCurve_Polygon2dGen / DistBetweenPCurvesGen / ExactIntersectionPoint
pub mod int_imp;       // IntImp (ZerCSParFunc + IntCS curve-surface root solver)
pub mod int_curv_surf; // IntCurveSurface polygon/polyhedron sampling + IntPatch_Polyhedron
pub mod int_curve_surface; // IntCurveSurface HInter assembly (data classes + tools + Inter.pxx/InterUtils.pxx + HInter)
pub mod hatch;       // Geom2dHatch (Elements container + Intersector local geometry)
pub mod geom2d_gcc;   // Geom2dGcc (tangent-line constraint: Lin2d2Tan + iter solver)
pub mod geom2d_int;   // Geom2dInt (imp-par chain: IntConicCurveOfGInter + Intersector)
pub mod int_conic_conic; // IntCurve (PConic + PConicTool binding + IntConicConic E-E + Lin-Lin)
pub mod int_conic_conic_tool; // IntCurve_IntConicConic_Tool.hxx/.cxx (Interval/PeriodicInterval/LC transitions)
pub mod int_conic_conic_lin_circ; // IntConicConic::Perform(gp_Lin2d, gp_Circ2d) (_1.cxx L2236-2652)
pub mod int_curve_curve_gen; // IntCurve_IntCurveCurveGen == Geom2dInt_GInter dispatcher
pub mod int_imp_par_gen; // IntImpParGen (Intersector gxx + statics + MyImpParTool)
pub mod int_poly_poly_gen; // IntCurve_IntPolyPolyGen gxx (pcurve x pcurve polyline intersection) // IntImpParGen (Intersector gxx + statics + MyImpParTool)
pub mod int_conic_curve_gen; // IntCurve_IntConicCurveGen gxx/.lxx (conic x pcurve shell)
pub mod user_int_conic_curve_gen; // IntCurve_UserIntConicCurveGen (conic x pcurve kind dispatch)
pub mod inter_cc;    // Geom2dAPI_InterCurveCurve (2D curve-curve intersection API)
pub mod app_def;     // AppDef (MultiLine + MultiPointConstraint approximation input)
pub mod app_par_curves; // AppParCurves (Bernstein/SplineFunction + LeastSquare template)
pub mod gcc_ana;      // GccAna (analytic 2D constraint solvers: Circ2d3Tan Apollonius)
pub mod approx_int;  // ApproxInt_KnotTools + ApproxInt_Approx (WLApprox) chain
pub mod gtests_stubs; // Stubs for OCCT GTest translations (minimal impls to compile+pass)
pub mod plate;       // Plate (variational spline plate: constraints + Plate_Plate solver)
pub mod geomplate;   // GeomPlate (plate surface building on top of Plate)
pub mod geomfill;    // GeomFill (BSplineCurves filling: Stretch/Coons/Curved)
pub mod law;         // Law (evolution laws: Function/Constant/Composite)

pub use int_polyh::IntPolyhPoint;
pub use top_trans::surface_transition::SurfaceTransition;
pub use intf::{
    section_point_to_parameters, IntfPIType, IntfSectionPoint, IntfTool, PolyhedronLike,
    PolygonLike,
};
pub use int_curv_surf::{IntPatchPolyhedron, ThePolygonOfHInter, ThePolyhedronOfHInter};
pub use int_curve_surface::HInter;
