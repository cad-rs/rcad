//! OCCT-aligned TKGeomBase GTest translations.
//!
//! This module translates all GTests from
//!   $OCCT_SRC/src/ModelingData/TKGeomBase/GTests/
//! into Rust unit tests and provides OCCT-aligned implementations.
//!
//! ## Completed translations
//!
//! | OCCT File | Rust Module | Status |
//! |-----------|-------------|--------|
//! | `BndLib_Test.cxx` | `bnd_lib.rs` / `bnd_lib_2d.rs` | ✅ Pre-existing |
//! | `GeomBndLib_Curve_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_Surface_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_OffsetCurve_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_OffsetSurface_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_SurfaceOfRevolution_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_SurfaceOfExtrusion_Test.cxx` | `bnd_lib.rs` | ✅ Pre-existing |
//! | `GeomBndLib_Curve2d_Test.cxx` | `bnd_lib_2d.rs` | ✅ Pre-existing |
//! | `GeomBndLib_OffsetCurve2d_Test.cxx` | `bnd_lib_2d.rs` | ✅ Pre-existing |
//! | `AdvApp2Var_*.cxx` | `adv_app2_var.rs` | ✅ Pre-existing |
//! | `AppCont_ContMatrices_Test.cxx` | `app_cont.rs` | ✅ Pre-existing |
//! | `Approx_BSplineApproxInterp_Test.cxx` | `bspline_approx_interp.rs` | ✅ Pre-existing |
//! | `GC_MakeSegment2d_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `GC_MakePlane_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `GC_MakeParabola2d_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `GC_MakeCircle2d_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `GC_MakeConicalSurface_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `GC_MakeArcOfCircle_Test.cxx` | `gc_make.rs` | ✅ Pre-existing |
//! | `gce_Make*.cxx` (5 files) | `gc_make.rs` | ✅ Pre-existing |
//! | `GCPnts_AbscissaPoint_Test.cxx` | `gcpnts.rs` | ✅ Fixed |
//! | `LProp_CurAndInf_Test.cxx` | `lprop_cur_and_inf.rs` | ✅ Fixed |
//! | `GeomLProp_CLProps2d_Test.cxx` | `curve2d_props` | ✅ New |
//! | `GeomLProp_CurAndInf2d_Test.cxx` | `curve2d_props` | ✅ New |
//! | `GProp_PGProps_Test.cxx` | `point_set_props` | ✅ New |
//! | `GProp_PEquation_Test.cxx` | `point_set_equation` | ✅ New |
//! | `ExtremaPC_Line_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_Circle_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_Ellipse_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_Parabola_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_Hyperbola_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_BezierCurve_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_BSplineCurve_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `ExtremaPC_OffsetCurve_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `GeomConvert_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `IntAna_IntQuadQuad_Test.cxx` | `extrema_pc_tests` | ✅ New |
//! | `Hermit_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `GeomConvert_CompCurveToBSplineCurve_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `ProjLib_Cone_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `ExtremaPC_SearchMode_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `ExtremaPC_Comparison_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `ExtremaPC_ExtendedGeometry_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `Extrema_ExtPC_Test.cxx` | `remaining_gtests` | ✅ New |
//! | `ProjLib_ComputeApproxOnPolarSurface_Test.cxx` | `remaining_gtests` | ✅ New (simplified) |
//!
//! All 43 OCCT GTest files in `TKGeomBase/GTests/` are now translated.

pub mod curve2d_props;
pub mod point_set_props;
pub mod point_set_equation;
pub mod extrema_pc_tests;
pub mod remaining_gtests;
