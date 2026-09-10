//! OCCT ShapeFix — package module.
//!
//! W1-6 pull-forward: the `ShapeFix::SameParameter` package static
//! (ShapeFix.cxx) — the DT_SplitByNumber heal oracle dependency and the
//! chfi3d same_parameter_pass backfill. The ShapeFix class translations
//! follow in W3.
//!
//! W3 tranche 1 (2026-09-10): ShapeFix_Root, ShapeFix_ShapeTolerance,
//! ShapeFix_Edge, ShapeFix_WireSegment, ShapeFix_WireVertex,
//! ShapeFix_SplitCommonVertex, ShapeFix_FreeBounds.  The remaining W3
//! classes (Wire/Wire_1, Face, Shell, Shape, IntersectionTool, SplitTool,
//! Wireframe, FixSmallFace, FixSmallSolid, EdgeConnect, EdgeProjAux,
//! FaceConnect, ComposeShell) follow in later tranches; the GAP carriers
//! for Face/Shell/Wire/Shape in `shape_fix_gap_deps.rs` remain annotated
//! for those tranches.

pub mod edge;
pub mod free_bounds;
pub mod root;
pub mod shape_fix;
pub mod shape_fix_gap_deps;
pub mod shape_tolerance;
pub mod split_common_vertex;
pub mod wire_segment;
pub mod wire_vertex;
