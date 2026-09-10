//! OCCT ShapeFix — package module.
//!
//! W1-6 pull-forward: the `ShapeFix::SameParameter` package static
//! (ShapeFix.cxx) — the DT_SplitByNumber heal oracle dependency and the
//! chfi3d same_parameter_pass backfill. The ShapeFix class translations
//! follow in W3.
//!
//! W3 tranche 1 (2026-09-10): ShapeFix_Root, ShapeFix_ShapeTolerance,
//! ShapeFix_Edge, ShapeFix_WireSegment, ShapeFix_WireVertex,
//! ShapeFix_SplitCommonVertex, ShapeFix_FreeBounds.
//!
//! W3 tranche 2 (2026-09-10): ShapeFix_Wire (+ Wire_1 + lxx) — the full
//! wire repair tool in `wire/` (mod + fix_api + fix_adv + fix_gaps +
//! wire_statics).
//!
//! W3 tranche 3 (2026-09-10): ShapeFix_Face (`face_a/face_b/face_c`),
//! ShapeFix_Shell (`shell` + `shell_statics`), ShapeFix_IntersectionTool
//! (`intersection_tool`) and ShapeFix_SplitTool (`split_tool`) — the
//! `ShapeFixIntersectionToolGap` / `ShapeFixSplitToolGap` /
//! `ShapeFixFaceGap` / `ShapeFixShellGap` carriers are retired (Rule 4).
//!
//! W3 tranche 4 (2026-09-10): ShapeFix_Shape (`shape_fix_shape`, cxx + lxx)
//! — the `ShapeFixShapeGap` carrier is retired (Rule 4) and the
//! `ShapeFix::RemoveSmallEdges` static drives the real class; the not-yet-
//! landed ShapeFix_Solid row is hosted by the `ShapeFixSolidGap` carrier in
//! `shape_fix_gap_deps.rs` (the real `ShapeFixShell` is embedded, so the
//! tool chain stays real).  The remaining W3 classes (Wireframe,
//! FixSmallFace, FixSmallSolid, EdgeConnect, EdgeProjAux, FaceConnect,
//! ComposeShell) follow in later tranches.

pub mod edge;
pub mod face_a;
pub mod face_b;
pub mod face_c;
pub mod free_bounds;
pub mod intersection_tool;
pub mod intersection_tool_fix;
pub mod root;
pub mod shape_fix;
pub mod shape_fix_gap_deps;
pub mod shape_fix_shape;
pub mod shape_tolerance;
pub mod shell;
pub mod shell_statics;
pub mod split_common_vertex;
pub mod split_tool;
pub mod wire;
pub mod wire_segment;
pub mod wire_vertex;
