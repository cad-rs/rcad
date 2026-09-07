//! OCCT TopOpeBRepBuild_HBuilder facade for the ChFi3d call surface.
//!
//! User ruling 2026-09-07: TKBool code (TopOpeBRepBuild / TopOpeBRepDS /
//! TopOpeBRepTool) is NOT translated 1:1.  Wherever TKFillet / TKOffset /
//! TKFeat depend on it, the rcad implementation goes through the aligned
//! TKBO pipeline (`bop::algo` / `bop::ds` / `bop::brep_algo_api`).
//!
//! The HBuilder surface consumed by ChFi3d (Builder.hxx) is the 7-method
//! set Perform / MergeSolid / IsSplit / Splits / Merged / NewEdges /
//! NewFaces plus the `Builder()` accessor (ChFi3d_Builder.hxx L180).  The
//! TKBO wiring of those methods lands in Stage 1g together with the ChFi3d
//! Perform flow that calls them; until then this file carries only the
//! member state the current ChFi3d call surface touches (construction in
//! `ChFi3d_Builder::ChFi3d_Builder`, return of `Builder()`).

/// OCCT TopOpeBRepBuild_HBuilder (TopOpeBRepBuild_Builder handle wrapper).
/// TKBO-backed facade; see the module doc for the ruling.
#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepBuildHBuilder {}
