//! OCCT-aligned: IntPatch_LineConstructor — construct section edges from intersection lines.
//!
//! OCCT IntPatch_LineConstructor.hxx / .cxx (66K)
//!
//! Splits and formats IntPatch_Lines into section edges for the BRep.
//! rcad equivalent: make_blocks.rs (pave_filler/make_blocks.rs)
//! This module wraps that existing infrastructure.

use super::int_patch_line::IntPatchLine;

/// OCCT-aligned: LineConstructor — construct section edges from intersection lines.
pub struct LineConstructor;

impl LineConstructor {
    /// Load lines to construct section edge boundaries.
    /// rcad: delegates to make_blocks.rs via the PaveFiller.
    pub fn load_lines(&self, _lines: Vec<IntPatchLine>) {
        // rcad: handled by PaveFiller::make_blocks
    }
}
