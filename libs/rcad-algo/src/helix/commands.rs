//! OCCT BRepTest_HelixCommands.cxx (Draw TKTopTest) — the DRAW command layer
//! over HelixBRep_BuilderHelix.
//!
//! 1:1 translation of the `helix` / `comphelix` / `spiral` / `comphelix2` /
//! `helix2` / `spiral2` commands (L117-547), including the static
//! `theHelixAxis = gp_Ax3(P0, DZ, OX)` default (L48).  Each function mirrors
//! the argument parsing and `SetParameters` overload used by the DRAW
//! command; `DisplayHelixResult` maps to the returned `Result` (non-zero
//! ErrorStatus = OCCT `catch` failure, WarningStatus/ToleranceReached stay
//! reachable on the builder).

use glam::DVec3;
use rcad_kernel::math::gp::Ax3;

use super::helix_brep::BuilderHelix;
use rcad_kernel::topo::topods::BRep;

/// OCCT static theHelixAxis (BRepTest_HelixCommands.cxx L48).
pub fn the_helix_axis() -> Ax3 {
    Ax3::from_pnt_n_vx(DVec3::ZERO, DVec3::Z, DVec3::X)
}

/// OCCT `comphelix name np D1 D2 [Di...] H1 [Hi...] P1 [Pi...] PF1 [PFi...]`
/// (L144-216).
#[allow(clippy::too_many_arguments)]
pub fn comphelix(
    diams: &[f64],
    heights: &[f64],
    pitches: &[f64],
    is_pitches: &[bool],
) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters(&the_helix_axis(), diams, heights, pitches, is_pitches);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT `helix name np D1 H1 [Hi...] P1 [Pi...] PF1 [PFi...]` (L220-289).
pub fn helix(diam: f64, heights: &[f64], pitches: &[f64], is_pitches: &[bool]) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters_helix(&the_helix_axis(), diam, heights, pitches, is_pitches);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT `spiral name np D1 D2 H1 [Hi...] P1 [Pi...] PF1 [PFi...]` (L293-365).
#[allow(clippy::too_many_arguments)]
pub fn spiral(
    diam1: f64,
    diam2: f64,
    heights: &[f64],
    pitches: &[f64],
    is_pitches: &[bool],
) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters_spiral(&the_helix_axis(), diam1, diam2, heights, pitches, is_pitches);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT `comphelix2 name np D1 D2 [Di...] P1 [Pi...] N1 [Ni...]` (L369-428).
pub fn comphelix2(diams: &[f64], pitches: &[f64], nb_turns: &[f64]) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters_turns(&the_helix_axis(), diams, pitches, nb_turns);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT `helix2 name np D1 P1 [Pi...] N1 [Ni...]` (L432-486).
pub fn helix2(diam: f64, pitches: &[f64], nb_turns: &[f64]) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters_helix_turns(&the_helix_axis(), diam, pitches, nb_turns);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT `spiral2 name np D1 D2 P1 [Pi...] N1 [Ni...]` (L490-547).
pub fn spiral2(diam1: f64, diam2: f64, pitches: &[f64], nb_turns: &[f64]) -> Result<BRep, i32> {
    let mut a_bh = BuilderHelix::new();
    a_bh.set_parameters_spiral_turns(&the_helix_axis(), diam1, diam2, pitches, nb_turns);
    a_bh.perform();
    finish(a_bh)
}

/// OCCT DisplayHelixResult (L59-75): on ErrorStatus == 0 the built shape is
/// registered; otherwise the error status is reported (`catch` fires).
fn finish(a_bh: BuilderHelix) -> Result<BRep, i32> {
    if a_bh.error_status() == 0 {
        Ok(a_bh.into_brep())
    } else {
        Err(a_bh.error_status())
    }
}
