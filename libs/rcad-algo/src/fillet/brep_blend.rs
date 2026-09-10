//! OCCT Blend package enums (TKFillet/Blend) — 1:1 port of
//! Blend_Status.hxx (L20-29) and Blend_DecrochStatus.hxx (L20-25).

/// OCCT Blend_Status (Blend_Status.hxx L20-29) — status of a walking step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendStatus {
    StepTooLarge,
    StepTooSmall,
    Backward,
    SamePoints,
    OnRst1,
    OnRst2,
    OnRst12,
    Ok,
}

/// OCCT Blend_DecrochStatus (Blend_DecrochStatus.hxx L20-25) — status of the
/// decrochement (shrinkage) test on a restriction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendDecrochStatus {
    NoDecroch,
    DecrochRst1,
    DecrochRst2,
    DecrochBoth,
}
