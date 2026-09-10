//! OCCT HLRAlgo_EdgeStatus (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_EdgeStatus.hxx` (L34-117) +
//! `HLRAlgo_EdgeStatus.cxx` (L25-123).

use crate::hlr::intrv::{Intervals, Interval};

/// OCCT HLRAlgo_EdgeStatus — the Hidden Line status of an edge: its bounds
/// with tolerances, full-visible / full-hidden flags, and the sequence of
/// visible intervals.
#[derive(Debug, Clone, Default)]
pub struct EdgeStatus {
    my_start: f64,
    my_end: f64,
    my_tol_start: f32,
    my_tol_end: f32,
    my_all_hidden: bool,
    my_all_visible: bool,
    my_visibles: Intervals,
}

impl EdgeStatus {
    /// OCCT HLRAlgo_EdgeStatus() — cxx L25-33.
    pub fn new() -> Self {
        EdgeStatus {
            my_start: 0.0,
            my_end: 0.0,
            my_tol_start: 0.0,
            my_tol_end: 0.0,
            my_all_hidden: false,
            my_all_visible: false,
            my_visibles: Intervals::new(),
        }
    }

    /// OCCT HLRAlgo_EdgeStatus(Start, TolStart, End, TolEnd) — cxx L37-49:
    /// default visible.
    pub fn from_bounds(start: f64, tol_start: f32, end: f64, tol_end: f32) -> Self {
        let mut status = EdgeStatus {
            my_start: start,
            my_end: end,
            my_tol_start: tol_start,
            my_tol_end: tol_end,
            my_all_hidden: false,
            my_all_visible: false,
            my_visibles: Intervals::new(),
        };
        status.show_all();
        status
    }

    /// OCCT Initialize(Start, TolStart, End, TolEnd) — cxx L53-63: default
    /// visible.
    pub fn initialize(&mut self, start: f64, tol_start: f32, end: f64, tol_end: f32) {
        self.my_start = start;
        self.my_tol_start = tol_start;
        self.my_end = end;
        self.my_tol_end = tol_end;
        self.show_all();
    }

    /// OCCT Bounds() — hxx L57-63.
    pub fn bounds(&self) -> (f64, f32, f64, f32) {
        (
            self.my_start,
            self.my_tol_start,
            self.my_end,
            self.my_tol_end,
        )
    }

    /// OCCT NbVisiblePart() — cxx L67-81.
    pub fn nb_visible_part(&self) -> usize {
        if self.all_hidden() {
            0
        } else if self.all_visible() {
            1
        } else {
            self.my_visibles.nb_intervals()
        }
    }

    /// OCCT VisiblePart(Index) — cxx L85-99.
    pub fn visible_part(&self, index: usize) -> (f64, f32, f64, f32) {
        if self.all_visible() {
            self.bounds()
        } else {
            self.my_visibles.value(index).bounds()
        }
    }

    /// OCCT Hide(Start, TolStart, End, TolEnd, OnFace, OnBoundary) —
    /// cxx L103-123: subtracts the interval from the visible parts (the
    /// OnBoundary flag is unused in OCCT, kept in the signature).
    #[allow(clippy::too_many_arguments)]
    pub fn hide(
        &mut self,
        start: f64,
        tol_start: f32,
        end: f64,
        tol_end: f32,
        on_face: bool,
        _on_boundary: bool,
    ) {
        if !on_face {
            if self.all_visible() {
                self.my_visibles =
                    Intervals::from_interval(Interval::from_bounds_tol(
                        self.my_start,
                        self.my_tol_start,
                        self.my_end,
                        self.my_tol_end,
                    ));
                self.set_all_visible(false);
            }
            self.my_visibles
                .subtract(&Interval::from_bounds_tol(start, tol_start, end, tol_end));
            if !self.all_hidden() {
                self.set_all_hidden(self.my_visibles.nb_intervals() == 0);
            }
        }
    }

    /// OCCT HideAll() — hxx L88-92.
    pub fn hide_all(&mut self) {
        self.set_all_visible(false);
        self.set_all_hidden(true);
    }

    /// OCCT ShowAll() — hxx L95-99.
    pub fn show_all(&mut self) {
        self.set_all_visible(true);
        self.set_all_hidden(false);
    }

    /// OCCT AllHidden() — hxx L101.
    pub fn all_hidden(&self) -> bool {
        self.my_all_hidden
    }

    /// OCCT AllHidden(B) — hxx L103.
    pub fn set_all_hidden(&mut self, b: bool) {
        self.my_all_hidden = b;
    }

    /// OCCT AllVisible() — hxx L105.
    pub fn all_visible(&self) -> bool {
        self.my_all_visible
    }

    /// OCCT AllVisible(B) — hxx L107.
    pub fn set_all_visible(&mut self, b: bool) {
        self.my_all_visible = b;
    }
}
