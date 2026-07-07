//! OCCT-aligned (placeholder): IntPatch_TheSearchInside
//!
//! WIP: full 1:1 translation of IntStart_SearchInside.gxx (10K)
//! Currently provides minimal interface that compiles.

use glam::{DVec2, DVec3};
use super::surf_function::SurfFunction;

/// Interior point (OCCT IntSurf_InteriorPoint)
#[derive(Clone, Debug)]
pub struct InteriorPoint {
    pub value: DVec3,
    pub u: f64,
    pub v: f64,
    pub direction: DVec3,
    pub direction_2d: DVec2,
}

/// Placeholder: IntPatch_TheSearchInside
pub struct SearchInside {
    done: bool,
    list: Vec<InteriorPoint>,
}

impl SearchInside {
    pub fn new() -> Self {
        Self { done: false, list: Vec::new() }
    }

    pub fn perform(
        &mut self,
        _func: &mut SurfFunction,
        _u_min: f64, _u_max: f64,
        _v_min: f64, _v_max: f64,
        _epsilon: f64,
    ) {
        self.list.clear();
        self.done = true;
    }

    pub fn perform_from_point(
        &mut self,
        _func: &mut SurfFunction,
        _u_start: f64,
        _v_start: f64,
    ) {
        self.list.clear();
        self.done = true;
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn nb_points(&self) -> usize { self.list.len() }
    pub fn value(&self, index: usize) -> &InteriorPoint { &self.list[index] }
}
