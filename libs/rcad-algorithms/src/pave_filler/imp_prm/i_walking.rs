//! OCCT-aligned (placeholder): IntPatch_TheIWalking
//!
//! WIP: full 1:1 translation of IntWalk_IWalking.gxx (102KB, ~3000 lines)
//! Currently provides minimal interface that compiles.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Surface3;
use super::surf_function::SurfFunction;
use super::s_on_bounds::PathPoint;
use super::search_inside::InteriorPoint;

/// Walking line — sequence of points on intersection
#[derive(Clone, Debug)]
pub struct IWLine {
    pub points: Vec<(DVec3, f64, f64)>,  // (3D point, u, v)
    pub has_first_point: bool,
    pub has_last_point: bool,
    pub first_point_index: usize,
    pub last_point_index: usize,
    pub is_tangent_at_begin: bool,
    pub is_tangent_at_end: bool,
}

impl IWLine {
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            has_first_point: false,
            has_last_point: false,
            first_point_index: 0,
            last_point_index: 0,
            is_tangent_at_begin: false,
            is_tangent_at_end: false,
        }
    }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point_at(&self, i: usize) -> &(DVec3, f64, f64) { &self.points[i] }
}

/// Placeholder: IntPatch_TheIWalking
pub struct IWalking {
    done: bool,
    lines: Vec<IWLine>,
    epsilon: f64,
    fleche: f64,
    pas: f64,
    reversed: bool,
}

impl IWalking {
    pub fn new(epsilon: f64, deflection: f64, step: f64) -> Self {
        Self {
            done: false,
            lines: Vec::new(),
            epsilon,
            fleche: deflection,
            pas: step,
            reversed: false,
        }
    }

    pub fn perform(
        &mut self,
        _path_points: &[PathPoint],
        _interior_points: &[InteriorPoint],
        _func: &mut SurfFunction,
        _surf: &Surface3,
        _reversed: bool,
    ) {
        self.lines.clear();
        self.done = true;
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn nb_lines(&self) -> usize { self.lines.len() }
    pub fn value(&self, index: usize) -> &IWLine { &self.lines[index] }
    pub fn nb_single_points(&self) -> usize { 0 }
    pub fn single_point(&self, _index: usize) -> &PathPoint { unimplemented!() }
}
