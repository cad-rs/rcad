//! OCCT-aligned (placeholder): IntPatch_TheSOnBounds
//!
//! WIP: full 1:1 translation of IntStart_SearchOnBoundaries.gxx (37K)
//! Currently provides minimal interface that compiles.

use glam::DVec3;
use super::arc_function::ArcFunction;

/// Boundary path point (OCCT IntPatch_ThePathPointOfTheSOnBounds)
#[derive(Clone, Debug)]
pub struct PathPoint {
    pub value: DVec3,
    pub parameter: f64,
    pub arc_index: usize,
    pub is_new: bool,
    pub tolerance: f64,
}

/// Boundary segment (OCCT IntPatch_TheSegmentOfTheSOnBounds)
#[derive(Clone, Debug)]
pub struct Segment {
    pub curve: rcad_kernel::geom::Curve2d,
    pub first_index: usize,
    pub last_index: usize,
}

/// Placeholder: IntPatch_TheSOnBounds
pub struct SOnBounds {
    done: bool,
    all: bool,
    points: Vec<PathPoint>,
    segments: Vec<Segment>,
}

impl SOnBounds {
    pub fn new() -> Self {
        Self { done: false, all: false, points: Vec::new(), segments: Vec::new() }
    }

    pub fn perform(
        &mut self,
        _func: &mut ArcFunction,
        _u_min: f64, _u_max: f64,
        _v_min: f64, _v_max: f64,
        _tol_boundary: f64,
        _tol_tangency: f64,
    ) {
        self.points.clear();
        self.segments.clear();
        self.done = true;
        self.all = false;
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn all_arc_solution(&self) -> bool { self.all }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point(&self, index: usize) -> &PathPoint { &self.points[index] }
    pub fn nb_segments(&self) -> usize { self.segments.len() }
    pub fn segment(&self, index: usize) -> &Segment { &self.segments[index] }
}
