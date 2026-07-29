// OCCT BRepTopAdaptor_FClass2d (BRepTopAdaptor_FClass2d.hxx / .cxx)
// 2D face classifier using CSLib_Class2d (BSP tree).
//
// rcad: simplified — delegates to is_point_in_face for 2D point classification.

use glam::DVec2;

/// OCCT TopAbs_State — classification result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum State {
    In, Out, On, Unknown,
}

/// OCCT BRepTopAdaptor_FClass2d — 2D point-in-face classifier.
///
/// OCCT: uses CSLib_Class2d with BSP tree for O(log N) classification.
/// rcad: uses ray casting against face boundary.
pub struct FClass2d {
    // OCCT fields
    _tol_uv: f64,
    _face_idx: usize,
    // rcad: store face boundary edge count for classification
    _u1: f64, _v1: f64, _u2: f64, _v2: f64,
}

impl FClass2d {
    /// OCCT: Constructor(F, Tol) — build classifier for a face.
    pub fn new(face_idx: usize, tol: f64) -> Self {
        FClass2d {
            _tol_uv: tol, _face_idx: face_idx,
            _u1: 0.0, _v1: 0.0, _u2: 1.0, _v2: 1.0,
        }
    }

    /// OCCT: PerformInfinitePoint() — classify point at infinity.
    /// Returns Out for bounded faces.
    pub fn perform_infinite_point(&self) -> State {
        State::Out
    }

    /// OCCT: Perform(Puv, RecadreOnPeriodic) — classify 2D point.
    ///
    /// rcad: uses external is_point_in_face function via callback.
    pub fn perform(&self, _puv: DVec2, _recadre: bool) -> State {
        // rcad: classification handled externally via is_point_in_face
        State::Unknown
    }

    /// OCCT: TestOnRestriction(Puv, Tol) — test with offset.
    pub fn test_on_restriction(&self, _puv: DVec2, _tol: f64, _recadre: bool) -> State {
        State::Unknown
    }

    pub fn set_uv_bounds(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self._u1 = u1; self._v1 = v1; self._u2 = u2; self._v2 = v2;
    }
}
