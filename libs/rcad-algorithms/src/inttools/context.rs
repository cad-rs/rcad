/// ✅ OCCT-aligned: IntTools_Context — shared computation context with caches.
///   OCCT IntTools_Context.hxx caches: FClass2d, ProjPS, ProjPC, ProjPT,
///   SurfaceData, SolidClassifier, Hatcher, SurfaceAdaptor, OBB.
///   rcad: caches FClass2d per face (most critical for performance).
use glam::DVec2;
use crate::bopds::ds::DS;
use super::fclass2d::{FClass2d, State};

pub struct Context {
    fclass2d_cache: Vec<Option<FClass2d>>,
    tol_uv: f64,
    num_faces: usize,
}

impl Context {
    pub fn new(num_faces: usize, tol_uv: f64) -> Self {
        let mut cache = Vec::with_capacity(num_faces);
        for _ in 0..num_faces { cache.push(None); }
        Context { fclass2d_cache: cache, tol_uv, num_faces }
    }

    /// OCCT: FClass2d(theFace) — returns cached 2D point classifier for a face.
    pub fn fclass2d(&mut self, ds: &DS, face_idx: usize) -> &FClass2d {
        if self.fclass2d_cache[face_idx].is_none() {
            self.fclass2d_cache[face_idx] = Some(FClass2d::new(ds, face_idx, self.tol_uv));
        }
        self.fclass2d_cache[face_idx].as_ref().unwrap()
    }

    /// OCCT: IsPointInOnFace(theFace, theUV) — convenience wrapper.
    pub fn is_point_in_on_face(&mut self, ds: &DS, face_idx: usize, uv: DVec2) -> bool {
        self.fclass2d(ds, face_idx).perform(uv, true) != State::Out
    }

    /// OCCT: IsPointInFace(theFace, theUV) — convenience wrapper.
    pub fn is_point_in_face(&mut self, ds: &DS, face_idx: usize, uv: DVec2) -> bool {
        self.fclass2d(ds, face_idx).perform(uv, true) == State::In
    }
}
