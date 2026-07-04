/// ✅ OCCT-aligned: IntTools_Context — shared computation context with caches.
///   OCCT IntTools_Context.hxx caches: FClass2d, ProjPS, ProjPC, ProjPT,
///   SurfaceData, SolidClassifier, Hatcher, SurfaceAdaptor, OBB.
///   rcad: caches FClass2d per face and surface/curve projectors.
use glam::DVec2;
use glam::DVec3;
use crate::bopds::ds::DS;
use super::fclass2d::{FClass2d, State};
use rcad_kernel::geom::{Curve3, Surface3, SurfaceEval};
use rcad_kernel::projection::{closest_point_on_surface, closest_point_on_curve,
    SurfaceProjection, CurveProjection};

pub struct Context {
    fclass2d_cache: Vec<Option<FClass2d>>,
    tol_uv: f64,
    num_faces: usize,
    /// OCCT: ProjPS — point-on-surface projector cache (latest result per face).
    proj_ps_latest: Vec<Option<SurfaceProjection>>,
    /// OCCT: ProjPC — point-on-curve projector cache (latest result per edge).
    proj_pc_latest: Vec<Option<CurveProjection>>,
    /// OCCT: ProjPT — single reusable point-on-curve projector for transient curves.
    proj_pt_latest: Option<CurveProjection>,
    /// OCCT: UVBounds cache — precomputed UV bounds per face.
    uv_bounds_cache: Vec<Option<[f64; 4]>>,
    /// OCCT: SurfaceAdaptor — precomputed surface references.
    surface_cache: Vec<Option<Surface3>>,
}

impl Context {
    pub fn new(num_faces: usize, tol_uv: f64) -> Self {
        let init_face: Vec<Option<FClass2d>> = (0..num_faces).map(|_| None).collect();
        Context {
            fclass2d_cache: init_face,
            tol_uv,
            num_faces,
            proj_ps_latest: (0..num_faces).map(|_| None).collect(),
            proj_pc_latest: Vec::new(),
            proj_pt_latest: None,
            uv_bounds_cache: (0..num_faces).map(|_| None).collect(),
            surface_cache: (0..num_faces).map(|_| None).collect(),
        }
    }

    pub fn tol_uv(&self) -> f64 { self.tol_uv }

    /// OCCT: FClass2d(theFace) — returns cached 2D point classifier for a face.
    pub fn fclass2d(&mut self, ds: &DS, face_idx: usize) -> &FClass2d {
        assert!(face_idx < self.num_faces, "Context: face_idx {} out of range ({})", face_idx, self.num_faces);
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

    /// OCCT: IsValidPointForFace(theP, theFace, theTol) — checks if a 3D point
    /// projects onto the face surface within tolerance.
    pub fn is_valid_point_for_face(&self, p: DVec3, face_idx: usize, tol: f64) -> bool {
        if face_idx >= self.num_faces { return false; }
        if let Some(surf) = &self.surface_cache[face_idx] {
            let proj = closest_point_on_surface(surf, p, 16);
            return proj.distance < tol && proj.distance.is_finite();
        }
        false
    }

    /// OCCT: ProjPS(theFace) — projects a 3D point onto the face surface.
    /// Returns (uv, 3d_point, distance) on success.
    pub fn proj_ps(&mut self, ds: &DS, face_idx: usize, p: DVec3) -> Option<(DVec2, DVec3, f64)> {
        if face_idx >= self.num_faces { return None; }
        let surf = &ds.faces[face_idx].surface;
        let proj = closest_point_on_surface(surf, p, 16);
        if proj.distance.is_finite() {
            self.proj_ps_latest[face_idx] = Some(proj);
            let cached = self.proj_ps_latest[face_idx].as_ref().unwrap();
            Some((glam::DVec2::new(cached.params.0, cached.params.1), cached.point, cached.distance))
        } else {
            None
        }
    }

    /// OCCT: ProjPC(theEdge) — projects a 3D point onto the edge's curve.
    /// Returns (param, 3d_point, distance) on success.
    pub fn proj_pc(&mut self, ds: &DS, edge_idx: usize, p: DVec3) -> Option<(f64, DVec3, f64)> {
        if edge_idx >= ds.edges.len() { return None; }
        let curve = &ds.edges[edge_idx].curve;
        let proj = closest_point_on_curve(curve, p, 16);
        if proj.distance.is_finite() && proj.param.is_finite() {
            if edge_idx >= self.proj_pc_latest.len() {
                self.proj_pc_latest.resize_with(edge_idx + 1, || None);
            }
            self.proj_pc_latest[edge_idx] = Some(proj);
            let cached = self.proj_pc_latest[edge_idx].as_ref().unwrap();
            Some((cached.param, cached.point, cached.distance))
        } else {
            None
        }
    }

    /// OCCT: ProjPT(theP, theC) — projects a 3D point onto a transient curve.
    /// Unlike ProjPC (which is keyed by edge index), this is a single reusable
    /// projector for one-off curve projections. Returns (param, 3d_point, distance).
    pub fn proj_pt(&mut self, curve: &Curve3, p: DVec3) -> Option<(f64, DVec3, f64)> {
        let proj = closest_point_on_curve(curve, p, 16);
        if proj.distance.is_finite() && proj.param.is_finite() {
            self.proj_pt_latest = Some(proj);
            let cached = self.proj_pt_latest.as_ref().unwrap();
            Some((cached.param, cached.point, cached.distance))
        } else {
            None
        }
    }

    /// OCCT: SurfaceAdaptor(theFace) — returns cached surface reference for a face.
    pub fn surface_adaptor(&mut self, ds: &DS, face_idx: usize) -> &Surface3
    where
        Self: 'static,
    {
        assert!(face_idx < self.num_faces, "surface_adaptor: face_idx {} >= {}", face_idx, self.num_faces);
        if self.surface_cache[face_idx].is_none() {
            self.surface_cache[face_idx] = Some(ds.faces[face_idx].surface.clone());
        }
        self.surface_cache[face_idx].as_ref().unwrap()
    }

    /// OCCT: UVBounds(theFace) — returns cached UV bounds [umin, umax, vmin, vmax].
    pub fn uv_bounds(&mut self, ds: &DS, face_idx: usize) -> [f64; 4] {
        if face_idx < self.num_faces {
            if let Some(bounds) = self.uv_bounds_cache[face_idx] {
                return bounds;
            }
            let bounds = ds.faces[face_idx].surface.default_domain();
            self.uv_bounds_cache[face_idx] = Some(bounds);
            return bounds;
        }
        [0.0, std::f64::consts::TAU, 0.0, std::f64::consts::PI]
    }
}
