// OCCT IntTools_EdgeFace — edge-face intersection.
//
// Finds intersection points between an edge curve and a face surface.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, CurveEval, Surface3};

/// Intersection hit between edge and face.
#[derive(Debug, Clone)]
pub struct EdgeFaceHit {
    pub point: DVec3,
    pub edge_param: f64,
    pub uv: DVec2,
    pub distance: f64,
}

/// IntTools_EdgeFace — edge-face intersection engine.
pub struct EdgeFace {
    curve: Curve3,
    range: [f64; 2],
    surface: Surface3,
    tol: f64,
    hits: Vec<EdgeFaceHit>,
    done: bool,
}

impl EdgeFace {
    pub fn new() -> Self {
        EdgeFace {
            curve: rcad_kernel::geom::Curve3::Line(
                rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }
            ),
            range: [0.0, 1.0],
            surface: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: glam::DVec3::ZERO, normal: glam::DVec3::Z,
                u_dir: glam::DVec3::X, v_dir: glam::DVec3::Y,
            }),
            tol: 1e-7,
            hits: Vec::new(),
            done: false,
        }
    }

    pub fn set_curve(&mut self, c: Curve3, r: [f64; 2]) { self.curve = c; self.range = r; }
    pub fn set_face(&mut self, s: Surface3) { self.surface = s; }
    pub fn set_tolerances(&mut self, t: f64) { self.tol = t.max(1e-7); }
    pub fn is_done(&self) -> bool { self.done }
    pub fn hits(&self) -> &[EdgeFaceHit] { &self.hits }
    pub fn common_parts(&self) -> &[EdgeFaceHit] { &self.hits }

    /// Perform edge-face intersection by sampling the curve and projecting to the surface.
    pub fn perform(&mut self) {
        self.hits.clear();
        let n_samples = 64usize;
        let tol = self.tol;
        let mut prev_uv: Option<DVec2> = None;

        for i in 0..=n_samples {
            let t = self.range[0] + (self.range[1] - self.range[0]) * i as f64 / n_samples as f64;
            let pt = self.curve.point_at(t);
            let proj = rcad_kernel::base::geom_api::project::closest_point_on_surface(&self.surface, pt, 64);
            let dist = proj.distance;

            if dist <= tol {
                // Check if this is a new hit (not consecutive with previous)
                let is_new = match prev_uv {
                    Some(puv) => (proj.params.0 - puv.x).abs() > 1e-6 || (proj.params.1 - puv.y).abs() > 1e-6,
                    None => true,
                };
                if is_new || self.hits.is_empty() {
                    self.hits.push(EdgeFaceHit {
                        point: proj.point,
                        edge_param: t,
                        uv: DVec2::new(proj.params.0, proj.params.1),
                        distance: dist,
                    });
                }
                prev_uv = Some(DVec2::new(proj.params.0, proj.params.1));
            }
        }
        self.done = true;
    }
}

impl Default for EdgeFace { fn default() -> Self { Self::new() } }
