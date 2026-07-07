//! OCCT-aligned: IntPatch_PrmPrmIntersection — parametric-parametric surface intersection.
//!
//! OCCT algorithm (137K, ~3500 lines):
//!   1. Build polyhedron approximations for both surfaces
//!   2. Find intersecting triangles → starting points
//!   3. Walk from each starting point along intersection
//!   4. Join/merge resulting Walking lines
//!
//! rcad: delegates to face_face::intersect_faces / marching / intss

use rcad_kernel::geom::{Surface3, SurfaceEval};
use super::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use super::int_patch_type::IntPatchIType;

pub struct PrmPrmIntersection {
    done: bool, empt: bool,
    spnt: Vec<super::int_patch_point::IntPatchPoint>,
    slin: Vec<IntPatchLine>,
}

impl PrmPrmIntersection {
    pub fn new() -> Self {
        Self { done: false, empt: true, spnt: Vec::new(), slin: Vec::new() }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn is_empty(&self) -> bool { self.empt }
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    pub fn slin_ref(&self) -> &[IntPatchLine] { &self.slin }

    /// OCCT: Perform — intersect two parametric surfaces.
    pub fn perform(&mut self, s1: &Surface3, s2: &Surface3,
                   tol_arc: f64, tol_tang: f64, _fleche: f64, _uv_max_step: f64) {
        self.done = false; self.empt = true; self.slin.clear(); self.spnt.clear();

        // Use marching/numeric intersection from intss via face_face
        let curves = crate::inttools::face_face::intersect_faces(s1, s2, tol_arc, tol_tang);
        self.slin = curves.into_iter().map(|c| IntPatchLine {
            line_type: IntPatchIType::Walking, curve: c.curve, t_range: c.t_range,
            pcurve1: c.pcurve1, pcurve2: c.pcurve2,
            tolerance: c.tolerance, tang_tolerance: c.tang_tolerance,
            wline_pnts: Vec::new(), is_purging_allowed: false,
            wl_type: crate::inttools::int_patch_line::WLineType::PrmPrm,
        }).collect();

        // If no analytic curves, try grid-based marching
        if self.slin.is_empty() {
            let pnts = self.sample_intersection(s1, s2, tol_arc);
            if pnts.len() >= 2 {
                self.slin.push(IntPatchLine::walking(pnts, WLineType::PrmPrm));
            }
        }

        self.empt = self.slin.is_empty();
        self.done = true;
    }

    /// Grid sampling to find zero crossings of distance function.
    fn sample_intersection(&self, s1: &Surface3, s2: &Surface3, tol: f64) -> Vec<WLinePnt> {
        let n_u = 40;
        let n_v = 40;
        let mut points = Vec::new();

        for i in 0..n_u {
            for j in 0..n_v-1 {
                let u = i as f64 / n_u as f64;
                let v1 = j as f64 / n_v as f64;
                let v2 = (j+1) as f64 / n_v as f64;
                let p1 = s1.point_at(u, v1);
                let p2 = s1.point_at(u, v2);
                // Distance between surfaces: project s1 point onto s2
                let d1 = surface_distance(s1, s2, u, v1);
                let d2 = surface_distance(s1, s2, u, v2);
                if d1 * d2 < 0.0 || d1.abs() < tol || d2.abs() < tol {
                    let t = d1.abs() / (d1.abs() + d2.abs()).max(1e-30);
                    let v = v1 + t * (v2 - v1);
                    points.push(WLinePnt { p3d: s1.point_at(u, v), u1: u, v1: v, u2: u, v2: v });
                }
            }
        }
        points
    }
}

/// Signed distance: positive if s1 point is outside s2's volume, negative if inside.
/// Approximation: project point onto s2 surface, check normal orientation.
fn surface_distance(s1: &Surface3, s2: &Surface3, u: f64, v: f64) -> f64 {
    let p1 = s1.point_at(u, v);
    let proj = rcad_kernel::projection::closest_point_on_surface(s2, p1, 16);
    let proj_pnt = proj.point;
    let d = (p1 - proj_pnt).length();
    let n = s2.normal_at(0.5, 0.5);
    let dir = p1 - proj_pnt;
    if dir.dot(n) >= 0.0 { d } else { -d }
}
