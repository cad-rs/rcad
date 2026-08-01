// OCCT Bnd_Box2d + BndLib_Add2dCurve (Bnd package) — 2D bounding box for the
// BRepClass_Intersector fast rejection tests.
//
// Used by BRepClass_Intersector::Perform (BRepClass_Intersector.cxx L348-356,
// L392-396): the edge pcurve's bbox is tested against the ray line/segment
// before the full intersection.

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};

/// OCCT Bnd_Box2d — axis-aligned bounding box in 2D with a gap (tolerance).
pub struct BndBox2d {
    x_min: f64,
    y_min: f64,
    x_max: f64,
    y_max: f64,
    gap: f64,
    void: bool,
}

impl BndBox2d {
    pub fn new() -> Self {
        BndBox2d {
            x_min: f64::INFINITY,
            y_min: f64::INFINITY,
            x_max: f64::NEG_INFINITY,
            y_max: f64::NEG_INFINITY,
            gap: 0.0,
            void: true,
        }
    }

    pub fn is_void(&self) -> bool {
        self.void
    }

    pub fn set_gap(&mut self, tol: f64) {
        self.gap = tol.abs();
    }

    pub fn add_point(&mut self, p: DVec2) {
        if self.void {
            self.x_min = p.x;
            self.x_max = p.x;
            self.y_min = p.y;
            self.y_max = p.y;
            self.void = false;
        } else {
            self.x_min = self.x_min.min(p.x);
            self.x_max = self.x_max.max(p.x);
            self.y_min = self.y_min.min(p.y);
            self.y_max = self.y_max.max(p.y);
        }
    }

    /// OCCT BndLib_Add2dCurve::Add(C, U1, U2, Tol, Box) — build the curve's
    /// 2D bounding box over [u1, u2]. Line/circle use exact bounds; other
    /// curves are sampled (conservative for the Intersector rejection test).
    pub fn add_curve(&mut self, curve: &Curve2d, u1: f64, u2: f64, _tol: f64) {
        match curve {
            Curve2d::Line(_) => {
                self.add_point(curve.point_at(u1));
                self.add_point(curve.point_at(u2));
            }
            Curve2d::Circle(c) => {
                // Exact bbox: endpoints plus the 4 cardinal frame points when
                // inside the arc.
                self.add_point(curve.point_at(u1));
                self.add_point(curve.point_at(u2));
                for i in 0..4 {
                    let ang = i as f64 * std::f64::consts::FRAC_PI_2;
                    // The point at `ang` in the circle frame (x_dir/y_dir).
                    let p = c.center + c.x_dir * (c.radius * ang.cos()) + c.y_dir * (c.radius * ang.sin());
                    // Include only when the angle lies inside the arc domain.
                    if ang >= u1 - 1e-12 && ang <= u2 + 1e-12 {
                        self.add_point(p);
                    }
                }
            }
            _ => {
                self.add_point(curve.point_at(u1));
                self.add_point(curve.point_at(u2));
                const N: usize = 64;
                for i in 1..N {
                    let t = u1 + (u2 - u1) * (i as f64) / (N as f64);
                    self.add_point(curve.point_at(t));
                }
            }
        }
    }

    /// OCCT Bnd_Box2d::IsOut(P) — point outside the box (with gap).
    pub fn is_out_point(&self, p: DVec2) -> bool {
        if self.void {
            return true;
        }
        p.x < self.x_min - self.gap
            || p.x > self.x_max + self.gap
            || p.y < self.y_min - self.gap
            || p.y > self.y_max + self.gap
    }

    /// OCCT Bnd_Box2d::IsOut(P1, P2) — segment outside the box.
    pub fn is_out_segment(&self, p1: DVec2, p2: DVec2) -> bool {
        if self.void {
            return true;
        }
        let (lo, hi) = (self.x_min - self.gap, self.x_max + self.gap);
        let (lo_y, hi_y) = (self.y_min - self.gap, self.y_max + self.gap);
        // Reject if both endpoints are on the same side of any slab.
        (p1.x < lo && p2.x < lo)
            || (p1.x > hi && p2.x > hi)
            || (p1.y < lo_y && p2.y < lo_y)
            || (p1.y > hi_y && p2.y > hi_y)
    }

    /// OCCT Bnd_Box2d::IsOut(L) — infinite line outside the box.
    pub fn is_out_line(&self, origin: DVec2, dir: DVec2) -> bool {
        if self.void {
            return true;
        }
        let (lo, hi) = (self.x_min - self.gap, self.x_max + self.gap);
        let (lo_y, hi_y) = (self.y_min - self.gap, self.y_max + self.gap);
        // Parametric t-intervals where the line is inside each slab.
        // x = ox + dx t ∈ [lo, hi].
        let t_range = |o: f64, d: f64, a: f64, b: f64| -> Option<(f64, f64)> {
            if d.abs() < 1e-30 {
                return if o >= a && o <= b { Some((f64::NEG_INFINITY, f64::INFINITY)) } else { None };
            }
            let mut t1 = (a - o) / d;
            let mut t2 = (b - o) / d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            Some((t1, t2))
        };
        let Some(rx) = t_range(origin.x, dir.x, lo, hi) else {
            return true;
        };
        let Some(ry) = t_range(origin.y, dir.y, lo_y, hi_y) else {
            return true;
        };
        let lo_t = rx.0.max(ry.0);
        let hi_t = rx.1.min(ry.1);
        lo_t <= hi_t
    }
}
