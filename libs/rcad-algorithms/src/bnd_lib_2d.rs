//! BndLib-style bounding algorithms for 2D geometry.
//!
//! Provides analytical bounding-box computation for 2D curves,
//! matching OCCT's `BndLib` (2D overloads) and `GeomBndLib_Curve2d` dispatch.
//!
//! ✅ OCCT-aligned: BndLib.hxx (2D), GeomBndLib_Curve2d
//!
//! | Type           | Method                                        |
//! |----------------|-----------------------------------------------|
//! | Line2d         | Endpoints (+ infinite handling)               |
//! | Circle2d       | Per-coordinate analytical (sqrt(R²·Xk²+R²·Yk²)) |
//! | Ellipse2d      | Per-coordinate analytical (major/minor axes)  |
//! | Parabola2d     | Endpoints + apex + infinite handling          |
//! | Hyperbola2d    | Endpoints + per-coordinate extremal via log   |
//! | BSplineCurve2  | Knot-interval sampling + control-point hull   |
//! | BezierCurve2   | Sampling + control-point hull                 |
//! | OffsetCurve2d  | Base curve box + |offset|                     |
//! | OtherCurve2d   | Adaptive sampling with deflection             |

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::geom::{
    BezierCurve2, BSplineCurve2, Circle2d, Ellipse2d, Hyperbola2d,
    Line2d, OffsetCurve2d, Parabola2d,
};

// =============================================================================
// BoundingBox2d — 2D Axis-Aligned Bounding Box
// =============================================================================
// ✅ OCCT-aligned: Bnd_Box2d (simplified — no open/whole-space flags)

/// A 2D axis-aligned bounding box.
///
/// OCCT-aligned: Bnd_Box2d.  Provides `Add(point)`, `Enlarge(tol)`,
/// `IsVoid()`, `Get()` accessors, and `from_corners`/`center`/etc.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox2d {
    /// Minimum corner.
    pub min: DVec2,
    /// Maximum corner.
    pub max: DVec2,
}

impl Default for BoundingBox2d {
    fn default() -> Self { Self::new() }
}

impl BoundingBox2d {
    /// Create a new (empty) bounding box.
    pub fn new() -> Self {
        Self {
            min: DVec2::splat(f64::INFINITY),
            max: DVec2::splat(f64::NEG_INFINITY),
        }
    }

    /// Create from a single point.
    pub fn from_point(p: DVec2) -> Self { Self { min: p, max: p } }

    /// Create from two corners (order agnostic).
    pub fn from_corners(p1: DVec2, p2: DVec2) -> Self {
        Self { min: p1.min(p2), max: p1.max(p2) }
    }

    /// Create from known-ordered min/max.
    pub fn from_min_max(min: DVec2, max: DVec2) -> Self { Self { min, max } }

    /// Whether the box contains at least one finite point.
    pub fn is_valid(&self) -> bool {
        self.min.x <= self.max.x
            && self.min.y <= self.max.y
            && self.min.x.is_finite()
            && self.max.x.is_finite()
    }

    /// Whether the box is empty.
    pub fn is_empty(&self) -> bool { !self.is_valid() }

    /// Expand to include a point.
    pub fn add_point(&mut self, p: DVec2) {
        self.min = self.min.min(p);
        self.max = self.max.max(p);
    }

    /// Expand to include another box.
    pub fn add_bbox(&mut self, other: &BoundingBox2d) {
        if other.is_valid() {
            self.min = self.min.min(other.min);
            self.max = self.max.max(other.max);
        }
    }

    /// Get corners as `(xmin, ymin, xmax, ymax)`.
    pub fn get(&self) -> (f64, f64, f64, f64) {
        (self.min.x, self.min.y, self.max.x, self.max.y)
    }

    /// Enlarge by tolerance on all sides.
    pub fn enlarge(&mut self, tol: f64) {
        if self.is_valid() {
            self.min -= DVec2::splat(tol);
            self.max += DVec2::splat(tol);
        }
    }

    /// Center point.
    pub fn center(&self) -> DVec2 { (self.min + self.max) * 0.5 }

    /// Size (extent) vector.
    pub fn size(&self) -> DVec2 { self.max - self.min }

    /// Width (X extent).
    pub fn width(&self) -> f64 { self.max.x - self.min.x }

    /// Height (Y extent).
    pub fn height(&self) -> f64 { self.max.y - self.min.y }

    /// Diagonal length.
    pub fn diagonal(&self) -> f64 { (self.max - self.min).length() }

    /// Area.
    pub fn area(&self) -> f64 {
        let s = self.size();
        if s.x > 0.0 && s.y > 0.0 { s.x * s.y } else { 0.0 }
    }
}

// =============================================================================
// Public API
// =============================================================================

/// Default domain for a Curve2d.
fn curve2d_default_domain(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Line(_) | Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => {
            [f64::NEG_INFINITY, f64::INFINITY]
        }
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
            [0.0, 2.0 * std::f64::consts::PI]
        }
        Curve2d::BSpline(c) => {
            let d = c.degree;
            let n = c.knots.len();
            if n > 2 * d { [c.knots[d], c.knots[n - d - 1]] }
            else { [c.knots[0], c.knots[n - 1]] }
        }
        Curve2d::Bezier(_) => [0.0, 1.0],
        Curve2d::Trimmed(c) => [c.t_min, c.t_max],
        // Offset delegates to a sampling fallback
        _ => [0.0, 1.0],
    }
}

/// Add a 2D curve to a bounding box over `[t1, t2]`.
pub fn add_curve_to_box(curve: &Curve2d, t1: f64, t2: f64, bbox: &mut BoundingBox2d, tol: f64) {
    let b = curve_box(curve, t1, t2, tol);
    if b.is_valid() {
        bbox.add_bbox(&b);
    }
}

/// Add optimal (tighter) bounding box for a 2D curve.
pub fn add_curve_optimal(curve: &Curve2d, t1: f64, t2: f64, bbox: &mut BoundingBox2d, tol: f64) {
    let b = curve_box_optimal(curve, t1, t2, tol);
    if b.is_valid() {
        bbox.add_bbox(&b);
    }
}

/// Compute bounding box for a 2D curve over its default domain.
pub fn curve_bounds(curve: &Curve2d, tol: f64) -> BoundingBox2d {
    let [t1, t2] = curve2d_default_domain(curve);
    curve_box(curve, t1, t2, tol)
}

/// Compute bounding box for a 2D curve over `[t1, t2]`.
pub fn curve_bounds_with_range(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    curve_box(curve, t1, t2, tol)
}

/// Compute optimal bounding box for a 2D curve over `[t1, t2]`.
pub fn curve_bounds_optimal(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    curve_box_optimal(curve, t1, t2, tol)
}

// =============================================================================
// Internal dispatch
// =============================================================================

fn curve_box(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    match curve {
        Curve2d::Line(c)    => line_box(c, t1, t2, tol),
        Curve2d::Circle(c)  => circle_box(c, t1, t2, tol),
        Curve2d::Ellipse(c) => ellipse_box(c, t1, t2, tol),
        Curve2d::Parabola(c)=> parabola_box(c, t1, t2, tol),
        Curve2d::Hyperbola(c)=> hyperbola_box(c, t1, t2, tol),
        Curve2d::BSpline(c) => bspline_curve_box(c, t1, t2, tol),
        Curve2d::Bezier(c)  => bezier_curve_box(c, t1, t2, tol),
        Curve2d::Offset(c)  => offset_curve_box(c, t1, t2, tol),
        _                   => other_curve_box(curve, t1, t2, tol),
    }
}

fn curve_box_optimal(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    match curve {
        Curve2d::Line(c)    => line_box(c, t1, t2, tol),
        Curve2d::Circle(c)  => circle_box(c, t1, t2, tol),
        Curve2d::Ellipse(c) => ellipse_box(c, t1, t2, tol),
        Curve2d::Parabola(c)=> parabola_box(c, t1, t2, tol),
        Curve2d::Hyperbola(c)=> hyperbola_box(c, t1, t2, tol),
        Curve2d::BSpline(c) => bspline_curve_box_optimal(c, t1, t2, tol),
        Curve2d::Bezier(c)  => bezier_curve_box_optimal(c, t1, t2, tol),
        Curve2d::Offset(c)  => offset_curve_box(c, t1, t2, tol),
        _                   => other_curve_box_optimal(curve, t1, t2, tol),
    }
}

// =============================================================================
// Analytical evaluators — 2D, matching GeomBndLib_Curve2d
// =============================================================================

fn eval_point(curve: &Curve2d, t: f64) -> DVec2 { curve.point_at(t) }

// ── Line2d ────────────────────────────────────────────────────────────────

fn line_box(line: &Line2d, mut t1: f64, mut t2: f64, tol: f64) -> BoundingBox2d {
    use rcad_kernel::geom::Curve2dEval;
    let mut bbox = BoundingBox2d::new();
    if t1.is_infinite() || t2.is_infinite() {
        // rcad BoundingBox2d has no open flags; add finite endpoints only.
        if t1.is_infinite() { t1 = 0.0; }
        if t2.is_infinite() { t2 = 0.0; }
    }
    let p1 = line.point_at(t1);
    let p2 = line.point_at(t2);
    bbox.add_point(p1);
    bbox.add_point(p2);
    bbox.enlarge(tol);
    bbox
}

// ── Circle2d ──────────────────────────────────────────────────────────────
// ✅ OCCT-aligned: GeomBndLib_Circle2d (same formula as 3D, 2 coordinates)

fn circle_box(circle: &Circle2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let r = circle.radius;
    let o = circle.center;
    let xd = circle.x_dir;
    let yd = circle.y_dir;

    let period = 2.0 * std::f64::consts::PI;
    if t2 - t1 >= period {
        return full_circle_box(o, xd, yd, r, tol);
    }

    let u1 = t1.rem_euclid(period);
    let u2 = t2.rem_euclid(period);
    let (a1, a2) = if u2 < u1 { (u2, u1) } else { (u1, u2) };

    let mut bbox = BoundingBox2d::new();

    let eval = |t: f64| -> DVec2 { o + r * (t.cos() * xd + t.sin() * yd) };
    bbox.add_point(eval(a1));
    bbox.add_point(eval(a2));

    for k in 0..2 {
        let xk = xd[k];
        let yk = yd[k];
        let t_extr = if xk.abs() > 1e-15 { yk.atan2(xk) }
                     else { std::f64::consts::PI / 2.0 };
        let t_extr2 = {
            let v = t_extr + std::f64::consts::PI;
            if v > period { v - period } else { v }
        };
        for &te in &[t_extr, t_extr2] {
            let tk = if te >= a1 { te } else { te + period };
            if tk >= a1 && tk <= a2 {
                bbox.add_point(eval(te));
            }
        }
    }
    bbox.enlarge(tol);
    bbox
}

fn full_circle_box(o: DVec2, xd: DVec2, yd: DVec2, r: f64, tol: f64) -> BoundingBox2d {
    let mut bmin = [f64::INFINITY; 2];
    let mut bmax = [f64::NEG_INFINITY; 2];
    for k in 0..2 {
        let xk = xd[k];
        let yk = yd[k];
        let amp = (r * r * xk * xk + r * r * yk * yk).sqrt();
        bmin[k] = o[k] - amp;
        bmax[k] = o[k] + amp;
    }
    let mut bbox = BoundingBox2d::from_corners(
        DVec2::new(bmin[0], bmin[1]),
        DVec2::new(bmax[0], bmax[1]),
    );
    bbox.enlarge(tol);
    bbox
}

// ── Ellipse2d ─────────────────────────────────────────────────────────────

fn ellipse_box(ell: &Ellipse2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let a_maj = ell.major_radius;
    let a_min = ell.minor_radius;
    let o = ell.center;
    // Perp: rotate_ccw_90(major_dir) = (-major_dir.y, major_dir.x)
    let xd = ell.major_dir;
    let yd = DVec2::new(-ell.major_dir.y, ell.major_dir.x);

    let period = 2.0 * std::f64::consts::PI;
    if t2 - t1 >= period {
        return full_ellipse_box(o, xd, yd, a_maj, a_min, tol);
    }

    let u1 = t1.rem_euclid(period);
    let u2 = t2.rem_euclid(period);
    let (a1, a2) = if u2 < u1 { (u2, u1) } else { (u1, u2) };

    let mut bbox = BoundingBox2d::new();
    let eval = |t: f64| -> DVec2 { o + a_maj * t.cos() * xd + a_min * t.sin() * yd };
    bbox.add_point(eval(a1));
    bbox.add_point(eval(a2));

    for k in 0..2 {
        let xk = xd[k];
        let yk = yd[k];
        let t_extr = if xk.abs() > 1e-15 { (a_min * yk).atan2(a_maj * xk) }
                     else { std::f64::consts::PI / 2.0 };
        let t_extr2 = if t_extr <= std::f64::consts::PI {
            t_extr + std::f64::consts::PI
        } else {
            t_extr - std::f64::consts::PI
        };
        for &te in &[t_extr, t_extr2] {
            let tk = if te >= a1 { te } else { te + period };
            if tk >= a1 && tk <= a2 {
                bbox.add_point(eval(te));
            }
        }
    }
    bbox.enlarge(tol);
    bbox
}

fn full_ellipse_box(o: DVec2, xd: DVec2, yd: DVec2, a_maj: f64, a_min: f64, tol: f64) -> BoundingBox2d {
    let mut bmin = [f64::INFINITY; 2];
    let mut bmax = [f64::NEG_INFINITY; 2];
    for k in 0..2 {
        let xk = xd[k];
        let yk = yd[k];
        let amp = (a_maj * a_maj * xk * xk + a_min * a_min * yk * yk).sqrt();
        bmin[k] = o[k] - amp;
        bmax[k] = o[k] + amp;
    }
    let mut bbox = BoundingBox2d::from_corners(
        DVec2::new(bmin[0], bmin[1]),
        DVec2::new(bmax[0], bmax[1]),
    );
    bbox.enlarge(tol);
    bbox
}

// ── Parabola2d ────────────────────────────────────────────────────────────

fn parabola_box(par: &Parabola2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    use rcad_kernel::geom::Curve2dEval;
    let mut bbox = BoundingBox2d::new();
    if t1.is_finite() { bbox.add_point(par.point_at(t1)); }
    if t2.is_finite() { bbox.add_point(par.point_at(t2)); }
    if t1 * t2 < 0.0 { bbox.add_point(par.point_at(0.0)); }
    bbox.enlarge(tol);
    bbox
}

// ── Hyperbola2d ───────────────────────────────────────────────────────────

fn hyperbola_box(hyp: &Hyperbola2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    use rcad_kernel::geom::Curve2dEval;
    let mut bbox = BoundingBox2d::new();
    if t1.is_finite() { bbox.add_point(hyp.point_at(t1)); }
    if t2.is_finite() { bbox.add_point(hyp.point_at(t2)); }
    if t1 * t2 < 0.0 { bbox.add_point(hyp.point_at(0.0)); }

    let a_maj = hyp.semi_major;
    let a_min = hyp.semi_minor;
    let xd = hyp.major_dir;
    let yd = DVec2::new(-hyp.major_dir.y, hyp.major_dir.x);
    let eps = f64::EPSILON.sqrt();

    for k in 0..2 {
        let a = a_min * yd[k];
        let b = a_maj * xd[k];
        let abp = (a + b).abs();
        let bam = (b - a).abs();
        if abp < eps || bam < eps { continue; }
        let cf = bam / abp;
        let t3 = 0.5 * cf.ln();
        if t3 >= t1.min(t2) && t3 <= t1.max(t2) {
            bbox.add_point(hyp.point_at(t3));
        }
    }
    bbox.enlarge(tol);
    bbox
}

// ── BSplineCurve2 ────────────────────────────────────────────────────────

fn bspline_curve_box(curve: &BSplineCurve2, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let a_weakness = 1.5;
    let (sampled_box, max_deflection) = fill_bspline_span_2d(&|t| curve.point_at(t), t1, t2);
    let mut bbox = if sampled_box.is_valid() {
        let mut enlarged = sampled_box;
        enlarged.enlarge(a_weakness * max_deflection);
        reduce_curve_box_2d(&curve.control_points, &enlarged)
    } else {
        BoundingBox2d::new()
    };
    bbox.enlarge(tol);
    bbox
}

fn bspline_curve_box_optimal(curve: &BSplineCurve2, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    other_curve_box_optimal_from_fn(&|t| curve.point_at(t), t1, t2, tol)
}

// ── BezierCurve2 ─────────────────────────────────────────────────────────

fn bezier_curve_box(curve: &BezierCurve2, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let a_weakness = 1.5;
    let (sampled_box, max_deflection) = fill_bspline_span_2d(&|t| curve.point_at(t), t1, t2);
    let mut bbox = if sampled_box.is_valid() {
        let mut enlarged = sampled_box;
        enlarged.enlarge(a_weakness * max_deflection);
        reduce_curve_box_2d(&curve.control_points, &enlarged)
    } else {
        BoundingBox2d::new()
    };
    bbox.enlarge(tol);
    bbox
}

fn bezier_curve_box_optimal(curve: &BezierCurve2, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    other_curve_box_optimal_from_fn(&|t| curve.point_at(t), t1, t2, tol)
}

// ── OffsetCurve2d ───────────────────────────────────────────────────────

fn offset_curve_box(curve: &OffsetCurve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let offset = curve.offset_distance.abs();
    let basis_box = curve_box(&curve.basis, t1, t2, 0.0);
    if !basis_box.is_valid() { return basis_box; }
    let mut bbox = basis_box;
    bbox.enlarge(offset + tol);
    bbox
}

// ── OtherCurve2d (sampling fallback) ─────────────────────────────────────

fn other_curve_box(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    let weakness = 1.5;
    let n = 33;
    let (mut bbox, max_deflection) = sample_curve_box_2d(&|t| curve.point_at(t), t1, t2, n);
    if bbox.is_valid() { bbox.enlarge(weakness * max_deflection); }
    bbox.enlarge(tol);
    bbox
}

fn other_curve_box_optimal(curve: &Curve2d, t1: f64, t2: f64, tol: f64) -> BoundingBox2d {
    other_curve_box_optimal_from_fn(&|t| curve.point_at(t), t1, t2, tol)
}

fn other_curve_box_optimal_from_fn<F: Fn(f64) -> DVec2>(
    eval: &F, t1: f64, t2: f64, tol: f64,
) -> BoundingBox2d {
    let n = 65usize.max(1);
    let dt = (t2 - t1) / (n - 1) as f64;
    let dt2 = dt / 2.0;

    let mut coord_min = [f64::INFINITY; 2];
    let mut coord_max = [f64::NEG_INFINITY; 2];
    let mut defl_max = [0.0f64; 2];

    let mut pts: Vec<DVec2> = Vec::with_capacity(n);

    for i in 0..n {
        let u = t1 + i as f64 * dt;
        let p = eval(u);
        pts.push(p);
        for k in 0..2 {
            coord_min[k] = coord_min[k].min(p[k]);
            coord_max[k] = coord_max[k].max(p[k]);
        }
        if i > 0 {
            let mid = eval(u - dt2);
            let chord_center = (pts[i - 1] + p) * 0.5;
            let d = mid - chord_center;
            for k in 0..2 {
                coord_min[k] = coord_min[k].min(mid[k]);
                coord_max[k] = coord_max[k].max(mid[k]);
                defl_max[k] = defl_max[k].max(d[k].abs());
            }
        }
    }

    let eps = if tol > 1e-7f64 { tol } else { 1e-7f64 };
    for k in 0..2 {
        let d = defl_max[k];
        if d > eps {
            coord_min[k] -= d;
            coord_max[k] += d;
        }
    }

    let mut bbox = BoundingBox2d::from_corners(
        DVec2::new(coord_min[0], coord_min[1]),
        DVec2::new(coord_max[0], coord_max[1]),
    );
    bbox.enlarge(eps);
    bbox
}

// =============================================================================
// Helpers
// =============================================================================

fn fill_bspline_span_2d<F: Fn(f64) -> DVec2>(eval: &F, a: f64, b: f64) -> (BoundingBox2d, f64) {
    let mut bbox = BoundingBox2d::new();
    let n = (8usize).max((b - a).abs().ceil() as usize * 4).min(64);
    let dt = (b - a) / n as f64;
    let mut prev = eval(a);
    bbox.add_point(prev);
    let mut max_deflection: f64 = 0.0;
    for i in 1..=n {
        let t = a + dt * i as f64;
        let p = eval(t);
        bbox.add_point(p);
        let mid = eval(t - dt * 0.5);
        let chord_center = (prev + p) * 0.5;
        let deflection = (mid - chord_center).length();
        if deflection > max_deflection { max_deflection = deflection; }
        prev = p;
    }
    (bbox, max_deflection)
}

fn reduce_curve_box_2d(poles: &[DVec2], sampled: &BoundingBox2d) -> BoundingBox2d {
    if poles.is_empty() { return *sampled; }
    let mut pole_min = DVec2::splat(f64::INFINITY);
    let mut pole_max = DVec2::splat(f64::NEG_INFINITY);
    for &p in poles {
        pole_min = pole_min.min(p);
        pole_max = pole_max.max(p);
    }
    BoundingBox2d::from_min_max(
        DVec2::new(sampled.min.x.max(pole_min.x), sampled.min.y.max(pole_min.y)),
        DVec2::new(sampled.max.x.min(pole_max.x), sampled.max.y.min(pole_max.y)),
    )
}

fn sample_curve_box_2d<F: Fn(f64) -> DVec2>(eval: &F, t1: f64, t2: f64, n: usize) -> (BoundingBox2d, f64) {
    let mut bbox = BoundingBox2d::new();
    if (t2 - t1).abs() < 1e-15 { return (bbox, 0.0f64); }
    let dt = (t2 - t1) / (2 * n) as f64;
    let mut p1 = eval(t1);
    bbox.add_point(p1);
    let mut max_tol: f64 = 0.0;
    let mut p = t1;
    for _ in 1..=n {
        p += dt;
        let p2 = eval(p);
        bbox.add_point(p2);
        p += dt;
        let p3 = eval(p);
        bbox.add_point(p3);
        let chord_center = (p1 + p3) * 0.5;
        let dist = (p2 - chord_center).length();
        if dist > max_tol { max_tol = dist; }
        p1 = p3;
    }
    (bbox, max_tol)
}

// =============================================================================
// Tests — translated from OCCT GTests
// =============================================================================
//
// OCCT source:
//   src/ModelingData/TKGeomBase/GTests/
//     BndLib_Test.cxx                          (2D gp-level tests)
//     GeomBndLib_Curve2d_Test.cxx              (GeomBndLib_Curve2d tests)
//     GeomBndLib_OffsetCurve2d_Test.cxx        (2D offset curve tests)

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::*;

    const TOL: f64 = 1e-7;

    // =====================================================================
    // BndLib 2D — gp-level (from BndLib_Test.cxx)
    // =====================================================================

    #[test]
    fn line2d_finite_segment() {
        let line = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::new(1.0, 1.0).normalize() });
        let b = curve_bounds_with_range(&line, 0.0, 10.0, 0.0);
        let expected = 10.0 / 2.0f64.sqrt();
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.x - expected).abs() < TOL);
        assert!((b.max.y - expected).abs() < TOL);
    }

    #[test]
    fn circle2d_full() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let b = curve_bounds(&c, 0.0);
        assert!((b.min.x - -5.0).abs() < TOL);
        assert!((b.max.x -  5.0).abs() < TOL);
        assert!((b.min.y - -5.0).abs() < TOL);
        assert!((b.max.y -  5.0).abs() < TOL);
    }

    #[test]
    fn circle2d_arc() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let b = curve_bounds_with_range(&c, 0.0, std::f64::consts::PI / 2.0, 0.0);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 5.0).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - 5.0).abs() < TOL);
    }

    #[test]
    fn ellipse2d_full() {
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X, major_radius: 8.0, minor_radius: 4.0,
        });
        let b = curve_bounds(&e, 0.0);
        assert!((b.min.x - -8.0).abs() < TOL);
        assert!((b.max.x -  8.0).abs() < TOL);
        assert!((b.min.y - -4.0).abs() < TOL);
        assert!((b.max.y -  4.0).abs() < TOL);
    }

    #[test]
    fn hyperbola2d_simple_arc() {
        let h = Curve2d::Hyperbola(Hyperbola2d {
            center: DVec2::ZERO, major_dir: DVec2::X, semi_major: 3.0, semi_minor: 2.0,
        });
        let b = curve_bounds_with_range(&h, -1.0, 1.0, 0.0);
        let xmax = 3.0 * 1.0f64.cosh();
        let yend = 2.0 * 1.0f64.sinh();
        assert!((b.min.x - 3.0).abs() < TOL);
        assert!((b.max.x - xmax).abs() < TOL);
        assert!((b.min.y - -yend).abs() < TOL);
        assert!((b.max.y -  yend).abs() < TOL);
    }

    #[test]
    fn parabola2d_finite_arc() {
        // OCCT gp_Parab2d with focal=2, but our Parabola2d has focal_param = 2*focal.
        // For expected values matching OCCT's test (aXmin=0, aXmax=2, aYmin=-4, aYmax=4)
        // we need focal_param = 4.
        let p = Curve2d::Parabola(Parabola2d {
            origin: DVec2::ZERO, axis_dir: DVec2::X, focal_param: 4.0,
        });
        let b = curve_bounds_with_range(&p, -4.0, 4.0, 0.0);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 2.0).abs() < TOL);
        assert!((b.min.y - -4.0).abs() < TOL);
        assert!((b.max.y -  4.0).abs() < TOL);
    }

    // =====================================================================
    // BndLib_Add2dCurveTest (from BndLib_Test.cxx)
    // =====================================================================

    #[test]
    fn add2d_circle_full() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let b = curve_bounds(&c, 0.0);
        assert!((b.min.x - -5.0).abs() < TOL);
        assert!((b.max.x -  5.0).abs() < TOL);
        assert!((b.min.y - -5.0).abs() < TOL);
        assert!((b.max.y -  5.0).abs() < TOL);
    }

    #[test]
    fn add2d_circle_arc() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let b = curve_bounds_with_range(&c, 0.0, std::f64::consts::PI / 2.0, 0.0);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 5.0).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - 5.0).abs() < TOL);
    }

    #[test]
    fn add2d_ellipse_full() {
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X, major_radius: 8.0, minor_radius: 4.0,
        });
        let b = curve_bounds(&e, 0.0);
        assert!((b.min.x - -8.0).abs() < TOL);
        assert!((b.max.x -  8.0).abs() < TOL);
        assert!((b.min.y - -4.0).abs() < TOL);
        assert!((b.max.y -  4.0).abs() < TOL);
    }

    #[test]
    fn add2d_line_segment() {
        let line = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::new(2.0, 1.0).normalize() });
        // Direction (2, 1) normalized: at t = sqrt(125) ≈ 11.18, P = (10, 5)
        let len = (125.0f64).sqrt();
        let b = curve_bounds_with_range(&line, 0.0, len, 0.0);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 10.0).abs() < 0.01);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - 5.0).abs() < 0.01);
    }

    #[test]
    fn add2d_bezier_curve() {
        let bezier = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![
                DVec2::new(0.0, 0.0), DVec2::new(1.0, 3.0),
                DVec2::new(3.0, 3.0), DVec2::new(4.0, 0.0),
            ],
            weights: vec![1.0; 4],
        });
        let b = curve_bounds(&bezier, 0.0);
        assert!(b.min.x <= 0.0);
        assert!(b.max.x >= 4.0);
        assert!(b.min.y <= 0.0);
        assert!(b.max.y > 1.0);
    }

    #[test]
    fn add2d_add_optimal_ellipse() {
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X, major_radius: 6.0, minor_radius: 3.0,
        });
        let b = curve_bounds_optimal(&e, 0.0, 2.0 * std::f64::consts::PI, 0.0);
        assert!((b.min.x - -6.0).abs() < TOL);
        assert!((b.max.x -  6.0).abs() < TOL);
        assert!((b.min.y - -3.0).abs() < TOL);
        assert!((b.max.y -  3.0).abs() < TOL);
    }

    #[test]
    fn add2d_adaptor_circle() {
        let c = Curve2d::Circle(Circle2d {
            center: DVec2::new(5.0, 5.0), x_dir: DVec2::X, y_dir: DVec2::Y, radius: 3.0,
        });
        let b = curve_bounds(&c, 0.0);
        assert!((b.min.x - 2.0).abs() < TOL);
        assert!((b.max.x - 8.0).abs() < TOL);
        assert!((b.min.y - 2.0).abs() < TOL);
        assert!((b.max.y - 8.0).abs() < TOL);
    }

    // =====================================================================
    // GeomBndLib_Curve2dTest (from GeomBndLib_Curve2d_Test.cxx)
    // =====================================================================

    #[test]
    fn curve2d_line_finite_segment() {
        let line = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::new(1.0, 1.0).normalize() });
        let b = curve_bounds_with_range(&line, 0.0, 10.0, TOL);
        let len = 10.0 / 2.0f64.sqrt();
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - len).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - len).abs() < TOL);
    }

    #[test]
    fn curve2d_circle_full() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let b = curve_bounds(&c, TOL);
        assert!((b.min.x - -5.0).abs() < TOL);
        assert!((b.max.x -  5.0).abs() < TOL);
        assert!((b.min.y - -5.0).abs() < TOL);
        assert!((b.max.y -  5.0).abs() < TOL);
    }

    #[test]
    fn curve2d_circle_arc() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 10.0));
        let b = curve_bounds_with_range(&c, 0.0, std::f64::consts::PI / 2.0, TOL);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 10.0).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - 10.0).abs() < TOL);
    }

    #[test]
    fn curve2d_ellipse_full() {
        let e = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X, major_radius: 10.0, minor_radius: 5.0,
        });
        let b = curve_bounds(&e, TOL);
        assert!((b.min.x - -10.0).abs() < TOL);
        assert!((b.max.x -  10.0).abs() < TOL);
        assert!((b.min.y -  -5.0).abs() < TOL);
        assert!((b.max.y -   5.0).abs() < TOL);
    }

    #[test]
    fn curve2d_hyperbola_finite_arc() {
        let h = Curve2d::Hyperbola(Hyperbola2d {
            center: DVec2::ZERO, major_dir: DVec2::X, semi_major: 5.0, semi_minor: 3.0,
        });
        let b = curve_bounds_with_range(&h, -1.0, 1.0, TOL);
        let xend = 5.0 * 1.0f64.cosh();
        let yend = 3.0 * 1.0f64.sinh();
        assert!(b.min.x <= 5.0 + TOL);
        assert!(b.max.x >= xend - TOL);
        assert!(b.min.y <= -yend + TOL);
        assert!(b.max.y >= yend - TOL);
    }

    #[test]
    fn curve2d_parabola_finite_arc() {
        let p = Curve2d::Parabola(Parabola2d {
            origin: DVec2::ZERO, axis_dir: DVec2::X, focal_param: 4.0,
        });
        let b = curve_bounds_with_range(&p, -5.0, 5.0, TOL);
        assert!(b.min.x <= TOL);
        assert!(b.max.x > 1.0);
        assert!(b.min.y <= -4.9);
        assert!(b.max.y >= 4.9);
    }

    #[test]
    fn curve2d_bezier_simple() {
        let bezier = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(5.0, 10.0), DVec2::new(10.0, 0.0)],
            weights: vec![1.0; 3],
        });
        let b = curve_bounds(&bezier, TOL);
        assert!(b.min.x <= TOL);
        assert!(b.max.x >= 10.0 - TOL);
        assert!(b.min.y <= TOL);
        assert!(b.max.y >= 0.0);
    }

    #[test]
    fn curve2d_bspline_simple() {
        let spline = Curve2d::BSpline(BSplineCurve2 {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.0, 0.0), DVec2::new(2.0, 5.0),
                DVec2::new(8.0, 5.0), DVec2::new(10.0, 0.0),
            ],
            weights: vec![1.0; 4],
        });
        let b = curve_bounds(&spline, TOL);
        assert!(b.min.x <= TOL);
        assert!(b.max.x >= 10.0 - TOL);
        assert!(b.min.y <= TOL);
        assert!(b.max.y >= 0.0);
    }

    #[test]
    fn curve2d_add_optimal_circle() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let standard = curve_bounds(&c, TOL);
        let optimal = curve_bounds_optimal(&c, 0.0, 2.0 * std::f64::consts::PI, TOL);
        assert!(optimal.min.x >= standard.min.x - TOL);
        assert!(optimal.min.y >= standard.min.y - TOL);
        assert!(optimal.max.x <= standard.max.x + TOL);
        assert!(optimal.max.y <= standard.max.y + TOL);
    }

    #[test]
    fn curve2d_trimmed_arc() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 10.0));
        let b = curve_bounds_with_range(&c, 0.0, std::f64::consts::PI / 2.0, TOL);
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 10.0).abs() < TOL);
        assert!((b.min.y - 0.0).abs() < TOL);
        assert!((b.max.y - 10.0).abs() < TOL);
    }

    // =====================================================================
    // OffsetCurve2d tests (from GeomBndLib_OffsetCurve2d_Test.cxx)
    // =====================================================================

    #[test]
    fn offset2d_circle_positive_full() {
        let basis = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 1.0,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
        // Circle r=5, offset +1 → effective r≈6
        assert!(b.min.x <= -5.9);
        assert!(b.max.x >= 5.9);
    }

    #[test]
    fn offset2d_circle_positive_arc() {
        let basis = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 1.0,
        });
        let b = curve_bounds_with_range(&off, 0.0, std::f64::consts::PI / 2.0, 0.0);
        assert!(b.is_valid());
    }

    #[test]
    fn offset2d_circle_negative_full() {
        let basis = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 5.0));
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: -1.0,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
    }

    #[test]
    fn offset2d_ellipse_full() {
        let basis = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO, major_dir: DVec2::X, major_radius: 8.0, minor_radius: 4.0,
        });
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 1.0,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
        assert!(b.min.x <= -8.9);
        assert!(b.max.x >= 8.9);
    }

    #[test]
    fn offset2d_line() {
        let basis = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 0.5,
        });
        let b = curve_bounds_with_range(&off, 0.0, 10.0, 0.0);
        assert!(b.is_valid());
    }

    #[test]
    fn offset2d_bspline() {
        let basis = Curve2d::BSpline(BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.0, 0.0), DVec2::new(2.0, 3.0),
                DVec2::new(5.0, 3.0), DVec2::new(7.0, 0.0),
            ],
            weights: vec![1.0; 4],
        });
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 0.5,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
    }

    #[test]
    fn offset2d_large_offset() {
        let basis = Curve2d::Circle(Circle2d::new(DVec2::new(100.0, 0.0), 3.0));
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 50.0,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
        assert!(b.min.x <= 47.0);
        assert!(b.max.x >= 153.0);
    }

    #[test]
    fn offset2d_off_center_circle() {
        let basis = Curve2d::Circle(Circle2d::new(DVec2::new(10.0, 20.0), 3.0));
        let off = Curve2d::Offset(OffsetCurve2d {
            basis: Box::new(basis),
            offset_distance: 0.5,
        });
        let b = curve_bounds(&off, 0.0);
        assert!(b.is_valid());
    }

    // =====================================================================
    // Edge case tests
    // =====================================================================

    #[test]
    fn circle2d_full_period() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 3.0));
        let full = curve_bounds(&c, 0.0);
        let arc = curve_bounds_with_range(&c, 0.0, 2.0 * std::f64::consts::PI, 0.0);
        assert!((full.min.x - arc.min.x).abs() < TOL);
        assert!((full.max.x - arc.max.x).abs() < TOL);
    }

    #[test]
    fn line2d_reversed_range() {
        let line = Curve2d::Line(Line2d { origin: DVec2::ZERO, direction: DVec2::X });
        let b = curve_bounds_with_range(&line, 10.0, 0.0, 0.0);
        assert!(b.is_valid());
        assert!((b.min.x - 0.0).abs() < TOL);
        assert!((b.max.x - 10.0).abs() < TOL);
    }

    #[test]
    fn zero_radius_circle2d() {
        let c = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 0.0));
        let b = curve_bounds(&c, 0.0);
        assert!(b.is_valid());
        assert!((b.center() - DVec2::ZERO).length() < TOL);
    }
}
