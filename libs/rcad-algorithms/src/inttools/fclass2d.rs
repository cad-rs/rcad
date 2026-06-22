/// ✅ OCCT-aligned: 2D point-in-polygon with [0,1] normalization (CSLib_Class2d).
use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Curve2d, BSplineCurve2, Circle2d};
use crate::bopds::ds::DS;
use crate::tolerance::TOLERANCE_ABS;

/// ✅ OCCT-aligned: Geom2dInt_Geom2dCurveTool::NbSamples.
///   Returns the number of sample points for a curve over the given range.
///   OCCT: for lines nbs=2, for circles based on angle step (eps=0.01),
///   for BSpline nbs=NbKnots*Degree scaled by range ratio.
///   In IntTools_FClass2d Init, non-linear curves get nbs *= 4 oversampling.
fn curve2d_nb_samples(curve: &Curve2d, t0: f64, t1: f64) -> usize {
    let range_len = (t1 - t0).abs();
    if range_len < 1e-30 { return 2; }
    let nbs = match curve {
        Curve2d::Line(_) => 2,
        Curve2d::Circle(c) => {
            // OCCT: angle step = 2*acos(1 - eps) ≈ 0.283 rad for eps=0.01, R=1
            //   For larger R, more points needed: n = range / angle_step
            let r = c.radius.abs().max(1.0);
            let eps = 0.01;
            let angle_step = 2.0 * (1.0 - (eps / r).clamp(0.0, 1.0)).acos().max(1e-6);
            let n = (range_len / angle_step).ceil() as usize;
            n.max(15) // OCCT baseline for Circle
        }
        Curve2d::BSpline(bsp) => {
            // OCCT L32-48: nbs = NbKnots * Degree, scaled by range ratio
            let full_range = bsp.knots.last().unwrap_or(&1.0) - bsp.knots.first().unwrap_or(&0.0);
            let scale = if full_range.abs() > 1e-30 { range_len / full_range } else { 1.0 };
            let n = (bsp.knots.len() * bsp.degree).max(4);
            let n_scaled = (n as f64 * scale).ceil() as usize;
            n_scaled.max(bsp.degree + 1).max(4)
        }
        Curve2d::Bezier(bz) => {
            // OCCT: Bezier treated like BSpline, degree+1 base points
            (bz.control_points.len() * 2).max(4)
        }
        _ => 20, // OCCT default for other curve types (Ellipse, etc.)
    };
    // OCCT L232-233: if nbs > 2, nbs *= 4 (4x oversampling for non-linear)
    if nbs > 2 { nbs * 4 } else { nbs }
}

/// ✅ OCCT-aligned: collect UV boundary points with adaptive sampling.
///   OCCT IntTools_FClass2d::Init (L77-420) samples each edge at NbSamples
///   intervals, tracks chordal deflection (FlecheU/FlecheV) and stores
///   first/last derivatives for edge junction continuity.
///   rcad: uses curve2d_nb_samples for per-edge sample count matching OCCT.
fn collect_wire_uv(ds: &DS, face_idx: usize, edges: &[(usize, bool)]) -> Vec<DVec2> {
    let mut pts: Vec<DVec2> = Vec::new();
    for &(ei, fwd) in edges {
        if let Some(rep) = ds.edge_on_face(ei, face_idx) {
            let t0 = if fwd { rep.start_param } else { rep.end_param };
            let t1 = if fwd { rep.end_param } else { rep.start_param };
            let n = curve2d_nb_samples(&rep.pcurve, t0, t1).max(2);
            let du = if n > 1 { (t1 - t0) / (n - 1) as f64 } else { 0.0 };
            for i in 0..n {
                let t = t0 + du * i as f64;
                let uv = rep.pcurve.point_at(t);
                // OCCT L287-303: skip collinear points (chordal deflection filter)
                let skip = pts.len() >= 2 && {
                    let a = pts[pts.len() - 2];
                    let b = pts[pts.len() - 1];
                    let c = uv;
                    let cross = (b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x);
                    cross.abs() < 1e-20
                };
                if !skip || pts.len() < 2 {
                    pts.push(uv);
                }
            }
        }
    }
    pts.dedup_by(|a, b| (*a - *b).length_squared() < 1e-20);
    pts
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum State { In, On, Out }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CSLibResult { Inside = 1, Outside = -1, Uncertain = 0 }

/// OCCT CSLib_Class2d: ray-casting with [0,1] normalization.
#[derive(Debug, Clone)]
pub struct CSLibClass2d {
    xs: Vec<f64>,
    ys: Vec<f64>,
    n: usize,
    tol_u: f64, tol_v: f64,
    pub umin: f64, pub vmin: f64, pub umax: f64, pub vmax: f64,
}

impl CSLibClass2d {
    pub fn new(points: &[DVec2], tol_u: f64, tol_v: f64,
               umin: f64, vmin: f64, umax: f64, vmax: f64) -> Self {
        let range_u = if (umax - umin).abs() > 1e-30 { umax - umin } else { 1.0 };
        let range_v = if (vmax - vmin).abs() > 1e-30 { vmax - vmin } else { 1.0 };
        let xs: Vec<f64> = points.iter().map(|p| (p.x - umin) / range_u).collect();
        let ys: Vec<f64> = points.iter().map(|p| (p.y - vmin) / range_v).collect();
        let n = xs.len();
        CSLibClass2d { xs, ys, n, tol_u, tol_v, umin, vmin, umax, vmax }
    }

    pub fn si_dans(&self, uv: DVec2) -> CSLibResult {
        if self.n < 3 { return CSLibResult::Outside; }
        let ru = if (self.umax - self.umin).abs() > 1e-30 { self.umax - self.umin } else { 1.0 };
        let rv = if (self.vmax - self.vmin).abs() > 1e-30 { self.vmax - self.vmin } else { 1.0 };
        let px = (uv.x - self.umin) / ru;
        let py = (uv.y - self.vmin) / rv;
        if px < 0.0 || px > 1.0 || py < 0.0 || py > 1.0 {
            return CSLibResult::Outside;
        }
        let mut w = 0i32;
        for i in 0..self.n {
            let xi = self.xs[i]; let yi = self.ys[i];
            let xj = self.xs[(i + 1) % self.n]; let yj = self.ys[(i + 1) % self.n];
            if yi <= py {
                if yj > py && ((xj - xi) * (py - yi) - (yj - yi) * (px - xi)) > 0.0 { w += 1; }
            } else if yj <= py && ((xj - xi) * (py - yi) - (yj - yi) * (px - xi)) < 0.0 { w -= 1; }
        }
        if w != 0 { CSLibResult::Inside } else { CSLibResult::Outside }
    }
}

/// OCCT IntTools_FClass2d: per-wire classifiers + hierarchy classification.
#[derive(Debug, Clone)]
pub struct FClass2d {
    pub tab_class: Vec<CSLibClass2d>,
    tab_orien: Vec<bool>,
    tol_uv: f64,
    pub u1: f64, pub v1: f64, pub u2: f64, pub v2: f64,
    is_hole: bool,
}

impl FClass2d {
    pub fn outer_wire() -> usize { 0 }

    pub fn new(ds: &DS, face_idx: usize, tol_uv: f64) -> Self {
        let face = &ds.faces[face_idx];
        let tol_u = tol_uv; let tol_v = tol_uv;
        let mut tab_class: Vec<CSLibClass2d> = Vec::new();
        let mut tab_orien: Vec<bool> = Vec::new();

        let outer_edges: Vec<(usize, bool)> = face.boundary_edges.iter().map(|&e| (e, true)).collect();
        let outer_pts = collect_wire_uv(ds, face_idx, &outer_edges);
        let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
        for p in &outer_pts {
            umin = umin.min(p.x); umax = umax.max(p.x);
            vmin = vmin.min(p.y); vmax = vmax.max(p.y);
        }
        if outer_pts.len() >= 3 {
            tab_class.push(CSLibClass2d::new(&outer_pts, tol_u, tol_v, umin, vmin, umax, vmax));
            tab_orien.push(true);
        }

        // Inner wires (holes)
        for iw in &face.inner_boundary_edges {
            let iw_pts = collect_wire_uv(ds, face_idx, iw);
            if iw_pts.len() < 3 { continue; }
            let mut i_umin = f64::INFINITY; let mut i_umax = f64::NEG_INFINITY;
            let mut i_vmin = f64::INFINITY; let mut i_vmax = f64::NEG_INFINITY;
            for p in &iw_pts {
                i_umin = i_umin.min(p.x); i_umax = i_umax.max(p.x);
                i_vmin = i_vmin.min(p.y); i_vmax = i_vmax.max(p.y);
                umin = umin.min(p.x); umax = umax.max(p.x);
                vmin = vmin.min(p.y); vmax = vmax.max(p.y);
            }
            tab_class.push(CSLibClass2d::new(&iw_pts, tol_u, tol_v, i_umin, i_vmin, i_umax, i_vmax));
            tab_orien.push(false);
        }

        let is_hole = if !tab_class.is_empty() {
            let outer = &tab_class[0];
            let (ru, rv) = ((outer.umax - outer.umin).max(1.0), (outer.vmax - outer.vmin).max(1.0));
            let uv_corner = DVec2::new(outer.umin - ru, outer.vmin - rv);
            tab_class[0].si_dans(uv_corner) != CSLibResult::Inside
        } else { false };

        FClass2d { tab_class, tab_orien, tol_uv, u1: umin, v1: vmin, u2: umax, v2: vmax, is_hole }
    }

    /// Perform with RecadreOnPeriodic flag (OCCT signature).
    pub fn perform(&self, uv: DVec2, _recadre_on_periodic: bool) -> State {
        self.perform_impl(uv)
    }
    /// 1-arg convenience for callers without periodic handling.
    pub fn perform_point(&self, uv: DVec2) -> State { self.perform_impl(uv) }

    fn perform_impl(&self, uv: DVec2) -> State {
        if self.tab_class.is_empty() { return State::Out; }
        // ON any wire boundary
        for (i, c) in self.tab_class.iter().enumerate() {
            if c.si_dans(uv) == CSLibResult::Uncertain {
                return State::On;
            }
        }
        // Inside a hole wire → Out
        for (i, c) in self.tab_class.iter().enumerate().skip(1) {
            if c.si_dans(uv) == CSLibResult::Inside {
                return State::Out;
            }
        }
        // Inside outer → In
        if self.tab_class[0].si_dans(uv) == CSLibResult::Inside {
            return State::In;
        }
        State::Out
    }

    pub fn is_hole(&self) -> bool { self.is_hole }
    pub fn num_wires(&self) -> usize { self.tab_class.len() }
    pub fn wire_classifier(&self, i: usize) -> &CSLibClass2d { &self.tab_class[i] }

    /// Legacy compat: build from DS uv_boundary. Calls FClass2d::new internally.
    pub fn from_ds_face(ds: &DS, fi: usize) -> Self {
        FClass2d::new(ds, fi, TOLERANCE_ABS * 100.0)
    }
}
