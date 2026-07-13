/// ✅ OCCT-aligned: 2D point-in-polygon with [0,1] normalization (CSLib_Class2d).
use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Curve2d, BSplineCurve2, Circle2d};
use crate::bopds::ds::DS;
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_LEN_SQ_DIV_SAFE};
use rcad_kernel::geom::Surface3;

/// ✅ OCCT-aligned: Geom2dInt_Geom2dCurveTool::NbSamples.
///   Returns the number of sample points for a curve over the given range.
///   OCCT: for lines nbs=2, for circles based on angle step (eps=0.01),
///   for BSpline nbs=NbKnots*Degree scaled by range ratio.
///   In IntTools_FClass2d Init, non-linear curves get nbs *= 4 oversampling.
pub(crate) fn curve2d_nb_samples(curve: &Curve2d, t0: f64, t1: f64) -> usize {
    let range_len = (t1 - t0).abs();
    if range_len < TOLERANCE_LEN_SQ_DIV_SAFE { return 2; }
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
            let scale = if full_range.abs() > TOLERANCE_LEN_SQ_DIV_SAFE { range_len / full_range } else { 1.0 };
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
///   The `sample_mult` parameter scales the sample count (default 1.0).
pub(crate) fn collect_wire_uv(ds: &DS, face_idx: usize, edges: &[(usize, bool)]) -> Vec<DVec2> {
    collect_wire_uv_with_mult(ds, face_idx, edges, 1.0)
}

fn collect_wire_uv_with_mult(ds: &DS, face_idx: usize, edges: &[(usize, bool)], sample_mult: f64) -> Vec<DVec2> {
    let mut pts: Vec<DVec2> = Vec::new();
    for &(ei, fwd) in edges {
        if let Some(rep) = ds.edge_on_face(ei, face_idx) {
            let t0 = if fwd { rep.start_param } else { rep.end_param };
            let t1 = if fwd { rep.end_param } else { rep.start_param };
            let n = ((curve2d_nb_samples(&rep.pcurve, t0, t1) as f64) * sample_mult).ceil() as usize;
            let n = n.max(2);
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
    tol_u: f64,
    tol_v: f64,
    umin: f64,
    vmin: f64,
    umax: f64,
    vmax: f64,
}

impl CSLibClass2d {
    pub fn new(points: &[DVec2], tol_u: f64, tol_v: f64,
               umin: f64, vmin: f64, umax: f64, vmax: f64) -> Self {
        let range_u = if (umax - umin).abs() > TOLERANCE_LEN_SQ_DIV_SAFE { umax - umin } else { 1.0 };
        let range_v = if (vmax - vmin).abs() > TOLERANCE_LEN_SQ_DIV_SAFE { vmax - vmin } else { 1.0 };
        let xs: Vec<f64> = points.iter().map(|p| (p.x - umin) / range_u).collect();
        let ys: Vec<f64> = points.iter().map(|p| (p.y - vmin) / range_v).collect();
        let n = xs.len();
        CSLibClass2d { xs, ys, n, tol_u, tol_v, umin, vmin, umax, vmax }
    }

    pub fn si_dans(&self, uv: DVec2) -> CSLibResult {
        if self.n < 3 { return CSLibResult::Outside; }
        let ru = if (self.umax - self.umin).abs() > TOLERANCE_LEN_SQ_DIV_SAFE { self.umax - self.umin } else { 1.0 };
        let rv = if (self.vmax - self.vmin).abs() > TOLERANCE_LEN_SQ_DIV_SAFE { self.vmax - self.vmin } else { 1.0 };

        // OCCT L155-160: quick rejection outside tolerance-expanded bounding box
        if uv.x < (self.umin - self.tol_u) || uv.x > (self.umax + self.tol_u)
            || uv.y < (self.vmin - self.tol_v) || uv.y > (self.vmax + self.tol_v)
        {
            return CSLibResult::Outside;
        }

        // Transform to normalized coordinates
        let px = (uv.x - self.umin) / ru;
        let py = (uv.y - self.vmin) / rv;

        // OCCT L166-171: internalSiDansOuOn returns Uncertain when point is ON vertex or edge
        let result = self.internal_si_dans_ou_on(px, py);
        if result == CSLibResult::Uncertain {
            return CSLibResult::Uncertain;
        }

        // OCCT L173-183: 4-corner tolerance check — shift by ±tol and re-classify
        if self.tol_u > 0.0 || self.tol_v > 0.0 {
            let tol_u_norm = self.tol_u / ru;
            let tol_v_norm = self.tol_v / rv;
            let is_inside = result == CSLibResult::Inside;
            if is_inside != self.internal_si_dans(px - tol_u_norm, py - tol_v_norm)
                || is_inside != self.internal_si_dans(px + tol_u_norm, py - tol_v_norm)
                || is_inside != self.internal_si_dans(px - tol_u_norm, py + tol_v_norm)
                || is_inside != self.internal_si_dans(px + tol_u_norm, py + tol_v_norm)
            {
                return CSLibResult::Uncertain;
            }
        }

        result
    }

    /// OCCT CSLib_Class2d::internalSiDansOuOn (L279-341): ray-casting with ON detection.
    /// Checks vertex proximity and edge ON before returning.
    fn internal_si_dans_ou_on(&self, px: f64, py: f64) -> CSLibResult {
        let mut a_nb_crossings: i32 = 0;
        let mut a_prev_dx = self.xs[0] - px;
        let mut a_prev_dy = self.ys[0] - py;
        let mut a_prev_y_is_neg = a_prev_dy < 0.0;

        for a_next_idx in 1..=self.n {
            let a_prev_idx = a_next_idx - 1;
            let a_curr_dx = self.xs[a_next_idx % self.n] - px;
            let a_curr_dy = self.ys[a_next_idx % self.n] - py;

            // OCCT L296-299: vertex proximity check
            if a_curr_dx < self.tol_u && a_curr_dx > -self.tol_u
                && a_curr_dy < self.tol_v && a_curr_dy > -self.tol_v
            {
                return CSLibResult::Uncertain;
            }

            // OCCT L301-316: edge ON detection — interpolate Y at test point's X
            let a_edge_dx = self.xs[a_next_idx % self.n] - self.xs[a_prev_idx];
            if (self.xs[a_prev_idx] - px) * a_curr_dx < 0.0 && a_edge_dx.abs() > 1e-12 {
                let a_interp_y = self.ys[a_next_idx % self.n]
                    - (self.ys[a_next_idx % self.n] - self.ys[a_prev_idx]) / a_edge_dx * a_curr_dx;
                let a_delta_y = a_interp_y - py;
                if a_delta_y >= -self.tol_v && a_delta_y <= self.tol_v {
                    return CSLibResult::Uncertain;
                }
            }

            // OCCT L318-338: ray-casting crossing count
            let a_curr_y_is_neg = a_curr_dy < 0.0;
            if a_curr_y_is_neg != a_prev_y_is_neg {
                if a_prev_dx > 0.0 && a_curr_dx > 0.0 {
                    a_nb_crossings += 1;
                } else if a_prev_dx > 0.0 || a_curr_dx > 0.0 {
                    let a_x_intersect = a_prev_dx - a_prev_dy * (a_curr_dx - a_prev_dx) / (a_curr_dy - a_prev_dy);
                    if a_x_intersect > 0.0 {
                        a_nb_crossings += 1;
                    }
                }
                a_prev_y_is_neg = a_curr_y_is_neg;
            }
            a_prev_dx = a_curr_dx;
            a_prev_dy = a_curr_dy;
        }

        if (a_nb_crossings & 1) != 0 { CSLibResult::Inside } else { CSLibResult::Outside }
    }

    /// OCCT CSLib_Class2d::internalSiDans (L234-275): pure ray-casting, no ON detection.
    fn internal_si_dans(&self, px: f64, py: f64) -> bool {
        let mut a_nb_crossings: i32 = 0;
        let mut a_prev_dx = self.xs[0] - px;
        let mut a_prev_dy = self.ys[0] - py;
        let mut a_prev_y_is_neg = a_prev_dy < 0.0;

        for a_next_idx in 1..=self.n {
            let a_curr_dx = self.xs[a_next_idx % self.n] - px;
            let a_curr_dy = self.ys[a_next_idx % self.n] - py;
            let a_curr_y_is_neg = a_curr_dy < 0.0;

            if a_curr_y_is_neg != a_prev_y_is_neg {
                if a_prev_dx > 0.0 && a_curr_dx > 0.0 {
                    a_nb_crossings += 1;
                } else if a_prev_dx > 0.0 || a_curr_dx > 0.0 {
                    let a_x_intersect = a_prev_dx - a_prev_dy * (a_curr_dx - a_prev_dx) / (a_curr_dy - a_prev_dy);
                    if a_x_intersect > 0.0 {
                        a_nb_crossings += 1;
                    }
                }
                a_prev_y_is_neg = a_curr_y_is_neg;
            }
            a_prev_dx = a_curr_dx;
            a_prev_dy = a_curr_dy;
        }

        (a_nb_crossings & 1) != 0
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
    /// ✅ OCCT-aligned: periodic surface tracking for RecadreOnPeriodic
    ///   (IntTools_FClass2d.cxx L655-678, L758-802).
    ///   U/V periodicity for sphere (U=2π, V=π) and cylinder (U=2π).
    is_u_periodic: bool,
    is_v_periodic: bool,
    u_period: f64,
    v_period: f64,
    /// ✅ OCCT-aligned: UV resolution (BRepAdaptor_Surface::UResolution/VResolution).
    ///   Minimum UV step corresponding to TolUV in 3D space.
    u_res: f64,
    v_res: f64,
}

impl FClass2d {
    pub fn outer_wire() -> usize { 0 }

    pub fn new(ds: &DS, face_idx: usize, tol_uv: f64) -> Self {
        let face = &ds.faces[face_idx];
        let tol_u = tol_uv; let tol_v = tol_uv;
        let mut tab_class: Vec<CSLibClass2d> = Vec::new();
        let mut tab_orien: Vec<bool> = Vec::new();

        let outer_edges: Vec<(usize, bool)> = face.boundary_edges.iter().map(|&e| (e, true)).collect();
        let mut outer_pts = collect_wire_uv(ds, face_idx, &outer_edges);
        // ✅ OCCT-aligned: for periodic surfaces with natural restriction (sphere, cylinder, cone),
        //   the uv_boundary provides a proper closed polygon, while edge-based sampling may
        //   produce a degenerate polygon when degenerate edges lack pcurves (sphere poles).
        if outer_pts.len() < 3 && face.natural_restriction {
            if let Some(ref uv_bnd) = face.uv_boundary {
                outer_pts = uv_bnd.clone();
            }
        }
        let mut fleche_u = tol_uv;
        let mut fleche_v = tol_uv;

        // OCCT L456-565: self-intersection polygon refinement.
        //   If chordal deflection is too large relative to area/perimeter ratio,
        //   the polygon may self-intersect.  Re-discretize with tighter tolerance.
        if outer_pts.len() >= 3 {
            let mut outer_area = polygon_signed_area(&outer_pts);
            let mut outer_perim = polygon_perimeter(&outer_pts);
            let mut an_exp_thick = (2.0 * outer_area.abs() / outer_perim).max(1e-7);
            let (mut fu, mut fv) = chordal_deflection(&outer_pts, tol_uv);
            let mut a_defl = fu.max(fv);
            let mut refine_iter = 0;
            while a_defl > an_exp_thick && a_defl > 1e-7 && refine_iter < 5 {
                let new_mult = ((a_defl / an_exp_thick) * 2.0).ceil().min(32.0);
                let new_pts = collect_wire_uv_with_mult(ds, face_idx, &outer_edges, new_mult);
                if new_pts.len() <= outer_pts.len() { break; }
                outer_pts = new_pts;
                outer_area = polygon_signed_area(&outer_pts);
                outer_perim = polygon_perimeter(&outer_pts);
                an_exp_thick = (2.0 * outer_area.abs() / outer_perim).max(1e-7);
                let (new_fu, new_fv) = chordal_deflection(&outer_pts, tol_uv);
                fu = new_fu; fv = new_fv;
                a_defl = fu.max(fv);
                if a_defl <= an_exp_thick { break; }
                refine_iter += 1;
            }
            fleche_u = fu; fleche_v = fv;
        }
        let mut umin = f64::INFINITY; let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY; let mut vmax = f64::NEG_INFINITY;
        for p in &outer_pts {
            umin = umin.min(p.x); umax = umax.max(p.x);
            vmin = vmin.min(p.y); vmax = vmax.max(p.y);
        }

        if outer_pts.len() >= 3 {
            // OCCT: TabOrien derived from polygon winding (CCW=FORWARD=true, CW=REVERSED=false).
            let outer_ccw = polygon_is_ccw(&outer_pts);
            tab_class.push(CSLibClass2d::new(&outer_pts, fleche_u, fleche_v, umin, vmin, umax, vmax));
            tab_orien.push(outer_ccw);
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
            let (fleche_u, fleche_v) = chordal_deflection(&iw_pts, tol_uv);
            let iw_ccw = polygon_is_ccw(&iw_pts);
            tab_class.push(CSLibClass2d::new(&iw_pts, fleche_u, fleche_v, i_umin, i_vmin, i_umax, i_vmax));
            tab_orien.push(iw_ccw);
        }

        let is_hole = if !tab_class.is_empty() {
            let outer = &tab_class[0];
            let (ru, rv) = ((outer.umax - outer.umin).max(1.0), (outer.vmax - outer.vmin).max(1.0));
            let uv_corner = DVec2::new(outer.umin - ru, outer.vmin - rv);
            tab_class[0].si_dans(uv_corner) != CSLibResult::Inside
        } else { false };

        // OCCT L655-678: detect periodic surfaces for RecadreOnPeriodic.
        let (is_u_periodic, is_v_periodic, u_period, v_period) = match &face.surface {
            Surface3::Sphere(_) => (true, true, std::f64::consts::TAU, std::f64::consts::PI),
            Surface3::Cylinder(_) => (true, false, std::f64::consts::TAU, 0.0),
            _ => (false, false, 0.0, 0.0),
        };

        // OCCT L586-619: periodic surface U1/U2 range expansion.
        let (u1, u2, v1, v2) = match &face.surface {
            Surface3::Sphere(_) => {
                let du = std::f64::consts::TAU - (umax - umin);
                let du = if du < 0.0 { 0.0 } else { du };
                (umin - du * 0.5, umin - du * 0.5 + std::f64::consts::TAU, vmin, vmax)
            }
            Surface3::Cylinder(_) => {
                let du = std::f64::consts::TAU - (umax - umin);
                let du = if du < 0.0 { 0.0 } else { du };
                (umin - du * 0.5, umin - du * 0.5 + std::f64::consts::TAU, vmin, vmax)
            }
            _ => (umin, umax, vmin, vmax),
        };

        // ✅ OCCT-aligned: compute UV resolution (UResolution/VResolution).
        let char_len = (umax - umin).max(vmax - vmin).max(1.0);
        let u_res = tol_uv / char_len;
        let v_res = tol_uv / char_len;

        FClass2d {
            tab_class, tab_orien, tol_uv,
            u1, v1, u2, v2,
            is_hole,
            is_u_periodic, is_v_periodic, u_period, v_period,
            u_res, v_res,
        }
    }

    /// ✅ OCCT-aligned: Perform with RecadreOnPeriodic (IntTools_FClass2d.cxx L637-803).
    ///   For periodic surfaces (sphere U/V, cylinder U), when the UV point
    ///   falls outside the cached [Umin,Umax]×[Vmin,Vmax] bounds, shift the
    ///   point by the period and reclassify.  OCCT L666-678 adjusts the UV
    ///   (GeomInt::AdjustPeriodic), L680-803 retries with period-shifted coordinates.
    pub fn perform(&self, uv: DVec2, recadre_on_periodic: bool) -> State {
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return State::In;
        }

        // OCCT L645-649: save original UV
        let mut u = uv.x;
        let mut v = uv.y;
        let mut uu = u;
        let mut vv = v;

        // OCCT L666-678: adjust periodic parameters (GeomInt::AdjustPeriodic)
        if recadre_on_periodic {
            if self.is_u_periodic {
                let (new_uu, _) = self.adjust_periodic(uu, self.u1, self.u2, self.u_period);
                uu = new_uu;
            }
            if self.is_v_periodic {
                let (new_vv, _) = self.adjust_periodic(vv, self.v1, self.v2, self.v_period);
                vv = new_vv;
            }
        }

        if !recadre_on_periodic || (!self.is_u_periodic && !self.is_v_periodic) {
            return self.perform_impl(DVec2::new(u, v));
        }

        // OCCT L680-803: main retry loop
        let mut urecadre = false;
        let mut vrecadre = false;
        loop {
            let result = self.perform_impl(DVec2::new(u, v));

            if result == State::In || result == State::On {
                return result;
            }

            // OCCT L768-779: U retry (reset to adjusted uu, or add period)
            if !urecadre {
                u = uu;
                urecadre = true;
            } else if self.is_u_periodic {
                u += self.u_period;
            }

            // OCCT L781-802: if U exhausted, try V shifts
            if u > self.u2 || !self.is_u_periodic {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if self.is_v_periodic {
                    v += self.v_period;
                }

                u = uu;

                // OCCT L798-801: V exhausted → return last result
                if v > self.v2 || !self.is_v_periodic {
                    return result;
                }
            }
        }
    }

    /// OCCT: GeomInt::AdjustPeriodic (GeomInt.cxx L21-48).
    /// Shifts `par` by multiples of `period` to bring it within [par_min, par_max].
    fn adjust_periodic(&self, par: f64, par_min: f64, par_max: f64, period: f64) -> (f64, f64) {
        let mut new_par = par;
        let mut offset = 0.0;
        let b_min = par_min - par > 1e-12;
        let b_max = par - par_max > 1e-12;
        if b_min || b_max {
            let dp = if b_min { par_max - par } else { par_min - par };
            let nb_per = (dp / period).trunc();
            offset = nb_per * period;
            new_par += offset;
        }
        (new_par, offset)
    }
    /// 1-arg convenience for callers without periodic handling.
    pub fn perform_point(&self, uv: DVec2) -> State { self.perform_impl(uv) }

    fn perform_impl(&self, uv: DVec2) -> State {
        if self.tab_class.is_empty() { return State::Out; }

        // OCCT L684-714: try per-wire classification with TabOrien.
        //   For each wire: if SiDans returns INSIDE and wire is FORWARD → IN.
        //   If SiDans returns INSIDE and wire is REVERSED → OUT.
        //   If SiDans returns UNCERTAIN → need fallback (BRepClass_FClassifier).
        let mut need_classifier = false;
        let mut dedans = 1i8;

        for (i, c) in self.tab_class.iter().enumerate() {
            let cur = c.si_dans(uv);
            let forw = self.tab_orien.get(i).copied().unwrap_or(true);

            if cur == CSLibResult::Inside {
                if !forw { dedans = -1; break; }
            } else if cur == CSLibResult::Outside {
                if forw { dedans = -1; break; }
            } else {
                need_classifier = true;
                break;
            }
        }

        if !need_classifier {
            return if dedans == 1 { State::In } else { State::Out };
        }

        // OCCT L726-756: BRepClass_FClassifier fallback.
        //   When SiDans is uncertain, OCCT uses BRepClass_FClassifier
        //   with tolerance from BRepAdaptor_Surface::UResolution.
        //   SiDans returns Uncertain only when the point is ON a wire
        //   boundary (within chordal tolerance), so BRepClass_FClassifier
        //   would return TopAbs_ON.  Equivalent: return State::On.
        State::On
    }

    pub fn is_hole(&self) -> bool { self.is_hole }
    pub fn num_wires(&self) -> usize { self.tab_class.len() }
    pub fn wire_classifier(&self, i: usize) -> &CSLibClass2d { &self.tab_class[i] }

    /// Legacy compat: build from DS uv_boundary. Calls FClass2d::new internally.
    pub fn from_ds_face(ds: &DS, fi: usize) -> Self {
        FClass2d::new(ds, fi, TOLERANCE_ABS * 100.0)
    }
}

/// Signed area of a 2D polygon (shoelace).  Positive = CCW (OCCT: aS).
fn polygon_signed_area(poly: &[DVec2]) -> f64 {
    let mut area = 0.0;
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    area * 0.5
}

/// Perimeter of a 2D polygon (sum of edge lengths).
fn polygon_perimeter(poly: &[DVec2]) -> f64 {
    let n = poly.len();
    if n < 2 { return 0.0; }
    let mut perim = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        perim += (poly[j] - poly[i]).length();
    }
    perim
}

/// Compute chordal deflection (FlecheU/FlecheV) for a UV polygon.
///   Maximum distance from each point to the chord connecting its neighbors.
fn chordal_deflection(pts: &[DVec2], tol_uv: f64) -> (f64, f64) {
    if pts.len() < 3 { return (tol_uv, tol_uv); }
    let mut fu = 0.0; let mut fv = 0.0;
    for i in 1..pts.len() - 1 {
        let a = pts[i - 1]; let b = pts[i]; let c = pts[i + 1];
        let ac = c - a; let len2 = ac.dot(ac);
        if len2 < TOLERANCE_LEN_SQ_DIV_SAFE { continue; }
        let t = ((b - a).dot(ac) / len2).clamp(0.0, 1.0);
        let proj = a + ac * t;
        let du = (b.x - proj.x).abs(); let dv = (b.y - proj.y).abs();
        if du > fu { fu = du; } if dv > fv { fv = dv; }
    }
    (fu.max(tol_uv), fv.max(tol_uv))
}

/// OCCT-aligned: detect polygon winding direction (signed area).
///   Positive area = counter-clockwise (FORWARD).  Used to set TabOrien.
fn polygon_is_ccw(poly: &[DVec2]) -> bool {
    polygon_signed_area(poly) > 0.0
}
