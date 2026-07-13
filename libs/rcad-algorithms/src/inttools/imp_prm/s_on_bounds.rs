//! ✅ OCCT-aligned: IntPatch_TheSOnBounds — boundary scanning for F(t)=0 on domain edges.
//!
//! OCCT IntStart_SearchOnBoundaries.gxx (1232 lines) — generic template instantiated
//! as IntPatch_TheSOnBounds (TheArc=Handle(Adaptor2d_Curve2d), TheFunction=ArcFunction).
//!
//! Algorithm (BoundedArc):
//!   1. Sample F(t)=Q(C(t)) along each boundary curve
//!   2. Use root-finding (OCCT: math_FunctionAllRoots / IntCurveSurface_HInter)
//!      rcad: sign-change detection with sample refinement
//!   3. Collect isolated points (NbPoints) and continuous segments (NbSegments)
//!   4. Special treatment for linear curves against quadric surfaces (TreatLC)

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use crate::tolerance::TOLERANCE_CLAMP_MIN;
use super::arc_function::ArcFunction;
use super::super::int_surf_quadric::Quadric;
use super::super::geom_abs_surface_type::GeomAbsSurfaceType;

// ── OCCT IntPatch_ThePathPointOfTheSOnBounds ──────────────────────────
#[derive(Clone, Debug)]
pub struct PathPoint {
    pub value: glam::DVec3,
    pub tolerance: f64,
    pub parameter: f64,      // curve parameter on boundary arc
    pub u: f64,              // U coordinate on the parametric surface at this point
    pub v: f64,              // V coordinate on the parametric surface at this point
    pub arc_index: usize,
    pub is_new: bool,
}

impl PathPoint {
    pub fn new(value: glam::DVec3, tol: f64, param: f64, u: f64, v: f64, arc_i: usize, is_new: bool) -> Self {
        Self { value, tolerance: tol, parameter: param, u, v, arc_index: arc_i, is_new }
    }
}

// ── OCCT IntPatch_TheSegmentOfTheSOnBounds ───────────────────────────
#[derive(Clone, Debug)]
pub struct Segment {
    pub curve: Curve2d,
    pub first_point_index: usize,
    pub last_point_index: usize,
}

impl Segment {
    pub fn new(curve: Curve2d) -> Self {
        Self { curve, first_point_index: 0, last_point_index: 0 }
    }
    pub fn has_first_point(&self) -> bool { self.first_point_index > 0 }
    pub fn has_last_point(&self) -> bool { self.last_point_index > 0 }
}

// ── OCCT IntPatch_TheSOnBounds ───────────────────────────────────────
pub struct SOnBounds {
    done: bool,
    all: bool,
    points: Vec<PathPoint>,
    segments: Vec<Segment>,
}

impl SOnBounds {
    // OCCT L44-46: constructor
    pub fn new() -> Self {
        Self { done: false, all: false, points: Vec::new(), segments: Vec::new() }
    }

    // OCCT L52-56: Perform(ArcFunction, Domain, TolBoundary, TolTangency, RecheckOnRegularity)
    // rcad: Domain = UV rectangle [u_min,u_max]x[v_min,v_max]
    pub fn perform(
        &mut self,
        func: &mut ArcFunction,
        u_min: f64, u_max: f64,
        v_min: f64, v_max: f64,
        tol_boundary: f64,
        tol_tangency: f64,
    ) {
        self.points.clear();
        self.segments.clear();
        self.done = false;
        self.all = false;

        // OCCT L52-56: four boundary arcs of the rectangular domain
        // Arc 0: V=const at Vmin (U varies)
        // Arc 1: V=const at Vmax (U varies)
        // Arc 2: U=const at Umin (V varies)
        // Arc 3: U=const at Umax (V varies)
        let arcs = [
            (Curve2d::Line(rcad_kernel::geom::Line2d { origin: DVec2::new(u_min, v_min), direction: DVec2::new(u_max - u_min, 0.0) }),
             u_min, u_max, 0usize),
            (Curve2d::Line(rcad_kernel::geom::Line2d { origin: DVec2::new(u_min, v_max), direction: DVec2::new(u_max - u_min, 0.0) }),
             u_min, u_max, 1usize),
            (Curve2d::Line(rcad_kernel::geom::Line2d { origin: DVec2::new(u_min, v_min), direction: DVec2::new(0.0, v_max - v_min) }),
             v_min, v_max, 2usize),
            (Curve2d::Line(rcad_kernel::geom::Line2d { origin: DVec2::new(u_max, v_min), direction: DVec2::new(0.0, v_max - v_min) }),
             v_min, v_max, 3usize),
        ];

        for (arc, p_deb, p_fin, arc_idx) in arcs {
            self.bounded_arc(func, &arc, p_deb, p_fin, arc_idx, tol_boundary, tol_tangency);
        }

        self.done = true;
        self.all = self.points.is_empty() && self.segments.is_empty();
    }

    // ── OCCT BoundedArc (L229-600+) ──────────────────────────────────
    // Analyze one boundary arc for F(t)=0 solutions.
    fn bounded_arc(
        &mut self,
        func: &mut ArcFunction,
        arc: &Curve2d,
        p_deb: f64,
        p_fin: f64,
        arc_idx: usize,
        tol_boundary: f64,
        tol_tangency: f64,
    ) {
        func.set_arc(arc.clone());
        let n_echant = func.nb_samples().max(100) as usize;

        // OCCT L273-281: adjust tolerances for short arcs
        let mut n_tol_tang = tol_tangency;
        if (p_fin - p_deb) < (tol_tangency * 10.0) {
            n_tol_tang = (p_fin - p_deb) * 0.1;
        }

        // ── Rejection test (OCCT L292-333) ──────────────────────────
        // Check if F(t) could cross zero on this arc
        let rejection = self.rejection_test(func, p_deb, p_fin, n_echant);

        if rejection {
            return; // No intersection on this arc
        }

        // ── Exact intersection for analytic quadric surfaces (OCCT L339-472) ──
        // rcad: skip exact intersection for now, use numerical root finding only
        // (full IntCurveSurface_HInter equivalent not yet available)

        // ── Numerical root finding (OCCT L475-489) ──────────────────
        // rcad: sign-change detection across samples (replaces math_FunctionAllRoots)
        let mut num_params: Vec<f64> = Vec::new();
        self.numerical_roots(func, p_deb, p_fin, n_echant, n_tol_tang, &mut num_params);

        // ── Process solutions (OCCT L568-600+) ──────────────────────
        let params = &num_params;

        // TreatLC: special handling for linear curves (OCCT L591) — skipped for now

        // Add solution points
        for &param in params {
            if let Some(f) = func.value(param) {
                if f.abs() <= n_tol_tang {
                    let uv = arc.point_at(param);
                    self.points.push(PathPoint::new(
                        *func.last_computed_point(),
                        n_tol_tang.max(tol_boundary),
                        param, uv.x, uv.y, arc_idx, true,
                    ));
                }
            }
        }
    }

    // ── OCCT Rejection test (L292-333) ──────────────────────────────
    fn rejection_test(&self, func: &mut ArcFunction, p_deb: f64, p_fin: f64, n_echant: usize) -> bool {
        let dur = (p_fin - p_deb) / 6.0f64.max(1.0);
        let mut minr = f64::MAX;
        let mut maxr = f64::MIN;
        let mut maxdr = f64::MIN;

        for i in 0..6 {
            let ur = p_deb + i as f64 * dur;
            if let Some((f, d)) = func.values(ur) {
                let d = d.abs();
                let dd = d * (dur * 2.0);
                if dd > maxdr { maxdr = dd; }
                let lminr = f - dd;
                let lmaxr = f + dd;
                if lminr < minr { minr = lminr; }
                if lmaxr > maxr { maxr = lmaxr; }
                if minr < 0.0 && maxr > 0.0 { return false; }
            }
        }

        // Soften bounds (L326-332)
        let soften = 0.001 + maxdr + (maxr - minr) * 0.1;
        minr -= soften;
        maxr += soften;
        !(minr < 0.0 && maxr > 0.0)
    }

    // ── OCCT math_FunctionAllRoots — robust 1D root finding (L475-489) ──
    //
    // Algorithm (matching OCCT IntStart_SearchOnBoundaries):
    //   1. Sample F(t) at NbEchant uniform points across [Pdeb, Pfin]
    //   2. For each interval [t_i, t_{i+1}]:
    //      a. Sign change → isolate root by linear interpolation
    //      b. Both near zero → segment (continuous region)
    //      c. One endpoint near zero → isolated point (tangent)
    //   3. Build both point list (spnt) and segment list (sseg)
    fn numerical_roots(
        &mut self,
        func: &mut ArcFunction,
        p_deb: f64,
        p_fin: f64,
        n_echant: usize,
        tol_tang: f64,
        params: &mut Vec<f64>,
    ) {
        let range = p_fin - p_deb;
        if range.abs() < TOLERANCE_CLAMP_MIN { return; }

        // OCCT: math_FunctionSample(Pdeb, Pfin, NbEchant)
        let dt = range / n_echant as f64;
        let mut prev_f: Option<f64> = None;

        // Collect sample values
        struct Sample { t: f64, f: f64 }
        let mut samples: Vec<Sample> = Vec::with_capacity(n_echant + 1);

        for i in 0..=n_echant {
            let t = p_deb + i as f64 * dt;
            if let Some(f) = func.value(t) {
                samples.push(Sample { t, f });
            }
        }

        if samples.is_empty() { return; }

        // OCCT L568-600: detect intervals and build points/segments
        let mut in_segment = false;
        let mut seg_start_t = 0.0;
        let mut seg_start_idx = 0usize;

        for i in 0..samples.len() - 1 {
            let s0 = &samples[i];
            let s1 = &samples[i + 1];

            let near0 = s0.f.abs() <= tol_tang;
            let near1 = s1.f.abs() <= tol_tang;
            let sign_change = s0.f * s1.f < 0.0;

            if near0 && near1 {
                // Both near zero → segment
                if !in_segment {
                    in_segment = true;
                    seg_start_t = s0.t;
                    seg_start_idx = params.len();
                    params.push(s0.t);
                }
                // Don't add s1 yet — wait for segment end
            } else {
                if in_segment {
                    // End of segment
                    in_segment = false;
                    // Create segment
                    let seg = Segment {
                        curve: Curve2d::Line(rcad_kernel::geom::Line2d {
                            origin: DVec2::new(p_deb, 0.0),
                            direction: DVec2::new(range, 0.0),
                        }),
                        first_point_index: seg_start_idx + 1,
                        last_point_index: params.len() + 1,
                    };
                    params.push(s0.t); // segment endpoint
                    self.segments.push(seg);
                }

                if sign_change || near0 || near1 {
                    // OCCT: root isolation via linear interpolation
                    let t_root = if near0 { s0.t }
                    else if near1 { s1.t }
                    else {
                        let alpha = s0.f.abs() / (s0.f.abs() + s1.f.abs());
                        s0.t + alpha * (s1.t - s0.t)
                    };
                    params.push(t_root);
                }
            }
        }

        // Close trailing segment
        if in_segment {
            let last_t = samples.last().unwrap().t;
            let seg = Segment {
                curve: Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: DVec2::new(p_deb, 0.0),
                    direction: DVec2::new(range, 0.0),
                }),
                first_point_index: seg_start_idx + 1,
                last_point_index: params.len() + 1,
            };
            params.push(last_t);
            self.segments.push(seg);
        }
    }

    // ── Public API ───────────────────────────────────────────────────
    pub fn is_done(&self) -> bool { self.done }
    pub fn all_arc_solution(&self) -> bool { self.all }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point(&self, index: usize) -> &PathPoint { &self.points[index] }
    pub fn nb_segments(&self) -> usize { self.segments.len() }
    pub fn segment(&self, index: usize) -> &Segment { &self.segments[index] }
}
