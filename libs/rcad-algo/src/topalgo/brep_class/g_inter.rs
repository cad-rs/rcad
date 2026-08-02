// OCCT Geom2dInt_GInter + IntImpParGen_Intersector — 2D intersection of a
// line (implicit curve) with a bounded parametric curve.
//
// For the FClass2d fallback this is always LINE vs edge-pcurve, so the
// implicit curve is the line: F(P) = cross(d, P - O) = 0. Roots of F(P(t))
// along the pcurve are found by sampling + Newton refinement (OCCT uses
// math_FunctionAllRoots + IntImpParGen_Intersector). Each root becomes an
// IntRes2d_IntersectionPoint with the transitions computed by
// IntImpParGen::DetermineTransition (IntImpParGen.cxx L86-206).

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};

use crate::geomalgo::int_res2d::{
    Domain, IntersectionPoint, IntersectionSegment, Position, Situation, Transition, TypeTrans,
};

// OCCT IntImpParGen.cxx L24-25.
const TOLERANCE_ANGULAIRE: f64 = 0.00000001;
const DERIVEE_PREMIERE_NULLE: f64 = 0.000000000001;

/// OCCT IntImpParGen::NormalizeOnDomain (IntImpParGen.cxx L28-46).
fn normalize_on_domain(param: f64, domain: &Domain) -> f64 {
    let mut mod_param = param;
    if domain.is_closed() {
        let (t, period_end) = domain.equivalent_parameters();
        let periode = period_end - t;
        while mod_param < domain.first_parameter()
            && mod_param + periode < domain.last_parameter()
        {
            mod_param += periode;
        }
        while mod_param > domain.last_parameter() && mod_param - periode > domain.first_parameter()
        {
            mod_param -= periode;
        }
    }
    mod_param
}

/// OCCT IntImpParGen::DeterminePosition (IntImpParGen.cxx L49-83).
fn determine_position(domain: &Domain, pnt: DVec2, param: f64) -> Position {
    let mut pos = Position::Middle;
    if domain.has_first_point() && pnt.distance(domain.first_point()) <= domain.first_tolerance() {
        pos = Position::Head;
    }
    if domain.has_last_point() && pnt.distance(domain.last_point()) <= domain.last_tolerance() {
        if pos == Position::Head {
            if (param - domain.last_parameter()).abs() < (param - domain.first_parameter()).abs() {
                pos = Position::End;
            }
        } else {
            pos = Position::End;
        }
    }
    pos
}

/// OCCT IntImpParGen::DetermineTransition (IntImpParGen.cxx L86-206).
/// `tan1/norm1` are the tangent/normal of the first curve (the line),
/// `tan2/norm2` of the second (the pcurve).
fn determine_transition(
    pos1: Position,
    mut tan1: DVec2,
    norm1: DVec2,
    pos2: Position,
    mut tan2: DVec2,
    norm2: DVec2,
) -> (Transition, Transition) {
    let mut courbure1 = true;
    let mut courbure2 = true;
    let mut decide = true;

    let mut t1 = Transition::empty();
    let mut t2 = Transition::empty();
    t1.set_position(pos1);
    t2.set_position(pos2);

    if tan1.length_squared() <= DERIVEE_PREMIERE_NULLE {
        tan1 = norm1;
        courbure1 = false;
        if tan1.length_squared() <= DERIVEE_PREMIERE_NULLE {
            decide = false;
        }
    }
    if tan2.length_squared() <= DERIVEE_PREMIERE_NULLE {
        tan2 = norm2;
        courbure2 = false;
        if tan2.length_squared() <= DERIVEE_PREMIERE_NULLE {
            decide = false;
        }
    }

    if !decide {
        t1 = Transition::undecided(pos1);
        t2 = Transition::undecided(pos2);
    } else {
        let sgn = tan1.x * tan2.y - tan1.y * tan2.x; // Crossed
        let norm = tan1.length() * tan2.length();
        if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
            // Transition TOUCH.
            let opos = tan1.dot(tan2) < 0.0;
            if !(courbure1 || courbure2) {
                t1 = Transition::touch(true, pos1, Situation::Unknown, opos);
                t2 = Transition::touch(true, pos2, Situation::Unknown, opos);
            } else {
                // Norm = (-Tan1.Y, Tan1.X) — rotate Tan1 by -90°.
                let nrm = DVec2::new(-tan1.y, tan1.x);
                let val1 = if !courbure1 {
                    0.0
                } else {
                    nrm.dot(norm1)
                };
                let val2 = if !courbure2 {
                    0.0
                } else {
                    nrm.dot(norm2)
                };
                if (val1 - val2).abs() <= TOLERANCE_ANGULAIRE {
                    t1 = Transition::touch(true, pos1, Situation::Unknown, opos);
                    t2 = Transition::touch(true, pos2, Situation::Unknown, opos);
                } else if val2 > val1 {
                    t2 = Transition::touch(true, pos2, Situation::Inside, opos);
                    let s1 = if opos { Situation::Inside } else { Situation::Outside };
                    t1 = Transition::touch(true, pos1, s1, opos);
                } else {
                    // Val1 > Val2.
                    t2 = Transition::touch(true, pos2, Situation::Outside, opos);
                    let s1 = if opos { Situation::Outside } else { Situation::Inside };
                    t1 = Transition::touch(true, pos1, s1, opos);
                }
            }
        } else if sgn < 0.0 {
            t1 = Transition::in_out(false, pos1, TypeTrans::In);
            t2 = Transition::in_out(false, pos2, TypeTrans::Out);
        } else {
            // sgn > 0.
            t1 = Transition::in_out(false, pos1, TypeTrans::Out);
            t2 = Transition::in_out(false, pos2, TypeTrans::In);
        }
    }
    (t1, t2)
}

/// OCCT Geom2dInt_GInter — result of intersecting a line with a pcurve.
pub struct GInter {
    points: Vec<IntersectionPoint>,
    segments: Vec<IntersectionSegment>,
    done: bool,
}

impl GInter {
    /// OCCT Geom2dInt_GInter(C1, D1, C2, D2, TolConf, Tol) with C1 the line
    /// (implicit) and C2 the pcurve.
    pub fn new(
        line_origin: DVec2,
        line_dir: DVec2,
        line_domain: &Domain,
        curve: &Curve2d,
        curve_domain: &Domain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut g = GInter {
            points: Vec::new(),
            segments: Vec::new(),
            done: false,
        };
        g.perform(line_origin, line_dir, line_domain, curve, curve_domain, tol_conf, tol);
        g
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty() && self.segments.is_empty()
    }

    pub fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// 1-indexed.
    pub fn point(&self, i: usize) -> &IntersectionPoint {
        &self.points[i - 1]
    }

    pub fn nb_segments(&self) -> usize {
        self.segments.len()
    }

    /// 1-indexed.
    pub fn segment(&self, i: usize) -> &IntersectionSegment {
        &self.segments[i - 1]
    }

    /// OCCT IntImpParGen_Intersector::Perform (IntImpParGen_Intersector.gxx
    /// L245-779), focused on the line-vs-curve point intersections.
    pub fn perform(
        &mut self,
        line_origin: DVec2,
        line_dir: DVec2,
        line_domain: &Domain,
        curve: &Curve2d,
        curve_domain: &Domain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.points.clear();
        self.segments.clear();
        self.done = false;

        let eps_x = 1.0e-10_f64.min(1.0e-10); // EPSX (curve-tool EpsX capped)
        let eps_nul = if tol_conf <= 1.0e-10 { 1.0e-10 } else { tol_conf };
        let eps_dist = if tol <= 1.0e-10 { 1.0e-10 } else { tol };

        if !(curve_domain.has_first_point() && curve_domain.has_last_point()) {
            self.done = false;
            return;
        }

        let (t_min, t_max) = (curve_domain.first_parameter(), curve_domain.last_parameter());
        if !t_min.is_finite() || !t_max.is_finite() {
            self.done = true;
            return;
        }

        // Implicit line function F(t) = cross(d, P(t) - O).
        let f_val = |t: f64| -> f64 {
            let p = curve.point_at(t);
            (p - line_origin).x * line_dir.y - (p - line_origin).y * line_dir.x
        };

        // Bracket roots by sampling on a grid.
        let nb_samples = 64usize;
        let mut roots: Vec<f64> = Vec::new();
        let mut prev_t = t_min;
        let mut prev_f = f_val(t_min);
        for i in 1..=nb_samples {
            let t = t_min + (t_max - t_min) * (i as f64) / (nb_samples as f64);
            let f = f_val(t);
            if prev_f == 0.0 {
                roots.push(prev_t);
            } else if f != 0.0 && prev_f.signum() != f.signum() {
                // Sign change in (prev_t, t): bisection + Newton refine.
                let root = self.refine_root(curve, line_origin, line_dir, prev_t, t, prev_f, f, eps_dist);
                if root.is_finite() {
                    roots.push(root);
                }
            }
            prev_t = t;
            prev_f = f;
        }
        if prev_f == 0.0 {
            roots.push(prev_t);
        }

        // Deduplicate roots closer than eps_x.
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots.dedup_by(|a, b| (*a - *b).abs() < eps_x);

        // Build an intersection point per root.
        for &t in &roots {
            let pt = curve.point_at(t);
            // Param on the line: (P - O)·d (d is unit).
            let param1 = (pt - line_origin).dot(line_dir);
            let param1 = normalize_on_domain(param1, line_domain);
            let param2 = t;

            let pos1 = determine_position(line_domain, pt, param1);
            let pos2 = determine_position(curve_domain, pt, param2);

            let tan1 = line_dir;
            let norm1 = DVec2::ZERO; // a line has zero second derivative
            let tan2 = curve.derivative_at(t);
            let norm2 = curve.derivative2_at(t);

            let (trans1, trans2) = determine_transition(pos1, tan1, norm1, pos2, tan2, norm2);

            let ip = IntersectionPoint::new(pt, param1, param2, trans1, trans2, false);
            self.insert(ip);
        }

        self.done = true;
    }

    /// Newton refinement of a sign-change bracket (prev_f * f < 0).
    fn refine_root(
        &self,
        curve: &Curve2d,
        line_origin: DVec2,
        line_dir: DVec2,
        t0: f64,
        t1: f64,
        mut f0: f64,
        mut f1: f64,
        _eps: f64,
    ) -> f64 {
        let mut lo = t0;
        let mut hi = t1;
        let mut t = 0.5 * (lo + hi);
        for _ in 0..30 {
            t = 0.5 * (lo + hi);
            let f = (curve.point_at(t) - line_origin).x * line_dir.y
                - (curve.point_at(t) - line_origin).y * line_dir.x;
            if f == 0.0 || (hi - lo) < 1e-13 {
                break;
            }
            // Newton step using dF/dt = cross(d, P'(t)).
            let dp = curve.derivative_at(t);
            let df = line_dir.x * dp.y - line_dir.y * dp.x;
            if df.abs() > 1e-30 {
                let tn = t - f / df;
                if tn > lo && tn < hi {
                    t = tn;
                    let fv = (curve.point_at(t) - line_origin).x * line_dir.y
                        - (curve.point_at(t) - line_origin).y * line_dir.x;
                    if fv.abs() < f.abs() && fv == 0.0 {
                        break;
                    }
                    // Re-bracket.
                    if f0.signum() != fv.signum() {
                        hi = t;
                        f1 = fv;
                    } else if f1.signum() != fv.signum() {
                        lo = t;
                        f0 = fv;
                    }
                    continue;
                }
            }
            // Bisection fallback.
            if f.signum() == f0.signum() {
                lo = t;
                f0 = f;
            } else {
                hi = t;
                f1 = f;
            }
        }
        t
    }

    /// OCCT IntRes2d_Intersection::Insert — insert keeping points sorted by
    /// ParamOnFirst, skipping duplicates.
    fn insert(&mut self, pnt: IntersectionPoint) {
        let n = self.points.len();
        if n == 0 {
            self.points.push(pnt);
            return;
        }
        let u = pnt.param_on_first();
        let mut b = n + 1;
        for (i, pnti) in self.points.iter().enumerate() {
            let ui = pnti.param_on_first();
            if ui >= u {
                b = i + 1;
                break;
            }
            if (ui - u).abs() < 1e-8 {
                if (pnt.param_on_second() - pnti.param_on_second()).abs() < 1e-8
                    && transition_equal(pnt.transition_of_first(), pnti.transition_of_first())
                    && transition_equal(pnt.transition_of_second(), pnti.transition_of_second())
                {
                    b = 0;
                    break;
                }
            }
        }
        if b > n {
            self.points.push(pnt);
        } else if b > 0 {
            self.points.insert(b - 1, pnt);
        }
    }
}

/// OCCT IntRes2d_Intersection.cxx L38-65: TransitionEqual.
fn transition_equal(t1: &Transition, t2: &Transition) -> bool {
    if t1.position_on_curve() == t2.position_on_curve()
        && t1.transition_type() == t2.transition_type()
    {
        if t1.transition_type() == TypeTrans::Touch {
            if t1.is_tangent() == t2.is_tangent()
                && t1.situation() == t2.situation()
                && t1.is_opposite() == t2.is_opposite()
            {
                return true;
            }
        } else {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle2d, Line2d};

    fn line_domain() -> Domain {
        // Semi-infinite line domain from the origin along +x: [0, inf].
        Domain::semi(DVec2::ZERO, 0.0, 1e-9, true)
    }

    fn circle_domain() -> Domain {
        let c = Circle2d::new(DVec2::new(2.0, 0.0), 1.0);
        let frame_point = |ang: f64| {
            c.center + c.x_dir * (c.radius * ang.cos()) + c.y_dir * (c.radius * ang.sin())
        };
        let p0 = frame_point(0.0);
        let p1 = frame_point(std::f64::consts::TAU);
        let mut d = Domain::bounded(p0, 0.0, 1e-9, p1, std::f64::consts::TAU, 1e-9);
        d.set_equivalent_parameters(0.0, std::f64::consts::TAU);
        d
    }

    #[test]
    fn line_x_circle_enter_exit() {
        // Line = x-axis (origin (0,0), dir (1,0)); circle center (2,0), r=1.
        // Crossings at line-param 1 (enter) and 3 (exit).
        let line = Line2d::new(DVec2::ZERO, DVec2::X);
        let circle = Curve2d::Circle(Circle2d::new(DVec2::new(2.0, 0.0), 1.0));
        let dl = line_domain();
        let de = circle_domain();
        let g = GInter::new(
            line.origin,
            line.direction,
            &dl,
            &circle,
            &de,
            1e-9,
            1e-9,
        );
        assert!(g.is_done(), "intersector should complete");
        assert_eq!(g.nb_points(), 2, "two crossings expected");
        // Sorted by ParamOnFirst.
        let p1 = g.point(1);
        let p2 = g.point(2);
        assert!((p1.param_on_first() - 1.0).abs() < 1e-6, "enter at line-param 1, got {}", p1.param_on_first());
        assert!((p2.param_on_first() - 3.0).abs() < 1e-6, "exit at line-param 3, got {}", p2.param_on_first());
        // First crossing: line enters the circle (transition In).
        assert_eq!(p1.transition_of_first().transition_type(), TypeTrans::In);
        // Second crossing: line exits (transition Out).
        assert_eq!(p2.transition_of_first().transition_type(), TypeTrans::Out);
        // The circle's transition is the opposite.
        assert_eq!(p1.transition_of_second().transition_type(), TypeTrans::Out);
        assert_eq!(p2.transition_of_second().transition_type(), TypeTrans::In);
    }

    #[test]
    fn line_x_circle_tangent() {
        // Line y=1 (origin (0,1), dir (1,0)); circle center (1,0), r=1. The
        // line touches the circle at exactly (1,1) (circle param π/2, a
        // sample point). The transition is Touch (the curves are parallel).
        let line = Line2d::new(DVec2::new(0.0, 1.0), DVec2::X);
        let circle = Curve2d::Circle(Circle2d::new(DVec2::new(1.0, 0.0), 1.0));
        let dl = Domain::semi(line.origin, 0.0, 1e-9, true);
        let de = circle_domain();
        let g = GInter::new(line.origin, line.direction, &dl, &circle, &de, 1e-9, 1e-9);
        assert!(g.is_done());
        // A tangency at a sampled zero is detected; the transition must be
        // Touch (or Undecided), never In/Out.
        for i in 1..=g.nb_points() {
            let t = g.point(i).transition_of_first().transition_type();
            assert!(
                t == TypeTrans::Touch || t == TypeTrans::Undecided,
                "tangent crossing should be Touch/Undecided, got {:?}",
                t
            );
        }
    }

    #[test]
    fn line_x_vertical_line() {
        // Line A = x-axis; line B = vertical x=1. Intersection at (1,0).
        let line = Line2d::new(DVec2::ZERO, DVec2::X);
        let vertical = Curve2d::Line(Line2d::new(DVec2::new(1.0, -5.0), DVec2::Y));
        let dl = line_domain();
        let dv = Domain::bounded(DVec2::new(1.0, -5.0), -5.0, 1e-9, DVec2::new(1.0, 5.0), 5.0, 1e-9);
        let g = GInter::new(line.origin, line.direction, &dl, &vertical, &dv, 1e-9, 1e-9);
        assert!(g.is_done());
        assert_eq!(g.nb_points(), 1);
        let p = g.point(1);
        assert!((p.param_on_first() - 1.0).abs() < 1e-6);
    }
}

