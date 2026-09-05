// OCCT IntCurve_IntConicConic_Tool (TKGeomAlgo IntCurve) — the shared
// interval / transition helpers of the IntConicConic Perform overloads.
//
// 1:1 translation of `IntCurve_IntConicConic_Tool.hxx` (L17-151) +
// `IntCurve_IntConicConic_Tool.cxx` (L17-294): the Interval and
// PeriodicInterval interval algebras, Determine_Transition_LC and
// NormalizeOnCircleDomain.

use glam::DVec2;

use super::int_res2d::{Domain as Res2dDomain, Position, Situation, Transition, TypeTrans};

/// OCCT Tool.hxx L33: `static double PIpPI = M_PI + M_PI;`
pub(crate) const PI_PPI: f64 = std::f64::consts::TAU;
/// OCCT Tool.cxx L20: `#define TOLERANCE_ANGULAIRE 0.00000001`.
const TOLERANCE_ANGULAIRE: f64 = 0.00000001;
/// OCCT gp::Resolution() (gp.hxx).
const GP_RESOLUTION: f64 = 1.0e-15;

/// OCCT IntCurve_IntConicConic_Tool.hxx L47-64 — a bounded-or-unbounded
/// real interval with domain-bound flags.
#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub binf: f64,
    pub bsup: f64,
    pub has_first_bound: bool,
    pub has_last_bound: bool,
    pub is_null: bool,
}

impl Interval {
    /// OCCT Interval() (Tool.cxx L215-224).
    pub fn new() -> Self {
        Interval {
            binf: 0.0,
            bsup: 0.0,
            has_first_bound: false,
            has_last_bound: false,
            is_null: true,
        }
    }

    /// OCCT Interval(a, b) (Tool.cxx L226-238).
    pub fn new_bounded(a: f64, b: f64) -> Self {
        let (binf, bsup) = if a < b { (a, b) } else { (b, a) };
        Interval {
            binf,
            bsup,
            has_first_bound: true,
            has_last_bound: true,
            is_null: false,
        }
    }

    /// OCCT Interval(Domain) (Tool.cxx L240-261) — the domain bounds widened
    /// by the endpoint tolerances.
    pub fn from_domain(domain: &Res2dDomain) -> Self {
        let mut r = Interval {
            binf: 0.0,
            bsup: 0.0,
            has_first_bound: false,
            has_last_bound: false,
            is_null: false,
        };
        if domain.has_first_point() {
            r.has_first_bound = true;
            r.binf = domain.first_parameter() - domain.first_tolerance();
        }
        if domain.has_last_point() {
            r.has_last_bound = true;
            r.bsup = domain.last_parameter() + domain.last_tolerance();
        }
        r
    }

    /// OCCT Interval(a, hf, b, hl) (Tool.cxx L263-271).
    pub fn new_half_bounded(a: f64, hf: bool, b: f64, hl: bool) -> Self {
        Interval {
            binf: a,
            bsup: b,
            is_null: false,
            has_first_bound: hf,
            has_last_bound: hl,
        }
    }

    /// OCCT Length() (Tool.cxx L273-276).
    pub fn length(&self) -> f64 {
        if self.is_null {
            -1.0
        } else {
            (self.bsup - self.binf).abs()
        }
    }

    /// OCCT IntersectionWithBounded(Inter) (Tool.cxx L278-294).
    pub fn intersection_with_bounded(&self, inter: &Interval) -> Interval {
        if self.is_null || inter.is_null {
            return Interval::new();
        }
        if !(self.has_first_bound || self.has_last_bound) {
            return Interval::new_bounded(inter.binf, inter.bsup);
        }
        let a;
        if self.has_first_bound {
            if inter.bsup < self.binf {
                return Interval::new();
            }
            a = if inter.binf < self.binf {
                self.binf
            } else {
                inter.binf
            };
        } else {
            a = inter.binf;
        }

        let b;
        if self.has_last_bound {
            if inter.binf > self.bsup {
                return Interval::new();
            }
            b = if inter.bsup > self.bsup {
                self.bsup
            } else {
                inter.bsup
            };
        } else {
            b = inter.bsup;
        }
        Interval::new_bounded(a, b)
    }
}

impl Default for Interval {
    fn default() -> Self {
        Interval::new()
    }
}

/// OCCT IntCurve_PeriodicInterval (Tool.hxx L66-151) — a periodic interval
/// of the [0, 2*PI) circle parameterization.
#[derive(Debug, Clone, Copy)]
pub struct PeriodicInterval {
    pub binf: f64,
    pub bsup: f64,
    pub isnull: bool,
}

impl PeriodicInterval {
    /// OCCT PeriodicInterval() (Tool.hxx L96-100).
    pub fn new() -> Self {
        PeriodicInterval {
            isnull: true,
            binf: 0.0,
            bsup: 0.0,
        }
    }

    /// OCCT PeriodicInterval(Domain) (Tool.hxx L84-94).
    pub fn from_domain(domain: &Res2dDomain) -> Self {
        let binf = if domain.has_first_point() {
            domain.first_parameter()
        } else {
            -1.0
        };
        let bsup = if domain.has_last_point() {
            domain.last_parameter()
        } else {
            20.0
        };
        PeriodicInterval {
            isnull: false,
            binf,
            bsup,
        }
    }

    /// OCCT PeriodicInterval(a, b) (Tool.hxx L102-110).
    pub fn new_ab(a: f64, b: f64) -> Self {
        let mut r = PeriodicInterval {
            isnull: false,
            binf: a,
            bsup: b,
        };
        if (b - a) < PI_PPI {
            r.normalize();
        }
        r
    }

    /// OCCT SetNull() (Tool.hxx L72-74).
    pub fn set_null(&mut self) {
        self.isnull = true;
    }

    /// OCCT IsNull() (Tool.hxx L76-78).
    pub fn is_null(&self) -> bool {
        self.isnull
    }

    /// OCCT SetValues(a, b) (Tool.hxx L112-120).
    pub fn set_values(&mut self, a: f64, b: f64) {
        self.isnull = false;
        self.binf = a;
        self.bsup = b;
        if (b - a) < PI_PPI {
            self.normalize();
        }
    }

    /// OCCT Complement() (Tool.hxx L80-90).
    pub fn complement(&mut self) {
        if !self.isnull {
            let t = self.binf;
            self.binf = self.bsup;
            self.bsup = t + PI_PPI;
            if self.binf > PI_PPI {
                self.binf -= PI_PPI;
                self.bsup -= PI_PPI;
            }
        }
    }

    /// OCCT Length() (Tool.hxx L92-93).
    pub fn length(&self) -> f64 {
        if self.isnull {
            -100.0
        } else {
            (self.bsup - self.binf).abs()
        }
    }

    /// OCCT Normalize() (Tool.hxx L122-136).
    pub fn normalize(&mut self) {
        if !self.isnull {
            while self.binf > PI_PPI {
                self.binf -= PI_PPI;
            }
            while self.binf < 0.0 {
                self.binf += PI_PPI;
            }
            while self.bsup < self.binf {
                self.bsup += PI_PPI;
            }
            while self.bsup >= (self.binf + PI_PPI) {
                self.bsup -= PI_PPI;
            }
        }
    }

    /// OCCT FirstIntersection(PInter) (Tool.cxx L127-172).
    pub fn first_intersection(&self, p_inter: &mut PeriodicInterval) -> PeriodicInterval {
        if p_inter.isnull || self.isnull {
            return PeriodicInterval::new();
        }
        if self.length() >= PI_PPI {
            return PeriodicInterval::new_ab(p_inter.binf, p_inter.bsup);
        }
        if p_inter.length() >= PI_PPI {
            return PeriodicInterval::new_ab(self.binf, self.bsup);
        }
        if p_inter.bsup <= self.binf {
            while p_inter.binf <= self.binf && p_inter.bsup <= self.binf {
                p_inter.binf += PI_PPI;
                p_inter.bsup += PI_PPI;
            }
        }
        if p_inter.binf >= self.bsup {
            while p_inter.binf >= self.bsup && p_inter.bsup >= self.bsup {
                p_inter.binf -= PI_PPI;
                p_inter.bsup -= PI_PPI;
            }
        }
        if (p_inter.bsup < self.binf) || (p_inter.binf > self.bsup) {
            return PeriodicInterval::new();
        }

        let a = if p_inter.binf > self.binf {
            p_inter.binf
        } else {
            self.binf
        };
        let b = if p_inter.bsup < self.bsup {
            p_inter.bsup
        } else {
            self.bsup
        };

        PeriodicInterval::new_ab(a, b)
    }

    /// OCCT SecondIntersection(PInter) (Tool.cxx L177-211).
    pub fn second_intersection(&self, p_inter: &PeriodicInterval) -> PeriodicInterval {
        if p_inter.isnull
            || self.isnull
            || self.length() >= PI_PPI
            || p_inter.length() >= PI_PPI
        {
            return PeriodicInterval::new();
        }

        let mut p_inter_inf = p_inter.binf + PI_PPI;
        let mut p_inter_sup = p_inter.bsup + PI_PPI;
        if p_inter_inf > self.bsup {
            p_inter_inf = p_inter.binf - PI_PPI;
            p_inter_sup = p_inter.bsup - PI_PPI;
        }
        if (p_inter_sup < self.binf) || (p_inter_inf > self.bsup) {
            return PeriodicInterval::new();
        }
        let a = if p_inter_inf > self.binf {
            p_inter_inf
        } else {
            self.binf
        };
        let b = if p_inter_sup < self.bsup {
            p_inter_sup
        } else {
            self.bsup
        };
        PeriodicInterval::new_ab(a, b)
    }
}

impl Default for PeriodicInterval {
    fn default() -> Self {
        PeriodicInterval::new()
    }
}

/// OCCT NormalizeOnCircleDomain (Tool.cxx L104-115).
pub(crate) fn normalize_on_circle_domain(param: f64, the_domain: &Res2dDomain) -> f64 {
    let mut param = param;
    while param < the_domain.first_parameter() {
        param += PI_PPI;
    }
    while param > the_domain.last_parameter() {
        param -= PI_PPI;
    }
    param
}

/// OCCT Determine_Transition_LC (Tool.cxx L26-100) — the line/conic
/// transition decision (TOUCH classification by the signed curvature
/// comparison, IN/OUT by the cross-product sign).
#[allow(clippy::too_many_arguments)]
pub(crate) fn determine_transition_lc(
    pos1: Position,
    tan1: &mut DVec2,
    norm1: DVec2,
    t1: &mut Transition,
    pos2: Position,
    tan2: &mut DVec2,
    norm2: DVec2,
    t2: &mut Transition,
    _tol: f64,
) {
    let sgn = tan1.x * tan2.y - tan1.y * tan2.x;
    let norm = tan1.length() * tan2.length();

    if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
        // Transition TOUCH #########
        let opos = tan1.dot(*tan2) < 0.0;

        //  Modified by Sergey KHROMOV - Thu Nov  2 17:57:15 2000 Begin
        *tan1 = tan1.normalize_or_zero();
        //  Modified by Sergey KHROMOV - Thu Nov  2 17:57:16 2000 End
        let norm_v = DVec2::new(-tan1.y, tan1.x);

        let val1 = norm_v.dot(norm1);
        let val2 = norm_v.dot(norm2);

        if (val1 - val2).abs() <= GP_RESOLUTION {
            t1.set_value_touch(true, pos1, Situation::Unknown, opos);
            t2.set_value_touch(true, pos2, Situation::Unknown, opos);
        } else if val2 > val1 {
            t2.set_value_touch(true, pos2, Situation::Inside, opos);
            if opos {
                t1.set_value_touch(true, pos1, Situation::Inside, opos);
            } else {
                t1.set_value_touch(true, pos1, Situation::Outside, opos);
            }
        } else {
            // Val1 > Val2
            t2.set_value_touch(true, pos2, Situation::Outside, opos);
            if opos {
                t1.set_value_touch(true, pos1, Situation::Outside, opos);
            } else {
                t1.set_value_touch(true, pos1, Situation::Inside, opos);
            }
        }
    } else if sgn < 0.0 {
        t1.set_value_in_out(false, pos1, TypeTrans::In);
        t2.set_value_in_out(false, pos2, TypeTrans::Out);
    } else {
        t1.set_value_in_out(false, pos1, TypeTrans::Out);
        t2.set_value_in_out(false, pos2, TypeTrans::In);
    }
}
