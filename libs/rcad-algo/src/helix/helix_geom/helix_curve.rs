//! OCCT HelixGeom_HelixCurve (TKHelix HelixGeom package).
//!
//! 1:1 translation of `HelixGeom_HelixCurve.hxx` + `.cxx` — adaptor for
//! analytic helix curves (cylindrical and tapered, CW/CCW).
//!
//! OCCT inherits from Adaptor3d_Curve; the adaptor surface here is the
//! inherent impl (the approximator consumes `eval_d0` / `eval_d1` / `eval_d2`
//! through the `ToolsEval` wrapper in `tools.rs`, standing in for
//! `HelixGeom_Tools_Eval`).

use glam::DVec3;
use rcad_kernel::math::GeomAbsShape;

/// OCCT HelixGeom_HelixCurve.
#[derive(Debug, Clone)]
pub struct HelixCurve {
    first: f64,
    last: f64,
    pitch: f64,
    r_start: f64,
    taper_angle: f64,
    is_clock_wise: bool,
    c1: f64,
    tg_beta: f64,
    tol_angle: f64,
}

impl Default for HelixCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl HelixCurve {
    /// OCCT HelixGeom_HelixCurve() — default constructor
    /// (HelixGeom_HelixCurve.cxx L30-42).
    pub fn new() -> Self {
        let last = 2.0 * std::f64::consts::PI;
        let pitch = 1.0;
        let c1 = pitch / last;
        HelixCurve {
            first: 0.0,
            last,
            pitch,
            r_start: 1.0,
            taper_angle: 0.0,
            is_clock_wise: true,
            c1,
            tg_beta: 0.0,
            tol_angle: 1.0e-4,
        }
    }

    /// OCCT HelixGeom_HelixCurve::Load() — sets default values for parameters
    /// (HelixGeom_HelixCurve.cxx L47-50).
    pub fn load_default(&mut self) {
        let (first, last, pitch, r_start, taper_angle, is_clock_wise) = (
            self.first,
            self.last,
            self.pitch,
            self.r_start,
            self.taper_angle,
            self.is_clock_wise,
        );
        self.load(first, last, pitch, r_start, taper_angle, is_clock_wise);
    }

    /// OCCT HelixGeom_HelixCurve::Load(aT1, aT2, aPitch, aRStart,
    /// aTaperAngle, aIsCW) — sets helix parameters
    /// (HelixGeom_HelixCurve.cxx L54-97).
    pub fn load(
        &mut self,
        a_t1: f64,
        a_t2: f64,
        a_pitch: f64,
        a_r_start: f64,
        a_taper_angle: f64,
        a_is_cw: bool,
    ) {
        let a_two_pi = 2.0 * std::f64::consts::PI;
        let a_half_pi = 0.5 * std::f64::consts::PI;
        // Store parameter values.
        self.first = a_t1;
        self.last = a_t2;
        self.pitch = a_pitch;
        self.r_start = a_r_start;
        self.taper_angle = a_taper_angle;
        self.is_clock_wise = a_is_cw;
        // Validate input parameters.
        if a_t1 >= a_t2 {
            panic!("HelixGeom_HelixCurve::Load");
        }
        if self.pitch < 0.0 {
            panic!("HelixGeom_HelixCurve::Load");
        }
        if self.r_start < 0.0 {
            panic!("HelixGeom_HelixCurve::Load");
        }
        if self.taper_angle <= -a_half_pi || self.taper_angle >= a_half_pi {
            panic!("HelixGeom_HelixCurve::Load");
        }
        // Calculate helix coefficient.
        self.c1 = self.pitch / a_two_pi;
        if self.taper_angle.abs() > self.tol_angle {
            self.tg_beta = self.taper_angle.tan();
        }
    }

    /// OCCT FirstParameter (HelixGeom_HelixCurve.cxx L100-102).
    pub fn first_parameter(&self) -> f64 {
        self.first
    }

    /// OCCT LastParameter (HelixGeom_HelixCurve.cxx L107-109).
    pub fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT Continuity (HelixGeom_HelixCurve.cxx L114-117).
    pub fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::CN
    }

    /// OCCT NbIntervals (HelixGeom_HelixCurve.cxx L121-124).
    pub fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Intervals (HelixGeom_HelixCurve.cxx L128-131).
    pub fn intervals(&self, t: &mut [f64], _s: GeomAbsShape) {
        t[0] = self.first;
        t[1] = self.last;
    }

    /// OCCT Resolution (HelixGeom_HelixCurve.cxx L136-139) — throws
    /// Standard_NotImplemented.
    pub fn resolution(&self, _r3d: f64) -> f64 {
        panic!("HelixGeom_HelixCurve::Resolution");
    }

    /// OCCT IsClosed (HelixGeom_HelixCurve.cxx L143-146) — throws
    /// Standard_NotImplemented.
    pub fn is_closed(&self) -> bool {
        panic!("HelixGeom_HelixCurve::IsClosed");
    }

    /// OCCT IsPeriodic (HelixGeom_HelixCurve.cxx L150-153) — throws
    /// Standard_NotImplemented.
    pub fn is_periodic(&self) -> bool {
        panic!("HelixGeom_HelixCurve::IsPeriodic");
    }

    /// OCCT Period (HelixGeom_HelixCurve.cxx L157-160) — throws
    /// Standard_DomainError.
    pub fn period(&self) -> f64 {
        panic!("HelixGeom_HelixCurve::Period");
    }

    /// OCCT EvalD0 (HelixGeom_HelixCurve.cxx L164-180).
    pub fn eval_d0(&self, the_t: f64) -> DVec3 {
        // Calculate trigonometric values and radius.
        let a_ct = the_t.cos();
        let a_st = the_t.sin();
        let a1 = self.r_start + self.c1 * self.tg_beta * the_t;
        // Calculate Cartesian coordinates.
        let a_x = a1 * a_ct;
        let mut a_y = a1 * a_st;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = self.c1 * the_t;
        DVec3::new(a_x, a_y, a_z)
    }

    /// OCCT EvalD1 (HelixGeom_HelixCurve.cxx L184-215) — point and first
    /// derivative.
    pub fn eval_d1(&self, the_t: f64) -> (DVec3, DVec3) {
        let a_ct = the_t.cos();
        let a_st = the_t.sin();
        // Calculate radius at parameter t.
        let a1 = self.r_start + self.c1 * self.tg_beta * the_t;
        // Calculate point coordinates.
        let a_x = a1 * a_ct;
        let mut a_y = a1 * a_st;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = self.c1 * the_t;
        let a_p = DVec3::new(a_x, a_y, a_z);
        // Calculate first derivative coefficients.
        let a1 = self.c1 * self.tg_beta;
        let a2 = self.r_start + a1 * the_t;
        // Calculate first derivative components.
        let a_x = a1 * a_ct - a2 * a_st;
        let mut a_y = a1 * a_st + a2 * a_ct;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = self.c1;
        let a_v1 = DVec3::new(a_x, a_y, a_z);
        (a_p, a_v1)
    }

    /// OCCT EvalD2 (HelixGeom_HelixCurve.cxx L218-259) — point, first and
    /// second derivatives.
    pub fn eval_d2(&self, the_t: f64) -> (DVec3, DVec3, DVec3) {
        let a_ct = the_t.cos();
        let a_st = the_t.sin();
        // Calculate radius at parameter t.
        let a1 = self.r_start + self.c1 * self.tg_beta * the_t;
        // Calculate point coordinates.
        let a_x = a1 * a_ct;
        let mut a_y = a1 * a_st;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = self.c1 * the_t;
        let a_p = DVec3::new(a_x, a_y, a_z);
        // Calculate first derivative coefficients.
        let a1 = self.c1 * self.tg_beta;
        let a2 = self.r_start + a1 * the_t;
        // Calculate first derivative components.
        let a_x = a1 * a_ct - a2 * a_st;
        let mut a_y = a1 * a_st + a2 * a_ct;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = self.c1;
        let a_v1 = DVec3::new(a_x, a_y, a_z);
        // Calculate second derivative.
        let a1 = 2.0 * a1;
        let a_x = -a2 * a_ct - a1 * a_st;
        let mut a_y = -a2 * a_st - a1 * a_ct;
        if !self.is_clock_wise {
            a_y = -a_y;
        }
        let a_z = 0.0;
        let a_v2 = DVec3::new(a_x, a_y, a_z);
        (a_p, a_v1, a_v2)
    }

    /// OCCT EvalDN (HelixGeom_HelixCurve.cxx L262-279).
    pub fn eval_dn(&self, the_t: f64, the_n: i32) -> DVec3 {
        match the_n {
            1 => {
                let (_p, d1) = self.eval_d1(the_t);
                d1
            }
            2 => {
                let (_p, _d1, d2) = self.eval_d2(the_t);
                d2
            }
            _ => panic!("HelixGeom_HelixCurve::EvalDN"),
        }
    }
}
