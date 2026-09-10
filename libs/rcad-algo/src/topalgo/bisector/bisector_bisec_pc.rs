//! OCCT Bisector_BisecPC — the bisector between a point and a curve,
//! 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_BisecPC.hxx (L35-216) / .cxx (L38-885).
//!
//! Architecture differences: `NCollection_Sequence<double>` maps to
//! `Vec<f64>` (OCCT 1-based Value(i) reads index i-1); gp_Trsf2d maps to
//! `glam::DAffine2`; Geom2dAPI_ProjectPointOnCurve is the GAP carrier in
//! super::deps_gap.

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::Curve2d;
use rcad_kernel::math::GeomAbsShape;

use super::bisector::is_convex;
use super::bisector_curve::{
    BisectorCurve, CurveKind, Geom2dCurveAdaptor, GP_RESOLUTION, ResD1, ResD2, ResD3,
};
use super::bisector_bisec_ana::adaptor_of;
use super::deps_gap::Geom2dAPIProjectPointOnCurve;

/// OCCT Bisector_BisecPC (Bisector_BisecPC.hxx L35-216).
pub struct BisectorBisecPC {
    /// OCCT curve (owned copy of the input handle: Cu->Copy()).
    pub(crate) curve: Arc<dyn BisectorCurve>,
    /// OCCT point.
    point: DVec2,
    /// OCCT sign.
    sign: f64,
    /// OCCT startIntervals.
    start_intervals: Vec<f64>,
    /// OCCT endIntervals.
    end_intervals: Vec<f64>,
    /// OCCT bisInterval.
    bis_interval: i32,
    /// OCCT currentInterval.
    current_interval: i32,
    /// OCCT shiftParameter.
    shift_parameter: f64,
    /// OCCT distMax.
    dist_max: f64,
    /// OCCT isEmpty.
    is_empty: bool,
    /// OCCT isConvex.
    is_convex: bool,
    /// OCCT extensionStart.
    extension_start: bool,
    /// OCCT extensionEnd.
    extension_end: bool,
    /// OCCT pointStartBis.
    point_start_bis: DVec2,
    /// OCCT pointEndBis.
    point_end_bis: DVec2,
}

impl Default for BisectorBisecPC {
    /// OCCT Bisector_BisecPC() (L38-48).
    fn default() -> Self {
        BisectorBisecPC {
            curve: Arc::new(super::bisector_curve::Geom2dCurveHandle::new(
                Curve2d::Line(rcad_kernel::geom::Line2d::new(DVec2::ZERO, DVec2::X)),
            )),
            point: DVec2::ZERO,
            sign: 0.0,
            start_intervals: Vec::new(),
            end_intervals: Vec::new(),
            bis_interval: 0,
            current_interval: 0,
            shift_parameter: 0.0,
            dist_max: 0.0,
            is_empty: true,
            is_convex: false,
            extension_start: false,
            extension_end: false,
            point_start_bis: DVec2::ZERO,
            point_end_bis: DVec2::ZERO,
        }
    }
}

/// OCCT `startIntervals.Value(i)` (1-based NCollection_Sequence).
fn seq_value(seq: &[f64], i: i32) -> f64 {
    seq[(i - 1) as usize]
}

impl BisectorBisecPC {
    /// OCCT Bisector_BisecPC() (L38-48).
    pub fn new() -> Self {
        BisectorBisecPC::default()
    }

    /// OCCT Bisector_BisecPC(Cu, P, Side, DistMax) (L53-59).
    pub fn new_full(cu: &Arc<dyn BisectorCurve>, p: DVec2, side: f64, dist_max: f64) -> Self {
        let mut s = BisectorBisecPC::new();
        s.perform(cu, p, side, dist_max);
        s
    }

    /// OCCT Bisector_BisecPC(Cu, P, Side, UMin, UMax) (L63-81).
    pub fn new_trimmed(
        cu: &Arc<dyn BisectorCurve>,
        p: DVec2,
        side: f64,
        u_min: f64,
        u_max: f64,
    ) -> Self {
        let mut s = BisectorBisecPC::new();
        s.curve = cu.copy_curve();
        s.point = p;
        s.sign = side;
        s.start_intervals.push(u_min);
        s.end_intervals.push(u_max);
        s.bis_interval = 1;
        s.extension_start = false;
        s.extension_end = false;
        // OCCT: pointStartBis = Value(UMin); pointEndBis = Value(UMax).
        s.point_start_bis = s.value(u_min);
        s.point_end_bis = s.value(u_max);
        s.is_convex = is_convex(&s.curve, s.sign);
        s
    }

    /// OCCT Perform(Cu, P, Side, DistMax) (L85-132).
    pub fn perform(&mut self, cu: &Arc<dyn BisectorCurve>, p: DVec2, side: f64, dist_max: f64) {
        self.curve = cu.copy_curve();
        self.point = p;
        self.dist_max = dist_max;
        self.sign = side;
        self.is_convex = is_convex(&self.curve, self.sign);
        //--------------------------------------------
        // Calculate interval of definition.
        //--------------------------------------------
        self.compute_intervals();
        if self.is_empty {
            return;
        }

        //-------------------------
        // Construction extensions.
        //-------------------------
        self.bis_interval = 1;
        self.extension_start = false;
        self.extension_end = false;
        self.point_start_bis = self.value(seq_value(&self.start_intervals, 1));
        self.point_end_bis = self.value(self.end_intervals.last().copied().unwrap_or(0.0));

        if !self.is_convex {
            let first_point = first_point_of(&self.curve);
            let last_point = last_point_of(&self.curve);
            if pnt_eq(self.point, first_point, CONFUSION) {
                self.extension_start = true;
                let u_first =
                    seq_value(&self.start_intervals, 1) - p.distance(self.point_start_bis);
                self.start_intervals.insert(0, u_first);
                let v2 = seq_value(&self.start_intervals, 2);
                self.end_intervals.insert(0, v2);
                self.bis_interval = 2;
            } else if pnt_eq(self.point, last_point, CONFUSION) {
                self.extension_end = true;
                let u_last = self.end_intervals.last().copied().unwrap_or(0.0)
                    + p.distance(self.point_end_bis);
                let last_val = self.end_intervals.last().copied().unwrap_or(0.0);
                self.start_intervals.push(last_val);
                self.end_intervals.push(u_last);
                self.bis_interval = 1;
            }
        }
    }

    /// OCCT IsExtendAtStart() (L136-139).
    pub fn is_extend_at_start(&self) -> bool {
        self.extension_start
    }

    /// OCCT IsExtendAtEnd() (L141-144).
    pub fn is_extend_at_end(&self) -> bool {
        self.extension_end
    }

    /// OCCT Reverse() (L150-153) — throws Standard_NotImplemented.
    pub fn reverse(&mut self) {
        panic!("Standard_NotImplemented: Bisector_BisecPC::Reverse");
    }

    /// OCCT ReversedParameter(U) (L157-160).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.last_parameter() + self.first_parameter() - u
    }

    /// OCCT Copy() (L164-185).
    pub fn copy_bisec_pc(&self) -> BisectorBisecPC {
        // OCCT deep-copies the curve then Init()s every member.
        let mut c = BisectorBisecPC::new();
        c.init(
            self.curve.copy_curve(),
            self.point,
            self.sign,
            self.start_intervals.clone(),
            self.end_intervals.clone(),
            self.bis_interval,
            self.current_interval,
            self.shift_parameter,
            self.dist_max,
            self.is_empty,
            self.is_convex,
            self.extension_start,
            self.extension_end,
            self.point_start_bis,
            self.point_end_bis,
        );
        c
    }

    /// OCCT Transform(T) (L189-195) — curve->Transform(T) plus the three
    /// point transforms.
    pub fn transform(&mut self, t: &glam::DAffine2) {
        // GAP: kernel Curve2d lacks a general gp_Trsf2d transform
        // (OCCT Geom2d_Geometry::Transform on the basis curve); the point
        // transforms mirror OCCT once the curve transform lands.
        unimplemented!("GAP: kernel Curve2d lacks a gp_Trsf2d transform");
    }

    /// OCCT IsCN(N) (L199-202).
    pub fn is_cn(&self, n: i32) -> bool {
        // OCCT: curve->IsCN(N + 1).
        self.curve.is_cn(n + 1)
    }

    /// OCCT FirstParameter() (L206-209).
    pub fn first_parameter(&self) -> f64 {
        self.start_intervals.first().copied().unwrap_or(0.0)
    }

    /// OCCT LastParameter() (L213-216).
    pub fn last_parameter(&self) -> f64 {
        self.end_intervals.last().copied().unwrap_or(0.0)
    }

    /// OCCT Continuity() (L220-237).
    pub fn continuity(&self) -> GeomAbsShape {
        // OCCT: Cont = curve->Continuity(); the C1->C0 ... demotion switch.
        let cont = curve2d_continuity(&self.curve);
        match cont {
            GeomAbsShape::C1 => GeomAbsShape::C0,
            GeomAbsShape::C2 => GeomAbsShape::C1,
            GeomAbsShape::C3 => GeomAbsShape::C2,
            GeomAbsShape::CN => GeomAbsShape::CN,
            _ => GeomAbsShape::C0,
        }
    }

    /// OCCT NbIntervals() (L241-244).
    pub fn nb_intervals(&self) -> i32 {
        self.start_intervals.len() as i32
    }

    /// OCCT IntervalFirst(I) (L248-251).
    pub fn interval_first(&self, i: i32) -> f64 {
        seq_value(&self.start_intervals, i)
    }

    /// OCCT IntervalLast(I) (L255-258).
    pub fn interval_last(&self, i: i32) -> f64 {
        seq_value(&self.end_intervals, i)
    }

    /// OCCT IntervalContinuity() (L262-279).
    pub fn interval_continuity(&self) -> GeomAbsShape {
        let cont = curve2d_continuity(&self.curve);
        match cont {
            GeomAbsShape::C1 => GeomAbsShape::C0,
            GeomAbsShape::C2 => GeomAbsShape::C1,
            GeomAbsShape::C3 => GeomAbsShape::C2,
            GeomAbsShape::CN => GeomAbsShape::CN,
            _ => GeomAbsShape::C0,
        }
    }

    /// OCCT IsClosed() (L283-298).
    pub fn is_closed(&self) -> bool {
        // OCCT quirk kept: endIntervals.First() (not Last()).
        if self.curve.is_closed() {
            if seq_value(&self.start_intervals, 1) == self.curve.first_parameter()
                && seq_value(&self.end_intervals, 1) == self.curve.last_parameter()
            {
                return true;
            }
        }
        false
    }

    /// OCCT IsPeriodic() (L302-305).
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT Extension(U, P, V1, V2, V3) (L309-348).
    fn extension(&self, u: f64, p: &mut DVec2, v1: &mut DVec2, v2: &mut DVec2, v3: &mut DVec2) {
        let d_u;

        *v1 = DVec2::ZERO;
        *v2 = DVec2::ZERO;
        *v3 = DVec2::ZERO;
        if u < seq_value(&self.start_intervals, self.bis_interval) {
            if pnt_eq(self.point_start_bis, self.point, PCONFUSION) {
                *p = self.point_start_bis;
            } else {
                let d_u_loc = u - seq_value(&self.start_intervals, self.bis_interval);
                d_u = d_u_loc;
                let dir_ext =
                    DVec2::new(self.point_start_bis.x - self.point.x, self.point_start_bis.y - self.point.y);
                *p = DVec2::new(
                    self.point_start_bis.x + d_u * dir_ext.x,
                    self.point_start_bis.y + d_u * dir_ext.y,
                );
                *v1 = dir_ext;
            }
        } else if u > seq_value(&self.end_intervals, self.bis_interval) {
            if pnt_eq(self.point_end_bis, self.point, PCONFUSION) {
                *p = self.point_end_bis;
            } else {
                let d_u_loc = u - seq_value(&self.end_intervals, self.bis_interval);
                d_u = d_u_loc;
                let dir_ext =
                    DVec2::new(self.point.x - self.point_end_bis.x, self.point.y - self.point_end_bis.y);
                *p = DVec2::new(
                    self.point_end_bis.x + d_u * dir_ext.x,
                    self.point_end_bis.y + d_u * dir_ext.y,
                );
                *v1 = dir_ext;
            }
        }
    }

    /// OCCT Values(U, N, P, V1, V2, V3) (L361-448).
    fn values(&self, u: f64, n: i32, p: &mut DVec2, v1: &mut DVec2, v2: &mut DVec2, v3: &mut DVec2) {
        if u < seq_value(&self.start_intervals, self.bis_interval) {
            self.extension(u, p, v1, v2, v3);
            return;
        } else if u > seq_value(&self.end_intervals, self.bis_interval) {
            self.extension(u, p, v1, v2, v3);
            return;
        }
        let u_on_curve = self.link_bis_curve(u);

        // OCCT switch (N): case 0 -> D1, case 1 -> D2, case 2 -> D3.
        // (N > 2 leaves the derivatives untouched — OCCT reads
        // uninitialized memory there; zeroed here as the defined path.)
        let (pc, tu, tuu, t3u) = match n {
            0 => {
                let (pc, tu) = self.curve.d1(u_on_curve);
                (pc, tu, DVec2::ZERO, DVec2::ZERO)
            }
            1 => {
                let (pc, tu, tuu) = self.curve.d2(u_on_curve);
                (pc, tu, tuu, DVec2::ZERO)
            }
            2 => {
                let (pc, tu, tuu, t3u) = self.curve.d3(u_on_curve);
                (pc, tu, tuu, t3u)
            }
            _ => (DVec2::ZERO, DVec2::ZERO, DVec2::ZERO, DVec2::ZERO),
        };

        let a_ppc = DVec2::new(pc.x - self.point.x, pc.y - self.point.y);
        let nor = DVec2::new(-tu.y, tu.x);

        let square_ppc = a_ppc.length_squared();
        let nor_ppc = nor.dot(a_ppc);
        let a1;

        if nor_ppc.abs() > GP_RESOLUTION && (nor_ppc * self.sign) < 0.0 {
            a1 = 0.5 * square_ppc / nor_ppc;
            *p = DVec2::new(pc.x - nor.x * a1, pc.y - nor.y * a1);
        } else {
            return;
        }

        if n == 0 {
            return; // End Calculation Point;
        }

        let nu = DVec2::new(-tuu.y, tuu.x); // derivative of the normal by U.
        let nu_ppc = nu.dot(a_ppc);
        let tu_ppc = tu.dot(a_ppc);
        let nor_ppce2 = nor_ppc * nor_ppc;
        let a2 = tu_ppc / nor_ppc - 0.5 * nu_ppc * square_ppc / nor_ppce2;

        //--------------------------
        *v1 = tu - a1 * nu - a2 * nor;
        //--------------------------
        if n == 1 {
            return; // End calculation D1.
        }

        let nuu = DVec2::new(-t3u.y, t3u.x);

        let nor_ppce4 = nor_ppce2 * nor_ppce2;
        let nuu_ppc = nuu.dot(a_ppc);
        let tuu_ppc = tuu.dot(a_ppc);

        let a21 = tuu_ppc / nor_ppc - tu_ppc * nu_ppc / nor_ppce2;
        let a22 = (0.5 * nuu_ppc * square_ppc + nu_ppc * tu_ppc) / nor_ppce2
            - nu_ppc * square_ppc * nor_ppc * nu_ppc / nor_ppce4;
        let a2u = a21 - a22;
        //----------------------------------------
        *v2 = tuu - 2.0 * a2 * nu - a1 * nuu - a2u * nor;
        //----------------------------------------
    }

    /// OCCT Distance(U) (L483-541) — squared distance to the curve/point.
    pub fn distance(&self, u: f64) -> f64 {
        let u_on_curve = self.link_bis_curve(u);

        let (pc, tan) = self.curve.d1(u_on_curve);
        let a_ppc = DVec2::new(pc.x - self.point.x, pc.y - self.point.y);
        let nor = DVec2::new(-tan.y, tan.x);

        let nor_nor = nor.length_squared();
        let square_mag_ppc = a_ppc.length_squared();
        let prosca = nor.dot(a_ppc);

        if pnt_eq(self.point, pc, CONFUSION) {
            if self.is_convex {
                return 0.0;
            } else {
                // the point is on a concave curve. The required point is not
                // the common point. This can avoid the discontinuity of the
                // bisectrice.
                return INFINITE_VALUE;
            }
        }

        if prosca.abs() < CONFUSION || (prosca * self.sign) > 0.0 {
            INFINITE_VALUE
        } else {
            let a = 0.5 * square_mag_ppc / prosca;
            let dist = a * a * nor_nor;
            dist
        }
    }

    /// OCCT EvalD0(U) (L545-551).
    pub fn eval_d0(&self, u: f64) -> DVec2 {
        let mut p = self.point;
        let mut v1 = DVec2::ZERO;
        let mut v2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        self.values(u, 0, &mut p, &mut v1, &mut v2, &mut v3);
        p
    }

    /// OCCT EvalD1(U) (L555-563).
    pub fn eval_d1(&self, u: f64) -> ResD1 {
        let mut point = self.point;
        let mut d1 = DVec2::ZERO;
        let mut v2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        self.values(u, 1, &mut point, &mut d1, &mut v2, &mut v3);
        ResD1 { point, d1 }
    }

    /// OCCT EvalD2(U) (L567-576).
    pub fn eval_d2(&self, u: f64) -> ResD2 {
        let mut point = self.point;
        let mut d1 = DVec2::ZERO;
        let mut d2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        self.values(u, 2, &mut point, &mut d1, &mut d2, &mut v3);
        ResD2 { point, d1, d2 }
    }

    /// OCCT EvalD3(U) (L580-589).
    pub fn eval_d3(&self, u: f64) -> ResD3 {
        let mut point = self.point;
        let mut d1 = DVec2::ZERO;
        let mut d2 = DVec2::ZERO;
        let mut d3 = DVec2::ZERO;
        // OCCT passes N=3 into a switch without case 3 (uninitialized reads
        // in C++); zeroed derivatives are the defined mirror here.
        self.values(u, 3, &mut point, &mut d1, &mut d2, &mut d3);
        ResD3 { point, d1, d2, d3 }
    }

    /// OCCT EvalDN(U, N) (L593-616).
    pub fn eval_dn(&self, u: f64, n: i32) -> DVec2 {
        if n < 1 {
            panic!("Geom2d_UndefinedDerivative: Bisector_BisecPC::EvalDN");
        }
        let mut p = self.point;
        let mut v1 = DVec2::ZERO;
        let mut v2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        self.values(u, n, &mut p, &mut v1, &mut v2, &mut v3);
        match n {
            1 => v1,
            2 => v2,
            3 => v3,
            _ => panic!("Geom2d_UndefinedDerivative: Bisector_BisecPC::EvalDN"),
        }
    }

    /// OCCT SearchBound(U1, U2) (L620-645).
    fn search_bound(&self, u1: f64, u2: f64) -> f64 {
        let mut u_mid = 0.0;
        let tol = PCONFUSION;
        let dist_max2 = self.dist_max * self.dist_max;
        let mut u11 = u1;
        let mut u22 = u2;
        let mut dist1 = self.distance(u11);

        while (u22 - u11) > tol {
            u_mid = 0.5 * (u22 + u11);
            let dist_mid = self.distance(u_mid);
            if (dist1 > dist_max2) == (dist_mid > dist_max2) {
                u11 = u_mid;
                dist1 = dist_mid;
            } else {
                u22 = u_mid;
            }
        }
        u_mid
    }

    /// OCCT CuspFilter() (L649-652) — throws Standard_NotImplemented.
    fn cusp_filter(&self) {
        panic!("Standard_NotImplemented: Bisector_BisecPC::CuspFilter");
    }

    /// OCCT ComputeIntervals() (L656-750).
    fn compute_intervals(&mut self) {
        let mut dist1;
        let mut dist2;
        let mut dist_proj;
        self.is_empty = false;
        self.shift_parameter = 0.0;
        let mut ya_proj = false;
        let dist_max2 = self.dist_max * self.dist_max;

        let u1 = self.curve.first_parameter();
        let u2 = self.curve.last_parameter();
        dist1 = self.distance(u1);
        dist2 = self.distance(u2);
        dist_proj = INFINITE_VALUE;

        let mut u_proj = 0.0;
        let mut u_start = 0.0;
        let mut u_end = 0.0;
        let proj = Geom2dAPIProjectPointOnCurve::new(self.point, &adaptor_of_kernel(&self.curve), u1, u2);
        if proj.nb_points() > 0 {
            u_proj = proj.lower_distance_parameter();
            dist_proj = self.distance(u_proj);
            ya_proj = true;
        }

        if dist1 < dist_max2 && dist2 < dist_max2 {
            if dist_proj > dist_max2 && ya_proj {
                self.is_empty = true;
            } else {
                self.start_intervals.push(u1);
                self.end_intervals.push(u2);
            }
            return;
        } else if dist1 > dist_max2 && dist2 > dist_max2 {
            if dist_proj < dist_max2 {
                u_start = self.search_bound(u1, u_proj);
                u_end = self.search_bound(u_proj, u2);
            } else {
                self.is_empty = true;
                return;
            }
        } else if dist1 < dist_max2 {
            u_start = u1;
            u_end = self.search_bound(u1, u2);
        } else if dist2 < dist_max2 {
            u_end = u2;
            u_start = self.search_bound(u1, u2);
        }
        self.start_intervals.push(u_start);
        self.end_intervals.push(u_end);

        //--------------------------------------------------------------------
        // Eventual offset of the parameter on the curve correspondingly to the
        // one on the curve. The offset can be done if the curve is periodical
        // and the point of initial parameter is less then the interval of
        // continuity.
        //--------------------------------------------------------------------
        if curve2d_is_periodic(&self.curve) {
            if self.start_intervals.len() > 1 {
                // Plusieurs intervals.
                if self.end_intervals.last().copied().unwrap_or(0.0) == self.curve.last_parameter()
                    && seq_value(&self.start_intervals, 1) == self.curve.first_parameter()
                {
                    //-----------------------------------------------------
                    // the bissectrice is defined at the origin.
                    // => Fusion of the first and the last interval.
                    //-----------------------------------------------------
                    self.start_intervals.remove(0);
                    self.end_intervals.pop();

                    self.shift_parameter =
                        curve2d_period(&self.curve) - seq_value(&self.start_intervals, 1);
                    for k in 0..self.start_intervals.len() {
                        self.end_intervals[k] += self.shift_parameter;
                        self.start_intervals[k] += self.shift_parameter;
                    }
                    self.start_intervals[0] = 0.0;
                }
            }
        }
        let _ = (&mut dist1, &mut dist2, &mut dist_proj);
    }

    /// OCCT LinkBisCurve(U) (L754-757).
    pub fn link_bis_curve(&self, u: f64) -> f64 {
        u - self.shift_parameter
    }

    /// OCCT LinkCurveBis(U) (L761-764).
    pub fn link_curve_bis(&self, u: f64) -> f64 {
        u + self.shift_parameter
    }

    /// OCCT IsEmpty() (L768-771).
    pub fn is_empty(&self) -> bool {
        self.is_empty
    }

    /// OCCT Parameter(P) (L775-815).
    pub fn parameter(&self, p: DVec2) -> f64 {
        let tol = CONFUSION;

        if pnt_eq(p, self.point_start_bis, tol) {
            return seq_value(&self.start_intervals, self.bis_interval);
        }
        if pnt_eq(p, self.point_end_bis, tol) {
            return seq_value(&self.end_intervals, self.bis_interval);
        }

        if self.extension_start {
            let axe_origin = self.point_start_bis;
            let axe_dir = DVec2::new(
                self.point_start_bis.x - p.x,
                self.point_start_bis.y - p.y,
            );
            let u = rcad_kernel::geom::Line2d::new(axe_origin, axe_dir);
            let u_par = crate::geomalgo::geom2d_int::elclib2d::line_parameter(
                axe_origin,
                axe_dir,
                p,
            );
            let proj = crate::geomalgo::geom2d_int::elclib2d::line_value(
                axe_origin,
                axe_dir,
                u_par,
            );
            let _ = u;
            if pnt_eq(proj, p, tol) && u_par < 0.0 {
                return u_par + seq_value(&self.start_intervals, self.bis_interval);
            }
        }
        if self.extension_end {
            let axe_origin = self.point_end_bis;
            let axe_dir = DVec2::new(p.x - self.point_end_bis.x, p.y - self.point_end_bis.y);
            let u_par =
                crate::geomalgo::geom2d_int::elclib2d::line_parameter(axe_origin, axe_dir, p);
            let proj = crate::geomalgo::geom2d_int::elclib2d::line_value(
                axe_origin,
                axe_dir,
                u_par,
            );
            if pnt_eq(proj, p, tol) && u_par > 0.0 {
                return u_par + seq_value(&self.end_intervals, self.bis_interval);
            }
        }
        let mut u_on_curve = 0.0;
        let proj = Geom2dAPIProjectPointOnCurve::new(
            p,
            &adaptor_of_kernel(&self.curve),
            self.curve.first_parameter(),
            self.curve.last_parameter(),
        );
        if proj.nb_points() > 0 {
            u_on_curve = proj.lower_distance_parameter();
        }
        self.link_curve_bis(u_on_curve)
    }

    /// OCCT Init(...) (L832-863).
    #[allow(clippy::too_many_arguments)]
    pub fn init(
        &mut self,
        curve: Arc<dyn BisectorCurve>,
        point: DVec2,
        sign: f64,
        start_intervals: Vec<f64>,
        end_intervals: Vec<f64>,
        bis_interval: i32,
        current_interval: i32,
        shift_parameter: f64,
        dist_max: f64,
        is_empty: bool,
        is_convex: bool,
        extension_start: bool,
        extension_end: bool,
        point_start_bis: DVec2,
        point_end_bis: DVec2,
    ) {
        self.curve = curve;
        self.point = point;
        self.sign = sign;
        self.start_intervals = start_intervals;
        self.end_intervals = end_intervals;
        self.bis_interval = bis_interval;
        self.current_interval = current_interval;
        self.shift_parameter = shift_parameter;
        self.dist_max = dist_max;
        self.is_empty = is_empty;
        self.is_convex = is_convex;
        self.extension_start = extension_start;
        self.extension_end = extension_end;
        self.point_start_bis = point_start_bis;
        self.point_end_bis = point_end_bis;
    }

    /// OCCT Dump(Deep, Offset) (L868-884).
    pub fn dump(&self, _deep: i32, _offset: i32) {
        println!("Bisector_BisecPC :");
        println!("Point :");
        println!(" X = {}", self.point.x);
        println!(" Y = {}", self.point.y);
        println!("Sign  :{}", self.sign);
        println!("Number Of Intervals :{}", self.start_intervals.len());
        for i in 1..=self.start_intervals.len() {
            println!(
                "Interval number :{}Start :{}  end :{}",
                i,
                self.start_intervals[i - 1],
                self.end_intervals[i - 1]
            );
        }
        println!("Index Current Interval :{}", self.current_interval);
    }
}

// ---------------------------------------------------------------------------
// Curve helpers (through the Geom2d_Curve handle trait).
// ---------------------------------------------------------------------------

fn first_point_of(curve: &Arc<dyn BisectorCurve>) -> DVec2 {
    curve.value(curve.first_parameter())
}

fn last_point_of(curve: &Arc<dyn BisectorCurve>) -> DVec2 {
    curve.value(curve.last_parameter())
}

fn pnt_eq(a: DVec2, b: DVec2, tol: f64) -> bool {
    a.distance(b) <= tol
}

/// OCCT curve->Period() — the conic period.
fn curve2d_period(curve: &Arc<dyn BisectorCurve>) -> f64 {
    curve.period()
}

fn curve2d_is_periodic(curve: &Arc<dyn BisectorCurve>) -> bool {
    curve.is_periodic()
}

/// OCCT curve->Continuity() — GAP: the kernel enum carries no continuity
/// field; the OCCT-typical value is mirrored (analytic = CN).
fn curve2d_continuity(curve: &Arc<dyn BisectorCurve>) -> GeomAbsShape {
    curve.continuity()
}

/// Geom2dAdaptor_Curve over the owned curve handle.
fn adaptor_of_kernel(curve: &Arc<dyn BisectorCurve>) -> Geom2dCurveAdaptor {
    adaptor_of(curve.clone())
}

// ---------------------------------------------------------------------------
// Bisector_Curve trait implementation for BisectorBisecPC
// (OCCT: class Bisector_BisecPC : public Bisector_Curve).
// ---------------------------------------------------------------------------

impl BisectorCurve for BisectorBisecPC {
    fn value(&self, u: f64) -> DVec2 {
        self.eval_d0(u)
    }

    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        let r = self.eval_d1(u);
        (r.point, r.d1)
    }

    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        let r = self.eval_d2(u);
        (r.point, r.d1, r.d2)
    }

    fn d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        let r = self.eval_d3(u);
        (r.point, r.d1, r.d2, r.d3)
    }

    fn dn(&self, u: f64, n: i32) -> DVec2 {
        self.eval_dn(u, n)
    }

    fn first_parameter(&self) -> f64 {
        BisectorBisecPC::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecPC::last_parameter(self)
    }

    fn is_closed(&self) -> bool {
        BisectorBisecPC::is_closed(self)
    }

    fn is_periodic(&self) -> bool {
        false
    }

    fn continuity(&self) -> GeomAbsShape {
        BisectorBisecPC::continuity(self)
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecPC::reversed_parameter(self, u)
    }

    fn reverse(&mut self) {
        BisectorBisecPC::reverse(self)
    }

    fn is_cn(&self, n: i32) -> bool {
        BisectorBisecPC::is_cn(self, n)
    }

    fn transform(&mut self, t: &glam::DAffine2) {
        BisectorBisecPC::transform(self, t)
    }

    fn copy_curve(&self) -> Arc<dyn BisectorCurve> {
        Arc::new(self.copy_bisec_pc())
    }

    fn parameter(&self, p: DVec2) -> f64 {
        BisectorBisecPC::parameter(self, p)
    }

    fn is_extend_at_start(&self) -> bool {
        BisectorBisecPC::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecPC::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecPC::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> f64 {
        BisectorBisecPC::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> f64 {
        BisectorBisecPC::interval_last(self, index)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::BisecPC
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync> {
        self
    }
}
