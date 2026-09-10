//! OCCT Bisector_BisecCC — the bisector between two curves, 1:1
//! translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_BisecCC.hxx (L35-239) / .cxx (L44-1959).
//!
//! Architecture differences: `handle<Geom2d_Curve>` maps to
//! `Arc<dyn BisectorCurve>`; `NCollection_Sequence<double>` maps to
//! `Vec<f64>` (1-based Value(i) reads index i-1); gp_Trsf2d maps to
//! `glam::DAffine2`.

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::Line2d;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::direct_polynomial_roots::epsilon;
use rcad_kernel::math::root::{FunctionRoots, FunctionValue};

use super::bisector::{cross2, dir2d_parallel, is_convex, pnt2d_equal};
use super::bisector_bisec_ana::adaptor_of;
use super::bisector_bisec_pc::BisectorBisecPC;
use super::bisector_curve::{BisectorCurve, CurveKind, GP_RESOLUTION, ResD1, ResD2, ResD3};
use super::bisector_function_h::FunctionH;
use super::bisector_point_on_bis::PointOnBis;
use super::bisector_poly_bis::PolyBis;
use super::deps_gap::{MathBissecNewton, MathFunctionRoot};
use crate::geomalgo::geom2d_int::GInter;

/// OCCT Bisector_BisecCC (Bisector_BisecCC.hxx L35-239).
pub struct BisectorBisecCC {
    /// OCCT curve1.
    curve1: Option<Arc<dyn BisectorCurve>>,
    /// OCCT curve2.
    curve2: Option<Arc<dyn BisectorCurve>>,
    /// OCCT sign1.
    sign1: f64,
    /// OCCT sign2.
    sign2: f64,
    /// OCCT startIntervals.
    start_intervals: Vec<f64>,
    /// OCCT endIntervals.
    end_intervals: Vec<f64>,
    /// OCCT currentInterval.
    current_interval: i32,
    /// OCCT myPolygon.
    my_polygon: PolyBis,
    /// OCCT shiftParameter.
    shift_parameter: f64,
    /// OCCT distMax.
    dist_max: f64,
    /// OCCT isEmpty.
    is_empty: bool,
    /// OCCT isConvex1.
    is_convex1: bool,
    /// OCCT isConvex2.
    is_convex2: bool,
    /// OCCT extensionStart.
    extension_start: bool,
    /// OCCT extensionEnd.
    extension_end: bool,
    /// OCCT pointStart.
    point_start: DVec2,
    /// OCCT pointEnd.
    point_end: DVec2,
}

impl Default for BisectorBisecCC {
    /// OCCT Bisector_BisecCC() (L64-76).
    fn default() -> Self {
        BisectorBisecCC {
            curve1: None,
            curve2: None,
            sign1: 0.0,
            sign2: 0.0,
            start_intervals: Vec::new(),
            end_intervals: Vec::new(),
            current_interval: 0,
            my_polygon: PolyBis::new(),
            shift_parameter: 0.0,
            dist_max: 0.0,
            is_empty: true,
            is_convex1: false,
            is_convex2: false,
            extension_start: false,
            extension_end: false,
            point_start: DVec2::ZERO,
            point_end: DVec2::ZERO,
        }
    }
}

/// OCCT `startIntervals.Value(i)` (1-based NCollection_Sequence).
fn seq_value(seq: &[f64], i: i32) -> f64 {
    seq[(i - 1) as usize]
}

/// OCCT Curvature(C, U, Tol) (L570-586).
fn curvature(c: &Arc<dyn BisectorCurve>, u: f64, tol: f64) -> f64 {
    let (_p, d1, d2) = c.d2(u);
    let norm2 = d1.length_squared();
    if norm2 < tol {
        0.0
    } else {
        cross2(d1, d2) / (norm2 * norm2.sqrt())
    }
}

/// OCCT DiscretPar(DU, EpsMin, EpsMax, NbMin, NbMax, Eps, Nb) (L1931-1959).
fn discret_par(
    du: f64,
    eps_min: f64,
    eps_max: f64,
    nb_min: i32,
    nb_max: i32,
    eps: &mut f64,
    nb: &mut i32,
) -> bool {
    if du <= nb_min as f64 * eps_min {
        *eps = du / (nb_min + 1) as f64;
        *nb = nb_min;
        return false;
    }

    *eps = eps_max.min(du / nb_max as f64);

    if *eps < eps_min {
        *eps = eps_min;
        *nb = (du / eps_min) as i32;
    } else {
        *nb = nb_max;
    }

    true
}

/// OCCT ProjOnCurve(P, C, theParam) (L1785-1833).
fn proj_on_curve(p: DVec2, c: &Arc<dyn BisectorCurve>, the_param: &mut f64) -> bool {
    *the_param = 0.0;

    let (pf, tf) = c.d1(c.first_parameter());
    let (pl, tl) = c.d1(c.last_parameter());

    if pnt2d_equal(p, pf, CONFUSION) {
        *the_param = c.first_parameter();
        return true;
    }

    if pnt2d_equal(p, pl, CONFUSION) {
        *the_param = c.last_parameter();
        return true;
    }

    let ppf = DVec2::new(pf.x - p.x, pf.y - p.y);
    let tf = tf.normalize_or_zero();

    if ppf.dot(tf).abs() < CONFUSION {
        *the_param = c.first_parameter();
        return true;
    }
    let ppl = DVec2::new(pl.x - p.x, pl.y - p.y);
    let tl = tl.normalize_or_zero();
    if ppl.dot(tl).abs() < CONFUSION {
        *the_param = c.last_parameter();
        return true;
    }
    let proj = super::deps_gap::Geom2dAPIProjectPointOnCurve::new(
        p,
        &adaptor_of(c.clone()),
        c.first_parameter(),
        c.last_parameter(),
    );
    if proj.nb_points() > 0 {
        *the_param = proj.lower_distance_parameter();
    } else {
        return false;
    }

    true
}

/// OCCT TestExtension(C1, C2, Start_End) (L1837-1875).
fn test_extension(c1: &Arc<dyn BisectorCurve>, c2: &Arc<dyn BisectorCurve>, start_end: i32) -> bool {
    let mut test = false;
    let (p1, mut t1) = if start_end == 1 {
        c1.d1(c1.first_parameter())
    } else {
        c1.d1(c1.last_parameter())
    };
    let (p2, mut t2) = c2.d1(c2.first_parameter());
    if pnt2d_equal(p1, p2, CONFUSION) {
        t1 = t1.normalize_or_zero();
        t2 = t2.normalize_or_zero();
        if t1.dot(t2) > 1.0 - CONFUSION {
            test = true;
        }
    } else {
        let (p2, t2b) = c2.d1(c2.last_parameter());
        t2 = t2b;
        if pnt2d_equal(p1, p2, CONFUSION) {
            t2 = t2.normalize_or_zero();
            if t1.dot(t2) > 1.0 - CONFUSION {
                test = true;
            }
        }
    }
    test
}

/// OCCT PointByInt(CA, CB, SignA, SignB, UOnA, UOnB, Dist) (L1328-1481).
fn point_by_int(
    ca: &Arc<dyn BisectorCurve>,
    cb: &Arc<dyn BisectorCurve>,
    sign_a: f64,
    sign_b: f64,
    u_on_a: f64,
    u_on_b: &mut f64,
    dist: &mut f64,
) -> bool {
    let is_convex_a = is_convex(ca, sign_a);
    let is_convex_b = is_convex(cb, sign_b);

    let (p1, tan1) = ca.d1(u_on_a);
    let n1 = DVec2::new(tan1.y, -tan1.x);

    //-------------------------------------------------------------------
    // test of confusion of P1 with extremity of curve2.
    //-------------------------------------------------------------------
    if p1.distance(cb.value(cb.first_parameter())) < CONFUSION {
        let u_on_b_loc = cb.first_parameter();
        let (_, tan2) = cb.d1(u_on_b_loc);
        if is_convex_a && is_convex_b {
            *dist = 0.0;
            *u_on_b = u_on_b_loc;
            return true;
        }
        if !dir2d_parallel(tan1, tan2, ANGULAR) {
            *dist = 0.0;
            return false;
        }
    }
    if p1.distance(cb.value(cb.last_parameter())) < CONFUSION {
        let u_on_b_loc = cb.last_parameter();
        let (_, tan2) = cb.d1(u_on_b_loc);
        if is_convex_a && is_convex_b {
            *dist = 0.0;
            *u_on_b = u_on_b_loc;
            return true;
        }
        if !dir2d_parallel(tan1, tan2, ANGULAR) {
            *dist = 0.0;
            return false;
        }
    }

    let mut d_min = INFINITE_VALUE;
    let mut ya_sol = false;
    let mut p_sol = DVec2::ZERO;
    //--------------------------------------------------------------------
    // Construction of the bisectrice point curve and of the straight line
    // passing through P1 and carried by the normal.
    //--------------------------------------------------------------------
    let bis_pc = Arc::new(BisectorBisecPC::new_full(cb, p1, sign_b, 500.0));
    //-------------------------------
    // Test if the bissectrice exists.
    //-------------------------------
    if bis_pc.is_empty() {
        *dist = INFINITE_VALUE;
        return false;
    }

    let nor_li = Line2d::new(p1, n1);

    let a_bis_pc = adaptor_of(bis_pc.clone());
    let a_nor_li = adaptor_of(Arc::new(super::bisector_curve::Geom2dCurveHandle::line(&nor_li)));
    let intersect = GInter::new_cc(&a_bis_pc, &a_nor_li, CONFUSION, CONFUSION);

    if intersect.is_done() && !intersect.base.is_empty() {
        for i in 1..=intersect.nb_points() {
            if intersect.point(i).param_on_second() * sign_a < PCONFUSION {
                let p = intersect.point(i).value();
                if p.distance_squared(p1) < d_min {
                    d_min = p.distance_squared(p1);
                    p_sol = p;
                    let upc = intersect.point(i).param_on_first();
                    *u_on_b = bis_pc.link_bis_curve(upc);
                    *dist = d_min;
                    ya_sol = true;
                }
            }
        }
    }
    if ya_sol {
        //--------------------------------------------------------------
        // Point found => Test distance curvature + Angular test
        //---------------------------------------------------------------
        let p2 = cb.value(*u_on_b);
        if p1.distance_squared(p_sol) < 1.0e-32 {
            return false;
        }
        if p2.distance_squared(p_sol) < 1.0e-32 {
            return false;
        }

        let pp1_unit = DVec2::new(p1.x - p_sol.x, p1.y - p_sol.y).normalize_or_zero();
        let pp2_unit = DVec2::new(p2.x - p_sol.x, p2.y - p_sol.y).normalize_or_zero();

        if pp1_unit.dot(pp2_unit) > 1.0 - ANGULAR {
            ya_sol = false;
        } else {
            *dist = dist.sqrt();
            if !is_convex_a {
                let k1 = curvature(ca, u_on_a, CONFUSION);
                if k1 != 0.0 && *dist > (1.0 / k1).abs() {
                    ya_sol = false;
                }
            }
            if ya_sol && !is_convex_b {
                let k2 = curvature(cb, *u_on_b, CONFUSION);
                if k2 != 0.0 && *dist > (1.0 / k2).abs() {
                    ya_sol = false;
                }
            }
        }
    }
    ya_sol
}

impl BisectorBisecCC {
    /// OCCT Bisector_BisecCC() (L64-76).
    pub fn new() -> Self {
        BisectorBisecCC::default()
    }

    /// OCCT Bisector_BisecCC(Cu1, Cu2, Side1, Side2, Origin, DistMax)
    /// (L80-88).
    pub fn new_full(
        cu1: &Arc<dyn BisectorCurve>,
        cu2: &Arc<dyn BisectorCurve>,
        side1: f64,
        side2: f64,
        origin: DVec2,
        dist_max: f64,
    ) -> Self {
        let mut s = BisectorBisecCC::new();
        s.perform(cu1, cu2, side1, side2, origin, dist_max);
        s
    }

    /// OCCT Perform(Cu1, Cu2, Side1, Side2, Origin, DistMax) (L92-327).
    pub fn perform(
        &mut self,
        cu1: &Arc<dyn BisectorCurve>,
        cu2: &Arc<dyn BisectorCurve>,
        side1: f64,
        side2: f64,
        origin: DVec2,
        dist_max: f64,
    ) {
        self.is_empty = false;
        self.dist_max = dist_max;

        self.curve1 = Some(cu1.copy_curve());
        self.curve2 = Some(cu2.copy_curve());

        let curve1 = self.curve1.clone().unwrap();
        let curve2 = self.curve2.clone().unwrap();

        self.sign1 = side1;
        self.sign2 = side2;

        self.is_convex1 = is_convex(&curve1, self.sign1);
        self.is_convex2 = is_convex(&curve2, self.sign2);

        let mut u;
        let mut uc1 = 0.0;
        let mut uc2 = 0.0;
        let mut dist = 0.0;
        let mut d_u;
        let mut u_sol;
        let mut p;
        let mut nb_pnts = 21;
        let eps_min = 10.0 * CONFUSION;
        let mut ya_poly = true;
        let mut ori_in_poly = false;
        //---------------------------------------------
        // Calculate first point of the polygon.
        //---------------------------------------------
        let mut u_proj = 0.0;
        let is_proj_done = proj_on_curve(origin, &curve1, &mut u_proj);
        u = u_proj;

        if !is_proj_done {
            self.is_empty = true;
            return;
        }

        p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
        if dist < CONFUSION {
            let a_p1 = curve1.value(curve1.last_parameter());
            let a_p2 = curve2.value(curve2.first_parameter());
            let dp = a_p1.distance(p) + a_p2.distance(p);
            let dorig = a_p1.distance(origin) + a_p2.distance(origin);
            if dp < dorig {
                self.is_empty = true;
                return;
            }
        }

        if dist < INFINITE_VALUE {
            //----------------------------------------------------
            // the parameter of the origin point gives a point
            // on the polygon.
            //----------------------------------------------------
            self.my_polygon
                .append(PointOnBis::new_full(uc1, uc2, u, dist, p));
            self.start_intervals.push(u);
            if pnt2d_equal(p, origin, CONFUSION) {
                //----------------------------------------
                // test if the first point is the origin.
                //----------------------------------------
                ori_in_poly = true;
            }
        } else {
            //-------------------------------------------------------
            // The origin point is on the extension.
            // Find the first point of the polygon by dichotomy.
            //-------------------------------------------------------
            d_u = (curve1.last_parameter() - u) / (nb_pnts - 1) as f64;
            u += d_u;
            for _i in 1..=nb_pnts - 1 {
                p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
                if dist < INFINITE_VALUE {
                    u_sol = self.search_bound(u - d_u, u);
                    p = self.value_by_int(u_sol, &mut uc1, &mut uc2, &mut dist);
                    self.start_intervals.push(u_sol);
                    self.my_polygon
                        .append(PointOnBis::new_full(uc1, uc2, u_sol, dist, p));
                    break;
                }
                u += d_u;
            }
        }

        if self.my_polygon.length() != 0 {
            self.sup_last_parameter();
            //----------------------------------------------
            // Construction of the polygon of the bissectrice.
            //---------------------------------------------
            u = self.first_parameter();
            let d_u_total = self.last_parameter() - u;

            if d_u_total < eps_min {
                nb_pnts = 3;
            }
            d_u = d_u_total / (nb_pnts - 1) as f64;

            u += d_u;
            // prevent addition of the same point.
            let mut prev_pnt = p;
            for _i in 1..=nb_pnts - 1 {
                p = self.value_by_int(u, &mut uc1, &mut uc2, &mut dist);
                if dist < INFINITE_VALUE {
                    if p.distance(prev_pnt) > CONFUSION {
                        self.my_polygon
                            .append(PointOnBis::new_full(uc1, uc2, u, dist, p));
                    }
                } else {
                    u_sol = self.search_bound(u - d_u, u);
                    p = self.value_by_int(u_sol, &mut uc1, &mut uc2, &mut dist);
                    self.end_intervals[0] = u_sol;
                    if p.distance(prev_pnt) > CONFUSION {
                        self.my_polygon
                            .append(PointOnBis::new_full(uc1, uc2, u_sol, dist, p));
                    }
                    break;
                }
                u += d_u;
                prev_pnt = p;
            }
        } else {
            //----------------
            // Empty Polygon.
            //----------------
            ya_poly = false;
        }

        self.extension_start = false;
        self.extension_end = false;
        self.point_start = origin;

        if self.is_convex1 && self.is_convex2 {
            if ya_poly {
                self.point_end = self.my_polygon.last().point();
            }
        } else {
            //---------------------------------------------------------------
            // Extension : The curve is extended at the beginning and/or the
            // end if - one of two curves is concave. - the curves have a
            // common point at the beginning and/or the end - the angle of
            // opening at the common point between two curves values M_PI.
            // the extension at the beginning is taken into account if the
            // origin is found above. ie : the origin is not the in the
            // polygon.
            //---------------------------------------------------------------

            //---------------------------------
            // Do the extensions exist ?
            //---------------------------------
            if ori_in_poly {
                self.extension_start = false;
            } else {
                self.extension_start = test_extension(&curve1, &curve2, 1);
            }
            self.extension_end = test_extension(&curve1, &curve2, 2);

            //-----------------
            // Calculate pointEnd.
            //-----------------
            if self.extension_end {
                self.point_end = curve1.value(curve1.last_parameter());
            } else if ya_poly {
                self.point_end = self.my_polygon.last().point();
            } else {
                self.compute_point_end();
            }
            //------------------------------------------------------
            // Update the Limits of intervals of definition.
            //------------------------------------------------------
            if ya_poly {
                if self.extension_start {
                    let p1 = self.my_polygon.first().point();
                    let u_first =
                        seq_value(&self.start_intervals, 1) - self.point_start.distance(p1);
                    self.start_intervals.insert(0, u_first);
                    let v2 = seq_value(&self.start_intervals, 2);
                    self.end_intervals.insert(0, v2);
                }
                if self.extension_end {
                    let p1 = self.my_polygon.last().point();
                    let u_first = *self.end_intervals.last().unwrap();
                    let u_last = u_first + self.point_end.distance(p1);
                    self.start_intervals.push(u_first);
                    self.end_intervals.push(u_last);
                }
            } else {
                //--------------------------------------------------
                // No polygon => the bissectrice is a segment.
                //--------------------------------------------------
                self.start_intervals.push(0.0);
                self.end_intervals
                    .push(self.point_end.distance(self.point_start));
            }
        }
        if !ya_poly && !self.extension_start && !self.extension_end {
            self.is_empty = true;
        }
        if self.my_polygon.length() <= 2 {
            self.is_empty = true;
        }
    }

    /// OCCT IsExtendAtStart() (L331-334).
    pub fn is_extend_at_start(&self) -> bool {
        self.extension_start
    }

    /// OCCT IsExtendAtEnd() (L338-341).
    pub fn is_extend_at_end(&self) -> bool {
        self.extension_end
    }

    /// OCCT IsEmpty() (L345-348).
    pub fn is_empty(&self) -> bool {
        self.is_empty
    }

    /// OCCT Reverse() (L352-355) — throws Standard_NotImplemented.
    pub fn reverse(&mut self) {
        panic!("Standard_NotImplemented: Bisector_BisecCC::Reverse");
    }

    /// OCCT ReversedParameter(U) (L359-362).
    pub fn reversed_parameter(&self, u: f64) -> f64 {
        self.last_parameter() + self.first_parameter() - u
    }

    /// OCCT Copy() (L366-390).
    pub fn copy_bisec_cc(&self) -> BisectorBisecCC {
        let mut c = BisectorBisecCC::new();

        c.set_curve(1, self.curve(1));
        c.set_curve(2, self.curve(2));
        c.set_sign(1, self.sign1);
        c.set_sign(2, self.sign2);
        c.set_is_convex(1, self.is_convex1);
        c.set_is_convex(2, self.is_convex2);
        c.set_polygon(self.my_polygon.clone());
        c.set_is_empty(self.is_empty);
        c.set_dist_max(self.dist_max);
        c.set_start_intervals(self.start_intervals.clone());
        c.set_end_intervals(self.end_intervals.clone());
        c.set_extension_start(self.extension_start);
        c.set_extension_end(self.extension_end);
        c.set_point_start(self.point_start);
        c.set_point_end(self.point_end);

        c
    }

    /// OCCT ChangeGuide() (L399-449) — same bisector with the curves
    /// inversed.
    pub fn change_guide(&self) -> Arc<BisectorBisecCC> {
        let mut c = BisectorBisecCC::new();

        c.set_curve(1, self.curve(2));
        c.set_curve(2, self.curve(1));
        c.set_sign(1, self.sign2);
        c.set_sign(2, self.sign1);
        c.set_is_convex(1, self.is_convex2);
        c.set_is_convex(2, self.is_convex1);

        //------------------------------------------------------------------
        // Construction of the new polygon from the initial one. inversion of
        // PointOnBis and Calculation of new parameters on the bissectrice.
        //------------------------------------------------------------------
        let mut poly = PolyBis::new();
        if self.sign1 == self.sign2 {
            //---------------------------------------------------------------
            // elements of the new polygon are ranked in the other direction.
            //---------------------------------------------------------------
            for i in (1..=self.my_polygon.length()).rev() {
                let p = self.my_polygon.value(i);
                let new_p = PointOnBis::new_full(
                    p.param_on_c2(),
                    p.param_on_c1(),
                    p.param_on_c2(),
                    p.distance(),
                    p.point(),
                );
                poly.append(new_p);
            }
        } else {
            for i in 1..=self.my_polygon.length() {
                let p = self.my_polygon.value(i);
                let new_p = PointOnBis::new_full(
                    p.param_on_c2(),
                    p.param_on_c1(),
                    p.param_on_c2(),
                    p.distance(),
                    p.point(),
                );
                poly.append(new_p);
            }
        }
        c.set_polygon(poly.clone());
        c.append_first_parameter(poly.first().param_on_bis());
        c.append_last_parameter(poly.last().param_on_bis());

        Arc::new(c)
    }

    /// OCCT Transform(T) (L453-460).
    pub fn transform(&mut self, t: &glam::DAffine2) {
        // GAP: kernel Curve2d lacks a general gp_Trsf2d transform (OCCT
        // Geom2d_Geometry::Transform on curve1/curve2).
        unimplemented!("GAP: kernel Curve2d lacks a gp_Trsf2d transform");
    }

    /// OCCT IsCN(N) (L464-467).
    pub fn is_cn(&self, n: i32) -> bool {
        self.curve(1).is_cn(n + 1) && self.curve(2).is_cn(n + 1)
    }

    /// OCCT FirstParameter() (L471-474).
    pub fn first_parameter(&self) -> f64 {
        self.start_intervals.first().copied().unwrap_or(0.0)
    }

    /// OCCT LastParameter() (L478-481).
    pub fn last_parameter(&self) -> f64 {
        self.end_intervals.last().copied().unwrap_or(0.0)
    }

    /// OCCT Continuity() (L485-502).
    pub fn continuity(&self) -> GeomAbsShape {
        let cont = self.curve(1).continuity();
        match cont {
            GeomAbsShape::C1 => GeomAbsShape::C0,
            GeomAbsShape::C2 => GeomAbsShape::C1,
            GeomAbsShape::C3 => GeomAbsShape::C2,
            GeomAbsShape::CN => GeomAbsShape::CN,
            _ => GeomAbsShape::C0,
        }
    }

    /// OCCT NbIntervals() (L506-509).
    pub fn nb_intervals(&self) -> i32 {
        self.start_intervals.len() as i32
    }

    /// OCCT IntervalFirst(Index) (L513-516).
    pub fn interval_first(&self, index: i32) -> f64 {
        seq_value(&self.start_intervals, index)
    }

    /// OCCT IntervalLast(Index) (L520-523).
    pub fn interval_last(&self, index: i32) -> f64 {
        seq_value(&self.end_intervals, index)
    }

    /// OCCT IntervalContinuity() (L527-544).
    pub fn interval_continuity(&self) -> GeomAbsShape {
        let cont = self.curve(1).continuity();
        match cont {
            GeomAbsShape::C1 => GeomAbsShape::C0,
            GeomAbsShape::C2 => GeomAbsShape::C1,
            GeomAbsShape::C3 => GeomAbsShape::C2,
            GeomAbsShape::CN => GeomAbsShape::CN,
            _ => GeomAbsShape::C0,
        }
    }

    /// OCCT IsClosed() (L548-559).
    pub fn is_closed(&self) -> bool {
        if self.curve(1).is_closed() {
            if self.start_intervals.first().copied().unwrap_or(0.0)
                == self.curve(1).first_parameter()
                && self.end_intervals.last().copied().unwrap_or(0.0)
                    == self.curve(1).last_parameter()
            {
                return true;
            }
        }
        false
    }

    /// OCCT IsPeriodic() (L563-566).
    pub fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT ValueAndDist(U, U1, U2, Dist) (L608-773) — the current point by
    /// iterative method.
    pub fn value_and_dist(&self, u: f64, u1: &mut f64, u2: &mut f64, dist: &mut f64) -> DVec2 {
        //-----------------------------------------------
        // is the polygon reduced to a point or empty?
        //-----------------------------------------------
        if self.my_polygon.length() <= 1 {
            // OCCT passes a local T by reference; the value is unused here.
            let mut _t = DVec2::ZERO;
            return self.extension(u, u1, u2, dist, &mut _t);
        }

        //-----------------------------------------------
        // test U out of the limits of the polygon.
        //-----------------------------------------------
        if u < self.my_polygon.first().param_on_bis() {
            let mut t = DVec2::ZERO;
            return self.extension(u, u1, u2, dist, &mut t);
        }
        if u > self.my_polygon.last().param_on_bis() {
            let mut t = DVec2::ZERO;
            return self.extension(u, u1, u2, dist, &mut t);
        }

        //-------------------------------------------------------
        // Find start parameter by using <myPolygon>.
        //-------------------------------------------------------
        let interval_index = self.my_polygon.interval(u);
        let u_min = self.my_polygon.value(interval_index).param_on_bis();
        let u_max = self.my_polygon.value(interval_index + 1).param_on_bis();
        let v_min = self.my_polygon.value(interval_index).param_on_c2();
        let v_max = self.my_polygon.value(interval_index + 1).param_on_c2();
        let v_init;

        if (u_max - u_min).abs() < GP_RESOLUTION {
            v_init = v_min;
        } else {
            let alpha = (u - u_min) / (u_max - u_min);
            v_init = v_min + alpha * (v_max - v_min);
        }

        *u1 = self.link_bis_curve(u);
        let curve1 = self.curve(1);
        let curve2 = self.curve(2);
        let (v_min, v_max) = if v_min <= v_max {
            (v_min, v_max)
        } else {
            (v_max, v_min)
        };
        let mut valid = true;
        //---------------------------------------------------------------
        // Calculate parameter U2 on curve C2 solution of H(u,v)=0
        //---------------------------------------------------------------
        let (p1, t1) = curve1.d1(*u1);
        let n1 = DVec2::new(t1.y, -t1.x);

        let mut u2_out;
        if (v_max - v_min) < PCONFUSION {
            u2_out = v_init;
        } else {
            let mut h = FunctionH::new(curve2.clone(), p1, self.sign1 * self.sign2 * t1);
            let eps_h = 1.0e-9;
            let eps_h100 = 1.0e-7;
            let f_init = h.value(v_init).unwrap_or(0.0);
            if f_init.abs() < eps_h {
                u2_out = v_init;
            } else {
                let mut a_new_solution = MathBissecNewton::new(eps_h);
                a_new_solution.perform(&mut h, v_min - eps_h100, v_max + eps_h100, 10);

                if a_new_solution.is_done() {
                    u2_out = a_new_solution.root();
                } else {
                    let mut sol_root =
                        MathFunctionRoot::new(&mut h, v_init, eps_h, v_min - eps_h100, v_max + eps_h100);

                    if sol_root.is_done() {
                        u2_out = sol_root.root();
                    } else {
                        valid = false;
                        u2_out = 0.0;
                    }
                }
            }
        }

        let mut p_bis = self.point_start;
        let mut dist_out = 0.0;
        //----------------
        // P(U) = F(U1,U2)
        //----------------
        if valid {
            let p2 = curve2.value(u2_out);
            let p2p1 = DVec2::new(p1.x - p2.x, p1.y - p2.y);
            let square_p2p1 = p2p1.length_squared();
            let n1p2p1 = n1.dot(p2p1);
            let an_eps = epsilon(1.0);

            if pnt2d_equal(p1, p2, CONFUSION) {
                p_bis = p1;
                dist_out = 0.0;
            } else if n1p2p1 * self.sign1 < an_eps {
                valid = false;
            } else {
                p_bis = p1 - n1 * (0.5 * square_p2p1 / n1p2p1);
                dist_out = p1.distance_squared(p_bis);
            }
        }

        //----------------------------------------------------------------
        // If the point is not valid calculate by intersection.
        //----------------------------------------------------------------
        if !valid {
            //--------------------------------------------------------------------
            // Construction of the bisectrice point curve and of the straight
            // line passing by P1 and carried by the normal. curve2 is
            // limited by VMin and VMax.
            //--------------------------------------------------------------------
            let mut d_min = INFINITE_VALUE;

            let bis_pc = Arc::new(BisectorBisecPC::new_trimmed(
                &curve2,
                p1,
                self.sign2,
                v_min,
                v_max,
            ));
            let nor_li = Line2d::new(p1, n1);

            let a_bis_pc = adaptor_of(bis_pc.clone());
            let a_nor_li =
                adaptor_of(Arc::new(super::bisector_curve::Geom2dCurveHandle::line(&nor_li)));
            let intersect = GInter::new_cc(&a_bis_pc, &a_nor_li, CONFUSION, CONFUSION);

            if intersect.is_done() && !intersect.base.is_empty() {
                for i in 1..=intersect.nb_points() {
                    if intersect.point(i).param_on_second() * self.sign1 < PCONFUSION {
                        let p = intersect.point(i).value();
                        if p.distance_squared(p1) < d_min {
                            d_min = p.distance_squared(p1);
                            p_bis = p;
                            u2_out = bis_pc
                                .link_bis_curve(intersect.point(i).param_on_first());
                            dist_out = d_min;
                        }
                    }
                }
            }
        }
        *u2 = u2_out;
        *dist = dist_out;
        p_bis
    }

    /// OCCT ValueByInt(U, U1, U2, Dist) (L788-988) — the current point by
    /// intersection.
    pub fn value_by_int(&self, u: f64, u1: &mut f64, u2: &mut f64, dist: &mut f64) -> DVec2 {
        //------------------------------------------------------------------
        // Return point, tangent, normal on C1 at parameter U.
        //-------------------------------------------------------------------
        *u1 = self.link_bis_curve(u);
        let curve1 = self.curve(1);
        let curve2 = self.curve(2);

        let (p1, tan1) = curve1.d1(*u1);
        let n1 = DVec2::new(tan1.y, -tan1.x);

        //--------------------------------------------------------------------------
        // test confusion of P1 with extremity of curve2.
        //--------------------------------------------------------------------------
        if p1.distance(curve2.value(curve2.first_parameter())) < CONFUSION {
            let u2_loc = curve2.first_parameter();
            let (_, tan2) = curve2.d1(u2_loc);
            if self.is_convex1 && self.is_convex2 {
                *dist = 0.0;
                *u2 = u2_loc;
                return p1;
            }
            if !dir2d_parallel(tan1, tan2, ANGULAR) {
                *dist = 0.0;
                return p1;
            }
        }
        if p1.distance(curve2.value(curve2.last_parameter())) < CONFUSION {
            let u2_loc = curve2.last_parameter();
            let (_, tan2) = curve2.d1(u2_loc);
            if self.is_convex1 && self.is_convex2 {
                *dist = 0.0;
                *u2 = u2_loc;
                return p1;
            }
            if !dir2d_parallel(tan1, tan2, ANGULAR) {
                *dist = 0.0;
                return p1;
            }
        }

        let mut ya_sol = false;
        let mut d_min = INFINITE_VALUE;
        let mut p_sol = DVec2::ZERO;
        let eps_max = 1.0e-6;
        let eps_x;
        let eps_h = 1.0e-8;
        let mut nb_samples = 20;
        let mut u_first_on_c2 = curve2.first_parameter();
        let mut u_last_on_c2 = curve2.last_parameter();

        if !self.my_polygon.is_empty() {
            if self.sign1 == self.sign2 {
                u_last_on_c2 = self.my_polygon.last().param_on_c2();
            } else {
                u_first_on_c2 = self.my_polygon.last().param_on_c2();
            }
        }

        if (u_last_on_c2 - u_first_on_c2).abs() < PCONFUSION / 100.0 {
            *dist = INFINITE_VALUE;
            return p1;
        }

        let mut eps_x_out = 0.0;
        discret_par(
            (u_last_on_c2 - u_first_on_c2).abs(),
            eps_h,
            eps_max,
            2,
            20,
            &mut eps_x_out,
            &mut nb_samples,
        );
        eps_x = eps_x_out;

        let mut h = FunctionH::new(curve2.clone(), p1, self.sign1 * self.sign2 * tan1);
        let mut sol_root = FunctionRoots::new(
            &mut h,
            u_first_on_c2,
            u_last_on_c2,
            nb_samples,
            eps_x,
            eps_h,
            eps_h,
            0.0,
        );
        if sol_root.is_done() {
            for j in 1..=sol_root.nb_solutions() {
                let u_sol = sol_root.value(j);
                let p2_curve2 = curve2.value(u_sol);
                let p2p1 = DVec2::new(p1.x - p2_curve2.x, p1.y - p2_curve2.y);
                let square_p2p1 = p2p1.length_squared();
                let n1p2p1 = n1.dot(p2p1);

                // Test if the solution is at the proper side of the curves.
                if n1p2p1 * self.sign1 > 0.0 {
                    let p = p1 - n1 * (0.5 * square_p2p1 / n1p2p1);
                    let dist_pp1 = p1.distance_squared(p);
                    if dist_pp1 < d_min {
                        d_min = dist_pp1;
                        p_sol = p;
                        *u2 = u_sol;
                        ya_sol = true;
                    }
                }
            }
        }

        if ya_sol {
            *dist = d_min;
            //--------------------------------------------------------------
            // Point found => Test curve distance + Angular Test
            //---------------------------------------------------------------
            let p2 = curve2.value(*u2);
            let pp1 = DVec2::new(p1.x - p_sol.x, p1.y - p_sol.y);
            let pp2 = DVec2::new(p2.x - p_sol.x, p2.y - p_sol.y);

            //-----------------------------------------------
            // Dist = product of norms = distance at the square.
            //-----------------------------------------------
            if pp1.dot(pp2) > (1.0 - ANGULAR) * *dist {
                ya_sol = false;
            } else {
                if !self.is_convex1 {
                    let k1 = curvature(&curve1, *u1, CONFUSION);
                    if k1 != 0.0 && *dist > 1.0 / (k1 * k1) {
                        ya_sol = false;
                    }
                }
                if ya_sol && !self.is_convex2 {
                    let k2 = curvature(&curve2, *u2, CONFUSION);
                    if k2 != 0.0 && *dist > 1.0 / (k2 * k2) {
                        ya_sol = false;
                    }
                }
            }
        }
        if !ya_sol {
            *dist = INFINITE_VALUE;
            p_sol = p1;
        }
        p_sol
    }

    /// OCCT EvalD0(U) (L992-997).
    pub fn eval_d0(&self, u: f64) -> DVec2 {
        let mut u1 = 0.0;
        let mut u2 = 0.0;
        let mut dist = 0.0;

        self.value_and_dist(u, &mut u1, &mut u2, &mut dist)
    }

    /// OCCT EvalD1(U) (L1001-1008).
    pub fn eval_d1(&self, u: f64) -> ResD1 {
        let mut v2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        let mut point = DVec2::ZERO;
        let mut d1 = DVec2::ZERO;
        self.values(u, 1, &mut point, &mut d1, &mut v2, &mut v3);
        ResD1 { point, d1 }
    }

    /// OCCT EvalD2(U) (L1012-1020).
    pub fn eval_d2(&self, u: f64) -> ResD2 {
        let mut v3 = DVec2::ZERO;
        let mut point = DVec2::ZERO;
        let mut d1 = DVec2::ZERO;
        let mut d2 = DVec2::ZERO;
        self.values(u, 2, &mut point, &mut d1, &mut d2, &mut v3);
        ResD2 { point, d1, d2 }
    }

    /// OCCT EvalD3(U) (L1024-1032).
    pub fn eval_d3(&self, u: f64) -> ResD3 {
        let mut point = DVec2::ZERO;
        let mut d1 = DVec2::ZERO;
        let mut d2 = DVec2::ZERO;
        let mut d3 = DVec2::ZERO;
        self.values(u, 3, &mut point, &mut d1, &mut d2, &mut d3);
        ResD3 { point, d1, d2, d3 }
    }

    /// OCCT EvalDN(U, N) (L1036-1059).
    pub fn eval_dn(&self, u: f64, n: i32) -> DVec2 {
        if n < 1 {
            panic!("Geom2d_UndefinedDerivative: Bisector_BisecCC::EvalDN");
        }
        let mut p = DVec2::ZERO;
        let mut v1 = DVec2::ZERO;
        let mut v2 = DVec2::ZERO;
        let mut v3 = DVec2::ZERO;
        self.values(u, n, &mut p, &mut v1, &mut v2, &mut v3);
        match n {
            1 => v1,
            2 => v2,
            3 => v3,
            _ => panic!("Geom2d_UndefinedDerivative: Bisector_BisecCC::EvalDN"),
        }
    }

    /// OCCT Values(U, N, P, V1, V2, V3) (L1081-1207).
    fn values(&self, u: f64, n: i32, p: &mut DVec2, v1: &mut DVec2, v2: &mut DVec2, v3: &mut DVec2) {
        *v1 = DVec2::ZERO;
        *v2 = DVec2::ZERO;
        *v3 = DVec2::ZERO;
        //------------------------------------------------------------------
        // Calculate the current point on the bisectrice and the parameters
        // on each curve.
        //------------------------------------------------------------------
        let mut u0 = 0.0;
        let mut v0 = 0.0;
        let mut dist = 0.0;

        //-----------------------------------------------
        // is the polygon reduced to a point or empty?
        //-----------------------------------------------
        if self.my_polygon.length() <= 1 {
            let mut t = DVec2::ZERO;
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, &mut t);
        }
        if u < self.my_polygon.first().param_on_bis() {
            let mut t = DVec2::ZERO;
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, &mut t);
            return;
        }
        if u > self.my_polygon.last().param_on_bis() {
            let mut t = DVec2::ZERO;
            *p = self.extension(u, &mut u0, &mut v0, &mut dist, &mut t);
            return;
        }
        *p = self.value_and_dist(u, &mut u0, &mut v0, &mut dist);

        if n == 0 {
            return;
        }
        //------------------------------------------------------------------
        // Return point, tangent, normal to C1 by parameter U0.
        //-------------------------------------------------------------------
        let curve1 = self.curve(1);
        let curve2 = self.curve(2);
        let (p1, tu, tuu) = curve1.d2(u0);
        let nor = DVec2::new(-tu.y, tu.x); // Normal by U0.
        let nu = DVec2::new(-tuu.y, tuu.x); // derivative of the normal by U0.

        //-------------------------------------------------------------------
        // Return point, tangent, normale to C2 by parameter V0.
        //-------------------------------------------------------------------
        let (p2, tv, tvv) = curve2.d2(v0);

        let pu_pv = DVec2::new(p2.x - p1.x, p2.y - p1.y);

        //-----------------------------
        // Calculate dH/du and dH/dv.
        //-----------------------------
        let tu_tu = tu.dot(tu);
        let tv_tv = tv.dot(tv);
        let tu_tv = tu.dot(tv);
        let tu_pu_pv = tu.dot(pu_pv);
        let tv_pu_pv = tv.dot(pu_pv);
        let tuu_pu_pv = tuu.dot(pu_pv);
        let tu_tuu = tu.dot(tuu);
        let tvv_pu_pv = tvv.dot(pu_pv);
        let tv_tvv = tv.dot(tvv);

        let d_hdu = 2.0
            * (tu_pu_pv * (tuu_pu_pv - tu_tu) * tv_tv
                + tv_pu_pv * tu_tv * tu_tu
                - tu_tuu * tv_pu_pv * tv_pu_pv);
        let d_hdv = 2.0
            * (tu_pu_pv * tu_tv * tv_tv + tv_tvv * tu_pu_pv * tu_pu_pv
                - tv_pu_pv * (tvv_pu_pv + tv_tv) * tu_tu);

        //-----------------------------
        // Calculate dF/du and dF/dv.
        //-----------------------------
        let nor_pu_pv = nor.dot(pu_pv);
        let nu_pu_pv = nu.dot(pu_pv);
        let nor_tv = nor.dot(tv);

        let a = 0.5 * pu_pv.length_squared();
        let b = -nor_pu_pv;
        let bb = b * b;
        let d_adu = -tu_pu_pv;
        let d_bdu = -nu_pu_pv;
        let d_adv = tv_pu_pv;
        let d_bdv = -nor_tv;

        //---------------------------------------
        // F(u,v) = Pu - (A(u,v)/B(u,v))*Nor(u)
        //----------------------------------------
        if bb < GP_RESOLUTION {
            *v1 = tu.normalize_or_zero() + tv.normalize_or_zero();
            *v1 = 0.5 * tu.length_squared() * *v1;
        } else {
            let d_fdu = tu - (d_adu / b - d_bdu * a / bb) * nor - (a / b) * nu;
            let d_fdv = (-d_adv / b + d_bdv * a / bb) * nor;

            if d_hdv.abs() > GP_RESOLUTION {
                *v1 = d_fdu + d_fdv * (-d_hdu / d_hdv);
            } else {
                *v1 = tu;
            }
        }
        if n == 1 {
            return;
        }
    }

    /// OCCT Extension(U, U1, U2, Dist, T) (L1214-1324) — the current point
    /// on the extensions by tangence of the curve.
    fn extension(
        &self,
        u: f64,
        u1: &mut f64,
        u2: &mut f64,
        dist: &mut f64,
        t: &mut DVec2,
    ) -> DVec2 {
        let mut p_ref = PointOnBis::new();
        let mut p;
        let mut p1;
        let mut d_u = 0.0;
        let mut extension_tangent = false;
        let mut tang;

        if self.my_polygon.length() == 0 {
            //---------------------------------------------
            // Empty Polygon => segment (pointStart,pointEnd)
            //---------------------------------------------
            d_u = u - self.start_intervals.first().copied().unwrap_or(0.0);
            p = self.point_start;
            p1 = self.point_end;
            *u1 = self.curve(1).last_parameter();
            if self.sign1 == self.sign2 {
                *u2 = self.curve(2).first_parameter();
            } else {
                *u2 = self.curve(2).last_parameter();
            }
            tang = DVec2::new(p1.x - p.x, p1.y - p.y);
        } else if u < self.my_polygon.first().param_on_bis() {
            p_ref = *self.my_polygon.first();
            p = p_ref.point();
            d_u = u - p_ref.param_on_bis();
            if self.extension_start {
                //------------------------------------------------------------
                // extension = segment (pointstart, first point of the
                // polygon.)
                //------------------------------------------------------------
                p1 = self.point_start;
                *u1 = self.curve(1).first_parameter();
                if self.sign1 == self.sign2 {
                    *u2 = self.curve(2).last_parameter();
                } else {
                    *u2 = self.curve(2).first_parameter();
                }
                tang = DVec2::new(p.x - p1.x, p.y - p1.y);
            } else {
                extension_tangent = true;
                tang = DVec2::ZERO;
                p1 = DVec2::ZERO;
            }
        } else if u > self.my_polygon.last().param_on_bis() {
            p_ref = *self.my_polygon.last();
            p = p_ref.point();
            d_u = u - p_ref.param_on_bis();
            if self.extension_end {
                //------------------------------------------------------------
                // extension = segment (last point of the polygon.pointEnd)
                //------------------------------------------------------------
                p1 = self.point_end;
                *u1 = self.curve(1).last_parameter();
                if self.sign1 == self.sign2 {
                    *u2 = self.curve(2).last_parameter();
                } else {
                    *u2 = self.curve(2).first_parameter();
                }
                tang = DVec2::new(p1.x - p.x, p1.y - p.y);
            } else {
                extension_tangent = true;
                tang = DVec2::ZERO;
                p1 = DVec2::ZERO;
            }
        } else {
            // Neither extension nor empty polygon (defensive; OCCT leaves
            // the locals default-initialized).
            p = DVec2::ZERO;
            p1 = DVec2::ZERO;
            tang = DVec2::ZERO;
        }

        if extension_tangent {
            //-----------------------------------------------------------
            // If the la curve has no a extension, it is extended by
            // tangency
            //------------------------------------------------------------
            *u1 = p_ref.param_on_c1();
            *u2 = p_ref.param_on_c2();
            let p2 = self.curve(2).value(*u2);
            let (p1_t, t1) = self.curve(1).d1(*u1);
            p1 = p1_t;
            tang = DVec2::new(
                2.0 * p.x - p1.x - p2.x,
                2.0 * p.y - p1.y - p2.y,
            );
            if tang.length() < CONFUSION {
                tang = t1;
            }
            if t1.dot(tang) < 0.0 {
                tang = -tang;
            }
        }

        *t = tang.normalize_or_zero();
        let p_bis = DVec2::new(p.x + d_u * t.x, p.y + d_u * t.y);
        *dist = p1.distance(p_bis);
        p_bis
    }

    /// OCCT SupLastParameter() (L1485-1510).
    fn sup_last_parameter(&mut self) {
        self.end_intervals.push(self.curve(1).last_parameter());
        //-------------------------------------------------------------------
        // Calculate parameter on curve1 associated to one or the other of
        // the extremities of curve2 following the values of sign1 and sign2.
        // the bissectrice is limited by the obtained parameters.
        //-------------------------------------------------------------------
        let mut u_on_c1 = 0.0;
        let mut dist = 0.0;
        let u_on_c2 = if self.sign1 == self.sign2 {
            self.curve(2).first_parameter()
        } else {
            self.curve(2).last_parameter()
        };
        let ya_sol = point_by_int(
            &self.curve(2),
            &self.curve(1),
            self.sign2,
            self.sign1,
            u_on_c2,
            &mut u_on_c1,
            &mut dist,
        );
        if ya_sol
            && u_on_c1 > self.start_intervals.first().copied().unwrap_or(0.0)
            && u_on_c1 < self.end_intervals.last().copied().unwrap_or(0.0)
        {
            self.end_intervals[0] = u_on_c1;
        }
    }

    /// OCCT Curve(I) const (L1514-1528) — the curve handle of index I.
    ///
    /// Architecture difference: the OCCT field is a handle (possibly null
    /// before Perform); a null deref is mirrored by the expect.
    pub fn curve(&self, i: i32) -> Arc<dyn BisectorCurve> {
        match i {
            1 => self.curve1.clone().expect("curve1"),
            2 => self.curve2.clone().expect("curve2"),
            _ => panic!("Standard_OutOfRange: Bisector_BisecCC::Curve"),
        }
    }

    /// OCCT LinkBisCurve(U) (L1532-1535).
    pub fn link_bis_curve(&self, u: f64) -> f64 {
        u - self.shift_parameter
    }

    /// OCCT LinkCurveBis(U) (L1539-1542).
    pub fn link_curve_bis(&self, u: f64) -> f64 {
        u + self.shift_parameter
    }

    /// OCCT Polygon() const (L1559-1562).
    pub fn polygon(&self) -> &PolyBis {
        &self.my_polygon
    }

    /// OCCT Parameter(P) (L1566-1584).
    pub fn parameter(&self, p: DVec2) -> f64 {
        let u_on_curve;

        if pnt2d_equal(p, self.value(self.first_parameter()), CONFUSION) {
            u_on_curve = self.first_parameter();
        } else if pnt2d_equal(p, self.value(self.last_parameter()), CONFUSION) {
            u_on_curve = self.last_parameter();
        } else {
            let mut u_proj = 0.0;
            proj_on_curve(p, &self.curve(1), &mut u_proj);
            u_on_curve = u_proj;
        }

        u_on_curve
    }

    /// OCCT Dump(Deep, Offset) (L1589-1606).
    pub fn dump(&self, _deep: i32, _offset: i32) {
        println!("Bisector_BisecCC :");
        println!("Sign1  :{}", self.sign1);
        println!("Sign2  :{}", self.sign2);

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

    // -------------------------------------------------------------------
    // Private setters (the OCCT private Set-Curves / Set-Sign / ...).
    // -------------------------------------------------------------------

    /// OCCT Curve(const int I, const handle(C)& C) (L1610-1624) — setter.
    fn set_curve(&mut self, i: i32, c: Arc<dyn BisectorCurve>) {
        match i {
            1 => self.curve1 = Some(c),
            2 => self.curve2 = Some(c),
            _ => panic!("Standard_OutOfRange: Bisector_BisecCC::Curve"),
        }
    }

    /// OCCT Sign(const int I, const double S) (L1628-1642) — setter.
    fn set_sign(&mut self, i: i32, s: f64) {
        match i {
            1 => self.sign1 = s,
            2 => self.sign2 = s,
            _ => panic!("Standard_OutOfRange: Bisector_BisecCC::Sign"),
        }
    }

    /// OCCT Polygon(const Bisector_PolyBis& P) (L1646-1649) — setter.
    fn set_polygon(&mut self, p: PolyBis) {
        self.my_polygon = p;
    }

    /// OCCT DistMax(const double D) (L1653-1656) — setter.
    fn set_dist_max(&mut self, d: f64) {
        self.dist_max = d;
    }

    /// OCCT IsConvex(const int I, const bool IsConvex) (L1660-1674) —
    /// setter.
    fn set_is_convex(&mut self, i: i32, is_convex: bool) {
        match i {
            1 => self.is_convex1 = is_convex,
            2 => self.is_convex2 = is_convex,
            _ => panic!("Standard_OutOfRange: Bisector_BisecCC::IsConvex"),
        }
    }

    /// OCCT IsEmpty(const bool IsEmpty) (L1678-1681) — setter.
    fn set_is_empty(&mut self, is_empty: bool) {
        self.is_empty = is_empty;
    }

    /// OCCT ExtensionStart(const bool) (L1685-1688) — setter.
    fn set_extension_start(&mut self, extension_start: bool) {
        self.extension_start = extension_start;
    }

    /// OCCT ExtensionEnd(const bool) (L1692-1695) — setter.
    fn set_extension_end(&mut self, extension_end: bool) {
        self.extension_end = extension_end;
    }

    /// OCCT PointStart(const gp_Pnt2d&) (L1699-1702) — setter.
    fn set_point_start(&mut self, point: DVec2) {
        self.point_start = point;
    }

    /// OCCT PointEnd(const gp_Pnt2d&) (L1706-1709) — setter.
    fn set_point_end(&mut self, point: DVec2) {
        self.point_end = point;
    }

    /// OCCT StartIntervals(const Sequence&) (L1713-1716) — setter.
    fn set_start_intervals(&mut self, start_intervals: Vec<f64>) {
        self.start_intervals = start_intervals;
    }

    /// OCCT EndIntervals(const Sequence&) (L1720-1723) — setter.
    fn set_end_intervals(&mut self, end_intervals: Vec<f64>) {
        self.end_intervals = end_intervals;
    }

    /// OCCT FirstParameter(const double U) (L1727-1730) — setter, appends to
    /// startIntervals.
    fn append_first_parameter(&mut self, u: f64) {
        self.start_intervals.push(u);
    }

    /// OCCT LastParameter(const double U) (L1734-1737) — setter, appends to
    /// endIntervals.
    fn append_last_parameter(&mut self, u: f64) {
        self.end_intervals.push(u);
    }

    /// OCCT SearchBound(U1, U2) (L1741-1781).
    fn search_bound(&self, u1: f64, u2: f64) -> f64 {
        let mut uc1 = 0.0;
        let mut uc2 = 0.0;
        let tol_pnt = CONFUSION;
        let tol_par = PCONFUSION;
        let mut u11 = u1;
        let mut u22 = u2;
        let mut dist1 = 0.0;
        let mut dist2 = 0.0;
        let mut dist_mid = 0.0;
        let mut p_bis_prec = self.value_by_int(u11, &mut uc1, &mut uc2, &mut dist1);
        let mut p_bis = self.value_by_int(u22, &mut uc1, &mut uc2, &mut dist2);

        while (u22 - u11) > tol_par
            || ((dist1 < INFINITE_VALUE
                && dist2 < INFINITE_VALUE
                && !pnt2d_equal(p_bis, p_bis_prec, tol_pnt)))
        {
            p_bis_prec = p_bis;
            let u_mid = 0.5 * (u22 + u11);
            p_bis = self.value_by_int(u_mid, &mut uc1, &mut uc2, &mut dist_mid);
            if (dist1 < INFINITE_VALUE) == (dist_mid < INFINITE_VALUE) {
                u11 = u_mid;
                dist1 = dist_mid;
            } else {
                u22 = u_mid;
                dist2 = dist_mid;
            }
        }
        p_bis = self.value_by_int(u11, &mut uc1, &mut uc2, &mut dist1);
        let u_mid = if dist1 < INFINITE_VALUE { u11 } else { u22 };
        u_mid
    }

    /// OCCT ComputePointEnd() (L1879-1927).
    fn compute_point_end(&mut self) {
        let curve1 = self.curve(1);
        let curve2 = self.curve(2);
        let u1 = curve1.first_parameter();
        let u2 = if self.sign1 == self.sign2 {
            curve2.last_parameter()
        } else {
            curve2.first_parameter()
        };
        let k1 = curvature(&curve1, u1, CONFUSION);
        let k2 = curvature(&curve2, u2, CONFUSION);
        let kc;
        if !self.is_convex1 && !self.is_convex2 {
            kc = k1.min(k2);
        } else if !self.is_convex1 {
            kc = k1;
        } else {
            kc = k2;
        }

        let (pf, tf_raw) = curve1.d1(u1);
        let tf = tf_raw.normalize_or_zero();
        let rc = if kc != 0.0 {
            (1.0 / kc).abs()
        } else {
            INFINITE_VALUE
        };
        self.point_end = DVec2::new(
            pf.x - self.sign1 * rc * tf.y,
            pf.y + self.sign1 * rc * tf.x,
        );
    }
}

// ---------------------------------------------------------------------------
// Bisector_Curve trait implementation for BisectorBisecCC
// (OCCT: class Bisector_BisecCC : public Bisector_Curve).
// ---------------------------------------------------------------------------

impl BisectorCurve for BisectorBisecCC {
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
        BisectorBisecCC::first_parameter(self)
    }

    fn last_parameter(&self) -> f64 {
        BisectorBisecCC::last_parameter(self)
    }

    fn is_closed(&self) -> bool {
        BisectorBisecCC::is_closed(self)
    }

    fn is_periodic(&self) -> bool {
        false
    }

    fn continuity(&self) -> GeomAbsShape {
        BisectorBisecCC::continuity(self)
    }

    fn reversed_parameter(&self, u: f64) -> f64 {
        BisectorBisecCC::reversed_parameter(self, u)
    }

    fn reverse(&mut self) {
        BisectorBisecCC::reverse(self)
    }

    fn is_cn(&self, n: i32) -> bool {
        BisectorBisecCC::is_cn(self, n)
    }

    fn transform(&mut self, t: &glam::DAffine2) {
        BisectorBisecCC::transform(self, t)
    }

    fn copy_curve(&self) -> Arc<dyn BisectorCurve> {
        Arc::new(self.copy_bisec_cc())
    }

    fn parameter(&self, p: DVec2) -> f64 {
        BisectorBisecCC::parameter(self, p)
    }

    fn is_extend_at_start(&self) -> bool {
        BisectorBisecCC::is_extend_at_start(self)
    }

    fn is_extend_at_end(&self) -> bool {
        BisectorBisecCC::is_extend_at_end(self)
    }

    fn nb_intervals(&self) -> i32 {
        BisectorBisecCC::nb_intervals(self)
    }

    fn interval_first(&self, index: i32) -> f64 {
        BisectorBisecCC::interval_first(self, index)
    }

    fn interval_last(&self, index: i32) -> f64 {
        BisectorBisecCC::interval_last(self, index)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::BisecCC
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn as_any_arc(self: Arc<Self>) -> Arc<dyn std::any::Any + Send + Sync> {
        self
    }
}

/// `down_cast<Bisector_BisecCC>(handle)` (reconciled mat2d surface).
pub fn cast_from(h: Arc<dyn BisectorCurve>) -> Option<Arc<BisectorBisecCC>> {
    h.as_any_arc().downcast::<BisectorBisecCC>().ok()
}
