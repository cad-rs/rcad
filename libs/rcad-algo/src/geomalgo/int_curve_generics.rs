//! OCCT IntCurve generic algorithm instantiations (TKGeomAlgo IntCurve
//! package) — the `TheCurve`/`TheCurveTool` template parameter pair maps to
//! the [`ProjPCurveTool`] trait.
//!
//! 1:1 translations:
//! - `IntCurve_Polygon2dGen.gxx` (L38-268) + `.lxx` (L20-93) — sampled 2D
//!   polygon over a curve with deflection-driven refinement.
//! - `IntCurve_DistBetweenPCurvesGen.gxx` (L33-105) — the
//!   math_FunctionSetWithDerivatives used by ExactIntersectionPoint.
//! - `IntCurve_ExactIntersectionPoint.gxx` (L25-271) — Newton solve for the
//!   exact intersection of two polygons' underlying curves, with the
//!   bound-widening retry loops kept verbatim.
//!
//! `IntCurve_IntConicCurveGen.gxx` / `IntCurve_UserIntConicCurveGen.gxx`
//! are deferred: their `Perform(IConicTool, ...)` engine bodies bind to
//! HLRBRep_Curve/HLRBRep_CurveTool (TKHLR Stage 3a) and the
//! IntImpParGen_Intersector instantiation (Stage 2a-2).

use glam::DVec2;
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

/// OCCT `#define MAJORATION_DEFLECTION 1.5` (Polygon2dGen.gxx L27).
const MAJORATION_DEFLECTION: f64 = 1.5;

/// OCCT `TheCurveTool` static tool for the IntCurve generics: curve
/// evaluation and parameter epsilon.  (The HLRBRep instantiation is
/// HLRBRep_CurveTool; Geom2dInt uses IntCurveCurveGen's own tool.)
pub trait ProjPCurveTool {
    /// OCCT `TheCurve` type.
    type Curve: ?Sized;
    /// OCCT TheCurveTool::Value(C, U) — point at parameter U.
    fn value(c: &Self::Curve, u: f64) -> DVec2;
    /// OCCT TheCurveTool::D1(C, U, P, T) — point and tangent.
    fn d1(c: &Self::Curve, u: f64) -> (DVec2, DVec2);
    /// OCCT TheCurveTool::EpsX(C) — parameter resolution.
    fn eps_x(c: &Self::Curve) -> f64;
}

/// OCCT IntCurve_Polygon2dGen — a polygon over a parametric curve: samples
/// at constant parameter steps, then refines where a mid-point deviates
/// from its chord by more than half the current deflection.
///
/// Arrays are 1-based in OCCT (`NCollection_Array1`, indices stored in
/// `TheIndex` are 1-based point indices); the Vec members here are 0-based
/// and every access subtracts 1.
#[derive(Debug, Clone)]
pub struct Polygon2dGen<C: ?Sized> {
    the_pnts: Vec<DVec2>,
    the_params: Vec<f64>,
    the_index: Vec<i32>,
    the_max_nb_points: usize,
    nb_pnt_in: usize,
    the_deflection: f64,
    binf: f64,
    bsup: f64,
    closed_polygon: bool,
    my_box: BndBox2d,
    _curve: std::marker::PhantomData<fn(&C)>,
}

impl<C: ?Sized> Default for Polygon2dGen<C> {
    fn default() -> Self {
        panic!("IntCurve_Polygon2dGen has no default constructor in OCCT")
    }
}

fn calcul_region(x: f64, y: f64, x1: f64, x2: f64, y1: f64, y2: f64) -> i32 {
    // OCCT IntCurve_Polygon2dGen::CalculRegion (lxx L58-93).
    let mut r;
    if x < x1 {
        r = 1;
    } else if x > x2 {
        r = 2;
    } else {
        r = 0;
    }
    if y < y1 {
        r |= 4;
    } else if y > y2 {
        r |= 8;
    }
    r
}

impl<C: ?Sized> Polygon2dGen<C> {
    /// OCCT DeflectionOverEstimation() — lxx L20-23.
    pub fn deflection_over_estimation(&self) -> f64 {
        self.the_deflection
    }

    /// OCCT SetDeflectionOverEstimation(x) — lxx L25-29.
    pub fn set_deflection_over_estimation(&mut self, x: f64) {
        self.the_deflection = x;
        self.my_box.enlarge(self.the_deflection);
    }

    /// OCCT Closed(flag) — lxx L32-35.
    pub fn set_closed(&mut self, flag: bool) {
        self.closed_polygon = flag;
    }

    /// OCCT NbSegments() — lxx L38-41.
    pub fn nb_segments(&self) -> i32 {
        if self.closed_polygon {
            self.nb_pnt_in as i32
        } else {
            self.nb_pnt_in as i32 - 1
        }
    }

    /// OCCT InfParameter() — lxx L44-47.
    pub fn inf_parameter(&self) -> f64 {
        self.the_params[self.the_index[0] as usize - 1]
    }

    /// OCCT SupParameter() — lxx L50-53.
    pub fn sup_parameter(&self) -> f64 {
        self.the_params[self.the_index[self.nb_pnt_in - 1] as usize - 1]
    }

    /// OCCT Bounding() (Intf_Polygon2d-style accessor used by callers).
    pub fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }
}

impl<C: ?Sized> Polygon2dGen<C> {
    /// OCCT IntCurve_Polygon2dGen(C, tNbPts, D, Tol) — gxx L38-117: uniform
    /// sampling and initial deflection estimate.
    pub fn new<T: ProjPCurveTool<Curve = C>>(
        c: &C,
        t_nb_pts: i32,
        d_first: f64,
        d_last: f64,
        tol: f64,
    ) -> Self {
        let size = if t_nb_pts < 3 { 6 } else { (t_nb_pts + t_nb_pts) as usize };
        let mut poly = Polygon2dGen {
            the_pnts: vec![DVec2::ZERO; size],
            the_params: vec![0.0; size],
            the_index: vec![0; size],
            the_max_nb_points: 0,
            nb_pnt_in: 0,
            the_deflection: 0.0,
            binf: 0.0,
            bsup: 0.0,
            closed_polygon: false,
            my_box: BndBox2d::new(),
            _curve: std::marker::PhantomData,
        };

        let nb_pts = if t_nb_pts < 3 { 3usize } else { t_nb_pts as usize };
        poly.the_max_nb_points = nb_pts + nb_pts;
        poly.nb_pnt_in = nb_pts;

        // Initialization of the breaking with constant parameter step.
        poly.binf = d_first;
        poly.bsup = d_last;
        let mut u = poly.binf;
        let u1 = poly.bsup;
        let du = (u1 - u) / (nb_pts as f64 - 1.0);
        let mut i = 1usize;

        loop {
            let p = T::value(c, u);
            poly.my_box.add_point(p);
            poly.the_index[i - 1] = i as i32;
            poly.the_pnts[i - 1] = p;
            poly.the_params[i - 1] = u;
            u += du;
            i += 1;
            if i > nb_pts {
                break;
            }
        }

        // Calculate a maximal deflection (gxx L82: min(1e-9, Tol/100)).
        poly.the_deflection = 0.000_000_001_f64.min(tol / 100.0);
        let mut i = 1usize;
        let mut u = d_first + du * 0.5;

        loop {
            let pm = T::value(c, u);
            let p1 = poly.the_pnts[i - 1];
            let p2 = poly.the_pnts[i];

            u += du;
            i += 1;

            let mut dx = p1.x - p2.x;
            if dx < 0.0 {
                dx = -dx;
            }
            let mut dy = p1.y - p2.y;
            if dy < 0.0 {
                dy = -dy;
            }
            if dx + dy > 1e-12 {
                // gp_Lin2d(P1, gp_Dir2d(gp_Vec2d(P1, P2))).Distance(Pm).
                let dir = DVec2::new(p2.x - p1.x, p2.y - p1.y).normalize_or_zero();
                let t = (DVec2::new(pm.x - p1.x, pm.y - p1.y)
                    - dir * (DVec2::new(pm.x - p1.x, pm.y - p1.y).dot(dir)))
                .length();
                if t > poly.the_deflection {
                    poly.the_deflection = t;
                }
            }
            if i >= nb_pts {
                break;
            }
        }

        poly.my_box.enlarge(poly.the_deflection * MAJORATION_DEFLECTION);
        poly.closed_polygon = false;
        poly
    }

    /// OCCT ComputeWithBox(C, BoxOtherPolygon) — gxx L121-266.
    pub fn compute_with_box<T: ProjPCurveTool<Curve = C>>(
        &mut self,
        c: &C,
        box_other_polygon: &BndBox2d,
    ) {
        if self.my_box.is_out_box(box_other_polygon) {
            self.nb_pnt_in = 2;
            self.my_box.set_void();
        } else {
            let (mut bx0, mut by0, mut bx1, mut by1) = box_other_polygon
                .get()
                .unwrap_or((-1.0, -1.0, 1.0, 1.0));

            bx0 -= self.the_deflection;
            by0 -= self.the_deflection;
            bx1 += self.the_deflection;
            by1 += self.the_deflection;
            let mut max_index_used = 1usize;
            let mut nbp: usize = 0;

            let mut x = self.the_pnts[self.the_index[0] as usize - 1].x;
            let mut y = self.the_pnts[self.the_index[0] as usize - 1].y;

            let mut rprec = calcul_region(x, y, bx0, bx1, by0, by1);
            for i in 2..=self.nb_pnt_in {
                let p2d = self.the_pnts[self.the_index[i - 1] as usize - 1];
                let ri = calcul_region(p2d.x, p2d.y, bx0, bx1, by0, by1);
                if (ri & rprec) == 0 {
                    if nbp != 0 {
                        if self.the_index[nbp - 1] != self.the_index[i - 2] {
                            nbp += 1;
                            self.the_index[nbp - 1] = self.the_index[i - 2];
                        }
                    } else {
                        nbp += 1;
                        self.the_index[nbp - 1] = self.the_index[i - 2];
                    }
                    nbp += 1;
                    self.the_index[nbp - 1] = self.the_index[i - 1];
                    if self.the_index[i - 1] as usize > max_index_used {
                        max_index_used = self.the_index[i - 1] as usize;
                    }
                    rprec = ri;
                }
                rprec = ri;
            }
            if nbp == 1 {
                self.nb_pnt_in = 2;
                self.my_box.set_void();
            } else {
                self.my_box.set_void();
                if nbp != 0 {
                    let p = self.the_pnts[self.the_index[0] as usize - 1];
                    self.my_box.add_point(p);
                }
                let mut ratio_deflection;
                let mut nbpassagedeflection = 0;
                let mut i;
                let mut nbp = nbp;
                loop {
                    nbpassagedeflection += 1;
                    let mut new_deflection = self.the_deflection;
                    i = 2usize;
                    while i <= nbp {
                        let i_i = self.the_index[i - 1] as usize;
                        let i_im1 = self.the_index[i - 2] as usize;
                        let pi = self.the_pnts[i_i - 1];
                        let pim1 = self.the_pnts[i_im1 - 1];
                        self.my_box.add_point(pi);
                        let reg_i = calcul_region(pi.x, pi.y, bx0, bx1, by0, by1);
                        let reg_im1 = calcul_region(pim1.x, pim1.y, bx0, bx1, by0, by1);
                        if (reg_i & reg_im1) == 0 {
                            let u = 0.5 * (self.the_params[i_i - 1] + self.the_params[i_im1 - 1]);
                            let pm = T::value(c, u);
                            let mut dx = pim1.x - pi.x;
                            if dx < 0.0 {
                                dx = -dx;
                            }
                            let mut dy = pim1.y - pi.y;
                            if dy < 0.0 {
                                dy = -dy;
                            }
                            // OCCT gxx L211: dx, dy, t are declared together
                            // with t = 0 BEFORE the branch.
                            let mut t = 0.0f64;
                            if dx + dy > 1e-12 {
                                let dir =
                                    DVec2::new(pi.x - pim1.x, pi.y - pim1.y).normalize_or_zero();
                                t = (DVec2::new(pm.x - pim1.x, pm.y - pim1.y)
                                    - dir * DVec2::new(pm.x - pim1.x, pm.y - pim1.y).dot(dir))
                                .length();
                                if (max_index_used < self.the_max_nb_points - 1)
                                    && (t > self.the_deflection * 0.5)
                                {
                                    let p1 = pim1;
                                    nbp += 1;
                                    let mut j = nbp;
                                    while j >= i + 1 {
                                        self.the_index[j - 1] = self.the_index[j - 2];
                                        j -= 1;
                                    }
                                    max_index_used += 1;
                                    self.the_index[i - 1] = max_index_used as i32;
                                    self.the_pnts[max_index_used - 1] = pm;
                                    self.the_params[max_index_used - 1] = u;

                                    let u1m = 0.5
                                        * (u + self.the_params[self.the_index[i - 2] as usize - 1]);
                                    let p1m = T::value(c, u1m);
                                    let dir1m =
                                        DVec2::new(pm.x - p1.x, pm.y - p1.y).normalize_or_zero();
                                    let t = (DVec2::new(p1m.x - p1.x, p1m.y - p1.y)
                                        - dir1m
                                            * DVec2::new(p1m.x - p1.x, p1m.y - p1.y).dot(dir1m))
                                    .length();
                                    i -= 1;
                                    let _ = t;
                                }
                            } else if t > new_deflection {
                                new_deflection = t;
                            }
                        }
                        i += 1;
                    }
                    if new_deflection != 0.0 {
                        ratio_deflection = self.the_deflection / new_deflection;
                    } else {
                        ratio_deflection = 10.0;
                    }
                    self.the_deflection = new_deflection;
                    self.nb_pnt_in = nbp;
                    let _ = new_deflection;
                    if !((ratio_deflection < 3.0)
                        && (nbpassagedeflection < 3)
                        && (max_index_used < self.the_max_nb_points - 2))
                    {
                        break;
                    }
                }
            }

            self.the_deflection *= MAJORATION_DEFLECTION;
            self.my_box.enlarge(self.the_deflection);
        }
        self.closed_polygon = false;
        // OCCT calls Dump() here (a no-op with static debug == 0).
    }

    /// OCCT AutoIntersectionIsPossible() — gxx L268-281.
    pub fn auto_intersection_is_possible(&self) -> bool {
        let p1 = self.the_pnts[self.the_index[0] as usize - 1];
        let p2 = self.the_pnts[self.the_index[1] as usize - 1];
        let v_ref = DVec2::new(p2.x - p1.x, p2.y - p1.y);
        for i in 3..=self.nb_pnt_in {
            let a = self.the_pnts[self.the_index[i - 2] as usize - 1];
            let b = self.the_pnts[self.the_index[i - 1] as usize - 1];
            let v = DVec2::new(b.x - a.x, b.y - a.y);
            if v.dot(v_ref) < 0.0 {
                return true;
            }
        }
        false
    }

    /// OCCT ApproxParamOnCurve(Aindex, TheParamOnLine) — gxx L285-310.
    pub fn approx_param_on_curve(&self, aindex: i32, the_param_on_line: f64) -> f64 {
        let mut index = aindex;
        let mut param_on_line = the_param_on_line;
        if index > self.nb_pnt_in as i32 {
            // OCCT prints "OutOfRange Polygon2d::ApproxParamOnCurve".
        }
        if (index == self.nb_pnt_in as i32) && (param_on_line == 0.0) {
            index -= 1;
            param_on_line = 1.0;
        }
        if index == 0 {
            index = 1;
            param_on_line = 0.0;
        }
        let indexp1 = self.the_index[index as usize];
        let index0 = self.the_index[index as usize - 1];

        let du = self.the_params[indexp1 as usize - 1] - self.the_params[index0 as usize - 1];
        self.the_params[index0 as usize - 1] + param_on_line * du
    }

    /// OCCT Segment(theIndex, theBegin, theEnd) — gxx L370-381.
    pub fn segment(&self, the_index: i32) -> (DVec2, DVec2) {
        let mut ind = the_index;
        let the_begin = self.the_pnts[self.the_index[the_index as usize - 1] as usize - 1];
        if the_index >= self.nb_pnt_in as i32 {
            if !self.closed_polygon {
                panic!("IntCurve_Polygon2dGen::Segment!");
            }
            ind = 0;
        }
        let the_end = self.the_pnts[self.the_index[ind as usize] as usize - 1];
        (the_begin, the_end)
    }
}

/// OCCT IntCurve_DistBetweenPCurvesGen — math_FunctionSetWithDerivatives
/// computing the vector difference of two parametrised curves
/// (DistBetweenPCurvesGen.gxx L33-105).  OCCT stores `void*` curve
/// pointers; the Rust struct borrows nothing and implements the function
/// set generically over the tool.
pub struct DistBetweenPCurvesGen<'a, C: ?Sized, T>
where
    T: ProjPCurveTool<Curve = C>,
{
    curve1: &'a C,
    curve2: &'a C,
    _tool: std::marker::PhantomData<fn(&C, &C) -> T>,
}

impl<'a, C: ?Sized, T> DistBetweenPCurvesGen<'a, C, T>
where
    T: ProjPCurveTool<Curve = C>,
{
    /// OCCT IntCurve_DistBetweenPCurvesGen(C1, C2) — gxx L33-38.
    pub fn new(curve1: &'a C, curve2: &'a C) -> Self {
        DistBetweenPCurvesGen {
            curve1,
            curve2,
            _tool: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, T> FunctionSetWithDerivatives for DistBetweenPCurvesGen<'_, C, T>
where
    T: ProjPCurveTool<Curve = C>,
{
    /// OCCT NbVariables() — gxx L42-45.
    fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT NbEquations() — gxx L49-52.
    fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, F) — gxx L56-64.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let p1 = T::value(self.curve1, x[0]);
        let p2 = T::value(self.curve2, x[1]);
        f[0] = p1.x - p2.x;
        f[1] = p1.y - p2.y;
        true
    }

    /// OCCT Derivatives(X, D) — gxx L68-81.
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        let (_p, t1) = T::d1(self.curve1, x[0]);
        df[0][0] = t1.x;
        df[1][0] = t1.y;
        let (_p, t2) = T::d1(self.curve2, x[1]);
        df[0][1] = -t2.x;
        df[1][1] = -t2.y;
        true
    }

    /// OCCT Values(X, F, D) — gxx L85-103.
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        let (p1, t1) = T::d1(self.curve1, x[0]);
        df[0][0] = t1.x;
        df[1][0] = t1.y;
        let (p2, t2) = T::d1(self.curve2, x[1]);
        df[0][1] = -t2.x;
        df[1][1] = -t2.y;
        f[0] = p1.x - p2.x;
        f[1] = p1.y - p2.y;
        true
    }
}

/// OCCT IntCurve_ExactIntersectionPoint — exact intersection of two curves
/// by Newton iteration seeded from polygon intersections, with the
/// bound-widening retry loops (ExactIntersectionPoint.gxx L25-271).
pub struct ExactIntersectionPoint<C: ?Sized, T>
where
    T: ProjPCurveTool<Curve = C>,
{
    done: bool,
    nbroots: i32,
    my_tol: f64,
    fct_dist: (),
    tolerance_vector: Vec<f64>,
    b_inf_vector: Vec<f64>,
    b_sup_vector: Vec<f64>,
    starting_point: Vec<f64>,
    root: Vec<f64>,
    an_error_occurred: bool,
    _tool: std::marker::PhantomData<fn(&C, &C) -> T>,
}

/// OCCT `IntCurve_ThePolygon2d` interface consumed by Perform — the
/// polygon accessors used by the retry loops.
pub trait Polygon2dLike {
    fn nb_segments(&self) -> i32;
    fn deflection_over_estimation(&self) -> f64;
    fn inf_parameter(&self) -> f64;
    fn sup_parameter(&self) -> f64;
    fn approx_param_on_curve(&self, index: i32, param_on_line: f64) -> f64;
}

impl<C: ?Sized> Polygon2dLike for Polygon2dGen<C> {
    fn nb_segments(&self) -> i32 {
        Polygon2dGen::nb_segments(self)
    }
    fn deflection_over_estimation(&self) -> f64 {
        Polygon2dGen::deflection_over_estimation(self)
    }
    fn inf_parameter(&self) -> f64 {
        Polygon2dGen::inf_parameter(self)
    }
    fn sup_parameter(&self) -> f64 {
        Polygon2dGen::sup_parameter(self)
    }
    fn approx_param_on_curve(&self, index: i32, param_on_line: f64) -> f64 {
        Polygon2dGen::approx_param_on_curve(self, index, param_on_line)
    }
}

impl<C: ?Sized, T> ExactIntersectionPoint<C, T>
where
    T: ProjPCurveTool<Curve = C>,
{
    /// OCCT IntCurve_ExactIntersectionPoint(C1, C2, Tol) — gxx L25-41.
    pub fn new(c1: &C, c2: &C, tol: f64) -> Self {
        ExactIntersectionPoint {
            done: false,
            nbroots: 0,
            my_tol: tol * tol,
            fct_dist: (),
            tolerance_vector: vec![T::eps_x(c1), T::eps_x(c2)],
            b_inf_vector: vec![0.0; 2],
            b_sup_vector: vec![0.0; 2],
            starting_point: vec![0.0; 2],
            root: vec![0.0; 2],
            an_error_occurred: false,
            _tool: std::marker::PhantomData,
        }
    }

    /// OCCT Perform(Poly1, Poly2, NumSegOn1, NumSegOn2, ParamOnSeg1,
    /// ParamOnSeg2) — gxx L45-200, including the bound-widening retry
    /// loops.
    pub fn perform<P1: Polygon2dLike, P2: Polygon2dLike>(
        &mut self,
        poly1: &P1,
        poly2: &P2,
        num_seg_on1: &mut i32,
        num_seg_on2: &mut i32,
        param_on_seg1: &mut f64,
        param_on_seg2: &mut f64,
        fct_dist: &mut dyn FunctionSetWithDerivatives,
    ) {
        // Search bounds: segment i-1 .. i+2 around the seed (gxx L52-79).
        if *num_seg_on1 >= poly1.nb_segments() && *param_on_seg1 == 0.0 {
            *num_seg_on1 -= 1;
            *param_on_seg1 = 1.0;
        }
        if *num_seg_on2 >= poly2.nb_segments() && *param_on_seg2 == 0.0 {
            *num_seg_on2 -= 1;
            *param_on_seg2 = 1.0;
        }
        if *num_seg_on1 <= 0 {
            *num_seg_on1 = 1;
            *param_on_seg1 = 0.0;
        }
        if *num_seg_on2 <= 0 {
            *num_seg_on2 = 1;
            *param_on_seg2 = 0.0;
        }

        self.starting_point[0] = poly1.approx_param_on_curve(*num_seg_on1, *param_on_seg1);
        if *num_seg_on1 <= 2 {
            self.b_inf_vector[0] = poly1.inf_parameter();
        } else {
            self.b_inf_vector[0] = poly1.approx_param_on_curve(*num_seg_on1 - 1, 0.0);
        }
        if *num_seg_on1 >= poly1.nb_segments() - 2 {
            self.b_sup_vector[0] = poly1.sup_parameter();
        } else {
            self.b_sup_vector[0] = poly1.approx_param_on_curve(*num_seg_on1 + 2, 0.0);
        }

        self.starting_point[1] = poly2.approx_param_on_curve(*num_seg_on2, *param_on_seg2);
        if *num_seg_on2 <= 2 {
            self.b_inf_vector[1] = poly2.inf_parameter();
        } else {
            self.b_inf_vector[1] = poly2.approx_param_on_curve(*num_seg_on2 - 1, 0.0);
        }
        if *num_seg_on2 >= poly2.nb_segments() - 2 {
            self.b_sup_vector[1] = poly2.sup_parameter();
        } else {
            self.b_sup_vector[1] = poly2.approx_param_on_curve(*num_seg_on2 + 2, 0.0);
        }

        self.math_perform(fct_dist);
        if self.nbroots == 0 {
            // OCCT L104-107: the deflection reads are discarded.
            let _ = poly1.deflection_over_estimation();
            let _ = poly2.deflection_over_estimation();
            {
                // The bounds on curve 1 may be too narrow.
                let mut diff = 1;
                let an_binf_vector = self.b_inf_vector[0];
                let an_bsup_vector = self.b_sup_vector[0];
                // Widen the bounds to the left (gxx L114-129).
                loop {
                    diff += 1;
                    if (*num_seg_on1 - diff) <= 1 {
                        self.b_inf_vector[0] = poly1.inf_parameter();
                        diff = 0;
                    } else {
                        self.b_inf_vector[0] = poly1.approx_param_on_curve(*num_seg_on1 - diff, 0.0);
                    }
                    self.math_perform(fct_dist);
                    // le 18 nov 97
                    if diff > 3 {
                        diff += *num_seg_on1 / 2;
                    }
                    if !(self.nbroots == 0 && diff != 0) {
                        break;
                    }
                }
                // Widen the bounds to the right (gxx L130-150).
                if self.nbroots == 0 {
                    self.b_inf_vector[0] = an_binf_vector;
                    let mut diff = 1;
                    loop {
                        diff += 1;
                        if (*num_seg_on1 + diff) >= (poly1.nb_segments() - 1) {
                            self.b_sup_vector[0] = poly1.sup_parameter();
                            diff = 0;
                        } else {
                            self.b_sup_vector[0] =
                                poly1.approx_param_on_curve(*num_seg_on1 + 1 + diff, 0.0);
                        }
                        self.math_perform(fct_dist);
                        // le 18 nov 97
                        if diff > 3 {
                            diff += 1 + (poly1.nb_segments() - *num_seg_on1) / 2;
                        }
                        if !(self.nbroots == 0 && diff != 0) {
                            break;
                        }
                    }
                }
                self.b_sup_vector[0] = an_bsup_vector;
            }

            if self.nbroots == 0 {
                // The bounds on curve 2 may be too narrow (gxx L154-198).
                let mut diff = 1;
                let an_binf_vector = self.b_inf_vector[1];
                let an_bsup_vector = self.b_sup_vector[1];
                loop {
                    diff += 1;
                    if (*num_seg_on2 - diff) <= 1 {
                        self.b_inf_vector[1] = poly2.inf_parameter();
                        diff = 0;
                    } else {
                        self.b_inf_vector[1] = poly2.approx_param_on_curve(*num_seg_on2 - diff, 0.0);
                    }
                    self.math_perform(fct_dist);
                    if diff > 3 {
                        diff += *num_seg_on2 / 2;
                    }
                    if !(self.nbroots == 0 && diff != 0) {
                        break;
                    }
                }
                if self.nbroots == 0 {
                    self.b_inf_vector[1] = an_binf_vector;
                    let mut diff = 1;
                    loop {
                        diff += 1;
                        if (*num_seg_on2 + diff) >= (poly2.nb_segments() - 1) {
                            self.b_sup_vector[1] = poly2.sup_parameter();
                            diff = 0;
                        } else {
                            self.b_sup_vector[1] =
                                poly2.approx_param_on_curve(*num_seg_on2 + 1 + diff, 0.0);
                        }
                        self.math_perform(fct_dist);
                        if diff > 3 {
                            diff += 1 + (poly2.nb_segments() - *num_seg_on2) / 2;
                        }
                        if !(self.nbroots == 0 && diff != 0) {
                            break;
                        }
                    }
                }
                self.b_sup_vector[1] = an_bsup_vector;
            }
        }
    }

    /// OCCT Perform(Uo, Vo, UInf, VInf, USup, VSup) — gxx L204-222.
    pub fn perform_bounds(
        &mut self,
        uo: f64,
        vo: f64,
        u_inf: f64,
        v_inf: f64,
        u_sup: f64,
        v_sup: f64,
        fct_dist: &mut dyn FunctionSetWithDerivatives,
    ) {
        self.done = true;
        self.b_inf_vector[0] = u_inf;
        self.b_inf_vector[1] = v_inf;
        self.b_sup_vector[0] = u_sup;
        self.b_sup_vector[1] = v_sup;
        self.starting_point[0] = uo;
        self.starting_point[1] = vo;
        self.math_perform(fct_dist);
    }

    /// OCCT NbRoots() — gxx L226-229.
    pub fn nb_roots(&self) -> i32 {
        self.nbroots
    }

    /// OCCT Roots(U, V) — gxx L233-237.
    pub fn roots(&self) -> (f64, f64) {
        (self.root[0], self.root[1])
    }

    /// OCCT MathPerform() — gxx L241-264.
    fn math_perform(&mut self, fct_dist: &mut dyn FunctionSetWithDerivatives) {
        let mut fct = FunctionSetRoot::new(fct_dist, &self.tolerance_vector, 60);
        fct.perform(
            fct_dist,
            &self.starting_point,
            &self.b_inf_vector,
            &self.b_sup_vector,
            false,
        );

        if fct.is_done() {
            self.root = fct.root();
            self.nbroots = 1;
            let mut xy = [0.0f64; 2];
            fct_dist.value(&self.root, &mut xy);
            let dist2 = xy[0] * xy[0] + xy[1] * xy[1];
            if dist2 > self.my_tol {
                self.nbroots = 0;
            }
        } else {
            self.an_error_occurred = true;
            self.nbroots = 0;
        }
    }

    /// OCCT AnErrorOccurred() — gxx L268-271.
    pub fn an_error_occurred(&self) -> bool {
        self.an_error_occurred
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Parabola y = x^2 as a TheCurve: value(u) = (u, u^2).
    struct ParaCurve;
    struct ParaTool;

    impl ProjPCurveTool for ParaTool {
        type Curve = ParaCurve;
        fn value(c: &ParaCurve, u: f64) -> DVec2 {
            let _ = c;
            DVec2::new(u, u * u)
        }
        fn d1(c: &ParaCurve, u: f64) -> (DVec2, DVec2) {
            let _ = c;
            (DVec2::new(u, u * u), DVec2::new(1.0, 2.0 * u))
        }
        fn eps_x(_c: &ParaCurve) -> f64 {
            1e-9
        }
    }

    /// OCCT Polygon2dGen ctor (gxx L38-117): uniform samples on [0,2] and a
    /// deflection consistent with the chord/mid-point deviation.
    #[test]
    fn polygon2d_gen_samples_and_deflection() {
        let c = ParaCurve;
        let poly = Polygon2dGen::<ParaCurve>::new::<ParaTool>(&c, 9, 0.0, 2.0, 1e-7);
        assert_eq!(poly.nb_segments(), 8);
        assert_eq!(poly.inf_parameter(), 0.0);
        assert_eq!(poly.sup_parameter(), 2.0);
        // TheDeflection >= min(1e-9, Tol/100) and the box is enlarged.
        assert!(poly.deflection_over_estimation() >= 1e-11);
        let (p0, p1) = poly.segment(1);
        assert_eq!(p0, DVec2::new(0.0, 0.0));
        assert_eq!(p1, DVec2::new(0.25, 0.0625));
        // ApproxParamOnCurve(1, 0.5) = midpoint parameter of segment 1.
        let u = poly.approx_param_on_curve(1, 0.5);
        assert!((u - 0.125).abs() < 1e-12);
    }

    /// OCCT ExactIntersectionPoint of two parabolas through the origin:
    /// C1: y = x, C2: y = -x intersect at (0, 0).
    #[test]
    fn exact_intersection_two_lines_as_curves() {
        struct LineCurve {
            a: f64,
            b: f64,
        };
        struct LineTool;
        impl ProjPCurveTool for LineTool {
            type Curve = LineCurve;
            fn value(c: &LineCurve, u: f64) -> DVec2 {
                DVec2::new(u, c.a * u + c.b)
            }
            fn d1(c: &LineCurve, u: f64) -> (DVec2, DVec2) {
                (Self::value(c, u), DVec2::new(1.0, c.a))
            }
            fn eps_x(_c: &LineCurve) -> f64 {
                1e-9
            }
        }

        let c1 = LineCurve { a: 1.0, b: 0.0 };
        let c2 = LineCurve { a: -1.0, b: 0.0 };
        let poly1 = Polygon2dGen::<LineCurve>::new::<LineTool>(&c1, 5, -1.0, 1.0, 1e-7);
        let poly2 = Polygon2dGen::<LineCurve>::new::<LineTool>(&c2, 5, -1.0, 1.0, 1e-7);

        let mut ex = ExactIntersectionPoint::<LineCurve, LineTool>::new(&c1, &c2, 1e-8);
        let mut fct = DistBetweenPCurvesGen::<LineCurve, LineTool>::new(&c1, &c2);
        let mut n1 = 2i32;
        let mut n2 = 2i32;
        let mut s1 = 0.5f64;
        let mut s2 = 0.5f64;
        ex.perform(&poly1, &poly2, &mut n1, &mut n2, &mut s1, &mut s2, &mut fct);
        assert_eq!(ex.nb_roots(), 1);
        let (u, v) = ex.roots();
        assert!((u.abs()) < 1e-5 && (v.abs()) < 1e-5, "u={u} v={v}");
        // On C1 the point is (u, u), on C2 the point is (v, -v); equal x.
        assert!((u - v).abs() < 1e-5);
        assert!(!ex.an_error_occurred());
    }
}

// ---------------------------------------------------------------------------
// OCCT IntCurve_Polygon2dGen : Intf_Polygon2d — the sampled polygon
// implements the Intf interference base-class interface (the hxx
// inheritance), consumed by Intf_InterferencePolygon2d.
// ---------------------------------------------------------------------------

impl<C: ?Sized> crate::geomalgo::intf_interference::IntfPolygon2d for Polygon2dGen<C> {
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }
    fn bounding_mut(&mut self) -> &mut BndBox2d {
        &mut self.my_box
    }
    fn closed(&self) -> bool {
        self.closed_polygon
    }
    fn deflection_over_estimation(&self) -> f64 {
        Polygon2dGen::deflection_over_estimation(self)
    }
    fn nb_segments(&self) -> i32 {
        Polygon2dGen::nb_segments(self)
    }
    fn segment(&self, the_index: i32) -> (DVec2, DVec2) {
        Polygon2dGen::segment(self, the_index)
    }
}
