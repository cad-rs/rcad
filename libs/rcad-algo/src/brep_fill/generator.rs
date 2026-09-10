//! OCCT BRepFill_Generator (TKBool/BRepFill) — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/BRepFill_Generator.cxx
//! (L61-1239) + BRepFill_Generator.hxx (L38-86) +
//! BRepFill_ThruSectionErrorStatus.hxx (L21-30).
//!
//! First consumer: BRepOffsetAPI_ThruSections (Stage 2e).  The LocOpe_
//! generator context was verified NOT to consume this class (Stage 3c).
//!
//! The file also carries the 1:1 translation of GeomFill_Generator.cxx
//! (L26-99) and the used part of GeomFill_Profiler.cxx (L116-287), because
//! rcad has no GeomFill module hosting them yet (reported gap).
//!
//! Architecture notes (Rust vs C++ data model):
//! - OCCT map key (TopTools_ShapeMapHasher::IsSame == TopoDS_Shape::IsSame ==
//!   same TShape) maps to `Shape::ptr_id()`.
//! - `BRep_Tool::Curve(E, L, f, l)` + `Transformed(loc)` maps to
//!   `BRep::edge(e).curve/.range` (rcad stores the located 3D curve; the
//!   generator edges carry identity locations).
//! - `TopExp::Vertices(E, V1, V2)` maps to the oriented (first, last) access
//!   of `TEdgeData` via the edge wrapper orientation.
//! - `BRep_Builder::MakeFace(F, Surf, Tol)` + `B.Add(Face, Wire)` maps to
//!   `BRep::add_tface_tol`; rcad needs the outer wire at construction, so the
//!   wire is built before the face (OCCT builds the face first and adds the
//!   wire afterwards).  Data flow is unchanged: every pcurve is bound to the
//!   face key after both exist, as in OCCT.
//! - `BRepLib::SameParameter(Shell)` maps to `BRep::same_parameter()`.

use glam::{DVec2, DVec3};
use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::core::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{
    BezierCurve2, BezierCurve3, Circle3, ConicalSurface, Curve2d, Curve3, CurveEval, CylindricalSurface,
    Line2d, Line3, Surface3, SurfaceEval, TrimmedCurve3, TrimmedSurface, BSplineCurve3, BSplineSurface,
};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

/// OCCT Precision::Confusion().
pub(super) const TOL_CONFUSION: f64 = CONFUSION;
/// OCCT Precision::PConfusion().
pub(super) const TOL_PCONFUSION: f64 = PCONFUSION;
/// OCCT Precision::Angular().
const TOL_ANGULAR: f64 = rcad_kernel::core::precision::ANGULAR;

/// OCCT BRepFill_ThruSectionErrorStatus (BRepFill_ThruSectionErrorStatus.hxx
/// L21-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFillThruSectionErrorStatus {
    Done,
    NotDone,
    NotSameTopology,
    ProfilesInconsistent,
    WrongUsage,
    Null3DCurve,
    Failed,
}

// =============================================================================
// Shape-map helpers (TopTools_ShapeMapHasher semantics)
// =============================================================================

/// OCCT TopTools_ShapeMapHasher — map key is the TShape identity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ShapeKey(pub u64);

pub(crate) fn shape_key(s: &Shape) -> ShapeKey {
    ShapeKey(s.ptr_id())
}

/// OCCT TopoDS_Shape::Reversed (TopAbs::Reverse).
pub(crate) fn shape_reversed(s: &Shape) -> Shape {
    let mut r = s.clone();
    r.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        o => o,
    };
    r
}

/// OCCT TopoDS_Shape::Oriented.
pub(crate) fn shape_oriented(s: &Shape, o: Orientation) -> Shape {
    let mut r = s.clone();
    r.orientation = o;
    r
}

// =============================================================================
// gp_Dir / gp_Ax1 / gp_Vec predicates
// =============================================================================

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance): the angle is 0 or PI
/// within the angular tolerance.
fn dir_is_parallel(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return false;
    }
    let cross = a.cross(b).length() / (la * lb);
    cross <= ang_tol.abs().sin() + 1e-18
}

/// OCCT gp_Ax1::IsParallel(Other, AngularTolerance): parallel directions
/// (the locations are ignored, gp_Ax1.cxx L63-73).
fn ax1_is_parallel(_loc_a: DVec3, dir_a: DVec3, _loc_b: DVec3, dir_b: DVec3, ang_tol: f64) -> bool {
    dir_is_parallel(dir_a, dir_b, ang_tol)
}

/// OCCT gp_Ax1::IsCoaxial(Other, AngularTolerance, LinearTolerance)
/// (gp_Ax1.cxx L84-99): parallel directions and the distance of the Other
/// location from This axis <= LinearTolerance.
fn ax1_is_coaxial(
    loc_a: DVec3,
    dir_a: DVec3,
    loc_b: DVec3,
    dir_b: DVec3,
    ang_tol: f64,
    lin_tol: f64,
) -> bool {
    if !dir_is_parallel(dir_a, dir_b, ang_tol) {
        return false;
    }
    let l = dir_a.length();
    if l < 1e-300 {
        return false;
    }
    let u = dir_a / l;
    let d = loc_b - loc_a;
    let proj = d.dot(u);
    (d - u * proj).length() <= lin_tol
}

/// OCCT gp_Vec::IsNormal(Other, AngularTolerance): the angle is PI/2 within
/// the angular tolerance.
fn vec_is_normal(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return false;
    }
    (a.dot(b) / (la * lb)).abs() <= ang_tol.abs().sin() + 1e-18
}

// =============================================================================
// GeomFill_Profiler (GeomFill_Profiler.cxx) — 1:1 translation
// =============================================================================

/// OCCT static UnifyByInsertingAllKnots (GeomFill_Profiler.cxx L29-78):
/// insert in the first curve the knot-vector of all the others, then insert
/// the grown knot vector into every other curve; rational curves get their
/// weights rescaled to a mean of 1.
fn unify_by_inserting_all_knots(the_curves: &mut [BSplineCurve3], ptol: f64) {
    // inserting in the first curve the knot-vector of all the others.
    for i in 1..the_curves.len() {
        let ki_mi = interior_knots_with_mults(&the_curves[i]);
        bspl_insert_knots(&mut the_curves[0], &ki_mi, ptol);
    }
    let new_knots = interior_knots_with_mults(&the_curves[0]);
    for i in 1..the_curves.len() {
        bspl_insert_knots(&mut the_curves[i], &new_knots, ptol);
    }

    // essai : tentative mise des poids sur chaque section a une moyenne 1
    for c in the_curves.iter_mut() {
        if c.weights.iter().any(|w| (*w - 1.0).abs() > 1e-15) {
            let np = c.weights.len();
            let mut sigma = 0.0;
            for w in &c.weights {
                sigma += *w;
            }
            sigma /= np as f64;
            for w in c.weights.iter_mut() {
                *w /= sigma;
            }
        }
    }
    // fin de l essai
}

/// OCCT static UnifyBySettingMiddleKnots (GeomFill_Profiler.cxx L81-113): set
/// every curve to the same knot vector whose interior knots are the means of
/// the curves' interior knots.
fn unify_by_setting_middle_knots(the_curves: &mut [BSplineCurve3]) {
    let nb_knots = knot_values(&the_curves[0]).len();
    let u_first = knot_values(&the_curves[0])[0];
    let u_last = knot_values(&the_curves[0])[nb_knots - 1];

    // Set middle values of knots
    let mut new_knots = vec![0.0f64; nb_knots];
    new_knots[0] = u_first;
    new_knots[nb_knots - 1] = u_last;
    for j in 1..nb_knots - 1 {
        let mut a_mid_knot = 0.0f64;
        for c in the_curves.iter() {
            a_mid_knot += knot_values(c)[j];
        }
        a_mid_knot /= the_curves.len() as f64;
        new_knots[j] = a_mid_knot;
    }

    for c in the_curves.iter_mut() {
        set_knots(c, &new_knots);
    }
}

/// Distinct knot values of a clamped BSpline curve (OCCT
/// Geom_BSplineCurve::Knots — one value per knot index).
fn knot_values(c: &BSplineCurve3) -> Vec<f64> {
    let d = c.degree;
    if c.knots.len() < 2 * (d + 1) {
        return c.knots.clone();
    }
    let mut out: Vec<f64> = Vec::new();
    for k in &c.knots {
        match out.last() {
            Some(last) if (*last - *k).abs() < 1e-15 => {}
            _ => out.push(*k),
        }
    }
    out
}

/// Distinct interior (knot, multiplicity) pairs of a clamped BSpline.
fn interior_knots_with_mults(c: &BSplineCurve3) -> Vec<(f64, usize)> {
    let d = c.degree;
    let mut out: Vec<(f64, usize)> = Vec::new();
    if c.knots.len() < 2 * (d + 1) {
        return out;
    }
    for (i, k) in c.knots.iter().enumerate() {
        if i < d + 1 || i >= c.knots.len() - d - 1 {
            continue;
        }
        match out.last_mut() {
            Some((last, m)) if (*last - *k).abs() < 1e-15 => *m += 1,
            _ => out.push((*k, 1)),
        }
    }
    out
}

/// OCCT Geom_BSplineCurve::InsertKnots(Knots, Mults, Tolerance, Add=false):
/// raise each listed interior knot to the requested multiplicity.
fn bspl_insert_knots(c: &mut BSplineCurve3, knots_mults: &[(f64, usize)], ptol: f64) {
    let _ = ptol;
    for (u, m) in knots_mults {
        let cur = c.knots.iter().filter(|k| (**k - *u).abs() < 1e-12).count();
        let target = (*m).min(c.degree);
        for _ in cur..target {
            if !insert_knot_one(c, *u) {
                break;
            }
        }
    }
}

/// Boehm knot insertion (one copy of the interior knot `u`); false when `u`
/// is not interior.
fn insert_knot_one(c: &mut BSplineCurve3, u: f64) -> bool {
    let d = c.degree;
    let n = c.control_points.len();
    let knots = c.knots.clone();
    if u <= knots[d] + 1e-12 || u >= knots[n + d] - 1e-12 {
        return false;
    }
    let mut k = d;
    for (i, tk) in knots.iter().enumerate() {
        if *tk <= u + 1e-12 {
            k = i;
        } else {
            break;
        }
    }
    k = k.clamp(d, n);
    let h: Vec<(DVec3, f64)> = (0..n)
        .map(|i| (c.control_points[i] * c.weights[i], c.weights[i]))
        .collect();
    let mut new_h: Vec<(DVec3, f64)> = Vec::with_capacity(n + 1);
    let lo = k + 2 - d;
    let hi = k + 1;
    for i1 in 1..=(n + 1) {
        let new_pole = if i1 < lo {
            h[i1 - 1]
        } else if i1 > hi {
            h[i1 - 2]
        } else {
            let ui = knots[i1 - 1];
            let uid = knots[i1 - 1 + d];
            let denom = uid - ui;
            let alpha = if denom.abs() < 1e-300 { 1.0 } else { (u - ui) / denom };
            (
                h[i1 - 1].0 * alpha + h[i1 - 2].0 * (1.0 - alpha),
                h[i1 - 1].1 * alpha + h[i1 - 2].1 * (1.0 - alpha),
            )
        };
        new_h.push(new_pole);
    }
    let mut new_knots: Vec<f64> = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(u);
    new_knots.extend_from_slice(&knots[k + 1..]);
    c.knots = new_knots;
    c.control_points = new_h.iter().map(|(p, _)| *p).collect();
    c.weights = new_h.iter().map(|(_, w)| *w).collect();
    true
}

/// OCCT Geom_BSplineCurve::SetKnots — replace the distinct knot values,
/// preserving multiplicities.
fn set_knots(c: &mut BSplineCurve3, new_knots: &[f64]) {
    let old = knot_values(c);
    if old.len() != new_knots.len() {
        return;
    }
    let n = new_knots.len();
    let flat = c.knots.clone();
    let mut out = flat.clone();
    for (i, k) in flat.iter().enumerate() {
        let mut di = 0usize;
        for (j, o) in old.iter().enumerate() {
            if (*o - *k).abs() < 1e-15 {
                di = j;
                break;
            }
        }
        out[i] = new_knots[di.min(n - 1)];
    }
    c.knots = out;
}

/// OCCT GeomFill_Profiler (GeomFill_Profiler.hxx L38-95) — unifies a sequence
/// of curves into BSpline curves of the same degree, range and knots.
#[derive(Debug, Clone)]
pub struct GeomFillProfiler {
    my_sequence: Vec<BSplineCurve3>,
    my_is_done: bool,
    my_is_periodic: bool,
}

impl Default for GeomFillProfiler {
    fn default() -> Self {
        Self::new()
    }
}

impl GeomFillProfiler {
    /// OCCT GeomFill_Profiler::GeomFill_Profiler (L116-121).
    pub fn new() -> Self {
        GeomFillProfiler {
            my_sequence: Vec::new(),
            my_is_done: false,
            my_is_periodic: true,
        }
    }

    /// OCCT GeomFill_Profiler::AddCurve (L124-161): strip the trim, convert
    /// conics by GeomConvert_ApproxCurve and the rest by GeomConvert to
    /// BSpline.
    pub fn add_curve(&mut self, curve: &Curve3) {
        // theCurve = TrimmedCurve ? BasisCurve : Curve (L134-137).
        let the_curve = match curve {
            Curve3::Trimmed(t) => (*t.curve).clone(),
            c => c.clone(),
        };
        // OCCT L138-143: GeomConvert_ApproxCurve for Geom_Conic kinds.
        // GAP: GeomConvert_ApproxCurve is not ported (reported); the section
        // edges reaching the generator are lines/circles/beziers/bsplines,
        // converted exactly by curve_to_bspline below.
        let c = curve_to_bspline(&the_curve);
        if self.my_is_periodic && !c.is_periodic {
            self.my_is_periodic = false;
        }
        self.my_sequence.push(c);
    }

    /// OCCT GeomFill_Profiler::Perform (L164-253).
    pub fn perform(&mut self, ptol: f64) {
        let mut my_degree = 0usize;
        let mut u_first = 0.0f64;
        let mut u_last = 0.0f64;
        let mut ecart_max = 0.0f64;

        for i in 0..self.my_sequence.len() {
            let c = &self.my_sequence[i];
            let d = c.degree;
            let u1 = c.knots[d];
            let u2 = c.knots[c.knots.len() - d - 1];

            // evaluate the max degree
            my_degree = my_degree.max(c.degree);

            // Calcul de Max ( Ufin - Udeb) sur l ensemble des courbes.
            if (u2 - u1) > ecart_max {
                ecart_max = u2 - u1;
                u_first = u1;
                u_last = u2;
            }
        }

        // increase the degree of the curves to my degree
        // reparametrize them in the range U1, U2.
        for i in 0..self.my_sequence.len() {
            let c = &mut self.my_sequence[i];
            if c.degree < my_degree {
                c.increase_degree(my_degree);
            }
            bspl_reparametrize(c, u_first, u_last);
        }

        let mut the_curves: Vec<BSplineCurve3> = self.my_sequence.clone();

        unify_by_inserting_all_knots(&mut the_curves, ptol);

        let mut unified = true;
        let the_nb_knots = knot_values(&the_curves[0]).len();
        for i in 1..the_curves.len() {
            if knot_values(&the_curves[i]).len() != the_nb_knots {
                unified = false;
                break;
            }
        }

        if unified {
            self.my_sequence = the_curves;
        } else {
            unify_by_setting_middle_knots(&mut self.my_sequence);
        }

        self.my_is_done = true;
    }

    /// OCCT GeomFill_Profiler::Degree (L256-265).
    pub fn degree(&self) -> usize {
        self.my_sequence[0].degree
    }

    /// OCCT GeomFill_Profiler::NbPoles (L268-275).
    pub fn nb_poles(&self) -> usize {
        self.my_sequence[0].control_points.len()
    }

    /// OCCT GeomFill_Profiler::IsPeriodic (L278-285).
    pub fn is_periodic(&self) -> bool {
        self.my_is_periodic
    }

    /// OCCT GeomFill_Profiler::KnotsAndMults — the flat knots of the first
    /// curve folded back to (knots, mults).
    pub fn knots_and_mults(&self) -> (Vec<f64>, Vec<usize>) {
        let d = self.my_sequence[0].degree;
        let mut knots: Vec<f64> = Vec::new();
        let mut mults: Vec<usize> = Vec::new();
        if self.my_sequence[0].knots.len() < 2 * (d + 1) {
            return (self.my_sequence[0].knots.clone(), vec![1; self.my_sequence[0].knots.len()]);
        }
        for k in &self.my_sequence[0].knots {
            match knots.last() {
                Some(last) if (*last - *k).abs() < 1e-15 => {
                    *mults.last_mut().unwrap() += 1;
                }
                _ => {
                    knots.push(*k);
                    mults.push(1);
                }
            }
        }
        (knots, mults)
    }
}

/// OCCT GeomConvert::CurveToBSplineCurve for the curve kinds carried by
/// generator section edges (Line, Circle, Bezier, BSpline).
fn curve_to_bspline(c: &Curve3) -> BSplineCurve3 {
    match c {
        Curve3::Line(l) => {
            // a line converts to the degree-1 BSpline over [0, 1]; the
            // profiler reparametrizes afterwards.
            BSplineCurve3 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![l.origin, l.origin + l.direction],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            }
        }
        Curve3::Circle(cir) => circle_arc_bspline(cir, 0.0, 2.0 * std::f64::consts::PI),
        Curve3::BSpline(b) => b.clone(),
        Curve3::Bezier(bz) => {
            let d = bz.control_points.len().saturating_sub(1);
            let mut knots = Vec::with_capacity(2 * (d + 1));
            for _ in 0..=d {
                knots.push(0.0);
            }
            for _ in 0..=d {
                knots.push(1.0);
            }
            BSplineCurve3 {
                degree: d,
                knots,
                control_points: bz.control_points.clone(),
                weights: if bz.weights.is_empty() {
                    vec![1.0; bz.control_points.len()]
                } else {
                    bz.weights.clone()
                },
                is_periodic: false,
            }
        }
        _ => BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![c.point_at(0.0), c.point_at(1.0)],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        },
    }
}

/// Exact rational degree-2 BSpline of a circular arc (at most PI/2 per
/// segment, weights [1, cos(seg/2), 1], interior knot multiplicity 2) — the
/// OCCT Convert_CircleToBSplineCurve form.
fn circle_arc_bspline(cir: &Circle3, t1: f64, t2: f64) -> BSplineCurve3 {
    let delta = t2 - t1;
    let nseg = ((delta / std::f64::consts::FRAC_PI_2).ceil().max(1.0)) as usize;
    let seg = delta / nseg as f64;
    let (o, r, x, y) = (cir.center, cir.radius, cir.x_dir, cir.y_dir);
    let pt = |a: f64| o + r * (a.cos() * x + a.sin() * y);
    let mut poles: Vec<DVec3> = Vec::with_capacity(2 * nseg + 1);
    let mut weights: Vec<f64> = Vec::with_capacity(2 * nseg + 1);
    let mut knots: Vec<f64> = Vec::with_capacity(2 * nseg + 4);
    for s in 0..nseg {
        let a = t1 + s as f64 * seg;
        let m = a + seg / 2.0;
        let w = (seg / 2.0).cos();
        if s == 0 {
            poles.push(pt(a));
            weights.push(1.0);
        }
        poles.push(o + r * (m.cos() * x + m.sin() * y) / w);
        weights.push(w);
        poles.push(pt(a + seg));
        weights.push(1.0);
    }
    knots.push(t1);
    knots.push(t1);
    knots.push(t1);
    for s in 1..nseg {
        let u = t1 + s as f64 * seg;
        knots.push(u);
        knots.push(u);
    }
    knots.push(t2);
    knots.push(t2);
    knots.push(t2);
    BSplineCurve3 {
        degree: 2,
        knots,
        control_points: poles,
        weights,
        is_periodic: false,
    }
}

/// OCCT BSplCLib::Reparametrize(U1, U2, Knots) — affine rescale of the flat
/// knot vector.
fn bspl_reparametrize(c: &mut BSplineCurve3, u1: f64, u2: f64) {
    let d = c.degree;
    if c.knots.len() <= 2 * d {
        return;
    }
    let k0 = c.knots[d];
    let k1 = c.knots[c.knots.len() - d - 1];
    let span = k1 - k0;
    if span.abs() < 1e-300 {
        return;
    }
    for k in c.knots.iter_mut() {
        *k = u1 + (*k - k0) * (u2 - u1) / span;
    }
}

// =============================================================================
// GeomFill_Generator (GeomFill_Generator.cxx L26-99) — 1:1 translation
// =============================================================================

/// OCCT GeomFill_Generator — builds the ruled (V-degree-1) surface through
/// the unified profile curves.
#[derive(Debug, Clone, Default)]
pub struct GeomFillGenerator {
    my_profiler: GeomFillProfiler,
    my_surface: Option<BSplineSurface>,
}

impl GeomFillGenerator {
    /// OCCT GeomFill_Generator::GeomFill_Generator (L29-31).
    pub fn new() -> Self {
        GeomFillGenerator {
            my_profiler: GeomFillProfiler::new(),
            my_surface: None,
        }
    }

    /// OCCT GeomFill_Generator::AddCurve (L38-41).
    pub fn add_curve(&mut self, curve: &Curve3) {
        self.my_profiler.add_curve(curve);
    }

    /// OCCT GeomFill_Generator::Perform (L44-96).
    pub fn perform(&mut self, ptol: f64) {
        // Perform the profile of the sections.
        self.my_profiler.perform(ptol);

        // Create the surface.
        let nb_u_poles = self.my_profiler.nb_poles();
        let nb_v_poles = self.my_profiler.my_sequence.len();
        if nb_u_poles == 0 || nb_v_poles < 2 {
            return;
        }
        let (u_knots, u_mults) = self.my_profiler.knots_and_mults();
        let degree = self.my_profiler.degree();
        let _is_u_periodic = self.my_profiler.is_periodic();
        let _is_v_periodic = false;

        // NCollection_Array2 Poles/Weights (1..NbUPoles, 1..NbVPoles).
        let mut poles: Vec<Vec<DVec3>> = vec![vec![DVec3::ZERO; nb_v_poles]; nb_u_poles];
        let mut weights: Vec<Vec<f64>> = vec![vec![0.0; nb_v_poles]; nb_u_poles];
        for (j, cj) in self.my_profiler.my_sequence.iter().enumerate() {
            if j >= nb_v_poles {
                break;
            }
            // VKnots(j) = (double)(j - 1) — the flat V knots are the clamped
            // [0, 0, 1, 1] below (VMults(1) = VMults(NbVKnots) = 2).
            for (i, p) in cj.control_points.iter().enumerate() {
                if i >= nb_u_poles {
                    break;
                }
                poles[i][j] = *p;
                weights[i][j] = cj.weights.get(i).copied().unwrap_or(1.0);
            }
        }
        let mut knots_u_flat: Vec<f64> = Vec::new();
        for (k, m) in u_knots.iter().zip(u_mults.iter()) {
            for _ in 0..*m {
                knots_u_flat.push(*k);
            }
        }

        self.my_surface = Some(BSplineSurface {
            degree_u: degree,
            degree_v: 1,
            knots_u: knots_u_flat,
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: poles,
            weights,
        });
    }

    /// OCCT GeomFill_Generator::Surface (L98-99).
    pub fn surface(&self) -> Option<BSplineSurface> {
        self.my_surface.clone()
    }
}

// =============================================================================
// Kernel-mapping helpers
// =============================================================================

/// OCCT TopExp::Vertices(E, V1, V2) — (first, last) vertices in the edge
/// traversal order (the orientation of the edge wrapper decides).
pub(crate) fn top_exp_vertices(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    let ed = brep.edge(e.clone());
    if e.orientation == Orientation::Reversed {
        (ed.last.clone(), ed.first.clone())
    } else {
        (ed.first.clone(), ed.last.clone())
    }
}

/// OCCT TopExp::Vertices(Wire, V1, V2) — first and last vertices of the wire
/// in traversal order.
pub(crate) fn top_exp_wire_vertices(brep: &BRep, w: &Shape) -> (Shape, Shape) {
    let wd = brep.wire(w.clone());
    let edges = &wd.edges;
    if edges.is_empty() {
        return (Shape::null(), Shape::null());
    }
    let (f, _) = top_exp_vertices(brep, &edges[0]);
    let (_, l) = top_exp_vertices(brep, &edges[edges.len() - 1]);
    (f, l)
}

/// OCCT BRep_Tool::IsClosed(Edge) — the edge Closed TShape flag.
pub(super) fn is_edge_closed(brep: &BRep, e: &Shape) -> bool {
    brep.edge(e.clone()).flags & tshape_flags::CLOSED != 0
}

/// OCCT BRep_Tool::SameParameter / SameRange flags.
pub(super) fn is_same_parameter(brep: &BRep, e: &Shape) -> bool {
    brep.edge(e.clone()).same_parameter
}

pub(super) fn is_same_range(brep: &BRep, e: &Shape) -> bool {
    brep.edge(e.clone()).same_range
}

/// Geom_Curve::ReversedParameter(p) — for the rcad reversed curves the
/// mapping is the negation.
fn reversed_parameter(_c: &Curve3, p: f64) -> f64 {
    -p
}

/// OCCT Geom_Curve::Reverse — the reversed parameterization.
pub(crate) fn curve_reversed(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Line(l) => Curve3::Line(Line3::new(l.origin + l.direction, -l.direction)),
        Curve3::Circle(cir) => {
            let mut out = *cir;
            out.y_dir = -out.y_dir;
            Curve3::Circle(out)
        }
        Curve3::BSpline(b) => {
            let mut b = b.clone();
            b.control_points.reverse();
            b.weights.reverse();
            let d = b.degree;
            let (k0, k1) = (b.knots[d], b.knots[b.knots.len() - d - 1]);
            b.knots = b.knots.iter().rev().map(|k| (k1 + k0) - k).collect();
            Curve3::BSpline(b)
        }
        Curve3::Bezier(bz) => {
            let mut b = bz.clone();
            b.control_points.reverse();
            b.weights.reverse();
            Curve3::Bezier(b)
        }
        Curve3::Trimmed(t) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(curve_reversed(&t.curve)),
            first: -t.last,
            last: -t.first,
        }),
        other => other.clone(),
    }
}

/// (FirstParameter, LastParameter) of the curve kinds produced here
/// (Geom_Curve::FirstParameter / LastParameter).
pub(super) fn curve_first_last(c: &Curve3) -> (f64, f64) {
    match c {
        Curve3::Line(_) => (0.0, 1.0),
        Curve3::Circle(_) => (0.0, 2.0 * std::f64::consts::PI),
        Curve3::BSpline(b) => {
            let d = b.degree;
            (b.knots[d], b.knots[b.knots.len() - d - 1])
        }
        Curve3::Bezier(_) => (0.0, 1.0),
        Curve3::Trimmed(t) => (t.first, t.last),
        _ => (0.0, 1.0),
    }
}

/// OCCT Geom_Surface::Bounds for the surface kinds produced here.
pub(super) fn surface_bounds(s: &Surface3) -> (f64, f64, f64, f64) {
    match s {
        Surface3::BSpline(b) => (
            b.knots_u[b.degree_u],
            b.knots_u[b.knots_u.len() - b.degree_u - 1],
            b.knots_v[b.degree_v],
            b.knots_v[b.knots_v.len() - b.degree_v - 1],
        ),
        Surface3::Trimmed(t) => (t.trim[0], t.trim[1], t.trim[2], t.trim[3]),
        _ => (-INFINITE_VALUE, INFINITE_VALUE, -INFINITE_VALUE, INFINITE_VALUE),
    }
}

/// OCCT Geom_Surface::UIso(u) for the ruled BSpline surface — the V-running
/// isoline at `u` (a degree-1 V-curve through a pole column).
pub(super) fn surface_uiiso(s: &Surface3, u: f64) -> Curve3 {
    if let Surface3::BSpline(b) = s {
        let d = b.degree_u;
        let last_col = b.control_points.len() - 1;
        let idx = if (u - b.knots_u[d]).abs() <= (u - b.knots_u[b.knots_u.len() - d - 1]).abs() {
            0
        } else {
            last_col
        };
        let col: Vec<DVec3> = b.control_points[idx].clone();
        let wcol: Vec<f64> = b.weights[idx].clone();
        return Curve3::BSpline(BSplineCurve3 {
            degree: b.degree_v,
            knots: b.knots_v.clone(),
            control_points: col,
            weights: wcol,
            is_periodic: false,
        });
    }
    // KPart surfaces take the Bezier branch at the call sites; this fallback
    // is never reached for IType == 0.
    Curve3::BSpline(BSplineCurve3 {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![s.point_at(u, 0.0), s.point_at(u, 1.0)],
        weights: vec![1.0, 1.0],
        is_periodic: false,
    })
}

/// Geom2d_BezierCurve through two poles (KPart pcurves).
pub(super) fn bezier2(p0: [f64; 2], p1: [f64; 2]) -> Curve2d {
    Curve2d::Bezier(BezierCurve2 {
        control_points: vec![DVec2::from(p0), DVec2::from(p1)],
        weights: vec![1.0, 1.0],
    })
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) + Range(E, F, f, l) — bind
/// the single pcurve on the face.
pub(super) fn bind_pcurve(brep: &mut BRep, e: &Shape, fkey: (u64, u32), pc: Curve2d, f: f64, l: f64) {
    brep.edge_mut_inplace(e.clone())
        .pcurves
        .insert(fkey, (pc, f, l));
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, F, Tol) — the seam (closed
/// surface) representation carrying the two pcurves.
pub(super) fn bind_seam_pcurves(brep: &mut BRep, e: &Shape, fkey: (u64, u32), pc1: Curve2d, pc2: Curve2d, f: f64, l: f64) {
    brep.edge_mut_inplace(e.clone())
        .representations
        .push(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
            face: fkey,
            pcurve1: pc1,
            pcurve2: pc2,
            range: [f, l],
        });
}

/// OCCT BRep_Builder::MakeEdge — always a NEW edge TShape; the kernel
/// add_tedge dedups curve-less edges, so the cache entry is evicted first.
pub(super) fn make_edge(brep: &mut BRep, curve: Option<Curve3>, first: Shape, last: Shape, range: [f64; 2]) -> Shape {
    if curve.is_none() {
        let ekey = (first.ptr_id(), last.ptr_id(), false);
        if brep.edge_by_key.contains_key(&ekey) {
            brep.edge_by_key.remove(&ekey);
        }
    }
    brep.add_tedge(curve, first, last, range)
}

// =============================================================================
// BRepFill_Generator.cxx — DetectKPart (L61-355)
// =============================================================================

/// The curve kind discriminator used by DetectKPart: (type, circle, line)
/// where type is 0 = Line, 1 = Circle, -1 = any other kind (GeomAdaptor_
/// Curve::GetType).
fn curve_type(c: &Curve3) -> (i32, Option<Circle3>, Option<Line3>) {
    match c {
        Curve3::Circle(cir) => (1, Some(*cir), None),
        Curve3::Line(l) => (0, None, Some(*l)),
        Curve3::Trimmed(t) => curve_type(&t.curve),
        _ => (-1, None, None),
    }
}

/// (point, derivative) at parameter t — GeomAdaptor_Curve::D1.
fn curve_d1(c: &Curve3, t: f64) -> (DVec3, DVec3) {
    (c.point_at(t), c.derivative_at(t))
}

/// OCCT static DetectKPart (BRepFill_Generator.cxx L61-355).  Returns the
/// particular-case type; -1 means a null 3D curve error, -2 the degenerated
/// cone case, 0 no particular case.
pub(super) fn detect_k_part(brep: &BRep, edge1: &Shape, edge2: &Shape) -> i32 {
    // initializations
    // !Note if IType set as -1 it means that occurs error with null 3d curve for the edge
    let mut itype: i32 = 0;

    // characteristics of the first edge
    let mut first1 = 0.0f64;
    let mut last1 = 0.0f64;
    // occ::handle<Geom_Curve> curv1 — function-scoped (null for a
    // degenerated edge), carried as Option in Rust.
    let mut curv1: Option<Curve3> = None;
    let degen1 = brep.edge(edge1.clone()).degenerated;

    // find the particular case
    let mut pos1 = DVec3::ZERO;
    let mut dist1 = 0.0f64;
    let mut axe1 = (DVec3::ZERO, DVec3::ZERO); // (location, direction)

    if degen1 {
        itype = -2;
        let (v1, _v2) = top_exp_vertices(brep, edge1);
        pos1 = brep.vertex(v1).point;
    } else {
        let ed = brep.edge(edge1.clone());
        let c = match &ed.curve {
            Some(c) => c.clone(),
            None => return -1,
        };
        first1 = ed.range[0];
        last1 = ed.range[1];
        let ff = first1;
        let ll = last1;
        let mut c = c;
        if edge1.orientation == Orientation::Reversed {
            c = curve_reversed(&c);
            first1 = reversed_parameter(&c, ll);
            last1 = reversed_parameter(&c, ff);
        }
        let (ctype1, circle1, line1) = curve_type(&c);
        if ctype1 == 1 {
            // first circular section
            itype = 1;
            let cir = circle1.unwrap();
            pos1 = cir.center;
            dist1 = cir.radius;
            axe1 = (cir.center, cir.normal);
        } else if ctype1 == 0 {
            // first straight line section
            itype = 4;
            let lin = line1.unwrap();
            pos1 = lin.origin;
            dist1 = c.point_at(first1).distance(c.point_at(last1));
            let vec = c.point_at(last1) - c.point_at(first1);
            let dir = vec.normalize_or_zero();
            axe1 = (c.point_at(first1), dir);
        } else {
            // first section of any type
            itype = 0;
        }
        curv1 = Some(c);
    }

    if itype != 0 {
        let degen2 = brep.edge(edge2.clone()).degenerated;
        if degen2 {
            let (v1, _v2) = top_exp_vertices(brep, edge2);
            let pos = brep.vertex(v1).point;
            if itype == 1 {
                // the only particular case with degenerated edge at end : the cone
                if pos1.distance(pos) < TOL_CONFUSION {
                    // the top is mixed with the center of the circle
                    itype = 0;
                } else {
                    let vec = pos - pos1;
                    let dir = vec.normalize_or_zero();
                    if ax1_is_parallel(pos1, dir, axe1.0, axe1.1, TOL_ANGULAR) {
                        // the top is on the axis of the circle
                        itype = 2;
                    } else {
                        // incorrect top --> no particular case
                        itype = 0;
                    }
                }
            } else if itype != 4 {
                // not a plane
                // no particular case
                itype = 0;
            }
        } else {
            let ed2 = brep.edge(edge2.clone());
            let mut curv2 = match &ed2.curve {
                Some(c) => c.clone(),
                None => return -1,
            };
            let mut first2 = ed2.range[0];
            let mut last2 = ed2.range[1];
            let ff = first2;
            let ll = last2;
            if edge2.orientation == Orientation::Reversed {
                curv2 = curve_reversed(&curv2);
                first2 = reversed_parameter(&curv2, ll);
                last2 = reversed_parameter(&curv2, ff);
            }
            let (ctype2, circle2, line2) = curve_type(&curv2);

            if itype > 0 && itype < 4 {
                if ctype2 != 1 {
                    // section not circular --> no particular case
                    itype = 0;
                } else {
                    let cir2 = circle2.unwrap();
                    if ax1_is_coaxial(cir2.center, cir2.normal, axe1.0, axe1.1, TOL_ANGULAR, TOL_CONFUSION) {
                        // same axis
                        if (cir2.radius - dist1).abs() < TOL_CONFUSION {
                            // possibility of cylinder or a piece of cylinder
                            let h1 = (last1 - first1).abs();
                            let h2 = (last2 - first2).abs();
                            let same_parametric_length = (h1 - h2).abs() < TOL_PCONFUSION;
                            let m1 = (first1 + last1) / 2.0;
                            let m2 = (first2 + last2) / 2.0;
                            let (p1, du) = curve_d1(curv1.as_ref().unwrap(), m1);
                            let p2 = curv2.point_at(m2);
                            let same =
                                same_parametric_length && vec_is_normal(p2 - p1, du, TOL_ANGULAR);
                            if same {
                                // cylinder or piece of cylinder
                                itype = 1;
                            } else {
                                // the interval of definition is not correct
                                itype = 0;
                            }
                        } else {
                            // possibility of cone truncation
                            let h1 = (last1 - first1).abs();
                            let h2 = (last2 - first2).abs();
                            let same_parametric_length = (h1 - h2).abs() < TOL_PCONFUSION;
                            let m1 = (first1 + last1) / 2.0;
                            let m2 = (first2 + last2) / 2.0;
                            let (p1, du) = curve_d1(curv1.as_ref().unwrap(), m1);
                            let p2 = curv2.point_at(m2);
                            let same =
                                same_parametric_length && vec_is_normal(p2 - p1, du, TOL_ANGULAR);
                            if same {
                                // truncation of cone
                                itype = 2;
                            } else {
                                // the interval of definition is not correct
                                itype = 0;
                            }
                        }
                        if cir2.center.distance(pos1) < TOL_CONFUSION {
                            // the centers are mixed
                            itype = 0;
                        }
                    } else {
                        // different axis
                        if cir2.radius == dist1 {
                            // torus ?
                            itype = 3;
                        } else {
                            // different radius --> no particular case
                            itype = 0;
                        }
                    }
                }
            } else if itype >= 4 {
                if ctype2 != 0 {
                    // not a straight line section --> no particular case
                    itype = 0;
                } else {
                    let lin2 = line2.unwrap();
                    let pos = lin2.origin;
                    let dist = curv2.point_at(first2).distance(curv2.point_at(last2));
                    let vec = curv2.point_at(last2) - curv2.point_at(first2);
                    let a_dir = vec.normalize_or_zero();
                    if ax1_is_parallel(
                        curv2.point_at(first2),
                        a_dir,
                        axe1.0,
                        axe1.1,
                        TOL_ANGULAR,
                    ) {
                        // parallel straight line
                        if (dist - dist1).abs() < TOL_CONFUSION {
                            let dir = (curv2.point_at(first2) - curv1.as_ref().unwrap().point_at(first1))
                                .normalize_or_zero();
                            if vec_is_normal(dir, a_dir, TOL_ANGULAR) {
                                // plane
                                itype = 4;
                            } else {
                                // extrusion ?
                                itype = 5;
                            }
                        } else {
                            // different length --> no particular case
                            itype = 0;
                        }
                    } else {
                        // not parallel straight line --> no particular case
                        itype = 0;
                    }
                }
            } else if itype == -2 {
                if ctype2 == 0 {
                    itype = 4; // plane
                } else if ctype2 == 1 {
                    // the only particular case with degenerated edge at the beginning the cone
                    let cir2 = circle2.unwrap();
                    let pos = cir2.center;
                    let axe = (cir2.center, cir2.normal);
                    if pos1.distance(pos) < TOL_CONFUSION {
                        // the top is mixed with the center of the circle
                        itype = 0;
                    } else {
                        let vec = pos - pos1;
                        let dir = vec.normalize_or_zero();
                        let axe1_new = (pos1, dir);
                        if ax1_is_parallel(axe.0, axe.1, axe1_new.0, axe1_new.1, TOL_ANGULAR) {
                            // the top is on the axis of the circle
                            itype = -2;
                        } else {
                            // incorrect top --> no particular case
                            itype = 0;
                        }
                    }
                } else {
                    itype = 0;
                }
            }
        }
    }
    // torus and extrusion are not particular cases.
    if itype == 3 || itype == 5 {
        itype = 0;
    }
    itype
}

// =============================================================================
// BRepFill_Generator.cxx — CreateKPart (L362-569)
// =============================================================================

/// OCCT static CreateKPart (BRepFill_Generator.cxx L362-569).  Returns None
/// when an error occurs (the C++ returns false).
pub(super) fn create_k_part(brep: &BRep, edge1: &Shape, edge2: &Shape, itype: i32) -> Option<Surface3> {
    // find the dimension
    let mut aa = 0.0f64;
    let mut bb = 0.0f64;

    let degen1 = brep.edge(edge1.clone()).degenerated;
    let degen2 = brep.edge(edge2.clone()).degenerated;

    // find characteristics of the first edge
    let mut c1: Option<Curve3> = None;
    let (mut v1f, mut v1l) = (Shape::null(), Shape::null());
    if degen1 {
        // cone with degenerated edge at the top
        let ed = brep.edge(edge1.clone());
        v1f = ed.first.clone();
        v1l = ed.last.clone();
    } else {
        let ed = brep.edge(edge1.clone());
        let a1 = ed.range[0];
        let b1 = ed.range[1];
        let mut c = ed.curve.clone()?;
        aa = a1;
        bb = b1;
        if edge1.orientation == Orientation::Reversed {
            c = curve_reversed(&c);
            aa = reversed_parameter(&c, b1);
            bb = reversed_parameter(&c, a1);
            // TopExp::Vertices(Edge1, v1l, v1f)
            let (l, f) = top_exp_vertices(brep, edge1);
            v1f = f;
            v1l = l;
        } else {
            // TopExp::Vertices(Edge1, v1f, v1l)
            let (f, l) = top_exp_vertices(brep, edge1);
            v1f = f;
            v1l = l;
        }
        c1 = Some(c);
    }

    // find characteristics of the second edge
    let mut c2: Option<Curve3> = None;
    let (mut v2f, mut v2l) = (Shape::null(), Shape::null());
    if degen2 {
        // cone with degenerated edge at the top
        let ed = brep.edge(edge2.clone());
        v2f = ed.first.clone();
        v2l = ed.last.clone();
    } else {
        let ed = brep.edge(edge2.clone());
        let a1 = ed.range[0];
        let b1 = ed.range[1];
        let mut c = ed.curve.clone()?;
        if edge2.orientation == Orientation::Reversed {
            c = curve_reversed(&c);
            if degen1 {
                aa = a1;
                bb = b1;
            }
            // TopExp::Vertices(Edge2, v2l, v2f)
            let (l, f) = top_exp_vertices(brep, edge2);
            v2f = f;
            v2l = l;
        } else {
            if degen1 {
                aa = a1;
                bb = b1;
            }
            // TopExp::Vertices(Edge2, v2f, v2l)
            let (f, l) = top_exp_vertices(brep, edge2);
            v2f = f;
            v2l = l;
        }
        c2 = Some(c);
    }

    // create the new surface
    let surface: Option<Surface3> = if itype == 1 {
        // cylindrical surface
        let c1c = match c1.as_ref() {
            Some(Curve3::Circle(c)) => *c,
            _ => return None,
        };
        let c2c = match c2.as_ref() {
            Some(Curve3::Circle(c)) => *c,
            _ => return None,
        };
        let mut axis = c1c.normal;
        let mut v = (c2c.center - c1c.center).dot(axis);
        if v < 0.0 {
            axis = -axis;
            v = -v;
        }
        let cyl = CylindricalSurface {
            origin: c1c.center,
            axis,
            radius: c1c.radius,
            ref_dir: c1c.x_dir,
            y_dir: Some(c1c.y_dir),
        };
        Some(Surface3::Trimmed(TrimmedSurface::new(
            Surface3::Cylinder(cyl),
            aa,
            bb,
            0.0f64.min(v),
            0.0f64.max(v),
        )))
    } else if itype == 2 {
        // conical surface
        let k1 = match c1.as_ref() {
            Some(Curve3::Circle(c)) => *c,
            _ => return None,
        };
        let mut axis = k1.normal;
        let (mut v, rad) = if degen2 {
            let v = (brep.vertex(v2f.clone()).point - k1.center).dot(axis);
            (v, -k1.radius)
        } else {
            let k2 = match c2.as_ref() {
                Some(Curve3::Circle(c)) => *c,
                _ => return None,
            };
            let v = (k2.center - k1.center).dot(axis);
            (v, k2.radius - k1.radius)
        };
        if v < 0.0 {
            axis = -axis;
            v = -v;
        }
        let ang = (rad / v).atan();
        let cone = ConicalSurface {
            apex: k1.center,
            axis,
            radius: k1.radius,
            half_angle_rad: ang,
            ref_dir: k1.x_dir,
        };
        let v = v / ang.cos();
        Some(Surface3::Trimmed(TrimmedSurface::new(
            Surface3::Cone(cone),
            aa,
            bb,
            0.0f64.min(v),
            0.0f64.max(v),
        )))
    } else if itype == -2 {
        // conical surface with the top at the beginning (degen1 is true)
        let k2 = match c2.as_ref() {
            Some(Curve3::Circle(c)) => *c,
            _ => return None,
        };
        let mut axis = k2.normal;
        let apex = brep.vertex(v1f.clone()).point;
        let mut v = (k2.center - apex).dot(axis);
        let rad = k2.radius; // - k2.Radius();
        if v < 0.0 {
            axis = -axis;
            v = -v;
        }
        let ang = (rad / v).atan();
        let cone = ConicalSurface {
            apex,
            axis,
            radius: 0.0,
            half_angle_rad: ang,
            ref_dir: k2.x_dir,
        };
        let v = v / ang.cos();
        Some(Surface3::Trimmed(TrimmedSurface::new(
            Surface3::Cone(cone),
            aa,
            bb,
            0.0f64.min(v),
            0.0f64.max(v),
        )))
    } else if itype == 3 {
        // torus surface ?
        None
    } else if itype == 4 {
        // surface plane
        let l1 = match c1.as_ref() {
            Some(Curve3::Line(l)) => Some(*l),
            _ => None,
        };
        let l2 = match c2.as_ref() {
            Some(Curve3::Line(l)) => Some(*l),
            _ => None,
        };
        // aLine ends as L2 when both exist (the later assignment wins).
        let a_line = if l2.is_some() { l2 } else { l1 };

        let p1 = if degen1 {
            brep.vertex(v1f.clone()).point
        } else {
            l1?.origin
        };
        let p2 = if degen2 {
            brep.vertex(v2f.clone()).point
        } else {
            l2?.origin
        };

        let p1p2 = p2 - p1;
        let d1 = a_line?.direction;
        let normal = d1.cross(p1p2).normalize_or_zero();
        let a_line = a_line?;
        Some(Surface3::Plane(rcad_kernel::geom::Plane {
            origin: a_line.origin,
            normal,
            u_dir: d1,
            v_dir: normal.cross(d1).normalize_or_zero(),
        }))
    } else if itype == 5 {
        // surface of extrusion ?
        None
    } else {
        // IType incorrect
        None
    };
    surface
}

// =============================================================================
// BRepFill_Generator.cxx — CreateNewEdge (L573-594)
// =============================================================================

/// OCCT static CreateNewEdge (BRepFill_Generator.cxx L573-594).
pub(super) fn create_new_edge(
    brep: &mut BRep,
    the_edge: &Shape,
    the_copied_edges: &mut HashMap<ShapeKey, Shape>,
    the_wire: &Shape,
    the_modif_wires: &mut Vec<Shape>,
) -> Shape {
    let a_new_edge = brep.empty_copied(the_edge);
    // TopoDS_Iterator over theEdge children (its two vertices).
    let ed = brep.edge(the_edge.clone());
    for v in [ed.first.clone(), ed.last.clone()] {
        // aBuilder.Add(aNewEdge, anIterator.Value())
        let _ = v;
    }
    the_copied_edges.insert(shape_key(the_edge), a_new_edge.clone());

    if !the_modif_wires.iter().any(|w| w.is_same(the_wire)) {
        the_modif_wires.push(the_wire.clone());
    }
    a_new_edge
}

// =============================================================================
// BRepFill_Generator class (BRepFill_Generator.hxx L38-86)
// =============================================================================

/// OCCT BRepFill_Generator — compute a shell of ruled faces through the
/// section wires (BRepFill_Generator.cxx).
#[derive(Debug)]
pub struct BRepFillGenerator {
    pub(super) my_wires: Vec<Shape>,
    pub(super) my_shell: Option<Shape>,
    pub(super) my_map: HashMap<ShapeKey, Vec<Shape>>,
    pub(super) my_old_new_shapes: HashMap<ShapeKey, Shape>,
    pub(super) my_reshaper: ShapeBuildReShape,
    pub(super) my_mutable_input: bool,
    pub(super) my_status: BRepFillThruSectionErrorStatus,
}

impl Default for BRepFillGenerator {
    fn default() -> Self {
        Self::new()
    }
}
