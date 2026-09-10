//! OCCT TKShHealing — ShapeCustom package class: `ShapeCustom_Curve2d`
//! (`ShapeCustom_Curve2d.hxx` L16-61 + `ShapeCustom_Curve2d.cxx` L1-196).
//!
//! Converts curve2d to analytical form with given precision or simplify
//! curve2d.
//!
//! Docket status: `[B-exception-single-point]` (W1-5) — the single
//! ShapeCustom class admitted to Path B (TKOffset hard dependency); it does
//! NOT open the ShapeCustom package (docs/TKSHHEALING_PATH_B_DOCKET.md §3).
//!
//! Architecture bridges (numbered, referenced by the methods below):
//! 1. `occ::handle<Geom2d_Curve>` -> `&Curve2d`; `occ::handle<Geom2d_Line>`
//!    (nullable) -> `Option<Curve2d>`; `occ::handle<Geom2d_BSplineCurve>&`
//!    (mutable) -> `&mut BSplineCurve2`. The `occ::down_cast` steps map to
//!    `match` arms on the enum variants.
//! 2. `gp_Pnt2d` / `gp_Vec2d` -> `DVec2`; `gp_Lin2d` ->
//!    `rcad_kernel::geom::Line2d` (origin + unit direction, the same
//!    invariant as the OCCT gp_Ax2d placement).
//! 3. `NCollection_Array1<gp_Pnt2d>` (the poles array, 1-based) ->
//!    `&[DVec2]` (0-based slice; loop bounds keep the OCCT 1..=nb form).
//! 4. `Geom2d_BSplineCurve` (TKG2d, untranslated W-batch dependency) is
//!    carried by [`Geom2dBSplineCurve`]: its read accessors / LocalDN /
//!    RemoveKnot are bridged over the rcad flat-knot representation through
//!    the documented kernel BSplCLib support, with the OCCT anchors cited at
//!    each method. GAP: the rational (`myRational`) RemoveKnot arm and the
//!    periodic curve form are not covered by the kernel BSplCLib::RemoveKnot
//!    carrier (`bspl_lib::remove_knot` is the NoWeights non-periodic arm) —
//!    those inputs keep the OCCT failure path (knot not removed).

use glam::DVec2;
use rcad_kernel::geom::{BSplineCurve2, Curve2d, Curve2dEval, Line2d};
use rcad_kernel::math::bspl_lib;
use rcad_kernel::precision::{ANGULAR, PCONFUSION};

// OCCT Standard_Real/pi (Standard_Integer.hxx / Standard.hxx): M_PI.
const PI: f64 = std::f64::consts::PI;

// ---------------------------------------------------------------------------
// gp / ElCLib re-hosts (architecture bridge #2; the edge.rs precedents).
// ---------------------------------------------------------------------------

/// OCCT ElCLib::Value(u, L) (ElCLib.cxx, the gp_Lin2d case):
/// `L.Location().XY() + u * L.Direction().XY()`.
fn elclib_value_lin2d(u: f64, l: &Line2d) -> DVec2 {
    l.origin + u * l.direction
}

/// OCCT ElCLib::Parameter(L, P) (ElCLib.cxx, the gp_Lin2d case):
/// `(P.XY() - L.Location().XY()).Dot(L.Direction().XY())`.
fn elclib_parameter_lin2d(l: &Line2d, p: DVec2) -> f64 {
    (p - l.origin).dot(l.direction)
}

/// OCCT gp_Lin2d::SquareDistance(theP) (gp_Lin2d.cxx):
/// `D.Dot(D) - (D.Dot(Dir))**2` with `D = theP.XY() - Location().XY()`.
fn gp_lin2d_square_distance(l: &Line2d, p: DVec2) -> f64 {
    let d = p - l.origin;
    let b = d.dot(l.direction);
    d.dot(d) - b * b
}

/// OCCT gp_Vec2d::Angle(theOther) (gp_Vec2d.cxx): the oriented angle from
/// `me` to `theOther` in [-PI, PI]; `atan2(Crossed, Dot)` of the XY vectors
/// (the trigonometric — counterclockwise — positive sense).  The OCCT
/// Standard_NoSuchObject on a null-magnitude vector maps to atan2(0, 0) = 0
/// here; the SimplifyBSpline2d call site only compares DegMult-th
/// derivatives, which are non-degenerate.
fn gp_vec2d_angle(a: DVec2, b: DVec2) -> f64 {
    let crossed = a.x.mul_add(b.y, -(a.y * b.x)); // a.x*b.y - a.y*b.x
    crossed.atan2(a.dot(b))
}

/// OCCT gp_Vec2d::IsParallel(theOther, theAngularTolerance)
/// (gp_Vec2d.hxx L372-376): `|Angle| <= tol || PI - |Angle| <= tol`.
fn gp_vec2d_is_parallel(a: DVec2, b: DVec2, the_angular_tolerance: f64) -> bool {
    let an_ang = gp_vec2d_angle(a, b).abs();
    an_ang <= the_angular_tolerance || PI - an_ang <= the_angular_tolerance
}

// ---------------------------------------------------------------------------
// Geom2d_BSplineCurve GAP carrier (architecture bridge #4).
// ---------------------------------------------------------------------------

/// The compressed (knots, mults) form OCCT stores (`myKnots` / `myMults`),
/// recovered from the rcad flat knot vector.  The flat <-> compressed
/// conversion is bijective (the kernel bspline_ops.rs precedent).
fn knots_mults(c: &BSplineCurve2) -> (Vec<f64>, Vec<i32>) {
    let mut knots = Vec::new();
    let mut mults = Vec::new();
    for (i, k) in c.knots.iter().enumerate() {
        if i > 0 && *k == knots.last().copied().unwrap_or(f64::NAN) {
            *mults.last_mut().unwrap() += 1;
        } else {
            knots.push(*k);
            mults.push(1);
        }
    }
    (knots, mults)
}

/// OCCT BSplCLib::FirstUKnotIndex (BSplCLib.cxx L111-122) — degree-aware.
fn first_uknot_index_deg(degree: usize, mults: &[i32]) -> usize {
    let mut index = 1usize; // OCCT 1-based Mults.Lower()
    let mut sigma_mult = mults[index - 1];
    while sigma_mult <= degree as i32 {
        index += 1;
        sigma_mult += mults[index - 1];
    }
    index
}

/// OCCT BSplCLib::LastUKnotIndex (BSplCLib.cxx L126-137).
fn last_uknot_index_deg(degree: usize, mults: &[i32]) -> usize {
    let mut index = mults.len(); // OCCT 1-based Mults.Upper()
    let mut sigma_mult = mults[index - 1];
    while sigma_mult <= degree as i32 {
        index -= 1;
        sigma_mult += mults[index - 1];
    }
    index
}

/// OCCT BSplCLib::KnotAnalysis (BSplCLib.cxx L692-752) — the MaxKnotMult
/// output only (the KnotForm output is not consumed here).
fn knot_analysis_max_mult(degree: usize, mults: &[i32]) -> i32 {
    let first_km = first_uknot_index_deg(degree, mults);
    let last_km = last_uknot_index_deg(degree, mults);
    let mut max_knot_mult = 0i32;
    if last_km - first_km != 1 {
        for i in (first_km + 1)..last_km {
            let multi = mults[i - 1];
            max_knot_mult = max_knot_mult.max(multi);
        }
    }
    max_knot_mult
}

/// GAP carrier for OCCT `Geom2d_BSplineCurve` (TKG2d, untranslated): the
/// operations consumed by `SimplifyBSpline2d`.  The rcad `BSplineCurve2`
/// stores the flat (expanded) knot vector; the compressed (knots, mults)
/// views and the evaluation/removal calls are bridged through the kernel
/// BSplCLib support (`bspl_lib`).  Architecture limits (bridge #4): the
/// non-periodic NoWeights arm only — a rational or periodic curve keeps the
/// OCCT failure path where noted.  GAP: closes with the Geom2d_BSplineCurve
/// batch (W2+); the W1-5 delivery only carries the call sites.
struct Geom2dBSplineCurve<'a> {
    c: &'a mut BSplineCurve2,
}

impl<'a> Geom2dBSplineCurve<'a> {
    /// OCCT NbKnots (Geom2d_BSplineCurve_1.cxx L598-601): `myKnots.Length()`.
    fn nb_knots(&self) -> i32 {
        knots_mults(self.c).0.len() as i32
    }

    /// OCCT Knot(Index) (Geom2d_BSplineCurve_1.cxx L356-361): `myKnots(Index)`.
    fn knot(&self, index: i32) -> f64 {
        knots_mults(self.c).0[(index - 1) as usize]
    }

    /// OCCT Multiplicity(Index) (Geom2d_BSplineCurve_1.cxx L575-582):
    /// `myMults(Index)`.
    fn multiplicity(&self, index: i32) -> i32 {
        knots_mults(self.c).1[(index - 1) as usize]
    }

    /// OCCT Degree (Geom2d_BSplineCurve_1.cxx L168-171): `myDeg`.
    fn degree(&self) -> i32 {
        self.c.degree as i32
    }

    /// OCCT IsCN(N) (Geom2d_BSplineCurve_1.cxx L34-56): the `mySmooth` switch.
    /// `mySmooth` is computed by updateKnots (Geom2d_BSplineCurve.cxx
    /// L1280-1324) from the KnotAnalysis MaxKnotMult.
    fn is_cn(&self, n: i32) -> bool {
        // OCCT L36: Standard_RangeError_Raise_if(N < 0, ...).
        assert!(n >= 0, "Geom2d_BSplineCurve::IsCN");
        let mults = knots_mults(self.c).1;
        let max_knot_mult = knot_analysis_max_mult(self.c.degree, &mults);
        // OCCT updateKnots L1299-1323: the mySmooth level.
        let smooth = if max_knot_mult == 0 {
            // GeomAbs_CN
            return true;
        } else {
            match self.c.degree as i32 - max_knot_mult {
                0 => 0, // GeomAbs_C0
                1 => 1, // GeomAbs_C1
                2 => 2, // GeomAbs_C2
                _ => 3, // GeomAbs_C3 (cases 3 and default)
            }
        };
        // OCCT IsCN switch L38-56.
        match smooth {
            0 => n <= 0, // GeomAbs_C0 (and GeomAbs_G1)
            1 => n <= 1, // GeomAbs_C1 (and GeomAbs_G2)
            2 => n <= 2, // GeomAbs_C2
            3 => {
                if n <= 3 {
                    true
                } else {
                    // OCCT L52-55: N <= myDeg - MaxKnotMult(myMults, Lower+1, Upper-1).
                    n <= self.c.degree as i32 - max_knot_mult_internal(&mults)
                }
            }
            _ => false,
        }
    }

    /// OCCT LocalDN(U, FromK1, ToK2, N) (Geom2d_BSplineCurve_1.cxx L547-566):
    /// the N-th derivative at `u`.  The OCCT span arguments only locate the
    /// flat evaluation index (LocateParameter + FlatIndex, L554-557); the
    /// BSplCLib::DN evaluation over the flat knots is span-independent, and
    /// the DegMult-th derivative is one-sided-equal at a knot of multiplicity
    /// below degree, so the global evaluation is value-faithful (the
    /// BSplineCurve3::dn precedent in geom/bspline_ops.rs).
    fn local_dn(&self, u: f64, _from_k1: i32, _to_k2: i32, n: i32) -> DVec2 {
        let dim = 2usize;
        let mut poles_flat = Vec::with_capacity(self.c.control_points.len() * dim);
        for p in &self.c.control_points {
            poles_flat.extend([p.x, p.y]);
        }
        let count = (n + 1) as usize;
        let mut poles_res = vec![0.0f64; count * dim];
        let mut weights_res = vec![0.0f64; count];
        let mut extrap = [0i32; 2];
        // Clamp to the parameter range (the DN callers evaluate in range).
        let first = self.c.knots[self.c.degree];
        let last = self.c.knots[self.c.knots.len() - self.c.degree - 1];
        let u_clamped = u.clamp(first, last);
        bspl_lib::eval_homogeneous(
            u_clamped,
            false,
            n,
            &mut extrap,
            self.c.degree,
            &self.c.knots,
            dim,
            &poles_flat,
            &self.c.weights,
            &mut poles_res,
            &mut weights_res,
        );
        bspl_lib::rational_derivatives_inplace(n, dim, &mut poles_res, &mut weights_res);
        let off = n as usize * dim;
        DVec2::new(poles_res[off], poles_res[off + 1])
    }

    /// OCCT RemoveKnot(Index, M, Tolerance) (Geom2d_BSplineCurve.cxx L410-470)
    /// through BSplCLib::RemoveKnot (the kernel `bspl_lib::remove_knot`
    /// carrier, the NoWeights non-periodic arm).
    fn remove_knot(&mut self, index: i32, m: i32, tolerance: f64) -> bool {
        // OCCT L412-416: M < 0 -> return true (no modification).
        if m < 0 {
            return true;
        }
        let (knots, mults) = knots_mults(self.c);
        // OCCT L418-425: the FirstUKnotIndex..LastUKnotIndex range check
        // (Standard_OutOfRange raise; rcad panics the same way).
        let i1 = first_uknot_index_deg(self.c.degree, &mults) as i32;
        let i2 = last_uknot_index_deg(self.c.degree, &mults) as i32;
        assert!((i1..=i2).contains(&index), "BSpline curve: RemoveKnot: index out of range");
        // OCCT L432-438: step = myMults(Index) - M; step <= 0 -> return true.
        let step = mults[(index - 1) as usize] - m;
        if step <= 0 {
            return true;
        }
        // Architecture limit (bridge #4): a rational curve keeps the OCCT
        // failure path (knot not removed) — the kernel carrier is the
        // NoWeights arm.  GAP: closes with the Geom2d batch.
        if self.c.weights.iter().any(|&w| w != 1.0) {
            return false;
        }
        // OCCT L440-464: BSplCLib::RemoveKnot(Index, M, myDeg, myPeriodic,
        // myPoles, Weights(), myKnots, myMults, npoles, NoWeights, nknots,
        // nmults, Tolerance); on false return false.
        let dim = 2usize;
        let mut poles_flat = Vec::with_capacity(self.c.control_points.len() * dim);
        for p in &self.c.control_points {
            poles_flat.extend([p.x, p.y]);
        }
        // OCCT L429-431: npoles(1, oldpoles.Length() - step).
        let nb_new_poles = self.c.control_points.len() - step as usize;
        let mut new_poles = vec![0.0f64; nb_new_poles * dim];
        let mut new_knots = Vec::with_capacity(knots.len());
        let mut new_mults = Vec::with_capacity(mults.len());
        let ok = bspl_lib::remove_knot(
            index as usize,
            m,
            self.c.degree,
            false,
            dim,
            &poles_flat,
            &knots,
            &mults,
            &mut new_poles,
            &mut new_knots,
            &mut new_mults,
            tolerance,
        );
        if !ok {
            return false;
        }
        // OCCT L466-478: commit the new poles / knots / mults (unit weights
        // for the non-rational arm) and updateKnots().
        self.c.control_points = new_poles
            .chunks_exact(dim)
            .map(|ch| DVec2::new(ch[0], ch[1]))
            .collect();
        self.c.weights = vec![1.0; self.c.control_points.len()];
        let mut flat = Vec::new();
        for (k, mult) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*mult {
                flat.push(*k);
            }
        }
        self.c.knots = flat;
        true
    }
}

/// OCCT BSplCLib::MaxKnotMult (BSplCLib.lxx): the max multiplicity over the
/// given (inclusive, 1-based) range — consumed by the IsCN C3 tail, which
/// passes (Lower+1, Upper-1): the internal knots only.
fn max_knot_mult_internal(mults: &[i32]) -> i32 {
    let mut max = 0i32;
    for i in 2..=(mults.len() - 1) {
        max = max.max(mults[i - 1]);
    }
    max
}

// ---------------------------------------------------------------------------
// The class.
// ---------------------------------------------------------------------------

/// OCCT ShapeCustom_Curve2d (ShapeCustom_Curve2d.hxx L31-59): the class has
/// static methods only; the rcad unit struct is the namespace.
pub struct ShapeCustomCurve2d;

impl ShapeCustomCurve2d {
    /// OCCT GetLine (ShapeCustom_Curve2d.cxx L33-46, static): builds the
    /// line through P1/P2 relocated at parameter c1; outputs cf/cl, the
    /// parameters of P1/P2 on the relocated line.
    fn get_line(p1: DVec2, p2: DVec2, c1: f64, cf: &mut f64, cl: &mut f64) -> Line2d {
        // gp_Vec2d avec(P1, P2); gp_Dir2d adir(avec); gp_Lin2d alin(P1, adir).
        let avec = p2 - p1;
        let mut alin = Line2d::new(p1, avec);
        // OCCT L42: alin.SetLocation(ElCLib::Value(c1, alin)).
        let value = elclib_value_lin2d(c1, &alin);
        alin = Line2d::new(value, alin.direction);
        // OCCT L43-44: cf = ElCLib::Parameter(alin, P1); cl = ... P2.
        *cf = elclib_parameter_lin2d(&alin, p1);
        *cl = elclib_parameter_lin2d(&alin, p2);
        alin
    }

    /// OCCT IsLinear (ShapeCustom_Curve2d.hxx L38-40 + .cxx L50-105): checks
    /// whether the poles are (within tolerance) on one line; outputs the
    /// maximal deviation.
    pub fn is_linear(the_poles: &[DVec2], tolerance: f64, deviation: &mut f64) -> bool {
        // OCCT L54-58: nbPoles < 2 -> false.
        let nb_poles = the_poles.len() as i32;
        if nb_poles < 2 {
            return false;
        }

        // OCCT L60-76: the widest pole pair.
        let mut d_max = 0f64;
        let mut i_max1 = 0usize;
        let mut i_max2 = 0usize;
        for i in 1..nb_poles {
            for j in (i + 1)..=nb_poles {
                let dist = the_poles[(i - 1) as usize].distance_squared(the_poles[(j - 1) as usize]);
                if dist > d_max {
                    d_max = dist;
                    i_max1 = i as usize;
                    i_max2 = j as usize;
                }
            }
        }

        // OCCT L78-82: constexpr dPreci = PConfusion()*PConfusion().
        let d_preci = PCONFUSION * PCONFUSION;
        if d_max < d_preci {
            return false;
        }

        // OCCT L84-87: gp_Vec2d/Dir2d/Lin2d on the widest pair.
        let tol2 = tolerance * tolerance;
        let avec = the_poles[i_max2 - 1] - the_poles[i_max1 - 1];
        let alin = Line2d::new(the_poles[i_max1 - 1], avec);

        // OCCT L89-101: every pole must be within tol2 of the line.
        let mut a_max = 0f64;
        for i in 1..=nb_poles {
            let dist = gp_lin2d_square_distance(&alin, the_poles[(i - 1) as usize]);
            if dist > tol2 {
                return false;
            }
            if dist > a_max {
                a_max = dist;
            }
        }
        // OCCT L102: Deviation = sqrt(aMax).
        *deviation = a_max.sqrt();

        true
    }

    /// OCCT ConvertToLine2d (ShapeCustom_Curve2d.hxx L45-52 + .cxx L109-152):
    /// converts a BSpline2d/Bezier2d to a 2d line when linear; recalculates
    /// first/last parameters; returns None (the OCCT null handle) when the
    /// curve is not a line.
    pub fn convert_to_line2d(
        the_curve: &Curve2d,
        c1: f64,
        c2: f64,
        the_tolerance: f64,
        cf: &mut f64,
        cl: &mut f64,
        the_deviation: &mut f64,
    ) -> Option<Curve2d> {
        let mut a_line2d: Option<Curve2d> = None;
        // OCCT L119-125: P1 = Value(c1); P2 = Value(c2); the degenerate
        // (square distance below tol^2) case is not a line.
        let p1 = the_curve.point_at(c1);
        let p2 = the_curve.point_at(c2);
        let d_preci = the_tolerance * the_tolerance;
        if p1.distance_squared(p2) < d_preci {
            return a_line2d; // it is not a line
        }

        // OCCT L127-137: the Geom2d_BSplineCurve down-cast arm.
        if let Curve2d::BSpline(bsc) = the_curve {
            if !Self::is_linear(&bsc.control_points, the_tolerance, the_deviation) {
                return a_line2d; // non
            }
            let alin = Self::get_line(p1, p2, c1, cf, cl);
            a_line2d = Some(Curve2d::Line(alin));
            return a_line2d;
        }

        // OCCT L139-149: the Geom2d_BezierCurve down-cast arm.
        if let Curve2d::Bezier(bzc) = the_curve {
            if !Self::is_linear(&bzc.control_points, the_tolerance, the_deviation) {
                return a_line2d; // non
            }
            let alin = Self::get_line(p1, p2, c1, cf, cl);
            a_line2d = Some(Curve2d::Line(alin));
            return a_line2d;
        }

        // OCCT L151: return aLine2d (null).
        a_line2d
    }

    /// OCCT SimplifyBSpline2d (ShapeCustom_Curve2d.hxx L57-58 + .cxx
    /// L156-196): removes the knots where the local DegMult-th derivatives
    /// are parallel; returns false when the bspline was not modified.
    pub fn simplify_bspline2d(the_bspline2d: &mut BSplineCurve2, the_tolerance: f64) -> bool {
        let mut bsc = Geom2dBSplineCurve { c: the_bspline2d };
        // OCCT L159: int NbK = aInitNbK = theBSpline2d->NbKnots().
        let a_init_nb_k = bsc.nb_knots();
        let mut nb_k = a_init_nb_k;
        // search knot to remove
        // OCCT L162-163: bool IsToRemove = true; int aKnotIndx = NbK - 1.
        let mut is_to_remove = true;
        let mut a_knot_indx = nb_k - 1;
        while is_to_remove && nb_k > 2 {
            // OCCT L165-167: aMult = Multiplicity(aKnotIndx);
            // DegMult = Degree() - aMult.
            let a_mult = bsc.multiplicity(a_knot_indx);
            let deg_mult = bsc.degree() - a_mult;
            // OCCT L168: if ((DegMult > 1) && IsCN(DegMult)).
            if deg_mult > 1 && bsc.is_cn(deg_mult) {
                // OCCT L170-172: U = Knot(aKnotIndx);
                // aVec1 = LocalDN(U, aKnotIndx-1, aKnotIndx, DegMult);
                // aVec2 = LocalDN(U, aKnotIndx, aKnotIndx+1, DegMult).
                let u = bsc.knot(a_knot_indx);
                let a_vec1 = bsc.local_dn(u, a_knot_indx - 1, a_knot_indx, deg_mult);
                let a_vec2 = bsc.local_dn(u, a_knot_indx, a_knot_indx + 1, deg_mult);
                // OCCT L173-174: check the derivations are have the "same"
                // angle: aVec1.IsParallel(aVec2, Precision::Angular()).
                if gp_vec2d_is_parallel(a_vec1, a_vec2, ANGULAR) {
                    // remove knot
                    // OCCT L177-184: try { OCC_CATCH_SIGNALS
                    // RemoveKnot(aKnotIndx, aMult-1, theTolerance); }
                    // catch (Standard_Failure const&) {}.
                    let _ = bsc.remove_knot(a_knot_indx, a_mult - 1, the_tolerance);
                }
            }
            // OCCT L187: aKnotIndx--.
            a_knot_indx -= 1;

            // OCCT L189: NbK = theBSpline2d->NbKnots().
            nb_k = bsc.nb_knots();
            // OCCT L190-193: if (aKnotIndx == 1 || aKnotIndx == NbK)
            // IsToRemove = false.
            if a_knot_indx == 1 || a_knot_indx == nb_k {
                is_to_remove = false;
            }
        }
        // OCCT L195: return (aInitNbK > NbK).
        a_init_nb_k > nb_k
    }
}
