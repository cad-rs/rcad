// OCCT Convert_CompBezierCurvesToBSplineCurve (TKMath/Convert) — 1:1 Rust
// translation of:
// - Convert_CompBezierCurvesToBSplineCurve.hxx L33-45 + .cxx L21-25 (the
//   gp_Pnt/gp_Vec instantiation)
// - Convert_CompBezierCurvesToBSplineCurveBase.hxx L29-222 (the generic
//   base class template; rcad carries the 3D instantiation directly)
//
// The degree elevation delegates to the same OCCT call chain:
//   Convert Base L82-86  BSplCLib::IncreaseDegree(NewDegree, Poles,
//                        NoWeights, aPoints, NoWeights)
//   -> BSplCLib_BzSyntaxes.cxx L32-43 (the gp_Pnt overload)
//   -> BSplCLib_CurveComputation.pxx L2080-2099 (BSplCLib_IncreaseDegree_
//      Bezier: the Bezier knots {0,1} / mults {Degree+1,Degree+1})
//   -> BSplCLib.cxx L2592 (the flat-pole engine == math::bspl_lib::
//      increase_degree, already translated).
//
// Architecture difference (the 2d instantiation): the OCCT base is a C++
// template over PointType/VecType (gp_Pnt/gp_Pnt2d); the `if constexpr
// (std::is_same_v<PointType, gp_Pnt>)` branch (the Base.hxx L120-135) is
// compile-time selected for the 3D instantiation carried here.  The 2D
// variant (Convert_CompBezierCurves2dToBSplineCurve) keeps the plain
// knot-value continuation and will re-use this file when its consumers land.

use glam::DVec3;

use crate::core::precision::REAL_LAST;
use crate::math::bspl_lib;

/// OCCT Standard::Epsilon(theValue) (Standard_Real.hxx L242-250) — the
/// distance to the nearest neighbouring double in the direction of infinity
/// with the same sign as theValue.
fn standard_epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        let next = if the_value == REAL_LAST {
            REAL_LAST
        } else {
            f64::from_bits(the_value.to_bits() + 1)
        };
        next - the_value
    } else {
        the_value - f64::from_bits(the_value.to_bits() - 1)
    }
}

/// OCCT gp::Resolution() == RealSmall() (gp.hxx L59-60) — the smallest
/// positive normalized double (the app_def_compute.rs precedent).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gp_Dir::Angle (gp_Dir.cxx L27-50) — the angle in [0, PI] between
/// the (unit) directions, computed with the acos/asin switch at 45 degrees.
fn gp_dir_angle(coord: DVec3, other: DVec3) -> f64 {
    let cosinus = coord.dot(other);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = coord.cross(other).length();
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-493) — raises
/// VectorWithNullMagnitude when a magnitude <= gp::Resolution(); the rcad
/// mapping of the raise is a panic.
fn gp_vec_angle(coord: DVec3, other: DVec3) -> f64 {
    assert!(
        coord.length() > GP_RESOLUTION && other.length() > GP_RESOLUTION,
        "gp_VectorWithNullMagnitude: gp_Vec::Angle"
    );
    // gp_Dir(coord) normalizes the coordinates.
    gp_dir_angle(coord.normalize(), other.normalize())
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L142-146).
fn gp_vec_is_parallel(coord: DVec3, other: DVec3, angular_tolerance: f64) -> bool {
    let an_ang = gp_vec_angle(coord, other);
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT BSplCLib_IncreaseDegree_Bezier (BSplCLib_CurveComputation.pxx
/// L2080-2099) — the Bezier (knots {0,1}, mults {Degree+1,Degree+1}) degree
/// elevation of a flat pole array through math::bspl_lib::increase_degree.
fn bspl_clib_increase_degree_bezier(
    new_degree: usize,
    poles: &[DVec3],
    new_poles: &mut [DVec3],
) {
    let a_degree = poles.len() - 1;
    // OCCT: BSplCLib_KnotArrays<2> aBezierKnots(aDegree).
    let knots = [0.0f64, 1.0];
    let mults = [(a_degree + 1) as i32, (a_degree + 1) as i32];

    // OCCT: PLib::SetPoles(Poles, poles) — the flat (x,y,z) pole buffer.
    let mut flat = vec![0.0f64; 3 * poles.len()];
    for (i, p) in poles.iter().enumerate() {
        flat[i * 3] = p.x;
        flat[i * 3 + 1] = p.y;
        flat[i * 3 + 2] = p.z;
    }
    let mut new_flat = vec![0.0f64; 3 * new_poles.len()];

    // OCCT: BSplCLib::IncreaseDegree(Degree, NewDegree, false, dim, poles,
    // Knot, Mult, newpoles, Knot, Mult).
    let mut new_knots = [0.0f64; 2];
    let mut new_mults = [0i32; 2];
    bspl_lib::increase_degree(
        a_degree,
        new_degree,
        false,
        3,
        &flat,
        &knots,
        &mults,
        &mut new_flat,
        &mut new_knots,
        &mut new_mults,
    );

    // OCCT: PLib::GetPoles(newpoles, NewPoles).
    for (i, p) in new_poles.iter_mut().enumerate() {
        *p = DVec3::new(new_flat[i * 3], new_flat[i * 3 + 1], new_flat[i * 3 + 2]);
    }
}

/// OCCT Convert_CompBezierCurvesToBSplineCurveBase<gp_Pnt, gp_Vec>
/// (the Base.hxx L29-222) — an algorithm to convert a sequence of adjacent
/// non-rational Bezier curves into a BSpline curve.
pub struct ConvertCompBezierCurvesToBSplineCurve {
    my_sequence: Vec<Vec<DVec3>>, // OCCT: mySequence (Base.hxx L214)
    my_curve_poles: Vec<DVec3>,   // OCCT: myCurvePoles (Base.hxx L215)
    my_curve_knots: Vec<f64>,     // OCCT: myCurveKnots (Base.hxx L216)
    my_knots_mults: Vec<i32>,     // OCCT: myKnotsMults (Base.hxx L217)
    my_degree: i32,               // OCCT: myDegree (Base.hxx L218)
    my_angular: f64,              // OCCT: myAngular (Base.hxx L219)
}

impl Default for ConvertCompBezierCurvesToBSplineCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl ConvertCompBezierCurvesToBSplineCurve {
    /// OCCT constructor (the Base.hxx L37-41; the .cxx L21-25) —
    /// theAngularTolerance defaults to 1.0e-4.
    pub fn new() -> Self {
        Self::new_with_angular_tolerance(1.0e-4)
    }

    /// OCCT constructor with an explicit angular tolerance.
    pub fn new_with_angular_tolerance(the_angular_tolerance: f64) -> Self {
        ConvertCompBezierCurvesToBSplineCurve {
            my_sequence: Vec::new(),
            my_curve_poles: Vec::new(),
            my_curve_knots: Vec::new(),
            my_knots_mults: Vec::new(),
            my_degree: 0,
            my_angular: the_angular_tolerance,
        }
    }

    /// OCCT AddCurve(thePoles) (the Base.hxx L46) — adds the Bezier curve
    /// defined by the table of poles to the sequence of adjacent Bezier
    /// curves to be converted.
    pub fn add_curve(&mut self, the_poles: &[DVec3]) {
        self.my_sequence.push(the_poles.to_vec());
    }

    /// OCCT Perform() (the Base.hxx L50-173) — computes all the data needed
    /// to build a BSpline curve equivalent to the adjacent Bezier curve
    /// sequence.
    pub fn perform(&mut self) {
        self.my_curve_poles.clear();
        self.my_curve_knots.clear();
        self.my_knots_mults.clear();
        if self.my_sequence.is_empty() {
            return;
        }
        let a_lower_i = 1usize; // OCCT: const int aLowerI = 1
        let an_upper_i = self.my_sequence.len(); // OCCT: mySequence.Length()
        let a_nbr_curv = an_upper_i - a_lower_i + 1;
        // OCCT: NCollection_Array1<double> aCurveKnVals(1, aNbrCurv).
        let mut a_curve_kn_vals = vec![0.0f64; a_nbr_curv];

        self.my_degree = 0;
        for i in 0..self.my_sequence.len() {
            self.my_degree = self.my_degree.max(self.my_sequence[i].len() as i32 - 1);
        }

        let mut a_det = 0.0f64;
        let mut a_p1 = DVec3::ZERO;
        // OCCT: PointType aP1, aP2, aP3 — aP2/aP3 are loop-local reads; aP1
        // carries the previous segment pole across iterations.
        let a_max_degree = self.my_degree;
        // OCCT: NCollection_Array1<PointType> aPoints(1, myDegree + 1).
        let mut a_points = vec![DVec3::ZERO; (self.my_degree + 1) as usize];

        for i in a_lower_i..=an_upper_i {
            // 1- Raise the Bezier curve to the maximum degree.
            // OCCT L78-79: aDeg = mySequence(i).Length() - 1; anInc =
            // myDegree - aDeg.
            let a_deg = self.my_sequence[i - 1].len() as i32 - 1;
            let an_inc = self.my_degree - a_deg;
            if an_inc > 0 {
                // OCCT L82-86: BSplCLib::IncreaseDegree(myDegree,
                // mySequence(i), NoWeights, aPoints, NoWeights).
                bspl_clib_increase_degree_bezier(
                    self.my_degree as usize,
                    &self.my_sequence[i - 1],
                    &mut a_points,
                );
            } else {
                // OCCT L90: aPoints = mySequence(i).
                a_points.copy_from_slice(&self.my_sequence[i - 1]);
            }

            // 2- Process the node of junction between 2 Bezier curves.
            if i == a_lower_i {
                // Processing of the initial node of the BSpline.
                // OCCT L97-100.
                for j in 1..=a_max_degree {
                    self.my_curve_poles.push(a_points[(j - 1) as usize]);
                }
                a_curve_kn_vals[1 - 1] = 1.0; // To begin the series.
                self.my_knots_mults.push(a_max_degree + 1);
                a_det = 1.0;
            }

            if i != a_lower_i {
                // OCCT L108-110.
                let a_p2 = a_points[0];
                let a_p3 = a_points[1];
                let a_v1 = a_p2 - a_p1;
                let a_v2 = a_p3 - a_p2;

                // Processing of the tangency between Bezier and the previous.
                // This allows to guarantee at least a C1 continuity if the
                // tangents are coherent.
                // OCCT L114-117.
                let a_d1 = a_v1.length_squared();
                let a_d2 = a_v2.length_squared();
                let is_parallel = a_max_degree > 1
                    && a_d1 > GP_RESOLUTION
                    && a_d2 > GP_RESOLUTION
                    && gp_vec_is_parallel(a_v1, a_v2, self.my_angular);
                if is_parallel {
                    // OCCT L119: aLambda = sqrt(aD2 / aD1).
                    let a_lambda = (a_d2 / a_d1).sqrt();
                    // The 3D instantiation (gp_Pnt) keeps the epsilon guard
                    // against numerically small knot values (the Base.hxx
                    // L120-135); the 2d branch skips the guard.
                    if a_curve_kn_vals[i - 2] * a_lambda > 10.0 * standard_epsilon(a_det) {
                        self.my_knots_mults.push(a_max_degree - 1);
                        a_curve_kn_vals[i - 1] = a_curve_kn_vals[i - 2] * a_lambda;
                    } else {
                        self.my_curve_poles.push(a_points[0]);
                        self.my_knots_mults.push(a_max_degree);
                        a_curve_kn_vals[i - 1] = 1.0;
                    }
                } else {
                    // OCCT L144-146.
                    self.my_curve_poles.push(a_points[0]);
                    self.my_knots_mults.push(a_max_degree);
                    a_curve_kn_vals[i - 1] = 1.0;
                }
                // OCCT L148: aDet += aCurveKnVals(i).
                a_det += a_curve_kn_vals[i - 1];

                // Store the poles.
                // OCCT L151-154.
                for j in 2..=a_max_degree {
                    self.my_curve_poles.push(a_points[(j - 1) as usize]);
                }
            }

            if i == an_upper_i {
                // Processing of the end node of the BSpline.
                // OCCT L160-161.
                self.my_curve_poles
                    .push(a_points[a_max_degree as usize]);
                self.my_knots_mults.push(a_max_degree + 1);
            }
            // OCCT L163: aP1 = aPoints(aMaxDegree).
            a_p1 = a_points[(a_max_degree - 1) as usize];
        }

        // Correct nodal values to make them variable within [0.,1.].
        // OCCT L167-172.
        self.my_curve_knots.push(0.0);
        for i in 2..=a_nbr_curv {
            let prev = self.my_curve_knots[i - 2];
            self.my_curve_knots.push(prev + (a_curve_kn_vals[i - 2] / a_det));
        }
        self.my_curve_knots.push(1.0);
    }

    /// OCCT Degree() (the Base.hxx L176) — returns the degree of the BSpline
    /// curve.
    pub fn degree(&self) -> i32 {
        self.my_degree
    }

    /// OCCT NbPoles() (the Base.hxx L179).
    pub fn nb_poles(&self) -> i32 {
        self.my_curve_poles.len() as i32
    }

    /// OCCT Poles(thePoles) (the Base.hxx L183-190) — loads the Poles table
    /// with the poles of the BSpline curve.
    pub fn poles(&self, the_poles: &mut Vec<DVec3>) {
        let mut k = 0usize;
        for p in the_poles.iter_mut() {
            *p = self.my_curve_poles[k];
            k += 1;
        }
    }

    /// OCCT NbKnots() (the Base.hxx L193).
    pub fn nb_knots(&self) -> i32 {
        self.my_curve_knots.len() as i32
    }

    /// OCCT KnotsAndMults(theKnots, theMults) (the Base.hxx L199-211).
    pub fn knots_and_mults(&self, the_knots: &mut Vec<f64>, the_mults: &mut Vec<i32>) {
        for (i, knot) in the_knots.iter_mut().enumerate() {
            *knot = self.my_curve_knots[i];
        }
        for (i, mult) in the_mults.iter_mut().enumerate() {
            *mult = self.my_knots_mults[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64) -> DVec3 {
        DVec3::new(x, y, 0.0)
    }

    /// Concat of a degree-2 Bezier and a non-parallel degree-3 Bezier
    /// (Base.hxx L106-155 junction processing, the `else` C0 branch):
    /// the result is a degree-3 BSpline with mults {4,3,4}, knots
    /// {0, 1/2, 1} and 7 poles (3 shared junction + segment 2 poles).
    #[test]
    fn concat_two_beziers_c0_junction() {
        let mut conv = ConvertCompBezierCurvesToBSplineCurve::new();
        // Degree-2 segment: poles (0,0), (1,1), (2,0).
        conv.add_curve(&[p(0.0, 0.0), p(1.0, 1.0), p(2.0, 0.0)]);
        // Degree-3 segment starting at (2,0) with a non-parallel tangent.
        conv.add_curve(&[p(2.0, 0.0), p(3.0, 2.0), p(4.0, -2.0), p(5.0, 0.0)]);
        conv.perform();

        assert_eq!(conv.degree(), 3);
        assert_eq!(conv.nb_poles(), 7);
        assert_eq!(conv.nb_knots(), 3);

        let mut knots = vec![0.0; conv.nb_knots() as usize];
        let mut mults = vec![0; conv.nb_knots() as usize];
        conv.knots_and_mults(&mut knots, &mut mults);
        assert!(
            knots.iter().zip([0.0, 0.5, 1.0]).all(|(a, b)| (a - b).abs() < 1e-15),
            "knots {:?}",
            knots
        );
        assert_eq!(mults, vec![4, 3, 4]);

        let mut poles = vec![DVec3::ZERO; conv.nb_poles() as usize];
        conv.poles(&mut poles);
        // Segment 1 raised 2 -> 3: Q1 = P0/3 + 2 P1/3, Q2 = 2 P1/3 + P2/3.
        let expected = [
            p(0.0, 0.0),
            p(2.0 / 3.0, 2.0 / 3.0),
            p(4.0 / 3.0, 2.0 / 3.0),
            p(2.0, 0.0),
            p(3.0, 2.0),
            p(4.0, -2.0),
            p(5.0, 0.0),
        ];
        for (i, want) in expected.iter().enumerate() {
            assert!(
                (poles[i] - *want).length() < 1e-12,
                "pole {}: {} vs {}",
                i + 1,
                poles[i],
                want
            );
        }
    }

    /// Concat of two Beziers with parallel end tangents (the C1 branch,
    /// Base.hxx L116-135 with the lambda guard taken): the junction poles
    /// are NOT appended (the C1 point is implied by the neighbouring poles
    /// and the knot ratio lambda = |V2| / |V1| = 3/2), the junction
    /// multiplicity drops to degree - 1 = 2, and the knot value is
    /// knVals(1) / (knVals(1) + knVals(2)) = 1 / 2.5 = 0.4.
    #[test]
    fn concat_two_beziers_c1_junction() {
        let mut conv = ConvertCompBezierCurvesToBSplineCurve::new();
        // Segment 1: end tangent Q3 - Q2 = (2/3, 2/3) (after the degree-2 ->
        // 3 elevation of poles (0,0), (1,-1), (2,0)).
        conv.add_curve(&[p(0.0, 0.0), p(1.0, -1.0), p(2.0, 0.0)]);
        // Segment 2: first tangent (1, 1) — parallel to (2/3, 2/3), length
        // ratio lambda = sqrt(2) / ((2/3) sqrt(2)) = 3/2.
        conv.add_curve(&[p(2.0, 0.0), p(3.0, 1.0), p(4.0, -1.0), p(5.0, 0.0)]);
        conv.perform();

        assert_eq!(conv.degree(), 3);
        // maxDegree leading poles + the (maxDegree - 1) non-junction poles
        // of segment 2 + the trailing pole = 3 + 2 + 1 = 6.
        assert_eq!(conv.nb_poles(), 6);
        assert_eq!(conv.nb_knots(), 3);

        let mut knots = vec![0.0; conv.nb_knots() as usize];
        let mut mults = vec![0; conv.nb_knots() as usize];
        conv.knots_and_mults(&mut knots, &mut mults);
        assert!(
            knots.iter().zip([0.0, 0.4, 1.0]).all(|(a, b)| (a - b).abs() < 1e-15),
            "knots {:?}",
            knots
        );
        assert_eq!(mults, vec![4, 2, 4]);

        let mut poles = vec![DVec3::ZERO; conv.nb_poles() as usize];
        conv.poles(&mut poles);
        let expected = [
            p(0.0, 0.0),
            p(2.0 / 3.0, -2.0 / 3.0),
            p(4.0 / 3.0, -2.0 / 3.0),
            p(3.0, 1.0),
            p(4.0, -1.0),
            p(5.0, 0.0),
        ];
        for (i, want) in expected.iter().enumerate() {
            assert!(
                (poles[i] - *want).length() < 1e-12,
                "pole {}: {} vs {}",
                i + 1,
                poles[i],
                want
            );
        }
    }

    /// OCCT Standard::Epsilon sanity (Standard_Real.hxx L242-250).
    #[test]
    fn standard_epsilon_matches_occt() {
        assert!((standard_epsilon(1.0) - f64::EPSILON).abs() < 1e-30);
        assert!(standard_epsilon(0.0) > 0.0);
    }

    /// Mirror of the OCCT GTest
    /// Convert_CompBezierCurvesToBSplineCurveTest.TwoAdjacentBeziers_C1
    /// (TKMath/GTests, L93-122): two cubic Beziers with parallel unit
    /// tangents at the junction — mults(2) == degree - 1.
    #[test]
    fn occt_gtest_two_adjacent_beziers_c1() {
        let mut conv = ConvertCompBezierCurvesToBSplineCurve::new();
        conv.add_curve(&[p(0.0, 0.0), p(1.0, 1.0), p(2.0, 1.0), p(3.0, 0.0)]);
        conv.add_curve(&[p(3.0, 0.0), p(4.0, -1.0), p(5.0, -1.0), p(6.0, 0.0)]);
        conv.perform();

        assert_eq!(conv.degree(), 3);
        assert_eq!(conv.nb_knots(), 3);

        let mut knots = vec![0.0; conv.nb_knots() as usize];
        let mut mults = vec![0; conv.nb_knots() as usize];
        conv.knots_and_mults(&mut knots, &mut mults);
        assert_eq!(mults[0], conv.degree() + 1);
        assert_eq!(mults[1], conv.degree() - 1);
        assert_eq!(mults[2], conv.degree() + 1);
        // lambda = 1: the junction knot value is 1/2.
        assert!((knots[1] - 0.5).abs() < 1e-15);
        // The output must interpolate the junction point (3,0) at t = 0.5.
        let poles_out = {
            let mut v = vec![DVec3::ZERO; conv.nb_poles() as usize];
            conv.poles(&mut v);
            v
        };
        // De Boor at the mult-2 knot: the junction point is the average of
        // the two adjacent poles (knot ratio 1).
        let junction = (poles_out[2] + poles_out[3]) * 0.5;
        assert!((junction - p(3.0, 0.0)).length() < 1e-12);
    }
}
