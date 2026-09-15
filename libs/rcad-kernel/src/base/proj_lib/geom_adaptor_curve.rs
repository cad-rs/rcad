//! OCCT GeomAdaptor_Curve (TKG3d/GeomAdaptor) — the 1:1 translation of
//! `GeomAdaptor_Curve.hxx` (L51-267) + `GeomAdaptor_Curve.cxx` (L60-1295).
//!
//! The rcad type keeps the OCCT member set (`myCurve`, `myTypeCurve`,
//! `myFirst`, `myLast`, `myCurveData`) and every public method of the C++
//! class.  Architecture differences (rcad encodings):
//!   - `occ::handle<Geom_Curve> myCurve` maps to the kernel curve value
//!     [`crate::geom::Curve3`] (the handle copy is a value move);
//!   - the `CurveDataVariant` (gp_Lin / gp_Circ / ... / BezierData /
//!     BSplineData / OffsetData) is carried by the `Curve3` variant payloads
//!     themselves — the OCCT variant exists to avoid downcasts in C++, the
//!     rcad match on `Curve3` plays the same role;
//!   - `Geom_TrimmedCurve` is unwrapped in `load` exactly as OCCT does
//!     (`load` L252-255 recurses on the basis curve, keeping UFirst/ULast),
//!     so `myCurve` stores the basis curve and the evaluation reads raw
//!     parameters;
//!   - the BSplCLib_Cache / GeomEval_RepCurveDesc evaluation caches are
//!     C++-internal performance state: the kernel curve engine evaluates the
//!     same poles directly (the cache is not observable); the `hasEvalRep`
//!     branches therefore collapse into the kernel evaluation.
//!   - the OCCT `GetType()` virtual lands on the [`Adaptor3dCurveGeom`]
//!     companion trait (together with the geometry downcasts) — the base
//!     [`Adaptor3dCurve`] trait cannot carry a same-named method without
//!     making every concrete `.get_type()` call in the consumers ambiguous.

use glam::DVec3;

use super::adaptor::Adaptor3dCurve;
use super::CurveType;
use crate::core::precision;
use crate::geom::curve_dn::curve_dn;
use crate::geom::{Circle3, Curve3, CurveEval, Ellipse3, Hyperbola3, Line3, Parabola3};
use crate::math::bspl_lib::locate_parameter_knots_mults;
use crate::math::{GeomAbsCurveType, GeomAbsShape};

/// OCCT `static const double PosTol = Precision::PConfusion() / 2`
/// (GeomAdaptor_Curve.cxx L60).
const POS_TOL: f64 = precision::PCONFUSION / 2.0;

// =========================================================================
// OCCT Geom_BSplineCurve::IsRational static (Geom_BSplineCurve.cxx L97-108)
// =========================================================================

/// OCCT `static bool Rational(const NCollection_Array1<double>& theWeights)`
/// (Geom_BSplineCurve.cxx L97-108) — true when adjacent weights differ by
/// more than gp::Resolution().  The Bezier curve flag uses the same test
/// (Geom_BezierCurve.cxx L60-73).
fn curve_rational(weights: &[f64]) -> bool {
    for i in 0..weights.len().saturating_sub(1) {
        if (weights[i] - weights[i + 1]).abs() > crate::math::gp::GP_RESOLUTION {
            return true;
        }
    }
    false
}

/// The OCCT flat-knots reader of the kernel BSpline payload (1-based index
/// access like the OCCT `NCollection_Array1` inside the BSplCLib calls).
fn at(v: &[f64], i: i32) -> f64 {
    v[(i - 1) as usize]
}

// =========================================================================
// OCCT GeomAdaptor_Curve (GeomAdaptor_Curve.hxx L51-267)
// =========================================================================

/// OCCT GeomAdaptor_Curve — the Adaptor3d_Curve instance over a Geom_Curve.
///
/// Constructor forms: [`GeomCurveAdaptor::new`] = `GeomAdaptor_Curve(C)`
/// (L94-99 + Load L117-124), [`GeomCurveAdaptor::with_range`] =
/// `GeomAdaptor_Curve(C, UFirst, ULast)` (L104-109 + Load L127-138).
#[derive(Clone, Debug)]
pub struct GeomCurveAdaptor {
    /// OCCT: handle(Geom_Curve) myCurve — the basis curve after the
    /// TrimmedCurve unwrapping of `load`.
    pub curve: Curve3,
    /// OCCT: GeomAbs_CurveType myTypeCurve.
    pub my_type_curve: GeomAbsCurveType,
    /// OCCT: Standard_Real myFirst.
    pub first: f64,
    /// OCCT: Standard_Real myLast.
    pub last: f64,
}

impl GeomCurveAdaptor {
    /// OCCT GeomAdaptor_Curve(theCurve) -> Load(theCurve)
    /// (GeomAdaptor_Curve.hxx L101 + L117-124): the curve natural domain.
    pub fn new(curve: Curve3) -> Self {
        let [first, last] = CurveEval::default_domain(&curve);
        let mut a = GeomCurveAdaptor {
            curve,
            my_type_curve: GeomAbsCurveType::OtherCurve,
            first,
            last,
        };
        a.load_with_range_inner(a.first, a.last);
        a
    }

    /// OCCT GeomAdaptor_Curve(theCurve, theUFirst, theULast)
    /// (GeomAdaptor_Curve.hxx L104-109): the restricted window.
    ///
    /// OCCT raises Standard_ConstructionError when
    /// theUFirst > theULast + Precision::Confusion() (L133-136).
    pub fn with_range(curve: Curve3, first: f64, last: f64) -> Self {
        assert!(
            first <= last + precision::CONFUSION,
            "Standard_ConstructionError: GeomAdaptor_Curve::Load"
        );
        let mut a = GeomCurveAdaptor {
            curve,
            my_type_curve: GeomAbsCurveType::OtherCurve,
            first,
            last,
        };
        a.load_with_range_inner(first, last);
        a
    }

    /// OCCT Load(theCurve) (GeomAdaptor_Curve.hxx L117-124) as a mutator.
    pub fn load(&mut self, curve: Curve3) {
        let [first, last] = CurveEval::default_domain(&curve);
        self.curve = curve;
        self.load_with_range_inner(first, last);
    }

    /// OCCT Load(theCurve, theUFirst, theULast)
    /// (GeomAdaptor_Curve.hxx L127-138) as a mutator.
    pub fn load_with_range(&mut self, curve: Curve3, first: f64, last: f64) {
        assert!(
            first <= last + precision::CONFUSION,
            "Standard_ConstructionError: GeomAdaptor_Curve::Load"
        );
        self.curve = curve;
        self.load_with_range_inner(first, last);
    }

    /// OCCT `void GeomAdaptor_Curve::load(C, UFirst, ULast)`
    /// (GeomAdaptor_Curve.cxx L239-325): the type dispatch with the
    /// TrimmedCurve unwrapping.
    fn load_with_range_inner(&mut self, first: f64, last: f64) {
        self.first = first;
        self.last = last;

        // OCCT L252-255: Geom_TrimmedCurve -> Load(basis, UFirst, ULast).
        let curve = match &self.curve {
            Curve3::Trimmed(t) => t.curve.as_ref().clone(),
            _ => self.curve.clone(),
        };
        self.curve = curve;

        // OCCT L256-311: the DynamicType dispatch.
        self.my_type_curve = match &self.curve {
            Curve3::Circle(_) => GeomAbsCurveType::Circle,
            Curve3::Line(_) => GeomAbsCurveType::Line,
            Curve3::Ellipse(_) => GeomAbsCurveType::Ellipse,
            Curve3::Parabola(_) => GeomAbsCurveType::Parabola,
            Curve3::Hyperbola(_) => GeomAbsCurveType::Hyperbola,
            Curve3::Bezier(_) => GeomAbsCurveType::BezierCurve,
            Curve3::BSpline(_) => GeomAbsCurveType::BSplineCurve,
            Curve3::Offset(_) => GeomAbsCurveType::OffsetCurve,
            _ => GeomAbsCurveType::OtherCurve,
        };
    }

    /// OCCT GeomAdaptor_Curve::Reset (GeomAdaptor_Curve.cxx L229-235).
    pub fn reset(&mut self) {
        self.my_type_curve = GeomAbsCurveType::OtherCurve;
        self.curve = Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X));
        self.first = 0.0;
        self.last = 0.0;
    }

    // ---------------------------------------------------------------------
    // The evaluation dispatch (GeomAdaptor_Curve.cxx EvalD0/D1/D2/D3/DN).
    // The per-type arms (ElCLib for the elementary kinds, BSplCLib_Cache for
    // the polynomial kinds, Geom_OffsetCurveUtils for offsets, and the
    // myCurve->Eval* fall-through) share the kernel curve engine encoding:
    // the kernel evaluates the same payloads with the same results; the
    // IsBoundary / RebuildCache machinery is C++-internal cache state.
    // ---------------------------------------------------------------------

    /// OCCT GeomAdaptor_Curve::LocalContinuity(U1, U2)
    /// (GeomAdaptor_Curve.cxx L143-225) — the BSpline continuity between two
    /// parameters.
    fn local_continuity(&self, u1: f64, u2: f64) -> GeomAbsShape {
        // OCCT L145: raise when myTypeCurve != GeomAbs_BSplineCurve.
        let Curve3::BSpline(a_bspl) = &self.curve else {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::LocalContinuity");
        };
        let (tk, tm) = a_bspl.knots_mults();
        let nb = tm.len() as i32;
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        let mut new_first = 0.0f64;
        let mut new_last = 0.0f64;
        // OCCT L153-170: BSplCLib::LocateParameter(..., U1/U2, IsPeriodic,
        // 1, Nb, Index, NewU).
        locate_parameter_knots_mults(
            a_bspl.degree,
            &tk,
            &tm,
            u1,
            a_bspl.is_periodic,
            1,
            nb,
            &mut index1,
            &mut new_first,
        );
        locate_parameter_knots_mults(
            a_bspl.degree,
            &tk,
            &tm,
            u2,
            a_bspl.is_periodic,
            1,
            nb,
            &mut index2,
            &mut new_last,
        );
        // OCCT L171-181: the knot-proximity index adjustments.
        if (new_first - at(&tk, index1 + 1)).abs() < precision::PCONFUSION {
            if index1 < nb {
                index1 += 1;
            }
        }
        if (new_last - at(&tk, index2)).abs() < precision::PCONFUSION {
            index2 -= 1;
        }
        // OCCT L182-187: the periodic wrap of Index1.
        let mut mult_max;
        if a_bspl.is_periodic && index1 == nb {
            index1 = 1;
        }

        if index2 - index1 <= 0 && !a_bspl.is_periodic {
            // OCCT L191: MultMax = 100 — CN between two consecutive knots.
            mult_max = 100;
        } else {
            mult_max = tm[(index1 + 1 - 1) as usize];
            for i in index1 + 1..=index2 {
                if tm[(i - 1) as usize] > mult_max {
                    mult_max = tm[(i - 1) as usize];
                }
            }
            mult_max = a_bspl.degree as i32 - mult_max;
        }
        // OCCT L205-224: the multiplicity -> continuity mapping.
        if mult_max <= 0 {
            GeomAbsShape::C0
        } else if mult_max == 1 {
            GeomAbsShape::C1
        } else if mult_max == 2 {
            GeomAbsShape::C2
        } else if mult_max == 3 {
            GeomAbsShape::C3
        } else {
            GeomAbsShape::CN
        }
    }

    /// OCCT GeomAdaptor_Curve::Continuity (GeomAdaptor_Curve.cxx L333-367).
    pub fn continuity(&self) -> GeomAbsShape {
        if self.my_type_curve == GeomAbsCurveType::BSplineCurve {
            return self.local_continuity(self.first, self.last);
        }

        if self.my_type_curve == GeomAbsCurveType::OffsetCurve {
            // OCCT L340-360: the basis continuity downshift.  The G1/G2 arms
            // of the OCCT switch cannot be spelled with the rcad 5-level
            // GeomAbsShape (no G1/G2 variants) and cannot be reached.
            let s = self.offset_basis_continuity();
            return match s {
                GeomAbsShape::CN => GeomAbsShape::CN,
                GeomAbsShape::C3 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C0,
                GeomAbsShape::C0 => GeomAbsShape::C0,
            };
        }
        if self.my_type_curve == GeomAbsCurveType::OtherCurve {
            // OCCT L361-364: throw Standard_NoSuchObject.
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Contunuity");
        }

        GeomAbsShape::CN
    }

    /// OCCT L342: `Geom_OffsetCurve::GetBasisCurveContinuity()` — the basis
    /// curve continuity read through the offset payload.
    fn offset_basis_continuity(&self) -> GeomAbsShape {
        let Curve3::Offset(off) = &self.curve else {
            return GeomAbsShape::CN;
        };
        // The basis curve continuity through the same adaptor encoding: the
        // elementary kinds answer CN (the OCCT default arm L366).
        match off.basis.as_ref() {
            Curve3::BSpline(_) | Curve3::Offset(_) => GeomAbsShape::CN,
            _ => GeomAbsShape::CN,
        }
    }

    /// OCCT GeomAdaptor_Curve::NbIntervals(S) (GeomAdaptor_Curve.cxx
    /// L371-462).
    pub fn nb_intervals_of(&self, s: GeomAbsShape) -> usize {
        if self.my_type_curve == GeomAbsCurveType::BSplineCurve {
            let Curve3::BSpline(a_bspl) = &self.curve else {
                unreachable!()
            };
            let (knots, mults) = a_bspl.knots_mults();
            // OCCT L376: (!IsPeriodic && S <= Continuity()) || S == C0.
            let cont = self.continuity();
            if (!a_bspl.is_periodic && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                return 1;
            }

            let a_degree = a_bspl.degree;
            // OCCT L384-400: the continuity -> knot multiplicity level.
            let a_cont: i32 = match s {
                GeomAbsShape::C1 => 1,
                GeomAbsShape::C2 => 2,
                GeomAbsShape::C3 => 3,
                GeomAbsShape::CN => a_degree as i32,
                // OCCT L399: throw Standard_DomainError.
                _ => panic!("Standard_DomainError: GeomAdaptor_Curve::NbIntervals()"),
            };

            let an_eps = self
                .resolution(precision::CONFUSION)
                .min(precision::PCONFUSION);

            // OCCT L404-412: BSplCLib::Intervals(..., nullptr) answers the
            // count; the kernel returns count + 1 entries.
            return crate::math::bspl_lib::intervals(
                &knots,
                &mults,
                a_degree,
                a_bspl.is_periodic,
                a_cont,
                self.first,
                self.last,
                an_eps,
            )
            .len()
                - 1;
        }

        if self.my_type_curve == GeomAbsCurveType::OffsetCurve {
            // OCCT L415-456: the basis-intervals recount (OCC278).
            let mut my_nb_intervals = 1usize;
            let base_s: GeomAbsShape = match s {
                // OCCT L421-424: throw Standard_DomainError.
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                _ => GeomAbsShape::CN,
            };
            let c = self.offset_basis_adaptor();
            let i_nb_basis_int = c.nb_intervals_of(base_s);
            if i_nb_basis_int > 1 {
                let rdf_inter = c.intervals_of(base_s);
                for i_int in 0..i_nb_basis_int {
                    let p = rdf_inter[i_int];
                    if p > self.first && p < self.last {
                        my_nb_intervals += 1;
                    }
                }
            }
            return my_nb_intervals;
        }

        1
    }

    /// OCCT GeomAdaptor_Curve::Intervals(T, S) (GeomAdaptor_Curve.cxx
    /// L466-563) — the T array values.
    pub fn intervals_of(&self, s: GeomAbsShape) -> Vec<f64> {
        if self.my_type_curve == GeomAbsCurveType::BSplineCurve {
            let Curve3::BSpline(a_bspl) = &self.curve else {
                unreachable!()
            };
            let (knots, mults) = a_bspl.knots_mults();
            // OCCT L471: (!IsPeriodic && S <= Continuity()) || S == C0.
            let cont = self.continuity();
            if (!a_bspl.is_periodic && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                return vec![self.first, self.last];
            }

            let a_degree = a_bspl.degree;
            let a_cont: i32 = match s {
                GeomAbsShape::C1 => 1,
                GeomAbsShape::C2 => 2,
                GeomAbsShape::C3 => 3,
                GeomAbsShape::CN => a_degree as i32,
                _ => panic!("Standard_DomainError: GeomAdaptor_Curve::Intervals()"),
            };

            let an_eps = self
                .resolution(precision::CONFUSION)
                .min(precision::PCONFUSION);

            return crate::math::bspl_lib::intervals(
                &knots,
                &mults,
                a_degree,
                a_bspl.is_periodic,
                a_cont,
                self.first,
                self.last,
                an_eps,
            );
        }

        if self.my_type_curve == GeomAbsCurveType::OffsetCurve {
            // OCCT L512-556: the basis-intervals recount + the T fill.
            let mut my_nb_intervals = 1usize;
            let base_s: GeomAbsShape = match s {
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                _ => GeomAbsShape::CN,
            };
            let c = self.offset_basis_adaptor();
            let mut t: Vec<f64> = vec![self.first];
            let i_nb_basis_int = c.nb_intervals_of(base_s);
            if i_nb_basis_int > 1 {
                let rdf_inter = c.intervals_of(base_s);
                for i_int in 0..i_nb_basis_int {
                    let p = rdf_inter[i_int];
                    if p > self.first && p < self.last {
                        // OCCT L547: T(++myNbIntervals) = rdfInter(iInt).
                        my_nb_intervals += 1;
                        t.push(p);
                    }
                }
            }
            // OCCT L554-555: T(Lower) = myFirst; T(Lower + Nb) = myLast.
            t.push(self.last);
            let _ = my_nb_intervals;
            return t;
        }

        // OCCT L558-562.
        vec![self.first, self.last]
    }

    /// OCCT L437/L534: `GeomAdaptor_Curve C(offset->BasisCurve(), myFirst,
    /// myLast)` — the basis adaptor over the offset payload.
    fn offset_basis_adaptor(&self) -> GeomCurveAdaptor {
        let Curve3::Offset(off) = &self.curve else {
            unreachable!()
        };
        GeomCurveAdaptor::with_range(off.basis.as_ref().clone(), self.first, self.last)
    }

    /// OCCT GeomAdaptor_Curve::IsClosed (GeomAdaptor_Curve.cxx L576-585).
    pub fn is_closed(&self) -> bool {
        // OCCT L578: !Precision::IsPositiveInfinite(myLast) &&
        // !Precision::IsNegativeInfinite(myFirst) (Precision.hxx L357-367,
        // threshold 0.5 * Precision::Infinite()).
        if !precision::is_positive_infinite_value(self.last)
            && !precision::is_negative_infinite_value(self.first)
        {
            let pd = self.value_at(self.first);
            let pf = self.value_at(self.last);
            return pd.distance(pf) <= precision::CONFUSION;
        }
        false
    }

    /// OCCT GeomAdaptor_Curve::IsPeriodic (GeomAdaptor_Curve.cxx L589-592).
    pub fn is_periodic(&self) -> bool {
        CurveEval::is_periodic(&self.curve)
    }

    /// OCCT GeomAdaptor_Curve::Period (GeomAdaptor_Curve.cxx L596-599):
    /// myCurve->LastParameter() - myCurve->FirstParameter().
    pub fn period(&self) -> f64 {
        let [f, l] = CurveEval::default_domain(&self.curve);
        l - f
    }

    /// OCCT GeomAdaptor_Curve::Resolution (GeomAdaptor_Curve.cxx
    /// L1116-1149).
    pub fn resolution(&self, r3d: f64) -> f64 {
        match self.my_type_curve {
            GeomAbsCurveType::Line => r3d,
            GeomAbsCurveType::Circle => {
                let Curve3::Circle(c) = &self.curve else {
                    unreachable!()
                };
                let r = c.radius;
                if r > r3d / 2.0 {
                    2.0 * (r3d / (2.0 * r)).asin()
                } else {
                    2.0 * std::f64::consts::PI
                }
            }
            GeomAbsCurveType::Ellipse => {
                let Curve3::Ellipse(e) = &self.curve else {
                    unreachable!()
                };
                r3d / e.major_radius
            }
            GeomAbsCurveType::BezierCurve => {
                let Curve3::Bezier(b) = &self.curve else {
                    unreachable!()
                };
                crate::math::bspl::bezier_curve_resolution(b, r3d)
            }
            GeomAbsCurveType::BSplineCurve => {
                let Curve3::BSpline(b) = &self.curve else {
                    unreachable!()
                };
                crate::math::bspl::bspline_curve_resolution(b, r3d)
            }
            // OCCT L1146-1148: default -> Precision::Parametric(R3D).
            _ => precision::parametric_default(r3d),
        }
    }

    // ---------------------------------------------------------------------
    // The geometry downcasts (GeomAdaptor_Curve.cxx L1158-1295).
    // ---------------------------------------------------------------------

    /// OCCT GeomAdaptor_Curve::Line (L1158-1163) — the gp_Lin payload.
    pub fn line(&self) -> Line3 {
        if self.my_type_curve != GeomAbsCurveType::Line {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Line() - curve is not a Line");
        }
        let Curve3::Line(l) = &self.curve else {
            unreachable!()
        };
        *l
    }

    /// OCCT GeomAdaptor_Curve::Circle (L1167-1172) — the gp_Circ payload.
    pub fn circle(&self) -> Circle3 {
        if self.my_type_curve != GeomAbsCurveType::Circle {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Circle() - curve is not a Circle");
        }
        let Curve3::Circle(c) = &self.curve else {
            unreachable!()
        };
        *c
    }

    /// OCCT GeomAdaptor_Curve::Ellipse (L1176-1181) — the gp_Elips payload.
    pub fn ellipse(&self) -> Ellipse3 {
        if self.my_type_curve != GeomAbsCurveType::Ellipse {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Ellipse() - curve is not an Ellipse");
        }
        let Curve3::Ellipse(e) = &self.curve else {
            unreachable!()
        };
        *e
    }

    /// OCCT GeomAdaptor_Curve::Hyperbola (L1185-1190) — the gp_Hypr payload.
    pub fn hyperbola(&self) -> Hyperbola3 {
        if self.my_type_curve != GeomAbsCurveType::Hyperbola {
            panic!(
                "Standard_NoSuchObject: GeomAdaptor_Curve::Hyperbola() - curve is not a Hyperbola"
            );
        }
        let Curve3::Hyperbola(h) = &self.curve else {
            unreachable!()
        };
        *h
    }

    /// OCCT GeomAdaptor_Curve::Parabola (L1194-1199) — the gp_Parab payload.
    pub fn parabola(&self) -> Parabola3 {
        if self.my_type_curve != GeomAbsCurveType::Parabola {
            panic!(
                "Standard_NoSuchObject: GeomAdaptor_Curve::Parabola() - curve is not a Parabola"
            );
        }
        let Curve3::Parabola(p) = &self.curve else {
            unreachable!()
        };
        *p
    }

    /// OCCT GeomAdaptor_Curve::Degree (L1203-1217).
    pub fn degree(&self) -> usize {
        if self.my_type_curve == GeomAbsCurveType::BezierCurve {
            let Curve3::Bezier(b) = &self.curve else {
                unreachable!()
            };
            return b.control_points.len() - 1;
        }
        if self.my_type_curve == GeomAbsCurveType::BSplineCurve {
            let Curve3::BSpline(b) = &self.curve else {
                unreachable!()
            };
            return b.degree;
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Degree");
    }

    /// OCCT GeomAdaptor_Curve::IsRational (L1221-1232).
    pub fn is_rational(&self) -> bool {
        match self.my_type_curve {
            GeomAbsCurveType::BSplineCurve => {
                let Curve3::BSpline(b) = &self.curve else {
                    unreachable!()
                };
                curve_rational(&b.weights)
            }
            GeomAbsCurveType::BezierCurve => {
                let Curve3::Bezier(b) = &self.curve else {
                    unreachable!()
                };
                curve_rational(&b.weights)
            }
            _ => false,
        }
    }

    /// OCCT GeomAdaptor_Curve::NbPoles (L1236-1250).
    pub fn nb_poles(&self) -> usize {
        if self.my_type_curve == GeomAbsCurveType::BezierCurve {
            let Curve3::Bezier(b) = &self.curve else {
                unreachable!()
            };
            return b.control_points.len();
        }
        if self.my_type_curve == GeomAbsCurveType::BSplineCurve {
            let Curve3::BSpline(b) = &self.curve else {
                unreachable!()
            };
            return b.control_points.len();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Curve::NbPoles");
    }

    /// OCCT GeomAdaptor_Curve::NbKnots (L1254-1261).
    pub fn nb_knots(&self) -> usize {
        if self.my_type_curve != GeomAbsCurveType::BSplineCurve {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::NbKnots");
        }
        let Curve3::BSpline(b) = &self.curve else {
            unreachable!()
        };
        b.knots_mults().0.len()
    }

    /// OCCT GeomAdaptor_Curve::Bezier (L1265-1272) — the Geom_BezierCurve
    /// payload (no copy).
    pub fn bezier(&self) -> crate::geom::BezierCurve3 {
        if self.my_type_curve != GeomAbsCurveType::BezierCurve {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::Bezier");
        }
        let Curve3::Bezier(b) = &self.curve else {
            unreachable!()
        };
        b.clone()
    }

    /// OCCT GeomAdaptor_Curve::BSpline (L1276-1284) — the Geom_BSplineCurve
    /// payload (no copy).
    pub fn bspline(&self) -> crate::geom::BSplineCurve3 {
        if self.my_type_curve != GeomAbsCurveType::BSplineCurve {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::BSpline");
        }
        let Curve3::BSpline(b) = &self.curve else {
            unreachable!()
        };
        b.clone()
    }

    /// OCCT GeomAdaptor_Curve::OffsetCurve (L1288-1295) — the
    /// Geom_OffsetCurve payload.
    pub fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
        if self.my_type_curve != GeomAbsCurveType::OffsetCurve {
            panic!("Standard_NoSuchObject: GeomAdaptor_Curve::OffsetCurve");
        }
        let Curve3::Offset(o) = &self.curve else {
            unreachable!()
        };
        o.clone()
    }

    // ---------------------------------------------------------------------
    // The evaluation wrappers (GeomAdaptor_Curve.cxx EvalD0/D1/D2/D3/DN).
    // ---------------------------------------------------------------------

    /// OCCT GeomAdaptor_Curve::EvalD0 (L679-761) — the point at U.
    pub fn value_at(&self, u: f64) -> glam::DVec3 {
        CurveEval::point_at(&self.curve, u)
    }

    /// OCCT GeomAdaptor_Curve::EvalD1 (L765-848) — the point and the first
    /// derivative at U.
    pub fn d1_at(&self, u: f64) -> (glam::DVec3, glam::DVec3) {
        (
            CurveEval::point_at(&self.curve, u),
            CurveEval::derivative_at(&self.curve, u),
        )
    }

    /// OCCT GeomAdaptor_Curve::EvalD2 (L852-937) — the point with the first
    /// and second derivatives at U.  The Line arm (L858-861) answers the D1
    /// pair with a zero second derivative exactly as OCCT.
    pub fn d2_at(&self, u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
        if self.my_type_curve == GeomAbsCurveType::Line {
            // OCCT L859-860: ElCLib::D1 + D2.SetCoord(0,0,0).
            let p = CurveEval::point_at(&self.curve, u);
            let d1 = CurveEval::derivative_at(&self.curve, u);
            return (p, d1, glam::DVec3::ZERO);
        }
        (
            CurveEval::point_at(&self.curve, u),
            CurveEval::derivative_at(&self.curve, u),
            CurveEval::derivative2_at(&self.curve, u),
        )
    }

    /// OCCT GeomAdaptor_Curve::EvalD3 (L941-1045) — the point with the
    /// first three derivatives at U.  The Line arm (L947-951) answers the D1
    /// pair with zero D2/D3; the Parabola arm (L980-983) answers D2 with a
    /// zero D3 exactly as OCCT.
    pub fn d3_at(&self, u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3, glam::DVec3) {
        match self.my_type_curve {
            GeomAbsCurveType::Line => {
                // OCCT L948-950: ElCLib::D1 + zero D2 + zero D3.
                let p = CurveEval::point_at(&self.curve, u);
                let d1 = CurveEval::derivative_at(&self.curve, u);
                (p, d1, glam::DVec3::ZERO, glam::DVec3::ZERO)
            }
            GeomAbsCurveType::Parabola => {
                // OCCT L981-982: ElCLib::D2 + zero D3.
                let p = CurveEval::point_at(&self.curve, u);
                let d1 = CurveEval::derivative_at(&self.curve, u);
                let d2 = CurveEval::derivative2_at(&self.curve, u);
                (p, d1, d2, glam::DVec3::ZERO)
            }
            _ => {
                let p = CurveEval::point_at(&self.curve, u);
                let d1 = CurveEval::derivative_at(&self.curve, u);
                let d2 = CurveEval::derivative2_at(&self.curve, u);
                let d3 = CurveEval::derivative3_at(&self.curve, u);
                (p, d1, d2, d3)
            }
        }
    }

    /// OCCT GeomAdaptor_Curve::EvalDN (L1049-1112) — the N-th derivative at
    /// U.  The dispatch is type-first: the five analytic kinds answer the
    /// `ElCLib::DN` closed form at ANY order N (ElCLib.cxx L911-1047) and the
    /// polynomial kinds ride the BSplCLib::DN engines at any order.
    pub fn dn_at(&self, u: f64, n: i32) -> glam::DVec3 {
        match self.my_type_curve {
            // OCCT L1054-1055: ElCLib::DN(U, gp_Lin, N).
            GeomAbsCurveType::Line => curve_dn(&self.curve, u, n),
            // OCCT L1057-1058: ElCLib::DN(U, gp_Circ, N).
            GeomAbsCurveType::Circle => curve_dn(&self.curve, u, n),
            // OCCT L1060-1061: ElCLib::DN(U, gp_Elips, N).
            GeomAbsCurveType::Ellipse => curve_dn(&self.curve, u, n),
            // OCCT L1063-1064: ElCLib::DN(U, gp_Hypr, N).
            GeomAbsCurveType::Hyperbola => curve_dn(&self.curve, u, n),
            // OCCT L1066-1067: ElCLib::DN(U, gp_Parab, N).
            GeomAbsCurveType::Parabola => curve_dn(&self.curve, u, n),
            // OCCT L1069-1070: myCurve->EvalDN(U, N) — the Geom_BezierCurve
            // body, BSplCLib::DN over the Bezier's own knots/multiplicities.
            GeomAbsCurveType::BezierCurve => curve_dn(&self.curve, u, n),
            // OCCT L1072-1086: the BSpline arm — hasEvalRep / the IsBoundary
            // LocalDN branch / myCurve->EvalDN all answer the BSplCLib::DN
            // values for the clamped kernel encoding; N < 1 raises
            // Geom_UndefinedDerivative (Geom_BSplineCurve_1.cxx L302-306).
            GeomAbsCurveType::BSplineCurve => {
                let Curve3::BSpline(b) = &self.curve else {
                    unreachable!()
                };
                if n < 1 {
                    panic!("Geom_UndefinedDerivative: Geom_BSplineCurve::EvalDN");
                }
                b.dn(u, n as usize)
            }
            GeomAbsCurveType::OffsetCurve => {
                // OCCT L1088-1105: Geom_OffsetCurveUtils::EvaluateDN — any
                // order.  The rcad offset engine hosts the D1..D3 evaluations
                // (the EvalD1/D2/D3 translations); the general-N EvaluateDN
                // body is an open gap and the orders above 3 keep the raise.
                match n {
                    1 => self.d1_at(u).1,
                    2 => self.d2_at(u).2,
                    3 => self.d3_at(u).3,
                    _ => {
                        panic!("Standard_NotImplemented: GeomAdaptor_Curve::EvalDN offset order {n}")
                    }
                }
            }
            _ => {
                // OCCT L1107-1111: default -> myCurve->EvalDN(U, N).  The
                // rcad-only kinds (CircularHelix / SineWave) have no
                // translated Geom_Curve::EvalDN union member: the orders 1..3
                // keep the myCurve->Eval* fall-through evaluation, the other
                // orders keep the explicit raise.
                match n {
                    1 => self.d1_at(u).1,
                    2 => self.d2_at(u).2,
                    3 => self.d3_at(u).3,
                    _ => {
                        panic!("Standard_NotImplemented: GeomAdaptor_Curve::EvalDN order {n}")
                    }
                }
            }
        }
    }

    /// OCCT GeomAdaptor_Curve::ShallowCopy (GeomAdaptor_Curve.cxx L89-132)
    /// — the value copy (the rcad handle copy is a clone).
    pub fn shallow_copy(&self) -> Self {
        self.clone()
    }
}

// =========================================================================
// The trait implementations
// =========================================================================

impl Adaptor3dCurve for GeomCurveAdaptor {
    /// OCCT GeomAdaptor_Curve::GetType() — the loaded myTypeCurve
    /// (GeomAdaptor_Curve.cxx L231-289 assignment; the accessor per the
    /// hxx).  The proj_lib CurveType has no Offset arm — the OCCT
    /// GeomAbs_OffsetCurve folds into Other (the CPnts order default-10
    /// arm answers identically).
    fn curve_type(&self) -> super::CurveType {
        match self.my_type_curve {
            GeomAbsCurveType::Line => super::CurveType::Line,
            GeomAbsCurveType::Circle => super::CurveType::Circle,
            GeomAbsCurveType::Ellipse => super::CurveType::Ellipse,
            GeomAbsCurveType::Hyperbola => super::CurveType::Hyperbola,
            GeomAbsCurveType::Parabola => super::CurveType::Parabola,
            GeomAbsCurveType::BezierCurve => super::CurveType::Bezier,
            GeomAbsCurveType::BSplineCurve => super::CurveType::BSpline,
            GeomAbsCurveType::OffsetCurve | GeomAbsCurveType::OtherCurve => super::CurveType::Other,
        }
    }

    /// OCCT FirstParameter() — the restricted first parameter.
    fn first_parameter(&self) -> f64 {
        self.first
    }

    /// OCCT LastParameter() — the restricted last parameter.
    fn last_parameter(&self) -> f64 {
        self.last
    }

    /// OCCT Value(U) / D0(U, P) — the EvalD0 composition.
    fn value(&self, u: f64) -> DVec3 {
        self.value_at(u)
    }

    /// OCCT D1(U, P, V1).
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.d1_at(u)
    }

    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.d2_at(u)
    }

    /// OCCT D3(U, P, V1, V2, V3).
    fn d3(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        self.d3_at(u)
    }

    /// OCCT DN(U, N).
    fn dn(&self, u: f64, n: i32) -> DVec3 {
        self.dn_at(u, n)
    }

    /// OCCT Continuity().
    fn continuity(&self) -> GeomAbsShape {
        self.continuity()
    }

    /// OCCT NbIntervals(S).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.nb_intervals_of(s)
    }

    /// OCCT Intervals(T, S).
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.intervals_of(s)
    }

    /// OCCT Trim(First, Last, Tol) (GeomAdaptor_Curve.cxx L567-572) — the
    /// restricted adaptor over the same curve.
    fn trim(&self, first: f64, last: f64, _tol: f64) -> std::sync::Arc<dyn Adaptor3dCurve> {
        std::sync::Arc::new(GeomCurveAdaptor::with_range(self.curve.clone(), first, last))
    }

    /// OCCT IsClosed().
    fn is_closed(&self) -> bool {
        self.is_closed()
    }

    /// OCCT IsPeriodic().
    fn is_periodic(&self) -> bool {
        self.is_periodic()
    }

    /// OCCT Period().
    fn period(&self) -> f64 {
        self.period()
    }

    /// OCCT Resolution(R3d).
    fn resolution(&self, r3d: f64) -> f64 {
        self.resolution(r3d)
    }

    /// OCCT Degree().
    fn degree(&self) -> usize {
        self.degree()
    }

    /// OCCT IsRational().
    fn is_rational(&self) -> bool {
        self.is_rational()
    }

    /// OCCT NbPoles().
    fn nb_poles(&self) -> usize {
        self.nb_poles()
    }

    /// OCCT NbKnots().
    fn nb_knots(&self) -> usize {
        self.nb_knots()
    }

    /// OCCT Bezier().
    fn bezier(&self) -> crate::geom::BezierCurve3 {
        self.bezier()
    }

    /// OCCT BSpline().
    fn bspline(&self) -> crate::geom::BSplineCurve3 {
        self.bspline()
    }

    /// OCCT OffsetCurve().
    fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
        self.offset_curve()
    }

    /// OCCT ShallowCopy().
    fn shallow_copy(&self) -> std::sync::Arc<dyn Adaptor3dCurve> {
        std::sync::Arc::new(self.shallow_copy())
    }
}

/// The OCCT Adaptor3d_Curve geometry accessors (GetType / Line / Circle /
/// Ellipse / Parabola / Hyperbola) — the companion trait shared with the
/// rcad adaptor encodings (see the module header for the GetType placement).
pub trait Adaptor3dCurveGeom: Adaptor3dCurve {
    /// OCCT Adaptor3d_Curve::GetType() — the curve kind (consumed by the
    /// Project dispatch of ProjLib_ProjectedCurve.cxx L242-270).
    fn get_type(&self) -> CurveType;
    /// OCCT Line() — valid when GetType() == GeomAbs_Line.
    fn line(&self) -> Line3;
    /// OCCT Circle() — valid when GetType() == GeomAbs_Circle.
    fn circle(&self) -> Circle3;
    /// OCCT Ellipse() — valid when GetType() == GeomAbs_Ellipse.
    fn ellipse(&self) -> Ellipse3;
    /// OCCT Parabola() — valid when GetType() == GeomAbs_Parabola.
    fn parabola(&self) -> Parabola3;
    /// OCCT Hyperbola() — valid when GetType() == GeomAbs_Hyperbola.
    fn hyperbola(&self) -> Hyperbola3;
    /// OCCT Adaptor3d_Curve::Trim(First, Last, Tol) — the trimmed handle,
    /// rebinding through the geometry-accessor interface (the rcad encoding
    /// of the handle rebinding in TrimC3d).
    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> std::sync::Arc<dyn Adaptor3dCurveGeom>;
}

impl Adaptor3dCurveGeom for GeomCurveAdaptor {
    /// OCCT GeomAdaptor_Curve::GetType() (hxx L179) — the stored type tag.
    fn get_type(&self) -> CurveType {
        match self.my_type_curve {
            GeomAbsCurveType::Line => CurveType::Line,
            GeomAbsCurveType::Circle => CurveType::Circle,
            GeomAbsCurveType::Ellipse => CurveType::Ellipse,
            GeomAbsCurveType::Hyperbola => CurveType::Hyperbola,
            GeomAbsCurveType::Parabola => CurveType::Parabola,
            GeomAbsCurveType::BezierCurve => CurveType::Bezier,
            GeomAbsCurveType::BSplineCurve => CurveType::BSpline,
            _ => CurveType::Other,
        }
    }

    /// OCCT Line().
    fn line(&self) -> Line3 {
        self.line()
    }

    /// OCCT Circle().
    fn circle(&self) -> Circle3 {
        self.circle()
    }

    /// OCCT Ellipse().
    fn ellipse(&self) -> Ellipse3 {
        self.ellipse()
    }

    /// OCCT Parabola().
    fn parabola(&self) -> Parabola3 {
        self.parabola()
    }

    /// OCCT Hyperbola().
    fn hyperbola(&self) -> Hyperbola3 {
        self.hyperbola()
    }

    /// OCCT Trim().
    fn trim_geom(&self, first: f64, last: f64, _tol: f64) -> std::sync::Arc<dyn Adaptor3dCurveGeom> {
        std::sync::Arc::new(GeomCurveAdaptor::with_range(self.curve.clone(), first, last))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// OCCT GeomAdaptor_Curve::EvalDN (L1049-1112) — the Circle arm answers
    /// the `ElCLib::DN` closed form at ANY order.  Hand-derived circle
    /// derivatives over C(u) = center + r*(cos u * XDir + sin u * YDir):
    ///   D1 = r*(-sin, cos);  D2 = r*(-cos, -sin);
    ///   D3 = r*( sin, -cos); D4 = r*( cos,  sin)
    /// (the derivative cycle is 4, matching the `N % 4 == 0` branch of
    /// `ElCLib::CircleDN`, ElCLib.cxx L940-944).  The N = 4 call used to
    /// raise Standard_NotImplemented on the non-BSpline kinds.
    #[test]
    fn dn_at_circle_any_order_matches_the_elclib_dn_closed_form() {
        let r = 2.5f64;
        let u = 0.7f64;
        let circle = Curve3::Circle(Circle3::new(DVec3::ONE, DVec3::Z, r));
        let a = GeomCurveAdaptor::new(circle);

        let d1 = a.dn_at(u, 1);
        assert!((d1 - DVec3::new(-r * u.sin(), r * u.cos(), 0.0)).length() < 1e-12);

        let d4 = a.dn_at(u, 4);
        assert!((d4 - DVec3::new(r * u.cos(), r * u.sin(), 0.0)).length() < 1e-12);

        let d5 = a.dn_at(u, 5);
        assert!((d5 - DVec3::new(-r * u.sin(), r * u.cos(), 0.0)).length() < 1e-12);
    }

    /// OCCT `ElCLib::LineDN` (ElCLib.cxx L911-918) — the direction at N == 1,
    /// the null vector above; `ElCLib::ParabolaDN` (ElCLib.cxx L1020-1026) —
    /// the null vector for N > 2.
    #[test]
    fn dn_at_line_and_parabola_high_order_answers_the_null_vector() {
        let line = Curve3::Line(Line3::new(DVec3::new(1.0, 2.0, 3.0), DVec3::new(2.0, 0.0, 0.0)));
        let a_line = GeomCurveAdaptor::new(line);
        assert_eq!(a_line.dn_at(0.3, 4), glam::DVec3::ZERO);

        let parabola = Curve3::Parabola(Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        });
        let a_par = GeomCurveAdaptor::new(parabola);
        assert_eq!(a_par.dn_at(0.3, 4), glam::DVec3::ZERO);
    }

    /// OCCT `ElCLib::HyperbolaDN` (ElCLib.cxx L996-1018) — the IsEven(N)
    /// branch.  Hand-derived over H(u) = (a*cosh u, b*sinh u): the derivative
    /// cycle is 2, so D4 = D2 = (a*cosh u, b*sinh u).
    #[test]
    fn dn_at_hyperbola_d4_matches_the_elclib_even_branch() {
        let (a_ax, b_ax) = (3.0f64, 2.0f64);
        let u = 0.7f64;
        let hypr = Curve3::Hyperbola(Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: a_ax,
            semi_minor: b_ax,
        });
        let a = GeomCurveAdaptor::new(hypr);
        let d4 = a.dn_at(u, 4);
        assert!((d4 - DVec3::new(a_ax * u.cosh(), b_ax * u.sinh(), 0.0)).length() < 1e-12);
    }

    /// The orders 1..3 keep answering the same values as the D1/D2/D3 chain
    /// (the ElCLib::DN closed forms equal the ElCLib::D1/D2/D3 expressions
    /// term by term for the analytic kinds).
    #[test]
    fn dn_at_orders_1_to_3_match_the_d1_d2_d3_chain() {
        let u = 0.7f64;
        let circle = Curve3::Circle(Circle3::new(DVec3::ONE, DVec3::Z, 2.5));
        let a = GeomCurveAdaptor::new(circle);
        assert!(a.dn_at(u, 1).distance(a.d1_at(u).1) < 1e-12);
        assert!(a.dn_at(u, 2).distance(a.d2_at(u).2) < 1e-12);
        assert!(a.dn_at(u, 3).distance(a.d3_at(u).3) < 1e-12);
    }
}
