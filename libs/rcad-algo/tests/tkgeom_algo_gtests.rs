//! TKGeomAlgo GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKGeomAlgo/GTests/
//!
//! Files translated:
//!   IntSurf_Quadric_Test.cxx — cone-apex gradient/distance finiteness
//!     (IntSurf_Quadric::Gradient / ValAndGrad).
//!   IntSurf_LineOn2S_Test.cxx — line-on-two-surfaces box caches, point
//!     replacement invalidation and Split semantics (IntSurf_LineOn2S +
//!     IntSurf_PntOn2S).
//!   TopTrans_SurfaceTransition_Test.cxx — instance-local undefined state and
//!     the Compare tolerance (TopTrans_SurfaceTransition).
//!   IntPolyh_Point_Test.cxx — the IntPolyh_Point arithmetic (Add/Sub/Divide/
//!     SquareModulus/SquareDistance/Dot).
//!   IntCurveSurface_IntersectionPoint_Test.cxx — the default/value
//!     constructors and SetValues of IntCurveSurface_IntersectionPoint.
//!   IntCurveSurface_InterUtils_Test.cxx — SectionPointToParameters with the
//!     degenerate-face fallback to the longest edge.
//!
//! 1:1 translations of the OCCT tests against the rcad OCCT-aligned APIs.

use rcad_algo::geomalgo::int_polyh::IntPolyhPoint;
use rcad_algo::geomalgo::int_patch::int_cs::{IntersectionPoint, TransitionOnCurve};
use rcad_algo::geomalgo::int_surf::{LineOn2S, PntOn2S, Quadric};
use rcad_algo::geomalgo::top_trans::surface_transition::SurfaceTransition;
use rcad_algo::geomalgo::intf::{
    section_point_to_parameters, IntfPIType, IntfSectionPoint, PolyhedronLike, PolygonLike,
};
use rcad_kernel::core::precision::{ANGULAR, COMPUTATIONAL, CONFUSION, SQUARE_CONFUSION};
use rcad_kernel::geom::ConicalSurface;
use rcad_kernel::topods::{Orientation, State};

// =============================================================================
// IntSurf_Quadric_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_surf_quadric_tests {
    use super::*;

    // TEST(IntSurf_Quadric_Test, ConeApexGradientRemainsFinite)
    #[test]
    fn cone_apex_gradient_remains_finite() {
        // gp_Cone(gp_Ax3(gp_Pnt(0,0,0), gp_Dir(0,0,1), gp_Dir(1,0,0)), 0.5, 0.0)
        let cone = ConicalSurface::new_with_ref_dir(
            glam::DVec3::ZERO,
            glam::DVec3::Z,
            0.0,
            0.5,
            glam::DVec3::X,
        );
        let quadric = Quadric::from_cone(&cone);
        let an_apex = cone.apex_point();

        // EXPECT_NO_THROW: the gradient at the apex must stay finite.
        let a_gradient = quadric.gradient(an_apex);
        assert!(
            a_gradient.length_squared() <= SQUARE_CONFUSION,
            "apex gradient {:?} must be null (SquareConfusion)",
            a_gradient
        );

        // ValAndGrad(apex, dist, grad): dist ~ 0, grad null.
        let (a_dist, a_grad) = quadric.val_and_grad(an_apex);
        assert!(
            (a_dist - 0.0).abs() <= CONFUSION,
            "distance at the apex should be 0, got {a_dist}"
        );
        assert!(
            a_grad.length_squared() <= SQUARE_CONFUSION,
            "ValAndGrad gradient {:?} must be null",
            a_grad
        );
    }
}

// =============================================================================
// IntSurf_LineOn2S_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_surf_line_on_2s_tests {
    use super::*;

    /// OCCT buildPoint: IntSurf_PntOn2S::SetValue(P, U1, V1, U2, V2).
    fn build_point(pt: glam::DVec3, u1: f64, v1: f64, u2: f64, v2: f64) -> PntOn2S {
        let mut p = PntOn2S::new();
        p.set_value_all(pt, u1, v1, u2, v2);
        p
    }

    // TEST(IntSurf_LineOn2S_Test, EmptyLineBoxesAreNotOut)
    #[test]
    fn empty_line_boxes_are_not_out() {
        let mut line = LineOn2S::new();
        assert!(!line.is_out_box(glam::DVec3::new(10.0, 10.0, 10.0)));
        assert!(!line.is_out_surf1_box(glam::DVec2::new(10.0, 10.0)));
        assert!(!line.is_out_surf2_box(glam::DVec2::new(10.0, 10.0)));
    }

    // TEST(IntSurf_LineOn2S_Test, PointReplacementInvalidatesCachedBoxes)
    #[test]
    fn point_replacement_invalidates_cached_boxes() {
        let mut line = LineOn2S::new();
        line.add(&build_point(
            glam::DVec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
            0.0,
            0.0,
        ));
        line.add(&build_point(
            glam::DVec3::new(1.0, 1.0, 1.0),
            1.0,
            1.0,
            1.0,
            1.0,
        ));

        assert!(line.is_out_box(glam::DVec3::new(100.0, 100.0, 100.0)));
        assert!(line.is_out_surf1_box(glam::DVec2::new(50.0, 50.0)));

        // SetPoint(2, ...) — 1-based index 2 -> rcad 0-based 1.
        line.set_point(1, glam::DVec3::new(100.0, 100.0, 100.0));
        assert!(!line.is_out_box(glam::DVec3::new(100.0, 100.0, 100.0)));

        // Value(2, P) — replace the point with full UV parameters.
        line.set_value(
            1,
            &build_point(
                glam::DVec3::new(100.0, 100.0, 100.0),
                50.0,
                50.0,
                60.0,
                60.0,
            ),
        );
        assert!(!line.is_out_surf1_box(glam::DVec2::new(50.0, 50.0)));
        assert!(!line.is_out_surf2_box(glam::DVec2::new(60.0, 60.0)));
    }

    // TEST(IntSurf_LineOn2S_Test, Split_DividesCorrectly)
    // OCCT NCollection_Sequence::Split(Index, SS) moves items Index..N into
    // SS, keeping 1..Index-1 in the original (1-based); rcad split uses a
    // 0-based index.
    #[test]
    fn split_divides_correctly() {
        let mut line = LineOn2S::new();
        line.add(&build_point(
            glam::DVec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
            0.0,
            0.0,
        ));
        line.add(&build_point(
            glam::DVec3::new(1.0, 0.0, 0.0),
            1.0,
            0.0,
            1.0,
            0.0,
        ));
        line.add(&build_point(
            glam::DVec3::new(2.0, 0.0, 0.0),
            2.0,
            0.0,
            2.0,
            0.0,
        ));
        line.add(&build_point(
            glam::DVec3::new(3.0, 0.0, 0.0),
            3.0,
            0.0,
            3.0,
            0.0,
        ));

        // Split(2) (1-based) -> keep 1 point, move 3 into the split line.
        let split = line.split(1);

        assert_eq!(line.nb_points(), 1);
        assert_eq!(split.nb_points(), 3);
        assert_eq!(line.value(0).value().x, 0.0);
        assert_eq!(split.value(0).value().x, 1.0);
    }
}

// =============================================================================
// TopTrans_SurfaceTransition_Test.cxx
// =============================================================================

#[cfg(test)]
mod top_trans_surface_transition_tests {
    use super::*;

    // TEST(TopTrans_SurfaceTransition_Test, UndefinedStateStaysInstanceLocal)
    #[test]
    fn undefined_state_stays_instance_local() {
        // aReferenceTransition.Reset(Dir(1,0,0), Dir(0,0,1));
        let mut reference = SurfaceTransition::new();
        reference.reset(glam::DVec3::X, glam::DVec3::Z);
        reference.compare(
            ANGULAR,
            glam::DVec3::Z,
            Orientation::Forward,
            Orientation::Forward,
        );

        let a_before = reference.state_before();
        let an_after = reference.state_after();
        assert_ne!(a_before, State::Unknown);
        assert_ne!(an_after, State::Unknown);

        // anInvalidTransition.Reset(Dir(1,0,0), Dir(0,0,1), Dir(1,0,0),
        //                          Dir(1,0,0), 1.0, 1.0);
        let mut invalid = SurfaceTransition::new();
        invalid.reset_full(
            glam::DVec3::X,
            glam::DVec3::Z,
            glam::DVec3::X,
            glam::DVec3::X,
            1.0,
            1.0,
        );

        assert_eq!(invalid.state_before(), State::Unknown);
        assert_eq!(invalid.state_after(), State::Unknown);
        assert_eq!(reference.state_before(), a_before);
        assert_eq!(reference.state_after(), an_after);
    }

    // TEST(TopTrans_SurfaceTransition_Test, CompareUsesProvidedTolerance)
    #[test]
    fn compare_uses_provided_tolerance() {
        let mut transition = SurfaceTransition::new();
        transition.reset(glam::DVec3::X, glam::DVec3::Z);
        transition.compare_full(
            1.0,
            glam::DVec3::Z,
            glam::DVec3::new(0.6, 0.8, 0.0),
            glam::DVec3::new(0.6, -0.8, 0.0),
            1.0,
            0.5,
            Orientation::Forward,
            Orientation::Forward,
        );

        assert_ne!(transition.state_before(), State::Unknown);
        assert_ne!(transition.state_after(), State::Unknown);
    }
}

// =============================================================================
// IntCurveSurface_IntersectionPoint_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_curve_surface_intersection_point_tests {
    use super::*;

    // TEST(IntCurveSurface_IntersectionPoint, DefaultConstructor_TransitionInitialized)
    #[test]
    fn default_constructor_transition_initialized() {
        let pt = IntersectionPoint::default();
        // Transition must be deterministic (initialized to Tangent).
        assert_eq!(pt.transition(), TransitionOnCurve::Tangent);
        // Scalar fields must be zero-initialized.
        assert_eq!(pt.u(), 0.0);
        assert_eq!(pt.v(), 0.0);
        assert_eq!(pt.w(), 0.0);
    }

    // TEST(IntCurveSurface_IntersectionPoint, ValueConstructor_AllFieldsSet)
    #[test]
    fn value_constructor_all_fields_set() {
        let a_p = glam::DVec3::new(1.0, 2.0, 3.0);
        let pt = IntersectionPoint::new(a_p, 0.5, 0.6, 0.7, TransitionOnCurve::In);

        assert!((pt.pnt().x - 1.0).abs() < 1e-15);
        assert!((pt.pnt().y - 2.0).abs() < 1e-15);
        assert!((pt.pnt().z - 3.0).abs() < 1e-15);
        assert_eq!(pt.u(), 0.5);
        assert_eq!(pt.v(), 0.6);
        assert_eq!(pt.w(), 0.7);
        assert_eq!(pt.transition(), TransitionOnCurve::In);
    }

    // TEST(IntCurveSurface_IntersectionPoint, SetValues_OverwritesTransition)
    #[test]
    fn set_values_overwrites_transition() {
        let mut pt = IntersectionPoint::default();
        assert_eq!(pt.transition(), TransitionOnCurve::Tangent);

        let a_p = glam::DVec3::new(5.0, 6.0, 7.0);
        pt.set_values(a_p, 1.0, 2.0, 3.0, TransitionOnCurve::Out);
        assert_eq!(pt.transition(), TransitionOnCurve::Out);
    }
}

// =============================================================================
// IntCurveSurface_InterUtils_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_curve_surface_inter_utils_tests {
    use super::*;

    /// OCCT TestPolyhedron — 3 points, 1 triangle (1-based indices).
    struct TestPolyhedron {
        points: [glam::DVec3; 3],
        u: [f64; 3],
        v: [f64; 3],
    }

    impl TestPolyhedron {
        fn new() -> Self {
            TestPolyhedron {
                points: [
                    glam::DVec3::new(0.0, 0.0, 0.0),
                    glam::DVec3::new(2.0, 0.0, 0.0),
                    glam::DVec3::new(2.0, 0.0, 0.0),
                ],
                u: [0.0, 1.0, 2.0],
                v: [0.0, 0.0, 0.0],
            }
        }
    }

    impl PolyhedronLike for TestPolyhedron {
        // Triangle(int, P1, P2, P3) { P1=1; P2=2; P3=3; } (1-based).
        fn triangle(&self, _t: usize) -> (usize, usize, usize) {
            (1, 2, 3)
        }
        // Point(theIndex) { return myPoints[theIndex - 1]; }
        fn point(&self, index: usize) -> glam::DVec3 {
            self.points[index - 1]
        }
        // Parameters(theIndex, U, V) { ... myU[theIndex - 1] ... }
        fn parameters(&self, index: usize) -> (f64, f64) {
            (self.u[index - 1], self.v[index - 1])
        }
    }

    /// OCCT TestPolygon.
    struct TestPolygon;

    impl PolygonLike for TestPolygon {
        // ApproxParamOnCurve(Index, ParamOnLine) = 10 * Index + ParamOnLine.
        fn approx_param_on_curve(&self, index: usize, param_on_line: f64) -> f64 {
            10.0 * index as f64 + param_on_line
        }
    }

    // TEST(IntCurveSurface_InterUtils, SectionPointToParameters_DegenerateFaceFallsBackToLongestEdge)
    #[test]
    fn section_point_to_parameters_degenerate_face_falls_back_to_longest_edge() {
        let polyhedron = TestPolyhedron::new();
        let polygon = TestPolygon;

        // Intf_SectionPoint(gp_Pnt(0.5,0,0), Intf_EDGE, 0, 2, 0.4, Intf_FACE, 1, 0, 0.0, 0.0)
        let section_point = IntfSectionPoint::new(
            glam::DVec3::new(0.5, 0.0, 0.0),
            IntfPIType::Edge,
            0,
            2,
            0.4,
            IntfPIType::Face,
            1,
            0,
            0.0,
            0.0,
        );

        let (a_u, a_v, a_w) = section_point_to_parameters(&section_point, &polyhedron, &polygon);

        assert!((a_u - 0.25).abs() < 1.0e-12);
        assert!((a_v - 0.0).abs() < 1.0e-12);
        assert!((a_w - 20.4).abs() < 1.0e-12);
    }
}

// =============================================================================
// IntPolyh_Point_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_polyh_point_tests {
    use super::*;

    // TEST(IntPolyh_Point, DefaultConstructor_AllZero)
    #[test]
    fn default_constructor_all_zero() {
        let pt = IntPolyhPoint::new();
        assert_eq!(pt.x(), 0.0);
        assert_eq!(pt.y(), 0.0);
        assert_eq!(pt.z(), 0.0);
        assert_eq!(pt.u(), 0.0);
        assert_eq!(pt.v(), 0.0);
    }

    // TEST(IntPolyh_Point, Divide_NormalDivisor_CorrectResult)
    #[test]
    fn divide_normal_divisor_correct_result() {
        let pt = IntPolyhPoint::new_uv(10.0, 20.0, 30.0, 0.5, 0.8);
        let res = pt.divide(2.0);
        assert_eq!(res.x(), 5.0);
        assert_eq!(res.y(), 10.0);
        assert_eq!(res.z(), 15.0);
        assert_eq!(res.u(), 0.25);
        assert_eq!(res.v(), 0.4);
    }

    // TEST(IntPolyh_Point, Divide_ZeroDivisor_ReturnsDefaultPoint)
    #[test]
    fn divide_zero_divisor_returns_default_point() {
        let pt = IntPolyhPoint::new_uv(10.0, 20.0, 30.0, 0.5, 0.8);
        let res = pt.divide(0.0);
        assert_eq!(res.x(), 0.0);
        assert_eq!(res.y(), 0.0);
        assert_eq!(res.z(), 0.0);
        assert_eq!(res.u(), 0.0);
        assert_eq!(res.v(), 0.0);
    }

    // TEST(IntPolyh_Point, Divide_NearZeroDivisor_ReturnsDefaultPoint)
    #[test]
    fn divide_near_zero_divisor_returns_default_point() {
        let pt = IntPolyhPoint::new_uv(1.0, 2.0, 3.0, 0.1, 0.2);
        // |1e-20| <= Precision::Computational() (machine epsilon).
        assert!(1e-20 <= COMPUTATIONAL);
        let res = pt.divide(1.0e-20);
        assert_eq!(res.x(), 0.0);
        assert_eq!(res.y(), 0.0);
        assert_eq!(res.z(), 0.0);
    }

    // TEST(IntPolyh_Point, Divide_NegativeDivisor_CorrectResult)
    #[test]
    fn divide_negative_divisor_correct_result() {
        let pt = IntPolyhPoint::new_uv(6.0, -9.0, 12.0, 0.3, 0.6);
        let res = pt.divide(-3.0);
        assert_eq!(res.x(), -2.0);
        assert_eq!(res.y(), 3.0);
        assert_eq!(res.z(), -4.0);
    }

    // TEST(IntPolyh_Point, Add_CorrectResult)
    #[test]
    fn add_correct_result() {
        let p1 = IntPolyhPoint::new_uv(1.0, 2.0, 3.0, 0.1, 0.2);
        let p2 = IntPolyhPoint::new_uv(4.0, 5.0, 6.0, 0.3, 0.4);
        let res = p1.add(&p2);
        assert_eq!(res.x(), 5.0);
        assert_eq!(res.y(), 7.0);
        assert_eq!(res.z(), 9.0);
    }

    // TEST(IntPolyh_Point, Sub_CorrectResult)
    #[test]
    fn sub_correct_result() {
        let p1 = IntPolyhPoint::new_uv(5.0, 7.0, 9.0, 0.5, 0.8);
        let p2 = IntPolyhPoint::new_uv(1.0, 2.0, 3.0, 0.1, 0.2);
        let res = p1.sub(&p2);
        assert_eq!(res.x(), 4.0);
        assert_eq!(res.y(), 5.0);
        assert_eq!(res.z(), 6.0);
    }

    // TEST(IntPolyh_Point, SquareModulus_CorrectResult)
    #[test]
    fn square_modulus_correct_result() {
        let pt = IntPolyhPoint::new_uv(1.0, 2.0, 3.0, 0.0, 0.0);
        assert_eq!(pt.square_modulus(), 14.0);
    }

    // TEST(IntPolyh_Point, SquareDistance_SamePoint_Zero)
    #[test]
    fn square_distance_same_point_zero() {
        let pt = IntPolyhPoint::new_uv(3.0, 4.0, 5.0, 0.0, 0.0);
        assert_eq!(pt.square_distance(&pt), 0.0);
    }

    // TEST(IntPolyh_Point, Dot_OrthogonalVectors_Zero)
    #[test]
    fn dot_orthogonal_vectors_zero() {
        let p1 = IntPolyhPoint::new_uv(1.0, 0.0, 0.0, 0.0, 0.0);
        let p2 = IntPolyhPoint::new_uv(0.0, 1.0, 0.0, 0.0, 0.0);
        assert_eq!(p1.dot(&p2), 0.0);
    }
}
