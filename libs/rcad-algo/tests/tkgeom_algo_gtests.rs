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
//!   IntCurveSurface_ThePolygonOfHInter_Test.cxx — the Closed flag of the
//!     curve-sampling polygon.
//!   IntCurveSurface_ThePolyhedronOfHInter_Test.cxx — singularity flags, the
//!     size/triangles/points of the surface-sampling polyhedron, the
//!     parameter-array constructor (min-size acceptance, single-value
//!     rejection) and PlaneEquation finiteness.
//!   IntPatch_Polyhedron_Test.cxx — the auto-computed subdivision, the
//!     zero-subdivision clamp and the TriConnex no-crash guarantees.
//!   Intf_Tool_Test.cxx — Hypr2dBox/Parab2dBox/ParabBox/HyprBox valid
//!     segments plus the no-intersection zero-segments case.
//!
//! 1:1 translations of the OCCT tests against the rcad OCCT-aligned APIs.

use glam::{DVec2, DVec3};
use rcad_algo::geomalgo::int_polyh::IntPolyhPoint;
use rcad_algo::geomalgo::int_patch::int_cs::{IntersectionPoint, TransitionOnCurve};
use rcad_algo::geomalgo::int_surf::{LineOn2S, PntOn2S, Quadric};
use rcad_algo::geomalgo::top_trans::surface_transition::SurfaceTransition;
use rcad_algo::geomalgo::intf::{
    section_point_to_parameters, IntfPIType, IntfSectionPoint, IntfTool, PolyhedronLike,
    PolygonLike,
};
use rcad_algo::geomalgo::{
    IntPatchPolyhedron, ThePolygonOfHInter, ThePolyhedronOfHInter,
};
use rcad_kernel::core::precision::{
    ANGULAR, COMPUTATIONAL, CONFUSION, INFINITE_VALUE, SQUARE_CONFUSION,
};
use rcad_kernel::geom::{
    ConicalSurface, Curve2d, Curve2dEval, Hyperbola2d, Hyperbola3, Line2d, Line3, Parabola2d,
    Parabola3, Plane, SphericalSurface, Surface3,
};
use rcad_kernel::math::bnd::{BndBox, BndBox2d};
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

// =============================================================================
// IntCurveSurface_ThePolygonOfHInter_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_curve_surface_the_polygon_of_hinter_tests {
    use super::*;

    // TEST(IntCurveSurface_ThePolygonOfHInter_Test, ClosedFlagReflectsStoredState)
    #[test]
    fn closed_flag_reflects_stored_state() {
        // Geom_Line(gp_Ax1(gp_Pnt(0,0,0), gp_Dir(1,0,0))).
        let line = Line3::new(DVec3::ZERO, DVec3::X);
        let mut polygon = ThePolygonOfHInter::new(&line, 6);

        assert!(!polygon.closed());
        polygon.set_closed(true);
        assert!(polygon.closed());
    }
}

// =============================================================================
// IntCurveSurface_ThePolyhedronOfHInter_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_curve_surface_the_polyhedron_of_hinter_tests {
    use super::*;

    /// OCCT makeTestPlane — gp_Ax3(Pnt(0,0,0), Dir(0,0,1)) adapted with the
    /// UV domain [-1, 1]².
    fn make_test_plane() -> Plane {
        Plane::new(DVec3::ZERO, DVec3::Z)
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, SingularityFlags_InitializedFalse)
    #[test]
    fn singularity_flags_initialized_false() {
        let surf = make_test_plane();
        let poly = ThePolyhedronOfHInter::new(&surf, 3, 3, -1.0, -1.0, 1.0, 1.0);

        assert!(!poly.has_u_min_singularity());
        assert!(!poly.has_u_max_singularity());
        assert!(!poly.has_v_min_singularity());
        assert!(!poly.has_v_max_singularity());
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, SingularitySetters_UpdateFlags)
    #[test]
    fn singularity_setters_update_flags() {
        let surf = make_test_plane();
        let mut poly = ThePolyhedronOfHInter::new(&surf, 3, 3, -1.0, -1.0, 1.0, 1.0);

        poly.set_u_min_singularity(true);
        assert!(poly.has_u_min_singularity());
        assert!(!poly.has_u_max_singularity());

        poly.set_v_max_singularity(true);
        assert!(poly.has_v_max_singularity());
        assert!(!poly.has_v_min_singularity());
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, BasicConstruction_ValidMesh)
    #[test]
    fn basic_construction_valid_mesh() {
        let surf = make_test_plane();
        let poly = ThePolyhedronOfHInter::new(&surf, 4, 4, -1.0, -1.0, 1.0, 1.0);

        let (nb_u, nb_v) = poly.size();
        assert_eq!(nb_u, 4);
        assert_eq!(nb_v, 4);
        assert_eq!(poly.nb_triangles(), 4 * 4 * 2);
        assert_eq!(poly.nb_points(), (4 + 1) * (4 + 1));
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, ParamArrayConstructor_MinimumSize)
    // Length-2 arrays give nbdelta = 1 — no division by zero.
    #[test]
    fn param_array_constructor_minimum_size() {
        let surf = make_test_plane();

        let u_pars = [-1.0, 1.0];
        let v_pars = [-1.0, 1.0];

        let poly = ThePolyhedronOfHInter::new_params(&surf, &u_pars, &v_pars);

        let (nb_u, nb_v) = poly.size();
        assert!(nb_u >= 1);
        assert!(nb_v >= 1);
        assert!(poly.nb_triangles() > 0);
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, ParamArrayConstructor_RejectsSingleValueArrays)
    // Constructor validation: a single-value array raises Standard_OutOfRange.
    #[test]
    fn param_array_constructor_rejects_single_value_arrays() {
        let surf = make_test_plane();

        // A single-value U array raises.
        assert!(std::panic::catch_unwind(|| {
            ThePolyhedronOfHInter::new_params(&surf, &[0.0], &[-1.0, 1.0]);
        })
        .is_err());

        // A single-value V array raises.
        assert!(std::panic::catch_unwind(|| {
            ThePolyhedronOfHInter::new_params(&surf, &[-1.0, 1.0], &[0.0]);
        })
        .is_err());
    }

    // TEST(IntCurveSurface_ThePolyhedronOfHInter, PlaneEquation_FiniteResults)
    #[test]
    fn plane_equation_finite_results() {
        let surf = make_test_plane();
        let poly = ThePolyhedronOfHInter::new(&surf, 3, 3, -1.0, -1.0, 1.0, 1.0);

        let (a_normal, a_polar_dist) = poly.plane_equation(1);
        assert!(a_normal.x.is_finite());
        assert!(a_normal.y.is_finite());
        assert!(a_normal.z.is_finite());
        assert!(a_polar_dist.is_finite());
    }
}

// =============================================================================
// IntPatch_Polyhedron_Test.cxx
// =============================================================================

#[cfg(test)]
mod int_patch_polyhedron_tests {
    use super::*;

    /// OCCT makePlane — gp_Ax3(Pnt(0,0,0), Dir(0,0,1)) adapted with [-1,1]².
    fn make_plane() -> Surface3 {
        Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))
    }

    /// OCCT makeSphere — gp_SphericalSurface(gp_Ax3(Pnt(0,0,0), Dir(0,0,1)), 1.0).
    fn make_sphere() -> Surface3 {
        Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
        })
    }

    // TEST(IntPatch_Polyhedron, DefaultConstructor_ProducesValidMesh)
    #[test]
    fn default_constructor_produces_valid_mesh() {
        let a_surf = make_sphere();
        let a_poly = IntPatchPolyhedron::new(&a_surf);

        let (nb_u, nb_v) = a_poly.size();
        assert!(nb_u > 0);
        assert!(nb_v > 0);
        assert!(a_poly.nb_triangles() > 0);
        assert!(a_poly.nb_points() > 0);
    }

    // TEST(IntPatch_Polyhedron, ZeroSubdivision_ClampedToMinimum)
    // Zero nbu/nbv must not cause division-by-zero — clamped to minimum 1.
    #[test]
    fn zero_subdivision_clamped_to_minimum() {
        let a_surf = make_plane();

        let a_poly = IntPatchPolyhedron::new_sub(&a_surf, 0, 0);
        let (nb_u, nb_v) = a_poly.size();
        assert!(nb_u >= 1);
        assert!(nb_v >= 1);
        assert!(a_poly.nb_triangles() > 0);
    }

    // TEST(IntPatch_Polyhedron, SmallSubdivision_ProducesValidMesh)
    #[test]
    fn small_subdivision_produces_valid_mesh() {
        let a_surf = make_plane();
        let a_poly = IntPatchPolyhedron::new_sub(&a_surf, 2, 2);

        let (nb_u, nb_v) = a_poly.size();
        assert_eq!(nb_u, 2);
        assert_eq!(nb_v, 2);
        assert_eq!(a_poly.nb_triangles(), 2 * 2 * 2);
    }

    // TEST(IntPatch_Polyhedron, TriConnex_PedgeZero_NoCrash)
    // TriConnex must not crash when Pedge == 0.
    #[test]
    fn tri_connex_pedge_zero_no_crash() {
        let a_surf = make_sphere();
        let a_poly = IntPatchPolyhedron::new_sub(&a_surf, 4, 4);

        let (a_p1, _a_p2, _a_p3) = a_poly.triangle(1);

        // Call TriConnex with Pedge = 0 (unknown edge mode).
        let (a_result, _a_tri_con, _an_other_p) = (a_poly.tri_connex(1, a_p1, 0).0, 0, 0);
        // Must not crash; result is a valid triangle index or 0.
        assert!(a_result >= 0);
    }

    // TEST(IntPatch_Polyhedron, TriConnex_AllVertices_NoCrash)
    #[test]
    fn tri_connex_all_vertices_no_crash() {
        let a_surf = make_sphere();
        let a_poly = IntPatchPolyhedron::new_sub(&a_surf, 3, 3);

        let (a_p1, a_p2, a_p3) = a_poly.triangle(1);

        // Call TriConnex with each vertex as Pedge and with Pedge=0.
        let _ = a_poly.tri_connex(1, a_p1, 0);
        let _ = a_poly.tri_connex(1, a_p1, a_p2);
        let _ = a_poly.tri_connex(1, a_p1, a_p3);
        let _ = a_poly.tri_connex(1, a_p2, a_p3);
        // All calls must complete without crash.
    }
}

// =============================================================================
// Intf_Tool_Test.cxx
// =============================================================================

#[cfg(test)]
mod intf_tool_tests {
    use super::*;

    // TEST(Intf_Tool, Hypr2dBox_ProducesValidSegments)
    #[test]
    fn hypr2d_box_produces_valid_segments() {
        // gp_Hypr2d(gp_Ax2d(gp_Pnt2d(0,0), gp_Dir2d(1,0)), 2.0, 1.0).
        let a_hypr = Hyperbola2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            semi_major: 2.0,
            semi_minor: 1.0,
        };
        let mut a_domain = BndBox2d::new();
        a_domain.update(-5.0, -5.0, 5.0, 5.0);

        let mut a_result = BndBox2d::new();
        let mut a_tool = IntfTool::new();
        a_tool.hypr_2d_box(&a_hypr, &a_domain, &mut a_result);

        let a_nb_seg = a_tool.nb_segments();
        assert!(a_nb_seg >= 0);
        assert!(a_nb_seg <= 6);

        for i in 1..=a_nb_seg {
            let a_begin = a_tool.begin_param(i);
            let a_end = a_tool.end_param(i);
            assert!(
                a_begin.is_finite() || a_begin == -INFINITE_VALUE,
                "Segment {i} begin={a_begin}"
            );
            assert!(
                a_end.is_finite() || a_end == INFINITE_VALUE,
                "Segment {i} end={a_end}"
            );
            assert!(a_begin <= a_end, "Segment {i} begin > end");
        }
    }

    // TEST(Intf_Tool, Parab2dBox_ProducesValidSegments)
    #[test]
    fn parab2d_box_produces_valid_segments() {
        // gp_Parab2d(gp_Ax2d(gp_Pnt2d(0,0), gp_Dir2d(1,0)), 1.0) — focal
        // length 1.0, so the rcad focal parameter p = 2.0.
        let a_parab = Parabola2d {
            origin: DVec2::ZERO,
            axis_dir: DVec2::X,
            focal_param: 2.0,
        };
        let mut a_domain = BndBox2d::new();
        a_domain.update(-5.0, -5.0, 5.0, 5.0);

        let mut a_result = BndBox2d::new();
        let mut a_tool = IntfTool::new();
        a_tool.parab_2d_box(&a_parab, &a_domain, &mut a_result);

        let a_nb_seg = a_tool.nb_segments();
        assert!(a_nb_seg >= 0);
        assert!(a_nb_seg <= 6);

        for i in 1..=a_nb_seg {
            let a_begin = a_tool.begin_param(i);
            let a_end = a_tool.end_param(i);
            assert!(
                a_begin.is_finite() || a_begin == -INFINITE_VALUE,
                "Segment {i} begin={a_begin}"
            );
            assert!(
                a_end.is_finite() || a_end == INFINITE_VALUE,
                "Segment {i} end={a_end}"
            );
            assert!(a_begin <= a_end, "Segment {i} begin > end");
        }
    }

    // TEST(Intf_Tool, ParabBox_ProducesValidSegments)
    #[test]
    fn parab_box_produces_valid_segments() {
        // gp_Parab(gp_Ax2(gp_Pnt(0,0,0), gp_Dir(0,0,1)), 1.0) — focal length
        // 1.0, rcad focal parameter p = 2.0, axis = X in the XY plane.
        let a_parab = Parabola3 {
            vertex: DVec3::ZERO,
            normal: DVec3::Z,
            axis_dir: DVec3::X,
            focal_param: 2.0,
        };
        let mut a_domain = BndBox::new();
        a_domain.update(-5.0, -5.0, -5.0, 5.0, 5.0, 5.0);

        let mut a_result = BndBox::new();
        let mut a_tool = IntfTool::new();
        a_tool.parab_box(&a_parab, &a_domain, &mut a_result);

        let a_nb_seg = a_tool.nb_segments();
        assert!(a_nb_seg >= 0);
        assert!(a_nb_seg <= 6);

        for i in 1..=a_nb_seg {
            let a_begin = a_tool.begin_param(i);
            let a_end = a_tool.end_param(i);
            assert!(a_begin <= a_end, "Segment {i} begin > end");
        }
    }

    // TEST(Intf_Tool, HyprBox_ProducesValidSegments)
    #[test]
    fn hypr_box_produces_valid_segments() {
        // gp_Hypr(gp_Ax2(gp_Pnt(0,0,0), gp_Dir(0,0,1)), 2.0, 1.0).
        let a_hypr = Hyperbola3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            major_dir: DVec3::X,
            semi_major: 2.0,
            semi_minor: 1.0,
        };
        let mut a_domain = BndBox::new();
        a_domain.update(-5.0, -5.0, -5.0, 5.0, 5.0, 5.0);

        let mut a_result = BndBox::new();
        let mut a_tool = IntfTool::new();
        a_tool.hypr_box(&a_hypr, &a_domain, &mut a_result);

        let a_nb_seg = a_tool.nb_segments();
        assert!(a_nb_seg >= 0);
        assert!(a_nb_seg <= 6);
    }

    // TEST(Intf_Tool, Hypr2dBox_NoIntersection_ZeroSegments)
    #[test]
    fn hypr2d_box_no_intersection_zero_segments() {
        // A small hyperbola far away from the small domain.
        let a_hypr = Hyperbola2d {
            center: DVec2::new(100.0, 100.0),
            major_dir: DVec2::X,
            semi_major: 0.1,
            semi_minor: 0.05,
        };
        let mut a_domain = BndBox2d::new();
        a_domain.update(-1.0, -1.0, 1.0, 1.0);

        let mut a_result = BndBox2d::new();
        let mut a_tool = IntfTool::new();
        a_tool.hypr_2d_box(&a_hypr, &a_domain, &mut a_result);

        assert_eq!(a_tool.nb_segments(), 0);
    }
}

// =============================================================================
// Geom2dAPI_Interpolate_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_api_interpolate_tests {
    use super::*;
    use rcad_kernel::base::geom_api::geom2d_interpolate::Geom2dInterpolate;
    use rcad_kernel::geom::BSplineCurve2;

    /// Distinct values of the expanded knot vector (OCCT BSplineCurve::Knot(i)).
    fn distinct_knots(curve: &BSplineCurve2) -> Vec<f64> {
        let mut out: Vec<f64> = Vec::new();
        for &k in &curve.knots {
            match out.last() {
                Some(&last) if (k - last).abs() < 1e-15 => {}
                _ => out.push(k),
            }
        }
        out
    }

    // TEST(Geom2dAPI_InterpolateTest, OCC28594_InterpolateWithAndWithoutTangentScale)
    #[test]
    fn occ28594_interpolate_with_and_without_tangent_scale() {
        let a_points = vec![
            DVec2::new(-30.4, 8.0),
            DVec2::new(-16.689912, 17.498217),
            DVec2::new(-23.803064, 24.748543),
            DVec2::new(-16.907466, 32.919615),
            DVec2::new(-8.543829, 26.549421),
            DVec2::new(0.0, 39.200000),
        ];
        let a_tangents = vec![
            DVec2::new(0.3, 0.4),
            DVec2::ZERO,
            DVec2::ZERO,
            DVec2::ZERO,
            DVec2::ZERO,
            DVec2::new(1.0, 0.0),
        ];
        let a_tangent_flags = vec![true, false, false, false, false, true];

        // Interpolation with tangent scale.
        let mut an_interp_with_scale = Geom2dInterpolate::new(a_points.clone(), false, CONFUSION);
        an_interp_with_scale.load(&a_tangents, &a_tangent_flags, true);
        an_interp_with_scale.perform();
        assert!(an_interp_with_scale.is_done());
        let a_curve_with_scale = an_interp_with_scale.curve();

        // Interpolation without tangent scale.
        let mut an_interp_without_scale = Geom2dInterpolate::new(a_points.clone(), false, CONFUSION);
        an_interp_without_scale.load(&a_tangents, &a_tangent_flags, false);
        an_interp_without_scale.perform();
        assert!(an_interp_without_scale.is_done());
        let a_curve_without_scale = an_interp_without_scale.curve();

        // Both curves must pass through all given points (the distinct knots
        // are the chord-length parameters).
        let a_tol = CONFUSION * 10.0;
        let a_knots = distinct_knots(&a_curve_with_scale);
        assert_eq!(a_knots.len(), a_points.len());
        for an_index in 0..a_points.len() {
            let a_pt_on_curve = a_curve_with_scale.point_at(a_knots[an_index]);
            let a_pt = a_points[an_index];
            assert!((a_pt.x - a_pt_on_curve.x).abs() < a_tol, " at point index {an_index}");
            assert!((a_pt.y - a_pt_on_curve.y).abs() < a_tol, " at point index {an_index}");
            // The without-scale curve passes through the same points.
            let a_pt_on_curve2 = a_curve_without_scale.point_at(a_knots[an_index]);
            assert!((a_pt.x - a_pt_on_curve2.x).abs() < a_tol, " at point index {an_index}");
            assert!((a_pt.y - a_pt_on_curve2.y).abs() < a_tol, " at point index {an_index}");
        }
    }
}

// =============================================================================
// GeomAPI_ProjectPointOnSurf_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom_api_project_point_on_surf_tests {
    use super::*;
    use rcad_kernel::base::geom_api::project_on_surf::ProjectPointOnSurf;
    use rcad_kernel::geom::{CylindricalSurface, TrimmedSurface};

    // TEST(GeomAPI_ProjectPointOnSurfTest, Bug867_InitWithTightBoundsNoException)
    // Bug OCC867: calling Init() with explicit UV bounds then Perform() must
    // not throw.  The bounds are tighter (U:[0,3]) than the trimmed surface
    // (U:[0,4]) to reproduce the original bug scenario.
    #[test]
    fn bug867_init_with_tight_bounds_no_exception() {
        // Geom_CylindricalSurface(gp::XOY(), 10.0).
        let a_cyl = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 10.0,
            ref_dir: DVec3::X,
            y_dir: None,
        });
        // Geom_RectangularTrimmedSurface(aCyl, 0.0, 4.0, 0.0, 2.0).
        let a_trim_surf = Surface3::Trimmed(TrimmedSurface::new(a_cyl, 0.0, 4.0, 0.0, 2.0));

        let a_point = DVec3::new(30.0, 30.0, 30.0);

        let mut a_projector = ProjectPointOnSurf::new();
        // EXPECT_NO_THROW — the original bug (OCC867) was that Init()+Perform()
        // threw an exception.
        a_projector.init_surface(&a_trim_surf, 0.0, 3.0, 0.0, 2.0);
        a_projector.perform(a_point);
        // Whether a projection exists within the given UV bounds is not the
        // concern here.
    }
}

// =============================================================================
// Geom2dHatch_Elements_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_hatch_elements_tests {
    use super::*;
    use rcad_algo::geomalgo::hatch::{HatchElement, HatchElements};
    use rcad_kernel::geom::Circle2d;

    /// OCCT makeCircleElement — gp_Ax2d(Pnt2d(0,0), Dir2d(1,0)) circle r=1,
    /// element with TopAbs_FORWARD.
    fn make_circle_element() -> HatchElement {
        let a_circle = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 1.0,
        };
        HatchElement::new(Curve2d::Circle(a_circle), Orientation::Forward)
    }

    // TEST(Geom2dHatch_Elements, CurrentEdge_ReturnsValidData)
    // Bug #38: CurrentEdge in const context must work; edge traversal returns
    // valid curve and orientation after Bind + InitWires/InitEdges.
    #[test]
    fn current_edge_returns_valid_data() {
        let mut an_elements = HatchElements::new();
        an_elements.bind(1, make_circle_element());

        an_elements.init_wires();
        assert!(an_elements.more_wires());

        an_elements.init_edges();
        assert!(an_elements.more_edges());

        let (a_edge, an_or) = an_elements.current_edge();

        assert_eq!(an_or, Orientation::Forward);
        // Circle has parameter range [0, 2*PI].
        let dom = a_edge.default_domain();
        assert!((dom[0] - 0.0).abs() < 1e-10);
        assert!(dom[1] > 6.0);
    }

    // TEST(Geom2dHatch_Elements, BindAndFind)
    #[test]
    fn bind_and_find() {
        let mut an_elements = HatchElements::new();
        assert!(!an_elements.is_bound(1));

        an_elements.bind(1, make_circle_element());
        assert!(an_elements.is_bound(1));

        let an_elem = an_elements.find(1);
        assert_eq!(an_elem.orientation(), Orientation::Forward);
    }

    // TEST(Geom2dHatch_Elements, Clear_RemovesAll)
    #[test]
    fn clear_removes_all() {
        let mut an_elements = HatchElements::new();
        an_elements.bind(1, make_circle_element());
        an_elements.bind(2, make_circle_element());
        assert!(an_elements.is_bound(1));

        an_elements.clear();
        assert!(!an_elements.is_bound(1));
        assert!(!an_elements.is_bound(2));
    }
}

// =============================================================================
// Geom2dHatch_Intersector_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_hatch_intersector_tests {
    use super::*;
    use rcad_algo::geomalgo::hatch::HatchIntersector;
    use rcad_kernel::geom::{BSplineCurve2, Circle2d, Line2d};

    // TEST(Geom2dHatch_Intersector, LocalGeometry_Circle_ValidOutputs)
    #[test]
    fn local_geometry_circle_valid_outputs() {
        // gp_Ax2d(Pnt2d(0,0), Dir2d(1,0)) circle radius 1.
        let a_circle = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 1.0,
        };
        let a_curve_adaptor = Curve2d::Circle(a_circle);

        let an_intersector = HatchIntersector::with_tolerances(1.0e-7, 1.0e-7);

        let (a_tang, _a_norm, a_curv) = an_intersector.local_geometry(&a_curve_adaptor, 0.0);

        // Circle tangent at param=0 is (0,1), curvature is 1.0.
        assert!((a_tang.x - 0.0).abs() < 1.0e-10);
        assert!((a_tang.y - 1.0).abs() < 1.0e-10);
        assert!((a_curv - 1.0).abs() < 1.0e-10);
    }

    // TEST(Geom2dHatch_Intersector, LocalGeometry_Line_ZeroCurvature)
    #[test]
    fn local_geometry_line_zero_curvature() {
        // gp_Ax2d(Pnt2d(0,0), Dir2d(1,0)) line.
        let a_line = Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        };
        let a_curve_adaptor = Curve2d::Line(a_line);

        let an_intersector = HatchIntersector::with_tolerances(1.0e-7, 1.0e-7);

        let (a_tang, a_norm, a_curv) = an_intersector.local_geometry(&a_curve_adaptor, 0.5);

        assert!((a_tang.x - 1.0).abs() < 1.0e-10);
        assert!((a_tang.y - 0.0).abs() < 1.0e-10);
        assert!((a_curv - 0.0).abs() < 1.0e-10);
        // Normal should be perpendicular to tangent.
        assert!((a_norm.x - 0.0).abs() < 1.0e-10);
        assert!(a_norm.y.abs() - 1.0 < 1.0e-10);
    }

    // TEST(Geom2dHatch_Intersector, LocalGeometry_DegenerateCurve_InitializedOutputs)
    // Bug #24: LocalGeometry with undefined tangent must still produce
    // initialized (deterministic) output values.
    #[test]
    fn local_geometry_degenerate_curve_initialized_outputs() {
        // A degenerate BSpline where all control points coincide.
        let a_degenerate_curve = BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![DVec2::ONE; 4],
            weights: vec![1.0; 4],
        };
        let a_curve_adaptor = Curve2d::BSpline(a_degenerate_curve);

        let an_intersector = HatchIntersector::with_tolerances(1.0e-7, 1.0e-7);

        // Must not crash and must produce initialized outputs.
        let (a_tang, a_norm, a_curv) = an_intersector.local_geometry(&a_curve_adaptor, 0.5);

        // Curvature should be initialized to 0 for degenerate case.
        assert!((a_curv - 0.0).abs() < 1.0e-10);
        // Tang and Norm should be valid unit directions (magnitude == 1).
        let a_tang_mag = (a_tang.x * a_tang.x + a_tang.y * a_tang.y).sqrt();
        assert!((a_tang_mag - 1.0).abs() < 1.0e-10);
        let a_norm_mag = (a_norm.x * a_norm.x + a_norm.y * a_norm.y).sqrt();
        assert!((a_norm_mag - 1.0).abs() < 1.0e-10);
    }
}

// =============================================================================
// Geom2dGcc_Lin2d2Tan_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_gcc_lin2d2tan_tests {
    use super::*;
    use rcad_algo::geomalgo::geom2d_gcc::{GccEntPosition, Lin2d2Tan, QualifiedCurve};
    use rcad_kernel::geom::{Circle2d, Ellipse2d};

    /// OCCT GeomAPI::To2d — project a point onto the plane's 2D frame.
    fn to2d_point(origin: DVec3, x_dir: DVec3, y_dir: DVec3, p: DVec3) -> DVec2 {
        DVec2::new((p - origin).dot(x_dir), (p - origin).dot(y_dir))
    }

    // TEST(Geom2dGcc_Lin2d2TanTest, OCC813_EllipseAndPoint)
    // Tangent line from a 2D point to a projected 2D ellipse.
    #[test]
    fn occ813_ellipse_and_point() {
        // gp_Ax2(Pnt(1262.224429, 425.040878, 363.609716),
        //        Dir(0.173648, 0.984808, 0), Dir(-0.932169, 0.164367, -0.322560)).
        let an_ax2_loc = DVec3::new(1262.224429, 425.040878, 363.609716);
        let an_ax2_z = DVec3::new(0.173648, 0.984808, 0.0);
        let an_ax2_x = DVec3::new(-0.932169, 0.164367, -0.322560);
        let an_ax2_y = an_ax2_z.cross(an_ax2_x);

        // GeomAPI::To2d(ellipse, plane): the ellipse lies in the plane, so the
        // 2D ellipse is centered at the origin with radii 150/100.
        let a_curve2d = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 150.0,
            minor_radius: 100.0,
        });
        let a_q_curve = QualifiedCurve::new(a_curve2d, GccEntPosition::Outside);

        // Query tangent line from the 2D point to the projected ellipse.
        let a_pnt2d = DVec2::new(200.0, 200.0);
        let a_lin_tan = Lin2d2Tan::curve_point(&a_q_curve, a_pnt2d, 0.1);

        assert!(a_lin_tan.nb_solutions() > 0, "Expected at least one tangent line solution");
    }

    // TEST(Geom2dGcc_Lin2d2TanTest, OCC814_CircleAndEllipse)
    // Tangent line between a 2D circle and a 2D ellipse.
    #[test]
    fn occ814_circle_and_ellipse() {
        let an_ax2_loc = DVec3::new(1262.224429, 425.040878, 363.609716);
        let an_ax2_z = DVec3::new(0.173648, 0.984808, 0.0);
        let an_ax2_x = DVec3::new(-0.932169, 0.164367, -0.322560);
        let an_ax2_y = an_ax2_z.cross(an_ax2_x);

        // Projected 2D ellipse (origin, radii 150/100).
        let a_curve2d = Curve2d::Ellipse(Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            major_radius: 150.0,
            minor_radius: 100.0,
        });
        let a_q_ell = QualifiedCurve::new(a_curve2d, GccEntPosition::Outside);

        // gp_Circle(gp_Ax2(Pnt(823.687192, 502.366825, 478.960440), same axes), 50)
        // projected onto the plane.
        let a_circ_center3d = DVec3::new(823.687192, 502.366825, 478.960440);
        let a_from_curve2d = Curve2d::Circle(Circle2d {
            center: to2d_point(an_ax2_loc, an_ax2_x, an_ax2_y, a_circ_center3d),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 50.0,
        });
        let a_q_cir = QualifiedCurve::new(a_from_curve2d, GccEntPosition::Outside);

        let a_lin_tan = Lin2d2Tan::curve_curve(&a_q_ell, &a_q_cir, 0.1);

        assert!(a_lin_tan.nb_solutions() > 0, "Expected at least one tangent line solution");
    }
}

// Geom2dGcc_Circ2d2TanRad_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_gcc_circ2d2tanrad_tests {
    use super::*;
    use rcad_algo::geomalgo::geom2d_gcc::{Circ2d2TanRad, GccEntPosition, QualifiedCurve};
    use rcad_kernel::geom::BezierCurve2;

    // TEST(Geom2dGcc_Circ2d2TanRadTest, BUC60897_TangentToLineAndBezier)
    // Circles of radius 10 tangent to a line and a Bezier curve; each tangency
    // point must lie at the circle radius from the circle center (within 1%).
    #[test]
    fn buc60897_tangent_to_line_and_bezier() {
        // Geom2d_Line(gp_Pnt2d(100, 0), gp_Dir2d(NX)).
        let a_line = Curve2d::Line(Line2d::new(DVec2::new(100.0, 0.0), DVec2::new(-1.0, 0.0)));

        // Geom2d_BezierCurve with three control points (non-rational, so the
        // homogeneous weights are 1.0 per the rcad BezierCurve2 convention).
        let a_curve = Curve2d::Bezier(BezierCurve2 {
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(50.0, 50.0),
                DVec2::new(0.0, 100.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });

        // Qualified curves with outside tangency.
        let a_qualif_curve1 = QualifiedCurve::new(a_line, GccEntPosition::Outside);
        let a_qualif_curve2 = QualifiedCurve::new(a_curve, GccEntPosition::Outside);

        // Find circles of radius 10 tangent to both curves.
        let a_radius = 10.0;
        let a_tolerance = 1e-7;
        let a_gcc_circ2d = Circ2d2TanRad::new_curve_curve(
            &a_qualif_curve1,
            &a_qualif_curve2,
            a_radius,
            a_tolerance,
        );

        assert!(a_gcc_circ2d.is_done(), "Geom2dGcc_Circ2d2TanRad failed to compute");
        assert!(
            a_gcc_circ2d.nb_solutions() > 0,
            "No tangent circles found"
        );

        let a_max_delta_percent = 1.0;
        for i in 1..=a_gcc_circ2d.nb_solutions() {
            let a_circ2d = a_gcc_circ2d.this_solution(i);
            let a_center = a_circ2d.center;
            let a_r = a_circ2d.radius;

            let (_a_par_sol1, _a_par_arg1, a_pnt_sol1) = a_gcc_circ2d.tangency1(i);
            let (_a_par_sol2, _a_par_arg2, a_pnt_sol2) = a_gcc_circ2d.tangency2(i);

            let a_d1 = a_pnt_sol1.distance(a_center);
            let a_delta1 = (a_d1 - a_r).abs() / a_r * 100.0;
            assert!(
                a_delta1 <= a_max_delta_percent,
                "Solution {}: tangency1 distance error {}% exceeds 1%",
                i,
                a_delta1
            );

            let a_d2 = a_pnt_sol2.distance(a_center);
            let a_delta2 = (a_d2 - a_r).abs() / a_r * 100.0;
            assert!(
                a_delta2 <= a_max_delta_percent,
                "Solution {}: tangency2 distance error {}% exceeds 1%",
                i,
                a_delta2
            );
        }
    }
}

// Geom2dGcc_Circ2d3Tan_Test.cxx
// =============================================================================

#[cfg(test)]
mod geom2d_gcc_circ2d3tan_tests {
    use super::*;
    use rcad_algo::geomalgo::geom2d_gcc::{Circ2d3Tan, GccEntPosition, QualifiedCurve};
    use rcad_kernel::geom::Circle2d;

    /// OCCT Geom2dGcc::Unqualified(adaptor) — a circle qualified as unqualified.
    fn create_qualified_circle(x: f64, y: f64, radius: f64) -> QualifiedCurve {
        let a_circle = Curve2d::Circle(Circle2d {
            center: DVec2::new(x, y),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius,
        });
        QualifiedCurve::new(a_circle, GccEntPosition::Unqualified)
    }

    /// OCCT verifySolutionValidity (Geom2dGcc_Circ2d3Tan_Test.cxx L82-104).
    fn verify_solution_validity(
        sol: &Circle2d,
        index: usize,
        min_x: f64,
        max_x: f64,
        min_y: f64,
        max_y: f64,
        max_radius: f64,
    ) {
        assert!(
            sol.radius > 0.0,
            "Solution {} should have positive radius",
            index
        );
        assert!(
            sol.radius < max_radius,
            "Solution {} should have reasonable radius",
            index
        );
        let center = sol.center;
        assert!(center.x > min_x, "Solution {} X coordinate should be reasonable", index);
        assert!(center.x < max_x, "Solution {} X coordinate should be reasonable", index);
        assert!(center.y > min_y, "Solution {} Y coordinate should be reasonable", index);
        assert!(center.y < max_y, "Solution {} Y coordinate should be reasonable", index);
    }

    /// OCCT verifyTangencyConstraints (L36-67) — the solution must be tangent
    /// to all three input circles (external or internal tangency).
    fn verify_tangency_constraints(
        sol: &Circle2d,
        c1: &Circle2d,
        c2: &Circle2d,
        c3: &Circle2d,
        index: usize,
        tolerance: f64,
    ) {
        let sol_center = sol.center;
        let sol_radius = sol.radius;
        let circles = [c1, c2, c3];
        for (i, c) in circles.iter().enumerate() {
            let dist = sol_center.distance(c.center);
            let expected_ext = sol_radius + c.radius;
            let expected_int = (sol_radius - c.radius).abs();
            let is_tangent = (dist - expected_ext).abs() <= tolerance
                || (dist - expected_int).abs() <= tolerance;
            assert!(
                is_tangent,
                "Solution {} should be tangent to circle {} (distance={}, expected external={}, expected internal={})",
                index,
                i + 1,
                dist,
                expected_ext,
                expected_int
            );
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, BUC60622_RegressionCase)
    #[test]
    fn buc60622_regression_case() {
        let my_tolerance = CONFUSION; // Precision::Confusion()
        let a_qual1 = create_qualified_circle(500.0, 1800.0, 500.0);
        let a_qual2 = create_qualified_circle(500.0, 1900.0, 400.0);
        let a_qual3 = create_qualified_circle(700.0, 1900.0, 200.0);

        let a_input_circ1 = Circle2d {
            center: DVec2::new(500.0, 1800.0),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 500.0,
        };
        let a_input_circ2 = Circle2d {
            center: DVec2::new(500.0, 1900.0),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 400.0,
        };
        let a_input_circ3 = Circle2d {
            center: DVec2::new(700.0, 1900.0),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 200.0,
        };

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Algorithm should succeed");

        let a_nb_solutions = a_solver.nb_solutions();
        assert_eq!(a_nb_solutions, 3, "BUC60622 case should find exactly 3 solutions");

        for i in 1..=a_nb_solutions {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -1000.0, 2000.0, 1000.0, 3000.0, 10000.0);
            verify_tangency_constraints(&a_sol, &a_input_circ1, &a_input_circ2, &a_input_circ3, i, 1e-6);
        }

        if a_nb_solutions >= 3 {
            let a_sol1 = a_solver.this_solution(1);
            assert!((a_sol1.center.x - 500.0).abs() <= 1.0);
            assert!((a_sol1.center.y - 1900.0).abs() <= 1.0);
            assert!((a_sol1.radius - 400.0).abs() <= 1.0);
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, ToleranceImpact_Analysis)
    #[test]
    fn tolerance_impact_analysis() {
        let a_qual1 = create_qualified_circle(500.0, 1800.0, 500.0);
        let a_qual2 = create_qualified_circle(500.0, 1900.0, 400.0);
        let a_qual3 = create_qualified_circle(700.0, 1900.0, 200.0);

        let a_test_tolerances = [CONFUSION, 1e-12, 1e-10, 1e-8];
        let mut a_default_sol_count = 0;

        for (idx, a_tol) in a_test_tolerances.iter().enumerate() {
            let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, *a_tol, 0.0, 0.0, 0.0);
            assert!(a_solver.is_done(), "Algorithm should succeed with tolerance {}", a_tol);

            let a_nb_sol = a_solver.nb_solutions();
            if idx == 0 {
                a_default_sol_count = a_nb_sol;
                assert!(a_nb_sol >= 1, "Should find at least 1 solution with default tolerance");
            } else {
                assert!(a_nb_sol >= 1, "Should find at least 1 solution with tolerance {}", a_tol);
                assert!(
                    a_nb_sol <= a_default_sol_count + 2,
                    "Solution count shouldn't increase dramatically with tolerance {}",
                    a_tol
                );
            }
            for i in 1..=a_nb_sol {
                let a_sol = a_solver.this_solution(i);
                verify_solution_validity(&a_sol, i, -1000.0, 2000.0, 1000.0, 3000.0, 10000.0);
            }
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, Simple_ThreeCircle_Case)
    #[test]
    fn simple_three_circle_case() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 2.0);
        let a_qual2 = create_qualified_circle(10.0, 0.0, 2.0);
        let a_qual3 = create_qualified_circle(5.0, 8.0, 2.0);
        let a_input1 = Circle2d::new(DVec2::new(0.0, 0.0), 2.0);
        let a_input2 = Circle2d::new(DVec2::new(10.0, 0.0), 2.0);
        let a_input3 = Circle2d::new(DVec2::new(5.0, 8.0), 2.0);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Simple case should always work");

        let a_nb_sol = a_solver.nb_solutions();
        assert!(a_nb_sol >= 1, "Should find at least one solution for simple case");
        for i in 1..=a_nb_sol {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -10000.0, 10000.0, -10000.0, 10000.0, 100000.0);
            verify_tangency_constraints(&a_sol, &a_input1, &a_input2, &a_input3, i, 1e-6);
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, Concentric_Circles_EdgeCase)
    #[test]
    fn concentric_circles_edge_case() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 1.0);
        let a_qual2 = create_qualified_circle(0.0, 0.0, 3.0);
        let a_qual3 = create_qualified_circle(10.0, 0.0, 2.0);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(
            a_solver.is_done(),
            "Algorithm should complete successfully even when no solutions exist"
        );
        let a_nb_sol = a_solver.nb_solutions();
        assert_eq!(a_nb_sol, 0, "This concentric configuration geometrically has no solutions");
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, SmallCircles_PrecisionTest)
    #[test]
    fn small_circles_precision_test() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 0.01);
        let a_qual2 = create_qualified_circle(0.1, 0.0, 0.01);
        let a_qual3 = create_qualified_circle(0.05, 0.08, 0.01);
        let a_input1 = Circle2d::new(DVec2::new(0.0, 0.0), 0.01);
        let a_input2 = Circle2d::new(DVec2::new(0.1, 0.0), 0.01);
        let a_input3 = Circle2d::new(DVec2::new(0.05, 0.08), 0.01);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Small circles should be handled correctly");

        let a_nb_sol = a_solver.nb_solutions();
        assert!(a_nb_sol >= 1, "Should find at least one solution for small circles");
        for i in 1..=a_nb_sol {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -10.0, 10.0, -10.0, 10.0, 10.0);
            verify_tangency_constraints(&a_sol, &a_input1, &a_input2, &a_input3, i, 1e-3);
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, LargeCircles_ScalingTest)
    #[test]
    fn large_circles_scaling_test() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 1000.0);
        let a_qual2 = create_qualified_circle(5000.0, 0.0, 1500.0);
        let a_qual3 = create_qualified_circle(2500.0, 4000.0, 800.0);
        let a_input1 = Circle2d::new(DVec2::new(0.0, 0.0), 1000.0);
        let a_input2 = Circle2d::new(DVec2::new(5000.0, 0.0), 1500.0);
        let a_input3 = Circle2d::new(DVec2::new(2500.0, 4000.0), 800.0);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Large circles should be handled correctly");

        let a_nb_sol = a_solver.nb_solutions();
        assert!(a_nb_sol >= 1, "Should find at least one solution for large circles");
        for i in 1..=a_nb_sol {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -10000.0, 10000.0, -10000.0, 10000.0, 50000.0);
            verify_tangency_constraints(&a_sol, &a_input1, &a_input2, &a_input3, i, 1.0);
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, LinearConfiguration_GeometricTest)
    #[test]
    fn linear_configuration_geometric_test() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 1.0);
        let a_qual2 = create_qualified_circle(5.0, 0.0, 1.5);
        let a_qual3 = create_qualified_circle(10.0, 0.0, 1.2);
        let a_input1 = Circle2d::new(DVec2::new(0.0, 0.0), 1.0);
        let a_input2 = Circle2d::new(DVec2::new(5.0, 0.0), 1.5);
        let a_input3 = Circle2d::new(DVec2::new(10.0, 0.0), 1.2);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Linear configuration should be solvable");

        let a_nb_sol = a_solver.nb_solutions();
        assert!(a_nb_sol >= 1, "Should find at least one solution for linear configuration");
        for i in 1..=a_nb_sol {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -10000.0, 10000.0, -10000.0, 10000.0, 100000.0);
            verify_tangency_constraints(&a_sol, &a_input1, &a_input2, &a_input3, i, 1e-6);
        }
    }

    // TEST_F(Geom2dGcc_Circ2d3TanTest, TouchingCircles_DegenerateCase)
    #[test]
    fn touching_circles_degenerate_case() {
        let my_tolerance = CONFUSION;
        let a_qual1 = create_qualified_circle(0.0, 0.0, 2.0);
        let a_qual2 = create_qualified_circle(4.0, 0.0, 2.0);
        let a_qual3 = create_qualified_circle(2.0, 5.0, 1.5);
        let a_input1 = Circle2d::new(DVec2::new(0.0, 0.0), 2.0);
        let a_input2 = Circle2d::new(DVec2::new(4.0, 0.0), 2.0);
        let a_input3 = Circle2d::new(DVec2::new(2.0, 5.0), 1.5);

        let a_solver = Circ2d3Tan::new_3_curves(&a_qual1, &a_qual2, &a_qual3, my_tolerance, 0.0, 0.0, 0.0);
        assert!(a_solver.is_done(), "Touching circles configuration should be solvable");

        let a_nb_sol = a_solver.nb_solutions();
        assert!(a_nb_sol >= 1, "Should find at least one solution for touching circles case");
        for i in 1..=a_nb_sol {
            let a_sol = a_solver.this_solution(i);
            verify_solution_validity(&a_sol, i, -10000.0, 10000.0, -10000.0, 10000.0, 100000.0);
            verify_tangency_constraints(&a_sol, &a_input1, &a_input2, &a_input3, i, 1e-6);
        }
    }
}
