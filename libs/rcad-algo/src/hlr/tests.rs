//! Stage unit tests for the TKHLR port (OCCT numeric anchors).

use super::intrv::{Intervals, Interval, Position, REAL_FIRST, REAL_LAST};
use super::top_bas::TestInterference;
use super::top_cnx::EdgeFaceTransition;
use super::algo::edges_block::{EdgesBlock, MinMaxIndices};
use super::algo::hlr_algo::HLRAlgo;
use super::algo::projector::Projector;
use glam::{DVec2, DVec3};
use rcad_kernel::math::gp::{Ax2, Trsf};
use rcad_kernel::topods::Orientation;

fn interval(s: f64, e: f64) -> Interval {
    Interval::from_bounds(s, e)
}

/// OCCT Intrv_Interval default ctor: the whole real line with zero
/// tolerances (Epsilon(RealFirst()) = RealFirst - nextafter(RealFirst,
/// RealFirst) = 0).
#[test]
fn intrv_default_interval_is_whole_line() {
    let i = Interval::new();
    assert_eq!(i.start(), REAL_FIRST);
    assert_eq!(i.end(), REAL_LAST);
    assert_eq!(i.tol_start(), 0.0);
    assert_eq!(i.tol_end(), 0.0);
    // Everything is inside an infinite-tolerance interval.
    assert_eq!(i.position(&interval(10.0, 20.0)), Position::Enclosing);
}

/// OCCT Intrv_Interval::Position — representative branches of the 13-way
/// classification (Intervals.cxx header comment table).
#[test]
fn intrv_position_classification() {
    let other = interval(10.0, 20.0);
    assert_eq!(interval(1.0, 5.0).position(&other), Position::Before);
    assert_eq!(interval(1.0, 10.0).position(&other), Position::JustBefore);
    assert_eq!(
        interval(1.0, 15.0).position(&other),
        Position::OverlappingAtStart
    );
    assert_eq!(
        interval(1.0, 20.0).position(&other),
        Position::JustEnclosingAtEnd
    );
    assert_eq!(interval(1.0, 30.0).position(&other), Position::Enclosing);
    assert_eq!(
        interval(10.0, 14.0).position(&other),
        Position::JustOverlappingAtStart
    );
    assert_eq!(interval(10.0, 20.0).position(&other), Position::Similar);
    assert_eq!(
        interval(10.0, 30.0).position(&other),
        Position::JustEnclosingAtStart
    );
    assert_eq!(interval(15.0, 17.0).position(&other), Position::Inside);
    assert_eq!(
        interval(15.0, 25.0).position(&other),
        Position::OverlappingAtEnd
    );
    assert_eq!(interval(20.0, 25.0).position(&other), Position::JustAfter);
    assert_eq!(interval(25.0, 30.0).position(&other), Position::After);
}

/// OCCT Intrv_Interval FuseAtStart / CutAtStart semantics
/// (Intrv_Interval.lxx L97-144).
#[test]
fn intrv_fuse_and_cut_at_start() {
    // Fuse: the start moves to the other start (min of the two starts).
    let mut i = interval(5.0, 30.0);
    i.fuse_at_start(2.0, 0.0);
    assert_eq!(i.start(), 2.0);
    // Cut: the start moves right to the max.
    let mut j = interval(5.0, 30.0);
    j.cut_at_start(8.0, 0.0);
    assert_eq!(j.start(), 8.0);
    // Fusing at the start of a RealFirst interval is a no-op.
    let mut k = Interval::new();
    k.fuse_at_start(2.0, 0.0);
    assert_eq!(k.start(), REAL_FIRST);
}

/// OCCT Intrv_Intervals Subtract / Unite / Intersect algebra
/// (Intrv_Intervals.cxx L68-248).
#[test]
fn intrv_intervals_algebra() {
    // Subtract: [0,30] minus [10,20] -> {[0,10], [20,30]}.
    let mut is = Intervals::from_interval(interval(0.0, 30.0));
    is.subtract(&interval(10.0, 20.0));
    assert_eq!(is.nb_intervals(), 2);
    assert_eq!(is.value(1).start(), 0.0);
    assert_eq!(is.value(1).end(), 10.0);
    assert_eq!(is.value(2).start(), 20.0);
    assert_eq!(is.value(2).end(), 30.0);

    // Unite: appending [25,40] merges nothing (kept sorted, non overlapping).
    is.unite(&interval(40.0, 50.0));
    assert_eq!(is.nb_intervals(), 3);

    // Intersect with [5, 45]: keeps {[5,10], [20,30], [40,45]}.
    is.intersect(&interval(5.0, 45.0));
    assert_eq!(is.nb_intervals(), 3);
    assert_eq!(is.value(1).start(), 5.0);
    assert_eq!(is.value(2).end(), 30.0);
    assert_eq!(is.value(3).start(), 40.0);
    assert_eq!(is.value(3).end(), 45.0);
}

/// OCCT TopBas_TestInterference accessors (TopBas_TestInterference.hxx).
#[test]
fn top_bas_test_interference_accessors() {
    let mut ti = TestInterference::from_parts(
        0.25,
        7,
        Orientation::Forward,
        Orientation::Reversed,
        Orientation::Internal,
    );
    assert_eq!(ti.intersection(), 0.25);
    assert_eq!(ti.boundary(), 7);
    assert_eq!(ti.orientation(), Orientation::Forward);
    assert_eq!(ti.transition(), Orientation::Reversed);
    assert_eq!(ti.boundary_transition(), Orientation::Internal);
    *ti.change_intersection() = 0.5;
    *ti.change_boundary() = 9;
    ti.set_boundary_transition(Orientation::External);
    assert_eq!(ti.intersection(), 0.5);
    assert_eq!(ti.boundary(), 9);
    assert_eq!(ti.boundary_transition(), Orientation::External);
}

/// OCCT TopCnx_EdgeFaceTransition::BoundaryTransition vote
/// (TopCnx_EdgeFaceTransition.cxx L127-134).
#[test]
fn top_cnx_boundary_transition_vote() {
    let mut t = EdgeFaceTransition::new();
    let x = DVec3::X;
    t.reset_linear(x);
    assert_eq!(t.boundary_transition(), Orientation::External);
    t.add_interference(1e-7, x, DVec3::Z, 0.0, Orientation::Forward, Orientation::Forward, Orientation::Forward);
    assert_eq!(t.boundary_transition(), Orientation::Forward);
    t.add_interference(1e-7, x, DVec3::Z, 0.0, Orientation::Reversed, Orientation::Forward, Orientation::Reversed);
    assert_eq!(t.boundary_transition(), Orientation::External);
    t.add_interference(1e-7, x, DVec3::Z, 0.0, Orientation::Reversed, Orientation::Forward, Orientation::Reversed);
    assert_eq!(t.boundary_transition(), Orientation::Reversed);
}

/// OCCT HLRAlgo_EdgesBlock flag bits (hxx L91-140).
#[test]
fn hlr_algo_edges_block_flags() {
    let mut b = EdgesBlock::new(3);
    b.set_edge(2, 41);
    assert_eq!(b.edge(2), 41);
    assert_eq!(b.nb_edges(), 3);
    b.set_orientation(1, Orientation::Reversed);
    assert_eq!(b.orientation(1), Orientation::Reversed);
    assert!(!b.out_line(1));
    b.set_out_line(1, true);
    b.set_internal(1, true);
    b.set_double(1, true);
    b.set_iso_line(1, true);
    assert!(b.out_line(1) && b.internal(1) && b.double(1) && b.iso_line(1));
    // Orientation survives the other flags being set.
    assert_eq!(b.orientation(1), Orientation::Reversed);
}

/// OCCT HLRAlgo EncodeMinMax / DecodeMinMax round-trip (HLRAlgo.cxx
/// L188-281) — values below 0x8000 pack exactly.
#[test]
fn hlr_algo_min_max_encode_decode_roundtrip() {
    let mut min = MinMaxIndices::default();
    let mut max = MinMaxIndices::default();
    for i in 0..8 {
        min.min[i] = 100 + i as i32;
        min.max[i] = 200 + i as i32;
        max.min[i] = 1000 + i as i32;
        max.max[i] = 2000 + i as i32;
    }
    let mut mm = MinMaxIndices::default();
    HLRAlgo::encode_min_max(&min, &max, &mut mm);
    let mut dmin = MinMaxIndices::default();
    let mut dmax = MinMaxIndices::default();
    HLRAlgo::decode_min_max(&mm, &mut dmin, &mut dmax);
    assert_eq!(dmin, min);
    assert_eq!(dmax, max);
    // AddMinMax: O grows to the union.
    let mut o_min = MinMaxIndices::default();
    let mut o_max = MinMaxIndices::default();
    o_min.min = [5; 8];
    o_max.max = [3000; 8];
    HLRAlgo::add_min_max(&min, &max, &mut o_min, &mut o_max);
    assert_eq!(o_min.min[0], 5);
    assert_eq!(o_max.max[7], 3000);
}

/// OCCT HLRAlgo::UpdateMinMax grows the 16-direction box (HLRAlgo.cxx
/// L123-158).
#[test]
fn hlr_algo_update_min_max() {
    let mut min = [f64::MAX; 16];
    let mut max = [-f64::MAX; 16];
    HLRAlgo::update_min_max(1.0, 2.0, 3.0, &mut min, &mut max);
    HLRAlgo::update_min_max(-1.0, -2.0, -3.0, &mut min, &mut max);
    // Direction 14/15 carry z directly.
    assert_eq!(min[14], -3.0);
    assert_eq!(max[14], 3.0);
    assert_eq!(max[15], 3.0);
    // cos/sin directions 0/1 span [-(x*|c|+y*|s|), +...] — just check the
    // box contains the projections of both points.
    let (c, s) = (
        (0.0f64 * std::f64::consts::PI / 14.0).cos(),
        (0.0f64 * std::f64::consts::PI / 14.0).sin(),
    );
    assert!(min[0] <= c * (-1.0) + s * (-2.0));
    assert!(max[0] >= c * 1.0 + s * 2.0);
}

/// OCCT HLRAlgo_Projector view recognition + projections
/// (HLRAlgo_Projector.cxx L112-317).
#[test]
fn hlr_algo_projector_standard_views() {
    // Top view: identity Ax2.
    let top = Projector::from_ax2(&Ax2::new(
        DVec3::ZERO,
        DVec3::Z,
        DVec3::X,
    ));
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    top.project_xyz(DVec3::new(1.0, 2.0, 3.0), &mut x, &mut y, &mut z);
    assert_eq!((x, y, z), (1.0, 2.0, 3.0));

    // Front view: main dir -Y (rows: X, Z, -Y per TrsfType type 2).
    let front_trsf = {
        let mut t = Trsf::identity();
        t.matrix = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, -1.0, 0.0],
        ];
        t
    };
    let front = Projector::from_trsf(&front_trsf, false, 0.0);
    front.project_xyz(DVec3::new(1.0, 2.0, 3.0), &mut x, &mut y, &mut z);
    assert_eq!((x, y, z), (1.0, 3.0, -2.0));

    // Unknown view falls back to the generic transform path.
    let mut generic_trsf = Trsf::identity();
    generic_trsf.matrix = [
        [0.6, 0.0, 0.8],
        [0.0, 1.0, 0.0],
        [-0.8, 0.0, 0.6],
    ];
    let generic = Projector::from_trsf(&generic_trsf, false, 0.0);
    let mut pout = DVec2::ZERO;
    generic.project_pnt(DVec3::new(1.0, 2.0, 3.0), &mut pout);
    let t = generic.transformation();
    let expected = t.apply(DVec3::new(1.0, 2.0, 3.0));
    assert_eq!(pout, DVec2::new(expected.x, expected.y));
}

/// OCCT HLRAlgo_Projector perspective formula (HLRAlgo_Projector.cxx
/// L228-236 and the D1 header comment L31-36).
#[test]
fn hlr_algo_projector_perspective() {
    let mut generic_trsf = Trsf::identity();
    generic_trsf.matrix = [
        [0.6, 0.0, 0.8],
        [0.0, 1.0, 0.0],
        [-0.8, 0.0, 0.6],
    ];
    let focus = 10.0;
    let p = Projector::from_trsf(&generic_trsf, true, focus);
    let mut pout = DVec2::ZERO;
    p.project_pnt(DVec3::new(1.0, 2.0, 3.0), &mut pout);
    let t = p.transformation();
    let p2 = t.apply(DVec3::new(1.0, 2.0, 3.0));
    let r = 1.0 - p2.z / focus;
    assert!((pout.x - p2.x / r).abs() < 1e-12);
    assert!((pout.y - p2.y / r).abs() < 1e-12);
    assert!(p.perspective());
    assert_eq!(p.focus(), focus);
}

/// OCCT HLRAlgo_Projector::Shoot — the eye ray through (X, Y), transformed
/// by the inverted transformation (cxx L346-359).
#[test]
fn hlr_algo_projector_shoot() {
    let top = Projector::from_ax2(&Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X));
    let l = top.shoot(0.5, -0.25);
    // Orthographic: the ray starts at (X, Y, 0) and runs along -Z.
    assert_eq!(l.pos, DVec3::new(0.5, -0.25, 0.0));
    assert!((l.dir + DVec3::Z).length() < 1e-12);
}
