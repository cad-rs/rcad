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

// ---- Stage 1: HLRAlgo data structures ----

use super::algo::bi_point::BiPoint;
use super::algo::coincidence::Coincidence;
use super::algo::edge_iterator::EdgeIterator;
use super::algo::edge_status::EdgeStatus;
use super::algo::interference::Interference;
use super::algo::intersection::Intersection;
use rcad_kernel::topods::State;

/// OCCT HLRAlgo_EdgeStatus hide/show state machine (EdgeStatus.cxx
/// L103-123): hiding on-face intervals subtracts from the visible sequence.
#[test]
fn hlr_algo_edge_status_state_machine() {
    // OCCT Interval(Start, TolStart, End, TolEnd) raises tolerances to the
    // representable epsilon of each bound (Interval.cxx L64-73).
    let eps30 = super::intrv::epsilon(30.0) as f32;
    let eps60 = super::intrv::epsilon(60.0) as f32;
    let eps100 = super::intrv::epsilon(100.0) as f32;
    let mut st = EdgeStatus::from_bounds(0.0, 0.0, 100.0, 0.0);
    assert!(st.all_visible());
    assert_eq!(st.nb_visible_part(), 1);
    assert_eq!(st.visible_part(1), (0.0, 0.0, 100.0, 0.0));

    // Hide the middle: two visible parts remain.
    st.hide(30.0, 0.0, 60.0, 0.0, false, false); // OCCT: subtract runs when OnFace == false
    assert!(!st.all_visible() && !st.all_hidden());
    assert_eq!(st.nb_visible_part(), 2);
    assert_eq!(st.visible_part(1), (0.0, 0.0, 30.0, eps30));
    assert_eq!(st.visible_part(2), (60.0, eps60, 100.0, eps100));

    // Hide everything: all hidden.
    st.hide(0.0, 0.0, 100.0, 0.0, false, false);
    assert!(st.all_hidden());
    assert_eq!(st.nb_visible_part(), 0);

    // ShowAll restores the full visible state.
    st.show_all();
    assert!(st.all_visible());
    assert_eq!(st.nb_visible_part(), 1);

    // HideAll flags hidden without touching the interval sequence.
    st.hide_all();
    assert!(st.all_hidden());

    // Hiding with OnFace == false is a no-op (cxx L110 guard).
    let mut st2 = EdgeStatus::from_bounds(0.0, 0.0, 50.0, 0.0);
    st2.hide(10.0, 0.0, 20.0, 0.0, true, false);
    assert!(st2.all_visible());
    assert_eq!(st2.nb_visible_part(), 1);
}

/// OCCT HLRAlgo_EdgeIterator over the visible parts (EdgeIterator.cxx
/// L42-94, lxx L41-70): the cached hidden interval runs between the end of
/// one visible part and the start of the next.
#[test]
fn hlr_algo_edge_iterator_visible_parts() {
    let eps30 = super::intrv::epsilon(30.0) as f32;
    let eps60 = super::intrv::epsilon(60.0) as f32;
    let eps100 = super::intrv::epsilon(100.0) as f32;
    let mut st = EdgeStatus::from_bounds(0.0, 0.0, 100.0, 0.0);
    st.hide(30.0, 0.0, 60.0, 0.0, false, false); // OCCT: subtract runs when OnFace == false

    let mut it = EdgeIterator::new();
    it.init_visible(&st);
    assert!(it.more_visible());
    assert_eq!(it.visible(), (0.0, 0.0, 30.0, eps30));
    it.next_visible();
    assert!(it.more_visible());
    assert_eq!(it.visible(), (60.0, eps60, 100.0, eps100));
    it.next_visible();
    assert!(!it.more_visible());

    // Hidden iterator: AllHidden == false -> 2 hidden intervals between the
    // visible parts (cached intervals) plus the outer ranges.
    it.init_hidden(&st);
    let mut seen = Vec::new();
    while it.more_hidden() {
        seen.push(it.hidden());
        it.next_hidden();
    }
    // The first cached interval = [end of part 1 = 30, start of part 2 = 60].
    assert_eq!(seen[0], (30.0, eps30, 60.0, eps60));
}

/// OCCT HLRAlgo_BiPoint flag bits (BiPoint.hxx L193-241) and the
/// Interference/Intersection/Coincidence accessors.
#[test]
fn hlr_algo_bi_point_flags_and_interference() {
    let mut bp = BiPoint::from_flags_bool(
        0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 5, true, false, true, false,
    );
    assert_eq!(bp.indices_ref().shape_index, 5);
    assert!(bp.rg1_line() && !bp.rgn_line() && bp.out_line() && !bp.int_line());
    assert!(!bp.hidden());
    bp.indices().seg_flags = 0;
    bp.set_hidden(true);
    assert!(bp.hidden());

    let inter = Interference::from_parts(
        Intersection::from_parts(
            Orientation::Reversed,
            2,
            3,
            4,
            0.5,
            1e-7,
            State::On,
        ),
        Coincidence::new(),
        Orientation::Forward,
        Orientation::Internal,
        Orientation::External,
    );
    assert_eq!(inter.intersection().parameter(), 0.5);
    assert_eq!(inter.intersection().state(), State::On);
    assert_eq!(inter.transition(), Orientation::Internal);
    assert_eq!(inter.boundary_transition(), Orientation::External);
}

// ---- Stage: the end-to-end HLR pipeline smoke (the OCCT VComputeHLR user
// path, ViewerTest_ObjectCommands.cxx L3299-3331) ----

use super::algo::projector::Projector as SmokeProjector;
use super::brep::algo::Algo as SmokeAlgo;
use super::brep::hider::Hider;
use super::brep::hlr_to_shape::HLRToShape as SmokeHLRToShape;
use rcad_kernel::geom::{Circle3, Curve3, Line3, Plane, Surface3};
use rcad_kernel::math::gp::Ax2 as SmokeAx2;
use rcad_kernel::topo::topods::{BRepBuilder, TShape};

/// A projected result edge summary: the curve kind, the trimmed length and
/// the two end points (the plane image of the 2D support: z = 0).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SegSummary {
    Line {
        len: f64,
        p1: glam::DVec3,
        p2: glam::DVec3,
    },
    Circle {
        len: f64,
        p1: glam::DVec3,
        p2: glam::DVec3,
        center: glam::DVec3,
        radius: f64,
    },
    Ellipse {
        len: f64,
        p1: glam::DVec3,
        p2: glam::DVec3,
        major: f64,
        minor: f64,
    },
    Other,
}

pub(crate) fn seg_summary(e: &rcad_kernel::topods::Shape) -> SegSummary {
    match &*e.data {
        TShape::Edge(ed) => {
            let v1 = ed
                .first
                .as_vertex()
                .map(|v| v.point)
                .unwrap_or(glam::DVec3::ZERO);
            let v2 = ed
                .last
                .as_vertex()
                .map(|v| v.point)
                .unwrap_or(glam::DVec3::ZERO);
            let [u1, u2] = ed.range;
            match &ed.curve {
                Some(Curve3::Line(l)) => {
                    let d = l.direction.length();
                    SegSummary::Line {
                        len: d * (u2 - u1),
                        p1: v1,
                        p2: v2,
                    }
                }
                Some(Curve3::Circle(c)) => SegSummary::Circle {
                    len: c.radius * (u2 - u1),
                    p1: v1,
                    p2: v2,
                    center: c.center,
                    radius: c.radius,
                },
                Some(Curve3::Ellipse(el)) => SegSummary::Ellipse {
                    len: f64::NAN,
                    p1: v1,
                    p2: v2,
                    major: el.major_radius,
                    minor: el.minor_radius,
                },
                _ => SegSummary::Other,
            }
        }
        _ => SegSummary::Other,
    }
}

/// The compound summary: the flattened edge summaries of a result compound.
pub(crate) fn compound_edges(
    s: &rcad_kernel::topods::Shape,
    out: &mut Vec<SegSummary>,
) {
    match &*s.data {
        TShape::Compound(c) => {
            for ch in c {
                compound_edges(ch, out);
            }
        }
        TShape::Edge(_) => out.push(seg_summary(s)),
        _ => {}
    }
}

/// The box fixture: the OCCT `box b 10 20 30` equivalent — six planar faces
/// with outward normals, twelve sharp line edges, one closed shell in a
/// solid.  Built with the kernel BRepBuilder (the plane_brep fixture form).
/// Returns the owning BRep (the kernel context the HLR loads read through)
/// and the solid shape.
pub(crate) fn smoke_box_solid() -> (rcad_kernel::BRep, rcad_kernel::topods::Shape) {
    let mut brep = rcad_kernel::BRep::new();
    let mut b = BRepBuilder::new();
    let p = |x: f64, y: f64, z: f64| glam::DVec3::new(x, y, z);
    let vs = [
        b.add_vertex(&mut brep, p(0.0, 0.0, 0.0), 1e-7),  // 0 A
        b.add_vertex(&mut brep, p(10.0, 0.0, 0.0), 1e-7), // 1 B
        b.add_vertex(&mut brep, p(10.0, 20.0, 0.0), 1e-7),// 2 C
        b.add_vertex(&mut brep, p(0.0, 20.0, 0.0), 1e-7), // 3 D
        b.add_vertex(&mut brep, p(0.0, 0.0, 30.0), 1e-7), // 4 E
        b.add_vertex(&mut brep, p(10.0, 0.0, 30.0), 1e-7),// 5 F
        b.add_vertex(&mut brep, p(10.0, 20.0, 30.0), 1e-7),// 6 G
        b.add_vertex(&mut brep, p(0.0, 20.0, 30.0), 1e-7),// 7 H
    ];
    // the 12 edges keyed by (i, j) corner indices.
    let mut em: std::collections::HashMap<(usize, usize), rcad_kernel::topods::Shape> =
        std::collections::HashMap::new();
    let corner = |i: usize| match i {
        0 => p(0.0, 0.0, 0.0),
        1 => p(10.0, 0.0, 0.0),
        2 => p(10.0, 20.0, 0.0),
        3 => p(0.0, 20.0, 0.0),
        4 => p(0.0, 0.0, 30.0),
        5 => p(10.0, 0.0, 30.0),
        6 => p(10.0, 20.0, 30.0),
        _ => p(0.0, 20.0, 30.0),
    };
    for i in 0..8 {
        for j in (i + 1)..8 {
            let is_edge = matches!(
                (i, j),
                (0, 1) | (1, 2) | (2, 3) | (0, 3) | (4, 5) | (5, 6) | (6, 7) | (4, 7)
                    | (0, 4) | (1, 5) | (2, 6) | (3, 7)
            );
            if !is_edge {
                continue;
            }
            let a = corner(i);
            let c = corner(j);
            // OCCT TopoDS_Edge: the last vertex child is stored REVERSED.
            let mut v_last = vs[j].clone();
            v_last.orientation = rcad_kernel::topods::Orientation::Reversed;
            let e = b.add_edge(
                &mut brep,
                Some(Curve3::Line(Line3 {
                    origin: a,
                    direction: (c - a).normalize(),
                })),
                vs[i].clone(),
                v_last,
                [0.0, (c - a).length()],
            );
            em.insert((i, j), e);
        }
    }
    // the oriented edge instance (i -> j): Forward along (i, j), Reversed
    // along (j, i).
    let dir_edge = |i: usize, j: usize| -> rcad_kernel::topods::Shape {
        let (k, l) = if i < j { (i, j) } else { (j, i) };
        let mut e = em[&(k, l)].clone();
        if i > j {
            e.orientation = rcad_kernel::topods::Orientation::Reversed;
        }
        e
    };
    // the six faces: (surface, wire corners CCW seen from outside).
    let plane = |o: glam::DVec3, n: glam::DVec3, u: glam::DVec3, v: glam::DVec3| {
        Surface3::Plane(Plane {
            origin: o,
            normal: n,
            u_dir: u,
            v_dir: v,
        })
    };
    let face_defs: [(Surface3, [usize; 4], [f64; 4]); 6] = [
        // bottom z=0, outward -Z: A D C B
        (
            plane(p(0.0, 0.0, 0.0), p(0.0, 0.0, -1.0), p(1.0, 0.0, 0.0), p(0.0, -1.0, 0.0)),
            [0, 3, 2, 1],
            [0.0, 20.0, -10.0, 0.0],
        ),
        // top z=30, outward +Z: E F G H
        (
            plane(p(0.0, 0.0, 30.0), p(0.0, 0.0, 1.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0)),
            [4, 5, 6, 7],
            [0.0, 10.0, 0.0, 20.0],
        ),
        // front y=0, outward -Y: A B F E
        (
            plane(p(0.0, 0.0, 0.0), p(0.0, -1.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 0.0, 1.0)),
            [0, 1, 5, 4],
            [0.0, 10.0, 0.0, 30.0],
        ),
        // back y=20, outward +Y: C D H G
        (
            plane(p(0.0, 20.0, 0.0), p(0.0, 1.0, 0.0), p(-1.0, 0.0, 0.0), p(0.0, 0.0, 1.0)),
            [2, 3, 7, 6],
            [0.0, 10.0, 0.0, 30.0],
        ),
        // left x=0, outward -X: A E H D
        (
            plane(p(0.0, 0.0, 0.0), p(-1.0, 0.0, 0.0), p(0.0, -1.0, 0.0), p(0.0, 0.0, 1.0)),
            [0, 4, 7, 3],
            [0.0, 20.0, 0.0, 30.0],
        ),
        // right x=10, outward +X: B C G F
        (
            plane(p(10.0, 0.0, 0.0), p(1.0, 0.0, 0.0), p(0.0, 1.0, 0.0), p(0.0, 0.0, 1.0)),
            [1, 2, 6, 5],
            [0.0, 20.0, 0.0, 30.0],
        ),
    ];
    let mut faces = Vec::new();
    for (surface, wire_corners, uv) in &face_defs {
        let mut wire_edges = Vec::new();
        for k in 0..4 {
            wire_edges.push(dir_edge(wire_corners[k], wire_corners[(k + 1) % 4]));
        }
        let wire = brep.add_twire(wire_edges);
        faces.push(brep.add_tface(
            Some(surface.clone()),
            wire,
            Vec::new(),
            None,
            Some(*uv),
            Vec::new(),
            true,
        ));
    }
    // the pcurves: the OCCT primitives carry a pcurve per (edge, face)
    // (BRep_Builder::UpdateEdge); the rcad FClass2d init early-returns
    // without one (the OCCT CurveOnPlane projection fallback is not landed,
    // fclass2d_topol.rs L322 — recorded as a pipeline gap), so the smoke
    // registers the isometric plane images of the 12 edges over their two
    // adjacent faces each.
    let plane_frame = |surface: &Surface3| -> (glam::DVec3, glam::DVec3, glam::DVec3) {
        match surface {
            Surface3::Plane(p) => (p.origin, p.u_dir, p.v_dir),
            _ => panic!("plane expected"),
        }
    };
    for (fi, (surface, wire_corners, _uv)) in face_defs.iter().enumerate() {
        let (o, u, v) = plane_frame(surface);
        for k in 0..4 {
            let (i, j) = (wire_corners[k], wire_corners[(k + 1) % 4]);
            let (kk, ll) = if i < j { (i, j) } else { (j, i) };
            let edge = em[&(kk, ll)].clone();
            let a3 = corner(i);
            let b3 = corner(j);
            let pa = glam::DVec2::new((a3 - o).dot(u), (a3 - o).dot(v));
            let pb = glam::DVec2::new((b3 - o).dot(u), (b3 - o).dot(v));
            let pc = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: pa,
                direction: (pb - pa).normalize(),
            });
            b.add_pcurve(
                &mut brep,
                edge,
                faces[fi].clone(),
                pc,
                0.0,
                (pb - pa).length(),
            );
        }
    }
    let shell = brep.add_tshell(faces.clone());
    let solid = brep.add_tsolid(vec![shell]);
    (brep, solid)
}

/// The VComputeHLR run (the exact-Algo branch): Add + Projector + Update +
/// Hide, then the HLRToShape extraction.  Returns (v, rg1v, rgnv, outline, h)
/// compounds — the sharp (typ 5), Rg1Line (typ 3), RgNLine (typ 4) visible
/// compounds, the outlines (typ 2) and the hidden sharp compound (OCCT
/// HLRBRep_HLRToShape.hxx L71-104).
///
/// The owning BRep rides with Add as the session kernel context — the
/// InternalAlgo Update builds every shape DS over it and feeds it to the
/// `myDS->Update(myProj)` phase (the OCCT global TShape graph stand-in the
/// loaded shapes live in).
pub(crate) fn smoke_run_hlr(
    solid: &rcad_kernel::topods::Shape,
    brep: rcad_kernel::BRep,
) -> (
    rcad_kernel::topods::Shape,
    rcad_kernel::topods::Shape,
    rcad_kernel::topods::Shape,
    rcad_kernel::topods::Shape,
    rcad_kernel::topods::Shape,
) {
    // OCCT L3303: aHlrAlgo->Projector(aProjector) — the V3d_XposYnegZpos
    // equivalent frame: dir (1,-1,1), up (-1,1,2); X = up ^ dir = (1,1,0).
    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));

    let mut algo = SmokeAlgo::new();
    // OCCT L3302: aHlrAlgo->Add(aSh, aNbIsolines);
    algo.add(&std::sync::Arc::new(brep), solid, 0);
    algo.set_projector(&proj);
    // OCCT L3304: aHlrAlgo->Update();
    algo.update();
    // OCCT L3305: aHlrAlgo->Hide();
    algo.hide();
    // OCCT L3307-3313: the HLRToShape filters.
    let mut hts = SmokeHLRToShape::new(&mut algo);
    let v = hts.v_compound();
    let rg1v = hts.rg1_line_v_compound();
    let rgnv = hts.rg_n_line_v_compound();
    let outline = hts.out_line_v_compound();
    let h = hts.h_compound();
    (v, rg1v, rgnv, outline, h)
}

/// The eight HLRToShape filter outputs of one exact-algorithm run
/// (ViewerTest_ObjectCommands.cxx L3308-3326).
#[allow(dead_code)]
pub(crate) struct VComputeFilters {
    /// aVisible[HLRBRep_Sharp] (L3308).
    pub v: rcad_kernel::topods::Shape,
    /// aVisible[HLRBRep_OutLine] (L3309).
    pub outline_v: rcad_kernel::topods::Shape,
    /// aVisible[HLRBRep_Rg1Line] (L3310).
    pub rg1_v: rcad_kernel::topods::Shape,
    /// aVisible[HLRBRep_IsoLine] (L3315).
    pub iso_v: rcad_kernel::topods::Shape,
    /// aHidden[HLRBRep_Sharp] (L3319).
    pub h: rcad_kernel::topods::Shape,
    /// aHidden[HLRBRep_OutLine] (L3320).
    pub outline_h: rcad_kernel::topods::Shape,
    /// aHidden[HLRBRep_Rg1Line] (L3321).
    pub rg1_h: rcad_kernel::topods::Shape,
    /// aHidden[HLRBRep_IsoLine] (L3326).
    pub iso_h: rcad_kernel::topods::Shape,
}

/// The VComputeHLR exact-algorithm arm (ViewerTest_ObjectCommands.cxx
/// L3301-3326): Add(aSh, aNbIsolines) / Projector / Update / Hide, then the
/// eight HLRToShape filter calls in the OCCT order.  The RgN filters stay
/// unextracted (toShowCNEdges == false, L3311-3313 / L3322-3324).
pub(crate) fn smoke_run_hlr_vcompute(
    solid: &rcad_kernel::topods::Shape,
    brep: rcad_kernel::BRep,
) -> VComputeFilters {
    // OCCT L3296-3300: the projector built from the DRAW view frame — the
    // V3d_XposYnegZpos equivalent: dir (1,-1,1), up (-1,1,2).
    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));

    let mut algo = SmokeAlgo::new();
    // OCCT L3302: aHlrAlgo->Add(aSh, aNbIsolines);
    algo.add(&std::sync::Arc::new(brep), solid, 0);
    algo.set_projector(&proj);
    // OCCT L3304: aHlrAlgo->Update();
    algo.update();
    // OCCT L3305: aHlrAlgo->Hide();
    algo.hide();

    let mut hts = SmokeHLRToShape::new(&mut algo);
    // OCCT L3308-3315: aVisible[...] = filter calls in order.
    let v = hts.v_compound();
    let outline_v = hts.out_line_v_compound();
    let rg1_v = hts.rg1_line_v_compound();
    let iso_v = hts.iso_line_v_compound();
    // OCCT L3319-3326: aHidden[...] = filter calls in order.
    let h = hts.h_compound();
    let outline_h = hts.out_line_h_compound();
    let rg1_h = hts.rg1_line_h_compound();
    let iso_h = hts.iso_line_h_compound();
    VComputeFilters {
        v,
        outline_v,
        rg1_v,
        iso_v,
        h,
        outline_h,
        rg1_h,
        iso_h,
    }
}

/// OCCT probe truth (session-18, bug25813_1): after the Load/ExploreFace
/// pass the small-cylinder lateral face carries a two-edge IntL wire (the
/// silhouette generator pieces) whose w_edge records hold the Internal
/// flag, and `Data::Update`'s OrientOutLine (Data.cxx L1746-1877) must have
/// rewritten every such w_edge orientation to FORWARD/REVERSED before Hide.
/// That rewrite is the precondition for the NextInterference F/R guard
/// (Data.cxx L1277) to intersect the IntL wedges at all - without it the
/// two IN interferences that hide the big-top rim middle (16.9334) never
/// come into existence.
#[test]
fn bug25813_1_intl_wedges_present_and_reoriented_after_update() {
    // pcylinder cc 10 30 / pcylinder cc2 8 50 / ttranslate cc2 0 0 2 / bfuse
    let cc = rcad_modeling::make_cylinder_brep(
        glam::DVec3::ZERO,
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        10.0,
        30.0,
    )
    .expect("pcylinder cc 10 30");
    let cc2 = rcad_modeling::make_cylinder_brep(
        glam::DVec3::new(0.0, 0.0, 2.0),
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        8.0,
        50.0,
    )
    .expect("pcylinder cc2 8 50");
    let fused = crate::fuse(&cc, &cc2).expect("bfuse a cc cc2");
    let solid = fused
        .tshapes
        .iter()
        .enumerate()
        .rev()
        .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
        .map(|(i, ts)| {
            rcad_kernel::topods::Shape::from_parts(
                ts.clone(),
                i,
                0,
                rcad_kernel::topods::Orientation::Forward,
            )
        })
        .expect("root solid");

    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));
    let mut algo = SmokeAlgo::new();
    algo.add(&std::sync::Arc::new(fused), &solid, 0);
    algo.set_projector(&proj);
    // OCCT L3304: aHlrAlgo->Update() — runs the withOutL/OrientOutLine pass.
    algo.update();

    let ds = algo.data_structure().expect("data structure");
    let mut intl_total = 0usize;
    let mut intl_unrewritten: Vec<(usize, usize, usize, Orientation)> = Vec::new();
    for (fi, fd) in ds.f_data_array().iter().enumerate() {
        let wb = fd.wires();
        let wb = unsafe { &mut *wb };
        for wi in 1..=wb.nb_wires() {
            let eb = wb.wire(wi);
            for ei in 1..=eb.nb_edges() {
                if eb.internal(ei) {
                    intl_total += 1;
                    let ori = eb.orientation(ei);
                    if !matches!(ori, Orientation::Forward | Orientation::Reversed) {
                        intl_unrewritten.push((fi + 1, wi, ei, ori));
                    }
                }
            }
        }
    }
    assert!(
        intl_total >= 2,
        "no IntL w_edges registered (OCCT X19: the outlined faces carry \
         two-edge IntL wires with the Internal flag)"
    );
    assert!(
        intl_unrewritten.is_empty(),
        "IntL w_edges still not reoriented to FORWARD/REVERSED after Update \
         (OCCT X20 requires F/R via OrientOutLine): {intl_unrewritten:?}"
    );

    // OCCT X17 truth: the rim-vs-silhouette crossings sit at 3D params
    // 1.4288992721907325 / 3.2834897081939571 on the big-top rim edge, i.e.
    // 2D params 0.6435/2.4981 after the myOX=-0.78539816339744828 offset.
    // The intersector only reports crossings inside the edge's 2D domain
    // [FirstParameter+myOX, LastParameter+myOX]; dump every analytic edge's
    // type and 3D range so the domain windows are visible in the test log.
    for (ei, ed) in ds.e_data_array().iter().enumerate() {
        let ec = ed.geometry();
        let t = ec.get_type();
        if matches!(
            t,
            rcad_kernel::base::proj_lib::CurveType::Circle
                | rcad_kernel::base::proj_lib::CurveType::Ellipse
        ) {
            let p1 = ec.parameter_3d(ec.first_parameter());
            let p2 = ec.parameter_3d(ec.last_parameter());
            println!(
                "INTLDS edge={} type={:?} range3d=({p1:.17}, {p2:.17})",
                ei + 1,
                t
            );
        }
    }
}

/// OCCT exact_hlr/bug25813_1 truth (session-18 probe X17): the full
/// pipeline hides the big-top rim far-upper middle (the two IN
/// interferences at 1.4288992721907325 / 3.2834897081939571 against the
/// small-cylinder lateral face's IntL silhouette wedges), giving a hidden
/// compound of exactly 3 arcs totalling 62.3361 projected length
/// (25.2237 big-bottom far half + 20.1790 z=30 inner-rim far half +
/// 16.9334 big-top far-upper middle) and a visible compound of 15 edges /
/// 204.1903.  Regression for the Classify LevelFlag Z-gate, which in OCCT
/// (Data.cxx L2026-2040, L2059-2073, L2087-2101, L2159-2173) has 15
/// comparisons - the symmetric MinMaxVert.Max[7]-iFaceMinMax.Min[7] term
/// does not exist.
#[test]
fn bug25813_1_hidden_compound_three_arcs() {
    let cc = rcad_modeling::make_cylinder_brep(
        glam::DVec3::ZERO,
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        10.0,
        30.0,
    )
    .expect("pcylinder cc 10 30");
    // the generated harness translates cc2 via apply_transform — measure
    // whether that fixture variant leaks the seam into the Iso compound.
    let cc2 = rcad_modeling::make_cylinder_brep(
        glam::DVec3::ZERO,
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        8.0,
        50.0,
    )
    .expect("pcylinder cc2 8 50");
    let cc2 = {
        let mut shape = cc2;
        shape.apply_transform(glam::DAffine3::from_translation(glam::DVec3::new(0.0, 0.0, 2.0)));
        shape
    };
    let fused = crate::fuse(&cc, &cc2).expect("bfuse a cc cc2");
    let solid = fused
        .tshapes
        .iter()
        .enumerate()
        .rev()
        .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
        .map(|(i, ts)| {
            rcad_kernel::topods::Shape::from_parts(
                ts.clone(),
                i,
                0,
                rcad_kernel::topods::Orientation::Forward,
            )
        })
        .expect("root solid");
    let brep = fused;
    // inline the pipeline (smoke_run_hlr's flow) so the Iso compound can be
    // measured too — the harness sums V + OutLine + Rg1 + Iso.
    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));
    let mut algo = SmokeAlgo::new();
    algo.add(&std::sync::Arc::new(brep), &solid, 0);
    algo.set_projector(&proj);
    algo.update();
    algo.hide();
    let mut hts = SmokeHLRToShape::new(&mut algo);
    let v = hts.v_compound();
    let rg1v = hts.rg1_line_v_compound();
    let rgnv = hts.rg_n_line_v_compound();
    let outline = hts.out_line_v_compound();
    let iso = hts.iso_line_v_compound();
    let h = hts.h_compound();

    let mut h_edges = Vec::new();
    assert!(!h.is_null(), "hidden compound is null");
    compound_edges(&h, &mut h_edges);
    assert_eq!(h_edges.len(), 3, "hidden edges: {h_edges:?}");
    let mass = super::acceptance::compound_mass(&h);
    assert!(
        (mass - 62.3360975).abs() < 1e-5,
        "hidden projected mass {mass} vs OCCT 62.3360975: {h_edges:?}"
    );

    let mut v_edges = Vec::new();
    assert!(!v.is_null(), "visible compound is null");
    compound_edges(&v, &mut v_edges);
    assert!(!v_edges.is_empty(), "visible compound empty");
    // OCCT's visible result = the 11 sharp V arcs + the 4 outline lines;
    // the RgN seam lines are NOT part of it (gap3: rcad still emits them).
    let vis_mass = super::acceptance::compound_mass(&v)
        + super::acceptance::compound_mass(&rg1v)
        + super::acceptance::compound_mass(&outline);
    assert!(
        (vis_mass - 204.1903269).abs() < 1e-5,
        "visible projected mass {vis_mass} vs OCCT 204.1903269"
    );
    // gap3 witness: OCCT's RgNLineVCompound is empty for this case; rcad
    // still emits the two seam lines (24.4949 + 17.9629 = 42.4578).  Flip
    // to 0 when the seam RgN routing is fixed.
    let rgn_mass = super::acceptance::compound_mass(&rgnv);
    assert!(
        (rgn_mass - 42.457822208241754).abs() < 1e-5,
        "RgN mass {rgn_mass}: expected the two seam lines 42.4578 (gap3)"
    );
    println!(
        "COMPOUND masses: v={:.6} rg1v={:.6} rgnv={:.6} outline={:.6} iso={:.6}",
        super::acceptance::compound_mass(&v),
        super::acceptance::compound_mass(&rg1v),
        rgn_mass,
        super::acceptance::compound_mass(&outline),
        super::acceptance::compound_mass(&iso),
    );
}

/// Diagnostic one level below [`bug25813_1_rim_middle_partially_hidden_by_outlined_faces`]:
/// drive the Data interference loop for the outlined faces and record every
/// surviving/rejected interference of the ellipse edges, to separate
/// "the intersector produces no crossings" from "RejectedPoint drops them".
#[test]
fn bug25813_1_interference_loop_over_outlined_faces() {
    let cc = rcad_modeling::make_cylinder_brep(
        glam::DVec3::ZERO,
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        10.0,
        30.0,
    )
    .expect("pcylinder cc 10 30");
    let cc2 = rcad_modeling::make_cylinder_brep(
        glam::DVec3::new(0.0, 0.0, 2.0),
        glam::DVec3::new(0.0, 0.0, 1.0),
        glam::DVec3::new(1.0, 0.0, 0.0),
        8.0,
        50.0,
    )
    .expect("pcylinder cc2 8 50");
    let fused = crate::fuse(&cc, &cc2).expect("bfuse a cc cc2");
    let solid = fused
        .tshapes
        .iter()
        .enumerate()
        .rev()
        .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
        .map(|(i, ts)| {
            rcad_kernel::topods::Shape::from_parts(
                ts.clone(),
                i,
                0,
                rcad_kernel::topods::Orientation::Forward,
            )
        })
        .expect("root solid");

    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));
    let mut algo = SmokeAlgo::new();
    algo.add(&std::sync::Arc::new(fused), &solid, 0);
    algo.set_projector(&proj);
    algo.update();
    // the real single-shape hide: runs Select + InitEdgeStatus +
    // InitBoundSort (which fills the sorted global edge list phase 2 of
    // MoreEdge iterates) + the per-face Hider passes.
    algo.hide_shape(1);

    for fi in [1usize, 4usize] {
        let ds = algo.data_structure_mut().expect("data structure");
        let mut mst: Vec<(
            rcad_kernel::topods::Shape,
            crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool,
        )> = Vec::new();
        ds.init_edge(fi, &mut mst);
        while ds.more_edge() {
            let e = ds.edge();
            let ec = ds.e_data_array()[e - 1].geometry();
            let t = ec.get_type();
            let analytic = matches!(
                t,
                rcad_kernel::base::proj_lib::CurveType::Circle
                    | rcad_kernel::base::proj_lib::CurveType::Ellipse
            );
            ds.init_interference();
            while ds.more_interference() {
                if ds.rejected_interference() {
                    let above = ds.above_interference();
                    if analytic && above {
                        println!("REJ-ABOVE face={fi} E={e} type={t:?}");
                    }
                } else if analytic {
                    let itf = ds.interference();
                    println!(
                        "HIT face={fi} E={e} type={t:?} st={:?} p={:?}",
                        itf.intersection().state(),
                        itf.intersection().parameter()
                    );
                }
                ds.next_interference();
            }
            // record which face-edges were visited for this analytic LE
            if analytic {
                println!(
                    "VISIT face={fi} E={e} lastFE={} ori={:?} int={}",
                    ds.fe(),
                    ds.fe_ori(),
                    ds.fe_internal()
                );
            }
            ds.next_edge(true);
        }

        // OCCT X17: the Hider classifies the rim midpoint (lf=1,
        // param=2.3561944901923448) against this face; the OCCT verdict is
        // IN (state=0, level=2). Probe the same verdict here.
        let ds2 = algo.data_structure().expect("data structure");
        let _ = ds2;
    }
    {
        let ds = algo.data_structure_mut().expect("data structure");
        // the big-top rim half candidates: ellipse edges with 3D range
        // [0.78539816339744828, 3.92699081698724139].
        for ei in 0..ds.e_data_array().len() {
            let (p1, p2, is_ell) = {
                let ed = &ds.e_data_array()[ei];
                let ec = ed.geometry();
                (
                    ec.parameter_3d(ec.first_parameter()),
                    ec.parameter_3d(ec.last_parameter()),
                    matches!(
                        ec.get_type(),
                        rcad_kernel::base::proj_lib::CurveType::Ellipse
                    ),
                )
            };
            if is_ell
                && (p1 - 0.78539816339744828).abs() < 1e-9
                && (p2 - 3.92699081698724139).abs() < 1e-9
            {
                let ed_ptr = &ds.e_data_array()[ei] as *const _;
                let mut lvl = 0i32;
                let st = ds.classify(
                    (ei + 1) as i32,
                    unsafe { &*ed_ptr },
                    true,
                    &mut lvl,
                    2.3561944901923448,
                );
                println!(
                    "CLASSIFY edge={ei} range=({p1:.6},{p2:.6}) -> state={st:?} level={lvl}"
                );
                if ei == 1 {
                    // replicate the LevelFlag gate (classify.rs L1842-1914)
                    // word by word for the failing edge.
                    let ec2 = unsafe { &*ed_ptr }.geometry();
                    let tol = unsafe { &*ed_ptr }.tolerance() as f64;
                    let p3 = ec2.value_3d(2.3561944901923448);
                    let mut x = 0.0f64;
                    let mut y = 0.0f64;
                    let mut z = 0.0f64;
                    ds.projector().project_xyz(p3, &mut x, &mut y, &mut z);
                    let mut tmin = [0.0f64; 16];
                    let mut tmax = [0.0f64; 16];
                    HLRAlgo::init_min_max(2.0e100, &mut tmin, &mut tmax);
                    HLRAlgo::update_min_max(x, y, z, &mut tmin, &mut tmax);
                    HLRAlgo::enlarge_min_max(tol, &mut tmin, &mut tmax);
                    let mut vmin = MinMaxIndices::default();
                    let mut vmax = MinMaxIndices::default();
                    crate::hlr::brep::data::reject1(
                        ds.deca(),
                        &tmin,
                        &tmax,
                        ds.sur_d(),
                        &mut vmin,
                        &mut vmax,
                    );
                    let mut mmv = MinMaxIndices::default();
                    HLRAlgo::encode_min_max(&vmin, &vmax, &mut mmv);
                    let mm = ds.i_face_min_max();
                    for k in 0..8usize {
                        let a = mm.max[k].wrapping_sub(mmv.min[k]) & 0x80008000u32 as i32;
                        let b = mmv.max[k].wrapping_sub(mm.min[k]) & 0x80008000u32 as i32;
                        println!(
                            "GATE k={k} a={a:08x} b={b:08x}{}",
                            if a | b != 0 { "  <-- REJECT" } else { "" }
                        );
                    }
                    println!(
                        "GATE z={z:.17} tol={tol:.17}"
                    );
                    println!(
                        "GATE face w6 min={:08x} max={:08x} | w7 min={:08x} max={:08x}",
                        mm.min[6], mm.max[6], mm.min[7], mm.max[7]
                    );
                    println!(
                        "GATE pt   w6 min={:08x} max={:08x} | w7 min={:08x} max={:08x}",
                        mmv.min[6], mmv.max[6], mmv.min[7], mmv.max[7]
                    );
                    println!(
                        "GATE raw z-lanes: deca14={:.17} surD14={:.17} tmin14={:.17} tmax14={:.17} | deca15={:.17} surD15={:.17} tmin15={:.17} tmax15={:.17}",
                        ds.deca()[14],
                        ds.sur_d()[14],
                        tmin[14],
                        tmax[14],
                        ds.deca()[15],
                        ds.sur_d()[15],
                        tmin[15],
                        tmax[15]
                    );
                }
            }
        }
    }
}

/// OCCT anchor (DRAWEXE `vcomputehlr b result -algoType algo 0 0 0 1 -1 1
/// -1 1 2` over `box b 10 20 30`): the visible compound keeps the 9 sharp
/// edges of the three visible faces (projected total 146.9694) and the
/// hidden compound the 3 far edges (48.9898); 9 + 3 = 12 = every box edge.
/// The projected segment endpoints are the analytic corner images
/// x' = (x+y)/sqrt(2), y' = (-x+y+2z)/sqrt(6).
#[test]
fn smoke_box_hlr_end_to_end() {
    let (brep, solid) = smoke_box_solid();
    let (v, _rg1v, _rgnv, outline, h) = smoke_run_hlr(&solid, brep);

    // no outlines on a pure-plane solid.
    let mut outline_edges = Vec::new();
    if !outline.is_null() {
        compound_edges(&outline, &mut outline_edges);
    }
    assert!(outline_edges.is_empty(), "plane solid has no outlines");

    let mut v_edges = Vec::new();
    assert!(!v.is_null(), "v_compound is null");
    compound_edges(&v, &mut v_edges);
    let mut h_edges = Vec::new();
    assert!(!h.is_null(), "h_compound is null");
    compound_edges(&h, &mut h_edges);
    println!("SMOKE v({}): {:?}", v_edges.len(), v_edges);
    println!("SMOKE h({}): {:?}", h_edges.len(), h_edges);

    let total_len = |segs: &[SegSummary]| -> f64 {
        segs
            .iter()
            .map(|s| match s {
                SegSummary::Line { len, .. } => *len,
                _ => f64::NAN,
            })
            .sum()
    };
    let sqrt_2_3 = (2.0f64 / 3.0).sqrt();
    // OCCT: 9 visible edges, mass 146.969; 3 hidden edges, mass 48.9898.
    assert_eq!(v_edges.len(), 9, "visible edges: {:?}", v_edges);
    assert_eq!(h_edges.len(), 3, "hidden edges: {:?}", h_edges);
    let v_len = total_len(&v_edges);
    let h_len = total_len(&h_edges);
    assert!(
        (v_len - 146.96939).abs() < 1e-4,
        "visible projected mass {} vs OCCT 146.969",
        v_len
    );
    assert!(
        (h_len - (60.0 * sqrt_2_3)).abs() < 1e-4,
        "hidden projected mass {} vs OCCT 48.9898",
        h_len
    );

    // the projected segment set: every segment matches the analytic corner
    // image pair (unordered within 1e-6); z = 0 (the projection plane image).
    let img = |x: f64, y: f64, z: f64| {
        glam::DVec3::new(
            (x + y) / (2.0f64).sqrt(),
            (-x + y + 2.0 * z) / (6.0f64).sqrt(),
            0.0,
        )
    };
    let corner_img = |i: usize| match i {
        0 => img(0.0, 0.0, 0.0),
        1 => img(10.0, 0.0, 0.0),
        2 => img(10.0, 20.0, 0.0),
        3 => img(0.0, 20.0, 0.0),
        4 => img(0.0, 0.0, 30.0),
        5 => img(10.0, 0.0, 30.0),
        6 => img(10.0, 20.0, 30.0),
        _ => img(0.0, 20.0, 30.0),
    };
    let mut expected: Vec<(glam::DVec3, glam::DVec3, bool)> = Vec::new();
    // the visible edge set: {01, 34, 45, 56, 67, 12, 26, 15, 03}-pairs per
    // the OCCT face walk: AB EF FG GH BC EH AE BF CG visible; DC AD DH hidden.
    for (i, j, vis) in [
        (0, 1, true),
        (5, 4, true),
        (5, 6, true),
        (6, 7, true),
        (1, 2, true),
        (4, 7, true),
        (0, 4, true),
        (1, 5, true),
        (2, 6, true),
        (2, 3, false),
        (0, 3, false),
        (3, 7, false),
    ] {
        expected.push((corner_img(i), corner_img(j), vis));
    }
    let close = |a: glam::DVec3, b: glam::DVec3| a.distance(b) < 1e-6;
    let seg_match = |segs: &[SegSummary], a: glam::DVec3, b: glam::DVec3| -> bool {
        segs.iter().any(|s| match s {
            SegSummary::Line { p1, p2, .. } => {
                (close(*p1, a) && close(*p2, b)) || (close(*p1, b) && close(*p2, a))
            }
            _ => false,
        })
    };
    for (a, b, vis) in &expected {
        let in_v = seg_match(&v_edges, *a, *b);
        let in_h = seg_match(&h_edges, *a, *b);
        if *vis {
            assert!(in_v && !in_h, "segment {:?}->{:?} must be visible", a, b);
        } else {
            assert!(in_h && !in_v, "segment {:?}->{:?} must be hidden", a, b);
        }
    }
}

/// Stability: three consecutive full runs give identical compounds (the
/// fclass2d seed-order regression guard).
#[test]
fn smoke_box_hlr_stable_over_three_runs() {
    let mut results = Vec::new();
    for _ in 0..3 {
        let (brep, solid) = smoke_box_solid();
        let (v, _rg1v, _rgnv, outline, h) = smoke_run_hlr(&solid, brep);
        let mut vs = Vec::new();
        if !v.is_null() {
            compound_edges(&v, &mut vs);
        }
        let mut os = Vec::new();
        if !outline.is_null() {
            compound_edges(&outline, &mut os);
        }
        let mut hs = Vec::new();
        if !h.is_null() {
            compound_edges(&h, &mut hs);
        }
        results.push((vs, os, hs));
    }
    assert_eq!(results[0], results[1]);
    assert_eq!(results[1], results[2]);
}


/// The cylinder fixture: the OCCT `pcylinder c 10 30` equivalent — a
/// cylindrical side face with a seam, two planar caps, radius 10, height
/// 30.  Built like the box fixture: outward surfaces, pcurves per
/// (edge, face) (the seam carries the closed pair u = 0 / u = 2*pi).
pub(crate) fn smoke_cylinder_solid() -> (rcad_kernel::BRep, rcad_kernel::topods::Shape) {
    use rcad_kernel::geom::{Circle2d, CylindricalSurface, Ellipse2d, Line2d};
    let r = 10.0f64;
    let h = 30.0f64;
    let tau = std::f64::consts::TAU;
    let mut brep = rcad_kernel::BRep::new();
    let mut b = BRepBuilder::new();
    // the seam vertices at (r, 0, z).
    let v_lo = b.add_vertex(&mut brep, glam::DVec3::new(r, 0.0, 0.0), 1e-7);
    let v_hi = b.add_vertex(&mut brep, glam::DVec3::new(r, 0.0, h), 1e-7);
    let oriented = |v: &rcad_kernel::topods::Shape, o: rcad_kernel::topods::Orientation| {
        let mut s = v.clone();
        s.orientation = o;
        s
    };
    // the two circle edges (closed: first == last seam vertex).
    let circ_lo = b.add_edge(
        &mut brep,
        Some(Curve3::Circle(Circle3 {
            center: glam::DVec3::new(0.0, 0.0, 0.0),
            normal: glam::DVec3::new(0.0, 0.0, 1.0),
            x_dir: glam::DVec3::new(1.0, 0.0, 0.0),
            y_dir: glam::DVec3::new(0.0, 1.0, 0.0),
            radius: r,
        })),
        v_lo.clone(),
        oriented(&v_lo, rcad_kernel::topods::Orientation::Reversed),
        [0.0, tau],
    );
    let circ_hi = b.add_edge(
        &mut brep,
        Some(Curve3::Circle(Circle3 {
            center: glam::DVec3::new(0.0, 0.0, h),
            normal: glam::DVec3::new(0.0, 0.0, 1.0),
            x_dir: glam::DVec3::new(1.0, 0.0, 0.0),
            y_dir: glam::DVec3::new(0.0, 1.0, 0.0),
            radius: r,
        })),
        v_hi.clone(),
        oriented(&v_hi, rcad_kernel::topods::Orientation::Reversed),
        [0.0, tau],
    );
    // the seam line (r, 0, 0) -> (r, 0, h).
    let seam = b.add_edge(
        &mut brep,
        Some(Curve3::Line(Line3 {
            origin: glam::DVec3::new(r, 0.0, 0.0),
            direction: glam::DVec3::new(0.0, 0.0, 1.0),
        })),
        v_lo.clone(),
        oriented(&v_hi, rcad_kernel::topods::Orientation::Reversed),
        [0.0, h],
    );
    let f = |e: &rcad_kernel::topods::Shape| oriented(e, rcad_kernel::topods::Orientation::Forward);
    let rv = |e: &rcad_kernel::topods::Shape| oriented(e, rcad_kernel::topods::Orientation::Reversed);
    // the side wire: the uv boundary (u, v) = bottom circle, seam, top
    // circle, seam.
    let side_wire = brep.add_twire(vec![
        f(&circ_lo),
        f(&seam),
        rv(&circ_hi),
        rv(&seam),
    ]);
    let top_wire = brep.add_twire(vec![f(&circ_hi)]);
    let bottom_wire = brep.add_twire(vec![rv(&circ_lo)]);
    let side = brep.add_tface(
        Some(Surface3::Cylinder(CylindricalSurface {
            origin: glam::DVec3::ZERO,
            axis: glam::DVec3::new(0.0, 0.0, 1.0),
            radius: r,
            ref_dir: glam::DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        })),
        side_wire,
        Vec::new(),
        None,
        Some([0.0, tau, 0.0, h]),
        Vec::new(),
        false,
    );
    let top_plane = Plane {
        origin: glam::DVec3::new(0.0, 0.0, h),
        normal: glam::DVec3::new(0.0, 0.0, 1.0),
        u_dir: glam::DVec3::new(1.0, 0.0, 0.0),
        v_dir: glam::DVec3::new(0.0, 1.0, 0.0),
    };
    let top = brep.add_tface(
        Some(Surface3::Plane(top_plane)),
        top_wire,
        Vec::new(),
        None,
        Some([-r, r, -r, r]),
        Vec::new(),
        true,
    );
    let bottom_plane = Plane {
        origin: glam::DVec3::ZERO,
        normal: glam::DVec3::new(0.0, 0.0, -1.0),
        u_dir: glam::DVec3::new(1.0, 0.0, 0.0),
        v_dir: glam::DVec3::new(0.0, -1.0, 0.0),
    };
    let bottom = brep.add_tface(
        Some(Surface3::Plane(bottom_plane)),
        bottom_wire,
        Vec::new(),
        None,
        Some([-r, r, -r, r]),
        Vec::new(),
        true,
    );
    // the pcurves.
    // side face: the circles are the v = 0 / v = h lines, the seam the
    // u = 0 / u = 2*pi closed pair.
    b.add_pcurve(
        &mut brep,
        circ_lo.clone(),
        side.clone(),
        rcad_kernel::geom::Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(1.0, 0.0),
        }),
        0.0,
        tau,
    );
    b.add_pcurve(
        &mut brep,
        circ_hi.clone(),
        side.clone(),
        rcad_kernel::geom::Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, h),
            direction: glam::DVec2::new(1.0, 0.0),
        }),
        0.0,
        tau,
    );
    b.update_edge_pcurve_closed(
        &mut brep,
        seam.clone(),
        rcad_kernel::geom::Curve2d::Line(Line2d {
            origin: glam::DVec2::new(0.0, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        }),
        rcad_kernel::geom::Curve2d::Line(Line2d {
            origin: glam::DVec2::new(tau, 0.0),
            direction: glam::DVec2::new(0.0, 1.0),
        }),
        side.clone(),
        0.0,
        h,
        1e-7,
    );
    // OCCT BRepPrim_Builder::SetPCurve(E, F, L1, L2) (BRepPrim_Builder.cxx
    // L107-118): after the closed pcurve pair UpdateEdge, the seam regularity
    // — myBuilder.Continuity(E, F, F, GeomAbs_CN).  This is the record the
    // HLRBRep_ShapeToHLR::Load reads back as rg = CN (reg1 AND regn true,
    // cxx L129-131).
    b.continuity(
        &mut brep,
        &seam,
        &side,
        &side,
        rcad_kernel::topods::GeomAbsShape::CN,
    );
    // caps: the circles project to full circles.
    b.add_pcurve(
        &mut brep,
        circ_hi.clone(),
        top.clone(),
        rcad_kernel::geom::Curve2d::Circle(Circle2d {
            center: glam::DVec2::new(0.0, 0.0),
            x_dir: glam::DVec2::new(1.0, 0.0),
            y_dir: glam::DVec2::new(0.0, 1.0),
            radius: r,
        }),
        0.0,
        tau,
    );
    b.add_pcurve(
        &mut brep,
        circ_lo.clone(),
        bottom.clone(),
        rcad_kernel::geom::Curve2d::Circle(Circle2d {
            center: glam::DVec2::new(0.0, 0.0),
            x_dir: glam::DVec2::new(1.0, 0.0),
            y_dir: glam::DVec2::new(0.0, -1.0),
            radius: r,
        }),
        0.0,
        tau,
    );
    let _ = Ellipse2d {
        center: glam::DVec2::ZERO,
        major_dir: glam::DVec2::X,
        minor_dir: glam::DVec2::Y,
        major_radius: r,
        minor_radius: r,
    };
    let shell = brep.add_tshell(vec![side, top, bottom]);
    let solid = brep.add_tsolid(vec![shell]);
    (brep, solid)
}

/// The projected ellipse-arc length over an angle span (Simpson).
pub(crate) fn ellipse_arc_len(a: f64, b: f64, u1: f64, u2: f64) -> f64 {
    let f = |t: f64| (a * a * t.sin() * t.sin() + b * b * t.cos() * t.cos()).sqrt();
    let n = 64;
    let hh = (u2 - u1) / n as f64;
    let mut s = f(u1) + f(u2);
    for k in 1..n {
        let w = if k % 2 == 1 { 4.0 } else { 2.0 };
        s += w * f(u1 + k as f64 * hh);
    }
    s * hh / 3.0
}

/// OCCT anchor (DRAWEXE `vcomputehlr c ... -algoType algo -showHiddenEdges`
/// over `pcylinder c 10 30`, dir (1,-1,1), up (-1,1,2)): the projected
/// top/bottom circles are the ellipses a = 10, b = 10/sqrt(3) centered at
/// (0, 24.49490) and (0, 0); OCCT's full reference is OutLineVCompound = the
/// two silhouette lines (projected length 24.49490 each, at x' = +/-10),
/// VCompound = the full top ellipse plus the near half of the bottom ellipse
/// (total 75.671, 5 arcs), HCompound = the far half of the bottom ellipse
/// (25.2237, 1 arc).
///
/// Pipeline gaps recorded against that reference (the anchors below pin the
/// currently matching subset):
/// - FIXED (agent T, the OCCT DomainIntersection rule in g_inter.rs): Contap
///   finds BOTH contour lines (the phi = 45 deg one at x' = +10 and the
///   phi = 225 deg one at x' = -10) and both cap circles split at the two
///   tangency points (ellipse params 0 and pi, plus the seam at 7*pi/4).
/// - FIXED (the hider layer): the bottom far half-ellipse (params [0, pi],
///   25.2237) is hidden by the cylindrical side face through the OCCT
///   Compare/Classify branch (no boundary interferences: the endpoints are
///   the silhouette tangency vertices and the seam is z-rejected).  This
///   needed two 1:1 corrections: HLRAlgo_Projector::Shoot transforms the
///   line by myInvTrsf where gp_Lin::Transform moves the gp_Ax1 position
///   AND its direction (gp_Lin.hxx L178 -> gp_Ax1.hxx L201-205), and
///   HLRBRep_Curve::First/LastParameter are Parameter2d of the raw curve
///   parameters (HLRBRep_Curve.lxx L66-76).
/// - CLOSED (gap 4, the seam regularity): the seam edge is a
///   BRepPrim_Builder::SetPCurve(E, F, L1, L2) closed pair, so OCCT records
///   Continuity(E, F, F, GeomAbs_CN) on it (BRepPrim_Builder.cxx L117);
///   HLRBRep_ShapeToHLR::Load reads rg = CN with the GeomAbs_Shape ranking
///   C0 < G1 < C1 < G2 < C2 < C3 < CN (GeomAbs_Shape.hxx L47-56), so
///   reg1 = (CN >= G1) AND regn = (CN >= G2) (cxx L129-131).  The
///   HLRBRep_HLRToShape typ filters (HLRBRep_HLRToShape.cxx L186-196) then
///   route the Rg1+RgN seam OUT of VCompound (typ 5) and Rg1LineVCompound
///   (typ 3: Rg1Line && !RgNLine) into RgNLineVCompound (typ 4: RgNLine).
/// - CLOSED (gap 5): each shared edge is drawn exactly once — the OCCT
///   Used-flag dedupe of the HLRBRep_HLRToShape face walk
///   (InternalCompound cxx L75-87 reset + DrawFace cxx L168-225 +
///   DrawEdge, aligned in brep/hlr_to_shape.rs and the EMaskUsed flag of
///   brep/edge_data.rs): the first face that reaches a circle edge draws
///   it and sets Used, the adjacent-face occurrence skips; the seam
///   closed pair (forward + reversed in the side wire) draws once.
#[test]
fn smoke_cylinder_hlr_end_to_end() {
    let (brep, solid) = smoke_cylinder_solid();
    let (v, rg1v, rgnv, outline, h) = smoke_run_hlr(&solid, brep);

    let b_minor = 10.0f64 / 3.0f64.sqrt();
    let per = ellipse_arc_len(10.0, b_minor, 0.0, std::f64::consts::TAU);

    // the outline: BOTH silhouette lines (OCCT-exact, the OutLineVCompound
    // reference) — projected length 24.49490 each, at x' = +10 (phi = 45 deg)
    // and x' = -10 (phi = 225 deg), y' from 0 to 24.49490.
    let mut o_edges = Vec::new();
    assert!(!outline.is_null(), "outline compound is null");
    compound_edges(&outline, &mut o_edges);
    assert!(!o_edges.is_empty(), "no outline edges");
    let mut plus_outline = false;
    let mut minus_outline = false;
    for o in &o_edges {
        match o {
            SegSummary::Line { len, p1, p2 } => {
                assert!((len - 30.0 * 2.0 / 6.0f64.sqrt()).abs() < 1e-9);
                assert!((p1.z).abs() < 1e-9 && (p2.z).abs() < 1e-9);
                // one end at the bottom plane image, the other at the top.
                let (lo, hi) = if p1.y < p2.y { (*p1, *p2) } else { (*p2, *p1) };
                assert!((lo.y).abs() < 1e-9, "silhouette start {:?}", lo);
                assert!((hi.y - 60.0 / 6.0f64.sqrt()).abs() < 1e-9);
                assert!((lo.x - hi.x).abs() < 1e-9);
                if (lo.x - 10.0).abs() < 1e-9 {
                    plus_outline = true;
                } else if (lo.x + 10.0).abs() < 1e-9 {
                    minus_outline = true;
                } else {
                    panic!("silhouette at unexpected x' = {:?}", lo.x);
                }
            }
            _ => panic!("outline edge is not a line: {:?}", o),
        }
    }
    assert!(
        plus_outline && minus_outline,
        "both phi = 45/225 deg silhouettes are required (OCCT truth): {:?}",
        o_edges
    );

    // the visible sharp edges: ellipse arcs of a = 10, b = 10/sqrt(3) — the
    // projected images of the cap circles ONLY.  The seam edge (reg1 AND
    // regn from the CN regularity) is excluded from VCompound by the typ 5
    // filter (!Rg1Line || OutLine) of HLRBRep_HLRToShape.cxx L196.
    let mut v_edges = Vec::new();
    assert!(!v.is_null(), "v_compound is null");
    compound_edges(&v, &mut v_edges);
    assert!(!v_edges.is_empty());
    let mut v_ellipse_arcs = 0;
    let mut v_seam = 0;
    for e in &v_edges {
        match e {
            SegSummary::Ellipse { major, minor, .. } => {
                assert!((major - 10.0).abs() < 1e-9 && (minor - b_minor).abs() < 1e-9);
                v_ellipse_arcs += 1;
            }
            SegSummary::Line { .. } => {
                v_seam += 1;
            }
            _ => panic!("unexpected visible edge: {:?}", e),
        }
    }
    assert_eq!(
        v_ellipse_arcs, 5,
        "OCCT VCompound = the full top ellipse (3 arcs) + the bottom near half (2 arcs): {:?}",
        v_edges
    );
    assert_eq!(v_seam, 0, "the seam must leave VCompound (typ 5): {:?}", v_edges);

    // the seam edge (rg1_line AND rg_n_line, HLRBRep_HLRToShape.cxx L188 vs
    // L192) is NOT in Rg1LineVCompound (typ 3 needs !RgNLine) — the Rg1
    // compound stays null.
    let mut rg1_edges = Vec::new();
    if !rg1v.is_null() {
        compound_edges(&rg1v, &mut rg1_edges);
    }
    assert!(rg1_edges.is_empty(), "no Rg1Line edge expected: {:?}", rg1_edges);

    // gap 4 pin: the seam lands in RgNLineVCompound (typ 4: RgNLine &&
    // !OutLine, cxx L190-192) — the projected (r, 0, z) line x' = r/sqrt(2),
    // y' from -r/sqrt(6) to (-r + 2h)/sqrt(6), projected length h*2/sqrt(6).
    assert!(!rgnv.is_null(), "rgn_line_v_compound is null");
    let mut rgn_edges = Vec::new();
    compound_edges(&rgnv, &mut rgn_edges);
    assert_eq!(rgn_edges.len(), 1, "exactly one seam edge drawn: {:?}", rgn_edges);
    match &rgn_edges[0] {
        SegSummary::Line { len, p1, p2 } => {
            assert!((len - 60.0 / 6.0f64.sqrt()).abs() < 1e-9);
            assert!((p1.z).abs() < 1e-9 && (p2.z).abs() < 1e-9);
            assert!((p1.x - p2.x).abs() < 1e-9);
            assert!((p1.x - 10.0 / 2.0f64.sqrt()).abs() < 1e-9);
            let (lo, hi) = if p1.y < p2.y { (*p1, *p2) } else { (*p2, *p1) };
            assert!((lo.y + 10.0 / 6.0f64.sqrt()).abs() < 1e-9, "seam lo {:?}", lo);
            assert!((hi.y - 50.0 / 6.0f64.sqrt()).abs() < 1e-9, "seam hi {:?}", hi);
        }
        _ => panic!("seam edge is not a line: {:?}", rgn_edges[0]),
    }

    // the drawn arc pieces: the TOP circle keeps [0, pi] + [pi, 7*pi/4] +
    // [7*pi/4, 2*pi] (the full ellipse), the BOTTOM circle keeps only the
    // near halves [pi, 7*pi/4] + [7*pi/4, 2*pi] — the far half [0, pi] went
    // to the hidden compound (the two tangency splits at the phi = 45/225
    // deg silhouettes and the seam split at 7*pi/4).
    let tau = std::f64::consts::TAU;
    let half = std::f64::consts::PI;
    let split = 7.0 * std::f64::consts::FRAC_PI_4;
    let mut ranges = edges_with_ranges(&v);
    ranges.retain(|(u1, u2)| *u2 - *u1 < tau); // drop the seam line
    ranges.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let expect = [
        (0.0, half),
        (half, split),
        (half, split),
        (split, tau),
        (split, tau),
    ];
    for (i, (u1, u2)) in ranges.iter().enumerate() {
        assert!(
            (u1 - expect[i].0).abs() < 1e-9 && (u2 - expect[i].1).abs() < 1e-9,
            "arc piece {}: ({},{}) vs {:?}",
            i,
            u1,
            u2,
            expect[i]
        );
    }

    // gap 5 pin: no arc segment is drawn twice — the ellipse arcs carry
    // pairwise distinct (unordered) endpoint pairs (the OCCT Used-flag
    // dedupe; a broken dedupe would repeat the first face's arcs verbatim).
    let arc_ends: Vec<(glam::DVec3, glam::DVec3)> = v_edges
        .iter()
        .filter_map(|e| match e {
            SegSummary::Ellipse { p1, p2, .. } => Some((*p1, *p2)),
            _ => None,
        })
        .collect();
    for i in 0..arc_ends.len() {
        for j in (i + 1)..arc_ends.len() {
            let (a1, a2) = arc_ends[i];
            let (b1, b2) = arc_ends[j];
            let same = |p: glam::DVec3, q: glam::DVec3| (p - q).length() < 1e-9;
            let dup = (same(a1, b1) && same(a2, b2)) || (same(a1, b2) && same(a2, b1));
            assert!(
                !dup,
                "shared edge drawn twice: arc[{i}] {:?} vs arc[{j}] {:?}",
                arc_ends[i],
                arc_ends[j]
            );
        }
    }

    // the arc total: the OCCT reference draws the top full ellipse + the
    // bottom near half (75.671 = 1.5 perimeters; the far half 25.2237 is in
    // the hidden compound).
    let mut v_arc = 0.0;
    for (u1, u2) in &ranges {
        v_arc += ellipse_arc_len(10.0, b_minor, *u1, *u2);
    }
    assert!((v_arc - 1.5 * per).abs() < 1e-2, "visible arc {}", v_arc);

    // the hidden compound: OCCT holds the far half of the bottom ellipse —
    // exactly one arc over [0, pi] (the projected 25.2237), endpoints at
    // the silhouette tangency images (x' = +/-10, y' = 0).
    assert!(!h.is_null(), "the hidden compound is null");
    let mut h_edges = Vec::new();
    compound_edges(&h, &mut h_edges);
    assert_eq!(h_edges.len(), 1, "one hidden arc: {:?}", h_edges);
    match &h_edges[0] {
        SegSummary::Ellipse { p1, p2, major, minor, .. } => {
            assert!((major - 10.0).abs() < 1e-9 && (minor - b_minor).abs() < 1e-9);
            let (lo, hi) = if p1.x < p2.x { (*p1, *p2) } else { (*p2, *p1) };
            assert!((lo.x + 10.0).abs() < 1e-9 && lo.y.abs() < 1e-9, "h lo {:?}", lo);
            assert!((hi.x - 10.0).abs() < 1e-9 && hi.y.abs() < 1e-9, "h hi {:?}", hi);
        }
        _ => panic!("hidden edge is not an ellipse arc: {:?}", h_edges[0]),
    }
    let mut h_ranges = edges_with_ranges(&h);
    assert_eq!(h_ranges.len(), 1);
    let (h1, h2) = h_ranges.remove(0);
    assert!(
        (h1 - 0.0).abs() < 1e-9 && (h2 - half).abs() < 1e-9,
        "hidden arc range ({},{})",
        h1,
        h2
    );
    assert!(
        (ellipse_arc_len(10.0, b_minor, h1, h2) - 25.2237).abs() < 1e-3,
        "hidden arc length {}",
        ellipse_arc_len(10.0, b_minor, h1, h2)
    );
}

/// The (range) pairs of the ellipse edges of a result compound.
pub(crate) fn edges_with_ranges(s: &rcad_kernel::topods::Shape) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    fn walk(s: &rcad_kernel::topods::Shape, out: &mut Vec<(f64, f64)>) {
        match &*s.data {
            TShape::Compound(c) => {
                for ch in c {
                    walk(ch, out);
                }
            }
            TShape::Edge(ed) => out.push((ed.range[0], ed.range[1])),
            _ => {}
        }
    }
    walk(s, &mut out);
    out
}

/// Stability: three consecutive full cylinder runs give identical compounds.
#[test]
fn smoke_cylinder_hlr_stable_over_three_runs() {
    let mut results = Vec::new();
    for _ in 0..3 {
        let (brep, solid) = smoke_cylinder_solid();
        let (v, _rg1v, _rgnv, outline, h) = smoke_run_hlr(&solid, brep);
        let mut vs = Vec::new();
        if !v.is_null() {
            compound_edges(&v, &mut vs);
        }
        let mut os = Vec::new();
        if !outline.is_null() {
            compound_edges(&outline, &mut os);
        }
        let mut hs = Vec::new();
        if !h.is_null() {
            compound_edges(&h, &mut hs);
        }
        results.push((vs, os, hs));
    }
    // NaN-tolerant compare (the Ellipse summaries carry a non-computed len).
    assert!(smoke_results_equal(&results[0], &results[1]));
    assert!(smoke_results_equal(&results[1], &results[2]));
}

/// NaN-tolerant equality of the smoke summaries (the Ellipse `len` is not
/// computed).
fn smoke_results_equal(
    a: &(
        Vec<SegSummary>,
        Vec<SegSummary>,
        Vec<SegSummary>,
    ),
    b: &(
        Vec<SegSummary>,
        Vec<SegSummary>,
        Vec<SegSummary>,
    ),
) -> bool {
    let seg_eq = |x: &SegSummary, y: &SegSummary| match (x, y) {
        (
            SegSummary::Line {
                len: l1,
                p1: a1,
                p2: b1,
            },
            SegSummary::Line {
                len: l2,
                p1: a2,
                p2: b2,
            },
        ) => (l1 - l2).abs() < 1e-9 && a1.distance(*a2) < 1e-9 && b1.distance(*b2) < 1e-9,
        _ => format!("{:?}", x) == format!("{:?}", y),
    };
    let vec_eq = |x: &Vec<SegSummary>, y: &Vec<SegSummary>| {
        x.len() == y.len() && x.iter().zip(y).all(|(p, q)| seg_eq(p, q))
    };
    let (av, ao, ah) = a;
    let (bv, bo, bh) = b;
    vec_eq(av, bv) && vec_eq(ao, bo) && vec_eq(ah, bh)
}



// ---- The cylinder Contap minimal repro (agent T) ----
// OCCT anchor: Contap_Contour::Perform over the `pcylinder c 10 30` side
// face with the smoke view direction (1,-1,1)/sqrt(3) — Contap_ContAna
// Perform(gp_Cylinder, gp_Dir) (Contap_ContAna.cxx L116-138) yields TWO
// silhouette generatrices (pt1 = loc + R*normale at u = 45 deg, pt2 = loc -
// R*normale at u = 225 deg); the restriction search finds both tangency
// points on each cap circle (u = pi/4 and u = 5*pi/4), so PerformAna
// (Contap_Contour.cxx L2156-2389) ends with one split Lin per generatrix,
// each carrying the two cap vertices.
#[test]
fn contour_cylinder_side_face_two_silhouettes() {
    use crate::hlr::contap::contour::Contour;
    use crate::hlr::contap::i_type::IType;
    use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;
    use std::sync::Arc;

    let (brep, solid) = smoke_cylinder_solid();

    // The direction OutLiner::fill computes: vecz transformed by the
    // inverted projection — the normalized (1,-1,1)/sqrt(3).
    let proj = SmokeProjector::from_ax2(&SmokeAx2::new(
        glam::DVec3::ZERO,
        glam::DVec3::new(1.0, -1.0, 1.0),
        glam::DVec3::new(1.0, 1.0, 0.0),
    ));
    let mut tr = *proj.transformation();
    tr.invert();
    let vecz = tr.transform_vec(glam::DVec3::new(0.0, 0.0, 1.0));
    let vecz = vecz.normalize();
    assert!((vecz.x - 1.0 / 3.0f64.sqrt()).abs() < 1e-12);
    assert!((vecz.y + 1.0 / 3.0f64.sqrt()).abs() < 1e-12);

    // The cylindrical side face of the fixture.
    let mut faces = Vec::new();
    fn walk(s: &rcad_kernel::topods::Shape, out: &mut Vec<rcad_kernel::topods::Shape>) {
        use rcad_kernel::topods::TShape;
        match &*s.data {
            TShape::Compound(c) => {
                for ch in c {
                    walk(ch, out);
                }
            }
            TShape::Solid(sl) => {
                for sh in &sl.shells {
                    walk(sh, out);
                }
            }
            TShape::Shell(sh) => {
                for f in &sh.faces {
                    walk(f, out);
                }
            }
            TShape::Face(_) => out.push(s.clone()),
            _ => {}
        }
    }
    walk(&solid, &mut faces);
    let tool_brep = Arc::new(brep);
    let mut silhouettes = Vec::new();
    for f in &faces {
        let mut brt = BRepTopAdaptorTool::new_face(tool_brep.clone(), f, 1e-7);
        let surface = brt.get_surface().unwrap().clone();
        if !matches!(
            surface.get_type(),
            crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
        ) {
            continue;
        }
        let mut fo = Contour::with_direction(vecz);
        fo.perform(&surface, brt.get_topol_tool());
        assert!(fo.is_done());
        for i in 1..=fo.nb_lines() {
            let l = fo.line(i);
            assert_eq!(l.type_contour(), IType::Lin);
            assert_eq!(l.nb_vertex(), 2);
            let lin = l.line();
            assert!((lin.direction.z - 1.0).abs() < 1e-12);
            silhouettes.push(lin.origin);
        }
    }

    // The two generatrices at u = 45 deg (x' = +10) and u = 225 deg
    // (x' = -10): origins (+-10/sqrt(2), +-10/sqrt(2), 0) on the base circle.
    assert_eq!(silhouettes.len(), 2, "silhouettes={:?}", silhouettes);
    let r2 = 10.0f64 / 2.0f64.sqrt();
    let mut matched_pos = false;
    let mut matched_neg = false;
    for o in &silhouettes {
        if (o.z).abs() < 1e-9 && (o.x - r2).abs() < 1e-9 && (o.y - r2).abs() < 1e-9 {
            matched_pos = true;
        }
        if (o.z).abs() < 1e-9 && (o.x + r2).abs() < 1e-9 && (o.y + r2).abs() < 1e-9 {
            matched_neg = true;
        }
    }
    assert!(
        matched_pos && matched_neg,
        "the phi=225 deg silhouette is missing: {:?}",
        silhouettes
    );

    // Both lines run the full height v in [0, 30] (the vertex parameters).
    // The cap-circle split parameters are checked through the smoke
    // end-to-end test below.
}
