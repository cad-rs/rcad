//! classify_in_restriction: the CSLib_Class2d ray-cast over a closed wire
//! polygon (one full-circle pcurve) must detect the boundary crossing even
//! when the query Y coincides with a sampled vertex Y (CSLib_Class2d.cxx
//! L235-268: the closing segment participates in the crossing count).

use rcad_algo::geomalgo::int_patch::so_on_bounds::{classify_in_restriction, BoundaryArc};
use rcad_kernel::geom::{Circle2d, Curve2d, Line2d};
use rcad_kernel::topods::State;
use glam::DVec2;

fn circle_segs(r: f64) -> Vec<BoundaryArc> {
    vec![BoundaryArc {
        arc: Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: r,
        }),
        first: 0.0,
        last: std::f64::consts::TAU,
        vtx_params: vec![],
    }]
}

#[test]
fn circle_polygon_interior_not_out() {
    let segs = circle_segs(4.0);
    let uv = [-4.0f64, 4.0, -4.0, 4.0];
    let tol = 3.5e-8;
    // Interior points (one at the exact last-sample Y of the polygon).
    assert_eq!(
        classify_in_restriction(2.0, -0.0, uv, &segs, tol, false, false),
        State::In
    );
    assert_eq!(
        classify_in_restriction(-2.0, 0.0, uv, &segs, tol, false, false),
        State::In
    );
    // Outside.
    assert_eq!(
        classify_in_restriction(5.0, 0.0, uv, &segs, tol, false, false),
        State::Out
    );
}

#[test]
fn two_seam_cylinder_polygon_interior() {
    // The closed cylinder lateral UV rectangle: bottom/top generatrix lines at
    // v = -20/0 spanning u in [3pi/2, 7pi/2] plus the doubled seam at u=3pi/2
    // and u=7pi/2 (the CurveOnClosedSurface pair).
    let two_pi = std::f64::consts::TAU;
    let mk_line = |o: (f64, f64), d: (f64, f64), len: f64| BoundaryArc {
        arc: Curve2d::Line(Line2d {
            origin: DVec2::new(o.0, o.1),
            direction: DVec2::new(d.0, d.1),
        }),
        first: 0.0,
        last: len,
        vtx_params: vec![],
    };
    let mut segs = vec![
        // Seam, FORWARD instance (west side of the UV region).
        mk_line((4.712389, 0.0), (0.0, -1.0), 20.0),
        // Seam, REVERSED instance shifted one period east (east side).
        mk_line((4.712389 + two_pi, 0.0), (0.0, -1.0), 20.0),
        // Bottom/top generating lines spanning the full period.
        mk_line((4.712389, 0.0), (1.0, 0.0), two_pi),
        mk_line((4.712389, -20.0), (1.0, 0.0), two_pi),
    ];
    let uv = [0.0f64, std::f64::consts::TAU, -20.0, 0.0];
    let state = classify_in_restriction(
        5.5,
        -10.0,
        uv,
        &segs,
        3.5e-8,
        true,
        false,
    );
    assert_eq!(state, State::In);
}
