//! Regression guard for the Geom_BSplineSurface::UIso / VIso arms
//! (Geom_BSplineSurface_1.cxx L598-630 / L775-807): the iso poles are
//! BSplSLib::Iso over the flat knot vector of the iso direction, and the
//! produced curve carries the opposite direction's knot vector.

use super::*;

/// A bilinear patch: `myUDeg = myVDeg = 1`, both knot vectors the flat
/// form of [0, 0, 1, 1].
fn bilinear_patch(weights: Vec<Vec<f64>>) -> BSplineSurface {
    BSplineSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 2.0)],
            vec![DVec3::new(3.0, 0.0, 1.0), DVec3::new(3.0, 1.0, -1.0)],
        ],
        weights,
        is_periodic_u: false,
        is_periodic_v: false,
    }
}

/// Expected values derived from the OCCT formula alone, for the uniform
/// 2.0 weight grid of [`bilinear_patch`] at V = 0.25:
///
/// - The static `Rational` (Geom_BSplineSurface.cxx L110-138) compares
///   weights adjacent in the pole grid; every pair is equal, so
///   `myURational || myVRational` is false and the `else` arm of
///   Geom_BSplineSurface::VIso (L792-800) runs `BSplSLib::Iso` with
///   `BSplSLib::NoWeights()`.
/// - A non-rational `Iso` writes `1.` into every `CWeights` entry
///   (BSplSLib.cxx L1732-1738): the uniform 2.0 must not reach the curve.
/// - The iso poles are the degree-1 corner cutting over the V window.
///   `LocateParameter(1, [0, 0, 1, 1], nullptr, 0.25, false, Index, u)`
///   gives `Index = 2` (1-based; the largest knot index with
///   `Knots(Index) <= 0.25`), so `BuildKnots` copies
///   `[Knots(Index), Knots(Index + 1)] = [0, 1]`
///   (BSplCLib.cxx L1570-1575), and `BSplCLib::Eval` combines the two V
///   poles with `X = (locknots[1] - U) / (locknots[1] - locknots[0])
///   = 0.75`, `Y = 1 - X` (BSplCLib.cxx L865-1006):
///   `CPoles(i) = 0.75 * Poles(i, 1) + 0.25 * Poles(i, 2)`.
/// - The curve carries the V data (L787): degree `myVDeg = 1`, knots
///   `myVFlatKnots`.
#[test]
fn viso_takes_the_non_rational_arm_for_a_uniform_weight_grid() {
    let surf = bilinear_patch(vec![vec![2.0, 2.0], vec![2.0, 2.0]]);
    let iso = bspline_surface_viso_full(&surf, 0.25);
    assert_eq!(iso.degree, 1);
    assert_eq!(iso.knots, vec![0.0, 0.0, 1.0, 1.0]);
    assert_eq!(iso.weights, vec![1.0, 1.0]);
    let expected = [DVec3::new(0.0, 0.25, 0.5), DVec3::new(3.0, 0.25, 0.5)];
    for (pole, want) in iso.control_points.iter().zip(expected.iter()) {
        assert!((*pole - *want).length() < 1e-15, "{pole:?} != {want:?}");
    }
}

/// The rational arm of the same patch: making one weight differ along the
/// U index sets `myVRational` (Geom_BSplineSurface.cxx L110-124), so
/// `BSplSLib::Iso` runs at dim = 4 with `Weights()` and returns the
/// evaluated homogeneous weights (BSplSLib.cxx L1705-1723) instead of the
/// non-rational 1.0.
#[test]
fn viso_takes_the_rational_arm_when_the_weights_differ() {
    let surf = bilinear_patch(vec![vec![1.0, 1.0], vec![1.0, 3.0]]);
    let iso = bspline_surface_viso_full(&surf, 0.25);
    assert!((iso.weights[0] - 1.0).abs() < 1e-15);
    assert!((iso.weights[1] - 1.5).abs() < 1e-15);
}
