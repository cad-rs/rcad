// Rule 5 file-size split: the per-kind evaluation impls moved to child
// modules; the evaluation traits and this module's `pub(crate)` basis
// functions keep their `crate::geom::eval::` paths through the re-exports.
#[path = "eval_basis.rs"]
mod eval_basis;
#[path = "eval_adaptor.rs"]
mod eval_adaptor;
#[path = "eval_conics.rs"]
mod eval_conics;
#[path = "eval_bspline.rs"]
mod eval_bspline;
#[path = "eval_surfaces.rs"]
mod eval_surfaces;
#[path = "eval_offset.rs"]
mod eval_offset;
#[path = "eval_enum_impls.rs"]
mod eval_enum_impls;

pub(crate) use eval_adaptor::*;
pub(crate) use eval_basis::*;

#[cfg(test)]
mod derivative_tests {
    use super::*;
    use crate::geom::*;

    /// A rational BSpline patch: degree 2 in U, degree 1 in V, with a weight
    /// grid that makes both directions rational.
    fn rational_bspline_patch() -> BSplineSurface {
        BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 0.0, 1.0)],
                vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 1.0)],
                vec![DVec3::new(3.0, 0.0, 0.0), DVec3::new(3.0, 0.0, 1.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![3.0, 3.0], vec![1.0, 1.0]],
            is_periodic_u: false,
            is_periodic_v: false,
        }
    }

    /// The non-rational counterpart of [`rational_bspline_patch`] (all weights
    /// equal, so `PrepareEval` keeps the `dim = 3` branch).
    fn plain_bspline_patch() -> BSplineSurface {
        BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 1.0, 0.0)],
                vec![DVec3::new(2.0, 0.0, 0.0), DVec3::new(2.0, 1.0, 0.0)],
            ],
            weights: vec![vec![1.0, 1.0], vec![1.0, 1.0]],
            is_periodic_u: false,
            is_periodic_v: false,
        }
    }

    /// Regression for the derivative-buffer defect: the buffers of
    /// `bspline_surface_dn` were sized for a single pole while the homogeneous
    /// evaluation writes `n + 1` poles, so any first-order `dn` on a BSpline
    /// surface indexed out of bounds (panic).
    #[test]
    fn bspline_surface_dn_first_order_no_panic() {
        let s = Surface3::BSpline(rational_bspline_patch());
        let d = s.dn(1.0, 0.5, 1, 0);
        assert!(d.x.is_finite() && d.y.is_finite() && d.z.is_finite(), "d={d:?}");
        assert!(d.length() > 0.0, "d={d:?}");
        let d = s.dn(1.0, 0.5, 0, 1);
        assert!(d.x.is_finite() && d.y.is_finite() && d.z.is_finite(), "d={d:?}");
    }

    /// The non-rational arm of the same defect (dim = 3 buffers).
    #[test]
    fn bspline_surface_dn_plain_first_order_no_panic() {
        let s = Surface3::BSpline(plain_bspline_patch());
        let d = s.dn(0.5, 0.5, 1, 0);
        assert!(d.x.is_finite() && d.y.is_finite() && d.z.is_finite(), "d={d:?}");
        assert!((d - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-12, "d={d:?}");
        let d = s.dn(0.5, 0.5, 0, 1);
        assert!((d - DVec3::new(0.0, 1.0, 0.0)).length() < 1e-12, "d={d:?}");
    }

    /// Second order terms (including the mixed one) must be finite for a
    /// rational patch.
    #[test]
    fn bspline_surface_dn_second_order_finite() {
        let s = Surface3::BSpline(rational_bspline_patch());
        for (nu, nv) in [(2, 0), (1, 1), (0, 1)] {
            let d = s.dn(1.0, 0.5, nu, nv);
            assert!(d.x.is_finite() && d.y.is_finite() && d.z.is_finite(), "d={d:?}");
        }
    }

    /// OCCT Geom_RectangularTrimmedSurface::EvalDN
    /// (Geom_RectangularTrimmedSurface.cxx L419-429): the range guard, then
    /// `basisSurf->EvalDN(U, V, Nu, Nv)` forwarded unchanged.  On a plane
    /// basis the forwarded value is the `ElSLib::PlaneDN` result (ElSLib.cxx
    /// L169-180): `Pos.XDirection()` for (1, 0), `Pos.YDirection()` for (0, 1)
    /// and zero for every other order.  Both expected vectors are read off the
    /// basis plane payload, so the assertion covers the delegation rather than
    /// the gp_Ax3 construction.
    #[test]
    fn trimmed_surface_dn_delegates_to_basis() {
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let trimmed =
            Surface3::Trimmed(TrimmedSurface::new(Surface3::Plane(plane), 0.5, 1.5, 0.5, 1.5));
        assert_eq!(trimmed.dn(1.0, 0.75, 1, 0), plane.u_dir);
        assert_eq!(trimmed.dn(1.0, 0.75, 0, 1), plane.v_dir);
        assert_eq!(trimmed.dn(1.0, 0.75, 2, 0), DVec3::ZERO);
        assert_eq!(trimmed.dn(1.0, 0.75, 1, 1), DVec3::ZERO);
    }

    /// The rational quarter circle in the XY plane as a degree-2 Bezier, with
    /// the arc extruded along Z (degree 1 in V): `S(u, v) = (cos(theta),
    /// sin(theta), v)` with `theta = u * pi / 2`.  All analytic derivatives
    /// below follow from that closed form, so this checks the rational
    /// machinery of `BSplSLib::DN` / `BSplSLib::D0` against exact values.
    ///
    /// The U weights are `{1, sqrt(2)/2, 1}` (the classic NURBS circle) and
    /// each weight is repeated over V, so the weight grid is
    /// `[w_i, w_i]` — V-constant columns.
    fn rational_arc_patch_control_points() -> (Vec<Vec<DVec3>>, Vec<Vec<f64>>) {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let pts = vec![
            vec![DVec3::new(1.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 1.0)],
            vec![DVec3::new(1.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0)],
            vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(0.0, 1.0, 1.0)],
        ];
        let weights = vec![vec![1.0, 1.0], vec![s, s], vec![1.0, 1.0]];
        (pts, weights)
    }

    fn rational_arc_bspline_patch() -> BSplineSurface {
        let (control_points, weights) = rational_arc_patch_control_points();
        BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points,
            weights,
            is_periodic_u: false,
            is_periodic_v: false,
        }
    }

    fn rational_arc_bezier_patch() -> BezierSurface {
        let (control_points, weights) = rational_arc_patch_control_points();
        BezierSurface { control_points, weights }
    }

    /// The exact arc values of the rational quadratic Bezier `N(u)/W(u)` with
    /// poles `{(1,0), (1,1), (0,1)}` and weights `{1, sqrt(2)/2, 1}` (the unit
    /// quarter circle), derived HERE from the Bernstein polynomials and the
    /// quotient rule `f = N/W`, `f' = (N' - f W') / W`,
    /// `f'' = (N'' - 2 f' W' - f W'') / W`.
    ///
    /// This is an independent oracle: it does not go through `PrepareEval` /
    /// `BSplCLib::Bohm` / `BSplSLib::RationalDerivative`.  NOTE the arc is NOT
    /// parameterised as `theta = u * pi / 2` — the rational quadratic has its
    /// own (non-uniform) parameterisation, e.g. `dP/du(0) = (0, sqrt(2), 0)`
    /// rather than `(0, pi/2, 0)`.
    fn rational_arc_exact(u: f64) -> (DVec3, DVec3, DVec3) {
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let b = [(1.0 - u) * (1.0 - u), 2.0 * u * (1.0 - u), u * u];
        let db = [-2.0 * (1.0 - u), 2.0 - 4.0 * u, 2.0 * u];
        let ddb = [2.0, -4.0, 2.0];
        let pts = [
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let w = [1.0, s, 1.0];
        let mut n = DVec3::ZERO;
        let mut dn = DVec3::ZERO;
        let mut ddn = DVec3::ZERO;
        let (mut wn, mut dwn, mut ddwn) = (0.0f64, 0.0f64, 0.0f64);
        for i in 0..3 {
            n += w[i] * b[i] * pts[i];
            dn += w[i] * db[i] * pts[i];
            ddn += w[i] * ddb[i] * pts[i];
            wn += w[i] * b[i];
            dwn += w[i] * db[i];
            ddwn += w[i] * ddb[i];
        }
        let f = n / wn;
        let df = (dn - f * dwn) / wn;
        let ddf = (ddn - 2.0 * df * dwn - f * ddwn) / wn;
        (f, df, ddf)
    }

    #[test]
    fn rational_bspline_surface_exact_derivatives() {
        let s = Surface3::BSpline(rational_arc_bspline_patch());
        for u in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let (f, df, _ddf) = rational_arc_exact(u);
            let (p, du, dv) = s.derivatives(u, 0.5);
            assert!((p - (f + DVec3::Z * 0.5)).length() < 1e-12, "u={u} p={p:?}");
            assert!((du - df).length() < 1e-12, "u={u} du={du:?} want={df:?}");
            assert!((dv - DVec3::Z).length() < 1e-12, "u={u} dv={dv:?}");
        }
    }

    #[test]
    fn rational_bspline_surface_exact_dn_second_order() {
        let s = Surface3::BSpline(rational_arc_bspline_patch());
        for u in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let (_f, df, ddf) = rational_arc_exact(u);
            let d1u = s.dn(u, 0.5, 1, 0);
            assert!((d1u - df).length() < 1e-12, "u={u} d1u={d1u:?} want={df:?}");
            let d2u = s.dn(u, 0.5, 2, 0);
            assert!(
                (d2u - ddf).length() < 1e-12,
                "u={u} d2u={d2u:?} want={ddf:?}"
            );
            // The mixed derivative vanishes: the V direction is a straight
            // line and the U weights do not depend on V.
            let duv = s.dn(u, 0.5, 1, 1);
            assert!(duv.length() < 1e-12, "u={u} duv={duv:?}");
        }
    }

    /// The Bezier arms must agree with the same closed form analytically (OCCT
    /// `Geom_BezierSurface::EvalD0/D1/D2`, all `BSplSLib::D0/D1/D2`).
    #[test]
    fn bezier_surface_exact_derivatives() {
        let bez = rational_arc_bezier_patch();
        for u in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let (f, df, ddf) = rational_arc_exact(u);
            let p = bez.point_at(u, 0.5);
            assert!((p - (f + DVec3::Z * 0.5)).length() < 1e-12, "u={u} p={p:?}");

            let (p2, du, dv) = bez.derivatives(u, 0.5);
            assert!((p2 - (f + DVec3::Z * 0.5)).length() < 1e-12);
            assert!((du - df).length() < 1e-12, "u={u} du={du:?} want={df:?}");
            assert!((dv - DVec3::Z).length() < 1e-12, "u={u} dv={dv:?}");

            let (_p3, _du3, _dv3, duu, duv, dvv) = bez.derivatives2(u, 0.5);
            assert!(
                (duu - ddf).length() < 1e-12,
                "u={u} duu={duu:?} want={ddf:?}"
            );
            assert!(duv.length() < 1e-12, "u={u} duv={duv:?}");
            assert!(dvv.length() < 1e-12, "u={u} dvv={dvv:?}");
        }
        // The `Surface3::dn` Bezier arm is the same leaf.
        let s = Surface3::Bezier(rational_arc_bezier_patch());
        let (_f, df, ddf) = rational_arc_exact(0.5);
        let d = s.dn(0.5, 0.5, 1, 0);
        assert!((d - df).length() < 1e-12, "d={d:?} want={df:?}");
        let d = s.dn(0.5, 0.5, 2, 0);
        assert!((d - ddf).length() < 1e-12, "d={d:?} want={ddf:?}");
    }

    /// `Geom_BezierSurface::EvalD0` is the rational TENSOR evaluation
    /// `N(u, v) / W(u, v)` (`BSplSLib::D0`), not the per-column rational
    /// u-evaluation followed by a unit-weight V-combine that this module used
    /// before: the two differ as soon as the weights vary in both directions.
    /// Here every weight is 1 except the centre one (`w[1][1] = 2`), so the
    /// tensor evaluation gives `(1, 1, 0)`, while the older scheme gave
    /// `(0.75, 1, 0)`.
    #[test]
    fn bezier_rational_point_is_the_tensor_evaluation() {
        let mut control_points = Vec::new();
        for i in 0..3 {
            let mut row = Vec::new();
            for j in 0..3 {
                row.push(DVec3::new(i as f64, j as f64, 0.0));
            }
            control_points.push(row);
        }
        let mut weights = vec![vec![1.0; 3]; 3];
        weights[1][1] = 2.0;
        let bez = BezierSurface { control_points, weights };
        let p = bez.point_at(0.5, 0.5);
        assert!((p - DVec3::new(1.0, 1.0, 0.0)).length() < 1e-12, "p={p:?}");
        // The same value must come out of the `Surface3` DN arm at (0, 0) --
        // no, at the point itself: `dn(0, 0)` is outside the OCCT contract, so
        // compare against the rational D0 of the BSpline form instead.
        let bs = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: bez.control_points.clone(),
            weights: bez.weights.clone(),
            is_periodic_u: false,
            is_periodic_v: false,
        };
        let p_bs = Surface3::BSpline(bs).point_at(0.5, 0.5);
        assert!((p_bs - p).length() < 1e-12, "p_bs={p_bs:?}");
    }
}
