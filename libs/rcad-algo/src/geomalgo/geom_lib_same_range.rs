//! OCCT GeomLib statics (TKGeomBase/GeomLib — GeomLib.cxx):
//! - SameRange (GeomLib.cxx L842-970) — rebuilds a 2d curve over a target
//!   parameter range (consumed by Filling cxx L1037 and UpdateEdge
//!   cxx L1767).  The 1:1 body lives in the kernel as
//!   `rcad_kernel::geom::same_range_2d`; this module is the OCCT-signature
//!   entry point (Tolerance first, `NewCurvePtr` out-parameter) used by the
//!   BRepFill/offset callers, so a single implementation serves both.
//! - ExtendSurfByLength (GeomLib.cxx L1485-1972) — extends a bounded surface
//!   along an iso (consumed by BuildShell cxx L2334/L2341). The 1:1 body is
//!   `fillet::chfi3d_builder_c2_geomlib::geom_lib_extend_surf_by_length`;
//!   this module is the OCCT-signature value-semantics entry point.

use rcad_kernel::geom::{Curve2d, Surface3};

/// OCCT GeomLib::SameRange(Tolerance, Curve2d, First, Last, NewFirst,
/// NewLast, NewCurve2d) (GeomLib.cxx L842-970).
///
/// Returns the OCCT `NewCurvePtr` result; the OCCT null-handle outcome (the
/// L902-922 guard leaving `NewCurvePtr` untouched, or a failed
/// `CurveToBSplineCurve`) is the failure the callers observe as a null
/// handle.
pub fn same_range(
    tolerance: f64,
    curve2d: &Curve2d,
    first: f64,
    last: f64,
    new_first: f64,
    new_last: f64,
) -> Curve2d {
    match rcad_kernel::geom::same_range_2d(
        tolerance,
        curve2d.clone(),
        first,
        last,
        new_first,
        new_last,
    ) {
        Some(c) => c,
        None => panic!(
            "Standard_NullObject: GeomLib::SameRange produced a null curve \
             (GeomLib.cxx L842-970)"
        ),
    }
}

/// OCCT GeomLib::ExtendSurfByLength(BoundedSurface, Length, Continuity,
/// InU, After) (GeomLib.cxx L1485-1972).
///
/// OCCT takes the surface by handle and extends it in place; rcad's value
/// semantics return the (possibly extended) surface. The 1:1 body is
/// `fillet::chfi3d_builder_c2_geomlib::geom_lib_extend_surf_by_length` (the
/// same OCCT static, first translated for the fillet BuildShell callers), so
/// a single implementation serves both; when it reports false the surface is
/// left exactly as OCCT leaves it.
pub fn extend_surf_by_length(
    surface: &Surface3,
    length: f64,
    continuity: i32,
    in_u: bool,
    after: bool,
) -> Surface3 {
    let mut extended = surface.clone();
    crate::fillet::chfi3d_builder_c2_geomlib::geom_lib_extend_surf_by_length(
        &mut extended,
        length,
        continuity,
        in_u,
        after,
    );
    extended
}

// ===========================================================================
// OCCT Geom2dConvert::C0BSplineToC1BSplineCurve (TKGeomBase/Geom2dConvert/
// Geom2dConvert.cxx L1459-1541) — reduces as far as possible the
// multiplicities of the knots of the BSpline BS (keeping the geometry) and
// returns a new BSpline which could still be C0.  The 1:1 body lives in the
// kernel (`geom2d_convert::c1_concat`); this module is the OCCT-signature
// entry point, with the OCCT handle out-parameter mapped to value semantics.
// ===========================================================================

use rcad_kernel::base::geom2d_convert::c1_concat::{
    comp_curve_add, concat_c1_default, curve_knots, curve_mults, gp_vec2d_is_parallel,
};
use rcad_kernel::base::geom2d_convert::Geom2dBSplineCurve;
use rcad_kernel::core::precision::PCONFUSION;

/// OCCT Geom2dConvert::C0BSplineToC1BSplineCurve(BS, tolerance)
/// (Geom2dConvert.cxx L1459-1541).  OCCT mutates the handle in place;
/// rcad's value semantics return the (possibly untouched) curve.
pub fn c0_bspline_to_c1_bspline_curve(
    bs: &Geom2dBSplineCurve,
    tolerance: f64,
) -> Geom2dBSplineCurve {
    let bs_mults = curve_mults(bs);
    let bs_knots = curve_knots(bs);
    let mut nbcurve_c1 = 1usize;
    // for (i = FirstUKnotIndex + 1; i <= LastUKnotIndex - 1; i++)
    for i in (bs.first_uknot_index() + 1)..=(bs.last_uknot_index() - 1) {
        if bs_mults[(i - 1) as usize] == bs.degree() as i32 {
            nbcurve_c1 += 1;
        }
    }

    nbcurve_c1 = nbcurve_c1.min((bs.nb_knots() - 1) as usize);

    if nbcurve_c1 > 1 {
        let mut array_of_curves: Vec<Geom2dBSplineCurve> = Vec::new();
        let array_of_toler = vec![tolerance; nbcurve_c1 - 1];

        let mut u2 = bs.first_parameter();
        let mut j = bs.first_uknot_index() + 1;
        for _i in 0..nbcurve_c1 {
            let u1 = u2;

            while j < bs.last_uknot_index()
                && bs_mults[(j - 1) as usize] < bs.degree() as i32
            {
                j += 1;
            }

            u2 = bs_knots[(j - 1) as usize];
            j += 1;
            // BSbis = BS->Copy(); BSbis->Segment(U1, U2) with the OCCT
            // default tolerance Precision::PConfusion().
            let mut bs_bis = bs.clone();
            bs_bis.segment(u1, u2, PCONFUSION);
            array_of_curves.push(bs_bis);
        }

        let an_angular_toler = 1.0e-7;

        let (point1, v1) = bs.eval_d1(bs.first_parameter()); // a verifier
        let (point2, v2) = bs.eval_d1(bs.last_parameter());

        let mut closed_flag = false;
        if (point1.distance_squared(point2) < tolerance * tolerance)
            && gp_vec2d_is_parallel(v1, v2, an_angular_toler)
        {
            closed_flag = true;
        }

        let (_array_of_indices, array_of_concatenated) =
            concat_c1_default(&mut array_of_curves, &array_of_toler, &mut closed_flag, tolerance);

        let mut result = array_of_concatenated[0].clone();
        if array_of_concatenated.len() >= 2 {
            for curve in array_of_concatenated.iter().skip(1) {
                // fusion = C.Add(Value(i), tolerance, true) — the after arm;
                // the OCCT failure throw mirrors comp_curve_add.
                result = comp_curve_add(&result, curve, tolerance, true);
            }
        }
        result
    } else {
        bs.clone()
    }
}

#[cfg(test)]
mod c0_to_c1_tests {
    use super::*;
    use rcad_kernel::base::geom2d_convert::c1_concat::curve_knots as knots_of;
    use glam::DVec2;

    /// Hand-derived C0-to-C1 reduction of the collinear cubic
    /// knots [0,1,2] x mults [4,3,4] with poles x = [0,1,2,3 | 4,5,6]
    /// (y = 0): the junction multiplicity 3 == Degree splits it into the
    /// Beziers [0,1,2,3] on [0,1] and [3,4,5,6] on [1,2] (equal end
    /// tangents -> C1, non-closed).  ConcatC1 merges them and the junction
    /// knot is removed down to 0 (every deviation is 0 on the straight
    /// line x = 3u), giving the single cubic on [0,2] with poles
    /// x = [0, 2, 4, 6], knots [0,2] x mults [4,4].
    #[test]
    fn collinear_cubic_reduces_to_single_bezier() {
        let poles: Vec<DVec2> = (0..7).map(|i| DVec2::new(i as f64, 0.0)).collect();
        let bs = Geom2dBSplineCurve::new(poles, vec![0.0, 1.0, 2.0], vec![4, 3, 4], 3, false);
        let out = c0_bspline_to_c1_bspline_curve(&bs, 1.0e-7);
        assert_eq!(out.degree(), 3);
        assert_eq!(out.nb_knots(), 2);
        assert_eq!((out.knot(1), out.knot(2)), (0.0, 2.0));
        assert_eq!((out.multiplicity(1), out.multiplicity(2)), (4, 4));
        assert_eq!(out.nb_poles_curve(), 4);
        for (i, expect) in [0.0f64, 2.0, 4.0, 6.0].iter().enumerate() {
            let p = out.pole(i as i32 + 1);
            assert!((p.x - expect).abs() < 1e-12, "pole{} x={}", i + 1, p.x);
            assert!(p.y.abs() < 1e-12, "pole{} y={}", i + 1, p.y);
        }
        // x = 3u on [0, 2].
        let p = out.eval_d0(1.0);
        assert!((p.x - 3.0).abs() < 1e-12 && p.y.abs() < 1e-12);
        let _ = knots_of(&out);
    }

    /// A curve with no interior knot of multiplicity == Degree
    /// (mults [4,1,4]) has nbcurveC1 == 1 and is returned untouched
    /// (OCCT L1481 guard: the body never runs).
    #[test]
    fn single_span_curve_untouched() {
        let poles: Vec<DVec2> = (0..5).map(|i| DVec2::new(i as f64, 1.0)).collect();
        let bs = Geom2dBSplineCurve::new(poles, vec![0.0, 1.0, 2.0], vec![4, 1, 4], 3, false);
        let out = c0_bspline_to_c1_bspline_curve(&bs, 1.0e-7);
        assert_eq!(out.nb_knots(), 3);
        assert_eq!(out.nb_poles_curve(), 5);
        assert_eq!(out.pole(3), DVec2::new(2.0, 1.0));
    }
}
