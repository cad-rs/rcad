// OCCT BndLib_Add2dCurve (BndLib package, TKGeomBase) — the 2D curve
// bounding-box filling used by BRepClass_Intersector::Perform
// (BRepClass_Intersector.cxx L348-356, L392-396): the edge pcurve's bbox is
// tested against the ray line/segment before the full intersection.
//
// Stand-in: the exact-bounds dispatch of BndLib_Add2dCurve::Add
// (BndLib_Add2dCurve.cxx — the per-curve-type ElCLib / Perel /
// Geom2dLProp_CurAndInf2d bounds) is not translated yet.  The rcad form
// keeps the rejection behaviour: Line and Circle use exact bounds, other
// curve types are sampled (conservative).
//
// The Bnd_Box2d data structure lives in its OCCT home,
// rcad-kernel/src/math/bnd (TKMath/Bnd); this module no longer declares a
// local copy of it.

use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::math::bnd::BndBox2d;

/// OCCT BndLib_Add2dCurve::Add(C, U1, U2, Tol, Box) — build the curve's 2D
/// bounding box over [the_u1, the_u2].
pub fn add_2d_curve(
    the_c: &Curve2d,
    the_u1: f64,
    the_u2: f64,
    _the_tol: f64,
    the_box: &mut BndBox2d,
) {
    match the_c {
        Curve2d::Line(_) => {
            the_box.add_point(the_c.point_at(the_u1));
            the_box.add_point(the_c.point_at(the_u2));
        }
        Curve2d::Circle(c) => {
            // Exact bbox: endpoints plus the 4 cardinal frame points when
            // inside the arc.
            the_box.add_point(the_c.point_at(the_u1));
            the_box.add_point(the_c.point_at(the_u2));
            for i in 0..4 {
                let ang = i as f64 * std::f64::consts::FRAC_PI_2;
                // The point at `ang` in the circle frame (x_dir/y_dir).
                let p = c.center
                    + c.x_dir * (c.radius * ang.cos())
                    + c.y_dir * (c.radius * ang.sin());
                // Include only when the angle lies inside the arc domain.
                if ang >= the_u1 - 1e-12 && ang <= the_u2 + 1e-12 {
                    the_box.add_point(p);
                }
            }
        }
        _ => {
            the_box.add_point(the_c.point_at(the_u1));
            the_box.add_point(the_c.point_at(the_u2));
            const N: usize = 64;
            for i in 1..N {
                let t = the_u1 + (the_u2 - the_u1) * (i as f64) / (N as f64);
                the_box.add_point(the_c.point_at(t));
            }
        }
    }
}
