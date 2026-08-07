//! Boolean plane-pcurve chain: project a 3D curve onto a plane and build the
//! 2D pcurve (OCCT BRepLib::BuildPCurveForEdgesOnPlane plane branch).
//!
//! Chain (all OCCT line refs below are from this translation's source):
//!
//!   BRep_Tool::CurveOnPlane (BRep_Tool.cxx L379-450)
//!     -> GeomProjLib::ProjectOnPlane (GeomProjLib.cxx L293-346)
//!     -> ProjLib_ProjectOnPlane::Load (ProjLib_ProjectOnPlane.cxx L513-968)
//!     -> ProjLib_ProjectedCurve::Perform, plane case (ProjLib_ProjectedCurve.cxx L391-395)
//!     -> ProjLib_Plane::Project (ProjLib_Plane.cxx L101-169)
//!     -> Geom2dAdaptor::MakeCurve (Geom2dAdaptor.cxx L33-117) + BasisCurve unwrap
//!        (BRep_Tool.cxx L443-447)
//!
//! Only projected curves of analytic type (Line / Circle / Ellipse / Parabola /
//! Hyperbola) yield a 2D pcurve: a projected BSpline / Bezier / Offset / Other
//! curve falls through ProjLib_ProjectedCurve::Project (L262-266) leaving the
//! result at GeomAbs_OtherCurve, and Geom2dAdaptor::MakeCurve throws
//! (Geom2dAdaptor.cxx L92-93).  Mirror here: those branches return None and no
//! pcurve is stored (OCCT keeps `bToUpdate = false`, BRepLib_1.cxx L338).
//!
//! Disclosed sub-steps, not translated (all unreachable in the boolean plane
//! path because BRep_Tool::CurveOnPlane passes KeepParametrization = true):
//! - PerformApprox (ProjLib_ProjectOnPlane.cxx L297-412, Approx_FitAndDivide)
//! - the non-keepParam eigen-axes path (L688-738, math_Jacobi)
//! - BuildParabolaByApex / BuildHyperbolaByApex (L1432-1542, LProp_CLProps3d)
//! - BuildByApprox (L1546-1567)

use glam::{DVec2, DVec3};

use crate::core::precision::{ANGULAR, APPROXIMATION, CONFUSION};
use crate::geom::{
    Circle2d, Circle3, Curve2d, Curve3, Ellipse2d, Ellipse3, Hyperbola2d, Hyperbola3, Line2d,
    Line3, Parabola2d, Parabola3, Plane,
};

// ============================================================================
// Projection primitives (ProjLib_ProjectOnPlane.cxx L486-509)
// ============================================================================

/// OCCT ProjLib_ProjectOnPlane::ProjectPnt (L486-497): project a point along
/// `dir` onto the plane (the plane is the target, `dir` the projection ray).
fn project_pnt(pl: &Plane, dir: DVec3, point: DVec3) -> DVec3 {
    let po = pl.origin - point; // PO(Point, Location) = Location - Point
    let alpha = po.dot(pl.normal) / dir.dot(pl.normal);
    point + alpha * dir
}

/// OCCT ProjLib_ProjectOnPlane::ProjectVec (L501-509): project a direction
/// along `dir` onto the plane.
fn project_vec(pl: &Plane, dir: DVec3, vec: DVec3) -> DVec3 {
    let z = pl.normal;
    vec - (vec.dot(z) / dir.dot(z)) * dir
}

// ============================================================================
// 2D evaluation helpers (ProjLib_Plane.cxx L87-97)
// ============================================================================

/// OCCT ProjLib_Plane::EvalPnt2d (L87-92): plane point -> (u, v).
fn eval_pnt2d(p: DVec3, pl: &Plane) -> DVec2 {
    let op = p - pl.origin;
    DVec2::new(op.dot(pl.u_dir), op.dot(pl.v_dir))
}

/// OCCT ProjLib_Plane::EvalDir2d (L94-97): direction -> (u, v), normalized to
/// unit (the OCCT return type is gp_Dir2d).
fn eval_dir2d(d: DVec3, pl: &Plane) -> DVec2 {
    DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir)).normalize_or_zero()
}

/// OCCT gp_Ax22d(P, Vx, Vy) constructor (gp_Ax22d.hxx L60-74): the Y direction
/// is re-oriented to be perpendicular to X; the input Vy only selects the side
/// (`Vx.Crossed(Vy) >= 0` -> R_ccw_90(X), else R_cw_90(X)).
fn ax22d_y_dir(x: DVec2, y: DVec2) -> DVec2 {
    let det = x.x * y.y - x.y * y.x;
    if det >= 0.0 {
        DVec2::new(-x.y, x.x) // R_ccw_90(x)
    } else {
        DVec2::new(x.y, -x.x) // R_cw_90(x)
    }
}

/// OCCT gp_Vec::IsNormal (gp_Vec.hxx L480-484): |PI/2 - Angle(a, b)| <= tol.
/// (Raises Standard_ConstructionError when |a| or |b| <= gp::Resolution; rcad
/// degrades that to false instead.)
fn vec_is_normal(a: DVec3, b: DVec3, tol: f64) -> bool {
    let ang = (a.dot(b) / (a.length() * b.length())).clamp(-1.0, 1.0).acos();
    (std::f64::consts::FRAC_PI_2 - ang).abs() <= tol
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L142-146): ang <= tol || PI - ang <= tol.
/// (Same zero-length raise divergence as `vec_is_normal`; degrades to false.)
fn vec_is_parallel(a: DVec3, b: DVec3, tol: f64) -> bool {
    let ang = (a.dot(b) / (a.length() * b.length())).clamp(-1.0, 1.0).acos();
    ang <= tol || std::f64::consts::PI - ang <= tol
}

// ============================================================================
// ProjLib_ProjectOnPlane::Load (ProjLib_ProjectOnPlane.cxx L513-968)
// ============================================================================

/// OCCT ProjLib_ProjectOnPlane::Load + GeomProjLib::ProjectOnPlane dispatch
/// (L513-968, GeomProjLib.cxx L306-337), restricted to the analytic cases that
/// yield a 2D pcurve.
///
/// `range` is the edge parameter range [f, l] — the trim bounds of the
/// Geom_TrimmedCurve wrapper (BRep_Tool::CurveOnPlane L432).  OCCT uses
/// `myCurve->FirstParameter()/LastParameter()` for it only to size the
/// transient TrimmedCurve wrapper (GeomProjLib.cxx L339-343) and the
/// degenerate-line BSpline knots; neither produces a 2D pcurve, so `range`
/// feeds no returned geometry (the caller stores the pcurve range separately).
///
/// `dir` is the projection direction (the plane normal on the boolean path).
///
/// Returns None where OCCT produces a BSpline / Bezier / Other result: those
/// make Geom2dAdaptor::MakeCurve throw, so no pcurve is stored.
pub fn project_on_plane(
    curve: &Curve3,
    range: [f64; 2],
    pl: &Plane,
    dir: DVec3,
) -> Option<Curve3> {
    // GeomAdaptor_Curve unwraps a Geom_TrimmedCurve to its basis curve
    // (GeomAdaptor_Curve.cxx Load), keeping the trim range as First/Last.
    let curve = match curve {
        Curve3::Trimmed(tc) => tc.basis_curve(),
        other => other,
    };
    match curve {
        Curve3::Line(l) => {
            // L545-635 (GeomAbs_Line case).
            let xc = project_vec(pl, dir, l.direction);
            if xc.length() < CONFUSION {
                // L553-570: line orthogonal to the plane -> degenerate degree-1
                // BSpline (2 identical poles, knots [f, l]).  BSpline -> no 2D
                // pcurve (Geom2dAdaptor::MakeCurve throws).
                let _ = range;
                None
            } else if (xc.length() - 1.0).abs() < CONFUSION {
                // L571-586: |Xc| ~ 1 -> the projected line stays a line.
                let p = project_pnt(pl, dir, l.origin);
                Some(Curve3::Line(Line3::new(p, xc)))
            } else {
                // L587-632: keepParam=true (CurveOnPlane L435) -> linear BSpline
                // between the projected endpoints, reparametrized to [f, l].
                // BSpline -> no 2D pcurve.
                let _ = range;
                None
            }
        }
        Curve3::Circle(c) => {
            // L636-802 (GeomAbs_Circle case; R1 = R2 = R, L645).
            let vdx = project_vec(pl, dir, c.x_dir);
            let vdy = project_vec(pl, dir, c.y_dir);
            let tol2 = APPROXIMATION * APPROXIMATION;
            let mut is_approx = vdx.length_squared() < tol2
                || vdy.length_squared() < tol2
                || vdx.cross(vdy).length_squared() < tol2; // L664-669
            if !is_approx {
                let dx = vdx.normalize();
                let dy = vdy.normalize();
                let o = c.center;
                let p = project_pnt(pl, dir, o);
                let px = project_pnt(pl, dir, o + c.radius * c.x_dir);
                let py = project_pnt(pl, dir, o + c.radius * c.y_dir);
                let major = p.distance(px);
                let minor = p.distance(py);
                // keepParam=true -> L682-685: isApprox = !IsNormal(Dx, Dy).
                is_approx = !vec_is_normal(dx, dy, ANGULAR);
                if !is_approx {
                    let normal = dx.cross(dy);
                    // L743-769: canonical Circle or Ellipse in the plane.
                    if (major - minor).abs() < CONFUSION {
                        return Some(Curve3::Circle(Circle3 {
                            center: p,
                            normal,
                            x_dir: dx,
                            y_dir: dy,
                            radius: major,
                        }));
                    } else if major > minor {
                        return Some(Curve3::Ellipse(Ellipse3 {
                            center: p,
                            normal,
                            major_dir: dx,
                            major_radius: major,
                            minor_radius: minor,
                        }));
                    } else {
                        is_approx = true; // L766-769
                    }
                }
            }
            // L773-782: no canonical curve -> PerformApprox BSpline -> None.
            let _ = is_approx;
            None
        }
        Curve3::Ellipse(e) => {
            // L636-802 (GeomAbs_Ellipse case; R1 = Major, R2 = Minor, L649-655).
            let vdx = project_vec(pl, dir, e.major_dir);
            let vdy = project_vec(pl, dir, e.normal.cross(e.major_dir).normalize());
            let tol2 = APPROXIMATION * APPROXIMATION;
            let mut is_approx = vdx.length_squared() < tol2
                || vdy.length_squared() < tol2
                || vdx.cross(vdy).length_squared() < tol2; // L664-669
            if !is_approx {
                let dx = vdx.normalize();
                let dy = vdy.normalize();
                let o = e.center;
                let p = project_pnt(pl, dir, o);
                let px = project_pnt(pl, dir, o + e.major_radius * e.major_dir);
                let py = project_pnt(
                    pl,
                    dir,
                    o + e.minor_radius * e.normal.cross(e.major_dir).normalize(),
                );
                let major = p.distance(px);
                let minor = p.distance(py);
                // keepParam=true -> L682-685.
                is_approx = !vec_is_normal(dx, dy, ANGULAR);
                if !is_approx {
                    let normal = dx.cross(dy);
                    if (major - minor).abs() < CONFUSION {
                        return Some(Curve3::Circle(Circle3 {
                            center: p,
                            normal,
                            x_dir: dx,
                            y_dir: dy,
                            radius: major,
                        }));
                    } else if major > minor {
                        return Some(Curve3::Ellipse(Ellipse3 {
                            center: p,
                            normal,
                            major_dir: dx,
                            major_radius: major,
                            minor_radius: minor,
                        }));
                    } else {
                        is_approx = true; // L766-769
                    }
                }
            }
            let _ = is_approx;
            None
        }
        Curve3::Parabola(p) => {
            // L804-853 (GeomAbs_Parabola case).
            let xc = project_vec(pl, dir, p.axis_dir);
            let yc = project_vec(pl, dir, p.normal.cross(p.axis_dir).normalize());
            let pp = project_pnt(pl, dir, p.vertex);
            // OCCT L815-843: pick the canonical projected curve, else approx.
            let mut is_approx = false;
            let result: Option<Curve3> =
                if (yc.length() - 1.0).abs() < CONFUSION && xc.length() < CONFUSION {
                    // L817-823: |Yc| ~ 1 and |Xc| ~ 0 -> projected line.
                    Some(Curve3::Line(Line3::new(pp, yc)))
                } else if vec_is_normal(xc, yc, ANGULAR) {
                    // L824-830: Xc perp Yc -> projected parabola.
                    // F = Focal / |Xc| (L827); rcad focal_param = 2 * OCCT Focal.
                    Some(Curve3::Parabola(Parabola3 {
                        vertex: pp,
                        normal: xc.cross(yc).normalize(),
                        axis_dir: xc.normalize(),
                        focal_param: p.focal_param / xc.length(),
                    }))
                } else if yc.length() < CONFUSION || vec_is_parallel(yc, xc, ANGULAR) {
                    // L831-834: |Yc| < Confusion or Yc parallel Xc -> approx.
                    is_approx = true;
                    None
                } else {
                    // L835-843: !keepParam -> BuildParabolaByApex (L1432-1542,
                    // disclosed; unreachable with KeepParametrization=true,
                    // CurveOnPlane L435); keepParam=true -> approx.
                    is_approx = true;
                    None
                };
            if is_approx {
                // L849-852: approx -> BuildByApprox BSpline -> no 2D pcurve.
                None
            } else {
                // L845-848: canonical -> GetTrimmedResult.
                result
            }
        }
        Curve3::Hyperbola(h) => {
            // L855-911 (GeomAbs_Hyperbola case).
            let xc = project_vec(pl, dir, h.major_dir);
            let yc = project_vec(pl, dir, h.normal.cross(h.major_dir).normalize());
            let pp = project_pnt(pl, dir, h.center);
            let a_r1 = h.semi_major;
            let a_r2 = h.semi_minor;
            let z = pl.normal;
            // OCCT L868-902: pick the canonical projected curve, else approx.
            let mut is_approx = false;
            let result: Option<Curve3> = if xc.length() < CONFUSION {
                // L870-876: |Xc| ~ 0 -> hyperbola with transverse axis Yc ^ Z.
                let x = yc.normalize_or_zero().cross(z);
                Some(Curve3::Hyperbola(Hyperbola3 {
                    center: pp,
                    normal: z,
                    major_dir: x,
                    semi_major: 0.0,
                    semi_minor: a_r2 * yc.length(),
                }))
            } else if yc.length() < CONFUSION {
                // L877-882: |Yc| ~ 0.
                Some(Curve3::Hyperbola(Hyperbola3 {
                    center: pp,
                    normal: z,
                    major_dir: xc.normalize(),
                    semi_major: a_r1 * xc.length(),
                    semi_minor: 0.0,
                }))
            } else if vec_is_normal(xc, yc, ANGULAR) {
                // L883-890: Xc perp Yc.
                Some(Curve3::Hyperbola(Hyperbola3 {
                    center: pp,
                    normal: xc.cross(yc).normalize(),
                    major_dir: xc.normalize(),
                    semi_major: a_r1 * xc.length(),
                    semi_minor: a_r2 * yc.length(),
                }))
            } else if yc.length() < CONFUSION || vec_is_parallel(yc, xc, ANGULAR) {
                // L891-894: |Yc| < Confusion or Yc parallel Xc -> approx.
                is_approx = true;
                None
            } else {
                // L895-902: !keepParam -> BuildHyperbolaByApex (L1491-1542,
                // disclosed; unreachable with KeepParametrization=true);
                // keepParam=true -> approx.
                is_approx = true;
                None
            };
            if is_approx {
                // L907-910: approx -> BuildByApprox BSpline -> no 2D pcurve.
                None
            } else {
                // L903-906: canonical -> GetTrimmedResult.
                result
            }
        }
        // L913-966: Bezier / BSpline (projected poles, keepParam forced true) and
        // the default case (PerformApprox) all give a BSpline / Bezier result ->
        // no 2D pcurve.
        _ => None,
    }
}

// ============================================================================
// 2D conversion (ProjLib_Plane::Project + Geom2dAdaptor::MakeCurve)
// ============================================================================

/// OCCT ProjLib_ProjectedCurve::Perform plane case (L391-395) +
/// ProjLib_Plane::Project (ProjLib_Plane.cxx L101-169) +
/// Geom2dAdaptor::MakeCurve (Geom2dAdaptor.cxx L33-117) + BasisCurve unwrap
/// (BRep_Tool.cxx L443-447).
///
/// Returns the bare 2D conic in the plane's (u, v) parameter space.  The
/// MakeCurve trim step (L97-113) is ignored: the result is unwrapped to the
/// basis conic by BRep_Tool::CurveOnPlane L443-447, and the caller stores the
/// pcurve range separately.
///
/// For Ellipse / Parabola / Hyperbola the 2D minor axis follows rcad's
/// `rotate_ccw_90(major)` convention, which matches OCCT's gp_Ax22d
/// re-orientation only when `det(EvalDir2d(Dx), EvalDir2d(Dy)) >= 0` (gp_Ax22d
/// re-orients Y to R_cw_90(X) when the determinant is negative, i.e. the
/// projected minor axis lands on the clockwise side).  rcad's 2D conics carry
/// no explicit minor-direction field, so the clockwise case cannot be
/// represented; it returns None (no pcurve) rather than a mirrored conic.
fn projected_curve_to_2d(proj: &Curve3, pl: &Plane) -> Option<Curve2d> {
    match proj {
        Curve3::Line(l) => {
            // Project(gp_Lin) L101-106: no Ax22d re-orientation.
            Some(Curve2d::Line(Line2d::new(
                eval_pnt2d(l.origin, pl),
                eval_dir2d(l.direction, pl),
            )))
        }
        Curve3::Circle(c) => {
            // Project(gp_Circ) L110-123: gp_Ax22d re-orients Y perp to X.
            let p2d = eval_pnt2d(c.center, pl);
            let x2d = eval_dir2d(c.x_dir, pl);
            let y2d = eval_dir2d(c.y_dir, pl);
            Some(Curve2d::Circle(Circle2d {
                center: p2d,
                x_dir: x2d,
                y_dir: ax22d_y_dir(x2d, y2d),
                radius: c.radius,
            }))
        }
        Curve3::Ellipse(e) => {
            // Project(gp_Elips) L127-139: axes EvalDir2d(X), EvalDir2d(Y).
            let p2d = eval_pnt2d(e.center, pl);
            let x2d = eval_dir2d(e.major_dir, pl);
            let y2d = eval_dir2d(e.normal.cross(e.major_dir).normalize(), pl);
            if x2d.x * y2d.y - x2d.y * y2d.x >= 0.0 {
                Some(Curve2d::Ellipse(Ellipse2d {
                    center: p2d,
                    major_dir: x2d,
                    major_radius: e.major_radius,
                    minor_radius: e.minor_radius,
                }))
            } else {
                None
            }
        }
        Curve3::Parabola(p) => {
            // Project(gp_Parab) L143-154: axes EvalDir2d(X), EvalDir2d(Y).
            let p2d = eval_pnt2d(p.vertex, pl);
            let x2d = eval_dir2d(p.axis_dir, pl);
            let y2d = eval_dir2d(p.normal.cross(p.axis_dir).normalize(), pl);
            if x2d.x * y2d.y - x2d.y * y2d.x >= 0.0 {
                Some(Curve2d::Parabola(Parabola2d {
                    origin: p2d,
                    axis_dir: x2d,
                    focal_param: p.focal_param,
                }))
            } else {
                None
            }
        }
        Curve3::Hyperbola(h) => {
            // Project(gp_Hypr) L158-169: axes EvalDir2d(X), EvalDir2d(Y).
            let p2d = eval_pnt2d(h.center, pl);
            let x2d = eval_dir2d(h.major_dir, pl);
            let y2d = eval_dir2d(h.normal.cross(h.major_dir).normalize(), pl);
            if x2d.x * y2d.y - x2d.y * y2d.x >= 0.0 {
                Some(Curve2d::Hyperbola(Hyperbola2d {
                    center: p2d,
                    major_dir: x2d,
                    semi_major: h.semi_major,
                    semi_minor: h.semi_minor,
                }))
            } else {
                None
            }
        }
        // BSpline / Bezier / Offset / Other -> Geom2dAdaptor::MakeCurve default
        // throws (L92-93): no pcurve.
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnPlane (BRep_Tool.cxx L379-450): project the edge's
/// 3D curve onto the plane and return the 2D pcurve in the plane's (u, v)
/// parameter space.
///
/// Equivalent to the full chain with Dir = plane normal and
/// KeepParametrization = true (L432-435), then ProjLib_ProjectedCurve +
/// Geom2dAdaptor::MakeCurve.
///
/// The location transform (BRep_Tool.cxx L421-428) is skipped: TEdgeData has
/// no TopLoc_Location (identity location in the boolean pipeline).
///
/// Returns None where OCCT stores no pcurve (Geom2dAdaptor::MakeCurve throws,
/// L92-93; the edge keeps `bToUpdate = false`).
pub fn curve_on_plane(curve: &Curve3, range: [f64; 2], pl: &Plane) -> Option<Curve2d> {
    let proj = project_on_plane(curve, range, pl, pl.normal)?;
    projected_curve_to_2d(&proj, pl)
}
