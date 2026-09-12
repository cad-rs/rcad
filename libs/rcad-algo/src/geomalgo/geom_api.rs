//! OCCT GeomAPI (TKGeomAlgo/GeomAPI) — the 1:1 translation of the two
//! conversion statics of GeomAPI.hxx:
//!
//! - `GeomAPI::To2d(C, P)`  (GeomAPI.cxx L37-52)
//! - `GeomAPI::To3d(C, P)`  (GeomAPI.cxx L56-65)
//!
//! `To2d` drives `ProjLib_ProjectedCurve` on the plane adaptor, so the plane
//! case of that request is translated here as the private helper
//! `proj_lib_projected_curve_plane`:
//!
//! - `ProjLib_ProjectedCurve::Perform`, `case GeomAbs_Plane`
//!   (ProjLib_ProjectedCurve.cxx L391-395): `ProjLib_Plane P(mySurface->
//!   Plane()); Project(P, myCurve); myResult = P;` — the `Project` dispatcher
//!   (L239-264) delegates the analytic types and breaks for
//!   BSpline / Bezier / Offset / Other ("try the approximation").
//! - the `!myResult.IsDone() && isAnalyticalSurf` fallback (L707-770) into
//!   `ProjLib_ComputeApprox::Perform`, whose plane×BSpline / plane×Bezier
//!   arms are the exact pole mapping (ProjLib_ComputeApprox.cxx L1197-1255).
//! - `Geom2dAdaptor::MakeCurve` (Geom2dAdaptor.cxx L33-117).
//!
//! `To3d` drives the `Adaptor3d_CurveOnSurface` plane branch of `EvalKPart`
//! (Adaptor3d_CurveOnSurface.cxx L1555-1563: for a `GeomAbs_Plane` surface
//! `myType = myCurve->GetType()`) plus the `to3d` statics of that file
//! (L57-99), then `GeomAdaptor::MakeCurve` (GeomAdaptor.cxx) with its
//! trailing trim step.
//!
//! Failure paths are the OCCT ones: `GeomAPI::To2d` leaves the result handle
//! null when the projected type is `GeomAbs_OffsetCurve` /
//! `GeomAbs_OtherCurve` (cxx L39 + L46-49) — the `None` below — and
//! `GeomAdaptor::MakeCurve` raises
//! `Standard_DomainError("GeomAdaptor::MakeCurve : OtherCurve")`.

use glam::DVec2;
use glam::DVec3;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::base::proj_lib::{CurveType, PlaneProjector};
use rcad_kernel::geom::{
    BezierCurve2, BezierCurve3, BSplineCurve2, BSplineCurve3, Circle2d, Circle3, Curve2d,
    Curve2dEval, Curve3, CurveEval, Ellipse2d, Ellipse3, Hyperbola2d, Hyperbola3, Line2d, Line3,
    Parabola2d, Parabola3, Plane, TrimmedCurve2, TrimmedCurve3,
};
use rcad_kernel::SurfaceEval;

// =========================================================================
// GeomAPI::To3d — Adaptor3d_CurveOnSurface statics (plane case)
// =========================================================================

/// OCCT Adaptor3d_CurveOnSurface::to3d(Pl, P) (Adaptor3d_CurveOnSurface.cxx
/// L57-60): `ElSLib::Value(P.X(), P.Y(), Pl)`.
fn plane_to3d_pnt(pl: &Plane, p: DVec2) -> DVec3 {
    pl.point_at(p.x, p.y)
}

/// OCCT Adaptor3d_CurveOnSurface::to3d(Pl, V) (cxx L62-71):
/// `V.X() * XAxis + V.Y() * YAxis`.
fn plane_to3d_vec(pl: &Plane, v: DVec2) -> DVec3 {
    pl.u_dir * v.x + pl.v_dir * v.y
}

/// OCCT Adaptor3d_CurveOnSurface::to3d(Pl, A) (cxx L72-78):
/// `gp_Ax2(P, VX.Crossed(VY), VX)` — returned as (origin, normal, x_dir).
fn plane_to3d_ax2(
    pl: &Plane,
    origin: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
) -> (DVec3, DVec3, DVec3) {
    let p = plane_to3d_pnt(pl, origin);
    let vx = plane_to3d_vec(pl, x_dir);
    let vy = plane_to3d_vec(pl, y_dir);
    (p, vx.cross(vy), vx)
}

/// OCCT GeomAdaptor::MakeCurve(Adaptor3d_Curve) (GeomAdaptor.cxx) — the
/// `switch (HC.GetType())` plus the trailing trim step.  `the_type` is
/// `HC.GetType()`, `the_first` / `the_last` are `HC.FirstParameter()` /
/// `HC.LastParameter()`, and `the_c` the already-built `Geom_*` of that type.
fn geom_adaptor_make_curve(
    the_c: Curve3,
    the_first: f64,
    the_last: f64,
    the_type: CurveType,
) -> Curve3 {
    match the_type {
        // case GeomAbs_Line: C = new Geom_Line(HC.Line()); ... case
        // GeomAbs_BSplineCurve: C = HC.BSpline()->Copy();
        CurveType::Line
        | CurveType::Circle
        | CurveType::Ellipse
        | CurveType::Parabola
        | CurveType::Hyperbola
        | CurveType::Bezier
        | CurveType::BSpline => {}
        // default: throw Standard_DomainError("GeomAdaptor::MakeCurve :
        // OtherCurve");
        _ => panic!("Standard_DomainError: GeomAdaptor::MakeCurve : OtherCurve"),
    }

    // trim the curve if necessary (GeomAdaptor.cxx):
    //   if ((!C.IsNull() && (HC.FirstParameter() != C->FirstParameter()))
    //       || (HC.LastParameter() != C->LastParameter()))
    //     C = new Geom_TrimmedCurve(C, HC.FirstParameter(), HC.LastParameter());
    let a_dom = CurveEval::default_domain(&the_c);
    if the_first != a_dom[0] || the_last != a_dom[1] {
        return Curve3::Trimmed(TrimmedCurve3::new(the_c, the_first, the_last));
    }
    the_c
}

/// OCCT GeomAPI::To3d(C, P) (GeomAPI.cxx L56-65).
pub fn to3d(the_c: &Curve2d, the_p: &Plane) -> Curve3 {
    // cxx L58: handle(Geom2dAdaptor_Curve) AHC = new Geom2dAdaptor_Curve(C); —
    // Geom2dAdaptor_Curve::Load unwraps a Geom2d_TrimmedCurve to its basis and
    // keeps the trim window as (myFirst, myLast).
    let (a_basis, a_first, a_last) = match the_c {
        Curve2d::Trimmed(t) => ((*t.curve).clone(), t.t_min, t.t_max),
        other => {
            let f = Curve2dEval::default_domain(other);
            (other.clone(), f[0], f[1])
        }
    };
    // cxx L60: handle(Geom_Plane) ThePlane = new Geom_Plane(P);
    // cxx L61: handle(GeomAdaptor_Surface) AHS = new
    // GeomAdaptor_Surface(ThePlane);
    //
    // cxx L63: Adaptor3d_CurveOnSurface COS(AHC, AHS); — EvalKPart
    // (Adaptor3d_CurveOnSurface.cxx L1555-1563): for a GeomAbs_Plane surface
    // myType = myCurve->GetType(), and the analytic forms are the plane
    // `to3d` images above.
    let a_type = curve2d_adaptor_type(&a_basis);
    let a_c = match a_type {
        // Adaptor3d_CurveOnSurface::Line() (cxx L1442-1461): the D1 image of
        // the 2D line at parameter 0 (`gp_Lin(P, V)` normalizes V).
        CurveType::Line => match &a_basis {
            Curve2d::Line(l) => Curve3::Line(Line3::new(
                plane_to3d_pnt(the_p, l.point_at(0.0)),
                plane_to3d_vec(the_p, l.direction),
            )),
            _ => unreachable!(),
        },
        // cxx L80-83: `gp_Circ(to3d(Pl, C.Axis()), C.Radius())`.
        CurveType::Circle => match &a_basis {
            Curve2d::Circle(c) => {
                let (o, n, x) = plane_to3d_ax2(the_p, c.center, c.x_dir, c.y_dir);
                let y = n.cross(x);
                Curve3::Circle(Circle3 {
                    center: o,
                    normal: n,
                    x_dir: x,
                    y_dir: y,
                    radius: c.radius,
                })
            }
            _ => unreachable!(),
        },
        // cxx L85-88: `gp_Elips(to3d(Pl, E.Axis()), MajorRadius, MinorRadius)`.
        CurveType::Ellipse => match &a_basis {
            Curve2d::Ellipse(e) => {
                let (o, n, x) = plane_to3d_ax2(the_p, e.center, e.major_dir, e.minor_dir);
                Curve3::Ellipse(Ellipse3 {
                    center: o,
                    normal: n,
                    major_dir: x,
                    major_radius: e.major_radius,
                    minor_radius: e.minor_radius,
                })
            }
            _ => unreachable!(),
        },
        // cxx L90-93: `gp_Hypr(to3d(Pl, H.Axis()), MajorRadius, MinorRadius)`.
        CurveType::Hyperbola => match &a_basis {
            Curve2d::Hyperbola(h) => {
                let (o, n, x) = plane_to3d_ax2(the_p, h.center, h.major_dir, DVec2::new(0.0, 1.0));
                Curve3::Hyperbola(Hyperbola3 {
                    center: o,
                    normal: n,
                    major_dir: x,
                    semi_major: h.semi_major,
                    semi_minor: h.semi_minor,
                })
            }
            _ => unreachable!(),
        },
        // cxx L95-98: `gp_Parab(to3d(Pl, P.Axis()), P.Focal())`.
        CurveType::Parabola => match &a_basis {
            Curve2d::Parabola(p) => {
                let (o, n, x) = plane_to3d_ax2(the_p, p.origin, p.axis_dir, DVec2::new(0.0, 1.0));
                Curve3::Parabola(Parabola3 {
                    vertex: o,
                    normal: n,
                    axis_dir: x,
                    focal_param: p.focal_param,
                })
            }
            _ => unreachable!(),
        },
        // Adaptor3d_CurveOnSurface::BSpline() plane form (cxx L1590+): the
        // poles are the plane `to3d` images, knots / degree / weights kept.
        CurveType::BSpline => match &a_basis {
            Curve2d::BSpline(b) => Curve3::BSpline(BSplineCurve3 {
                degree: b.degree,
                knots: b.knots.clone(),
                control_points: b
                    .control_points
                    .iter()
                    .map(|p| plane_to3d_pnt(the_p, *p))
                    .collect(),
                weights: b.weights.clone(),
                is_periodic: false,
            }),
            _ => unreachable!(),
        },
        // Adaptor3d_CurveOnSurface::Bezier() plane form (cxx L1458+): the
        // same pole mapping.
        CurveType::Bezier => match &a_basis {
            Curve2d::Bezier(b) => Curve3::Bezier(BezierCurve3 {
                control_points: b
                    .control_points
                    .iter()
                    .map(|p| plane_to3d_pnt(the_p, *p))
                    .collect(),
                weights: b.weights.clone(),
            }),
            _ => unreachable!(),
        },
        // Every other 2D type reports GeomAbs_OtherCurve from EvalKPart, which
        // is the Standard_DomainError of GeomAdaptor::MakeCurve.
        CurveType::Other => {
            return geom_adaptor_make_curve(
                Curve3::Line(Line3::new(DVec3::ZERO, DVec3::X)),
                a_first,
                a_last,
                a_type,
            )
        }
    };
    // cxx L64: return GeomAdaptor::MakeCurve(COS);
    geom_adaptor_make_curve(a_c, a_first, a_last, a_type)
}

/// `Geom2dAdaptor_Curve::GetType()` over the rcad `Curve2d` — the
/// `CurveType` value `Adaptor3d_CurveOnSurface::EvalKPart` copies for a plane
/// surface.
fn curve2d_adaptor_type(the_c: &Curve2d) -> CurveType {
    match the_c {
        Curve2d::Line(_) => CurveType::Line,
        Curve2d::Circle(_) => CurveType::Circle,
        Curve2d::Ellipse(_) => CurveType::Ellipse,
        Curve2d::Parabola(_) => CurveType::Parabola,
        Curve2d::Hyperbola(_) => CurveType::Hyperbola,
        Curve2d::Bezier(_) => CurveType::Bezier,
        Curve2d::BSpline(_) => CurveType::BSpline,
        _ => CurveType::Other,
    }
}

// =========================================================================
// GeomAPI::To2d — ProjLib_ProjectedCurve (plane case) + Geom2dAdaptor
// =========================================================================

/// The `GeomAbs_Plane` arm of OCCT `ProjLib_ProjectedCurve::Perform`
/// (ProjLib_ProjectedCurve.cxx L391-395) plus the `ProjLib_ComputeApprox`
/// fallback for BSpline / Bezier curves (L707-770 -> the plane arms at
/// ProjLib_ComputeApprox.cxx L1197-1255).  Returns the `myResult` projected
/// curve, or `None` where OCCT leaves `myResult` not done (which is the
/// `GeomAbs_OtherCurve` the `GeomAPI::To2d` gate rejects).
fn proj_lib_projected_curve_plane(the_hc: &GeomCurveAdaptor, the_p: &Plane) -> Option<Curve2d> {
    // --- ProjLib_ProjectedCurve::Perform, case GeomAbs_Plane (L391-395):
    //     ProjLib_Plane P(mySurface->Plane()); Project(P, myCurve);
    //     myResult = P; ---
    // The `Project(ProjLib_Projector&, Adaptor3d_Curve&)` dispatcher
    // (L239-264) switches on `C->GetType()`: the analytic types call
    // `P.Project(...)`, and GeomAbs_BSplineCurve / GeomAbs_BezierCurve /
    // GeomAbs_OffsetCurve / GeomAbs_OtherCurve `break` ("try the
    // approximation") leaving myResult not done.
    let mut a_plane_proj = PlaneProjector::with_plane(the_p);
    let a_analytic = match &the_hc.curve {
        Curve3::Line(l) => {
            a_plane_proj.project_line(l);
            true
        }
        Curve3::Circle(c) => {
            a_plane_proj.project_circle(c);
            true
        }
        Curve3::Ellipse(e) => {
            a_plane_proj.project_ellipse(e);
            true
        }
        Curve3::Hyperbola(h) => {
            a_plane_proj.project_hyperbola(h);
            true
        }
        Curve3::Parabola(p) => {
            a_plane_proj.project_parabola(p);
            true
        }
        _ => false,
    };
    if a_analytic {
        // ProjLib_Plane::Project(gp_Lin/gp_Circ/gp_Elips/gp_Parab/gp_Hypr)
        // (ProjLib_Plane.cxx L101-169) — the gp_*2d results, which
        // Geom2dAdaptor::MakeCurve turns into the Geom2d_* classes below.
        let a_inner = a_plane_proj.projector();
        return Some(match a_inner.get_type() {
            CurveType::Line => {
                let l = a_inner.line();
                Curve2d::Line(Line2d::new(
                    DVec2::new(l.origin.x, l.origin.y),
                    DVec2::new(l.direction.x, l.direction.y),
                ))
            }
            CurveType::Circle => {
                let c = a_inner.circle();
                Curve2d::Circle(Circle2d {
                    center: DVec2::new(c.center.x, c.center.y),
                    x_dir: DVec2::new(c.x_dir.x, c.x_dir.y),
                    y_dir: DVec2::new(c.y_dir.x, c.y_dir.y),
                    radius: c.radius,
                })
            }
            CurveType::Ellipse => {
                let e = a_inner.ellipse();
                Curve2d::Ellipse(Ellipse2d {
                    center: DVec2::new(e.center.x, e.center.y),
                    major_dir: DVec2::new(e.major_dir.x, e.major_dir.y),
                    // OCCT gp_Ax22d keeps the minor axis as the YDirection of
                    // the projected frame (ProjLib_Plane::Project(gp_Elips),
                    // cxx L126-131).
                    minor_dir: DVec2::new(-e.major_dir.y, e.major_dir.x),
                    major_radius: e.major_radius,
                    minor_radius: e.minor_radius,
                })
            }
            CurveType::Hyperbola => {
                let h = a_inner.hyperbola();
                Curve2d::Hyperbola(Hyperbola2d {
                    center: DVec2::new(h.center.x, h.center.y),
                    major_dir: DVec2::new(h.major_dir.x, h.major_dir.y),
                    semi_major: h.semi_major,
                    semi_minor: h.semi_minor,
                })
            }
            CurveType::Parabola => {
                let p = a_inner.parabola();
                Curve2d::Parabola(Parabola2d {
                    origin: DVec2::new(p.vertex.x, p.vertex.y),
                    axis_dir: DVec2::new(p.axis_dir.x, p.axis_dir.y),
                    focal_param: p.focal_param,
                })
            }
            // The plane analytic projector never produces these.
            _ => return None,
        });
    }

    // --- if (!myResult.IsDone() && isAnalyticalSurf) (L707-770):
    //     ProjLib_ComputeApprox Comp; Comp.Perform(myCurve, mySurface); ---
    match &the_hc.curve {
        // ProjLib_ComputeApprox.cxx L1197-1222:
        // if (CType == GeomAbs_BSplineCurve && SType == GeomAbs_Plane) — the
        // poles through the plane frame, knots / mults / weights kept.
        Curve3::BSpline(b) => Some(Curve2d::BSpline(BSplineCurve2 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b
                .control_points
                .iter()
                .map(|p| plane_projector_project(the_p, *p))
                .collect(),
            weights: b.weights.clone(),
        })),
        // ProjLib_ComputeApprox.cxx L1229-1255:
        // else if (CType == GeomAbs_BezierCurve && SType == GeomAbs_Plane).
        Curve3::Bezier(b) => Some(Curve2d::Bezier(BezierCurve2 {
            control_points: b
                .control_points
                .iter()
                .map(|p| plane_projector_project(the_p, *p))
                .collect(),
            weights: b.weights.clone(),
        })),
        // The remaining types run the ProjLib_Function / AppParCurves
        // approximation body (ProjLib_ComputeApprox.cxx) — not translated.
        // OCCT's own failure path is preserved: `Comp.Bezier().IsNull() &&
        // Comp.BSpline().IsNull()` -> `return;` (ProjLib_ProjectedCurve.cxx
        // L716-719) leaves myResult not done, which the To2d gate below
        // rejects exactly as OCCT does.
        _ => None,
    }
}

/// `PlaneProjector::Project(P)` (ProjLib_ComputeApprox.cxx L1207/L1237): the
/// pole written in the plane frame.
fn plane_projector_project(the_p: &Plane, the_pnt: DVec3) -> DVec2 {
    let v = the_pnt - the_p.origin;
    DVec2::new(v.dot(the_p.u_dir), v.dot(the_p.v_dir))
}

/// OCCT Geom2dAdaptor::MakeCurve(Adaptor2d_Curve2d) (Geom2dAdaptor.cxx
/// L33-117) — the `switch (HC.GetType())` plus the trim step.
fn geom2d_adaptor_make_curve(
    the_type: CurveType,
    the_c: &Curve2d,
    the_first: f64,
    the_last: f64,
) -> Curve2d {
    let c2d = match the_type {
        // case GeomAbs_Line: C2D = new Geom2d_Line(HC.Line());
        CurveType::Line => match the_c {
            Curve2d::Line(_) => the_c.clone(),
            _ => unreachable!(),
        },
        // case GeomAbs_Circle / Ellipse / Parabola / Hyperbola.
        CurveType::Circle => match the_c {
            Curve2d::Circle(_) => the_c.clone(),
            _ => unreachable!(),
        },
        CurveType::Ellipse => match the_c {
            Curve2d::Ellipse(_) => the_c.clone(),
            _ => unreachable!(),
        },
        CurveType::Parabola => match the_c {
            Curve2d::Parabola(_) => the_c.clone(),
            _ => unreachable!(),
        },
        CurveType::Hyperbola => match the_c {
            Curve2d::Hyperbola(_) => the_c.clone(),
            _ => unreachable!(),
        },
        // case GeomAbs_BezierCurve: C2D = HC.Bezier();
        CurveType::Bezier => match the_c {
            Curve2d::Bezier(_) => the_c.clone(),
            _ => unreachable!(),
        },
        // case GeomAbs_BSplineCurve: C2D = HC.BSpline();
        CurveType::BSpline => match the_c {
            Curve2d::BSpline(_) => the_c.clone(),
            _ => unreachable!(),
        },
        // default: throw Standard_DomainError("Geom2dAdaptor::MakeCurve,
        // OtherCurve");
        CurveType::Other => panic!("Standard_DomainError: Geom2dAdaptor::MakeCurve, OtherCurve"),
    };

    // trim the curve if necessary (cxx L97-113).
    let c_dom = Curve2dEval::default_domain(&c2d);
    if the_first != c_dom[0] || the_last != c_dom[1] {
        // if (C2D->IsPeriodic() || (HC.FirstParameter() >= C2D->FirstParameter()
        //     && HC.LastParameter() <= C2D->LastParameter()))
        let is_periodic = matches!(c2d, Curve2d::Circle(_) | Curve2d::Ellipse(_));
        if is_periodic || (the_first >= c_dom[0] && the_last <= c_dom[1]) {
            return Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(c2d),
                t_min: the_first,
                t_max: the_last,
            });
        }
        // else { tf = max(HC.FirstParameter(), C2D->FirstParameter());
        //        tl = min(HC.LastParameter(), C2D->LastParameter()); }
        let tf = the_first.max(c_dom[0]);
        let tl = the_last.min(c_dom[1]);
        return Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c2d),
            t_min: tf,
            t_max: tl,
        });
    }
    c2d
}

/// OCCT GeomAPI::To2d(C, P) (GeomAPI.cxx L37-52).
pub fn to2d(the_c: &Curve3, the_p: &Plane) -> Option<Curve2d> {
    // cxx L40: handle(GeomAdaptor_Curve) HC = new GeomAdaptor_Curve(C);
    let a_hc = GeomCurveAdaptor::new(the_c.clone());
    // cxx L41: handle(Geom_Plane) Plane = new Geom_Plane(P);
    // cxx L42: handle(GeomAdaptor_Surface) HS = new
    // GeomAdaptor_Surface(Plane);
    //
    // cxx L44: ProjLib_ProjectedCurve Proj(HS, HC);
    let a_curve = proj_lib_projected_curve_plane(&a_hc, the_p)?;
    // cxx L46: if (Proj.GetType() != GeomAbs_OffsetCurve && Proj.GetType() !=
    // GeomAbs_OtherCurve) { result = Geom2dAdaptor::MakeCurve(Proj); } — the
    // null `result` of the other cases is the None above.
    //
    // cxx L48: result = Geom2dAdaptor::MakeCurve(Proj);
    Some(geom2d_adaptor_make_curve(
        curve2d_adaptor_type(&a_curve),
        &a_curve,
        a_hc.first,
        a_hc.last,
    ))
}
