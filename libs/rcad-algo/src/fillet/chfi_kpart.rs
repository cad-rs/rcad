//! OCCT ChFiKPart package — 1:1 translation of the analytic ("particular
//! case") fillet/chamfer surface computation.
//!
//! Sources:
//!   - ChFiKPart_ComputeData_Fcts.cxx (InPeriod, IndexCurveInDS,
//!     IndexSurfaceInDS, PCurve)
//!   - ChFiKPart_ComputeData.cxx L51-641 (Compute dispatch)
//!   - ChFiKPart_ComputeData_FilPlnPln.cxx L42-174 (MakeFillet plane-plane)
//!
//! Pending branches of the dispatch (Plane/Cylinder, Plane/Cone, Sphere,
//! Rotule, chamfer PlnPln/PlnCyl/PlnCon) carry their OCCT file references
//! and report failure exactly like the OCCT `return false` paths.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::topods;

use super::chfi3d::TopOpeBRepDSHDataStructure;
use super::chfi_kpart_gp::{surface3_ax3, surface3_d0, surface3_d1};
use super::chfi_ds::{ChFiDSSpineHandle, ChFiDSSurfData};
use super::chfi3d_ds::{TopOpeBRepDSCurve, TopOpeBRepDSSurface};

// =========================================================================
// OCCT ChFiKPart_ComputeData_Fcts.cxx L27-45 — ChFiKPart_InPeriod.
// =========================================================================
pub fn chfi_kpart_in_period(u: f64, ufirst: f64, ulast: f64, eps: f64) -> f64 {
    let mut u = u;
    let period = ulast - ufirst;
    while eps < ufirst - u {
        u += period;
    }
    while eps > ulast - u {
        u -= period;
    }
    if u < ufirst {
        u = ufirst;
    }
    u
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_Fcts.cxx L60-75 — ChFiKPart_PCurve.
// Calculate a straight line in form of BSpline to guarantee the parameters.
// =========================================================================
pub fn chfi_kpart_pcurve(
    uv1: DVec2,
    uv2: DVec2,
    pardeb: f64,
    parfin: f64,
) -> rcad_kernel::geom::Curve2d {
    // OCCT L65-73: Geom2d_BSplineCurve(p, k, m, 1) with poles (UV1, UV2),
    // knots (Pardeb, Parfin), multiplicities (2, 2), degree 1.
    rcad_kernel::geom::Curve2d::BSpline(rcad_kernel::geom::BSplineCurve2 {
        degree: 1,
        knots: vec![pardeb, parfin],
        control_points: vec![uv1, uv2],
        weights: vec![1.0, 1.0],
    })
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_Fcts.cxx L83-135 — ChFiKPart_ProjPC.
// For spherical corners the contours which of are not isos the circle is
// projected.
// =========================================================================
pub fn chfi_kpart_proj_pc(
    cg: &rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor,
    sg: &rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor,
    pcurv: &mut rcad_kernel::geom::Curve2d,
) {
    use rcad_kernel::base::proj_lib::proj_lib_projected_curve_b::ProjLibProjectedCurve;
    use rcad_kernel::base::proj_lib::{Adaptor3dSurface, GeomAbsSurfaceType};
    use rcad_kernel::geom::{BezierCurve2, BSplineCurve2, Curve2d};

    // OCCT L87: if (Sg.GetType() < GeomAbs_BezierSurface) — the five
    // analytic kinds precede BezierSurface in the GeomAbs_SurfaceType
    // ordering.
    let is_analytic = matches!(
        sg.get_type(),
        GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Torus
    );
    if is_analytic {
        // OCCT L88-91: HCg = new GeomAdaptor_Curve(Cg);
        //              HSg = new GeomAdaptor_Surface(Sg);
        //              ProjLib_ProjectedCurve Projc(HSg, HCg);
        let hc = std::sync::Arc::new(cg.clone());
        let hs = std::sync::Arc::new(sg.clone());
        let projc = ProjLibProjectedCurve::with_surface_curve(hs, hc);
        // OCCT L92-129: the Line / Bezier / BSpline switch over the
        // projection result.
        match projc.get_type() {
            rcad_kernel::base::proj_lib::CurveType::Line => {
                // OCCT L95-97: Pcurv = new Geom2d_Line(Projc.Line()).
                *pcurv = Curve2d::Line(projc.line());
            }
            rcad_kernel::base::proj_lib::CurveType::Bezier => {
                // OCCT L99-111: the poles (+ the weights when rational)
                // rebuild the Geom2d_BezierCurve.
                if let Some(bez) = projc.bezier() {
                    *pcurv = Curve2d::Bezier(BezierCurve2 {
                        control_points: bez.control_points,
                        weights: bez.weights,
                    });
                } else {
                    panic!("Standard_NotImplemented: failed approximation of the pcurve ");
                }
            }
            rcad_kernel::base::proj_lib::CurveType::BSpline => {
                // OCCT L113-129: the poles / knots / multiplicities (+ the
                // weights when rational) rebuild the Geom2d_BSplineCurve.
                // The rcad BSplineCurve2 carries the same data with the knot
                // vector flat-expanded.
                if let Some(b) = projc.bspline() {
                    *pcurv = Curve2d::BSpline(BSplineCurve2 {
                        degree: b.degree,
                        knots: b.knots,
                        control_points: b.control_points,
                        weights: b.weights,
                    });
                } else {
                    panic!("Standard_NotImplemented: failed approximation of the pcurve ");
                }
            }
            // OCCT L130-133: default: throw
            // Standard_NotImplemented("failed approximation of the pcurve ").
            _ => panic!("Standard_NotImplemented: failed approximation of the pcurve "),
        }
    } else {
        // OCCT L134-137: throw
        // Standard_NotImplemented("approximate pcurve on the left surface").
        panic!("Standard_NotImplemented: approximate pcurve on the left surface");
    }
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_Fcts.cxx L142-155 — IndexCurveInDS /
// IndexSurfaceInDS (DStr.AddCurve / DStr.AddSurface).
// =========================================================================
pub fn chfi_kpart_index_curve_in_ds(c: Curve3, dstr: &mut TopOpeBRepDSHDataStructure) -> i32 {
    // OCCT Fcts.cxx L142-145: DStr.AddCurve(TopOpeBRepDS_Curve(C, 0.)).
    dstr.add_curve(TopOpeBRepDSCurve::new(Some(c), 0.0))
}

/// OCCT Fcts.cxx L142-145 with a null Geom_Curve handle (the Rotule /
/// Sphere pointed-side interferences index a null curve).
pub fn chfi_kpart_index_curve_in_ds_option(
    c: Option<Curve3>,
    dstr: &mut TopOpeBRepDSHDataStructure,
) -> i32 {
    dstr.add_curve(TopOpeBRepDSCurve::new(c, 0.0))
}

pub fn chfi_kpart_index_surface_in_ds(s: Surface3, dstr: &mut TopOpeBRepDSHDataStructure) -> i32 {
    // OCCT Fcts.cxx L152-155: DStr.AddSurface(TopOpeBRepDS_Surface(S, 0.)).
    dstr.add_surface(TopOpeBRepDSSurface::new(s, 0.0))
}

// =========================================================================
// OCCT ElCLib / ElSLib analytic kernels used by FilPlnPln.
// =========================================================================

/// OCCT ElCLib::Value(U, L) — point at parameter U on a line.
pub fn elclib_line_value(u: f64, line: &Line3) -> DVec3 {
    line.origin + line.direction * u
}

/// OCCT ElCLib::Parameter(L, P) — parameter of the projection of P on a
/// line.
pub fn elclib_line_parameter(line: &Line3, p: DVec3) -> f64 {
    (p - line.origin).dot(line.direction)
}

/// OCCT ElSLib::PlaneParameters(Pos, P, u, v) — UV of P in the plane frame.
/// rcad Plane carries origin+normal; the gp_Ax3 X/Y directions are derived
/// deterministically (x = any perpendicular of the normal, y = n ^ x).
pub fn elslib_plane_parameters(plane: &Plane, p: DVec3) -> DVec2 {
    let xdir = plane.u_dir.normalize();
    let ydir = plane.v_dir.normalize();
    let d = p - plane.origin;
    DVec2::new(d.dot(xdir), d.dot(ydir))
}

/// OCCT ElSLib::CylinderD1(u, v, Pos, R, P, du, dv) — point and first
/// derivatives on a cylinder with the frame (x, y = axis ^ x, axis).
pub fn elslib_cylinder_d1(
    u: f64,
    v: f64,
    origin: DVec3,
    xdir: DVec3,
    axis: DVec3,
    radius: f64,
) -> (DVec3, DVec3, DVec3) {
    let ydir = axis.cross(xdir).normalize();
    let p = origin + (xdir * (radius * u.cos())) + (ydir * (radius * u.sin())) + (axis * v);
    let du = (xdir * (-u.sin())) + (ydir * u.cos());
    (p, du * radius, dv_of(axis))
}

fn dv_of(axis: DVec3) -> DVec3 {
    axis
}

// =========================================================================
// OCCT ChFiKPart_ComputeData.cxx L51-641 — ChFiKPart_ComputeData::Compute.
//
// The fillet branch is translated for the analytic combinations; the
// Plane/Cylinder, Plane/Cone, Sphere and Rotule combinations
// (ChFiKPart_ComputeData_FilPlnCyl.cxx L42-599, FilPlnCon.cxx, Sphere.cxx,
// Rotule.cxx) and the chamfer branch (ChPlnPln/ChPlnCyl/ChPlnCon/
// ChAsymPln*) are pending translations and report the OCCT failure path.
// =========================================================================
pub fn compute_data_compute(
    brep: &topods::BRep,
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    s1: &Shape,
    s2: &Shape,
    or1: Orientation,
    or2: Orientation,
    spine: &ChFiDSSpineHandle,
    iedge: usize,
) -> bool {
    let wref = 0.0f64;
    let _ = wref;
    let _ = brep;

    let surf_type = |s: &Shape| -> Option<&'static str> {
        let fd = s.as_face()?;
        let surf = fd.surface.as_ref()?;
        Some(match surf {
            Surface3::Plane(_) => "Plane",
            Surface3::Cylinder(_) => "Cylinder",
            Surface3::Cone(_) => "Cone",
            Surface3::Sphere(_) => "Sphere",
            Surface3::Torus(_) => "Torus",
            _ => "Other",
        })
    };
    let typ1 = surf_type(s1).unwrap_or("Other");
    let typ2 = surf_type(s2).unwrap_or("Other");

    // OCCT: the elementary spine curve type (Line or Circle).
    let ctyp: Option<&'static str> = {
        let base = spine.base();
        let e = base.edges(iedge);
        e.as_edge()
            .and_then(|ed| ed.curve.as_ref())
            .map(|c| match c {
                Curve3::Line(_) => "Line",
                Curve3::Circle(_) => "Circle",
                _ => "Other",
            })
    };
    let Some(ctyp) = ctyp else {
        return false;
    };

    // Return orientations.
    let or_face1 = s1.orientation;
    let or_face2 = s2.orientation;

    match spine {
        ChFiDSSpineHandle::Fil(sp) => {
            let radius = sp.radius_on(iedge);

            if typ1 == "Plane" && typ2 == "Plane" {
                // OCCT: ChFiKPart_MakeFillet(DStr, Data, S1->Plane(),
                // S2->Plane(), Or1, Or2, Radius, Spine->Line(), Wref,
                // OrFace1) — the gp_Lin overload.
                let Some(pl1) = face_plane(s1) else {
                    return false;
                };
                let Some(pl2) = face_plane(s2) else {
                    return false;
                };
                let Some(line) = spine_line(&sp.base, iedge) else {
                    return false;
                };
                make_fillet_plane_plane_lin(
                    dstr,
                    data,
                    &pl1,
                    &pl2,
                    or1,
                    or2,
                    radius,
                    &line,
                    wref,
                    or_face1,
                )
            } else if typ1 == "Plane" && typ2 == "Cylinder" {
                // OCCT L107-141: ChFiKPart_MakeFillet plane/cylinder — the
                // elementary spine type selects the gp_Lin / gp_Circ overload.
                let Some(pl1) = face_plane(s1) else {
                    return false;
                };
                let Some(cyl2) = face_cylinder(s2) else {
                    return false;
                };
                let (fu, lu) = face_u_range(s2);
                if ctyp == "Line" {
                    let Some(line) = spine_line(&sp.base, iedge) else {
                        return false;
                    };
                    super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_cyl_lin(
                        dstr,
                        data,
                        &pl1,
                        &cyl2,
                        fu,
                        lu,
                        or1,
                        or2,
                        radius,
                        &line,
                        wref,
                        or_face1,
                        true,
                    )
                } else {
                    let Some(circ) = spine_circle_gp(&sp.base, iedge) else {
                        return false;
                    };
                    super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_cyl_circ(
                        dstr,
                        data,
                        &pl1,
                        &cyl2,
                        fu,
                        lu,
                        or1,
                        or2,
                        radius,
                        &circ,
                        wref,
                        or_face1,
                        true,
                    )
                }
            } else if typ1 == "Cylinder" && typ2 == "Plane" {
                // OCCT L142-176: the swapped plane/cylinder case
                // (S2 plane first, plandab = false, OrFace2).
                let Some(pl2) = face_plane(s2) else {
                    return false;
                };
                let Some(cyl1) = face_cylinder(s1) else {
                    return false;
                };
                let (fu, lu) = face_u_range(s1);
                if ctyp == "Line" {
                    let Some(line) = spine_line(&sp.base, iedge) else {
                        return false;
                    };
                    super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_cyl_lin(
                        dstr,
                        data,
                        &pl2,
                        &cyl1,
                        fu,
                        lu,
                        or2,
                        or1,
                        radius,
                        &line,
                        wref,
                        or_face2,
                        false,
                    )
                } else {
                    let Some(circ) = spine_circle_gp(&sp.base, iedge) else {
                        return false;
                    };
                    super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_cyl_circ(
                        dstr,
                        data,
                        &pl2,
                        &cyl1,
                        fu,
                        lu,
                        or2,
                        or1,
                        radius,
                        &circ,
                        wref,
                        or_face2,
                        false,
                    )
                }
            } else if typ1 == "Plane" && typ2 == "Cone" {
                // OCCT L177-192: FilPlnCon with the gp_Circ spine, plandab.
                let Some(pl1) = face_plane(s1) else {
                    return false;
                };
                let Some(con2) = face_cone(s2) else {
                    return false;
                };
                let (fu, lu) = face_u_range(s2);
                let Some(circ) = spine_circle_gp(&sp.base, iedge) else {
                    return false;
                };
                super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_con_circ(
                    dstr,
                    data,
                    &pl1,
                    &con2,
                    fu,
                    lu,
                    or1,
                    or2,
                    radius,
                    &circ,
                    wref,
                    or_face1,
                    true,
                )
            } else if typ1 == "Cone" && typ2 == "Plane" {
                // OCCT L193-208: the swapped plane/cone case.
                let Some(pl2) = face_plane(s2) else {
                    return false;
                };
                let Some(con1) = face_cone(s1) else {
                    return false;
                };
                let (fu, lu) = face_u_range(s1);
                let Some(circ) = spine_circle_gp(&sp.base, iedge) else {
                    return false;
                };
                super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_con_circ(
                    dstr,
                    data,
                    &pl2,
                    &con1,
                    fu,
                    lu,
                    or2,
                    or1,
                    radius,
                    &circ,
                    wref,
                    or_face2,
                    false,
                )
            } else {
                // OCCT L209-212: throw Standard_NotImplemented.
                panic!("Standard_NotImplemented: particular case not written");
            }
        }
        ChFiDSSpineHandle::Chamf(csp) => {
            // OCCT L214-635: the chamfer branch dispatches on IsChamfer()
            // and Mode() into ChPlnPln / ChPlnCyl / ChPlnCon / ChAsym*.
            let a_mode = csp.base.my_mode;
            if csp.is_chamfer() == super::chfi_ds::ChFiDS_ChamfMethod::Sym {
                let dis = csp.get_dist();
                chamfer_dispatch_sym(
                    dstr, data, a_mode, csp, s1, s2, iedge, typ1, typ2, ctyp, or1, or2, dis,
                    wref, or_face1, or_face2,
                )
            } else if csp.is_chamfer() == super::chfi_ds::ChFiDS_ChamfMethod::TwoDist {
                let (dis1, dis2) = csp.dists();
                chamfer_dispatch_two_dist(
                    dstr, data, a_mode, csp, s1, s2, iedge, typ1, typ2, ctyp, or1, or2,
                    dis1, dis2, wref, or_face1, or_face2,
                )
            } else {
                let (dis, angle) = csp.get_dist_angle();
                let dis_on_p = true;
                chamfer_dispatch_dist_angle(
                    dstr, data, csp, s1, s2, iedge, typ1, typ2, ctyp, or1, or2, dis, angle,
                    wref, or_face1, or_face2, dis_on_p,
                )
            }
        }
    }
}

/// OCCT ChFiDS_Spine::Line() — the current elementary spine gp_Lin.
pub(crate) fn spine_line(sp: &super::chfi_ds::ChFiDSSpine, iedge: usize) -> Option<Line3> {
    let e = sp.edges(iedge).clone();
    let ed = e.as_edge()?;
    match ed.curve.as_ref()? {
        Curve3::Line(l) => Some(l.clone()),
        Curve3::Circle(_) => None,
        _ => None,
    }
}

/// OCCT ChFiDS_Spine::Circle() — the current elementary spine gp_Circ
/// (the gp_Circ carrier is [`GpCirc`]).
pub(crate) fn spine_circle_gp(
    sp: &super::chfi_ds::ChFiDSSpine,
    iedge: usize,
) -> Option<super::chfi_kpart_gp::GpCirc> {
    let e = sp.edges(iedge).clone();
    let ed = e.as_edge()?;
    match ed.curve.as_ref()? {
        Curve3::Circle(c) => Some(super::chfi_kpart_gp::GpCirc::from_circle3(c)),
        _ => None,
    }
}

#[allow(dead_code)]
fn spine_circle(sp: &super::chfi_ds::ChFiDSSpine, iedge: usize) -> Option<Circle3> {
    let e = sp.edges(iedge).clone();
    let ed = e.as_edge()?;
    match ed.curve.as_ref()? {
        Curve3::Circle(c) => Some(c.clone()),
        _ => None,
    }
}

/// OCCT S->Plane() — the plane payload of the support face.
pub(crate) fn face_plane(s: &Shape) -> Option<Plane> {
    match s.as_face()?.surface.as_ref()? {
        Surface3::Plane(p) => Some(*p),
        _ => None,
    }
}

/// OCCT S->Cylinder() — the cylinder payload as a gp_Cylinder carrier.
pub(crate) fn face_cylinder(s: &Shape) -> Option<super::chfi_kpart_gp::GpCylindricalSurface> {
    let c = match s.as_face()?.surface.as_ref()? {
        Surface3::Cylinder(c) => c,
        _ => return None,
    };
    Some(super::chfi_kpart_gp::GpCylindricalSurface {
        pos: super::chfi_kpart_gp::GpAx3 {
            location: c.origin,
            vxdir: c.ref_dir,
            vydir: c.y_axis(),
            vzdir: c.axis,
        },
        radius: c.radius,
    })
}

/// OCCT S->Cone() — the cone payload as a gp_Cone carrier.
pub(crate) fn face_cone(s: &Shape) -> Option<super::chfi_kpart_gp::GpConicalSurface> {
    let c = match s.as_face()?.surface.as_ref()? {
        Surface3::Cone(c) => *c,
        _ => return None,
    };
    Some(super::chfi_kpart_gp::GpConicalSurface {
        pos: super::chfi_kpart_gp::GpAx3 {
            location: c.apex,
            vxdir: c.ref_dir,
            vydir: c.axis.cross(c.ref_dir).normalize(),
            vzdir: c.axis,
        },
        semi_angle: c.half_angle_rad,
        ref_radius: c.radius,
    })
}

/// OCCT ChFiKPart_ComputeData.cxx L115-116 (and its L131/L150/L165
/// siblings) — `S2->FirstUParameter(), S2->LastUParameter()` of the
/// BRepAdaptor_Surface built by ChFi3d_Builder::ConexFaces
/// (ChFi3d_Builder_2.cxx L884: `Sb.Initialize(F)` with the default
/// Restriction = true, BRepAdaptor_Surface.hxx L65).
/// BRepAdaptor_Surface::Initialize (BRepAdaptor_Surface.cxx L71-76) loads
/// `BRepTools::UVBounds(F, umin, umax, vmin, vmax)` as the adaptor domain,
/// so the first/last U parameters are the FACE's parameter-space box.  The
/// rcad TFaceData::uv_domain cache carries that box; when it is absent the
/// OCCT AddUVBounds empty-box fallback applies — the surface natural
/// domain (BRepTools.cxx L139-153), and a null surface leaves the box void
/// (L76-79).  The former rcad (0, 2*pi) constant is gone.
pub(crate) fn face_u_range(s: &Shape) -> (f64, f64) {
    if let Some([umin, umax, ..]) = s.as_face().and_then(|f| f.uv_domain) {
        (umin, umax)
    } else {
        // OCCT BRepTools::AddUVBounds L139-153: the empty box takes the
        // surface natural bounds; a null surface leaves it void (zeros).
        match s.as_face().and_then(|f| f.surface.as_ref()) {
            Some(surf) => {
                let [u1, u2, _, _] = rcad_kernel::geom::SurfaceEval::default_domain(surf);
                (u1, u2)
            }
            None => (0.0, 0.0),
        }
    }
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_FilPlnPln.cxx L42-174 — ChFiKPart_MakeFillet
// (plane-plane, gp_Lin spine overload).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn make_fillet_plane_plane_lin(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pl1: &Plane,
    pl2: &Plane,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &Line3,
    first: f64,
    of1: Orientation,
) -> bool {
    // calcul du cylindre
    // OCCT: D1 = Pos1.XDirection().Crossed(Pos1.YDirection()) — the plane
    // normal; rcad Plane stores the normal directly.
    let mut d1 = pl1.normal.normalize();
    if or1 == Orientation::Reversed {
        d1 = -d1;
    }
    let mut d2 = pl2.normal.normalize();
    if or2 == Orientation::Reversed {
        d2 = -d2;
    }

    // OCCT: IntAna_QuadQuadGeo LInt(Pl1, Pl2, Angular, Confusion).
    let lint = rcad_kernel::base::int_ana::intersect_plane_plane_intana(pl1, pl2);
    let pv;
    match &lint {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(lint_line) => {
            // On met l origine du cylindre au point de depart fourni sur la
            // ligne guide: ElCLib::Value(Parameter(LIntLine, ElCLib::Value(First, Spine)), LIntLine).
            let p0 = elclib_line_value(first, spine);
            let par = elclib_line_parameter(lint_line, p0);
            pv = elclib_line_value(par, lint_line);
        }
        _ => return false,
    }

    let axis_cylinder = spine.direction.normalize();
    let ang = {
        let dot = d1.dot(d2).clamp(-1.0, 1.0);
        dot.acos()
    };
    let v = d1 + d2;
    let sdir = v.normalize();
    let fac = radius / (ang / 2.0).cos();
    let c = pv + sdir * fac;
    let xdir = -d1;
    // OCCT: gp_Ax3 CylAx3(C, AxisCylinder, xdir); if (YDirection().Dot(D2) >= 0) YReverse.
    let mut ydir = axis_cylinder.cross(xdir).normalize();
    if ydir.dot(d2) >= 0.0 {
        ydir = -ydir;
    }
    let gcyl = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
        origin: c,
        axis: axis_cylinder,
        radius,
        ref_dir: xdir,
        y_dir: Some(ydir),
    });
    let surf_index = chfi_kpart_index_surface_in_ds(gcyl.clone(), dstr);
    data.change_surf(surf_index);

    // On regarde si l orientation du cylindre est la meme que celle des faces.
    let (_p, deru, derv) = elslib_cylinder_d1(0.0, 0.0, c, xdir, axis_cylinder, radius);
    let norcyl = deru.cross(derv).normalize();
    let norpl = pl1.normal.normalize();
    let mut norface = norpl;
    if of1 == Orientation::Reversed {
        norface = -norface;
    }
    let toreverse = norcyl.dot(norface) <= 0.0;
    *data.change_orientation() = if toreverse {
        Orientation::Reversed
    } else {
        Orientation::Forward
    };

    // On charge les FaceInterferences avec les pcurves et courbes 3d.
    // La face 1.
    let mut p2dpln = elslib_plane_parameters(pl1, pv);
    let dir2dpln = DVec2::new(axis_cylinder.dot(xdir_of(pl1)), axis_cylinder.dot(ydir_of(pl1)));
    let mut lin2dpln = (p2dpln, dir2dpln);
    let linpln = (pv, axis_cylinder);
    let lin2dcyl = DVec2::new(0.0, 0.0);
    let trans;
    let mut toreverse2 = norcyl.dot(norpl) <= 0.0;
    if toreverse2 {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }
    let glin2dpln1 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: lin2dpln.0,
        direction: lin2dpln.1,
    });
    let glinpln1 = Curve3::Line(Line3::new(linpln.0, linpln.1));
    let glin2dcyl1 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: lin2dcyl,
        direction: DVec2::new(0.0, 1.0),
    });
    data.change_interference_on_s1().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln1, dstr),
        trans,
        Some(glin2dpln1),
        Some(glin2dcyl1),
    );

    // La face 2.
    let (_p, deru, derv) = elslib_cylinder_d1(ang, 0.0, c, xdir, axis_cylinder, radius);
    let norcyl = deru.cross(derv).normalize();
    let norpl2 = pl2.normal.normalize();
    toreverse2 = norcyl.dot(norpl2) <= 0.0;
    p2dpln = elslib_plane_parameters(pl2, pv);
    lin2dpln = (
        p2dpln,
        DVec2::new(axis_cylinder.dot(xdir_of(pl2)), axis_cylinder.dot(ydir_of(pl2))),
    );
    let glin2dpln2 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: lin2dpln.0,
        direction: lin2dpln.1,
    });
    let glinpln2 = Curve3::Line(Line3::new(pv, axis_cylinder));
    let glin2dcyl2 = rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
        origin: DVec2::new(ang, 0.0),
        direction: DVec2::new(0.0, 1.0),
    });
    let trans2 = if toreverse2 {
        Orientation::Forward
    } else {
        Orientation::Reversed
    };
    data.change_interference_on_s2().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln2, dstr),
        trans2,
        Some(glin2dpln2),
        Some(glin2dcyl2),
    );
    true
}

fn xdir_of(plane: &Plane) -> DVec3 {
    plane.u_dir.normalize()
}

fn ydir_of(plane: &Plane) -> DVec3 {
    let x = xdir_of(plane);
    plane.normal.cross(x).normalize()
}

// =========================================================================
// Chamfer dispatch helpers — the OCCT ChFiKPart_ComputeData::Compute
// chamfer branch (ChFiKPart_ComputeData.cxx L214-635) split by
// IsChamfer() value.  Each branch extracts the spine curve
// (CSpine->Line() / CSpine->Circle()) and the support payloads, then calls
// the translated case function with the OCCT argument order.
// =========================================================================

/// OCCT ChFiKPart_ComputeData.cxx L219-357 — IsChamfer() == ChFiDS_Sym.
#[allow(clippy::too_many_arguments)]
fn chamfer_dispatch_sym(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    a_mode: super::chfi_ds::ChFiDS_ChamfMode,
    csp: &super::chfi_ds::ChFiDSChamfSpine,
    s1: &Shape,
    s2: &Shape,
    iedge: usize,
    typ1: &str,
    typ2: &str,
    ctyp: &str,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    wref: f64,
    or_face1: Orientation,
    or_face2: Orientation,
) -> bool {
    if typ1 == "Plane" && typ2 == "Plane" {
        // OCCT L224-238: MakeChamfer(aMode, S1->Plane(), S2->Plane(), Or1,
        // Or2, dis, dis, CSpine->Line(), Wref, OrFace1).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(line) = spine_line(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch::chfi_kpart_make_chamfer_pln_pln(
            dstr, data, a_mode, &pl1, &pl2, or1, or2, dis, dis, &line, wref, or_face1,
        )
    } else if typ1 == "Plane" && typ2 == "Cylinder" {
        // OCCT L239-277: the plane is S1, plandab = true, OrFace1.
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(cyl2) = face_cylinder(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        chamfer_pln_cyl(
            dstr,
            data,
            a_mode,
            csp,
            iedge,
            ctyp,
            &pl1,
            &cyl2,
            fu,
            lu,
            or1,
            or2,
            dis,
            dis,
            wref,
            or_face1,
            true,
        )
    } else if typ1 == "Cylinder" && typ2 == "Plane" {
        // OCCT L278-316: the plane is S2, plandab = false, OrFace2.
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(cyl1) = face_cylinder(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        chamfer_pln_cyl(
            dstr,
            data,
            a_mode,
            csp,
            iedge,
            ctyp,
            &pl2,
            &cyl1,
            fu,
            lu,
            or2,
            or1,
            dis,
            dis,
            wref,
            or_face2,
            false,
        )
    } else if typ1 == "Plane" && typ2 == "Cone" {
        // OCCT L317-334: MakeChamfer(aMode, S1->Plane(), S2->Cone(), ...,
        // dis, dis, CSpine->Circle(), Wref, OrFace1, true).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(con2) = face_cone(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_con_circ(
            dstr, data, a_mode, &pl1, &con2, fu, lu, or1, or2, dis, dis, &circ, wref,
            or_face1, true,
        )
    } else if typ1 == "Cone" && typ2 == "Plane" {
        // OCCT L335-352: the swapped cone/plane case.
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(con1) = face_cone(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_con_circ(
            dstr, data, a_mode, &pl2, &con1, fu, lu, or2, or1, dis, dis, &circ, wref,
            or_face2, false,
        )
    } else {
        // OCCT L353-356: throw Standard_NotImplemented.
        panic!("Standard_NotImplemented: particular case not written");
    }
}

/// The plane/cylinder chamfer sub-dispatch — OCCT selects the gp_Circ or
/// gp_Lin spine overload from the elementary spine type (L241-276 with
/// ctyp == Circle / else).
#[allow(clippy::too_many_arguments)]
fn chamfer_pln_cyl(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    a_mode: super::chfi_ds::ChFiDS_ChamfMode,
    csp: &super::chfi_ds::ChFiDSChamfSpine,
    iedge: usize,
    ctyp: &str,
    pln: &Plane,
    cyl: &super::chfi_kpart_gp::GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis1: f64,
    dis2: f64,
    wref: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    if ctyp == "Circle" {
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_cyl_circ(
            dstr, data, a_mode, pln, cyl, fu, lu, or1, or2, dis1, dis2, &circ, wref, ofpl,
            plandab,
        )
    } else {
        let Some(line) = spine_line(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_cyl_lin(
            dstr, data, a_mode, pln, cyl, fu, lu, or1, or2, dis1, dis2, &line, wref, ofpl,
            plandab,
        )
    }
}

/// OCCT ChFiKPart_ComputeData.cxx L358-495 — IsChamfer() == ChFiDS_TwoDist.
#[allow(clippy::too_many_arguments)]
fn chamfer_dispatch_two_dist(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    a_mode: super::chfi_ds::ChFiDS_ChamfMode,
    csp: &super::chfi_ds::ChFiDSChamfSpine,
    s1: &Shape,
    s2: &Shape,
    iedge: usize,
    typ1: &str,
    typ2: &str,
    ctyp: &str,
    or1: Orientation,
    or2: Orientation,
    dis1: f64,
    dis2: f64,
    wref: f64,
    or_face1: Orientation,
    or_face2: Orientation,
) -> bool {
    if typ1 == "Plane" && typ2 == "Plane" {
        // OCCT L362-376: MakeChamfer(aMode, Pl1, Pl2, Or1, Or2, dis1, dis2,
        // CSpine->Line(), Wref, OrFace1).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(line) = spine_line(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch::chfi_kpart_make_chamfer_pln_pln(
            dstr, data, a_mode, &pl1, &pl2, or1, or2, dis1, dis2, &line, wref, or_face1,
        )
    } else if typ1 == "Plane" && typ2 == "Cylinder" {
        // OCCT L377-415: (dis1, dis2), plandab = true.
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(cyl2) = face_cylinder(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        chamfer_pln_cyl(
            dstr, data, a_mode, csp, iedge, ctyp, &pl1, &cyl2, fu, lu, or1, or2, dis1,
            dis2, wref, or_face1, true,
        )
    } else if typ1 == "Cylinder" && typ2 == "Plane" {
        // OCCT L416-454: the distances are swapped (dis2, dis1) and the
        // surfaces reversed; plandab = false, OrFace2.
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(cyl1) = face_cylinder(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        chamfer_pln_cyl(
            dstr, data, a_mode, csp, iedge, ctyp, &pl2, &cyl1, fu, lu, or2, or1, dis2,
            dis1, wref, or_face2, false,
        )
    } else if typ1 == "Plane" && typ2 == "Cone" {
        // OCCT L455-472: MakeChamfer(aMode, S1->Plane(), S2->Cone(), ...,
        // dis1, dis2, CSpine->Circle(), Wref, OrFace1, true).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(con2) = face_cone(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_con_circ(
            dstr, data, a_mode, &pl1, &con2, fu, lu, or1, or2, dis1, dis2, &circ, wref,
            or_face1, true,
        )
    } else if typ1 == "Cone" && typ2 == "Plane" {
        // OCCT L473-490: the swapped cone/plane case (dis1, dis2).
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(con1) = face_cone(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch_plncyl::chfi_kpart_make_chamfer_pln_con_circ(
            dstr, data, a_mode, &pl2, &con1, fu, lu, or2, or1, dis1, dis2, &circ, wref,
            or_face2, false,
        )
    } else {
        // OCCT L491-494: throw Standard_NotImplemented.
        panic!("Standard_NotImplemented: particular case not written");
    }
}

/// OCCT ChFiKPart_ComputeData.cxx L496-634 — the DistAngle (ChAsym) branch.
#[allow(clippy::too_many_arguments)]
#[allow(unused_variables)]
fn chamfer_dispatch_dist_angle(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    csp: &super::chfi_ds::ChFiDSChamfSpine,
    s1: &Shape,
    s2: &Shape,
    iedge: usize,
    typ1: &str,
    typ2: &str,
    ctyp: &str,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    wref: f64,
    or_face1: Orientation,
    or_face2: Orientation,
    dis_on_p: bool,
) -> bool {
    if typ1 == "Plane" && typ2 == "Plane" {
        // OCCT L501-515: MakeChAsym(Pl1, Pl2, Or1, Or2, dis, Angle,
        // CSpine->Line(), Wref, OrFace1, DisOnP).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(line) = spine_line(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_ch::chfi_kpart_make_ch_asym_pln_pln(
            dstr, data, &pl1, &pl2, or1, or2, dis, angle, &line, wref, or_face1, dis_on_p,
        )
    } else if typ1 == "Plane" && typ2 == "Cylinder" {
        // OCCT L516-554: plandab = true.
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(cyl2) = face_cylinder(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        ch_asym_pln_cyl(
            dstr, data, csp, iedge, ctyp, &pl1, &cyl2, fu, lu, or1, or2, dis, angle,
            wref, or_face1, true, dis_on_p,
        )
    } else if typ1 == "Cylinder" && typ2 == "Plane" {
        // OCCT L555-593: the swapped cylinder/plane case; plandab = false.
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(cyl1) = face_cylinder(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        ch_asym_pln_cyl(
            dstr, data, csp, iedge, ctyp, &pl2, &cyl1, fu, lu, or2, or1, dis, angle,
            wref, or_face2, false, dis_on_p,
        )
    } else if typ1 == "Plane" && typ2 == "Cone" {
        // OCCT L594-611: MakeChAsym(S1->Plane(), S2->Cone(), ..., Circle,
        // Wref, OrFace1, true, DisOnP).
        let Some(pl1) = face_plane(s1) else {
            return false;
        };
        let Some(con2) = face_cone(s2) else {
            return false;
        };
        let (fu, lu) = face_u_range(s2);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_chasym::chfi_kpart_make_ch_asym_pln_con_circ(
            dstr, data, &pl1, &con2, fu, lu, or1, or2, dis, angle, &circ, wref, or_face1,
            true, dis_on_p,
        )
    } else if typ1 == "Cone" && typ2 == "Plane" {
        // OCCT L612-629: the swapped cone/plane case.
        let Some(pl2) = face_plane(s2) else {
            return false;
        };
        let Some(con1) = face_cone(s1) else {
            return false;
        };
        let (fu, lu) = face_u_range(s1);
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_chasym::chfi_kpart_make_ch_asym_pln_con_circ(
            dstr, data, &pl2, &con1, fu, lu, or2, or1, dis, angle, &circ, wref, or_face2,
            false, dis_on_p,
        )
    } else {
        // OCCT L630-633: throw Standard_NotImplemented.
        panic!("Standard_NotImplemented: particular case not written");
    }
}

/// The plane/cylinder ChAsym sub-dispatch — OCCT selects the gp_Circ or
/// gp_Lin spine overload from the elementary spine type (L518-553).
#[allow(clippy::too_many_arguments)]
fn ch_asym_pln_cyl(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    csp: &super::chfi_ds::ChFiDSChamfSpine,
    iedge: usize,
    ctyp: &str,
    pln: &Plane,
    cyl: &super::chfi_kpart_gp::GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    wref: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> bool {
    if ctyp == "Circle" {
        let Some(circ) = spine_circle_gp(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_chasym::chfi_kpart_make_ch_asym_pln_cyl_circ(
            dstr, data, pln, cyl, fu, lu, or1, or2, dis, angle, &circ, wref, ofpl, plandab,
            dis_on_p,
        )
    } else {
        let Some(line) = spine_line(&csp.base, iedge) else {
            return false;
        };
        super::chfi_kpart_chasym::chfi_kpart_make_ch_asym_pln_cyl_lin(
            dstr, data, pln, cyl, fu, lu, or1, or2, dis, angle, &line, wref, ofpl, plandab,
            dis_on_p,
        )
    }
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_CS.cxx L24-73 — ChFiKPart_CornerSpine.
// The guideline is the circle corresponding to the section of S2, and other
// construction elements.
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_corner_spine(
    s1: &Surface3,
    s2: &Surface3,
    p1s1: DVec2,
    _p2s1: DVec2,
    p1s2: DVec2,
    p2s2: DVec2,
    r: f64,
) -> (
    super::chfi_kpart_gp::GpCylindricalSurface,
    super::chfi_kpart_gp::GpCirc,
    f64,
    f64,
) {
    use super::chfi_kpart_gp::{elclib_circle_parameter, GpAx3, GpCirc, GpCylindricalSurface};
    // OCCT L37-39: ax = S1->Plane().Position(); V1 = XDirection; V2 = YDirection.
    let ax = surface3_ax3(s1);
    let v1 = ax.vxdir;
    let v2 = ax.vydir;
    // OCCT L42: S2->D1(P1S2.X(), P1S2.Y(), P, du, dv).
    let (p, du, dv) = surface3_d1(s2, p1s2.x, p1s2.y);
    // OCCT L43-45: V = gp_Vec(P, S1->Value(P1S1)); V = (V.V1)*V1 + (V.V2)*V2.
    let value_p1s1 = surface3_d0(s1, p1s1.x, p1s1.y);
    let mut v = value_p1s1 - p;
    v = v1 * v.dot(v1) + v2 * v.dot(v2);
    v = v.normalize();
    // OCCT L46-47: P2 = S2->Value(P2S2); Vorien = gp_Vec(P, P2).
    let p2 = surface3_d0(s2, p2s2.x, p2s2.y);
    let vorien = p2 - p;
    let mut dx;
    let cent;
    if v.dot(vorien) >= 0.0 {
        // OCCT L50-54: cent = P + R*V; dx.Reverse().
        cent = p + v * r;
        dx = -v;
    } else {
        // OCCT L55-58: cent = P - R*V.
        cent = p - v * r;
        dx = v;
    }
    // OCCT L59-60: dy = (dx ^ dy) ^ dx (gp_Dir ops normalize).
    let dy = dx.cross(p2 - cent).normalize();
    let dy = dx.cross(dy).normalize();
    // OCCT L61-62: circax2 = gp_Ax2(cent, dx ^ dy, dx); cylax3 = gp_Ax3(circax2).
    let dxdy = dx.cross(dy).normalize();
    let mut circax2 = GpAx3::new_pn_vx(cent, dxdy, dx);
    let mut cylax3 = circax2;
    // OCCT L63-66: if ((du ^ dv).Dot(dx) < 0.) cylax3.ZReverse().
    if du.cross(dv).dot(dx) < 0.0 {
        cylax3.z_reverse();
    }
    // OCCT L67-68: First = 0; Last = ElCLib::CircleParameter(circax2, P2).
    let first = 0.0;
    let last = elclib_circle_parameter(&circax2, p2);
    let circ = GpCirc::new(circax2, r);
    let cyl = GpCylindricalSurface::new(cylax3, r);
    let _ = &mut circax2;
    (cyl, circ, first, last)
}

// =========================================================================
// OCCT ChFiKPart_ComputeData.cxx L641-712 — ComputeCorner (the toric corner
// against a cylindrical or other support S2).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn compute_data_compute_corner_cyl(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    s1: &Shape,
    s2: &Shape,
    orface1: Orientation,
    _orface2: Orientation,
    or1: Orientation,
    or2: Orientation,
    minrad: f64,
    majrad: f64,
    p1s1: DVec2,
    p2s1: DVec2,
    p1s2: DVec2,
    p2s2: DVec2,
) -> bool {
    let surf = |s: &Shape| -> Surface3 {
        s.as_face()
            .and_then(|f| f.surface.clone())
            .expect("ComputeCorner: no surface")
    };
    let surf1 = surf(s1);
    let surf2 = surf(s2);
    let typ1 = super::chfi_kpart_gp::surface3_kind(&surf1);
    let typ2 = super::chfi_kpart_gp::surface3_kind(&surf2);
    if typ1 != super::chfi_kpart_gp::SurfaceKind::Plane {
        // OCCT L659-662: throw Standard_ConstructionError.
        panic!("Standard_ConstructionError: la face du conge torique doit etre plane");
    }
    // The guideline is the circle corresponding
    // to the section of S2, and other construction elements.
    let (mut cyl, circ, first, last) = chfi_kpart_corner_spine(
        &surf1, &surf2, p1s1, p2s1, p1s2, p2s2, majrad,
    );
    let fu;
    let lu;
    if typ2 == super::chfi_kpart_gp::SurfaceKind::Cylinder {
        // OCCT L670-675: cyl = S2->Cylinder(); fu = P1S2.X(); lu = P2S2.X().
        cyl = face_cylinder(s2).expect("ComputeCorner: cylinder expected");
        fu = p1s2.x;
        lu = p2s2.x;
    } else {
        fu = first;
        lu = last;
    }
    let Some(pl1) = face_plane(s1) else {
        return false;
    };
    let surfok = super::chfi_kpart_fil::chfi_kpart_make_fillet_pln_cyl_circ(
        dstr, data, &pl1, &cyl, fu, lu, or1, or2, minrad, &circ, first, orface1, true,
    );
    if surfok {
        if typ2 != super::chfi_kpart_gp::SurfaceKind::Cylinder {
            // OCCT L698-700: ChangePCurveOnFace() = ChFiKPart_PCurve(...).
            let pcurv = chfi_kpart_pcurve(p1s2, p2s2, first, last);
            *data
                .change_interference_on_s2()
                .change_pcurve_on_face() = Some(pcurv);
        }
        // OCCT L701-704: the four vertex points.
        data.change_vertex_first_on_s1()
            .set_point(surface3_d0(&surf1, p1s1.x, p1s1.y));
        data.change_vertex_last_on_s1()
            .set_point(surface3_d0(&surf1, p2s1.x, p2s1.y));
        data.change_vertex_first_on_s2()
            .set_point(surface3_d0(&surf2, p1s2.x, p1s2.y));
        data.change_vertex_last_on_s2()
            .set_point(surface3_d0(&surf2, p2s2.x, p2s2.y));
        // OCCT L705-708: the interference parameter ranges.
        data.change_interference_on_s1().set_first_parameter(first);
        data.change_interference_on_s1().set_last_parameter(last);
        data.change_interference_on_s2().set_first_parameter(first);
        data.change_interference_on_s2().set_last_parameter(last);
        return true;
    }
    false
}

// =========================================================================
// OCCT ChFiKPart_ComputeData.cxx L716-730 — ComputeCorner (the spherical
// corner from three vertices).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn compute_data_compute_corner_sphere(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    s1: &Shape,
    s2: &Shape,
    orface1: Orientation,
    orface2: Orientation,
    or1: Orientation,
    or2: Orientation,
    rad: f64,
    ps1: DVec2,
    p1s2: DVec2,
    p2s2: DVec2,
) -> bool {
    let surf = |s: &Shape| -> Surface3 {
        s.as_face()
            .and_then(|f| f.surface.clone())
            .expect("ComputeCorner: no surface")
    };
    let surf1 = surf(s1);
    let surf2 = surf(s2);
    super::chfi_kpart_fil::chfi_kpart_sphere(
        dstr,
        data,
        &surf1,
        &surf2,
        orface1,
        orface2,
        or1,
        or2,
        rad,
        ps1,
        p1s2,
        p2s2,
    )
}

// =========================================================================
// OCCT ChFiKPart_ComputeData.cxx L734-762 — ComputeCorner (the toric joint
// between three planes, ChFiKPart_MakeRotule).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn compute_data_compute_corner_rotule(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    s: &Shape,
    s1: &Shape,
    s2: &Shape,
    ofs: Orientation,
    os: Orientation,
    os1: Orientation,
    os2: Orientation,
    radius: f64,
) -> bool {
    let surf = |sh: &Shape| -> Surface3 {
        sh.as_face()
            .and_then(|f| f.surface.clone())
            .expect("ComputeCorner: no surface")
    };
    let typ = super::chfi_kpart_gp::surface3_kind(&surf(s));
    let typ1 = super::chfi_kpart_gp::surface3_kind(&surf(s1));
    let typ2 = super::chfi_kpart_gp::surface3_kind(&surf(s2));
    if typ != super::chfi_kpart_gp::SurfaceKind::Plane
        || typ1 != super::chfi_kpart_gp::SurfaceKind::Plane
        || typ2 != super::chfi_kpart_gp::SurfaceKind::Plane
    {
        // OCCT L748-751: throw Standard_ConstructionError.
        panic!("Standard_ConstructionError: torus joint only between the planes");
    }
    let Some(pl) = face_plane(s) else {
        return false;
    };
    let Some(pl1) = face_plane(s1) else {
        return false;
    };
    let Some(pl2) = face_plane(s2) else {
        return false;
    };
    super::chfi_kpart_fil::chfi_kpart_make_rotule(
        dstr, data, &pl, &pl1, &pl2, os, os1, os2, radius, ofs,
    )
}
