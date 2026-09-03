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
use super::chfi_ds::{ChFiDSSpineHandle, ChFiDSSurfData};
use super::topopebrepds::{TopOpeBRepDSCurve, TopOpeBRepDSSurface};

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
// OCCT ChFiKPart_ComputeData_Fcts.cxx L142-155 — IndexCurveInDS /
// IndexSurfaceInDS (DStr.AddCurve / DStr.AddSurface).
// =========================================================================
pub fn chfi_kpart_index_curve_in_ds(c: Curve3, dstr: &mut TopOpeBRepDSHDataStructure) -> i32 {
    // OCCT Fcts.cxx L142-145: DStr.AddCurve(TopOpeBRepDS_Curve(C, 0.)).
    dstr.add_curve(TopOpeBRepDSCurve::new(Some(c), 0.0))
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
                let Some(line) = spine_line(sp, iedge) else {
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
            } else if (typ1 == "Plane" && typ2 == "Cylinder")
                || (typ1 == "Cylinder" && typ2 == "Plane")
            {
                // OCCT L107-176: ChFiKPart_MakeFillet plane/cylinder
                // overloads (FilPlnCyl.cxx) — pending translation.
                false
            } else if (typ1 == "Plane" && typ2 == "Cone") || (typ1 == "Cone" && typ2 == "Plane") {
                // OCCT L177-208: FilPlnCon.cxx — pending translation.
                false
            } else {
                // OCCT L209-212: throw Standard_NotImplemented.
                panic!("Standard_NotImplemented: particular case not written");
            }
        }
        ChFiDSSpineHandle::Chamf(csp) => {
            // OCCT L214-...: the chamfer branch dispatches on IsChamfer()
            // and Mode() into ChPlnPln / ChPlnCyl / ChPlnCon / ChAsym*
            // (ChFiKPart_ComputeData_Ch*.cxx) — pending translations.
            let _ = (csp, ctyp, brep);
            false
        }
    }
}

fn face_plane(s: &Shape) -> Option<Plane> {
    match s.as_face()?.surface.as_ref()? {
        Surface3::Plane(p) => Some(p.clone()),
        _ => None,
    }
}

fn spine_line(sp: &super::chfi_ds::ChFiDSFilSpine, iedge: usize) -> Option<Line3> {
    let e = sp.base.edges(iedge).clone();
    let ed = e.as_edge()?;
    match ed.curve.as_ref()? {
        Curve3::Line(l) => Some(l.clone()),
        Curve3::Circle(_) => None,
        _ => None,
    }
}

#[allow(dead_code)]
fn spine_circle(sp: &super::chfi_ds::ChFiDSFilSpine, iedge: usize) -> Option<Circle3> {
    let e = sp.base.edges(iedge).clone();
    let ed = e.as_edge()?;
    match ed.curve.as_ref()? {
        Curve3::Circle(c) => Some(c.clone()),
        _ => None,
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
    let glinpln1 = Curve3::Line(Line3 {
        origin: linpln.0,
        direction: linpln.1,
    });
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
    let glinpln2 = Curve3::Line(Line3 {
        origin: pv,
        direction: axis_cylinder,
    });
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
