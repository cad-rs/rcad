//! OCCT ChFiKPart — fillet particular cases:
//!   - ChFiKPart_ComputeData_Sphere.cxx L47-207 (spherical corner)
//!   - ChFiKPart_ComputeData_Rotule.cxx L43-172 (toric corner plane/plane)
//!   - ChFiKPart_ComputeData_FilPlnCyl.cxx L52-276 (gp_Lin spine overload)
//!   - ChFiKPart_ComputeData_FilPlnCyl.cxx L283-599 (gp_Circ spine overload)
//!   - ChFiKPart_ComputeData_FilPlnCon.cxx L48-341 (gp_Circ spine overload)
//!
//! 1:1 translation; the OCCT gp/ElSLib carriers come from
//! [`super::chfi_kpart_gp`].

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Line2d, Line3, Plane, Surface3};
use rcad_kernel::topo::topods::Orientation;

use super::chfi_kpart::{
    chfi_kpart_index_curve_in_ds_option, chfi_kpart_index_surface_in_ds, chfi_kpart_in_period,
    chfi_kpart_pcurve, chfi_kpart_proj_pc,
};
use super::chfi_kpart_gp::{
    curve2d_line, dir_angle, elslib_cone_d1, elslib_cone_parameters, elslib_cylinder_d1,
    elslib_cylinder_parameters, elslib_plane_d0, elslib_plane_parameters, elslib_sphere_d1,
    elslib_sphere_parameters, elslib_torus_d0, elslib_torus_d1, elslib_torus_parameters,
    gce_make_circ, surface3_ax3, surface3_d0, surface3_d1, surface3_sphere, surface3_torus,
    topabs_reverse, GpAx22d, GpAx3, GpCirc, GpConicalSurface, GpCylindricalSurface,
    CONFUSION, P_CONFUSION,
};
use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::ChFiDSSurfData;
use rcad_kernel::geom::CylindricalSurface;

// =========================================================================
// OCCT ChFiKPart_ComputeData_FilPlnCyl.cxx L52-276 — ChFiKPart_MakeFillet
// (plane/cylinder, gp_Lin spine overload).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_fillet_pln_cyl_lin(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &Line3,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    // calculate the cylinder fillet.

    // plane deviated from radius
    // OCCT L69-72: AxPln = Pln.Position(); NorPln = X ^ Y; NorF = NorPln.
    let axpln = surface3_ax3(&Surface3::Plane(*pln));
    let norpln = axpln.vxdir.cross(axpln.vydir).normalize();
    let mut norf = norpln;
    let axcyl = cyl.pos;

    if or1 == Orientation::Reversed {
        norf = -norf;
    }
    // OCCT L78: PlanOffset = Pln.Translated(Radius * gp_Vec(NorF)).
    let plan_offset = Plane {
        origin: pln.origin + norf * radius,
        ..*pln
    };

    // Parallel cylinder
    let mut roff = cyl.radius;

    if (or2 == Orientation::Forward && axcyl.direct())
        || (or2 == Orientation::Reversed && !axcyl.direct())
    {
        roff += radius;
    } else if radius < roff {
        roff -= radius;
    } else {
        // OCCT L93-96: "the fillet does not pass".
        return false;
    }
    // intersection of the parallel plane and of the parallel cylinder.
    // OCCT L99-100: gp_Cylinder CylOffset(Cyl.Position(), ROff);
    // IntAna_QuadQuadGeo LInt(PlanOffset, CylOffset, Angular, Confusion).
    let cyl_offset = CylindricalSurface {
        origin: axcyl.location,
        axis: axcyl.vzdir,
        radius: roff,
        ref_dir: axcyl.vxdir,
        y_dir: Some(axcyl.vydir),
    };
    let lint = rcad_kernel::base::int_ana::intersect_plane_cylinder_intana(&plan_offset, &cyl_offset);
    let orspine = spine.origin + spine.direction * first;
    let orfillet;
    let mut dirfillet;
    match &lint {
        rcad_kernel::base::int_ana::PlnCylResult::TwoLines(l1, l2) => {
            dirfillet = l1.direction;
            let p1 = l1.origin + l1.direction * (orspine - l1.origin).dot(l1.direction);
            // OCCT L108-119: NbSolutions() == 2 — keep the line whose point
            // is the closest to OrSpine.
            let p2 = l2.origin + l2.direction * (orspine - l2.origin).dot(l2.direction);
            if p1.distance_squared(orspine) < p2.distance_squared(orspine) {
                orfillet = p1;
            } else {
                orfillet = p2;
            }
        }
        rcad_kernel::base::int_ana::PlnCylResult::TangentLine(l1) => {
            dirfillet = l1.direction;
            orfillet = l1.origin + l1.direction * (orspine - l1.origin).dot(l1.direction);
        }
        // OCCT: the circle/ellipse solutions carry no line; LInt.Line(1)
        // raises Standard_OutOfRange there.
        _ => return false,
    }

    // Construction fillet
    if dirfillet.dot(spine.direction) < 0.0 {
        dirfillet = -dirfillet;
    }
    dirfillet = dirfillet.normalize();

    let mut uoncyl;
    let voncyl;
    let uonpln;
    let vonpln;
    let (u_cyl, v_cyl) = elslib_cylinder_parameters(&axcyl, cyl.radius, orfillet);
    uoncyl = u_cyl;
    voncyl = v_cyl;
    let tesp = CONFUSION;
    if uoncyl < fu - tesp || uoncyl > lu + tesp {
        uoncyl = chfi_kpart_in_period(uoncyl, fu, fu + 2.0 * std::f64::consts::PI, tesp);
    }
    let (u_pln, v_pln) = elslib_plane_parameters(&axpln, orfillet);
    uonpln = u_pln;
    vonpln = v_pln;

    let mut xdir;
    let mut otherdir;
    if plandab {
        xdir = -norf;
        // OCCT L149: gp_Vec(OrFillet, ElSLib::Value(UOnCyl, VOnCyl, Cyl)).
        let val = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius).0;
        otherdir = (val - orfillet).normalize();
    } else {
        otherdir = -norf;
        let val = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius).0;
        xdir = (val - orfillet).normalize();
    }
    xdir = xdir.normalize();
    otherdir = otherdir.normalize();

    // OCCT L159-164: gp_Ax3 AxFil(OrFillet, DirFillet, XDir); YReverse test.
    let mut axfil = GpAx3::new_pn_vx(orfillet, dirfillet, xdir);
    let aprod = xdir.cross(otherdir);
    if aprod.dot(dirfillet) < 0.0 {
        axfil.y_reverse();
    }

    // OCCT L166-167: Geom_CylindricalSurface(AxFil, Radius); ChangeSurf.
    let fillet = GpCylindricalSurface {
        pos: axfil,
        radius,
    };
    data.change_surf(chfi_kpart_index_surface_in_ds(fillet.to_surface3(), dstr));

    // FaceInterferences are loaded with pcurves and curves 3D.
    // edge plane-Fillet
    let ppln2d = DVec2::new(uonpln, vonpln);
    let vpln2d = DVec2::new(dirfillet.dot(axpln.vxdir), dirfillet.dot(axpln.vydir));
    let _ = ppln2d;
    let lin2dpln = (ppln2d, vpln2d);
    let ponpln = elslib_plane_d0(uonpln, vonpln, &axpln);
    let c3d = Line3 {
        origin: ponpln,
        direction: dirfillet,
    };

    let mut uonfillet;
    let v;
    let (uf, vf) =
        elslib_cylinder_parameters(&axfil, radius, ponpln);
    uonfillet = uf;
    v = vf;

    if uonfillet > std::f64::consts::PI {
        uonfillet = 0.0;
    }
    let lonfillet = (DVec2::new(uonfillet, v), DVec2::new(0.0, 1.0));
    let l3d = Curve3::Line(c3d);
    let lfac = curve2d_line(lin2dpln.0, lin2dpln.1);
    let lfil = curve2d_line(lonfillet.0, lonfillet.1);
    let (_p, deru, derv) = elslib_cylinder_d1(uonfillet, v, &axfil, radius);
    let norfil = deru.cross(derv).normalize();
    let toreverse = norfil.dot(norpln) <= 0.0;
    // It is checked if the orientation of the cylinder is the same as of the plane.
    if toreverse {
        *data.change_orientation() = topabs_reverse(ofpl);
    } else {
        *data.change_orientation() = ofpl;
    }

    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }
    if plandab {
        data.change_interference_on_s1().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    } else {
        data.change_interference_on_s2().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    }

    // edge cylinder-Fillet.
    let pcyl2d = DVec2::new(uoncyl, voncyl);
    let mut dpc = DVec2::new(0.0, 1.0);
    if dirfillet.dot(axcyl.vzdir) < 0.0 {
        dpc = -dpc;
    }
    let lin2dcyl = (pcyl2d, dpc);
    let poncyl = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius).0;
    let c3d = Line3 {
        origin: poncyl,
        direction: dirfillet,
    };
    let (uf, vf) = elslib_cylinder_parameters(&axfil, radius, poncyl);
    uonfillet = uf;
    let v = vf;
    if uonfillet > std::f64::consts::PI {
        uonfillet = 0.0;
    }
    let lonfillet = (DVec2::new(uonfillet, v), DVec2::new(0.0, 1.0));
    let l3d = Curve3::Line(c3d);
    let lfac = curve2d_line(lin2dcyl.0, lin2dcyl.1);
    let lfil = curve2d_line(lonfillet.0, lonfillet.1);

    let (_p, deru, derv) = elslib_cylinder_d1(uonfillet, v, &axfil, radius);
    let norfil = deru.cross(derv).normalize();
    let (_p, deru, derv) = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius);
    let norcyl = deru.cross(derv).normalize();

    let toreverse = norfil.dot(norcyl) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }
    if plandab {
        data.change_interference_on_s2().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    } else {
        data.change_interference_on_s1().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    }
    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_FilPlnCyl.cxx L283-599 — ChFiKPart_MakeFillet
// (plane/cylinder, gp_Circ spine overload; torus or sphere fillet).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_fillet_pln_cyl_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    // calculation of the fillet (torus or sphere).
    let mut c1sphere = false;
    let pospl = surface3_ax3(&Surface3::Plane(*pln));
    let dpnat = pospl.vxdir.cross(pospl.vydir).normalize();
    let mut dp = dpnat;
    let mut df = dp;
    if or1 == Orientation::Reversed {
        dp = -dp;
    }
    if ofpl == Orientation::Reversed {
        df = -df;
    }

    let mut or = cyl.pos.location;
    let (u0, v0) = elslib_plane_parameters(&pospl, or);
    let c2dpln = DVec2::new(u0, v0);
    // OCCT L317-318: ElSLib::PlaneD0(u, v, PosPl, Or); cPln = Or.
    let orproj = elslib_plane_d0(u0, v0, &pospl);
    let cpln = orproj;
    // OCCT ChFiKPart_ComputeData_FilPlnCyl.cxx L317-319: ElSLib::PlaneD0
    // REASSIGNS Or to its projection on the plane, then the Radius*Dp offset
    // starts from that projected point — the torus/sphere centre circle.
    or = cpln + dp * radius;
    let mut ptsp;
    let mut dsp;
    // Modification for the PtSp found at the wrong side of the sewing edge.
    let mut ptsp2;
    let mut dsp2;
    let acote = 1e-7;
    let (ps, ds) =
        super::chfi_kpart_gp::elclib_circle_d1(first, &spine.pos, spine.radius);
    ptsp = ps;
    dsp = ds.normalize();
    let (mut u, mut v) = elslib_cylinder_parameters(&cyl.pos, cyl.radius, ptsp);
    if u.abs() < acote || (u - 2.0 * std::f64::consts::PI).abs() < acote {
        let (ps2, ds2) =
            super::chfi_kpart_gp::elclib_circle_d1(first + 0.2, &spine.pos, spine.radius);
        ptsp2 = ps2;
        dsp2 = ds2;
        let (u2, _v2) = elslib_cylinder_parameters(&cyl.pos, cyl.radius, ptsp2);
        if (u2 - u).abs() > std::f64::consts::PI {
            u = 2.0 * std::f64::consts::PI - u;
            ptsp = elslib_cylinder_d1(u, v, &cyl.pos, cyl.radius).0;
            // OCCT L338-339: PR = ElCLib::Parameter(Spine, PtSp);
            // ElCLib::D1(PR, Spine, PtSp2, DSp);
            let pr = super::chfi_kpart_gp::elclib_circle_parameter(&spine.pos, ptsp);
            let (ps2, ds2) =
                super::chfi_kpart_gp::elclib_circle_d1(pr, &spine.pos, spine.radius);
            ptsp2 = ps2;
            dsp = ds2.normalize();
        }
    }
    // end of modif
    let mut dx = (ptsp - or).normalize();
    dx = dp.cross(dx.cross(dp)).normalize();
    let dy = dsp;
    let (ptcyl, vu, vv) = elslib_cylinder_d1(u, v, &cyl.pos, cyl.radius);
    let _ = ptcyl;
    let mut dc = vu.cross(vv).normalize();
    if or2 == Orientation::Reversed {
        dc = -dc;
    }
    let mut dz = dp;
    let cylrad = cyl.radius;
    let mut rad;
    let dedans = dx.dot(dc) <= 0.0;
    if dedans {
        if !plandab {
            dz = -dz;
        }
        rad = cylrad - radius;
        if rad.abs() <= CONFUSION {
            c1sphere = true;
        } else if rad < 0.0 {
            // OCCT L370-373: "the fillet can't pass".
            return false;
        }
    } else {
        if plandab {
            dz = -dz;
        }
        rad = cylrad + radius;
    }
    let mut filax3 = GpAx3::new_pn_vx(or, dz, dx);
    if filax3.vydir.dot(dy) <= 0.0 {
        filax3.y_reverse();
    }

    let filsurf;
    if c1sphere {
        // OCCT L392-393: Geom_SphericalSurface(FilAx3, Radius).
        filsurf = surface3_sphere(&filax3, radius);
        data.change_surf(chfi_kpart_index_surface_in_ds(filsurf.clone(), dstr));
    } else {
        // OCCT L397-398: Geom_ToroidalSurface(FilAx3, Rad, Radius).
        filsurf = surface3_torus(&filax3, rad, radius);
        data.change_surf(chfi_kpart_index_surface_in_ds(filsurf.clone(), dstr));
    }

    // It is checked if the orientation of the fillet is the same as of faces.
    let p;
    let mut pp = DVec3::ZERO;
    let mut deru;
    let mut derv;
    let pset = cpln + dx * rad;
    p = pset;
    let mut u;
    let mut v;
    u = 0.0;
    if (dedans && plandab) || (!dedans && !plandab) {
        if c1sphere {
            v = -std::f64::consts::FRAC_PI_2;
        } else {
            v = 3.0 * std::f64::consts::PI / 2.0;
        }
    } else {
        v = std::f64::consts::FRAC_PI_2;
    }
    let mut norfil;
    if c1sphere {
        let (_pp, du, dv) = elslib_sphere_d1(u, v, &filax3, cylrad);
        deru = du;
        derv = dv;
        pp = _pp;
        norfil = filax3.vxdir.cross(filax3.vydir).normalize();
        if v < 0.0 {
            norfil = -norfil;
        }
    } else {
        let (_pp, du, dv) = elslib_torus_d1(u, v, &filax3, rad, radius);
        deru = du;
        derv = dv;
        pp = _pp;
        norfil = deru.cross(derv).normalize();
    }
    let mut p2dfil = DVec2::new(0.0, v);
    let mut toreverse = norfil.dot(df) <= 0.0;
    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // FaceInterferences are loaded with pcurves and curves 3d.

    // The plane face.
    let gcirc2dpln: Option<rcad_kernel::geom::Curve2d>;
    let gcircpln: Option<Curve3>;
    let mut circax2 = filax3.ax2();
    if !c1sphere {
        let (up, vp) = elslib_plane_parameters(&pospl, p);
        u = up;
        v = vp;
        let p2dpln = DVec2::new(u, v);
        let d2d = DVec2::new(dsp.dot(pospl.vxdir), dsp.dot(pospl.vydir));
        let ax2dpln = GpAx22d::new(
            c2dpln,
            (p2dpln - c2dpln).normalize(),
            d2d.normalize(),
        );
        gcirc2dpln = Some(ax2dpln.to_circle2d(rad));
        circax2.set_location(cpln);
        let circpln = GpCirc::new(circax2, rad);
        gcircpln = Some(Curve3::Circle(circpln.to_circle3()));
    } else {
        let (up, vp) = elslib_plane_parameters(&pospl, p);
        u = up;
        v = vp;
        let p2dpln = DVec2::new(u, v);
        let pbid = DVec2::ZERO;
        if plandab {
            data.set_2d_points(p2dpln, p2dpln, pbid, pbid);
        } else {
            data.set_2d_points(pbid, pbid, p2dpln, p2dpln);
        }
        gcirc2dpln = None;
        gcircpln = None;
    }
    let lin2dfil = (p2dfil, DVec2::new(1.0, 0.0));
    let glin2dfil1 = curve2d_line(lin2dfil.0, lin2dfil.1);
    toreverse = norfil.dot(dpnat) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }
    if plandab {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcircpln, dstr),
            trans,
            gcirc2dpln,
            Some(glin2dfil1),
        );
    } else {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcircpln, dstr),
            trans,
            gcirc2dpln,
            Some(glin2dfil1),
        );
    }

    // The cylindrical face.
    let p = or + dx * cylrad;
    u = 0.0;
    if dedans {
        if plandab && !c1sphere {
            v = 2.0 * std::f64::consts::PI;
        } else {
            v = 0.0;
        }
    } else {
        v = std::f64::consts::PI;
    }
    p2dfil = DVec2::new(u, v);
    if c1sphere {
        let (_pp, du, dv) = elslib_sphere_d1(u, v, &filax3, cylrad);
        deru = du;
        derv = dv;
        pp = _pp;
    } else {
        let (_pp, du, dv) = elslib_torus_d1(u, v, &filax3, rad, radius);
        deru = du;
        derv = dv;
        pp = _pp;
    }
    let _ = (pp, p);
    let norfil = deru.cross(derv).normalize();
    let lin2dfil = (p2dfil, lin2dfil.1);
    let glin2dfil2 = curve2d_line(lin2dfil.0, lin2dfil.1);
    let (mut u, v) = elslib_cylinder_parameters(&cyl.pos, cyl.radius, p);
    let tol = P_CONFUSION;
    let mut careaboutsens = false;
    if (lu - fu - 2.0 * std::f64::consts::PI).abs() < tol {
        careaboutsens = true;
    }
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    if u < fu || u > lu {
        u = chfi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, tol);
    }
    let (_pp, deru, derv) = elslib_cylinder_d1(u, v, &cyl.pos, cyl.radius);
    let norcyl = deru.cross(derv).normalize();
    let mut d2dcyl = DVec2::new(1.0, 0.0);
    if deru.dot(dy) < 0.0 {
        d2dcyl = -d2dcyl;
        if careaboutsens && (fu - u).abs() < tol {
            u = lu;
        }
    } else if careaboutsens && (lu - u).abs() < tol {
        u = fu;
    }
    let p2dcyl = DVec2::new(u, v);
    let lin2dcyl = (p2dcyl, d2dcyl);
    let glin2dcyl = curve2d_line(lin2dcyl.0, lin2dcyl.1);
    circax2.set_location(or);
    let circcyl = GpCirc::new(circax2, cylrad);
    let gcirccyl = Curve3::Circle(circcyl.to_circle3());
    let toreverse = norfil.dot(norcyl) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }
    if plandab {
        data.change_interference_on_s2().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(gcirccyl, dstr),
            trans,
            Some(glin2dcyl),
            Some(glin2dfil2),
        );
    } else {
        data.change_interference_on_s1().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(gcirccyl, dstr),
            trans,
            Some(glin2dcyl),
            Some(glin2dfil2),
        );
    }
    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_FilPlnCon.cxx L48-341 — ChFiKPart_MakeFillet
// (cone/plane or plane/cone, gp_Circ spine overload; torus or sphere fillet).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_fillet_pln_con_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    con: &GpConicalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    radius: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    // calculate the fillet (torus or sphere).
    let mut c1sphere = false;
    let pospl = surface3_ax3(&Surface3::Plane(*pln));
    let dpnat = pospl.vxdir.cross(pospl.vydir).normalize();
    let mut dp = dpnat;
    let mut df = dp;
    if or1 == Orientation::Reversed {
        dp = -dp;
    }
    if ofpl == Orientation::Reversed {
        df = -df;
    }

    let mut or = con.pos.location;
    let (u0, v0) = elslib_plane_parameters(&pospl, or);
    let c2dpln = DVec2::new(u0, v0);
    // OCCT L81-82: ElSLib::PlaneD0(u, v, PosPl, Or); cPln = Or.
    let cpln = elslib_plane_d0(u0, v0, &pospl);
    // OCCT ChFiKPart_ComputeData_FilPlnCon.cxx L81-83: ElSLib::PlaneD0
    // REASSIGNS Or to its projection on the plane, then the Radius*Dp offset
    // starts from that projected point.
    or = cpln + dp * radius;

    let (ptsp, dsp) =
        super::chfi_kpart_gp::elclib_circle_d1(first, &spine.pos, spine.radius);
    let dsp = dsp.normalize();
    // OCCT L88-88: IntAna_QuadQuadGeo CInt(Pln, Con, Angular, Confusion).
    let con_payload = con.to_surface3();
    let cint = rcad_kernel::base::int_ana::intersect_plane_cone_intana(pln, match &con_payload {
        Surface3::Cone(c) => c,
        _ => unreachable!(),
    });
    let pv;
    match &cint {
        rcad_kernel::base::int_ana::PlnConResult::Circle(c) => {
            // The origin of the fillet is set at the start point on the
            // guideline.  OCCT L94: Pv = ElCLib::Value(Parameter(Circ, PtSp), Circ).
            let gc = GpCirc::from_circle3(c);
            let cpar = super::chfi_kpart_gp::elclib_circle_parameter(&gc.pos, ptsp);
            pv = gc.pos.location
                + gc.pos.vxdir * (cpar.cos() * gc.radius)
                + gc.pos.vydir * (cpar.sin() * gc.radius);
        }
        // OCCT: CInt.Circle(1) on a non-circle solution raises.
        _ => return false,
    }
    let dx = (pv - cpln).normalize();
    let dy = dsp;
    let (u0c, v0c) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, pv);
    let u = u0c;
    let v = v0c;
    let (_ptcon, vu, vv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
    let _ = _ptcon;
    let mut dc = vu.cross(vv).normalize();
    if or2 == Orientation::Reversed {
        dc = -dc;
    }
    let mut dz = dp;

    // OCCT L113-124: the plane-normal step point pp and the cone V direction.
    let pp0 = pv + dc;
    let (up, vp) = elslib_plane_parameters(&pospl, pp0);
    let pp_proj = elslib_plane_d0(up, vp, &pospl);
    let ddp = (pp_proj - pv).normalize();
    let (u2, v2) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, pv);
    let u = u2;
    let v = v2;
    let (pp_cone, dcu, dcv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
    let pp = pp_cone;
    let _ = dcu;
    let mut ddc = dcv.normalize();
    if ddc.dot(dp) < 0.0 {
        ddc = -ddc;
    }
    let ang = dir_angle(ddp, ddc);
    let rabio = radius / (ang / 2.0).tan();
    let maxrad = cpln.distance(pv);
    let mut rad;
    let dedans = dx.dot(dc) <= 0.0;
    if dedans {
        if !plandab {
            dz = -dz;
        }
        rad = maxrad - rabio;
        if rad.abs() <= CONFUSION {
            c1sphere = true;
        } else if rad < 0.0 {
            // OCCT L143-146: "the fillet does not pass".
            return false;
        }
    } else {
        if plandab {
            dz = -dz;
        }
        rad = maxrad + rabio;
    }
    let mut filax3 = GpAx3::new_pn_vx(or, dz, dx);
    if filax3.vydir.dot(dy) <= 0.0 {
        filax3.y_reverse();
    }

    let filsurf;
    if c1sphere {
        filsurf = surface3_sphere(&filax3, radius);
        data.change_surf(chfi_kpart_index_surface_in_ds(filsurf.clone(), dstr));
    } else {
        filsurf = surface3_torus(&filax3, rad, radius);
        data.change_surf(chfi_kpart_index_surface_in_ds(filsurf.clone(), dstr));
    }

    // It is checked if the orientation of the fillet is the same
    // as of the faces.
    let mut u;
    let mut v;
    let p = cpln + dx * rad;
    let mut norfil;
    let mut deru;
    let mut derv;
    let mut pp = DVec3::ZERO;
    if c1sphere {
        let (us, vs) = elslib_sphere_parameters(&filax3, rad, p);
        u = us;
        v = vs;
        let (_pp, du, dv) = elslib_sphere_d1(u, v, &filax3, rad);
        deru = du;
        derv = dv;
        pp = _pp;
    } else {
        let (ut, vt) =
            elslib_torus_parameters(&filax3, rad, radius, p);
        u = ut;
        v = vt;
        let (_pp, du, dv) = elslib_torus_d1(u, v, &filax3, rad, radius);
        deru = du;
        derv = dv;
        pp = _pp;
        if !plandab && ang < std::f64::consts::FRAC_PI_2 && dedans {
            v = v + 2.0 * std::f64::consts::PI;
        }
    }
    let mut p2dfil = DVec2::new(0.0, v);
    norfil = deru.cross(derv).normalize();
    let mut toreverse = norfil.dot(df) <= 0.0;
    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // FaceInterferences are loaded with pcurves and curves 3d.

    // The plane face.
    let gcirc2dpln: Option<rcad_kernel::geom::Curve2d>;
    let gcircpln: Option<Curve3>;
    let mut circax2 = filax3.ax2();
    if !c1sphere {
        let (up, vp) = elslib_plane_parameters(&pospl, p);
        u = up;
        v = vp;
        let p2dpln = DVec2::new(u, v);
        let d2d = DVec2::new(dsp.dot(pospl.vxdir), dsp.dot(pospl.vydir));
        let ax2dpln = GpAx22d::new(
            c2dpln,
            (p2dpln - c2dpln).normalize(),
            d2d.normalize(),
        );
        gcirc2dpln = Some(ax2dpln.to_circle2d(rad));
        circax2.set_location(cpln);
        let circpln = GpCirc::new(circax2, rad);
        gcircpln = Some(Curve3::Circle(circpln.to_circle3()));
    } else {
        gcirc2dpln = None;
        gcircpln = None;
    }
    let lin2dfil = (p2dfil, DVec2::new(1.0, 0.0));
    let glin2dfil1 = curve2d_line(lin2dfil.0, lin2dfil.1);
    toreverse = norfil.dot(dpnat) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }
    if plandab {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcircpln, dstr),
            trans,
            gcirc2dpln,
            Some(glin2dfil1),
        );
    } else {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcircpln, dstr),
            trans,
            gcirc2dpln,
            Some(glin2dfil1),
        );
    }

    // The conic face.
    let p = pv + ddc * rabio;
    if c1sphere {
        let (us, vs) = elslib_sphere_parameters(&filax3, radius, p);
        u = us;
        v = vs;
        let (_pp, du, dv) = elslib_sphere_d1(u, v, &filax3, radius);
        deru = du;
        derv = dv;
        pp = _pp;
    } else {
        let (ut, vt) =
            elslib_torus_parameters(&filax3, rad, radius, p);
        u = ut;
        v = vt;
        let (_pp, du, dv) = elslib_torus_d1(u, v, &filax3, rad, radius);
        deru = du;
        derv = dv;
        pp = _pp;
        if plandab && ang < std::f64::consts::FRAC_PI_2 && dedans {
            v = v + 2.0 * std::f64::consts::PI;
        }
    }
    let _ = pp;
    let norfil = deru.cross(derv).normalize();
    p2dfil = DVec2::new(0.0, v);
    let lin2dfil = (p2dfil, lin2dfil.1);
    let glin2dfil2 = curve2d_line(lin2dfil.0, lin2dfil.1);
    let (mut u, v) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, p);
    let tol = P_CONFUSION;
    let mut careaboutsens = false;
    if (lu - fu - 2.0 * std::f64::consts::PI).abs() < tol {
        careaboutsens = true;
    }
    if u >= fu - tol && u < fu {
        u = fu;
    }
    if u <= lu + tol && u > lu {
        u = lu;
    }
    if u < fu || u > lu {
        // OCCT L292: ElCLib::InPeriod(u, fu, fu + 2*PI) (default Eps).
        u = chfi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, P_CONFUSION);
    }
    let (_pp, deru, derv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
    let norcon = deru.cross(derv).normalize();
    let mut d2dcon = DVec2::new(1.0, 0.0);
    if deru.dot(dy) < 0.0 {
        d2dcon = -d2dcon;
        if careaboutsens && (fu - u).abs() < tol {
            u = lu;
        }
    } else if careaboutsens && (lu - u).abs() < tol {
        u = fu;
    }
    let p2dcon = DVec2::new(u, v);
    let lin2dcon = (p2dcon, d2dcon);
    let glin2dcon = curve2d_line(lin2dcon.0, lin2dcon.1);
    let scal = dp.dot(pv - p);
    let ppc = cpln + dp * scal;
    circax2.set_location(ppc);
    let circcon = GpCirc::new(circax2, p.distance(ppc));
    let gcirccon = Curve3::Circle(circcon.to_circle3());
    let toreverse = norfil.dot(norcon) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }
    if plandab {
        data.change_interference_on_s2().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(gcirccon, dstr),
            trans,
            Some(glin2dcon),
            Some(glin2dfil2),
        );
    } else {
        data.change_interference_on_s1().set_interference(
            super::chfi_kpart::chfi_kpart_index_curve_in_ds(gcirccon, dstr),
            trans,
            Some(glin2dcon),
            Some(glin2dfil2),
        );
    }
    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_Sphere.cxx L47-207 — ChFiKPart_Sphere (the
// spherical corner whose contours are not all isos, from three vertices).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_sphere(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    s1: &Surface3,
    s2: &Surface3,
    orface1: Orientation,
    _orface2: Orientation,
    or1: Orientation,
    _or2: Orientation,
    rad: f64,
    ps1: DVec2,
    p1s2: DVec2,
    p2s2: DVec2,
) -> bool {
    // Construction of the sphere :
    // - pole south on PS1
    // - origine of u given by P1S2
    // - u+ to P2S2

    let ptol = CONFUSION;
    let (p1, v1, v2) = surface3_d1(s1, ps1.x, ps1.y);
    let mut ds1 = v1.cross(v2).normalize();
    let dnat1 = ds1;
    let mut df1 = ds1;
    if or1 == Orientation::Reversed {
        ds1 = -ds1;
    }
    if orface1 == Orientation::Reversed {
        df1 = -df1;
    }
    let p2 = surface3_d0(s2, p1s2.x, p1s2.y);
    let p3 = surface3_d0(s2, p2s2.x, p2s2.y);
    let ci = gce_make_circ(p1, p2, p3).expect("StdFail_NotDone: gce_MakeCirc");
    let di = ci.axis_direction();
    let pp = ci.location();
    let rr = ci.radius();
    let delta = (rad * rad - rr * rr).sqrt();
    let mut cen = pp + di * delta;
    let mut dz = (cen - p1).normalize();
    if (ds1.dot(dz) - 1.0).abs() > ptol {
        cen = pp - di * delta;
        dz = (cen - p1).normalize();
        if (ds1.dot(dz) - 1.0).abs() > ptol {
            // OCCT L95-98: "center of the spherical corner not found".
            return false;
        }
    }
    let ddx = (p2 - cen).normalize();
    let mut dddx = ddx;
    let ddy = (p3 - cen).normalize();
    let dx = dz.cross(ddx.cross(dz)).normalize();
    let mut ddz = -dz;
    let mut filax3 = GpAx3::new_pn_vx(cen, dz, dx);
    if filax3.vydir.dot(ddy) <= 0.0 {
        filax3.y_reverse();
        ddz = -ddz;
        dddx = -dddx;
    }
    // OCCT L113-114: Geom_SphericalSurface(FilAx3, Rad); ChangeSurf.
    let gsph = surface3_sphere(&filax3, rad);
    data.change_surf(chfi_kpart_index_surface_in_ds(gsph.clone(), dstr));

    // the normal of the sphere is compared to the normal of the face
    // oriented to determine the final orientation of the fillet.
    let mut toreverse = ddz.dot(df1) <= 0.0;
    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // Parameters of p2 and p3 are calculated on the Sphere to have
    // ranges of curves.
    let mut uu1;
    let vv1;
    let uu2;
    let _vv2;
    let (u1, w1) = elslib_sphere_parameters(&filax3, rad, p2);
    uu1 = u1;
    vv1 = w1;
    uu1 = 0.0;
    let (u2, w2) = elslib_sphere_parameters(&filax3, rad, p3);
    uu2 = u2;
    _vv2 = w2;

    // FaceInterferences are loaded with pcurves and curves 3d.

    // Pointed side.

    let c: Option<Curve3> = None;
    let c2d: Option<rcad_kernel::geom::Curve2d> = None;
    let p2dfil = DVec2::new(0.0, -std::f64::consts::FRAC_PI_2);
    // gp::DX2d()
    let c2dfil0 = curve2d_line(p2dfil, DVec2::new(1.0, 0.0));
    toreverse = ddz.dot(dnat1) <= 0.0;
    let mut trans = Orientation::Reversed;
    if toreverse {
        trans = Orientation::Forward;
    }
    let mut c2dfil = c2dfil0;
    data.change_interference_on_s1().set_interference(
        chfi_kpart_index_curve_in_ds_option(c, dstr),
        trans,
        c2d,
        Some(c2dfil.clone()),
    );

    // The other side.

    let ang = dir_angle(ddx, ddy);
    let dci = ddx.cross(ddy).normalize();
    let axci = GpAx3::new_pn_vx(cen, dci, ddx);
    let ci2 = GpCirc::new(axci, rad);
    let c = Curve3::Circle(ci2.to_circle3());
    // OCCT L162-164: GeomAdaptor_Surface AS(gsph); GeomAdaptor_Curve AC(C, 0., ang);
    // ChFiKPart_ProjPC(AC, AS, C2dFil) — the pcurve of the circle on the
    // sphere support surface via ProjLib_ProjectedCurve.
    let as_ = rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor::new(
        gsph.clone(),
    );
    let ac = rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor::with_range(
        c.clone(),
        0.0,
        ang,
    );
    chfi_kpart_proj_pc(&ac, &as_, &mut c2dfil);
    // OCCT L165-166: p2dbid = C2dFil->Value(0.); pp2dbid = (uu1, vv1).
    let p2dbid = curve2d_value_zero(&c2dfil);
    let pp2dbid = DVec2::new(uu1, vv1);
    if !pp2dbid.abs_diff_eq(p2dbid, ptol) {
        // OCCT L169-170: C2dFil->Translate(v2dbid).
        c2dfil = translate_curve2d(&c2dfil, pp2dbid - p2dbid);
        // OCCT: C2dFil is a shared handle, so the Translate above is visible
        // through the S1 interference registered earlier; rcad copies values,
        // so the S1 stored pcurve is refreshed to keep the end-state equal.
        if let Some(pc) = data.change_interference_on_s1().change_pcurve_on_surf() {
            *pc = c2dfil.clone();
        }
    }
    let v2d = p2s2 - p1s2;
    let d2d = v2d.normalize();
    let c2d;
    if (v2d.length() - ang).abs() <= ptol {
        c2d = Some(curve2d_line(p1s2, d2d));
    } else {
        c2d = Some(chfi_kpart_pcurve(p1s2, p2s2, 0.0, ang));
    }
    let (pp1, v1, v2) = surface3_d1(s2, p1s2.x, p1s2.y);
    let _ = pp1;
    let ds2 = v1.cross(v2).normalize();
    toreverse = ds2.dot(dddx) <= 0.0;
    trans = Orientation::Reversed;
    if !toreverse {
        trans = Orientation::Forward;
    }
    data.change_interference_on_s2().set_interference(
        super::chfi_kpart::chfi_kpart_index_curve_in_ds(c, dstr),
        trans,
        c2d,
        Some(c2dfil),
    );

    data.change_vertex_first_on_s1().set_point(p1);
    data.change_vertex_last_on_s1().set_point(p1);
    data.change_vertex_first_on_s2().set_point(p2);
    data.change_vertex_last_on_s2().set_point(p3);
    data.change_interference_on_s1().set_first_parameter(0.0);
    data.change_interference_on_s1().set_last_parameter(uu2);
    data.change_interference_on_s2().set_first_parameter(0.0);
    data.change_interference_on_s2().set_last_parameter(ang);

    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_Rotule.cxx L43-172 — ChFiKPart_MakeRotule (the
// toric corner joint between three planes).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_rotule(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pl: &Plane,
    pl1: &Plane,
    pl2: &Plane,
    opl: Orientation,
    opl1: Orientation,
    opl2: Orientation,
    r: f64,
    ofpl: Orientation,
) -> bool {
    // calcul du tore.
    //---------------
    let pos = surface3_ax3(&Surface3::Plane(*pl));
    let mut dpl = pos.vxdir.cross(pos.vydir).normalize();
    let dfpl;
    let dplnat = dpl;
    if opl == Orientation::Reversed {
        dpl = -dpl;
    }
    if ofpl == Orientation::Reversed {
        dfpl = -dpl;
    } else {
        dfpl = dpl;
    }
    let pos1 = surface3_ax3(&Surface3::Plane(*pl1));
    let mut dpl1 = pos1.vxdir.cross(pos1.vydir).normalize();
    if opl1 == Orientation::Reversed {
        dpl1 = -dpl1;
    }
    let pos2 = surface3_ax3(&Surface3::Plane(*pl2));
    let mut dpl2 = pos2.vxdir.cross(pos2.vydir).normalize();
    if opl2 == Orientation::Reversed {
        dpl2 = -dpl2;
    }

    let alpha = dir_angle(dpl1, dpl2);

    // OCCT L84: IntAna_QuadQuadGeo LInt(pl1, pl2, Angular, Confusion).
    let lint = rcad_kernel::base::int_ana::intersect_plane_plane_intana(pl1, pl2);
    let ptor;
    let pcirc;
    match &lint {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(lint_line) => {
            let par = (pl.origin - lint_line.origin).dot(lint_line.direction);
            pcirc = lint_line.origin + lint_line.direction * par;
            ptor = pcirc + dpl * r;
        }
        _ => {
            return false;
        }
    }

    let mut ppos = GpAx3::new_pn_vx(ptor, -dpl, dpl1);
    if ppos.vydir.dot(dpl2) < 0.0 {
        ppos.y_reverse();
    }
    // OCCT L102-103: Geom_ToroidalSurface(ppos, r, r); ChangeSurf.
    let gtor = surface3_torus(&ppos, r, r);
    data.change_surf(chfi_kpart_index_surface_in_ds(gtor, dstr));

    // on compare l orientation du tore a celle de la face en bout.
    //------------------------------------------------------------
    let (_pp, du, dv) = elslib_torus_d1(
        0.0,
        std::f64::consts::FRAC_PI_2,
        &ppos,
        r,
        r,
    );
    let mut pp = _pp;
    let drot = du.cross(dv).normalize();
    let reversecur = drot.dot(dplnat) <= 0.0;
    let reversefil = drot.dot(dfpl) <= 0.0;
    if reversefil {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // on charge les FaceInterferences avec les pcurves et courbes 3d.
    //-----------------------------------------------------------------

    // du cote du plan
    //---------------
    let mut circax2 = ppos.ax2();
    circax2.set_location(pcirc);
    let gc = GpCirc::new(circax2, r);
    let (u, v) = elslib_plane_parameters(&pos, pcirc);
    let p2dcirc = DVec2::new(u, v);
    let dx2d = DVec2::new(dpl1.dot(pl.u_dir), dpl1.dot(pl.v_dir));
    let dy2d = DVec2::new(ppos.vydir.dot(pl.u_dir), ppos.vydir.dot(pl.v_dir));
    let circ2dax = GpAx22d::new(p2dcirc, dx2d.normalize(), dy2d.normalize());
    let gc2d = circ2dax.to_circle2d(r);
    let p2dlin = DVec2::new(0.0, std::f64::consts::FRAC_PI_2);
    let gl2d = curve2d_line(p2dlin, DVec2::new(1.0, 0.0));
    let mut trans = Orientation::Reversed;
    if reversecur {
        trans = Orientation::Forward;
    }
    data.change_interference_on_s1().set_interference(
        super::chfi_kpart::chfi_kpart_index_curve_in_ds(
            Curve3::Circle(gc.to_circle3()),
            dstr,
        ),
        trans,
        Some(gc2d),
        Some(gl2d),
    );

    // du cote pointu
    //--------------
    let bid: Option<Curve3> = None;
    let bid2d: Option<rcad_kernel::geom::Curve2d> = None;
    let p2dlin = DVec2::new(0.0, std::f64::consts::PI);
    let gl2dcoin = curve2d_line(p2dlin, DVec2::new(1.0, 0.0));
    data.change_interference_on_s2().set_interference(
        chfi_kpart_index_curve_in_ds_option(bid, dstr),
        trans,
        bid2d,
        Some(gl2dcoin),
    );

    // et les points
    //-------------
    data.change_vertex_first_on_s1().set_point(pp);
    pp = elslib_torus_d0(alpha, std::f64::consts::FRAC_PI_2, &ppos, r, r);
    data.change_vertex_last_on_s1().set_point(pp);
    data.change_interference_on_s1().set_first_parameter(0.0);
    data.change_interference_on_s1()
        .set_last_parameter(alpha);
    data.change_interference_on_s2().set_first_parameter(0.0);
    data.change_interference_on_s2()
        .set_last_parameter(alpha);

    true
}

// =========================================================================
// rcad-side Curve2d helpers for the OCCT Geom2d_Line operations used by the
// Sphere case (Value(0) / Translate are Geom2d_Line bodies).
// =========================================================================

/// OCCT Geom2d_Line::Value(0.) — the line origin for a straight pcurve.
/// Pending boundary: only Line pcurves are supported until ProjLib exists.
fn curve2d_value_zero(c: &rcad_kernel::geom::Curve2d) -> DVec2 {
    match c {
        rcad_kernel::geom::Curve2d::Line(l) => l.origin,
        _ => panic!("Standard_NotImplemented: pcurve value pending ProjLib boundary"),
    }
}

/// OCCT Geom2d_Line::Translate(V) — shifts the line origin.
fn translate_curve2d(c: &rcad_kernel::geom::Curve2d, v: DVec2) -> rcad_kernel::geom::Curve2d {
    match c {
        rcad_kernel::geom::Curve2d::Line(l) => {
            rcad_kernel::geom::Curve2d::Line(Line2d {
                origin: l.origin + v,
                direction: l.direction,
            })
        }
        _ => panic!("Standard_NotImplemented: pcurve translate pending ProjLib boundary"),
    }
}

/// OCCT ElCLib::Value(First, Spine) for a gp_Circ spine (line form).
#[allow(dead_code)]
fn elclib_circ_value(u: f64, spine: &GpCirc) -> DVec3 {
    spine.pos.location
        + spine.pos.vxdir * (u.cos() * spine.radius)
        + spine.pos.vydir * (u.sin() * spine.radius)
}
