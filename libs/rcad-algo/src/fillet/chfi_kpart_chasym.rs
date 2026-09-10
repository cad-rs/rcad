//! OCCT ChFiKPart — asymmetric chamfer particular cases:
//!   - ChFiKPart_ComputeData_ChAsymPlnCyl.cxx L65-369 (gp_Circ spine overload)
//!   - ChFiKPart_ComputeData_ChAsymPlnCyl.cxx L376-656 (gp_Lin spine overload)
//!   - ChFiKPart_ComputeData_ChAsymPlnCon.cxx L57-708 (gp_Circ spine overload)
//!
//! 1:1 translation; the OCCT gp/ElSLib carriers come from
//! [`super::chfi_kpart_gp`].

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topo::topods::Orientation;

use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::ChFiDSSurfData;
use super::chfi_kpart::{
    chfi_kpart_index_curve_in_ds, chfi_kpart_index_curve_in_ds_option,
    chfi_kpart_index_surface_in_ds, chfi_kpart_in_period,
};
use super::chfi_kpart_gp::{
    curve2d_line, dir_angle, elslib_cone_d1, elslib_cone_parameters, elslib_cylinder_d1,
    elslib_cylinder_parameters, elslib_plane_d0, elslib_plane_parameters, rotate_vec_around_ax1,
    surface3_ax3, GpAx22d, GpAx3, GpCirc, GpConicalSurface, GpCylindricalSurface, CONFUSION,
    P_CONFUSION,
};

/// OCCT gp_Dir(DVec3) helper — normalized direction.
fn gp_dir(v: DVec3) -> DVec3 {
    v.normalize()
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_ChAsymPlnCyl.cxx L65-369 — ChFiKPart_MakeChAsym
// (plane/cylinder, gp_Circ spine overload; conical chamfer surface).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_ch_asym_pln_cyl_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> bool {
    // compute the chamfer surface(cone)

    // compute the normals to the plane surface & to the plane face
    let pospl = surface3_ax3(&Surface3::Plane(*pln));
    let mut dpl = pospl.vxdir.cross(pospl.vydir).normalize();
    let mut norf = dpl;
    if ofpl == Orientation::Reversed {
        norf = -norf;
    }
    if or1 == Orientation::Reversed {
        dpl = -dpl;
    }

    // compute the origin Or of the cone
    let or = cyl.pos.location;
    let (u0, v0) = elslib_plane_parameters(&pospl, or);
    let pt2dpln = DVec2::new(u0, v0);
    // OCCT L101: ElSLib::PlaneD0(u, v, PosPl, Or) — the projection PtPl.
    let or_proj = elslib_plane_d0(u0, v0, &pospl);
    let ptpl = or_proj;

    let (ptsp_raw, dsp_raw) =
        super::chfi_kpart_gp::elclib_circle_d1(first, &spine.pos, spine.radius);
    let ptsp = ptsp_raw;
    let dsp = dsp_raw.normalize();
    let dx = gp_dir(ptsp - or_proj);
    let dy = dsp;
    let (u1, v1) = elslib_cylinder_parameters(&cyl.pos, cyl.radius, ptsp);
    let u = u1;
    let v = v1;
    let (_ptcyl, vu, vv) = elslib_cylinder_d1(u, v, &cyl.pos, cyl.radius);
    let _ = _ptcyl;
    let mut dcyl = vu.cross(vv).normalize(); // normal to the cylinder in PtSp
    if or2 == Orientation::Reversed {
        dcyl = -dcyl;
    }
    let dedans = dcyl.dot(dx) <= 0.0;

    let mut pointu = false;
    let conrad;
    let mut rad;
    let mut semiangl;

    // Calculation of distance
    let dis1;
    let dis2;
    let cosnpcyl;
    let sinnpcyl;

    if (plandab && dis_on_p) || (!plandab && !dis_on_p) {
        dis1 = dis;
        cosnpcyl = dpl.dot(dcyl);
        sinnpcyl = (1.0 - cosnpcyl * cosnpcyl).sqrt();
        dis2 = dis / (sinnpcyl / angle.tan() - cosnpcyl);
    } else {
        dis2 = dis;
        cosnpcyl = dpl.dot(dcyl);
        sinnpcyl = (1.0 - cosnpcyl * cosnpcyl).sqrt();
        dis1 = dis / (sinnpcyl / angle.tan() - cosnpcyl);
    }

    let mut or = or_proj + dpl * dis2;

    // variables used to compute the semiangle of the cone
    let vec1 = gp_dir(or - ptpl);
    let pt = or + pospl.vxdir * dis1;
    let vec2 = gp_dir(pt - ptpl);

    // compute the parameters of the conical surface
    if dedans {
        rad = cyl.radius - dis1;
        if rad.abs() <= CONFUSION {
            pointu = true;
        }
        if rad < 0.0 {
            // OCCT L161-164: "the chamfer can't pass".
            return false;
        }
    } else {
        rad = cyl.radius + dis1;
        dpl = -dpl;
    }
    conrad = cyl.radius;
    semiangl = dir_angle(vec1, vec2);
    let mut conax3 = GpAx3::new_pn_vx(or, dpl, dx);

    // OCCT L177: Geom_ConicalSurface(ConAx3, SemiAngl, ConRad).
    let mut gcon = GpConicalSurface::new(conax3, semiangl, conrad);

    // changes due to the fact the parameters of the chamfer must go increasing
    // from surface S1 to surface S2
    if (dedans && !plandab) || (!dedans && plandab) {
        gcon.v_reverse(); // be careful : the SemiAngle was changed
        conax3 = gcon.position();
        semiangl = gcon.semi_angle();
    }

    // changes due to the fact we have reversed the V direction of
    // parametrization
    if conax3.vydir.dot(dsp) <= 0.0 {
        conax3.y_reverse();
        gcon.set_position(conax3);
    }

    data.change_surf(chfi_kpart_index_surface_in_ds(gcon.to_surface3(), dstr));

    // compute the chamfer's orientation according to the orientation
    // of the faces

    // search the normal to the cone
    let (_pt, deru0, derv0) = elslib_cone_d1(0.0, 0.0, &conax3, conrad, semiangl);
    let _ = _pt;

    let norcon = deru0.cross(derv0).normalize();

    let toreverse = norcon.dot(norf) <= 0.0;
    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // we load of the faceInterference with the pcurves and
    //  the 3d curves

    // Case of the plane face
    // NB: in the case 'pointu', no pcurve on the plane surface
    // and no intersection plane-chamfer are needed
    let gcir2dpln: Option<rcad_kernel::geom::Curve2d>;
    let gcirpln: Option<Curve3>;
    let mut cirax2 = conax3.ax2();
    cirax2.set_location(ptpl);

    if !pointu {
        // intersection plane-chamfer
        let cirpln = GpCirc::new(cirax2, rad);
        gcirpln = Some(Curve3::Circle(cirpln.to_circle3()));

        // pcurve on the plane
        let (up, vp) = elslib_plane_parameters(&pospl, pt);
        let u = up;
        let v = vp;
        let p2dpln = DVec2::new(u, v);
        let d2d = DVec2::new(dsp.dot(pospl.vxdir), dsp.dot(pospl.vydir));
        let ax2dpln = GpAx22d::new(
            pt2dpln,
            (p2dpln - pt2dpln).normalize(),
            d2d.normalize(),
        );
        gcir2dpln = Some(ax2dpln.to_circle2d(rad));
    } else {
        gcir2dpln = None;
        gcirpln = None;
    }

    // pcurve on the chamfer
    let vch;
    if plandab {
        vch = -(dis1 * dis1 + dis2 * dis2).sqrt();
    } else {
        vch = (dis1 * dis1 + dis2 * dis2).sqrt();
    }
    let p2dch = DVec2::new(0.0, vch);
    let (_pt, _deru, _derv) = elslib_cone_d1(0.0, vch, &conax3, conrad, semiangl);
    let lin2dch = (p2dch, DVec2::new(1.0, 0.0));
    let glin2dch1 = curve2d_line(lin2dch.0, lin2dch.1);

    // orientation
    let trans;
    let norpl = pospl.vxdir.cross(pospl.vydir).normalize();
    let toreverse = norcon.dot(norpl) <= 0.0;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }

    if plandab {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcirpln, dstr),
            trans,
            gcir2dpln,
            Some(glin2dch1),
        );
    } else {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds_option(gcirpln, dstr),
            trans,
            gcir2dpln,
            Some(glin2dch1),
        );
    }

    // Case of the cylindrical face

    // intersection cylinder-chamfer
    cirax2.set_location(or);
    let circyl = GpCirc::new(cirax2, conrad);
    let gcircyl = Curve3::Circle(circyl.to_circle3());

    // pcurve on the chamfer
    let p2dch = DVec2::new(0.0, 0.0);
    let (_pt, deru, derv) = elslib_cone_d1(0.0, 0.0, &conax3, conrad, semiangl);
    let lin2dch = (p2dch, lin2dch.1);
    let glin2dch2 = curve2d_line(lin2dch.0, lin2dch.1);

    // pcurve on the cylinder
    // OCCT L300: norCon.SetXYZ(deru.Crossed(derv).XYZ());
    let norcon = deru.cross(derv).normalize();

    let pt = or + conrad * dx;
    let (u3, v3) = elslib_cylinder_parameters(&cyl.pos, cyl.radius, pt);
    let mut u = u3;
    let v = v3;
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

    let (_pt, deru, derv) = elslib_cylinder_d1(u, v, &cyl.pos, cyl.radius);
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

    // orientation
    let toreverse = norcon.dot(norcyl) <= 0.0;
    let trans;
    if (toreverse && plandab) || (!toreverse && !plandab) {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }

    if plandab {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds(gcircyl, dstr),
            trans,
            Some(glin2dcyl),
            Some(glin2dch2),
        );
    } else {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds(gcircyl, dstr),
            trans,
            Some(glin2dcyl),
            Some(glin2dch2),
        );
    }

    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_ChAsymPlnCyl.cxx L376-656 — ChFiKPart_MakeChAsym
// (plane/cylinder, gp_Lin spine overload; plane chamfer surface).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_ch_asym_pln_cyl_lin(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    _fu: f64,
    _lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &Line3,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> bool {
    // calculation of the fillet plane.
    // or1 and or2 permit to determine in which of four sides created by
    // intersection of 2 surfaces we are
    //        _|_          Ofpl is orientation of the plane face allowing
    //         |4          to determine the side of the material

    let orspine = spine.origin + spine.direction * first;

    let xdir = gp_dir(spine.direction);
    let axpln = surface3_ax3(&Surface3::Plane(*pln));
    let norpln = axpln.vxdir.cross(axpln.vydir).normalize();
    let mut norf = norpln;
    if or1 == Orientation::Reversed {
        norf = -norf;
    }

    let axcyl = cyl.pos;
    // OrCyl is the point on axis of cylinder in the plane normal to the
    // axis containing OrSpine (OCCT L413-417).
    let loc = axcyl.location;
    let locsp = orspine - loc;
    let axdir = axcyl.vzdir;
    let temp = axdir * locsp.dot(axdir);
    let orcyl = loc + temp;

    // construction of POnPln
    let mut tmp = orcyl - orspine;
    if (or2 == Orientation::Forward && axcyl.direct())
        || (or2 == Orientation::Reversed && !axcyl.direct())
    {
        tmp = -tmp;
    }

    let mut vectranslpln = xdir.cross(norpln).normalize();
    if vectranslpln.dot(tmp) <= 0.0 {
        vectranslpln = -vectranslpln;
    }

    let vectranslcyl = orspine - orcyl;

    // Calculation of distances dis1 and dis2, depending on Dis and Angle
    let dirsorc = gp_dir(vectranslcyl);
    let cos_a1 = dirsorc.dot(vectranslpln.normalize());
    let sin_a1 = (1.0 - cos_a1 * cos_a1).sqrt();
    let mut dis1 = 0.0;
    let mut dis2;
    let ray = cyl.radius;
    let is_dis_on_p = (plandab && dis_on_p) || (!plandab && !dis_on_p);

    if is_dis_on_p {
        dis1 = dis;
        let sin_al = angle.sin();
        let cos_al = angle.cos();
        let h = dis1 * sin_al;
        let cos_ahoc = cos_a1 * sin_al + sin_a1 * cos_al;
        let mut sin_ahoc = sin_a1 * sin_al - cos_a1 * cos_al;
        if cos_a1 > 0.0 {
            sin_ahoc = -sin_ahoc;
        }
        let temp1 = h / ray;
        let temp2 = sin_ahoc * sin_ahoc + temp1 * cos_ahoc;

        dis2 = temp2 + temp1 * (cos_ahoc - temp1);

        if dis2 < -1.0e-9 {
            // OCCT L462-465: "too great angle of chamfer".
            return false;
        } else if dis2 < 1.0e-9 {
            dis2 = ray * (2.0 * temp2).sqrt();
        } else {
            dis2 = ray * (2.0 * (temp2 - sin_ahoc * dis2.sqrt())).sqrt();
        }
    } else {
        dis2 = dis;
    }

    // construction of POnCyl
    let alpha = 2.0 * (dis2 * 0.5 / ray).asin();
    let mut vectemp = -vectranslcyl;

    if xdir.cross(gp_dir(vectranslcyl)).normalize().dot(norf) < 0.0 {
        vectemp = rotate_vec_around_ax1(vectemp, xdir, alpha);
    } else {
        vectemp = rotate_vec_around_ax1(vectemp, -xdir, alpha);
    }

    let poncyl0 = orcyl + vectemp;
    let poncyl = poncyl0;
    let (uoncyl, voncyl) = elslib_cylinder_parameters(&axcyl, cyl.radius, poncyl);
    let (poncyl_d1, duoncyl, dvoncyl) =
        elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius);
    let poncyl = poncyl_d1;

    // Construction of the point on the plane
    if !is_dis_on_p {
        let mut corde = orspine - poncyl;
        let mut tcyl = duoncyl.cross(dvoncyl);

        tcyl = xdir.cross(tcyl);
        tcyl = tcyl.normalize();
        corde = corde.normalize();
        let cos_cortan = tcyl.dot(corde);
        let mut tg_cortan = 1.0 / (cos_cortan * cos_cortan);
        tg_cortan = (tg_cortan - 1.0).sqrt();

        let mut tg_ang = angle.tan();
        tg_ang = (tg_ang + tg_cortan) / (1.0 - tg_ang * tg_cortan);

        let cos_a11 = dis2 / (2.0 * ray);
        let mut sin_a11 = (1.0 - cos_a11 * cos_a11).sqrt();
        if cos_a1 > 0.0 {
            sin_a11 = -sin_a11;
        }
        dis1 = (sin_a1 + cos_a1 * tg_ang) * cos_a11;
        dis1 -= (cos_a1 - sin_a1 * tg_ang) * sin_a11;
        dis1 = (dis2 * tg_ang) / dis1;
    }

    vectranslpln = vectranslpln * dis1;

    let ponpln0 = orspine + vectranslpln;
    let ponpln = ponpln0;

    // construction of the chamfer
    let (uonpln, vonpln) = elslib_plane_parameters(&axpln, ponpln);
    let ponpln = elslib_plane_d0(uonpln, vonpln, &axpln);

    // construction of YDir to go from face1 to face2.
    let mut ydir = poncyl - ponpln;
    if !plandab {
        ydir = -ydir;
    }
    let axch = GpAx3::new_pn_vx(ponpln, xdir.cross(ydir).normalize(), xdir);

    // OCCT L543-544: Geom_Plane(AxCh); ChangeSurf.
    let chamfer = super::chfi_kpart_gp::surface3_plane(&axch);
    data.change_surf(chfi_kpart_index_surface_in_ds(chamfer, dstr));

    // FaceInterferences are loaded with pcurves and curves 3d.
    //----------- edge plan-Chamfer
    let ppln2d = DVec2::new(uonpln, vonpln);
    let vpln2d = DVec2::new(xdir.dot(axpln.vxdir), xdir.dot(axpln.vydir));
    let lin2dpln = (ppln2d, vpln2d);

    let ponpln = elslib_plane_d0(uonpln, vonpln, &axpln);
    let c3d = Line3 {
        origin: ponpln,
        direction: xdir,
    };

    let (u, vonchamfer) = elslib_plane_parameters(&axch, ponpln);
    let lonchamfer = (DVec2::new(u, vonchamfer), DVec2::new(1.0, 0.0));

    let l3d = Curve3::Line(c3d);
    let lfac = curve2d_line(lin2dpln.0, lin2dpln.1);
    let lfil = curve2d_line(lonchamfer.0, lonchamfer.1);

    let norfil = axch.vzdir;
    let mut toreverse = norfil.dot(norpln) <= 0.0;
    let dirplncyl = gp_dir(poncyl - ponpln);
    let dirspln = gp_dir(ponpln - orspine);
    let poschamfpln = dirplncyl.dot(dirspln) > 0.0;

    if !is_dis_on_p && poschamfpln {
        toreverse = !toreverse;
    }
    // It is checked if the orientation of the Chamfer is the same as of the plane
    if toreverse {
        *data.change_orientation() = super::chfi_kpart_gp::topabs_reverse(ofpl);
    } else {
        *data.change_orientation() = ofpl;
    }

    let mut trans = Orientation::Forward;
    if (!plandab && toreverse) || (plandab && !toreverse) {
        trans = Orientation::Reversed;
    }

    // trans allows to determine the "material" side on S1(2) limited by L3d
    if plandab {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    } else {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    }

    //------------edge cylindre-Chamfer
    let pcyl2d = DVec2::new(uoncyl, voncyl);
    let mut vcyl2d = DVec2::new(0.0, 1.0);
    if xdir.dot(axcyl.vzdir) < 0.0 {
        vcyl2d = -vcyl2d;
    }
    let lin2dcyl = (pcyl2d, vcyl2d);

    let poncyl = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius).0;
    let c3d = Line3 {
        origin: poncyl,
        direction: xdir,
    };

    let (u, vonchamfer) = elslib_plane_parameters(&axch, poncyl);
    let lonchamfer = (DVec2::new(u, vonchamfer), DVec2::new(1.0, 0.0));

    let l3d = Curve3::Line(c3d);
    let lfac = curve2d_line(lin2dcyl.0, lin2dcyl.1);
    let lfil = curve2d_line(lonchamfer.0, lonchamfer.1);

    let (_p, deru, derv) = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius);
    let norcyl = deru.cross(derv).normalize();
    let dirscyl = gp_dir(poncyl - orspine);
    let poschamfcyl = dirplncyl.dot(dirscyl) < 0.0;
    let mut toreverse = norfil.dot(norcyl) <= 0.0;

    if is_dis_on_p && poschamfcyl {
        toreverse = !toreverse;
    }
    let mut trans = Orientation::Reversed;
    if (!plandab && toreverse) || (plandab && !toreverse) {
        trans = Orientation::Forward;
    }

    if plandab {
        data.change_interference_on_s2().set_interference(
            chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    } else {
        data.change_interference_on_s1().set_interference(
            chfi_kpart_index_curve_in_ds(l3d, dstr),
            trans,
            Some(lfac),
            Some(lfil),
        );
    }
    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_ChAsymPlnCon.cxx L57-708 — ChFiKPart_MakeChAsym
// (plane/cone, gp_Circ spine overload; conical or cylindrical chamfer).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_ch_asym_pln_con_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pln: &Plane,
    con: &GpConicalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
    dis_on_p: bool,
) -> bool {
    // Compute the chamfer surface(cone)
    let pospl = surface3_ax3(&Surface3::Plane(*pln));
    let mut dpl = pospl.vxdir.cross(pospl.vydir).normalize();
    let mut norf = dpl;
    if ofpl == Orientation::Reversed {
        norf = -norf;
    }
    if or1 == Orientation::Reversed {
        dpl = -dpl;
    }

    // compute the origin of the conical chamfer PtPl
    let or0 = con.pos.location;
    let (u0, v0) = elslib_plane_parameters(&pospl, or0);
    let pt2dpln = DVec2::new(u0, v0);
    // OCCT L91: ElSLib::PlaneD0(u, v, PosPl, Or) — the projection PtPl.
    let ptpl = elslib_plane_d0(u0, v0, &pospl);

    let (ptsp_raw, dsp_raw) =
        super::chfi_kpart_gp::elclib_circle_d1(first, &spine.pos, spine.radius);
    let ptsp = ptsp_raw;
    let dsp = dsp_raw.normalize();
    let dx = gp_dir(ptsp - ptpl);

    // compute the normal to the cone in PtSp
    let (u1, v1) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, ptsp);
    let u = u1;
    let v = v1;
    let (_ptcon, deru, derv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
    let _ = _ptcon;
    let mut dcon = deru.cross(derv).normalize();
    if or2 == Orientation::Reversed {
        dcon = -dcon;
    }

    let dedans = dx.dot(dcon) <= 0.0;
    let ouvert = dpl.dot(dcon) >= 0.0;

    // variables used to compute the semiangle of the chamfer
    let angcon = con.semi_angle;
    // OCCT keeps the Dis parameter and a separate local `dis`; Rust cannot
    // distinguish them by case, so the parameter is captured here.
    let dis_param = dis;
    let mut move_;
    let mut chamfrad;
    let mut semiangl;
    let mut pointu = false;
    let mut dis;
    let mut iscylinder = false;
    let mut is_con_par = false;

    if (plandab && dis_on_p) || (!plandab && !dis_on_p) {
        let tgang = angle.tan();
        let tgcon = angcon.tan().abs();
        let dis11;
        if ouvert {
            move_ = dis_param * tgang / (1.0 - tgcon * tgcon);
            dis11 = move_ * tgcon;
            dis = dis_param + dis11;
        } else {
            move_ = dis_param * tgang / (1.0 + tgcon * tgcon);
            dis11 = move_ * tgcon;
            dis = dis_param - dis11;
        }

        // compute the parameters of the conical chamfer
        if dedans {
            chamfrad = spine.radius - dis;
            if chamfrad.abs() < CONFUSION {
                pointu = true;
            }
            if chamfrad < 0.0 {
                // OCCT L149-152: "the chamfer can't pass".
                return false;
            }
            semiangl = std::f64::consts::FRAC_PI_2 - angle;
        } else {
            chamfrad = spine.radius + dis;
            semiangl = angle - std::f64::consts::FRAC_PI_2;
        }

        if ouvert {
            if angcon.abs() - semiangl.abs() > -CONFUSION {
                // OCCT L166-169: "wrong choice of angle for the chamfer".
                return false;
            }
        }
    } else {
        let dis1;
        move_ = dis_param * angcon.cos();
        if ouvert {
            semiangl = angcon.abs() + angle;

            if (std::f64::consts::FRAC_PI_2 - semiangl) < CONFUSION {
                // OCCT L183-186: "wrong choice of angle for the chamfer".
                return false;
            }
            dis1 = move_ * semiangl.tan() - dis_param * angcon.sin().abs();

            if !dedans {
                semiangl = -semiangl;
            }
        } else {
            semiangl = angcon.abs() - angle;

            if semiangl.abs() < CONFUSION {
                iscylinder = true;
                dis1 = dis_param * angcon.sin().abs();
            } else {
                dis1 = dis_param * angcon.sin().abs() - move_ * semiangl.tan();
            }

            if semiangl > CONFUSION {
                is_con_par = true;
            }

            if dedans {
                semiangl = -semiangl;
            }
        }

        // compute the parameters of the conical chamfer
        if dedans {
            chamfrad = spine.radius - dis1;

            if chamfrad.abs() < CONFUSION {
                pointu = true;
            }
            if chamfrad < 0.0 {
                // OCCT L231-234: "the chamfer can't pass".
                return false;
            }
        } else {
            chamfrad = spine.radius + dis1;
        }

        if ouvert {
            dis = dis1 + dis_param * angcon.sin().abs();
        } else {
            dis = dis1 - dis_param * angcon.sin().abs();
        }
    }

    let or = ptpl + dpl * move_;
    let pt0 = or + pospl.vxdir * dis;
    let pt = pt0;

    let mut chamfax3 = GpAx3::new_pn_vx(ptpl, dpl, dx);

    if iscylinder {
        // OCCT L262: Geom_CylindricalSurface(ChamfAx3, ChamfRad).
        let mut gcyl = GpCylindricalSurface::new(chamfax3, chamfrad);

        // changes due to the fact the parameters of the chamfer must go
        // increasing from surface S1 to surface S2
        if !plandab {
            gcyl.v_reverse(); // be careful : the SemiAngle was changed
            chamfax3 = gcyl.position();
        }

        // changes due to the fact we have reversed the V direction of
        // parametrization
        if chamfax3.vydir.dot(dsp) <= 0.0 {
            chamfax3.y_reverse();
            gcyl.set_position(chamfax3);
        }

        data.change_surf(chfi_kpart_index_surface_in_ds(gcyl.to_surface3(), dstr));

        let mut torevcha = !chamfax3.direct();
        let mut cylaxe = chamfax3.vzdir;
        torevcha = (torevcha && !plandab) || (!torevcha && plandab);

        if torevcha {
            cylaxe = -cylaxe;
        }
        let toreverse = norf.dot(cylaxe) < 0.0;

        if (toreverse && dedans) || (!toreverse && !dedans) {
            *data.change_orientation() = Orientation::Reversed;
        } else {
            *data.change_orientation() = Orientation::Forward;
        }

        // we load the faceInterference with the pcurves and
        //  the 3d curves

        // Case of the plane face
        // NB: in the case 'pointu', no pcurve on the plane surface
        // and no intersection plane-chamfer are needed

        // intersection plane-chamfer
        let mut cirax2 = chamfax3.ax2();
        cirax2.set_location(ptpl);

        let pt = ptpl + dx * chamfrad;
        let cirpln = GpCirc::new(cirax2, chamfrad);
        let gcirpln = Curve3::Circle(cirpln.to_circle3());

        // pcurve on the plane
        let (up, vp) = elslib_plane_parameters(&pospl, pt);
        let u = up;
        let v = vp;
        let p2dpln = DVec2::new(u, v);
        let d2d = DVec2::new(dsp.dot(pospl.vxdir), dsp.dot(pospl.vydir));
        let ax2dpln = GpAx22d::new(
            pt2dpln,
            (p2dpln - pt2dpln).normalize(),
            d2d.normalize(),
        );
        let gcir2dpln = ax2dpln.to_circle2d(chamfrad);

        // pcurve on chamfer
        let p2dch = DVec2::new(0.0, 0.0);
        let lin2dch = (p2dch, DVec2::new(1.0, 0.0));
        let glin2dch1 = curve2d_line(lin2dch.0, lin2dch.1);

        // orientation
        let trans;
        let norpl = pospl.vxdir.cross(pospl.vydir).normalize();
        let mut toreverse = norpl.dot(cylaxe) < 0.0;

        toreverse = (toreverse && plandab) || (!toreverse && !plandab);

        if (toreverse && dedans) || (!toreverse && !dedans) {
            trans = Orientation::Forward;
        } else {
            trans = Orientation::Reversed;
        }

        if plandab {
            data.change_interference_on_s1().set_interference(
                chfi_kpart_index_curve_in_ds(gcirpln, dstr),
                trans,
                Some(gcir2dpln),
                Some(glin2dch1),
            );
        } else {
            data.change_interference_on_s2().set_interference(
                chfi_kpart_index_curve_in_ds(gcirpln, dstr),
                trans,
                Some(gcir2dpln),
                Some(glin2dch1),
            );
        }

        // Case of the conical face

        // intersection cone-chamfer
        let rad;
        if dedans {
            rad = chamfrad + dis;
        } else {
            rad = chamfrad - dis;
        }

        cirax2.set_location(or);
        let circon = GpCirc::new(cirax2, rad);
        let gcircon = Curve3::Circle(circon.to_circle3());

        // pcurve on chamfer
        let vch;
        if plandab {
            vch = (dis * dis + move_ * move_).sqrt();
        } else {
            vch = -(dis * dis + move_ * move_).sqrt();
        }
        let p2dch = DVec2::new(0.0, vch);
        let (_pt, _deru, _derv) =
            elslib_cylinder_d1(0.0, vch, &chamfax3, chamfrad);
        let lin2dch = (p2dch, lin2dch.1);
        let glin2dch2 = curve2d_line(lin2dch.0, lin2dch.1);

        // pcurve on cone
        let pt = or + rad * dx;
        let (u3, v3) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, pt);
        let mut u = u3;
        let v = v3;
        let tol = P_CONFUSION;
        if u >= 2.0 * std::f64::consts::PI - tol && u <= 2.0 * std::f64::consts::PI {
            u = 0.0;
        }
        if u >= fu - tol && u < fu {
            u = fu;
        }
        if u <= lu + tol && u > lu {
            u = lu;
        }
        if u < fu || u > lu {
            u = chfi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, P_CONFUSION);
        }
        let (_pt, deru, derv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
        let p2dcon = DVec2::new(u, v);
        let d2dcon;
        if deru.dot(dsp) <= 0.0 {
            d2dcon = -DVec2::new(1.0, 0.0);
        } else {
            d2dcon = DVec2::new(1.0, 0.0);
        }
        let lin2dcon = (p2dcon, d2dcon);
        let glin2dcon = curve2d_line(lin2dcon.0, lin2dcon.1);

        // orientation
        let norcon = deru.cross(derv).normalize();

        let mut dircon = con.pos.vzdir;
        if angcon > CONFUSION {
            dircon = -dircon;
        }
        let torevcon = norcon.dot(dircon) < 0.0;

        let trans;
        if (torevcon && dedans) || (!torevcon && !dedans) {
            trans = Orientation::Reversed;
        } else {
            trans = Orientation::Forward;
        }

        if plandab {
            data.change_interference_on_s2().set_interference(
                chfi_kpart_index_curve_in_ds(gcircon, dstr),
                trans,
                Some(glin2dcon),
                Some(glin2dch2),
            );
        } else {
            data.change_interference_on_s1().set_interference(
                chfi_kpart_index_curve_in_ds(gcircon, dstr),
                trans,
                Some(glin2dcon),
                Some(glin2dch2),
            );
        }
    } else {
        // OCCT L467: Geom_ConicalSurface(ChamfAx3, SemiAngl, ChamfRad).
        let mut gcon = GpConicalSurface::new(chamfax3, semiangl, chamfrad);

        // changes due to the fact the parameters of the chamfer must go
        // increasing from surface S1 to surface S2
        if !plandab {
            gcon.v_reverse(); // be careful : the SemiAngle was changed
            chamfax3 = gcon.position();
            semiangl = gcon.semi_angle();
        }

        // changes due to the fact we have reversed the V direction of
        // parametrization
        if chamfax3.vydir.dot(dsp) <= 0.0 {
            chamfax3.y_reverse();
            gcon.set_position(chamfax3);
        }

        data.change_surf(chfi_kpart_index_surface_in_ds(gcon.to_surface3(), dstr));

        // compute the chamfer's orientation according to the orientation
        //  of the faces

        // search the normal to the conical chamfer
        let mut u;
        let v;
        u = 0.0;
        if plandab {
            v = (dis * dis + move_ * move_).sqrt();
        } else {
            v = -(dis * dis + move_ * move_).sqrt();
        }

        let (_p, deru, derv) = elslib_cone_d1(u, v, &chamfax3, chamfrad, semiangl);
        let mut norchamf = deru.cross(derv).normalize();

        let mut toreverse = norf.dot(norchamf) < 0.0;

        if is_con_par {
            toreverse = !toreverse;
        }

        if toreverse {
            *data.change_orientation() = Orientation::Reversed;
        } else {
            *data.change_orientation() = Orientation::Forward;
        }

        // we load the faceInterference with the pcurves and
        //  the 3d curves

        // Case of the plane face
        // NB: in the case 'pointu', no pcurve on the plane surface
        // and no intersection plane-chamfer are needed

        // intersection plane-chamfer
        let gcir2dpln: Option<rcad_kernel::geom::Curve2d>;
        let gcirpln: Option<Curve3>;
        let mut cirax2 = chamfax3.ax2();
        cirax2.set_location(ptpl);

        if !pointu {
            let pt = ptpl + dx * chamfrad;
            let cirpln = GpCirc::new(cirax2, chamfrad);
            gcirpln = Some(Curve3::Circle(cirpln.to_circle3()));

            // pcurve on the plane
            let (up, vp) = elslib_plane_parameters(&pospl, pt);
            let u = up;
            let v = vp;
            let p2dpln = DVec2::new(u, v);
            let d2d = DVec2::new(dsp.dot(pospl.vxdir), dsp.dot(pospl.vydir));
            let ax2dpln = GpAx22d::new(
                pt2dpln,
                (p2dpln - pt2dpln).normalize(),
                d2d.normalize(),
            );
            gcir2dpln = Some(ax2dpln.to_circle2d(chamfrad));
        } else {
            gcir2dpln = None;
            gcirpln = None;
        }

        // pcurve on chamfer
        let p2dch = DVec2::new(0.0, 0.0);
        let (_pt, deru, derv) = elslib_cone_d1(0.0, 0.0, &chamfax3, chamfrad, semiangl);
        let lin2dch = (p2dch, DVec2::new(1.0, 0.0));
        let glin2dch1 = curve2d_line(lin2dch.0, lin2dch.1);

        // orientation
        let trans;
        let norpl = pospl.vxdir.cross(pospl.vydir).normalize();
        if !pointu {
            norchamf = deru.cross(derv).normalize();
        }
        let mut toreverse = norchamf.dot(norpl) <= 0.0;

        if is_con_par {
            toreverse = !toreverse;
        }

        if (toreverse && plandab) || (!toreverse && !plandab) {
            trans = Orientation::Forward;
        } else {
            trans = Orientation::Reversed;
        }

        if plandab {
            data.change_interference_on_s1().set_interference(
                chfi_kpart_index_curve_in_ds_option(gcirpln, dstr),
                trans,
                gcir2dpln,
                Some(glin2dch1),
            );
        } else {
            data.change_interference_on_s2().set_interference(
                chfi_kpart_index_curve_in_ds_option(gcirpln, dstr),
                trans,
                gcir2dpln,
                Some(glin2dch1),
            );
        }

        // Case of the conical face

        // intersection cone-chamfer
        let rad;
        if dedans {
            rad = chamfrad + dis;
        } else {
            rad = chamfrad - dis;
        }

        cirax2.set_location(or);
        let circon = GpCirc::new(cirax2, rad);
        let gcircon = Curve3::Circle(circon.to_circle3());

        // pcurve on chamfer
        let vch;
        if plandab {
            vch = (dis * dis + move_ * move_).sqrt();
        } else {
            vch = -(dis * dis + move_ * move_).sqrt();
        }
        let p2dch = DVec2::new(0.0, vch);
        let (_pt, deru, derv) = elslib_cone_d1(0.0, vch, &chamfax3, chamfrad, semiangl);
        let lin2dch = (p2dch, lin2dch.1);
        let glin2dch2 = curve2d_line(lin2dch.0, lin2dch.1);

        // pcurve on cone
        // OCCT L629: norchamf.SetXYZ(deru.Crossed(derv).XYZ());
        let norchamf = deru.cross(derv).normalize();

        let pt = or + rad * dx;
        let (u3, v3) = elslib_cone_parameters(&con.pos, con.ref_radius, con.semi_angle, pt);
        let mut u = u3;
        let v = v3;
        let tol = P_CONFUSION;
        if u >= 2.0 * std::f64::consts::PI - tol && u <= 2.0 * std::f64::consts::PI {
            u = 0.0;
        }
        if u >= fu - tol && u < fu {
            u = fu;
        }
        if u <= lu + tol && u > lu {
            u = lu;
        }
        if u < fu || u > lu {
            u = chfi_kpart_in_period(u, fu, fu + 2.0 * std::f64::consts::PI, P_CONFUSION);
        }
        let (_pt, deru, derv) = elslib_cone_d1(u, v, &con.pos, con.ref_radius, con.semi_angle);
        let p2dcon = DVec2::new(u, v);
        let d2dcon;
        if deru.dot(dsp) <= 0.0 {
            d2dcon = -DVec2::new(1.0, 0.0);
        } else {
            d2dcon = DVec2::new(1.0, 0.0);
        }
        let lin2dcon = (p2dcon, d2dcon);
        let glin2dcon = curve2d_line(lin2dcon.0, lin2dcon.1);

        // orientation
        let norcon = deru.cross(derv).normalize();

        let mut dircon = con.pos.vzdir;
        let mut dirchamf = gcon.axis_direction();
        if angcon > CONFUSION {
            dircon = -dircon;
        }
        if semiangl > CONFUSION {
            dirchamf = -dirchamf;
        }

        let torevcon = norcon.dot(dircon) > 0.0;
        let torevcha = norchamf.dot(dirchamf) > 0.0;

        let toreverse = (torevcon && !torevcha) || (!torevcon && torevcha);

        let trans;
        if (toreverse && plandab) || (!toreverse && !plandab) {
            trans = Orientation::Reversed;
        } else {
            trans = Orientation::Forward;
        }

        if plandab {
            data.change_interference_on_s2().set_interference(
                chfi_kpart_index_curve_in_ds(gcircon, dstr),
                trans,
                Some(glin2dcon),
                Some(glin2dch2),
            );
        } else {
            data.change_interference_on_s1().set_interference(
                chfi_kpart_index_curve_in_ds(gcircon, dstr),
                trans,
                Some(glin2dcon),
                Some(glin2dch2),
            );
        }
    }
    true
}
