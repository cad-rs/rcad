//! OCCT ChFiKPart — chamfer particular cases plane/cylinder and plane/cone:
//!   - ChFiKPart_ComputeData_ChPlnCyl.cxx L64-365 (gp_Circ spine overload)
//!   - ChFiKPart_ComputeData_ChPlnCyl.cxx L372-595 (gp_Lin spine overload)
//!   - ChFiKPart_ComputeData_ChPlnCon.cxx L49-164 (gp_Circ spine overload;
//!     delegates to ChFiKPart_MakeChAsym plane/cone)
//!
//! 1:1 translation; the OCCT gp/ElSLib carriers come from
//! [`super::chfi_kpart_gp`].

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topo::topods::Orientation;

use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::{ChFiDS_ChamfMode, ChFiDSSurfData};
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
// OCCT ChFiKPart_ComputeData_ChPlnCyl.cxx L64-365 — ChFiKPart_MakeChamfer
// (plane/cylinder, gp_Circ spine overload; conical chamfer surface).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_chamfer_pln_cyl_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    the_mode: ChFiDS_ChamfMode,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    the_dis1: f64,
    the_dis2: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    // compute the chamfer surface(cone)

    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    if the_mode == ChFiDS_ChamfMode::ConstThroatChamfer {
        dis1 = the_dis1 * 2.0f64.sqrt();
        dis2 = dis1;
    } else if the_mode == ChFiDS_ChamfMode::ConstThroatWithPenetrationChamfer {
        let a_dis2 = the_dis1.min(the_dis2);
        let a_dis1 = the_dis1.max(the_dis2);
        dis2 = (a_dis1 * a_dis1 - a_dis2 * a_dis2).sqrt();
        dis1 = a_dis1 * a_dis1 / a_dis2 - a_dis2;
    }

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
    // OCCT L114: ElSLib::PlaneD0(u, v, PosPl, Or) — Or becomes the
    // projection of the cylinder origin on the plane.
    let or_proj = elslib_plane_d0(u0, v0, &pospl);
    let ptpl = or_proj;

    // start 3d point on the Spine + tangent vector
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
    let mut or = or_proj + dpl * dis2;

    // variables used to compute the semiangle of the cone
    let vec1 = gp_dir(or - ptpl);
    let pt0 = or + pospl.vxdir * dis1;
    let pt = pt0;
    let vec2 = gp_dir(pt - ptpl);

    // compute the parameters of the conical surface
    if dedans {
        rad = cyl.radius - dis1;
        if rad.abs() <= CONFUSION {
            pointu = true;
        }
        if rad < 0.0 {
            // OCCT L155-158: "the chamfer can't pass".
            return false;
        }
    } else {
        rad = cyl.radius + dis1;
        dpl = -dpl;
    }
    conrad = cyl.radius;
    semiangl = dir_angle(vec1, vec2);
    let mut conax3 = GpAx3::new_pn_vx(or, dpl, dx);

    // OCCT L171: Geom_ConicalSurface(ConAx3, SemiAngl, ConRad).
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
    let mut deru;
    let mut derv;
    let (_pt, du, dv) = elslib_cone_d1(0.0, 0.0, &conax3, conrad, semiangl);
    deru = du;
    derv = dv;

    let norcon = deru.cross(derv).normalize();

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
    let (_pt, du, dv) = elslib_cone_d1(0.0, vch, &conax3, conrad, semiangl);
    deru = du;
    derv = dv;
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
    let (_pt, du, dv) = elslib_cone_d1(0.0, 0.0, &conax3, conrad, semiangl);
    deru = du;
    derv = dv;
    let lin2dch = (p2dch, lin2dch.1);
    let glin2dch2 = curve2d_line(lin2dch.0, lin2dch.1);

    // pcurve on the cylinder
    // OCCT L296: norCon.SetXYZ(deru.Crossed(derv).XYZ());
    let norcon = deru.cross(derv).normalize();

    let pt = or + dx * conrad;
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
// OCCT ChFiKPart_ComputeData_ChPlnCyl.cxx L372-595 — ChFiKPart_MakeChamfer
// (plane/cylinder, gp_Lin spine overload; plane chamfer surface).
// =========================================================================
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
pub fn chfi_kpart_make_chamfer_pln_cyl_lin(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    the_mode: ChFiDS_ChamfMode,
    pln: &Plane,
    cyl: &GpCylindricalSurface,
    _fu: f64,
    _lu: f64,
    or1: Orientation,
    or2: Orientation,
    dis1: f64,
    dis2: f64,
    spine: &Line3,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    // calculation of the fillet plane.
    // or1 and or2 permit to determine in which of four sides created by
    // intersection of 2 surfaces we are
    //        _|_          Ofpl is orientation of the plane face allowing
    //         |4          to determine the side of the material

    if the_mode != ChFiDS_ChamfMode::ClassicChamfer {
        return false;
    }

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
    // axis containing OrSpine
    // Project <OrSpine> onto <AxCyl> (OCCT L415-418).
    let axloc = axcyl.location;
    let axdir = axcyl.vzdir;
    let parameter = (orspine - axloc).dot(axdir);
    let orcyl = axloc + axdir * parameter;

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
    vectranslpln = vectranslpln * dis1;

    let ponpln0 = orspine + vectranslpln;
    let ponpln = ponpln0;

    // construction of POnCyl
    let alpha = 2.0 * (dis2 * 0.5 / cyl.radius).asin();

    let mut veccyltransl = orspine - orcyl;

    if xdir.cross(gp_dir(veccyltransl)).normalize().dot(norf) > 0.0 {
        veccyltransl = rotate_vec_around_ax1(veccyltransl, xdir, alpha);
    } else {
        veccyltransl = rotate_vec_around_ax1(veccyltransl, -xdir, alpha);
    }

    let poncyl0 = orcyl + veccyltransl;
    let poncyl = poncyl0;

    // construction of chamfer
    let (uoncyl, voncyl) = elslib_cylinder_parameters(&axcyl, cyl.radius, poncyl);
    let poncyl = elslib_cylinder_d1(uoncyl, voncyl, &axcyl, cyl.radius).0;
    let (uonpln, vonpln) = elslib_plane_parameters(&axpln, ponpln);
    let ponpln = elslib_plane_d0(uonpln, vonpln, &axpln);

    // construction of YDir to go to face1 from face2.
    let mut ydir = poncyl - ponpln;
    if !plandab {
        ydir = -ydir;
    }
    let axch = GpAx3::new_pn_vx(ponpln, xdir.cross(ydir).normalize(), xdir);

    // OCCT L478-479: Geom_Plane(AxCh); ChangeSurf.
    let chamfer = super::chfi_kpart_gp::surface3_plane(&axch);
    data.change_surf(chfi_kpart_index_surface_in_ds(chamfer, dstr));

    // FaceInterferences are loaded with pcurves and curves 3d.
    //----------- edge plane-Chamfer
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

    if poschamfpln {
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

    // trans permits to determine the "material" side on S1(2) limited by L3d
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

    //------------edge cylinder-Chamfer
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

    let mut toreverse = norfil.dot(norcyl) <= 0.0;

    let dirscyl = gp_dir(poncyl - orspine);
    let poschamfcyl = dirplncyl.dot(dirscyl) < 0.0;

    if poschamfcyl {
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
// OCCT ChFiKPart_ComputeData_ChPlnCon.cxx L49-164 — ChFiKPart_MakeChamfer
// (plane/cone, gp_Circ spine overload; delegates to ChFiKPart_MakeChAsym).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_chamfer_pln_con_circ(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    the_mode: ChFiDS_ChamfMode,
    pln: &Plane,
    con: &GpConicalSurface,
    fu: f64,
    lu: f64,
    or1: Orientation,
    or2: Orientation,
    the_dis1: f64,
    the_dis2: f64,
    spine: &GpCirc,
    first: f64,
    ofpl: Orientation,
    plandab: bool,
) -> bool {
    let angcon = con.semi_angle;

    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    let alpha = std::f64::consts::FRAC_PI_2 - angcon;
    let cos_half_alpha = (alpha / 2.0).cos();
    if the_mode == ChFiDS_ChamfMode::ConstThroatChamfer {
        dis1 = the_dis1 / cos_half_alpha;
        dis2 = dis1;
    } else if the_mode == ChFiDS_ChamfMode::ConstThroatWithPenetrationChamfer {
        let a_dis1 = the_dis1.min(the_dis2);
        let a_dis2 = the_dis1.max(the_dis2);
        let dis1dis1 = a_dis1 * a_dis1;
        let dis2dis2 = a_dis2 * a_dis2;
        let sin_alpha = alpha.sin();
        let cos_alpha = alpha.cos();
        let cotan_alpha = cos_alpha / sin_alpha;
        dis1 = (dis2dis2 - dis1dis1).sqrt() - a_dis1 * cotan_alpha;
        let cos_beta =
            (1.0 - dis1dis1 / dis2dis2).sqrt() * cos_alpha + a_dis1 / a_dis2 * sin_alpha;
        let full_dist1 = a_dis2 / cos_beta;
        dis2 = full_dist1 - a_dis1 / sin_alpha;
    }

    let sincon = angcon.sin().abs();
    let angle;

    let pospl = surface3_ax3(&Surface3::Plane(*pln));
    let mut dpl = pospl.vxdir.cross(pospl.vydir).normalize();
    if or1 == Orientation::Reversed {
        dpl = -dpl;
    }

    // compute the origin of the conical chamfer PtPl
    let or = con.pos.location;
    let (u0, v0) = elslib_plane_parameters(&pospl, or);
    let _ = (u0, v0);
    let _or_proj = elslib_plane_d0(u0, v0, &pospl);

    let (ptsp_raw, dsp_raw) =
        super::chfi_kpart_gp::elclib_circle_d1(first, &spine.pos, spine.radius);
    let ptsp = ptsp_raw;
    let dsp = dsp_raw.normalize();
    let _ = dsp;

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

    let ouvert = dpl.dot(dcon) >= 0.0;

    if !ouvert {
        if (dis1 - dis2 * sincon).abs() > CONFUSION {
            let abscos = (dis2 - dis1 * sincon).abs();
            angle = (dis1 * angcon.cos() / abscos).atan();
        } else {
            angle = angcon;
        }
    } else {
        angle = (dis1 * angcon.cos() / (dis2 + dis1 * sincon)).atan();
    }

    let dis_on_p = false;

    let is_resol = super::chfi_kpart_chasym::chfi_kpart_make_ch_asym_pln_con_circ(
        dstr,
        data,
        pln,
        con,
        fu,
        lu,
        or1,
        or2,
        dis2,
        angle,
        spine,
        first,
        ofpl,
        plandab,
        dis_on_p,
    );

    is_resol
}
