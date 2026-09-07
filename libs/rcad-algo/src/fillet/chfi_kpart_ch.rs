//! OCCT ChFiKPart — chamfer particular cases on plane/plane:
//!   - ChFiKPart_ComputeData_ChPlnPln.cxx L53-272 (ChFiKPart_MakeChamfer)
//!   - ChFiKPart_ComputeData_ChAsymPlnPln.cxx L53-266 (ChFiKPart_MakeChAsym)
//!
//! 1:1 translation; the OCCT gp/ElSLib carriers come from
//! [`super::chfi_kpart_gp`].

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topo::topods::Orientation;

use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::{ChFiDS_ChamfMode, ChFiDSSurfData};
use super::chfi_kpart::chfi_kpart_index_curve_in_ds;
use super::chfi_kpart_gp::{
    curve2d_line, dir_angle, elslib_plane_parameters, surface3_ax3, surface3_plane, GpAx3,
};

/// OCCT gp_Dir(DVec3) helper — normalized direction.
fn gp_dir(v: DVec3) -> DVec3 {
    v.normalize()
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_ChPlnPln.cxx L53-272 — ChFiKPart_MakeChamfer
// (plane/plane, all symmetric chamfer modes).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_chamfer_pln_pln(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    the_mode: ChFiDS_ChamfMode,
    pl1: &Plane,
    pl2: &Plane,
    or1: Orientation,
    or2: Orientation,
    the_dis1: f64,
    the_dis2: f64,
    spine: &Line3,
    first: f64,
    of1: Orientation,
) -> bool {
    // Creation of the plane which carry the chamfer

    // compute the normals to the planes Pl1 and Pl2
    let pos1 = surface3_ax3(&Surface3::Plane(*pl1));
    let mut d1 = pos1.vxdir.cross(pos1.vydir).normalize();
    if or1 == Orientation::Reversed {
        d1 = -d1;
    }
    let pos2 = surface3_ax3(&Surface3::Plane(*pl2));
    let mut d2 = pos2.vxdir.cross(pos2.vydir).normalize();
    if or2 == Orientation::Reversed {
        d2 = -d2;
    }

    // compute the intersection line of Pl1 and Pl2
    let lint = rcad_kernel::base::int_ana::intersect_plane_plane_intana(pl1, pl2);

    let p;
    let fint;
    match &lint {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(lint_line) => {
            let orspine = spine.origin + spine.direction * first;
            fint = (orspine - lint_line.origin).dot(lint_line.direction);
            p = lint_line.origin + lint_line.direction * fint;
        }
        _ => {
            return false;
        }
    }

    let linax1 = gp_dir(spine.direction);
    let mut vectransl1 = gp_dir(linax1.cross(d1));
    if vectransl1.dot(d2) <= 0.0 {
        vectransl1 = -vectransl1;
    }

    let mut vectransl2 = gp_dir(linax1.cross(d2));
    if vectransl2.dot(d1) <= 0.0 {
        vectransl2 = -vectransl2;
    }

    let mut dis1 = the_dis1;
    let mut dis2 = the_dis2;
    let alpha = dir_angle(vectransl1, vectransl2);
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

    // Compute a point on the plane Pl1 and on the chamfer
    let p1 = p + vectransl1 * dis1;

    // Point on the plane Pl2 and on the chamfer
    let p2 = p + vectransl2 * dis2;

    // the middle point of P1 P2 is the origin of the chamfer
    let po = (p1 + p2) / 2.0;

    // compute a second point on the plane Pl2
    let pp = {
        let lint_line = match &lint {
            rcad_kernel::base::int_ana::PlnPlnResult::Line(l) => l,
            _ => unreachable!(),
        };
        lint_line.origin + lint_line.direction * (fint + 10.0)
    };
    let p22 = pp + vectransl2 * dis2;

    // Compute the normal vector <AxisPlan> to the chamfer's plane
    let v1 = gp_dir(p2 - p1);
    let v2 = gp_dir(p22 - p1);
    let axisplan = gp_dir(v1.cross(v2));

    let xdir = linax1; // u axis
    let mut planax3 = GpAx3::new_pn_vx(po, axisplan, xdir);
    if planax3.vydir.dot(d2) >= 0.0 {
        planax3.y_reverse();
    }

    // OCCT L163-164: Geom_Plane(PlanAx3); ChangeSurf.
    let gpl = surface3_plane(&planax3);
    data.change_surf(super::chfi_kpart::chfi_kpart_index_surface_in_ds(gpl, dstr));

    // About the orientation of the chamfer plane
    // Compute the normal to the face 1
    let norpl = pos1.vxdir.cross(pos1.vydir).normalize();
    let mut norface1 = norpl;
    if of1 == Orientation::Reversed {
        norface1 = -norface1;
    }

    // Compute the orientation of the chamfer plane
    let norplch = gp_dir(planax3.vxdir.cross(planax3.vydir));

    let dirch12 = gp_dir(p2 - p1);
    let mut toreverse = norplch.dot(norface1) <= 0.0;
    if vectransl1.dot(dirch12) > 0.0 {
        toreverse = !toreverse;
    }

    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // Loading of the FaceInterferences with pcurves & 3d curves.

    // case face 1
    let linpln1 = Line3 {
        origin: p1,
        direction: xdir,
    };
    let glinpln1 = Curve3::Line(linpln1);

    let (u1, v1p) = elslib_plane_parameters(&pos1, p1);
    let u = u1;
    let v = v1p;
    let p2dpln = DVec2::new(u, v);
    let dir2dpln = DVec2::new(xdir.dot(pos1.vxdir), xdir.dot(pos1.vydir));
    let glin2dpln1 = curve2d_line(p2dpln, dir2dpln);

    let (u2, v2p) = elslib_plane_parameters(&planax3, p1);
    let u = u2;
    let v = v2p;
    let glin2dplnch1 = curve2d_line(DVec2::new(u, v), DVec2::new(1.0, 0.0));

    let trans;
    let mut toreverse = norplch.dot(norpl) <= 0.0;
    if vectransl1.dot(dirch12) > 0.0 {
        toreverse = !toreverse;
    }
    if toreverse {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }

    data.change_interference_on_s1().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln1, dstr),
        trans,
        Some(glin2dpln1),
        Some(glin2dplnch1),
    );

    // case face 2

    let linpln2 = Line3 {
        origin: p2,
        direction: xdir,
    };
    let glinpln2 = Curve3::Line(linpln2);

    let (u3, v3p) = elslib_plane_parameters(&pos2, p2);
    let u = u3;
    let v = v3p;
    let dir2dpln = DVec2::new(xdir.dot(pos2.vxdir), xdir.dot(pos2.vydir));
    let glin2dpln2 = curve2d_line(DVec2::new(u, v), dir2dpln);

    let (u4, v4p) = elslib_plane_parameters(&planax3, p2);
    let u = u4;
    let v = v4p;
    let glin2dplnch2 = curve2d_line(DVec2::new(u, v), DVec2::new(1.0, 0.0));

    let norpl = pos2.vxdir.cross(pos2.vydir).normalize();
    let mut toreverse = norplch.dot(norpl) <= 0.0;
    if vectransl2.dot(dirch12) < 0.0 {
        toreverse = !toreverse;
    }
    let trans;
    if toreverse {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }

    data.change_interference_on_s2().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln2, dstr),
        trans,
        Some(glin2dpln2),
        Some(glin2dplnch2),
    );

    true
}

// =========================================================================
// OCCT ChFiKPart_ComputeData_ChAsymPlnPln.cxx L53-266 — ChFiKPart_MakeChAsym
// (plane/plane, distance + angle chamfer).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi_kpart_make_ch_asym_pln_pln(
    dstr: &mut TopOpeBRepDSHDataStructure,
    data: &mut ChFiDSSurfData,
    pl1: &Plane,
    pl2: &Plane,
    or1: Orientation,
    or2: Orientation,
    dis: f64,
    angle: f64,
    spine: &Line3,
    first: f64,
    of1: Orientation,
    dis_on_p1: bool,
) -> bool {
    // Creation of the plane which carry the chamfer

    // compute the normals to the planes Pl1 and Pl2
    let pos1 = surface3_ax3(&Surface3::Plane(*pl1));
    let mut d1 = pos1.vxdir.cross(pos1.vydir).normalize();
    if or1 == Orientation::Reversed {
        d1 = -d1;
    }

    let pos2 = surface3_ax3(&Surface3::Plane(*pl2));
    let mut d2 = pos2.vxdir.cross(pos2.vydir).normalize();
    if or2 == Orientation::Reversed {
        d2 = -d2;
    }

    // compute the intersection line of Pl1 and Pl2
    let lint = rcad_kernel::base::int_ana::intersect_plane_plane_intana(pl1, pl2);

    let p;
    let fint;
    match &lint {
        rcad_kernel::base::int_ana::PlnPlnResult::Line(lint_line) => {
            let orspine = spine.origin + spine.direction * first;
            fint = (orspine - lint_line.origin).dot(lint_line.direction);
            p = lint_line.origin + lint_line.direction * fint;
        }
        _ => {
            return false;
        }
    }

    let linax1 = gp_dir(spine.direction);
    let mut vectransl1 = gp_dir(linax1.cross(d1));
    if vectransl1.dot(d2) < 0.0 {
        vectransl1 = -vectransl1;
    }

    let mut vectransl2 = gp_dir(linax1.cross(d2));
    if vectransl2.dot(d1) < 0.0 {
        vectransl2 = -vectransl2;
    }

    let cosp = vectransl1.dot(vectransl2);
    let sinp = (1.0 - cosp * cosp).sqrt();

    let dis1;
    let dis2;
    if dis_on_p1 {
        dis1 = dis;
        dis2 = dis / (cosp + sinp / angle.tan());
    } else {
        dis1 = dis / (cosp + sinp / angle.tan());
        dis2 = dis;
    }
    // Compute a point on the plane Pl1 and on the chamfer
    let p1 = p + vectransl1 * dis1;

    // Point on the plane Pl2 and on the chamfer
    let p2 = p + vectransl2 * dis2;

    // the middle point of P1 P2 is the origin of the chamfer
    let po = (p1 + p2) / 2.0;

    // compute a second point on the plane Pl2
    let pp = {
        let lint_line = match &lint {
            rcad_kernel::base::int_ana::PlnPlnResult::Line(l) => l,
            _ => unreachable!(),
        };
        lint_line.origin + lint_line.direction * (fint + 10.0)
    };
    let p22 = pp + vectransl2 * dis2;

    // Compute the normal vector <AxisPlan> to the chamfer's plane
    let v1 = gp_dir(p2 - p1);
    let v2 = gp_dir(p22 - p1);
    let axisplan = gp_dir(v1.cross(v2));

    let xdir = linax1; // u axis
    let mut planax3 = GpAx3::new_pn_vx(po, axisplan, xdir);
    if planax3.vydir.dot(d2) >= 0.0 {
        planax3.y_reverse();
    }

    // OCCT L157-158: Geom_Plane(PlanAx3); ChangeSurf.
    let gpl = surface3_plane(&planax3);
    data.change_surf(super::chfi_kpart::chfi_kpart_index_surface_in_ds(gpl, dstr));

    // About the orientation of the chamfer plane
    // Compute the normal to the face 1
    let norpl = pos1.vxdir.cross(pos1.vydir).normalize();
    let mut norface1 = norpl;
    if of1 == Orientation::Reversed {
        norface1 = -norface1;
    }

    // Compute the orientation of the chamfer plane
    let norplch = gp_dir(planax3.vxdir.cross(planax3.vydir));

    let dirch12 = gp_dir(p2 - p1);
    let mut toreverse = norplch.dot(norface1) <= 0.0;
    if vectransl1.dot(dirch12) > 0.0 {
        toreverse = !toreverse;
    }

    if toreverse {
        *data.change_orientation() = Orientation::Reversed;
    } else {
        *data.change_orientation() = Orientation::Forward;
    }

    // Loading of the FaceInterferences with pcurves & 3d curves.

    // case face 1
    let linpln1 = Line3 {
        origin: p1,
        direction: xdir,
    };
    let glinpln1 = Curve3::Line(linpln1);

    let (u1, v1p) = elslib_plane_parameters(&pos1, p1);
    let u = u1;
    let v = v1p;
    let p2dpln = DVec2::new(u, v);
    let dir2dpln = DVec2::new(xdir.dot(pos1.vxdir), xdir.dot(pos1.vydir));
    let glin2dpln1 = curve2d_line(p2dpln, dir2dpln);

    let (u2, v2p) = elslib_plane_parameters(&planax3, p1);
    let u = u2;
    let v = v2p;
    let glin2dplnch1 = curve2d_line(DVec2::new(u, v), DVec2::new(1.0, 0.0));

    let trans;
    let mut toreverse = norplch.dot(norpl) <= 0.0;
    if vectransl1.dot(dirch12) > 0.0 {
        toreverse = !toreverse;
    }
    if toreverse {
        trans = Orientation::Forward;
    } else {
        trans = Orientation::Reversed;
    }

    data.change_interference_on_s1().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln1, dstr),
        trans,
        Some(glin2dpln1),
        Some(glin2dplnch1),
    );

    // case face 2

    let linpln2 = Line3 {
        origin: p2,
        direction: xdir,
    };
    let glinpln2 = Curve3::Line(linpln2);

    let (u3, v3p) = elslib_plane_parameters(&pos2, p2);
    let u = u3;
    let v = v3p;
    let dir2dpln = DVec2::new(xdir.dot(pos2.vxdir), xdir.dot(pos2.vydir));
    let glin2dpln2 = curve2d_line(DVec2::new(u, v), dir2dpln);

    let (u4, v4p) = elslib_plane_parameters(&planax3, p2);
    let u = u4;
    let v = v4p;
    let glin2dplnch2 = curve2d_line(DVec2::new(u, v), DVec2::new(1.0, 0.0));

    let norpl = pos2.vxdir.cross(pos2.vydir).normalize();
    let mut toreverse = norplch.dot(norpl) <= 0.0;
    if vectransl2.dot(dirch12) < 0.0 {
        toreverse = !toreverse;
    }
    let trans;
    if toreverse {
        trans = Orientation::Reversed;
    } else {
        trans = Orientation::Forward;
    }

    data.change_interference_on_s2().set_interference(
        chfi_kpart_index_curve_in_ds(glinpln2, dstr),
        trans,
        Some(glin2dpln2),
        Some(glin2dplnch2),
    );

    true
}
