//! OCCT ChFi3d_Builder_C2.cxx (TKFillet/ChFi3d) — ChFi3d_Builder::
//! PerformTwoCornerbyInter (L108-1043) and the file statics Reduce
//! (L78-101), plus the ChFi3d_Builder_0.cxx helpers it consumes
//! (ChFi3d_InPeriod L92-108, ChFi3d_IntCS L4302-4385,
//! ChFi3d_ComputesIntPC L4389-4487).
//!
//! Architecture mappings:
//!   - OCCT `GeomAdaptor_Surface` -> the fillet-module `GeomAdaptorSurface`
//!     (surface + trimmed UV bounds); the HInter entry takes the bop
//!     `BRepAdaptorSurface` built from the same surface + bounds.
//!   - OCCT `GeomAdaptor_Curve` -> the bop `BRepAdaptorCurve` (curve +
//!     trimmed parameter range; Load(C, f, l) = with_range).
//!   - OCCT `Adaptor3d_CurveOnSurface` (a pcurve image on a face) has no
//!     Curve3 carrier; its plane image is exact for every supported family
//!     (`Adaptor3d_CurveOnSurface::Value` = the surface image of the 2d
//!     point), and on other face kinds the HInter keeps its documented
//!     not-done path (the OCCT polyhedron branch over the same image).
//!   - `Extrema_ExtPC(P, C, Tol)` -> the real kernel engine
//!     `rcad_kernel::base::extrema_ext_pc::ExtremaExtPC`.
//!   - `GeomLib::ExtendCurveToPoint` -> the translated body in
//!     `chfi3d_geom_lib` (BSpline poles; the OCCT body concatenates a
//!     GeomConvert_CompCurveToBSplineCurve, so the extension narrows to the
//!     BSpline case).

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{
    BezierCurve3, BSplineCurve3, Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _,
    Ellipse3, Hyperbola3, Line3, Parabola3, Surface3, SurfaceEval as _,
};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape, TShape};

use super::chfi3d::{chfi3d_index_of_surf_data, chfi3d_index_point_in_ds, topabs_compose, topabs_reverse, ChFi3dBuilder};
use super::chfi3d_builder_0::{
    chfi3d_bound_fac, chfi3d_bound_surf, chfi3d_boite, chfi3d_compute_curves,
    chfi3d_enlarge_box_curve, chfi3d_enlarge_box_dstr, chfi3d_enlarge_box_edge_faces,
    chfi3d_enlarge_box_surf_pc, chfi3d_fil_curve_in_ds, chfi3d_fil_point_in_ds,
    chfi3d_set_point_tolerance, BndBox, BRepAdaptorSurface, GeomAdaptorSurface,
};
use super::chfi3d_builder_2::BRepAdaptorCurve2d;
use super::chfi3d_builder_cncrn::chfi3d_is_in_front;
use super::chfi3d_ds::TopOpeBRepDSCurve;
use super::chfi3d_geom_lib::geom_lib_extend_curve_to_point;
use super::chfi_ds::ChFiDSSurfData;

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L92-108 — ChFi3d_InPeriod.
// =========================================================================
pub(crate) fn chfi3d_in_period(u: f64, ufirst: f64, ulast: f64, eps: f64) -> f64 {
    // OCCT L95-96.
    let period = ulast - ufirst;
    let mut u = u;
    // OCCT L97-100.
    while eps < ufirst - u {
        u += period;
    }
    // OCCT L101-104.
    while eps > ulast - u {
        u -= period;
    }
    // OCCT L105-107.
    if u < ufirst {
        u = ufirst;
    }
    u
}

// =========================================================================
// OCCT ChFi3d_Builder_C2.cxx L78-92 — static Reduce (surface overload).
// =========================================================================
fn reduce_surf(p1: f64, p2: f64, hs1: &mut GeomAdaptorSurface, hs2: &mut GeomAdaptorSurface) {
    // OCCT L79-84: the bounds of the shared underlying surface.
    let surf = hs1.surface.clone();
    let (ud, uf, vd, vf) = {
        let b = surf.default_domain();
        (b[0], b[1], b[2], b[3])
    };
    let milmoins = 0.51 * vd + 0.49 * vf;
    let milplus = 0.49 * vd + 0.51 * vf;
    // OCCT L85-91.
    if p1 < p2 {
        hs1.load_bounded(surf.clone(), ud, uf, vd, milmoins);
        hs2.load_bounded(surf, ud, uf, milplus, vf);
    } else {
        hs1.load_bounded(surf.clone(), ud, uf, milplus, vf);
        hs2.load_bounded(surf, ud, uf, vd, milmoins);
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C2.cxx L94-101 — static Reduce (curve overload).  The
// OCCT GeomAdaptor_Curve carries (curve, first, last); the bounds triple is
// mutated in place (the bop BRepAdaptorCurve is the rcad carrier).
// =========================================================================
fn reduce_curv(p1: f64, p2: f64, first: &mut f64, last: &mut f64) {
    // OCCT L95-98.
    let f = *first;
    let l = *last;
    let milmoins = 0.51 * f + 0.49 * l;
    let milplus = 0.49 * f + 0.51 * l;
    // OCCT L99-100.
    if p1 < p2 {
        *first = f;
        *last = milmoins;
    } else {
        *first = milplus;
        *last = l;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L4302-4385 — ChFi3d_IntCS (fast calculation of
// the intersection curve surface).  The OCCT curve adaptor (GeomAdaptor_Curve
// bounds [c_first, c_last], or the restricted image of a CurveOnSurface) is
// passed as (curve, first, last); the surface adaptor keeps its UV window.
// =========================================================================
pub(crate) fn chfi3d_int_cs(
    s: &GeomAdaptorSurface,
    c_curve: &Curve3,
    c_first: f64,
    c_last: f64,
    p2ds: &mut DVec2,
    wc: &mut f64,
) -> bool {
    // OCCT L4310-4313.
    let uf = c_first;
    let ul = c_last;
    let u1 = s.first_u_parameter();
    let u2 = s.last_u_parameter();
    let v1 = s.first_v_parameter();
    let v2 = s.last_v_parameter();
    // OCCT L4315-4319.
    let keepfirst = *wc < -1.0e100;
    let keeplast = *wc > 1.0e100;
    let mut temp = 0.0f64;
    if keepfirst {
        temp = 1.0e100;
    }
    if keeplast {
        temp = -1.0e100;
    }
    // OCCT L4320-4321.
    let mut dist = 2.0e100;

    // OCCT L4322-4323: IntCurveSurface_HInter Intersection;
    // Intersection.Perform(C, S).
    let ba_curve = crate::bop::int_tools::bean_face_intersector::BRepAdaptorCurve::with_range(
        c_curve.clone(),
        uf,
        ul,
    );
    let ba_surface =
        crate::bop::int_tools::bean_face_intersector::BRepAdaptorSurface::with_uv_bounds(
            s.surface.clone(),
            [
                s.first_u_parameter(),
                s.last_u_parameter(),
                s.first_v_parameter(),
                s.last_v_parameter(),
            ],
        );
    let mut inter = crate::bop::int_tools::bean_face_intersector::IntCurveSurfaceHInter::new();
    inter.perform(&ba_curve, &ba_surface);

    // OCCT L4324-4327.
    let mut isol: usize = 0;
    if inter.is_done() {
        // OCCT L4329-4350.
        for i in 1..=inter.nb_points() {
            let pint = inter.point(i);
            let mut up = pint.u();
            let mut vp = pint.v();
            if s.is_u_periodic() {
                up = chfi3d_in_period(up, u1, u1 + s.u_period(), 1.0e-8);
            }
            if s.is_v_periodic() {
                vp = chfi3d_in_period(vp, v1, v1 + s.v_period(), 1.0e-8);
            }
            if uf <= pint.w() && ul >= pint.w() && u1 <= up && u2 >= up && v1 <= vp && v2 >= vp {
                if keepfirst && pint.w() < temp {
                    temp = pint.w();
                    isol = i;
                } else if keeplast && pint.w() > temp {
                    temp = pint.w();
                    isol = i;
                } else if (pint.w() - *wc).abs() < dist {
                    dist = (pint.w() - *wc).abs();
                    isol = i;
                }
            }
        }
        // OCCT L4352-4355.
        if isol == 0 {
            return false;
        }
        // OCCT L4356-4368.
        let pint = inter.point(isol);
        let mut up = pint.u();
        let mut vp = pint.v();
        if s.is_u_periodic() {
            up = chfi3d_in_period(up, u1, u1 + s.u_period(), 1.0e-8);
        }
        if s.is_v_periodic() {
            vp = chfi3d_in_period(vp, v1, v1 + s.v_period(), 1.0e-8);
        }
        *p2ds = DVec2::new(up, vp);
        *wc = pint.w();
        return true;
    }
    false
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L4389-4399 — ChFi3d_ComputesIntPC (5-arg
// overload; marshals to the 6-arg one with a bidon point).  The
// FaceInterference pairs are marshaled as (pcurve_on_surf, first, last).
// =========================================================================
pub(crate) fn chfi3d_computes_int_pc(
    fi1: (&Curve2d, f64, f64),
    hs1: &GeomAdaptorSurface,
    fi2: (&Curve2d, f64, f64),
    hs2: &GeomAdaptorSurface,
    u_int1: &mut f64,
    u_int2: &mut f64,
) {
    // OCCT L4391-4397: gp_Pnt bid; ChFi3d_ComputesIntPC(Fi1, Fi2, HS1, HS2,
    // UInt1, UInt2, bid).
    let mut bid = glam::DVec3::ZERO;
    chfi3d_computes_int_pc_p(fi1, hs1, fi2, hs2, u_int1, u_int2, &mut bid);
}

/// OCCT ChFi3d_Builder_0.cxx L4402-4487 — ChFi3d_ComputesIntPC (6-arg).
#[allow(clippy::too_many_arguments)]
pub(crate) fn chfi3d_computes_int_pc_p(
    fi1: (&Curve2d, f64, f64),
    hs1: &GeomAdaptorSurface,
    fi2: (&Curve2d, f64, f64),
    hs2: &GeomAdaptorSurface,
    u_int1: &mut f64,
    u_int2: &mut f64,
    p: &mut glam::DVec3,
) {
    let (pc1, f1_first, f1_last) = fi1;
    let (pc2, f2_first, f2_last) = fi2;
    // OCCT L4412-4417: distref2 from the two seeded pcurve points.
    let xy = pc1.point_at(*u_int1);
    let p3d1 = hs1.value(xy.x, xy.y);
    let xy = pc2.point_at(*u_int2);
    let p3d2 = hs2.value(xy.x, xy.y);
    let distref2 = p3d1.distance_squared(p3d2);
    // OCCT L4418-4419.
    *p = 0.5 * (p3d1 + p3d2);
    // OCCT L4421-4423: delt1 = Min(0.1, 0.05*(Last-First)); hc2d1 restricted
    // to [UInt1-delt1, UInt1+delt1]; cons1(hc2d1, HS1).
    let delt1 = 0.05_f64 * (f1_last - f1_first);
    let delt1 = 0.1_f64.min(delt1);
    let im1 = plane_image_or_none(pc1, hs1);
    // OCCT L4425-4428: delt2 and cons2.
    let delt2 = 0.05_f64 * (f2_last - f2_first);
    let delt2 = 0.1_f64.min(delt2);
    let im2 = plane_image_or_none(pc2, hs2);
    // OCCT L4429-4447: Extrema_LocateExtCC ext(cons1, cons2, UInt1, UInt2).
    if let (Some(i1), Some(i2)) = (im1, im2) {
        let cons1 = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(i1),
            first: *u_int1 - delt1,
            last: *u_int1 + delt1,
        });
        let cons2 = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(i2),
            first: *u_int2 - delt2,
            last: *u_int2 + delt2,
        });
        if let Some((par1, par2, pnt1, pnt2)) =
            rcad_kernel::base::extrema::extrema_locate_ext_cc(&cons1, &cons2, *u_int1, *u_int2)
        {
            let dist2 = pnt1.distance_squared(pnt2);
            if dist2 < distref2 {
                *u_int1 = par1;
                *u_int2 = par2;
                *p = 0.5 * (pnt1 + pnt2);
            }
        }
    }
}

/// The 3d image of an interference pcurve on its carrying surface — the
/// `Adaptor3d_CurveOnSurface` value semantics.  On a plane the image of
/// every supported 2d family is the same canonic family through the plane
/// placement; on other surfaces the image is not a canonic Curve3 and None
/// keeps the OCCT downstream failure path.
fn plane_image_or_none(pc: &Curve2d, hs: &GeomAdaptorSurface) -> Option<Curve3> {
    let pln = match &hs.surface {
        Surface3::Plane(p) => *p,
        _ => return None,
    };
    let p3 = |p: DVec2| pln.point_at(p.x, p.y);
    let d3 = |d: DVec2| pln.u_dir * d.x + pln.v_dir * d.y;
    match pc {
        Curve2d::Line(l) => Some(Curve3::Line(Line3::new(p3(l.origin), d3(l.direction)))),
        Curve2d::Circle(c2) => Some(Curve3::Circle(Circle3 {
            center: p3(c2.center),
            normal: pln.normal,
            x_dir: d3(c2.x_dir),
            y_dir: d3(c2.y_dir),
            radius: c2.radius,
        })),
        Curve2d::Ellipse(e2) => Some(Curve3::Ellipse(Ellipse3 {
            center: p3(e2.center),
            normal: pln.normal,
            major_dir: d3(e2.major_dir),
            major_radius: e2.major_radius,
            minor_radius: e2.minor_radius,
        })),
        Curve2d::Hyperbola(h2) => Some(Curve3::Hyperbola(Hyperbola3 {
            center: p3(h2.center),
            normal: pln.normal,
            major_dir: d3(h2.major_dir),
            semi_major: h2.semi_major,
            semi_minor: h2.semi_minor,
        })),
        Curve2d::Parabola(p2) => Some(Curve3::Parabola(Parabola3 {
            vertex: p3(p2.origin),
            normal: pln.normal,
            axis_dir: d3(p2.axis_dir),
            focal_param: p2.focal_param,
        })),
        Curve2d::Bezier(b) => Some(Curve3::Bezier(BezierCurve3 {
            control_points: b.control_points.iter().map(|p| p3(*p)).collect(),
            weights: b.weights.clone(),
        })),
        Curve2d::BSpline(b) => Some(Curve3::BSpline(BSplineCurve3 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b.control_points.iter().map(|p| p3(*p)).collect(),
            weights: b.weights.clone(),
            is_periodic: false,
        })),
        _ => None,
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C2.cxx L108-1043 — ChFi3d_Builder::
// PerformTwoCornerbyInter.  Performs PerformTwoCorner by intersection.  In
// case of Biseau for all cases the path is used; 3D curve and 2 pcurves are
// approximated.
// =========================================================================
impl ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_C2.cxx L108 — the member entry (this).
    #[allow(clippy::too_many_lines)]
    pub fn perform_two_cornerby_inter(&mut self, index: usize) -> bool {
        let fb = self;
    let brep = fb.my_brep.clone();
    // OCCT L110-112.
    let mut done = false;
    let vtx = fb.my_vdata_map.find_key(index).clone();
    let mut dstr = fb.my_ds.take().expect("DS");

    // Information on fillets is extracted
    //------------------------------------------------------

    // the first
    //---------- (OCCT L115-124)
    let vdata: Vec<super::chfi_ds::SharedStripe> = fb.my_vdata_map.find_from_index(index).clone();
    let corner1 = vdata[0].clone();
    let mut sens1 = 0i32;
    let ifd1 = {
        let g = corner1.read().expect("stripe lock");
        chfi3d_index_of_surf_data(&vtx, &g, &mut sens1)
    };
    let seqfil1: Vec<super::chfi_ds::SharedSurfData> = {
        let g = corner1.read().expect("stripe lock");
        g.set_of_surf_data().clone()
    };
    let fd1 = seqfil1[(ifd1 - 1) as usize].clone();

    // the second
    //---------- (OCCT L126-140)
    let corner2 = vdata[1].clone();
    let mut sens2 = 0i32;
    let ifd2 = if Arc::ptr_eq(&corner2, &corner1) {
        // OCCT L133-136.
        let g = corner2.read().expect("stripe lock");
        sens2 = -1;
        g.set_of_surf_data().len() as i32
    } else {
        let g = corner2.read().expect("stripe lock");
        chfi3d_index_of_surf_data(&vtx, &g, &mut sens2)
    };
    let seqfil2: Vec<super::chfi_ds::SharedSurfData> = {
        let g = corner2.read().expect("stripe lock");
        g.set_of_surf_data().clone()
    };
    let fd2 = seqfil2[(ifd2 - 1) as usize].clone();

    // The concavities are analysed in case of different concavities,
    // preview an evolutionary connection of type ThreeCorner of R to 0.
    // Otherwise the opposite face
    // and the eventual intersection of 2 pcurves on this face are found.
    // (OCCT L142-170)
    let isfirst1 = sens1 == 1;
    let isfirst2 = sens2 == 1;
    let mut okintercc;
    let mut okvisavis = false;
    let mut same_side = false;
    let mut ifa_co1 = 0i32;
    let mut ifa_co2 = 0i32;
    let mut u_intpc1 = 0.0f64;
    let mut u_intpc2 = 0.0f64;
    let mut faco = Shape::null();
    okintercc = chfi3d_is_in_front(
        &brep, &dstr, &corner1, &corner2, ifd1, ifd2, sens1, sens2, &mut u_intpc1, &mut u_intpc2,
        &mut faco, &mut same_side, &mut ifa_co1, &mut ifa_co2, &mut okvisavis, &vtx, true, false,
    );
    if !okvisavis {
        // OCCT L171-181: "TwoCorner : pas de face commune".
        done = false;
        fb.my_ds = Some(dstr);
        return done;
    }
    if !okintercc {
        // OCCT L183-206: the intersection of pcurves is calculated without
        // restricting them by common points (enlarge = true).
        okintercc = chfi3d_is_in_front(
            &brep, &dstr, &corner1, &corner2, ifd1, ifd2, sens1, sens2, &mut u_intpc1,
            &mut u_intpc2, &mut faco, &mut same_side, &mut ifa_co1, &mut ifa_co2, &mut okvisavis,
            &vtx, true, true,
        );
    }

    if !okvisavis {
        // OCCT L208-215: "TwoCorner : no common face".
        done = false;
        fb.my_ds = Some(dstr);
        return done;
    }
    if !okintercc {
        // OCCT L216-225: "biseau : failed intersection of tangency lines on
        // common face".
        done = false;
        fb.my_ds = Some(dstr);
        return done;
    }
    // OCCT L227.
    let ifa_arc1 = 3 - ifa_co1;
    let ifa_arc2 = 3 - ifa_co2;

    // It is checked if the fillets have a commonpoint on a common arc.
    // This edge is the pivot of the bevel or of the kneecap.
    // (OCCT L230-247)
    let (cp1, cp2): (super::chfi_ds::ChFiDS_CommonPoint, super::chfi_ds::ChFiDS_CommonPoint) = {
        let mut fd1g = fd1.write().expect("surfdata lock");
        let mut fd2g = fd2.write().expect("surfdata lock");
        (
            fd1g.change_vertex(isfirst1, ifa_arc1).clone(),
            fd2g.change_vertex(isfirst2, ifa_arc2).clone(),
        )
    };
    if !cp1.is_on_arc() || !cp2.is_on_arc() {
        // OCCT L234-242: "fail 1 of 2 fillets are not on arc".
        done = false;
        fb.my_ds = Some(dstr);
        return done;
    }
    if !cp1.arc().is_same(cp2.arc()) {
        // OCCT L243-256: look like OnSame + OnDiff case (eap, Arp 9 2002,
        // occ266).  OCCT reads the same myDS reference; the rcad take/put-back
        // hands the DS over for the nested call.
        done = true;
        fb.my_ds = Some(dstr);
        fb.perform_more_three_corner(index, 2);
        return done;
    }

    // OCCT L258-276: pivot + parCP1/parCP2 + Hpivot + the ExtendCurveToPoint
    // block (unreachable after the L242 early return — translated for form).
    let pivot = cp1.arc().clone();
    let parcp1 = cp1.parameter_on_arc();
    let mut parcp2 = cp2.parameter_on_arc();
    let (hpivot, _, _) =
        crate::brep_algo::tool::brep_tool_curve(&pivot).expect("Standard_NullObject: pivot curve");
    if !pivot.is_same(cp2.arc()) {
        // OCCT L259-261: csau = BRep_Tool::Curve(pivot, ubid, vbid).
        if let Some((mut c1, _, _)) = crate::brep_algo::tool::brep_tool_curve(&pivot) {
            // OCCT L262-263: C1 = down_cast<Geom_BoundedCurve>(csau) — every
            // rcad Curve3 is bounded, so the downcast mirrors the null test.
            // OCCT L265: GeomLib::ExtendCurveToPoint(C1, CP2.Point(), 1,
            // false) — the translated body extends the BSpline poles (the
            // OCCT body concatenates through
            // GeomConvert_CompCurveToBSplineCurve).
            if let Curve3::BSpline(bs) = &mut c1 {
                geom_lib_extend_curve_to_point(bs, cp2.point(), 1, false);
            }
            // OCCT L264/L266-267: cad.Load(C1); Extrema_ExtPC ext(CP2.Point(),
            // cad, 1.e-4); parCP2 = ext.Point(1).Parameter().
            let a_adaptor = GeomCurveAdaptor::new(c1.clone());
            let a_tool = CurveToolHandle::for_curve3(&c1, &a_adaptor, &a_adaptor);
            let ext = ExtremaExtPC::new_point_curve(cp2.point(), &a_tool, 1.0e-4);
            parcp2 = ext.point(1).param;
        }
    }
    // OCCT L268-270.
    let psp1 = hpivot.point_at(parcp1);
    let psp2 = hpivot.point_at(parcp2);
    let sameparam = psp1.distance(psp2) < 10.0 * fb.tolapp3d;

    // OCCT L272-296: FF1/FF2 + the myEFMap ok1/ok2 check.
    let ff1 = dstr
        .shape(fd1.read().expect("surfdata lock").index_of(ifa_arc1))
        .clone();
    let ff2 = dstr
        .shape(fd2.read().expect("surfdata lock").index_of(ifa_arc2))
        .clone();
    {
        let mut ok1 = false;
        let mut ok2 = false;
        let faces = if fb.my_ef_map.contains(&pivot) {
            fb.my_ef_map.find(&pivot).clone()
        } else {
            Vec::new()
        };
        for kt in &faces {
            if !ok1 && ff1.is_same(kt) {
                ok1 = true;
            }
            if !ok2 && ff2.is_same(kt) {
                ok2 = true;
            }
        }
        if !ok1 || !ok2 {
            // OCCT L289-295: "fail one of surfaces has no common base face
            // with the pivot edge".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        }
    }

    // OCCT L298-301.
    let mut hs1 = chfi3d_bound_surf(&dstr, &fd1.read().expect("surfdata lock"), ifa_co1, ifa_arc1);
    let mut hs2 = chfi3d_bound_surf(&dstr, &fd2.read().expect("surfdata lock"), ifa_co2, ifa_arc2);

    // OCCT L303-309: Pardeb/Parfin(1,4) + PGc1/PGc2/Gc.
    let mut pardeb = [0.0f64; 4];
    let mut parfin = [0.0f64; 4];

    if sameparam {
        // Side common face, calculation of Pardeb. (OCCT L313-326)
        {
            let fd1g = fd1.read().expect("surfdata lock");
            let fd2g = fd2.read().expect("surfdata lock");
            chfi3d_computes_int_pc(
                interference_bundle(&fd1g, ifa_co1),
                &hs1,
                interference_bundle(&fd2g, ifa_co2),
                &hs2,
                &mut u_intpc1,
                &mut u_intpc2,
            );
        }
        {
            let fd1g = fd1.read().expect("surfdata lock");
            let uv = fd1g
                .interference(ifa_co1)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intpc1);
            pardeb[0] = uv.x;
            pardeb[1] = uv.y;
        }
        {
            let fd2g = fd2.read().expect("surfdata lock");
            let uv = fd2g
                .interference(ifa_co2)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intpc2);
            pardeb[2] = uv.x;
            pardeb[3] = uv.y;
        }
        // OCCT L325-326.
        let pfa_co = hs1.value(pardeb[0], pardeb[1]);

        // Side arc, calculation of Parfin. (OCCT L328-344)
        let mut u_intarc1 = {
            let fd1g = fd1.read().expect("surfdata lock");
            fd1g.interference(ifa_arc1).parameter(isfirst1)
        };
        let mut u_intarc2 = {
            let fd2g = fd2.read().expect("surfdata lock");
            fd2g.interference(ifa_arc2).parameter(isfirst2)
        };
        {
            let fd1g = fd1.read().expect("surfdata lock");
            let fd2g = fd2.read().expect("surfdata lock");
            chfi3d_computes_int_pc(
                interference_bundle(&fd1g, ifa_arc1),
                &hs1,
                interference_bundle(&fd2g, ifa_arc2),
                &hs2,
                &mut u_intarc1,
                &mut u_intarc2,
            );
        }
        {
            let fd1g = fd1.read().expect("surfdata lock");
            let uv = fd1g
                .interference(ifa_arc1)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intarc1);
            parfin[0] = uv.x;
            parfin[1] = uv.y;
        }
        {
            let fd2g = fd2.read().expect("surfdata lock");
            let uv = fd2g
                .interference(ifa_arc2)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intarc2);
            parfin[2] = uv.x;
            parfin[3] = uv.y;
        }

        // OCCT L346-348.
        let same_surf =
            fd1.read().expect("surfdata lock").surf() == fd2.read().expect("surfdata lock").surf();
        if same_surf {
            reduce_surf(u_intpc1, u_intpc2, &mut hs1, &mut hs2);
        }

        // OCCT L350-380: tolreached + ChFi3d_ComputeCurves on the ordered
        // parameter bundles.
        let mut tolreached = fb.tolapp3d;
        let computed = if ifa_co1 == 1 {
            chfi3d_compute_curves(&hs1, &hs2, pardeb, parfin, fb.tolapp3d, fb.tol2d, &mut tolreached)
        } else if ifa_co1 == 2 {
            chfi3d_compute_curves(&hs1, &hs2, parfin, pardeb, fb.tolapp3d, fb.tol2d, &mut tolreached)
        } else {
            // OCCT: neither branch runs and the null Gc would raise on use;
            // the not-done return keeps the OCCT failure outcome.
            None
        };
        let Some(cc) = computed else {
            // OCCT L355-365/L366-380: "failed to calculate bevel error
            // interSS".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        };
        let gc = cc.c3d.clone();
        let pgc1 = cc.pc1.clone();
        let pgc2 = cc.pc2.clone();

        // CornerData are updated with results of the intersection.
        // (OCCT L382-430)
        let (wfirst, wlast) = {
            let d = gc.default_domain();
            (d[0], d[1])
        };
        // OCCT L386-390: tolerances from the pre-mutation CommonPoints.
        let tolpco = {
            let fd1g = fd1.read().expect("surfdata lock");
            let fd2g = fd2.read().expect("surfdata lock");
            let cpco1 = fd1g.vertex(isfirst1, ifa_co1);
            let cpco2 = fd2g.vertex(isfirst2, ifa_co2);
            cpco1.tolerance().max(cpco2.tolerance())
        };
        let mut tolparc = {
            let fd1g = fd1.read().expect("surfdata lock");
            let fd2g = fd2.read().expect("surfdata lock");
            let cparc1 = fd1g.vertex(isfirst1, ifa_arc1);
            let cparc2 = fd2g.vertex(isfirst2, ifa_arc2);
            cparc1.tolerance().max(cparc2.tolerance())
        };
        // OCCT L391.
        let icurv =
            dstr.add_curve(TopOpeBRepDSCurve::new(Some(gc.clone()), tolreached));
        // OCCT L392-395: Corner1.
        {
            let mut c1g = corner1.write().expect("stripe lock");
            c1g.set_parameters(isfirst1, wfirst, wlast);
            c1g.set_curve(icurv, isfirst1);
            c1g.change_pcurve(isfirst1, pgc1.clone());
        }
        {
            let mut fd1g = fd1.write().expect("surfdata lock");
            let cpco1 = fd1g.change_vertex(isfirst1, ifa_co1);
            cpco1.reset();
            cpco1.set_point(pfa_co);
            cpco1.set_tolerance(tolreached.max(tolpco));
            fd1g.change_interference(ifa_co1)
                .set_parameter(isfirst1, u_intpc1);
            tolparc = tolparc.max(tolreached);
            let cparc1 = fd1g.change_vertex(isfirst1, ifa_arc1);
            cparc1.set_tolerance(tolparc.max(tolreached));
        }
        // OCCT L397-400.
        let ipoin1 = {
            let fd1g = fd1.read().expect("surfdata lock");
            chfi3d_index_point_in_ds(fd1g.vertex(isfirst1, 1), &mut dstr)
        };
        corner1
            .write()
            .expect("stripe lock")
            .set_index_point(ipoin1, isfirst1, 1);
        let ipoin2 = {
            let fd1g = fd1.read().expect("surfdata lock");
            chfi3d_index_point_in_ds(fd1g.vertex(isfirst1, 2), &mut dstr)
        };
        corner1
            .write()
            .expect("stripe lock")
            .set_index_point(ipoin2, isfirst1, 2);
        // OCCT L402-405: Corner2.
        {
            let mut c2g = corner2.write().expect("stripe lock");
            c2g.set_parameters(isfirst2, wfirst, wlast);
            c2g.set_curve(icurv, isfirst2);
            c2g.change_pcurve(isfirst2, pgc2.clone());
        }
        {
            let mut fd2g = fd2.write().expect("surfdata lock");
            fd2g.change_interference(ifa_co2)
                .set_parameter(isfirst2, u_intpc2);
            let cpco1_clone = fd1
                .read()
                .expect("surfdata lock")
                .vertex(isfirst1, ifa_co1)
                .clone();
            let cparc1_clone = fd1
                .read()
                .expect("surfdata lock")
                .vertex(isfirst1, ifa_arc1)
                .clone();
            *fd2g.change_vertex(isfirst2, ifa_co2) = cpco1_clone;
            *fd2g.change_vertex(isfirst2, ifa_arc2) = cparc1_clone;
        }
        // OCCT L406-412.
        if ifa_co1 != ifa_co2 {
            corner2
                .write()
                .expect("stripe lock")
                .set_orientation(Orientation::Reversed, isfirst2);
        }
        {
            let idx_co = corner1
                .read()
                .expect("stripe lock")
                .index_point(isfirst1, ifa_co1);
            corner2
                .write()
                .expect("stripe lock")
                .set_index_point(idx_co, isfirst2, ifa_co2);
            let idx_arc = corner1
                .read()
                .expect("stripe lock")
                .index_point(isfirst1, ifa_arc1);
            corner2
                .write()
                .expect("stripe lock")
                .set_index_point(idx_arc, isfirst2, ifa_arc2);
        }
        // The tolerances of points are updated. (OCCT L414-430)
        let mut bco = BndBox::default();
        let mut barc = BndBox::default();
        if ifa_co1 == 1 {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&corner1.read().expect("stripe lock")),
                &fd1.read().expect("surfdata lock"),
                &mut bco,
                &mut barc,
                isfirst1,
            );
        } else {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&corner1.read().expect("stripe lock")),
                &fd1.read().expect("surfdata lock"),
                &mut barc,
                &mut bco,
                isfirst1,
            );
        }
        if ifa_co2 == 1 {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&corner2.read().expect("stripe lock")),
                &fd2.read().expect("surfdata lock"),
                &mut bco,
                &mut barc,
                isfirst2,
            );
        } else {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&corner2.read().expect("stripe lock")),
                &fd2.read().expect("surfdata lock"),
                &mut barc,
                &mut bco,
                isfirst2,
            );
        }
        let cparc = fd1
            .read()
            .expect("surfdata lock")
            .vertex(isfirst1, ifa_arc1)
            .clone();
        {
            let lf = if fb.my_ef_map.contains(cparc.arc()) {
                fb.my_ef_map.find(cparc.arc()).clone()
            } else {
                Vec::new()
            };
            chfi3d_enlarge_box_edge_faces(
                &brep,
                cparc.arc(),
                &lf,
                cparc.parameter_on_arc(),
                &mut barc,
            );
        }
        chfi3d_set_point_tolerance(
            &mut dstr,
            &barc,
            corner1
                .read()
                .expect("stripe lock")
                .index_point(isfirst1, ifa_arc1),
        );
        chfi3d_set_point_tolerance(
            &mut dstr,
            &bco,
            corner1
                .read()
                .expect("stripe lock")
                .index_point(isfirst1, ifa_co1),
        );
    } else {
        // It is necessary to identify the border surface,
        // find the end point of the intersection Surf/Surf
        // by the intersection of the tangency line of the small
        // on the opposing face with the surface of the big,
        // and finally intersect the big with the face at end
        // between this point and the point on arc. (OCCT L433-474)
        let mut parcrois = false;
        {
            // OCCT L435-444: explorer over the FORWARD-oriented pivot's
            // vertices; parcrois = (matched orientation == TopAbs_FORWARD).
            let verts: [(Shape, Orientation); 2] = match pivot.data.as_ref() {
                TShape::Edge(ed) => [
                    (ed.first.clone(), Orientation::Forward),
                    (ed.last.clone(), Orientation::Reversed),
                ],
                _ => panic!("pivot is not an edge"),
            };
            for (cur, cur_o) in verts {
                if cur.is_same(&vtx) {
                    parcrois = cur_o == Orientation::Forward;
                    break;
                }
            }
        }
        let big_cd;
        let sma_cd;
        let big_fd;
        let sma_fd;
        let mut big_hs;
        let mut sma_hs;
        let ifa_co_big;
        let ifa_co_sma;
        let ifa_arc_big;
        let ifa_arc_sma;
        let isfirst_big;
        let isfirst_sma;
        let mut u_intpc_big;
        let mut u_intpc_sma;

        if (parcrois && parcp2 > parcp1) || (!parcrois && parcp2 < parcp1) {
            // OCCT L446-462.
            u_intpc_big = u_intpc2;
            u_intpc_sma = u_intpc1;
            big_hs = hs2.clone();
            sma_hs = hs1.clone();
            big_cd = corner2.clone();
            sma_cd = corner1.clone();
            big_fd = fd2.clone();
            sma_fd = fd1.clone();
            ifa_co_big = ifa_co2;
            ifa_co_sma = ifa_co1;
            ifa_arc_big = ifa_arc2;
            ifa_arc_sma = ifa_arc1;
            isfirst_big = isfirst2;
            isfirst_sma = isfirst1;
        } else {
            // OCCT L463-479.
            u_intpc_big = u_intpc1;
            u_intpc_sma = u_intpc2;
            big_hs = hs1.clone();
            sma_hs = hs2.clone();
            big_cd = corner1.clone();
            sma_cd = corner2.clone();
            big_fd = fd1.clone();
            sma_fd = fd2.clone();
            ifa_co_big = ifa_co1;
            ifa_co_sma = ifa_co2;
            ifa_arc_big = ifa_arc1;
            ifa_arc_sma = ifa_arc2;
            isfirst_big = isfirst1;
            isfirst_sma = isfirst2;
        }

        // Intersection of the big with the small :
        //------------------------------------

        // Pardeb (parameters of point PFaCo)
        // the intersection is checked (OCCT L484-504)
        {
            let smag = sma_fd.read().expect("surfdata lock");
            let bigg = big_fd.read().expect("surfdata lock");
            chfi3d_computes_int_pc(
                interference_bundle(&smag, ifa_co_sma),
                &sma_hs,
                interference_bundle(&bigg, ifa_co_big),
                &big_hs,
                &mut u_intpc_sma,
                &mut u_intpc_big,
            );
        }
        let uvi_big = {
            let bigg = big_fd.read().expect("surfdata lock");
            bigg
                .interference(ifa_co_big)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intpc_big)
        };
        pardeb[2] = uvi_big.x;
        pardeb[3] = uvi_big.y;
        let uvi_sma = {
            let smag = sma_fd.read().expect("surfdata lock");
            smag
                .interference(ifa_co_sma)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(u_intpc_sma)
        };
        pardeb[0] = uvi_sma.x;
        pardeb[1] = uvi_sma.y;
        let pfa_co = sma_hs.value(uvi_sma.x, uvi_sma.y);

        // Parfin (parameters of point PMil) (OCCT L506-535)
        let fi_arc_sma_first;
        let fi_arc_sma_last;
        {
            let smag = sma_fd.read().expect("surfdata lock");
            let fi_arc_sma = smag.interference(ifa_arc_sma);
            fi_arc_sma_first = fi_arc_sma.parameter_first();
            fi_arc_sma_last = fi_arc_sma.parameter_last();
        }
        let ctg = {
            let smag = sma_fd.read().expect("surfdata lock");
            dstr.curve(smag.interference(ifa_arc_sma).line_index())
                .curve
                .clone()
                .expect("DS curve null")
        };
        // OCCT L508-525: Hctg + wi/temp.
        let mut wi;
        let mut hctg_first;
        let mut hctg_last;
        if isfirst_sma {
            wi = fi_arc_sma_first;
            let mut temp_v = fi_arc_sma_first;
            if u_intpc_sma < temp_v {
                temp_v = u_intpc_sma;
            }
            hctg_first = temp_v;
            hctg_last = fi_arc_sma_last;
        } else {
            wi = fi_arc_sma_last;
            let mut temp_v = fi_arc_sma_last;
            if u_intpc_sma > temp_v {
                temp_v = u_intpc_sma;
            }
            hctg_first = fi_arc_sma_first;
            hctg_last = temp_v;
        }
        // OCCT L526-529: Reduce over the shared surface and the tangency
        // curve.
        let same_surf = sma_fd.read().expect("surfdata lock").surf()
            == big_fd.read().expect("surfdata lock").surf();
        if same_surf {
            reduce_surf(u_intpc_sma, u_intpc_big, &mut sma_hs, &mut big_hs);
            reduce_curv(u_intpc_sma, u_intpc_big, &mut hctg_first, &mut hctg_last);
        }
        // OCCT L531-538: ChFi3d_IntCS(BigHS, Hctg, UVi, wi).
        let mut uvi = DVec2::ZERO;
        if !chfi3d_int_cs(&big_hs, &ctg, hctg_first, hctg_last, &mut uvi, &mut wi) {
            // OCCT L533-538: "bevel : failed inter C S".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        }
        // OCCT L540-546.
        parfin[2] = uvi.x;
        parfin[3] = uvi.y;
        uvi = {
            let smag = sma_fd.read().expect("surfdata lock");
            smag
                .interference(ifa_arc_sma)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(wi)
        };
        parfin[0] = uvi.x;
        parfin[1] = uvi.y;
        let p_mil = sma_hs.value(parfin[0], parfin[1]);

        // OCCT L548-559: ChFi3d_ComputeCurves(SmaHS, BigHS, Pardeb, Parfin,
        // ...).
        let mut tolreached = 0.0f64;
        let computed =
            chfi3d_compute_curves(&sma_hs, &big_hs, pardeb, parfin, fb.tolapp3d, fb.tol2d, &mut tolreached);
        let Some(cc) = computed else {
            // OCCT L552-558: "failed to calculate bevel failed interSS".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        };
        let gc = cc.c3d.clone();
        let pgc1 = cc.pc1.clone();
        let pgc2 = cc.pc2.clone();

        // SmaCD is updated, for it this is all. (OCCT L561-589)
        let (wfirst, wlast) = {
            let d = gc.default_domain();
            (d[0], d[1])
        };
        let tolpco = {
            let smag = sma_fd.read().expect("surfdata lock");
            let bigg = big_fd.read().expect("surfdata lock");
            let psmaco = smag.vertex(isfirst_sma, ifa_co_sma);
            let pbigco = bigg.vertex(isfirst_big, ifa_co_big);
            psmaco.tolerance().max(pbigco.tolerance())
        };
        let tolpmil = {
            let smag = sma_fd.read().expect("surfdata lock");
            let psmamil = smag.vertex(isfirst_sma, ifa_arc_sma);
            psmamil.tolerance()
        };
        // OCCT L568.
        let icurv = dstr.add_curve(TopOpeBRepDSCurve::new(Some(gc.clone()), tolreached));
        {
            let mut scg = sma_cd.write().expect("stripe lock");
            scg.set_parameters(isfirst_sma, wfirst, wlast);
            scg.set_curve(icurv, isfirst_sma);
            scg.change_pcurve(isfirst_sma, pgc1.clone());
        }
        {
            let mut smag = sma_fd.write().expect("surfdata lock");
            let psmaco = smag.change_vertex(isfirst_sma, ifa_co_sma);
            psmaco.reset();
            psmaco.set_point(pfa_co);
            psmaco.set_tolerance(tolpco.max(tolreached));
            smag.change_interference(ifa_co_sma)
                .set_parameter(isfirst_sma, u_intpc_sma);
            let psmamil = smag.change_vertex(isfirst_sma, ifa_arc_sma);
            psmamil.reset();
            psmamil.set_point(p_mil);
            psmamil.set_tolerance(tolpmil.max(tolreached));
            smag.change_interference(ifa_arc_sma)
                .set_parameter(isfirst_sma, wi);
        }
        // OCCT L580-585.
        let ipoint_co = {
            let smag = sma_fd.read().expect("surfdata lock");
            let psmaco = smag.vertex(isfirst_sma, ifa_co_sma).clone();
            chfi3d_index_point_in_ds(&psmaco, &mut dstr)
        };
        sma_cd
            .write()
            .expect("stripe lock")
            .set_index_point(ipoint_co, isfirst_sma, ifa_co_sma);
        let ipoint_mil = {
            let smag = sma_fd.read().expect("surfdata lock");
            let psmamil = smag.vertex(isfirst_sma, ifa_arc_sma).clone();
            chfi3d_index_point_in_ds(&psmamil, &mut dstr)
        };
        sma_cd
            .write()
            .expect("stripe lock")
            .set_index_point(ipoint_mil, isfirst_sma, ifa_arc_sma);
        // OCCT L586-589.
        if ifa_co_sma == 2 {
            sma_cd
                .write()
                .expect("stripe lock")
                .set_orientation(Orientation::Reversed, isfirst_sma);
        }
        // For BigCD the first results are met in the DS. (OCCT L591-613)
        {
            big_cd
                .write()
                .expect("stripe lock")
                .set_index_point(ipoint_co, isfirst_big, ifa_co_big);
        }
        {
            let psmaco = sma_fd
                .read()
                .expect("surfdata lock")
                .vertex(isfirst_sma, ifa_co_sma)
                .clone();
            let mut bigg = big_fd.write().expect("surfdata lock");
            *bigg.change_vertex(isfirst_big, ifa_co_big) = psmaco;
            bigg.change_interference(ifa_co_big)
                .set_parameter(isfirst_big, u_intpc_big);
        }
        // OCCT L606-613: the point interferences on the curve.
        {
            let li = dstr.change_curve_interferences(icurv);
            li.push(chfi3d_fil_point_in_ds(
                Orientation::Forward,
                icurv,
                ipoint_co,
                wfirst,
                false,
            ));
            li.push(chfi3d_fil_point_in_ds(
                Orientation::Reversed,
                icurv,
                ipoint_mil,
                wlast,
                false,
            ));
        }
        // the transition of curves of intersection on the Big
        // (OCCT L615-632)
        let tra = big_fd
            .read()
            .expect("surfdata lock")
            .interference_on_s1()
            .transition();
        let ofac = dstr
            .shape(big_fd.read().expect("surfdata lock").index_of_s1)
            .orientation;
        let ofil = big_fd.read().expect("surfdata lock").orientation();
        let mut tracurv = topabs_compose(ofac, ofil);
        tracurv = topabs_compose(tracurv, tra);
        if !isfirst_big {
            tracurv = topabs_reverse(tracurv);
        }
        if ifa_co_big != 1 {
            tracurv = topabs_reverse(tracurv);
        }
        let isurf = big_fd.read().expect("surfdata lock").surf();
        let interfc = chfi3d_fil_curve_in_ds(icurv, isurf, Some(pgc2.clone()), tracurv);
        dstr.change_surface_interferences(isurf).push(interfc);

        // The tolerances of points are updated (beginning).
        // (OCCT L634-641)
        let mut bco = BndBox::default();
        let mut bmil = BndBox::default();
        let mut barc = BndBox::default();
        if ifa_co_sma == 1 {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&sma_cd.read().expect("stripe lock")),
                &sma_fd.read().expect("surfdata lock"),
                &mut bco,
                &mut bmil,
                isfirst_sma,
            );
        } else {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                Some(&sma_cd.read().expect("stripe lock")),
                &sma_fd.read().expect("surfdata lock"),
                &mut bmil,
                &mut bco,
                isfirst_sma,
            );
        }
        chfi3d_enlarge_box_surf_pc(&big_hs, &pgc2, wfirst, wlast, &mut bco, &mut bmil);

        // Intersection of the big with the face at end :
        // -------------------------------------------

        // Pardeb (parameters of PMil)
        // The intersection curve surface is tried again, now with
        // representation pcurve on face of the curve to be sure.
        // (OCCT L648-668)
        let f = dstr
            .shape(sma_fd.read().expect("surfdata lock").index_of(ifa_arc_sma))
            .clone();
        let mut hf = BRepAdaptorSurface::initialize(&brep, &f);
        let (fsma_raw, lsma_raw) = {
            let smag = sma_fd.read().expect("surfdata lock");
            let fi_arc_sma = smag.interference(ifa_arc_sma);
            (fi_arc_sma.parameter_first(), fi_arc_sma.parameter_last())
        };
        let deltsma = 0.05 * (lsma_raw - fsma_raw);
        let pcpc = sma_fd
            .read()
            .expect("surfdata lock")
            .interference(ifa_arc_sma)
            .pcurve_on_face()
            .expect("PCurveOnFace")
            .clone();
        let (pcpc_first, pcpc_last) = {
            let d = pcpc.default_domain();
            (d[0], d[1])
        };
        let fsma = pcpc_first.max(wi - deltsma);
        let lsma = pcpc_last.min(wi + deltsma);
        if lsma < fsma {
            // OCCT L659-661.
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        }
        // OCCT L652-653: c2df over the restricted pcurve; consf =
        // Adaptor3d_CurveOnSurface(c2df, HF); Hconsf wrapper.  The rcad
        // HInter consumes the 3d image curve: on a plane face the image is
        // exact for every supported family; on other faces the HInter keeps
        // its documented not-done path (the failure branch below).
        let hconsf_ok = match plane_image_or_none(&pcpc, &hf_adaptor_geom(&hf)) {
            Some(image) => chfi3d_int_cs(&big_hs, &image, fsma, lsma, &mut uvi, &mut wi),
            None => false,
        };
        if !hconsf_ok {
            // OCCT L687-692: "bevel : failed inter C S".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        }
        // OCCT L663-668.
        pardeb[2] = uvi.x;
        pardeb[3] = uvi.y;
        uvi = pcpc.point_at(wi);
        pardeb[0] = uvi.x;
        pardeb[1] = uvi.y;
        let ppff1 = uvi;

        // Parfin (parameters of the point cpend) (OCCT L670-682)
        let ptg = big_fd
            .read()
            .expect("surfdata lock")
            .interference(ifa_arc_big)
            .parameter(isfirst_big);
        let uvi_big_arc = {
            let bigg = big_fd.read().expect("surfdata lock");
            bigg
                .interference(ifa_arc_big)
                .pcurve_on_surf()
                .expect("PCurveOnSurf")
                .point_at(ptg)
        };
        parfin[2] = uvi_big_arc.x;
        parfin[3] = uvi_big_arc.y;
        let cpend = {
            let bigg = big_fd.read().expect("surfdata lock");
            bigg.vertex(isfirst_big, ifa_arc_big).clone()
        };
        let mut etest = cpend.arc().clone();
        // OCCT L677-679: BRep_Tool::IsClosed(etest, F) -> etest.Reverse().
        if brep.is_edge_closed_on_face(&etest, &f) {
            etest.orientation = Orientation::Reversed;
        }
        let mut arc = BRepAdaptorCurve2d::new();
        arc.initialize(&brep, &etest, &f);
        let uvi_arc = arc.value(cpend.parameter_on_arc());
        parfin[0] = uvi_arc.x;
        parfin[1] = uvi_arc.y;
        let ppff2 = uvi_arc;

        // Intersection. (OCCT L684-693)
        let (uu1, uu2, vv1, vv2) = chfi3d_boite(ppff1, ppff2);
        // for the case when two chamfers are on two edges OnSame,
        // it is necessary to extend the surface carrying F, or at least
        // not to limit it.
        chfi3d_bound_fac(&mut hf, uu1, uu2, vv1, vv2, true);

        let computed =
            chfi3d_compute_curves(&hf_as_geom(&hf), &big_hs, pardeb, parfin, fb.tolapp3d, fb.tol2d, &mut tolreached);
        let Some(cc) = computed else {
            // OCCT L698-707: "fail calculation bevel fail interSS".
            done = false;
            fb.my_ds = Some(dstr);
            return done;
        };
        let gc = cc.c3d.clone();
        let pgc1 = cc.pc1.clone();
        let pgc2 = cc.pc2.clone();

        // End of update of the BigCD and the DS. (OCCT L709-727)
        let (wfirst, wlast) = {
            let d = gc.default_domain();
            (d[0], d[1])
        };
        let icurv = dstr.add_curve(TopOpeBRepDSCurve::new(Some(gc.clone()), tolreached));
        {
            let mut bigg = big_fd.write().expect("surfdata lock");
            let cpendg = bigg.change_vertex(isfirst_big, ifa_arc_big);
            let sav = cpendg.tolerance();
            cpendg.set_tolerance(sav.max(tolreached));
        }
        let cpend = {
            let bigg = big_fd.read().expect("surfdata lock");
            bigg.vertex(isfirst_big, ifa_arc_big).clone()
        };
        let ipoint_arc = chfi3d_index_point_in_ds(&cpend, &mut dstr);
        big_cd
            .write()
            .expect("stripe lock")
            .set_index_point(ipoint_arc, isfirst_big, ifa_arc_big);
        {
            let li = dstr.change_curve_interferences(icurv);
            li.push(chfi3d_fil_point_in_ds(
                Orientation::Forward,
                icurv,
                ipoint_mil,
                wfirst,
                false,
            ));
            li.push(chfi3d_fil_point_in_ds(
                Orientation::Reversed,
                icurv,
                ipoint_arc,
                wlast,
                false,
            ));
        }
        let interfc = chfi3d_fil_curve_in_ds(icurv, isurf, Some(pgc2.clone()), tracurv);
        dstr.change_surface_interferences(isurf).push(interfc);
        big_cd.write().expect("stripe lock").in_ds(isfirst_big, 1);

        // Finally the information on faces is placed in the DS.
        // (OCCT L729-745)
        let ishape = dstr.add_shape(&f);
        if sma_fd.read().expect("surfdata lock").surf()
            == big_fd.read().expect("surfdata lock").surf()
        {
            tracurv = topabs_compose(etest.orientation, cpend.transition_on_arc());
        } else {
            let f_edges = super::chfi3d_filbuilder_c2::face_edges(&brep, &f);
            for cur in &f_edges {
                if cur.is_same(&etest) {
                    tracurv = topabs_compose(cur.orientation, cpend.transition_on_arc());
                    break;
                }
            }
        }
        let interfc = chfi3d_fil_curve_in_ds(icurv, ishape, Some(pgc1.clone()), tracurv);
        dstr.change_shape_interferences(ishape).push(interfc);

        // The tolerances of points are updated (end). (OCCT L747-765)
        if ifa_co_big == 1 {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                None,
                &big_fd.read().expect("surfdata lock"),
                &mut bco,
                &mut barc,
                isfirst_big,
            );
        } else {
            chfi3d_enlarge_box_dstr(
                &brep,
                &dstr,
                None,
                &big_fd.read().expect("surfdata lock"),
                &mut barc,
                &mut bco,
                isfirst_big,
            );
        }
        chfi3d_enlarge_box_surf_pc(&big_hs, &pgc2, wfirst, wlast, &mut bmil, &mut barc);
        chfi3d_enlarge_box_surf_pc(&hf_as_geom(&hf), &pgc1, wfirst, wlast, &mut bmil, &mut barc);
        chfi3d_enlarge_box_curve(&gc, wfirst, wlast, &mut bmil, &mut barc);
        {
            let cparc = big_fd
                .read()
                .expect("surfdata lock")
                .vertex(isfirst_big, ifa_arc_big)
                .clone();
            let lf = if fb.my_ef_map.contains(cparc.arc()) {
                fb.my_ef_map.find(cparc.arc()).clone()
            } else {
                Vec::new()
            };
            chfi3d_enlarge_box_edge_faces(&brep, cparc.arc(), &lf, cparc.parameter_on_arc(), &mut barc);
        }
        chfi3d_set_point_tolerance(
            &mut dstr,
            &bco,
            sma_cd
                .read()
                .expect("stripe lock")
                .index_point(isfirst_sma, ifa_co_sma),
        );
        chfi3d_set_point_tolerance(
            &mut dstr,
            &bmil,
            sma_cd
                .read()
                .expect("stripe lock")
                .index_point(isfirst_sma, ifa_arc_sma),
        );
        chfi3d_set_point_tolerance(
            &mut dstr,
            &barc,
            big_cd
                .read()
                .expect("stripe lock")
                .index_point(isfirst_big, ifa_arc_big),
        );
    }
    // OCCT L767-769.
    done = true;
    fb.my_ds = Some(dstr);
    done
    }
}

/// The interference bundle (pcurve-on-surf + parameter range) marshaling a
/// ChFiDS_FaceInterference for ChFi3d_ComputesIntPC.
fn interference_bundle(fd: &ChFiDSSurfData, on_s: i32) -> (&Curve2d, f64, f64) {
    let fi = fd.interference(on_s);
    (
        fi.pcurve_on_surf().expect("PCurveOnSurf"),
        fi.parameter_first(),
        fi.parameter_last(),
    )
}

/// The GeomAdaptor view of a face adaptor (Value-only; OCCT passes the
/// BRepAdaptor_Surface HF itself where rcad needs the GeomAdaptor kind).
fn hf_adaptor_geom(hf: &BRepAdaptorSurface) -> GeomAdaptorSurface {
    GeomAdaptorSurface::new(hf.surface.clone())
}

/// The GeomAdaptor view keeping the (possibly BoundFac-adjusted) UV window.
fn hf_as_geom(hf: &BRepAdaptorSurface) -> GeomAdaptorSurface {
    let mut gs = GeomAdaptorSurface::new(hf.surface.clone());
    gs.load_bounded(hf.surface.clone(), hf.ufirst, hf.ulast, hf.vfirst, hf.vlast);
    gs
}
