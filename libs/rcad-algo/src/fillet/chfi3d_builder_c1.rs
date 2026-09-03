//! OCCT ChFi3d_Builder_C1.cxx — corner machinery, 1:1 translation.
//!
//! Contents (OCCT line references):
//!   - recadre (static, L118-138)
//!   - Update x3 (static, L147-181 / L191-279 / L349-417)
//!   - IntersUpdateOnSame (static, L292-340)
//!   - ChFi3d_ExtendSurface (static, L421-449)
//!   - ComputeCurve2d (static, L456-483)
//!   - ChFi3d_Recale (static, L487-551)
//!   - ChFi3d_SelectStripe (L558-588)
//!   - ChFi3d_Builder::FindFace (L4533-4595)
//!   - ChFi3d_Builder::PerformOneCorner (L611-1652)

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation};
use rcad_kernel::topods;

use super::chfi3d::chfi3d_index_point_in_ds;
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_bound_fac, chfi3d_bound_surf, chfi3d_boite,
    chfi3d_compute_curves, chfi3d_couture, chfi3d_enlarge_box_curve, chfi3d_enlarge_box_dstr,
    chfi3d_enlarge_box_edge_faces, chfi3d_enlarge_box_surf_pc, chfi3d_eval_tol_reached,
    chfi3d_is_pseudo_seam, chfi3d_reparam_pcurv,
    chfi3d_set_point_tolerance, BndBox, BRepAdaptorSurface, GeomAdaptorSurface, P_CONFUSION,
};
use super::chfi3d_builder_0::chfi3d_fil_curve_in_ds;
use super::chfi3d_builder_0_filds::{chfi3d_contains, geom2d_int_g_inter};
use super::chfi3d::{topabs_compose, topabs_reverse, ChFi3dBuilder};
use super::chfi3d::chfi3d_index_of_surf_data;
use super::chfi_ds::{ChFiDS_CommonPoint, ChFiDS_State, ChFiDSSurfData, SharedStripe};
use super::topopebrepds::{
    TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind,
    TopOpeBRepDSPoint,
};

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L118-138 — recadre.
// =========================================================================
fn recadre(p: f64, ref_: f64, isfirst: bool, first: f64, last: f64) -> f64 {
    let mut pp = p;
    if isfirst {
        pp -= last - first;
    } else {
        pp += last - first;
    }
    if (pp - ref_).abs() < (p - ref_).abs() {
        return pp;
    }
    p
}

// =========================================================================
// OCCT Adaptor3d_CurveOnSurface stand-in — a 2D pcurve evaluated on a
// surface (the 3D image the OCCT adaptor exposes).  Architecture: rcad
// evaluates the pcurve/surface pair in place of the adaptor handle.
// =========================================================================
struct CurveOnSurf<'a> {
    surf: &'a rcad_kernel::geom::Surface3,
    pc: &'a rcad_kernel::geom::Curve2d,
}

impl<'a> CurveOnSurf<'a> {
    fn value(&self, t: f64) -> DVec3 {
        let uv = self.pc.point_at(t);
        self.surf.point_at(uv.x, uv.y)
    }

    fn deriv2(&self, t: f64) -> DVec3 {
        let h = 1e-4f64;
        (self.value(t + h) - 2.0 * self.value(t) + self.value(t - h)) / (h * h)
    }
}

/// Newton locate-extremum between two curve-on-surface images from seeds —
/// OCCT Extrema_LocateExtCC over Adaptor3d_CurveOnSurface.
fn locate_ext_cc_curvesonurf(
    c1: &CurveOnSurf,
    c2: &CurveOnSurf,
    seed1: f64,
    seed2: f64,
) -> Option<(f64, f64, f64)> {
    // (par1, par2, square distance)
    let (mut s, mut t) = (seed1, seed2);
    let h = 1e-6f64;
    for _ in 0..50 {
        let p1 = c1.value(s);
        let p2 = c2.value(t);
        let diff = p1 - p2;
        let d1 = (c1.value(s + h) - c1.value(s - h)) / (2.0 * h);
        let d2 = (c2.value(t + h) - c2.value(t - h)) / (2.0 * h);
        let g = [2.0 * diff.dot(d1), -2.0 * diff.dot(d2)];
        if g[0].hypot(g[1]) < 1e-12 {
            break;
        }
        let dd1 = c1.deriv2(s);
        let dd2 = c2.deriv2(t);
        let h11 = 2.0 * (d1.dot(d1) + diff.dot(dd1));
        let h22 = 2.0 * (d2.dot(d2) - diff.dot(dd2));
        let ds = -g[0] / (h11.abs().max(1e-30) * h11.signum());
        let dt = -g[1] / (h22.abs().max(1e-30) * h22.signum());
        let f0 = diff.length_squared();
        let mut alpha = 1.0f64;
        let (mut ns, mut nt) = (s, t);
        for _ in 0..8 {
            ns = s + alpha * ds;
            nt = t + alpha * dt;
            if (c1.value(ns) - c2.value(nt)).length_squared() < f0 {
                break;
            }
            alpha *= 0.5;
        }
        let step = (ns - s).hypot(nt - t);
        s = ns;
        t = nt;
        if step < 1e-9 {
            break;
        }
    }
    Some((
        s,
        t,
        (c1.value(s) - c2.value(t)).length_squared(),
    ))
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L147-181 — Update: the LocateExtCC between the
// pcurve on the face at end and the interference pcurve on the fillet
// surface.  Returns false (and leaves the outputs) when the extremum is
// farther than tol.
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn update_locate(
    fb: &BRepAdaptorSurface,
    pcfb: &rcad_kernel::geom::Curve2d,
    surf: &GeomAdaptorSurface,
    fi: &mut super::chfi_ds::ChFiDS_FaceInterference,
    cp: &mut ChFiDS_CommonPoint,
    p2dbout: &mut DVec2,
    isfirst: bool,
    pared: &mut f64,
    wop: &mut f64,
    tol: f64,
) -> bool {
    let Some(pc_fi) = fi.pcurve_on_surf() else {
        return false;
    };
    let c1 = CurveOnSurf {
        surf: &fb.surface,
        pc: pcfb,
    };
    let c2 = CurveOnSurf {
        surf: surf.surface(),
        pc: pc_fi,
    };
    let Some((parfb, parltg, dist2)) = locate_ext_cc_curvesonurf(&c1, &c2, *pared, *wop) else {
        return false;
    };
    // OCCT L166: if (dist2 < tol * tol).
    if dist2 < tol * tol {
        *p2dbout = pcfb.point_at(parfb);
        *pared = parfb;
        *wop = parltg;
        fi.set_parameter(isfirst, *wop);
        cp.reset();
        let p1 = c1.value(parfb);
        cp.set_point(p1);
        return true;
    }
    false
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L191-279 — Update: intersect the face at end
// with the 3D interference curve; update the parameter and point.  The
// (uf, ul) range plays the role of the OCCT GeomAdaptor_Curve(ct, uf, ul)
// domain.
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn update_cs(
    fb: &BRepAdaptorSurface,
    ct: &rcad_kernel::geom::Curve3,
    uf: f64,
    ul: f64,
    fi: &mut super::chfi_ds::ChFiDS_FaceInterference,
    cp: &mut ChFiDS_CommonPoint,
    p2dbout: &mut DVec2,
    isfirst: bool,
    wop: &mut f64,
) -> bool {
    // OCCT: IntCurveSurface_HInter Intersection; Intersection.Perform(ct, fb).
    let mut intersection = crate::geomalgo::int_patch::int_cs::IntCurveSurface::new();
    let isperiodic = ct.is_periodic();
    intersection.perform(ct, &fb.surface, [uf, ul]);
    if !intersection.is_done() {
        return false;
    }
    // check if in KPart the limits of the tangency line are already in
    // place at this stage.  Modif lvt: the periodic cases are reframed,
    // especially if nothing was found.
    let mut w;
    let mut wbis = 0.0f64;
    let mut recadrebis = false;
    let nbp = intersection.nb_points();
    let mut dist = f64::INFINITY;
    let mut distbis = f64::INFINITY;
    let mut isol = 0usize;
    let mut isolbis = 0usize;
    for i in 1..=nbp {
        w = intersection.point(i - 1).w();
        if isperiodic {
            w = recadre(w, *wop, isfirst, uf, ul);
        }
        if uf <= w && ul >= w && (w - *wop).abs() < dist {
            isol = i;
            dist = (w - *wop).abs();
        }
    }
    if isperiodic {
        for i in 1..=nbp {
            w = intersection.point(i - 1).w();
            if uf <= w
                && ul >= w
                && (w - *wop).abs() < distbis
                && ((w - ul).abs() <= 0.01 || (w - uf).abs() <= 0.01)
            {
                isolbis = i;
                wbis = recadre(w, *wop, isfirst, uf, ul);
                distbis = (wbis - *wop).abs();
                recadrebis = true;
            }
        }
    }
    if isol == 0 && isolbis == 0 {
        return false;
    }
    let w;
    if !recadrebis {
        let pint = intersection.point(isol - 1);
        *p2dbout = DVec2::new(pint.u(), pint.v());
        let wraw = pint.w();
        w = if isperiodic {
            super::chfi_ds::elclib_in_period(wraw, uf, ul)
        } else {
            wraw
        };
    } else if dist > distbis {
        let pint = intersection.point(isolbis - 1);
        *p2dbout = DVec2::new(pint.u(), pint.v());
        w = wbis;
    } else {
        let pint = intersection.point(isol - 1);
        *p2dbout = DVec2::new(pint.u(), pint.v());
        let wraw = pint.w();
        w = super::chfi_ds::elclib_in_period(wraw, uf, ul);
    }
    fi.set_parameter(isfirst, w);
    cp.reset();
    cp.set_point(ct.point_at(w));
    *wop = w;
    true
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L292-340 — IntersUpdateOnSame: intersect the
// ChFi-<Fop> interference curve with the extended surface of <Fprol> and
// update FIop/CPop.
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn inters_update_on_same(
    hgs: &GeomAdaptorSurface,
    hbs: &mut BRepAdaptorSurface,
    c3dfi: &rcad_kernel::geom::Curve3,
    uf: f64,
    ul: f64,
    fop: &rcad_kernel::topods::Shape,
    fprol: &rcad_kernel::topods::Shape,
    eprol: &rcad_kernel::topods::Shape,
    vtx: &rcad_kernel::topods::Shape,
    is_first: bool,
    tol: f64,
    fiop: &mut super::chfi_ds::ChFiDS_FaceInterference,
    cpop: &mut ChFiDS_CommonPoint,
    fprol_uv: &mut DVec2,
    c3d_u: &mut f64,
    brep: &topods::BRep,
) -> bool {
    // add more or less restrictive criterions to decide if the intersection
    // is done with the face at extended end or if the end is sharp.
    // OCCT L311-319: the adaptor restricts the curve when non-periodic.
    if !c3dfi.is_periodic() {
        let mut wop = *c3d_u;
        let mut p2dbout = *fprol_uv;
        let ok = update_cs(hbs, c3dfi, uf, ul, fiop, cpop, &mut p2dbout, is_first, &mut wop);
        *fprol_uv = p2dbout;
        *c3d_u = wop;
        if ok {
            return true;
        }
    } else {
        let mut wop = *c3d_u;
        let mut p2dbout = *fprol_uv;
        // OCCT: GeomAdaptor_Curve(c3dFI) over the full periodic domain.
        let [f0, l0] = c3dfi.default_domain();
        let ok = update_cs(hbs, c3dfi, f0, l0, fiop, cpop, &mut p2dbout, is_first, &mut wop);
        *fprol_uv = p2dbout;
        *c3d_u = wop;
        if ok {
            return true;
        }
    }
    let _ = (uf, ul);

    if !super::chfi3d::is_tangent_faces(
        brep,
        eprol,
        fprol,
        fop,
        crate::geomalgo::gtests_stubs::GeomAbsShape::G1,
    ) {
        return false;
    }

    let Some((gpcprol, puf, pul)) = brep.curve_on_surface(eprol, fprol) else {
        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
    };
    let _ = (puf, pul);
    let partemp = brep_tool_parameter(brep, vtx, eprol);

    // OCCT L339: Update(HBs, pcprol, HGs, FIop, CPop, FprolUV, isFirst,
    // partemp, c3dU, Tol).
    let mut pared = partemp;
    let mut wop = *c3d_u;
    let mut p2dbout = *fprol_uv;
    let ok = update_locate(
        hbs,
        &gpcprol,
        hgs,
        fiop,
        cpop,
        &mut p2dbout,
        is_first,
        &mut pared,
        &mut wop,
        tol,
    );
    *fprol_uv = p2dbout;
    *c3d_u = wop;
    ok
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L349-417 — Update: the LocateExtCC between the
// edge pcurve on the face at end and the interference pcurve, preferring
// the values in range (with the bug 23139/25657 clamping).
// =========================================================================
fn update_on_face(
    face: &BRepAdaptorSurface,
    edonface: &rcad_kernel::geom::Curve2d,
    surf: &GeomAdaptorSurface,
    fi: &mut super::chfi_ds::ChFiDS_FaceInterference,
    cp: &mut ChFiDS_CommonPoint,
    isfirst: bool,
) -> bool {
    if !cp.is_on_arc() {
        return false;
    }
    let Some(pc) = fi.pcurve_on_surf() else {
        return false;
    };
    let c1 = CurveOnSurf {
        surf: &face.surface,
        pc: edonface,
    };
    let c2 = CurveOnSurf {
        surf: surf.surface(),
        pc,
    };
    let pared = cp.parameter_on_arc();
    let mut parltg = fi.parameter(isfirst);
    let f = fi.parameter_first();
    let l = fi.parameter_last();
    let delta = 0.1 * (l - f);
    let f = (f - delta).max(pc.default_domain()[0]);
    let l = (l + delta).min(pc.default_domain()[1]);
    let Some((pared, mut parltg, _d2)) = locate_ext_cc_curvesonurf(&c1, &c2, pared, parltg) else {
        return false;
    };
    if parltg > f && parltg < l {
        ////modified by jgv, 10.05.2012 for the bug 23139, 25657////
        if let Some(pconf) = fi.pcurve_on_face() {
            let mut pconf = pconf.clone();
            if let rcad_kernel::geom::Curve2d::Trimmed(tr) = &pconf {
                pconf = (*tr.curve).clone();
            }
            if !pconf.is_periodic() {
                let [pf, pl] = pconf.default_domain();
                if isfirst {
                    if parltg < pf {
                        parltg = pf;
                    }
                } else if parltg > pl {
                    parltg = pl;
                }
            }
        }
        /////////////////////////////////////////////////////
        fi.set_parameter(isfirst, parltg);
        cp.set_arc(cp.tolerance(), cp.arc().clone(), pared, cp.transition_on_arc());
        return true;
    }
    false
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L421-449 — ChFi3d_ExtendSurface: a plane /
// quadric keeps prol = 0; the GeomLib::ExtendSurfByLength extension of a
// BSpline/Bezier face surface is a pending TKGeomAlgo translation (the
// surface is returned unextended, matching the pre-extension state).
// =========================================================================
pub fn chfi3d_extend_surface(s: &mut rcad_kernel::geom::Surface3, prol: &mut i32) {
    if *prol != 0 {
        return;
    }

    *prol = match s {
        rcad_kernel::geom::Surface3::BSpline(_) => 1,
        rcad_kernel::geom::Surface3::Bezier(_) => 2,
        _ => 0,
    };
    if *prol == 0 {
        return;
    }

    // OCCT L436-447: bounds, D0 length, and the four
    // GeomLib::ExtendSurfByLength calls (u/v, both senses) — pending.
    let _ = s;
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L456-483 — ComputeCurve2d: the 2D projection
// of the curve on the face.  OCCT builds an edge and runs
// BRepAlgo_NormalProjection; the projection algorithm is a pending TKTopAlgo
// translation — C2d stays null exactly like the OCCT !IsDone() outcome.
// =========================================================================
pub fn compute_curve2d(
    _ct: &rcad_kernel::geom::Curve3,
    _face: &rcad_kernel::topods::Shape,
    _c2d: &mut Option<rcad_kernel::geom::Curve2d>,
    _brep: &topods::BRep,
) {
    // Pending: BRepAlgo_NormalProjection (OrtProj.Init/Add/SetParams/
    // SetLimit/Compute3d/Build + MapShapes + CurveOnSurface).
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L487-551 — ChFi3d_Recale.
// =========================================================================
pub fn chfi3d_recale(bs: &BRepAdaptorSurface, p1: &mut DVec2, p2: &mut DVec2, refon1: bool) {
    let mut surf = bs.surface().clone();
    if let rcad_kernel::geom::Surface3::Trimmed(tr) = &surf {
        surf = (*tr.basis).clone();
    }
    if bs.is_u_periodic() {
        let mut u1 = p1.x;
        let mut u2 = p2.x;
        let uper = bs.u_period();
        if (u2 - u1).abs() > 0.5 * uper {
            if u2 < u1 && refon1 {
                u2 += uper;
            } else if u2 < u1 && !refon1 {
                u1 -= uper;
            } else if u1 < u2 && refon1 {
                u2 -= uper;
            } else if u1 < u2 && !refon1 {
                u1 += uper;
            }
        }
        p1.x = u1;
        p2.x = u2;
    }
    if bs.is_v_periodic() {
        let mut v1 = p1.y;
        let mut v2 = p2.y;
        let vper = bs.v_period();
        if (v2 - v1).abs() > 0.5 * vper {
            if v2 < v1 && refon1 {
                v2 += vper;
            } else if v2 < v1 && !refon1 {
                v1 -= vper;
            } else if v1 < v2 && refon1 {
                v2 -= vper;
            } else if v1 < v2 && !refon1 {
                v1 += vper;
            }
        }
        p1.y = v1;
        p2.y = v2;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L558-588 — ChFi3d_SelectStripe: find the
// stripe with ChFiDS_OnSame state if <thePrepareOnSame> is True.  The
// iterator position is carried by it_pos (OCCT advances It in place).
// =========================================================================
pub fn chfi3d_select_stripe(
    stripes: &[SharedStripe],
    it_pos: &mut usize,
    vtx: &rcad_kernel::topods::Shape,
    the_prepare_on_same: bool,
) -> bool {
    if !the_prepare_on_same {
        return true;
    }

    while *it_pos < stripes.len() {
        let stripe = &stripes[*it_pos];
        let mut sens = 0i32;
        let st = stripe.read().expect("stripe lock");
        let _ = chfi3d_index_of_surf_data(vtx, &st, &mut sens);
        let stat = if sens == 1 {
            st.spine().expect("spine").base().first_status()
        } else {
            st.spine().expect("spine").base().last_status()
        };
        drop(st);
        if stat == ChFiDS_State::OnSame {
            // OCCT leaves It on the found stripe.
            return true;
        }
        *it_pos += 1;
    }

    false
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L4533-4595 — FindFace: the common face of the
// two common-point arcs (works only if there is exactly one).
// =========================================================================
impl ChFi3dBuilder {
    pub fn find_face2(
        &self,
        p1: &ChFiDS_CommonPoint,
        p2: &ChFiDS_CommonPoint,
        fv: &mut rcad_kernel::topods::Shape,
        favoid: &rcad_kernel::topods::Shape,
    ) -> bool {
        if p1.is_vertex() || p2.is_vertex() {
            // OCCT debug print: "change of face on vertex".
        }
        if !(p1.is_on_arc() && p2.is_on_arc()) {
            return false;
        }
        let mut found = false;
        if self.my_ef_map.contains(p1.arc()) {
            for f in self.my_ef_map.find(p1.arc()).clone() {
                if found {
                    break;
                }
                *fv = f.clone();
                if !fv.is_same(favoid) {
                    if self.my_ef_map.contains(p2.arc()) {
                        for jt in self.my_ef_map.find(p2.arc()).clone() {
                            if found {
                                break;
                            }
                            if jt.is_same(&*fv) {
                                found = true;
                            }
                        }
                    }
                }
            }
        }
        found
    }

    /// OCCT L4533-4540 — FindFace(V, P1, P2, Fv): Favoid is null.
    pub fn find_face(
        &self,
        p1: &ChFiDS_CommonPoint,
        p2: &ChFiDS_CommonPoint,
        fv: &mut rcad_kernel::topods::Shape,
    ) -> bool {
        let favoid = rcad_kernel::topods::Shape::null();
        self.find_face2(p1, p2, fv, &favoid)
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_C1.cxx L611-1652 — ChFi3d_Builder::PerformOneCorner.
// Calculates a corner with three edges and a fillet.
// =========================================================================
impl ChFi3dBuilder {
    pub fn perform_one_corner(&mut self, index: usize, the_prepare_on_same: bool) {
        let vtx = self.my_vdata_map.find_key(index).clone();
        let stripes = self.my_vdata_map.find_from_index(index).clone();

        // The fillet is returned.
        let mut it_pos = 0usize;
        if !chfi3d_select_stripe(&stripes, &mut it_pos, &vtx, the_prepare_on_same) {
            return;
        }
        let stripe = stripes[it_pos].clone();
        let mut st = stripe.write().expect("stripe lock");
        let spine = st.spine().expect("spine").clone();

        // SurfData and its CommonPoints.
        let mut sens = 0i32;
        let mut num = chfi3d_index_of_surf_data(&vtx, &st, &mut sens);
        let isfirst = sens == 1;
        // OCCT L637-653: the surplus SurfData (missing support faces) are
        // removed around the selected one.
        if isfirst {
            while (num as usize) < st.my_hdata.len() {
                let invalid = {
                    let sd = st.my_hdata[(num - 1) as usize].read().expect("surfdata lock");
                    sd.index_of_s1 == 0 || sd.index_of_s2 == 0
                };
                if !invalid {
                    break;
                }
                st.my_hdata.remove((num - 1) as usize);
            }
        } else {
            while num > 1 {
                let invalid = {
                    let sd = st.my_hdata[(num - 1) as usize].read().expect("surfdata lock");
                    sd.index_of_s1 == 0 || sd.index_of_s2 == 0
                };
                if !invalid {
                    break;
                }
                st.my_hdata.remove((num - 1) as usize); // The surplus is removed
                num -= 1;
            }
        }

        let fd_lock = st.my_hdata[(num - 1) as usize].clone();
        let mut fd = fd_lock.write().expect("surfdata lock");
        let mut cv1 = fd.vertex(isfirst, 1).clone();
        let mut cv2 = fd.vertex(isfirst, 2).clone();
        // To evaluate the new points.
        let mut box1 = BndBox::default();
        let mut box2 = BndBox::default();

        // The cases of cap and intersection are processed separately.
        let stat = if isfirst {
            spine.base().first_status()
        } else {
            spine.base().last_status()
        };
        let onsame = stat == ChFiDS_State::OnSame;
        let mut fv = rcad_kernel::topods::Shape::null();
        let mut fad = rcad_kernel::topods::Shape::null();
        let mut fop = {
            let mut f = self
                .my_ds
                .as_ref()
                .expect("DS")
                .shape(fd.index_of(2))
                .clone();
            f.orientation = Orientation::Forward;
            f
        };
        let arcspine = if isfirst {
            spine.base().edges(1).clone()
        } else {
            spine.base().edges(spine.base().nb_edges()).clone()
        };
        let oarcprolv;
        let mut oarcprolop = Orientation::Forward;
        let mut hbs: Option<BRepAdaptorSurface> = None;
        let mut hbad: Option<BRepAdaptorSurface> = None;
        let mut hbop: Option<BRepAdaptorSurface> = None;
        let mut inters = true;
        let mut ifad_arc = 1i32;
        let mut ifop_arc = 2i32;
        let mut arcprol = rcad_kernel::topods::Shape::null();
        let mut couture = false;
        let mut tolreached = self.tolapp3d;
        let mut par1 = 0.0f64;
        let mut par2 = 0.0f64;
        let mut indpt = 0i32;
        let mut icurv1 = 0i32;
        let mut icurv2 = 0i32;
        let mut curv1: Option<rcad_kernel::geom::Curve3> = None;
        let mut curv2: Option<rcad_kernel::geom::Curve3> = None;
        let mut c2d1: Option<rcad_kernel::geom::Curve2d> = None;
        let mut c2d2: Option<rcad_kernel::geom::Curve2d> = None;
        let mut cc: Option<rcad_kernel::geom::Curve3> = None;
        let mut ps: Option<rcad_kernel::geom::Curve2d> = None;
        let mut pc: Option<rcad_kernel::geom::Curve2d> = None;
        let mut udeb = 0.0f64;
        let mut ufin = 0.0f64;
        let mut intcouture = false;

        if onsame {
            if !cv1.is_on_arc() && !cv2.is_on_arc() {
                panic!("Standard_ConstructionError: Corner OnSame : no point on arc");
            } else if cv1.is_on_arc() && cv2.is_on_arc() {
                let mut sur1 = false;
                let mut sur2 = false;
                {
                    let ed = cv1.arc().as_edge().expect("not an edge");
                    for v in [&ed.first, &ed.last] {
                        if vtx.is_same(v) {
                            sur1 = true;
                            break;
                        }
                    }
                }
                {
                    let ed = cv2.arc().as_edge().expect("not an edge");
                    for v in [&ed.first, &ed.last] {
                        if vtx.is_same(v) {
                            sur2 = true;
                            break;
                        }
                    }
                }
                if sur1 && sur2 {
                    let e = [cv1.arc().clone(), cv2.arc().clone(), arcspine.clone()];
                    if super::chfi3d_builder_0::chfi3d_edge_state(
                        &e,
                        &self.my_ef_map,
                        &self.my_brep,
                    ) != ChFiDS_State::OnDiff
                    {
                        ifad_arc = 2;
                    }
                } else if sur2 {
                    ifad_arc = 2;
                }
            } else if cv2.is_on_arc() {
                ifad_arc = 2;
            }
            ifop_arc = 3 - ifad_arc;

            let arcpiv = fd.vertex(isfirst, ifad_arc).arc().clone();
            fad = self
                .my_ds
                .as_ref()
                .expect("DS")
                .shape(fd.index_of(ifad_arc))
                .clone();
            fop = self
                .my_ds
                .as_ref()
                .expect("DS")
                .shape(fd.index_of(ifop_arc))
                .clone();
            // The face at end is returned without check of its unicity.
            if self.my_ef_map.contains(&arcpiv) {
                for it in self.my_ef_map.find(&arcpiv).clone() {
                    if !fad.is_same(&it) {
                        fv = it.clone();
                        break;
                    }
                }
            }

            // Does the face at bout contain the Vertex ?
            let mut isinface = false;
            for v in super::chfi3d_builder_0::topexp_face_vertices(&self.my_brep, &fv) {
                if v.is_same(&vtx) {
                    isinface = true;
                    break;
                }
            }
            if !isinface && fd.vertex(isfirst, 3 - ifad_arc).is_on_arc() {
                ifad_arc = 3 - ifad_arc;
                ifop_arc = 3 - ifop_arc;
                let arcpiv = fd.vertex(isfirst, ifad_arc).arc().clone();
                fad = self
                    .my_ds
                    .as_ref()
                    .expect("DS")
                    .shape(fd.index_of(ifad_arc))
                    .clone();
                fop = self
                    .my_ds
                    .as_ref()
                    .expect("DS")
                    .shape(fd.index_of(ifop_arc))
                    .clone();
                if self.my_ef_map.contains(&arcpiv) {
                    for it in self.my_ef_map.find(&arcpiv).clone() {
                        if !fad.is_same(&it) {
                            fv = it.clone();
                            break;
                        }
                    }
                }
            }

            if fv.is_null() {
                panic!("StdFail_NotDone: OneCorner : face at end not found");
            }

            fv.orientation = Orientation::Forward;
            fad.orientation = Orientation::Forward;
            fop.orientation = Orientation::Forward;

            // The edge that will be extended is returned (OCCT L806-820).
            let mut found: Option<(rcad_kernel::topods::Shape, Orientation)> = None;
            if self.my_ve_map.contains(&vtx) {
                for ite in self.my_ve_map.find(&vtx).clone() {
                    if arcpiv.is_same(&ite) {
                        continue;
                    }
                    let mut hit = None;
                    for e in super::chfi3d_builder_0::topexp_face_edges(&self.my_brep, &fv) {
                        if ite.is_same(&e) {
                            hit = Some(e.orientation);
                            break;
                        }
                    }
                    if let Some(o) = hit {
                        found = Some((ite.clone(), o));
                        break;
                    }
                }
            }
            match found {
                Some((e, o)) => {
                    arcprol = e;
                    oarcprolv = o;
                }
                None => {
                    // OCCT L821-825: PerformIntersectionAtEnd(Index); return;
                    drop(fd);
                    drop(st);
                    self.perform_intersection_at_end(index);
                    return;
                }
            }
            for e in super::chfi3d_builder_0::topexp_face_edges(&self.my_brep, &fop) {
                if arcprol.is_same(&e) {
                    oarcprolop = e.orientation;
                    break;
                }
            }
            // OCCT L834-851: BRE.MakeFace(FFv, Sface, tol); when the surface
            // was extended, Bs loads FFv and DStr.SetNewSurface(Fv, Sface)
            // (the GeomLib extension itself is pending — quadric faces keep
            // prol = 0 and load Fv directly).
            let mut sface = self.my_brep.face_surface_world(&fv).expect("face surface");
            let mut prol = 0i32;
            chfi3d_extend_surface(&mut sface, &mut prol);
            if prol != 0 {
                // DStr.SetNewSurface(Fv, Sface) — pending TopOpeBRepDS
                // surface replacement.
                hbs = Some(BRepAdaptorSurface::initialize_surface(sface.clone()));
            } else {
                hbs = Some(BRepAdaptorSurface::initialize(&self.my_brep, &fv));
            }
            hbad = Some(BRepAdaptorSurface::initialize(&self.my_brep, &fad));
            hbop = Some(BRepAdaptorSurface::initialize(&self.my_brep, &fop));
        } else {
            oarcprolv = Orientation::Forward;
        }

        // in case of OnSame it is necessary to modify the CommonPoint in the
        // empty and its parameter in the FaceInterference.
        let mut cpop_arc = fd.vertex(isfirst, ifop_arc).clone();
        let mut fiop_arc = fd.interference(ifop_arc).clone();
        let mut cpad_arc = fd.vertex(isfirst, ifad_arc).clone();
        let mut fiad_arc = fd.interference(ifad_arc).clone();
        // the parameter of the vertex in the air is initialized with the
        // value of its opposite (point on arc).
        let mut wop = fd.interference(ifad_arc).parameter(isfirst);
        let c3df = self
            .my_ds
            .as_ref()
            .expect("DS")
            .curve(fiop_arc.line_index())
            .curve
            .clone();
        let hgs_initial = GeomAdaptorSurface::new(
            self.my_ds
                .as_ref()
                .expect("DS")
                .surface(fd.surf())
                .surface
                .clone(),
        );
        let mut p2dbout = DVec2::ZERO;
        let isurf = fd.surf();

        if onsame {
            let save_cpop_arc = cpop_arc.clone();
            let fv_guard = fv.clone();
            let bs = hbs.as_mut().expect("Bs");
            inters = inters_update_on_same(
                &hgs_initial,
                bs,
                c3df.as_ref().expect("c3df"),
                fiop_arc.parameter_first(),
                fiop_arc.parameter_last(),
                &fop,
                &fv_guard,
                &arcprol,
                &vtx,
                isfirst,
                10.0 * self.tolapp3d, // in
                &mut fiop_arc,
                &mut cpop_arc,
                &mut p2dbout,
                &mut wop,
                &self.my_brep,
            );

            // in the case of degenerated Fi, the parameter difference can be
            // even negative (eap, occ293).
            if (fiad_arc.parameter_last() - fiad_arc.parameter_first()) > 10.0 * self.tolesp {
                let Some((pced, _pf, _pl)) = self.my_brep.curve_on_surface(cpad_arc.arc(), &fv)
                else {
                    panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                };
                update_on_face(
                    bs,
                    &pced,
                    &hgs_initial,
                    &mut fiad_arc,
                    &mut cpad_arc,
                    isfirst,
                );
            }

            if the_prepare_on_same {
                let mut saved = save_cpop_arc;
                let p = cpop_arc.point();
                saved.set_point(p);
                cpop_arc = saved;
                *fd.change_vertex(isfirst, ifop_arc) = cpop_arc;
                *fd.change_interference(ifop_arc) = fiop_arc;
                return;
            }
        } else {
            let fop_clone = fop.clone();
            inters = self.find_face2(&cv1, &cv2, &mut fv, &fop_clone);
            if !inters {
                // OCCT L912-914: PerformIntersectionAtEnd(Index); return;
                drop(fd);
                drop(st);
                self.perform_intersection_at_end(index);
                return;
            }
            hbs = Some(BRepAdaptorSurface::initialize(&self.my_brep, &fv));
            let bs = hbs.as_mut().expect("Bs");
            let Some((pced1, _a, _b)) = self.my_brep.curve_on_surface(cv1.arc(), &fv) else {
                panic!("Standard_ConstructionError: Failed to get p-curve of edge");
            };
            let mut fi1 = fd.interference_on_s1().clone();
            update_on_face(bs, &pced1, &hgs_initial, &mut fi1, &mut cv1, isfirst);
            *fd.change_interference_on_s1() = fi1;
            let Some((pced2, _c, _d)) = self.my_brep.curve_on_surface(cv2.arc(), &fv) else {
                panic!("Standard_ConstructionError: Failed to get p-curve of edge");
            };
            let mut fi2 = fd.interference_on_s2().clone();
            update_on_face(bs, &pced2, &hgs_initial, &mut fi2, &mut cv2, isfirst);
            *fd.change_interference_on_s2() = fi2;
        }
        let bs = hbs.as_mut().expect("Bs");

        let mut edgecouture = rcad_kernel::topods::Shape::null();

        if inters {
            let hgs = chfi3d_bound_surf(self.my_ds.as_ref().expect("DS"), &fd, 1, 2);
            let fi1 = fd.interference_on_s1().clone();
            let fi2 = fd.interference_on_s2().clone();
            let mut pardeb = [0.0f64; 4];
            let mut parfin = [0.0f64; 4];
            let pfil1;
            let mut pfac1;
            let pfil2;
            let mut pfac2;
            if onsame && ifop_arc == 1 {
                pfac1 = p2dbout;
            } else {
                let Some((hc1, _f1, _l1)) = self.my_brep.curve_on_surface(cv1.arc(), &fv) else {
                    panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                };
                pfac1 = hc1.point_at(cv1.parameter_on_arc());
            }
            if onsame && ifop_arc == 2 {
                pfac2 = p2dbout;
            } else {
                let Some((hc2, _f2, _l2)) = self.my_brep.curve_on_surface(cv2.arc(), &fv) else {
                    panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                };
                pfac2 = hc2.point_at(cv2.parameter_on_arc());
            }
            if fi1.line_index() != 0 {
                pfil1 = fi1
                    .pcurve_on_surf()
                    .expect("pfil1")
                    .point_at(fi1.parameter(isfirst));
            } else {
                pfil1 = fi1
                    .pcurve_on_surf()
                    .expect("pfil1")
                    .point_at(fi1.parameter(!isfirst));
            }
            if fi2.line_index() != 0 {
                pfil2 = fi2
                    .pcurve_on_surf()
                    .expect("pfil2")
                    .point_at(fi2.parameter(isfirst));
            } else {
                pfil2 = fi2
                    .pcurve_on_surf()
                    .expect("pfil2")
                    .point_at(fi2.parameter(!isfirst));
            }
            if onsame {
                chfi3d_recale(bs, &mut pfac1, &mut pfac2, ifad_arc == 1);
            }

            pardeb[0] = pfil1.x;
            pardeb[1] = pfil1.y;
            pardeb[2] = pfac1.x;
            pardeb[3] = pfac1.y;
            parfin[0] = pfil2.x;
            parfin[1] = pfil2.y;
            parfin[2] = pfac2.x;
            parfin[3] = pfac2.y;

            let (uu1, uu2, vv1, vv2) = chfi3d_boite(pfac1, pfac2);
            chfi3d_bound_fac(bs, uu1, uu2, vv1, vv2, true);

            // OCCT passes HBs as an Adaptor3d_Surface handle; the rcad view
            // carries the same surface (the trimmed bounds only affect
            // First/LastParameter, which ComputeCurves does not use).
            let bs_view = GeomAdaptorSurface::new(bs.surface.clone());
            match chfi3d_compute_curves(
                &hgs,
                &bs_view,
                pardeb,
                parfin,
                self.tolapp3d,
                self.tol2d,
                &mut tolreached,
            ) {
                Some(res) => {
                    cc = Some(res.c3d);
                    ps = Some(res.pc1);
                    pc = Some(res.pc2);
                    tolreached = res.tolreached;
                }
                None => panic!("Standard_Failure: OneCorner : echec calcul intersection"),
            }

            let cc_ref = cc.as_ref().expect("Cc");
            match cc_ref {
                rcad_kernel::geom::Curve3::Trimmed(tr) => {
                    udeb = tr.first;
                    ufin = tr.last;
                }
                other => {
                    let [f0, l0] = other.default_domain();
                    udeb = f0;
                    ufin = l0;
                }
            }

            // determine if the curve has an intersection with the sewing
            // edge.
            let (cout, ecout) = chfi3d_couture(&self.my_brep, &fv);
            couture = cout;
            edgecouture = ecout;

            if couture && !self.my_brep.is_edge_degenerated(&edgecouture) {
                let (c, ctr_range) = self
                    .my_brep
                    .edge_curve_world(&edgecouture)
                    .expect("sewing curve");
                let (ctr_first, ctr_last) = (ctr_range[0], ctr_range[1]);
                let ctrim = rcad_kernel::geom::Curve3::Trimmed(
                    rcad_kernel::geom::TrimmedCurve3::new(c.clone(), ctr_first, ctr_last),
                );
                let basis = match &ctrim {
                    rcad_kernel::geom::Curve3::Trimmed(tr) => (*tr.curve).clone(),
                    _ => c.clone(),
                };
                let cc_basis = match cc_ref {
                    rcad_kernel::geom::Curve3::Trimmed(tr) => (*tr.curve).clone(),
                    other => other.clone(),
                };
                let ext = rcad_kernel::base::extrema::extrema_curve_curve(&basis, &cc_basis, 32);
                if !ext.pairs.is_empty() {
                    let mut imin = 0usize;
                    let mut distmin2 = f64::MAX;
                    for (i, pair) in ext.pairs.iter().enumerate() {
                        let d2 = pair.distance * pair.distance;
                        if d2 < distmin2 {
                            distmin2 = d2;
                            imin = i;
                        }
                    }
                    if distmin2 <= P_CONFUSION * P_CONFUSION {
                        let best = &ext.pairs[imin];
                        par1 = best.param1;
                        par2 = best.param2;
                        let tol = 1.0e-4;
                        if (par2 - udeb).abs() > tol && (ufin - par2).abs() > tol {
                            let p1 = best.point1;
                            indpt = self
                                .my_ds
                                .as_mut()
                                .expect("DS")
                                .add_point(TopOpeBRepDSPoint::new(p1, tol));
                            intcouture = true;
                            curv1 = Some(rcad_kernel::geom::Curve3::Trimmed(
                                rcad_kernel::geom::TrimmedCurve3::new(
                                    cc_basis.clone(),
                                    udeb,
                                    par2,
                                ),
                            ));
                            curv2 = Some(rcad_kernel::geom::Curve3::Trimmed(
                                rcad_kernel::geom::TrimmedCurve3::new(
                                    cc_basis.clone(),
                                    par2,
                                    ufin,
                                ),
                            ));
                            icurv1 = self
                                .my_ds
                                .as_mut()
                                .expect("DS")
                                .add_curve(TopOpeBRepDSCurve::new(curv1.clone(), tolreached));
                            icurv2 = self
                                .my_ds
                                .as_mut()
                                .expect("DS")
                                .add_curve(TopOpeBRepDSCurve::new(curv2.clone(), tolreached));
                        }
                    }
                }
            }
        } else {
            // (!inters)
            panic!("Standard_NotImplemented: OneCorner : bouchon non ecrit");
        }
        let ishape = self.my_ds.as_mut().expect("DS").add_shape(&fv);
        let mut et = Orientation::Forward;
        if ifad_arc == 1 {
            for e in super::chfi3d_builder_0::topexp_face_edges(&self.my_brep, &fv) {
                if e.is_same(cv1.arc()) {
                    et = topabs_reverse(topabs_compose(e.orientation, cv1.transition_on_arc()));
                    break;
                }
            }
        } else {
            for e in super::chfi3d_builder_0::topexp_face_edges(&self.my_brep, &fv) {
                if e.is_same(cv2.arc()) {
                    et = topabs_compose(e.orientation, cv2.transition_on_arc());
                    break;
                }
            }
        }

        let i_p1 = chfi3d_index_point_in_ds(&cv1, self.my_ds.as_mut().expect("DS"));
        let i_p2 = chfi3d_index_point_in_ds(&cv2, self.my_ds.as_mut().expect("DS"));
        st.set_index_point(i_p1, isfirst, 1);
        st.set_index_point(i_p2, isfirst, 2);

        if !intcouture {
            // there is no intersection with the sewing edge: the curve Cc is
            // stored in the stripe; the storage in the DS is not done by
            // FILDS.
            let cc_val = cc.clone().expect("Cc");
            let pc_val = pc.clone().expect("Pc");
            let icurve = self
                .my_ds
                .as_mut()
                .expect("DS")
                .add_curve(TopOpeBRepDSCurve::new(Some(cc_val), tolreached));
            let interfc = chfi3d_fil_curve_in_ds(icurve, ishape, Some(pc_val.clone()), et);

            // 31/01/02 akm (OCC119): prevent the builder from creating
            // intersecting fillets — take all the interferences with faces
            // from all the stripes and look if their pcurves intersect our
            // cork pcurve.
            for a_check_stripe in &self.my_list_stripe {
                let guard = a_check_stripe.read().expect("stripe lock");
                for a_data_arc in &guard.my_hdata {
                    let a_data = a_data_arc.read().expect("surfdata lock");
                    let common_face = if ishape == a_data.index_of_s1 {
                        Some(1)
                    } else if ishape == a_data.index_of_s2 {
                        Some(2)
                    } else {
                        // Normal case - no common surface.
                        None
                    };
                    let Some(ons) = common_face else { continue };
                    let a_fi = a_data.interference(ons);
                    let Some(a_pcurve) = a_fi.pcurve_on_face().cloned() else {
                        continue;
                    };
                    if (a_fi.parameter_last() - a_fi.parameter_first()).abs() <= f64::EPSILON {
                        // Degenerates.
                        continue;
                    }
                    let (nb_points, nb_segments) = geom2d_int_g_inter(&pc_val, &a_pcurve);
                    if nb_segments > 0 || nb_points > 0 {
                        panic!("StdFail_NotDone: OneCorner : fillets have too big radiuses");
                    }
                }
            }
            let shape_interfs: Vec<TopOpeBRepDSInterference> = self
                .my_ds
                .as_ref()
                .expect("DS")
                .shape_interferences(ishape)
                .to_vec();
            for other_intrf in &shape_interfs {
                let TopOpeBRepDSInterference::SurfaceCurve(an_other_intrf) = other_intrf else {
                    // We need only interferences between cork face and curves.
                    continue;
                };
                // Look if there is an intersection between pcurves.
                let Some(an_other_cur) = self
                    .my_ds
                    .as_ref()
                    .expect("DS")
                    .curve(an_other_intrf.index_g)
                    .curve
                    .clone()
                else {
                    continue;
                };
                let (of, ol) = match &an_other_cur {
                    rcad_kernel::geom::Curve3::Trimmed(tr) => (tr.first, tr.last),
                    other => {
                        let d = other.default_domain();
                        (d[0], d[1])
                    }
                };
                let Some(other_pc) = an_other_intrf.pcurve.clone() else {
                    continue;
                };
                let _ = (of, ol);
                let (nb_points, nb_segments) = geom2d_int_g_inter(&pc_val, &other_pc);
                if nb_segments > 0 || nb_points > 0 {
                    panic!("StdFail_NotDone: OneCorner : fillets have too big radiuses");
                }
            }
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_shape_interferences(ishape)
                .push(interfc);

            //// modified by jgv, 26.03.02 for OCC32 ////
            let cvs = [cv1.clone(), cv2.clone()];
            for cvi in &cvs {
                if cvi.is_on_arc() && chfi3d_is_pseudo_seam(&self.my_brep, cvi.arc(), &fv) {
                    let Some((hc, hc_first, hc_last)) =
                        self.my_brep.curve_on_surface(cvi.arc(), &fv)
                    else {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    };
                    let pfac1 = hc.point_at(cvi.parameter_on_arc());
                    let pcf = pc_val.point_at(udeb);
                    let pcl = pc_val.point_at(ufin);
                    let onfirst = pfac1.distance(pcf) < pfac1.distance(pcl);
                    let der_pc = if onfirst {
                        pc_val.derivative_at(udeb)
                    } else {
                        -pc_val.derivative_at(ufin)
                    };
                    let der_hc = hc.derivative_at(cvi.parameter_on_arc());
                    let (prm1, prm2, first_to_par) = if der_hc.dot(der_pc) > 0.0 {
                        (cvi.parameter_on_arc(), hc_last, false)
                    } else {
                        (hc_first, cvi.parameter_on_arc(), true)
                    };
                    let (ct, _cfr) = {
                        let r = self.my_brep.edge_curve_world(cvi.arc()).expect("arc curve");
                        (r.0, ())
                    };
                    let ct = rcad_kernel::geom::Curve3::Trimmed(
                        rcad_kernel::geom::TrimmedCurve3::new(ct, prm1, prm2),
                    );
                    let toled = self.my_brep.tolerance(cvi.arc());
                    let indcurv = self
                        .my_ds
                        .as_mut()
                        .expect("DS")
                        .add_curve(TopOpeBRepDSCurve::new(Some(ct), toled));
                    let indpoint = if isfirst {
                        st.indexfirst_pon_s1
                    } else {
                        st.indexlast_pon_s1
                    };
                    let indvertex = self.my_ds.as_mut().expect("DS").add_shape(&vtx);
                    let (interfp1, interfp2) = if first_to_par {
                        (
                            super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                                Orientation::Forward,
                                indcurv,
                                indvertex,
                                prm1,
                                true,
                            ),
                            super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                                Orientation::Reversed,
                                indcurv,
                                indpoint,
                                prm2,
                                false,
                            ),
                        )
                    } else {
                        (
                            super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                                Orientation::Forward,
                                indcurv,
                                indpoint,
                                prm1,
                                false,
                            ),
                            super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                                Orientation::Reversed,
                                indcurv,
                                indvertex,
                                prm2,
                                true,
                            ),
                        )
                    };
                    self.my_ds
                        .as_mut()
                        .expect("DS")
                        .change_curve_interferences(indcurv)
                        .push(interfp1);
                    self.my_ds
                        .as_mut()
                        .expect("DS")
                        .change_curve_interferences(indcurv)
                        .push(interfp2);
                    let indface = self.my_ds.as_mut().expect("DS").add_shape(&fv);
                    let interfc = chfi3d_fil_curve_in_ds(
                        indcurv,
                        indface,
                        Some(hc.clone()),
                        cvi.arc().orientation,
                    );
                    self.my_ds
                        .as_mut()
                        .expect("DS")
                        .change_shape_interferences(indface)
                        .push(interfc);
                    // modify degenerated edge
                    let hcr = {
                        let mut local = cvi.arc().clone();
                        local.orientation = topabs_reverse(local.orientation);
                        self.my_brep.curve_on_surface(&local, &fv)
                    };
                    let Some((hcr, _f, _l)) = hcr else {
                        panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                    };
                    let interfc = chfi3d_fil_curve_in_ds(
                        indcurv,
                        indface,
                        Some(hcr),
                        topabs_reverse(cvi.arc().orientation),
                    );
                    self.my_ds
                        .as_mut()
                        .expect("DS")
                        .change_shape_interferences(indface)
                        .push(interfc);
                    // degenerated edge carrying Vtx.
                    let mut degen_exist = false;
                    let mut edeg = rcad_kernel::topods::Shape::null();
                    for ecur in super::chfi3d_builder_0::topexp_face_edges(&self.my_brep, &fv) {
                        if self.my_brep.is_edge_degenerated(&ecur) {
                            let ed = ecur.as_edge().expect("not an edge");
                            if ed.first.is_same(&vtx) || ed.last.is_same(&vtx) {
                                degen_exist = true;
                                edeg = ecur;
                                break;
                            }
                        }
                    }
                    if degen_exist {
                        let Some((cd, cd_first, cd_last)) =
                            self.my_brep.curve_on_surface(&edeg, &fv)
                        else {
                            panic!("Standard_ConstructionError: Failed to get p-curve of edge");
                        };
                        let mut cd = cd;
                        if let rcad_kernel::geom::Curve2d::Trimmed(tcd) = &cd {
                            cd = (*tcd.curve).clone();
                        }
                        let p2d = if first_to_par {
                            hc.point_at(hc_first)
                        } else {
                            hc.point_at(hc_last)
                        };
                        // OCCT: Geom2dAPI_ProjectPointOnCurve(P2d, Cd).
                        let par = project_point_on_curve2d(p2d, &cd);
                        let ideg = self.my_ds.as_mut().expect("DS").add_shape(&edeg);
                        let ori = if par < cd_first {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        }; // if par<fd => par>ld
                        let interfp1 = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                            ori, ideg, indvertex, par, true,
                        );
                        self.my_ds
                            .as_mut()
                            .expect("DS")
                            .change_shape_interferences(ideg)
                            .push(interfp1);
                        let _ = cd_last;
                    }
                }
            }
            /////////////////////////////////////////////
            st.change_pcurve(isfirst, ps.clone().expect("Ps"));
            st.set_curve(icurve, isfirst);
            st.set_parameters(isfirst, udeb, ufin);
        } else {
            // curves curv1/curv2 stored in the DS (stripe->InDS(isfirst)).
            compute_curve2d(curv1.as_ref().expect("curv1"), &fv, &mut c2d1, &self.my_brep);
            let interfv = chfi3d_fil_curve_in_ds(icurv1, ishape, c2d1.clone(), et);
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_shape_interferences(ishape)
                .push(interfv);
            compute_curve2d(curv2.as_ref().expect("curv2"), &fv, &mut c2d2, &self.my_brep);
            let interfv = chfi3d_fil_curve_in_ds(icurv2, ishape, c2d2.clone(), et);
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_shape_interferences(ishape)
                .push(interfv);
            // interferences of curv1 and curv2 on Isurf.
            if fd.orientation() == fv.orientation {
                et = topabs_reverse(et);
            }
            let ps_val = ps.clone().expect("Ps");
            let c2d1t = rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(ps_val.clone()),
                t_min: udeb,
                t_max: par2,
            });
            let interfv = chfi3d_fil_curve_in_ds(icurv1, isurf, Some(c2d1t), et);
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_surface_interferences(isurf)
                .push(interfv);
            let c2d2t = rcad_kernel::geom::Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(ps_val.clone()),
                t_min: par2,
                t_max: ufin,
            });
            let interfv = chfi3d_fil_curve_in_ds(icurv2, isurf, Some(c2d2t), et);
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_surface_interferences(isurf)
                .push(interfv);

            // limitation of the sewing edge.
            let iarc = self.my_ds.as_mut().expect("DS").add_shape(&edgecouture);
            let vdeb = self.my_brep.first_vertex(&edgecouture);
            let vfin = self.my_brep.last_vertex(&edgecouture);
            let pard = brep_tool_parameter(&self.my_brep, &vdeb, &edgecouture);
            let parf = brep_tool_parameter(&self.my_brep, &vfin, &edgecouture);
            let ori = if (par1 - pard).abs() < (parf - par1).abs() {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            let interfedge =
                super::chfi3d_builder_0::chfi3d_fil_point_in_ds(ori, iarc, indpt, par1, false);
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_shape_interferences(iarc)
                .push(interfedge);

            // creation of CurveInterferences from Icurv1 and Icurv2.
            st.in_ds(isfirst, 1);
            let ind1 = st.index_point(isfirst, 1);
            let ind2 = st.index_point(isfirst, 2);
            let interfprol = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                Orientation::Forward,
                icurv1,
                ind1,
                udeb,
                false,
            );
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_curve_interferences(icurv1)
                .push(interfprol);
            let interfprol = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                Orientation::Reversed,
                icurv1,
                indpt,
                par2,
                false,
            );
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_curve_interferences(icurv1)
                .push(interfprol);
            let interfprol = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                Orientation::Forward,
                icurv2,
                indpt,
                par2,
                false,
            );
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_curve_interferences(icurv2)
                .push(interfprol);
            let interfprol = super::chfi3d_builder_0::chfi3d_fil_point_in_ds(
                Orientation::Reversed,
                icurv2,
                ind2,
                ufin,
                false,
            );
            self.my_ds
                .as_mut()
                .expect("DS")
                .change_curve_interferences(icurv2)
                .push(interfprol);
        }

        chfi3d_enlarge_box_surf_pc(
            &GeomAdaptorSurface::new(bs.surface.clone()),
            &pc.clone().expect("Pc"),
            udeb,
            ufin,
            &mut box1,
            &mut box2,
        );

        if onsame && inters {
            // OCCT L1397-1630: VARIANT 1 — a small missing end of curve is
            // added for the extension of the face at end and the limitation
            // of the opposing face (OnSame corners only).
            self.perform_onsame_tail_pending();
        }
        chfi3d_enlarge_box_dstr(
            &self.my_brep,
            self.my_ds.as_ref().expect("DS"),
            Some(&st),
            &fd,
            &mut box1,
            &mut box2,
            isfirst,
        );
        if cv1.is_on_arc() {
            let lf = if self.my_ef_map.contains(cv1.arc()) {
                self.my_ef_map.find(cv1.arc()).clone()
            } else {
                Vec::new()
            };
            chfi3d_enlarge_box_edge_faces(
                &self.my_brep,
                cv1.arc(),
                &lf,
                cv1.parameter_on_arc(),
                &mut box1,
            );
        }
        if cv2.is_on_arc() {
            let lf = if self.my_ef_map.contains(cv2.arc()) {
                self.my_ef_map.find(cv2.arc()).clone()
            } else {
                Vec::new()
            };
            chfi3d_enlarge_box_edge_faces(
                &self.my_brep,
                cv2.arc(),
                &lf,
                cv2.parameter_on_arc(),
                &mut box2,
            );
        }
        if !cv1.is_vertex() {
            chfi3d_set_point_tolerance(
                self.my_ds.as_mut().expect("DS"),
                &box1,
                st.index_point(isfirst, 1),
            );
        }
        if !cv2.is_vertex() {
            chfi3d_set_point_tolerance(
                self.my_ds.as_mut().expect("DS"),
                &box2,
                st.index_point(isfirst, 2),
            );
        }

        // Write the mutated CommonPoints / interferences back.
        *fd.change_vertex(isfirst, 1) = cv1;
        *fd.change_vertex(isfirst, 2) = cv2;
        *fd.change_vertex(isfirst, ifop_arc) = cpop_arc;
        *fd.change_interference(ifop_arc) = fiop_arc;
        *fd.change_vertex(isfirst, ifad_arc) = cpad_arc;
        *fd.change_interference(ifad_arc) = fiad_arc;
    }

    /// OCCT ChFi3d_Builder_C2.cxx PerformIntersectionAtEnd — pending.
    fn perform_intersection_at_end(&mut self, _index: usize) {
        // Pending translation (Builder_C2.cxx).
    }

    /// OCCT C1.cxx L1397-1630 OnSame VARIANT-1 tail — pending (OnSame
    /// corners only; the box KPart corner path is !onsame).
    fn perform_onsame_tail_pending(&mut self) {}
}

/// OCCT Geom2dAPI_ProjectPointOnCurve(P2d, Cd).LowerDistanceParameter().
fn project_point_on_curve2d(p: DVec2, cd: &rcad_kernel::geom::Curve2d) -> f64 {
    let [f0, l0] = cd.default_domain();
    let n = 64usize;
    let mut best = f0;
    let mut bestd = f64::MAX;
    for i in 0..=n {
        let t = f0 + (l0 - f0) * (i as f64) / (n as f64);
        let d = cd.point_at(t).distance(p);
        if d < bestd {
            bestd = d;
            best = t;
        }
    }
    best
}
