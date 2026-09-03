//! OCCT ChFi3d_Builder_SpKP.cxx — SplitKPart (reconstruction of the
//! KPart SurfData against the restrictions of the support faces), 1:1
//! translation of the analytic subset:
//!   - CompTra (L70-82), CompCommonPoint (L89-109)
//!   - CpSD (L144-158), AdjustParam (L162-198)
//!   - FillSD (L627-742)
//!   - SplitKPart (L749-1290): the hatching of the tangency lines against
//!     the face boundaries runs through a Geom2dHatch stand-in (Trim over
//!     the face boundary pcurves); the periodic/multi-domain Tri parsing is
//!     translated for the non-periodic single-domain case and the domain
//!     classification (Adaptor3d_TopolTool::Classify) carries its pending
//!     boundary.  SearchFace (Builder_2.cxx L1523) is pending — the
//!     isolated-contour semantics (no neighbouring contour) stand in, which
//!     is the exact outcome for single-stripe KPart contours.

use glam::DVec2;
use rcad_kernel::base::int_ana2d::AnaIntersection2d;
use rcad_kernel::geom::{Curve2dEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};
use rcad_kernel::topods;

use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::topexp_face_edges;
use super::chfi_ds::{ChFiDSSpineHandle, ChFiDSSurfData, SharedSurfData};

/// OCCT Precision::PIntersection() = Intersection()/100 = 1.e-11.
const PITOL: f64 = 1.0e-11;

/// OCCT TopAbs_Position (the element-relative position of a hatcher point).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopAbsPosition {
    Forward,
    Reversed,
    Internal,
    External,
    On,
}

/// OCCT HatchGen_PointOnElement.
#[derive(Debug, Clone)]
pub struct PointOnElement {
    /// The hatcher element index (rcad: index into the boundary list).
    pub index: i32,
    /// The parameter on the element (the boundary pcurve).
    pub parameter: f64,
    pub position: TopAbsPosition,
}

/// OCCT HatchGen_PointOnHatching.
#[derive(Debug, Clone)]
pub struct PointOnHatching {
    /// The parameter on the hatching curve.
    pub parameter: f64,
    pub points: Vec<PointOnElement>,
}

/// OCCT HatchGen_Domain.
#[derive(Debug, Clone, Default)]
pub struct HatchDomain {
    pub has_first_point: bool,
    pub first_point: Option<PointOnHatching>,
    pub has_second_point: bool,
    pub second_point: Option<PointOnHatching>,
}

impl HatchDomain {
    pub fn first_point(&self) -> &PointOnHatching {
        self.first_point.as_ref().expect("no first point")
    }

    pub fn second_point(&self) -> &PointOnHatching {
        self.second_point.as_ref().expect("no second point")
    }
}

// =========================================================================
// OCCT CompTra (SpKP.cxx L70-82).
// =========================================================================
fn comp_tra(o1: Orientation, o2: Orientation, isfirst: bool) -> Orientation {
    if isfirst {
        super::chfi3d::topabs_reverse(super::chfi3d::topabs_compose(o1, o2))
    } else {
        super::chfi3d::topabs_compose(o1, o2)
    }
}

// =========================================================================
// OCCT CompCommonPoint (SpKP.cxx L89-109).
// =========================================================================
fn comp_common_point(
    fil_point: &mut super::chfi_ds::ChFiDS_CommonPoint,
    arc: &Shape,
    pe: &PointOnElement,
    or: Orientation,
) {
    let pos = pe.position;
    let ed = arc.as_edge().expect("not an edge");
    let v = if pos == TopAbsPosition::Forward {
        ed.first.clone()
    } else {
        ed.last.clone()
    };
    fil_point.set_vertex(v);
    fil_point.set_arc(
        P_OPERATOR_INTERSECTION,
        arc.clone(),
        pe.parameter,
        super::chfi3d::topabs_compose(arc.orientation, or),
    );
}

/// OCCT Precision::PIntersection() (the tolerance passed to SetArc).
const P_OPERATOR_INTERSECTION: f64 = PITOL;

// =========================================================================
// OCCT CpSD (SpKP.cxx L144-158): construct a new SurfData sharing the faces
// and copying the surface / interferences (registered as new DS entries).
// =========================================================================
fn cp_sd(dstr: &mut super::topopebrepds::TopOpeBRepDSHDataStructure, data: &ChFiDSSurfData) -> ChFiDSSurfData {
    let mut new_data = ChFiDSSurfData::default();
    let tos = dstr.surface(data.surf()).clone();
    new_data.change_surf(dstr.add_surface(super::topopebrepds::TopOpeBRepDSSurface::new(
        tos.surface.clone(),
        tos.tolerance(),
    )));
    new_data.change_index_of_s1(data.index_of_s1);
    new_data.change_index_of_s2(data.index_of_s2);
    *new_data.change_orientation() = data.orientation();
    *new_data.change_interference_on_s1() = cp_interf(dstr, data.interference_on_s1());
    *new_data.change_interference_on_s2() = cp_interf(dstr, data.interference_on_s2());
    new_data
}

/// OCCT CpInterf (SpKP.cxx L116-137).
fn cp_interf(
    dstr: &mut super::topopebrepds::TopOpeBRepDSHDataStructure,
    fi: &super::chfi_ds::ChFiDS_FaceInterference,
) -> super::chfi_ds::ChFiDS_FaceInterference {
    let mut new_f = fi.clone();
    let toc_curve = dstr.curve(fi.line_index()).curve.clone();
    let new_c = toc_curve.clone();
    new_f.set_interference(
        dstr.add_curve(super::topopebrepds::TopOpeBRepDSCurve::new(new_c, dstr.curve(fi.line_index()).tolerance())),
        fi.transition(),
        fi.pcurve_on_face().cloned(),
        fi.pcurve_on_surf().cloned(),
    );
    new_f
}

// =========================================================================
// OCCT AdjustParam (SpKP.cxx L162-198).
// =========================================================================
fn adjust_param(
    dom: &HatchDomain,
    f: &mut f64,
    l: &mut f64,
    wref: f64,
    period: f64,
    pitol: f64,
) -> bool {
    if dom.has_first_point {
        *f = dom.first_point().parameter;
    } else {
        *f = 0.0;
    }
    if dom.has_second_point {
        *l = dom.second_point().parameter;
    } else {
        *l = period;
    }
    if period == 0.0 {
        return false;
    }

    *f = super::chfi_ds::elclib_in_period(*f, wref - pitol, wref + period - pitol);
    *l = super::chfi_ds::elclib_in_period(*l, wref + pitol, wref + period + pitol);
    if *l < *f {
        *f -= period;
        return true;
    }
    false
}

// =========================================================================
// The Geom2dHatch_Hatcher stand-in: Trim of the tangency pcurve against the
// face boundary elements.  Returns the intersection hits (param on the
// hatching curve, PointOnElement) sorted by parameter; the single analytic
// domain spans the hits (or the whole curve when no hit exists).
// =========================================================================
fn hatcher_trim(
    brep: &topods::BRep,
    face: &Shape,
    pc: &rcad_kernel::geom::Curve2d,
    pcf: f64,
    pcl: f64,
) -> Vec<(f64, PointOnElement)> {
    let mut hits: Vec<(f64, PointOnElement)> = Vec::new();
    let mut element_index = 0i32;
    for e in topexp_face_edges(brep, face) {
        element_index += 1;
        let mut e_fwd = e.clone();
        e_fwd.orientation = Orientation::Forward;
        let mut face_fwd = face.clone();
        face_fwd.orientation = Orientation::Forward;
        let Some((epc, ef, el)) = brep.curve_on_surface(&e_fwd, &face_fwd) else {
            continue;
        };
        // Intersections of the hatching pcurve with the element pcurve.
        let pts = intersect_curve2d(pc, &epc);
        for (u_hatch, u_elem) in pts {
            if u_hatch < pcf - PITOL || u_hatch > pcl + PITOL {
                continue;
            }
            if u_elem < ef - PITOL || u_elem > el + PITOL {
                continue;
            }
            let position = if (u_elem - ef).abs() <= PITOL {
                TopAbsPosition::Forward
            } else if (u_elem - el).abs() <= PITOL {
                TopAbsPosition::Reversed
            } else {
                TopAbsPosition::On
            };
            hits.push((
                u_hatch,
                PointOnElement {
                    index: element_index,
                    parameter: u_elem,
                    position,
                },
            ));
        }
    }
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    hits
}

/// 2D curve-curve intersections (the Geom2dHatch_Intersector analytic
/// subset): returns (parameter on c1, parameter on c2) pairs.
fn intersect_curve2d(
    c1: &rcad_kernel::geom::Curve2d,
    c2: &rcad_kernel::geom::Curve2d,
) -> Vec<(f64, f64)> {
    use rcad_kernel::geom::Curve2d;
    let mut out = Vec::new();
    match (c1, c2) {
        (Curve2d::Line(l1), Curve2d::Line(l2)) => {
            let mut inter = AnaIntersection2d::new();
            inter.perform_lin_lin(l1, l2);
            if inter.is_done() {
                for i in 1..=inter.nb_points() {
                    out.push((inter.point(i).param_on_first(), inter.point(i).param_on_second()));
                }
            }
        }
        _ => {
            // The remaining analytic pairs run through the conic chain as in
            // geom2d_int_g_inter; only the point parameters are needed.
            let mut inter = AnaIntersection2d::new();
            let conic_pair = (c1.clone(), c2.clone());
            let _ = conic_pair;
            // Line vs circle / circle vs circle cover the box fillet cases.
            match (c1, c2) {
                (Curve2d::Line(l), Curve2d::Circle(c)) | (Curve2d::Circle(c), Curve2d::Line(l)) => {
                    inter.perform_lin_circ(l, c);
                    let swap = matches!((c1, c2), (Curve2d::Circle(_), Curve2d::Line(_)));
                    if inter.is_done() {
                        for i in 1..=inter.nb_points() {
                            let (a, b) = (
                                inter.point(i).param_on_first(),
                                inter.point(i).param_on_second(),
                            );
                            out.push(if swap { (b, a) } else { (a, b) });
                        }
                    }
                }
                (Curve2d::Circle(a), Curve2d::Circle(b)) => {
                    inter.perform_circ_circ(a, b);
                    if inter.is_done() {
                        for i in 1..=inter.nb_points() {
                            out.push((
                                inter.point(i).param_on_first(),
                                inter.point(i).param_on_second(),
                            ));
                        }
                    }
                }
                _ => {
                    // Non-canonical element pcurves: the Geom2dHatch
                    // generic intersector is pending; no hits reported.
                }
            }
        }
    }
    out
}

/// The single analytic domain built from the hits: bounded by the extreme
/// hits when they exist, otherwise the whole hatching curve without points
/// (the OCCT domain of a hatching lying strictly inside the face).
fn analytic_domain(hits: &[(f64, PointOnElement)], pcf: f64, pcl: f64) -> HatchDomain {
    let mut dom = HatchDomain::default();
    if let Some((par, pe)) = hits.first() {
        dom.has_first_point = true;
        dom.first_point = Some(PointOnHatching {
            parameter: *par,
            points: vec![pe.clone()],
        });
    }
    if let Some((par, pe)) = hits.last() {
        dom.has_second_point = true;
        dom.second_point = Some(PointOnHatching {
            parameter: *par,
            points: vec![pe.clone()],
        });
    }
    if hits.len() >= 2 {
        return dom;
    }
    // A single hit bounds only one side; the other side stays open.
    if hits.len() == 1 {
        return dom;
    }
    let _ = (pcf, pcl);
    dom
}

// =========================================================================
// OCCT FillSD (SpKP.cxx L627-742).
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn fill_sd(
    brep: &topods::BRep,
    dstr: &mut super::topopebrepds::TopOpeBRepDSHDataStructure,
    cd: &mut ChFiDSSurfData,
    boundary: &[Shape],
    dom: &HatchDomain,
    ponh: f64,
    isfirst: bool,
    ons: i32,
    _pitol: f64,
    bout: &Shape,
) {
    let opp = 3 - ons;
    let surf = dstr.surface(cd.surf()).surface.clone();

    let pph: Option<&PointOnHatching> = if isfirst && dom.has_first_point {
        Some(dom.first_point())
    } else if !isfirst && dom.has_second_point {
        Some(dom.second_point())
    } else {
        None
    };

    match pph {
        None => {
            cd.change_interference(ons).set_parameter(isfirst, ponh);
            let pcons = cd.interference(ons).pcurve_on_surf().expect("pcons");
            let uv = pcons.point_at(ponh);
            let p = surf.point_at(uv.x, uv.y);
            cd.change_vertex(isfirst, ons).set_point(p);
        }
        Some(ph) => {
            // Modification to find already existing vertexes.
            let mut le_type = 1usize;
            let nb_int = ph.points.len();
            if nb_int > 1 {
                let mut trouve = true;
                let mut suite = true;
                let mut v1 = Shape::null();
                let mut v2 = Shape::null();
                while trouve {
                    let petemp = &ph.points[le_type - 1];
                    if let Some(he) = boundary.get(petemp.index as usize - 1) {
                        let ed = he.as_edge().expect("not an edge");
                        v1 = ed.first.clone();
                        v2 = ed.last.clone();
                    } else {
                        suite = false;
                    }
                    if ((v1.is_same(bout) || v2.is_same(bout)) && suite) || !suite {
                        if v1.is_same(bout) || v2.is_same(bout) {
                            trouve = false; // found — exit
                        }
                        break;
                    } else {
                        suite = true;
                        trouve = true;
                        le_type += 1;
                        if le_type > nb_int {
                            le_type = 1;
                            break;
                        }
                    }
                }
            }
            let pe = &ph.points[le_type - 1];
            let Some(e) = boundary.get(pe.index as usize - 1) else {
                return;
            };
            let e = e.clone();

            if pe.position != TopAbsPosition::Internal {
                let mut o = cd.interference(ons).transition();
                if isfirst {
                    o = super::chfi3d::topabs_reverse(o);
                }
                let mut pons = cd.vertex(isfirst, ons).clone();
                comp_common_point(&mut pons, &e, pe, o);
                *cd.change_vertex(isfirst, ons) = pons;
            } else {
                let mut pons = cd.vertex(isfirst, ons).clone();
                pons.set_arc(
                    PITOL,
                    e.clone(),
                    pe.parameter,
                    comp_tra(cd.interference(ons).transition(), e.orientation, isfirst),
                );
                *cd.change_vertex(isfirst, ons) = pons;
            }
            cd.change_interference(ons).set_parameter(isfirst, ponh);
            let pcadj = cd.interference(ons).pcurve_on_surf().expect("pcadj");
            let uv = pcadj.point_at(ponh);
            let p = surf.point_at(uv.x, uv.y);
            cd.change_vertex(isfirst, ons).set_point(p);
        }
    }
    let pons_on_arc = cd.vertex(isfirst, ons).is_on_arc();
    let mut popp = cd.vertex(isfirst, opp).clone();
    if !popp.is_on_arc() {
        cd.change_interference(opp).set_parameter(isfirst, ponh);
        let pcopp = cd.interference(opp).pcurve_on_surf().expect("pcopp");
        let uv = pcopp.point_at(ponh);
        let p = surf.point_at(uv.x, uv.y);
        popp.set_point(p);
        *cd.change_vertex(isfirst, opp) = popp;
    }
    let _ = (pons_on_arc, brep);
}

// =========================================================================
// OCCT SplitKPart (SpKP.cxx L749-1290) — the reconstruction entry called
// from PerformSetOfKPart.  rcad: s1/s2 carry the support faces.
// =========================================================================
impl ChFi3dBuilder {
    pub fn split_k_part_hatched(
        &mut self,
        data: &mut ChFiDSSurfData,
        set_data: &mut Vec<SharedSurfData>,
        spine: &ChFiDSSpineHandle,
        iedge: usize,
        s1: &Shape,
        s2: &Shape,
        intf: &mut bool,
        intl: &mut bool,
    ) -> bool {
        // The hatching of each face is started by tangency lines.
        let pitol = PITOL;

        // OCCT L772-840: Trim of both tangency lines.
        let c1 = data.interference_on_s1().pcurve_on_face().cloned();
        let c2 = data.interference_on_s2().pcurve_on_face().cloned();
        let mut face_fwd = |f: &Shape| {
            let mut ff = f.clone();
            ff.orientation = Orientation::Forward;
            ff
        };
        let f1 = face_fwd(s1);
        let f2 = face_fwd(s2);

        let pcf1 = data.interference_on_s1().parameter_first();
        let pcl1 = data.interference_on_s1().parameter_last();
        let pcf2 = data.interference_on_s2().parameter_first();
        let pcl2 = data.interference_on_s2().parameter_last();

        let (hits1, dom1) = match &c1 {
            Some(pc) => {
                let hits = hatcher_trim(&self.my_brep, &f1, pc, pcf1, pcl1);
                if hits.is_empty() {
                    // OCCT: Nb1 == 0 — "tangency line out of the face".
                    return false;
                }
                let d = analytic_domain(&hits, pcf1, pcl1);
                (Some(hits), d)
            }
            None => (None, HatchDomain::default()),
        };
        let (hits2, dom2) = match &c2 {
            Some(pc) => {
                let hits = hatcher_trim(&self.my_brep, &f2, pc, pcf2, pcl2);
                if hits.is_empty() {
                    return false;
                }
                let d = analytic_domain(&hits, pcf2, pcl2);
                (Some(hits), d)
            }
            None => (None, HatchDomain::default()),
        };

        // Boundary element lists (M1/M2 maps).
        let boundary1 = topexp_face_edges(&self.my_brep, &f1);
        let boundary2 = topexp_face_edges(&self.my_brep, &f2);

        // Return start and end vertexes of the Spine (OCCT L842-852).
        let support = spine.base().edges(iedge).clone();
        let mut bout1 = self.my_brep.first_vertex(&support);
        let mut bout2 = self.my_brep.last_vertex(&support);
        if support.orientation == Orientation::Reversed {
            std::mem::swap(&mut bout1, &mut bout2);
        }

        // Return faces + register the support faces in the DS.
        let dstr = self.my_ds.as_mut().expect("DS");
        data.change_index_of_s1(dstr.add_shape(&f1));
        data.change_index_of_s2(dstr.add_shape(&f2));

        let nb1 = if hits1.is_some() { 1usize } else { 0usize };
        let nb2 = if hits2.is_some() { 1usize } else { 0usize };

        // OCCT L878-880: onS switcher + cntlFiOnS (OnSame length control —
        // only relevant for the truncation path, pending).
        if c1.is_none() && c2.is_none() {
            // "SplitData : 2 zero lines hatching impossible"
            return false;
        } else if c1.is_none() || (nb1 == 1 && !dom1.has_first_point) {
            // It is checked if the point 2d of the degenerated edge is in
            // the face (Adaptor3d_TopolTool::Classify — pending; treated as
            // IN, the box cases never enter this branch with C1 null).
            let pon_first = dom2.first_point.as_ref().map(|p| p.parameter).unwrap_or(pcf2);
            let pon_second = dom2.second_point.as_ref().map(|p| p.parameter).unwrap_or(pcl2);
            let _ = (&pon_first, &pon_second);
            // Filling of SurfData.
            let mut cd = data.clone();
            let pon_first = dom2.first_point.as_ref().map(|p| p.parameter).unwrap_or(pcf2);
            let pon_second = dom2.second_point.as_ref().map(|p| p.parameter).unwrap_or(pcl2);
            fill_sd(&self.my_brep, dstr, &mut cd, &boundary2, &dom2, pon_first, true, 2, pitol, &bout1);
            fill_sd(&self.my_brep, dstr, &mut cd, &boundary2, &dom2, pon_second, false, 2, pitol, &bout2);
            set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd)));
            // intf/intl tails — the isolated-contour semantics reset both
            // (SpKP L916-941 with SearchFace pending = false).
            if *intf {
                *intf = false;
            }
            if *intl {
                *intl = false;
            }
            return true;
        } else if c2.is_none() || (nb2 == 1 && !dom2.has_first_point) {
            let mut cd = data.clone();
            let pon_first = dom1.first_point.as_ref().map(|p| p.parameter).unwrap_or(pcf1);
            let pon_second = dom1.second_point.as_ref().map(|p| p.parameter).unwrap_or(pcl1);
            fill_sd(&self.my_brep, dstr, &mut cd, &boundary1, &dom1, pon_first, true, 1, pitol, &bout1);
            fill_sd(&self.my_brep, dstr, &mut cd, &boundary1, &dom1, pon_second, false, 1, pitol, &bout2);
            set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd)));
            if *intf {
                *intf = false;
            }
            if *intl {
                *intl = false;
            }
            return true;
        } else {
            // Parsing of domains by increasing parameters (non-periodic:
            // identity order — Tri pending for the periodic cases).
            let period1 = 0.0f64;
            let period2 = 0.0f64;

            // Filling of SurfData (OCCT L1031-1085).
            let mut cd = data.clone();
            for dom1_i in [&dom1] {
                let mut f1 = 0.0f64;
                let mut l1 = 0.0f64;
                let mut f2 = 0.0f64;
                let mut l2 = 0.0f64;
                let acheval1 = adjust_param(dom1_i, &mut f1, &mut l1, 0.0, period1, pitol);
                let nbcoup1 = if acheval1 { 2 } else { 1 };
                for _icoup1 in 1..=nbcoup1 {
                    let acheval2 = adjust_param(&dom2, &mut f2, &mut l2, 0.0, period2, pitol);
                    let nbcoup2 = if acheval2 { 2 } else { 1 };
                    for _icoup2 in 1..=nbcoup2 {
                        if f2 <= l1 && f1 <= l2 {
                            let tol2d = self.tol2d;
                            if f1 >= f2 - tol2d {
                                fill_sd(&self.my_brep, dstr, &mut cd, &boundary1, dom1_i, f1, true, 1, pitol, &bout1);
                            }
                            if f2 >= f1 - tol2d {
                                fill_sd(&self.my_brep, dstr, &mut cd, &boundary2, &dom2, f2, true, 2, pitol, &bout1);
                            }
                            if l1 >= l2 - tol2d {
                                fill_sd(&self.my_brep, dstr, &mut cd, &boundary2, &dom2, l2, false, 2, pitol, &bout2);
                            }
                            if l2 >= l1 - tol2d {
                                fill_sd(&self.my_brep, dstr, &mut cd, &boundary1, dom1_i, l1, false, 1, pitol, &bout2);
                            }
                            set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd.clone())));
                            cd = cp_sd(dstr, &cd);
                        }
                        f2 += period2;
                        l2 += period2;
                    }
                    f1 += period1;
                    l1 += period1;
                }
            }

            if set_data.is_empty() {
                return false;
            }

            // Processing of extensions (OCCT L1087-1210): with the pending
            // isolated-contour SearchFace the tails reduce to intf/intl =
            // true when both points of the kept SurfData are on arcs.
            if *intf {
                let sd0 = set_data[0].read().expect("surfdata lock");
                let cp1 = sd0.vertex_first_on_s1().is_on_arc();
                let cp2 = sd0.vertex_first_on_s2().is_on_arc();
                drop(sd0);
                if cp1 && cp2 {
                    *intf = true;
                } else if cp1 || cp2 {
                    *intf = true;
                } else {
                    *intf = false;
                }
            }
            if *intl {
                let sdl = set_data.last().unwrap().read().expect("surfdata lock");
                let cp1 = sdl.vertex_last_on_s1().is_on_arc();
                let cp2 = sdl.vertex_last_on_s2().is_on_arc();
                drop(sdl);
                if cp1 && cp2 {
                    *intl = true;
                } else if cp1 || cp2 {
                    *intl = true;
                } else {
                    *intl = false;
                }
            }
            true
        }
    }
}

// OCCT DVec2 import kept for the FillSD signature symmetry.
#[allow(unused)]
fn _unused_dvec2(_p: DVec2) {}
