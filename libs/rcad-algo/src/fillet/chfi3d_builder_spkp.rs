//! OCCT ChFi3d_Builder_SpKP.cxx — SplitKPart (reconstruction of the
//! KPart SurfData against the restrictions of the support faces), 1:1
//! translation of the analytic subset:
//!   - CompTra (L70-82), CompCommonPoint (L89-109)
//!   - CpSD (L144-158), AdjustParam (L162-198)
//!   - Tri (L548-623), FillSD (L627-742)
//!   - SplitKPart (L749-1290): the hatching of the tangency lines against
//!     the face boundaries runs through the landed 1:1 Geom2dHatch_Hatcher
//!     (geomalgo/hatch/hatcher.rs) as the OCCT statements do: AddElement per
//!     face restriction, Trim, ComputeDomains, IsDone, NbDomains, then
//!     Domain(iH, Ind(i)) with the Tri ordering.
//!     The tails of the extension processing (SearchFace, ChFi3d_EdgeState)
//!     and the Adaptor3d_TopolTool::Classify of the degenerate-edge 2d point
//!     are still pending — see the branch comments.

use rcad_kernel::geom::{Curve2d, Curve2dEval as _, SurfaceEval as _, TrimmedCurve2};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::hatch::hatch_gen::{Domain as HatchDomain, PointOnElement, PointOnHatching};
use crate::geomalgo::hatch::hatcher::Hatcher;
use crate::geomalgo::hatch::intersector::HatchIntersector;

use super::chfi3d::ChFi3dBuilder;
use super::chfi3d_builder_0::topexp_face_edges;
use super::chfi_ds::{ChFiDSSpineHandle, ChFiDSSurfData, SharedSurfData};

/// OCCT Precision::PIntersection() = Intersection()/100 = 1.e-11.
const PITOL: f64 = 1.0e-11;

// The hatching types are the landed 1:1 HatchGen ones (geomalgo/hatch/
// hatch_gen.rs): HatchGen_Domain (HatchDomain), HatchGen_PointOnHatching
// (PointOnHatching) and HatchGen_PointOnElement (PointOnElement).

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
    let pos = pe.position();
    let ed = arc.as_edge().expect("not an edge");
    let v = if pos == Orientation::Forward {
        ed.first.clone()
    } else {
        ed.last.clone()
    };
    fil_point.set_vertex(v);
    fil_point.set_arc(
        P_OPERATOR_INTERSECTION,
        arc.clone(),
        pe.parameter(),
        super::chfi3d::topabs_compose(arc.orientation, or),
    );
}

/// OCCT Precision::PIntersection() (the tolerance passed to SetArc).
const P_OPERATOR_INTERSECTION: f64 = PITOL;

// =========================================================================
// OCCT CpSD (SpKP.cxx L144-158): construct a new SurfData sharing the faces
// and copying the surface / interferences (registered as new DS entries).
// =========================================================================
fn cp_sd(dstr: &mut super::chfi3d_ds::TopOpeBRepDSHDataStructure, data: &ChFiDSSurfData) -> ChFiDSSurfData {
    let mut new_data = ChFiDSSurfData::default();
    let tos = dstr.surface(data.surf()).clone();
    new_data.change_surf(dstr.add_surface(super::chfi3d_ds::TopOpeBRepDSSurface::new(
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
    dstr: &mut super::chfi3d_ds::TopOpeBRepDSHDataStructure,
    fi: &super::chfi_ds::ChFiDS_FaceInterference,
) -> super::chfi_ds::ChFiDS_FaceInterference {
    let mut new_f = fi.clone();
    let toc_curve = dstr.curve(fi.line_index()).curve.clone();
    let new_c = toc_curve.clone();
    new_f.set_interference(
        dstr.add_curve(super::chfi3d_ds::TopOpeBRepDSCurve::new(new_c, dstr.curve(fi.line_index()).tolerance())),
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
    if dom.has_first_point() {
        *f = dom.first_point().parameter();
    } else {
        *f = 0.0;
    }
    if dom.has_second_point() {
        *l = dom.second_point().parameter();
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
// OCCT Tri (SpKP.cxx L548-623): Ind holds the ranks of the domains of the
// IndH-th hatching sorted by increasing adjusted parameter; the domain
// without a first point is closed on the first point of the domain without a
// second point, shifted one period back.
// =========================================================================
fn tri(
    h: &mut Hatcher,
    ih: usize,
    ind: &mut [usize],
    wref: f64,
    period: f64,
    pitol: f64,
    nbdom: &mut usize,
) -> bool {
    // OCCT L556-561: Ind(i) = i.
    for i in 1..=*nbdom {
        ind[i - 1] = i;
    }
    let mut f1 = 0.0f64;
    let mut f2 = 0.0f64;
    let mut l = 0.0f64;
    // bool Invert = true; while (Invert) { Invert = false;
    //                                    for (i = 1; i < Nbdom; i++) {...} }
    let mut invert = true;
    while invert {
        invert = false;
        for i in 1..*nbdom {
            adjust_param(h.domain(ih, ind[i - 1]), &mut f1, &mut l, wref, period, pitol);
            adjust_param(h.domain(ih, ind[i]), &mut f2, &mut l, wref, period, pitol);
            if f2 < f1 {
                ind.swap(i - 1, i);
                invert = true;
            }
        }
    }

    // OCCT L583-598.
    let mut i_sans_first = 0usize;
    let mut i_sans_last = 0usize;
    if *nbdom != 1 {
        for i in 1..=*nbdom {
            if !h.domain(ih, ind[i - 1]).has_first_point() {
                i_sans_first = i;
            }
            if !h.domain(ih, ind[i - 1]).has_second_point() {
                i_sans_last = i;
            }
        }
    }
    if i_sans_first != 0 {
        if i_sans_last == 0 {
            // OCCT L603-606: "Parsing : Pb of Hatcher".
            return false;
        }
        // OCCT L608-614: the first point of the domain without a first point
        // is taken from the domain without a second point, one period back.
        // OCCT reaches the two domains through a const_cast of H.Domain();
        // rcad takes the mutable accessor of the hatching for the same
        // write-back.
        let new_par = h.hatching_curve(ih).first_parameter() - period
            + h.domain(ih, ind[i_sans_last - 1]).first_point().parameter();
        let mut ph = h.domain(ih, ind[i_sans_last - 1]).first_point().clone();
        ph.set_parameter(new_par);
        h.hatching(ih)
            .change_domain(ind[i_sans_last - 1])
            .set_first_point(&ph);
        h.hatching(ih)
            .change_domain(ind[i_sans_first - 1])
            .set_first_point(&ph);

        // OCCT L616-620: Ind(k) = Ind(k + 1) for k = iSansLast..Nbdom-1; Nbdom--.
        for k in i_sans_last..*nbdom {
            ind[k - 1] = ind[k];
        }
        *nbdom -= 1;
    }
    true
}

// =========================================================================
// OCCT FillSD (SpKP.cxx L627-742).
// =========================================================================
#[allow(clippy::too_many_arguments)]
fn fill_sd(
    dstr: &mut super::chfi3d_ds::TopOpeBRepDSHDataStructure,
    cd: &mut ChFiDSSurfData,
    m: &[Shape],
    dom: &HatchDomain,
    ponh: f64,
    isfirst: bool,
    ons: i32,
    pitol: f64,
    bout: &Shape,
) {
    let opp = 3 - ons;
    let surf = dstr.surface(cd.surf()).surface.clone();

    // OCCT L639-651: the point on hatching of the requested end of the domain.
    let pph: Option<&PointOnHatching> = if isfirst && dom.has_first_point() {
        Some(dom.first_point())
    } else if !isfirst && dom.has_second_point() {
        Some(dom.second_point())
    } else {
        None
    };

    match pph {
        None => {
            // OCCT L656-661.
            cd.change_interference(ons).set_parameter(isfirst, ponh);
            let pcons = cd.interference(ons).pcurve_on_surf().expect("pcons");
            let uv = pcons.point_at(ponh);
            let p = surf.point_at(uv.x, uv.y);
            cd.change_vertex(isfirst, ons).set_point(p);
        }
        Some(ph) => {
            // Modification to find already existing vertexes (L664-704).
            let mut le_type = 1usize;
            let nb_int = ph.nb_points();
            if nb_int > 1 {
                let mut trouve = true;
                let mut suite = true;
                let mut v1 = Shape::null();
                let mut v2 = Shape::null();
                while trouve {
                    let petemp = ph.point(le_type);
                    if let Some(he) = m.get(petemp.index() as usize - 1) {
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
            let pe = ph.point(le_type);
            let Some(e) = m.get(pe.index() as usize - 1) else {
                return;
            };
            let e = e.clone();

            if pe.position() != Orientation::Internal {
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
                    pitol,
                    e.clone(),
                    pe.parameter(),
                    comp_tra(cd.interference(ons).transition(), e.orientation, isfirst),
                );
                *cd.change_vertex(isfirst, ons) = pons;
            }
            let pcadj = cd.interference(ons).pcurve_on_surf().expect("pcadj");
            let uv = pcadj.point_at(ponh);
            let p = surf.point_at(uv.x, uv.y);
            cd.change_interference(ons).set_parameter(isfirst, ponh);
            cd.change_vertex(isfirst, ons).set_point(p);
        }
    }
    // OCCT L735-741.
    let mut popp = cd.vertex(isfirst, opp).clone();
    if !popp.is_on_arc() {
        cd.change_interference(opp).set_parameter(isfirst, ponh);
        let pcopp = cd.interference(opp).pcurve_on_surf().expect("pcopp");
        let uv = pcopp.point_at(ponh);
        let p = surf.point_at(uv.x, uv.y);
        popp.set_point(p);
        *cd.change_vertex(isfirst, opp) = popp;
    }
}

// =========================================================================
// OCCT SplitKPart (SpKP.cxx L749-1290) — the reconstruction entry called
// from PerformSetOfKPart.  rcad: s1/s2 carry the support faces (the OCCT
// S1/I1 and S2/I2 pairs — the topol tool of a BRep face only ever yields
// BRepAdaptor_Curve2d over the face edges, BRepTopAdaptor_TopolTool.cxx
// L86-94).
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
        // The hatching of each faces is started by tangency lines.

        let pitol = PITOL;

        // OCCT L764-771: the M1/M2 maps (element index -> restriction), the
        // hatching indices, the domain counts and the shared intersector.
        let mut m1: Vec<Shape> = Vec::new();
        let mut m2: Vec<Shape> = Vec::new();
        let mut ih1 = 0usize;
        let mut ih2 = 0usize;
        let mut nb1 = 1usize;
        let mut nb2 = 1usize;

        // Cutting of tangency lines (hatching).
        let mut h1 = Hatcher::new(
            HatchIntersector::with_tolerances(pitol, pitol),
            self.tol2d,
            self.tolapp3d,
            false,
            false,
        );
        let mut h2 = Hatcher::new(
            HatchIntersector::with_tolerances(pitol, pitol),
            self.tol2d,
            self.tolapp3d,
            false,
            false,
        );

        let face_fwd = |f: &Shape| {
            let mut ff = f.clone();
            ff.orientation = Orientation::Forward;
            ff
        };
        let face1 = face_fwd(s1);
        let face2 = face_fwd(s2);

        // OCCT L772-773: C1 = Data->InterferenceOnS1().PCurveOnFace().
        let c1 = data.interference_on_s1().pcurve_on_face().cloned();
        if let Some(c1) = &c1 {
            // OCCT L777-790: the elements are the restrictions of I1 (the
            // face edges), each one bound to the handle it comes from.
            for e in topexp_face_edges(&self.my_brep, &face1) {
                // OCCT L779-780: Bc = down_cast<BRepAdaptor_Curve2d>
                // (I1->Value()), Gc = down_cast<Geom2dAdaptor_Curve>
                // (I1->Value()).
                let Some((pc, first, last)) = self.my_brep.curve_on_surface(&e, &face1) else {
                    // OCCT BRepAdaptor_Curve2d leaves the adaptor unloaded when
                    // BRep_Tool::CurveOnSurface finds no pcurve (cxx L50-62);
                    // such a restriction carries no hatching element.
                    continue;
                };
                // OCCT Geom2dAdaptor_Curve::Load(C, First, Last) (cxx
                // L267-273) — the loaded range of the adaptor.  rcad's
                // Curve2d carries no range of its own, so the pcurve is
                // restricted to the loaded one.
                let bc = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(pc),
                    t_min: first,
                    t_max: last,
                });
                // OCCT L783-788: ie = H1.AddElement(*Bc, Bc->Edge().Orientation())
                // (the Gc branch is unreachable for a BRep face).
                let ie = h1.add_element(&bc, e.orientation);
                // OCCT L789: M1.Bind(ie, I1->Value()).
                if m1.len() < ie {
                    m1.resize(ie, Shape::null());
                }
                m1[ie - 1] = e.clone();
            }
            // OCCT L791-799: iH1 = H1.Trim(ll1); H1.ComputeDomains(iH1);
            // if (!H1.IsDone(iH1)) return false; Nb1 = H1.NbDomains(iH1);
            // if (Nb1 == 0) return false.
            ih1 = h1.trim_curve(c1);
            h1.compute_domains_hatching(ih1);
            if !h1.is_done(ih1) {
                return false;
            }
            nb1 = h1.nb_domains(ih1);
            if nb1 == 0 {
                // "SplitKPart : tangency line out of the face"
                return false;
            }
        }

        // OCCT L807-840: the same block for the second face.
        let c2 = data.interference_on_s2().pcurve_on_face().cloned();
        if let Some(c2) = &c2 {
            for e in topexp_face_edges(&self.my_brep, &face2) {
                let Some((pc, first, last)) = self.my_brep.curve_on_surface(&e, &face2) else {
                    continue;
                };
                let bc = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(pc),
                    t_min: first,
                    t_max: last,
                });
                let ie = h2.add_element(&bc, e.orientation);
                // OCCT L824: M2.Bind(ie, I2->Value()).
                if m2.len() < ie {
                    m2.resize(ie, Shape::null());
                }
                m2[ie - 1] = e.clone();
            }
            ih2 = h2.trim_curve(c2);
            h2.compute_domains_hatching(ih2);
            if !h2.is_done(ih2) {
                return false;
            }
            nb2 = h2.nb_domains(ih2);
            if nb2 == 0 {
                // "SplitKPart : tangency line out of the face"
                return false;
            }
        }

        // Return start and end vertexes of the Spine (OCCT L842-852).
        let support = spine.base().edges(iedge).clone();
        let mut bout1 = self.my_brep.first_vertex(&support);
        let mut bout2 = self.my_brep.last_vertex(&support);
        if support.orientation == Orientation::Reversed {
            std::mem::swap(&mut bout1, &mut bout2);
        }

        // Return faces + register the support faces in the DS.
        let dstr = self.my_ds.as_mut().expect("DS");
        // D6 routing: shape registry -> BOPDS DS (AppendShape); OCCT ChFi3d_Builder_SpKP.cxx L871
        data.change_index_of_s1(dstr.add_shape(&face1));
        // D6 routing: shape registry -> BOPDS DS (AppendShape); OCCT ChFi3d_Builder_SpKP.cxx L872
        data.change_index_of_s2(dstr.add_shape(&face2));

        // OCCT L874-876: the domain parameter cursors, the Ind arrays and wref.
        let mut f1 = 0.0f64;
        let mut l1 = 0.0f64;
        let mut f2 = 0.0f64;
        let mut l2 = 0.0f64;
        let mut ind1: Vec<usize> = vec![0; nb1];
        let mut ind2: Vec<usize> = vec![0; nb2];
        let wref = 0.0f64;

        // OCCT L878-880: onS switcher + cntlFiOnS (OnSame length control —
        // only relevant for the pending truncation path).
        if c1.is_none() && c2.is_none() {
            // "SplitData : 2 zero lines hatching impossible"
            return false;
        } else if c1.is_none() || (nb1 == 1 && !h1.domain(ih1, 1).has_first_point()) {
            // OCCT L891-900: "It is checked if the point 2d of the degenerated
            // edge is in the face" — CD->Get2dPoints(false, 1) plus
            // I1->Classify(p2d1, 1.e-8, false) (Adaptor3d_TopolTool::Classify
            // over BRepTopAdaptor_FClass2d) is still pending; the C1 non-null
            // case never enters the classification.
            //
            // Parsing of domains by increasing parameters (OCCT L902-906).
            if !tri(&mut h2, ih2, &mut ind2, wref, 0.0, pitol, &mut nb2) {
                return false;
            }
            // Filling of SurfData (OCCT L907-915).
            let mut cd = data.clone();
            for i in 1..=nb2 {
                let dom2 = h2.domain(ih2, ind2[i - 1]);
                fill_sd(dstr, &mut cd, &m2, dom2, dom2.first_point().parameter(), true, 2, pitol, &bout1);
                fill_sd(dstr, &mut cd, &m2, dom2, dom2.second_point().parameter(), false, 2, pitol, &bout2);
                set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd.clone())));
                cd = cp_sd(dstr, &cd);
            }
            // OCCT L916-941: the intf/intl tails run through SearchFace, which
            // is still pending — the rcad tail keeps the previous stand-in
            // (both flags reset).
            if *intf {
                *intf = false;
            }
            if *intl {
                *intl = false;
            }
            return true;
        } else if c2.is_none() || (nb2 == 1 && !h2.domain(ih2, 1).has_first_point()) {
            // OCCT L945-954: the same degenerate-edge check for the second
            // face (I2->Classify — pending).
            //
            // Parsing of domains by increasing parameters (OCCT L956-960).
            if !tri(&mut h1, ih1, &mut ind1, wref, 0.0, pitol, &mut nb1) {
                return false;
            }
            // Filling of SurfData (OCCT L961-969).
            let mut cd = data.clone();
            for i in 1..=nb1 {
                let dom1 = h1.domain(ih1, ind1[i - 1]);
                fill_sd(dstr, &mut cd, &m1, dom1, dom1.first_point().parameter(), true, 1, pitol, &bout1);
                fill_sd(dstr, &mut cd, &m1, dom1, dom1.second_point().parameter(), false, 1, pitol, &bout2);
                set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd.clone())));
                cd = cp_sd(dstr, &cd);
            }
            // OCCT L970-995: the intf/intl tails (SearchFace pending — same
            // stand-in as the first branch).
            if *intf {
                *intf = false;
            }
            if *intl {
                *intl = false;
            }
            return true;
        } else {
            // The else branch runs with both tangency lines: OCCT L997-1029.
            let (Some(ll1), Some(ll2)) = (c1.as_ref(), c2.as_ref()) else {
                return false;
            };

            // Parsing of domains by increasing parameters,
            // if there is a 2d circle on a plane, one goes on 2D line of
            // opposite face.
            let mut period1 = 0.0f64;
            let mut period2 = 0.0f64;
            // OCCT Geom2dAdaptor_Curve::IsPeriodic / Period of ll1, ll2.
            if Curve2dAdaptor::is_periodic(ll1) {
                if !tri(&mut h2, ih2, &mut ind2, wref, 0.0, pitol, &mut nb2) {
                    return false;
                }
                period1 = Curve2dAdaptor::period(ll1);
                if !tri(&mut h1, ih1, &mut ind1, wref, period1, pitol, &mut nb1) {
                    return false;
                }
            } else {
                if !tri(&mut h1, ih1, &mut ind1, wref, 0.0, pitol, &mut nb1) {
                    return false;
                }
                if Curve2dAdaptor::is_periodic(ll2) {
                    period2 = Curve2dAdaptor::period(ll2);
                }
                if !tri(&mut h2, ih2, &mut ind2, wref, period2, pitol, &mut nb2) {
                    return false;
                }
            }

            // Filling of SurfData (OCCT L1031-1085).
            let mut cd = data.clone();
            for i in 1..=nb1 {
                let dom1 = h1.domain(ih1, ind1[i - 1]);
                let mut nbcoup1 = 1usize;
                let acheval1 = adjust_param(dom1, &mut f1, &mut l1, wref, period1, pitol);
                if acheval1 {
                    nbcoup1 = 2;
                }
                for _icoup1 in 1..=nbcoup1 {
                    for j in 1..=nb2 {
                        // OCCT L1046: the second domain is read with the loop
                        // index (not with Ind2(j)).
                        let dom2 = h2.domain(ih2, j);
                        let mut nbcoup2 = 1usize;
                        let acheval2 = adjust_param(dom2, &mut f2, &mut l2, wref, period2, pitol);
                        if acheval2 {
                            nbcoup2 = 2;
                        }
                        for _icoup2 in 1..=nbcoup2 {
                            if f2 <= l1 && f1 <= l2 {
                                let tol2d = self.tol2d;
                                if f1 >= f2 - tol2d {
                                    fill_sd(dstr, &mut cd, &m1, dom1, f1, true, 1, pitol, &bout1);
                                }
                                if f2 >= f1 - tol2d {
                                    fill_sd(dstr, &mut cd, &m2, dom2, f2, true, 2, pitol, &bout1);
                                }
                                if l1 >= l2 - tol2d {
                                    fill_sd(dstr, &mut cd, &m2, dom2, l2, false, 2, pitol, &bout2);
                                }
                                if l2 >= l1 - tol2d {
                                    fill_sd(dstr, &mut cd, &m1, dom1, l1, false, 1, pitol, &bout2);
                                }
                                set_data.push(std::sync::Arc::new(std::sync::RwLock::new(cd.clone())));
                                cd = cp_sd(dstr, &cd);
                            }
                            f2 += period2;
                            l2 += period2;
                        }
                    }
                    f1 += period1;
                    l1 += period1;
                }
            }

            if set_data.is_empty() {
                return false;
            }

            // OCCT L1094-1227: the extension processing of the beginning of
            // the spine (SearchFace, ChFi3d_cherche_element, ChFi3d_EdgeState)
            // is pending; with no face search the tails reduce to the OnArc
            // state of the kept SurfData ends.
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
            // OCCT L1228-1330: the extension processing of the end of the
            // spine (pending, same reduction).
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
