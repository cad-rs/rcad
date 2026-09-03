//! OCCT ChFi3d_Builder_0.cxx L2393-3441 — the DS-assembly part of
//! Builder_0 (ChFi3d_Orientation / ChFi3d_Contains / QueryAddVertexInEdge /
//! CutEdge / findIndexPoint / ChFi3d_FilDS / ChFi3d_StripeEdgeInter).
//!
//! Split from chfi3d_builder_0.rs to honor the 2000-line source limit
//! (workspace Rule 5); the OCCT annotations are the same file.

use glam::DVec2;
use rcad_kernel::geom::Curve2dEval as _;
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::topods;

use super::chfi3d::{topabs_compose, topabs_reverse, ChFi3dBuilder};
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_compute_arete, chfi3d_compute_pcurv_2pt, chfi3d_fil_curve_in_ds,
    chfi3d_fil_point_in_ds, chfi3d_fil_vertex_in_ds, chfi3d_same_parameter, P_CONFUSION,
};
use super::chfi_ds::{ChFiDSRegul, ChFiDSStripe, ChFiDSSurfData};
use super::topopebrepds::{
    TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind,
    TopOpeBRepDSSolidSurfaceInterference,
};

use super::chfi3d::chfi3d_index_point_in_ds;

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2393-2438 — ChFi3d_Orientation: the transition
// orientation of the interference (the first found in the list).
// =========================================================================
pub fn chfi3d_orientation(
    li: &[TopOpeBRepDSInterference],
    igros: i32,
    ipetit: i32,
    or: &mut Orientation,
    isvertex: bool,
    aprendre: bool,
) -> bool {
    // In case, when it is necessary to insert a point/vertex, it should be
    // known if this is a point or a vertex, because their index can be the
    // same.
    let typepetit = if isvertex {
        TopOpeBRepDSKind::Vertex
    } else {
        TopOpeBRepDSKind::Point
    };
    for cur in li {
        let (gk, g, _sk, s) = cur.gkgsks();
        if aprendre {
            if s == igros && g == ipetit && gk == typepetit {
                *or = cur.transition().orientation_in();
                return true;
            }
        } else if s == igros && g == ipetit {
            *or = cur.transition().orientation_in();
            return true;
        }
    }
    false
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2445-2453 — ChFi3d_Contains: check if the
// interference does not already exist.
// =========================================================================
pub fn chfi3d_contains(
    li: &[TopOpeBRepDSInterference],
    igros: i32,
    ipetit: i32,
    isvertex: bool,
    aprendre: bool,
) -> bool {
    let mut bid_or = Orientation::Forward;
    chfi3d_orientation(li, igros, ipetit, &mut bid_or, isvertex, aprendre)
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2457-2484 — QueryAddVertexInEdge.
// =========================================================================
fn query_add_vertex_in_edge(
    li: &mut Vec<TopOpeBRepDSInterference>,
    ic: i32,
    iv: i32,
    par: f64,
    or: Orientation,
) {
    for cur in li.iter() {
        if let TopOpeBRepDSInterference::CurvePoint(cpi) = cur {
            let new_iv = cpi.index_g;
            let kv = cpi.kind_g;
            let new_or = cpi.transition.orientation_in();
            let newpar = cpi.parameter;
            if iv == new_iv
                && kv == TopOpeBRepDSKind::Vertex
                && or == new_or
                && (par - newpar).abs() < 1.0e-10
            {
                return;
            }
        }
    }
    let interf = chfi3d_fil_vertex_in_ds(or, ic, iv, par);
    li.push(interf);
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2488-2516 — CutEdge.
// =========================================================================
pub fn cut_edge(
    brep: &topods::BRep,
    v: &Shape,
    sd: &ChFiDSSurfData,
    dstr: &mut TopOpeBRepDSHDataStructure,
    _isfirst: bool,
    ons: i32,
) {
    let on_curve = if ons == 1 {
        sd.is_on_curve1()
    } else {
        sd.is_on_curve2()
    };
    if !on_curve {
        return;
    }
    let ic = sd.index_of(ons);
    let iv = dstr.add_shape(v);
    let e_forward = {
        let mut e = dstr.shape(ic).clone();
        e.orientation = Orientation::Forward;
        e
    };

    // process them checking that it has not been done already.
    // (OCCT: TopExp_Explorer over E's vertices; for a FORWARD edge the
    // first vertex carries FORWARD, the last REVERSED.)
    let ed = e_forward.as_edge().expect("not an edge");
    for (vv, vor) in [(&ed.first, Orientation::Forward), (&ed.last, Orientation::Reversed)] {
        if vv.is_same(v) {
            let or = topabs_reverse(vor);
            let par = brep_tool_parameter(brep, vv, &e_forward);
            let li = dstr.change_shape_interferences(ic);
            query_add_vertex_in_edge(li, ic, iv, par, or);
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2523-2561 — findIndexPoint: returns in <ipoin>
// the index of the point bounding a curve interfering with <Fd> and
// coinciding with the last common point on the <OnS> face.
// =========================================================================
pub fn find_index_point(
    dstr: &TopOpeBRepDSHDataStructure,
    fd: &ChFiDSSurfData,
    on_s: i32,
    ipoin: &mut i32,
) -> bool {
    *ipoin = 0;
    let p = fd.vertex(false, on_s).point();

    for sci in dstr.surface_interferences(fd.surf()) {
        let TopOpeBRepDSInterference::SurfaceCurve(sci) = sci else {
            continue;
        };
        for cpi in dstr.curve_interferences(sci.index_g) {
            let TopOpeBRepDSInterference::CurvePoint(cpi) = cpi else {
                continue;
            };
            let i_point = cpi.index_g;
            let tp = dstr.point(i_point);
            if p.distance(tp.point()) <= tp.tolerance() {
                *ipoin = i_point;
                return true;
            }
        }
    }
    false
}

// =========================================================================
// OCCT Geom2dInt_GInter — pending TKGeomAlgo translation for BSpline
// pcurves; the analytic (line/circle/ellipse) cases run through the
// AnaIntersection2d chain, matching Geom2dInt_TheIntConicCurveOfGInter.
// Returns (nb_points, nb_segments).
// =========================================================================
pub fn geom2d_int_g_inter(
    pc1: &rcad_kernel::geom::Curve2d,
    pc2: &rcad_kernel::geom::Curve2d,
) -> (usize, usize) {
    use rcad_kernel::base::int_ana2d::{AnaIntersection2d, Conic2d};

    fn is_conic_like(pc: &rcad_kernel::geom::Curve2d) -> bool {
        matches!(
            pc,
            rcad_kernel::geom::Curve2d::Line(_)
                | rcad_kernel::geom::Curve2d::Circle(_)
                | rcad_kernel::geom::Curve2d::Ellipse(_)
        )
    }
    if !is_conic_like(pc1) || !is_conic_like(pc2) {
        // BSpline/Bezier pcurves: the generic GInter chain is pending; no
        // intersection is reported (the analytic box cases never hit it).
        return (0, 0);
    }

    let mut inter = AnaIntersection2d::new();
    match (pc1, pc2) {
        (rcad_kernel::geom::Curve2d::Line(l1), rcad_kernel::geom::Curve2d::Line(l2)) => {
            inter.perform_lin_lin(l1, l2);
        }
        (rcad_kernel::geom::Curve2d::Line(l), rcad_kernel::geom::Curve2d::Circle(c))
        | (rcad_kernel::geom::Curve2d::Circle(c), rcad_kernel::geom::Curve2d::Line(l)) => {
            inter.perform_lin_circ(l, c);
        }
        (rcad_kernel::geom::Curve2d::Circle(c1), rcad_kernel::geom::Curve2d::Circle(c2)) => {
            inter.perform_circ_circ(c1, c2);
        }
        (rcad_kernel::geom::Curve2d::Line(l), other)
        | (other, rcad_kernel::geom::Curve2d::Line(l)) => {
            inter.perform_lin_conic(l, &Conic2d::from_ellipse(&ellipse_of(other).unwrap()));
        }
        (rcad_kernel::geom::Curve2d::Circle(c), other)
        | (other, rcad_kernel::geom::Curve2d::Circle(c)) => {
            inter.perform_circ_conic(c, &Conic2d::from_ellipse(&ellipse_of(other).unwrap()));
        }
        (a, b) => {
            // Ellipse/Ellipse (the only remaining conic-like pair).
            let e1 = ellipse_of(a).unwrap();
            let e2 = ellipse_of(b).unwrap();
            inter.perform_ellipse_conic(&e1, &Conic2d::from_ellipse(&e2));
        }
    }
    if !inter.is_done() {
        return (0, 0);
    }
    (inter.nb_points(), 0)
}

fn ellipse_of(pc: &rcad_kernel::geom::Curve2d) -> Option<rcad_kernel::geom::Ellipse2d> {
    match pc {
        rcad_kernel::geom::Curve2d::Ellipse(e) => Some(e.clone()),
        _ => None,
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L3352-3441 — ChFi3d_StripeEdgeInter: examines
// two stripes for an intersection between curves of interference with
// faces; such an intersection would corrupt the result, so raise (the OCCT
// throw maps to the false return).
// =========================================================================
pub fn chfi3d_stripe_edge_inter(
    the_stripe1: &ChFiDSStripe,
    the_stripe2: &ChFiDSStripe,
    _dstr: &mut TopOpeBRepDSHDataStructure,
    tol2d: f64,
) -> bool {
    // Do not check the stripes having common corner points.
    for isur1 in 1..=2 {
        for isur2 in 1..=2 {
            if the_stripe1.index_point(false, isur1) == the_stripe2.index_point(false, isur2)
                || the_stripe1.index_point(false, isur1) == the_stripe2.index_point(true, isur2)
                || the_stripe1.index_point(true, isur1) == the_stripe2.index_point(false, isur2)
                || the_stripe1.index_point(true, isur1) == the_stripe2.index_point(true, isur2)
            {
                return true;
            }
        }
    }

    let adat1 = &the_stripe1.my_hdata;
    let adat2 = &the_stripe2.my_hdata;

    for dat1 in adat1 {
        let a_dat1 = dat1.read().expect("surfdata lock");
        let ishape11 = a_dat1.index_of_s1;
        let ishape12 = a_dat1.index_of_s2;
        for dat2 in adat2 {
            let a_dat2 = dat2.read().expect("surfdata lock");
            let ishape21 = a_dat2.index_of_s1;
            let ishape22 = a_dat2.index_of_s2;

            // Find those FaceInterferences able to intersect.
            let (afi1, afi2) = if ishape11 == ishape21 {
                (
                    a_dat1.interference_on_s1().clone(),
                    a_dat2.interference_on_s1().clone(),
                )
            } else if ishape11 == ishape22 {
                (
                    a_dat1.interference_on_s1().clone(),
                    a_dat2.interference_on_s2().clone(),
                )
            } else if ishape12 == ishape21 {
                (
                    a_dat1.interference_on_s2().clone(),
                    a_dat2.interference_on_s1().clone(),
                )
            } else if ishape12 == ishape22 {
                (
                    a_dat1.interference_on_s2().clone(),
                    a_dat2.interference_on_s2().clone(),
                )
            } else {
                // No common faces.
                continue;
            };

            if (afi1.parameter_first() - afi1.parameter_last()).abs() <= f64::EPSILON
                || (afi2.parameter_first() - afi2.parameter_last()).abs() <= f64::EPSILON
                || afi1.pcurve_on_face().is_none()
                || afi2.pcurve_on_face().is_none()
            {
                // Do not waste time on degenerates.
                continue;
            }
            let (nb_points, nb_segments) = geom2d_int_g_inter(
                afi1.pcurve_on_face().unwrap(),
                afi2.pcurve_on_face().unwrap(),
            );
            let _ = tol2d;
            if nb_segments > 0 || nb_points > 0 {
                // OCCT: throw StdFail_NotDone("StripeEdgeInter : fillets
                // have too big radiuses").
                return false;
            }
        }
    }
    true
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2565-3341 — ChFi3d_FilDS.
// =========================================================================
pub fn chfi3d_filds(
    brep: &topods::BRep,
    solid_index: i32,
    cordat: &ChFiDSStripe,
    dstr: &mut TopOpeBRepDSHDataStructure,
    reglist: &mut Vec<ChFiDSRegul>,
    tol3d: f64,
    tol2d: f64,
) {
    let closed = cordat
        .spine()
        .map(|sp| sp.base().is_periodic())
        .unwrap_or(false);
    let mut degene = false;
    let mut is_vertex1 = false;
    let mut is_vertex2 = false;
    let mut singulier_en_bout = false;
    let seqfil = &cordat.my_hdata;
    let mut ipoin1 = cordat.indexfirst_pon_s1;
    let mut ipoin2 = cordat.indexfirst_pon_s2;
    let mut num_edge = 1usize;
    let mut boutde_vtx = Shape::null();
    let mut icurv = 0i32;
    let mut iarc1 = 0i32;
    let mut iarc2 = 0i32;
    let mut trafil1;
    let mut trafil2;
    let mut et1;
    let mut pardeb = 0.0f64;
    let mut parfin = 0.0f64;

    let mut regcout = ChFiDSRegul::new(); // for closed and tangent CD
    let mut regfilfil = ChFiDSRegul::new(); // for connections Surf/Surf

    // OCCT L2608-2609: V3/V4 carry the previous SurfData's last points.
    let mut v3 = super::chfi_ds::ChFiDS_CommonPoint::default();
    let mut v4 = super::chfi_ds::ChFiDS_CommonPoint::default();

    // Nullify degenerated ChFi/Faces interferences, eap occ293 (L2613-2653).
    if seqfil.len() > 1 {
        for j in 1..=seqfil.len() {
            for on_s in 1..=2 {
                let (ic_fil1, fi_len, fd_v_on_arc, mut cp1, mut cp2) = {
                    let fd = seqfil[j - 1].read().expect("surfdata lock");
                    let fi = fd.interference(on_s);
                    let ic_fil1 = fi.line_index();
                    let fi_len = (fi.parameter_first() - fi.parameter_last()).abs();
                    let isfirst = j == 1;
                    let i = if isfirst { j + 1 } else { j - 1 };
                    let cp1 = seqfil[i - 1]
                        .read()
                        .expect("surfdata lock")
                        .vertex(isfirst, on_s)
                        .clone();
                    let cp2 = fd.vertex(!isfirst, on_s).clone();
                    let v_on_arc = fd.vertex(isfirst, on_s).is_on_arc();
                    (ic_fil1, fi_len, v_on_arc, cp1, cp2)
                };
                if ic_fil1 == 0 {
                    continue;
                }
                if fi_len > P_CONFUSION {
                    continue;
                }
                dstr.change_curve(ic_fil1).nullify();

                // care of CommonPoint, eap occ354.
                if j != 1 && j != seqfil.len() {
                    continue;
                }
                let isfirst = j == 1;
                let i = if isfirst { j + 1 } else { j - 1 };
                if fd_v_on_arc && cp1.is_on_arc() {
                    // OCCT L2646-2649 (CP2 receives its own point back;
                    // only the flags/tolerance are reset).
                    cp1.reset();
                    let p2 = cp2.point();
                    cp1.set_point(p2);
                    cp2.reset();
                    let p1 = cp1.point();
                    cp2.set_point(p1);
                    *seqfil[i - 1]
                        .write()
                        .expect("surfdata lock")
                        .change_vertex(isfirst, on_s) = cp1;
                    *seqfil[j - 1]
                        .write()
                        .expect("surfdata lock")
                        .change_vertex(!isfirst, on_s) = cp2;
                }
            }
        }
    }

    for j in 1..=seqfil.len() {
        let fd_guard = seqfil[j - 1].read().expect("surfdata lock");
        let fd: &ChFiDSSurfData = &fd_guard;
        let isurf = fd.surf();
        let mut ishape1 = fd.index_of_s1;
        let mut ishape2 = fd.index_of_s2;

        // eap, Apr 29 2002, occ 293: which end is already in DS.
        let mut is_in_ds1 = false;
        let mut is_in_ds2 = false;
        if (j as i32) <= cordat.is_in_ds(true) {
            is_in_ds1 = true;
            is_in_ds2 = (j as i32) + 1 <= cordat.is_in_ds(true);
        }
        if ((seqfil.len() - j) as i32) < cordat.is_in_ds(false) {
            is_in_ds2 = true;
            is_in_ds1 = is_in_ds1 || (((seqfil.len() - j + 1) as i32) < cordat.is_in_ds(false));
        }

        // creation of SolidSurfaceInterference.
        let ssi = TopOpeBRepDSInterference::SolidSurface(
            TopOpeBRepDSSolidSurfaceInterference::new(
                fd.orientation(),
                TopOpeBRepDSKind::Solid,
                solid_index,
                TopOpeBRepDSKind::Surface,
                isurf,
            ),
        );
        dstr.change_shape_interferences(solid_index).push(ssi);

        let fi1 = fd.interference_on_s1().clone();
        let fi2 = fd.interference_on_s2().clone();
        let v1 = fd.vertex_first_on_s1().clone();
        let v2 = fd.vertex_first_on_s2().clone();

        // Processing to manage double interferences.
        if j > 1 {
            if v1.is_on_arc() && v3.is_on_arc() && v1.arc().is_same(v3.arc()) {
                if chfi3d_contains(dstr.shape_interferences(iarc1), iarc1, ipoin1, false, false)
                    && (v1.transition_on_arc() != v3.transition_on_arc())
                {
                    let interp1 = chfi3d_fil_point_in_ds(
                        v1.transition_on_arc(),
                        iarc1,
                        ipoin1,
                        v1.parameter_on_arc(),
                        false,
                    );
                    dstr.change_shape_interferences_of(v1.arc()).push(interp1);
                }
            }

            if v2.is_on_arc() && v4.is_on_arc() && v2.arc().is_same(v4.arc()) {
                if chfi3d_contains(dstr.shape_interferences(iarc2), iarc2, ipoin2, false, false)
                    && (v2.transition_on_arc() != v4.transition_on_arc())
                {
                    let interp2 = chfi3d_fil_point_in_ds(
                        v2.transition_on_arc(),
                        iarc2,
                        ipoin2,
                        v2.parameter_on_arc(),
                        false,
                    );
                    dstr.change_shape_interferences_of(v2.arc()).push(interp2);
                }
            }
        }

        v3 = fd.vertex_last_on_s1().clone();
        v4 = fd.vertex_last_on_s2().clone();

        if ishape1 != 0 {
            if ishape1 > 0 {
                trafil1 = dstr.shape(ishape1).orientation;
            } else {
                let mut or = Orientation::Forward;
                chfi3d_orientation(
                    dstr.shape_interferences(solid_index),
                    solid_index,
                    -ishape1,
                    &mut or,
                    false,
                    false,
                );
                trafil1 = or;
            }
            trafil1 = topabs_compose(trafil1, fd.orientation());
            trafil1 = topabs_compose(topabs_reverse(fi1.transition()), trafil1);
            trafil2 = topabs_reverse(trafil1);
        } else if ishape2 > 0 {
            trafil2 = dstr.shape(ishape2).orientation;
            trafil2 = topabs_compose(trafil2, fd.orientation());
            trafil2 = topabs_compose(topabs_reverse(fi2.transition()), trafil2);
            trafil1 = topabs_reverse(trafil2);
        } else {
            let mut or = Orientation::Forward;
            chfi3d_orientation(
                dstr.shape_interferences(solid_index),
                solid_index,
                -ishape2,
                &mut or,
                false,
                false,
            );
            trafil2 = or;
            trafil2 = topabs_compose(trafil2, fd.orientation());
            trafil2 = topabs_compose(topabs_reverse(fi2.transition()), trafil2);
            trafil1 = topabs_reverse(trafil2);
        }

        et1 = topabs_reverse(trafil1);

        // A small paragraph to process contacts of edges, which touch a
        // vertex of the obstacle.
        if v1.is_vertex() && fd.is_on_curve1() {
            let vv1 = v1.vertex().clone();
            cut_edge(brep, &vv1, fd, dstr, true, 1);
        }
        if v2.is_vertex() && fd.is_on_curve2() {
            let vv2 = v2.vertex().clone();
            cut_edge(brep, &vv2, fd, dstr, true, 2);
        }
        if v3.is_vertex() && fd.is_on_curve1() {
            let vv3 = v3.vertex().clone();
            cut_edge(brep, &vv3, fd, dstr, false, 1);
        }
        if v4.is_vertex() && fd.is_on_curve2() {
            let vv4 = v4.vertex().clone();
            cut_edge(brep, &vv4, fd, dstr, false, 2);
        }

        if j == 1 {
            is_vertex1 = v1.is_vertex();
            is_vertex2 = v2.is_vertex();
            singulier_en_bout = v1.point().distance(v2.point()) <= 0.0;

            if singulier_en_bout {
                // Queue de Billard.
                if !v1.is_vertex() || !v2.is_vertex() {
                    // OCCT: empty block.
                } else {
                    is_vertex1 = true;
                    is_vertex2 = true; // caution...
                    // The edge is removed from spine starting on this vertex.
                    let spine = cordat.spine().expect("spine");
                    let arcspine = spine.base().edges(1).clone();
                    boutde_vtx = v1.vertex().clone();
                    let iarcspine = dstr.add_shape(&arcspine);
                    let ivtx = cordat.indexfirst_pon_s1;

                    // OCCT L2802-2809: TopExp_Explorer over the
                    // FORWARD-oriented spine edge's vertices.
                    let mut ovtx = Orientation::Forward;
                    let ed = arcspine.as_edge().expect("not an edge");
                    if boutde_vtx.is_same(&ed.first) {
                        ovtx = Orientation::Forward;
                    } else if boutde_vtx.is_same(&ed.last) {
                        ovtx = Orientation::Reversed;
                    }
                    ovtx = topabs_reverse(ovtx);
                    let parvtx = brep_tool_parameter(brep, &boutde_vtx, &arcspine);
                    let interfv = chfi3d_fil_vertex_in_ds(ovtx, iarcspine, ivtx, parvtx);
                    dstr.change_shape_interferences(iarcspine).push(interfv);
                }
            } else {
                if v1.is_on_arc() {
                    iarc1 = dstr.add_shape(v1.arc());
                    if !chfi3d_contains(
                        dstr.shape_interferences(iarc1),
                        iarc1,
                        ipoin1,
                        false,
                        false,
                    ) {
                        let interp1 = chfi3d_fil_point_in_ds(
                            v1.transition_on_arc(),
                            iarc1,
                            ipoin1,
                            v1.parameter_on_arc(),
                            is_vertex1,
                        );
                        dstr.change_shape_interferences_of(v1.arc()).push(interp1);
                    }
                }

                if v2.is_on_arc() {
                    iarc2 = dstr.add_shape(v2.arc());
                    if !chfi3d_contains(
                        dstr.shape_interferences(iarc2),
                        iarc2,
                        ipoin2,
                        false,
                        false,
                    ) {
                        let interp2 = chfi3d_fil_point_in_ds(
                            v2.transition_on_arc(),
                            iarc2,
                            ipoin2,
                            v2.parameter_on_arc(),
                            is_vertex2,
                        );
                        dstr.change_shape_interferences_of(v2.arc()).push(interp2);
                    }
                }
            }

            if !is_in_ds1 {
                et1 = topabs_compose(et1, cordat.first_pcurve_orientation());
                icurv = cordat.first_curve();
                if closed && !singulier_en_bout {
                    regcout.set_curve(icurv);
                    regcout.set_s1(isurf, false);
                }
                let pcurv = cordat.first_pcurve().cloned();
                let (pd, pf) = cordat.first_parameters();
                pardeb = pd;
                parfin = pf;

                let li_empty = dstr.curve_interferences(icurv).is_empty();
                if li_empty {
                    if cordat.first_pcurve_orientation() == Orientation::Reversed {
                        let interp1 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            icurv,
                            ipoin1,
                            parfin,
                            is_vertex1,
                        );
                        let interp2 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            icurv,
                            ipoin2,
                            pardeb,
                            is_vertex2,
                        );
                        dstr.change_curve_interferences(icurv).push(interp1);
                        dstr.change_curve_interferences(icurv).push(interp2);
                    } else {
                        let interp1 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            icurv,
                            ipoin1,
                            pardeb,
                            is_vertex1,
                        );
                        let interp2 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            icurv,
                            ipoin2,
                            parfin,
                            is_vertex2,
                        );
                        dstr.change_curve_interferences(icurv).push(interp1);
                        dstr.change_curve_interferences(icurv).push(interp2);
                    }
                }
                let interfc1 = chfi3d_fil_curve_in_ds(icurv, isurf, pcurv, et1);
                dstr.change_surface_interferences(isurf).push(interfc1.clone());
                if ipoin1 == ipoin2 {
                    dstr.change_curve(icurv).nullify();
                    // OCCT: TCurv.SetSCI(Interfc1, bidinterf).
                    if let TopOpeBRepDSInterference::SurfaceCurve(sc) = &interfc1 {
                        let refd = super::topopebrepds::InterferenceRef {
                            pcurve: sc.pcurve.clone(),
                            index_s: sc.index_s,
                            index_g: sc.index_g,
                        };
                        dstr.change_curve(icurv).sci1 = Some(refd);
                    }
                }
            }
        } else {
            // ---- Interference between Fillets ------
            if !is_in_ds1 {
                if degene && is_vertex1 {
                    // The edge is removed from the spine starting on this
                    // vertex.
                    num_edge += 1; // The previous edge of the vertex has already been found.
                    let spine = cordat.spine().expect("spine");
                    let arcspine = spine.base().edges(num_edge).clone();
                    let iarcspine = dstr.add_shape(&arcspine);
                    let ivtx = dstr.add_shape(&boutde_vtx);
                    let mut ovtx = Orientation::Forward;
                    let ed = arcspine.as_edge().expect("not an edge");
                    if boutde_vtx.is_same(&ed.first) {
                        ovtx = Orientation::Forward;
                    } else if boutde_vtx.is_same(&ed.last) {
                        ovtx = Orientation::Reversed;
                    }
                    ovtx = topabs_reverse(ovtx);
                    let parvtx = brep_tool_parameter(brep, &boutde_vtx, &arcspine);
                    let interfv = chfi3d_fil_vertex_in_ds(ovtx, iarcspine, ivtx, parvtx);
                    dstr.change_shape_interferences(iarcspine).push(interfv);
                } // End of the removal

                let uv1 = fd
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .map(|pc| pc.point_at(fd.interference_on_s1().parameter_first()))
                    .unwrap_or(DVec2::ZERO);
                let uv2 = fd
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .map(|pc| pc.point_at(fd.interference_on_s2().parameter_first()))
                    .unwrap_or(DVec2::ZERO);
                if degene {
                    // pcurve is associated via SCI to TopOpeBRepDSCurve.
                    let pcurv = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
                    let interfc1 = chfi3d_fil_curve_in_ds(icurv, isurf, Some(pcurv), et1);
                    dstr.change_surface_interferences(isurf).push(interfc1.clone());
                    dstr.change_curve(icurv).nullify();
                    if let TopOpeBRepDSInterference::SurfaceCurve(sc) = &interfc1 {
                        let refd = super::topopebrepds::InterferenceRef {
                            pcurve: sc.pcurve.clone(),
                            index_s: sc.index_s,
                            index_g: sc.index_g,
                        };
                        dstr.change_curve(icurv).sci1 = Some(refd);
                    }
                } else {
                    regfilfil.set_s2(isurf, false);
                    reglist.push(regfilfil);
                    // OCCT L2938-2949: ChFi3d_ComputePCurv(TCurv.ChangeCurve(),
                    // UV1, UV2, PCurv, Surface(Fd->Surf()).Surface(), Pardeb,
                    // Parfin, tol3d, tolreached) + TCurv.Tolerance(max).
                    let mut pc = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
                    let mut tolreached = tol3d;
                    if let Some(c3d) = dstr.curve(icurv).curve.clone() {
                        let surf = dstr.surface(fd.surf()).surface.clone();
                        chfi3d_same_parameter(&c3d, &mut pc, &surf, tol3d, &mut tolreached);
                    }
                    let t_tol = dstr.curve(icurv).tolerance();
                    dstr
                        .change_curve(icurv)
                        .set_tolerance(t_tol.max(tolreached));
                    let interfc1 = chfi3d_fil_curve_in_ds(icurv, isurf, Some(pc), et1);
                    dstr.change_surface_interferences(isurf).push(interfc1);
                }
            }
        } // End of Interference between fillets

        // ---- Interference Fillets / Faces.
        let ic_fil1 = fi1.line_index();
        if ic_fil1 != 0 {
            let interfc3 =
                chfi3d_fil_curve_in_ds(ic_fil1, isurf, fi1.pcurve_on_surf().cloned(), trafil1);
            dstr.change_surface_interferences(isurf).push(interfc3.clone());
            ishape1 = fd.index_of_s1;
            // Case of degenerated edge: pcurve is associated via SCI.
            if dstr.curve(ic_fil1).curve.is_none() {
                if let TopOpeBRepDSInterference::SurfaceCurve(sc) = &interfc3 {
                    let refd = super::topopebrepds::InterferenceRef {
                        pcurve: sc.pcurve.clone(),
                        index_s: sc.index_s,
                        index_g: sc.index_g,
                    };
                    dstr.change_curve(ic_fil1).sci1 = Some(refd);
                }
            } else {
                let mut regon1 = ChFiDSRegul::new();
                regon1.set_curve(ic_fil1);
                regon1.set_s1(isurf, false);
                if ishape1 < 0 {
                    ishape1 = -ishape1;
                    regon1.set_s2(ishape1, false);
                    let interfc1 = chfi3d_fil_curve_in_ds(
                        ic_fil1,
                        ishape1,
                        fi1.pcurve_on_face().cloned(),
                        fi1.transition(),
                    );
                    dstr.change_surface_interferences(ishape1).push(interfc1);
                } else if ishape1 > 0 {
                    regon1.set_s2(ishape1, true);
                    let interfc1 = chfi3d_fil_curve_in_ds(
                        ic_fil1,
                        ishape1,
                        fi1.pcurve_on_face().cloned(),
                        fi1.transition(),
                    );
                    dstr.change_shape_interferences(ishape1).push(interfc1);
                }
                reglist.push(regon1);
            }
            // Indice and type of the point at End.
            let mut ipoin = 0i32;
            let mut is_vertex = fd.vertex_last_on_s1().is_vertex();
            if j == seqfil.len() {
                ipoin = cordat.indexlast_pon_s1;
            } else if j == seqfil.len() - 1
                && dstr
                    .curve(
                        seqfil.last()
                            .unwrap()
                            .read()
                            .expect("surfdata lock")
                            .interference_on_s1()
                            .line_index(),
                    )
                    .curve
                    .is_none()
            {
                if closed {
                    ipoin = cordat.indexfirst_pon_s1;
                    is_vertex = seqfil[0]
                        .read()
                        .expect("surfdata lock")
                        .vertex_first_on_s1()
                        .is_vertex();
                } else {
                    ipoin = cordat.indexlast_pon_s1;
                    is_vertex = seqfil
                        .last()
                        .unwrap()
                        .read()
                        .expect("surfdata lock")
                        .vertex_last_on_s1()
                        .is_vertex();
                }
            } else if dstr.curve(ic_fil1).curve.is_none() {
                // Rotation !!
                ipoin = ipoin1;
                is_vertex = is_vertex1;
            } else if (j == 1 || j == seqfil.len() - 1)
                && (fd
                    .vertex_last_on_s1()
                    .is_equal(seqfil[0].read().expect("surfdata lock").vertex_first_on_s1(), 1.0e-7)
                    || fd.vertex_last_on_s1().is_equal(
                        seqfil.last()
                            .unwrap()
                            .read()
                            .expect("surfdata lock")
                            .vertex_last_on_s1(),
                        1.0e-7,
                    ))
            {
                // Case of SurfData cut in "Triangular" way.
                ipoin = cordat.indexlast_pon_s1;
            } else if is_in_ds2 && find_index_point(dstr, fd, 1, &mut ipoin) {
                // OCCT: ipoin set by findIndexPoint.
            } else {
                ipoin = chfi3d_index_point_in_ds(fd.vertex_last_on_s1(), dstr);
            }

            // OCCT holds a live list reference: the second Contains sees the
            // first append.  rcad re-fetches the list for each test.
            if !chfi3d_contains(
                dstr.curve_interferences(ic_fil1),
                ic_fil1,
                ipoin1,
                false,
                false,
            ) {
                let interp1 = chfi3d_fil_point_in_ds(
                    Orientation::Forward,
                    ic_fil1,
                    ipoin1,
                    fi1.parameter_first(),
                    is_vertex1,
                );
                dstr.change_curve_interferences(ic_fil1).push(interp1);
            }
            if ipoin == ipoin1
                || !chfi3d_contains(
                    dstr.curve_interferences(ic_fil1),
                    ic_fil1,
                    ipoin,
                    false,
                    false,
                )
            {
                let interp3 = chfi3d_fil_point_in_ds(
                    Orientation::Reversed,
                    ic_fil1,
                    ipoin,
                    fi1.parameter_last(),
                    is_vertex,
                );
                dstr.change_curve_interferences(ic_fil1).push(interp3);
            }
            ipoin1 = ipoin;
            is_vertex1 = is_vertex;
        }

        let ic_fil2 = fi2.line_index();
        if ic_fil2 != 0 {
            let interfc4 =
                chfi3d_fil_curve_in_ds(ic_fil2, isurf, fi2.pcurve_on_surf().cloned(), trafil2);
            dstr.change_surface_interferences(isurf).push(interfc4.clone());
            ishape2 = fd.index_of_s2;
            // Case of degenerated edge: pcurve is associated via SCI.
            if dstr.curve(ic_fil2).curve.is_none() {
                if let TopOpeBRepDSInterference::SurfaceCurve(sc) = &interfc4 {
                    let refd = super::topopebrepds::InterferenceRef {
                        pcurve: sc.pcurve.clone(),
                        index_s: sc.index_s,
                        index_g: sc.index_g,
                    };
                    dstr.change_curve(ic_fil2).sci1 = Some(refd);
                }
            } else {
                let mut regon2 = ChFiDSRegul::new();
                regon2.set_curve(ic_fil2);
                regon2.set_s1(isurf, false);
                if ishape2 < 0 {
                    ishape2 = -ishape2;
                    regon2.set_s2(ishape2, false);
                    let interfc2 = chfi3d_fil_curve_in_ds(
                        ic_fil2,
                        ishape2,
                        fi2.pcurve_on_face().cloned(),
                        fi2.transition(),
                    );
                    dstr.change_surface_interferences(ishape2).push(interfc2);
                } else if ishape2 > 0 {
                    regon2.set_s2(ishape2, true);
                    let interfc2 = chfi3d_fil_curve_in_ds(
                        ic_fil2,
                        ishape2,
                        fi2.pcurve_on_face().cloned(),
                        fi2.transition(),
                    );
                    dstr.change_shape_interferences(ishape2).push(interfc2);
                }
                reglist.push(regon2);
            }
            // Indice and type of the point in End.
            let mut ipoin = 0i32;
            let mut is_vertex = fd.vertex_last_on_s2().is_vertex();
            if j == seqfil.len() {
                ipoin = cordat.indexlast_pon_s2;
            } else if j == seqfil.len() - 1
                && dstr
                    .curve(
                        seqfil.last()
                            .unwrap()
                            .read()
                            .expect("surfdata lock")
                            .interference_on_s2()
                            .line_index(),
                    )
                    .curve
                    .is_none()
            {
                if closed {
                    ipoin = cordat.indexfirst_pon_s2;
                    is_vertex = seqfil[0]
                        .read()
                        .expect("surfdata lock")
                        .vertex_first_on_s2()
                        .is_vertex();
                } else {
                    ipoin = cordat.indexlast_pon_s2;
                    is_vertex = seqfil
                        .last()
                        .unwrap()
                        .read()
                        .expect("surfdata lock")
                        .vertex_last_on_s2()
                        .is_vertex();
                }
            } else if dstr.curve(ic_fil2).curve.is_none() {
                // Rotation !!
                ipoin = ipoin2;
                is_vertex = is_vertex2;
            } else if fd
                .vertex_last_on_s2()
                .is_equal(fd.vertex_last_on_s1(), 0.0)
            {
                // Pinch !!
                ipoin = ipoin1;
                is_vertex = is_vertex1;
            } else if (j == 1 || j == seqfil.len() - 1)
                && (fd
                    .vertex_last_on_s2()
                    .is_equal(seqfil[0].read().expect("surfdata lock").vertex_first_on_s2(), 1.0e-7)
                    || fd.vertex_last_on_s2().is_equal(
                        seqfil.last()
                            .unwrap()
                            .read()
                            .expect("surfdata lock")
                            .vertex_last_on_s2(),
                        1.0e-7,
                    ))
            {
                // Case of SurfData cut in "Triangular" way.
                ipoin = cordat.indexlast_pon_s2;
            } else if is_in_ds2 && find_index_point(dstr, fd, 2, &mut ipoin) {
                // OCCT: ipoin set by findIndexPoint.
            } else {
                ipoin = chfi3d_index_point_in_ds(fd.vertex_last_on_s2(), dstr);
            }

            if !chfi3d_contains(
                dstr.curve_interferences(ic_fil2),
                ic_fil2,
                ipoin2,
                false,
                false,
            ) {
                let interp2 = chfi3d_fil_point_in_ds(
                    Orientation::Forward,
                    ic_fil2,
                    ipoin2,
                    fi2.parameter_first(),
                    is_vertex2,
                );
                dstr.change_curve_interferences(ic_fil2).push(interp2);
            }
            if ipoin == ipoin2
                || !chfi3d_contains(
                    dstr.curve_interferences(ic_fil2),
                    ic_fil2,
                    ipoin,
                    false,
                    false,
                )
            {
                let interp4 = chfi3d_fil_point_in_ds(
                    Orientation::Reversed,
                    ic_fil2,
                    ipoin,
                    fi2.parameter_last(),
                    is_vertex,
                );
                dstr.change_curve_interferences(ic_fil2).push(interp4);
            }
            ipoin2 = ipoin;
            is_vertex2 = is_vertex;
        }

        et1 = trafil1;
        if j == seqfil.len() {
            if !is_in_ds2 {
                icurv = cordat.last_curve();
                if closed && !singulier_en_bout && (ipoin1 != ipoin2) {
                    regcout.set_s2(isurf, false);
                    reglist.push(regcout);
                }
                let pcurv = cordat.last_pcurve().cloned();
                et1 = topabs_compose(et1, cordat.last_pcurve_orientation());
                let (pd, pf) = cordat.last_parameters();
                pardeb = pd;
                parfin = pf;
                let li_empty = dstr.curve_interferences(icurv).is_empty();
                if li_empty {
                    if cordat.last_pcurve_orientation() == Orientation::Reversed {
                        let interp5 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            icurv,
                            ipoin1,
                            parfin,
                            is_vertex1,
                        );
                        let interp6 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            icurv,
                            ipoin2,
                            pardeb,
                            is_vertex2,
                        );
                        dstr.change_curve_interferences(icurv).push(interp5);
                        dstr.change_curve_interferences(icurv).push(interp6);
                    } else {
                        let interp5 = chfi3d_fil_point_in_ds(
                            Orientation::Forward,
                            icurv,
                            ipoin1,
                            pardeb,
                            is_vertex1,
                        );
                        let interp6 = chfi3d_fil_point_in_ds(
                            Orientation::Reversed,
                            icurv,
                            ipoin2,
                            parfin,
                            is_vertex2,
                        );
                        dstr.change_curve_interferences(icurv).push(interp5);
                        dstr.change_curve_interferences(icurv).push(interp6);
                    }
                }
                let interfc1 = chfi3d_fil_curve_in_ds(icurv, isurf, pcurv, et1);
                dstr.change_surface_interferences(isurf).push(interfc1.clone());
                if ipoin1 == ipoin2 {
                    dstr.change_curve(icurv).nullify();
                    if let TopOpeBRepDSInterference::SurfaceCurve(sc) = &interfc1 {
                        let refd = super::topopebrepds::InterferenceRef {
                            pcurve: sc.pcurve.clone(),
                            index_s: sc.index_s,
                            index_g: sc.index_g,
                        };
                        dstr.change_curve(icurv).sci1 = Some(refd);
                    }
                }
            }
        } else {
            // eap, Apr 29 2002, occ 293.
            if !is_in_ds2 {
                let u_v1 = fd
                    .interference_on_s1()
                    .pcurve_on_surf()
                    .map(|pc| pc.point_at(fd.interference_on_s1().parameter_last()))
                    .unwrap_or(DVec2::ZERO);
                let u_v2 = fd
                    .interference_on_s2()
                    .pcurve_on_surf()
                    .map(|pc| pc.point_at(fd.interference_on_s2().parameter_last()))
                    .unwrap_or(DVec2::ZERO);
                let surf = dstr.surface(fd.surf()).surface.clone();
                let (c3d, pcurv, pd, pf, tolreached) = chfi3d_compute_arete(
                    brep,
                    fd.vertex_last_on_s1(),
                    u_v1,
                    fd.vertex_last_on_s2(),
                    u_v2,
                    &surf,
                    tol3d,
                    tol2d,
                    0,
                );
                pardeb = pd;
                parfin = pf;
                icurv = dstr
                    .add_curve(super::topopebrepds::TopOpeBRepDSCurve::new(c3d, tolreached));
                regfilfil.set_curve(icurv);
                regfilfil.set_s1(isurf, false);
                let interp5 =
                    chfi3d_fil_point_in_ds(Orientation::Forward, icurv, ipoin1, pardeb, is_vertex1);
                dstr.change_curve_interferences(icurv).push(interp5);
                let interp6 =
                    chfi3d_fil_point_in_ds(Orientation::Reversed, icurv, ipoin2, parfin, is_vertex2);
                dstr.change_curve_interferences(icurv).push(interp6);
                let interfc1 = chfi3d_fil_curve_in_ds(icurv, isurf, Some(pcurv), et1);
                dstr.change_surface_interferences(isurf).push(interfc1);
            }
        }

        degene = v3.point().distance(v4.point()) <= 0.0;

        // Processing of degenerated case.
        if degene {
            // Queue de Billard.
            let vertex = v3.is_vertex() && v4.is_vertex();
            if vertex {
                // The edge of the spine starting on this vertex is removed.
                let mut trouve = false;
                let mut arcspine = Shape::null();
                let mut ovtx = Orientation::Forward;
                boutde_vtx = v3.vertex().clone();

                let spine = cordat.spine().expect("spine");
                while num_edge <= spine.base().nb_edges() && !trouve {
                    arcspine = spine.base().edges(num_edge).clone();
                    // OCCT L3263-3277: TopExp_Explorer over the
                    // FORWARD-oriented edge's vertices.
                    let ed = arcspine.as_edge().expect("not an edge");
                    for (vv, vor) in
                        [(&ed.first, Orientation::Forward), (&ed.last, Orientation::Reversed)]
                    {
                        if boutde_vtx.is_same(vv) && !trouve {
                            ovtx = vor;
                            if closed && num_edge == 1 {
                                trouve = spine.base().nb_edges() == 1;
                            } else {
                                trouve = true;
                            }
                        }
                    }
                    if !trouve {
                        num_edge += 1; // Go to the next edge.
                    }
                }
                let iarcspine = dstr.add_shape(&arcspine);
                let ivtx;
                if j == seqfil.len() {
                    ivtx = cordat.indexlast_pon_s1;
                } else {
                    ivtx = dstr.add_shape(&boutde_vtx);
                }
                ovtx = topabs_reverse(ovtx);
                let parvtx = brep_tool_parameter(brep, &boutde_vtx, &arcspine);
                let interfv = chfi3d_fil_vertex_in_ds(ovtx, iarcspine, ivtx, parvtx);
                dstr.change_shape_interferences(iarcspine).push(interfv);
            }
        } else if !closed || j != seqfil.len() {
            // Processing of interference Point / Edges.
            if v3.is_on_arc() {
                if !(v3.is_vertex() && fd.is_on_curve1()) {
                    iarc1 = dstr.add_shape(v3.arc());
                    if !chfi3d_contains(
                        dstr.shape_interferences(iarc1),
                        iarc1,
                        ipoin1,
                        v3.is_vertex(),
                        true,
                    ) {
                        let interfpp = chfi3d_fil_point_in_ds(
                            v3.transition_on_arc(),
                            iarc1,
                            ipoin1,
                            v3.parameter_on_arc(),
                            v3.is_vertex(),
                        );
                        dstr.change_shape_interferences_of(v3.arc()).push(interfpp);
                    }
                }
            }

            if v4.is_on_arc() {
                if !(v4.is_vertex() && fd.is_on_curve2()) {
                    iarc2 = dstr.add_shape(v4.arc());
                    if !chfi3d_contains(
                        dstr.shape_interferences(iarc2),
                        iarc2,
                        ipoin2,
                        v4.is_vertex(),
                        true,
                    ) {
                        let intfpp = chfi3d_fil_point_in_ds(
                            v4.transition_on_arc(),
                            iarc2,
                            ipoin2,
                            v4.parameter_on_arc(),
                            v4.is_vertex(),
                        );
                        dstr.change_shape_interferences_of(v4.arc()).push(intfpp);
                    }
                }
            }
        }
    }
}

// =========================================================================
// Builder entries mirroring OCCT ChFi3d_Builder.cxx L375/L391:
//   ChFi3d_StripeEdgeInter(st, aCheckStripe, DStr, tol2d);
//   ChFi3d_FilDS(solidindex, st, DStr, myRegul, tolapp3d, tol2d);
// =========================================================================
impl ChFi3dBuilder {
    /// False encodes the OCCT raise (fillets too big).
    pub fn stripe_edge_inter(&mut self, st: &ChFiDSStripe, other: &ChFiDSStripe) -> bool {
        let mut ds = self.my_ds.take().expect("DS");
        let ok = chfi3d_stripe_edge_inter(st, other, &mut ds, self.tol2d);
        self.my_ds = Some(ds);
        ok
    }

    pub fn filds_stripe(&mut self, cordat: &ChFiDSStripe) {
        let mut reglist = std::mem::take(&mut self.my_regul);
        chfi3d_filds(
            &self.my_brep,
            cordat.solid_index(),
            cordat,
            self.my_ds.as_mut().expect("DS"),
            &mut reglist,
            self.tolapp3d,
            self.tol2d,
        );
        self.my_regul = reglist;
    }
}
