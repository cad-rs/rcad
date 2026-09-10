//! OCCT Bisector_Bisec — the trimmed-bisector container, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_Bisec.hxx (L51-115) / .cxx (L49-682).
//!
//! Architecture differences: `handle<Geom2d_Curve>` maps to
//! `Arc<dyn BisectorCurve>`; `Geom2d_Point` maps to `DVec2`; the result
//! `Handle(Geom2d_TrimmedCurve)` maps to [`HandleGeom2dTrimmedCurve`]
//! (Arc<RwLock<TrimmedCurve>>).

use std::sync::{Arc, RwLock};

use glam::DVec2;
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::geom::{Line2d, TrimmedCurve2};

use super::bisector::{cross2, is_convex};
use super::bisector_bisec_ana::BisectorBisecAna;
use super::bisector_bisec_cc::BisectorBisecCC;
use super::bisector_bisec_pc::BisectorBisecPC;
use super::bisector_curve::{BisectorCurve, Geom2dCurveHandle, CurveKind, TrimmedCurve};
use super::bisector::GeomAbsJoinType;

/// OCCT `Handle(Geom2d_TrimmedCurve)` — the shared trimmed result.
pub type HandleGeom2dTrimmedCurve = Arc<RwLock<TrimmedCurve>>;

/// OCCT Bisector_Bisec (Bisector_Bisec.hxx L51-115).
#[derive(Default)]
pub struct BisectorBisec {
    /// OCCT thebisector (null handle until a Perform succeeds).
    thebisector: Option<HandleGeom2dTrimmedCurve>,
}

/// OCCT `GC_MakeSegment2d(P1, P2)` — a trimmed line over [P1,P2]
/// (GC_MakeSegment2d.cxx: parameter range [0, P1.Distance(P2)]).
fn make_segment_2d(p1: DVec2, p2: DVec2) -> Arc<dyn BisectorCurve> {
    let line = Line2d::new(p1, p2 - p1);
    Arc::new(Geom2dCurveHandle::new(rcad_kernel::geom::Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(rcad_kernel::geom::Curve2d::Line(line)),
        t_min: 0.0,
        t_max: p1.distance(p2),
    })))
}

/// OCCT `IsMaxRC(C, U, R)` (L616-681) — tests whether the curvature radius
/// at the U extremity is the maximal one of the curve ends.
fn is_max_rc(c: &Arc<dyn BisectorCurve>, u: f64, r: &mut f64) -> bool {
    let us = c.first_parameter();
    let ul = c.last_parameter();

    let (_p, d1, d2) = c.d2(us);
    let norm2 = d1.length_squared();
    let kf = if norm2 < crate::topalgo::bisector::bisector_curve::GP_RESOLUTION {
        0.0
    } else {
        cross2(d1, d2).abs() / (norm2 * norm2.sqrt())
    };

    let (_p, d1, d2) = c.d2(ul);
    let norm2 = d1.length_squared();
    let kl = if norm2 < crate::topalgo::bisector::bisector_curve::GP_RESOLUTION {
        0.0
    } else {
        cross2(d1, d2).abs() / (norm2 * norm2.sqrt())
    };

    let mut is_max = false;

    if u == ul {
        if kl < kf {
            if kl == 0.0 {
                *r = INFINITE_VALUE;
            } else {
                *r = 1.0 / kl;
            }
            is_max = true;
        }
    } else if kf < kl {
        if kf == 0.0 {
            *r = INFINITE_VALUE;
        } else {
            *r = 1.0 / kf;
        }
        is_max = true;
    }
    is_max
}

/// OCCT `ReplaceByLineIfIsToSmall(Bis, UFirst, ULast)` (L584-612) — if the
/// size of an algorithmic bissectrice is negligible, it is replaced by a
/// half-straight.
fn replace_by_line_if_is_to_small(
    bis: &mut Arc<dyn BisectorCurve>,
    u_first: &mut f64,
    u_last: &mut f64,
) {
    if (*u_last - *u_first).abs() > 2.0 * PCONFUSION * 10.0 {
        return; // patch
    }

    let pf = bis.value(*u_first);
    let pl = bis.value(*u_last);

    if pf.distance(pl) > CONFUSION * 10.0 {
        return;
    }

    let t1 = bis.dn(*u_first, 1);

    let line = Line2d::new(pf, t1);
    let bis_l = super::bisector_curve::TrimmedCurve::new(
        Arc::new(Geom2dCurveHandle::line(&line)),
        0.0,
        INFINITE_VALUE,
    );
    let mut bis_ana = BisectorBisecAna::new();
    bis_ana.init(bis_l);
    *u_first = bis_ana.parameter_of_start_point();
    *u_last = bis_ana.parameter_of_end_point();
    *bis = Arc::new(bis_ana);
}

impl BisectorBisec {
    /// OCCT Bisector_Bisec() (L49).
    pub fn new() -> Self {
        BisectorBisec::default()
    }

    /// OCCT Perform(Cu1, Cu2, P, V1, V2, Sense, ajointype, Tolerance,
    /// oncurve) (L63-246) — the bissectrice between two curves.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_curve(
        &mut self,
        afirstcurve: &Arc<dyn BisectorCurve>,
        asecondcurve: &Arc<dyn BisectorCurve>,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        ajointype: GeomAbsJoinType,
        tolerance: f64,
        oncurve: bool,
    ) {
        // OCCT: Type1/Type2 are the DynamicTypes, looked through a trimmed
        // curve for its BasisCurve type (used for the BSpline / analytic
        // dispatch below).
        let _type1 = curve_kind_through_trim(afirstcurve);
        let _type2 = curve_kind_through_trim(asecondcurve);
        let mut bis: Arc<dyn BisectorCurve>;
        let mut u_first;
        let mut u_last;

        let mut afirstcurve1 = afirstcurve.clone();
        let mut asecondcurve1 = asecondcurve.clone();

        if let Some(a_bs) = bspline_through_trim(afirstcurve) {
            if a_bs.degree == 1 && a_bs.control_points.len() == 2 {
                let p1 = a_bs.control_points[0];
                let p2 = a_bs.control_points[1];
                if p1.distance(p2) < 1.0e-4 {
                    afirstcurve1 = make_segment_2d(p1, p2);
                }
            }
        }

        if let Some(a_bs) = bspline_through_trim(asecondcurve) {
            if a_bs.degree == 1 && a_bs.control_points.len() == 2 {
                let p1 = a_bs.control_points[0];
                let p2 = a_bs.control_points[1];
                if p1.distance(p2) < 1.0e-4 {
                    asecondcurve1 = make_segment_2d(p1, p2);
                }
            }
        }

        let kind1 = curve_kind_through_trim(&afirstcurve1);
        let kind2 = curve_kind_through_trim(&asecondcurve1);

        if (kind1 == CurveKind::Circle || kind1 == CurveKind::Line)
            && (kind2 == CurveKind::Circle || kind2 == CurveKind::Line)
        {
            //------------------------------------------------------------------
            // Analytic Bissectrice.
            //------------------------------------------------------------------
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform(
                &afirstcurve1,
                &asecondcurve1,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                ajointype,
                tolerance,
                oncurve,
            );
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::new(bis_ana);
        } else {
            let mut is_line = false;

            if oncurve {
                let fd = afirstvector.normalize_or_zero();
                let sd = asecondvector.normalize_or_zero();
                // if (Fd.Dot(Sd) < std::sqrt(2. * Precision::Angular()) - 1.)
                if fd.dot(sd) < (2.0 * ANGULAR).sqrt() - 1.0 {
                    is_line = true;
                }
            }
            if is_line {
                //------------------------------------------------------------------
                // Half-Staight.
                //------------------------------------------------------------------
                let n = DVec2::new(
                    -adirection * afirstvector.y,
                    adirection * afirstvector.x,
                );
                let l = Line2d::new(apoint, n);
                let bis_l = super::bisector_curve::TrimmedCurve::new(
                    Arc::new(Geom2dCurveHandle::line(&l)),
                    0.0,
                    INFINITE_VALUE,
                );
                let mut bis_ana = BisectorBisecAna::new();
                bis_ana.init(bis_l);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::new(bis_ana);
            } else {
                //------------------------------------------------------------------
                // Bissectrice algo.
                //------------------------------------------------------------------
                let mut bis_cc = BisectorBisecCC::new();
                bis_cc.perform(&asecondcurve1, &afirstcurve1, adirection, adirection, apoint, 500.0);

                if bis_cc.is_empty() {
                    // bissectrice is empty. a point is projected at the end of
                    // the guide curve. Construction of a false bissectrice.
                    // (modified by NIZHNY-EAP Mon Feb 21 12:00:13 2000)
                    let a_p1 = afirstcurve1.value(afirstcurve1.last_parameter());
                    let a_p2 = asecondcurve1.value(asecondcurve1.first_parameter());
                    let a_pm = DVec2::new(0.5 * (a_p1.x + a_p2.x), 0.5 * (a_p1.y + a_p2.y));
                    let nx;
                    let ny;
                    if a_pm.distance(apoint) > 10.0 * CONFUSION {
                        let mut nx_loc = apoint.x - a_pm.x;
                        let mut ny_loc = apoint.y - a_pm.y;
                        if adirection < 0.0 {
                            nx_loc = -nx_loc;
                            ny_loc = -ny_loc;
                        }
                        nx = nx_loc;
                        ny = ny_loc;
                    } else {
                        let dir1 = afirstvector.normalize_or_zero();
                        let dir2 = asecondvector.normalize_or_zero();
                        let mut nx_loc = -dir1.x - dir2.x;
                        let mut ny_loc = -dir1.y - dir2.y;
                        if nx_loc.abs() <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                            && ny_loc.abs()
                                <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                        {
                            nx_loc = -afirstvector.y;
                            ny_loc = afirstvector.x;
                        }
                        nx = nx_loc;
                        ny = ny_loc;
                    }
                    let n = DVec2::new(adirection * nx, adirection * ny);

                    let l = Line2d::new(apoint, n);
                    let bis_l = super::bisector_curve::TrimmedCurve::new(
                        Arc::new(Geom2dCurveHandle::line(&l)),
                        0.0,
                        INFINITE_VALUE,
                    );
                    let mut bis_ana = BisectorBisecAna::new();
                    bis_ana.init(bis_l);
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::new(bis_ana);
                } else {
                    u_first = bis_cc.first_parameter();
                    u_last = bis_cc.last_parameter();
                    bis = Arc::new(bis_cc);
                    replace_by_line_if_is_to_small(&mut bis, &mut u_first, &mut u_last);
                }
            }
        }
        u_first = u_first.max(bis.first_parameter());
        u_last = u_last.min(bis.last_parameter());
        self.thebisector = Some(Arc::new(RwLock::new(TrimmedCurve::new(bis, u_first, u_last),
        )));
    }

    /// OCCT Perform(Cu, Pnt, P, V1, V2, Sense, Tolerance, oncurve) (L260-391)
    /// — the bissectrice between a curve and a point.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_curve_point(
        &mut self,
        afirstcurve: &Arc<dyn BisectorCurve>,
        asecondpoint: DVec2,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut bis: Arc<dyn BisectorCurve>;
        let mut u_first;
        let mut u_last;

        let type1 = curve_kind_through_trim(afirstcurve);

        if type1 == CurveKind::Circle || type1 == CurveKind::Line {
            //------------------------------------------------------------------
            // Analytic Bissectrice.
            //------------------------------------------------------------------
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform_curve_point(
                afirstcurve,
                asecondpoint,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                tolerance,
                oncurve,
            );
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::new(bis_ana);
        } else {
            let mut is_line = false;
            let mut rc = INFINITE_VALUE;

            if oncurve
                && (is_convex(afirstcurve, adirection)
                    || is_max_rc(afirstcurve, afirstcurve.last_parameter(), &mut rc))
            {
                is_line = true;
            }
            if is_line {
                //------------------------------------------------------------------
                // Half-Right.
                //------------------------------------------------------------------
                let n = DVec2::new(
                    -adirection * afirstvector.y,
                    adirection * afirstvector.x,
                );
                let l = Line2d::new(apoint, n);
                let bis_l = super::bisector_curve::TrimmedCurve::new(
                    Arc::new(Geom2dCurveHandle::line(&l)),
                    0.0,
                    rc,
                );
                let mut bis_ana = BisectorBisecAna::new();
                bis_ana.init(bis_l);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::new(bis_ana);
            } else {
                //------------------------------------------------------------------
                // Bissectrice algo.
                //------------------------------------------------------------------
                let mut bis_pc = BisectorBisecPC::new();
                let afirstcurvereverse = afirstcurve.reversed();

                bis_pc.perform(&afirstcurvereverse, asecondpoint, -adirection, 500.0);
                // Modified by Sergey KHROMOV - Thu Feb 21 16:49:54 2002 Begin
                if bis_pc.is_empty() {
                    let dir1 = afirstvector.normalize_or_zero();
                    let dir2 = asecondvector.normalize_or_zero();
                    let mut nx = -dir1.x - dir2.x;
                    let mut ny = -dir1.y - dir2.y;
                    if nx.abs() <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                        && ny.abs() <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                    {
                        nx = -afirstvector.y;
                        ny = afirstvector.x;
                    }
                    let n = DVec2::new(adirection * nx, adirection * ny);
                    let l = Line2d::new(apoint, n);
                    let bis_l = super::bisector_curve::TrimmedCurve::new(
                        Arc::new(Geom2dCurveHandle::line(&l)),
                        0.0,
                        rc,
                    );
                    let mut bis_ana = BisectorBisecAna::new();
                    bis_ana.init(bis_l);
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::new(bis_ana);
                } else {
                    // Modified by Sergey KHROMOV - Wed Mar  6 17:01:08 2002 End
                    u_first = bis_pc.parameter(apoint);
                    u_last = bis_pc.last_parameter();
                    if u_first >= u_last {
                        // Extrapolate by line.
                        let v = bis_pc.value(bis_pc.first_parameter())
                            - bis_pc.value(u_last);
                        let l = Line2d::new(apoint, v);
                        let bis_l = super::bisector_curve::TrimmedCurve::new(
                            Arc::new(Geom2dCurveHandle::line(&l)),
                            0.0,
                            rc,
                        );
                        let mut bis_ana = BisectorBisecAna::new();
                        bis_ana.init(bis_l);
                        u_first = bis_ana.parameter_of_start_point();
                        u_last = bis_ana.parameter_of_end_point();
                        bis = Arc::new(bis_ana);
                    } else {
                        bis = Arc::new(bis_pc);
                    }
                }
            }
        }
        if u_first < bis.first_parameter() {
            u_first = bis.first_parameter();
        }
        if u_last > bis.last_parameter() {
            u_last = bis.last_parameter();
        }
        self.thebisector = Some(Arc::new(RwLock::new(TrimmedCurve::new(bis, u_first, u_last),
        )));
    }

    /// OCCT Perform(Pnt, Cu, P, V1, V2, Sense, Tolerance, oncurve) (L405-529)
    /// — the bissectrice between a point and a curve.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_curve(
        &mut self,
        afirstpoint: DVec2,
        asecondcurve: &Arc<dyn BisectorCurve>,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut bis: Arc<dyn BisectorCurve>;
        let mut u_first;
        let mut u_last;

        let type1 = curve_kind_through_trim(asecondcurve);

        if type1 == CurveKind::Circle || type1 == CurveKind::Line {
            //------------------------------------------------------------------
            // Analytic Bissectrice.
            //------------------------------------------------------------------
            let mut bis_ana = BisectorBisecAna::new();
            bis_ana.perform_point_curve(
                afirstpoint,
                asecondcurve,
                apoint,
                afirstvector,
                asecondvector,
                adirection,
                tolerance,
                oncurve,
            );
            u_first = bis_ana.parameter_of_start_point();
            u_last = bis_ana.parameter_of_end_point();
            bis = Arc::new(bis_ana);
        } else {
            let mut is_line = false;
            let mut rc = INFINITE_VALUE;

            if oncurve
                && (is_convex(asecondcurve, adirection)
                    || is_max_rc(asecondcurve, asecondcurve.first_parameter(), &mut rc))
            {
                is_line = true;
            }
            if is_line {
                //------------------------------------------------------------------
                // Half-Staight.
                //------------------------------------------------------------------
                let n = DVec2::new(
                    -adirection * afirstvector.y,
                    adirection * afirstvector.x,
                );
                let l = Line2d::new(apoint, n);
                let bis_l = super::bisector_curve::TrimmedCurve::new(
                    Arc::new(Geom2dCurveHandle::line(&l)),
                    0.0,
                    rc,
                );
                let mut bis_ana = BisectorBisecAna::new();
                bis_ana.init(bis_l);
                u_first = bis_ana.parameter_of_start_point();
                u_last = bis_ana.parameter_of_end_point();
                bis = Arc::new(bis_ana);
            } else {
                //------------------------------------------------------------------
                // Bissectrice algo.
                //------------------------------------------------------------------
                let mut bis_pc = BisectorBisecPC::new();
                bis_pc.perform(asecondcurve, afirstpoint, adirection, 500.0);
                // Modified by Sergey KHROMOV - Thu Feb 21 16:49:54 2002 Begin
                if bis_pc.is_empty() {
                    let dir1 = afirstvector.normalize_or_zero();
                    let dir2 = asecondvector.normalize_or_zero();
                    let mut nx = -dir1.x - dir2.x;
                    let mut ny = -dir1.y - dir2.y;
                    if nx.abs() <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                        && ny.abs() <= crate::topalgo::bisector::bisector_curve::GP_RESOLUTION
                    {
                        nx = -afirstvector.y;
                        ny = afirstvector.x;
                    }
                    let n = DVec2::new(adirection * nx, adirection * ny);
                    let l = Line2d::new(apoint, n);
                    let bis_l = super::bisector_curve::TrimmedCurve::new(
                        Arc::new(Geom2dCurveHandle::line(&l)),
                        0.0,
                        rc,
                    );
                    let mut bis_ana = BisectorBisecAna::new();
                    bis_ana.init(bis_l);
                    u_first = bis_ana.parameter_of_start_point();
                    u_last = bis_ana.parameter_of_end_point();
                    bis = Arc::new(bis_ana);
                } else {
                    // Modified by Sergey KHROMOV - Thu Feb 21 16:49:58 2002 End
                    u_first = bis_pc.parameter(apoint);
                    u_last = bis_pc.last_parameter();
                    if u_first >= u_last {
                        // Extrapolate by line.
                        let v = bis_pc.value(bis_pc.first_parameter())
                            - bis_pc.value(u_last);
                        let l = Line2d::new(apoint, v);
                        let bis_l = super::bisector_curve::TrimmedCurve::new(
                            Arc::new(Geom2dCurveHandle::line(&l)),
                            0.0,
                            rc,
                        );
                        let mut bis_ana = BisectorBisecAna::new();
                        bis_ana.init(bis_l);
                        u_first = bis_ana.parameter_of_start_point();
                        u_last = bis_ana.parameter_of_end_point();
                        bis = Arc::new(bis_ana);
                    } else {
                        bis = Arc::new(bis_pc);
                    }
                }
            }
        }

        u_first = u_first.max(bis.first_parameter());
        u_last = u_last.min(bis.last_parameter());
        self.thebisector = Some(Arc::new(RwLock::new(TrimmedCurve::new(bis, u_first, u_last),
        )));
    }

    /// OCCT Perform(Pnt1, Pnt2, P, V1, V2, Sense, Tolerance, oncurve)
    /// (L542-563) — the bissectrice between two points.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_point_point(
        &mut self,
        afirstpoint: DVec2,
        asecondpoint: DVec2,
        apoint: DVec2,
        afirstvector: DVec2,
        asecondvector: DVec2,
        adirection: f64,
        tolerance: f64,
        oncurve: bool,
    ) {
        let mut bis = BisectorBisecAna::new();

        bis.perform_point_point(
            afirstpoint,
            asecondpoint,
            apoint,
            afirstvector,
            asecondvector,
            adirection,
            tolerance,
            oncurve,
        );
        let u_first = bis.parameter_of_start_point();
        let u_last = bis.parameter_of_end_point();
        self.thebisector = Some(Arc::new(RwLock::new(TrimmedCurve::new(
            Arc::new(bis),
            u_first,
            u_last,
        ))));
    }

    /// OCCT Value() const (L567-570) — the Curve of <me>.
    pub fn value(&self) -> HandleGeom2dTrimmedCurve {
        self.thebisector
            .as_ref()
            .expect("Bisector_Bisec::Value on a null handle")
            .clone()
    }

    /// OCCT ChangeValue() (L574-577) — the Curve of <me>.
    pub fn change_value(&self) -> HandleGeom2dTrimmedCurve {
        self.value()
    }
}

/// OCCT `Type1 == STANDARD_TYPE(Geom2d_TrimmedCurve) ? BasisCurve type :
/// DynamicType` — the kind of a handle, looked through a trimmed curve.
fn curve_kind_through_trim(c: &Arc<dyn BisectorCurve>) -> CurveKind {
    match c.as_any().downcast_ref::<super::bisector_curve::TrimmedCurve>() {
        Some(t) => t.basis().kind(),
        None => c.kind(),
    }
}

/// OCCT `Type1 == STANDARD_TYPE(Geom2d_BSplineCurve)` with
/// `aBS = down_cast<Geom2d_BSplineCurve>(trimmed ? BasisCurve : itself)` —
/// the kernel BSpline payload of the handle, looked through a trimmed curve.
fn bspline_through_trim(
    c: &Arc<dyn BisectorCurve>,
) -> Option<&rcad_kernel::geom::BSplineCurve2> {
    let basis: &Arc<dyn BisectorCurve> = match c.as_any().downcast_ref::<super::bisector_curve::TrimmedCurve>() {
        Some(t) => t.basis(),
        None => c,
    };
    let handle = basis.as_any().downcast_ref::<Geom2dCurveHandle>()?;
    match &handle.curve {
        rcad_kernel::geom::Curve2d::BSpline(b) => Some(b),
        _ => None,
    }
}
