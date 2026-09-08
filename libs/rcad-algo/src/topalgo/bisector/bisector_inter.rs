//! OCCT Bisector_Inter — intersection between two Bisec from Bisector,
//! 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_Inter.hxx (L32-89) / .cxx (L44-463).
//!
//! Architecture differences: `class Bisector_Inter : public
//! IntRes2d_Intersection` maps to composition over
//! `geomalgo::int_res2d::IntersectionBase`; the `done` flag is mirrored
//! locally because NeighbourPerform writes the base-class member directly.

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::core::precision::PCONFUSION;
use rcad_kernel::geom::Line2d;

use super::bisector_bisec::BisectorBisec;
use super::bisector_bisec_cc::BisectorBisecCC;
use super::bisector_curve::{
    BisectorCurve, CurveKind, Geom2dCurveAdaptor, Geom2dCurveHandle,
};
use super::bisector_function_inter::FunctionInter;
use super::deps_gap::MathBissecNewton;
use crate::geomalgo::geom2d_int::{elclib2d, GInter};
use crate::geomalgo::int_res2d::{
    Domain, IntersectionBase, IntersectionPoint, Transition,
};

/// OCCT RealFirst()/RealLast().
const REAL_FIRST: f64 = -f64::MAX;
const REAL_LAST: f64 = f64::MAX;

/// OCCT static ConstructSegment(PMin, PMax, UMin, UMax) (L61-71).
fn construct_segment(
    p_min: DVec2,
    p_max: DVec2,
    u_min: f64,
    _u_max: f64,
) -> Arc<dyn BisectorCurve> {
    let dir = DVec2::new(p_max.x - p_min.x, p_max.y - p_min.y);
    let l = Line2d::new(
        DVec2::new(p_min.x - u_min * dir.x, p_min.y - u_min * dir.y),
        dir,
    );
    Arc::new(Geom2dCurveHandle::line(&l))
}

/// OCCT Bisector_Inter (Bisector_Inter.hxx L32-89).
pub struct BisectorInter {
    /// OCCT base class IntRes2d_Intersection sub-object.
    base: IntersectionBase,
    /// OCCT IntRes2d_Intersection::done — mirrored locally (see module
    /// note); NeighbourPerform writes it directly.
    done: bool,
}

impl Default for BisectorInter {
    /// OCCT Bisector_Inter() (L44).
    fn default() -> Self {
        BisectorInter { base: IntersectionBase::new(), done: false }
    }
}

impl BisectorInter {
    /// OCCT Bisector_Inter() (L44).
    pub fn new() -> Self {
        BisectorInter::default()
    }

    /// OCCT Bisector_Inter(C1, D1, C2, D2, TolConf, Tol, ComunElement)
    /// (L48-57).
    #[allow(clippy::too_many_arguments)]
    pub fn with_curves(
        c1: &BisectorBisec,
        d1: &Domain,
        c2: &BisectorBisec,
        d2: &Domain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) -> Self {
        let mut i = BisectorInter::new();
        i.perform(c1, d1, c2, d2, tol_conf, tol, comun_element);
        i
    }

    /// OCCT Perform(C1, D1, C2, D2, TolConf, Tol, ComunElement) (L75-221).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        c1: &BisectorBisec,
        d1: &Domain,
        c2: &BisectorBisec,
        d2: &Domain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) {
        let bis1 = c1.value().read().unwrap().basis().clone();
        let bis2 = c2.value().read().unwrap().basis().clone();

        let size1 = (bis1.nb_intervals() + 1) as usize;
        let size2 = (bis2.nb_intervals() + 1) as usize;
        let mut sbis1: Vec<Option<Arc<dyn BisectorCurve>>> =
            (0..size1).map(|_| None).collect();
        let mut sbis2: Vec<Option<Arc<dyn BisectorCurve>>> =
            (0..size2).map(|_| None).collect();
        let mut sd1: Vec<Domain> = (0..size1).map(|_| Domain::infinite()).collect();
        let mut sd2: Vec<Domain> = (0..size2).map(|_| Domain::infinite()).collect();

        let mut nb1: i32 = 0;
        let mut nb2: i32 = 0;
        let mut min_domain;
        let mut max_domain;
        let mut u_min;
        let mut u_max;
        let mut p_min;
        let mut p_max;

        //------------------------------------------------------
        // Return Min Max domain1.
        //------------------------------------------------------
        if d1.has_first_point() {
            min_domain = d1.first_parameter();
        } else {
            min_domain = REAL_FIRST;
        }

        if d1.has_last_point() {
            max_domain = d1.last_parameter();
        } else {
            max_domain = REAL_LAST;
        }

        //----------------------------------------------------------
        // Cutting the first curve by the intervals of
        // continuity taking account of D1
        //----------------------------------------------------------
        for ib1 in 1..=bis1.nb_intervals() {
            u_min = bis1.interval_first(ib1);
            u_max = bis1.interval_last(ib1);
            if u_max > min_domain && u_min < max_domain {
                u_min = u_min.max(min_domain);
                u_max = u_max.min(max_domain);
                p_min = bis1.value(u_min);
                p_max = bis1.value(u_max);
                sd1[ib1 as usize].set_values_bounded(
                    p_min,
                    u_min,
                    d1.first_tolerance(),
                    p_max,
                    u_max,
                    d1.last_tolerance(),
                );

                if (ib1 == 1 && bis1.is_extend_at_start())
                    || (ib1 == bis1.nb_intervals() && bis1.is_extend_at_end())
                {
                    //--------------------------------------------------------
                    // Part corresponding to an extension is a segment.
                    //--------------------------------------------------------
                    sbis1[ib1 as usize] = Some(construct_segment(p_min, p_max, u_min, u_max));
                } else {
                    sbis1[ib1 as usize] = Some(bis1.clone());
                }
                nb1 += 1;
            }
        }

        //------------------------------------------------------
        // Return Min Max domain2.
        //------------------------------------------------------
        if d2.has_first_point() {
            min_domain = d2.first_parameter();
        } else {
            min_domain = REAL_FIRST;
        }

        if d2.has_last_point() {
            max_domain = d2.last_parameter();
        } else {
            max_domain = REAL_LAST;
        }

        //----------------------------------------------------------
        // Cut the second curve following the intervals of
        // continuity taking account of D2
        //----------------------------------------------------------
        for ib2 in 1..=bis2.nb_intervals() {
            u_min = bis2.interval_first(ib2);
            u_max = bis2.interval_last(ib2);
            if u_max > min_domain && u_min < max_domain {
                u_min = u_min.max(min_domain);
                u_max = u_max.min(max_domain);
                p_min = bis2.value(u_min);
                p_max = bis2.value(u_max);
                sd2[ib2 as usize].set_values_bounded(
                    p_min,
                    u_min,
                    d2.first_tolerance(),
                    p_max,
                    u_max,
                    d2.last_tolerance(),
                );

                // OCCT quirk kept: Bis1->NbIntervals() in the second clause.
                if (ib2 == 1 && bis2.is_extend_at_start())
                    || (ib2 == bis1.nb_intervals() && bis2.is_extend_at_end())
                {
                    //--------------------------------------------------------
                    // Part corresponding to an extension is a segment.
                    //--------------------------------------------------------
                    sbis2[ib2 as usize] = Some(construct_segment(p_min, p_max, u_min, u_max));
                } else {
                    sbis2[ib2 as usize] = Some(bis2.clone());
                }
                nb2 += 1;
            }
        }

        //--------------------------------------------------------------
        // Loop on the intersections of parts of each curve.
        //--------------------------------------------------------------
        for ib1 in 1..=nb1 {
            for ib2 in 1..=nb2 {
                let s1 = sbis1[ib1 as usize].clone().expect("SBis1");
                let s2 = sbis2[ib2 as usize].clone().expect("SBis2");
                self.single_perform(
                    &s1,
                    &sd1[ib1 as usize],
                    &s2,
                    &sd2[ib2 as usize],
                    tol_conf,
                    tol,
                    comun_element,
                );
            }
        }
    }

    /// OCCT SinglePerform(CBis1, D1, CBis2, D2, TolConf, Tol, ComunElement)
    /// (L225-323).
    #[allow(clippy::too_many_arguments)]
    fn single_perform(
        &mut self,
        cbis1: &Arc<dyn BisectorCurve>,
        d1: &Domain,
        cbis2: &Arc<dyn BisectorCurve>,
        d2: &Domain,
        tol_conf: f64,
        tol: f64,
        comun_element: bool,
    ) {
        let bis1 = cbis1;
        let bis2 = cbis2;

        let type1 = bis1.kind();
        let type2 = bis2.kind();

        if type1 == CurveKind::BisecAna || type2 == CurveKind::BisecAna {
            let c2bis1: Arc<dyn BisectorCurve> = if type1 == CurveKind::BisecAna {
                super::bisector_bisec_ana::downcast_to_ana(bis1).geom2d_curve()
            } else {
                bis1.clone()
            };
            let c2bis2: Arc<dyn BisectorCurve> = if type2 == CurveKind::BisecAna {
                super::bisector_bisec_ana::downcast_to_ana(bis2).geom2d_curve()
            } else {
                bis2.clone()
            };
            let type1 = c2bis1.kind();
            let type2 = c2bis2.kind();
            if type1 == CurveKind::Line && type2 != CurveKind::Line {
                self.test_bound(&c2bis1, d1, &c2bis2, d2, tol_conf, false);
            } else if type2 == CurveKind::Line && type1 != CurveKind::Line {
                self.test_bound(&c2bis2, d2, &c2bis1, d1, tol_conf, true);
            }
            let ac2bis1 = Geom2dCurveAdaptor::new(c2bis1);
            let ac2bis2 = Geom2dCurveAdaptor::new(c2bis2);
            let mut intersect = GInter::new_cd_cd(&ac2bis1, d1, &ac2bis2, d2, tol_conf, tol);
            self.append_result(
                &mut intersect,
                d1.first_parameter(),
                d1.last_parameter(),
                d2.first_parameter(),
                d2.last_parameter(),
            );
        } else if type1 == CurveKind::BisecPC || type2 == CurveKind::BisecPC {
            let abis1 = Geom2dCurveAdaptor::new(bis1.clone());
            let abis2 = Geom2dCurveAdaptor::new(bis2.clone());
            let mut intersect = GInter::new_cd_cd(&abis1, d1, &abis2, d2, tol_conf, tol);
            self.append_result(
                &mut intersect,
                d1.first_parameter(),
                d1.last_parameter(),
                d2.first_parameter(),
                d2.last_parameter(),
            );
        } else if comun_element && type1 == CurveKind::BisecCC && type2 == CurveKind::BisecCC {
            let bcc1 = super::bisector_bisec_cc::cast_from(bis1.clone()).expect("down_cast");
            let bcc2 = super::bisector_bisec_cc::cast_from(bis2.clone()).expect("down_cast");
            self.neighbour_perform(&bcc1, d1, &bcc2, d2, tol);
        } else {
            // If we are here one of two bissectrices is a segment.
            // If one of bissectrices is not a segment, it is tested if
            // its extremities are on the straight line.

            if type1 == CurveKind::Line && type2 != CurveKind::Line {
                self.test_bound(bis1, d1, bis2, d2, tol_conf, false);
            } else if type2 == CurveKind::Line && type1 != CurveKind::Line {
                self.test_bound(bis2, d2, bis1, d1, tol_conf, true);
            }
            let abis1 = Geom2dCurveAdaptor::new(bis1.clone());
            let abis2 = Geom2dCurveAdaptor::new(bis2.clone());
            let mut intersect = GInter::new_cd_cd(&abis1, d1, &abis2, d2, tol_conf, tol);
            self.append_result(
                &mut intersect,
                d1.first_parameter(),
                d1.last_parameter(),
                d2.first_parameter(),
                d2.last_parameter(),
            );
        }
    }

    /// OCCT NeighbourPerform(Bis1, D1, Bis2, D2, Tol) (L336-391) — the
    /// intersection of 2 neighbor bissectrices curve/curve.
    fn neighbour_perform(
        &mut self,
        bis1: &Arc<BisectorBisecCC>,
        d1: &Domain,
        bis2: &Arc<BisectorBisecCC>,
        d2: &Domain,
        tol: f64,
    ) {
        let eps = PCONFUSION;

        // Change guideline on Bis2.
        let bis_temp = bis2.change_guide();
        let guide = bis2.curve(2);

        // note: returned points are not used in the code, but can be useful
        // for consulting in debugger.
        let mut u1 = 0.0;
        let mut u_max = 0.0;
        let mut u_min = 0.0;
        let mut dist = 0.0;
        let _p2s = bis2.value_and_dist(d2.first_parameter(), &mut u1, &mut u_max, &mut dist);
        let _p2e = bis2.value_and_dist(d2.last_parameter(), &mut u1, &mut u_min, &mut dist);

        // Calculate the domain of intersection on the guideline.
        let u_min = d1.first_parameter().max(u_min);
        let u_max = d1.last_parameter().min(u_max);

        self.done = true;

        if u_min - eps > u_max + eps {
            return;
        }

        // Solution F = 0 to find the common point.
        let mut f_int = FunctionInter::with_curves(guide, bis1.clone(), bis_temp.clone());

        let mut a_solution = MathBissecNewton::new(tol);
        a_solution.perform(&mut f_int, u_min, u_max, 20);

        if !a_solution.is_done() {
            return;
        }
        let u_sol = a_solution.root();

        let mut u2 = 0.0;
        let p_sol = bis_temp.value_and_dist(u_sol, &mut u1, &mut u2, &mut dist);

        let trans1 = Transition::empty();
        let trans2 = Transition::empty();
        let point_inter_sol = IntersectionPoint::new(
            p_sol,
            u_sol,
            u2,
            trans1,
            trans2,
            false,
        );
        self.base.append_point(&point_inter_sol);
    }

    /// OCCT TestBound(Bis1, D1, Bis2, D2, TolConf, Reverse) (L397-463) —
    /// tests if the extremities of Bis2 are on the segment of Bis1.
    fn test_bound(
        &mut self,
        bis1: &Arc<dyn BisectorCurve>,
        d1: &Domain,
        bis2: &Arc<dyn BisectorCurve>,
        d2: &Domain,
        tol_conf: f64,
        reverse: bool,
    ) {
        let trans1 = Transition::empty();
        let trans2 = Transition::empty();

        // OCCT: gp_Lin2d L1 = Bis1->Lin2d() — the down-cast to Geom2d_Line.
        let handle = bis1
            .as_any()
            .downcast_ref::<Geom2dCurveHandle>()
            .expect("Geom2d_Line down-cast");
        let l1 = super::bisector_bisec_ana::line_of_kernel(&handle.curve);
        let mut pf = bis2.value(d2.first_parameter());
        let mut pl = bis2.value(d2.last_parameter());
        // Modified by skv - Mon May 5 14:43:28 2003 OCC616 Begin
        let tol = tol_conf;
        // Modified by skv - Mon May 5 14:43:30 2003 OCC616 End

        let mut bisec_algo = false;
        if bis2.kind() == CurveKind::BisecCC {
            bisec_algo = true;
        }

        if l1.distance(pf) < tol {
            let u1 = elclib2d::line_parameter(l1.origin, l1.direction, pf);
            // Modified by skv - Mon May 5 14:48:12 2003 OCC616 Begin
            if d1.first_parameter() - d1.first_tolerance() < u1
                && d1.last_parameter() + d1.last_tolerance() > u1
            {
                // Modified by skv - Mon May 5 14:48:14 2003 OCC616 End
                // PF est sur L1
                if bisec_algo {
                    pf = elclib2d::line_value(l1.origin, l1.direction, u1);
                }
                let mut point_inter_sol = IntersectionPoint::empty();
                point_inter_sol.set_values(
                    pf,
                    u1,
                    d2.first_parameter(),
                    trans1.clone(),
                    trans2.clone(),
                    reverse,
                );
                self.base.append_point(&point_inter_sol);
            }
        }

        if l1.distance(pl) < tol {
            let u1 = elclib2d::line_parameter(l1.origin, l1.direction, pl);
            // Modified by skv - Mon May 5 15:05:48 2003 OCC616 Begin
            if d1.first_parameter() - d1.first_tolerance() < u1
                && d1.last_parameter() + d1.last_tolerance() > u1
            {
                // Modified by skv - Mon May 5 15:05:49 2003 OCC616 End
                if bisec_algo {
                    pl = elclib2d::line_value(l1.origin, l1.direction, u1);
                }
                let mut point_inter_sol = IntersectionPoint::empty();
                point_inter_sol.set_values(
                    pl,
                    u1,
                    d2.last_parameter(),
                    trans1.clone(),
                    trans2.clone(),
                    reverse,
                );
                self.base.append_point(&point_inter_sol);
            }
        }
    }

    /// OCCT Append(PointInterSol) — the IntRes2d_Intersection base method.
    pub fn append_point(&mut self, point: &IntersectionPoint) {
        self.base.append_point(point);
    }

    /// OCCT Append(Intersect, FirstParam1, LastParam1, FirstParam2,
    /// LastParam2) — the IntRes2d_Intersection base method.
    pub fn append_result(
        &mut self,
        intersect: &mut GInter,
        first_param1: f64,
        last_param1: f64,
        first_param2: f64,
        last_param2: f64,
    ) {
        self.base.append_intersector(
            &intersect.base,
            first_param1,
            last_param1,
            first_param2,
            last_param2,
        );
    }

    // -------------------------------------------------------------------
    // IntRes2d_Intersection result surface.
    // -------------------------------------------------------------------

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsEmpty().
    pub fn is_empty(&self) -> bool {
        self.base.is_empty()
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// OCCT Point(N) — 1-based.
    pub fn point(&self, n: usize) -> &IntersectionPoint {
        self.base.point(n)
    }

    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// OCCT Segment(N) — 1-based.
    pub fn segment(&self, n: usize) -> &crate::geomalgo::int_res2d::IntersectionSegment {
        self.base.segment(n)
    }
}

