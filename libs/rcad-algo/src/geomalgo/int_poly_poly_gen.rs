//! OCCT IntCurve_IntPolyPolyGen (TKGeomAlgo IntCurve package) — polyline /
//! polyline intersection of two parametrised curves (and the
//! auto-intersection of one), refined to exact points.
//!
//! 1:1 translation of `IntCurve_IntPolyPolyGen.gxx` (L1-1797) over the
//! template parameters `TheCurve`/`TheCurveTool`/`TheProjPCur`, mapped to
//! the curve type `C` plus the tool trait bundle `T`
//! ([`ProjPCurveTool`] + [`ParTool`] + [`ProjectOnPCurveTool`]); the
//! polygon side keeps the concrete IntCurve_ThePolygon2d ==
//! [`Polygon2dGen`]:
//! - Geom2dInt_TheIntPCurvePCurveOfGInter — TheCurve = Adaptor2d_Curve2d,
//!   TheCurveTool = Geom2dInt_Geom2dCurveTool, TheProjPCur =
//!   Geom2dInt_TheProjPCurOfGInter;
//! - HLRBRep_TheIntPCurvePCurveOfCInter — TheCurve = HLRBRep_Curve (Stage 3b
//!   supplies the tool impls over the same traits).
//!
//! The static helpers `HeadOrEndPoint` (gxx L872-1031) and `GetIntersection`
//! (gxx L1575-1783) are free functions here.

use glam::DVec2;
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::direct_polynomial_roots::epsilon;

use super::int_curve_generics::{ExactIntersectionPoint, Polygon2dGen};
use super::int_imp_par_gen::{determine_transition, determine_transition_in_out, ParTool, ProjectOnPCurveTool};
use super::int_curve_generics::ProjPCurveTool;
use super::intf_interference_polygon2d::InterferencePolygon2d;
use super::int_res2d::{
    Domain as Res2dDomain, IntersectionBase, IntersectionPoint, IntersectionSegment, Position,
    Transition,
};

/// OCCT gxx L50: NBITER_MAX_POLYGON.
const NBITER_MAX_POLYGON: i32 = 10;
/// OCCT gxx L51: TOL_CONF_MINI.
const TOL_CONF_MINI: f64 = 0.0000000001;
/// OCCT gxx L52: TOL_MINI.
const TOL_MINI: f64 = 0.0000000001;
/// OCCT Precision::PConfusion() (Precision.hxx).
const PRECISION_P_CONFUSION: f64 = 1.0e-9;

/// OCCT IntCurve_IntPolyPolyGen — the polyline-based pcurve x pcurve
/// intersection engine, generic over the curve type and the tool triple.
pub struct IntPolyPolyGen<C: ?Sized, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    pub base: IntersectionBase,
    /// OCCT DomainOnCurve1.
    pub domain_on_curve1: Res2dDomain,
    /// OCCT DomainOnCurve2.
    pub domain_on_curve2: Res2dDomain,
    /// OCCT myMinPntNb.
    my_min_pnt_nb: i32,
    _tool: std::marker::PhantomData<fn(&T)>,
}

impl<C: ?Sized, T> Clone for IntPolyPolyGen<C, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    fn clone(&self) -> Self {
        IntPolyPolyGen {
            base: self.base.clone(),
            domain_on_curve1: self.domain_on_curve1.clone(),
            domain_on_curve2: self.domain_on_curve2.clone(),
            my_min_pnt_nb: self.my_min_pnt_nb,
            _tool: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, T> std::fmt::Debug for IntPolyPolyGen<C, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IntPolyPolyGen").field("base", &self.base).finish()
    }
}

impl<C: ?Sized, T> IntPolyPolyGen<C, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    /// OCCT IntCurve_IntPolyPolyGen() (gxx L84-89).
    pub fn new() -> Self {
        // const int aMinPntNb = 20; // Minimum number of samples.
        let a_min_pnt_nb = 20;
        IntPolyPolyGen {
            base: IntersectionBase::new(),
            domain_on_curve1: Res2dDomain::infinite(),
            domain_on_curve2: Res2dDomain::infinite(),
            my_min_pnt_nb: a_min_pnt_nb,
            _tool: std::marker::PhantomData,
        }
    }

    /// OCCT GetMinNbSamples() (gxx L1787-1790).
    pub fn get_min_nb_samples(&self) -> i32 {
        self.my_min_pnt_nb
    }

    /// OCCT SetMinNbSamples(theMinNbSamples) (gxx L1794-1797).
    pub fn set_min_nb_samples(&mut self, the_min_nb_samples: i32) {
        self.my_min_pnt_nb = the_min_nb_samples;
    }

    /// OCCT Perform(C1, D1, C2, D2, TheTolConf, TheTol) (gxx L93-292) — the
    /// public pair Perform with the end-point processing.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        c1: &C,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        the_tol_conf: f64,
        the_tol: f64,
    ) {
        self.base.reset_fields();
        self.domain_on_curve1 = d1.clone();
        self.domain_on_curve2 = d2.clone();
        let du = d1.last_parameter() - d1.first_parameter();
        let dv = d2.last_parameter() - d2.first_parameter();
        let tl = if the_tol < TOL_MINI { TOL_MINI } else { the_tol };
        let tl_conf = if the_tol_conf < TOL_CONF_MINI {
            TOL_CONF_MINI
        } else {
            the_tol_conf
        };
        self.perform_pair_internal(c1, d1, c2, d2, tl_conf, tl, 0, du, dv);
        //----------------------------------------------------------------------
        //-- Processing of end points
        //----------------------------------------------------------------------
        let mut head_on1 = false;
        let mut head_on2 = false;
        let mut end_on1 = false;
        let mut end_on2 = false;

        //--------------------------------------------------------------------
        //-- The points Head Head ... End End are not rejected if
        //-- they are already present at the end of segment
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        let mut pos_segment = 0i32;

        let n = self.base.nb_points();
        for i in 1..=n {
            let pos1 = self.base.point(i).transition_of_first().position_on_curve();
            if pos1 == Position::Head {
                head_on1 = true;
            } else if pos1 == Position::End {
                end_on1 = true;
            }

            let pos2 = self.base.point(i).transition_of_second().position_on_curve();
            if pos2 == Position::Head {
                head_on2 = true;
            } else if pos2 == Position::End {
                end_on2 = true;
            }

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }
        }

        let n = self.base.nb_segments();
        for i in 1..=n {
            let seg = self.base.segment(i).clone();
            let pos1 = seg.first_point().transition_of_first().position_on_curve();
            if pos1 == Position::Head {
                head_on1 = true;
            } else if pos1 == Position::End {
                end_on1 = true;
            }

            let pos2 = seg.first_point().transition_of_second().position_on_curve();
            if pos2 == Position::Head {
                head_on2 = true;
            } else if pos2 == Position::End {
                end_on2 = true;
            }

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }

            let pos1 = seg.last_point().transition_of_first().position_on_curve();
            if pos1 == Position::Head {
                head_on1 = true;
            } else if pos1 == Position::End {
                end_on1 = true;
            }

            let pos2 = seg.last_point().transition_of_second().position_on_curve();
            if pos2 == Position::Head {
                head_on2 = true;
            } else if pos2 == Position::End {
                end_on2 = true;
            }

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }
        }

        let u0 = d1.first_parameter();
        let u1 = d1.last_parameter();
        let v0 = d2.first_parameter();
        let v1 = d2.last_parameter();
        let mut int_pt = IntersectionPoint::empty();

        if d1.first_tolerance() != 0.0 || d2.first_tolerance() != 0.0 {
            if head_or_end_point::<C, T>(
                d1, c1, u0, d2, c2, v0, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.first_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point::<C, T>(
                d1, c1, u0, d2, c2, v1, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.first_tolerance() != 0.0 {
            if head_or_end_point::<C, T>(
                d1, c1, u1, d2, c2, v0, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
        if d1.last_tolerance() != 0.0 || d2.last_tolerance() != 0.0 {
            if head_or_end_point::<C, T>(
                d1, c1, u1, d2, c2, v1, the_tol_conf, &mut int_pt, &mut head_on1, &mut head_on2,
                &mut end_on1, &mut end_on2, pos_segment,
            ) {
                self.base.insert(&int_pt);
            }
        }
    }

    /// OCCT Perform(C1, D1, TheTolConf, TheTol) (gxx L297-386) — the public
    /// auto-intersection Perform.
    pub fn perform_auto(&mut self, c1: &C, d1: &Res2dDomain, the_tol_conf: f64, the_tol: f64) {
        self.base.reset_fields();
        self.domain_on_curve1 = d1.clone();
        self.domain_on_curve2 = d1.clone();
        let du = d1.last_parameter() - d1.first_parameter();
        let tl = if the_tol < TOL_MINI { TOL_MINI } else { the_tol };
        let tl_conf = if the_tol_conf < TOL_CONF_MINI {
            TOL_CONF_MINI
        } else {
            the_tol_conf
        };
        self.perform_auto_internal(c1, d1, tl_conf, tl, 0, du, du);
        //--------------------------------------------------------------------
        //-- The points Head Head ... End End are not rejected if
        //-- they are already present at the end of segment
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        let mut pos_segment = 0i32;

        let n = self.base.nb_points();
        for i in 1..=n {
            let pos1 = self.base.point(i).transition_of_first().position_on_curve();
            let pos2 = self.base.point(i).transition_of_second().position_on_curve();

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }
        }

        let n = self.base.nb_segments();
        for i in 1..=n {
            let seg = self.base.segment(i).clone();
            let pos1 = seg.first_point().transition_of_first().position_on_curve();
            let pos2 = seg.first_point().transition_of_second().position_on_curve();

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }

            let pos1 = seg.last_point().transition_of_first().position_on_curve();
            let pos2 = seg.last_point().transition_of_second().position_on_curve();

            if pos1 == Position::Head {
                if pos2 == Position::Head {
                    pos_segment |= 1;
                } else if pos2 == Position::End {
                    pos_segment |= 2;
                }
            } else if pos1 == Position::End {
                if pos2 == Position::Head {
                    pos_segment |= 4;
                } else if pos2 == Position::End {
                    pos_segment |= 8;
                }
            }
        }
        let _ = pos_segment;
    }

    /// OCCT Perform(C1, D1, TolConf, Tol, NbIter, DeltaU, DeltaV)
    /// (gxx L390-870) — the auto-intersection engine.
    #[allow(clippy::too_many_arguments)]
    fn perform_auto_internal(
        &mut self,
        c1: &C,
        d1: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        nb_iter: i32,
        _delta_u: f64,
        _delta_v: f64,
    ) {
        let mut tan1;
        let mut tan2;
        let mut norm1 = DVec2::ZERO;
        let mut norm2 = DVec2::ZERO;
        let mut p1;
        let mut p2;
        self.base.done = false;

        let mut nbsamples = <T as ParTool<C>>::nb_samples_uv(c1, d1.first_parameter(), d1.last_parameter());

        if nb_iter > 3 || (nb_iter > 2 && nbsamples > 100) {
            return;
        }

        //---  We take systematically two times more points
        //--   than on a normal curve.
        //--   Auto-intersecting curves often produce
        //--   polygons rather far from the curve with parameter ct.
        nbsamples *= 2;

        if nb_iter > 0 {
            nbsamples = (3 * (nbsamples * nb_iter)) / 2;
        }
        let poly1 = Polygon2dGen::<C>::new::<T>(
            c1,
            nbsamples,
            d1.first_parameter(),
            d1.last_parameter(),
            tol,
        );
        if !poly1.auto_intersection_is_possible() {
            self.base.done = true;
            return;
        }
        //----------------------------------------------------------------------
        //-- If the deflection is less than the Tolerance of Confusion
        //-- then the deflection of the polygon is set in TolConf
        //-- (Detection of Tangency Zones)
        //----------------------------------------------------------------------
        let mut poly1 = poly1;
        if poly1.deflection_over_estimation() < tol_conf {
            poly1.set_deflection_over_estimation(tol_conf);
        }

        let inter_pp = InterferencePolygon2d::new_polygon_auto(&poly1);
        let mut eip = ExactIntersectionPoint::<C, T>::new(c1, c1, tol_conf);

        //----------------------------------------------------------------------
        //-- Processing of SectionPoint
        //----------------------------------------------------------------------
        let nbsp = inter_pp.interf.nb_section_points();
        if nbsp >= 1 {
            //-- -------------------------------------------------------------------
            //-- filtering, filtering, filtering ...
            //--
            let mut tri_index: Vec<i32> = vec![0; nbsp + 1];
            let mut ptr_seg_index1: Vec<i32> = vec![0; nbsp + 1];
            let mut ptr_seg_index2: Vec<i32> = vec![0; nbsp + 1];
            for i in 1..=nbsp {
                tri_index[i] = i as i32;
                let spnt1 = inter_pp.interf.pnt_value(i);
                let (_typ, seg_index1, _param_on1) = spnt1.info_first();
                ptr_seg_index1[i] = seg_index1;
                let (_typ, seg_index2, _param_on2) = spnt1.info_second();
                ptr_seg_index2[i] = seg_index2;
            }

            loop {
                let mut triok = true;

                for tr in 1..nbsp {
                    let seg_index1 = ptr_seg_index1[tri_index[tr] as usize];
                    let seg_index_1 = ptr_seg_index1[tri_index[tr + 1] as usize];

                    let seg_index2 = ptr_seg_index2[tri_index[tr] as usize];
                    let seg_index_2 = ptr_seg_index2[tri_index[tr + 1] as usize];

                    if seg_index1 > seg_index_1 {
                        let q = tri_index[tr];
                        tri_index[tr] = tri_index[tr + 1];
                        tri_index[tr + 1] = q;
                        triok = false;
                    } else if seg_index1 == seg_index_1 && seg_index2 > seg_index_2 {
                        let q = tri_index[tr];
                        tri_index[tr] = tri_index[tr + 1];
                        tri_index[tr + 1] = q;
                        triok = false;
                    }
                }
                if triok {
                    break;
                }
            }

            //-- supression des doublons Si Si !
            for i in 1..nbsp {
                if ptr_seg_index1[tri_index[i] as usize] == ptr_seg_index1[tri_index[i + 1] as usize]
                    && ptr_seg_index2[tri_index[i] as usize]
                        == ptr_seg_index2[tri_index[i + 1] as usize]
                {
                    tri_index[i] = -(i as i32);
                }
            }

            for sp in 1..=nbsp {
                if tri_index[sp] > 0 {
                    let spnt = inter_pp.interf.pnt_value(tri_index[sp] as usize);
                    let (_typ, seg_index1, param_on1) = spnt.info_first();
                    let (_typ, seg_index2, param_on2) = spnt.info_second();

                    if (seg_index1 - seg_index2).abs() > 1 {
                        let mut num_seg_on1 = seg_index1;
                        let mut num_seg_on2 = seg_index2;
                        let mut param_on_seg1 = param_on1;
                        let mut param_on_seg2 = param_on2;
                        let mut fct = super::int_curve_generics::DistBetweenPCurvesGen::<C, T>::new(c1, c1);
                        eip.perform(
                            &poly1,
                            &poly1,
                            &mut num_seg_on1,
                            &mut num_seg_on2,
                            &mut param_on_seg1,
                            &mut param_on_seg2,
                            &mut fct,
                        );
                        if eip.nb_roots() >= 1 {
                            //----------------------------------------------------------------
                            //-- It is checked if the found point is a root
                            //----------------------------------------------------------------
                            let (u, v) = eip.roots();

                            let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, u);
                            p1 = p1_r;
                            tan1 = tan1_r;
                            let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c1, v);
                            p2 = p2_r;
                            tan2 = tan2_r;
                            let mut dist = p1.distance(p2);
                            let eps_x1 = 10.0 * <T as ParTool<C>>::eps_x(c1);

                            if (u - v).abs() <= eps_x1 {
                                //---------------------------------------------
                                //-- Solution not valid
                                //-- The maths should have converged in a
                                //-- trivial solution  ( point U = V )
                                //---------------------------------------------
                                dist = tol_conf + 1.0;
                            }

                            //---------------------------------------------------------
                            //-- It is checked if the point (u,v) already exists
                            //--
                            self.base.done = true;
                            let nbp = self.base.nb_points();

                            for p in 1..=nbp {
                                let pt = self.base.point(p);
                                if (u - pt.param_on_first()).abs() <= eps_x1
                                    && (v - pt.param_on_second()).abs() <= eps_x1
                                {
                                    dist = tol_conf + 1.0;
                                    break;
                                }
                            }

                            if dist <= tol_conf {
                                //-- Or the point is already present
                                let mut pos1 = Position::Middle;
                                let mut pos2 = Position::Middle;
                                let mut trans1 = Transition::empty();
                                let mut trans2 = Transition::empty();
                                //---------------------------------------------------------
                                //-- Calculate Positions of Points on the curve
                                //--
                                if p1.distance(self.domain_on_curve1.first_point())
                                    <= self.domain_on_curve1.first_tolerance()
                                {
                                    pos1 = Position::Head;
                                } else if p1.distance(self.domain_on_curve1.last_point())
                                    <= self.domain_on_curve1.last_tolerance()
                                {
                                    pos1 = Position::End;
                                }

                                if p2.distance(self.domain_on_curve2.first_point())
                                    <= self.domain_on_curve2.first_tolerance()
                                {
                                    pos2 = Position::Head;
                                } else if p2.distance(self.domain_on_curve2.last_point())
                                    <= self.domain_on_curve2.last_tolerance()
                                {
                                    pos2 = Position::End;
                                }
                                //---------------------------------------------------------
                                if !determine_transition_in_out(
                                    pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2,
                                    tol_conf,
                                ) {
                                    let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, u);
                                    p1 = p1_r;
                                    tan1 = tan1_r;
                                    norm1 = norm1_r;
                                    let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c1, v);
                                    p2 = p2_r;
                                    tan2 = tan2_r;
                                    norm2 = norm2_r;
                                    determine_transition(
                                        pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2,
                                        norm2, &mut trans2, tol_conf,
                                    );
                                }
                                let ip = IntersectionPoint::new(p1, u, v, trans1, trans2, false);
                                self.base.insert(&ip);
                            }
                        }
                    }
                }
            }
        }

        //----------------------------------------------------------------------
        //-- Processing of TangentZone
        //----------------------------------------------------------------------
        let nbtz = inter_pp.interf.nb_tangent_zones();
        for tz in 1..=nbtz {
            let nb_pnts = inter_pp.interf.zone_value(tz).number_of_points();
            //==================================================================
            //== Find the first and the last point in the tangency zone.
            //==================================================================
            let mut param_sup_on_curve1 = -f64::MAX;
            let mut param_sup_on_curve2 = -f64::MAX;
            let mut param_inf_on_curve1 = f64::MAX;
            let mut param_inf_on_curve2 = f64::MAX;
            for qq in 1..=nb_pnts {
                let spnt1 = inter_pp.interf.zone_value(tz).get_point(qq);
                //================================================================
                //== The zones of tangency are discretized
                //== Test of stop : Check if
                //==     (Deflection  < Tolerance)
                //==  Or (Sample < EpsX)   (normally the first condition is
                //==                           more strict)
                //================================================================
                let (_typ, mut seg_index1on_p1, mut param_on_line) = spnt1.info_first();
                if seg_index1on_p1 > poly1.nb_segments() {
                    seg_index1on_p1 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p1 <= 0 {
                    seg_index1on_p1 = 1;
                    param_on_line = 0.0;
                }
                let poly_u_inf_zone = poly1.approx_param_on_curve(seg_index1on_p1, param_on_line);

                let (_typ, mut seg_index1on_p2, mut param_on_line2) = spnt1.info_second();
                if seg_index1on_p2 > poly1.nb_segments() {
                    seg_index1on_p2 -= 1;
                    param_on_line2 = 1.0;
                }
                if seg_index1on_p2 <= 0 {
                    seg_index1on_p2 = 1;
                    param_on_line2 = 0.0;
                }
                let poly_v_inf_zone = poly1.approx_param_on_curve(seg_index1on_p2, param_on_line2);

                //----------------------------------------------------------------

                if param_inf_on_curve1 > poly_u_inf_zone {
                    param_inf_on_curve1 = poly_u_inf_zone;
                }
                if param_inf_on_curve2 > poly_v_inf_zone {
                    param_inf_on_curve2 = poly_v_inf_zone;
                }

                if param_sup_on_curve1 < poly_u_inf_zone {
                    param_sup_on_curve1 = poly_u_inf_zone;
                }
                if param_sup_on_curve2 < poly_v_inf_zone {
                    param_sup_on_curve2 = poly_v_inf_zone;
                }
            }

            let mut poly_u_inf = param_inf_on_curve1;
            let mut poly_u_sup = param_sup_on_curve1;
            let mut poly_v_inf = param_inf_on_curve2;
            let mut poly_v_sup = param_sup_on_curve2;

            p1 = <T as ParTool<C>>::value(c1, poly_u_inf);
            p2 = <T as ParTool<C>>::value(c1, poly_v_inf);
            let distmemesens = p1.distance_squared(p2);
            p2 = <T as ParTool<C>>::value(c1, poly_v_sup);
            let distdiffsens = p1.distance_squared(p2);
            if distmemesens > distdiffsens {
                let qwerty = poly_v_inf;
                poly_v_inf = poly_v_sup;
                poly_v_sup = qwerty;
            }

            //-----------------------------------------------------------------
            //-- Calculate Positions of Points on the curve and
            //-- Transitions on each limit of the segment

            let mut pos1 = Position::Middle;
            let mut pos2 = Position::Middle;
            let mut trans1 = Transition::empty();
            let mut trans2 = Transition::empty();

            let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, poly_u_inf);
            p1 = p1_r;
            tan1 = tan1_r;
            let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c1, poly_v_inf);
            p2 = p2_r;
            tan2 = tan2_r;

            if p1.distance(self.domain_on_curve1.first_point())
                <= self.domain_on_curve1.first_tolerance()
            {
                pos1 = Position::Head;
            } else if p1.distance(self.domain_on_curve1.last_point())
                <= self.domain_on_curve1.last_tolerance()
            {
                pos1 = Position::End;
            }
            if p2.distance(self.domain_on_curve2.first_point())
                <= self.domain_on_curve2.first_tolerance()
            {
                pos2 = Position::Head;
            } else if p2.distance(self.domain_on_curve2.last_point())
                <= self.domain_on_curve2.last_tolerance()
            {
                pos2 = Position::End;
            }

            if pos1 == Position::Middle && pos2 != Position::Middle {
                poly_u_inf = T::find_parameter_between(
                    c1,
                    p2,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    <T as ParTool<C>>::eps_x(c1),
                );
            } else if pos1 != Position::Middle && pos2 == Position::Middle {
                poly_v_inf = T::find_parameter_between(
                    c1,
                    p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    <T as ParTool<C>>::eps_x(c1),
                );
            } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                > (param_inf_on_curve2 - param_sup_on_curve2).abs()
            {
                poly_v_inf = T::find_parameter_between(
                    c1,
                    p1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    <T as ParTool<C>>::eps_x(c1),
                );
            } else {
                poly_u_inf = T::find_parameter_between(
                    c1,
                    p2,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    <T as ParTool<C>>::eps_x(c1),
                );
            }

            if !determine_transition_in_out(
                pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2, tol_conf,
            ) {
                let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, poly_u_inf);
                p1 = p1_r;
                tan1 = tan1_r;
                norm1 = norm1_r;
                let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c1, poly_v_inf);
                p2 = p2_r;
                tan2 = tan2_r;
                norm2 = norm2_r;
                determine_transition(
                    pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2, norm2, &mut trans2,
                    tol_conf,
                );
            }
            let pt_seg1 =
                IntersectionPoint::new(p1, poly_u_inf, poly_v_inf, trans1.clone(), trans2.clone(), false);
            //----------------------------------------------------------------------

            if (poly_u_inf - poly_u_sup).abs() <= <T as ParTool<C>>::eps_x(c1)
                || (poly_v_inf - poly_v_sup).abs() <= <T as ParTool<C>>::eps_x(c1)
            {
                // bad segment
            } else {
                let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, poly_u_sup);
                p1 = p1_r;
                tan1 = tan1_r;
                let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c1, poly_v_sup);
                p2 = p2_r;
                tan2 = tan2_r;
                pos1 = Position::Middle;
                pos2 = Position::Middle;

                if p1.distance(self.domain_on_curve1.first_point())
                    <= self.domain_on_curve1.first_tolerance()
                {
                    pos1 = Position::Head;
                } else if p1.distance(self.domain_on_curve1.last_point())
                    <= self.domain_on_curve1.last_tolerance()
                {
                    pos1 = Position::End;
                }
                if p2.distance(self.domain_on_curve2.first_point())
                    <= self.domain_on_curve2.first_tolerance()
                {
                    pos2 = Position::Head;
                } else if p2.distance(self.domain_on_curve2.last_point())
                    <= self.domain_on_curve2.last_tolerance()
                {
                    pos2 = Position::End;
                }

                if pos1 == Position::Middle && pos2 != Position::Middle {
                    poly_u_sup = T::find_parameter_between(
                        c1,
                        p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                } else if pos1 != Position::Middle && pos2 == Position::Middle {
                    poly_v_sup = T::find_parameter_between(
                        c1,
                        p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_sup = T::find_parameter_between(
                        c1,
                        p1,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                } else {
                    poly_u_sup = T::find_parameter_between(
                        c1,
                        p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                }

                if !determine_transition_in_out(
                    pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2, tol_conf,
                ) {
                    let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, poly_u_sup);
                    p1 = p1_r;
                    tan1 = tan1_r;
                    norm1 = norm1_r;
                    let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c1, poly_v_sup);
                    p2 = p2_r;
                    tan2 = tan2_r;
                    norm2 = norm2_r;
                    determine_transition(
                        pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2, norm2, &mut trans2,
                        tol_conf,
                    );
                }
                let pt_seg2 = IntersectionPoint::new(
                    p1,
                    poly_u_sup,
                    poly_v_sup,
                    trans1.clone(),
                    trans2.clone(),
                    false,
                );

                let oppos = !(tan1.dot(tan2) > 0.0);
                if param_inf_on_curve1 > param_sup_on_curve1 {
                    let seg = IntersectionSegment::with_points(&pt_seg2, &pt_seg1, oppos, false);
                    self.base.append_segment(&seg);
                } else {
                    let seg = IntersectionSegment::with_points(&pt_seg1, &pt_seg2, oppos, false);
                    self.base.append_segment(&seg);
                }
            }
        } // end of processing of TangentZone

        self.base.done = true;
    }
}

impl<C: ?Sized, T> IntPolyPolyGen<C, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    /// OCCT Perform(C1, D1, C2, D2, TolConf, Tol, NbIter, DeltaU, DeltaV)
    /// (gxx L1038-1146) — base method to perform polyline / polyline
    /// intersection for a pair of curves.
    #[allow(clippy::too_many_arguments)]
    fn perform_pair_internal(
        &mut self,
        c1: &C,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        nb_iter: i32,
        delta_u: f64,
        delta_v: f64,
    ) {
        self.base.done = false;

        if nb_iter > NBITER_MAX_POLYGON {
            return;
        }

        // Number of samples running.
        let mut nbsamples_on_c1 = <T as ParTool<C>>::nb_samples_uv(c1, d1.first_parameter(), d1.last_parameter());
        let mut nbsamples_on_c2 = <T as ParTool<C>>::nb_samples_uv(c2, d2.first_parameter(), d2.last_parameter());

        if nb_iter == 0 {
            // Minimal number of points.
            nbsamples_on_c1 = nbsamples_on_c1.max(self.my_min_pnt_nb);
            nbsamples_on_c2 = nbsamples_on_c2.max(self.my_min_pnt_nb);
        } else {
            // Increase number of samples in second and next iterations.
            nbsamples_on_c1 = (5 * (nbsamples_on_c1 * nb_iter)) / 4;
            nbsamples_on_c2 = (5 * (nbsamples_on_c2 * nb_iter)) / 4;
        }

        let mut a_poly1 = Polygon2dGen::<C>::new::<T>(
            c1,
            nbsamples_on_c1,
            d1.first_parameter(),
            d1.last_parameter(),
            tol,
        );
        let mut a_poly2 = Polygon2dGen::<C>::new::<T>(
            c2,
            nbsamples_on_c2,
            d2.first_parameter(),
            d2.last_parameter(),
            tol,
        );

        if (a_poly1.deflection_over_estimation() > tol_conf)
            && (a_poly2.deflection_over_estimation() > tol_conf)
        {
            let a_deflection_sum = a_poly1.deflection_over_estimation().max(tol_conf)
                + a_poly2.deflection_over_estimation().max(tol_conf);

            if nbsamples_on_c2 > nbsamples_on_c1 {
                a_poly2.compute_with_box::<T>(c2, a_poly1.bounding());
                a_poly1.set_deflection_over_estimation(a_deflection_sum);
                a_poly1.compute_with_box::<T>(c1, a_poly2.bounding());
            } else {
                a_poly1.compute_with_box::<T>(c1, a_poly2.bounding());
                a_poly2.set_deflection_over_estimation(a_deflection_sum);
                a_poly2.compute_with_box::<T>(c2, a_poly1.bounding());
            }
        }

        //----------------------------------------------------------------------
        //-- if the deflection less then the Tolerance of Confusion
        //-- Then the deflection of the polygon is set in TolConf
        //-- (Detection of Tangency Zones)
        //----------------------------------------------------------------------

        if a_poly1.deflection_over_estimation() < tol_conf {
            a_poly1.set_deflection_over_estimation(tol_conf);
        }
        if a_poly2.deflection_over_estimation() < tol_conf {
            a_poly2.set_deflection_over_estimation(tol_conf);
        }
        // for case when a few polygon points were replaced by line
        // if exact solution was not found
        // then search of precise solution will be repeated
        // for polygon contains all initial points
        // secondary search will be performed only for case when initial points
        // were dropped
        let is_full_representation = a_poly1.nb_segments() == nbsamples_on_c1
            && a_poly2.nb_segments() == nbsamples_on_c2;

        if !self.find_intersect(
            c1, d1, c2, d2, tol_conf, tol, nb_iter, delta_u, delta_v, &a_poly1, &a_poly2,
            is_full_representation,
        ) && !is_full_representation
        {
            if a_poly1.nb_segments() < nbsamples_on_c1 {
                a_poly1 = Polygon2dGen::<C>::new::<T>(
                    c1,
                    nbsamples_on_c1,
                    d1.first_parameter(),
                    d1.last_parameter(),
                    tol,
                );
            }
            if a_poly2.nb_segments() < nbsamples_on_c2 {
                a_poly2 = Polygon2dGen::<C>::new::<T>(
                    c2,
                    nbsamples_on_c2,
                    d2.first_parameter(),
                    d2.last_parameter(),
                    tol,
                );
            }

            self.find_intersect(
                c1, d1, c2, d2, tol_conf, tol, nb_iter, delta_u, delta_v, &a_poly1, &a_poly2,
                true,
            );
        }

        self.base.done = true;
    }

    /// OCCT findIntersect (gxx L1152-1569).
    #[allow(clippy::too_many_arguments)]
    fn find_intersect(
        &mut self,
        c1: &C,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        nb_iter: i32,
        delta_u: f64,
        delta_v: f64,
        the_poly1: &Polygon2dGen<C>,
        the_poly2: &Polygon2dGen<C>,
        is_full_polygon: bool,
    ) -> bool {
        let mut tan1;
        let mut tan2;
        let mut norm1 = DVec2::ZERO;
        let mut norm2 = DVec2::ZERO;
        let mut p1;
        let mut p2;
        let inter_pp = InterferencePolygon2d::new_polygon_polygon(the_poly1, the_poly2);
        let mut eip = ExactIntersectionPoint::<C, T>::new(c1, c2, tol_conf);
        // double U = 0., V = 0.;
        let mut an_error_occurred = false;
        self.base.done = true; // To prevent exception in nbp=NbPoints();
        //----------------------------------------------------------------------
        //-- Processing of SectionPoint
        //----------------------------------------------------------------------
        let nbsp = inter_pp.interf.nb_section_points();
        for sp in 1..=nbsp {
            let spnt = inter_pp.interf.pnt_value(sp);
            let (_typ, mut seg_index1, mut param_on1) = spnt.info_first();
            let (_typ2, mut seg_index2, mut param_on2) = spnt.info_second();
            let mut u = 0.0;
            let mut v = 0.0;
            let mut fct = super::int_curve_generics::DistBetweenPCurvesGen::<C, T>::new(c1, c2);
            eip.perform(
                the_poly1,
                the_poly2,
                &mut seg_index1,
                &mut seg_index2,
                &mut param_on1,
                &mut param_on2,
                &mut fct,
            );
            an_error_occurred = eip.an_error_occurred();

            if eip.nb_roots() == 0 && !is_full_polygon {
                return false;
            }

            if an_error_occurred {
                continue;
            }

            //--------------------------------------------------------------------
            //-- It is checked if the found point is really a root
            //--------------------------------------------------------------------

            let (u_r, v_r) = eip.roots();
            u = u_r;
            v = v_r;
            let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, u);
            p1 = p1_r;
            tan1 = tan1_r;
            let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c2, v);
            p2 = p2_r;
            tan2 = tan2_r;
            let mut dist = p1.distance(p2);
            if eip.nb_roots() == 0 && dist > tol_conf {
                let a_trans = Transition::empty();
                let mut a_p_int =
                    IntersectionPoint::new(p1, u, v, a_trans.clone(), a_trans.clone(), false);
                let a_t1f = the_poly1.approx_param_on_curve(seg_index1, 0.0);
                let a_t1l = the_poly1.approx_param_on_curve(seg_index1, 1.0);
                let a_t2f = the_poly2.approx_param_on_curve(seg_index2, 0.0);
                let a_t2l = the_poly2.approx_param_on_curve(seg_index2, 1.0);
                //
                let a_max_count = 16;
                let mut a_count = 0;
                get_intersection::<C, T>(
                    c1, a_t1f, a_t1l, c2, a_t2f, a_t2l, tol_conf, a_max_count, &mut a_p_int,
                    &mut dist, &mut a_count,
                );
                u = a_p_int.param_on_first();
                v = a_p_int.param_on_second();
                let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, u);
                p1 = p1_r;
                tan1 = tan1_r;
                let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c2, v);
                p2 = p2_r;
                tan2 = tan2_r;
                dist = p1.distance(p2);
            }
            //-----------------------------------------------------------------
            //-- It is checked if the point (u,v) does not exist already
            //--
            let nbp = self.base.nb_points();
            let eps_x1 = 10.0 * <T as ParTool<C>>::eps_x(c1);
            let eps_x2 = 10.0 * <T as ParTool<C>>::eps_x(c2);
            for p in 1..=nbp {
                let pt = self.base.point(p);
                if (u - pt.param_on_first()).abs() <= eps_x1
                    && (v - pt.param_on_second()).abs() <= eps_x2
                {
                    dist = tol_conf + 1.0;
                    break;
                }
            }

            if dist <= tol_conf {
                //-- Or the point is already present
                let mut pos1 = Position::Middle;
                let mut pos2 = Position::Middle;
                let mut trans1 = Transition::empty();
                let mut trans2 = Transition::empty();
                //-----------------------------------------------------------------
                //-- Calculate the Positions of Points on the curve
                //--
                if p1.distance(self.domain_on_curve1.first_point())
                    <= self.domain_on_curve1.first_tolerance()
                {
                    pos1 = Position::Head;
                } else if p1.distance(self.domain_on_curve1.last_point())
                    <= self.domain_on_curve1.last_tolerance()
                {
                    pos1 = Position::End;
                }

                if p2.distance(self.domain_on_curve2.first_point())
                    <= self.domain_on_curve2.first_tolerance()
                {
                    pos2 = Position::Head;
                } else if p2.distance(self.domain_on_curve2.last_point())
                    <= self.domain_on_curve2.last_tolerance()
                {
                    pos2 = Position::End;
                }
                //-----------------------------------------------------------------
                //-- Calculate the Transitions (see IntImpParGen.cxx)
                //--
                if !determine_transition_in_out(
                    pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2, tol_conf,
                ) {
                    let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, u);
                    p1 = p1_r;
                    tan1 = tan1_r;
                    norm1 = norm1_r;
                    let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c2, v);
                    p2 = p2_r;
                    tan2 = tan2_r;
                    norm2 = norm2_r;
                    determine_transition(
                        pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2, norm2, &mut trans2,
                        tol_conf,
                    );
                }
                let ip = IntersectionPoint::new(p1, u, v, trans1, trans2, false);
                self.base.insert(&ip);
            }
        }

        //----------------------------------------------------------------------
        //-- Processing of TangentZone
        //----------------------------------------------------------------------
        let nbtz = inter_pp.interf.nb_tangent_zones();
        for tz in 1..=nbtz {
            let nb_pnts = inter_pp.interf.zone_value(tz).number_of_points();
            //==================================================================
            //== Find the first and the last point in the tangency zone.
            //==================================================================
            let mut param_sup_on_curve1 = -f64::MAX;
            let mut param_sup_on_curve2 = -f64::MAX;
            let mut param_inf_on_curve1 = f64::MAX;
            let mut param_inf_on_curve2 = f64::MAX;
            for qq in 1..=nb_pnts {
                let spnt1 = inter_pp.interf.zone_value(tz).get_point(qq);
                //== The zones of tangency are discretized
                let (_typ, mut seg_index1on_p1, mut param_on_line) = spnt1.info_first();
                if seg_index1on_p1 > the_poly1.nb_segments() {
                    seg_index1on_p1 -= 1;
                    param_on_line = 1.0;
                }
                if seg_index1on_p1 <= 0 {
                    seg_index1on_p1 = 1;
                    param_on_line = 0.0;
                }
                let poly_u_inf_zone =
                    the_poly1.approx_param_on_curve(seg_index1on_p1, param_on_line);

                let (_typ, mut seg_index1on_p2, mut param_on_line2) = spnt1.info_second();
                if seg_index1on_p2 > the_poly2.nb_segments() {
                    seg_index1on_p2 -= 1;
                    param_on_line2 = 1.0;
                }
                if seg_index1on_p2 <= 0 {
                    seg_index1on_p2 = 1;
                    param_on_line2 = 0.0;
                }
                let poly_v_inf_zone =
                    the_poly2.approx_param_on_curve(seg_index1on_p2, param_on_line2);

                //----------------------------------------------------------------

                if param_inf_on_curve1 > poly_u_inf_zone {
                    param_inf_on_curve1 = poly_u_inf_zone;
                }
                if param_inf_on_curve2 > poly_v_inf_zone {
                    param_inf_on_curve2 = poly_v_inf_zone;
                }

                if param_sup_on_curve1 < poly_u_inf_zone {
                    param_sup_on_curve1 = poly_u_inf_zone;
                }
                if param_sup_on_curve2 < poly_v_inf_zone {
                    param_sup_on_curve2 = poly_v_inf_zone;
                }
            }

            let mut poly_u_inf = param_inf_on_curve1;
            let mut poly_u_sup = param_sup_on_curve1;
            let mut poly_v_inf = param_inf_on_curve2;
            let mut poly_v_sup = param_sup_on_curve2;

            p1 = <T as ParTool<C>>::value(c1, poly_u_inf);
            p2 = <T as ParTool<C>>::value(c2, poly_v_inf);
            let distmemesens = p1.distance_squared(p2);
            p2 = <T as ParTool<C>>::value(c2, poly_v_sup);
            let distdiffsens = p1.distance_squared(p2);
            if distmemesens > distdiffsens {
                let qwerty = poly_v_inf;
                poly_v_inf = poly_v_sup;
                poly_v_sup = qwerty;
            }

            if ((the_poly1.deflection_over_estimation() > tol_conf)
                || (the_poly2.deflection_over_estimation() > tol_conf))
                && (nb_iter < NBITER_MAX_POLYGON)
            {
                let recurs_d1 = Res2dDomain::bounded(
                    <T as ParTool<C>>::value(c1, param_inf_on_curve1),
                    param_inf_on_curve1,
                    tol_conf,
                    <T as ParTool<C>>::value(c1, param_sup_on_curve1),
                    param_sup_on_curve1,
                    tol_conf,
                );
                let recurs_d2 = Res2dDomain::bounded(
                    <T as ParTool<C>>::value(c2, param_inf_on_curve2),
                    param_inf_on_curve2,
                    tol_conf,
                    <T as ParTool<C>>::value(c2, param_sup_on_curve2),
                    param_sup_on_curve2,
                    tol_conf,
                );
                //-- thePoly1(2) are not deleted,
                //-- finally they are destroyed.
                //-- !! No untimely return !!
                // OCCT calls Perform(C1, RecursD1, C2, RecursD2, Tol, TolConf,
                // ...) — the Tol/TolConf swap is in the source (L1390).
                self.perform_pair_internal(
                    c1,
                    &recurs_d1,
                    c2,
                    &recurs_d2,
                    tol,      // Tol  -> the TolConf slot
                    tol_conf, // TolConf -> the Tol slot
                    nb_iter + 1,
                    delta_u,
                    delta_v,
                );
            } else {
                //-----------------------------------------------------------------
                //-- Calculate Positions of Points on the curve and
                //-- Transitions on each limit of the segment

                let mut pos1 = Position::Middle;
                let mut pos2 = Position::Middle;
                let mut trans1 = Transition::empty();
                let mut trans2 = Transition::empty();

                let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, poly_u_inf);
                p1 = p1_r;
                tan1 = tan1_r;
                let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c2, poly_v_inf);
                p2 = p2_r;
                tan2 = tan2_r;

                if p1.distance(self.domain_on_curve1.first_point())
                    <= self.domain_on_curve1.first_tolerance()
                {
                    pos1 = Position::Head;
                } else if p1.distance(self.domain_on_curve1.last_point())
                    <= self.domain_on_curve1.last_tolerance()
                {
                    pos1 = Position::End;
                }
                if p2.distance(self.domain_on_curve2.first_point())
                    <= self.domain_on_curve2.first_tolerance()
                {
                    pos2 = Position::Head;
                } else if p2.distance(self.domain_on_curve2.last_point())
                    <= self.domain_on_curve2.last_tolerance()
                {
                    pos2 = Position::End;
                }

                if pos1 == Position::Middle && pos2 != Position::Middle {
                    poly_u_inf = T::find_parameter_between(
                        c1,
                        p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                } else if pos1 != Position::Middle && pos2 == Position::Middle {
                    poly_v_inf = T::find_parameter_between(
                        c2,
                        p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
                        <T as ParTool<C>>::eps_x(c2),
                    );
                } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                    > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                {
                    poly_v_inf = T::find_parameter_between(
                        c2,
                        p1,
                        d2.first_parameter(),
                        d2.last_parameter(),
                        <T as ParTool<C>>::eps_x(c2),
                    );
                } else {
                    poly_u_inf = T::find_parameter_between(
                        c1,
                        p2,
                        d1.first_parameter(),
                        d1.last_parameter(),
                        <T as ParTool<C>>::eps_x(c1),
                    );
                }

                if !determine_transition_in_out(
                    pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2, tol_conf,
                ) {
                    let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, poly_u_inf);
                    p1 = p1_r;
                    tan1 = tan1_r;
                    norm1 = norm1_r;
                    let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c2, poly_v_inf);
                    p2 = p2_r;
                    tan2 = tan2_r;
                    norm2 = norm2_r;
                    determine_transition(
                        pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2, norm2, &mut trans2,
                        tol_conf,
                    );
                }
                let pt_seg1 = IntersectionPoint::new(
                    p1,
                    poly_u_inf,
                    poly_v_inf,
                    trans1.clone(),
                    trans2.clone(),
                    false,
                );
                //----------------------------------------------------------------------

                if (poly_u_inf - poly_u_sup).abs() <= <T as ParTool<C>>::eps_x(c1)
                    || (poly_v_inf - poly_v_sup).abs() <= <T as ParTool<C>>::eps_x(c2)
                {
                    self.base.insert(&pt_seg1);
                } else {
                    let (p1_r, tan1_r) = <T as ParTool<C>>::d1(c1, poly_u_sup);
                    p1 = p1_r;
                    tan1 = tan1_r;
                    let (p2_r, tan2_r) = <T as ParTool<C>>::d1(c2, poly_v_sup);
                    p2 = p2_r;
                    tan2 = tan2_r;
                    pos1 = Position::Middle;
                    pos2 = Position::Middle;

                    if p1.distance(self.domain_on_curve1.first_point())
                        <= self.domain_on_curve1.first_tolerance()
                    {
                        pos1 = Position::Head;
                    } else if p1.distance(self.domain_on_curve1.last_point())
                        <= self.domain_on_curve1.last_tolerance()
                    {
                        pos1 = Position::End;
                    }
                    if p2.distance(self.domain_on_curve2.first_point())
                        <= self.domain_on_curve2.first_tolerance()
                    {
                        pos2 = Position::Head;
                    } else if p2.distance(self.domain_on_curve2.last_point())
                        <= self.domain_on_curve2.last_tolerance()
                    {
                        pos2 = Position::End;
                    }

                    if pos1 == Position::Middle && pos2 != Position::Middle {
                        poly_u_sup = T::find_parameter_between(
                            c1,
                            p2,
                            d1.first_parameter(),
                            d1.last_parameter(),
                            <T as ParTool<C>>::eps_x(c1),
                        );
                    } else if pos1 != Position::Middle && pos2 == Position::Middle {
                        poly_v_sup = T::find_parameter_between(
                            c2,
                            p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
                            <T as ParTool<C>>::eps_x(c2),
                        );
                    } else if (param_inf_on_curve1 - param_sup_on_curve1).abs()
                        > (param_inf_on_curve2 - param_sup_on_curve2).abs()
                    {
                        poly_v_sup = T::find_parameter_between(
                            c2,
                            p1,
                            d2.first_parameter(),
                            d2.last_parameter(),
                            <T as ParTool<C>>::eps_x(c2),
                        );
                    } else {
                        poly_u_sup = T::find_parameter_between(
                            c1,
                            p2,
                            d1.first_parameter(),
                            d1.last_parameter(),
                            <T as ParTool<C>>::eps_x(c1),
                        );
                    }

                    if !determine_transition_in_out(
                        pos1, &mut tan1, &mut trans1, pos2, &mut tan2, &mut trans2, tol_conf,
                    ) {
                        let (p1_r, tan1_r, norm1_r) = <T as ParTool<C>>::d2(c1, poly_u_sup);
                        p1 = p1_r;
                        tan1 = tan1_r;
                        norm1 = norm1_r;
                        let (p2_r, tan2_r, norm2_r) = <T as ParTool<C>>::d2(c2, poly_v_sup);
                        p2 = p2_r;
                        tan2 = tan2_r;
                        norm2 = norm2_r;
                        determine_transition(
                            pos1, &mut tan1, norm1, &mut trans1, pos2, &mut tan2, norm2,
                            &mut trans2, tol_conf,
                        );
                    }
                    let pt_seg2 = IntersectionPoint::new(
                        p1,
                        poly_u_sup,
                        poly_v_sup,
                        trans1.clone(),
                        trans2.clone(),
                        false,
                    );

                    let oppos = !(tan1.dot(tan2) > 0.0);
                    if param_inf_on_curve1 > param_sup_on_curve1 {
                        let seg = IntersectionSegment::with_points(&pt_seg2, &pt_seg1, oppos, false);
                        self.base.append_segment(&seg);
                    } else {
                        let seg = IntersectionSegment::with_points(&pt_seg1, &pt_seg2, oppos, false);
                        self.base.append_segment(&seg);
                    }
                }
            }
        }
        true
    }
}

/// OCCT HeadOrEndPoint (gxx L872-1031) — the domain end-point intersection
/// test.  Returns true when a point has been produced in `int_pt`.
#[allow(clippy::too_many_arguments)]
fn head_or_end_point<C: ?Sized, T>(
    d1: &Res2dDomain,
    c1: &C,
    tu: f64,
    d2: &Res2dDomain,
    c2: &C,
    tv: f64,
    tol_conf: f64,
    int_pt: &mut IntersectionPoint,
    head_on1: &mut bool,
    head_on2: &mut bool,
    end_on1: &mut bool,
    end_on2: &mut bool,
    pos_segment: i32,
) -> bool
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    let mut p1;
    let mut p2;
    let mut sp1 = DVec2::ZERO;
    let mut sp2 = DVec2::ZERO;
    let mut t1;
    let mut t2;
    let mut n1 = DVec2::ZERO;
    let mut n2 = DVec2::ZERO;
    let mut u = tu;
    let mut v = tv;
    let svu = u;
    let svv = v;

    let (p1_r, t1_r) = <T as ParTool<C>>::d1(c1, u);
    p1 = p1_r;
    t1 = t1_r;
    let (p2_r, t2_r) = <T as ParTool<C>>::d1(c2, v);
    p2 = p2_r;
    t2 = t2_r;

    let mut pos1 = Position::Middle;
    let mut pos2 = Position::Middle;
    let mut trans1 = Transition::empty();
    let mut trans2 = Transition::empty();

    //----------------------------------------------------------------------
    //-- Head On 1   :        Head1 <-> P2
    if p2.distance(d1.first_point()) <= d1.first_tolerance() {
        pos1 = Position::Head;
        *head_on1 = true;
        sp1 = d1.first_point();
        u = d1.first_parameter();
    }
    //----------------------------------------------------------------------
    //-- End On 1   :         End1 <-> P2
    else if p2.distance(d1.last_point()) <= d1.last_tolerance() {
        pos1 = Position::End;
        *end_on1 = true;
        sp1 = d1.last_point();
        u = d1.last_parameter();
    }
    //----------------------------------------------------------------------
    //-- Head On 2   :        Head2 <-> P1
    else if p1.distance(d2.first_point()) <= d2.first_tolerance() {
        pos2 = Position::Head;
        *head_on2 = true;
        sp2 = d2.first_point();
        v = d2.first_parameter();
    }
    //----------------------------------------------------------------------
    //-- End On 2   :        End2 <-> P1
    else if p1.distance(d2.last_point()) <= d2.last_tolerance() {
        pos2 = Position::End;
        *end_on2 = true;
        sp2 = d2.last_point();
        v = d2.last_parameter();
    }

    let eps_x1 = <T as ParTool<C>>::eps_x(c1);
    let eps_x2 = <T as ParTool<C>>::eps_x(c2);

    if pos1 != Position::Middle || pos2 != Position::Middle {
        if pos1 == Position::Middle {
            if (u - d1.first_parameter()).abs() <= eps_x1 {
                pos1 = Position::Head;
                p1 = d1.first_point();
                *head_on1 = true;
            } else if (u - d1.last_parameter()).abs() <= eps_x1 {
                pos1 = Position::End;
                p1 = d1.last_point();
                *end_on1 = true;
            }
        } else if u != tu {
            p1 = sp1;
        }

        if pos2 == Position::Middle {
            if (v - d2.first_parameter()).abs() <= eps_x2 {
                pos2 = Position::Head;
                *head_on2 = true;
                p2 = d2.first_point();
                if pos1 != Position::Middle {
                    p1 = DVec2::new(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y));
                } else {
                    p2 = p1;
                }
            } else if (v - d2.last_parameter()).abs() <= eps_x2 {
                pos2 = Position::End;
                *end_on2 = true;
                p2 = d2.last_point();
                if pos1 != Position::Middle {
                    p1 = DVec2::new(0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y));
                } else {
                    p2 = p1;
                }
            }
        }

        //--------------------------------------------------------------------
        //-- It is tested if a point at the end of segment already has its
        //-- transitions.  If Yes, the new point is not created.
        //-- PosSegment =            1    if Head Head
        //--                       2      if Head End
        //--                     4        if End  Head
        //--                   8          if End  End
        //--------------------------------------------------------------------
        if pos1 == Position::Head {
            if pos2 == Position::Head && (pos_segment & 1) != 0 {
                return false;
            }
            if pos2 == Position::End && (pos_segment & 2) != 0 {
                return false;
            }
        } else if pos1 == Position::End {
            if pos2 == Position::Head && (pos_segment & 4) != 0 {
                return false;
            }
            if pos2 == Position::End && (pos_segment & 8) != 0 {
                return false;
            }
        }

        if !determine_transition_in_out(
            pos1, &mut t1, &mut trans1, pos2, &mut t2, &mut trans2, tol_conf,
        ) {
            let (p1_r, t1_r, n1_r) = <T as ParTool<C>>::d2(c1, svu);
            p1 = p1_r;
            t1 = t1_r;
            n1 = n1_r;
            let (p2_r, t2_r, n2_r) = <T as ParTool<C>>::d2(c2, svv);
            p2 = p2_r;
            t2 = t2_r;
            n2 = n2_r;
            determine_transition(pos1, &mut t1, n1, &mut trans1, pos2, &mut t2, n2, &mut trans2, tol_conf);
        }
        int_pt.set_values(p1, u, v, trans1, trans2, false);
        true
    } else {
        false
    }
}

/// OCCT GetIntersection (gxx L1575-1783) — the recursive interval
/// bisection used when the exact solver found no root: narrows the pair of
/// parameter brackets until the closest sample pair is inside tolerance.
#[allow(clippy::too_many_arguments)]
fn get_intersection<C: ?Sized, T>(
    the_c1: &C,
    the_t1f: f64,
    the_t1l: f64,
    the_c2: &C,
    the_t2f: f64,
    the_t2l: f64,
    the_tol_conf: f64,
    the_max_count: i32,
    the_p_int: &mut IntersectionPoint,
    the_dist: &mut f64,
    the_count: &mut i32,
) where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    *the_count += 1;
    //
    let a_tol2 = the_tol_conf * the_tol_conf;
    let a_ptol1 = (100.0 * epsilon((the_t1f.abs()).max(the_t1l.abs()))).max(PRECISION_P_CONFUSION);
    let a_ptol2 = (100.0 * epsilon((the_t2f.abs()).max(the_t2l.abs()))).max(PRECISION_P_CONFUSION);
    let a_p1f = <T as ParTool<C>>::value(the_c1, the_t1f);
    let a_p1l = <T as ParTool<C>>::value(the_c1, the_t1l);
    let mut a_b1 = BndBox2d::new();
    a_b1.add_point(a_p1f);
    a_b1.add_point(a_p1l);
    a_b1.enlarge(the_tol_conf);
    //
    let a_p2f = <T as ParTool<C>>::value(the_c2, the_t2f);
    let a_p2l = <T as ParTool<C>>::value(the_c2, the_t2l);
    let mut a_b2 = BndBox2d::new();
    a_b2.add_point(a_p2f);
    a_b2.add_point(a_p2l);
    a_b2.enlarge(the_tol_conf);
    //
    if a_b1.is_out_box(&a_b2) {
        *the_count -= 1;
        return;
    }
    //
    let is_small1 = (the_t1l - the_t1f) <= a_ptol1 || a_p1f.distance_squared(a_p1l) / 4.0 <= a_tol2;
    let is_small2 = (the_t2l - the_t2f) <= a_ptol2 || a_p2f.distance_squared(a_p2l) / 4.0 <= a_tol2;

    if (is_small1 && is_small2) || *the_count > the_max_count {
        // Seems to be intersection
        // Simple treatment of segment intersection
        let a_pnts1 = [
            a_p1f,
            (a_p1f + a_p1l) * 0.5,
            a_p1l,
        ];
        let a_pnts2 = [
            a_p2f,
            (a_p2f + a_p2l) * 0.5,
            a_p2l,
        ];
        let mut imin = 0usize;
        let mut jmin = 0usize;
        let mut dmin = f64::MAX;
        for (i, a) in a_pnts1.iter().enumerate() {
            for (j, b) in a_pnts2.iter().enumerate() {
                let d = (*a - *b).length_squared();
                if d < dmin {
                    dmin = d;
                    imin = i;
                    jmin = j;
                }
            }
        }
        //
        dmin = dmin.sqrt();
        if *the_dist > dmin {
            *the_dist = dmin;
            //
            let t1 = match imin {
                0 => the_t1f,
                1 => (the_t1f + the_t1l) / 2.0,
                _ => the_t1l,
            };
            //
            let t2 = match jmin {
                0 => the_t2f,
                1 => (the_t2f + the_t2l) / 2.0,
                _ => the_t2l,
            };
            //
            let a_pint = (a_pnts1[imin] + a_pnts2[jmin]) * 0.5;
            //
            let a_trans1 = Transition::empty();
            let a_trans2 = Transition::empty();
            the_p_int.set_values(a_pint, t1, t2, a_trans1, a_trans2, false);
        }
        *the_count -= 1;
        return;
    }

    if is_small1 {
        let a_t2m = (the_t2l + the_t2f) / 2.0;
get_intersection::<C, T>(
            the_c1, the_t1f, the_t1l, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
get_intersection::<C, T>(
            the_c1, the_t1f, the_t1l, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
    } else if is_small2 {
        let a_t1m = (the_t1l + the_t1f) / 2.0;
get_intersection::<C, T>(
            the_c1, the_t1f, a_t1m, the_c2, the_t2f, the_t2l, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
get_intersection::<C, T>(
            the_c1, a_t1m, the_t1l, the_c2, the_t2f, the_t2l, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
    } else {
        let a_t1m = (the_t1l + the_t1f) / 2.0;
        let a_t2m = (the_t2l + the_t2f) / 2.0;
get_intersection::<C, T>(
            the_c1, the_t1f, a_t1m, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
get_intersection::<C, T>(
            the_c1, the_t1f, a_t1m, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
get_intersection::<C, T>(
            the_c1, a_t1m, the_t1l, the_c2, the_t2f, a_t2m, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
get_intersection::<C, T>(
            the_c1, a_t1m, the_t1l, the_c2, a_t2m, the_t2l, the_tol_conf, the_max_count,
            the_p_int, the_dist, the_count,
        );
    }
}

impl<C: ?Sized, T> Default for IntPolyPolyGen<C, T>
where
    T: ProjPCurveTool<Curve = C> + ParTool<C> + ProjectOnPCurveTool<C>,
{
    /// OCCT IntCurve_IntPolyPolyGen() (gxx L84-89).
    fn default() -> Self {
        IntPolyPolyGen::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::geom2d_int::{Curve2dAdaptor, TheIntPCurvePCurveOfGInter};
    use rcad_kernel::geom::{BezierCurve2, Curve2d, Line2d};

    const TOL: f64 = 1.0e-6;

    /// OCCT Geom2dInt_TheIntPCurvePCurveOfGInter::Perform(C1, D1, C2, D2,
    /// TolConf, Tol): two straight parametrised segments crossing at (1, 1)
    /// give one exact transversal point.
    #[test]
    fn pair_crossing_lines() {
        let c1 = Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::new(1.0, 1.0)));
        let c2 = Curve2d::Line(Line2d::new(DVec2::new(2.0, 0.0), DVec2::new(-1.0, 1.0)));
        let d1 = Res2dDomain::bounded(c1.value(0.0), 0.0, TOL, c1.value(2.0 * 2.0f64.sqrt()), 2.0 * 2.0f64.sqrt(), TOL);
        let d2 = Res2dDomain::bounded(c2.value(0.0), 0.0, TOL, c2.value(2.0 * 2.0f64.sqrt()), 2.0 * 2.0f64.sqrt(), TOL);
        let mut itp = TheIntPCurvePCurveOfGInter::new();
        itp.perform(&c1, &d1, &c2, &d2, TOL, TOL);
        assert!(itp.base.is_done());
        assert_eq!(itp.base.nb_points(), 1, "got {}", itp.base.nb_points());
        let pt = itp.base.point(1);
        let p = pt.value();
        assert!((p.x - 1.0).abs() < 1.0e-7 && (p.y - 1.0).abs() < 1.0e-7, "got {p:?}");
        let s = 2.0f64.sqrt();
        assert!((pt.param_on_first() - s).abs() < 1.0e-6, "u={}", pt.param_on_first());
        assert!((pt.param_on_second() - s).abs() < 1.0e-6, "v={}", pt.param_on_second());
        // Transversal crossing: In on one curve, Out on the other.
        let t1 = pt.transition_of_first();
        let t2 = pt.transition_of_second();
        assert!((t1.is_tangent() == t2.is_tangent()) && !t1.is_tangent());
    }

    /// OCCT Perform(C1, D1, TheTolConf, TheTol): a degree-1 open figure-eight
    /// BSpline self-intersects at (1, 1) with parameters 1/6 and 5/6.
    #[test]
    fn auto_intersection_figure_eight() {
        let third = 1.0 / 3.0;
        let bs = rcad_kernel::geom::BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, third, 2.0 * third, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(2.0, 2.0),
                DVec2::new(2.0, 0.0),
                DVec2::new(0.0, 2.0),
            ],
            weights: vec![1.0; 4],
        };
        let c1 = Curve2d::BSpline(bs);
        let d1 = Res2dDomain::bounded(c1.value(0.0), 0.0, TOL, c1.value(1.0), 1.0, TOL);
        let mut itp = TheIntPCurvePCurveOfGInter::new();
        itp.perform_auto(&c1, &d1, TOL, TOL);
        assert!(itp.base.is_done());
        assert_eq!(itp.base.nb_points(), 1, "got {}", itp.base.nb_points());
        let pt = itp.base.point(1);
        let p = pt.value();
        assert!((p.x - 1.0).abs() < 1.0e-6 && (p.y - 1.0).abs() < 1.0e-6, "got {p:?}");
        assert!((pt.param_on_first() - 1.0 / 6.0).abs() < 1.0e-5, "u={}", pt.param_on_first());
        assert!((pt.param_on_second() - 5.0 / 6.0).abs() < 1.0e-5, "v={}", pt.param_on_second());
    }

    /// Disjoint curves produce an empty (but done) result.
    #[test]
    fn pair_disjoint() {
        let c1 = Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X));
        let c2 = Curve2d::Line(Line2d::new(DVec2::new(0.0, 5.0), DVec2::X));
        let d1 = Res2dDomain::bounded(DVec2::ZERO, 0.0, TOL, DVec2::new(1.0, 0.0), 1.0, TOL);
        let d2 = Res2dDomain::bounded(DVec2::new(0.0, 5.0), 0.0, TOL, DVec2::new(1.0, 5.0), 1.0, TOL);
        let mut itp = TheIntPCurvePCurveOfGInter::new();
        itp.perform(&c1, &d1, &c2, &d2, TOL, TOL);
        assert!(itp.base.is_done());
        assert_eq!(itp.base.nb_points(), 0);
        assert_eq!(itp.base.nb_segments(), 0);
    }

    /// A parabola arc and its mirror cross twice at exact points.
    #[test]
    fn pair_crossing_beziers() {
        // Two quadratic Beziers forming an X: C1 from (-2,-2) to (2,2),
        // C2 from (-2,2) to (2,-2); both pass through the origin.
        let bez1 = BezierCurve2 {
            control_points: vec![DVec2::new(-2.0, -2.0), DVec2::ZERO, DVec2::new(2.0, 2.0)],
            weights: vec![1.0; 3],
        };
        let bez2 = BezierCurve2 {
            control_points: vec![DVec2::new(-2.0, 2.0), DVec2::ZERO, DVec2::new(2.0, -2.0)],
            weights: vec![1.0; 3],
        };
        let c1 = Curve2d::Bezier(bez1);
        let c2 = Curve2d::Bezier(bez2);
        let d1 = Res2dDomain::bounded(c1.value(0.0), 0.0, TOL, c1.value(1.0), 1.0, TOL);
        let d2 = Res2dDomain::bounded(c2.value(0.0), 0.0, TOL, c2.value(1.0), 1.0, TOL);
        let mut itp = TheIntPCurvePCurveOfGInter::new();
        itp.perform(&c1, &d1, &c2, &d2, TOL, TOL);
        assert!(itp.base.is_done());
        assert_eq!(itp.base.nb_points(), 1, "got {}", itp.base.nb_points());
        let pt = itp.base.point(1);
        let p = pt.value();
        assert!(p.length() < 1.0e-6, "got {p:?}");
        assert!((pt.param_on_first() - 0.5).abs() < 1.0e-6);
        assert!((pt.param_on_second() - 0.5).abs() < 1.0e-6);
    }

    /// GetMinNbSamples / SetMinNbSamples (gxx L1787-1797).
    #[test]
    fn min_nb_samples_accessors() {
        let itp = TheIntPCurvePCurveOfGInter::new();
        assert_eq!(itp.get_min_nb_samples(), 20);
        let mut itp2 = TheIntPCurvePCurveOfGInter::new();
        itp2.set_min_nb_samples(50);
        assert_eq!(itp2.get_min_nb_samples(), 50);
    }
}
