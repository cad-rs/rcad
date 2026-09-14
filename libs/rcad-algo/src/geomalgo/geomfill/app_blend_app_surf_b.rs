//! OCCT AppBlend_AppSurf — the InternalPerform body
//! (AppBlend_AppSurf.gxx L188-593) of the engine declared in the parent
//! module `app_blend_app_surf`.  GAP note: the UseSmoothing branch drives
//! AppDef_Variational, kept as the [`AppDefVariational`] failure-path
//! carrier (see the parent module header).

use glam::{DVec2, DVec3};
use rcad_kernel::math::math_matrix::Vector as RVector;
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::app_def::{MultiLine, MultiPointConstraint};
use crate::geomalgo::app_def_compute::Compute;
use crate::geomalgo::approx_int::{
    ApproxParamType, AppParConstraint, ConstraintCouple, MultiBSpCurve,
};
use crate::geomalgo::bspl_compute_line::BSplineCompute;

use super::{AppBlendAppSurf, TheSectionGenerator, GP_RESOLUTION};
use super::super::line::Line;

/// GAP carrier: OCCT AppDef_Variational (TKGeomBase/AppDef) — the
/// variational-smoothing approximation has no rcad port yet.  The carrier
/// preserves the OCCT failure path of the gxx flow: IsCreated() = false
/// makes the translated body return early with done = false, exactly like a
/// non-created Variational run (the points_to_bspline.rs convention).  The
/// OCCT call shape is recorded at the use site.
pub(super) struct AppDefVariational;

impl AppDefVariational {
    /// OCCT AppDef_Variational(MultiLine, FirstPoint, LastPoint,
    /// Constraints) (AppDef_Variational.cxx).
    #[allow(unused_variables)]
    pub(super) fn new(
        line: &MultiLine,
        first_point: i32,
        last_point: i32,
        constraints: &[ConstraintCouple],
    ) -> Self {
        AppDefVariational
    }

    /// OCCT SetMaxDegree(Degree).
    #[allow(unused_variables)]
    pub(super) fn set_max_degree(&mut self, degree: i32) {}

    /// OCCT SetContinuity(Continuity).
    #[allow(unused_variables)]
    pub(super) fn set_continuity(&mut self, continuity: GeomAbsShape) {}

    /// OCCT SetMaxSegment(MaxSegment).
    #[allow(unused_variables)]
    pub(super) fn set_max_segment(&mut self, max_segment: i32) {}

    /// OCCT SetTolerance(Tol).
    #[allow(unused_variables)]
    pub(super) fn set_tolerance(&mut self, tol: f64) {}

    /// OCCT SetWithMinMax(WithMinMax).
    #[allow(unused_variables)]
    pub(super) fn set_with_min_max(&mut self, with_min_max: bool) {}

    /// OCCT SetWithCutting(WithCutting).
    #[allow(unused_variables)]
    pub(super) fn set_with_cutting(&mut self, with_cutting: bool) {}

    /// OCCT SetNbIterations(NbIterations).
    #[allow(unused_variables)]
    pub(super) fn set_nb_iterations(&mut self, nb_iterations: i32) {}

    /// OCCT SetCriteriumWeight(W1, W2, W3).
    #[allow(unused_variables)]
    pub(super) fn set_criterium_weight(&mut self, w1: f64, w2: f64, w3: f64) {}

    /// OCCT IsCreated() — the GAP failure path (false).
    pub(super) fn is_created(&self) -> bool {
        false
    }

    /// OCCT IsOverConstrained() — behind the IsCreated() early return.
    pub(super) fn is_over_constrained(&self) -> bool {
        false
    }

    /// OCCT Approximate() — behind the IsCreated() early return.
    pub(super) fn approximate(&mut self) {
        panic!("GAP: AppDef_Variational (TKGeomBase/AppDef) is not translated — see file header")
    }

    /// OCCT IsDone() — behind the IsCreated() early return.
    pub(super) fn is_done(&self) -> bool {
        panic!("GAP: AppDef_Variational (TKGeomBase/AppDef) is not translated — see file header")
    }

    /// OCCT MaxError() — behind the IsCreated() early return.
    pub(super) fn max_error(&self) -> f64 {
        panic!("GAP: AppDef_Variational (TKGeomBase/AppDef) is not translated — see file header")
    }

    /// OCCT Value() — behind the IsCreated() early return.
    pub(super) fn value(&self) -> MultiBSpCurve {
        panic!("GAP: AppDef_Variational (TKGeomBase/AppDef) is not translated — see file header")
    }
}

impl AppBlendAppSurf {
    /// OCCT AppBlend_AppSurf::InternalPerform (gxx L188-593).
    pub(super) fn internal_perform<G: TheSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &mut G,
        sp_approx: bool,
        use_smoothing: bool,
    ) {
        self.done = false;
        // OCCT L195-198: if (Lin.IsNull()) { return; } — the handle null
        // check is unrepresentable over the rcad value Line.

        let nb_point;
        let mut mult_p;
        let mut mult_l;
        let mut withderiv;
        let cfirst: AppParConstraint;
        let clast: AppParConstraint;

        let mut mytol3d: f64;
        let mut mytol2d: f64;
        let mut new_dv: DVec3;

        self.seq_poles2d.clear();

        nb_point = lin.nb_points();
        // OCCT L210: AppDef_MultiPointConstraint multP; — the default
        // ctor; rcad defers the init (always overwritten before the first
        // read at SetValue(1)).
        mult_l = MultiLine::new_nb_mult(nb_point as usize);

        let (mut nb_u_poles, mut nb_u_knots, mut nb_poles2d) = (0i32, 0i32, 0i32);
        f.get_shape(&mut nb_u_poles, &mut nb_u_knots, &mut self.udeg, &mut nb_poles2d);

        self.tab_u_knots = Some(vec![0.0; nb_u_knots as usize]);
        self.tab_u_mults = Some(vec![0; nb_u_knots as usize]);

        f.knots(self.tab_u_knots.as_mut().unwrap());
        f.mults(self.tab_u_mults.as_mut().unwrap());

        let mut tab_app_p = vec![DVec3::ZERO; nb_u_poles as usize];
        let mut tab_app_v = vec![DVec3::ZERO; nb_u_poles as usize];

        let mut tab_p2d = vec![DVec2::ZERO; nb_poles2d.max(1) as usize];
        let mut tab_v2d = vec![DVec2::ZERO; nb_poles2d.max(1) as usize];

        let mut tab_w = vec![0.0f64; nb_u_poles as usize];
        let mut tab_dw = vec![0.0f64; nb_u_poles as usize];

        let mut tab_app_p2d = vec![DVec2::ZERO; (nb_poles2d + nb_u_poles) as usize];
        let mut tab_app_v2d = vec![DVec2::ZERO; (nb_poles2d + nb_u_poles) as usize];

        // OCCT L232: AppParCurves_MultiBSpCurve multC; — the default
        // ctor; rcad defers the init (each branch assigns before use).
        let mult_c: MultiBSpCurve;

        // OCCT gxx L234: //  bool SpApprox = false;

        withderiv = f.section_d1(
            lin.point(1),
            &mut tab_app_p,
            &mut tab_app_v,
            &mut tab_p2d,
            &mut tab_v2d,
            &mut tab_w,
            &mut tab_dw,
        );

        if super::app_blend_get_context_approx_with_no_tgt() {
            withderiv = false;
        }

        for j in 1..=nb_poles2d {
            tab_app_p2d[(j - 1) as usize] = tab_p2d[(j - 1) as usize];
            if withderiv {
                tab_app_v2d[(j - 1) as usize] = tab_v2d[(j - 1) as usize];
            }
        }
        for j in 1..=nb_u_poles {
            // pour les courbes rationnelles il faut multiplier les poles par
            // leurs poids respectifs
            if withderiv {
                tab_app_v2d[(nb_poles2d + j - 1) as usize] =
                    DVec2::new(tab_dw[(j - 1) as usize], 0.0);
                new_dv = tab_app_p[(j - 1) as usize] * tab_dw[(j - 1) as usize]
                    + tab_app_v[(j - 1) as usize] * tab_w[(j - 1) as usize];
                tab_app_v[(j - 1) as usize] = new_dv;
            }
            tab_app_p[(j - 1) as usize] *= tab_w[(j - 1) as usize];
            tab_app_p2d[(nb_poles2d + j - 1) as usize] =
                DVec2::new(tab_w[(j - 1) as usize], 0.0);
        }

        if withderiv {
            mult_p = MultiPointConstraint::new_tangency(
                &tab_app_p,
                &tab_app_p2d,
                &tab_app_v,
                &tab_app_v2d,
            );
            cfirst = AppParConstraint::TangencyPoint;
        } else {
            mult_p = MultiPointConstraint::new_tab_p_p2d(&tab_app_p, &tab_app_p2d);
            cfirst = AppParConstraint::PassPoint;
        }
        mult_l.set_value(1, &mult_p);

        for i in 2..=nb_point - 1 {
            if sp_approx {
                f.section(lin.point(i), &mut tab_app_p, &mut tab_p2d, &mut tab_w);
                for j in 1..=nb_poles2d {
                    tab_app_p2d[(j - 1) as usize] = tab_p2d[(j - 1) as usize];
                }
                for j in 1..=nb_u_poles {
                    // pour les courbes rationnelles il faut multiplier les poles par
                    // leurs poids respectifs
                    tab_app_p[(j - 1) as usize] *= tab_w[(j - 1) as usize];
                    tab_app_p2d[(nb_poles2d + j - 1) as usize] =
                        DVec2::new(tab_w[(j - 1) as usize], 0.0);
                }
                mult_p = MultiPointConstraint::new_tab_p_p2d(&tab_app_p, &tab_app_p2d);
                mult_l.set_value(i as usize, &mult_p);
            }
            // ***********************
            else {
                withderiv = f.section_d1(
                    lin.point(i),
                    &mut tab_app_p,
                    &mut tab_app_v,
                    &mut tab_p2d,
                    &mut tab_v2d,
                    &mut tab_w,
                    &mut tab_dw,
                );
                if super::app_blend_get_context_approx_with_no_tgt() {
                    withderiv = false;
                }

                for j in 1..=nb_poles2d {
                    tab_app_p2d[(j - 1) as usize] = tab_p2d[(j - 1) as usize];
                    if withderiv {
                        tab_app_v2d[(j - 1) as usize] = tab_v2d[(j - 1) as usize];
                    }
                }
                for j in 1..=nb_u_poles {
                    // pour les courbes rationnelles il faut multiplier les poles par
                    // leurs poids respectifs
                    if withderiv {
                        tab_app_v2d[(nb_poles2d + j - 1) as usize] =
                            DVec2::new(tab_dw[(j - 1) as usize], 0.0);
                        new_dv = tab_app_p[(j - 1) as usize] * tab_dw[(j - 1) as usize]
                            + tab_app_v[(j - 1) as usize] * tab_w[(j - 1) as usize];
                        tab_app_v[(j - 1) as usize] = new_dv;
                    }
                    tab_app_p[(j - 1) as usize] *= tab_w[(j - 1) as usize];
                    tab_app_p2d[(nb_poles2d + j - 1) as usize] =
                        DVec2::new(tab_w[(j - 1) as usize], 0.0);
                }
                if withderiv {
                    mult_p = MultiPointConstraint::new_tangency(
                        &tab_app_p,
                        &tab_app_p2d,
                        &tab_app_v,
                        &tab_app_v2d,
                    );
                } else {
                    mult_p = MultiPointConstraint::new_tab_p_p2d(&tab_app_p, &tab_app_p2d);
                }
                mult_l.set_value(i as usize, &mult_p);
            }
            // ******************************
        }

        withderiv = f.section_d1(
            lin.point(nb_point),
            &mut tab_app_p,
            &mut tab_app_v,
            &mut tab_p2d,
            &mut tab_v2d,
            &mut tab_w,
            &mut tab_dw,
        );
        if super::app_blend_get_context_approx_with_no_tgt() {
            withderiv = false;
        }

        for j in 1..=nb_poles2d {
            tab_app_p2d[(j - 1) as usize] = tab_p2d[(j - 1) as usize];
            if withderiv {
                tab_app_v2d[(j - 1) as usize] = tab_v2d[(j - 1) as usize];
            }
        }
        for j in 1..=nb_u_poles {
            // pour les courbes rationnelles il faut multiplier les poles par
            // leurs poids respectifs
            if withderiv {
                tab_app_v2d[(nb_poles2d + j - 1) as usize] =
                    DVec2::new(tab_dw[(j - 1) as usize], 0.0);
                new_dv = tab_app_p[(j - 1) as usize] * tab_dw[(j - 1) as usize]
                    + tab_app_v[(j - 1) as usize] * tab_w[(j - 1) as usize];
                tab_app_v[(j - 1) as usize] = new_dv;
            }
            tab_app_p[(j - 1) as usize] *= tab_w[(j - 1) as usize];
            tab_app_p2d[(nb_poles2d + j - 1) as usize] = DVec2::new(tab_w[(j - 1) as usize], 0.0);
        }

        if withderiv {
            mult_p = MultiPointConstraint::new_tangency(
                &tab_app_p,
                &tab_app_p2d,
                &tab_app_v,
                &tab_app_v2d,
            );
            clast = AppParConstraint::TangencyPoint;
        } else {
            mult_p = MultiPointConstraint::new_tab_p_p2d(&tab_app_p, &tab_app_p2d);
            clast = AppParConstraint::PassPoint;
        }
        mult_l.set_value(nb_point as usize, &mult_p);

        // IFV 04.06.07 occ13904
        if nb_point == 2 {
            self.dmin = 1;
            if cfirst == AppParConstraint::PassPoint && clast == AppParConstraint::PassPoint {
                self.dmax = 1;
            }
        }

        if !sp_approx {
            // OCCT L385: AppDef_Compute theapprox(dmin, dmax, tol3d, tol2d,
            // nbit, true, paramtype) — the trailing `false` is the OCCT
            // Squares default.
            let mut theapprox = Compute::new(
                self.dmin,
                self.dmax,
                self.tol3d,
                self.tol2d,
                self.nbit,
                true,
                self.paramtype,
                false,
            );
            if self.knownp {
                let mut the_params = RVector::new(1, nb_point);

                // On recale les parametres entre 0 et 1.
                the_params.set(1, 0.0);
                the_params.set(nb_point, 1.0);
                let u_f = f.parameter(lin.point(1));
                let u_l = f.parameter(lin.point(nb_point)) - u_f;
                for i in 2..nb_point {
                    the_params.set(i, (f.parameter(lin.point(i)) - u_f) / u_l);
                }
                let the_app_def = Compute::new_with_parameters(
                    &the_params,
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    true,
                );
                theapprox = the_app_def;
            }
            theapprox.set_constraints(cfirst, clast);
            theapprox.perform(&mult_l);

            let mut the_tol3d;
            let mut the_tol2d;
            mytol3d = 0.0;
            mytol2d = 0.0;
            for index in 1..=theapprox.nb_multi_curves() {
                (the_tol3d, the_tol2d) = theapprox.error(index);
                mytol3d = the_tol3d.max(mytol3d);
                mytol2d = the_tol2d.max(mytol2d);
            }
            mult_c = theapprox.spline_value().clone();
        } else {
            if !use_smoothing {
                let mut use_squares = false;
                if self.nbit == 0 {
                    use_squares = true;
                }
                let mut theapprox = BSplineCompute::new(
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    self.paramtype,
                    use_squares,
                );
                // OCCT L428-443 quirk kept literally: the first `if` is not
                // part of the following if/else-if chain, so a C0 continuity
                // is set to 0 and then overwritten by the chain's `else`
                // arm (SetContinuity(3)).
                if self.continuity == GeomAbsShape::C0 {
                    theapprox.set_continuity(0);
                }
                if self.continuity == GeomAbsShape::C1 {
                    theapprox.set_continuity(1);
                } else if self.continuity == GeomAbsShape::C2 {
                    theapprox.set_continuity(2);
                } else {
                    theapprox.set_continuity(3);
                }

                theapprox.set_constraints(cfirst, clast);

                if self.knownp {
                    let mut the_params = RVector::new(1, nb_point);
                    // On recale les parametres entre 0 et 1.
                    the_params.set(1, 0.0);
                    the_params.set(nb_point, 1.0);
                    let u_f = f.parameter(lin.point(1));
                    let u_l = f.parameter(lin.point(nb_point)) - u_f;
                    for i in 2..nb_point {
                        the_params.set(i, (f.parameter(lin.point(i)) - u_f) / u_l);
                    }

                    theapprox.init(
                        self.dmin,
                        self.dmax,
                        self.tol3d,
                        self.tol2d,
                        self.nbit,
                        true,
                        ApproxParamType::IsoParametric,
                        true,
                    );
                    theapprox.set_parameters(&the_params);
                }
                theapprox.perform(&mult_l);
                let (the_tol3d, the_tol2d) = theapprox.error();
                mytol3d = the_tol3d;
                mytol2d = the_tol2d;
                self.tol3dreached = mytol3d;
                self.tol2dreached = mytol2d;
                mult_c = theapprox.value().clone();
            } else {
                // Variational algo
                // OCCT L476-477: the (1, NbPoint) array of constraint
                // couples.
                let mut tabofcc: Vec<ConstraintCouple> =
                    vec![ConstraintCouple { index: 0, constraint: AppParConstraint::NoConstraint };
                        nb_point as usize];
                let constraint = AppParConstraint::NoConstraint;

                for i in 1..=nb_point {
                    let acc = ConstraintCouple {
                        index: i,
                        constraint,
                    };
                    tabofcc[(i - 1) as usize] = acc;
                }

                tabofcc[0].constraint = cfirst;
                tabofcc[(nb_point - 1) as usize].constraint = clast;

                // GAP: AppDef_Variational not ported — the translated call
                // shape below keeps the OCCT flow; IsCreated() = false makes
                // every arm behave like the OCCT failure path (done = false).
                let mut variation =
                    AppDefVariational::new(&mult_l, 1, nb_point, &tabofcc);

                //===================================
                let the_max_segments = 1000i32;
                let the_with_min_max = false;
                let the_with_cutting = true;
                //===================================

                variation.set_max_degree(self.dmax);
                variation.set_continuity(self.continuity);
                variation.set_max_segment(the_max_segments);

                variation.set_tolerance(self.tol3d);
                variation.set_with_min_max(the_with_min_max);
                variation.set_with_cutting(the_with_cutting);
                variation.set_nb_iterations(self.nbit);

                variation.set_criterium_weight(
                    self.critweights[0],
                    self.critweights[1],
                    self.critweights[2],
                );

                if !variation.is_created() {
                    return;
                }

                if variation.is_over_constrained() {
                    return;
                }

                // OCCT L518-525: try { Variation.Approximate(); }
                // catch (Standard_Failure const&) { return; }
                variation.approximate();

                if !variation.is_done() {
                    return;
                }

                mytol3d = variation.max_error();
                mytol2d = 0.0;
                self.tol3dreached = mytol3d;
                self.tol2dreached = mytol2d;
                mult_c = variation.value();
            }
        }

        self.vdeg = mult_c.degree as i32;
        let nb_v_poles = mult_c.nb_poles();

        self.tab_poles = Some(vec![vec![DVec3::ZERO; nb_v_poles]; nb_u_poles as usize]);
        self.tab_weights = Some(vec![vec![0.0; nb_v_poles]; nb_u_poles as usize]);
        self.tab_v_knots = Some(mult_c.knots.clone());

        if self.knownp && !use_smoothing {
            let u1 = f.parameter(lin.point(1));
            let u2 = f.parameter(lin.point(nb_point));
            Self::bspl_reparametrize(u1, u2, self.tab_v_knots.as_mut().unwrap());
        }

        self.tab_v_mults = Some(mult_c.mults.iter().map(|m| *m as i32).collect());

        let mut newtab_p = vec![DVec3::ZERO; nb_v_poles];
        let mut newtab_p2d: Vec<DVec2> = Vec::new();
        for j in 1..=nb_u_poles {
            mult_c.curve(j as usize, &mut newtab_p);
            mult_c.curve2d((j + nb_u_poles + nb_poles2d) as usize, &mut newtab_p2d);
            for k in 1..=nb_v_poles {
                // pour les courbes rationnelles il faut maintenant diviser
                // les poles par leurs poids respectifs
                let a_weight = newtab_p2d[(k - 1) as usize].x;
                self.tab_poles.as_mut().unwrap()[(j - 1) as usize][(k - 1) as usize] =
                    newtab_p[(k - 1) as usize] / a_weight;
                if a_weight < GP_RESOLUTION {
                    self.done = false;
                    return;
                }
                self.tab_weights.as_mut().unwrap()[(j - 1) as usize][(k - 1) as usize] = a_weight;
            }
        }

        for j in 1..=nb_poles2d {
            let mut newtab_p2d = vec![DVec2::ZERO; nb_v_poles];
            mult_c.curve2d((nb_u_poles + j) as usize, &mut newtab_p2d);
            self.seq_poles2d.push(newtab_p2d);
        }

        self.done = true;
    }
}
