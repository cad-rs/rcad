//! OCCT IntImpParGen package (TKGeomAlgo) — generic algorithm intersecting an
//! implicit curve with a bounded parametric curve.
//!
//! 1:1 translations:
//! - [`Intersector`] — `IntImpParGen_Intersector.gxx` (L1-824): the walking
//!   intersection engine (point solutions + segment solutions + domain
//!   clipping + transitions + endpoint tests).  The OCCT template parameters
//!   map to the [`ParTool`]/[`ProjectOnPCurveTool`] traits (the ImpTool
//!   parameter is [`IConicTool`] in every OCCT instantiation):
//!   - Geom2dInt_TheIntersectorOfTheIntConicCurveOfGInter —
//!     ParTool = Geom2dInt_Geom2dCurveTool, ProjectOnPCurveTool =
//!     Geom2dInt_TheProjPCurOfGInter, ParCurve = Adaptor2d_Curve2d;
//!   - HLRBRep_TheIntersectorOfTheIntConicCurveOfCInter —
//!     ParTool = HLRBRep_CurveTool, ProjectOnPCurveTool =
//!     HLRBRep_TheProjPCurOfCInter, ParCurve = HLRBRep_CurvePtr (Stage 3a).
//! - The `IntImpParGen.cxx` package statics (NormalizeOnDomain /
//!   DeterminePosition / DetermineTransition x2).
//! - `IntImpParGen_MyImpParTool` — the signed-distance function
//!   F(u) = Dist(ImpCurve, P(u)) consumed by math_FunctionAllRoots; the
//!   per-instantiation bodies (Geom2dInt_MyImpParToolOfTheIntersectorOf…,
//!   HLRBRep_MyImpParToolOfTheIntersectorOf…) are textually identical
//!   modulo the tool types, hence one generic [`MyImpParTool`].

use glam::DVec2;
use rcad_kernel::math::root::{FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative};

use super::geom2d_int::IConicTool;
use super::int_res2d::{
    Domain as Res2dDomain, IntersectionBase, IntersectionPoint, IntersectionSegment, Position,
    Situation, Transition, TypeTrans,
};

/// gp::Resolution() — OCCT gp.hxx (compared against squared magnitudes).
const GP_RESOLUTION: f64 = 1e-15;

/// OCCT `ParTool` template parameter of IntImpParGen_Intersector.gxx — the
/// static accessors over the parametric curve (Geom2dInt_Geom2dCurveTool for
/// the GInter instantiation, HLRBRep_CurveTool for the CInter one).  The
/// curve type is a trait parameter (`ParTool<C>`), mirroring the OCCT
/// template signature where the ParCurve type is fixed by the
/// instantiation; `C = dyn Curve2dAdaptor + 'a` keeps reference signatures
/// lifetime-elastic like the C++ `const Adaptor2d_Curve2d&`.
pub trait ParTool<C: ?Sized> {
    /// OCCT ParTool::NbSamples(C, U1, U2).
    fn nb_samples_uv(c: &C, u1: f64, u2: f64) -> i32;
    /// OCCT ParTool::EpsX(C).
    fn eps_x(c: &C) -> f64;
    /// OCCT ParTool::Value(C, U).
    fn value(c: &C, u: f64) -> DVec2;
    /// OCCT ParTool::D1(C, U, P, T).
    fn d1(c: &C, u: f64) -> (DVec2, DVec2);
    /// OCCT ParTool::D2(C, U, P, T, N).
    fn d2(c: &C, u: f64) -> (DVec2, DVec2, DVec2);
}

/// OCCT `ProjectOnPCurveTool` template parameter of IntImpParGen_Intersector.gxx —
/// the point projection on the parametric curve (Geom2dInt_TheProjPCurOfGInter
/// for the GInter instantiation, HLRBRep_TheProjPCurOfCInter for the CInter one).
pub trait ProjectOnPCurveTool<C: ?Sized> {
    /// OCCT ProjectOnPCurveTool::FindParameter(C, P, Tol).
    fn find_parameter(c: &C, p: DVec2, tol: f64) -> f64;
    /// OCCT ProjectOnPCurveTool::FindParameter(C, P, Low, High, Tol).
    fn find_parameter_between(c: &C, p: DVec2, low: f64, high: f64, tol: f64) -> f64;
}

/// OCCT IntImpParGen::NormalizeOnDomain (IntImpParGen.cxx L28-46).
pub fn normalize_on_domain(param: f64, domain: &Res2dDomain) -> f64 {
    let mut mod_param = param;
    if domain.is_closed() {
        let (t, mut periode) = domain.equivalent_parameters();
        periode -= t;
        while mod_param < domain.first_parameter() && mod_param + periode < domain.last_parameter() {
            mod_param += periode;
        }
        while mod_param > domain.last_parameter() && mod_param - periode > domain.first_parameter() {
            mod_param -= periode;
        }
    }
    mod_param
}

/// OCCT IntImpParGen::DeterminePosition (IntImpParGen.cxx L49-83).
pub fn determine_position(pos: &mut Position, domain: &Res2dDomain, pnt: DVec2, param: f64) {
    *pos = Position::Middle;
    if domain.has_first_point() {
        if pnt.distance(domain.first_point()) <= domain.first_tolerance() {
            *pos = Position::Head;
        }
    }
    if domain.has_last_point() {
        if pnt.distance(domain.last_point()) <= domain.last_tolerance() {
            if *pos == Position::Head {
                if (param - domain.last_parameter()).abs() < (param - domain.first_parameter()).abs() {
                    *pos = Position::End;
                }
            } else {
                *pos = Position::End;
            }
        }
    }
}

/// OCCT IntImpParGen::DetermineTransition (IntImpParGen.cxx L86-206) — the
/// full overload with the second derivatives (TOUCH classification).
#[allow(clippy::too_many_arguments)]
pub fn determine_transition(
    pos1: Position,
    tan1: &mut DVec2,
    norm1: DVec2,
    t1: &mut Transition,
    pos2: Position,
    tan2: &mut DVec2,
    norm2: DVec2,
    t2: &mut Transition,
    _tol: f64,
) {
    let mut courbure1 = true;
    let mut courbure2 = true;
    let mut decide = true;

    t1.set_position(pos1);
    t2.set_position(pos2);

    if tan1.length_squared() <= DERIVEE_PREMIERE_NULLE {
        *tan1 = norm1;
        courbure1 = false;
        if tan1.length_squared() <= DERIVEE_PREMIERE_NULLE {
            decide = false;
        }
    }
    if tan2.length_squared() <= DERIVEE_PREMIERE_NULLE {
        *tan2 = norm2;
        courbure2 = false;
        if tan2.length_squared() <= DERIVEE_PREMIERE_NULLE {
            decide = false;
        }
    }

    if !decide {
        t1.set_value_undecided(pos1);
        t2.set_value_undecided(pos2);
    } else {
        let sgn = tan1.x * tan2.y - tan1.y * tan2.x;
        let norm = tan1.length() * tan2.length();
        if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
            // Transition TOUCH.
            let opos = tan1.dot(*tan2) < 0.0;
            if !(courbure1 || courbure2) {
                t1.set_value_touch(true, pos1, Situation::Unknown, opos);
                t2.set_value_touch(true, pos2, Situation::Unknown, opos);
            } else {
                let norm_v = DVec2::new(-tan1.y, tan1.x);
                let val1 = if !courbure1 { 0.0 } else { norm_v.dot(norm1) };
                let val2 = if !courbure2 { 0.0 } else { norm_v.dot(norm2) };
                if (val1 - val2).abs() <= TOLERANCE_ANGULAIRE {
                    t1.set_value_touch(true, pos1, Situation::Unknown, opos);
                    t2.set_value_touch(true, pos2, Situation::Unknown, opos);
                } else if val2 > val1 {
                    t2.set_value_touch(true, pos2, Situation::Inside, opos);
                    if opos {
                        t1.set_value_touch(true, pos1, Situation::Inside, opos);
                    } else {
                        t1.set_value_touch(true, pos1, Situation::Outside, opos);
                    }
                } else {
                    // val1 > val2
                    t2.set_value_touch(true, pos2, Situation::Outside, opos);
                    if opos {
                        t1.set_value_touch(true, pos1, Situation::Outside, opos);
                    } else {
                        t1.set_value_touch(true, pos1, Situation::Inside, opos);
                    }
                }
            }
        } else if sgn < 0.0 {
            t1.set_value_in_out(false, pos1, TypeTrans::In);
            t2.set_value_in_out(false, pos2, TypeTrans::Out);
        } else {
            // sgn > 0
            t1.set_value_in_out(false, pos1, TypeTrans::Out);
            t2.set_value_in_out(false, pos2, TypeTrans::In);
        }
    }
}

/// OCCT IntImpParGen::DetermineTransition (IntImpParGen.cxx L209-251) — the
/// IN/OUT-only overload (returns false when the transition cannot be decided).
#[allow(clippy::too_many_arguments)]
pub fn determine_transition_in_out(
    pos1: Position,
    tan1: &mut DVec2,
    t1: &mut Transition,
    pos2: Position,
    tan2: &mut DVec2,
    t2: &mut Transition,
    _tol: f64,
) -> bool {
    t1.set_position(pos1);
    t2.set_position(pos2);

    let tan1_mag = tan1.length();
    if tan1_mag <= DERIVEE_PREMIERE_NULLE {
        return false;
    }
    let tan2_mag = tan2.length();
    if tan2_mag <= DERIVEE_PREMIERE_NULLE {
        return false;
    }

    let sgn = tan1.x * tan2.y - tan1.y * tan2.x;
    let norm = tan1_mag * tan2_mag;
    if sgn.abs() <= TOLERANCE_ANGULAIRE * norm {
        return false;
    } else if sgn < 0.0 {
        t1.set_value_in_out(false, pos1, TypeTrans::In);
        t2.set_value_in_out(false, pos2, TypeTrans::Out);
    } else {
        t1.set_value_in_out(false, pos1, TypeTrans::Out);
        t2.set_value_in_out(false, pos2, TypeTrans::In);
    }
    true
}

/// OCCT IntImpParGen.cxx L24: TOLERANCE_ANGULAIRE.
const TOLERANCE_ANGULAIRE: f64 = 0.00000001;
/// OCCT IntImpParGen.cxx L25: DERIVEE_PREMIERE_NULLE.
const DERIVEE_PREMIERE_NULLE: f64 = 0.000000000001;

/// OCCT Geom2dInt_MyImpParToolOfTheIntersectorOfTheIntConicCurveOfGInter /
/// HLRBRep_MyImpParToolOfTheIntersectorOfTheIntConicCurveOfCInter — the
/// signed-distance function F(u) = Dist(ImpCurve, P(u)) with derivative
/// F'(u) = GradDist(P(u)).P'(u).  The per-instantiation bodies are
/// textually identical modulo the tool types.
pub struct MyImpParTool<'a, C: ?Sized, PT: ParTool<C>> {
    imp_tool: &'a IConicTool,
    par_curve: &'a C,
    _tools: std::marker::PhantomData<fn(&PT)>,
}

impl<'a, C: ?Sized, PT: ParTool<C>> MyImpParTool<'a, C, PT> {
    pub fn new(imp_tool: &'a IConicTool, par_curve: &'a C) -> Self {
        MyImpParTool {
            imp_tool,
            par_curve,
            _tools: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, PT: ParTool<C>> FunctionValue for MyImpParTool<'_, C, PT> {
    fn value(&mut self, param: f64) -> Option<f64> {
        Some(self.imp_tool.distance(PT::value(self.par_curve, param)))
    }
}

impl<C: ?Sized, PT: ParTool<C>> FunctionWithDerivative for MyImpParTool<'_, C, PT> {
    fn derivative(&mut self, param: f64) -> Option<f64> {
        let pt = PT::value(self.par_curve, param);
        let grad = self.imp_tool.grad_distance(pt);
        let (_, tan) = PT::d1(self.par_curve, param);
        Some(grad.dot(tan))
    }
    fn values(&mut self, param: f64) -> Option<(f64, f64)> {
        let v = self.value(param)?;
        let d = self.derivative(param)?;
        Some((v, d))
    }
}

/// OCCT IntImpParGen_Intersector (IntImpParGen_Intersector.gxx) — the
/// walking intersection of an implicit conic with a parametric curve,
/// generic over the template parameters `ParTool` and
/// `ProjectOnPCurveTool` (the ImpTool parameter is `IntCurve_IConicTool`
/// in every OCCT instantiation); the parametric-curve type is the `C`
/// parameter shared by both tools.
pub struct Intersector<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> {
    pub base: IntersectionBase,
    _tools: std::marker::PhantomData<fn(&C, &PT, &JT)>,
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> Clone for Intersector<C, PT, JT> {
    fn clone(&self) -> Self {
        Intersector {
            base: self.base.clone(),
            _tools: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> std::fmt::Debug
    for Intersector<C, PT, JT>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Intersector").field("base", &self.base).finish()
    }
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> Intersector<C, PT, JT> {
    /// OCCT IntImpParGen_Intersector::IntImpParGen_Intersector() (L226-229).
    pub fn new() -> Self {
        Intersector {
            base: IntersectionBase::new(),
            _tools: std::marker::PhantomData,
        }
    }

    /// OCCT IntImpParGen_Intersector::FindU (L781-788).
    fn find_u(&self, parameter: f64, par_curve: &C, imp_tool: &IConicTool) -> (DVec2, f64) {
        let point = PT::value(par_curve, parameter);
        (point, imp_tool.find_parameter(point))
    }

    /// OCCT IntImpParGen_Intersector::FindV (L790-824).
    #[allow(clippy::too_many_arguments)]
    fn find_v(
        &self,
        parameter: f64,
        imp_tool: &IConicTool,
        par_curve: &C,
        par_domain: &Res2dDomain,
        v0: f64,
        v1: f64,
        tolerance: f64,
    ) -> f64 {
        let point = imp_tool.value(parameter);
        if par_domain.is_closed() {
            let v = JT::find_parameter(par_curve, point, tolerance);
            normalize_on_domain(v, par_domain)
        } else {
            let (mut vv0, mut vv1) = (v0, v1);
            if v1 < v0 {
                vv0 = v1;
                vv1 = v0;
            }
            // Modif le 15 Septembre 1992: the returned parameter is tested
            // against the bracket (gxx L814-822).
            let mut x = JT::find_parameter_between(par_curve, point, vv0, vv1, tolerance);
            if x > vv1 {
                x = vv1;
            } else if x < vv0 {
                x = vv0;
            }
            x
        }
    }

    /// OCCT IntImpParGen_Intersector::And_Domaine_Objet1_Intersections
    /// (L42-222).
    #[allow(clippy::too_many_arguments)]
    fn and_domaine_objet1_intersections(
        &self,
        imp_tool: &IConicTool,
        par_curve: &C,
        imp_domain: &Res2dDomain,
        par_domain: &Res2dDomain,
        nb_resultats: &mut usize,
        inter2_and_domain2: &[f64],
        inter1: &[f64],
        resultat1: &mut [f64],
        resultat2: &mut [f64],
        eps_nul: f64,
    ) {
        let nb_bornes_intersection = *nb_resultats;
        *nb_resultats = 0;

        let mut i = 0usize;
        while i < nb_bornes_intersection {
            let mut param1 = inter1[i];
            let mut param2 = inter1[i + 1];
            let (mut indice_1, mut indice_2) = (i, i + 1);
            if param1 > param2 {
                let t = param1;
                param1 = param2;
                param2 = t;
                indice_1 = i + 1;
                indice_2 = i;
            }

            let pt1 = imp_tool.value(param1);
            let pt2 = imp_tool.value(param2);

            let mut is_on_the_imp_curve_domain1 = true;
            let mut is_on_the_imp_curve_domain2 = true;
            if imp_domain.has_first_point() {
                if param1 < imp_domain.first_parameter() {
                    if pt1.distance(imp_domain.first_point()) > imp_domain.first_tolerance() {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain1 && imp_domain.has_last_point() {
                if param1 > imp_domain.last_parameter() {
                    if pt1.distance(imp_domain.last_point()) > imp_domain.last_tolerance() {
                        is_on_the_imp_curve_domain1 = false;
                    }
                }
            }
            if imp_domain.has_first_point() {
                if param2 < imp_domain.first_parameter() {
                    if pt2.distance(imp_domain.first_point()) > imp_domain.first_tolerance() {
                        is_on_the_imp_curve_domain2 = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain2 && imp_domain.has_last_point() {
                if param2 > imp_domain.last_parameter() {
                    if pt2.distance(imp_domain.last_point()) > imp_domain.last_tolerance() {
                        is_on_the_imp_curve_domain2 = false;
                    }
                }
            }

            if is_on_the_imp_curve_domain1 {
                // Bound 1 is on the domain.
                *nb_resultats += 1;
                resultat1[*nb_resultats - 1] = inter1[indice_1];
                resultat2[*nb_resultats - 1] = inter2_and_domain2[indice_1];
                // Bound 2 is also on the domain.
                if is_on_the_imp_curve_domain2 {
                    *nb_resultats += 1;
                    resultat1[*nb_resultats - 1] = inter1[indice_2];
                    resultat2[*nb_resultats - 1] = inter2_and_domain2[indice_2];
                } else {
                    // Bound 1 on the domain and bound 2 outside.
                    let t = imp_domain.last_parameter();
                    *nb_resultats += 1;
                    resultat1[*nb_resultats - 1] = t;
                    resultat2[*nb_resultats - 1] = self.find_v(
                        t,
                        imp_tool,
                        par_curve,
                        par_domain,
                        inter2_and_domain2[indice_1],
                        inter2_and_domain2[indice_2],
                        eps_nul,
                    );
                }
            } else if is_on_the_imp_curve_domain2 {
                // Bound 1 is not on the domain.
                let t = imp_domain.first_parameter();
                *nb_resultats += 1;
                resultat1[*nb_resultats - 1] = t;
                resultat2[*nb_resultats - 1] = self.find_v(
                    t,
                    imp_tool,
                    par_curve,
                    par_domain,
                    inter2_and_domain2[indice_1],
                    inter2_and_domain2[indice_2],
                    eps_nul,
                );
                *nb_resultats += 1;
                resultat1[*nb_resultats - 1] = inter1[indice_2];
                resultat2[*nb_resultats - 1] = inter2_and_domain2[indice_2];
            } else if param1 < imp_domain.first_parameter() && param2 > imp_domain.last_parameter() {
                // Both bounds are outside the domain.
                let t = imp_domain.first_parameter();
                *nb_resultats += 1;
                resultat1[*nb_resultats - 1] = t;
                resultat2[*nb_resultats - 1] = self.find_v(
                    t,
                    imp_tool,
                    par_curve,
                    par_domain,
                    inter2_and_domain2[indice_1],
                    inter2_and_domain2[indice_2],
                    eps_nul,
                );
                let t = imp_domain.last_parameter();
                *nb_resultats += 1;
                resultat1[*nb_resultats - 1] = t;
                resultat2[*nb_resultats - 1] = self.find_v(
                    t,
                    imp_tool,
                    par_curve,
                    par_domain,
                    inter2_and_domain2[indice_1],
                    inter2_and_domain2[indice_2],
                    eps_nul,
                );
            }
            i += 2;
        }
    }

    /// OCCT IntImpParGen_Intersector::Perform (L245-779).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        imp_tool: &IConicTool,
        imp_domain: &Res2dDomain,
        par_curve: &C,
        par_domain: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut head_on_imp = false;
        let mut head_on_par = false;
        let mut end_on_imp = false;
        let mut end_on_par = false;

        self.base.reset_fields();

        let mut imp_par_tool = MyImpParTool::<C, PT>::new(imp_tool, par_curve);

        if !(par_domain.has_first_point() && par_domain.has_last_point()) {
            panic!("Standard_ConstructionError: Domaine sur courbe incorrect");
        }

        let nb_echantillons = PT::nb_samples_uv(
            par_curve,
            par_domain.first_parameter(),
            par_domain.last_parameter(),
        );

        let mut eps_x = PT::eps_x(par_curve);
        if eps_x > 1.0e-10 {
            eps_x = 1.0e-10;
        }
        let eps_nul = if tol_conf <= 1.0e-10 { 1.0e-10 } else { tol_conf };
        let eps_dist = if tol <= 1.0e-10 { 1.0e-10 } else { tol };

        let tolerance_angulaire = eps_dist;

        if (par_domain.last_parameter() - par_domain.first_parameter()) < 100.0 * eps_x {
            eps_x = (par_domain.last_parameter() - par_domain.first_parameter()) * 0.01;
        }

        let sample2 = FunctionSample::new(
            par_domain.first_parameter(),
            par_domain.last_parameter(),
            nb_echantillons,
        );
        let mut sol = FunctionAllRoots::new(&mut imp_par_tool, &sample2, eps_x, eps_dist, eps_nul);

        if !sol.is_done() {
            self.base.done = false;
            return;
        }

        let nb_segments_solution = sol.nb_intervals();
        let nb_points_solution = sol.nb_points();

        // ---- Treatment of the point solutions (L313-376) ----
        for i in 1..=nb_points_solution {
            let param2 = sol.get_point(i);
            let (pt, mut param1) = self.find_u(param2, par_curve, imp_tool);

            if imp_domain.is_closed() {
                param1 = normalize_on_domain(param1, imp_domain);
            }

            let mut is_on_the_imp_curve_domain = true;
            if imp_domain.has_first_point() {
                if param1 < imp_domain.first_parameter() {
                    if pt.distance(imp_domain.first_point()) > imp_domain.first_tolerance() {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }
            if is_on_the_imp_curve_domain && imp_domain.has_last_point() {
                if param1 > imp_domain.last_parameter() {
                    if pt.distance(imp_domain.last_point()) > imp_domain.last_tolerance() {
                        is_on_the_imp_curve_domain = false;
                    }
                }
            }

            if is_on_the_imp_curve_domain {
                let (pt1, mut tan1, norm1) = imp_tool.d2(param1);
                let (pt2, mut tan2, norm2) = PT::d2(par_curve, param2);

                let mut pos1 = Position::Middle;
                let mut pos2 = Position::Middle;
                determine_position(&mut pos1, imp_domain, pt1, param1);
                determine_position(&mut pos2, par_domain, pt2, param2);

                if pos1 == Position::End {
                    end_on_imp = true;
                } else if pos1 == Position::Head {
                    head_on_imp = true;
                }
                if pos2 == Position::End {
                    end_on_par = true;
                } else if pos2 == Position::Head {
                    head_on_par = true;
                }

                let mut trans1 = Transition::empty();
                let mut trans2 = Transition::empty();
                determine_transition(
                    pos1,
                    &mut tan1,
                    norm1,
                    &mut trans1,
                    pos2,
                    &mut tan2,
                    norm2,
                    &mut trans2,
                    tolerance_angulaire,
                );

                let ip = IntersectionPoint::new(
                    pt1,
                    param1,
                    param2,
                    trans1,
                    trans2,
                    self.base.reversed_parameters(),
                );
                self.base.insert(&ip);
            }
        }
        // ---- End of the treatment of the point solutions ----

        // ---- Treatment of the segments (L379-491) ----
        // The worst case is 8 created segments per solution segment (the
        // closed-implicit-curve shift creates at most 2 extra pairs), hence
        // the OCCT array sizing 2 + 8*nb_segments_solution.
        let mut inter2_and_domaine2: Vec<f64> = vec![0.0; 2 + 8 * nb_segments_solution];
        let mut inter1: Vec<f64> = vec![0.0; 2 + 8 * nb_segments_solution];
        let mut nb_segments_crees = 0usize;

        let mut j2 = 0usize;
        for j in 1..=nb_segments_solution {
            let (param2_inf, param2_sup) = sol.get_interval(j);
            let (_, mut param1_inf) = self.find_u(param2_inf, par_curve, imp_tool);
            let (_, mut param1_sup) = self.find_u(param2_sup, par_curve, imp_tool);

            // ---- Closed implicit curve ----
            if imp_domain.is_closed() {
                let (param1_origine, param1_fin) = imp_domain.equivalent_parameters();
                let periode = param1_fin - param1_origine;

                while param1_inf < param1_origine {
                    param1_inf += periode;
                }
                while param1_sup < param1_origine {
                    param1_sup += periode;
                }

                let (_, mut t2, n2) = PT::d2(par_curve, param2_inf);
                let (_, mut t1, n1) = imp_tool.d2(param1_inf);
                if t1.length_squared() <= GP_RESOLUTION {
                    t1 = n1;
                }
                if t2.length_squared() <= GP_RESOLUTION {
                    t2 = n2;
                }

                if t1.dot(t2) >= 0.0 {
                    // param1_inf designates an entering point (and T1 points
                    // towards the matter).
                    if param1_inf >= param1_sup {
                        param1_sup += periode;
                    }
                } else if param1_inf <= param1_sup {
                    // param1_inf: outgoing point (and T1 points outside).
                    param1_inf += periode;
                }

                // Create a new segment shifted by the period (gxx L445-455).
                let decal1 = if param1_inf > param1_sup {
                    param1_sup + periode
                } else {
                    param1_inf + periode
                };
                if imp_domain.last_parameter() > decal1 {
                    inter2_and_domaine2[j2] = param2_inf;
                    inter1[j2] = param1_inf + periode;
                    inter2_and_domaine2[j2 + 1] = param2_sup;
                    inter1[j2 + 1] = param1_sup + periode;
                    j2 += 2;
                    nb_segments_crees += 1;
                }

                let decal2 = if param1_inf < param1_sup {
                    param1_sup - periode
                } else {
                    param1_inf - periode
                };
                if imp_domain.first_parameter() < decal2 {
                    inter2_and_domaine2[j2] = param2_inf;
                    inter1[j2] = param1_inf - periode;
                    inter2_and_domaine2[j2 + 1] = param2_sup;
                    inter1[j2 + 1] = param1_sup - periode;
                    j2 += 2;
                    nb_segments_crees += 1;
                }
            }

            inter2_and_domaine2[j2] = param2_inf;
            inter1[j2] = param1_inf;
            inter2_and_domaine2[j2 + 1] = param2_sup;
            inter1[j2 + 1] = param1_sup;
        }

        // INTER2_DOMAINE2 : intersection AND curve domain as a function of
        // PARAM2; INTER1 : intersection AND curve domain as a function of
        // PARAM1 (L493-512).
        let nb_segments_solution_total = nb_segments_solution + nb_segments_crees;
        let mut resultat1: Vec<f64> = vec![0.0; 2 + (1 + nb_segments_solution_total) * 2];
        let mut resultat2: Vec<f64> = vec![0.0; 2 + (1 + nb_segments_solution_total) * 2];
        let mut nb_resultats = nb_segments_solution_total * 2;

        self.and_domaine_objet1_intersections(
            imp_tool,
            par_curve,
            imp_domain,
            par_domain,
            &mut nb_resultats,
            &inter2_and_domaine2,
            &inter1,
            &mut resultat1,
            &mut resultat2,
            eps_nul,
        );

        // Inlined Calcule_Toutes_Transitions (L514-658).
        {
            let dist_mini_imp_curve = eps_nul;
            let tolerance_angulaire_dist_mini = dist_mini_imp_curve;

            let mut k = 0usize;
            while k < nb_resultats {
                let ip1 = k + 1;
                let mut only_one_point = false;

                let mut param1_on1 = resultat1[k];
                let mut param1_on2 = resultat2[k];
                let mut param2_on1 = resultat1[ip1];
                let mut param2_on2 = resultat2[ip1];

                let pt1_on1 = imp_tool.value(param1_on1);
                let pt2_on1 = imp_tool.value(param2_on1);
                let pt1_on2 = PT::value(par_curve, param1_on2);
                let pt2_on2 = PT::value(par_curve, param2_on2);

                if !imp_domain.is_closed() {
                    if pt1_on1.distance(pt2_on1) <= dist_mini_imp_curve {
                        if pt1_on2.distance(pt2_on2) <= dist_mini_imp_curve {
                            only_one_point = true;
                        }
                    }
                }

                param1_on1 = normalize_on_domain(param1_on1, imp_domain);
                param1_on2 = normalize_on_domain(param1_on2, par_domain);

                let (mut pt1_on1_2, mut tan1, norm1) = imp_tool.d2(param1_on1);
                let (mut pt1_on2_2, mut tan2, norm2) = PT::d2(par_curve, param1_on2);

                let mut pos1 = Position::Middle;
                let mut pos2 = Position::Middle;
                determine_position(&mut pos1, imp_domain, pt1_on1_2, param1_on1);
                determine_position(&mut pos2, par_domain, pt1_on2_2, param1_on2);

                if pos1 == Position::End {
                    end_on_imp = true;
                } else if pos1 == Position::Head {
                    head_on_imp = true;
                }
                if pos2 == Position::End {
                    end_on_par = true;
                } else if pos2 == Position::Head {
                    head_on_par = true;
                }

                let mut trans1 = Transition::empty();
                let mut trans2 = Transition::empty();
                determine_transition(
                    pos1,
                    &mut tan1,
                    norm1,
                    &mut trans1,
                    pos2,
                    &mut tan2,
                    norm2,
                    &mut trans2,
                    tolerance_angulaire_dist_mini,
                );

                // Detection of the case: the intersection is at the end of
                // both domains.
                if pos1 != Position::Middle && pos2 != Position::Middle {
                    let m = 0.5 * (pt1_on1_2.x + pt1_on2_2.x);
                    pt1_on1_2.x = m;
                    let m = 0.5 * (pt1_on1_2.y + pt1_on2_2.y);
                    pt1_on1_2.y = m;
                }

                let new_p1 = IntersectionPoint::new(
                    pt1_on1_2,
                    param1_on1,
                    param1_on2,
                    trans1,
                    trans2,
                    self.base.reversed_parameters(),
                );
                if !only_one_point {
                    let mut new_p2 = IntersectionPoint::empty();

                    param2_on1 = normalize_on_domain(param2_on1, imp_domain);
                    param2_on2 = normalize_on_domain(param2_on2, par_domain);

                    let (mut pt2_on1_2, mut tan1b, norm1b) = imp_tool.d2(param2_on1);
                    let (mut pt2_on2_2, mut tan2b, norm2b) = PT::d2(par_curve, param2_on2);

                    let mut pos1b = Position::Middle;
                    let mut pos2b = Position::Middle;
                    determine_position(&mut pos1b, imp_domain, pt2_on1_2, param2_on1);
                    determine_position(&mut pos2b, par_domain, pt2_on2_2, param2_on2);

                    if pos1b == Position::End {
                        end_on_imp = true;
                    } else if pos1b == Position::Head {
                        head_on_imp = true;
                    }
                    if pos2b == Position::End {
                        end_on_par = true;
                    } else if pos2b == Position::Head {
                        head_on_par = true;
                    }

                    let mut trans1b = Transition::empty();
                    let mut trans2b = Transition::empty();
                    determine_transition(
                        pos1b,
                        &mut tan1b,
                        norm1b,
                        &mut trans1b,
                        pos2b,
                        &mut tan2b,
                        norm2b,
                        &mut trans2b,
                        tolerance_angulaire_dist_mini,
                    );

                    // Detection of the case: the intersection is at the end
                    // of both domains.
                    if pos1b != Position::Middle && pos2b != Position::Middle {
                        let m = 0.5 * (pt2_on1_2.x + pt2_on2_2.x);
                        pt2_on1_2.x = m;
                        let m = 0.5 * (pt2_on1_2.y + pt2_on2_2.y);
                        pt2_on1_2.y = m;
                    }

                    new_p2.set_values(
                        pt2_on1_2,
                        param2_on1,
                        param2_on2,
                        trans1b,
                        trans2b,
                        self.base.reversed_parameters(),
                    );

                    let segopposite = tan1b.dot(tan2b) < 0.0;

                    let new_seg = IntersectionSegment::with_points(
                        &new_p1,
                        &new_p2,
                        segopposite,
                        self.base.reversed_parameters(),
                    );
                    self.base.append_segment(&new_seg);
                } else {
                    self.base.insert(&new_p1);
                }
                k += 2;
            }
        }

        // ---- The boundary points are tested as solutions (L660-777) ----
        if !head_on_imp && imp_domain.has_first_point() {
            if !head_on_par {
                if imp_domain.first_point().distance(par_domain.first_point())
                    <= imp_domain.first_tolerance().max(par_domain.first_tolerance())
                {
                    let param1 = imp_domain.first_parameter();
                    let param2 = par_domain.first_parameter();
                    let (pt1, mut tan1, norm1) = imp_tool.d2(param1);
                    let (pt2, mut tan2, norm2) = PT::d2(par_curve, param2);
                    let mut trans1 = Transition::empty();
                    let mut trans2 = Transition::empty();
                    determine_transition(
                        Position::Head,
                        &mut tan1,
                        norm1,
                        &mut trans1,
                        Position::Head,
                        &mut tan2,
                        norm2,
                        &mut trans2,
                        tolerance_angulaire,
                    );
                    let ip = IntersectionPoint::new(
                        imp_domain.first_point(),
                        param1,
                        param2,
                        trans1,
                        trans2,
                        self.base.reversed_parameters(),
                    );
                    let _ = pt1;
                    let _ = pt2;
                    self.base.insert(&ip);
                }
            }
            if !end_on_par {
                if imp_domain.first_point().distance(par_domain.last_point())
                    <= imp_domain.first_tolerance().max(par_domain.last_tolerance())
                {
                    let param1 = imp_domain.first_parameter();
                    let param2 = par_domain.last_parameter();
                    let (_, mut tan1, norm1) = imp_tool.d2(param1);
                    let (_, mut tan2, norm2) = PT::d2(par_curve, param2);
                    let mut trans1 = Transition::empty();
                    let mut trans2 = Transition::empty();
                    determine_transition(
                        Position::Head,
                        &mut tan1,
                        norm1,
                        &mut trans1,
                        Position::End,
                        &mut tan2,
                        norm2,
                        &mut trans2,
                        tolerance_angulaire,
                    );
                    let ip = IntersectionPoint::new(
                        imp_domain.first_point(),
                        param1,
                        param2,
                        trans1,
                        trans2,
                        self.base.reversed_parameters(),
                    );
                    self.base.insert(&ip);
                }
            }
        }

        if !end_on_imp && imp_domain.has_last_point() {
            if !head_on_par {
                if imp_domain.last_point().distance(par_domain.first_point())
                    <= imp_domain.last_tolerance().max(par_domain.first_tolerance())
                {
                    let param1 = imp_domain.last_parameter();
                    let param2 = par_domain.first_parameter();
                    let (_, mut tan1, norm1) = imp_tool.d2(param1);
                    let (_, mut tan2, norm2) = PT::d2(par_curve, param2);
                    let mut trans1 = Transition::empty();
                    let mut trans2 = Transition::empty();
                    determine_transition(
                        Position::End,
                        &mut tan1,
                        norm1,
                        &mut trans1,
                        Position::Head,
                        &mut tan2,
                        norm2,
                        &mut trans2,
                        tolerance_angulaire,
                    );
                    let ip = IntersectionPoint::new(
                        imp_domain.last_point(),
                        param1,
                        param2,
                        trans1,
                        trans2,
                        self.base.reversed_parameters(),
                    );
                    self.base.insert(&ip);
                }
            }
            if !end_on_par {
                if imp_domain.last_point().distance(par_domain.last_point())
                    <= imp_domain.last_tolerance().max(par_domain.last_tolerance())
                {
                    let param1 = imp_domain.last_parameter();
                    let param2 = par_domain.last_parameter();
                    let (_, mut tan1, norm1) = imp_tool.d2(param1);
                    let (_, mut tan2, norm2) = PT::d2(par_curve, param2);
                    let mut trans1 = Transition::empty();
                    let mut trans2 = Transition::empty();
                    determine_transition(
                        Position::End,
                        &mut tan1,
                        norm1,
                        &mut trans1,
                        Position::End,
                        &mut tan2,
                        norm2,
                        &mut trans2,
                        tolerance_angulaire,
                    );
                    let ip = IntersectionPoint::new(
                        imp_domain.last_point(),
                        param1,
                        param2,
                        trans1,
                        trans2,
                        self.base.reversed_parameters(),
                    );
                    self.base.insert(&ip);
                }
            }
        }
        self.base.done = true;
    }
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> Default for Intersector<C, PT, JT> {
    fn default() -> Self {
        Intersector::new()
    }
}
