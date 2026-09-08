//! OCCT GeomFill_GuideTrihedronPlan (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_GuideTrihedronPlan.hxx (members) +
//! GeomFill_GuideTrihedronPlan.cxx (whole file L51-580, including the
//! file-local InGoodPeriod static).
//!
//! Architecture differences:
//! - OCCT declares private `myTrimmed` / `myCurve` handles that SHADOW the
//!   inherited GeomFill_TrihedronLaw members (C++ name hiding, hxx L120-121):
//!   the rcad struct carries both shadow fields and the trait base.
//! - `IntCurveSurface_HInter` over generic Curve3/Surface3 needs the
//!   production host-tool markers (see [`int_curve_surface_h_inter`]) — the
//!   carrier preserves the OCCT failure path.
//! - `Pole` is the `NCollection_HArray2<gp_Pnt2d>(1, 1, 1, myNbPts)` array;
//!   `X` / `XTol` / `Inf` / `Sup` are `math_Vector(1, 1)`.
//! - `math_FunctionRoot(F, Guess, Tol, A, B, NbIter)` is routed through
//!   [`FunctionSetRoot`] exactly as in OCCT (math_FunctionRoot.cxx L95-118
//!   wraps math_FunctionSetRoot::Perform(F, V, Aa, Bb)).
//! - The OCCT D0/D1 bodies write `X` (through InitX) and `myStatus` through
//!   `&self`; both caches are `Cell` (Rust trait methods take &self).

use std::cell::Cell;

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve3, CurveEval, Plane, Surface3, TrimmedCurve3};
use rcad_kernel::math::el::in_period;
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};
use rcad_kernel::math::GeomAbsShape;

use super::frenet::Frenet;
use super::int_curve_surface_h_inter::IntCurveSurfaceHInter;
use super::plan_func::PlanFunc;
use super::trihedron_law::{
    curve_first_parameter, curve_last_parameter, PipeError, TrihedronLaw, TrihedronLawBase,
};
use super::trihedron_with_guide::{TrihedronWithGuide, TrihedronWithGuideBase};

/// OCCT GeomFill_GuideTrihedronPlan (GeomFill_GuideTrihedronPlan.hxx
/// L115-127).
#[derive(Debug, Clone)]
pub struct GuideTrihedronPlan {
    /// OCCT GeomFill_TrihedronLaw protected base (myCurve / myTrimmed) —
    /// unused because of the private shadow members below.
    pub(crate) base: TrihedronLawBase,
    /// OCCT GeomFill_TrihedronWithGuide protected base
    /// (myGuide / myTrimG / myCurPointOnGuide).
    pub(crate) guide_base: TrihedronWithGuideBase,
    /// OCCT private handle(Adaptor3d_Curve) myTrimmed (hxx L120 — SHADOW).
    my_trimmed: Option<Curve3>,
    /// OCCT private handle(Adaptor3d_Curve) myCurve (hxx L121 — SHADOW).
    my_curve: Option<Curve3>,
    /// OCCT handle(NCollection_HArray2<gp_Pnt2d>) Pole — (1, 1, 1, myNbPts);
    /// stored as the flat row Value(1, ii).
    pole: Vec<DVec2>,
    /// OCCT int myNbPts.
    my_nb_pts: usize,
    /// OCCT handle(GeomFill_Frenet) frenet.
    frenet: Frenet,
    /// OCCT math_Vector X(1, 1) — InitX writes it from D0/D1 (&self).
    x: Cell<[f64; 1]>,
    /// OCCT math_Vector XTol(1, 1).
    x_tol: [f64; 1],
    /// OCCT math_Vector Inf(1, 1).
    inf: [f64; 1],
    /// OCCT math_Vector Sup(1, 1).
    sup: [f64; 1],
    /// OCCT GeomFill_PipeError myStatus — D0/D1 write it (&self).
    my_status: Cell<PipeError>,
}

/// OCCT static InGoodPeriod (L51-63).
fn in_good_period(prec: f64, period: f64, current: &mut f64) {
    let mut diff = *current - prec;
    let nb = (diff / period).trunc() as i32;
    *current -= nb as f64 * period;
    diff = *current - prec;
    if diff > period / 2.0 {
        *current -= period;
    } else if diff < -period / 2.0 {
        *current += period;
    }
}

/// OCCT math_FunctionRoot(F, Guess, Tolerance, A, B, NbIterations)
/// (math_FunctionRoot.cxx L95-118) — routes through math_FunctionSetRoot
/// exactly as the OCCT wrapper does.  Returns (IsDone, Root()).
fn function_root_bounded(
    f: &mut PlanFunc,
    guess: f64,
    tolerance: f64,
    a: f64,
    b: f64,
    nb_iterations: i32,
) -> (bool, f64) {
    let v = [guess];
    let tol = [tolerance];
    // OCCT: math_FunctionSetRoot Sol(Ff, Tol, NbIterations);
    // Sol.Perform(Ff, V, Aa, Bb).
    let mut sol = FunctionSetRoot::new(f, &tol, nb_iterations);
    // OCCT Perform(F, V, Aa, Bb) — no stop-on-divergent switch in the OCCT
    // signature (the rcad port exposes it; false matches the OCCT behavior).
    sol.perform(f, &v, &[a], &[b], false);
    let done = sol.is_done();
    let root = if done { sol.root()[0] } else { f64::MAX };
    (done, root)
}

impl GuideTrihedronPlan {
    /// OCCT GeomFill_GuideTrihedronPlan::GeomFill_GuideTrihedronPlan
    /// (L69-88).
    pub fn new(the_guide: &Curve3) -> Self {
        // OCCT: X(1,1), XTol(1,1), Inf(1,1), Sup(1,1), myStatus(GeomFill_PipeOk).
        let my_nb_pts = 20usize; // nb points pour calculs
        // OCCT: Pole = new NCollection_HArray2<gp_Pnt2d>(1, 1, 1, myNbPts);
        // frenet = new GeomFill_Frenet(); XTol.Init(1.e-6);
        // XTol(1) = myGuide->Resolution(1.e-6).
        let mut x_tol = [1.0e-6f64; 1];
        x_tol[0] = the_guide.resolution(1.0e-6);
        GuideTrihedronPlan {
            base: TrihedronLawBase::default(),
            guide_base: TrihedronWithGuideBase {
                my_guide: Some(the_guide.clone()), // guide
                my_trim_g: Some(the_guide.clone()),
                my_cur_point_on_guide: Cell::new(DVec3::ZERO),
            },
            my_trimmed: None,
            my_curve: None,
            pole: vec![DVec2::ZERO; my_nb_pts],
            my_nb_pts,
            frenet: Frenet::new(),
            x: Cell::new([0.0; 1]),
            x_tol,
            inf: [0.0; 1],
            sup: [0.0; 1],
            my_status: Cell::new(PipeError::PipeOk),
        }
    }

    /// OCCT Init (L92-193) — calcule myNbPts points sur la courbe guide
    /// (<=> normale).
    pub fn init(&mut self) {
        self.my_status.set(PipeError::PipeOk);
        // OCCT: f = myCurve->FirstParameter(); l = myCurve->LastParameter()
        // (the SHADOW myCurve).
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        let f = curve_first_parameter(my_curve);
        let l = curve_last_parameter(my_curve);

        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        let mut w = 0.0f64;
        let mut delta_g =
            (curve_last_parameter(my_guide) - curve_first_parameter(my_guide)) / 2.0;

        self.inf[0] = curve_first_parameter(my_guide) - delta_g;
        self.sup[0] = curve_last_parameter(my_guide) + delta_g;

        if !my_guide.is_periodic() {
            // OCCT: myTrimG = myGuide->Trim(First - DeltaG/100,
            // Last + DeltaG/100, DeltaG * 1.e-7).
            self.guide_base.my_trim_g = Some(Curve3::Trimmed(TrimmedCurve3::new(
                my_guide.clone(),
                curve_first_parameter(my_guide) - delta_g / 100.0,
                curve_last_parameter(my_guide) + delta_g / 100.0,
            )));
        } else {
            self.guide_base.my_trim_g = Some(my_guide.clone());
        }
        //  double Step = DeltaG/100;
        delta_g /= 3.0;
        for ii in 1..=self.my_nb_pts {
            let mut t = ((self.my_nb_pts - ii) as f64) * f + ((ii - 1) as f64) * l;
            t /= (self.my_nb_pts - 1) as f64;
            let p = my_curve.point_at(t);
            let mut tangent = DVec3::ZERO;
            let mut normal = DVec3::ZERO;
            let mut bi_normal = DVec3::ZERO;
            self.frenet.d0(t, &mut tangent, &mut normal, &mut bi_normal);
            // OCCT: Plan = new Geom_Plane(P, Tangent);
            // Pl = new GeomAdaptor_Surface(Plan).
            let plan = Surface3::Plane(Plane::new(p, tangent));

            let mut int = IntCurveSurfaceHInter::default();
            let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
            int.perform(my_trim_g, &plan); // intersection plan / guide
            if int.nb_points() == 0 {
                w = if (curve_last_parameter(my_guide) - w).abs()
                    > (curve_first_parameter(my_guide) - w).abs()
                {
                    curve_first_parameter(my_guide)
                } else {
                    curve_last_parameter(my_guide)
                };

                self.my_status.set(PipeError::PlaneNotIntersectGuide);
                // return;
            } else {
                let mut p_int = int.point(1);
                let mut pmin = p_int.pnt();
                let mut dmin = p.distance(pmin);
                for jj in 2..=int.nb_points() {
                    pmin = int.point(jj).pnt();
                    if p.distance(pmin) < dmin {
                        p_int = int.point(jj);
                        dmin = p.distance(pmin);
                    }
                } // for_jj

                w = p_int.w();
            }
            if ii > 1 {
                let mut diff = w - self.pole[ii - 2].y;
                if diff.abs() > delta_g {
                    if my_guide.is_periodic() {
                        // OCCT: InGoodPeriod(Pole->Value(1, ii-1).Y(),
                        // myGuide->Period(), w).
                        let period =
                            curve_last_parameter(my_guide) - curve_first_parameter(my_guide);
                        in_good_period(self.pole[ii - 2].y, period, &mut w);

                        diff = w - self.pole[ii - 2].y;
                    }
                }
            }

            // OCCT: gp_Pnt2d p1(t, w); Pole->SetValue(1, ii, p1) — on stocke
            // les parametres.
            let p1 = DVec2::new(t, w);
            self.pole[ii - 1] = p1;
        } // for_ii
    }

    /// OCCT InitX (L525-580) — recherche par interpolation d'une valeur
    /// initiale.
    fn init_x(&self, param: f64) {
        let mut ideb = 1usize;
        let mut ifin = self.pole.len();
        let mut valeur;

        valeur = self.pole[ideb - 1].x;
        if param == valeur {
            ifin = ideb + 1;
        }

        valeur = self.pole[ifin - 1].x;
        if param == valeur {
            ideb = ifin - 1;
        }

        while ideb + 1 != ifin {
            let idemi = (ideb + ifin) / 2;
            valeur = self.pole[idemi - 1].x;
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                ideb = idemi;
                ifin = ideb + 1;
            }
        }

        let t1 = self.pole[ideb - 1].x;
        let t2 = self.pole[ifin - 1].x;
        let diff = t2 - t1;
        if diff > 1.0e-7 {
            let b = (param - t1) / diff;
            let a = (t2 - param) / diff;

            // OCCT: X(1) = Pole->Value(1, Ideb).Coord(2) * a
            //            + Pole->Value(1, Ifin).Coord(2) * b; // param guide
            self.x.set([self.pole[ideb - 1].y * a + self.pole[ifin - 1].y * b]);
        } else {
            self.x.set([(self.pole[ideb - 1].y + self.pole[ifin - 1].y) / 2.0]);
        }
        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        if my_guide.is_periodic() {
            let x0 = self.x.get()[0];
            self.x.set([in_period(
                x0,
                curve_first_parameter(my_guide),
                curve_last_parameter(my_guide),
            )]);
        }
    }
}

impl TrihedronLaw for GuideTrihedronPlan {
    fn my_curve(&self) -> &Option<Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<Curve3>) {
        self.base.my_trimmed = c;
    }

    /// OCCT SetCurve (L195-204).
    fn set_curve(&mut self, c: Curve3) -> bool {
        // OCCT: myCurve = C — the SHADOW member.
        self.my_curve = Some(c);
        // OCCT: if (!myCurve.IsNull()) { Init(); }.
        self.init();
        true
    }

    /// OCCT D0 (L214-266).
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        // myCurve->D0(Param, P).
        let p = my_curve.point_at(param);

        self.frenet.d0(param, tangent, normal, binormal);

        // initialisation de la recherche
        self.init_x(param);

        let iter = 50;

        // fonction dont il faut trouver la racine : G(W)-Pl(U,V)=0
        let mut e = PlanFunc::new(p, *tangent, my_curve);

        // resolution
        let (is_done, res) =
            function_root_bounded(&mut e, self.x.get()[0], self.x_tol[0], self.inf[0], self.sup[0], iter);

        if is_done {
            // OCCT: Pprime = myTrimG->Value(Res) — pt sur courbe guide.
            let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
            let pprime = my_trim_g.point_at(res);
            // gp_Vec n(P, Pprime) — vecteur definissant la normale du triedre.
            let n = pprime - p;

            *normal = n.normalize_or_zero();
            *binormal = tangent.cross(*normal);
            *binormal = binormal.normalize_or_zero();
        } else {
            // Erreur...
            self.my_status.set(PipeError::PlaneNotIntersectGuide);
            return false;
        }

        true
    }

    /// OCCT D1 (L268-356).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
    ) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        // triedre de frenet sur la trajectoire
        // myCurve->D1(Param, P, To).
        let p = my_curve.point_at(param);
        let to = my_curve.derivative_at(param);
        self.frenet
            .d1(param, tangent, dtangent, normal, dnormal, binormal, dbinormal);

        // tolerance sur E
        let iter = 50;

        // fonction dont il faut trouver la racine : G(W)-Pl(U,V)=0
        self.init_x(param);
        let mut e = PlanFunc::new(p, *tangent, my_curve);

        // resolution
        let (is_done, res) =
            function_root_bounded(&mut e, self.x.get()[0], self.x_tol[0], self.inf[0], self.sup[0], iter);

        if is_done {
            let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
            // myTrimG->D1(Res, PG, TG).
            let pg = my_trim_g.point_at(res);
            let tg = my_trim_g.derivative_at(res);
            let mut n = pg - p;
            let mut norm = n.length();
            if norm < 1.0e-12 {
                norm = 1.0;
            }
            n /= norm;

            *normal = n;
            *binormal = tangent.cross(*normal);

            // derivee premiere du triedre
            let mut dedx = 0.0;
            let mut dedt = 0.0;
            e.derivative(res, &mut dedx);
            e.dedt(res, to, *dtangent, &mut dedt);
            let dtg_dt = -dedt / dedx;

            // OCCT: dn.SetLinearForm(dtg_dt, TG, -1, To).
            let dn = dtg_dt * tg - to;

            // OCCT: DNormal.SetLinearForm(-(n * dn), n, dn); DNormal /= Norm.
            *dnormal = -(n.dot(dn)) * n + dn;
            *dnormal /= norm;
            // OCCT: DBiNormal.SetLinearForm(Tangent.Crossed(DNormal),
            // DTangent.Crossed(Normal)).
            *dbinormal = tangent.cross(*dnormal) + dtangent.cross(*normal);
        } else {
            // Erreur...
            self.my_status.set(PipeError::PlaneNotIntersectGuide);
            return false;
        }

        true
    }

    /// OCCT D2 (L358-391) — the OCCT body computes the Frenet D2 and then
    /// returns false verbatim.
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        d2tangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        d2normal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        d2binormal: &mut DVec3,
    ) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        // myCurve->D2(Param, P, To, DTo).
        let _p = my_curve.point_at(param);
        let _to = my_curve.derivative_at(param);
        let _d_to = my_curve.derivative2_at(param);

        // triedre de Frenet sur la trajectoire
        self.frenet.d2(
            param,
            tangent,
            dtangent,
            d2tangent,
            normal,
            dnormal,
            d2normal,
            binormal,
            dbinormal,
            d2binormal,
        );

        false
    }

    /// OCCT Copy (L393-400).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let guide = self
            .guide_base
            .my_guide
            .clone()
            .expect("null myGuide in GuideTrihedronPlan::Copy");
        let mut copy = GuideTrihedronPlan::new(&guide);
        if let Some(curve) = self.my_curve.clone() {
            TrihedronLaw::set_curve(&mut copy, curve);
        }
        Box::new(copy)
    }

    /// OCCT ErrorStatus (L402-409).
    fn error_status(&self) -> PipeError {
        self.my_status.get()
    }

    /// OCCT NbIntervals (L411-434) — Version provisoire : Il faut tenir
    /// compte du guide.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let tmp_s = lifted_continuity(s);

        // OCCT: Nb = myCurve->NbIntervals(tmpS) — the SHADOW myCurve.
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        super::frenet::curve_nb_intervals(my_curve, tmp_s)
    }

    /// OCCT Intervals (L436-457).
    fn intervals(&self, tt: &mut Vec<f64>, s: GeomAbsShape) {
        let tmp_s = lifted_continuity(s);
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        let disc = super::frenet::curve_intervals(my_curve, tmp_s);
        for (i, v) in disc.iter().enumerate() {
            tt[i] = *v;
        }
    }

    /// OCCT SetInterval (L459-463).
    fn set_interval(&mut self, first: f64, last: f64) {
        // OCCT: myTrimmed = myCurve->Trim(First, Last, Precision::Confusion())
        // — the SHADOW members.
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        self.my_trimmed = Some(Curve3::Trimmed(TrimmedCurve3::new(
            my_curve.clone(),
            first,
            last,
        )));
    }

    /// OCCT GetAverageLaw (L466-491).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        let delta = (curve_last_parameter(my_curve) - curve_first_parameter(my_curve)) / 20.001;

        *atangent = DVec3::ZERO;
        *anormal = DVec3::ZERO;
        *abinormal = DVec3::ZERO;
        let mut t = DVec3::ZERO;
        let mut n = DVec3::ZERO;
        let mut b = DVec3::ZERO;

        for ii in 1..=20 {
            let tt = curve_first_parameter(my_curve) + (ii - 1) as f64 * delta;
            TrihedronLaw::d0(self, tt, &mut t, &mut n, &mut b);
            *atangent += t;
            *anormal += n;
            *abinormal += b;
        }
        *atangent /= 20.0;
        *anormal /= 20.0;
        *abinormal /= 20.0;
    }

    /// OCCT IsConstant (L493-508).
    fn is_constant(&self) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        if let (Curve3::Line(c), Curve3::Line(g)) = (my_curve, my_guide) {
            // OCCT: Angle = myCurve->Line().Angle(myGuide->Line()) — the
            // angle of the two directions in [0, Pi] (gp_Dir::Angle:
            // ATan2(|D1 x D2|, D1 . D2); for unit directions
            // acos(dot) is the same value).
            let dot = c.direction.dot(g.direction).clamp(-1.0, 1.0);
            let angle = dot.acos();
            if angle < 1.0e-12 || (2.0 * std::f64::consts::PI - angle) < 1.0e-12 {
                return true;
            }
        }

        false
    }

    /// OCCT IsOnlyBy3dCurve (L510-513).
    fn is_only_by3d_curve(&self) -> bool {
        false
    }
}

/// OCCT switch (S) { case C0: tmpS = C1; ... default: tmpS = CN } — shared
/// by NbIntervals / Intervals (pure math/tool mapping).
fn lifted_continuity(s: GeomAbsShape) -> GeomAbsShape {
    match s {
        GeomAbsShape::C0 => GeomAbsShape::C1,
        GeomAbsShape::C1 => GeomAbsShape::C2,
        GeomAbsShape::C2 => GeomAbsShape::C3,
        _ => GeomAbsShape::CN,
    }
}

impl TrihedronWithGuide for GuideTrihedronPlan {
    /// OCCT Guide (L207-211).
    fn guide(&self) -> Option<Curve3> {
        self.guide_base.my_guide.clone()
    }

    /// OCCT Origine (L519-520) — Nothing!!
    fn origine(&mut self, _param1: f64, _param2: f64) {}

    fn current_point_on_guide(&self) -> DVec3 {
        self.guide_base.my_cur_point_on_guide.get()
    }

    fn with_guide_base(&self) -> &TrihedronWithGuideBase {
        &self.guide_base
    }

    fn with_guide_base_mut(&mut self) -> &mut TrihedronWithGuideBase {
        &mut self.guide_base
    }
}
