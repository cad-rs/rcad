//! OCCT GeomFill_GuideTrihedronAC (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_GuideTrihedronAC.hxx (members) + GeomFill_GuideTrihedronAC.cxx
//! (whole file L35-386).
//!
//! Architecture differences:
//! - OCCT declares a private `myCurve` that SHADOWS the inherited
//!   GeomFill_TrihedronLaw::myCurve (C++ name hiding): SetCurve assigns the
//!   shadow while the base member stays null.  The rcad struct carries both.
//! - `Approx_CurvlinFunc` (TKGeomAlgo/Approx) is not translated — the
//!   [`ApproxCurvlinFunc`] carrier preserves the OCCT failure path.
//! - OCCT `myTrimmed->D1/D2/D3` evaluations map to the rcad `CurveEval`
//!   trait (same precedent as the Frenet law's `law_d2` helper).
//! - `myCurPointOnGuide` is a D0/D1/D2 output cache; the OCCT write through
//!   `&self` maps to a `Cell` (Rust trait methods take &self).

use std::cell::Cell;

use glam::DVec3;

use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::geom::{Curve3, CurveEval, TrimmedCurve3};
use rcad_kernel::math::GeomAbsShape;

use super::approx_curvlin_func::ApproxCurvlinFunc;
use super::trihedron_law::{
    curve_first_parameter, curve_last_parameter, TrihedronLaw, TrihedronLawBase,
};
use super::trihedron_with_guide::{TrihedronWithGuide, TrihedronWithGuideBase};

/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-14;

/// OCCT GeomFill_GuideTrihedronAC (GeomFill_GuideTrihedronAC.hxx L110-125).
#[derive(Debug, Clone)]
pub struct GuideTrihedronAC {
    /// OCCT GeomFill_TrihedronLaw protected base (myCurve / myTrimmed).
    pub(crate) base: TrihedronLawBase,
    /// OCCT GeomFill_TrihedronWithGuide protected base
    /// (myGuide / myTrimG / myCurPointOnGuide).
    pub(crate) guide_base: TrihedronWithGuideBase,
    /// OCCT handle(Approx_CurvlinFunc) myGuideAC.
    my_guide_ac: ApproxCurvlinFunc,
    /// OCCT double Lguide.
    lguide: f64,
    /// OCCT handle(Approx_CurvlinFunc) myCurveAC — the OCCT null handle
    /// until SetCurve.
    my_curve_ac: Option<ApproxCurvlinFunc>,
    /// OCCT double L.
    l: f64,
    /// OCCT private handle(Adaptor3d_Curve) myCurve — SHADOWS the base
    /// member (C++ name hiding; see module doc).
    my_curve: Option<Curve3>,
    /// OCCT double UTol.
    utol: f64,
    /// OCCT double STol.
    #[allow(dead_code)]
    stol: f64,
    /// OCCT double Orig1.
    orig1: f64,
    /// OCCT double Orig2.
    orig2: f64,
}

/// OCCT myTrimmed->D1(Param, P, To).
fn curve_d1(c: &Curve3, u: f64) -> (DVec3, DVec3) {
    (c.point_at(u), c.derivative_at(u))
}

/// OCCT myTrimmed->D2(Param, P, To, DTo).
fn curve_d2(c: &Curve3, u: f64) -> (DVec3, DVec3, DVec3) {
    (c.point_at(u), c.derivative_at(u), c.derivative2_at(u))
}

/// OCCT myTrimmed->D3(Param, P, To, DTo, D2To).
fn curve_d3(c: &Curve3, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
    (
        c.point_at(u),
        c.derivative_at(u),
        c.derivative2_at(u),
        c.derivative3_at(u),
    )
}

impl GuideTrihedronAC {
    /// OCCT GeomFill_GuideTrihedronAC::GeomFill_GuideTrihedronAC (L35-47).
    pub fn new(guide: &Curve3) -> Self {
        // OCCT ctor: myCurve.Nullify() (the SHADOW member); myGuide = guide;
        // myTrimG = guide; myGuideAC = new Approx_CurvlinFunc(myGuide,
        // 1.e-7); Lguide = myGuideAC->GetLength(); UTol = STol =
        // Precision::PConfusion(); Orig1 = 0; Orig2 = 1.
        let my_guide_ac = ApproxCurvlinFunc::new(guide, 1.0e-7);
        let lguide = my_guide_ac.get_length();
        GuideTrihedronAC {
            base: TrihedronLawBase::default(),
            guide_base: TrihedronWithGuideBase {
                my_guide: Some(guide.clone()),
                my_trim_g: Some(guide.clone()),
                my_cur_point_on_guide: Cell::new(DVec3::ZERO),
            },
            my_guide_ac,
            lguide,
            my_curve_ac: None,
            l: 0.0,
            my_curve: None,
            utol: P_CONFUSION,
            stol: P_CONFUSION,
            orig1: 0.0, // origines pour le cas path multi-edges
            orig2: 1.0,
        }
    }

    /// OCCT Origine (L381-386).
    pub fn set_origine(&mut self, or_acr1: f64, or_acr2: f64) {
        self.orig1 = or_acr1;
        self.orig2 = or_acr2;
    }
}

impl TrihedronLaw for GuideTrihedronAC {
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

    /// OCCT Copy (L267-275).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let guide = self
            .guide_base
            .my_guide
            .clone()
            .expect("null myGuide in GuideTrihedronAC::Copy");
        let mut copy = GuideTrihedronAC::new(&guide);
        if let Some(curve) = self.my_curve.clone() {
            TrihedronLaw::set_curve(&mut copy, curve);
        }
        copy.set_origine(self.orig1, self.orig2);
        Box::new(copy)
    }

    /// OCCT SetCurve (L277-288).
    fn set_curve(&mut self, c: Curve3) -> bool {
        // OCCT: myCurve = C — the SHADOW member; myTrimmed = C — the base
        // member (no shadow there).
        self.my_curve = Some(c.clone());
        self.base.my_trimmed = Some(c.clone());
        // OCCT: if (!myCurve.IsNull()) { myCurveAC = new
        // Approx_CurvlinFunc(C, 1.e-7); L = myCurveAC->GetLength();
        // //    CorrectOrient(myGuide); }
        let my_curve_ac = ApproxCurvlinFunc::new(&c, 1.0e-7);
        self.l = my_curve_ac.get_length();
        self.my_curve_ac = Some(my_curve_ac);
        true
    }

    /// OCCT D0 (L56-87).
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        // double s = myCurveAC->GetSParameter(Param); // abscisse curviligne <=> Param
        let s = self
            .my_curve_ac
            .as_ref()
            .expect("null myCurveAC")
            .get_s_parameter(param);
        // double OrigG = Orig1 + s * (Orig2 - Orig1); // abscisse curv sur le guide (cas multi-edges)
        let orig_g = self.orig1 + s * (self.orig2 - self.orig1);
        // double tG = myGuideAC->GetUParameter(*myGuide, OrigG, 1); // param <=> s sur theGuide
        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        let tg = self.my_guide_ac.get_u_parameter(my_guide, orig_g, 1);

        let my_trimmed = self.base.my_trimmed.as_ref().expect("null myTrimmed");
        let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
        // myTrimmed->D1(Param, P, To); // point et derivee au parametre Param sur myCurve
        let (p, to) = curve_d1(my_trimmed, param);
        // myTrimG->D0(tG, PG); // point au parametre tG sur myGuide
        let pg = my_trim_g.point_at(tg);
        // myCurPointOnGuide = PG.
        self.guide_base.my_cur_point_on_guide.set(pg);

        // gp_Vec n(P, PG); // vecteur definissant la normale
        let n = pg - p;

        // TODO: finding #8 - no zero-magnitude guard before
        // Normalized()/Magnitude().  Adding guards (return false on
        // near-zero) caused blend regressions.  Needs investigation to find
        // a safe fallback strategy.
        *normal = n.normalize_or_zero();
        let b = to.cross(*normal);
        *binormal = b / b.length();
        *tangent = normal.cross(*binormal);
        *tangent = tangent.normalize_or_zero();

        true
    }

    /// OCCT D1 (L89-156).
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
        // triedre
        // abscisse curviligne <=> Param
        let s = self
            .my_curve_ac
            .as_ref()
            .expect("null myCurveAC")
            .get_s_parameter(param);
        // parametre <=> s sur theGuide
        let orig_g = self.orig1 + s * (self.orig2 - self.orig1);
        // parametre <=> s sur  theGuide
        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        let tg = self.my_guide_ac.get_u_parameter(my_guide, orig_g, 1);

        let my_trimmed = self.base.my_trimmed.as_ref().expect("null myTrimmed");
        let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
        // myTrimmed->D2(Param, P, To, DTo).
        let (p, to, d_to) = curve_d2(my_trimmed, param);
        // myTrimG->D1(tG, PG, TG).
        let (pg, tg_vec) = curve_d1(my_trim_g, tg);
        // myCurPointOnGuide = PG.
        self.guide_base.my_cur_point_on_guide.set(pg);

        let mut n = pg - p;
        let mut norm = n.length();
        if norm < 1.0e-12 {
            norm = 1.0;
        }

        n /= norm;
        // TODO: finding #8 - no zero-magnitude guard before TG.Magnitude()/L
        // division.  Adding guards caused blend regressions.  Needs safe
        // fallback strategy.
        let dtg = (self.orig2 - self.orig1)
            * (to.length() / tg_vec.length())
            * (self.lguide / self.l);
        // OCCT: dn.SetLinearForm(dtg, TG, -1, To).
        let mut dn = dtg * tg_vec - to;
        dn /= norm;

        // triedre
        *normal = n;
        let mut b = to.cross(*normal);
        let norm_b = b.length();
        b /= norm_b;

        *binormal = b;

        *tangent = normal.cross(*binormal);
        *tangent = tangent.normalize_or_zero();

        // derivee premiere
        // OCCT: DNormal.SetLinearForm(-(n.Dot(dn)), n, dn).
        *dnormal = -(n.dot(dn)) * n + dn;

        // OCCT: BPrim.SetLinearForm(DTo.Crossed(Normal), To.Crossed(DNormal)).
        let b_prim = d_to.cross(*normal) + to.cross(*dnormal);

        // OCCT: DBiNormal.SetLinearForm(-(B.Dot(BPrim)), B, BPrim).
        *dbinormal = -(b.dot(b_prim)) * b + b_prim;
        *dbinormal /= norm_b;

        // OCCT: DTangent.SetLinearForm(Normal.Crossed(DBiNormal),
        // DNormal.Crossed(BiNormal)).
        *dtangent = normal.cross(*dbinormal) + dnormal.cross(*binormal);

        true
    }

    /// OCCT D2 (L158-265).
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
        // abscisse curviligne <=> Param
        let s = self
            .my_curve_ac
            .as_ref()
            .expect("null myCurveAC")
            .get_s_parameter(param);
        // parametre <=> s sur theGuide
        let orig_g = self.orig1 + s * (self.orig2 - self.orig1);
        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        let tg = self.my_guide_ac.get_u_parameter(my_guide, orig_g, 1);

        let my_trimmed = self.base.my_trimmed.as_ref().expect("null myTrimmed");
        let my_trim_g = self.guide_base.my_trim_g.as_ref().expect("null myTrimG");
        // myTrimmed->D3(Param, P, To, DTo, D2To).
        let (p, to, d_to, d2_to) = curve_d3(my_trimmed, param);
        // myTrimG->D2(tG, PG, TG, DTG).
        let (pg, tg_vec, dtg_vec) = curve_d2(my_trim_g, tg);
        // myCurPointOnGuide = PG.
        self.guide_base.my_cur_point_on_guide.set(pg);

        let n_to = to.length();
        let n2_to = to.dot(to);
        let n_tg = tg_vec.length();
        let n2_tp = tg_vec.dot(tg_vec);
        let dtg_dt = (self.orig2 - self.orig1) * (n_to / n_tg) * (self.lguide / self.l);

        // gp_Vec n(P, PG); // vecteur definissant la normale
        let mut n = pg - p;
        let mut norm = n.length();
        // derivee de n par rapport a Param
        // OCCT: dn.SetLinearForm(dtg_dt, TG, -1, To).
        let mut dn = dtg_dt * tg_vec - to;

        // derivee seconde de tG par rapport a Param
        let d2tp_dt2 = (self.orig2 - self.orig1)
            * (self.lguide / self.l)
            * (d_to.dot(to) / (n_to * n_tg)
                - n2_to * tg_vec.dot(dtg_vec) * (self.lguide / self.l) / (n2_tp * n2_tp));
        // derivee seconde de n par rapport a Param
        // OCCT: d2n.SetLinearForm(dtg_dt * dtg_dt, DTG, d2tp_dt2, TG, -1, DTo).
        let mut d2n = (dtg_dt * dtg_dt) * dtg_vec + d2tp_dt2 * tg_vec - d_to;

        if norm > 1.0e-9 {
            n /= norm;
            dn /= norm;
            d2n /= norm;
        }
        // triedre
        *normal = n;

        let mut tn = to.cross(*normal);

        let norma = tn.length();
        if norma > 1.0e-9 {
            tn /= norma;
        }

        *binormal = tn;

        *tangent = normal.cross(*binormal);
        //  Tangent.Normalize();

        // derivee premiere du triedre
        let ndn = n.dot(dn);
        // OCCT: DNormal.SetLinearForm(-ndn, n, dn).
        *dnormal = -ndn * n + dn;

        // OCCT: DTN.SetLinearForm(DTo.Crossed(Normal), To.Crossed(DNormal)).
        let mut dtn = d_to.cross(*normal) + to.cross(*dnormal);
        dtn /= norma;
        let tn_dtn = tn.dot(dtn);

        // OCCT: DBiNormal.SetLinearForm(-TNDTN, TN, DTN).
        *dbinormal = -tn_dtn * tn + dtn;

        // OCCT: DTangent.SetLinearForm(Normal.Crossed(DBiNormal),
        // DNormal.Crossed(BiNormal)).
        *dtangent = normal.cross(*dbinormal) + dnormal.cross(*binormal);

        // derivee seconde du triedre
        let tn2 = tn.dot(tn);

        // OCCT: D2Normal.SetLinearForm(-2 * ndn, dn,
        // 3 * ndn * ndn - (dn.SquareMagnitude() + n.Dot(d2n)), n, d2n).
        *d2normal = (-2.0 * ndn) * dn + (3.0 * ndn * ndn - (dn.dot(dn) + n.dot(d2n))) * n + d2n;

        // OCCT: D2TN.SetLinearForm(1, D2To.Crossed(Normal), 2,
        // DTo.Crossed(DNormal), To.Crossed(D2Normal)).
        let mut d2tn = d2_to.cross(*normal) + 2.0 * d_to.cross(*dnormal) + to.cross(*d2normal);
        d2tn /= norma;

        // OCCT: D2BiNormal.SetLinearForm(-2 * TNDTN, DTN,
        // 3 * TNDTN * TNDTN - (TN2 + TN.Dot(D2TN)), TN, D2TN).
        *d2binormal = (-2.0 * tn_dtn) * dtn + (3.0 * tn_dtn * tn_dtn - (tn2 + tn.dot(d2tn))) * tn
            + d2tn;

        // OCCT: D2Tangent.SetLinearForm(1, D2Normal.Crossed(BiNormal), 2,
        // DNormal.Crossed(DBiNormal), Normal.Crossed(D2BiNormal)).
        *d2tangent = d2normal.cross(*binormal)
            + 2.0 * dnormal.cross(*dbinormal)
            + normal.cross(*d2binormal);

        true
    }

    /// OCCT NbIntervals (L290-305).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let my_curve_ac = self.my_curve_ac.as_ref().expect("null myCurveAC");
        let nb_c = my_curve_ac.nb_intervals(s);
        let mut disc_c = vec![0.0; nb_c + 1];
        my_curve_ac.intervals(&mut disc_c, s);
        let nb_g = self.my_guide_ac.nb_intervals(s);
        let mut disc_g = vec![0.0; nb_g + 1];
        self.my_guide_ac.intervals(&mut disc_g, s);

        let mut seq: Vec<f64> = Vec::new();
        fuse_intervals(&disc_c, &disc_g, &mut seq, P_CONFUSION, true);

        seq.len() - 1
    }

    /// OCCT Intervals (L307-324).
    fn intervals(&self, tt: &mut Vec<f64>, s: GeomAbsShape) {
        let my_curve_ac = self.my_curve_ac.as_ref().expect("null myCurveAC");
        let nb_c = my_curve_ac.nb_intervals(s);
        let mut disc_c = vec![0.0; nb_c + 1];
        my_curve_ac.intervals(&mut disc_c, s);
        let nb_g = self.my_guide_ac.nb_intervals(s);
        let mut disc_g = vec![0.0; nb_g + 1];
        self.my_guide_ac.intervals(&mut disc_g, s);

        let mut seq: Vec<f64> = Vec::new();
        fuse_intervals(&disc_c, &disc_g, &mut seq, P_CONFUSION, true);
        let nb = seq.len();

        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        for ii in 0..nb {
            tt[ii] = my_curve_ac.get_u_parameter(my_curve, seq[ii], 1);
        }
    }

    /// OCCT SetInterval (L326-343).
    fn set_interval(&mut self, first: f64, last: f64) {
        // OCCT: myTrimmed = myCurve->Trim(First, Last, UTol).
        let my_curve = self.my_curve.as_ref().expect("null myCurve (shadow)");
        self.base.my_trimmed = Some(Curve3::Trimmed(TrimmedCurve3::new(
            my_curve.clone(),
            first,
            last,
        )));

        let my_curve_ac = self.my_curve_ac.as_ref().expect("null myCurveAC");
        let mut sf;
        let mut sl;

        sf = my_curve_ac.get_s_parameter(first);
        sl = my_curve_ac.get_s_parameter(last);
        //  if (Sl>1) Sl=1;
        //  myCurveAC->Trim(Sf, Sl, UTol);

        let my_guide = self.guide_base.my_guide.as_ref().expect("null myGuide");
        let mut u = self.orig1 + sf * (self.orig2 - self.orig1);
        sf = self.my_guide_ac.get_u_parameter(my_guide, u, 1);
        u = self.orig1 + sl * (self.orig2 - self.orig1);
        sl = self.my_guide_ac.get_u_parameter(my_guide, u, 1);
        // OCCT: myTrimG = myGuide->Trim(Sf, Sl, UTol).
        self.guide_base.my_trim_g = Some(Curve3::Trimmed(TrimmedCurve3::new(
            my_guide.clone(),
            sf,
            sl,
        )));
    }

    /// OCCT GetAverageLaw (L345-369).
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

    /// OCCT IsConstant (L371-374).
    fn is_constant(&self) -> bool {
        false
    }

    /// OCCT IsOnlyBy3dCurve (L376-379).
    fn is_only_by3d_curve(&self) -> bool {
        false
    }
}

impl TrihedronWithGuide for GuideTrihedronAC {
    /// OCCT Guide (L49-53).
    fn guide(&self) -> Option<Curve3> {
        self.guide_base.my_guide.clone()
    }

    /// OCCT Origine (L381-386).
    fn origine(&mut self, param1: f64, param2: f64) {
        self.set_origine(param1, param2);
    }

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
