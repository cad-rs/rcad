//! OCCT GeomFill_CircularBlendFunc (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_CircularBlendFunc.hxx (members) + GeomFill_CircularBlendFunc.cxx
//! (whole file L24-676, including the file-local TolAng constant,
//! GeomFillNextShape and GeomFillFusInt statics).
//!
//! The consumed GeomFill package statics (`GetShape` / `Knots` / `Mults` /
//! `GetCircle` x3 / `GetTolerance` / `GetMinimalWeights`) are already
//! translated in [`super::geom_fill`] (earlier round).  GAP carrier: the 3D
//! `GCPnts_QuasiUniformDeflection` ([`QuasiUniformDeflection3d`]) is not
//! translated (rcad hosts the 2D variant only) — its use inside Discret
//! keeps the OCCT failure path.

use glam::{DVec2, DVec3};

use rcad_kernel::base::convert::ConvertParameterisation;
use rcad_kernel::geom::{Curve3, CurveEval, TrimmedCurve3};
use rcad_kernel::math::GeomAbsShape;

use super::approx_sweep_function::ApproxSweepFunction;
use super::frenet::{curve_intervals, curve_nb_intervals};
use super::geom_fill as geom_fill_statics;
use super::sweep_section_generator::QuasiUniformDeflection3d;

/// OCCT static const double TolAng = 1.e-6 (L26).
const TOL_ANG: f64 = 1.0e-6;
/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-14;

/// OCCT static GeomFillNextShape (L28-40).
fn geom_fill_next_shape(s: GeomAbsShape) -> GeomAbsShape {
    match s {
        GeomAbsShape::C0 => GeomAbsShape::C1,
        GeomAbsShape::C1 => GeomAbsShape::C2,
        GeomAbsShape::C2 => GeomAbsShape::C3,
        GeomAbsShape::C3 => GeomAbsShape::CN,
        _ => GeomAbsShape::CN,
    }
}

/// OCCT static GeomFillFusInt (L42-106) — the interval merge (Epspar =
/// Precision::PConfusion() * 0.99, "en suposant que le positionement
/// fonctionne a PConfusion()/2").
fn geom_fill_fus_int(i1: &[f64], i2: &[f64], seq: &mut Vec<f64>) {
    let mut ind1 = 1usize;
    let mut ind2 = 1usize;
    let epspar: f64 = P_CONFUSION * 0.99;
    let mut v1;
    let mut v2;
    // Initialisations : les IND1 et IND2 pointent sur le 1er element
    // de chacune des 2 tables a traiter.

    //--- On remplit TABSOR en parcourant TABLE1 et TABLE2 simultanement ---
    //------------------ en eliminant les occurrences multiples ------------

    while ind1 <= i1.len() && ind2 <= i2.len() {
        v1 = i1[ind1 - 1];
        v2 = i2[ind2 - 1];
        if (v1 - v2).abs() <= epspar {
            // Ici les elements de I1 et I2 conviennent .
            seq.push((v1 + v2) / 2.0);
            ind1 += 1;
            ind2 += 1;
        } else if v1 < v2 {
            // Ici l' element de I1 convient.
            seq.push(v1);
            ind1 += 1;
        } else {
            // Ici l' element de TABLE2 convient.
            seq.push(v2);
            ind2 += 1;
        }
    }

    if ind1 > i1.len() {
        //----- Ici I1 est epuise, on complete avec la fin de TABLE2 -------
        while ind2 <= i2.len() {
            seq.push(i2[ind2 - 1]);
            ind2 += 1;
        }
    }

    if ind2 > i2.len() {
        //----- Ici I2 est epuise, on complete avec la fin de I1 -------
        while ind1 <= i1.len() {
            seq.push(i1[ind1 - 1]);
            ind1 += 1;
        }
    }
}

/// OCCT GeomFill_CircularBlendFunc (hxx L153-172).
#[derive(Debug, Clone)]
pub struct CircularBlendFunc {
    /// OCCT gp_Pnt myBary.
    my_bary: DVec3,
    /// OCCT double myRadius.
    my_radius: f64,
    /// OCCT double maxang.
    maxang: f64,
    /// OCCT double minang.
    minang: f64,
    /// OCCT double distmin.
    distmin: f64,
    /// OCCT handle(Adaptor3d_Curve) myPath.
    my_path: Curve3,
    /// OCCT handle(Adaptor3d_Curve) myCurve1.
    my_curve1: Curve3,
    /// OCCT handle(Adaptor3d_Curve) myCurve2.
    my_curve2: Curve3,
    /// OCCT handle(Adaptor3d_Curve) myTPath.
    my_tpath: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myTCurve1.
    my_tcurve1: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myTCurve2.
    my_tcurve2: Option<Curve3>,
    /// OCCT int myDegree.
    my_degree: i32,
    /// OCCT int myNbKnots.
    my_nb_knots: i32,
    /// OCCT int myNbPoles.
    my_nb_poles: i32,
    /// OCCT Convert_ParameterisationType myTConv.
    my_tconv: ConvertParameterisation,
    /// OCCT bool myreverse.
    myreverse: bool,
}

impl CircularBlendFunc {
    /// OCCT GeomFill_CircularBlendFunc::GeomFill_CircularBlendFunc
    /// (L108-143).
    pub fn new(
        path: &Curve3,
        curve1: &Curve3,
        curve2: &Curve3,
        radius: f64,
        polynomial: bool,
    ) -> Self {
        // OCCT: maxang(RealFirst()), minang(RealLast()), distmin(RealLast()).
        let mut blend = CircularBlendFunc {
            my_bary: DVec3::ZERO,
            my_radius: radius,
            maxang: f64::MIN,
            minang: f64::MAX,
            distmin: f64::MAX,
            my_path: path.clone(),
            my_curve1: curve1.clone(),
            my_curve2: curve2.clone(),
            my_tpath: Some(path.clone()),
            my_tcurve1: Some(curve1.clone()),
            my_tcurve2: Some(curve2.clone()),
            my_degree: 0,
            my_nb_knots: 0,
            my_nb_poles: 0,
            my_tconv: ConvertParameterisation::QuasiAngular,
            myreverse: false,
        };

        // Recopie des arguments — done above.

        // Estimations numeriques
        blend.discret();

        // Type de convertion ?
        if polynomial {
            blend.my_tconv = ConvertParameterisation::Polynomial;
        } else if blend.maxang > 0.65 * std::f64::consts::PI {
            blend.my_tconv = ConvertParameterisation::QuasiAngular; // car c'est Continue
        } else {
            blend.my_tconv = ConvertParameterisation::TgtThetaOver2;
        }
        // car c'est le plus performant

        // On en deduit la structure
        // OCCT: GeomFill::GetShape(maxang, myNbPoles, myNbKnots, myDegree,
        // myTConv).
        let mut nb_poles = 0i32;
        let mut nb_knots = 0i32;
        let mut degree = 0i32;
        geom_fill_statics::get_shape(
            blend.maxang,
            &mut nb_poles,
            &mut nb_knots,
            &mut degree,
            &mut blend.my_tconv,
        );
        blend.my_nb_poles = nb_poles;
        blend.my_nb_knots = nb_knots;
        blend.my_degree = degree;
        blend
    }

    /// OCCT Discret (L145-266).
    fn discret(&mut self) {
        let t_first = self.my_path.default_domain()[0];
        let t_last = self.my_path.default_domain()[1];
        let mut t;
        let l1;
        let l2;
        let mut l;
        let percent;

        let mut p1 = self.my_curve1.point_at(t_first);
        let mut p2 = self.my_curve1.point_at((t_first + t_last) / 2.0);
        let mut p3 = self.my_curve1.point_at(t_last);
        l1 = p1.distance(p2) + p2.distance(p3);

        p1 = self.my_curve2.point_at(t_first);
        p2 = self.my_curve2.point_at((t_first + t_last) / 2.0);
        p3 = self.my_curve2.point_at(t_last);
        l2 = p1.distance(p2) + p2.distance(p3);

        let c: &Curve3;
        if l1 > l2 {
            l = l1;
            c = &self.my_curve1;
        } else {
            l = l2;
            c = &self.my_curve2;
        }
        let _ = &mut l;

        let fleche = 1.0e-2 * l;
        let mut angle;
        let mut cosa;
        // OCCT: GCPnts_QuasiUniformDeflection Samp; Samp.Initialize(*C,
        // Fleche) — the 3D QuasiUniformDeflection is the GAP carrier.
        let samp = QuasiUniformDeflection3d::initialize(c, fleche);
        self.my_bary = DVec3::ZERO;
        let mut ns1;
        let mut ns2;

        if samp.is_done() {
            percent = 1.0 / (2.0 * samp.nb_points() as f64);
            //    char name[100];
            for ii in 1..=samp.nb_points() {
                t = samp.parameter(ii);
                p1 = self.my_curve1.point_at(t);
                p2 = self.my_curve2.point_at(t);
                let center = self.my_path.point_at(t);
                ns1 = center - p1;
                ns2 = center - p2;
                ns1 = ns1.normalize_or_zero();
                ns2 = ns2.normalize_or_zero();
                cosa = ns1.dot(ns2);
                if cosa > 1.0 {
                    cosa = 1.0;
                }
                angle = cosa.acos().abs();
                if angle > self.maxang {
                    self.maxang = angle;
                }
                if angle < self.minang {
                    self.minang = angle;
                }
                self.distmin = self.distmin.min(p1.distance(p2));
                self.my_bary += p1 + p2;
            }
        } else {
            percent = 1.0 / 42.0;
            let delta = (t_last - t_first) / 20.0;
            let mut ii = 0usize;
            t = t_first;
            while ii <= 20 {
                p1 = self.my_curve1.point_at(t);
                p2 = self.my_curve2.point_at(t);
                let center = self.my_path.point_at(t);

                ns1 = center - p1;
                ns2 = center - p2;
                ns1 = ns1.normalize_or_zero();
                ns2 = ns2.normalize_or_zero();
                cosa = ns1.dot(ns2);
                if cosa > 1.0 {
                    cosa = 1.0;
                }
                angle = cosa.acos().abs();

                if angle > self.maxang {
                    self.maxang = angle;
                }
                if angle < self.minang {
                    self.minang = angle;
                }
                self.distmin = self.distmin.min(p1.distance(p2));
                self.my_bary += p1 + p2;

                ii += 1;
                t += delta;
            }
        }
        self.my_bary *= percent;

        // Faut il inverser la trajectoire ?
        t = (t_first + t_last) / 2.0;
        p1 = self.my_curve1.point_at(t);
        p2 = self.my_curve2.point_at(t);
        let center = self.my_path.point_at(t);
        let d_center = self.my_path.derivative_at(t);

        ns1 = center - p1;
        ns2 = center - p2;

        // myreverse = (DCenter.Dot(ns1.Crossed(ns2)) < 0);
        self.myreverse = false;
        let _ = d_center;
    }
}

impl ApproxSweepFunction for CircularBlendFunc {
    /// OCCT D0 (L268-310).
    fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        _poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        let my_tpath = self.my_tpath.as_ref().expect("null myTPath");
        let my_tcurve1 = self.my_tcurve1.as_ref().expect("null myTCurve1");
        let my_tcurve2 = self.my_tcurve2.as_ref().expect("null myTCurve2");

        // Positionnement
        let mut center = my_tpath.point_at(param);
        let p1 = my_tcurve1.point_at(param);
        let p2 = my_tcurve2.point_at(param);
        let mut ns1 = center - p1;
        let mut ns2 = center - p2;
        let mut nplan;
        if !gp_vec_is_parallel(ns1, ns2, TOL_ANG) {
            nplan = ns1.cross(ns2);
        } else {
            nplan = my_tpath.derivative_at(param);
            if self.myreverse {
                nplan = -nplan;
            }
        }

        // Normalisation
        ns1 = ns1.normalize_or_zero();
        ns2 = ns2.normalize_or_zero();
        nplan = nplan.normalize_or_zero();

        // OCCT: temp.SetLinearForm(myRadius, ns1, myRadius, ns2, 1, P1, P2).
        let temp = self.my_radius * ns1 + self.my_radius * ns2 + p1 + p2;
        center = 0.5 * temp;

        // Section — OCCT: GeomFill::GetCircle(myTConv, ns1, ns2, nplan, P1,
        // P2, myRadius, Center, Poles, Weigths).
        geom_fill_statics::get_circle(
            self.my_tconv,
            ns1,
            ns2,
            nplan,
            p1,
            p2,
            self.my_radius,
            center,
            poles,
            weigths,
        );

        true
    }

    /// OCCT D1 (L312-401).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        _poles2d: &mut [DVec2],
        _dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        let my_tpath = self.my_tpath.as_ref().expect("null myTPath");
        let my_tcurve1 = self.my_tcurve1.as_ref().expect("null myTCurve1");
        let my_tcurve2 = self.my_tcurve2.as_ref().expect("null myTCurve2");

        // Positionemment
        let mut center = my_tpath.point_at(param);
        let d_center0 = my_tpath.derivative_at(param);
        let p1 = my_tcurve1.point_at(param);
        let dp1 = my_tcurve1.derivative_at(param);
        let p2 = my_tcurve2.point_at(param);
        let dp2 = my_tcurve2.derivative_at(param);

        let mut ns1 = center - p1;
        let mut ns2 = center - p2;
        let mut dns1 = d_center0 - dp1;
        let mut dns2 = d_center0 - dp2;

        let mut nplan;
        let mut dnplan;
        if !gp_vec_is_parallel(ns1, ns2, TOL_ANG) {
            nplan = ns1.cross(ns2);
            dnplan = dns1.cross(ns2) + ns1.cross(dns2);
        } else {
            nplan = my_tpath.derivative_at(param);
            dnplan = my_tpath.derivative2_at(param);
            if self.myreverse {
                nplan = -nplan;
                dnplan = -dnplan;
            }
        }

        // Normalisation
        let invnorm1 = 1.0 / ns1.length();
        let invnorm2 = 1.0 / ns2.length();

        ns1 *= invnorm1;
        dns1 = -(dns1.dot(ns1)) * ns1 + dns1;
        dns1 *= invnorm1;

        ns2 *= invnorm2;
        dns2 = -(dns2.dot(ns2)) * ns2 + dns2;
        dns2 *= invnorm2;

        let temp = self.my_radius * ns1 + self.my_radius * ns2 + p1 + p2;
        center = 0.5 * temp;
        let mut d_center = self.my_radius * dns1 + self.my_radius * dns2 + dp1 + dp2;
        d_center *= 0.5;

        let invnormp = 1.0 / nplan.length();
        nplan *= invnormp;
        dnplan = -(dnplan.dot(nplan)) * nplan + dnplan;
        dnplan *= invnormp;

        // OCCT: GeomFill::GetCircle(myTConv, ns1, ns2, Dns1, Dns2, nplan,
        // dnplan, P1, P2, DP1, DP2, myRadius, 0, Center, DCenter, Poles,
        // DPoles, Weigths, DWeigths).
        geom_fill_statics::get_circle_d1(
            self.my_tconv,
            ns1,
            ns2,
            dns1,
            dns2,
            nplan,
            dnplan,
            p1,
            p2,
            dp1,
            dp2,
            self.my_radius,
            0.0,
            center,
            d_center,
            poles,
            dpoles,
            weigths,
            dweigths,
        );
        true
    }

    /// OCCT D2 (L403-532).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        _poles2d: &mut [DVec2],
        _dpoles2d: &mut [DVec2],
        _d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        let my_tpath = self.my_tpath.as_ref().expect("null myTPath");
        let my_tcurve1 = self.my_tcurve1.as_ref().expect("null myTCurve1");
        let my_tcurve2 = self.my_tcurve2.as_ref().expect("null myTCurve2");

        // Positionement
        let mut center = my_tpath.point_at(param);
        let mut d_center = my_tpath.derivative_at(param);
        let d2_center0 = my_tpath.derivative2_at(param);
        let p1 = my_tcurve1.point_at(param);
        let dp1 = my_tcurve1.derivative_at(param);
        let d2p1 = my_tcurve1.derivative2_at(param);
        let p2 = my_tcurve2.point_at(param);
        let dp2 = my_tcurve2.derivative_at(param);
        let d2p2 = my_tcurve2.derivative2_at(param);

        let mut ns1 = center - p1;
        let mut dns1 = d_center - dp1;
        let mut d2ns1 = d2_center0 - d2p1;
        let mut ns2 = center - p2;
        let mut dns2 = d_center - dp2;
        let mut d2ns2 = d2_center0 - d2p2;

        let mut nplan;
        let mut dnplan;
        let mut d2nplan;
        if !gp_vec_is_parallel(ns1, ns2, TOL_ANG) {
            nplan = ns1.cross(ns2);
            dnplan = dns1.cross(ns2) + ns1.cross(dns2);
            // OCCT: d2nplan.SetLinearForm(1, D2ns1.Crossed(ns2), 2,
            // Dns1.Crossed(Dns2), ns1.Crossed(D2ns2)).
            d2nplan = d2ns1.cross(ns2) + 2.0 * dns1.cross(dns2) + ns1.cross(d2ns2);
        } else {
            nplan = my_tpath.derivative_at(param);
            dnplan = my_tpath.derivative2_at(param);
            d2nplan = my_tpath.derivative3_at(param);
            if self.myreverse {
                nplan = -nplan;
                dnplan = -dnplan;
                d2nplan = -d2nplan;
            }
        }

        // Normalisation
        let invnorm1 = 1.0 / ns1.length();
        let invnorm2 = 1.0 / ns2.length();

        ns1 *= invnorm1;
        let mut sc = dns1.dot(ns1);
        // OCCT: D2ns1.SetLinearForm(3*sc*sc*invnorm1 - D2ns1.Dot(ns1)
        // - invnorm1*Dns1.SquareMagnitude(), ns1, -2*sc*invnorm1, Dns1, D2ns1).
        d2ns1 = (3.0 * sc * sc * invnorm1 - d2ns1.dot(ns1) - invnorm1 * dns1.dot(dns1)) * ns1
            + (-2.0 * sc * invnorm1) * dns1
            + d2ns1;
        dns1 = -(dns1.dot(ns1)) * ns1 + dns1;
        d2ns1 *= invnorm1;
        dns1 *= invnorm1;

        ns2 *= invnorm2;
        sc = dns2.dot(ns2);
        d2ns2 = (3.0 * sc * sc * invnorm2 - d2ns2.dot(ns2) - invnorm2 * dns2.dot(dns2)) * ns2
            + (-2.0 * sc * invnorm2) * dns2
            + d2ns2;
        // OCCT quirk: the D2 branch reuses `sc` in the Dns2 formula instead
        // of Dns2.Dot(ns2) — reproduced verbatim.
        dns2 = -sc * ns2 + dns2;
        d2ns2 *= invnorm2;
        dns2 *= invnorm2;

        let temp = self.my_radius * ns1 + self.my_radius * ns2 + p1 + p2;
        center = 0.5 * temp;
        d_center = self.my_radius * dns1 + self.my_radius * dns2 + dp1 + dp2;
        d_center *= 0.5;
        let mut d2_center = self.my_radius * d2ns1 + self.my_radius * d2ns2 + d2p1 + d2p2;
        d2_center *= 0.5;

        let invnormp = 1.0 / nplan.length();
        nplan *= invnormp;
        sc = dnplan.dot(nplan);
        d2nplan = (3.0 * sc * sc * invnormp - d2nplan.dot(nplan) - invnormp * dnplan.dot(dnplan))
            * nplan
            + (-2.0 * sc * invnormp) * dnplan
            + d2nplan;
        dnplan = -sc * nplan + dnplan;
        dnplan *= invnormp;
        d2nplan *= invnormp;

        // OCCT: GeomFill::GetCircle(myTConv, ns1, ns2, Dns1, Dns2, D2ns1,
        // D2ns2, nplan, dnplan, d2nplan, P1, P2, DP1, DP2, D2P1, D2P2,
        // myRadius, 0, 0, Center, DCenter, D2Center, Poles, DPoles, D2Poles,
        // Weigths, DWeigths, D2Weigths).
        geom_fill_statics::get_circle_d2(
            self.my_tconv,
            ns1,
            ns2,
            dns1,
            dns2,
            d2ns1,
            d2ns2,
            nplan,
            dnplan,
            d2nplan,
            p1,
            p2,
            dp1,
            dp2,
            d2p1,
            d2p2,
            self.my_radius,
            0.0,
            0.0,
            center,
            d_center,
            d2_center,
            poles,
            dpoles,
            d2poles,
            weigths,
            dweigths,
            d2weigths,
        );
        true
    }

    /// OCCT Nb2dCurves (L534-537).
    fn nb_2d_curves(&self) -> usize {
        0
    }

    /// OCCT SectionShape (L539-544).
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize) {
        *nb_poles = self.my_nb_poles as usize;
        *nb_knots = self.my_nb_knots as usize;
        *degree = self.my_degree as usize;
    }

    /// OCCT Knots (L546-549).
    fn knots(&self, t_knots: &mut [f64]) {
        geom_fill_statics::knots(self.my_tconv, t_knots);
    }

    /// OCCT Mults (L551-554).
    fn mults(&self, t_mults: &mut [i32]) {
        geom_fill_statics::mults(self.my_tconv, t_mults);
    }

    /// OCCT IsRational (L556-559).
    fn is_rational(&self) -> bool {
        self.my_tconv != ConvertParameterisation::Polynomial
    }

    /// OCCT NbIntervals (L561-599).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let nb_i_center = curve_nb_intervals(&self.my_path, geom_fill_next_shape(s));
        let nb_i_cb1 = curve_nb_intervals(&self.my_curve1, s);
        let nb_i_cb2 = curve_nb_intervals(&self.my_curve2, s);

        let i_center = curve_intervals(&self.my_path, geom_fill_next_shape(s));
        let i_cb1 = curve_intervals(&self.my_curve1, s);
        let i_cb2 = curve_intervals(&self.my_curve2, s);

        let mut inter: Vec<f64> = Vec::new();
        geom_fill_fus_int(&i_cb1, &i_cb2, &mut inter);

        let icbs = inter.clone();

        let mut inter: Vec<f64> = Vec::new();
        geom_fill_fus_int(&i_center, &icbs, &mut inter);

        inter.len() - 1
    }

    /// OCCT Intervals (L601-634).
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let i_center = curve_intervals(&self.my_path, geom_fill_next_shape(s));
        let i_cb1 = curve_intervals(&self.my_curve1, s);
        let i_cb2 = curve_intervals(&self.my_curve2, s);

        let mut inter: Vec<f64> = Vec::new();
        geom_fill_fus_int(&i_cb1, &i_cb2, &mut inter);

        let icbs = inter.clone();

        let mut inter: Vec<f64> = Vec::new();
        geom_fill_fus_int(&i_center, &icbs, &mut inter);

        // Recopie du resultat
        for (ii, value) in inter.iter().enumerate() {
            t[ii] = *value;
        }
    }

    /// OCCT SetInterval (L636-642).
    fn set_interval(&mut self, first: f64, last: f64) {
        let eps = P_CONFUSION;
        self.my_tpath = Some(Curve3::Trimmed(TrimmedCurve3::new(
            self.my_path.clone(),
            first,
            last,
        )));
        self.my_tcurve1 = Some(Curve3::Trimmed(TrimmedCurve3::new(
            self.my_curve1.clone(),
            first,
            last,
        )));
        self.my_tcurve2 = Some(Curve3::Trimmed(TrimmedCurve3::new(
            self.my_curve2.clone(),
            first,
            last,
        )));
        let _ = eps;
    }

    /// OCCT GetTolerance (L644-656).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        let low = 0usize;
        let up = tol3d.len() - 1;

        // OCCT: Tol = GeomFill::GetTolerance(myTConv, minang, myRadius,
        // AngleTol, SurfTol).
        let tol = geom_fill_statics::get_tolerance(
            self.my_tconv,
            self.minang,
            self.my_radius,
            angle_tol,
            surf_tol,
        );
        for value in tol3d.iter_mut() {
            *value = surf_tol;
        }
        tol3d[low + 1] = tol.min(surf_tol);
        tol3d[up - 1] = tol.min(surf_tol);
        tol3d[low] = tol.min(bound_tol);
        tol3d[up] = tol.min(bound_tol);
    }

    /// OCCT SetTolerance (L658-661) — "y rien a faire !".
    fn set_tolerance(&mut self, _tol3d: f64, _tol2d: f64) {}

    /// OCCT BarycentreOfSurf (L663-666).
    fn barycentre_of_surf(&self) -> DVec3 {
        self.my_bary
    }

    /// OCCT MaximalSection (L668-671).
    fn maximal_section(&self) -> f64 {
        self.maxang * self.my_radius
    }

    /// OCCT GetMinimalWeight (L673-676).
    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        geom_fill_statics::get_minimal_weights(self.my_tconv, self.minang, self.maxang, weigths);
    }
}

// Re-export the gp helper used above (shared with the sweep generator).
use super::sweep_section_generator::gp_vec_is_parallel;
