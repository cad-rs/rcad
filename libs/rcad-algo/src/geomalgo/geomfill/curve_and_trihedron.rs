//! OCCT GeomFill_CurveAndTrihedron (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_CurveAndTrihedron.hxx (members) + GeomFill_CurveAndTrihedron.cxx
//! (whole file L35-344).
//!
//! Architecture differences:
//! - The `handle(GeomFill_TrihedronLaw)` member maps to
//!   `RefCell<Box<dyn TrihedronLaw>>` (SetCurve mutates the law state while
//!   D0/D1/D2 only read it through `&self`).
//! - The OCCT member caches `Point` / `V1` / `V2` / `V3` are written inside
//!   the const evaluators (`&self`); they are `Cell<DVec3>` — same cache,
//!   no cross-call reads in OCCT either.

use std::cell::{Cell, RefCell};

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::geom::{Curve3, CurveEval, TrimmedCurve3};
use rcad_kernel::math::GeomAbsShape;

use super::gp_mat::GpMat;
use super::location_law::LocationLaw;
use super::trihedron_law::{curve_first_parameter, curve_last_parameter, TrihedronLaw};

/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-14;

/// OCCT GeomFill_CurveAndTrihedron (GeomFill_CurveAndTrihedron.hxx
/// L100-110).
pub struct CurveAndTrihedron {
    /// OCCT handle(GeomFill_TrihedronLaw) myLaw.
    my_law: RefCell<Box<dyn TrihedronLaw>>,
    /// OCCT handle(Adaptor3d_Curve) myTrimmed.
    my_trimmed: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myCurve.
    my_curve: Option<Curve3>,
    /// OCCT gp_Pnt Point (D0/D1/D2 cache).
    point: Cell<DVec3>,
    /// OCCT gp_Vec V1 (D0 cache).
    v1: Cell<DVec3>,
    /// OCCT gp_Vec V2 (D0 cache).
    v2: Cell<DVec3>,
    /// OCCT gp_Vec V3 (D0 cache).
    v3: Cell<DVec3>,
    /// OCCT gp_Mat Trans.
    trans: GpMat,
    /// OCCT bool WithTrans.
    with_trans: bool,
}

impl CurveAndTrihedron {
    /// OCCT GeomFill_CurveAndTrihedron::GeomFill_CurveAndTrihedron
    /// (L35-42).
    pub fn new(trihedron: Box<dyn TrihedronLaw>) -> Self {
        CurveAndTrihedron {
            my_law: RefCell::new(trihedron),
            my_trimmed: None,
            my_curve: None,
            point: Cell::new(DVec3::ZERO),
            v1: Cell::new(DVec3::ZERO),
            v2: Cell::new(DVec3::ZERO),
            v3: Cell::new(DVec3::ZERO),
            trans: GpMat::identity(),
            with_trans: false,
        }
    }

    /// OCCT SetTrsf (L72-90).
    pub fn set_trsf(&mut self, transfo: GpMat) {
        self.trans = transfo;
        // OCCT: gp_Mat Aux; Aux.SetIdentity(); Aux -= Trans.
        let mut aux = GpMat::identity();
        aux = aux.added(&self.trans.multiplied_scalar(-1.0));
        self.with_trans = false; // Au cas ou Trans = I
        let mut ii = 1;
        while ii <= 3 && !self.with_trans {
            let mut jj = 1;
            while jj <= 3 && !self.with_trans {
                if (aux.mat[ii - 1][jj - 1]).abs() > 1.0e-14 {
                    self.with_trans = true;
                }
                jj += 1;
            }
            ii += 1;
        }
    }
}

impl LocationLaw for CurveAndTrihedron {
    /// OCCT Copy (L46-55).
    fn copy_law(&self) -> Box<dyn LocationLaw> {
        let law_copy = self.my_law.borrow().copy_law();
        let mut copy = CurveAndTrihedron::new(law_copy);
        if let Some(curve) = self.my_curve.clone() {
            LocationLaw::set_curve(&mut copy, curve);
        }
        copy.set_trsf(self.trans);
        Box::new(copy)
    }

    /// OCCT SetCurve (L58-63).
    fn set_curve(&mut self, c: Curve3) -> bool {
        self.my_curve = Some(c.clone());
        self.my_trimmed = Some(c.clone());
        self.my_law.borrow_mut().set_curve(c)
    }

    /// OCCT GetCurve (L65-68).
    fn get_curve(&self) -> Option<Curve3> {
        self.my_curve.clone()
    }

    /// OCCT SetTrsf — the trait entry; the body lives in the inherent
    /// method above.
    fn set_trsf(&mut self, transfo: GpMat) {
        CurveAndTrihedron::set_trsf(self, transfo);
    }

    /// OCCT D0(Param, M, V) (L93-108).
    fn d0(&self, param: f64, m: &mut GpMat, v: &mut DVec3) -> bool {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        // myTrimmed->D0(Param, Point).
        self.point.set(my_trimmed.point_at(param));
        *v = self.point.get();

        let mut v1 = self.v1.get();
        let mut v2 = self.v2.get();
        let mut v3 = self.v3.get();
        let ok = self
            .my_law
            .borrow()
            .d0(param, &mut v1, &mut v2, &mut v3);
        self.v1.set(v1);
        self.v2.set(v2);
        self.v3.set(v3);
        m.set_cols(self.v2.get(), self.v3.get(), self.v1.get());

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
        }
        ok
    }

    /// OCCT D0(Param, M, V, Pnts2d) (L111-130).
    fn d0_2d(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut DVec3,
        _pnts2d: &mut [DVec2],
    ) -> bool {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        // myTrimmed->D0(Param, Point).
        self.point.set(my_trimmed.point_at(param));
        *v = self.point.get();

        let mut v1 = self.v1.get();
        let mut v2 = self.v2.get();
        let mut v3 = self.v3.get();
        let ok = self
            .my_law
            .borrow()
            .d0(param, &mut v1, &mut v2, &mut v3);
        self.v1.set(v1);
        self.v2.set(v2);
        self.v3.set(v3);
        m.set_cols(self.v2.get(), self.v3.get(), self.v1.get());

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
        }
        ok
    }

    /// OCCT D1 (L132-158).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut DVec3,
        dm: &mut GpMat,
        dv: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _vecs2d: &mut [DVec2],
    ) -> bool {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        // myTrimmed->D1(Param, Point, DV).
        self.point.set(my_trimmed.point_at(param));
        let d_point = my_trimmed.derivative_at(param);
        *v = self.point.get();
        *dv = d_point;

        // OCCT: gp_Vec DV1, DV2, DV3; myLaw->D1(Param, V1, DV1, V2, DV2, V3, DV3).
        let mut v1 = self.v1.get();
        let mut v2 = self.v2.get();
        let mut v3 = self.v3.get();
        let (mut dv1, mut dv2, mut dv3) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let ok = self.my_law.borrow().d1(
            param,
            &mut v1,
            &mut dv1,
            &mut v2,
            &mut dv2,
            &mut v3,
            &mut dv3,
        );
        self.v1.set(v1);
        self.v2.set(v2);
        self.v3.set(v3);
        m.set_cols(self.v2.get(), self.v3.get(), self.v1.get());
        dm.set_cols(dv2, dv3, dv1);

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
            *dm = dm.multiplied_mat(&self.trans);
        }

        ok
    }

    /// OCCT D2 (L160-193).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut DVec3,
        dm: &mut GpMat,
        dv: &mut DVec3,
        d2m: &mut GpMat,
        d2v: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _d1vecs2d: &mut [DVec2],
        _d2vecs2d: &mut [DVec2],
    ) -> bool {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        // myTrimmed->D2(Param, Point, DV, D2V).
        self.point.set(my_trimmed.point_at(param));
        let d_point = my_trimmed.derivative_at(param);
        let d2_point = my_trimmed.derivative2_at(param);
        *v = self.point.get();
        *dv = d_point;
        *d2v = d2_point;

        // OCCT: gp_Vec DV1..D2V3; myLaw->D2(Param, V1, DV1, D2V1, V2, DV2,
        // D2V2, V3, DV3, D2V3).
        let mut v1 = self.v1.get();
        let mut v2 = self.v2.get();
        let mut v3 = self.v3.get();
        let (mut dv1, mut dv2, mut dv3) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let (mut d2v1, mut d2v2, mut d2v3) = (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        let ok = self.my_law.borrow().d2(
            param,
            &mut v1,
            &mut dv1,
            &mut d2v1,
            &mut v2,
            &mut dv2,
            &mut d2v2,
            &mut v3,
            &mut dv3,
            &mut d2v3,
        );
        self.v1.set(v1);
        self.v2.set(v2);
        self.v3.set(v3);

        m.set_cols(self.v2.get(), self.v3.get(), self.v1.get());
        dm.set_cols(dv2, dv3, dv1);
        d2m.set_cols(d2v2, d2v3, d2v1);

        if self.with_trans {
            *m = m.multiplied_mat(&self.trans);
            *dm = dm.multiplied_mat(&self.trans);
            *d2m = d2m.multiplied_mat(&self.trans);
        }

        ok
    }

    /// OCCT NbIntervals (L195-220).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        let nb_sec = super::frenet::curve_nb_intervals(my_trimmed, s);
        let nb_law = TrihedronLaw::nb_intervals(self.my_law.borrow().as_ref(), s);

        if nb_sec == 1 {
            return nb_law;
        } else if nb_law == 1 {
            return nb_sec;
        }

        let int_c = super::frenet::curve_intervals(my_trimmed, s);
        let mut int_l = vec![0.0; nb_law + 1];
        TrihedronLaw::intervals(self.my_law.borrow().as_ref(), &mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        fuse_intervals(&int_c, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        inter.len() - 1
    }

    /// OCCT Intervals (L222-253).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        let nb_sec = super::frenet::curve_nb_intervals(my_trimmed, s);
        let nb_law = TrihedronLaw::nb_intervals(self.my_law.borrow().as_ref(), s);

        if nb_sec == 1 {
            TrihedronLaw::intervals(self.my_law.borrow().as_ref(), t, s);
            return;
        } else if nb_law == 1 {
            let disc = super::frenet::curve_intervals(my_trimmed, s);
            t[..disc.len()].copy_from_slice(&disc);
            return;
        }

        let int_c = super::frenet::curve_intervals(my_trimmed, s);
        let mut int_l = vec![0.0; nb_law + 1];
        TrihedronLaw::intervals(self.my_law.borrow().as_ref(), &mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        fuse_intervals(&int_c, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        // OCCT: for (ii = 1; ii <= Inter.Length(); ii++) T(ii) = Inter(ii).
        for (ii, value) in inter.iter().enumerate() {
            t[ii] = *value;
        }
    }

    /// OCCT SetInterval (L255-261).
    fn set_interval(&mut self, first: f64, last: f64) {
        self.my_law.borrow_mut().set_interval(first, last);
        // OCCT: myTrimmed = myCurve->Trim(First, Last, 0).
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.my_trimmed = Some(Curve3::Trimmed(TrimmedCurve3::new(
            my_curve.clone(),
            first,
            last,
        )));
    }

    /// OCCT GetInterval (L263-269).
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        *first = curve_first_parameter(my_trimmed);
        *last = curve_last_parameter(my_trimmed);
    }

    /// OCCT GetDomain (L271-280).
    fn get_domain(&self, first: &mut f64, last: &mut f64) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        *first = curve_first_parameter(my_curve);
        *last = curve_last_parameter(my_curve);
    }

    /// OCCT GetMaximalNorm (L282-287) — On suppose les triedre normee
    /// => return 1.
    fn get_maximal_norm(&self) -> f64 {
        1.0
    }

    /// OCCT GetAverageLaw (L289-309).
    fn get_average_law(&self, am: &mut GpMat, av: &mut DVec3) {
        let my_trimmed = self.my_trimmed.as_ref().expect("null myTrimmed");
        // OCCT: myLaw->GetAverageLaw(V1, V2, V3) — fills the member caches.
        let mut v1 = self.v1.get();
        let mut v2 = self.v2.get();
        let mut v3 = self.v3.get();
        self.my_law
            .borrow()
            .get_average_law(&mut v1, &mut v2, &mut v3);
        self.v1.set(v1);
        self.v2.set(v2);
        self.v3.set(v3);
        am.set_cols(self.v1.get(), self.v2.get(), self.v3.get());

        *av = DVec3::ZERO;
        let delta = (curve_last_parameter(my_trimmed) - curve_first_parameter(my_trimmed)) / 10.0;
        let mut u = curve_first_parameter(my_trimmed);
        for _ii in 0..=10 {
            // OCCT: V.SetXYZ(myTrimmed->Value(U).XYZ()); AV += V.
            *av += my_trimmed.point_at(u);
            u += delta;
        }
        *av /= 11.0;
    }

    /// OCCT IsTranslation (L311-323).
    fn is_translation(&self, error: &mut f64) -> bool {
        *error = 0.0;
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // OCCT: Type = myCurve->GetType(); if (Type == GeomAbs_Line) ...
        let is_line = matches!(my_curve, Curve3::Line(_));
        if is_line {
            return TrihedronLaw::is_constant(self.my_law.borrow().as_ref())
                || TrihedronLaw::is_only_by3d_curve(self.my_law.borrow().as_ref());
        }
        false
    }

    /// OCCT IsRotation (L325-337).
    fn is_rotation(&self, error: &mut f64) -> bool {
        *error = 0.0;
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        // OCCT: Type = myCurve->GetType(); if (Type == GeomAbs_Circle) ...
        let is_circle = matches!(my_curve, Curve3::Circle(_));
        if is_circle {
            return TrihedronLaw::is_only_by3d_curve(self.my_law.borrow().as_ref());
        }
        false
    }

    /// OCCT Rotation (L339-344).
    fn rotation(&self, centre: &mut DVec3) {
        // OCCT: Centre = myCurve->Circle().Location().
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        match my_curve {
            Curve3::Circle(circle) => *centre = circle.center,
            _ => panic!("Standard_NoSuchObject: GeomFill_CurveAndTrihedron::Rotation"),
        }
    }
}
