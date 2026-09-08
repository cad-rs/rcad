//! OCCT GeomFill_SweepFunction (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SweepFunction.hxx (members) + GeomFill_SweepFunction.cxx (whole
//! file L36-371).
//!
//! Architecture differences:
//! - The OCCT `handle(GeomFill_SectionLaw)` / `handle(GeomFill_LocationLaw)`
//!   members are shared, and the laws are later mutated in place through the
//!   same handles (SetInterval / SetTolerance) — the rcad carriers are
//!   `Rc<RefCell<...>>`.
//! - The scratch members M / V / DM / DV / D2M / D2V are written inside the
//!   const D0/D1/D2 evaluators — `Cell`s (same cache pattern as batch 1).
//! - The class derives Approx_SweepFunction; the rcad form implements the
//!   [`ApproxSweepFunction`] trait on the same struct.

use std::cell::{Cell, RefCell};
use std::rc::Rc;

use glam::{DVec2, DVec3};

use rcad_kernel::math::GeomAbsShape;

use super::gp_mat::GpMat;
use super::location_law::LocationLaw;
use super::section_law::SectionLaw;
use crate::geomalgo::geomfill::approx_sweep_function::ApproxSweepFunction;

/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-9;

/// OCCT GeomFill_SweepFunction (GeomFill_SweepFunction.hxx).
pub struct SweepFunction {
    /// OCCT handle(GeomFill_LocationLaw) myLoc.
    my_loc: Rc<RefCell<dyn LocationLaw>>,
    /// OCCT handle(GeomFill_SectionLaw) mySec.
    my_sec: Rc<RefCell<dyn SectionLaw>>,
    /// OCCT double myf.
    myf: f64,
    /// OCCT double myfOnS.
    myfons: f64,
    /// OCCT double myRatio.
    my_ratio: f64,
    /// OCCT gp_Mat M (D0/D1/D2 scratch).
    m: Cell<GpMat>,
    /// OCCT gp_Vec V (D0/D1/D2 scratch).
    v: Cell<DVec3>,
    /// OCCT gp_Mat DM (D1/D2 scratch).
    dm: Cell<GpMat>,
    /// OCCT gp_Vec DV (D1/D2 scratch).
    dv: Cell<DVec3>,
    /// OCCT gp_Mat D2M (D2 scratch).
    d2m: Cell<GpMat>,
    /// OCCT gp_Vec D2V (D2 scratch).
    d2v: Cell<DVec3>,
}

impl SweepFunction {
    /// OCCT GeomFill_SweepFunction::GeomFill_SweepFunction (L36-47).
    pub fn new(
        section: Rc<RefCell<dyn SectionLaw>>,
        location: Rc<RefCell<dyn LocationLaw>>,
        first_parameter: f64,
        first_parameter_on_s: f64,
        ratio_parameter_on_s: f64,
    ) -> Self {
        SweepFunction {
            my_loc: location,
            my_sec: section,
            myf: first_parameter,
            myfons: first_parameter_on_s,
            my_ratio: ratio_parameter_on_s,
            m: Cell::new(GpMat::identity()),
            v: Cell::new(DVec3::ZERO),
            dm: Cell::new(GpMat::identity()),
            dv: Cell::new(DVec3::ZERO),
            d2m: Cell::new(GpMat::identity()),
            d2v: Cell::new(DVec3::ZERO),
        }
    }

    /// OCCT GeomFill_SweepFunction::D0 (L51-81).
    pub fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        let t = self.myfons + (param - self.myf) * self.my_ratio;
        let l = poles.len();

        let mut m = self.m.get();
        let mut v = self.v.get();
        let ok = self
            .my_loc
            .borrow()
            .d0_2d(param, &mut m, &mut v, poles2d);
        self.m.set(m);
        self.v.set(v);
        if !ok {
            return ok;
        }
        let ok = self.my_sec.borrow().d0(t, poles, weigths);
        if !ok {
            return ok;
        }
        let m = self.m.get();
        let v = self.v.get();

        for ii in 1..=l {
            // OCCT: aux *= M; aux += V.XYZ() — the gp_XYZ row-vector form.
            poles[ii - 1] = GpMat::multiply_xyz_row(poles[ii - 1], &m) + v;
        }
        true
    }

    /// OCCT GeomFill_SweepFunction::D1 (L85-127).
    #[allow(clippy::too_many_arguments)]
    pub fn d1(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        let t = self.myfons + (param - self.myf) * self.my_ratio;
        let l = poles.len();

        let mut m = self.m.get();
        let mut v = self.v.get();
        let mut dm = self.dm.get();
        let mut dv = self.dv.get();
        let ok = self.my_loc.borrow().d1(
            param,
            &mut m,
            &mut v,
            &mut dm,
            &mut dv,
            poles2d,
            dpoles2d,
        );
        self.m.set(m);
        self.v.set(v);
        self.dm.set(dm);
        self.dv.set(dv);
        if !ok {
            return ok;
        }
        let ok = self
            .my_sec
            .borrow()
            .d1(t, poles, dpoles, weigths, dweigths);
        if !ok {
            return ok;
        }
        let m = self.m.get();
        let v = self.v.get();
        let dm = self.dm.get();
        let dv = self.dv.get();

        for ii in 1..=l {
            let mut pprim = dpoles[ii - 1];
            let mut p = poles[ii - 1];
            // PPrim *= myRatio; DWeigths(ii) *= myRatio.
            pprim *= self.my_ratio;
            dweigths[ii - 1] *= self.my_ratio;
            // PPrim *= M; PPrim += DM * P; PPrim += DV.XYZ().
            pprim = GpMat::multiply_xyz_row(pprim, &m);
            pprim += dm.multiplied_xyz(p);
            pprim += dv;
            dpoles[ii - 1] = pprim;

            // P *= M; P += V.XYZ().
            p = GpMat::multiply_xyz_row(p, &m) + v;
            poles[ii - 1] = p;
        }
        true
    }

    /// OCCT GeomFill_SweepFunction::D2 (L131-186).
    #[allow(clippy::too_many_arguments)]
    pub fn d2(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        let t = self.myfons + (param - self.myf) * self.my_ratio;
        let squareratio = self.my_ratio * self.my_ratio;
        let l = poles.len();

        let mut m = self.m.get();
        let mut v = self.v.get();
        let mut dm = self.dm.get();
        let mut dv = self.dv.get();
        let mut d2m = self.d2m.get();
        let mut d2v = self.d2v.get();
        let ok = self.my_loc.borrow().d2(
            param,
            &mut m,
            &mut v,
            &mut dm,
            &mut dv,
            &mut d2m,
            &mut d2v,
            poles2d,
            dpoles2d,
            d2poles2d,
        );
        self.m.set(m);
        self.v.set(v);
        self.dm.set(dm);
        self.dv.set(dv);
        self.d2m.set(d2m);
        self.d2v.set(d2v);
        if !ok {
            return ok;
        }
        let ok = self.my_sec.borrow().d2(
            t,
            poles,
            dpoles,
            d2poles,
            weigths,
            dweigths,
            d2weigths,
        );
        if !ok {
            return ok;
        }
        let m = self.m.get();
        let v = self.v.get();
        let dm = self.dm.get();
        let dv = self.dv.get();
        let d2m = self.d2m.get();
        let d2v = self.d2v.get();

        for ii in 1..=l {
            let mut psecn = d2poles[ii - 1];
            let mut pprim = dpoles[ii - 1];
            let mut p = poles[ii - 1];
            // PPrim *= myRatio; DWeigths(ii) *= myRatio; PSecn *=
            // squareratio; D2Weigths(ii) *= squareratio.
            pprim *= self.my_ratio;
            dweigths[ii - 1] *= self.my_ratio;
            psecn *= squareratio;
            d2weigths[ii - 1] *= squareratio;

            // PSecn *= M; PSecn += 2 * (DM * PPrim); PSecn += D2M * P;
            // PSecn += D2V.XYZ().
            psecn = GpMat::multiply_xyz_row(psecn, &m);
            psecn += 2.0 * dm.multiplied_xyz(pprim);
            psecn += d2m.multiplied_xyz(p);
            psecn += d2v;
            d2poles[ii - 1] = psecn;

            // PPrim *= M; PPrim += DM * P; PPrim += DV.XYZ().
            pprim = GpMat::multiply_xyz_row(pprim, &m);
            pprim += dm.multiplied_xyz(p);
            pprim += dv;
            dpoles[ii - 1] = pprim;

            // P *= M; P += V.XYZ().
            p = GpMat::multiply_xyz_row(p, &m) + v;
            poles[ii - 1] = p;
        }
        true
    }

    /// OCCT Nb2dCurves (L190-193).
    pub fn nb_2d_curves(&self) -> usize {
        self.my_loc.borrow().nb_2d_curves()
    }

    /// OCCT SectionShape (L197-200).
    pub fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize) {
        self.my_sec.borrow().section_shape(nb_poles, nb_knots, degree);
    }

    /// OCCT Knots (L204-207).
    pub fn knots(&self, t_knots: &mut [f64]) {
        self.my_sec.borrow().knots(t_knots);
    }

    /// OCCT Mults (L211-214).
    pub fn mults(&self, t_mults: &mut [i32]) {
        self.my_sec.borrow().mults(t_mults);
    }

    /// OCCT IsRational (L218-221).
    pub fn is_rational(&self) -> bool {
        self.my_sec.borrow().is_rational()
    }

    /// OCCT NbIntervals (L225-255).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let nb_sec = self.my_sec.borrow().nb_intervals(s);
        let nb_loc = self.my_loc.borrow().nb_intervals(s);

        if nb_sec == 1 {
            return nb_loc;
        } else if nb_loc == 1 {
            return nb_sec;
        }

        let mut int_s = vec![0.0; nb_sec + 1];
        let mut int_l = vec![0.0; nb_loc + 1];
        self.my_sec.borrow().intervals(&mut int_s, s);
        for ii in 1..=nb_sec + 1 {
            // T = (IntS(ii) - myfOnS) / myRatio + myf.
            let t = (int_s[ii - 1] - self.myfons) / self.my_ratio + self.myf;
            int_s[ii - 1] = t;
        }
        self.my_loc.borrow().intervals(&mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        rcad_kernel::base::geom_lib::fuse_intervals(&int_s, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        inter.len() - 1
    }

    /// OCCT Intervals (L259-300) — writes the caller-allocated T array in
    /// place.
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let nb_sec = self.my_sec.borrow().nb_intervals(s);
        let nb_loc = self.my_loc.borrow().nb_intervals(s);

        if nb_sec == 1 {
            // OCCT: myLoc->Intervals(T, S) — the rcad LocationLaw form is
            // Vec-based; the values are copied back into T.
            let mut tv = vec![0.0; nb_loc + 1];
            self.my_loc.borrow().intervals(&mut tv, s);
            t.copy_from_slice(&tv);
            return;
        } else if nb_loc == 1 {
            // OCCT: mySec->Intervals(T, S); the section bounds are mapped
            // back to the path parameter in place.
            self.my_sec.borrow().intervals(t, s);
            for ii in 1..=nb_sec + 1 {
                let tv = (t[ii - 1] - self.myfons) / self.my_ratio + self.myf;
                t[ii - 1] = tv;
            }
            return;
        }

        let mut int_s = vec![0.0; nb_sec + 1];
        let mut int_l = vec![0.0; nb_loc + 1];

        self.my_sec.borrow().intervals(&mut int_s, s);
        for ii in 1..=nb_sec + 1 {
            let tv = (int_s[ii - 1] - self.myfons) / self.my_ratio + self.myf;
            int_s[ii - 1] = tv;
        }
        self.my_loc.borrow().intervals(&mut int_l, s);

        let mut inter: Vec<f64> = Vec::new();
        rcad_kernel::base::geom_lib::fuse_intervals(&int_s, &int_l, &mut inter, P_CONFUSION * 0.99, true);
        // OCCT: for (ii = 1; ii <= Inter.Length(); ii++) T(ii) = Inter(ii).
        for (ii, value) in inter.iter().enumerate() {
            t[ii] = *value;
        }
    }

    /// OCCT SetInterval (L304-311).
    pub fn set_interval(&mut self, first: f64, last: f64) {
        self.my_loc.borrow_mut().set_interval(first, last);
        // OCCT literal: uf = myf + (First - myf) * myRatio (the myfOnS
        // offset is not applied — the OCCT form is kept).
        let uf = self.myf + (first - self.myf) * self.my_ratio;
        let ul = self.myf + (last - self.myf) * self.my_ratio;
        self.my_sec.borrow_mut().set_interval(uf, ul);
    }

    /// OCCT GetTolerance (L315-321).
    pub fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        self.my_sec
            .borrow()
            .get_tolerance(bound_tol, surf_tol, angle_tol, tol3d);
    }

    /// OCCT Resolution (L325-331).
    pub fn resolution(&self, index: usize, tol: f64, tolu: &mut f64, tolv: &mut f64) {
        self.my_loc.borrow().resolution(index, tol, tolu, tolv);
    }

    /// OCCT SetTolerance (L335-339).
    pub fn set_tolerance(&mut self, tol3d: f64, tol2d: f64) {
        self.my_sec.borrow_mut().set_tolerance(tol3d, tol2d);
        self.my_loc.borrow_mut().set_tolerance(tol3d, tol2d);
    }

    /// OCCT BarycentreOfSurf (L343-355).
    pub fn barycentre_of_surf(&self) -> DVec3 {
        let mut a_m = GpMat::identity();
        let mut translate = DVec3::ZERO;

        let mut bary = self.my_sec.borrow().barycentre_of_surf();
        self.my_loc.borrow().get_average_law(&mut a_m, &mut translate);
        // OCCT: Bary.ChangeCoord() *= aM; += Translate.XYZ().
        bary = GpMat::multiply_xyz_row(bary, &a_m) + translate;

        bary
    }

    /// OCCT MaximalSection (L359-364).
    pub fn maximal_section(&self) -> f64 {
        let l = self.my_sec.borrow().maximal_section();
        l * self.my_loc.borrow().get_maximal_norm()
    }

    /// OCCT GetMinimalWeight (L368-371).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        self.my_sec.borrow().get_minimal_weight(weigths);
    }
}

impl ApproxSweepFunction for SweepFunction {
    /// OCCT D0 — the trait entry over the inherent method.
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        SweepFunction::d0(self, param, first, last, poles, poles2d, weigths)
    }

    /// OCCT D1 — the trait entry over the inherent method.
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        SweepFunction::d1(
            self,
            param,
            first,
            last,
            poles,
            dpoles,
            poles2d,
            dpoles2d,
            weigths,
            dweigths,
        )
    }

    /// OCCT D2 — the trait entry over the inherent method.
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        SweepFunction::d2(
            self,
            param,
            first,
            last,
            poles,
            dpoles,
            d2poles,
            poles2d,
            dpoles2d,
            d2poles2d,
            weigths,
            dweigths,
            d2weigths,
        )
    }

    /// OCCT Nb2dCurves.
    fn nb_2d_curves(&self) -> usize {
        SweepFunction::nb_2d_curves(self)
    }

    /// OCCT SectionShape.
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize) {
        SweepFunction::section_shape(self, nb_poles, nb_knots, degree)
    }

    /// OCCT Knots.
    fn knots(&self, t_knots: &mut [f64]) {
        SweepFunction::knots(self, t_knots)
    }

    /// OCCT Mults.
    fn mults(&self, t_mults: &mut [i32]) {
        SweepFunction::mults(self, t_mults)
    }

    /// OCCT IsRational.
    fn is_rational(&self) -> bool {
        SweepFunction::is_rational(self)
    }

    /// OCCT NbIntervals.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        SweepFunction::nb_intervals(self, s)
    }

    /// OCCT Intervals.
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        SweepFunction::intervals(self, t, s)
    }

    /// OCCT SetInterval.
    fn set_interval(&mut self, first: f64, last: f64) {
        SweepFunction::set_interval(self, first, last)
    }

    /// OCCT Resolution.
    fn resolution(&self, index: usize, tol: f64, tolu: &mut f64, tolv: &mut f64) {
        SweepFunction::resolution(self, index, tol, tolu, tolv)
    }

    /// OCCT GetTolerance.
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        SweepFunction::get_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d)
    }

    /// OCCT SetTolerance.
    fn set_tolerance(&mut self, tol3d: f64, tol2d: f64) {
        SweepFunction::set_tolerance(self, tol3d, tol2d)
    }

    /// OCCT BarycentreOfSurf.
    fn barycentre_of_surf(&self) -> DVec3 {
        SweepFunction::barycentre_of_surf(self)
    }

    /// OCCT MaximalSection.
    fn maximal_section(&self) -> f64 {
        SweepFunction::maximal_section(self)
    }

    /// OCCT GetMinimalWeight.
    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        SweepFunction::get_minimal_weight(self, weigths)
    }
}
