//! OCCT AppDef_LinearCriteria (ModelingData/TKGeomBase/AppDef).
//!
//! 1:1 translation of `AppDef_LinearCriteria.hxx` (L17-118) and
//! `AppDef_LinearCriteria.cxx` (L38-857).
//!
//! Handle / architecture mapping:
//! - `occ::handle<FEmTool_Curve> myCurve` -> [`CurveHandle`]
//!   (`Rc<RefCell<Curve>>`), preserving the OCCT aliasing (the owner keeps
//!   mutating the curve between criteria calls);
//! - `occ::handle<FEmTool_ElementaryCriterion> myCriteria[3]` ->
//!   `[Option<Box<dyn ElementaryCriterion>>; 3]` (null handle == None,
//!   created by `SetCurve` exactly like OCCT);
//! - `occ::handle<NCollection_HArray1<double>> myParameters` ->
//!   `Option<RealArray1>` by value (stored once by `SetParameters`, only
//!   read afterwards);
//! - `occ::handle<NCollection_HArray1<double>> myCache` ->
//!   `Option<Vec<f64>>` (the OCCT `BasicValue` array view over the cache
//!   cells becomes a subslice of the same storage);
//! - `occ::handle<NCollection_HArray2<...>>` (assembly table) ->
//!   `fem_tool::AssemblingTable` (freshly built by `AssemblyTable()`).
//!
//! The OCCT static `order` (cxx L38-41) becomes the module-private
//! [`order`] function.  `AppDef_MyLineTool` statics map to the free
//! functions of `geomalgo::app_def::my_line_tool`.

use std::cell::RefCell;
use std::rc::Rc;

use glam::{DVec2, DVec3};

use rcad_kernel::math::fem_tool::elementary_criterion::{CoeffHandle, ElementaryCriterion};
use rcad_kernel::math::fem_tool::{
    AssemblingTable, IntArray2, LinearFlexion, LinearJerk, LinearTension,
};
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::GeomAbsShape;

use super::app_def::my_line_tool;
use super::app_def::MultiLine;
use super::app_def_smooth_criterion::{CurveHandle, RealArray1, SmoothCriterion};

/// OCCT static order (AppDef_LinearCriteria.cxx L38-41).
fn order(b: &rcad_kernel::math::plib::hermit_jacobi::HermitJacobi) -> i32 {
    b.niv_constr()
}

/// OCCT AppDef_LinearCriteria - defined an Linear Criteria to used in
/// variational Smoothing of points (AppDef_LinearCriteria.hxx L37-116).
pub struct LinearCriteria {
    /// OCCT mySSP.
    my_ssp: MultiLine,
    /// OCCT myParameters (null until SetParameters).
    my_parameters: Option<RealArray1>,
    /// OCCT myCache (null until first BuildCache).
    my_cache: Option<Vec<f64>>,
    /// OCCT myCriteria[3] (null handles until SetCurve).
    my_criteria: [Option<Box<dyn ElementaryCriterion>>; 3],
    /// OCCT myEstimation[3].
    my_estimation: [f64; 3],
    /// OCCT myQuadraticWeight.
    my_quadratic_weight: f64,
    /// OCCT myQualityWeight.
    my_quality_weight: f64,
    /// OCCT myPercent[3].
    my_percent: [f64; 3],
    /// OCCT myPntWeight.
    my_pnt_weight: RealArray1,
    /// OCCT myCurve (shared handle).
    my_curve: Option<CurveHandle>,
    /// OCCT myLength.
    my_length: f64,
    /// OCCT myE.
    my_e: i32,
    /// OCCT IF (`if` is a Rust keyword).
    if_: i32,
    /// OCCT IL.
    il: i32,
}

impl LinearCriteria {
    /// OCCT AppDef_LinearCriteria::AppDef_LinearCriteria (cxx L45-60).
    pub fn new(ssp: &MultiLine, first_point: i32, last_point: i32) -> Self {
        // memset(myEstimation, 0, sizeof(myEstimation));
        // memset(myPercent, 0, sizeof(myPercent));
        // myPntWeight.Init(1.);
        let mut my_pnt_weight = RealArray1::new(first_point, last_point);
        my_pnt_weight.init(1.0);

        LinearCriteria {
            my_ssp: ssp.clone(),
            my_parameters: None,
            my_cache: None,
            my_criteria: [None, None, None],
            my_estimation: [0.0; 3],
            my_quadratic_weight: 0.0,
            my_quality_weight: 0.0,
            my_percent: [0.0; 3],
            my_pnt_weight,
            my_curve: None,
            my_length: 0.0,
            my_e: 0,
            if_: 0,
            il: 0,
        }
    }

    /// OCCT AppDef_LinearCriteria::BuildCache (cxx L809-857).
    fn build_cache(&mut self, element: i32) {
        let curve = self
            .my_curve
            .as_ref()
            .expect("AppDef_LinearCriteria: myCurve is null")
            .clone();
        let params = self
            .my_parameters
            .as_ref()
            .expect("AppDef_LinearCriteria: myParameters is null")
            .clone();

        // UFirst = myCurve->Knots()(Element); ULast = myCurve->Knots()(Element + 1);
        // (the Rust Knots() accessor requires &mut self, hence RefMut).
        let (ufirst, ulast) = {
            let mut mc = curve.borrow_mut();
            (mc.knots()[(element - 1) as usize], mc.knots()[element as usize])
        };

        self.if_ = 0;
        for ipnt in params.lower()..=params.upper() {
            let t = params.value(ipnt);
            if (t > ufirst && t <= ulast) || (element == 1 && t == ufirst) {
                if self.if_ == 0 {
                    self.if_ = ipnt;
                }
                self.il = ipnt;
            } else if t > ulast {
                break;
            }
        }

        if self.if_ != 0 {
            let order = curve.borrow().base().work_degree() + 1;
            self.my_cache = Some(vec![0.0; ((self.il - self.if_ + 1) * order) as usize]);

            // const PLib_HermitJacobi& myBase = myCurve->Base(); (cxx L837)
            let mc = curve.borrow();
            let my_base = mc.base();

            let mut ii = 1i32;
            let mut ipnt = self.if_;
            while ipnt <= self.il {
                // OCCT: cache = &myCache->ChangeValue(ii);
                //       NCollection_Array1<double> BasicValue(cache[0], 0, order - 1);
                // (a view over `order` cache cells starting at ii).
                let t = params.value(ipnt);
                let coeff = 2.0 / (ulast - ufirst);
                let c0 = -(ulast + ufirst) / 2.0;
                let s = (t + c0) * coeff;
                let cache = self.my_cache.as_mut().unwrap();
                my_base.d0(s, &mut cache[(ii - 1) as usize..((ii - 1) + order) as usize]);
                ii += order;
                ipnt += 1;
            }
        } else {
            // no points in the interval.
            self.if_ = self.il;
            self.il -= 1;
        }
        self.my_e = element;
    }
}

impl SmoothCriterion for LinearCriteria {
    /// OCCT AppDef_LinearCriteria::SetParameters (cxx L64-69).
    fn set_parameters(&mut self, parameters: &RealArray1) {
        // myParameters = Parameters; (handle store -> value copy; OCCT does
        // not alias-mutate the array afterwards).
        self.my_parameters = Some(parameters.clone());
        self.my_e = 0; // Cache become invalid.
    }

    /// OCCT AppDef_LinearCriteria::SetCurve (cxx L73-156).
    fn set_curve(&mut self, c: &CurveHandle) {
        if self.my_curve.is_none() {
            self.my_curve = Some(c.clone());

            let curve = self.my_curve.clone().unwrap();
            let (mx_deg, nb_dim, order_v) = {
                let mc = curve.borrow();
                (mc.base().work_degree(), mc.dimension(), order(mc.base()))
            };

            let mut constraint_order = GeomAbsShape::C0;
            match order_v {
                0 => constraint_order = GeomAbsShape::C0,
                1 => constraint_order = GeomAbsShape::C1,
                2 => constraint_order = GeomAbsShape::C2,
                _ => {} // OCCT: no default case, ConstraintOrder stays C0.
            }

            self.my_criteria[0] = Some(Box::new(LinearTension::new(mx_deg, constraint_order)));
            self.my_criteria[1] = Some(Box::new(LinearFlexion::new(mx_deg, constraint_order)));
            self.my_criteria[2] = Some(Box::new(LinearJerk::new(mx_deg, constraint_order)));

            let coeff: CoeffHandle = Rc::new(RefCell::new(Matrix::new(0, 0, 1, nb_dim)));

            self.my_criteria[0].as_mut().unwrap().set_coeff(&coeff);
            self.my_criteria[1].as_mut().unwrap().set_coeff(&coeff);
            self.my_criteria[2].as_mut().unwrap().set_coeff(&coeff);
        } else if !Rc::ptr_eq(self.my_curve.as_ref().unwrap(), c) {
            let (old_mx_deg, old_nb_dim, old_order) = {
                let mc = self.my_curve.as_ref().unwrap().borrow();
                (mc.base().work_degree(), mc.dimension(), order(mc.base()))
            };

            self.my_curve = Some(c.clone());

            let (mx_deg, nb_dim, order_v) = {
                let mc = self.my_curve.as_ref().unwrap().borrow();
                (mc.base().work_degree(), mc.dimension(), order(mc.base()))
            };

            if mx_deg != old_mx_deg || order_v != old_order {
                let mut constraint_order = GeomAbsShape::C0;
                match order_v {
                    0 => constraint_order = GeomAbsShape::C0,
                    1 => constraint_order = GeomAbsShape::C1,
                    2 => constraint_order = GeomAbsShape::C2,
                    _ => {} // OCCT: no default case, ConstraintOrder stays C0.
                }

                self.my_criteria[0] = Some(Box::new(LinearTension::new(mx_deg, constraint_order)));
                self.my_criteria[1] = Some(Box::new(LinearFlexion::new(mx_deg, constraint_order)));
                self.my_criteria[2] = Some(Box::new(LinearJerk::new(mx_deg, constraint_order)));

                let coeff: CoeffHandle = Rc::new(RefCell::new(Matrix::new(0, 0, 1, nb_dim)));

                self.my_criteria[0].as_mut().unwrap().set_coeff(&coeff);
                self.my_criteria[1].as_mut().unwrap().set_coeff(&coeff);
                self.my_criteria[2].as_mut().unwrap().set_coeff(&coeff);
            } else if nb_dim != old_nb_dim {
                let coeff: CoeffHandle = Rc::new(RefCell::new(Matrix::new(0, 0, 1, nb_dim)));

                self.my_criteria[0].as_mut().unwrap().set_coeff(&coeff);
                self.my_criteria[1].as_mut().unwrap().set_coeff(&coeff);
                self.my_criteria[2].as_mut().unwrap().set_coeff(&coeff);
            }
        }
    }

    /// OCCT AppDef_LinearCriteria::GetCurve (cxx L160-163).
    fn get_curve(&self, c: &mut Option<CurveHandle>) {
        *c = self.my_curve.clone();
    }

    /// OCCT AppDef_LinearCriteria::SetEstimation (cxx L167-172).
    fn set_estimation(&mut self, e1: f64, e2: f64, e3: f64) {
        self.my_estimation[0] = e1;
        self.my_estimation[1] = e2;
        self.my_estimation[2] = e3;
    }

    /// OCCT AppDef_LinearCriteria::EstLength (cxx L174-177).
    fn est_length(&mut self) -> &mut f64 {
        &mut self.my_length
    }

    /// OCCT AppDef_LinearCriteria::GetEstimation (cxx L181-186).
    fn get_estimation(&self, e1: &mut f64, e2: &mut f64, e3: &mut f64) {
        *e1 = self.my_estimation[0];
        *e2 = self.my_estimation[1];
        *e3 = self.my_estimation[2];
    }

    /// OCCT AppDef_LinearCriteria::AssemblyTable (cxx L190-273).
    fn assembly_table(&self) -> AssemblingTable {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::AssemblyTable")
            .clone();

        let (nb_dim, nb_elm, nc1, mx_deg) = {
            let mc = curve.borrow();
            (
                mc.dimension(),
                mc.nb_elements(),
                order(mc.base()) + 1,
                mc.base().work_degree(),
            )
        };

        let mut ass_table = AssemblingTable::new(1, nb_dim, 1, nb_elm, 0, mx_deg);

        let mut nb_glob_var = 0i32;
        let mut gi0: i32;

        // For dim = 1
        // For first element (el = 1)
        {
            let glob_index = ass_table.change_value(1, 1);

            for i in 0..nc1 {
                nb_glob_var += 1;
                glob_index.set_value(i, nb_glob_var);
            }
            gi0 = mx_deg - 2 * nc1 + 1;
            for i in nc1..(2 * nc1) {
                nb_glob_var += 1;
                glob_index.set_value(i, nb_glob_var + gi0);
            }
            for i in (2 * nc1)..=mx_deg {
                nb_glob_var += 1;
                glob_index.set_value(i, nb_glob_var - nc1);
            }
            gi0 = nb_glob_var - nc1 + 1;
        }

        // For rest elements
        for el in 2..=nb_elm {
            let glob_index = ass_table.change_value(1, el);

            for i in 0..nc1 {
                glob_index.set_value(i, gi0 + i);
            }

            gi0 = mx_deg - 2 * nc1 + 1;
            for i in nc1..(2 * nc1) {
                nb_glob_var += 1;
                glob_index.set_value(i, nb_glob_var + gi0);
            }
            for i in (2 * nc1)..=mx_deg {
                nb_glob_var += 1;
                glob_index.set_value(i, nb_glob_var - nc1);
            }
            gi0 = nb_glob_var - nc1 + 1;
        }

        // For other dimensions
        let gi0 = nb_glob_var;
        for dim in 2..=nb_dim {
            for el in 1..=nb_elm {
                // Aux = AssTable->Value(1, el);
                let aux: Vec<i32> = {
                    let a = ass_table.value(1, el);
                    (0..=mx_deg).map(|i| a.value(i)).collect()
                };
                let glob_index = ass_table.change_value(dim, el);
                for i in 0..=mx_deg {
                    glob_index.set_value(i, aux[(i) as usize] + nb_glob_var);
                }
            }
            nb_glob_var += gi0;
        }

        ass_table
    }

    /// OCCT AppDef_LinearCriteria::DependenceTable (cxx L277-294).
    fn dependence_table(&self) -> IntArray2 {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::DependenceTable")
            .clone();

        let dim = curve.borrow().dimension();

        let mut dep_tab = IntArray2::new_init(1, dim, 1, dim, 0);
        for i in 1..=dim {
            dep_tab.set_value(i, i, 1);
        }

        dep_tab
    }

    /// OCCT AppDef_LinearCriteria::QualityValues (cxx L298-440).
    fn quality_values(
        &mut self,
        j1min: f64,
        j2min: f64,
        j3min: f64,
        j1: &mut f64,
        j2: &mut f64,
        j3: &mut f64,
    ) -> i32 {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::QualityValues")
            .clone();

        let (nb_dim, nb_elm) = {
            let mc = curve.borrow();
            (mc.dimension(), mc.nb_elements())
        };

        // NCollection_Array1<double>& Knots = myCurve->Knots(); - read
        // element-wise below (the array is not mutated in the loop).

        let mut deg = 0i32;
        let mut coeff: Option<CoeffHandle> = None;

        *j1 = 0.0;
        *j2 = 0.0;
        *j3 = 0.0;
        for el in 1..=nb_elm {
            let (curdeg, ufirst, ulast) = {
                let mut mc = curve.borrow_mut();
                let deg_i = mc.degree(el);
                let ufirst = mc.knots()[(el - 1) as usize];
                let ulast = mc.knots()[el as usize];
                (deg_i, ufirst, ulast)
            };

            if deg != curdeg {
                deg = curdeg;
                coeff = Some(Rc::new(RefCell::new(Matrix::new(0, deg, 1, nb_dim))));
            }
            let coeff = coeff.as_ref().unwrap();

            {
                let mut cm = coeff.borrow_mut();
                curve.borrow_mut().get_element(el, &mut cm);
            }

            {
                let c = self.my_criteria[0].as_mut().unwrap();
                c.set_coeff(coeff);
                c.set_knots(ufirst, ulast);
                *j1 += c.value();
            }

            {
                let c = self.my_criteria[1].as_mut().unwrap();
                c.set_coeff(coeff);
                c.set_knots(ufirst, ulast);
                *j2 += c.value();
            }

            {
                let c = self.my_criteria[2].as_mut().unwrap();
                c.set_coeff(coeff);
                c.set_knots(ufirst, ulast);
                *j3 += c.value();
            }
        }

        // Calculation of ICDANA - see MOTEST.f
        //  double JEsMin[3] = {.01, .001, .001}; // from MOTLIS.f
        let jes_min = [j1min, j2min, j3min];
        let val_cri = [*j1, *j2, *j3];

        let mut icdana = 0i32;

        //   (2) Test the improvement of estimates
        //       (overestimated criterion => No minimization)

        for i in 0..=2usize {
            if (val_cri[i] < 0.8 * self.my_estimation[i]) && (self.my_estimation[i] > jes_min[i]) {
                if icdana < 1 {
                    icdana = 1;
                }
                if val_cri[i] < 0.1 * self.my_estimation[i] {
                    icdana = 2;
                }
                self.my_estimation[i] = (1.05 * val_cri[i]).max(jes_min[i]);
            }
        }

        //  (3) Update the Estimates
        //     (underestimated criterion => poor conditioning)
        if val_cri[0] > self.my_estimation[0] * 2.0 {
            self.my_estimation[0] += val_cri[0] * 0.1;
            if icdana == 0 {
                if val_cri[0] > self.my_estimation[0] * 10.0 {
                    icdana = 2;
                } else {
                    icdana = 1;
                }
            } else {
                icdana = 2;
            }
        }
        if val_cri[1] > self.my_estimation[1] * 20.0 {
            self.my_estimation[1] += val_cri[1] * 0.1;
            if icdana == 0 {
                if val_cri[1] > self.my_estimation[1] * 100.0 {
                    icdana = 2;
                } else {
                    icdana = 1;
                }
            } else {
                icdana = 2;
            }
        }
        if val_cri[2] > self.my_estimation[2] * 20.0 {
            self.my_estimation[2] += val_cri[2] * 0.05;
            if icdana == 0 {
                if val_cri[2] > self.my_estimation[2] * 100.0 {
                    icdana = 2;
                } else {
                    icdana = 1;
                }
            } else {
                icdana = 2;
            }
        }

        icdana
    }

    /// OCCT AppDef_LinearCriteria::ErrorValues (cxx L444-510).
    fn error_values(
        &mut self,
        max_error: &mut f64,
        quadratic_error: &mut f64,
        average_error: &mut f64,
    ) {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::ErrorValues")
            .clone();

        let nb_dim = curve.borrow().dimension();

        let my_nb_p2d = my_line_tool::nb_p2d(&self.my_ssp) as i32;
        let my_nb_p3d = my_line_tool::nb_p3d(&self.my_ssp) as i32;

        if nb_dim != (2 * my_nb_p2d + 3 * my_nb_p3d) {
            panic!("Standard_DomainError: AppDef_LinearCriteria::ErrorValues");
        }

        let mut tab_p3d = vec![DVec3::ZERO; 1.max(my_nb_p3d) as usize];
        let mut tab_p2d = vec![DVec2::ZERO; 1.max(my_nb_p2d) as usize];
        let mut base_point = vec![0.0f64; nb_dim as usize];

        let params = self
            .my_parameters
            .as_ref()
            .expect("AppDef_LinearCriteria: myParameters is null");

        let mut c0;
        *max_error = 0.0;
        *quadratic_error = 0.0;
        *average_error = 0.0;

        for i in params.lower()..=params.upper() {
            curve.borrow_mut().d0(params.value(i), &mut base_point);

            c0 = 0;
            my_line_tool::value_3d(&self.my_ssp, i as usize, &mut tab_p3d);
            for ipnt in 1..=my_nb_p3d {
                // P3d.SetCoord(BasePoint(c0 + 1), BasePoint(c0 + 2), BasePoint(c0 + 3));
                let p3d = DVec3::new(
                    base_point[c0 as usize],
                    base_point[(c0 + 1) as usize],
                    base_point[(c0 + 2) as usize],
                );
                let sqr_dist = p3d.distance_squared(tab_p3d[(ipnt - 1) as usize]);
                let dist = sqr_dist.sqrt();
                *max_error = (*max_error).max(dist);
                *quadratic_error += sqr_dist;
                *average_error += dist;
                c0 += 3;
            }

            if my_nb_p3d == 0 {
                my_line_tool::value_2d(&self.my_ssp, i as usize, &mut tab_p2d);
            } else {
                my_line_tool::value_3d_2d(&self.my_ssp, i as usize, &mut tab_p3d, &mut tab_p2d);
            }
            for ipnt in 1..=my_nb_p2d {
                // P2d.SetCoord(BasePoint(c0 + 1), BasePoint(c0 + 2));
                let p2d = DVec2::new(base_point[c0 as usize], base_point[(c0 + 1) as usize]);
                let sqr_dist = p2d.distance_squared(tab_p2d[(ipnt - 1) as usize]);
                let dist = sqr_dist.sqrt();
                *max_error = (*max_error).max(dist);
                *quadratic_error += sqr_dist;
                *average_error += dist;
                c0 += 2;
            }
        }
    }

    /// OCCT AppDef_LinearCriteria::Hessian (cxx L514-619).
    fn hessian(&mut self, element: i32, dimension1: i32, dimension2: i32, h: &mut Matrix) {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::Hessian")
            .clone();

        if self.dependence_table().value(dimension1, dimension2) == 0 {
            panic!("Standard_DomainError: AppDef_LinearCriteria::Hessian");
        }

        let (mx_deg, order_v, ufirst, ulast) = {
            let mut mc = curve.borrow_mut();
            let mx_deg = mc.base().work_degree();
            let order_v = order(mc.base());
            let ufirst = mc.knots()[(element - 1) as usize];
            let ulast = mc.knots()[element as usize];
            (mx_deg, order_v, ufirst, ulast)
        };

        // math_Matrix AuxH(0, H.RowNumber() - 1, 0, H.ColNumber() - 1, 0.);
        let mut aux_h = Matrix::new(0, h.row_number() - 1, 0, h.col_number() - 1);
        aux_h.init(0.0);

        // Quality criterion part of Hessian

        h.init(0.0);

        for icrit in 0..=2usize {
            // (myQualityWeight * myPercent[icrit] / myEstimation[icrit])
            // is read before the criteria borrow (the criteria calls do not
            // touch the weights).
            let scale =
                self.my_quality_weight * self.my_percent[icrit] / self.my_estimation[icrit];

            let c = self.my_criteria[icrit].as_mut().unwrap();
            c.set_knots(ufirst, ulast);
            c.hessian(dimension1, dimension2, &mut aux_h);
            // H += scale * AuxH
            for i in h.lower_row()..=h.upper_row() {
                for j in h.lower_col()..=h.upper_col() {
                    let v = h.get(i, j) + scale * aux_h.get(i, j);
                    h.set(i, j, v);
                }
            }
        }

        // Least square part of Hessian

        aux_h.init(0.0);

        let coeff = (ulast - ufirst) / 2.0;
        let degh = 2 * order_v + 1;

        let i0 = h.lower_row();
        let j0 = h.lower_col();
        let params = self
            .my_parameters
            .as_ref()
            .expect("AppDef_LinearCriteria: myParameters is null");
        let di = self.my_pnt_weight.lower() - params.lower();

        // BuilCache
        if self.my_e != element {
            self.build_cache(element);
        }

        // Compute the least square Hessian
        let mut ii = 1i32;
        let mut ipnt = self.if_;
        while ipnt <= self.il {
            let poid = self.my_pnt_weight.value(di + ipnt) * 2.0;
            // const double* BV = &myCache->Value(ii); (BV[i] is cache cell ii-1+i)
            let bv = &self.my_cache.as_ref().unwrap()[(ii - 1) as usize..];

            // Hermite*Hermite part of matrix
            for i in 0..=degh {
                let k1 = if i <= order_v { i } else { i - order_v - 1 };
                let curcoeff = coeff.powi(k1) * poid * bv[i as usize];

                // Hermite*Hermite part of matrix
                for j in i..=degh {
                    let k2 = if j <= order_v { j } else { j - order_v - 1 };
                    let v = aux_h.get(i, j) + curcoeff * coeff.powi(k2) * bv[j as usize];
                    aux_h.set(i, j, v);
                }
                // Hermite*Jacobi part of matrix
                for j in (degh + 1)..=mx_deg {
                    let v = aux_h.get(i, j) + curcoeff * bv[j as usize];
                    aux_h.set(i, j, v);
                }
            }

            // Jacoby*Jacobi part of matrix
            for i in (degh + 1)..=mx_deg {
                let curcoeff = bv[i as usize] * poid;
                for j in i..=mx_deg {
                    let v = aux_h.get(i, j) + curcoeff * bv[j as usize];
                    aux_h.set(i, j, v);
                }
            }

            ipnt += 1;
            ii += mx_deg + 1;
        }

        let mut i1 = i0;
        for i in 0..=mx_deg {
            let mut j1 = j0 + i;
            for j in i..=mx_deg {
                let v = h.get(i1, j1) + self.my_quadratic_weight * aux_h.get(i, j);
                h.set(i1, j1, v);
                let v = h.get(i1, j1);
                h.set(j1, i1, v);
                j1 += 1;
            }
            i1 += 1;
        }
    }

    /// OCCT AppDef_LinearCriteria::Gradient (cxx L623-733).
    fn gradient(&mut self, element: i32, dimension: i32, g: &mut Vector) {
        let curve = self
            .my_curve
            .as_ref()
            .expect("Standard_DomainError: AppDef_LinearCriteria::ErrorValues")
            .clone();

        let my_nb_p2d = my_line_tool::nb_p2d(&self.my_ssp) as i32;
        let my_nb_p3d = my_line_tool::nb_p3d(&self.my_ssp) as i32;

        if dimension > (2 * my_nb_p2d + 3 * my_nb_p3d) {
            panic!("Standard_DomainError: AppDef_LinearCriteria::ErrorValues");
        }

        let mut tab_p3d = vec![DVec3::ZERO; 1.max(my_nb_p3d) as usize];
        let mut tab_p2d = vec![DVec2::ZERO; 1.max(my_nb_p2d) as usize];

        let in3d;
        let mut ind_pnt;
        let mut ind_crd;

        if dimension <= 3 * my_nb_p3d {
            in3d = true;
            ind_crd = dimension % 3;
            ind_pnt = dimension / 3;
            if ind_crd == 0 {
                ind_crd = 3;
            } else {
                ind_pnt += 1;
            }
        } else {
            in3d = false;
            ind_crd = (dimension - 3 * my_nb_p3d) % 2;
            ind_pnt = (dimension - 3 * my_nb_p3d) / 2;
            if ind_crd == 0 {
                ind_crd = 2;
            } else {
                ind_pnt += 1;
            }
        }

        let (ufirst, ulast, order_v, mx_deg) = {
            let mut mc = curve.borrow_mut();
            let ufirst = mc.knots()[(element - 1) as usize];
            let ulast = mc.knots()[element as usize];
            let order_v = order(mc.base());
            let mx_deg = mc.base().work_degree();
            (ufirst, ulast, order_v, mx_deg)
        };
        let coeff = (ulast - ufirst) / 2.0;

        let degh = 2 * order_v + 1;
        let i0 = g.lower();
        let params = self
            .my_parameters
            .as_ref()
            .expect("AppDef_LinearCriteria: myParameters is null");
        let di = self.my_pnt_weight.lower() - params.lower();

        if self.my_e != element {
            self.build_cache(element);
        }

        // G.Init(0.)
        for r in g.lower()..=g.upper() {
            g.set(r, 0.0);
        }

        let bv = &self.my_cache.as_ref().unwrap()[..];

        let mut ii = 1i32;
        let mut ipnt = self.if_;
        while ipnt <= self.il {
            let pnt;
            if in3d {
                my_line_tool::value_3d(&self.my_ssp, ipnt as usize, &mut tab_p3d);
                pnt = tab_p3d[(ind_pnt - 1) as usize][(ind_crd - 1) as usize];
            } else {
                if my_nb_p3d == 0 {
                    my_line_tool::value_2d(&self.my_ssp, ipnt as usize, &mut tab_p2d);
                } else {
                    my_line_tool::value_3d_2d(
                        &self.my_ssp,
                        ipnt as usize,
                        &mut tab_p3d,
                        &mut tab_p2d,
                    );
                }
                pnt = tab_p2d[(ind_pnt - 1) as usize][(ind_crd - 1) as usize];
            }

            let curcoeff = pnt * self.my_pnt_weight.value(di + ipnt);
            for i in 0..=mx_deg {
                let v = g.get(i0 + i) + bv[(ii - 1) as usize] * curcoeff;
                g.set(i0 + i, v);
                ii += 1;
            }
            ipnt += 1;
        }

        // G *= 2. * myQuadraticWeight
        for r in g.lower()..=g.upper() {
            let v = g.get(r) * 2.0 * self.my_quadratic_weight;
            g.set(r, v);
        }

        for i in 0..=degh {
            let k = if i <= order_v { i } else { i - order_v - 1 };
            let curcoeff = coeff.powi(k);
            let v = g.get(i0 + i) * curcoeff;
            g.set(i0 + i, v);
        }
    }

    /// OCCT AppDef_LinearCriteria::InputVector (cxx L737-763).
    fn input_vector(&mut self, x: &Vector, ass_table: &AssemblingTable) {
        let curve = self.my_curve.as_ref().unwrap().clone();

        let (nb_dim, nb_elm, mx_deg) = {
            let mc = curve.borrow();
            (mc.dimension(), mc.nb_elements(), mc.base().work_degree())
        };
        let mut coeff_el = Matrix::new(0, mx_deg, 1, nb_dim);

        let i0 = x.lower() - 1;

        for el in 1..=nb_elm {
            for dim in 1..=nb_dim {
                let glob_index = ass_table.value(dim, el);
                for i in 0..=mx_deg {
                    coeff_el.set(i, dim, x.get(i0 + glob_index.value(i)));
                }
            }
            curve.borrow_mut().set_degree(el, mx_deg);
            curve.borrow_mut().set_element(el, &coeff_el);
        }
    }

    /// OCCT AppDef_LinearCriteria::SetWeight (cxx L767-789).
    fn set_weight(
        &mut self,
        quadratic_weight: f64,
        quality_weight: f64,
        percent_j1: f64,
        percent_j2: f64,
        percent_j3: f64,
    ) {
        if quadratic_weight < 0.0 || quality_weight < 0.0 {
            panic!("Standard_DomainError: AppDef_LinearCriteria::SetWeight");
        }
        if percent_j1 < 0.0 || percent_j2 < 0.0 || percent_j3 < 0.0 {
            panic!("Standard_DomainError: AppDef_LinearCriteria::SetWeight");
        }

        self.my_quadratic_weight = quadratic_weight;
        self.my_quality_weight = quality_weight;

        let total = percent_j1 + percent_j2 + percent_j3;
        self.my_percent[0] = percent_j1 / total;
        self.my_percent[1] = percent_j2 / total;
        self.my_percent[2] = percent_j3 / total;
    }

    /// OCCT AppDef_LinearCriteria::GetWeight (cxx L793-798).
    fn get_weight(&self, quadratic_weight: &mut f64, quality_weight: &mut f64) {
        *quadratic_weight = self.my_quadratic_weight;
        *quality_weight = self.my_quality_weight;
    }

    /// OCCT AppDef_LinearCriteria::SetWeight(const NCollection_Array1<double>&)
    /// (cxx L802-805).  The OCCT Array1 assignment copies bounds and data
    /// (NCollection_Array1.hxx L273/L325) -> the Rust clone.
    fn set_weight_array(&mut self, weight: &RealArray1) {
        self.my_pnt_weight = weight.clone();
    }
}

#[cfg(test)]
mod tests {
    use super::super::app_def::{MultiLine, MultiPointConstraint};
    use rcad_kernel::math::fem_tool::Curve;
    use super::*;

    /// Shared fixture: a 1-dimensional, 1-element FEmTool curve on the
    /// HermitJacobi(2, C0) base (WorkDegree = 2, NivConstr = 0), knots
    /// [First, Last], one line of 3d points.
    fn make_criteria(
        first_knot: f64,
        last_knot: f64,
        params: &[f64],
        points: &[DVec3],
    ) -> LinearCriteria {
        let ml = MultiLine::new_tab_p3d(points);
        let mut crit = LinearCriteria::new(&ml, 1, points.len() as i32);

        let base = rcad_kernel::math::plib::hermit_jacobi::HermitJacobi::new(
            2,
            GeomAbsShape::C0,
        );
        let mut curve = Curve::new(1, 1, &base, 0.0);
        curve.knots()[0] = first_knot;
        curve.knots()[1] = last_knot;
        let curve: CurveHandle = Rc::new(RefCell::new(curve));

        let mut parameters = RealArray1::new(1, params.len() as i32);
        for (i, p) in params.iter().enumerate() {
            parameters.set_value(1 + i as i32, *p);
        }

        crit.set_parameters(&parameters);
        crit.set_curve(&curve);
        crit
    }

    /// AssemblyTable global index map, hand-derived for a 2-element,
    /// 1-dimension curve on HermitJacobi(2, C0): nc1 = Order+1 = 1,
    /// MxDeg = 2.
    ///
    /// OCCT cxx L206-270 by hand:
    /// el=1: block1 (i=0): G(0)=1; gi0 = 2-2+1 = 1;
    ///       block2 (i=1): G(1) = 2 + 1 = 3;
    ///       block3 (i=2): G(2) = 3 - 1 = 2;  gi0 = 3-1+1 = 3.
    /// el=2: block1 (i=0): G(0)=3; gi0 = 1;
    ///       block2 (i=1): G(1) = 4 + 1 = 5;
    ///       block3 (i=2): G(2) = 5 - 1 = 4;  gi0 = 5.
    /// i.e. the shared end value of el1 (index 3) equals the start value of
    /// el2, the two interior Jacobi dofs (2, 4) are private, NbGlobVar = 5.
    #[test]
    fn linear_criteria_assembly_table_hand() {
        let crit = make_criteria(0.0, 1.0, &[0.0], &[DVec3::ZERO, DVec3::ZERO]);
        // Replace the curve by a 2-element one (SetCurve first-call path
        // already created the criteria; the second call just rebinds).
        let base = rcad_kernel::math::plib::hermit_jacobi::HermitJacobi::new(2, GeomAbsShape::C0);
        let mut curve2 = Curve::new(1, 2, &base, 0.0);
        curve2.knots()[0] = 0.0;
        curve2.knots()[1] = 1.0;
        curve2.knots()[2] = 2.0;
        let curve2: CurveHandle = Rc::new(RefCell::new(curve2));
        let mut crit = crit;
        crit.set_curve(&curve2);

        let table = crit.assembly_table();
        let g1 = table.value(1, 1);
        let g2 = table.value(1, 2);
        assert_eq!(
            [g1.value(0), g1.value(1), g1.value(2)],
            [1, 3, 2]
        );
        assert_eq!(
            [g2.value(0), g2.value(1), g2.value(2)],
            [3, 5, 4]
        );

        // DependenceTable of a 1-dimension curve is the 1x1 identity.
        let dep = crit.dependence_table();
        assert_eq!(dep.value(1, 1), 1);
        assert_eq!(dep.col_length(), 1);
    }

    /// Gradient of the least-square part on a degenerate configuration,
    /// hand-derived.  One element on [0, 2] (HermitJacobi(2, C0): basis =
    /// [H00, H01, J0], J0 = TN0*(1-s^2) vanishes at s = +-1), pass points
    /// at u = 0 and u = 2 with x-coordinates 1 and 3, point weight 1,
    /// QuadraticWeight 1.  BuildCache maps u=0 -> s=-1, u=2 -> s=+1, so the
    /// basis rows are BV1 = [1,0,0] and BV2 = [0,1,0].  The OCCT statements
    /// give G_j = 2*QuadWeight*sum_i w_i * Pnt_i * BV_i[j], and the final
    /// Hermite scaling multiplies by coeff^k with coeff = (2-0)/2 = 1 and
    /// k = 0 (C0), i.e. G = [2*1, 2*3, 0].
    #[test]
    fn linear_criteria_gradient_hand() {
        let mut crit = make_criteria(
            0.0,
            2.0,
            &[0.0, 2.0],
            &[DVec3::new(1.0, 7.0, 8.0), DVec3::new(3.0, 7.0, 8.0)],
        );
        crit.set_weight(1.0, 0.0, 1.0, 1.0, 1.0);

        let mut g = Vector::new(0, 2);
        crit.gradient(1, 1, &mut g);

        assert!((g.get(0) - 2.0).abs() < 1e-12, "g0 = {}", g.get(0));
        assert!((g.get(1) - 6.0).abs() < 1e-12, "g1 = {}", g.get(1));
        assert!(g.get(2).abs() < 1e-12, "g2 = {}", g.get(2));

        // Same configuration with knots [0, 4]: coeff = (ULast-UFirst)/2 = 2
        // and the u->s map keeps BV = [1,0,0] / [0,1,0]; the k-mapping for
        // C0 gives k(0) = k(1) = 0, so the values are unchanged.
        let mut crit = make_criteria(
            0.0,
            4.0,
            &[0.0, 4.0],
            &[DVec3::new(1.0, 7.0, 8.0), DVec3::new(3.0, 7.0, 8.0)],
        );
        crit.set_weight(1.0, 0.0, 1.0, 1.0, 1.0);
        let mut g = Vector::new(0, 2);
        crit.gradient(1, 1, &mut g);
        assert!((g.get(0) - 2.0).abs() < 1e-12, "g0 = {}", g.get(0));
        assert!((g.get(1) - 6.0).abs() < 1e-12, "g1 = {}", g.get(1));
    }

    /// Hessian (quality + least-square parts) hand-derived.  Same fixture as
    /// `linear_criteria_gradient_hand` with knots [0, 2], pass points at
    /// u = 0, 2, weights: set_estimation(1, 1, 1), set_weight(QuadW = 1,
    /// QualW = 2, 1/1/1) -> myPercent = 1/3 each, so the quality part is
    /// (2/3) * sum_icrit Hessian_icrit.
    ///
    /// Hand-derived elementary criteria on [0,2] (coeff_half = 1, cteh3 = 2,
    /// R = int B_i' B_j' for B = [H00, H01, J0], J0 = TN0*(1-s^2),
    /// TN0 = sqrt(15/16)):
    ///   H00' = -1/2, H01' = +1/2, J0' = -2s*TN0
    ///   R00 = 1/2, R01 = -1/2, R11 = 1/2, R02 = R12 = 0 (odd integrand),
    ///   R22 = 4*TN0^2*int s^2 = 4*(15/16)*(2/3) = 5/2.
    ///   Tension Hessian = (2/h)   * R   = [[1,-1,0],[-1,1,0],[0,0,5]]
    ///   Flexion Hessian = (2/h^3) * R'' = [[0,0,0],[0,0,0],[0,0,2*(15/2)]]
    ///     (H00'' = H01'' = 0, J0'' = -2*TN0, R22'' = 4*TN0^2*int 1 = 15/2)
    ///   Jerk Hessian    = 0              (third derivative of a degree-2
    ///     basis vanishes).
    /// Least-square part: poid = 2 per point, BV1 = [1,0,0], BV2 = [0,1,0]
    ///   -> AuxH = diag(2, 2, 0), symmetrized into H with QuadW = 1.
    /// Total:
    ///   H(0,0) = 2 + (2/3)*1            = 8/3
    ///   H(0,1) =     (2/3)*(-1)         = -2/3
    ///   H(1,1) = 2 + (2/3)*1            = 8/3
    ///   H(2,2) =     (2/3)*(5 + 15 + 0) = 40/3
    #[test]
    fn linear_criteria_hessian_hand() {
        let mut crit = make_criteria(
            0.0,
            2.0,
            &[0.0, 2.0],
            &[DVec3::new(1.0, 7.0, 8.0), DVec3::new(3.0, 7.0, 8.0)],
        );
        crit.set_estimation(1.0, 1.0, 1.0);
        crit.set_weight(1.0, 2.0, 1.0, 1.0, 1.0);

        let mut h = Matrix::new(0, 2, 0, 2);
        crit.hessian(1, 1, 1, &mut h);

        assert!((h.get(0, 0) - 8.0 / 3.0).abs() < 1e-9, "h00 = {}", h.get(0, 0));
        assert!((h.get(0, 1) - (-2.0 / 3.0)).abs() < 1e-9, "h01 = {}", h.get(0, 1));
        assert!((h.get(1, 0) - (-2.0 / 3.0)).abs() < 1e-9, "h10 = {}", h.get(1, 0));
        assert!((h.get(1, 1) - 8.0 / 3.0).abs() < 1e-9, "h11 = {}", h.get(1, 1));
        assert!((h.get(2, 2) - 40.0 / 3.0).abs() < 1e-9, "h22 = {}", h.get(2, 2));
        assert!(h.get(0, 2).abs() < 1e-9, "h02 = {}", h.get(0, 2));
        assert!(h.get(1, 2).abs() < 1e-9, "h12 = {}", h.get(1, 2));
    }

    /// QualityValues on a linear element, hand-derived.  Element on [0, 2]
    /// with coefficients (c0, c1) = (0, 1): the curve is the straight line
    /// P(s) = H01(s) (Jacobi coefficient 0).  Tension (hand-derived above):
    /// J1 = cteh3 * (0.5*R11*c1^2) = 2 * 0.5 * 1/2 = 1/2.  Flexion / Jerk:
    /// every R''/R''' entry reached by a nonzero coefficient vanishes
    /// (H01'' = H01''' = 0), so J2 = J3 = 0 exactly.  ICDANA state machine
    /// (cxx L363-439) with myEstimation = [0,0,0]: only step (3) fires for
    /// J1 (0.5 > 0*2 -> est0 = 0.05, 0.5 > 0.5*10 is false -> ICDANA = 1).
    #[test]
    fn linear_criteria_quality_values_hand() {
        let mut crit = make_criteria(0.0, 2.0, &[0.0], &[DVec3::ZERO, DVec3::ZERO]);

        // Fill the single element with the linear coefficients (0, 1).  The
        // table is sized by the curve degree (Degree(1) = WorkDegree = 2),
        // matching OCCT `Coeff = new HArray2<double>(0, deg, 1, NbDim)`.
        let curve = crit.curve().unwrap();
        let mut coeffs = Matrix::new(0, 2, 1, 1);
        coeffs.set(0, 1, 0.0);
        coeffs.set(1, 1, 1.0);
        curve.borrow_mut().set_element(1, &coeffs);

        let mut j1 = 0.0;
        let mut j2 = 0.0;
        let mut j3 = 0.0;
        let icdana = crit.quality_values(0.01, 0.001, 0.001, &mut j1, &mut j2, &mut j3);

        assert!((j1 - 0.5).abs() < 1e-9, "j1 = {}", j1);
        assert!(j2.abs() < 1e-9, "j2 = {}", j2);
        assert!(j3.abs() < 1e-9, "j3 = {}", j3);
        assert_eq!(icdana, 1);

        // Step (3) also updated myEstimation[0] by J1 * 0.1 = 0.05.
        let mut e1 = 0.0;
        let mut e2 = 0.0;
        let mut e3 = 0.0;
        crit.get_estimation(&mut e1, &mut e2, &mut e3);
        assert!((e1 - 0.05).abs() < 1e-12, "e1 = {}", e1);
        assert!(e2.abs() < 1e-12 && e3.abs() < 1e-12);
    }

    /// InputVector round-trip through the hand-derived AssemblyTable of
    /// `linear_criteria_assembly_table_hand`: X = [10, 20, 30, 40, 50] with
    /// the global map el1 = (1,3,2), el2 = (3,5,4) must produce element
    /// coefficients el1 = (10, 30, 20), el2 = (30, 50, 40)
    /// (CoeffEl(i, dim) = X(i0 + GlobIndex(i)), cxx L750-762).
    #[test]
    fn linear_criteria_input_vector_hand() {
        let crit = make_criteria(0.0, 1.0, &[0.0], &[DVec3::ZERO, DVec3::ZERO]);
        let base = rcad_kernel::math::plib::hermit_jacobi::HermitJacobi::new(2, GeomAbsShape::C0);
        let mut curve2 = Curve::new(1, 2, &base, 0.0);
        curve2.knots()[0] = 0.0;
        curve2.knots()[1] = 1.0;
        curve2.knots()[2] = 2.0;
        let curve2: CurveHandle = Rc::new(RefCell::new(curve2));
        let mut crit = crit;
        crit.set_curve(&curve2);

        let table = crit.assembly_table();
        let mut x = Vector::new(1, 5);
        for (i, v) in [10.0, 20.0, 30.0, 40.0, 50.0].iter().enumerate() {
            x.set(1 + i as i32, *v);
        }
        crit.input_vector(&x, &table);

        let mut cm = Matrix::new(0, 2, 1, 1);
        let mut c2 = curve2.borrow_mut();
        c2.get_element(1, &mut cm);
        assert_eq!(cm.get(0, 1), 10.0);
        assert_eq!(cm.get(1, 1), 30.0);
        assert_eq!(cm.get(2, 1), 20.0);
        c2.get_element(2, &mut cm);
        assert_eq!(cm.get(0, 1), 30.0);
        assert_eq!(cm.get(1, 1), 50.0);
        assert_eq!(cm.get(2, 1), 40.0);
    }

    /// EstLength returns the mutable `double&` (cxx L174-177): writing
    /// through it must be observable via GetEstimation-free direct reads
    /// and via the SmoothCriterion trait object dispatch.
    #[test]
    fn linear_criteria_est_length_and_trait_dispatch() {
        let crit = make_criteria(0.0, 2.0, &[0.0], &[DVec3::ZERO, DVec3::ZERO]);
        let mut crit: Box<dyn SmoothCriterion> = Box::new(crit);

        *crit.est_length() = 3.5;
        assert!((*crit.est_length() - 3.5).abs() < 1e-15);

        // GetCurve / Curve() convenience on the trait object.
        let c = crit.curve().expect("curve handle");
        assert_eq!(c.borrow().nb_elements(), 1);
        assert_eq!(c.borrow().dimension(), 1);

        // SetEstimation / GetEstimation round trip.
        crit.set_estimation(1.0, 2.0, 3.0);
        let (mut e1, mut e2, mut e3) = (0.0, 0.0, 0.0);
        crit.get_estimation(&mut e1, &mut e2, &mut e3);
        assert_eq!((e1, e2, e3), (1.0, 2.0, 3.0));

        // SetWeight / GetWeight round trip and percent normalization.
        let (mut qw, mut jw) = (0.0, 0.0);
        crit.set_weight(4.0, 6.0, 1.0, 1.0, 2.0);
        crit.get_weight(&mut qw, &mut jw);
        assert_eq!((qw, jw), (4.0, 6.0));
        crit.set_weight_array(&RealArray1::new(1, 2));
    }

    /// Unused-import guard: MultiPointConstraint is re-exported here only to
    /// keep the fixture imports honest (the MultiLine constructors already
    /// build the constraints internally).
    #[allow(dead_code)]
    fn _mpc_type_witness() -> MultiPointConstraint {
        MultiPointConstraint::new()
    }
}
