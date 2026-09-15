//! OCCT AppDef_Variational - private engine methods, part 3
//! (AppDef_Variational.cxx L2129-3403): InitSmoothCriterion,
//! InitParameters, InitCriterionEstimations, EstTangent, EstSecnd,
//! InitCutting, Adjusting, AssemblingConstraints, InitTthetaF.
//! Also hosts the discriminating unit tests.
//!
//! Split off from `app_def_variational.rs` via `#[path]`; the parent's
//! private fields and helpers are visible here (child module).

use std::cell::RefCell;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::math::fem_tool::{Assembly, Curve};
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::math::plib::hermit_jacobi::HermitJacobi;
use rcad_kernel::core::precision::CONFUSION;

use crate::geomalgo::app_def_linear_criteria::LinearCriteria;
use super::variational_b::near_index;
use super::{
    not_parallel, v_add, v_div_assign, v_mul, v_norm, v_norm2, v_sub, AppDefVariational,
    CriterionHandle, CurveHandle, RealArray1,
};

/// OCCT math_VectorBase::Init (math_VectorBase.lxx) - fills with a value.
fn v_init(v: &mut Vector, value: f64) {
    for i in v.lower()..=v.upper() {
        v.set(i, value);
    }
}

impl AppDefVariational {
    /// OCCT AppDef_Variational::InitSmoothCriterion (cxx L2129-2204).
    pub(super) fn init_smooth_criterion(&mut self) {

        let eps2 = 1.0e-6;
        let eps3 = 1.0e-9;
        //  const double J1 = .01, J2 = .001, J3 = .001;

        let mut length = 0.0;

        self.init_parameters(&mut length);

        self.my_smooth_criterion
            .borrow_mut()
            .set_parameters(&self.my_parameters);

        let (mut e1, mut e2, mut e3) = (0.0f64, 0.0f64, 0.0f64);

        self.init_criterion_estimations(length, &mut e1, &mut e2, &mut e3);
        /*
        J1 = 1.e-8; J2 = J3 = (E1 + 1.e-8) * 1.e-6;

        if(E1 < J1) E1 = J1;
        if(E2 < J2) E2 = J2;
        if(E3 < J3) E3 = J3;
        */
        *self.my_smooth_criterion.borrow_mut().est_length() = length;
        self.my_smooth_criterion.borrow_mut().set_estimation(e1, e2, e3);

        let mut w_quadratic;
        let w_quality;

        if !self.my_with_min_max && self.my_tolerance != 0.0 {
            w_quality = self.my_tolerance;
        } else if self.my_tolerance == 0.0 {
            w_quality = 1.0;
        } else {
            w_quality = self.my_tolerance.max(eps2 * length);
        }

        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        w_quadratic = ((self.my_nb_points - nb_constr) as f64).sqrt() * w_quality;
        if w_quadratic > eps3 {
            w_quadratic = 1.0 / w_quadratic;
        }

        if w_quadratic == 0.0 {
            w_quadratic = e1.sqrt().max(1.0);
        }

        self.my_smooth_criterion.borrow_mut().set_weight(
            w_quadratic,
            w_quality,
            self.my_percent[0],
            self.my_percent[1],
            self.my_percent[2],
        );

        let the_base = HermitJacobi::new(self.my_max_degree, self.my_continuity);
        let mut the_curve: Option<CurveHandle> = None;
        let nb_elem;
        let curv_tol = eps2 * length / self.my_nb_points as f64;

        // Subdivision of the interval according to constraints
        if self.my_with_cutting && nb_constr != 0 {

            self.init_cutting(&the_base, curv_tol, &mut the_curve);
        } else {

            nb_elem = 1;
            the_curve = Some(Rc::new(RefCell::new(Curve::new(
                self.my_dimension,
                nb_elem,
                &the_base,
                curv_tol,
            ))));
            // TheCurve->Knots().SetValue(TheCurve->Knots().Lower(), ...);
            // TheCurve->Knots().SetValue(TheCurve->Knots().Upper(), ...);
            let nb_knots = the_curve.as_ref().unwrap().borrow().nb_elements() + 1;
            the_curve.as_ref().unwrap().borrow_mut().knots()[0] =
                self.my_parameters.value(self.my_first_point);
            the_curve.as_ref().unwrap().borrow_mut().knots()[(nb_knots - 1) as usize] =
                self.my_parameters.value(self.my_last_point);
        }

        self.my_smooth_criterion
            .borrow_mut()
            .set_curve(the_curve.as_ref().unwrap());
    }

    /// OCCT AppDef_Variational::InitParameters (cxx L2209-2251).
    fn init_parameters(&mut self, length: &mut f64) {

        // constexpr double Eps1 = Precision::Confusion() * .01;
        let eps1 = CONFUSION * 0.01;

        let mut aux;
        let mut dist;
        let mut i1 = 0;
        let mut i0;

        *length = 0.0;
        self.my_parameters.set_value(self.my_first_point, *length);

        for ipoint in (self.my_first_point + 1)..=self.my_last_point {
            i0 = i1;
            i1 += self.my_dimension;
            dist = 0.0;
            for i in 1..=self.my_dimension {
                aux = self.my_tab_points.value(i1 + i) - self.my_tab_points.value(i0 + i);
                dist += aux * aux;
            }
            *length += dist.sqrt();
            self.my_parameters.set_value(ipoint, *length);
        }

        if *length <= eps1 {
            panic!("Standard_ConstructionError: AppDef_Variational::InitParameters");
        }

        for ipoint in (self.my_first_point + 1)..=(self.my_last_point - 1) {
            let v = self.my_parameters.value(ipoint) / *length;
            self.my_parameters.set_value(ipoint, v);
        }

        self.my_parameters.set_value(self.my_last_point, 1.0);

        // With few points there is likely underestimation ...
        if self.my_nb_points < 10 {
            *length *= 1.0 + 0.1 / (self.my_nb_points - 1) as f64;
        }
    }

    /// OCCT AppDef_Variational::InitCriterionEstimations (cxx L2256-2395).
    fn init_criterion_estimations(&self, length: f64, e1: &mut f64, e2: &mut f64, e3: &mut f64) {
        *e1 = length * length;

        // constexpr double Eps1 = Precision::Confusion() * .01;
        let eps1 = CONFUSION * 0.01;

        let mut vtang1 = Vector::new(1, self.my_dimension);
        let mut vtang2 = Vector::new(1, self.my_dimension);
        let mut vtang3 = Vector::new(1, self.my_dimension);
        let mut vscnd1 = Vector::new(1, self.my_dimension);
        let mut vscnd2 = Vector::new(1, self.my_dimension);
        let mut vscnd3 = Vector::new(1, self.my_dimension);

        // ========== Treatment of first point =================

        let mut ipnt = self.my_first_point;

        self.est_tangent(ipnt, &mut vtang1);
        ipnt += 1;
        self.est_tangent(ipnt, &mut vtang2);
        ipnt += 1;
        self.est_tangent(ipnt, &mut vtang3);

        let mut ipnt = self.my_first_point;
        self.est_secnd(ipnt, &vtang1, &vtang2, length, &mut vscnd1);
        ipnt += 1;
        self.est_secnd(ipnt, &vtang1, &vtang3, length, &mut vscnd2);

        //  Modified by skv - Fri Apr  8 14:58:12 2005 OCC8559 Begin
        //   double Delta = .5 * (myParameters->Value(ipnt) - myParameters->Value(--ipnt));
        let an_ind = ipnt;
        ipnt -= 1; // --ipnt
        let mut delta = 0.5 * (self.my_parameters.value(an_ind) - self.my_parameters.value(ipnt));
        //  Modified by skv - Fri Apr  8 14:58:12 2005 OCC8559 End

        if delta <= eps1 {
            delta = 1.0;
        }

        *e2 = v_norm2(&vscnd1) * delta;

        *e3 = if delta > eps1 {
            v_norm2(&v_sub(&vscnd2, &vscnd1)) / (4.0 * delta)
        } else {
            0.0
        };
        // ========== Treatment of internal points =================

        let mut curr_point = 2;

        let mut ipnt = self.my_first_point + 1;
        while ipnt < self.my_last_point {

            delta = 0.5 * (self.my_parameters.value(ipnt + 1) - self.my_parameters.value(ipnt - 1));

            if curr_point == 1 {
                if ipnt + 1 != self.my_last_point {
                    self.est_tangent(ipnt + 2, &mut vtang3);
                    self.est_secnd(ipnt + 1, &vtang1, &vtang3, length, &mut vscnd2);
                } else {
                    self.est_secnd(ipnt + 1, &vtang1, &vtang2, length, &mut vscnd2);
                }

                *e2 += v_norm2(&vscnd1) * delta;
                *e3 += if delta > eps1 {
                    v_norm2(&v_sub(&vscnd2, &vscnd3)) / (4.0 * delta)
                } else {
                    0.0
                };
            } else if curr_point == 2 {
                if ipnt + 1 != self.my_last_point {
                    self.est_tangent(ipnt + 2, &mut vtang1);
                    self.est_secnd(ipnt + 1, &vtang2, &vtang1, length, &mut vscnd3);
                } else {
                    self.est_secnd(ipnt + 1, &vtang2, &vtang3, length, &mut vscnd3);
                }

                *e2 += v_norm2(&vscnd2) * delta;
                *e3 += if delta > eps1 {
                    v_norm2(&v_sub(&vscnd3, &vscnd1)) / (4.0 * delta)
                } else {
                    0.0
                };
            } else {
                if ipnt + 1 != self.my_last_point {
                    self.est_tangent(ipnt + 2, &mut vtang2);
                    self.est_secnd(ipnt + 1, &vtang3, &vtang2, length, &mut vscnd1);
                } else {
                    self.est_secnd(ipnt + 1, &vtang3, &vtang1, length, &mut vscnd1);
                }

                *e2 += v_norm2(&vscnd3) * delta;
                *e3 += if delta > eps1 {
                    v_norm2(&v_sub(&vscnd1, &vscnd2)) / (4.0 * delta)
                } else {
                    0.0
                };
            }

            curr_point += 1;
            if curr_point == 4 {
                curr_point = 1;
            }
            ipnt += 1;
        }

        // ========== Treatment of last point =================

        let mut delta =
            0.5 * (self.my_parameters.value(self.my_last_point) - self.my_parameters.value(self.my_last_point - 1));
        if delta <= eps1 {
            delta = 1.0;
        }

        let aux;

        if curr_point == 1 {

            *e2 += v_norm2(&vscnd1) * delta;
            aux = v_norm2(&v_sub(&vscnd1, &vscnd3));
            *e3 += if delta > eps1 { aux / (4.0 * delta) } else { aux };
        } else if curr_point == 2 {

            *e2 += v_norm2(&vscnd2) * delta;
            aux = v_norm2(&v_sub(&vscnd2, &vscnd1));
            *e3 += if delta > eps1 { aux / (4.0 * delta) } else { aux };
        } else {

            *e2 += v_norm2(&vscnd3) * delta;
            aux = v_norm2(&v_sub(&vscnd3, &vscnd2));
            *e3 += if delta > eps1 { aux / (4.0 * delta) } else { aux };
        }

        let aux = length * length;

        *e2 *= aux;
        *e3 *= aux;
    }

    /// OCCT AppDef_Variational::EstTangent (cxx L2401-2580).
    fn est_tangent(&self, ipnt: i32, vtang: &mut Vector) {

        // constexpr double Eps1 = Precision::Confusion() * .01;
        let eps1 = CONFUSION * 0.01;
        let eps_norm = 1.0e-9;

        let mut wpnt = 1.0;

        if ipnt == self.my_first_point {
            // Estimation at first point
            if self.my_nb_points < 3 {
                wpnt = 0.0;
            } else {

                let adr1 = 1;
                let adr2 = adr1 + self.my_dimension;
                let adr3 = adr2 + self.my_dimension;

                let pnt1 = self.tab_points_vector(adr1);
                let pnt2 = self.tab_points_vector(adr2);
                let pnt3 = self.tab_points_vector(adr3);

                // Parabolic interpolation
                // if we have parabolic interpolation: F(t) = A0 + A1*t + A2*t*t,
                // first derivative for t=0 is A1 = ((d2-1)*P1 + P2 - d2*P3)/(d*(1-d))
                //       d= |P2-P1|/(|P2-P1|+|P3-P2|), d2 = d*d
                let v1 = v_norm(&v_sub(&pnt2, &pnt1));
                let mut v2 = 0.0;
                if v1 > eps1 {
                    v2 = v_norm(&v_sub(&pnt3, &pnt2));
                }
                if v2 > eps1 {
                    let mut d = v1 / (v1 + v2);
                    let d1;
                    d1 = 1.0 / (d * (1.0 - d));
                    d *= d;
                    // VTang = ((d - 1.) * Pnt1 + Pnt2 - d * Pnt3) * d1;
                    *vtang = v_mul(
                        &v_sub(&v_add(&v_mul(&pnt1, d - 1.0), &pnt2), &v_mul(&pnt3, d)),
                        d1,
                    );
                } else {
                    // Simple 2-point estimation

                    *vtang = v_sub(&pnt2, &pnt1);
                }
            }
        } else if ipnt == self.my_last_point {
            // Estimation at last point
            if self.my_nb_points < 3 {
                wpnt = 0.0;
            } else {

                let adr1 = (self.my_last_point - 3) * self.my_dimension + 1;
                let adr2 = adr1 + self.my_dimension;
                let adr3 = adr2 + self.my_dimension;

                let pnt1 = self.tab_points_vector(adr1);
                let pnt2 = self.tab_points_vector(adr2);
                let pnt3 = self.tab_points_vector(adr3);

                // Parabolic interpolation
                // if we have parabolic interpolation: F(t) = A0 + A1*t + A2*t*t,
                // first derivative for t=1 is 2*A2 + A1 = ((d2+1)*P1 - P2 - d2*P3)/(d*(1-d))
                //       d= |P2-P1|/(|P2-P1|+|P3-P2|), d2 = d*(d-2)
                let v1 = v_norm(&v_sub(&pnt2, &pnt1));
                let mut v2 = 0.0;
                if v1 > eps1 {
                    v2 = v_norm(&v_sub(&pnt3, &pnt2));
                }
                if v2 > eps1 {
                    let mut d = v1 / (v1 + v2);
                    let d1;
                    d1 = 1.0 / (d * (1.0 - d));
                    d *= d - 2.0;
                    // VTang = ((d + 1.) * Pnt1 - Pnt2 - d * Pnt3) * d1;
                    *vtang = v_mul(
                        &v_sub(&v_sub(&v_mul(&pnt1, d + 1.0), &pnt2), &v_mul(&pnt3, d)),
                        d1,
                    );
                } else {
                    // Simple 2-point estimation

                    *vtang = v_sub(&pnt3, &pnt2);
                }
            }
        } else {

            let adr1 = (ipnt - self.my_first_point - 1) * self.my_dimension + 1;
            let adr2 = adr1 + 2 * self.my_dimension;

            let pnt1 = self.tab_points_vector(adr1);
            let pnt2 = self.tab_points_vector(adr2);

            *vtang = v_sub(&pnt2, &pnt1);
        }

        let mut vnorm = v_norm(vtang);

        if vnorm <= eps_norm {
            v_init(vtang, 0.0);
        } else {
            v_div_assign(vtang, vnorm);
        }

        // Estimation with constraints

        let mut wcnt = 0.0;
        let mut id_cnt = 1;

        // Warning!! Here it is suppoused that all points are in range [myFirstPoint, myLastPoint]

        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        let mut vcnt = Vector::new_init(1, self.my_dimension, 0.0);

        if nb_constr > 0 {

            while self.my_typ_constraints.value(2 * id_cnt - 1) < ipnt && id_cnt <= nb_constr {
                id_cnt += 1;
            }
            if (self.my_typ_constraints.value(2 * id_cnt - 1) == ipnt)
                && (self.my_typ_constraints.value(2 * id_cnt) >= 1)
            {
                wcnt = 1.0;
                let mut i0 = 2 * self.my_dimension * (id_cnt - 1);
                let mut k = 0;
                for _i in 1..=self.my_nb_p3d {
                    for _j in 1..=3 {
                        k += 1;
                        i0 += 1;
                        vcnt.set(k, self.my_tab_constraints.value(i0));
                    }
                    i0 += 3;
                }
                for _i in 1..=self.my_nb_p2d {
                    for _j in 1..=2 {
                        k += 1;
                        i0 += 1;
                        vcnt.set(k, self.my_tab_constraints.value(i0));
                    }
                    i0 += 2;
                }
            }
        }

        // Averaging of estimation

        let mut denom = wpnt + wcnt;
        if denom == 0.0 {
            denom = 1.0;
        } else {
            denom = 1.0 / denom;
        }

        // VTang = (Wpnt * VTang + Wcnt * VCnt) * Denom;
        *vtang = v_mul(&v_add(&v_mul(vtang, wpnt), &v_mul(&vcnt, wcnt)), denom);

        vnorm = v_norm(vtang);

        if vnorm <= eps_norm {
            v_init(vtang, 0.0);
        } else {
            v_div_assign(vtang, vnorm);
        }
    }

    /// OCCT AppDef_Variational::EstSecnd (cxx L2585-2679).
    fn est_secnd(
        &self,
        ipnt: i32,
        vtang1: &Vector,
        vtang2: &Vector,
        length: f64,
        vscnd: &mut Vector,
    ) {
        // const double Eps = 1.e-9;
        let eps = 1.0e-9;

        let wpnt = 1.0;

        let mut aux;

        if ipnt == self.my_first_point {
            aux = self.my_parameters.value(ipnt + 1) - self.my_parameters.value(ipnt);
        } else if ipnt == self.my_last_point {
            aux = self.my_parameters.value(ipnt) - self.my_parameters.value(ipnt - 1);
        } else {
            aux = self.my_parameters.value(ipnt + 1) - self.my_parameters.value(ipnt - 1);
        }

        if aux <= eps {
            aux = 1.0;
        } else {
            aux = 1.0 / aux;
        }

        // VScnd = (VTang2 - VTang1) * aux;
        *vscnd = v_mul(&v_sub(vtang2, vtang1), aux);

        // Estimation with constraints

        let mut wcnt = 0.0;
        let mut id_cnt = 1;

        // Warning!! Here it is suppoused that all points are in range [myFirstPoint, myLastPoint]

        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        let mut vcnt = Vector::new_init(1, self.my_dimension, 0.0);

        if nb_constr > 0 {

            while self.my_typ_constraints.value(2 * id_cnt - 1) < ipnt && id_cnt <= nb_constr {
                id_cnt += 1;
            }

            if (self.my_typ_constraints.value(2 * id_cnt - 1) == ipnt)
                && (self.my_typ_constraints.value(2 * id_cnt) >= 2)
            {
                wcnt = 1.0;
                let mut i0 = 2 * self.my_dimension * (id_cnt - 1) + 3;
                let mut k = 0;
                for _i in 1..=self.my_nb_p3d {
                    for _j in 1..=3 {
                        k += 1;
                        i0 += 1;
                        vcnt.set(k, self.my_tab_constraints.value(i0));
                    }
                    i0 += 3;
                }
                i0 -= 1;
                for _i in 1..=self.my_nb_p2d {
                    for _j in 1..=2 {
                        k += 1;
                        i0 += 1;
                        vcnt.set(k, self.my_tab_constraints.value(i0));
                    }
                    i0 += 2;
                }
            }
        }

        // Averaging of estimation

        let mut denom = wpnt + wcnt;
        if denom == 0.0 {
            denom = 1.0;
        } else {
            denom = 1.0 / denom;
        }

        // VScnd = (Wpnt * VScnd + (Wcnt * Length) * VCnt) * Denom;
        *vscnd = v_mul(
            &v_add(&v_mul(vscnd, wpnt), &v_mul(&vcnt, wcnt * length)),
            denom,
        );
    }

    /// OCCT AppDef_Variational::InitCutting (cxx L2684-2813).
    fn init_cutting(&self, a_base: &HermitJacobi, curv_tol: f64, a_curve: &mut Option<CurveHandle>) {

        // Definition of number of elements
        let mut orcmx = -1;
        let mut ncont = 0;
        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;

        for i in 1..=nb_constr {
            let kk = self.my_typ_constraints.value(2 * i).abs() + 1;
            orcmx = orcmx.max(kk);
            ncont += kk;
        }

        if orcmx > self.my_max_degree - self.my_niv_cont {
            panic!("Standard_ConstructionError: AppDef_Variational::InitCutting");
        }

        let mut nlibre = (self.my_max_degree - self.my_niv_cont - (self.my_max_degree + 1) / 4)
            .max(self.my_niv_cont + 1);

        let mut nb_elem = if ncont % nlibre == 0 {
            ncont / nlibre
        } else {
            ncont / nlibre + 1
        };

        while (nb_elem > self.my_max_segment) && (nlibre < self.my_max_degree - self.my_niv_cont) {

            nlibre += 1;
            nb_elem = if ncont % nlibre == 0 {
                ncont / nlibre
            } else {
                ncont / nlibre + 1
            };
        }

        if nb_elem > self.my_max_segment {
            panic!("Standard_ConstructionError: AppDef_Variational::InitCutting");
        }

        let curve = Rc::new(RefCell::new(Curve::new(self.my_dimension, nb_elem, a_base, curv_tol)));

        let mut n_cnt = (ncont - 1) / nb_elem + 1;
        let mut nplus = nb_elem - (n_cnt * nb_elem - ncont);

        // NCollection_Array1<double>& Knot = aCurve->Knots();
        let mut i_deb = 0;
        let mut i_fin = nb_constr + 1;
        let mut n_deb = 0;
        let mut n_fin = 0;
        let mut ind_el = 1; // Knot.Lower()
        let i_upper = nb_elem + 1; // Knot.Upper()
        let mut nb_el = 0;

        curve.borrow_mut().knots()[(ind_el - 1) as usize] =
            self.my_parameters.value(self.my_first_point);
        curve.borrow_mut().knots()[(i_upper - 1) as usize] =
            self.my_parameters.value(self.my_last_point);

        while nb_elem - nb_el > 1 {

            ind_el += 1;
            nb_el += 1;
            if nplus == 0 {
                n_cnt -= 1;
            }

            while n_deb < n_cnt && i_deb < i_fin {
                i_deb += 1;
                n_deb += self.my_typ_constraints.value(2 * i_deb).abs() + 1;
            }

            if n_deb == n_cnt {
                n_deb = 0;
                if nplus == 1
                    && self.my_parameters.value(self.my_typ_constraints.value(2 * i_deb - 1))
                        > curve.borrow_mut().knots()[(ind_el - 1 - 1) as usize]
                {

                    curve.borrow_mut().knots()[(ind_el - 1) as usize] =
                        self.my_parameters.value(self.my_typ_constraints.value(2 * i_deb - 1));
                } else {
                    curve.borrow_mut().knots()[(ind_el - 1) as usize] = (self
                        .my_parameters
                        .value(self.my_typ_constraints.value(2 * i_deb - 1))
                        + self
                            .my_parameters
                            .value(self.my_typ_constraints.value(2 * i_deb + 1)))
                        / 2.0;
                }
            } else {
                n_deb -= n_cnt;
                curve.borrow_mut().knots()[(ind_el - 1) as usize] =
                    self.my_parameters.value(self.my_typ_constraints.value(2 * i_deb - 1));
            }

            nplus -= 1;
            if nplus == 0 {
                n_cnt -= 1;
            }

            if nb_elem - nb_el == 1 {
                break;
            }

            nb_el += 1;

            while n_fin < n_cnt && i_deb < i_fin {
                i_fin -= 1;
                n_fin += self.my_typ_constraints.value(2 * i_fin).abs() + 1;
            }

            if n_fin == n_cnt {
                n_fin = 0;
                let v = (self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 1))
                    + self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 3)))
                    / 2.0;
                curve.borrow_mut().knots()[((i_upper + 1 - ind_el) - 1) as usize] = v;
            } else {
                n_fin -= n_cnt;
                if self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 1))
                    < curve.borrow_mut().knots()[((i_upper - ind_el + 1) - 1) as usize]
                {
                    curve.borrow_mut().knots()[((i_upper + 1 - ind_el) - 1) as usize] =
                        self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 1));
                } else {
                    let v = (self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 1))
                        + self.my_parameters.value(self.my_typ_constraints.value(2 * i_fin - 3)))
                        / 2.0;
                    curve.borrow_mut().knots()[((i_upper + 1 - ind_el) - 1) as usize] = v;
                }
            }
        }

        *a_curve = Some(curve);
    }

    /// OCCT AppDef_Variational::Adjusting (cxx L2816-2974).
    pub(super) fn adjusting(
        &mut self,
        j: &mut CriterionHandle,
        w_quadratic: &mut f64,
        w_quality: &mut f64,
        the_curve: &mut CurveHandle,
        ecarts: &mut RealArray1,
    ) {

        //  std::cout << "=========== Adjusting =============" << std::endl;

        /* Initialized data */

        let mxiter = 2;
        let eps1 = 1.0e-6;
        let nbr_pnt = self.my_last_point - self.my_first_point + 1;
        let nbr_constraint = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        let curv_tol = eps1 * *j.borrow_mut().est_length() / nbr_pnt as f64;

        /* Local variables */
        let mut iter;
        let mut ecart;
        let mut erold;
        let mut emold;
        let mut tpara;
        let mut vocri = [0.0f64; 4];
        let mut j1cibl;
        let vtest;
        let mut vseuil;
        let mut numint;
        let mut flag = 0;
        let mut tbpoid = RealArray1::new(self.my_first_point, self.my_last_point);
        let mut loptim;
        let mut lrejet;
        let mut j_new: CriterionHandle;
        let mut c_new: CurveHandle;
        let (mut ge1, mut ge2, mut ge3) = (0.0f64, 0.0f64, 0.0f64);

        /* (0.b) Initialisations */

        loptim = true;
        iter = 0;
        tbpoid.init(1.0);

        /* ============   loop on the smoothing engine  ============== */

        vtest = *w_quality * 0.9;
        j1cibl = (self.my_criterium[0] / (nbr_pnt - nbr_constraint) as f64).sqrt();

        while loptim {

            iter += 1;

            /*     (1) Save the previous state */

            vocri[0] = self.my_criterium[0];
            vocri[1] = self.my_criterium[1];
            vocri[2] = self.my_criterium[2];
            vocri[3] = self.my_criterium[3];
            erold = self.my_max_error;
            emold = self.my_average_error;

            /*     (2) Increase the least squares weight */

            if j1cibl > vtest {
                *w_quadratic = j1cibl / vtest * *w_quadratic;
            }

            /*     (3) Increase the weight associated with problematic points */

            vseuil = *w_quality * 0.88;

            for ipnt in self.my_first_point..=self.my_last_point {
                if ecarts.value(ipnt) > vtest {
                    ecart = (ecarts.value(ipnt) - vseuil) / *w_quality;
                    let v = (ecart * 3.0 + 1.0) * tbpoid.value(ipnt);
                    tbpoid.set_value(ipnt, v);
                }
            }

            /*     (4) Forced splitting */

            if the_curve.borrow().nb_elements() < self.my_max_segment && self.my_with_cutting {

                // numint = NearIndex(myParameters->Value(myMaxErrorIndex),
                //                    TheCurve->Knots(), 0, flag); (the knots are
                // only read here - snapshot with identical semantics)
                let knots_snapshot: Vec<f64> = the_curve.borrow_mut().knots().to_vec();
                numint = near_index(
                    self.my_parameters.value(self.my_max_error_index),
                    &knots_snapshot,
                    0.0,
                    &mut flag,
                );

                tpara = (knots_snapshot[(numint - 1) as usize]
                    + knots_snapshot[numint as usize]
                    + self.my_parameters.value(self.my_max_error_index) * 2.0)
                    / 4.0;

                c_new = Rc::new(RefCell::new(Curve::new(
                    self.my_dimension,
                    the_curve.borrow().nb_elements() + 1,
                    the_curve.borrow().base(),
                    curv_tol,
                )));

                let nb_knots = the_curve.borrow().nb_elements() + 1;
                for i in 1..=numint {
                    let v = knots_snapshot[(i - 1) as usize];
                    c_new.borrow_mut().knots()[(i - 1) as usize] = v;
                }
                for i in (numint + 1)..=nb_knots {
                    let v = knots_snapshot[(i - 1) as usize];
                    c_new.borrow_mut().knots()[i as usize] = v;
                }

                c_new.borrow_mut().knots()[numint as usize] = tpara;
            } else {

                c_new = Rc::new(RefCell::new(Curve::new(
                    self.my_dimension,
                    the_curve.borrow().nb_elements(),
                    the_curve.borrow().base(),
                    curv_tol,
                )));

                // CNew->Knots() = TheCurve->Knots();
                {
                    let src = the_curve.borrow_mut().knots().to_vec();
                    c_new.borrow_mut().knots().copy_from_slice(&src);
                }
            }

            j_new = Rc::new(RefCell::new(LinearCriteria::new(
                &self.my_ssp,
                self.my_first_point,
                self.my_last_point,
            )));

            *j_new.borrow_mut().est_length() = *j.borrow_mut().est_length();

            j.borrow().get_estimation(&mut ge1, &mut ge2, &mut ge3);

            j_new.borrow_mut().set_estimation(ge1, ge2, ge3);

            j_new.borrow_mut().set_parameters(&self.my_parameters);

            j_new.borrow_mut().set_weight(
                *w_quadratic,
                *w_quality,
                self.my_percent[0],
                self.my_percent[1],
                self.my_percent[2],
            );

            j_new.borrow_mut().set_weight_array(&tbpoid);

            j_new.borrow_mut().set_curve(&c_new);

            /*     (5) Relissage */

            self.the_motor(&j_new, *w_quadratic, *w_quality, &mut c_new, ecarts);

            /*     (6) Tests de rejet */

            j1cibl = (self.my_criterium[0] / (nbr_pnt - nbr_constraint) as f64).sqrt();
            vseuil = vocri[1].sqrt() + (erold - self.my_max_error) * 4.0;

            lrejet = (self.my_max_error > *w_quality && self.my_max_error > erold * 1.01)
                || (self.my_criterium[1].sqrt() > vseuil * 1.05);

            if lrejet {
                self.my_criterium[0] = vocri[0];
                self.my_criterium[1] = vocri[1];
                self.my_criterium[2] = vocri[2];
                self.my_criterium[3] = vocri[3];
                self.my_max_error = erold;
                self.my_average_error = emold;

                loptim = false;
            } else {
                *j = j_new.clone();
                *the_curve = c_new.clone();
                j.borrow_mut().set_curve(the_curve);
            }

            /*     (7) Test de convergence */

            if ((iter >= mxiter) && (self.my_max_segment == c_new.borrow().nb_elements()))
                || self.my_max_error < *w_quality
            {
                loptim = false;
            }
        }
    }

    /// OCCT AppDef_Variational::AssemblingConstraints (cxx L2993-3313).
    pub(super) fn assembling_constraints(
        &self,
        curve: &CurveHandle,
        parameters: &RealArray1,
        cblong: f64,
        a: &mut Assembly,
    ) {

        let (mx_deg, nb_elm, nb_dim) = {
            let c = curve.borrow();
            (c.base().work_degree(), c.nb_elements(), c.dimension())
        };

        // OCCT: NCollection_Array1<double> G0(0, MxDeg), G1(0, MxDeg), G2(0,
        // MxDeg); math_Vector V0((double*)&G0(0), 0, MxDeg), V1(...), V2(...);
        // (the math_Vector is an aliased view over the Array1 storage; the
        // Rust Vector is the single storage written by D0/D1/D2 and read by
        // AddConstraint.)
        let mut v0 = Vector::new(0, mx_deg);
        let mut v1 = Vector::new(0, mx_deg);
        let mut v2 = Vector::new(0, mx_deg);

        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        let p0 = parameters.lower() - self.my_first_point;
        let mut curel = 1;
        let ntheta = 6 * self.my_nb_p3d + 2 * self.my_nb_p2d;

        //  Ng3d = 3 * NbConstr + 2 * myNbTangPoints + 5 * myNbCurvPoints;
        //  Ng2d = 2 * NbConstr + 1 * myNbTangPoints + 3 * myNbCurvPoints;
        let ng3d = 3 * nb_constr + 3 * self.my_nb_tang_points + 5 * self.my_nb_curv_points;
        let ng2d = 2 * nb_constr + 2 * self.my_nb_tang_points + 3 * self.my_nb_curv_points;
        let nbeg2d = ng3d * self.my_nb_p3d;
        //  NgPC1 = NbConstr + myNbCurvPoints;
        let ngpc1 = nb_constr + self.my_nb_tang_points + self.my_nb_curv_points;
        let mut npass = 0;
        let mut ntang3d = 3 * ngpc1;
        let mut ntang2d = 2 * ngpc1;

        // OCCT: NCollection_Array1<double>& Intervals = Curve->Knots(); - a
        // read-only view in this scope; snapshot (identical semantics - the
        // knots are not modified here).
        let intervals: Vec<f64> = curve.borrow_mut().knots().to_vec();

        // OCCT: const PLib_HermitJacobi& myHermitJacobi = Curve->Base();
        // (cloned so the Ref is not held across the AddConstraint calls).
        let my_hermit_jacobi = curve.borrow().base().clone();
        let order = my_hermit_jacobi.niv_constr() + 1;

        a.nullify_constraint();

        let mut ipnt = -1;
        let mut ityp = 0;
        for i in 1..=nb_constr {

            ipnt += 2;
            ityp += 2;

            let point = self.my_typ_constraints.value(ipnt);
            let typ_of_constr = self.my_typ_constraints.value(ityp);

            let mut t = parameters.value(p0 + point);

            let mut el = curel;
            while el <= nb_elm {
                el += 1;
                if t <= intervals[(el - 1) as usize] {
                    curel = el - 1;
                    break;
                }
            }

            let ufirst = intervals[(curel - 1) as usize];
            let ulast = intervals[curel as usize];
            let coeff = (ulast - ufirst) / 2.0;
            let c0 = (ulast + ufirst) / 2.0;

            t = (t - c0) / coeff;

            if typ_of_constr == 0 {
                my_hermit_jacobi.d0(t, v0.data.v.as_mut_slice());
                for k in 1..order {
                    let mfact = coeff.powi(k);
                    let g = v0.get(k) * mfact;
                    v0.set(k, g);
                    let g = v0.get(k + order) * mfact;
                    v0.set(k + order, g);
                }
            } else if typ_of_constr == 1 {
                my_hermit_jacobi.d1(t, v0.data.v.as_mut_slice(), v1.data.v.as_mut_slice());
                for k in 1..order {
                    let mfact = coeff.powi(k);
                    let g = v0.get(k) * mfact;
                    v0.set(k, g);
                    let g = v0.get(k + order) * mfact;
                    v0.set(k + order, g);
                    let g = v1.get(k) * mfact;
                    v1.set(k, g);
                    let g = v1.get(k + order) * mfact;
                    v1.set(k + order, g);
                }
                let mfact = 1.0 / coeff;
                for k in 0..=mx_deg {
                    let g = v1.get(k) * mfact;
                    v1.set(k, g);
                }
            } else {
                my_hermit_jacobi.d2(
                    t,
                    v0.data.v.as_mut_slice(),
                    v1.data.v.as_mut_slice(),
                    v2.data.v.as_mut_slice(),
                );
                for k in 1..order {
                    let mfact = coeff.powi(k);
                    let g = v0.get(k) * mfact;
                    v0.set(k, g);
                    let g = v0.get(k + order) * mfact;
                    v0.set(k + order, g);
                    let g = v1.get(k) * mfact;
                    v1.set(k, g);
                    let g = v1.get(k + order) * mfact;
                    v1.set(k + order, g);
                    let g = v2.get(k) * mfact;
                    v2.set(k, g);
                    let g = v2.get(k + order) * mfact;
                    v2.set(k + order, g);
                }
                let mfact = 1.0 / coeff;
                let mfact1 = mfact / coeff;
                for k in 0..=mx_deg {
                    let g = v1.get(k) * mfact;
                    v1.set(k, g);
                    let g = v2.get(k) * mfact1;
                    v2.set(k, g);
                }
            }

            npass += 1; // NPass++

            let mut j = nb_dim * (point - self.my_first_point);
            let mut n0 = npass;
            let mut curdim = 0;
            for _pnt in 1..=self.my_nb_p3d {
                let mut index_of_constraint = n0;
                for k in 1..=3 {
                    curdim += 1;
                    a.add_constraint(
                        index_of_constraint,
                        curel,
                        curdim,
                        &v0,
                        self.my_tab_points.value(j + k),
                    );
                    index_of_constraint += ngpc1;
                }
                j += 3;
                n0 += ng3d;
            }

            n0 = npass + nbeg2d;
            for _pnt in 1..=self.my_nb_p2d {
                let mut index_of_constraint = n0;
                for k in 1..=2 {
                    curdim += 1;
                    a.add_constraint(
                        index_of_constraint,
                        curel,
                        curdim,
                        &v0,
                        self.my_tab_points.value(j + k),
                    );
                    index_of_constraint += ngpc1;
                }
                j += 2;
                n0 += ng2d;
            }

            /*    if(TypOfConstr == 1) { ... } (disabled OCCT block) */
            if typ_of_constr == 1 {

                npass += 1;
                n0 = npass;
                j = 2 * nb_dim * (i - 1);
                curdim = 0;
                for _pnt in 1..=self.my_nb_p3d {
                    let mut index_of_constraint = n0;
                    for k in 1..=3 {
                        curdim += 1;
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v1,
                            cblong * self.my_tab_constraints.value(j + k),
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng3d;
                    j += 6;
                }

                n0 = npass + nbeg2d;
                for _pnt in 1..=self.my_nb_p2d {
                    let mut index_of_constraint = n0;
                    for k in 1..=2 {
                        curdim += 1;
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v1,
                            cblong * self.my_tab_constraints.value(j + k),
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng2d;
                    j += 4;
                }
            }
            if typ_of_constr == 2 {

                npass += 1;
                n0 = npass;
                j = 2 * nb_dim * (i - 1);
                curdim = 0;
                for _pnt in 1..=self.my_nb_p3d {
                    let mut index_of_constraint = n0;
                    for k in 1..=3 {
                        curdim += 1;
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v1,
                            cblong * self.my_tab_constraints.value(j + k),
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng3d;
                    j += 6;
                }

                n0 = npass + nbeg2d;
                for _pnt in 1..=self.my_nb_p2d {
                    let mut index_of_constraint = n0;
                    for k in 1..=2 {
                        curdim += 1;
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v1,
                            cblong * self.my_tab_constraints.value(j + k),
                        );
                        index_of_constraint += ngpc1;
                    }
                    n0 += ng2d;
                    j += 4;
                }

                j = 2 * nb_dim * (i - 1) + 3;
                let mut jt = ntheta * (i - 1);
                let mut index_of_constraint = ntang3d + 1;
                curdim = 0;
                for _pnt in 1..=self.my_nb_p3d {
                    let mut r1 = 0.0;
                    let mut r2 = 0.0;
                    for k in 1..=3 {
                        r1 += self.my_tab_constraints.value(j + k) * self.my_ttheta.value(jt + k);
                        r2 += self.my_tab_constraints.value(j + k)
                            * self.my_ttheta.value(jt + 3 + k);
                    }
                    r1 *= cblong * cblong;
                    r2 *= cblong * cblong;
                    for k in 1..=3 {
                        curdim += 1;
                        if k > 1 {
                            r1 = 0.0;
                            r2 = 0.0;
                        }
                        // A.AddConstraint(..., myTfthet->Value(jt + k) * V2, R1);
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v_mul(&v2, self.my_tfthet.value(jt + k)),
                            r1,
                        );
                        a.add_constraint(
                            index_of_constraint + 1,
                            curel,
                            curdim,
                            &v_mul(&v2, self.my_tfthet.value(jt + 3 + k)),
                            r2,
                        );
                    }
                    index_of_constraint += ng3d;
                    j += 6;
                    jt += 6;
                }

                j -= 1;
                let mut index_of_constraint = nbeg2d + ntang2d + 1;
                for _pnt in 1..=self.my_nb_p2d {
                    let mut r1 = 0.0;
                    for k in 1..=2 {
                        r1 += self.my_tab_constraints.value(j + k) * self.my_ttheta.value(jt + k);
                    }
                    r1 *= cblong * cblong;
                    for k in 1..=2 {
                        curdim += 1;
                        if k > 1 {
                            r1 = 0.0;
                        }
                        a.add_constraint(
                            index_of_constraint,
                            curel,
                            curdim,
                            &v_mul(&v2, self.my_tfthet.value(jt + k)),
                            r1,
                        );
                    }
                    index_of_constraint += ng2d;
                    j += 4;
                    jt += 2;
                }

                ntang3d += 2;
                ntang2d += 1;
            }
        }
    }

    /// OCCT AppDef_Variational::InitTthetaF (cxx L3315-3403).
    pub(super) fn init_ttheta_f(
        &mut self,
        ndimen: i32,
        typcon: crate::geomalgo::approx_int::AppParConstraint,
        begin: i32,
        jndex: i32,
    ) -> bool {
        if !(2..=3).contains(&ndimen) {
            return false;
        }
        let mut t = DVec3::ZERO;
        let mut v = DVec3::ZERO;
        // OCCT: gp_Vec theta2 default-constructed; assigned in the 3d branch.
        let mut theta2 = DVec3::ZERO;

        if typcon == crate::geomalgo::approx_int::AppParConstraint::TangencyPoint
            || typcon == crate::geomalgo::approx_int::AppParConstraint::CurvaturePoint
        {
            t.x = self.my_tab_constraints.value(jndex);
            t.y = self.my_tab_constraints.value(jndex + 1);
            if ndimen == 3 {
                t.z = self.my_tab_constraints.value(jndex + 2);
            } else {
                t.z = 0.0;
            }
            if ndimen == 2 {
                v.x = 0.0;
                v.y = 0.0;
                v.z = 1.0;
            }
            if ndimen == 3 {
                if !not_parallel(t, &mut v) {
                    return false;
                }
            }
            let theta1 = v.cross(t).normalize();
            self.my_ttheta.set_value(begin, theta1.x);
            self.my_ttheta.set_value(begin + 1, theta1.y);
            if ndimen == 3 {
                theta2 = t.cross(theta1).normalize();
                self.my_ttheta.set_value(begin + 2, theta1.z);
                self.my_ttheta.set_value(begin + 3, theta2.x);
                self.my_ttheta.set_value(begin + 4, theta2.y);
                self.my_ttheta.set_value(begin + 5, theta2.z);
            }

            // Calculation of myTfthet
            if typcon == crate::geomalgo::approx_int::AppParConstraint::CurvaturePoint {
                let xx = t.x * t.x; // std::pow(T.X(), 2)
                let xy = t.x * t.y;
                let yy = t.y * t.y;
                if ndimen == 2 {
                    let f = DVec3::new(
                        yy * theta1.x - xy * theta1.y,
                        xx * theta1.y - xy * theta1.x,
                        0.0,
                    );
                    self.my_tfthet.set_value(begin, f.x);
                    self.my_tfthet.set_value(begin + 1, f.y);
                }
                if ndimen == 3 {
                    let xz = t.x * t.z;
                    let yz = t.y * t.z;
                    let zz = t.z * t.z;

                    let mut f = DVec3::new(
                        (zz + yy) * theta1.x - xy * theta1.y - xz * theta1.z,
                        (xx + zz) * theta1.y - xy * theta1.x - yz * theta1.z,
                        (xx + yy) * theta1.z - xz * theta1.x - yz * theta1.y,
                    );
                    self.my_tfthet.set_value(begin, f.x);
                    self.my_tfthet.set_value(begin + 1, f.y);
                    self.my_tfthet.set_value(begin + 2, f.z);
                    f = DVec3::new(
                        (zz + yy) * theta2.x - xy * theta2.y - xz * theta2.z,
                        (xx + zz) * theta2.y - xy * theta2.x - yz * theta2.z,
                        (xx + yy) * theta2.z - xz * theta2.x - yz * theta2.y,
                    );
                    self.my_tfthet.set_value(begin + 3, f.x);
                    self.my_tfthet.set_value(begin + 4, f.y);
                    self.my_tfthet.set_value(begin + 5, f.z);
                }
            }
        }
        true
    }

    /// OCCT math_Vector Pnt((double*)&myTabPoints->Value(adr), 1,
    /// myDimension) - a view over myDimension cells of myTabPoints starting
    /// at adr; the Rust mirror copies the cells into a local Vector (the
    /// OCCT code only reads the view).
    fn tab_points_vector(&self, adr: i32) -> Vector {
        let mut v = Vector::new(1, self.my_dimension);
        for i in 1..=self.my_dimension {
            v.set(i, self.my_tab_points.value(adr + i - 1));
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::super::{AppDefVariational, CoupleArray1};
    use crate::geomalgo::app_def::MultiLine;
    use crate::geomalgo::approx_int::{AppParConstraint, ConstraintCouple};
    use glam::DVec3;
    use rcad_kernel::math::math_matrix::Matrix;
    use rcad_kernel::math::GeomAbsShape;

    fn pass_point_couples(nb_point: i32) -> CoupleArray1 {
        // The AppBlend_AppSurf.gxx construction shape (L475-485) with every
        // constraint set to PassPoint.
        let mut cc = CoupleArray1::new(1, nb_point);
        for i in 1..=nb_point {
            cc.set_value(
                i,
                ConstraintCouple {
                    index: i,
                    constraint: AppParConstraint::PassPoint,
                },
            );
        }
        cc
    }

    /// Over-constraint detection, hand-derived from the OCCT Init formula
    /// (AppDef_Variational.cxx L462-478).  4 pass points, MaxDegree = 1,
    /// Continuity C2 (NivCont = 2), WithCutting = true (MaxSeg = MaxSegment
    /// = 100):
    ///   (MaxDegree - NivCont) * MaxSeg - NbPass - 2*Tang - 3*Curv
    ///   = (1 - 2) * 100 - 4 - 0 - 0 = -104 < 0
    ///   -> IsOverConstrained = true, IsCreated = false.
    /// With MaxDegree = 2: (2 - 2) * 100 - 4 = -4 < 0, so SetMaxDegree(2)
    /// must be rejected (returns false and keeps MaxDegree = 1 - observable
    /// through max_degree()).
    #[test]
    fn variational_overconstraint_detection_hand() {
        let pts = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
            DVec3::new(6.0, 0.0, 0.0),
        ];
        let ml = MultiLine::new_tab_p3d(&pts);
        let cc = pass_point_couples(4);
        let mut variation = AppDefVariational::new_full(
            &ml, 1, 4, cc, 1, // MaxDegree
            100, // MaxSegment
            GeomAbsShape::C2,
            false, true, 1.0, 2,
        );

        assert!(variation.is_over_constrained());
        assert!(!variation.is_created());

        // SetMaxDegree(2): ((2 - 2) * 100 - 4) < 0 -> rejected.
        assert!(!variation.set_max_degree(2));
        assert_eq!(variation.max_degree(), 1);
    }

    /// Degenerate configuration in which the variational smoothing reduces
    /// to interpolation, hand-derived.  4 collinear points (non-uniform
    /// chord spacing 1, 2, 3), every point a PassPoint (a hard constraint
    /// assembled by AssemblingConstraints + Assembly::Solve), OCCT default
    /// fields (MaxDegree = 14, MaxSegment = 100, C2, Tolerance = 1).
    ///   - OverConstraint formula: (14 - 2) * 100 - 4 = 1196 >= 0 -> created;
    ///   - InitParameters (cxx L2209-2251): total chord = 6, so the
    ///     parameters are exactly (0, 1/6, 3/6, 1);
    ///   - the degree-14 HermitJacobi element contains the straight line, so
    ///     the constrained least-squares minimizer interpolates every data
    ///     point: MaxError / QuadraticError / AverageError / Distance ~ 0 up
    ///     to the linear-solver accuracy.
    #[test]
    fn variational_collinear_pass_points_interpolate() {
        let pts = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
            DVec3::new(6.0, 0.0, 0.0),
        ];
        let ml = MultiLine::new_tab_p3d(&pts);
        let cc = pass_point_couples(4);
        let mut variation = AppDefVariational::new(&ml, 1, 4, cc);

        assert!(variation.is_created());
        assert!(!variation.is_over_constrained());

        variation.approximate();
        assert!(variation.is_done());

        assert!(
            variation.max_error() < 1.0e-6,
            "max error {}",
            variation.max_error()
        );
        assert!(
            variation.quadratic_error() < 1.0e-12,
            "quadratic error {}",
            variation.quadratic_error()
        );
        assert!(
            variation.average_error() < 1.0e-6,
            "average error {}",
            variation.average_error()
        );

        // Chord-length parameters: (0, 1/6, 1/2, 1).
        let params = variation.parameters();
        assert!((params.value(1) - 0.0).abs() < 1.0e-12);
        assert!((params.value(2) - 1.0 / 6.0).abs() < 1.0e-9);
        assert!((params.value(3) - 0.5).abs() < 1.0e-9);
        assert!((params.value(4) - 1.0).abs() < 1.0e-12);

        // The result knots span [0, 1] (OCCT ctor L141-143 + Approximate).
        let knots = variation.knots();
        assert!((knots.value(knots.lower()) - 0.0).abs() < 1.0e-12);
        assert!((knots.value(knots.upper()) - 1.0).abs() < 1.0e-12);

        // Distance matrix (cxx L735-787): every entry must vanish.
        let mut mat = Matrix::new(1, 3, 1, 4);
        variation.distance(&mut mat);
        for i in 1..=3 {
            for j in 1..=4 {
                assert!(
                    mat.get(i, j) < 1.0e-6,
                    "distance({}, {}) = {}",
                    i,
                    j,
                    mat.get(i, j)
                );
            }
        }

        // MaxErrorIndex is a point index in [1, 4].
        let idx = variation.max_error_index();
        assert!((1..=4).contains(&idx), "max error index {}", idx);
    }

    // NOTE on the deliberately absent two-point case: with NbPoints = 2 the
    // OCCT InitCriterionEstimations (cxx L2270-2276) unconditionally calls
    // EstTangent(myFirstPoint + 2), whose middle-point branch reads
    // myTabPoints cells beyond the storage end - AppDef_Variational.cxx
    // disables the range checks (`No_Standard_OutOfRange`, cxx L27), so the
    // OCCT behavior there is undefined-memory-read-dependent and cannot be
    // reproduced in safe Rust (nor hand-derived).  See the batch report.
}
