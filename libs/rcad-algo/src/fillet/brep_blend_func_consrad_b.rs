//! OCCT BlendFunc_ConstRad — part 2: ComputeValues (BlendFunc_ConstRad.cxx
//! L136-776) and the four Section overloads (L1141-1199 obsolete,
//! L1280-1351, L1355-1524, L1528-1903).  The struct and the remaining
//! methods live in [`super::brep_blend_func_consrad`].

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::geom::{Circle3, CurveEval as _, SurfaceEval as _};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{MatD, VecD};

use super::brep_blend_func::{blend_func_compute_dnormal, blend_func_compute_normal};
use super::brep_blend_func_consrad::EPS;
use super::brep_blend_func_consrad::{geomfill_get_circle_pending, BlendFuncConstRad};
use super::brep_blend_point::BlendPoint;

impl<'a> BlendFuncConstRad<'a> {
    /// OCCT ComputeValues(X, Order, ByParam, Param)
    /// (BlendFunc_ConstRad.cxx L136-776) — OBLIGATORY passage for all
    /// calculations: positions on Surfaces and Curve, computes the equations
    /// and their partial derivatives, stocks intermediate results in fields.
    pub(crate) fn compute_values(&mut self, x: &[f64], order: i32, by_param: bool, param: f64) -> bool {
        // OCCT L143-146: the `static` locals (d3u1... d3gui, ptgui,
        // invnormtg, dinvnormtg) are a reallocation optimisation; the cache
        // subset (ptgui, d1gui, d2gui, invnormtg) lives on the struct, the
        // rest are plain locals below.
        let mut t = param;
        let aux;

        // Case of implicite parameter
        if !by_param {
            t = self.param;
        }

        // Is the work already done ?
        let mut my_x_ok = order <= self.my_x_order;
        let mut ii = 1;
        while (ii <= 4) && my_x_ok {
            my_x_ok = x[ii - 1] == self.xval[ii - 1];
            ii += 1;
        }

        let t_ok = (t == self.tval) && ((order <= self.my_t_order) || (!by_param));

        if my_x_ok && t_ok {
            return true;
        }

        // Processing of t
        if !t_ok {
            self.tval = t;
            if by_param {
                self.my_t_order = order;
            } else {
                self.my_t_order = 0;
            }
            //----- Positioning on the curve ----------------
            match self.my_t_order {
                0 => {
                    // OCCT: tcurv->D1(T, ptgui, d1gui); nplan = d1gui.Normalized();
                    self.ptgui = self.tcurv().point_at(t);
                    self.d1gui = self.tcurv().derivative_at(t);
                    self.nplan = self.d1gui.normalize_or_zero();
                }
                1 => {
                    // OCCT: tcurv->D2(T, ptgui, d1gui, d2gui);
                    self.ptgui = self.tcurv().point_at(t);
                    self.d1gui = self.tcurv().derivative_at(t);
                    self.d2gui = self.tcurv().derivative2_at(t);
                    self.nplan = self.d1gui.normalize_or_zero();
                    self.invnormtg = 1.0 / self.d1gui.length();
                    // OCCT: dnplan.SetLinearForm(invnormtg, d2gui,
                    //                           -invnormtg * (nplan.Dot(d2gui)), nplan);
                    self.dnplan = self.invnormtg * self.d2gui
                        + (-self.invnormtg * self.nplan.dot(self.d2gui)) * self.nplan;
                }
                2 => {
                    // OCCT: tcurv->D3(T, ptgui, d1gui, d2gui, d3gui);
                    self.ptgui = self.tcurv().point_at(t);
                    self.d1gui = self.tcurv().derivative_at(t);
                    self.d2gui = self.tcurv().derivative2_at(t);
                    let d3gui = self.tcurv().derivative3_at(t);
                    self.nplan = self.d1gui.normalize_or_zero();
                    self.invnormtg = 1.0 / self.d1gui.length();
                    self.dnplan = self.invnormtg * self.d2gui
                        + (-self.invnormtg * self.nplan.dot(self.d2gui)) * self.nplan;
                    // OCCT: dinvnormtg = -nplan.Dot(d2gui) * invnormtg * invnormtg;
                    let dinvnormtg = -self.nplan.dot(self.d2gui) * self.invnormtg * self.invnormtg;
                    // OCCT: d2nplan.SetLinearForm(invnormtg, d3gui, dinvnormtg, d2gui);
                    let mut d2nplan = self.invnormtg * d3gui + dinvnormtg * self.d2gui;
                    aux = dinvnormtg * (self.nplan.dot(self.d2gui))
                        + self.invnormtg * (self.dnplan.dot(self.d2gui) + self.nplan.dot(d3gui));
                    // OCCT: d2nplan.SetLinearForm(-invnormtg * (nplan.Dot(d2gui)), dnplan,
                    //                             -aux, nplan, d2nplan);
                    d2nplan = (-self.invnormtg * self.nplan.dot(self.d2gui)) * self.dnplan
                        + (-aux) * self.nplan
                        + d2nplan;
                    self.d2nplan = d2nplan;
                }
                _ => return false,
            }
        }

        // Processing of X
        if !my_x_ok {
            self.xval = [x[0], x[1], x[2], x[3]];
            self.my_x_order = order;
            //-------------- Positioning on surfaces -----------------
            match self.my_x_order {
                0 => {
                    // OCCT: surf1->D1(X(1), X(2), pts1, d1u1, d1v1);
                    let (p, du, dv) = self.surf1.derivatives(x[0], x[1]);
                    self.pts1 = p;
                    self.d1u1 = du;
                    self.d1v1 = dv;
                    self.nsurf1 = self.d1u1.cross(self.d1v1);
                    // OCCT: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
                    let (p, du, dv) = self.surf2.derivatives(x[2], x[3]);
                    self.pts2 = p;
                    self.d1u2 = du;
                    self.d1v2 = dv;
                    self.nsurf2 = self.d1u2.cross(self.d1v2);
                }
                1 => {
                    // OCCT: surf1->D2(X(1), X(2), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf1.derivatives2(x[0], x[1]);
                    self.pts1 = p;
                    self.d1u1 = du;
                    self.d1v1 = dv;
                    self.d2u1 = d2u;
                    self.d2uv1 = d2uv;
                    self.d2v1 = d2v;
                    self.nsurf1 = self.d1u1.cross(self.d1v1);
                    // OCCT: dns1u1 = d2u1.Crossed(d1v1).Added(d1u1.Crossed(d2uv1));
                    self.dns1u1 = self.d2u1.cross(self.d1v1) + self.d1u1.cross(self.d2uv1);
                    self.dns1v1 = self.d2uv1.cross(self.d1v1) + self.d1u1.cross(self.d2v1);
                    // OCCT: surf2->D2(X(3), X(4), pts2, d1u2, d1v2, d2u2, d2v2, d2uv2);
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf2.derivatives2(x[2], x[3]);
                    self.pts2 = p;
                    self.d1u2 = du;
                    self.d1v2 = dv;
                    self.d2u2 = d2u;
                    self.d2uv2 = d2uv;
                    self.d2v2 = d2v;
                    self.nsurf2 = self.d1u2.cross(self.d1v2);
                    self.dns1u2 = self.d2u2.cross(self.d1v2) + self.d1u2.cross(self.d2uv2);
                    self.dns1v2 = self.d2uv2.cross(self.d1v2) + self.d1u2.cross(self.d2v2);
                }
                2 => {
                    // OCCT: surf1->D3(X(1), X(2), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1,
                    //                  d3u1, d3v1, d3uuv1, d3uvv1);
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf1.derivatives2(x[0], x[1]);
                    self.pts1 = p;
                    self.d1u1 = du;
                    self.d1v1 = dv;
                    self.d2u1 = d2u;
                    self.d2uv1 = d2uv;
                    self.d2v1 = d2v;
                    self.nsurf1 = self.d1u1.cross(self.d1v1);
                    self.dns1u1 = self.d2u1.cross(self.d1v1) + self.d1u1.cross(self.d2uv1);
                    self.dns1v1 = self.d2uv1.cross(self.d1v1) + self.d1u1.cross(self.d2v1);

                    // OCCT: surf2->D3(X(3), X(4), pts2, d1u2, d1v2, d2u2, d2v2, d2uv2,
                    //                  d3u2, d3v2, d3uuv2, d3uvv2);
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf2.derivatives2(x[2], x[3]);
                    self.pts2 = p;
                    self.d1u2 = du;
                    self.d1v2 = dv;
                    self.d2u2 = d2u;
                    self.d2uv2 = d2uv;
                    self.d2v2 = d2v;
                    self.nsurf2 = self.d1u2.cross(self.d1v2);
                    self.dns1u2 = self.d2u2.cross(self.d1v2) + self.d1u2.cross(self.d2uv2);
                    self.dns1v2 = self.d2uv2.cross(self.d1v2) + self.d1u2.cross(self.d2v2);
                }
                _ => return false,
            }
            // Case of degenerated surfaces
            if self.nsurf1.length() < EPS {
                let p = (x[0], x[1]);
                if order == 0 {
                    // OCCT: BlendFunc::ComputeNormal(surf1, P, nsurf1);
                    let mut n = DVec3::ZERO;
                    blend_func_compute_normal(self.surf1, p, &mut n);
                    self.nsurf1 = n;
                } else {
                    // OCCT: BlendFunc::ComputeDNormal(surf1, P, nsurf1, dns1u1, dns1v1);
                    let mut n = DVec3::ZERO;
                    let mut dnu = DVec3::ZERO;
                    let mut dnv = DVec3::ZERO;
                    blend_func_compute_dnormal(self.surf1, p, &mut n, &mut dnu, &mut dnv);
                    self.nsurf1 = n;
                    self.dns1u1 = dnu;
                    self.dns1v1 = dnv;
                }
            }
            if self.nsurf2.length() < EPS {
                let p = (x[2], x[3]);
                if order == 0 {
                    // OCCT: BlendFunc::ComputeNormal(surf2, P, nsurf2);
                    let mut n = DVec3::ZERO;
                    blend_func_compute_normal(self.surf2, p, &mut n);
                    self.nsurf2 = n;
                } else {
                    // OCCT: BlendFunc::ComputeDNormal(surf2, P, nsurf2, dns1u2, dns1v2);
                    let mut n = DVec3::ZERO;
                    let mut dnu = DVec3::ZERO;
                    let mut dnv = DVec3::ZERO;
                    blend_func_compute_dnormal(self.surf2, p, &mut n, &mut dnu, &mut dnv);
                    self.nsurf2 = n;
                    self.dns1u2 = dnu;
                    self.dns1v2 = dnv;
                }
            }
        }

        // -------------------- Positioning of order 0 ---------------------
        let the_d = -self.nplan.dot(self.ptgui);

        // OCCT L289-292: E(1) is computed componentwise — the exact operation
        // order is preserved.
        self.e[0] = (self.nplan.x * (self.pts1.x + self.pts2.x)
            + self.nplan.y * (self.pts1.y + self.pts2.y)
            + self.nplan.z * (self.pts1.z + self.pts2.z))
            / 2.0
            + the_d;

        let ncrossns1 = self.nplan.cross(self.nsurf1);
        let ncrossns2 = self.nplan.cross(self.nsurf2);

        let mut invnorm1 = ncrossns1.length();
        let mut invnorm2 = ncrossns2.length();

        if invnorm1 > EPS {
            invnorm1 = 1.0 / invnorm1;
        } else {
            invnorm1 = 1.0; // Unsatisfactory, but it is not necessary to crash
        }
        if invnorm2 > EPS {
            invnorm2 = 1.0 / invnorm2;
        } else {
            invnorm2 = 1.0; // Unsatisfactory, but it is not necessary to crash
        }

        let ndotns1 = self.nplan.dot(self.nsurf1);
        let ndotns2 = self.nplan.dot(self.nsurf2);

        // OCCT: temp.SetLinearForm(ndotns1, nplan, -1., nsurf1); temp.Multiply(invnorm1);
        let mut temp = (ndotns1 * self.nplan + (-1.0) * self.nsurf1) * invnorm1;

        // OCCT: resul.SetLinearForm(ray1, temp, gp_Vec(pts2, pts1));
        let mut resul = self.ray1 * temp + (self.pts1 - self.pts2);
        // OCCT: temp.SetLinearForm(ndotns2, nplan, -1., nsurf2); temp.Multiply(invnorm2);
        temp = (ndotns2 * self.nplan + (-1.0) * self.nsurf2) * invnorm2;
        // OCCT: resul.Subtract(ray2 * temp);
        resul -= self.ray2 * temp;

        self.e[1] = resul.x;
        self.e[2] = resul.y;
        self.e[3] = resul.z;

        // -------------------- Positioning of order 1 ---------------------
        if order >= 1 {
            self.dedx[0][0] = self.nplan.dot(self.d1u1) / 2.0;
            self.dedx[0][1] = self.nplan.dot(self.d1v1) / 2.0;
            self.dedx[0][2] = self.nplan.dot(self.d1u2) / 2.0;
            self.dedx[0][3] = self.nplan.dot(self.d1v2) / 2.0;

            let mut cube = invnorm1 * invnorm1 * invnorm1;
            // Derived in relation to u1
            let mut grosterme = -ncrossns1.dot(self.nplan.cross(self.dns1u1)) * cube;
            // OCCT: dndu1.SetLinearForm(grosterme * ndotns1 + invnorm1 * nplan.Dot(dns1u1),
            //                           nplan, -grosterme, nsurf1, -invnorm1, dns1u1);
            self.dndu1 = (grosterme * ndotns1 + invnorm1 * self.nplan.dot(self.dns1u1))
                * self.nplan
                + (-grosterme) * self.nsurf1
                + (-invnorm1) * self.dns1u1;

            // OCCT: resul.SetLinearForm(ray1, dndu1, d1u1);
            let resul = self.ray1 * self.dndu1 + self.d1u1;
            self.dedx[1][0] = resul.x;
            self.dedx[2][0] = resul.y;
            self.dedx[3][0] = resul.z;

            // Derived in relation to v1

            grosterme = -ncrossns1.dot(self.nplan.cross(self.dns1v1)) * cube;
            self.dndv1 = (grosterme * ndotns1 + invnorm1 * self.nplan.dot(self.dns1v1))
                * self.nplan
                + (-grosterme) * self.nsurf1
                + (-invnorm1) * self.dns1v1;

            let resul = self.ray1 * self.dndv1 + self.d1v1;
            self.dedx[1][1] = resul.x;
            self.dedx[2][1] = resul.y;
            self.dedx[3][1] = resul.z;

            cube = invnorm2 * invnorm2 * invnorm2;
            // Derived in relation to u2
            grosterme = -ncrossns2.dot(self.nplan.cross(self.dns1u2)) * cube;
            self.dndu2 = (grosterme * ndotns2 + invnorm2 * self.nplan.dot(self.dns1u2))
                * self.nplan
                + (-grosterme) * self.nsurf2
                + (-invnorm2) * self.dns1u2;

            // OCCT: resul.SetLinearForm(-ray2, dndu2, -1, d1u2);
            let resul = -self.ray2 * self.dndu2 + (-1.0) * self.d1u2;
            self.dedx[1][2] = resul.x;
            self.dedx[2][2] = resul.y;
            self.dedx[3][2] = resul.z;

            // Derived in relation to v2
            grosterme = -ncrossns2.dot(self.nplan.cross(self.dns1v2)) * cube;
            self.dndv2 = (grosterme * ndotns2 + invnorm2 * self.nplan.dot(self.dns1v2))
                * self.nplan
                + (-grosterme) * self.nsurf2
                + (-invnorm2) * self.dns1v2;

            let resul = -self.ray2 * self.dndv2 + (-1.0) * self.d1v2;
            self.dedx[1][3] = resul.x;
            self.dedx[2][3] = resul.y;
            self.dedx[3][3] = resul.z;

            if by_param {
                // OCCT: temp.SetXYZ((pts1.XYZ() + pts2.XYZ()) / 2 - ptgui.XYZ());
                let temp = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
                // Derived from n1 in relation to w
                grosterme = ncrossns1.dot(self.dnplan.cross(self.nsurf1)) * invnorm1 * invnorm1;
                // OCCT: dn1w.SetLinearForm((dnplan.Dot(nsurf1) - grosterme * ndotns1) * invnorm1,
                //                           nplan, ndotns1 * invnorm1, dnplan,
                //                           grosterme * invnorm1, nsurf1);
                self.dn1w = ((self.dnplan.dot(self.nsurf1) - grosterme * ndotns1) * invnorm1)
                    * self.nplan
                    + (ndotns1 * invnorm1) * self.dnplan
                    + (grosterme * invnorm1) * self.nsurf1;

                // Derivee from n2 in relation to w
                grosterme = ncrossns2.dot(self.dnplan.cross(self.nsurf2)) * invnorm2 * invnorm2;
                self.dn2w = ((self.dnplan.dot(self.nsurf2) - grosterme * ndotns2) * invnorm2)
                    * self.nplan
                    + (ndotns2 * invnorm2) * self.dnplan
                    + (grosterme * invnorm2) * self.nsurf2;

                self.dedt[0] = self.dnplan.dot(temp) - 1.0 / self.invnormtg;
                self.dedt[1] = self.ray1 * self.dn1w.x - self.ray2 * self.dn2w.x;
                self.dedt[2] = self.ray1 * self.dn1w.y - self.ray2 * self.dn2w.y;
                self.dedt[3] = self.ray1 * self.dn1w.z - self.ray2 * self.dn2w.z;
            }
            // ------   Positioning of order 2  -----------------------------
            if order == 2 {
                self.compute_d2(
                    invnorm1, invnorm2, ndotns1, ndotns2, ncrossns1, ncrossns2, by_param,
                );
            }
        }
        true
    }

    /// OCCT BlendFunc_ConstRad::ComputeValues order-2 block
    /// (BlendFunc_ConstRad.cxx L434-772) — extracted verbatim; the
    /// invnorm / ndotns / ncrossns values are the order-1 locals of the
    /// caller.
    #[allow(clippy::too_many_arguments)]
    fn compute_d2(
        &mut self,
        invnorm1: f64,
        invnorm2: f64,
        ndotns1: f64,
        ndotns2: f64,
        ncrossns1: DVec3,
        ncrossns2: DVec3,
        by_param: bool,
    ) {
        //     gp_Vec d2ndu1,  d2ndu2, d2ndv1, d2ndv2, d2nduv1, d2nduv2;
        self.d2edx2.init(0.0);

        *self.d2edx2.change_value(1, 1, 1) = self.nplan.dot(self.d2u1) / 2.0;
        *self.d2edx2.change_value(1, 2, 1) = self.nplan.dot(self.d2uv1) / 2.0;
        *self.d2edx2.change_value(1, 1, 2) = self.nplan.dot(self.d2uv1) / 2.0;
        *self.d2edx2.change_value(1, 2, 2) = self.nplan.dot(self.d2v1) / 2.0;

        *self.d2edx2.change_value(1, 3, 3) = self.nplan.dot(self.d2u2) / 2.0;
        *self.d2edx2.change_value(1, 4, 3) = self.nplan.dot(self.d2uv2) / 2.0;
        *self.d2edx2.change_value(1, 3, 4) = self.nplan.dot(self.d2uv2) / 2.0;
        *self.d2edx2.change_value(1, 4, 4) = self.nplan.dot(self.d2v2) / 2.0;
        // ================
        // ==  Surface 1 ==
        // ================
        let mut carre = invnorm1 * invnorm1;
        let mut cube = carre * invnorm1;
        // Derived double compared to u1
        // Derived from the norm
        // OCCT: d2ns1u1.SetLinearForm(1, d3u1.Crossed(d1v1), 2, d2u1.Crossed(d2uv1),
        //                             1, d1u1.Crossed(d3uuv1));
        let (d3u1, d3v1, d3uuv1, d3uvv1) = self.surf1_d3();
        let d2ns1u1 = d3u1.cross(self.d1v1)
            + 2.0 * (self.d2u1.cross(self.d2uv1))
            + self.d1u1.cross(d3uuv1);
        let dprim = ncrossns1.dot(self.nplan.cross(self.dns1u1));
        let smallterm = -2.0 * dprim * cube;
        // OCCT: DSecn = ncrossns1.Dot(nplan.Crossed(d2ns1u1))
        //       + (nplan.Crossed(dns1u1)).SquareMagnitude();
        let dsecn = ncrossns1.dot(self.nplan.cross(d2ns1u1))
            + self.nplan.cross(self.dns1u1).length_squared();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
        let p1 = self.nplan.dot(self.dns1u1);
        let p2 = self.nplan.dot(d2ns1u1);
        // OCCT: d2ndu1.SetLinearForm(invnorm1 * p2 + smallterm * p1, nplan,
        //                           -smallterm, dns1u1, -invnorm1, d2ns1u1);
        let mut d2ndu1 = (invnorm1 * p2 + smallterm * p1) * self.nplan
            + (-smallterm) * self.dns1u1
            + (-invnorm1) * d2ns1u1;
        d2ndu1 += temp;
        // OCCT: resul.SetLinearForm(ray1, d2ndu1, d2u1);
        let resul = self.ray1 * d2ndu1 + self.d2u1;
        *self.d2edx2.change_value(2, 1, 1) = resul.x;
        *self.d2edx2.change_value(3, 1, 1) = resul.y;
        *self.d2edx2.change_value(4, 1, 1) = resul.z;

        // Derived double compared to u1, v1
        // Derived from the norm
        // OCCT: d2ns1uv1 = d3uuv1.Crossed(d1v1) + d2u1.Crossed(d2v1) + d1u1.Crossed(d3uvv1);
        let d2ns1uv1 =
            d3uuv1.cross(self.d1v1) + self.d2u1.cross(self.d2v1) + self.d1u1.cross(d3uvv1);
        let mut uterm = ncrossns1.dot(self.nplan.cross(self.dns1u1));
        let mut vterm = ncrossns1.dot(self.nplan.cross(self.dns1v1));
        let dsecn = self.nplan.cross(self.dns1v1).dot(self.nplan.cross(self.dns1u1))
            + ncrossns1.dot(self.nplan.cross(d2ns1uv1));
        let grosterme = (3.0 * uterm * vterm * carre - dsecn) * cube;
        uterm *= -cube; // and only now
        vterm *= -cube;

        let p1 = self.nplan.dot(self.dns1u1);
        let p2 = self.nplan.dot(self.dns1v1);
        // OCCT: temp.SetLinearForm(grosterme * ndotns1, nplan, -grosterme, nsurf1,
        //                         -invnorm1, d2ns1uv1);
        let temp = (grosterme * ndotns1) * self.nplan
            + (-grosterme) * self.nsurf1
            + (-invnorm1) * d2ns1uv1;
        // OCCT: d2nduv1.SetLinearForm(invnorm1 * nplan.Dot(d2ns1uv1) + uterm * p2
        //                             + vterm * p1, nplan, -uterm, dns1v1, -vterm, dns1u1);
        let mut d2nduv1 = (invnorm1 * self.nplan.dot(d2ns1uv1) + uterm * p2 + vterm * p1)
            * self.nplan
            + (-uterm) * self.dns1v1
            + (-vterm) * self.dns1u1;

        d2nduv1 += temp;
        // OCCT: resul.SetLinearForm(ray1, d2nduv1, d2uv1);
        let resul = self.ray1 * d2nduv1 + self.d2uv1;

        *self.d2edx2.change_value(2, 2, 1) = resul.x;
        *self.d2edx2.change_value(2, 1, 2) = resul.x;
        *self.d2edx2.change_value(3, 2, 1) = resul.y;
        *self.d2edx2.change_value(3, 1, 2) = resul.y;
        *self.d2edx2.change_value(4, 2, 1) = resul.z;
        *self.d2edx2.change_value(4, 1, 2) = resul.z;

        // Derived double compared to v1
        // Derived from the norm
        // OCCT: d2ns1v1.SetLinearForm(1, d1u1.Crossed(d3v1), 2, d2uv1.Crossed(d2v1),
        //                             1, d3uvv1.Crossed(d1v1));
        let d2ns1v1 = self.d1u1.cross(d3v1)
            + 2.0 * (self.d2uv1.cross(self.d2v1))
            + d3uvv1.cross(self.d1v1);
        let dprim = ncrossns1.dot(self.nplan.cross(self.dns1v1));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = ncrossns1.dot(self.nplan.cross(d2ns1v1))
            + self.nplan.cross(self.dns1v1).length_squared();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.nplan.dot(self.dns1v1);
        let p2 = self.nplan.dot(d2ns1v1);
        let temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
        let mut d2ndv1 = (invnorm1 * p2 + smallterm * p1) * self.nplan
            + (-smallterm) * self.dns1v1
            + (-invnorm1) * d2ns1v1;
        d2ndv1 += temp;
        // OCCT: resul.SetLinearForm(ray1, d2ndv1, d2v1);
        let resul = self.ray1 * d2ndv1 + self.d2v1;

        *self.d2edx2.change_value(2, 2, 2) = resul.x;
        *self.d2edx2.change_value(3, 2, 2) = resul.y;
        *self.d2edx2.change_value(4, 2, 2) = resul.z;
        // ================
        // ==  Surface 2 ==
        // ================
        carre = invnorm2 * invnorm2;
        cube = carre * invnorm2;
        // Derived double compared to u2
        // Derived from the norm
        let (d3u2, d3v2, d3uuv2, d3uvv2) = self.surf2_d3();
        let d2ns1u2 = d3u2.cross(self.d1v2)
            + 2.0 * (self.d2u2.cross(self.d2uv2))
            + self.d1u2.cross(d3uuv2);
        let dprim = ncrossns2.dot(self.nplan.cross(self.dns1u2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = ncrossns2.dot(self.nplan.cross(d2ns1u2))
            + self.nplan.cross(self.dns1u2).length_squared();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
        let p1 = self.nplan.dot(self.dns1u2);
        let p2 = self.nplan.dot(d2ns1u2);
        let mut d2ndu2 = (invnorm2 * p2 + smallterm * p1) * self.nplan
            + (-smallterm) * self.dns1u2
            + (-invnorm2) * d2ns1u2;
        d2ndu2 += temp;
        // OCCT: resul.SetLinearForm(-ray2, d2ndu2, -1, d2u2);
        let resul = -self.ray2 * d2ndu2 + (-1.0) * self.d2u2;
        *self.d2edx2.change_value(2, 3, 3) = resul.x;
        *self.d2edx2.change_value(3, 3, 3) = resul.y;
        *self.d2edx2.change_value(4, 3, 3) = resul.z;

        // Derived double compared to u2, v2
        // Derived from the norm
        let d2ns1uv2 =
            d3uuv2.cross(self.d1v2) + self.d2u2.cross(self.d2v2) + self.d1u2.cross(d3uvv2);
        let mut uterm = ncrossns2.dot(self.nplan.cross(self.dns1u2));
        let mut vterm = ncrossns2.dot(self.nplan.cross(self.dns1v2));
        let dsecn = self.nplan.cross(self.dns1v2).dot(self.nplan.cross(self.dns1u2))
            + ncrossns2.dot(self.nplan.cross(d2ns1uv2));
        let grosterme = (3.0 * uterm * vterm * carre - dsecn) * cube;
        uterm *= -cube; // and only now
        vterm *= -cube;

        let p1 = self.nplan.dot(self.dns1u2);
        let p2 = self.nplan.dot(self.dns1v2);
        let temp = (grosterme * ndotns2) * self.nplan
            + (-grosterme) * self.nsurf2
            + (-invnorm2) * d2ns1uv2;
        let mut d2nduv2 = (invnorm2 * self.nplan.dot(d2ns1uv2) + uterm * p2 + vterm * p1)
            * self.nplan
            + (-uterm) * self.dns1v2
            + (-vterm) * self.dns1u2;

        d2nduv2 += temp;
        // OCCT: resul.SetLinearForm(-ray2, d2nduv2, -1, d2uv2);
        let resul = -self.ray2 * d2nduv2 + (-1.0) * self.d2uv2;

        *self.d2edx2.change_value(2, 4, 3) = resul.x;
        *self.d2edx2.change_value(2, 3, 4) = resul.x;
        *self.d2edx2.change_value(3, 4, 3) = resul.y;
        *self.d2edx2.change_value(3, 3, 4) = resul.y;
        *self.d2edx2.change_value(4, 4, 3) = resul.z;
        *self.d2edx2.change_value(4, 3, 4) = resul.z;

        // Derived double compared to v2
        // Derived from the norm
        let d2ns1v2 = self.d1u2.cross(d3v2)
            + 2.0 * (self.d2uv2.cross(self.d2v2))
            + d3uvv2.cross(self.d1v2);
        let dprim = ncrossns2.dot(self.nplan.cross(self.dns1v2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = ncrossns2.dot(self.nplan.cross(d2ns1v2))
            + self.nplan.cross(self.dns1v2).length_squared();
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.nplan.dot(self.dns1v2);
        let p2 = self.nplan.dot(d2ns1v2);
        let temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
        let mut d2ndv2 = (invnorm2 * p2 + smallterm * p1) * self.nplan
            + (-smallterm) * self.dns1v2
            + (-invnorm2) * d2ns1v2;
        d2ndv2 += temp;
        // OCCT: resul.SetLinearForm(-ray2, d2ndv2, -1, d2v2);
        let resul = -self.ray2 * d2ndv2 + (-1.0) * self.d2v2;

        *self.d2edx2.change_value(2, 4, 4) = resul.x;
        *self.d2edx2.change_value(3, 4, 4) = resul.y;
        *self.d2edx2.change_value(4, 4, 4) = resul.z;

        if by_param {
            self.compute_d2_t(invnorm1, invnorm2, ndotns1, ndotns2, ncrossns1, ncrossns2);
        }
    }

    /// OCCT BlendFunc_ConstRad::ComputeValues byParam order-2 block
    /// (BlendFunc_ConstRad.cxx L608-772).
    #[allow(clippy::too_many_arguments)]
    fn compute_d2_t(
        &mut self,
        invnorm1: f64,
        invnorm2: f64,
        ndotns1: f64,
        ndotns2: f64,
        ncrossns1: DVec3,
        ncrossns2: DVec3,
    ) {
        //  ---------- Derivation double in t, X --------------------------
        self.d2edxdt[0][0] = self.dnplan.dot(self.d1u1) / 2.0;
        self.d2edxdt[0][1] = self.dnplan.dot(self.d1v1) / 2.0;
        self.d2edxdt[0][2] = self.dnplan.dot(self.d1u2) / 2.0;
        self.d2edxdt[0][3] = self.dnplan.dot(self.d1v2) / 2.0;

        let carre = invnorm1 * invnorm1;
        let cube = carre * invnorm1;
        //--> Derived compared to u1 and t
        let tterm = ncrossns1.dot(self.dnplan.cross(self.nsurf1));
        let smallterm = -tterm * cube;
        // Derived from the norm
        let mut uterm = ncrossns1.dot(self.nplan.cross(self.dns1u1));
        let dsecn = self.nplan.cross(self.dns1u1).dot(self.dnplan.cross(self.nsurf1))
            + ncrossns1.dot(self.dnplan.cross(self.dns1u1));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(self.nsurf1);
        let p2 = self.nplan.dot(self.dns1u1);
        let p12 = self.dnplan.dot(self.dns1u1);

        // OCCT: d2ndtu1.SetLinearForm(invnorm1 * p12 + smallterm * p2 + uterm * p1
        //                             + grosterme * ndotns1, nplan,
        //                             invnorm1 * p2 + uterm * ndotns1, dnplan,
        //                             -smallterm, dns1u1);
        let mut d2ndtu1 = (invnorm1 * p12 + smallterm * p2 + uterm * p1 + grosterme * ndotns1)
            * self.nplan
            + (invnorm1 * p2 + uterm * ndotns1) * self.dnplan
            + (-smallterm) * self.dns1u1;
        // OCCT: d2ndtu1 -= grosterme * nsurf1;
        d2ndtu1 -= grosterme * self.nsurf1;

        // OCCT: resul = ray1 * d2ndtu1;
        let resul = self.ray1 * d2ndtu1;
        self.d2edxdt[1][0] = resul.x;
        self.d2edxdt[2][0] = resul.y;
        self.d2edxdt[3][0] = resul.z;

        //--> Derived compared to v1 and t
        // Derived from the norm
        uterm = ncrossns1.dot(self.nplan.cross(self.dns1v1));
        let dsecn = self.nplan.cross(self.dns1v1).dot(self.dnplan.cross(self.nsurf1))
            + ncrossns1.dot(self.dnplan.cross(self.dns1v1));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(self.nsurf1);
        let p2 = self.nplan.dot(self.dns1v1);
        let p12 = self.dnplan.dot(self.dns1v1);
        let mut d2ndtv1 = (invnorm1 * p12 + uterm * p1 + smallterm * p2 + grosterme * ndotns1)
            * self.nplan
            + (invnorm1 * p2 + uterm * ndotns1) * self.dnplan
            + (-smallterm) * self.dns1v1;
        d2ndtv1 -= grosterme * self.nsurf1;

        let resul = self.ray1 * d2ndtv1;
        self.d2edxdt[1][1] = resul.x;
        self.d2edxdt[2][1] = resul.y;
        self.d2edxdt[3][1] = resul.z;

        let carre = invnorm2 * invnorm2;
        let cube = carre * invnorm2;
        //--> Derived compared to u2 and t
        let tterm = ncrossns2.dot(self.dnplan.cross(self.nsurf2));
        let smallterm = -tterm * cube;
        // Derived from the norm
        let mut uterm = ncrossns2.dot(self.nplan.cross(self.dns1u2));
        let dsecn = self.nplan.cross(self.dns1u2).dot(self.dnplan.cross(self.nsurf2))
            + ncrossns2.dot(self.dnplan.cross(self.dns1u2));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(self.nsurf2);
        let p2 = self.nplan.dot(self.dns1u2);
        let p12 = self.dnplan.dot(self.dns1u2);

        let mut d2ndtu2 = (invnorm2 * p12 + smallterm * p2 + uterm * p1 + grosterme * ndotns2)
            * self.nplan
            + (invnorm2 * p2 + uterm * ndotns2) * self.dnplan
            + (-smallterm) * self.dns1u2;
        d2ndtu2 -= grosterme * self.nsurf2;

        let resul = -self.ray2 * d2ndtu2;
        self.d2edxdt[1][2] = resul.x;
        self.d2edxdt[2][2] = resul.y;
        self.d2edxdt[3][2] = resul.z;

        //--> Derived compared to v2 and t
        // Derived from the norm
        uterm = ncrossns2.dot(self.nplan.cross(self.dns1v2));
        let dsecn = self.nplan.cross(self.dns1v2).dot(self.dnplan.cross(self.nsurf2))
            + ncrossns2.dot(self.dnplan.cross(self.dns1v2));
        let grosterme = (3.0 * uterm * tterm * carre - dsecn) * cube;
        uterm *= -cube;

        let p1 = self.dnplan.dot(self.nsurf2);
        let p2 = self.nplan.dot(self.dns1v2);
        let p12 = self.dnplan.dot(self.dns1v2);

        let mut d2ndtv2 = (invnorm2 * p12 + smallterm * p2 + uterm * p1 + grosterme * ndotns2)
            * self.nplan
            + (invnorm2 * p2 + uterm * ndotns2) * self.dnplan
            + (-smallterm) * self.dns1v2;
        d2ndtv2 -= grosterme * self.nsurf2;

        let resul = -self.ray2 * d2ndtv2;
        self.d2edxdt[1][3] = resul.x;
        self.d2edxdt[2][3] = resul.y;
        self.d2edxdt[3][3] = resul.z;

        //  ---------- Derivation double in t -----------------------------
        // Derived from n1 compared to w
        let carre = invnorm1 * invnorm1;
        let cube = carre * invnorm1;
        // Derived from the norm
        let dprim = ncrossns1.dot(self.dnplan.cross(self.nsurf1));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = self.dnplan.cross(self.nsurf1).length_squared()
            + ncrossns1.dot(self.d2nplan.cross(self.nsurf1));
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.dnplan.dot(self.nsurf1);
        let p2 = self.d2nplan.dot(self.nsurf1);

        let temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
        // OCCT: d2n1w.SetLinearForm(smallterm * p1 + invnorm1 * p2, nplan,
        //                           smallterm * ndotns1 + 2 * invnorm1 * p1, dnplan,
        //                           ndotns1 * invnorm1, d2nplan);
        let mut d2n1w = (smallterm * p1 + invnorm1 * p2) * self.nplan
            + (smallterm * ndotns1 + 2.0 * invnorm1 * p1) * self.dnplan
            + (ndotns1 * invnorm1) * self.d2nplan;
        d2n1w += temp;

        // Derived from n2 compared to w
        let carre = invnorm2 * invnorm2;
        let cube = carre * invnorm2;
        // Derived from the norm
        let dprim = ncrossns2.dot(self.dnplan.cross(self.nsurf2));
        let smallterm = -2.0 * dprim * cube;
        let dsecn = self.dnplan.cross(self.nsurf2).length_squared()
            + ncrossns2.dot(self.d2nplan.cross(self.nsurf2));
        let grosterme = (3.0 * dprim * dprim * carre - dsecn) * cube;

        let p1 = self.dnplan.dot(self.nsurf2);
        let p2 = self.d2nplan.dot(self.nsurf2);

        let temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
        let mut d2n2w = (smallterm * p1 + invnorm2 * p2) * self.nplan
            + (smallterm * ndotns2 + 2.0 * invnorm2 * p1) * self.dnplan
            + (ndotns2 * invnorm2) * self.d2nplan;
        d2n2w += temp;

        let temp = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
        // OCCT: D2EDT2(1) = d2nplan.Dot(temp) - 2 * dnplan.Dot(d1gui) - nplan.Dot(d2gui);
        self.d2edt2[0] =
            self.d2nplan.dot(temp) - 2.0 * self.dnplan.dot(self.d1gui) - self.nplan.dot(self.d2gui);
        self.d2edt2[1] = self.ray1 * d2n1w.x - self.ray2 * d2n2w.x;
        self.d2edt2[2] = self.ray1 * d2n1w.y - self.ray2 * d2n2w.y;
        self.d2edt2[3] = self.ray1 * d2n1w.z - self.ray2 * d2n2w.z;
    }

    /// The third-order surface derivatives needed by the order-2 block
    /// (OCCT: surf1->D3 outputs d3u1, d3v1, d3uuv1, d3uvv1), taken with the
    /// generic partial-derivative accessor at the cached X.
    fn surf1_d3(&self) -> (DVec3, DVec3, DVec3, DVec3) {
        let (u, v) = (self.xval[0], self.xval[1]);
        (
            self.surf1.dn(u, v, 3, 0),
            self.surf1.dn(u, v, 0, 3),
            self.surf1.dn(u, v, 2, 1),
            self.surf1.dn(u, v, 1, 2),
        )
    }

    /// OCCT surf2->D3 outputs d3u2, d3v2, d3uuv2, d3uvv2.
    fn surf2_d3(&self) -> (DVec3, DVec3, DVec3, DVec3) {
        let (u, v) = (self.xval[2], self.xval[3]);
        (
            self.surf2.dn(u, v, 3, 0),
            self.surf2.dn(u, v, 0, 3),
            self.surf2.dn(u, v, 2, 1),
            self.surf2.dn(u, v, 1, 2),
        )
    }

    /// OCCT Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)
    /// (BlendFunc_ConstRad.cxx L1141-1199) — useful for a quick and
    /// approximate visualization of the surface area.
    #[allow(clippy::too_many_arguments)]
    pub fn section(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pdeb: &mut f64,
        pfin: &mut f64,
        c: &mut Circle3,
    ) {
        // OCCT: math_Vector X(1, 4);
        let x = [u1, v1, u2, v2];
        let prm = param;

        self.compute_values(&x, 0, true, prm);

        let mut ns1 = self.nsurf1;
        let mut np = self.nplan;

        let norm1 = self.nplan.cross(ns1).length();
        let norm1 = if norm1 < EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm1
        };
        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        // OCCT: Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());
        let center = self.pts1 + self.ray1 * ns1;

        // ns1 is oriented from the center to pts1,

        if self.ray1 > 0.0 {
            ns1 = -ns1;
        }
        if self.choix % 2 != 0 {
            np = -np;
        }
        // OCCT: C.SetRadius(std::abs(ray1)); C.SetPosition(gp_Ax2(Center, np, ns1));
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx.)
        c.radius = self.ray1.abs();
        c.center = center;
        c.normal = np;
        c.x_dir = ns1;
        c.y_dir = np.cross(ns1);
        *pdeb = 0.0;
        // OCCT: Pfin = ElCLib::Parameter(C, pts2);
        *pfin = super::brep_blend_func_consrad::elclib_circle_parameter(c, self.pts2);
        // Test negative and almost null angles : Singular Case
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns1);
            *pfin = super::brep_blend_func_consrad::elclib_circle_parameter(c, self.pts2);
        }
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT Section(P, Poles, Poles2d, Weights) (BlendFunc_ConstRad.cxx
    /// L1280-1351).
    pub fn section_simple(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        // OCCT: math_Vector X(1, 4);
        let mut x = [0.0f64; 4];
        let prm = p.parameter();

        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()

        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        x[0] = u1;
        x[1] = v1;
        x[2] = u2;
        x[3] = v2;

        self.compute_values(&x, 0, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avoid crushing the fields !
        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;

        poles_2d[0] = DVec2::new(x[0], x[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(x[2], x[3]);

        if self.my_s_shape == super::brep_blend_func::BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            return;
        }

        let mut norm1 = self.nplan.cross(ns1).length();
        let mut norm2 = self.nplan.cross(ns2).length();
        if norm1 < EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());
        let center = self.pts1 + self.ray1 * ns1;

        // ns1 (resp. ns2) is oriented from center to pts1 (resp. pts2),
        // and the triedron ns1,ns2,nplan is made direct.

        if self.ray1 > 0.0 {
            ns1 = -ns1;
        }
        if self.ray2 > 0.0 {
            ns2 = -ns2;
        }
        if self.choix % 2 != 0 {
            np = -np;
        }

        // OCCT L1350: GeomFill::GetCircle(myTConv, ns1, ns2, np, pts1, pts2,
        // std::abs(ray1), Center, Poles, Weights).
        let _ = (ns1, ns2, np, center);
        geomfill_get_circle_pending();
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BlendFunc_ConstRad.cxx L1355-1524) — used for the first and last
    /// section.
    #[allow(clippy::too_many_arguments)]
    #[allow(unreachable_code)] // the math_SVD / GetCircle fallbacks are pending
    pub fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        let mut ns1: DVec3;
        let mut ns2: DVec3;
        let mut np: DVec3;
        let mut dnp: DVec3;
        let mut dnorm1w: DVec3;
        let mut dnorm2w: DVec3;
        let tgc: DVec3;
        let norm1: f64;
        let norm2: f64;

        // OCCT: math_Vector sol(1, 4), secmember(1, 4);
        let mut sol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];

        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()
        let mut istgt = true;

        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        sol[0] = u1;
        sol[1] = v1;
        sol[2] = u2;
        sol[3] = v2;

        // Calculation of equations
        self.compute_values(&sol, 1, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avoid crushing the fields !
        ns1 = self.nsurf1;
        ns2 = self.nsurf2;
        np = self.nplan;
        dnp = self.dnplan;

        if !self.pts1.distance(self.pts2).le(&(1.0e-4)) {
            // OCCT L1387: if (!pts1.IsEqual(pts2, 1.e-4)) — the negation
            // drives the Gauss solve.
            // Calculation of derivates Processing Normal
            // OCCT: math_Gauss Resol(DEDX, 1.e-9);
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, self.dedx[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);

            if resol.is_done() {
                // OCCT: Resol.Solve(-DEDT, secmember);
                let mut sv = VecD::new(4);
                for i in 1..=4 {
                    sv.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut sv);
                for i in 1..=4 {
                    secmember[i - 1] = sv.get(i);
                }
                istgt = false;
            }
        }

        if istgt {
            // OCCT L1400-1408: math_SVD SingRS(DEDX); if (SingRS.IsDone()) {
            // SingRS.Solve(-DEDT, secmember, 1.e-6); istgt = false; }.
            // GAP (plan 0.6): math_SVD pending in rcad-kernel — see IsSolution.
        }

        if !istgt {
            // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
            self.tg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1;
            self.tg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2;

            // OCCT: dnorm1w.SetLinearForm(secmember(1), dndu1, secmember(2), dndv1, dn1w);
            dnorm1w = secmember[0] * self.dndu1 + secmember[1] * self.dndv1 + self.dn1w;
            dnorm2w = secmember[2] * self.dndu2 + secmember[3] * self.dndv2 + self.dn2w;
        } else {
            dnorm1w = DVec3::ZERO;
            dnorm2w = DVec3::ZERO;
        }

        // Tops 2d
        poles_2d[0] = DVec2::new(sol[0], sol[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
        if !istgt {
            d_poles_2d[0] = DVec2::new(secmember[0], secmember[1]);
            d_poles_2d[poles_2d.len() - 1] = DVec2::new(secmember[2], secmember[3]);
        }

        // the linear case is processed...
        if self.my_s_shape == super::brep_blend_func::BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            if !istgt {
                d_poles[low] = self.tg1;
                d_poles[upp] = self.tg2;
                d_weigths[low] = 0.0;
                d_weigths[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle
        norm1 = self.nplan.cross(ns1).length();
        norm2 = self.nplan.cross(ns2).length();
        let norm1 = if norm1 < EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm1
        };
        let norm2 = if norm2 < EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm2
        };

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());
        let center = self.pts1 + self.ray1 * ns1;
        if !istgt {
            // OCCT: tgc.SetLinearForm(ray1, dnorm1w, tg1);  = tg1.Added(ray1*dn1w);
            tgc = self.ray1 * dnorm1w + self.tg1;
        } else {
            tgc = DVec3::ZERO;
        }

        // ns1 is oriented from the center to pts1, and ns2 from the center to pts2
        // and the trihedron ns1,ns2,nplan is made direct

        if self.ray1 > 0.0 {
            ns1 = -ns1;
            if !istgt {
                dnorm1w = -dnorm1w;
            }
        }
        if self.ray2 > 0.0 {
            ns2 = -ns2;
            if !istgt {
                dnorm2w = -dnorm2w;
            }
        }
        if self.choix % 2 != 0 {
            np = -np;
            dnp = -dnp;
        }

        if !istgt {
            // OCCT L1499-1517: GeomFill::GetCircle(myTConv, ns1, ns2, dnorm1w,
            // dnorm2w, np, dnp, pts1, pts2, tg1, tg2, std::abs(ray1), 0,
            // Center, tgc, Poles, DPoles, Weights, DWeights).
            let _ = (ns1, ns2, dnorm1w, dnorm2w, np, dnp, center, tgc);
            geomfill_get_circle_pending()
        } else {
            // OCCT L1521: GeomFill::GetCircle(myTConv, ns1, ns2, np, pts1,
            // pts2, std::abs(ray1), Center, Poles, Weights).
            let _ = (ns1, ns2, np, center);
            geomfill_get_circle_pending()
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BlendFunc_ConstRad.cxx L1528-1903) —
    /// used for the first and last section.
    #[allow(clippy::too_many_arguments)]
    #[allow(unreachable_code)] // the math_SVD / GetCircle fallbacks are pending
    pub fn section_d2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        _d2_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        d2_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
        d2_weigths: &mut [f64],
    ) -> bool {
        let mut ns1: DVec3;
        let mut ns2: DVec3;
        let mut np: DVec3;
        let mut dnp: DVec3;
        let mut d2np: DVec3;
        let mut dnorm1w: DVec3;
        let mut dnorm2w: DVec3;
        let mut d2norm1w: DVec3;
        let mut d2norm2w: DVec3;
        let tgc: DVec3;
        let dtgc: DVec3;
        let dtg1: DVec3;
        let dtg2: DVec3;
        let norm1: f64;
        let norm2: f64;

        // OCCT: math_Vector X(1, 4), sol(1, 4), secmember(1, 4);
        //       math_Matrix D2DXdSdt(1, 4, 1, 4);
        let mut x = [0.0f64; 4];
        let mut sol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut d2dxdsdt = vec![vec![0.0f64; 4]; 4];

        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()
        let mut istgt = true;

        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        x[0] = u1;
        x[1] = v1;
        x[2] = u2;
        x[3] = v2;

        // Calculation of equations
        self.compute_values(&x, 2, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avois crushing the fields !
        ns1 = self.nsurf1;
        ns2 = self.nsurf2;
        np = self.nplan;
        dnp = self.dnplan;
        d2np = self.d2nplan;

        // Calculation of derivatives

        if self.pts1.distance(self.pts2) > 1.0e-4 {
            // OCCT L1706: if (!pts1.IsEqual(pts2, 1.e-4)).
            // OCCT: math_Gauss Resol(DEDX, 1.e-9); // Precise tolerance !!!!!
            // Calculation of derivatives Processing Normal
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, self.dedx[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                // OCCT: Resol.Solve(-DEDT, sol);
                let mut sv = VecD::new(4);
                for i in 1..=4 {
                    sv.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut sv);
                for i in 1..=4 {
                    sol[i - 1] = sv.get(i);
                }
                // OCCT: D2EDX2.Multiply(sol, D2DXdSdt);
                self.d2edx2.multiply(&sol, &mut d2dxdsdt);
                // OCCT: secmember = -(D2EDT2 + (2 * D2EDXDT + D2DXdSdt) * sol);
                for i in 0..4 {
                    let mut somme = 0.0;
                    for j in 0..4 {
                        somme += (2.0 * self.d2edxdt[i][j] + d2dxdsdt[i][j]) * sol[j];
                    }
                    secmember[i] = -(self.d2edt2[i] + somme);
                }
                // OCCT: Resol.Solve(secmember);
                let mut sv = VecD::new(4);
                for i in 1..=4 {
                    sv.set(i, secmember[i - 1]);
                }
                resol.solve(&mut sv);
                for i in 1..=4 {
                    secmember[i - 1] = sv.get(i);
                }
                istgt = false;
            }
        }

        if istgt {
            // OCCT L1720-1732: math_SVD SingRS(DEDX); math_Vector Vbis(1, 4);
            // if (SingRS.IsDone()) { SingRS.Solve(-DEDT, sol, 1.e-6);
            // D2EDX2.Multiply(sol, D2DXdSdt);
            // Vbis = -(D2EDT2 + (2 * D2EDXDT + D2DXdSdt) * sol);
            // SingRS.Solve(Vbis, secmember, 1.e-6); istgt = false; }.
            // GAP (plan 0.6): math_SVD pending in rcad-kernel — see IsSolution.
        }

        if !istgt {
            // OCCT: tg1.SetLinearForm(sol(1), d1u1, sol(2), d1v1);
            self.tg1 = sol[0] * self.d1u1 + sol[1] * self.d1v1;
            self.tg2 = sol[2] * self.d1u2 + sol[3] * self.d1v2;

            // OCCT: dnorm1w.SetLinearForm(sol(1), dndu1, sol(2), dndv1, dn1w);
            dnorm1w = sol[0] * self.dndu1 + sol[1] * self.dndv1 + self.dn1w;
            dnorm2w = sol[2] * self.dndu2 + sol[3] * self.dndv2 + self.dn2w;
            // OCCT: temp.SetLinearForm(sol(1) * sol(1), d2u1, 2 * sol(1) * sol(2), d2uv1,
            //                         sol(2) * sol(2), d2v1);
            let temp = (sol[0] * sol[0]) * self.d2u1
                + (2.0 * sol[0] * sol[1]) * self.d2uv1
                + (sol[1] * sol[1]) * self.d2v1;

            // OCCT: dtg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1, temp);
            dtg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1 + temp;

            let temp = (sol[2] * sol[2]) * self.d2u2
                + (2.0 * sol[2] * sol[3]) * self.d2uv2
                + (sol[3] * sol[3]) * self.d2v2;
            dtg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2 + temp;

            let mut temp = (sol[0] * sol[0]) * self.d2ndu1
                + (2.0 * sol[0] * sol[1]) * self.d2nduv1
                + (sol[1] * sol[1]) * self.d2ndv1;

            let tempbis = (2.0 * sol[0]) * self.d2ndtu1
                + (2.0 * sol[1]) * self.d2ndtv1
                + self.d2n1w;
            temp += tempbis;
            // OCCT: d2norm1w.SetLinearForm(secmember(1), dndu1, secmember(2), dndv1, temp);
            d2norm1w = secmember[0] * self.dndu1 + secmember[1] * self.dndv1 + temp;

            let mut temp = (sol[2] * sol[2]) * self.d2ndu2
                + (2.0 * sol[2] * sol[3]) * self.d2nduv2
                + (sol[3] * sol[3]) * self.d2ndv2;
            let tempbis = (2.0 * sol[2]) * self.d2ndtu2
                + (2.0 * sol[3]) * self.d2ndtv2
                + self.d2n2w;
            temp += tempbis;
            d2norm2w = secmember[2] * self.dndu2 + secmember[3] * self.dndv2 + temp;
        } else {
            dtg1 = DVec3::ZERO;
            dtg2 = DVec3::ZERO;
            dnorm1w = DVec3::ZERO;
            dnorm2w = DVec3::ZERO;
            d2norm1w = DVec3::ZERO;
            d2norm2w = DVec3::ZERO;
        }

        // Tops 2d
        poles_2d[0] = DVec2::new(x[0], x[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(x[2], x[3]);
        if !istgt {
            d_poles_2d[0] = DVec2::new(sol[0], sol[1]);
            d_poles_2d[poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
            d2_poles_2d[0] = DVec2::new(secmember[0], secmember[1]);
            d2_poles_2d[poles_2d.len() - 1] = DVec2::new(secmember[2], secmember[3]);
        }

        // linear case is processed...
        if self.my_s_shape == super::brep_blend_func::BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            if !istgt {
                // NOTE: the double assignment DPoles(low) = tg1; DPoles(low) = dtg1;
                // is translated literally from OCCT (L1790-1793).
                d_poles[low] = self.tg1;
                d_poles[upp] = self.tg2;
                d_poles[low] = dtg1;
                d_poles[upp] = dtg2;
                d_weigths[low] = 0.0;
                d_weigths[upp] = 0.0;
                d2_weigths[low] = 0.0;
                d2_weigths[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of circle
        norm1 = self.nplan.cross(ns1).length();
        norm2 = self.nplan.cross(ns2).length();
        let norm1 = if norm1 < EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm1
        };
        let norm2 = if norm2 < EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm2
        };

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + ray1 * ns1.XYZ());
        let center = self.pts1 + self.ray1 * ns1;
        if !istgt {
            // OCCT: tgc.SetLinearForm(ray1, dnorm1w, tg1);
            tgc = self.ray1 * dnorm1w + self.tg1;
            // OCCT: dtgc.SetLinearForm(ray1, d2norm1w, dtg1);
            dtgc = self.ray1 * d2norm1w + dtg1;
        } else {
            tgc = DVec3::ZERO;
            dtgc = DVec3::ZERO;
        }

        // ns1 is oriented from the center to pts1 and ns2 from the center to pts2
        // trihedron ns1,ns2,nplan is made direct

        if self.ray1 > 0.0 {
            ns1 = -ns1;
            if !istgt {
                dnorm1w = -dnorm1w;
                d2norm1w = -d2norm1w;
            }
        }
        if self.ray2 > 0.0 {
            ns2 = -ns2;
            if !istgt {
                dnorm2w = -dnorm2w;
                d2norm2w = -d2norm2w;
            }
        }
        if self.choix % 2 != 0 {
            np = -np;
            dnp = -dnp;
            d2np = -d2np;
        }

        if !istgt {
            // OCCT L1860-1887: GeomFill::GetCircle(myTConv, ns1, ns2, dnorm1w,
            // dnorm2w, d2norm1w, d2norm2w, np, dnp, d2np, pts1, pts2, tg1, tg2,
            // dtg1, dtg2, std::abs(ray1), 0, 0, Center, tgc, dtgc, Poles,
            // DPoles, D2Poles, Weights, DWeights, D2Weights).
            let _ = (ns1, ns2, dnorm1w, dnorm2w, d2norm1w, d2norm2w, np, dnp, d2np, center, tgc, dtgc);
            geomfill_get_circle_pending()
        } else {
            // OCCT L1891-1900: GeomFill::GetCircle(myTConv, ns1, ns2, nplan,
            // pts1, pts2, std::abs(ray1), Center, Poles, Weights).
            let _ = (ns1, ns2, self.nplan, center);
            geomfill_get_circle_pending()
        }
    }
}
