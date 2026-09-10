//! OCCT BlendFunc_EvolRad — part 2: ComputeValues (BlendFunc_EvolRad.cxx
//! L194-843), Tangent (L1058-1109) and the four Section overloads
//! (L1135-1193, L1376-1454, L1458-1634, L1638-2019).  The struct and the
//! remaining methods live in [`super::brep_blend_func_evolrad`].

use glam::{DVec2, DVec3};

use rcad_kernel::base::convert::ConvertParameterisation;
use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::geom::{Circle3, CurveEval as _, SurfaceEval as _};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{MatD, VecD};

use crate::geomalgo::geomfill::geom_fill::{get_circle, get_circle_d1, get_circle_d2};

use super::brep_blend_func::{blend_func_compute_dnormal, blend_func_compute_normal};
use super::brep_blend_func_consrad::elclib_circle_parameter;
use super::brep_blend_func_evolrad::BlendFuncEvolRad;
use super::brep_blend_point::BlendPoint;

/// Architecture mapping: OCCT consumes GeomFill's
/// `Convert_ParameterisationType` through both BlendFunc and GeomFill; the
/// rcad kernel exposes the same enumeration as
/// [`ConvertParameterisation`].  This is the identity mapping between the
/// two rcad spellings of the OCCT enum.
fn tconv(t_conv: super::brep_blend_func::ConvertParameterisationType) -> ConvertParameterisation {
    match t_conv {
        super::brep_blend_func::ConvertParameterisationType::TgtThetaOver2 => {
            ConvertParameterisation::TgtThetaOver2
        }
        super::brep_blend_func::ConvertParameterisationType::TgtThetaOver2_1 => {
            ConvertParameterisation::TgtThetaOver2_1
        }
        super::brep_blend_func::ConvertParameterisationType::TgtThetaOver2_2 => {
            ConvertParameterisation::TgtThetaOver2_2
        }
        super::brep_blend_func::ConvertParameterisationType::TgtThetaOver2_3 => {
            ConvertParameterisation::TgtThetaOver2_3
        }
        super::brep_blend_func::ConvertParameterisationType::TgtThetaOver2_4 => {
            ConvertParameterisation::TgtThetaOver2_4
        }
        super::brep_blend_func::ConvertParameterisationType::QuasiAngular => {
            ConvertParameterisation::QuasiAngular
        }
        super::brep_blend_func::ConvertParameterisationType::RationalC1 => {
            ConvertParameterisation::RationalC1
        }
        super::brep_blend_func::ConvertParameterisationType::Polynomial => {
            ConvertParameterisation::Polynomial
        }
    }
}

impl<'a> BlendFuncEvolRad<'a> {
    /// OCCT ComputeValues(X, Order, ByParam, Param)
    /// (BlendFunc_EvolRad.cxx L194-843) — OBLIGATORY passage for all
    /// computations: positions on Surfaces and Curve and the law, partial
    /// calculation of the equations and their derivatives, storage of
    /// intermediate results in fields.
    pub(crate) fn compute_values(
        &mut self,
        x: &[f64],
        order: i32,
        by_param: bool,
        param: f64,
    ) -> bool {
        // OCCT L201-205: the `static` locals (d3u1... d3gui, ptgui,
        // invnormtg, dinvnormtg) are a reallocation optimisation; the subset
        // read across calls on cached states (ptgui, d1gui, d2gui, invnormtg,
        // dinvnormtg) lives on the struct, the rest are plain locals below.
        let mut t = param;

        // Case of implicit parameter
        if !by_param {
            t = self.param;
        }

        // The work is done already?
        let mut l_x_ok = order <= self.my_x_order;
        let mut ii = 1usize;
        while ii <= x.len() && l_x_ok {
            l_x_ok = x[ii - 1] == self.xval[ii - 1];
            ii += 1;
        }

        let t_ok = (t == self.tval) && ((order <= self.my_t_order) || (!by_param));

        if l_x_ok && t_ok {
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
            //----- Positioning on the curve and the law----------------
            match self.my_t_order {
                0 => {
                    // OCCT: tcurv->D1(T, ptgui, d1gui);
                    let (p, d1) = {
                        let c = self.tcurv();
                        (c.point_at(t), c.derivative_at(t))
                    };
                    self.ptgui = p;
                    self.d1gui = d1;
                    // OCCT: nplan = d1gui.Normalized();
                    self.nplan = self.d1gui.normalize();
                    // OCCT: ray = tevol->Value(T);
                    let law = self.tevol().clone();
                    self.ray = law.borrow_mut().value(t);
                }

                1 => {
                    // OCCT: tcurv->D2(T, ptgui, d1gui, d2gui);
                    let (p, d1, d2) = {
                        let c = self.tcurv();
                        (c.point_at(t), c.derivative_at(t), c.derivative2_at(t))
                    };
                    self.ptgui = p;
                    self.d1gui = d1;
                    self.d2gui = d2;
                    let d2gui = self.d2gui;
                    // OCCT: nplan = d1gui.Normalized();
                    self.nplan = self.d1gui.normalize();
                    // OCCT: invnormtg = ((double)1) / d1gui.Magnitude();
                    self.invnormtg = 1.0 / self.d1gui.length();
                    // OCCT: dnplan.SetLinearForm(invnormtg, d2gui,
                    //           -invnormtg * (nplan.Dot(d2gui)), nplan);
                    self.dnplan = self.invnormtg * d2gui
                        + (-self.invnormtg * self.nplan.dot(d2gui)) * self.nplan;

                    // OCCT: tevol->D1(T, ray, dray);
                    let law = self.tevol().clone();
                    let mut f = 0.0;
                    let mut d = 0.0;
                    law.borrow_mut().d1(t, &mut f, &mut d);
                    self.ray = f;
                    self.dray = d;
                }
                2 => {
                    // OCCT: tcurv->D3(T, ptgui, d1gui, d2gui, d3gui);
                    // (d3gui is written by OCCT but never read — the write is
                    // kept as a discarded local for statement parity.)
                    let (p, d1, d2, d3) = {
                        let c = self.tcurv();
                        (
                            c.point_at(t),
                            c.derivative_at(t),
                            c.derivative2_at(t),
                            c.derivative3_at(t),
                        )
                    };
                    self.ptgui = p;
                    self.d1gui = d1;
                    self.d2gui = d2;
                    let d3gui = d3;
                    let d2gui = self.d2gui;
                    // OCCT: nplan = d1gui.Normalized();
                    self.nplan = self.d1gui.normalize();
                    self.invnormtg = 1.0 / self.d1gui.length();
                    // OCCT: dnplan.SetLinearForm(invnormtg, d2gui,
                    //           -invnormtg * (nplan.Dot(d2gui)), nplan);
                    self.dnplan = self.invnormtg * d2gui
                        + (-self.invnormtg * self.nplan.dot(d2gui)) * self.nplan;
                    // OCCT: dinvnormtg = -nplan.Dot(d2gui) * invnormtg * invnormtg;
                    self.dinvnormtg =
                        -self.nplan.dot(d2gui) * self.invnormtg * self.invnormtg;
                    // OCCT: d2nplan.SetLinearForm(invnormtg, d3gui, dinvnormtg, d2gui);
                    let d2nplan = self.invnormtg * d3gui + self.dinvnormtg * d2gui;
                    // OCCT: aux = dinvnormtg * (nplan.Dot(d2gui))
                    //           + invnormtg * (dnplan.Dot(d2gui) + nplan.Dot(d3gui));
                    let aux = self.dinvnormtg * self.nplan.dot(d2gui)
                        + self.invnormtg
                            * (self.dnplan.dot(d2gui) + self.nplan.dot(d3gui));
                    // OCCT: d2nplan.SetLinearForm(-invnormtg * (nplan.Dot(d2gui)),
                    //           dnplan, -aux, nplan, d2nplan);
                    self.d2nplan = (-self.invnormtg * self.nplan.dot(d2gui)) * self.dnplan
                        + (-aux) * self.nplan
                        + d2nplan;

                    // OCCT: tevol->D2(T, ray, dray, d2ray);
                    let law = self.tevol().clone();
                    let mut f = 0.0;
                    let mut d = 0.0;
                    let mut d2 = 0.0;
                    law.borrow_mut().d2(t, &mut f, &mut d, &mut d2);
                    self.ray = f;
                    self.dray = d;
                    self.d2ray = d2;
                }
                _ => return false,
            }
        }

        // Processing of X
        if !l_x_ok {
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
                    // OCCT: nsurf1 = d1u1.Crossed(d1v1);
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
                    self.d2v1 = d2v;
                    self.d2uv1 = d2uv;
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
                    self.d2v2 = d2v;
                    self.d2uv2 = d2uv;
                    self.nsurf2 = self.d1u2.cross(self.d1v2);
                    self.dns1u2 = self.d2u2.cross(self.d1v2) + self.d1u2.cross(self.d2uv2);
                    self.dns1v2 = self.d2uv2.cross(self.d1v2) + self.d1u2.cross(self.d2v2);
                }
                2 => {
                    // OCCT: surf1->D3(X(1), X(2), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1,
                    //                  d3u1, d3v1, d3uuv1, d3uvv1);
                    // (the D3 outputs are read later in the order-2 block via
                    // the xval-cached surf1_d3()/surf2_d3() accessors.)
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf1.derivatives2(x[0], x[1]);
                    self.pts1 = p;
                    self.d1u1 = du;
                    self.d1v1 = dv;
                    self.d2u1 = d2u;
                    self.d2v1 = d2v;
                    self.d2uv1 = d2uv;
                    self.nsurf1 = self.d1u1.cross(self.d1v1);
                    // OCCT: surf2->D3(X(3), X(4), pts2, d1u2, d1v2, d2u2, d2v2, d2uv2,
                    //                  d3u2, d3v2, d3uuv2, d3uvv2);
                    let (p, du, dv, d2u, d2uv, d2v) = self.surf2.derivatives2(x[2], x[3]);
                    self.pts2 = p;
                    self.d1u2 = du;
                    self.d1v2 = dv;
                    self.d2u2 = d2u;
                    self.d2v2 = d2v;
                    self.d2uv2 = d2uv;
                    self.nsurf2 = self.d1u2.cross(self.d1v2);
                }
                _ => return false,
            }
            // Case of degenerated surfaces
            if self.nsurf1.length() < super::brep_blend_func_evolrad::EPS {
                // OCCT: gp_Pnt2d P(X(1), X(2));
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
            if self.nsurf2.length() < super::brep_blend_func_evolrad::EPS {
                // OCCT: gp_Pnt2d P(X(3), X(4));
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
        let ray1 = self.sg1 * self.ray;
        let ray2 = self.sg2 * self.ray;

        // OCCT: theD = -(nplan.XYZ().Dot(ptgui.XYZ()));
        let the_d = -self.nplan.dot(self.ptgui);

        // OCCT: E(1) = (nplan.X() * (pts1.X() + pts2.X()) + ... ) / 2 + theD;
        self.e[0] = (self.nplan.x * (self.pts1.x + self.pts2.x)
            + self.nplan.y * (self.pts1.y + self.pts2.y)
            + self.nplan.z * (self.pts1.z + self.pts2.z))
            / 2.0
            + the_d;

        let ncrossns1 = self.nplan.cross(self.nsurf1);
        let ncrossns2 = self.nplan.cross(self.nsurf2);
        let mut invnorm1 = ncrossns1.length();
        let mut invnorm2 = ncrossns2.length();

        invnorm1 = if invnorm1 > super::brep_blend_func_evolrad::EPS {
            1.0 / invnorm1
        } else {
            1.0 // Unsatisfactory, but it is not necessary to stop
        };
        invnorm2 = if invnorm2 > super::brep_blend_func_evolrad::EPS {
            1.0 / invnorm2
        } else {
            1.0 // Unsatisfactory, but it is not necessary to stop
        };

        let ndotns1 = self.nplan.dot(self.nsurf1);
        let ndotns2 = self.nplan.dot(self.nsurf2);

        // OCCT: n1.SetLinearForm(ndotns1, nplan, -1., nsurf1); n1.Multiply(invnorm1);
        let n1 = (ndotns1 * self.nplan - self.nsurf1) * invnorm1;
        let n2 = (ndotns2 * self.nplan - self.nsurf2) * invnorm2;

        // OCCT: resul.SetLinearForm(ray1, n1, -ray2, n2, gp_Vec(pts2, pts1));
        let resul = ray1 * n1 - ray2 * n2 + (self.pts1 - self.pts2);

        self.e[1] = resul.x;
        self.e[2] = resul.y;
        self.e[3] = resul.z;

        // -------------------- Positioning of order 1 ---------------------
        if order >= 1 {
            // OCCT: DEDX(1, 1..4) = nplan.Dot(d1u1|d1v1|d1u2|d1v2) / 2
            self.dedx[0][0] = self.nplan.dot(self.d1u1) / 2.0;
            self.dedx[0][1] = self.nplan.dot(self.d1v1) / 2.0;
            self.dedx[0][2] = self.nplan.dot(self.d1u2) / 2.0;
            self.dedx[0][3] = self.nplan.dot(self.d1v2) / 2.0;

            let mut cube = invnorm1 * invnorm1 * invnorm1;
            // Derived compared to u1
            // OCCT: grosterme = -ncrossns1.Dot(nplan.Crossed(dns1u1)) * cube;
            let mut grosterme = -ncrossns1.dot(self.nplan.cross(self.dns1u1)) * cube;
            // OCCT: dndu1.SetLinearForm(grosterme * ndotns1 + invnorm1 * nplan.Dot(dns1u1),
            //           nplan, -grosterme, nsurf1, -invnorm1, dns1u1);
            self.dndu1 = (grosterme * ndotns1 + invnorm1 * self.nplan.dot(self.dns1u1))
                * self.nplan
                + (-grosterme) * self.nsurf1
                + (-invnorm1) * self.dns1u1;

            // OCCT: resul.SetLinearForm(ray1, dndu1, d1u1);
            let mut resul = ray1 * self.dndu1 + self.d1u1;
            self.dedx[1][0] = resul.x;
            self.dedx[2][0] = resul.y;
            self.dedx[3][0] = resul.z;

            // Derived compared to v1
            grosterme = -ncrossns1.dot(self.nplan.cross(self.dns1v1)) * cube;
            self.dndv1 = (grosterme * ndotns1 + invnorm1 * self.nplan.dot(self.dns1v1))
                * self.nplan
                + (-grosterme) * self.nsurf1
                + (-invnorm1) * self.dns1v1;

            resul = ray1 * self.dndv1 + self.d1v1;
            self.dedx[1][1] = resul.x;
            self.dedx[2][1] = resul.y;
            self.dedx[3][1] = resul.z;

            cube = invnorm2 * invnorm2 * invnorm2;
            // Derivee par rapport a u2
            grosterme = -ncrossns2.dot(self.nplan.cross(self.dns1u2)) * cube;
            self.dndu2 = (grosterme * ndotns2 + invnorm2 * self.nplan.dot(self.dns1u2))
                * self.nplan
                + (-grosterme) * self.nsurf2
                + (-invnorm2) * self.dns1u2;

            // OCCT: resul.SetLinearForm(-ray2, dndu2, -1, d1u2);
            resul = -ray2 * self.dndu2 - self.d1u2;
            self.dedx[1][2] = resul.x;
            self.dedx[2][2] = resul.y;
            self.dedx[3][2] = resul.z;

            // Derived compared to v2
            grosterme = -ncrossns2.dot(self.nplan.cross(self.dns1v2)) * cube;
            self.dndv2 = (grosterme * ndotns2 + invnorm2 * self.nplan.dot(self.dns1v2))
                * self.nplan
                + (-grosterme) * self.nsurf2
                + (-invnorm2) * self.dns1v2;

            resul = -ray2 * self.dndv2 - self.d1v2;
            self.dedx[1][3] = resul.x;
            self.dedx[2][3] = resul.y;
            self.dedx[3][3] = resul.z;

            if by_param {
                // OCCT: temp.SetXYZ((pts1.XYZ() + pts2.XYZ()) / 2 - ptgui.XYZ());
                let mut temp = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
                // Derived from n1 compared to w
                // OCCT: grosterme = ncrossns1.Dot(dnplan.Crossed(nsurf1)) * invnorm1 * invnorm1;
                let grosterme = ncrossns1.dot(self.dnplan.cross(self.nsurf1)) * invnorm1 * invnorm1;
                // OCCT: dn1w.SetLinearForm((dnplan.Dot(nsurf1) - grosterme * ndotns1) * invnorm1,
                //           nplan, ndotns1 * invnorm1, dnplan, grosterme * invnorm1, nsurf1);
                self.dn1w = ((self.dnplan.dot(self.nsurf1) - grosterme * ndotns1) * invnorm1)
                    * self.nplan
                    + (ndotns1 * invnorm1) * self.dnplan
                    + (grosterme * invnorm1) * self.nsurf1;

                // Derived from n2 compared to w
                let grosterme = ncrossns2.dot(self.dnplan.cross(self.nsurf2)) * invnorm2 * invnorm2;
                self.dn2w = ((self.dnplan.dot(self.nsurf2) - grosterme * ndotns2) * invnorm2)
                    * self.nplan
                    + (ndotns2 * invnorm2) * self.dnplan
                    + (grosterme * invnorm2) * self.nsurf2;

                // OCCT: DEDT(1) = dnplan.Dot(temp) - 1. / invnormtg;
                self.dedt[0] = self.dnplan.dot(temp) - 1.0 / self.invnormtg;
                // OCCT: temp.SetLinearForm(ray2, dn2w, sg2 * dray, n2);
                temp = ray2 * self.dn2w + (self.sg2 * self.dray) * n2;
                // OCCT: resul.SetLinearForm(ray1, dn1w, sg1 * dray, n1, -1, temp);
                let resul = ray1 * self.dn1w + (self.sg1 * self.dray) * n1 - temp;
                self.dedt[1] = resul.x;
                self.dedt[2] = resul.y;
                self.dedt[3] = resul.z;
            }
            // ------   Positioning of order 2  -----------------------------
            if order == 2 {
                // OCCT L499: the d3 surface outputs of the order-2 block
                // (statics d3u1, d3v1, d3uuv1, d3uvv1 / d3u2...), read from
                // the xval cache exactly as the OCCT statics.
                let (d3u1, d3v1, d3uuv1, d3uvv1) = self.surf1_d3();
                let (d3u2, d3v2, d3uuv2, d3uvv2) = self.surf2_d3();
                let mut uterm: f64;
                let mut vterm: f64;
                let mut smallterm: f64;
                let mut p1: f64;
                let mut p2: f64;
                let mut p12: f64;
                let mut d_prim: f64;
                let mut d_secn: f64;
                // OCCT: D2EDX2.Init(0);
                self.d2edx2.init(0.0);

                // OCCT: D2EDX2(1, 1, 1) = nplan.Dot(d2u1) / 2;
                *self.d2edx2.change_value(1, 1, 1) = self.nplan.dot(self.d2u1) / 2.0;
                // OCCT: D2EDX2(1, 2, 1) = D2EDX2(1, 1, 2) = nplan.Dot(d2uv1) / 2;
                *self.d2edx2.change_value(1, 2, 1) = self.nplan.dot(self.d2uv1) / 2.0;
                *self.d2edx2.change_value(1, 1, 2) = self.nplan.dot(self.d2uv1) / 2.0;
                // OCCT: D2EDX2(1, 2, 2) = nplan.Dot(d2v1) / 2;
                *self.d2edx2.change_value(1, 2, 2) = self.nplan.dot(self.d2v1) / 2.0;

                // OCCT: D2EDX2(1, 3, 3) = nplan.Dot(d2u2) / 2;
                *self.d2edx2.change_value(1, 3, 3) = self.nplan.dot(self.d2u2) / 2.0;
                // OCCT: D2EDX2(1, 4, 3) = D2EDX2(1, 3, 4) = nplan.Dot(d2uv2) / 2;
                *self.d2edx2.change_value(1, 4, 3) = self.nplan.dot(self.d2uv2) / 2.0;
                *self.d2edx2.change_value(1, 3, 4) = self.nplan.dot(self.d2uv2) / 2.0;
                // OCCT: D2EDX2(1, 4, 4) = nplan.Dot(d2v2) / 2;
                *self.d2edx2.change_value(1, 4, 4) = self.nplan.dot(self.d2v2) / 2.0;
                // ================
                // ==  Surface 1 ==
                // ================
                let mut carre = invnorm1 * invnorm1;
                let mut cube = carre * invnorm1;
                // Derived double compared to u1
                // Derived from the norm
                // OCCT: d2ns1u1.SetLinearForm(1, d3u1.Crossed(d1v1),
                //           2, d2u1.Crossed(d2uv1), 1, d1u1.Crossed(d3uuv1));
                let d2ns1u1 = d3u1.cross(self.d1v1)
                    + 2.0 * (self.d2u1.cross(self.d2uv1))
                    + self.d1u1.cross(d3uuv1);
                d_prim = ncrossns1.dot(self.nplan.cross(self.dns1u1));
                smallterm = -2.0 * d_prim * cube;
                // OCCT: DSecn = ncrossns1.Dot(nplan.Crossed(d2ns1u1))
                //           + (nplan.Crossed(dns1u1)).SquareMagnitude();
                d_secn = ncrossns1.dot(self.nplan.cross(d2ns1u1))
                    + self.nplan.cross(self.dns1u1).dot(self.nplan.cross(self.dns1u1));
                // OCCT: grosterme = (3 * DPrim * DPrim * carre - DSecn) * cube;
                let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                // OCCT: temp.SetLinearForm(grosterme * ndotns1, nplan, -grosterme, nsurf1);
                let mut temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
                p1 = self.nplan.dot(self.dns1u1);
                p2 = self.nplan.dot(d2ns1u1);
                // OCCT: d2ndu1.SetLinearForm(invnorm1 * p2 + smallterm * p1,
                //           nplan, -smallterm, dns1u1, -invnorm1, d2ns1u1);
                self.d2ndu1 = (invnorm1 * p2 + smallterm * p1) * self.nplan
                    + (-smallterm) * self.dns1u1
                    + (-invnorm1) * d2ns1u1;
                // OCCT: d2ndu1 += temp;
                self.d2ndu1 += temp;
                // OCCT: resul.SetLinearForm(ray1, d2ndu1, d2u1);
                let mut resul = ray1 * self.d2ndu1 + self.d2u1;
                *self.d2edx2.change_value(2, 1, 1) = resul.x;
                *self.d2edx2.change_value(3, 1, 1) = resul.y;
                *self.d2edx2.change_value(4, 1, 1) = resul.z;

                // Derived double compared to u1, v1
                // Derived from the norm
                // OCCT: d2ns1uv1 = (d3uuv1.Crossed(d1v1)) + (d2u1.Crossed(d2v1))
                //           + (d1u1.Crossed(d3uvv1));
                let d2ns1uv1 =
                    d3uuv1.cross(self.d1v1) + self.d2u1.cross(self.d2v1) + self.d1u1.cross(d3uvv1);
                uterm = ncrossns1.dot(self.nplan.cross(self.dns1u1));
                vterm = ncrossns1.dot(self.nplan.cross(self.dns1v1));
                // OCCT: DSecn = (nplan.Crossed(dns1v1)).Dot(nplan.Crossed(dns1u1))
                //           + ncrossns1.Dot(nplan.Crossed(d2ns1uv1));
                d_secn = self.nplan.cross(self.dns1v1).dot(self.nplan.cross(self.dns1u1))
                    + ncrossns1.dot(self.nplan.cross(d2ns1uv1));
                // OCCT: grosterme = (3 * uterm * vterm * carre - DSecn) * cube;
                let grosterme = (3.0 * uterm * vterm * carre - d_secn) * cube;
                // OCCT: uterm *= -cube; // and only now
                uterm *= -cube;
                // OCCT: vterm *= -cube;
                vterm *= -cube;

                p1 = self.nplan.dot(self.dns1u1);
                p2 = self.nplan.dot(self.dns1v1);
                // OCCT: temp.SetLinearForm(grosterme * ndotns1, nplan, -grosterme, nsurf1,
                //           -invnorm1, d2ns1uv1);
                temp = (grosterme * ndotns1) * self.nplan
                    + (-grosterme) * self.nsurf1
                    + (-invnorm1) * d2ns1uv1;
                // OCCT: d2nduv1.SetLinearForm(invnorm1 * nplan.Dot(d2ns1uv1) + uterm * p2
                //           + vterm * p1, nplan, -uterm, dns1v1, -vterm, dns1u1);
                self.d2nduv1 = (invnorm1 * self.nplan.dot(d2ns1uv1) + uterm * p2 + vterm * p1)
                    * self.nplan
                    + (-uterm) * self.dns1v1
                    + (-vterm) * self.dns1u1;

                // OCCT: d2nduv1 += temp;
                self.d2nduv1 += temp;
                // OCCT: resul.SetLinearForm(ray1, d2nduv1, d2uv1);
                resul = ray1 * self.d2nduv1 + self.d2uv1;

                // OCCT: D2EDX2(2, 2, 1) = D2EDX2(2, 1, 2) = resul.X(); ...
                *self.d2edx2.change_value(2, 2, 1) = resul.x;
                *self.d2edx2.change_value(2, 1, 2) = resul.x;
                *self.d2edx2.change_value(3, 2, 1) = resul.y;
                *self.d2edx2.change_value(3, 1, 2) = resul.y;
                *self.d2edx2.change_value(4, 2, 1) = resul.z;
                *self.d2edx2.change_value(4, 1, 2) = resul.z;

                // Derived double compared to v1
                // Derived from the norm
                // OCCT: d2ns1v1.SetLinearForm(1, d1u1.Crossed(d3v1),
                //           2, d2uv1.Crossed(d2v1), 1, d3uvv1.Crossed(d1v1));
                let d2ns1v1 = self.d1u1.cross(d3v1)
                    + 2.0 * (self.d2uv1.cross(self.d2v1))
                    + d3uvv1.cross(self.d1v1);
                d_prim = ncrossns1.dot(self.nplan.cross(self.dns1v1));
                smallterm = -2.0 * d_prim * cube;
                d_secn = ncrossns1.dot(self.nplan.cross(d2ns1v1))
                    + self.nplan.cross(self.dns1v1).dot(self.nplan.cross(self.dns1v1));
                let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                p1 = self.nplan.dot(self.dns1v1);
                p2 = self.nplan.dot(d2ns1v1);
                temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
                self.d2ndv1 = (invnorm1 * p2 + smallterm * p1) * self.nplan
                    + (-smallterm) * self.dns1v1
                    + (-invnorm1) * d2ns1v1;
                self.d2ndv1 += temp;
                resul = ray1 * self.d2ndv1 + self.d2v1;

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
                let d2ns1u2 = d3u2.cross(self.d1v2)
                    + 2.0 * (self.d2u2.cross(self.d2uv2))
                    + self.d1u2.cross(d3uuv2);
                d_prim = ncrossns2.dot(self.nplan.cross(self.dns1u2));
                smallterm = -2.0 * d_prim * cube;
                d_secn = ncrossns2.dot(self.nplan.cross(d2ns1u2))
                    + self.nplan.cross(self.dns1u2).dot(self.nplan.cross(self.dns1u2));
                let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
                p1 = self.nplan.dot(self.dns1u2);
                p2 = self.nplan.dot(d2ns1u2);
                self.d2ndu2 = (invnorm2 * p2 + smallterm * p1) * self.nplan
                    + (-smallterm) * self.dns1u2
                    + (-invnorm2) * d2ns1u2;
                self.d2ndu2 += temp;
                // OCCT: resul.SetLinearForm(-ray2, d2ndu2, -1, d2u2);
                resul = -ray2 * self.d2ndu2 - self.d2u2;
                *self.d2edx2.change_value(2, 3, 3) = resul.x;
                *self.d2edx2.change_value(3, 3, 3) = resul.y;
                *self.d2edx2.change_value(4, 3, 3) = resul.z;

                // Derived double compared to u2, v2
                // Derived from the norm
                let d2ns1uv2 =
                    d3uuv2.cross(self.d1v2) + self.d2u2.cross(self.d2v2) + self.d1u2.cross(d3uvv2);
                uterm = ncrossns2.dot(self.nplan.cross(self.dns1u2));
                vterm = ncrossns2.dot(self.nplan.cross(self.dns1v2));
                d_secn = self.nplan.cross(self.dns1v2).dot(self.nplan.cross(self.dns1u2))
                    + ncrossns2.dot(self.nplan.cross(d2ns1uv2));
                let grosterme = (3.0 * uterm * vterm * carre - d_secn) * cube;
                uterm *= -cube; // and only now
                vterm *= -cube;

                p1 = self.nplan.dot(self.dns1u2);
                p2 = self.nplan.dot(self.dns1v2);
                temp = (grosterme * ndotns2) * self.nplan
                    + (-grosterme) * self.nsurf2
                    + (-invnorm2) * d2ns1uv2;
                self.d2nduv2 = (invnorm2 * self.nplan.dot(d2ns1uv2) + uterm * p2 + vterm * p1)
                    * self.nplan
                    + (-uterm) * self.dns1v2
                    + (-vterm) * self.dns1u2;

                self.d2nduv2 += temp;
                resul = -ray2 * self.d2nduv2 - self.d2uv2;

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
                d_prim = ncrossns2.dot(self.nplan.cross(self.dns1v2));
                smallterm = -2.0 * d_prim * cube;
                d_secn = ncrossns2.dot(self.nplan.cross(d2ns1v2))
                    + self.nplan.cross(self.dns1v2).dot(self.nplan.cross(self.dns1v2));
                let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                p1 = self.nplan.dot(self.dns1v2);
                p2 = self.nplan.dot(d2ns1v2);
                temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
                self.d2ndv2 = (invnorm2 * p2 + smallterm * p1) * self.nplan
                    + (-smallterm) * self.dns1v2
                    + (-invnorm2) * d2ns1v2;
                self.d2ndv2 += temp;
                resul = -ray2 * self.d2ndv2 - self.d2v2;

                *self.d2edx2.change_value(2, 4, 4) = resul.x;
                *self.d2edx2.change_value(3, 4, 4) = resul.y;
                *self.d2edx2.change_value(4, 4, 4) = resul.z;

                if by_param {
                    let mut tterm: f64;
                    //  ---------- Double Derivation on t, X --------------------------
                    // OCCT: D2EDXDT(1, 1..4) = dnplan.Dot(d1u1|d1v1|d1u2|d1v2) / 2
                    self.d2edxdt[0][0] = self.dnplan.dot(self.d1u1) / 2.0;
                    self.d2edxdt[0][1] = self.dnplan.dot(self.d1v1) / 2.0;
                    self.d2edxdt[0][2] = self.dnplan.dot(self.d1u2) / 2.0;
                    self.d2edxdt[0][3] = self.dnplan.dot(self.d1v2) / 2.0;

                    carre = invnorm1 * invnorm1;
                    cube = carre * invnorm1;
                    //--> Derived compared to u1 and t
                    tterm = ncrossns1.dot(self.dnplan.cross(self.nsurf1));
                    smallterm = -tterm * cube;
                    // Derived from the norm
                    uterm = ncrossns1.dot(self.nplan.cross(self.dns1u1));
                    d_secn = self.nplan.cross(self.dns1u1).dot(self.dnplan.cross(self.nsurf1))
                        + ncrossns1.dot(self.dnplan.cross(self.dns1u1));
                    let mut grosterme = (3.0 * uterm * tterm * carre - d_secn) * cube;
                    uterm *= -cube;

                    p1 = self.dnplan.dot(self.nsurf1);
                    p2 = self.nplan.dot(self.dns1u1);
                    p12 = self.dnplan.dot(self.dns1u1);

                    // OCCT: d2ndtu1.SetLinearForm(invnorm1 * p12 + smallterm * p2
                    //           + uterm * p1 + grosterme * ndotns1, nplan,
                    //           invnorm1 * p2 + uterm * ndotns1, dnplan,
                    //           -smallterm, dns1u1);
                    self.d2ndtu1 = (invnorm1 * p12 + smallterm * p2 + uterm * p1
                        + grosterme * ndotns1)
                        * self.nplan
                        + (invnorm1 * p2 + uterm * ndotns1) * self.dnplan
                        + (-smallterm) * self.dns1u1;
                    // OCCT: d2ndtu1 -= grosterme * nsurf1;
                    self.d2ndtu1 -= grosterme * self.nsurf1;

                    // OCCT: resul.SetLinearForm(ray1, d2ndtu1, sg1 * dray, dndu1);
                    resul = ray1 * self.d2ndtu1 + (self.sg1 * self.dray) * self.dndu1;
                    self.d2edxdt[1][0] = resul.x;
                    self.d2edxdt[2][0] = resul.y;
                    self.d2edxdt[3][0] = resul.z;

                    //--> Derived compared to v1 and t
                    // Derived from the norm
                    uterm = ncrossns1.dot(self.nplan.cross(self.dns1v1));
                    d_secn = self.nplan.cross(self.dns1v1).dot(self.dnplan.cross(self.nsurf1))
                        + ncrossns1.dot(self.dnplan.cross(self.dns1v1));
                    grosterme = (3.0 * uterm * tterm * carre - d_secn) * cube;
                    uterm *= -cube;

                    p1 = self.dnplan.dot(self.nsurf1);
                    p2 = self.nplan.dot(self.dns1v1);
                    p12 = self.dnplan.dot(self.dns1v1);
                    self.d2ndtv1 = (invnorm1 * p12 + uterm * p1 + smallterm * p2
                        + grosterme * ndotns1)
                        * self.nplan
                        + (invnorm1 * p2 + uterm * ndotns1) * self.dnplan
                        + (-smallterm) * self.dns1v1;
                    self.d2ndtv1 -= grosterme * self.nsurf1;

                    resul = ray1 * self.d2ndtv1 + (self.sg1 * self.dray) * self.dndv1;
                    self.d2edxdt[1][1] = resul.x;
                    self.d2edxdt[2][1] = resul.y;
                    self.d2edxdt[3][1] = resul.z;

                    carre = invnorm2 * invnorm2;
                    cube = carre * invnorm2;
                    //--> Derived compared to u2 and t
                    tterm = ncrossns2.dot(self.dnplan.cross(self.nsurf2));
                    smallterm = -tterm * cube;
                    // Derived from the norm
                    uterm = ncrossns2.dot(self.nplan.cross(self.dns1u2));
                    d_secn = self.nplan.cross(self.dns1u2).dot(self.dnplan.cross(self.nsurf2))
                        + ncrossns2.dot(self.dnplan.cross(self.dns1u2));
                    grosterme = (3.0 * uterm * tterm * carre - d_secn) * cube;
                    uterm *= -cube;

                    p1 = self.dnplan.dot(self.nsurf2);
                    p2 = self.nplan.dot(self.dns1u2);
                    p12 = self.dnplan.dot(self.dns1u2);

                    self.d2ndtu2 = (invnorm2 * p12 + smallterm * p2 + uterm * p1
                        + grosterme * ndotns2)
                        * self.nplan
                        + (invnorm2 * p2 + uterm * ndotns2) * self.dnplan
                        + (-smallterm) * self.dns1u2;
                    self.d2ndtu2 -= grosterme * self.nsurf2;

                    // OCCT: resul.SetLinearForm(-ray2, d2ndtu2, -sg2 * dray, dndu2);
                    resul = -ray2 * self.d2ndtu2 - (self.sg2 * self.dray) * self.dndu2;
                    self.d2edxdt[1][2] = resul.x;
                    self.d2edxdt[2][2] = resul.y;
                    self.d2edxdt[3][2] = resul.z;

                    //--> Derived compared to v2 and t
                    // Derived from the norm
                    uterm = ncrossns2.dot(self.nplan.cross(self.dns1v2));
                    d_secn = self.nplan.cross(self.dns1v2).dot(self.dnplan.cross(self.nsurf2))
                        + ncrossns2.dot(self.dnplan.cross(self.dns1v2));
                    grosterme = (3.0 * uterm * tterm * carre - d_secn) * cube;
                    uterm *= -cube;

                    p1 = self.dnplan.dot(self.nsurf2);
                    p2 = self.nplan.dot(self.dns1v2);
                    p12 = self.dnplan.dot(self.dns1v2);

                    self.d2ndtv2 = (invnorm2 * p12 + smallterm * p2 + uterm * p1
                        + grosterme * ndotns2)
                        * self.nplan
                        + (invnorm2 * p2 + uterm * ndotns2) * self.dnplan
                        + (-smallterm) * self.dns1v2;
                    self.d2ndtv2 -= grosterme * self.nsurf2;

                    resul = -ray2 * self.d2ndtv2 - (self.sg2 * self.dray) * self.dndv2;
                    self.d2edxdt[1][3] = resul.x;
                    self.d2edxdt[2][3] = resul.y;
                    self.d2edxdt[3][3] = resul.z;

                    //  ---------- Double derivation on t -----------------------------
                    // Derived from n1 compared to w
                    carre = invnorm1 * invnorm1;
                    cube = carre * invnorm1;
                    // Derived from the norm
                    d_prim = ncrossns1.dot(self.dnplan.cross(self.nsurf1));
                    smallterm = -2.0 * d_prim * cube;
                    // OCCT: DSecn = (dnplan.Crossed(nsurf1)).SquareMagnitude()
                    //           + ncrossns1.Dot(d2nplan.Crossed(nsurf1));
                    d_secn = self.dnplan.cross(self.nsurf1).dot(self.dnplan.cross(self.nsurf1))
                        + ncrossns1.dot(self.d2nplan.cross(self.nsurf1));
                    let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                    p1 = self.dnplan.dot(self.nsurf1);
                    p2 = self.d2nplan.dot(self.nsurf1);

                    temp = (grosterme * ndotns1) * self.nplan + (-grosterme) * self.nsurf1;
                    // OCCT: d2n1w.SetLinearForm(smallterm * p1 + invnorm1 * p2, nplan,
                    //           smallterm * ndotns1 + 2 * invnorm1 * p1, dnplan,
                    //           ndotns1 * invnorm1, d2nplan);
                    self.d2n1w = (smallterm * p1 + invnorm1 * p2) * self.nplan
                        + (smallterm * ndotns1 + 2.0 * invnorm1 * p1) * self.dnplan
                        + (ndotns1 * invnorm1) * self.d2nplan;
                    // OCCT: d2n1w += temp;
                    self.d2n1w += temp;

                    // Derived from n2 compared to w
                    carre = invnorm2 * invnorm2;
                    cube = carre * invnorm2;
                    // Derived from the norm
                    d_prim = ncrossns2.dot(self.dnplan.cross(self.nsurf2));
                    smallterm = -2.0 * d_prim * cube;
                    d_secn = self.dnplan.cross(self.nsurf2).dot(self.dnplan.cross(self.nsurf2))
                        + ncrossns2.dot(self.d2nplan.cross(self.nsurf2));
                    let grosterme = (3.0 * d_prim * d_prim * carre - d_secn) * cube;

                    p1 = self.dnplan.dot(self.nsurf2);
                    p2 = self.d2nplan.dot(self.nsurf2);

                    temp = (grosterme * ndotns2) * self.nplan + (-grosterme) * self.nsurf2;
                    self.d2n2w = (smallterm * p1 + invnorm2 * p2) * self.nplan
                        + (smallterm * ndotns2 + 2.0 * invnorm2 * p1) * self.dnplan
                        + (ndotns2 * invnorm2) * self.d2nplan;
                    self.d2n2w += temp;

                    // OCCT: temp.SetXYZ((pts1.XYZ() + pts2.XYZ()) / 2 - ptgui.XYZ());
                    temp = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
                    // OCCT: D2EDT2(1) = d2nplan.Dot(temp) - 2 * dnplan.Dot(d1gui)
                    //           - nplan.Dot(d2gui);
                    self.d2edt2[0] = self.d2nplan.dot(temp)
                        - 2.0 * self.dnplan.dot(self.d1gui)
                        - self.nplan.dot(self.d2gui);

                    // OCCT: resul.SetLinearForm(ray1, d2n1w, 2 * sg1 * dray, dn1w,
                    //           sg1 * d2ray, n1);
                    resul = ray1 * self.d2n1w
                        + (2.0 * self.sg1 * self.dray) * self.dn1w
                        + (self.sg1 * self.d2ray) * n1;
                    // OCCT: temp.SetLinearForm(ray2, d2n2w, 2 * sg2 * dray, dn2w,
                    //           sg2 * d2ray, n2);
                    temp = ray2 * self.d2n2w
                        + (2.0 * self.sg2 * self.dray) * self.dn2w
                        + (self.sg2 * self.d2ray) * n2;
                    // OCCT: resul -= temp;
                    resul -= temp;

                    self.d2edt2[1] = resul.x;
                    self.d2edt2[2] = resul.y;
                    self.d2edt2[3] = resul.z;
                }
            }
        }
        true
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

    /// OCCT Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)
    /// (BlendFunc_EvolRad.cxx L1058-1109).
    #[allow(clippy::too_many_arguments)]
    pub fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut DVec3,
        tg_last: &mut DVec3,
        norm_first: &mut DVec3,
        norm_last: &mut DVec3,
    ) {
        let mut ns1: DVec3;

        if (u1 != self.xval[0]) || (v1 != self.xval[1]) || (u2 != self.xval[2]) || (v2 != self.xval[3])
        {
            // OCCT: surf1->D1(U1, V1, bid, d1u, d1v);
            //       NmF = ns1 = d1u.Crossed(d1v);
            let (_, du, dv) = self.surf1.derivatives(u1, v1);
            ns1 = du.cross(dv);
            *norm_first = ns1;
            // OCCT: surf2->D1(U2, V2, bid, d1u, d1v); NmL = d1u.Crossed(d1v);
            let (_, du, dv) = self.surf2.derivatives(u2, v2);
            *norm_last = du.cross(dv);
        } else {
            // OCCT: NmF = ns1 = nsurf1; NmL = nsurf2;
            ns1 = self.nsurf1;
            *norm_first = ns1;
            *norm_last = self.nsurf2;
        }

        // OCCT: invnorm1 = nplan.Crossed(ns1).Magnitude();
        let mut invnorm1 = self.nplan.cross(ns1).length();
        if invnorm1 < super::brep_blend_func_evolrad::EPS {
            invnorm1 = 1.0;
        } else {
            invnorm1 = 1.0 / invnorm1;
        }
        // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) * invnorm1, nplan, -invnorm1, ns1);
        ns1 = (self.nplan.dot(ns1) * invnorm1) * self.nplan + (-invnorm1) * ns1;

        // OCCT: Center.SetXYZ(pts1.XYZ() + sg1 * ray * ns1.XYZ());
        let center = self.pts1 + self.sg1 * self.ray * ns1;

        // OCCT: TgF = nplan.Crossed(gp_Vec(Center, pts1));
        *tg_first = self.nplan.cross(self.pts1 - center);
        // OCCT: TgL = nplan.Crossed(gp_Vec(Center, pts2));
        *tg_last = self.nplan.cross(self.pts2 - center);
        if self.choix % 2 == 1 {
            *tg_first = -*tg_first;
            *tg_last = -*tg_last;
        }
    }

    /// OCCT Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)
    /// (BlendFunc_EvolRad.cxx L1135-1193) — useful for a quick and
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
        let norm1 = if norm1 < super::brep_blend_func_evolrad::EPS {
            1.0 // Unsatisfactory, but it is not necessary to stop
        } else {
            norm1
        };
        // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) / norm1, nplan, -1. / norm1, ns1);
        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;

        // OCCT: Center.SetXYZ(pts1.XYZ() + sg1 * ray * ns1.XYZ());
        let center = self.pts1 + self.sg1 * self.ray * ns1;

        // ns1 is oriented from the center to pts1
        if self.sg1 > 0.0 {
            ns1 = -ns1;
        }
        if self.choix % 2 != 0 {
            np = -np;
        }
        // OCCT: C.SetRadius(std::abs(ray)); C.SetPosition(gp_Ax2(Center, np, ns1));
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx.)
        c.radius = self.ray.abs();
        c.center = center;
        c.normal = np;
        c.x_dir = ns1;
        c.y_dir = np.cross(ns1);
        *pdeb = 0.0;
        // OCCT: Pfin = ElCLib::Parameter(C, pts2);
        *pfin = elclib_circle_parameter(c, self.pts2);
        // Test of negative and almost null angles : Single Case
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns1);
            *pfin = elclib_circle_parameter(c, self.pts2);
        }
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT Section(P, Poles, Poles2d, Weights) (BlendFunc_EvolRad.cxx
    /// L1376-1454).
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

        // Calculation and storage of distmin
        self.compute_values(&x, 0, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avoid crashing the fields !
        let mut ns1 = self.nsurf1;
        let mut ns2 = self.nsurf2;
        let mut np = self.nplan;

        // OCCT: Poles2d(Poles2d.Lower()).SetCoord(X(1), X(2)); ...
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
        if norm1 < super::brep_blend_func_evolrad::EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < super::brep_blend_func_evolrad::EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + sg1 * ray * ns1.XYZ());
        let center = self.pts1 + self.sg1 * self.ray * ns1;

        // ns1 (resp. ns2) is oriented from center to pts1 (resp. pts2),
        // and the trihedron ns1,ns2,nplan is made direct.

        if self.sg1 > 0.0 {
            ns1 = -ns1;
        }
        if self.sg2 > 0.0 {
            ns2 = -ns2;
        }
        if self.choix % 2 != 0 {
            np = -np;
        }

        // OCCT L1453: GeomFill::GetCircle(myTConv, ns1, ns2, np, pts1, pts2,
        // std::abs(ray), Center, Poles, Weights).
        get_circle(
            tconv(self.my_t_conv),
            ns1,
            ns2,
            np,
            self.pts1,
            self.pts2,
            self.ray.abs(),
            center,
            poles,
            weigths,
        );
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BlendFunc_EvolRad.cxx L1458-1634) — used for the first and last
    /// section.
    #[allow(clippy::too_many_arguments)]
    #[allow(unreachable_code)] // the math_SVD fallback is a pending kernel gap
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
        // OCCT locals default-constructed (null) — rcad binds them to ZERO.
        let mut dnorm1w = DVec3::ZERO;
        let mut dnorm2w = DVec3::ZERO;

        // OCCT: gp_Pnt Center; (assigned in the circle case below)
        let center;
        // OCCT: math_Vector sol(1, 4), secmember(1, 4);
        let mut sol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];

        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()
        let mut istgt = true;

        // OCCT: P.ParametersOnS1(sol(1), sol(2)); P.ParametersOnS2(sol(3), sol(4));
        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        sol[0] = u1;
        sol[1] = v1;
        sol[2] = u2;
        sol[3] = v2;

        // Calculation of equations
        self.compute_values(&sol, 1, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avoid crashing fields !
        ns1 = self.nsurf1;
        ns2 = self.nsurf2;
        np = self.nplan;
        let mut dnp = self.dnplan;
        let mut rayprim = self.dray;

        // OCCT: if (!pts1.IsEqual(pts2, 1.e-4))
        if self.pts1.distance(self.pts2) > 1.0e-4 {
            // Calculation of derived  Normal processing
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
                let mut x = VecD::new(4);
                for i in 1..=4 {
                    x.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=4 {
                    secmember[i - 1] = x.get(i);
                }
                istgt = false;
            }
        }

        if istgt {
            // OCCT L1503-1511: math_SVD SingRS(DEDX);
            // if (SingRS.IsDone()) { SingRS.Solve(-DEDT, secmember, 1.e-6);
            //                        istgt = false; }
            // GAP (plan 0.6): math_SVD is pending in rcad-kernel; the OCCT
            // !IsDone() path (istgt stays true) is preserved.
        }

        if !istgt {
            // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
            self.tg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1;
            self.tg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2;

            // OCCT: dnorm1w.SetLinearForm(secmember(1), dndu1, secmember(2), dndv1, dn1w);
            dnorm1w = secmember[0] * self.dndu1 + secmember[1] * self.dndv1 + self.dn1w;
            dnorm2w = secmember[2] * self.dndu2 + secmember[3] * self.dndv2 + self.dn2w;

            istgt = false;
        }

        // Tops 2D
        // OCCT: Poles2d(Poles2d.Lower()).SetCoord(sol(1), sol(2)); ...
        poles_2d[0] = DVec2::new(sol[0], sol[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
        if !istgt {
            d_poles_2d[0] = DVec2::new(secmember[0], secmember[1]);
            d_poles_2d[d_poles_2d.len() - 1] = DVec2::new(secmember[2], secmember[3]);
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
        let mut norm1 = self.nplan.cross(ns1).length();
        let mut norm2 = self.nplan.cross(ns2).length();
        if norm1 < super::brep_blend_func_evolrad::EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < super::brep_blend_func_evolrad::EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + sg1 * ray * ns1.XYZ());
        center = self.pts1 + self.sg1 * self.ray * ns1;
        let mut tgc = DVec3::ZERO;
        if !istgt {
            // OCCT: tgc.SetLinearForm(sg1 * ray, dnorm1w, sg1 * dray, ns1, tg1);
            tgc = self.sg1 * self.ray * dnorm1w + (self.sg1 * self.dray) * ns1 + self.tg1;
        }

        // ns1 is oriented from center to pts1, and  ns2 from center to pts2
        // and the trihedron ns1,ns2,nplan is made direct

        if self.sg1 > 0.0 {
            ns1 = -ns1;
            if !istgt {
                dnorm1w = -dnorm1w;
            }
        }
        if self.sg2 > 0.0 {
            ns2 = -ns2;
            if !istgt {
                dnorm2w = -dnorm2w;
            }
        }
        if self.choix % 2 != 0 {
            np = -np;
            dnp = -dnp;
        }

        // OCCT: if (ray < 0.) { rayprim = -rayprim; } — to avoid
        // std::abs(dray) some lines below
        if self.ray < 0.0 {
            rayprim = -rayprim;
        }

        if !istgt {
            // OCCT L1609-1627: GeomFill::GetCircle(myTConv, ns1, ns2, dnorm1w,
            // dnorm2w, np, dnp, pts1, pts2, tg1, tg2, std::abs(ray), rayprim,
            // Center, tgc, Poles, DPoles, Weights, DWeights).
            get_circle_d1(
                tconv(self.my_t_conv),
                ns1,
                ns2,
                dnorm1w,
                dnorm2w,
                np,
                dnp,
                self.pts1,
                self.pts2,
                self.tg1,
                self.tg2,
                self.ray.abs(),
                rayprim,
                center,
                tgc,
                poles,
                d_poles,
                weigths,
                d_weigths,
            )
        } else {
            // OCCT L1631: GeomFill::GetCircle(myTConv, ns1, ns2, np, pts1,
            // pts2, std::abs(ray), Center, Poles, Weights);
            get_circle(
                tconv(self.my_t_conv),
                ns1,
                ns2,
                np,
                self.pts1,
                self.pts2,
                self.ray.abs(),
                center,
                poles,
                weigths,
            );
            false
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BlendFunc_EvolRad.cxx L1638-2019) —
    /// used for the first and last section.
    #[allow(clippy::too_many_arguments)]
    #[allow(unreachable_code)] // the math_SVD fallback is a pending kernel gap
    pub fn section_d2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        d2_poles: &mut [DVec3],
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
        // OCCT locals default-constructed (null) — rcad binds them to ZERO.
        let mut dnorm1w = DVec3::ZERO;
        let mut dnorm2w = DVec3::ZERO;
        let mut d2norm1w = DVec3::ZERO;
        let mut d2norm2w = DVec3::ZERO;
        let mut dtg1 = DVec3::ZERO;
        let mut dtg2 = DVec3::ZERO;

        // OCCT: gp_Pnt Center; (assigned in the circle case below)
        let center;
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

        // OCCT: P.ParametersOnS1(X(1), X(2)); P.ParametersOnS2(X(3), X(4));
        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        x[0] = u1;
        x[1] = v1;
        x[2] = u2;
        x[3] = v2;

        // Calculs des equations
        // OCCT: ComputeValues(X, 2, true, prm);
        self.compute_values(&x, 2, true, prm);
        self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

        // ns1, ns2, np are copied locally to avoid crashing the fields
        ns1 = self.nsurf1;
        ns2 = self.nsurf2;
        np = self.nplan;
        let mut dnp = self.dnplan;
        let mut d2np = self.d2nplan;
        let mut rayprim = self.dray;
        let mut raysecn = self.d2ray;

        // OCCT: if (!pts1.IsEqual(pts2, 1.e-4))
        if self.pts1.distance(self.pts2) > 1.0e-4 {
            // OCCT: math_Gauss Resol(DEDX, 1.e-9); // Tolerance to precise
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, self.dedx[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            // Calculation of derived Normal Processing
            if resol.is_done() {
                // OCCT: Resol.Solve(-DEDT, sol);
                let mut xg = VecD::new(4);
                for i in 1..=4 {
                    xg.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut xg);
                for i in 1..=4 {
                    sol[i - 1] = xg.get(i);
                }
                // OCCT: D2EDX2.Multiply(sol, D2DXdSdt);
                self.d2edx2.multiply(&sol, &mut d2dxdsdt);
                // OCCT: secmember = -(D2EDT2 + (2 * D2EDXDT + D2DXdSdt) * sol);
                for i in 1..=4 {
                    let mut somme = 0.0;
                    for j in 1..=4 {
                        somme +=
                            (2.0 * self.d2edxdt[i - 1][j - 1] + d2dxdsdt[i - 1][j - 1]) * sol[j - 1];
                    }
                    secmember[i - 1] = -(self.d2edt2[i - 1] + somme);
                }
                // OCCT: Resol.Solve(secmember);
                let mut xb = VecD::new(4);
                for i in 1..=4 {
                    xb.set(i, secmember[i - 1]);
                }
                resol.solve(&mut xb);
                for i in 1..=4 {
                    secmember[i - 1] = xb.get(i);
                }
                istgt = false;
            }
        }

        if istgt {
            // OCCT L1829-1841: math_SVD SingRS(DEDX); math_Vector Vbis(1, 4);
            // if (SingRS.IsDone()) { SingRS.Solve(-DEDT, sol, 1.e-6);
            //   D2EDX2.Multiply(sol, D2DXdSdt);
            //   Vbis = -(D2EDT2 + (2 * D2EDXDT + D2DXdSdt) * sol);
            //   SingRS.Solve(Vbis, secmember, 1.e-6); istgt = false; }
            // GAP (plan 0.6): math_SVD is pending in rcad-kernel; the OCCT
            // !IsDone() path (istgt stays true) is preserved.
        }

        if !istgt {
            // OCCT: tg1.SetLinearForm(sol(1), d1u1, sol(2), d1v1);
            self.tg1 = sol[0] * self.d1u1 + sol[1] * self.d1v1;
            self.tg2 = sol[2] * self.d1u2 + sol[3] * self.d1v2;

            // OCCT: dnorm1w.SetLinearForm(sol(1), dndu1, sol(2), dndv1, dn1w);
            dnorm1w = sol[0] * self.dndu1 + sol[1] * self.dndv1 + self.dn1w;
            dnorm2w = sol[2] * self.dndu2 + sol[3] * self.dndv2 + self.dn2w;
            // OCCT: temp.SetLinearForm(sol(1) * sol(1), d2u1, 2 * sol(1) * sol(2), d2uv1,
            //           sol(2) * sol(2), d2v1);
            let mut temp = (sol[0] * sol[0]) * self.d2u1
                + (2.0 * sol[0] * sol[1]) * self.d2uv1
                + (sol[1] * sol[1]) * self.d2v1;

            // OCCT: dtg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1, temp);
            dtg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1 + temp;

            // OCCT: temp.SetLinearForm(sol(3) * sol(3), d2u2, 2 * sol(3) * sol(4), d2uv2,
            //           sol(4) * sol(4), d2v2);
            temp = (sol[2] * sol[2]) * self.d2u2
                + (2.0 * sol[2] * sol[3]) * self.d2uv2
                + (sol[3] * sol[3]) * self.d2v2;
            // OCCT: dtg2.SetLinearForm(secmember(3), d1u2, secmember(4), d1v2, temp);
            dtg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2 + temp;

            // OCCT: temp.SetLinearForm(sol(1) * sol(1), d2ndu1, 2 * sol(1) * sol(2),
            //           d2nduv1, sol(2) * sol(2), d2ndv1);
            temp = (sol[0] * sol[0]) * self.d2ndu1
                + (2.0 * sol[0] * sol[1]) * self.d2nduv1
                + (sol[1] * sol[1]) * self.d2ndv1;

            // OCCT: tempbis.SetLinearForm(2 * sol(1), d2ndtu1, 2 * sol(2), d2ndtv1, d2n1w);
            let tempbis = (2.0 * sol[0]) * self.d2ndtu1
                + (2.0 * sol[1]) * self.d2ndtv1
                + self.d2n1w;
            // OCCT: temp += tempbis;
            temp += tempbis;
            // OCCT: d2norm1w.SetLinearForm(secmember(1), dndu1, secmember(2), dndv1, temp);
            d2norm1w = secmember[0] * self.dndu1 + secmember[1] * self.dndv1 + temp;

            // OCCT: temp.SetLinearForm(sol(3) * sol(3), d2ndu2, 2 * sol(3) * sol(4),
            //           d2nduv2, sol(4) * sol(4), d2ndv2);
            temp = (sol[2] * sol[2]) * self.d2ndu2
                + (2.0 * sol[2] * sol[3]) * self.d2nduv2
                + (sol[3] * sol[3]) * self.d2ndv2;
            // OCCT: tempbis.SetLinearForm(2 * sol(3), d2ndtu2, 2 * sol(4), d2ndtv2, d2n2w);
            let tempbis = (2.0 * sol[2]) * self.d2ndtu2
                + (2.0 * sol[3]) * self.d2ndtv2
                + self.d2n2w;
            // OCCT: temp += tempbis;
            temp += tempbis;
            // OCCT: d2norm2w.SetLinearForm(secmember(3), dndu2, secmember(4), dndv2, temp);
            d2norm2w = secmember[2] * self.dndu2 + secmember[3] * self.dndv2 + temp;
        }

        // Tops 2d
        // OCCT: Poles2d(Poles2d.Lower()).SetCoord(X(1), X(2)); ...
        poles_2d[0] = DVec2::new(x[0], x[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(x[2], x[3]);
        if !istgt {
            d_poles_2d[0] = DVec2::new(sol[0], sol[1]);
            d_poles_2d[d_poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
            d2_poles_2d[0] = DVec2::new(secmember[0], secmember[1]);
            d2_poles_2d[d2_poles_2d.len() - 1] = DVec2::new(secmember[2], secmember[3]);
        }

        // the linear is processed...
        if self.my_s_shape == super::brep_blend_func::BlendFuncSectionShape::Linear {
            poles[low] = self.pts1;
            poles[upp] = self.pts2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            if !istgt {
                // OCCT L1899-1902: DPoles(low) = tg1; DPoles(upp) = tg2;
                // DPoles(low) = dtg1; DPoles(upp) = dtg2; — the first two
                // assignments are overwritten by the second two in the OCCT
                // source (kept literally, bug-compatible).
                d_poles[low] = dtg1;
                d_poles[upp] = dtg2;
                d_weigths[low] = 0.0;
                d_weigths[upp] = 0.0;
                d2_weigths[low] = 0.0;
                d2_weigths[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle
        let mut norm1 = self.nplan.cross(ns1).length();
        let mut norm2 = self.nplan.cross(ns2).length();
        if norm1 < super::brep_blend_func_evolrad::EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < super::brep_blend_func_evolrad::EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        // OCCT: Center.SetXYZ(pts1.XYZ() + sg1 * ray * ns1.XYZ());
        center = self.pts1 + self.sg1 * self.ray * ns1;
        let mut tgc = DVec3::ZERO;
        let mut dtgc = DVec3::ZERO;
        if !istgt {
            // OCCT: tgc.SetLinearForm(sg1 * ray, dnorm1w, sg1 * dray, ns1, tg1);
            tgc = self.sg1 * self.ray * dnorm1w + (self.sg1 * self.dray) * ns1 + self.tg1;
            // OCCT: dtgc.SetLinearForm(sg1 * ray, d2norm1w, 2 * sg1 * dray, dnorm1w,
            //           sg1 * d2ray, ns1);
            dtgc = self.sg1 * self.ray * d2norm1w
                + (2.0 * self.sg1 * self.dray) * dnorm1w
                + (self.sg1 * self.d2ray) * ns1;
            // OCCT: dtgc += dtg1;
            dtgc += dtg1;
        }

        // ns1 is oriented from the center to pts1, and ns2 from the center to pts2
        // and the trihedron ns1,ns2,nplan is made direct

        if self.sg1 > 0.0 {
            ns1 = -ns1;
            if !istgt {
                dnorm1w = -dnorm1w;
                d2norm1w = -d2norm1w;
            }
        }
        if self.sg2 > 0.0 {
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

        // OCCT: if (ray < 0.) { rayprim = -rayprim; raysecn = -raysecn; } —
        // to avoid std::abs(dray) several lines below
        if self.ray < 0.0 {
            rayprim = -rayprim;
            raysecn = -raysecn;
        }

        if !istgt {
            // OCCT L1976-2003: GeomFill::GetCircle(myTConv, ns1, ns2, dnorm1w,
            // dnorm2w, d2norm1w, d2norm2w, np, dnp, d2np, pts1, pts2, tg1,
            // tg2, dtg1, dtg2, std::abs(ray), rayprim, raysecn, Center, tgc,
            // dtgc, Poles, DPoles, D2Poles, Weights, DWeights, D2Weights).
            get_circle_d2(
                tconv(self.my_t_conv),
                ns1,
                ns2,
                dnorm1w,
                dnorm2w,
                d2norm1w,
                d2norm2w,
                np,
                dnp,
                d2np,
                self.pts1,
                self.pts2,
                self.tg1,
                self.tg2,
                dtg1,
                dtg2,
                self.ray.abs(),
                rayprim,
                raysecn,
                center,
                tgc,
                dtgc,
                poles,
                d_poles,
                d2_poles,
                weigths,
                d_weigths,
                d2_weigths,
            )
        } else {
            // OCCT L2007-2016: GeomFill::GetCircle(myTConv, ns1, ns2, nplan,
            // pts1, pts2, std::abs(ray), Center, Poles, Weights).
            get_circle(
                tconv(self.my_t_conv),
                ns1,
                ns2,
                self.nplan,
                self.pts1,
                self.pts2,
                self.ray.abs(),
                center,
                poles,
                weigths,
            );
            false
        }
    }
}
