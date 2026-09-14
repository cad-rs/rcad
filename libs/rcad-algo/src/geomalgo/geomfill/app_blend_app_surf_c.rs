//! OCCT AppBlend_AppSurf — the Perform(Lin, F, NbMaxP) body
//! (AppBlend_AppSurf.gxx L597-1049) of the engine declared in the parent
//! module `app_blend_app_surf`, plus the AppParCurves_MultiCurve re-host
//! helpers consumed by the body.

use glam::{DVec2, DVec3};
use rcad_kernel::math::math_matrix::Vector as RVector;

use crate::geomalgo::app_def::{MultiLine, MultiPointConstraint};
use crate::geomalgo::app_def_compute::Compute;
use crate::geomalgo::approx_int::{
    ApproxParamType, AppParConstraint, MCurvesToBSpCurve, MultiCurve,
};

use super::{
    AppBlendAppSurf, TheSectionGenerator, GP_RESOLUTION, REAL_FIRST, REAL_LAST, SCAL,
};
use super::super::line::Line;

/// OCCT AppParCurves_MultiCurve::Pole(CuIndex, Nieme)
/// (AppParCurves_MultiCurve.cxx L119-127) —
/// `return tabPoint->Value(Nieme).Point(CuIndex);`
/// (the rcad MultiCurve re-host over the public pole storage).
pub(super) fn multi_curve_pole(mucu: &MultiCurve, cu_index: usize, nieme: usize) -> DVec3 {
    mucu.value(nieme).point(cu_index)
}

/// OCCT AppParCurves_MultiCurve::Transform(CuIndex, x, dx, y, dy, z, dz)
/// (AppParCurves_MultiCurve.cxx L146-163) —
/// `for (i = 1; i <= tabPoint->Length(); i++)
///    (tabPoint->ChangeValue(i)).Transform(CuIndex, x, dx, y, dy, z, dz);`
/// (the Dimension != 3 range check is carried by MultiPoint::transform).
pub(super) fn multi_curve_transform(
    mucu: &mut MultiCurve,
    cu_index: usize,
    x: f64,
    dx: f64,
    y: f64,
    dy: f64,
    z: f64,
    dz: f64,
) {
    for point in mucu.poles.iter_mut() {
        point.transform(cu_index, x, dx, y, dy, z, dz);
    }
}

/// OCCT AppParCurves_MultiCurve::Transform2d(CuIndex, x, dx, y, dy)
/// (AppParCurves_MultiCurve.cxx L165-179).
pub(super) fn multi_curve_transform2d(
    mucu: &mut MultiCurve,
    cu_index: usize,
    x: f64,
    dx: f64,
    y: f64,
    dy: f64,
) {
    for point in mucu.poles.iter_mut() {
        point.transform2d(cu_index, x, dx, y, dy);
    }
}

impl AppBlendAppSurf {
    /// OCCT AppBlend_AppSurf::Perform(Lin, F, NbMaxP) (gxx L597-1049).
    pub fn perform_nb_max_p<G: TheSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &mut G,
        nb_max_p: i32,
    ) {
        self.done = false;
        // OCCT L602-605: the handle null check is unrepresentable over the
        // rcad value Line.

        let mut withderiv;
        let mut cfirst = AppParConstraint::NoConstraint;
        let mut clast = AppParConstraint::NoConstraint;

        let mut mytol3d = 0.0f64;
        let mut mytol2d = 0.0f64;
        let mut new_dv: DVec3;

        self.seq_poles2d.clear();

        let nb_point_tot = lin.nb_points();

        let (mut nb_u_poles, mut nb_u_knots, mut nb_poles2d) = (0i32, 0i32, 0i32);
        f.get_shape(&mut nb_u_poles, &mut nb_u_knots, &mut self.udeg, &mut nb_poles2d);

        self.tab_u_knots = Some(vec![0.0; nb_u_knots as usize]);
        self.tab_u_mults = Some(vec![0; nb_u_knots as usize]);

        f.knots(self.tab_u_knots.as_mut().unwrap());
        f.mults(self.tab_u_mults.as_mut().unwrap());

        let mut tab_app_p = vec![DVec3::ZERO; nb_u_poles as usize];
        let mut tab_app_v = vec![DVec3::ZERO; nb_u_poles as usize];
        let mut x = REAL_LAST;
        let mut y = REAL_LAST;
        let mut z = REAL_LAST;
        let mut dx = REAL_FIRST;
        let mut dy = REAL_FIRST;
        let mut dz = REAL_FIRST;

        let mut tab_p2d = vec![DVec2::ZERO; nb_poles2d.max(1) as usize];
        let mut tab_v2d = vec![DVec2::ZERO; nb_poles2d.max(1) as usize];
        let mut x2d = vec![REAL_LAST; nb_poles2d.max(1) as usize];
        let mut y2d = vec![REAL_LAST; nb_poles2d.max(1) as usize];
        let mut dx2d = vec![REAL_FIRST; nb_poles2d.max(1) as usize];
        let mut dy2d = vec![REAL_FIRST; nb_poles2d.max(1) as usize];

        let mut tab_w = vec![0.0f64; nb_u_poles as usize];
        let mut tab_dw = vec![0.0f64; nb_u_poles as usize];

        let mut tab_app_p2d = vec![DVec2::ZERO; (nb_poles2d + nb_u_poles) as usize];
        let mut tab_app_v2d = vec![DVec2::ZERO; (nb_poles2d + nb_u_poles) as usize];

        // On calcule les boites de chaque ligne (box for all lines)
        let (mut xc, mut yc, mut zc): (f64, f64, f64);
        for i in 1..=nb_point_tot {
            f.section_d1(
                lin.point(i),
                &mut tab_app_p,
                &mut tab_app_v,
                &mut tab_p2d,
                &mut tab_v2d,
                &mut tab_w,
                &mut tab_dw,
            );
            for j in 1..=nb_u_poles {
                (xc, yc, zc) = (
                    tab_app_p[(j - 1) as usize].x,
                    tab_app_p[(j - 1) as usize].y,
                    tab_app_p[(j - 1) as usize].z,
                );
                if xc < x {
                    x = xc;
                }
                if xc > dx {
                    dx = xc;
                }
                if yc < y {
                    y = yc;
                }
                if yc > dy {
                    dy = yc;
                }
                if zc < z {
                    z = zc;
                }
                if zc > dz {
                    dz = zc;
                }
            }
            for j in 1..=nb_poles2d {
                (xc, yc) = (tab_p2d[(j - 1) as usize].x, tab_p2d[(j - 1) as usize].y);
                if xc < x2d[(j - 1) as usize] {
                    x2d[(j - 1) as usize] = xc;
                }
                if xc > dx2d[(j - 1) as usize] {
                    dx2d[(j - 1) as usize] = xc;
                }
                if yc < y2d[(j - 1) as usize] {
                    y2d[(j - 1) as usize] = yc;
                }
                if yc > dy2d[(j - 1) as usize] {
                    dy2d[(j - 1) as usize] = yc;
                }
            }
        }
        // On calcule pour chaque ligne la transformation vers 0 1.
        let seuil = 1000.0 * self.tol3d;
        let seuil2d = 1000.0 * self.tol2d;
        if (dx - x) < seuil {
            dx = 1.0;
            x = 0.0;
        } else {
            dx = 1.0 / (dx - x);
            x *= -dx;
        }
        if (dy - y) < seuil {
            dy = 1.0;
            y = 0.0;
        } else {
            dy = 1.0 / (dy - y);
            y *= -dy;
        }
        if (dz - z) < seuil {
            dz = 1.0;
            z = 0.0;
        } else {
            dz = 1.0 / (dz - z);
            z *= -dz;
        }
        for j in 1..=nb_poles2d {
            let jx = (j - 1) as usize;
            if (dx2d[jx] - x2d[jx]) < seuil2d {
                dx2d[jx] = 1.0;
                x2d[jx] = 0.0;
            } else {
                dx2d[jx] = 1.0 / (dx2d[jx] - x2d[jx]);
                x2d[jx] *= -dx2d[jx];
            }
            if (dy2d[jx] - y2d[jx]) < seuil2d {
                dy2d[jx] = 1.0;
                y2d[jx] = 0.0;
            } else {
                dy2d[jx] = 1.0 / (dy2d[jx] - y2d[jx]);
                y2d[jx] *= -dy2d[jx];
            }
        }
        if !SCAL {
            dx = 1.0;
            x = 0.0;
            dy = 1.0;
            y = 0.0;
            dz = 1.0;
            z = 0.0;
            for j in 1..=nb_poles2d {
                let jx = (j - 1) as usize;
                dx2d[jx] = 1.0;
                x2d[jx] = 0.0;
                dy2d[jx] = 1.0;
                y2d[jx] = 0.0;
            }
        }
        //  modified by eap Thu Jan  3 14:45:22 2002 ___BEGIN___
        // Keep "inter-troncons" parameters, not only first and last
        let mut a_param_seq: Vec<f64> = Vec::new();
        if self.knownp {
            a_param_seq.push(f.parameter(lin.point(1)));
        }
        //  modified by EAP Thu Jan  3 14:45:41 2002 ___END___

        let mut concat = MCurvesToBSpCurve::new();

        // On calcule le nombre de troncons.
        let mut nbtronc = nb_point_tot / nb_max_p;
        let mut reste = nb_point_tot - (nbtronc * nb_max_p);
        // On regarde si il faut prendre un troncon de plus.
        let mut nmax = nb_max_p;
        if nbtronc > 0 && reste > 0 {
            nmax = nb_point_tot / (nbtronc + 1);
            if nmax > (2 * nb_max_p) / 3 {
                nbtronc += 1;
                reste = nb_point_tot - (nbtronc * nmax);
            } else {
                nmax = nb_max_p;
            }
        } else if nbtronc == 0 {
            nbtronc = 1;
            nmax = reste;
            reste = 0;
        }

        // Approximate each "troncon" with nb of Bezier's using AppDef_Compute
        // and concat them into BSpline with Approx_MCurvesToBSpCurve

        let mut troncsize = vec![0i32; nbtronc as usize];
        let mut troncstart = vec![0i32; nbtronc as usize];

        let rab = reste / nbtronc + 1;
        let mut start = 1i32;
        for itronc in 1..=nbtronc {
            troncstart[(itronc - 1) as usize] = start;
            let rabrab = rab.min(reste);
            if reste > 0 {
                reste -= rabrab;
            }
            troncsize[(itronc - 1) as usize] = nmax + rabrab + 1;
            start += nmax + rabrab;
        }
        let last_tronc = (nbtronc - 1) as usize;
        troncsize[last_tronc] = troncsize[last_tronc] - 1;
        for itronc in 1..=nbtronc {
            let nb_point = troncsize[(itronc - 1) as usize];
            let st_point = troncstart[(itronc - 1) as usize];
            let mut mult_p;
            let mut mult_l = MultiLine::new_nb_mult(nb_point as usize);

            for i in 1..=nb_point {
                let i_lin = st_point + i - 1;
                withderiv = f.section_d1(
                    lin.point(i_lin),
                    &mut tab_app_p,
                    &mut tab_app_v,
                    &mut tab_p2d,
                    &mut tab_v2d,
                    &mut tab_w,
                    &mut tab_dw,
                );
                if super::app_blend_get_context_approx_with_no_tgt() {
                    withderiv = false;
                }

                for j in 1..=nb_poles2d {
                    let (p2x, p2y) =
                        (tab_p2d[(j - 1) as usize].x, tab_p2d[(j - 1) as usize].y);
                    tab_app_p2d[(j - 1) as usize] = DVec2::new(
                        dx2d[(j - 1) as usize] * p2x + x2d[(j - 1) as usize],
                        dy2d[(j - 1) as usize] * p2y + y2d[(j - 1) as usize],
                    );
                    if withderiv {
                        let (v2x, v2y) =
                            (tab_v2d[(j - 1) as usize].x, tab_v2d[(j - 1) as usize].y);
                        tab_app_v2d[(j - 1) as usize] = DVec2::new(
                            dx2d[(j - 1) as usize] * v2x,
                            dy2d[(j - 1) as usize] * v2y,
                        );
                    }
                }
                for j in 1..=nb_u_poles {
                    // pour les courbes rationnelles il faut multiplier les poles par
                    // leurs poids respectifs
                    if withderiv {
                        tab_app_v2d[(nb_poles2d + j - 1) as usize] =
                            DVec2::new(tab_dw[(j - 1) as usize], 0.0);
                        new_dv = tab_app_p[(j - 1) as usize] * tab_dw[(j - 1) as usize]
                            + tab_app_v[(j - 1) as usize] * tab_w[(j - 1) as usize];
                        tab_app_v[(j - 1) as usize] =
                            DVec3::new(dx * new_dv.x, dy * new_dv.y, dz * new_dv.z);
                    }
                    tab_app_p[(j - 1) as usize] *= tab_w[(j - 1) as usize];
                    tab_app_p2d[(nb_poles2d + j - 1) as usize] =
                        DVec2::new(tab_w[(j - 1) as usize], 0.0);
                    let (px, py, pz) = (
                        tab_app_p[(j - 1) as usize].x,
                        tab_app_p[(j - 1) as usize].y,
                        tab_app_p[(j - 1) as usize].z,
                    );
                    tab_app_p[(j - 1) as usize] =
                        DVec3::new(dx * px + x, dy * py + y, dz * pz + z);
                }
                if withderiv {
                    mult_p = MultiPointConstraint::new_tangency(
                        &tab_app_p,
                        &tab_app_p2d,
                        &tab_app_v,
                        &tab_app_v2d,
                    );
                    if i == 1 {
                        cfirst = AppParConstraint::TangencyPoint;
                    } else if i == nb_point {
                        clast = AppParConstraint::TangencyPoint;
                    }
                } else {
                    mult_p = MultiPointConstraint::new_tab_p_p2d(&tab_app_p, &tab_app_p2d);
                    if i == 1 {
                        cfirst = AppParConstraint::PassPoint;
                    } else if i == nb_point {
                        clast = AppParConstraint::PassPoint;
                    }
                }
                mult_l.set_value(i as usize, &mult_p);
            }

            // IFV 04.06.07 occ13904
            if nb_point == 2 {
                self.dmin = 1;
                if cfirst == AppParConstraint::PassPoint
                    && clast == AppParConstraint::PassPoint
                {
                    self.dmax = 1;
                }
            }

            //  modified by EAP Thu Jan  3 15:44:13 2002 ___BEGIN___
            let mut u_floc = 0.0f64;
            let mut u_lloc = 0.0f64;
            // OCCT L903: AppDef_Compute theapprox(dmin, dmax, tol3d, tol2d,
            // nbit) — the OCCT cutting/parametrization/Squares defaults
            // (true / ChordLength / false).
            let mut theapprox = Compute::new(
                self.dmin,
                self.dmax,
                self.tol3d,
                self.tol2d,
                self.nbit,
                true,
                ApproxParamType::ChordLength,
                false,
            );
            if self.knownp {
                let mut the_params = RVector::new(1, nb_point);
                // On recale les parametres entre 0 et 1.
                u_floc = f.parameter(lin.point(st_point));
                u_lloc = f.parameter(lin.point(st_point + nb_point - 1));
                //  modified by EAP Thu Jan  3 15:45:17 2002 ___END___
                for i in 1..=nb_point {
                    let i_lin = st_point + i - 1;
                    the_params.set(i, (f.parameter(lin.point(i_lin)) - u_floc) / (u_lloc - u_floc));
                }
                let the_app_def1 = Compute::new_with_parameters(
                    &the_params,
                    self.dmin,
                    self.dmax,
                    self.tol3d,
                    self.tol2d,
                    self.nbit,
                    true,
                    true,
                );
                theapprox = the_app_def1;
            }
            theapprox.set_constraints(cfirst, clast);
            theapprox.perform(&mult_l);

            //  modified by EAP Thu Jan  3 16:00:43 2002 ___BEGIN___
            // To know internal parameters if multicurve is approximated by several Bezier's
            let mut a_pole_dist_seq: Vec<f64> = Vec::new();
            let mut a_whole_dist = 0.0f64;
            //  modified by EAP Thu Jan  3 16:45:48 2002 ___END___
            // OCCT re-reads theapprox.NbMultiCurves() inside the loop; the
            // count is loop-invariant (the transformations below only touch
            // pole values), so the rcad form hoists it before the mutable
            // borrow (Rust borrow discipline).
            let nb_multi_curves = theapprox.nb_multi_curves();
            for index in 1..=nb_multi_curves {
                // OCCT L930: AppParCurves_MultiCurve& mucu =
                // theapprox.ChangeValue(Index); — the rcad error() query is
                // evaluated before the mutable borrow (Rust borrow
                // discipline; OCCT aliases both over the same object).
                let (the_tol3d, the_tol2d) = theapprox.error(index);
                let mucu = theapprox.change_value(index);
                mytol3d = (the_tol3d / dx).max(mytol3d);
                mytol3d = (the_tol3d / dy).max(mytol3d);
                mytol3d = (the_tol3d / dz).max(mytol3d);
                for j in 1..=nb_u_poles {
                    multi_curve_transform(
                        mucu,
                        j as usize,
                        -x / dx,
                        1.0 / dx,
                        -y / dy,
                        1.0 / dy,
                        -z / dz,
                        1.0 / dz,
                    );
                }
                for j in 1..=nb_poles2d {
                    multi_curve_transform2d(
                        mucu,
                        (j + nb_u_poles) as usize,
                        -x2d[(j - 1) as usize] / dx2d[(j - 1) as usize],
                        1.0 / dx2d[(j - 1) as usize],
                        -y2d[(j - 1) as usize] / dy2d[(j - 1) as usize],
                        1.0 / dy2d[(j - 1) as usize],
                    );
                    mytol2d = (the_tol2d / dx2d[(j - 1) as usize]).max(mytol2d);
                    mytol2d = (the_tol2d / dy2d[(j - 1) as usize]).max(mytol2d);
                }
                concat.append(mucu.clone());

                //  modified by EAP Thu Jan  3 15:45:23 2002 ___BEGIN___
                if self.knownp && nb_multi_curves > 1 {
                    let a_first_pole = multi_curve_pole(mucu, index, 1);
                    let a_last_pole = multi_curve_pole(mucu, index, mucu.nb_poles());
                    a_pole_dist_seq.push(a_first_pole.distance(a_last_pole));
                    a_whole_dist += a_pole_dist_seq[a_pole_dist_seq.len() - 1];
                }
            }
            if self.knownp {
                let mut i_u = u_floc;
                for i_dist in 1..a_pole_dist_seq.len() {
                    i_u += a_pole_dist_seq[i_dist - 1] / a_whole_dist * (u_lloc - u_floc);
                    a_param_seq.push(i_u);
                }
                a_param_seq.push(u_lloc);
            }
            //  modified by EAP Thu Jan  3 15:45:27 2002 ___END___
        }
        self.tol3dreached = mytol3d;
        self.tol2dreached = mytol2d;
        concat.perform();
        let mult_c = concat.value();
        self.vdeg = mult_c.degree as i32;
        let nb_v_poles = mult_c.nb_poles();

        self.tab_poles = Some(vec![vec![DVec3::ZERO; nb_v_poles]; nb_u_poles as usize]);
        self.tab_weights = Some(vec![vec![0.0; nb_v_poles]; nb_u_poles as usize]);
        self.tab_v_knots = Some(mult_c.knots.clone());

        if self.knownp {
            //  modified by EAP Fri Jan  4 12:07:30 2002 ___BEGIN___
            if a_param_seq.len() != self.tab_v_knots.as_ref().unwrap().len() {
                Self::bspl_reparametrize(
                    f.parameter(lin.point(1)),
                    f.parameter(lin.point(lin.nb_points())),
                    self.tab_v_knots.as_mut().unwrap(),
                );
            } else {
                let mut i_tab_knot = 0usize; // tabVKnots->Lower()
                for i_knot in 1..=a_param_seq.len() {
                    self.tab_v_knots.as_mut().unwrap()[i_tab_knot] = a_param_seq[i_knot - 1];
                    i_tab_knot += 1;
                }
            }
            //  modified by EAP Fri Jan  4 12:07:35 2002 ___END___
        }

        self.tab_v_mults = Some(mult_c.mults.iter().map(|m| *m as i32).collect());

        let mut newtab_p = vec![DVec3::ZERO; nb_v_poles];
        let mut newtab_p2d: Vec<DVec2> = Vec::new();
        for j in 1..=nb_u_poles {
            mult_c.curve(j as usize, &mut newtab_p);
            mult_c.curve2d((j + nb_u_poles + nb_poles2d) as usize, &mut newtab_p2d);
            for k in 1..=nb_v_poles {
                // pour les courbes rationnelles il faut maintenant diviser
                // les poles par leurs poids respectifs
                let a_weight = newtab_p2d[(k - 1) as usize].x;
                self.tab_poles.as_mut().unwrap()[(j - 1) as usize][(k - 1) as usize] =
                    newtab_p[(k - 1) as usize] / a_weight;
                if a_weight < GP_RESOLUTION {
                    self.done = false;
                    return;
                }
                self.tab_weights.as_mut().unwrap()[(j - 1) as usize][(k - 1) as usize] = a_weight;
            }
        }

        for j in 1..=nb_poles2d {
            let mut newtab_p2d = vec![DVec2::ZERO; nb_v_poles];
            mult_c.curve2d((nb_u_poles + j) as usize, &mut newtab_p2d);
            self.seq_poles2d.push(newtab_p2d);
        }

        self.done = true;
    }
}
