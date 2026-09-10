//! OCCT BlendFunc package services (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_SectionShape.hxx (L19-24), Convert_ParameterisationType.hxx
//! (TKMath/Convert, L79-89, required by BlendFunc::GetShape's out parameter)
//! and BlendFunc.cxx (L40-280: GetShape, GetMinimalWeights, NextShape,
//! ComputeNormal, ComputeDNormal).

use glam::DVec3;

use rcad_kernel::geom::{Circle3, Curve3, Surface3, SurfaceEval as _, TrimmedCurve3};
use rcad_kernel::math::cs_lib::{
    dnnormal, dnnuv_array2, normal_max_order, NormalStatus,
};
use rcad_kernel::math::GeomAbsShape;

/// OCCT Convert_ParameterisationType (TKMath/Convert,
/// Convert_ParameterisationType.hxx L79-89).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConvertParameterisationType {
    TgtThetaOver2,
    TgtThetaOver2_1,
    TgtThetaOver2_2,
    TgtThetaOver2_3,
    TgtThetaOver2_4,
    QuasiAngular,
    RationalC1,
    Polynomial,
}

/// OCCT BlendFunc_SectionShape (BlendFunc_SectionShape.hxx L19-24).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendFuncSectionShape {
    Rational,
    QuasiAngular,
    Polynomial,
    Linear,
}

/// OCCT BlendFunc::GetShape (BlendFunc.cxx L40-88).
pub fn blend_func_get_shape(
    s_shape: BlendFuncSectionShape,
    max_ang: f64,
    nb_poles: &mut i32,
    nb_knots: &mut i32,
    degree: &mut i32,
    t_conv: &mut ConvertParameterisationType,
) {
    match s_shape {
        BlendFuncSectionShape::Rational => {
            let nb_span = (3.0 * max_ang.abs() / 2.0 / std::f64::consts::PI).ceil() as i32;
            *nb_poles = 2 * nb_span + 1;
            *nb_knots = nb_span + 1;
            *degree = 2;
            if nb_span == 1 {
                *t_conv = ConvertParameterisationType::TgtThetaOver2_1;
            } else {
                // QuasiAngular affin d'etre C1 (et meme beaucoup plus)
                *nb_poles = 7;
                *nb_knots = 2;
                *degree = 6;
                *t_conv = ConvertParameterisationType::QuasiAngular;
            }
        }
        BlendFuncSectionShape::QuasiAngular => {
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
            *t_conv = ConvertParameterisationType::QuasiAngular;
        }
        BlendFuncSectionShape::Polynomial => {
            *nb_poles = 8;
            *nb_knots = 2;
            *degree = 7;
            *t_conv = ConvertParameterisationType::Polynomial;
        }
        BlendFuncSectionShape::Linear => {
            *nb_poles = 2;
            *nb_knots = 2;
            *degree = 1;
        }
    }
}

/// OCCT BlendFunc::GetMinimalWeights (BlendFunc.cxx L96-134) — on suppose
/// les extremum de poids sont obtenus pour les extremums d'angles.
pub fn blend_func_get_minimal_weights(
    s_shape: BlendFuncSectionShape,
    t_conv: ConvertParameterisationType,
    min_ang: f64,
    max_ang: f64,
    weights: &mut [f64],
) {
    use rcad_kernel::base::convert::{
        geom_convert_curve_to_bspline_curve, ConvertParameterisation,
    };
    match s_shape {
        BlendFuncSectionShape::Polynomial | BlendFuncSectionShape::Linear => {
            for w in weights.iter_mut() {
                *w = 1.0;
            }
        }
        BlendFuncSectionShape::Rational | BlendFuncSectionShape::QuasiAngular => {
            // OCCT L110-131: builds Geom_TrimmedCurve(Geom_Circle, 0, MaxAng),
            // converts with GeomConvert::CurveToBSplineCurve(Sect, TConv) and
            // takes the pole-wise minimum against the MinAng conversion.
            let tconv = match t_conv {
                ConvertParameterisationType::TgtThetaOver2 => ConvertParameterisation::TgtThetaOver2,
                ConvertParameterisationType::TgtThetaOver2_1 => {
                    ConvertParameterisation::TgtThetaOver2_1
                }
                ConvertParameterisationType::TgtThetaOver2_2 => {
                    ConvertParameterisation::TgtThetaOver2_2
                }
                ConvertParameterisationType::TgtThetaOver2_3 => {
                    ConvertParameterisation::TgtThetaOver2_3
                }
                ConvertParameterisationType::TgtThetaOver2_4 => {
                    ConvertParameterisation::TgtThetaOver2_4
                }
                ConvertParameterisationType::QuasiAngular => ConvertParameterisation::QuasiAngular,
                ConvertParameterisationType::RationalC1 => ConvertParameterisation::RationalC1,
                ConvertParameterisationType::Polynomial => ConvertParameterisation::Polynomial,
            };
            // OCCT: gp_Circ C(gp_Ax2(Origin, Z), 1).
            let circle = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
            let sect1 = Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(Curve3::Circle(circle)),
                first: 0.0,
                last: max_ang,
            });
            let ctobspl = geom_convert_curve_to_bspline_curve(&sect1, tconv);
            // OCCT: Weights.Assign(CtoBspl->WeightsArray()).
            for (w, sw) in weights.iter_mut().zip(ctobspl.weights.iter()) {
                *w = *sw;
            }

            // OCCT: angle_min = max(Precision::PConfusion(), MinAng).
            let angle_min = 1e-9f64.max(min_ang);
            let sect2 = Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(Curve3::Circle(circle)),
                first: 0.0,
                last: angle_min,
            });
            let ctobspl = geom_convert_curve_to_bspline_curve(&sect2, tconv);
            for (w, poids) in weights.iter_mut().zip(ctobspl.weights.iter()) {
                if *poids < *w {
                    *w = *poids;
                }
            }
        }
    }
}

/// OCCT BlendFunc::NextShape (BlendFunc.cxx L138-152) — the next level of
/// continuity.
pub fn blend_func_next_shape(s: GeomAbsShape) -> GeomAbsShape {
    match s {
        GeomAbsShape::C0 => GeomAbsShape::C1,
        GeomAbsShape::C1 => GeomAbsShape::C2,
        GeomAbsShape::C2 => GeomAbsShape::C3,
        _ => GeomAbsShape::CN,
    }
}

/// OCCT BlendFunc::ComputeNormal(Surf, p2d, Normal) (BlendFunc.cxx L156-213).
pub fn blend_func_compute_normal(
    surf: &Surface3,
    p2d: (f64, f64),
    normal: &mut glam::DVec3,
) -> bool {
    let max_order = 3i32;
    let u = p2d.0;
    let v = p2d.1;

    // OCCT: NCollection_Array2<gp_Vec> DerSurf(0, MaxOrder + 1, 0, MaxOrder + 1).
    let mut der_surf = vec![vec![glam::DVec3::ZERO; (max_order + 2) as usize]; (max_order + 2) as usize];
    for i in 1..=max_order + 1 {
        der_surf[i as usize][0] = surf.dn(u, v, i, 0);
    }
    for i in 0..=max_order + 1 {
        for j in 1..=max_order + 1 {
            der_surf[i as usize][j as usize] = surf.dn(u, v, i, j);
        }
    }

    // OCCT: NCollection_Array2<gp_Vec> DerNUV(0, MaxOrder, 0, MaxOrder).
    let mut der_nuv = vec![vec![glam::DVec3::ZERO; (max_order + 1) as usize]; (max_order + 1) as usize];
    for i in 0..=max_order {
        for j in 0..=max_order {
            der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, &der_surf);
        }
    }

    let (umin, umax, vmin, vmax) = {
        let domain = surf.default_domain();
        (domain[0], domain[1], domain[2], domain[3])
    };
    let (stat, thenormal, _order_u, _order_v) = normal_max_order(
        max_order,
        &der_nuv,
        1e-9,
        u,
        v,
        umin,
        umax,
        vmin,
        vmax,
    );
    if stat == NormalStatus::Defined {
        *normal = thenormal.expect("CSLib_Defined without a direction");
        true
    } else {
        false
    }
}

/// OCCT BlendFunc::ComputeDNormal(Surf, p2d, Normal, DNu, DNv)
/// (BlendFunc.cxx L216-280).
#[allow(clippy::too_many_arguments)]
pub fn blend_func_compute_dnormal(
    surf: &Surface3,
    p2d: (f64, f64),
    normal: &mut glam::DVec3,
    d_nu: &mut glam::DVec3,
    d_nv: &mut glam::DVec3,
) -> bool {
    let max_order = 3i32;
    let u = p2d.0;
    let v = p2d.1;

    let mut der_surf = vec![vec![glam::DVec3::ZERO; (max_order + 2) as usize]; (max_order + 2) as usize];
    for i in 1..=max_order + 1 {
        der_surf[i as usize][0] = surf.dn(u, v, i, 0);
    }
    for i in 0..=max_order + 1 {
        for j in 1..=max_order + 1 {
            der_surf[i as usize][j as usize] = surf.dn(u, v, i, j);
        }
    }

    let mut der_nuv = vec![vec![glam::DVec3::ZERO; (max_order + 1) as usize]; (max_order + 1) as usize];
    for i in 0..=max_order {
        for j in 0..=max_order {
            der_nuv[i as usize][j as usize] = dnnuv_array2(i, j, &der_surf);
        }
    }

    let (umin, umax, vmin, vmax) = {
        let domain = surf.default_domain();
        (domain[0], domain[1], domain[2], domain[3])
    };
    let (stat, thenormal, order_u, order_v) = normal_max_order(
        max_order,
        &der_nuv,
        1e-9,
        u,
        v,
        umin,
        umax,
        vmin,
        vmax,
    );
    if stat == NormalStatus::Defined {
        *normal = thenormal.expect("CSLib_Defined without a direction");
        *d_nu = dnnormal(1, 0, &der_nuv, order_u, order_v);
        *d_nv = dnnormal(0, 1, &der_nuv, order_u, order_v);
        true
    } else {
        false
    }
}
