//! OCCT BlendFunc package services (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_SectionShape.hxx (L19-24), Convert_ParameterisationType.hxx
//! (TKMath/Convert, L79-89, required by BlendFunc::GetShape's out parameter)
//! and BlendFunc.cxx (L40-152: GetShape, GetMinimalWeights, NextShape).
//!
//! Pending (kernel gaps, not translated here):
//! - BlendFunc::ComputeNormal / ComputeDNormal (BlendFunc.cxx L156-280)
//!   need `Adaptor3d_Surface::DN(U, V, Nu, Nv)` (arbitrary-order
//!   derivatives) and the `CSLib::DNNUV` / `CSLib::Normal` Array2 forms;
//!   rcad `SurfaceEval` / `cs_lib` do not expose them yet.
//! - The `Rational`/`QuasiAngular` branch of GetMinimalWeights needs
//!   `GeomConvert::CurveToBSplineCurve` (Convert package).

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
/// The `Rational`/`QuasiAngular` branch is pending GeomConvert (see the
/// module header); it is unreachable until the constant-radius function
/// family (second batch) lands.
pub fn blend_func_get_minimal_weights(
    s_shape: BlendFuncSectionShape,
    _t_conv: ConvertParameterisationType,
    _min_ang: f64,
    _max_ang: f64,
    weights: &mut [f64],
) {
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
            // Pending: GeomConvert::CurveToBSplineCurve is not translated.
            unimplemented!("BlendFunc::GetMinimalWeights Rational/QuasiAngular: pending GeomConvert");
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
