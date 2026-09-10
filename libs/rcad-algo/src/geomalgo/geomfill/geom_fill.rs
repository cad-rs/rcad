//! OCCT GeomFill class statics (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill.hxx (L41-127 declarations) + GeomFill.cxx (whole file L42-827):
//! Surface, GetShape, GetMinimalWeights, Knots, Mults, GetTolerance and the
//! three GetCircle overloads (0th / 1st / 2nd order).
//!
//! GAP notes:
//! - `GeomFill::Surface`'s fallback reuses the already-translated
//!   GeomFill_Generator port hosted in `crate::brep_fill::generator`
//!   (BRepFill_Generator.cxx carries it; relocation pending).
//! - GetMinimalWeights / GetTolerance use the rcad-kernel
//!   GeomConvert::CurveToBSplineCurve port
//!   (`rcad_kernel::base::convert::geom_convert_curve_to_bspline_curve`).

use glam::DVec3;

use rcad_kernel::base::convert::{
    geom_convert_curve_to_bspline_curve, ConvertParameterisation,
};
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve3, CylindricalSurface, Line3, Plane, Surface3, TrimmedCurve3,
    TrimmedSurface,
};
use rcad_kernel::math::gp::{Ax2, Ax3};

use crate::brep_fill::generator::GeomFillGenerator;

use super::polynomial_convertor::PolynomialConvertor;
use super::quasi_angular_convertor::QuasiAngularConvertor;

/// OCCT GeomFill::Surface (GeomFill.cxx L42-182) — builds a ruled surface
/// between the two curves. Special-cased outputs: plane (2 trimmed lines),
/// cylinder / cone (2 coaxial circles), otherwise the GeomFill_Generator
/// ruled BSpline surface.
pub fn surface(curve1: &Curve3, curve2: &Curve3) -> Surface3 {
    // OCCT: occ::handle<Geom_Curve> TheCurve1, TheCurve2; the trimmed-case
    // basis curve extraction maps to cloned Curve3 handles.
    let the_curve1: Curve3;
    let the_curve2: Curve3;
    let mut a1 = 0.0f64;
    let mut a2 = 0.0f64;
    let mut b1 = 0.0f64;
    let mut b2 = 0.0f64;
    let mut trim1 = false;
    let mut trim2 = false;

    // recherche du type de la surface resultat (trimmed curve unwrapping).
    match curve1 {
        Curve3::Trimmed(ctrim) => {
            the_curve1 = ctrim.basis_curve().clone();
            a1 = ctrim.first;
            b1 = ctrim.last;
            trim1 = true;
        }
        _ => {
            the_curve1 = curve1.clone();
        }
    }
    match curve2 {
        Curve3::Trimmed(ctrim) => {
            the_curve2 = ctrim.basis_curve().clone();
            a2 = ctrim.first;
            b2 = ctrim.last;
            trim2 = true;
        }
        _ => {
            the_curve2 = curve2.clone();
        }
    }

    let mut surf: Option<Surface3> = None;
    let mut is_done = false;

    let is_line1 = matches!(the_curve1, Curve3::Line(_));
    let is_line2 = matches!(the_curve2, Curve3::Line(_));
    let is_circle1 = matches!(the_curve1, Curve3::Circle(_));
    let is_circle2 = matches!(the_curve2, Curve3::Circle(_));

    // Les deux courbes sont des droites.
    if is_line1 && is_line2 && trim1 && trim2 {
        let line1 = line_of(&the_curve1);
        let line2 = line_of(&the_curve2);
        // OCCT gp_Lin::Location / Direction.
        let loc1 = line1.origin;
        let d1 = line1.direction;
        let d2 = line2.direction;

        if dir_is_parallel(d1, d2, ANGULAR) {
            let p1p2 = loc1 - line2.origin; // OCCT gp_Vec P1P2(L1.Location(), L2.Location())
            let proj = p1p2.dot(d1);

            if dir_is_equal(d1, d2, ANGULAR) {
                if (a1 - proj - a2).abs() <= CONFUSION && (b1 - proj - b2).abs() <= CONFUSION {
                    // OCCT: gp_Ax3 Ax(L1.Location(), gp_Dir(D1.Crossed(P1P2)), D1).
                    let ax = Ax3::from_pnt_n_vx(loc1, d1.cross(p1p2), d1);
                    // OCCT: new Geom_Plane(Ax) — keeps the Ax3 frame.
                    let p = Plane {
                        origin: ax.location(),
                        normal: ax.direction(),
                        u_dir: ax.x_direction,
                        v_dir: ax.y_direction,
                    };
                    // OCCT: double V = P1P2.Dot(Ax.YDirection()).
                    let v = p1p2.dot(ax.y_direction);
                    surf = Some(Surface3::Trimmed(TrimmedSurface::new(
                        Surface3::Plane(p),
                        a1,
                        b1,
                        0.0_f64.min(v),
                        0.0_f64.max(v),
                    )));
                    is_done = true;
                }
            }
            if dir_is_opposite(d1, d2, ANGULAR) {
                if (a1 - proj + b2).abs() <= CONFUSION && (b1 - proj + a2).abs() <= CONFUSION {
                    let ax = Ax3::from_pnt_n_vx(loc1, d1.cross(p1p2), d1);
                    let p = Plane {
                        origin: ax.location(),
                        normal: ax.direction(),
                        u_dir: ax.x_direction,
                        v_dir: ax.y_direction,
                    };
                    let v = p1p2.dot(ax.y_direction);
                    surf = Some(Surface3::Trimmed(TrimmedSurface::new(
                        Surface3::Plane(p),
                        a1,
                        b1,
                        0.0_f64.min(v),
                        0.0_f64.max(v),
                    )));
                    is_done = true;
                }
            }
        }
    }
    // Les deux courbes sont des cercles.
    else if is_circle1 && is_circle2 {
        let circ1 = circle_of(&the_curve1);
        let circ2 = circle_of(&the_curve2);

        // OCCT: gp_Ax3 A1 = C1.Position(); gp_Ax3 A2 = C2.Position().
        let a1_ax3 = Ax3::from_ax2(&circle_ax2(&circ1));
        let a2_ax3 = Ax3::from_ax2(&circle_ax2(&circ2));

        // first, A1 & A2 must be coaxials
        if ax1_is_coaxial(a1_ax3.axis, a2_ax3.axis, ANGULAR, CONFUSION) {
            let mut v = (a2_ax3.location() - a1_ax3.location()).dot(a1_ax3.direction());
            if !trim1 && !trim2 {
                if (circ1.radius - circ2.radius).abs() < CONFUSION {
                    // OCCT: new Geom_CylindricalSurface(A1, C1.Radius()).
                    let c = CylindricalSurface {
                        origin: a1_ax3.location(),
                        axis: a1_ax3.direction(),
                        radius: circ1.radius,
                        ref_dir: a1_ax3.x_direction,
                        y_dir: Some(a1_ax3.y_direction),
                    };
                    // OCCT: Geom_RectangularTrimmedSurface(C, min(0,V),
                    // max(0,V), false) — UTrim = false: the [min,max] pair
                    // bounds V, U stays untrimmed (rcad's 4-field trim model
                    // cannot express an untrimmed direction; U keeps 0..0).
                    let c = Surface3::Cylinder(c);
                    surf = Some(Surface3::Trimmed(TrimmedSurface::new(
                        c,
                        0.0,
                        0.0,
                        0.0_f64.min(v),
                        0.0_f64.max(v),
                    )));
                } else {
                    let rad = circ2.radius - circ1.radius;
                    let mut ang = (rad / v).atan();
                    let mut axis1 = a1_ax3.direction();
                    if ang < 0.0 {
                        // OCCT: A1.ZReverse() — reverses ONLY the main
                        // direction (gp_Ax3.hxx L131).
                        axis1 = -axis1;
                        v = -v;
                        ang = -ang;
                    }
                    // OCCT: new Geom_ConicalSurface(A1, Ang, C1.Radius()).
                    let c = ConicalSurface::new_with_ref_dir(
                        a1_ax3.location(),
                        axis1,
                        circ1.radius,
                        ang,
                        a1_ax3.x_direction,
                    );
                    v /= ang.cos();
                    // OCCT: Geom_RectangularTrimmedSurface(C, min(0,V),
                    // max(0,V), false) — same UTrim mapping note as above.
                    let c = Surface3::Cone(c);
                    surf = Some(Surface3::Trimmed(TrimmedSurface::new(
                        c,
                        0.0,
                        0.0,
                        0.0_f64.min(v),
                        0.0_f64.max(v),
                    )));
                }
                is_done = true;
            } else if trim1 && trim2 {
                // OCCT L166-168: empty branch.
            }
        }
    }

    // OCCT L172-179: GeomFill_Generator fallback.
    if !is_done {
        let mut generator = GeomFillGenerator::new();
        generator.add_curve(curve1);
        generator.add_curve(curve2);
        generator.perform(PCONFUSION);
        if let Some(bspl) = generator.surface() {
            surf = Some(Surface3::BSpline(bspl));
        }
    }

    surf.expect("GeomFill::Surface: generator produced no surface")
}

/// OCCT GeomFill::GetShape (GeomFill.cxx L186-225).
pub fn get_shape(
    max_ang: f64,
    nb_poles: &mut i32,
    nb_knots: &mut i32,
    degree: &mut i32,
    t_conv: &mut ConvertParameterisation,
) {
    match t_conv {
        ConvertParameterisation::QuasiAngular => {
            *nb_poles = 7;
            *nb_knots = 2;
            *degree = 6;
        }
        ConvertParameterisation::Polynomial => {
            *nb_poles = 8;
            *nb_knots = 2;
            *degree = 7;
        }
        _ => {
            let nb_span = (3.0 * max_ang.abs() / 2.0 / std::f64::consts::PI).ceil() as i32;
            *nb_poles = 2 * nb_span + 1;
            *nb_knots = nb_span + 1;
            *degree = 2;
            if nb_span == 1 {
                *t_conv = ConvertParameterisation::TgtThetaOver2_1;
            } else if nb_span == 2 {
                *t_conv = ConvertParameterisation::TgtThetaOver2_2;
            } else if nb_span == 3 {
                *t_conv = ConvertParameterisation::TgtThetaOver2_3;
            }
        }
    }
}

/// OCCT GeomFill::GetMinimalWeights (GeomFill.cxx L228-264) — on suppose les
/// extremum de poids sont obtenus pour les extremums d'angles.
pub fn get_minimal_weights(
    t_conv: ConvertParameterisation,
    min_ang: f64,
    max_ang: f64,
    weights: &mut [f64],
) {
    if t_conv == ConvertParameterisation::Polynomial {
        // OCCT: Weights.Init(1).
        for w in weights.iter_mut() {
            *w = 1.0;
        }
    } else {
        // OCCT: gp_Ax2 popAx2(gp_Pnt(0,0,0), gp_Dir(Z)); gp_Circ C(popAx2, 1).
        let pop_circ = Circle3::new(DVec3::ZERO, DVec3::Z, 1.0);
        // OCCT: Sect1 = new Geom_TrimmedCurve(new Geom_Circle(C), 0., MaxAng).
        let sect1 = Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Circle(pop_circ),
            0.0,
            max_ang,
        ));
        // OCCT: CtoBspl = GeomConvert::CurveToBSplineCurve(Sect1, TConv).
        let ctobspl = geom_convert_curve_to_bspline_curve(&sect1, t_conv);
        // OCCT: Weights.Assign(CtoBspl->WeightsArray()).
        let arr = &ctobspl.weights;
        for (i, w) in weights.iter_mut().enumerate() {
            *w = arr[i];
        }

        // OCCT: angle_min = max(PConfusion, MinAng); Sect2 = trimmed(0,
        // angle_min); CtoBspl = convert; poids = WeightsArray().
        let angle_min = PCONFUSION.max(min_ang);
        let sect2 = Curve3::Trimmed(TrimmedCurve3::new(
            Curve3::Circle(pop_circ),
            0.0,
            angle_min,
        ));
        let ctobspl = geom_convert_curve_to_bspline_curve(&sect2, t_conv);
        let poids = &ctobspl.weights;

        // OCCT: for ii: if (poids(ii) < Weights(ii)) Weights(ii) = poids(ii).
        for ii in 0..weights.len() {
            if poids[ii] < weights[ii] {
                weights[ii] = poids[ii];
            }
        }
    }
}

/// OCCT GeomFill::Knots (GeomFill.cxx L268-285).
pub fn knots(t_conv: ConvertParameterisation, t_knots: &mut [f64]) {
    if (t_conv != ConvertParameterisation::QuasiAngular)
        && (t_conv != ConvertParameterisation::Polynomial)
    {
        // Cas rational classique (flat 0..n knots).
        let mut val = 0.0f64;
        for k in t_knots.iter_mut() {
            *k = val;
            val += 1.0;
        }
    } else {
        // OCCT: TKnots(1) = 0.; TKnots(2) = 1.
        t_knots[0] = 0.0;
        t_knots[1] = 1.0;
    }
}

/// OCCT GeomFill::Mults (GeomFill.cxx L289-315).
pub fn mults(t_conv: ConvertParameterisation, t_mults: &mut [i32]) {
    match t_conv {
        ConvertParameterisation::QuasiAngular => {
            // OCCT: TMults(1) = 7; TMults(2) = 7.
            t_mults[0] = 7;
            t_mults[1] = 7;
        }
        ConvertParameterisation::Polynomial => {
            // OCCT: TMults(1) = 8; TMults(2) = 8.
            t_mults[0] = 8;
            t_mults[1] = 8;
        }
        _ => {
            // Cas rational classsique: TMults(Lower)=3, inner = 2,
            // TMults(Upper)=3.
            let low = 0usize;
            let upp = t_mults.len() - 1;
            t_mults[low] = 3;
            for i in (low + 1)..upp {
                t_mults[i] = 2;
            }
            t_mults[upp] = 3;
        }
    }
}

/// OCCT GeomFill::GetTolerance (GeomFill.cxx L318-340) — determiner la
/// tolerance 3d permetant de respecter la tolerance de continuite G1.
pub fn get_tolerance(
    t_conv: ConvertParameterisation,
    angle_min: f64,
    radius: f64,
    angular_tol: f64,
    spatial_tol: f64,
) -> f64 {
    // OCCT: gp_Circ C(popAx2, Radius); Sect = trimmed(0, max(AngleMin, 0.02))
    // — 0.02 est proche d'1 degree, en desous on ne se preocupe pas de la
    // tngence afin d'eviter des tolerances d'approximation tendant vers 0 !
    let pop_circ = Circle3::new(DVec3::ZERO, DVec3::Z, radius);
    let sect = Curve3::Trimmed(TrimmedCurve3::new(
        Curve3::Circle(pop_circ),
        0.0,
        angle_min.max(0.02),
    ));
    let ctobspl = geom_convert_curve_to_bspline_curve(&sect, t_conv);
    // OCCT: Dist = CtoBspl->Pole(1).Distance(CtoBspl->Pole(2)) + SpatialTol.
    let dist = ctobspl.control_points[0].distance(ctobspl.control_points[1]) + spatial_tol;
    dist * angular_tol / 2.0
}

/// OCCT GeomFill::GetCircle (GeomFill.cxx L342-445) — calculs les poles et
/// poids d'un cercle definie par ses extremites et son rayon. On evite (si
/// possible) de passer par les convertions pour 1) des problemes de
/// performances, 2) assurer la coherance entre cette methode est celle qui
/// donne la derive.
///
/// `ns1`/`ns2`: normal rentrente au premier / second point; `nplan`: normal
/// au plan; `rayon`: rayon (doit etre positif).
///
/// Array mapping: OCCT `NCollection_Array1` is 1-based (`Low..Upp`); rcad
/// slices are 0-based (`0..n-1`), so `Poles(jj)` maps to `poles[jj - 1]`.
#[allow(clippy::too_many_arguments)]
pub fn get_circle(
    t_conv: ConvertParameterisation,
    ns1: DVec3,
    ns2: DVec3,
    nplan: DVec3,
    pts1: DVec3,
    pts2: DVec3,
    rayon: f64,
    center: DVec3,
    poles: &mut [DVec3],
    weights: &mut [f64],
) {
    // La classe de convertion
    let mut cosa: f64;
    let mut sina: f64;
    let low = 0usize; // OCCT: int low = Poles.Lower().
    let upp = poles.len() - 1; // OCCT: int upp = Poles.Upper().

    cosa = ns1.dot(ns2);
    sina = nplan.dot(ns1.cross(ns2));

    if cosa < -1.0 {
        cosa = -1.0;
        sina = 0.0;
    }
    if cosa > 1.0 {
        cosa = 1.0;
        sina = 0.0;
    }
    let mut angle = cosa.acos();
    // Recadrage sur ]-pi/2, 3pi/2]
    if sina < 0.0 {
        if cosa > 0.0 {
            angle = -angle;
        } else {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
    }

    match t_conv {
        ConvertParameterisation::QuasiAngular => {
            let mut q_convertor = QuasiAngularConvertor::new();
            q_convertor.init();
            q_convertor.section(pts1, center, nplan, angle, poles, weights);
        }
        ConvertParameterisation::Polynomial => {
            let mut p_convertor = PolynomialConvertor::new();
            p_convertor.init();
            p_convertor.section(pts1, center, nplan, angle, poles);
            for w in weights.iter_mut() {
                *w = 1.0;
            }
        }
        _ => {
            // Cas Rational, on utilise une expression directe beaucoup plus
            // performente que GeomConvert
            let nb_span = (poles.len() as i32 - 1) / 2;

            poles[low] = pts1;
            poles[upp] = pts2;
            weights[low] = 1.0;
            weights[upp] = 1.0;

            // OCCT: np2 = nplan.Crossed(ns1).
            let np2 = nplan.cross(ns1);

            let alpha = angle / (nb_span as f64);
            let cosas2 = (alpha / 2.0).cos();

            // OCCT: for (i = 1, jj = low + 2; i <= NbSpan - 1; i++, jj += 2).
            let mut jj = (low + 2) as i32;
            for i in 1..=nb_span - 1 {
                let lambda = (i as f64) * alpha;
                cosa = lambda.cos();
                sina = lambda.sin();
                // OCCT: temp.SetLinearForm(Cosa - 1, ns1, Sina, np2).
                let temp = (cosa - 1.0) * ns1 + sina * np2;
                poles[(jj - 1) as usize] = pts1 + rayon * temp;
                weights[(jj - 1) as usize] = 1.0;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            // OCCT: for (i = 1, jj = low + 1; i <= NbSpan; i++, jj += 2).
            let mut jj = (low + 1) as i32;
            for _ in 1..=nb_span {
                let temp = poles[(jj - 1) as usize]
                    + poles[(jj + 1) as usize]
                    - 2.0 * center;
                poles[jj as usize] = center + lambda * temp;
                weights[jj as usize] = cosas2;
                jj += 2;
            }
        }
    }
}

/// OCCT GeomFill::GetCircle first-derivative overload (GeomFill.cxx L447-603).
#[allow(clippy::too_many_arguments)]
pub fn get_circle_d1(
    t_conv: ConvertParameterisation,
    ns1: DVec3,
    ns2: DVec3,
    dn1w: DVec3,
    dn2w: DVec3,
    nplan: DVec3,
    dnplan: DVec3,
    pts1: DVec3,
    pts2: DVec3,
    tang1: DVec3,
    tang2: DVec3,
    rayon: f64,
    drayon: f64,
    center: DVec3,
    dcenter: DVec3,
    poles: &mut [DVec3],
    dpoles: &mut [DVec3],
    weights: &mut [f64],
    dweights: &mut [f64],
) -> bool {
    let mut cosa: f64;
    let mut sina: f64;
    let low = 0usize;
    let upp = poles.len() - 1;
    let nb_span = (poles.len() as i32 - 1) / 2;

    cosa = ns1.dot(ns2);
    sina = nplan.dot(ns1.cross(ns2));

    if cosa < -1.0 {
        cosa = -1.0;
        sina = 0.0;
    }
    if cosa > 1.0 {
        cosa = 1.0;
        sina = 0.0;
    }
    let mut angle = cosa.acos();
    // Recadrage sur ]-pi/2, 3pi/2]
    if sina < 0.0 {
        if cosa > 0.0 {
            angle = -angle;
        } else {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
    }

    let dangle: f64;
    if sina.abs() > cosa.abs() {
        dangle = -(dn1w.dot(ns2) + ns1.dot(dn2w)) / sina;
    } else {
        dangle = (dnplan.dot(ns1.cross(ns2)) + nplan.dot(dn1w.cross(ns2) + ns1.cross(dn2w))) / cosa;
    }

    // Aux Extremites.
    poles[low] = pts1;
    poles[upp] = pts2;
    weights[low] = 1.0;
    weights[upp] = 1.0;

    dpoles[low] = tang1;
    dpoles[upp] = tang2;
    dweights[low] = 0.0;
    dweights[upp] = 0.0;

    match t_conv {
        ConvertParameterisation::QuasiAngular => {
            let mut q_convertor = QuasiAngularConvertor::new();
            q_convertor.init();
            q_convertor.section_d1(
                pts1,
                tang1,
                center,
                dcenter,
                nplan,
                dnplan,
                angle,
                dangle,
                poles,
                dpoles,
                weights,
                dweights,
            );
            return true;
        }
        ConvertParameterisation::Polynomial => {
            let mut p_convertor = PolynomialConvertor::new();
            p_convertor.init();
            p_convertor.section_d1(
                pts1,
                tang1,
                center,
                dcenter,
                nplan,
                dnplan,
                angle,
                dangle,
                poles,
                dpoles,
            );
            for w in weights.iter_mut() {
                *w = 1.0;
            }
            for w in dweights.iter_mut() {
                *w = 0.0;
            }
            return true;
        }
        _ => {
            // Cas rationel classique
            // OCCT: np2 = nplan.Crossed(ns1);
            let np2 = nplan.cross(ns1);
            // OCCT: dnp2 = dnplan.Crossed(ns1).Added(nplan.Crossed(dn1w)).
            let dnp2 = dnplan.cross(ns1) + nplan.cross(dn1w);

            let alpha = angle / (nb_span as f64);
            let cosas2 = (alpha / 2.0).cos();
            let sinas2 = (alpha / 2.0).sin();

            // OCCT: for (i = 1, jj = low + 2; i <= NbSpan - 1; i++, jj += 2).
            let mut jj = (low + 2) as i32;
            for i in 1..=nb_span - 1 {
                let lambda = (i as f64) * alpha;
                cosa = lambda.cos();
                sina = lambda.sin();
                // OCCT: temp.SetLinearForm(Cosa - 1, ns1, Sina, np2).
                let temp = (cosa - 1.0) * ns1 + sina * np2;
                poles[(jj - 1) as usize] = pts1 + rayon * temp;

                // OCCT: DPoles(jj).SetLinearForm(DRayon, temp, tang1).
                let mut dpol = drayon * temp + tang1;
                // OCCT: temp.SetLinearForm(-Sina, ns1, Cosa, np2);
                // temp.Multiply(i/NbSpan * DAngle);
                let mut temp = -sina * ns1 + cosa * np2;
                temp *= (i as f64) / (nb_span as f64) * dangle;
                // OCCT: temp.Add(((Cosa - 1) * dn1w).Added(Sina * dnp2)).
                temp += (cosa - 1.0) * dn1w + sina * dnp2;
                // OCCT: DPoles(jj) += Rayon * temp.
                dpol += rayon * temp;
                dpoles[(jj - 1) as usize] = dpol;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            let mut dlambda = (lambda * sinas2 * dangle) / (cosas2 * nb_span as f64);

            // OCCT: for (i = 1, jj = low; i <= NbSpan; i++, jj += 2).
            let mut jj = low as i32;
            for _ in 1..=nb_span {
                // OCCT: temp.SetXYZ(Poles(jj).XYZ() + Poles(jj + 2).XYZ()
                // - 2. * Center.XYZ()).
                let temp = poles[jj as usize] + poles[(jj + 2) as usize] - 2.0 * center;
                poles[(jj + 1) as usize] = center + lambda * temp;
                // OCCT: DPoles(jj + 1).SetLinearForm(Dlambda, temp,
                // 1. - 2 * lambda, DCenter, lambda, (DPoles(jj)+DPoles(jj+2))).
                dpoles[(jj + 1) as usize] = dlambda * temp
                    + (1.0 - 2.0 * lambda) * dcenter
                    + lambda * (dpoles[jj as usize] + dpoles[(jj + 2) as usize]);
                jj += 2;
            }

            // Les poids
            // OCCT: Dlambda = -Sinas2 * DAngle / (2 * NbSpan).
            dlambda = -sinas2 * dangle / (2.0 * nb_span as f64);
            let mut i = low as i32;
            while i < upp as i32 {
                weights[i as usize] = 1.0;
                weights[(i + 1) as usize] = cosas2;
                dweights[i as usize] = 0.0;
                dweights[(i + 1) as usize] = dlambda;
                i += 2;
            }
            true
        }
    }
}

/// OCCT GeomFill::GetCircle second-derivative overload (GeomFill.cxx L605-827).
#[allow(clippy::too_many_arguments)]
pub fn get_circle_d2(
    t_conv: ConvertParameterisation,
    ns1: DVec3,
    ns2: DVec3,
    dn1w: DVec3,
    dn2w: DVec3,
    d2n1w: DVec3,
    d2n2w: DVec3,
    nplan: DVec3,
    dnplan: DVec3,
    d2nplan: DVec3,
    pts1: DVec3,
    pts2: DVec3,
    tang1: DVec3,
    tang2: DVec3,
    dtang1: DVec3,
    dtang2: DVec3,
    rayon: f64,
    drayon: f64,
    d2rayon: f64,
    center: DVec3,
    dcenter: DVec3,
    d2center: DVec3,
    poles: &mut [DVec3],
    dpoles: &mut [DVec3],
    d2poles: &mut [DVec3],
    weights: &mut [f64],
    dweights: &mut [f64],
    d2weights: &mut [f64],
) -> bool {
    let mut cosa: f64;
    let mut sina: f64;
    let aux: f64;
    let low = 0usize;
    let upp = poles.len() - 1;
    let nb_span = (poles.len() as i32 - 1) / 2;

    cosa = ns1.dot(ns2);
    sina = nplan.dot(ns1.cross(ns2));

    if cosa < -1.0 {
        cosa = -1.0;
        sina = 0.0;
    }
    if cosa > 1.0 {
        cosa = 1.0;
        sina = 0.0;
    }
    let mut angle = cosa.acos();
    // Recadrage sur ]-pi/2, 3pi/2]
    if sina < 0.0 {
        if cosa > 0.0 {
            angle = -angle;
        } else {
            angle = 2.0 * std::f64::consts::PI - angle;
        }
    }

    let dangle: f64;
    let d2angle: f64;
    if sina.abs() > cosa.abs() {
        aux = dn1w.dot(ns2) + ns1.dot(dn2w);
        dangle = -aux / sina;
        d2angle = -(d2n1w.dot(ns2) + 2.0 * dn1w.dot(dn2w) + ns1.dot(d2n2w)) / sina
            + aux
                * (dnplan.dot(ns1.cross(ns2)) + nplan.dot(dn1w.cross(ns2) + ns1.cross(dn2w)))
                / (sina * sina);
    } else {
        // OCCT: temp = dn1w.Crossed(ns2) + ns1.Crossed(dn2w).
        let temp = dn1w.cross(ns2) + ns1.cross(dn2w);
        dangle = (dnplan.dot(ns1.cross(ns2)) + nplan.dot(temp)) / cosa;
        d2angle = (d2nplan.dot(ns1.cross(ns2))
            + 2.0 * dnplan.dot(temp)
            + nplan.dot(d2n1w.cross(ns2) + 2.0 * dn1w.cross(dn2w) + ns1.cross(d2n2w)))
            / cosa
            - (dn1w.dot(ns2) + ns1.dot(dn2w)) * (dnplan.dot(ns1.cross(ns2)) + nplan.dot(temp))
                / (cosa * cosa);
    }

    // Aux Extremites.
    poles[low] = pts1;
    poles[upp] = pts2;
    weights[low] = 1.0;
    weights[upp] = 1.0;

    dpoles[low] = tang1;
    dpoles[upp] = tang2;
    dweights[low] = 0.0;
    dweights[upp] = 0.0;

    d2poles[low] = dtang1;
    d2poles[upp] = dtang2;
    d2weights[low] = 0.0;
    d2weights[upp] = 0.0;

    match t_conv {
        ConvertParameterisation::QuasiAngular => {
            let mut q_convertor = QuasiAngularConvertor::new();
            q_convertor.init();
            q_convertor.section_d2(
                pts1,
                tang1,
                dtang1,
                center,
                dcenter,
                d2center,
                nplan,
                dnplan,
                d2nplan,
                angle,
                dangle,
                d2angle,
                poles,
                dpoles,
                d2poles,
                weights,
                dweights,
                d2weights,
            );
            return true;
        }
        ConvertParameterisation::Polynomial => {
            let mut p_convertor = PolynomialConvertor::new();
            p_convertor.init();
            p_convertor.section_d2(
                pts1,
                tang1,
                dtang1,
                center,
                dcenter,
                d2center,
                nplan,
                dnplan,
                d2nplan,
                angle,
                dangle,
                d2angle,
                poles,
                dpoles,
                d2poles,
            );
            for w in weights.iter_mut() {
                *w = 1.0;
            }
            for w in dweights.iter_mut() {
                *w = 0.0;
            }
            for w in d2weights.iter_mut() {
                *w = 0.0;
            }
            return true;
        }
        _ => {
            // OCCT: np2 = nplan.Crossed(ns1).
            let np2 = nplan.cross(ns1);
            // OCCT: dnp2 = dnplan.Crossed(ns1).Added(nplan.Crossed(dn1w)).
            let dnp2 = dnplan.cross(ns1) + nplan.cross(dn1w);
            // OCCT: d2np2 = d2nplan.Crossed(ns1).Added(nplan.Crossed(dn2w));
            // d2np2 += 2 * dnplan.Crossed(dn1w).
            let mut d2np2 = d2nplan.cross(ns1) + nplan.cross(dn2w);
            d2np2 += 2.0 * dnplan.cross(dn1w);

            let alpha = angle / (nb_span as f64);
            let cosas2 = (alpha / 2.0).cos();
            let sinas2 = (alpha / 2.0).sin();

            // OCCT: for (i = 1, jj = low + 2; i <= NbSpan - 1; i++, jj += 2).
            let mut jj = (low + 2) as i32;
            for i in 1..=nb_span - 1 {
                let lambda = (i as f64) * alpha;
                cosa = lambda.cos();
                sina = lambda.sin();
                // OCCT: temp.SetLinearForm(Cosa - 1, ns1, Sina, np2).
                let temp = (cosa - 1.0) * ns1 + sina * np2;
                poles[(jj - 1) as usize] = pts1 + rayon * temp;

                // OCCT: DPoles(jj).SetLinearForm(DRayon, temp, tang1).
                let mut dpol = drayon * temp + tang1;
                // OCCT: dtemp.SetLinearForm(-Sina, ns1, Cosa, np2);
                // aux = i / NbSpan; dtemp.Multiply(aux * DAngle);
                let mut dtemp = -sina * ns1 + cosa * np2;
                let aux = (i as f64) / (nb_span as f64);
                dtemp *= aux * dangle;
                dtemp += (cosa - 1.0) * dn1w + sina * dnp2;
                dpol += rayon * dtemp;
                dpoles[(jj - 1) as usize] = dpol;

                // OCCT: D2Poles(jj).SetLinearForm(D2Rayon, temp, 2 * DRayon,
                // dtemp, Dtang1).
                let mut d2pol = d2rayon * temp + 2.0 * drayon * dtemp + dtang1;
                // OCCT: temp.SetLinearForm(Cosa - 1, dn2w, Sina, d2np2).
                let mut temp = (cosa - 1.0) * dn2w + sina * d2np2;
                // OCCT: dtemp.SetLinearForm(-Sina, ns1, Cosa, np2).
                dtemp = -sina * ns1 + cosa * np2;
                // OCCT: temp += (aux * aux * D2Angle) * dtemp.
                temp += (aux * aux * d2angle) * dtemp;
                // OCCT: dtemp.SetLinearForm(-Sina, dn1w + np2, Cosa, dnp2,
                // -Cosa, ns1).
                dtemp = -sina * (dn1w + np2) + cosa * dnp2 + (-cosa) * ns1;
                // OCCT: temp += (aux * DAngle) * dtemp.
                temp += (aux * dangle) * dtemp;
                // OCCT: D2Poles(jj) += Rayon * temp.
                d2pol += rayon * temp;
                d2poles[(jj - 1) as usize] = d2pol;
                jj += 2;
            }

            let lambda = 1.0 / (2.0 * cosas2 * cosas2);
            let mut dlambda = (lambda * sinas2 * dangle) / (cosas2 * nb_span as f64);
            let aux = sinas2 / cosas2;
            let d2lambda = (dlambda * aux * dangle
                + d2angle * aux * lambda
                + (1.0 + aux * aux) * (dangle / (2.0 * nb_span as f64)) * dangle * lambda)
                / nb_span as f64;

            // OCCT: for (i = 1, jj = low; i <= NbSpan; i++, jj += 2).
            let mut jj = low as i32;
            for _ in 1..=nb_span {
                // OCCT: temp = Poles(jj) + Poles(jj+2) - 2*Center.
                let temp = poles[jj as usize] + poles[(jj + 2) as usize] - 2.0 * center;
                poles[(jj + 1) as usize] = center + lambda * temp;

                // OCCT: dtemp = DPoles(jj) + DPoles(jj+2) - 2*DCenter.
                let dtemp = dpoles[jj as usize] + dpoles[(jj + 2) as usize] - 2.0 * dcenter;
                // OCCT: DPoles(jj + 1).SetLinearForm(Dlambda, temp, lambda,
                // dtemp, DCenter).
                dpoles[(jj + 1) as usize] = dlambda * temp + lambda * dtemp + dcenter;
                // OCCT: D2Poles(jj + 1).SetLinearForm(D2lambda, temp,
                // 2 * Dlambda, dtemp, lambda, (D2Poles(jj) + D2Poles(jj+2)));
                // D2Poles(jj + 1) += (1 - 2 * lambda) * D2Center.
                let mut d2pol = d2lambda * temp
                    + 2.0 * dlambda * dtemp
                    + lambda * (d2poles[jj as usize] + d2poles[(jj + 2) as usize]);
                d2pol += (1.0 - 2.0 * lambda) * d2center;
                d2poles[(jj + 1) as usize] = d2pol;
                jj += 2;
            }

            // Les poids
            // OCCT: Dlambda = -Sinas2 * DAngle / (2 * NbSpan);
            // D2lambda = -Sinas2 * D2Angle / (2 * NbSpan)
            //          - Cosas2 * pow(DAngle / (2 * NbSpan), 2).
            dlambda = -sinas2 * dangle / (2.0 * nb_span as f64);
            let d2lambda_w = -sinas2 * d2angle / (2.0 * nb_span as f64)
                - cosas2 * (dangle / (2.0 * nb_span as f64) * (dangle / (2.0 * nb_span as f64)));
            let mut i = low as i32;
            while i < upp as i32 {
                weights[i as usize] = 1.0;
                weights[(i + 1) as usize] = cosas2;
                dweights[i as usize] = 0.0;
                dweights[(i + 1) as usize] = dlambda;
                d2weights[i as usize] = 0.0;
                d2weights[(i + 1) as usize] = d2lambda_w;
                i += 2;
            }
            true
        }
    }
}

// =============================================================================
// Kernel-mapping helpers (gp module gaps carried here, annotated)
// =============================================================================

/// OCCT Geom_Line::Lin() mapping — rcad stores the line directly.
fn line_of(c: &Curve3) -> Line3 {
    match c {
        Curve3::Line(l) => *l,
        _ => unreachable!("line_of: not a line"),
    }
}

/// OCCT Geom_Circle::Circ() mapping — rcad stores the circle directly.
fn circle_of(c: &Curve3) -> Circle3 {
    match c {
        Curve3::Circle(circ) => *circ,
        _ => unreachable!("circle_of: not a circle"),
    }
}

/// OCCT gp_Circ::Position() (gp_Ax2) mapping.
fn circle_ax2(c: &Circle3) -> Ax2 {
    Ax2::new(c.center, c.normal, c.x_dir)
}

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance) — gp_Dir.lxx:
/// `Us.IsParallel(Vs, AngularTolerance)` is
/// `std::abs(X*Vx + Y*Vy + Z*Vz) >= Cos(AngularTolerance)` (dot of units).
fn dir_is_parallel(d1: DVec3, d2: DVec3, angular_tolerance: f64) -> bool {
    d1.dot(d2).abs() >= angular_tolerance.cos()
}

/// OCCT gp_Dir::IsEqual(Other, AngularTolerance) — dot >= cos(tolerance).
fn dir_is_equal(d1: DVec3, d2: DVec3, angular_tolerance: f64) -> bool {
    d1.dot(d2) >= angular_tolerance.cos()
}

/// OCCT gp_Dir::IsOpposite(Other, AngularTolerance) — dot <= -cos(tolerance).
fn dir_is_opposite(d1: DVec3, d2: DVec3, angular_tolerance: f64) -> bool {
    d1.dot(d2) <= -(angular_tolerance.cos())
}

/// GAP carrier: OCCT gp_Ax1::IsCoaxial (TKMath/gp/gp_Ax1.cxx) — not yet
/// translated in rcad-kernel's gp module; carried here pending relocation.
fn ax1_is_coaxial(
    this: rcad_kernel::math::gp::Ax1,
    other: rcad_kernel::math::gp::Ax1,
    angular_tolerance: f64,
    linear_tolerance: f64,
) -> bool {
    let xyz1 = this.location - other.location;
    let xyz1 = xyz1.cross(other.direction);
    let d1 = xyz1.length();
    let xyz2 = other.location - this.location;
    let xyz2 = xyz2.cross(this.direction);
    let d2 = xyz2.length();
    dir_is_equal(this.direction, other.direction, angular_tolerance)
        && d1 <= linear_tolerance
        && d2 <= linear_tolerance
}
