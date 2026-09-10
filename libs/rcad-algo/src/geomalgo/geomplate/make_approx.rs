//! OCCT GeomPlate_MakeApprox (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_MakeApprox.cxx (whole file L37-558).
//!
//! Complete: the GeomPlate_MakeApprox_Eval evaluator functor (L37-265,
//! served to the approximation driver through the AdvApp2Var_EvaluatorFunc2Var
//! interface), both constructors including the Seq2d/Seq3d constraint
//! collection, the RealBounds/EnlargeCoeff scaling and the seuil clamping
//! (L269-537), and the result reads `mySurface = AppPlate.Surface(1)` /
//! `myAppError = AppPlate.MaxError(3, 1)` / `myCritError =
//! AppPlate.CritError(3, 1)` (L327-329, L451-453, L489-491, L528-530) over
//! the landed AdvApp2Var class layer
//! ([`crate::geomalgo::adv_app2_var::approx_afunc2var::ApproxAFunc2Var`])
//! with the kernel `AdvApprox_DichoCutting` cutting tool.
//!
//! `handle(Geom_BSplineSurface)` -> `Option<BSplineSurface>` (kernel type;
//! None models the null handle).

use glam::{DVec2, DVec3};

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::IsoType;
use rcad_kernel::geom::BSplineSurface;
use rcad_kernel::math::adv_approx::{Cutting, DichoCutting};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::adv_app2_var::approx_afunc2var::ApproxAFunc2Var;
use crate::geomalgo::adv_app2_var::approxf2var_c::EvaluatorFunc2Var;
use crate::geomalgo::adv_app2_var::criterion::Criterion;
use crate::geomalgo::adv_app2_var::nc_array::Array2;

use super::plate_g0_criterion::{
    AdvApp2VarCriterionRepartition, AdvApp2VarCriterionType, PlateG0Criterion,
};
use super::plate_g1_criterion::PlateG1Criterion;
use super::surface::GeomPlateSurface;

// ---------------------------------------------------------------------------
// GeomPlate_MakeApprox_Eval (GeomPlate_MakeApprox.cxx L37-265)
// ---------------------------------------------------------------------------

/// OCCT GeomPlate_MakeApprox_Eval : AdvApp2Var_EvaluatorFunc2Var — the plate
/// surface evaluator served to the approximation driver.
pub struct GeomPlateMakeApproxEval<'a> {
    /// L59: handle(Geom_Surface) mySurf.
    my_surf: &'a GeomPlateSurface,
}

impl<'a> GeomPlateMakeApproxEval<'a> {
    /// OCCT ctor (L41-44).
    pub fn new(surf: &'a GeomPlateSurface) -> Self {
        GeomPlateMakeApproxEval { my_surf: surf }
    }
}

impl EvaluatorFunc2Var for GeomPlateMakeApproxEval<'_> {
    /// OCCT Evaluate (L62-265) — the raw-buffer calling convention of
    /// AdvApp2Var_EvaluatorFunc2Var is preserved (all parameters in OCCT
    /// declaration order; the Result layout is Result[Dimension, N]).
    #[allow(clippy::too_many_arguments)]
    fn evaluate(
        &self,
        dimension: &i32,
        u_start_end: &[f64],
        v_start_end: &[f64],
        favor_iso: &i32,
        const_param: &f64,
        nb_params: &i32,
        parameters: &[f64],
        u_order: &i32,
        v_order: &i32,
        result: &mut [f64],
        error_code: &mut i32,
    ) {
        // *ErrorCode = 0;
        *error_code = 0;
        let mut upar = 0.0f64;
        let mut vpar = 0.0f64;

        // Dimension incorrecte
        if *dimension != 3 {
            *error_code = 1;
        }

        // Parametres incorrects
        if *favor_iso == 1 {
            upar = *const_param;
            if upar < u_start_end[0] || upar > u_start_end[1] {
                *error_code = 2;
            }
            for jpar in 1..=*nb_params {
                vpar = parameters[jpar as usize - 1];
                if vpar < v_start_end[0] || vpar > v_start_end[1] {
                    *error_code = 2;
                }
            }
        } else {
            vpar = *const_param;
            if vpar < v_start_end[0] || vpar > v_start_end[1] {
                *error_code = 2;
            }
            for jpar in 1..=*nb_params {
                upar = parameters[jpar as usize - 1];
                if upar < u_start_end[0] || upar > u_start_end[1] {
                    *error_code = 2;
                }
            }
        }

        // Initialisation
        for idim in 1..=*dimension {
            for jpar in 1..=*nb_params {
                result[(idim - 1 + (jpar - 1) * *dimension) as usize] = 0.0;
            }
        }

        let order = *u_order + *v_order;

        if *favor_iso == 1 {
            upar = *const_param;
            match order {
                0 => {
                    for jpar in 1..=*nb_params {
                        vpar = parameters[jpar as usize - 1];
                        let pnt = self.my_surf.eval_d0(upar, vpar);
                        result[((jpar - 1) * dimension) as usize] = pnt.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pnt.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pnt.z;
                    }
                }
                1 => {
                    for jpar in 1..=*nb_params {
                        vpar = parameters[jpar as usize - 1];
                        let (pnt, v1, v2) = self.my_surf.eval_d1(upar, vpar);
                        let pick = if *u_order == 1 { v1 } else { v2 };
                        let _ = pnt;
                        result[((jpar - 1) * dimension) as usize] = pick.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pick.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pick.z;
                    }
                }
                2 => {
                    for jpar in 1..=*nb_params {
                        vpar = parameters[jpar as usize - 1];
                        let (_p, v1, v2, v3, v4, v5) = self.my_surf.eval_d2(upar, vpar);
                        let pick = if *u_order == 2 {
                            v3
                        } else if *u_order == 1 {
                            v5
                        } else if *u_order == 0 {
                            v4
                        } else {
                            DVec3::ZERO
                        };
                        result[((jpar - 1) * dimension) as usize] = pick.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pick.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pick.z;
                    }
                }
                _ => {}
            }
        } else {
            vpar = *const_param;
            match order {
                0 => {
                    for jpar in 1..=*nb_params {
                        upar = parameters[jpar as usize - 1];
                        let pnt = self.my_surf.eval_d0(upar, vpar);
                        result[((jpar - 1) * dimension) as usize] = pnt.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pnt.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pnt.z;
                    }
                }
                1 => {
                    for jpar in 1..=*nb_params {
                        upar = parameters[jpar as usize - 1];
                        let (pnt, v1, v2) = self.my_surf.eval_d1(upar, vpar);
                        let pick = if *u_order == 1 { v1 } else { v2 };
                        let _ = pnt;
                        result[((jpar - 1) * dimension) as usize] = pick.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pick.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pick.z;
                    }
                }
                2 => {
                    for jpar in 1..=*nb_params {
                        upar = parameters[jpar as usize - 1];
                        let (_p, v1, v2, v3, v4, v5) = self.my_surf.eval_d2(upar, vpar);
                        let pick = if *u_order == 2 {
                            v3
                        } else if *u_order == 1 {
                            v5
                        } else if *u_order == 0 {
                            v4
                        } else {
                            DVec3::ZERO
                        };
                        result[((jpar - 1) * dimension) as usize] = pick.x;
                        result[(1 + (jpar - 1) * dimension) as usize] = pick.y;
                        result[(2 + (jpar - 1) * dimension) as usize] = pick.z;
                    }
                }
                _ => {}
            }
        }
    }
}

// ---------------------------------------------------------------------------
// GeomPlate_MakeApprox (GeomPlate_MakeApprox.cxx L269-558)
// ---------------------------------------------------------------------------

/// OCCT GeomPlate_MakeApprox — the conversion of the plate surface into a
/// BSpline surface.
pub struct MakeApprox {
    /// hxx: handle(GeomPlate_Surface) myPlate.
    my_plate: GeomPlateSurface,
    /// hxx: handle(Geom_BSplineSurface) mySurface.
    my_surface: Option<BSplineSurface>,
    /// hxx: double myAppError.
    my_app_error: f64,
    /// hxx: double myCritError.
    my_crit_error: f64,
}

impl MakeApprox {
    /// OCCT ctor with a criterion (GeomPlate_MakeApprox.cxx L269-335) —
    /// the criterion flavor (PlateCrit carries the AdvApp2Var_Criterion
    /// interface).
    pub fn new_with_criterion(
        surf_plate: &GeomPlateSurface,
        plate_crit: &dyn Criterion,
        tol3d: f64,
        nbmax: i32,
        dgmax: i32,
        continuity: GeomAbsShape,
        enlarge_coeff: f64,
    ) -> Self {
        // myPlate = SurfPlate;
        let mut r = MakeApprox {
            my_plate: surf_plate.clone(),
            my_surface: None,
            my_app_error: 0.0,
            my_crit_error: 0.0,
        };

        // double U0 = 0., U1 = 0., V0 = 0., V1 = 0.;
        // myPlate->RealBounds(U0, U1, V0, V1);           (L279-280)
        // U0 = EnlargeCoeff * U0; ...                    (L281-284)
        let (mut u0, mut u1, mut v0, mut v1) = surf_plate.real_bounds();
        u0 = enlarge_coeff * u0;
        u1 = enlarge_coeff * u1;
        v0 = enlarge_coeff * v0;
        v1 = enlarge_coeff * v1;

        // int nb1 = 0, nb2 = 0, nb3 = 1;                 (L286)
        let nb1 = 0;
        let nb2 = 0;
        let nb3 = 1;
        // nul1 = new NCollection_HArray1<double>(1, 1); nul1->Init(0.);   (L287-288)
        let mut nul1 = vec![0.0f64; 1];
        nul1.fill(0.0);
        // nul2 = new NCollection_HArray2<double>(1, 1, 1, 4); nul2->Init(0.); (L289-290)
        let mut nul2 = Array2::<f64>::new(1, 1, 1, 4);
        nul2.init(0.0);
        // eps3D = new NCollection_HArray1<double>(1, 1); eps3D->Init(Tol3d); (L291-292)
        let mut eps3d = vec![0.0f64; 1];
        eps3d.fill(tol3d);
        // epsfr = new NCollection_HArray2<double>(1, 1, 1, 4); epsfr->Init(Tol3d); (L293-294)
        let mut epsfr = Array2::<f64>::new(1, 1, 1, 4);
        epsfr.init(tol3d);
        // GeomAbs_IsoType myType = GeomAbs_IsoV;         (L295)
        let my_type = IsoType::IsoV;
        // int myPrec = 0;                                (L296)
        let my_prec = 0;

        // AdvApprox_DichoCutting myDec;                  (L298)
        let my_dec = DichoCutting;

        // POP pour WNT
        // GeomPlate_MakeApprox_Eval ev(myPlate);         (L301)
        let ev = GeomPlateMakeApproxEval::new(surf_plate);
        // AdvApp2Var_ApproxAFunc2Var AppPlate(...)       (L302-326)
        let app_plate = ApproxAFunc2Var::new_with_criterion(
            nb1,
            nb2,
            nb3,
            &nul1,
            &nul1,
            &eps3d,
            &nul2,
            &nul2,
            &epsfr,
            u0,
            u1,
            v0,
            v1,
            my_type,
            continuity,
            continuity,
            my_prec,
            dgmax,
            dgmax,
            nbmax,
            &ev,
            plate_crit,
            &my_dec,
            &my_dec,
        );
        // mySurface   = AppPlate.Surface(1);             (L327)
        r.my_surface = app_plate.surface(1).cloned();
        // myAppError  = AppPlate.MaxError(3, 1);         (L328)
        r.my_app_error = app_plate.max_error_at(3, 1);
        // myCritError = AppPlate.CritError(3, 1);        (L329)
        r.my_crit_error = app_plate.crit_error(3, 1);
        r
    }

    /// OCCT ctor with CritOrder (GeomPlate_MakeApprox.cxx L339-537) —
    /// the CritOrder == -1 / 0 / 1 branches; OCCT default arguments:
    /// CritOrder = -1, Continuity = GeomAbs_C0, EnlargeCoeff = 1.3.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surf_plate: &GeomPlateSurface,
        tol3d: f64,
        nbmax: i32,
        dgmax: i32,
        dmax: f64,
        crit_order: i32,
        continuity: GeomAbsShape,
        enlarge_coeff: f64,
    ) -> Self {
        // myPlate = SurfPlate;
        let mut r = MakeApprox {
            my_plate: surf_plate.clone(),
            my_surface: None,
            my_app_error: 0.0,
            my_crit_error: 0.0,
        };

        // NCollection_Sequence<gp_XY>  Seq2d;            (L350)
        let mut seq2d: Vec<DVec2> = Vec::new();
        // NCollection_Sequence<gp_XYZ> Seq3d;            (L351)
        let mut seq3d: Vec<DVec3> = Vec::new();

        // if (CritOrder >= 0)                            (L353)
        if crit_order >= 0 {
            // contraintes 2d d'ordre 0
            // myPlate->Constraints(Seq2d);                (L357)
            surf_plate.constraints(&mut seq2d);

            // contraintes 3d correspondantes sur plate
            // int i, nbp = Seq2d.Length();                (L360)
            let nbp = seq2d.len();
            // for (i = 1; i <= nbp; i++)                  (L361)
            for i in 1..=nbp {
                // gp_XY  P2d = Seq2d.Value(i);
                let p2d = seq2d[i - 1];
                if crit_order == 0 {
                    // a l'ordre 0
                    // myPlate->D0(P2d.X(), P2d.Y(), PP);   (L369)
                    let pp = surf_plate.eval_d0(p2d.x, p2d.y);
                    // gp_XYZ P3d(PP.X(), PP.Y(), PP.Z()); Seq3d.Append(P3d);
                    seq3d.push(DVec3::new(pp.x, pp.y, pp.z));
                } else {
                    // a l'ordre 1
                    // myPlate->D1(P2d.X(), P2d.Y(), PP, v1h, v2h);  (L376)
                    let (_pp, v1h, v2h) = surf_plate.eval_d1(p2d.x, p2d.y);
                    // v3h = v1h ^ v2h;
                    let v3h = v1h.cross(v2h);
                    // gp_XYZ P3d(v3h.X(), v3h.Y(), v3h.Z()); Seq3d.Append(P3d);
                    seq3d.push(DVec3::new(v3h.x, v3h.y, v3h.z));
                }
            }
        }

        // double U0 = 0., U1 = 0., V0 = 0., V1 = 0.;
        // myPlate->RealBounds(U0, U1, V0, V1);           (L384-385)
        // U0 = EnlargeCoeff * U0; ...                    (L386-389)
        let (mut u0, mut u1, mut v0, mut v1) = surf_plate.real_bounds();
        u0 = enlarge_coeff * u0;
        u1 = enlarge_coeff * u1;
        v0 = enlarge_coeff * v0;
        v1 = enlarge_coeff * v1;

        // double seuil = Tol3d;                          (L391)
        let mut seuil = tol3d;
        // if (CritOrder == 0 && Tol3d < 10 * dmax)       (L392-399)
        if crit_order == 0 && tol3d < 10.0 * dmax {
            seuil = 10.0 * dmax;
        }
        // if (CritOrder == 1 && Tol3d < 10 * dmax)       (L400-407)
        if crit_order == 1 && tol3d < 10.0 * dmax {
            seuil = 10.0 * dmax;
        }

        // int nb1 = 0, nb2 = 0, nb3 = 1;                 (L408)
        let nb1 = 0;
        let nb2 = 0;
        let nb3 = 1;
        // nul1 = new NCollection_HArray1<double>(1, 1); nul1->Init(0.);   (L409-410)
        let mut nul1 = vec![0.0f64; 1];
        nul1.fill(0.0);
        // nul2 = new NCollection_HArray2<double>(1, 1, 1, 4); nul2->Init(0.); (L411-412)
        let mut nul2 = Array2::<f64>::new(1, 1, 1, 4);
        nul2.init(0.0);
        // eps3D = new NCollection_HArray1<double>(1, 1); eps3D->Init(Tol3d); (L413-414)
        let mut eps3d = vec![0.0f64; 1];
        eps3d.fill(tol3d);
        // epsfr = new NCollection_HArray2<double>(1, 1, 1, 4); epsfr->Init(Tol3d); (L415-416)
        let mut epsfr = Array2::<f64>::new(1, 1, 1, 4);
        epsfr.init(tol3d);

        // GeomAbs_IsoType myType = GeomAbs_IsoV;         (L418)
        let my_type = IsoType::IsoV;
        // int myPrec = 0;                                (L419)
        let mut my_prec = 0;

        // AdvApprox_DichoCutting myDec;                  (L421)
        let my_dec = DichoCutting;

        // if (CritOrder == -1)                           (L423-458)
        if crit_order == -1 {
            my_prec = 1;
            // POP pour NT
            // GeomPlate_MakeApprox_Eval ev(myPlate);      (L427)
            let ev = GeomPlateMakeApproxEval::new(surf_plate);
            // AdvApp2Var_ApproxAFunc2Var AppPlate(...)    (L428-450)
            let app_plate = ApproxAFunc2Var::new(
                nb1,
                nb2,
                nb3,
                &nul1,
                &nul1,
                &eps3d,
                &nul2,
                &nul2,
                &epsfr,
                u0,
                u1,
                v0,
                v1,
                my_type,
                continuity,
                continuity,
                my_prec,
                dgmax,
                dgmax,
                nbmax,
                &ev,
                &my_dec,
                &my_dec,
            );
            // mySurface   = AppPlate.Surface(1);          (L451)
            r.my_surface = app_plate.surface(1).cloned();
            // myAppError  = AppPlate.MaxError(3, 1);      (L452)
            r.my_app_error = app_plate.max_error_at(3, 1);
            // myCritError = 0.;                           (L453)
            r.my_crit_error = 0.0;
        }
        // else if (CritOrder == 0)                       (L459-497)
        else if crit_order == 0 {
            // GeomPlate_PlateG0Criterion Crit0(Seq2d, Seq3d, seuil);  (L461)
            // (default args: Type = Absolute, Repart = Regular)
            let crit0 = PlateG0Criterion::new(
                &seq2d,
                &seq3d,
                seuil,
                AdvApp2VarCriterionType::Absolute,
                AdvApp2VarCriterionRepartition::Regular,
            );
            // POP pour NT
            // GeomPlate_MakeApprox_Eval ev(myPlate);      (L463)
            let ev = GeomPlateMakeApproxEval::new(surf_plate);
            // AdvApp2Var_ApproxAFunc2Var AppPlate(...)    (L464-488)
            let app_plate = ApproxAFunc2Var::new_with_criterion(
                nb1,
                nb2,
                nb3,
                &nul1,
                &nul1,
                &eps3d,
                &nul2,
                &nul2,
                &epsfr,
                u0,
                u1,
                v0,
                v1,
                my_type,
                continuity,
                continuity,
                my_prec,
                dgmax,
                dgmax,
                nbmax,
                &ev,
                &crit0,
                &my_dec,
                &my_dec,
            );
            // mySurface   = AppPlate.Surface(1);          (L489)
            r.my_surface = app_plate.surface(1).cloned();
            // myAppError  = AppPlate.MaxError(3, 1);      (L490)
            r.my_app_error = app_plate.max_error_at(3, 1);
            // myCritError = AppPlate.CritError(3, 1);     (L491)
            r.my_crit_error = app_plate.crit_error(3, 1);
        }
        // else if (CritOrder == 1)                       (L498-536)
        else if crit_order == 1 {
            // GeomPlate_PlateG1Criterion Crit1(Seq2d, Seq3d, seuil);  (L500)
            // (default args: Type = Absolute, Repart = Regular)
            let crit1 = PlateG1Criterion::new(
                &seq2d,
                &seq3d,
                seuil,
                AdvApp2VarCriterionType::Absolute,
                AdvApp2VarCriterionRepartition::Regular,
            );
            // POP pour NT
            // GeomPlate_MakeApprox_Eval ev(myPlate);      (L502)
            let ev = GeomPlateMakeApproxEval::new(surf_plate);
            // AdvApp2Var_ApproxAFunc2Var AppPlate(...)    (L503-527)
            let app_plate = ApproxAFunc2Var::new_with_criterion(
                nb1,
                nb2,
                nb3,
                &nul1,
                &nul1,
                &eps3d,
                &nul2,
                &nul2,
                &epsfr,
                u0,
                u1,
                v0,
                v1,
                my_type,
                continuity,
                continuity,
                my_prec,
                dgmax,
                dgmax,
                nbmax,
                &ev,
                &crit1,
                &my_dec,
                &my_dec,
            );
            // mySurface   = AppPlate.Surface(1);          (L528)
            r.my_surface = app_plate.surface(1).cloned();
            // myAppError  = AppPlate.MaxError(3, 1);      (L529)
            r.my_app_error = app_plate.max_error_at(3, 1);
            // myCritError = AppPlate.CritError(3, 1);     (L530)
            r.my_crit_error = app_plate.crit_error(3, 1);
        }
        r
    }

    /// OCCT Surface() (GeomPlate_MakeApprox.cxx L541-544).
    pub fn surface(&self) -> Option<&BSplineSurface> {
        self.my_surface.as_ref()
    }

    /// OCCT ApproxError() (GeomPlate_MakeApprox.cxx L548-551).
    pub fn approx_error(&self) -> f64 {
        self.my_app_error
    }

    /// OCCT CriterionError() (GeomPlate_MakeApprox.cxx L555-558).
    pub fn criterion_error(&self) -> f64 {
        self.my_crit_error
    }

    /// OCCT myPlate member access (hxx) — kept for the AdvApp2Var wiring.
    pub fn plate(&self) -> &GeomPlateSurface {
        &self.my_plate
    }
}
