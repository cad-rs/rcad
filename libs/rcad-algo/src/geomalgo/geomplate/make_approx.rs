//! OCCT GeomPlate_MakeApprox (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_MakeApprox.cxx.
//!
//! Complete: the GeomPlate_MakeApprox_Eval evaluator functor (L37-265) and
//! both constructors up to the AdvApp2Var_ApproxAFunc2Var call, including
//! the Seq2d/Seq3d constraint collection, the RealBounds/EnlargeCoeff
//! scaling and the seuil clamping.
//!
//! GAP leaf: `AdvApp2Var_ApproxAFunc2Var` (with `AdvApprox_DichoCutting`,
//! `AdvApp2Var_Criterion` dispatch and `AdvApp2Var_Patch` evaluation) is the
//! untranslated AdvApp2Var package.  Each OCCT AppPlate construction site is
//! preserved at its exact anchor; until it lands the approximation result
//! stays null (mySurface = AppPlate.Surface(1) never executes), which keeps
//! the OCCT failure path where a failed approximation leaves the surface
//! null and the caller falls back (GeomPlate_BuildPlateSurface Perform
//! L537-580 chains on the MakeApprox result).
//!
//! `handle(Geom_BSplineSurface)` -> `Option<BSplineSurface>` (kernel type;
//! None models the null handle).

use glam::{DVec2, DVec3};

use rcad_kernel::geom::BSplineSurface;

use super::plate_g0_criterion::PlateG0Criterion;
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

    /// OCCT Evaluate (L62-265) — the raw-buffer calling convention of
    /// AdvApp2Var_EvaluatorFunc2Var is preserved (all parameters in OCCT
    /// declaration order; the Result layout is Result[Dimension, N]).
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
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
        _plate_crit: &PlateG0Criterion,
        tol3d: f64,
        nbmax: i32,
        dgmax: i32,
        continuity: i32,
        enlarge_coeff: f64,
    ) -> Self {
        let _ = (tol3d, nbmax, dgmax, continuity, enlarge_coeff);
        // myPlate = SurfPlate;
        let mut r = MakeApprox {
            my_plate: surf_plate.clone(),
            my_surface: None,
            my_app_error: 0.0,
            my_crit_error: 0.0,
        };

        // double U0 = 0., U1 = 0., V0 = 0., V1 = 0.;
        // myPlate->RealBounds(U0, U1, V0, V1);
        // U0 = EnlargeCoeff * U0; ... (cxx L279-284)
        let (mut u0, mut u1, mut v0, mut v1) = surf_plate.real_bounds();
        u0 = enlarge_coeff * u0;
        u1 = enlarge_coeff * u1;
        v0 = enlarge_coeff * v0;
        v1 = enlarge_coeff * v1;
        let _ = (u0, u1, v0, v1);

        // The AdvApp2Var_ApproxAFunc2Var construction (cxx L286-326) —
        // GAP leaf (untranslated AdvApp2Var package); mySurface stays null,
        // preserving the failed-approximation path.
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
        continuity: i32,
        enlarge_coeff: f64,
    ) -> Self {
        let _ = (nbmax, dgmax, continuity);
        // myPlate = SurfPlate;
        let mut r = MakeApprox {
            my_plate: surf_plate.clone(),
            my_surface: None,
            my_app_error: 0.0,
            my_crit_error: 0.0,
        };

        // NCollection_Sequence<gp_XY> Seq2d; NCollection_Sequence<gp_XYZ> Seq3d;
        let mut seq2d: Vec<DVec2> = Vec::new();
        let mut seq3d: Vec<DVec3> = Vec::new();

        if crit_order >= 0 {
            // contraintes 2d d'ordre 0
            surf_plate.constraints(&mut seq2d);

            // contraintes 3d correspondantes sur plate
            let nbp = seq2d.len();
            for p2d in &seq2d[..nbp] {
                if crit_order == 0 {
                    // a l'ordre 0
                    let pp = surf_plate.eval_d0(p2d.x, p2d.y);
                    seq3d.push(DVec3::new(pp.x, pp.y, pp.z));
                } else {
                    // a l'ordre 1
                    let (_pp, v1h, v2h) = surf_plate.eval_d1(p2d.x, p2d.y);
                    let v3h = v1h.cross(v2h);
                    seq3d.push(DVec3::new(v3h.x, v3h.y, v3h.z));
                }
            }
        }

        // double U0 = 0., U1 = 0., V0 = 0., V1 = 0.;
        // myPlate->RealBounds(U0, U1, V0, V1); U0 = EnlargeCoeff * U0; ...
        // (cxx L384-389)
        let (u0, u1, v0, v1) = surf_plate.real_bounds();
        let (u0, u1, v0, v1) = (
            enlarge_coeff * u0,
            enlarge_coeff * u1,
            enlarge_coeff * v0,
            enlarge_coeff * v1,
        );

        // double seuil = Tol3d; if (CritOrder == 0 && Tol3d < 10 * dmax)
        // seuil = 10 * dmax; ... (cxx L391-407)
        let mut seuil = tol3d;
        if crit_order == 0 && tol3d < 10.0 * dmax {
            seuil = 10.0 * dmax;
        }
        if crit_order == 1 && tol3d < 10.0 * dmax {
            seuil = 10.0 * dmax;
        }

        // int nb1 = 0, nb2 = 0, nb3 = 1; ... nul1/nul2/eps3D/epsfr setup
        // (cxx L408-419).
        let _ = (u0, u1, v0, v1, seuil);

        // if (CritOrder == -1) { ... AppPlate ... }        (cxx L423-458)
        // else if (CritOrder == 0) { Crit0; AppPlate }     (cxx L459-497)
        // else if (CritOrder == 1) { Crit1; AppPlate }     (cxx L498-536)
        //
        // GAP leaf: every branch constructs AdvApp2Var_ApproxAFunc2Var over
        // the untranslated AdvApp2Var package; the three construction
        // anchors are preserved above and mySurface stays null (the OCCT
        // failed-approximation path), until AdvApp2Var lands.
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
