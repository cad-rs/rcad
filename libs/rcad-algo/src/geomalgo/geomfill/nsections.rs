//! OCCT GeomFill_NSections (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_NSections.hxx (L35-218) + GeomFill_NSections.cxx (whole file
//! L26-1047): a GeomFill_SectionLaw built from N section curves.
//!
//! Architecture mappings: `NCollection_Sequence<Geom_Curve>` ->
//! `Vec<Curve3>`; `Geom_BSplineSurface` -> rcad `BSplineSurface` (flat
//! knots, no periodic flag — the surface built by the constructor is
//! non-periodic, so `IsVPeriodic`/`IsUPeriodic` report false and the
//! `SetVNotPeriodic` branches are anchor-out-of-scope).  The multi-curve
//! ComputeSurface path depends on the GeomFill_AppSurf / SectionGenerator
//! approximation machinery (not yet ported) — it follows the ThruSections
//! precedent as an `unimplemented!()` skeleton; the single-curve path
//! early-returns in OCCT itself (the OCC27875 GTest anchor).
//! `GCPnts_AbscissaPoint::Length` -> rcad `gcpnts::abscissa_point::
//! arc_length` over the full parameter range.

use glam::DVec3;

use rcad_kernel::base::gcpnts::abscissa_point::arc_length;
use rcad_kernel::math::bspl_lib::{
    increase_degree as bspl_increase_degree, intervals as bspl_intervals,
};
use rcad_kernel::math::bspl::de_boor_homo;
use rcad_kernel::math::bspl_lib::eval_flat;
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::geom::{BSplineCurve3, BSplineSurface, Circle3, Curve3, CurveEval, SurfaceEval};

use super::section_law::SectionLaw;

// OCCT Precision / local constants.
const P_CONFUSION: f64 = 1e-12;

/// OCCT GeomAbs_Shape -> the integer continuity rank consumed by
/// BSplCLib::Intervals (C0 = 0 ... C3 = 3, CN = a value above any degree).
fn continuity_rank(s: GeomAbsShape) -> i32 {
    match s {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::C1 => 1,
        GeomAbsShape::C2 => 2,
        GeomAbsShape::C3 => 3,
        GeomAbsShape::CN => 1000,
    }
}

/// The (knots, mults) form of a flat knot vector.
fn knots_mults_of(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for &k in flat {
        if knots.last() == Some(&k) {
            *mults.last_mut().unwrap() += 1;
        } else {
            knots.push(k);
            mults.push(1);
        }
    }
    (knots, mults)
}

/// OCCT static ResultEval (L142-201) — evaluates the poles/weights
/// derivatives of `surf` in V at `v` into the flat `result` array.
fn result_eval(surf: &BSplineSurface, v: f64, deriv: usize, result: &mut [f64]) {
    let rational = surf.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
    let gap = if rational { 4 } else { 3 };
    let cdeg = surf.degree_v;
    let cdim = surf.control_points.len() * gap;
    let nb_p = surf.control_points[0].len();
    // les noeuds plats
    let fknots = &surf.knots_v;
    // les poles (flattened, homogeneous when rational)
    let mut spoles = vec![0.0f64; cdim * nb_p];
    let mut ipole = 0usize;
    for jj in 0..nb_p {
        for ii in 0..surf.control_points.len() {
            let p = surf.control_points[ii][jj];
            spoles[ipole] = p.x;
            spoles[ipole + 1] = p.y;
            spoles[ipole + 2] = p.z;
            if rational {
                let w = surf.weights[ii][jj];
                spoles[ipole + 3] = w;
                spoles[ipole] *= w;
                spoles[ipole + 1] *= w;
                spoles[ipole + 2] *= w;
            }
            ipole += gap;
        }
    }
    let mut extrap_mode = [cdeg as i32, cdeg as i32];
    let mut eval_bs = vec![0.0f64; cdim * (deriv + 1)];
    eval_flat(
        v,
        false,
        deriv as i32,
        &mut extrap_mode,
        cdeg,
        fknots,
        cdim,
        &spoles,
        &mut eval_bs,
    );
    for ii in 0..cdim {
        result[ii] = eval_bs[ii + deriv * cdim];
    }
}

/// OCCT GeomFill_NSections.
#[derive(Debug, Clone)]
pub struct NSections {
    u_first: f64,
    u_last: f64,
    v_first: f64,
    v_last: f64,
    my_sections: Vec<Curve3>,
    my_trsfs: Vec<Trsf>,
    my_params: Vec<f64>,
    my_surface: Option<BSplineSurface>,
    my_ref_surf: Option<BSplineSurface>,
}

impl NSections {
    /// OCCT ctor (NC) (L203-212).
    pub fn new(nc: Vec<Curve3>) -> Self {
        let mut law = NSections {
            u_first: 0.0,
            u_last: 1.0,
            v_first: 0.0,
            v_last: 1.0,
            my_sections: nc,
            my_trsfs: Vec::new(),
            my_params: Vec::new(),
            my_surface: None,
            my_ref_surf: None,
        };
        law.compute_surface();
        law
    }

    /// OCCT ctor (NC, NP) (L214-224).
    pub fn new_with_params(nc: Vec<Curve3>, np: Vec<f64>) -> Self {
        let mut law = NSections {
            u_first: 0.0,
            u_last: 1.0,
            v_first: 0.0,
            v_last: 1.0,
            my_sections: nc,
            my_trsfs: Vec::new(),
            my_params: np,
            my_surface: None,
            my_ref_surf: None,
        };
        law.compute_surface();
        law
    }

    /// OCCT ctor (NC, NP, UF, UL) (L226-236).
    pub fn new_with_params_bounds(nc: Vec<Curve3>, np: Vec<f64>, uf: f64, ul: f64) -> Self {
        let mut law = NSections {
            u_first: uf,
            u_last: ul,
            v_first: 0.0,
            v_last: 1.0,
            my_sections: nc,
            my_trsfs: Vec::new(),
            my_params: np,
            my_surface: None,
            my_ref_surf: None,
        };
        law.compute_surface();
        law
    }

    /// OCCT ctor (NC, NP, UF, UL, VF, VL) (L238-250).
    pub fn new_with_bounds(
        nc: Vec<Curve3>,
        np: Vec<f64>,
        uf: f64,
        ul: f64,
        vf: f64,
        vl: f64,
    ) -> Self {
        let mut law = NSections {
            u_first: uf,
            u_last: ul,
            v_first: vf,
            v_last: vl,
            my_sections: nc,
            my_trsfs: Vec::new(),
            my_params: np,
            my_surface: None,
            my_ref_surf: None,
        };
        law.compute_surface();
        law
    }

    /// OCCT ctor (NC, Trsfs, NP, UF, UL, VF, VL, Surf) (L252-268).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_reference(
        nc: Vec<Curve3>,
        trsfs: Vec<Trsf>,
        np: Vec<f64>,
        uf: f64,
        ul: f64,
        vf: f64,
        vl: f64,
        surf: Option<BSplineSurface>,
    ) -> Self {
        let mut law = NSections {
            u_first: uf,
            u_last: ul,
            v_first: vf,
            v_last: vl,
            my_sections: nc,
            my_trsfs: trsfs,
            my_params: np,
            my_surface: None,
            my_ref_surf: surf,
        };
        law.compute_surface();
        law
    }

    /// OCCT ComputeSurface (L502-617).
    ///
    /// The single-section path early-returns in OCCT ("We will not be able
    /// to create surface from single curve").  The multi-section path uses
    /// the GeomFill_SectionGenerator + GeomFill_AppSurf approximation
    /// machinery, which is not ported yet (ThruSections precedent).
    pub fn compute_surface(&mut self) {
        let bs: Option<BSplineSurface>;
        if self.my_ref_surf.is_none() {
            let jfin = self.my_sections.len();
            if jfin <= 1 {
                // We will not be able to create surface from single curve.
                // (mySurface stays null; D0/D1/D2 report false — the
                // OCC27875 anchor exercises exactly this path.)
                bs = None;
            } else {
                let my_pres3d = 1.0e-06f64;
                let _ = my_pres3d;
                // OCCT: SectionGenerator + AppSurf approximation to build BS.
                unimplemented!(
                    "GeomFill_NSections multi-curve ComputeSurface needs the                      GeomFill_AppSurf/SectionGenerator approximation machinery (not ported)"
                );
            }
        } else {
            // segmentation de myRefSurf
            let ref_surf = self.my_ref_surf.as_ref().unwrap().clone();
            let uk = flat_knots_u(&ref_surf);
            let (uk_distinct, _um) = knots_mults_of(&uk);
            let mut ui1 = snap_to_knot(self.u_first, &uk_distinct, P_CONFUSION);
            let mut ui2 = snap_to_knot(self.u_last, &uk_distinct, P_CONFUSION);
            let v0 = ref_surf.knots_v[ref_surf.degree_v];
            let v1 = ref_surf.knots_v[ref_surf.knots_v.len() - ref_surf.degree_v - 1];
            let _ = (&mut ui1, &mut ui2, v0, v1);
            // OCCT: BS = myRefSurf->Copy(); BS->CheckAndSegment(Ui1, Ui2, V0, V1)
            // — the surface segmentation kernel is not ported yet.
            unimplemented!(
                "GeomFill_NSections reference-surface segmentation                  (Geom_BSplineSurface::CheckAndSegment) is not ported"
            );
        }
        self.my_surface = bs;
        // On augmente le degre pour que le positionnement D2 soit correct
        if let Some(surface) = &mut self.my_surface {
            if surface.degree_v < 2 {
                // OCCT: IncreaseDegree(UDegree, 2) in V.
                increase_degree_v(surface, 2);
            }
        }
    }

    /// OCCT SetSurface (L? hxx) — sets the reference surface.
    pub fn set_surface(&mut self, ref_surf: BSplineSurface) {
        self.my_ref_surf = Some(ref_surf);
    }

    /// OCCT D0 helper — the V-iso curve of the surface at `v`.
    fn v_iso(&self, v: f64) -> BSplineCurve3 {
        let surf = self.my_surface.as_ref().unwrap();
        let rational = surf.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
        let nb_u = surf.control_points.len();
        let nb_v = surf.control_points[0].len();
        let mut poles = Vec::with_capacity(nb_u);
        let mut weights = Vec::with_capacity(nb_u);
        for ii in 0..nb_u {
            let column: Vec<DVec3> = (0..nb_v).map(|jj| surf.control_points[ii][jj]).collect();
            let column_w: Vec<f64> = (0..nb_v)
                .map(|jj| if rational { surf.weights[ii][jj] } else { 1.0 })
                .collect();
            let h = de_boor_homo(surf.degree_v, &surf.knots_v, &column, &column_w, v);
            weights.push(h[3]);
            poles.push(DVec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
        }
        BSplineCurve3 {
            degree: surf.degree_u,
            knots: surf.knots_u.clone(),
            control_points: poles,
            weights,
            is_periodic: false,
        }
    }
}

fn flat_knots_u(surf: &BSplineSurface) -> Vec<f64> {
    surf.knots_u.clone()
}

/// Snap `u` to a knot within `tol` (the OCCT LocateU + knot-equality idiom).
fn snap_to_knot(u: f64, distinct_knots: &[f64], tol: f64) -> f64 {
    for &k in distinct_knots {
        if (u - k).abs() <= tol {
            return k;
        }
    }
    u
}

/// OCCT Geom_BSplineSurface::IncreaseDegree in the V direction only —
/// degree elevation per U row through the 1D kernel.
fn increase_degree_v(surf: &mut BSplineSurface, new_degree: usize) {
    let nb_u = surf.control_points.len();
    let nb_v = surf.control_points[0].len();
    let rational = surf.weights.iter().any(|col| col.iter().any(|&w| w != 1.0));
    let (vk, vm) = knots_mults_of(&surf.knots_v);
    let mut new_poles = vec![vec![DVec3::ZERO; nb_v + (new_degree - surf.degree_v)]; nb_u];
    let mut new_weights = vec![vec![1.0f64; new_poles[0].len()]; nb_u];
    let mut new_knots = vec![0.0f64; vk.len()];
    let mut new_mults = vm.clone();
    for ii in 0..nb_u {
        let column: Vec<DVec3> = (0..nb_v).map(|jj| surf.control_points[ii][jj]).collect();
        let column_w: Vec<f64> = (0..nb_v)
            .map(|jj| if rational { surf.weights[ii][jj] } else { 1.0 })
            .collect();
        let (flat_p, dim) = flatten_column(&column, &column_w);
        let mut new_flat = vec![0.0f64; new_poles[0].len() * dim];
        bspl_increase_degree(
            surf.degree_v,
            new_degree,
            false,
            dim,
            &flat_p,
            &vk,
            &vm,
            &mut new_flat,
            &mut new_knots,
            &mut new_mults,
        );
        for jj in 0..new_poles[0].len() {
            if dim == 4 {
                let w = new_flat[jj * 4 + 3];
                new_poles[ii][jj] =
                    DVec3::new(new_flat[jj * 4] / w, new_flat[jj * 4 + 1] / w, new_flat[jj * 4 + 2] / w);
                new_weights[ii][jj] = w;
            } else {
                new_poles[ii][jj] =
                    DVec3::new(new_flat[jj * 3], new_flat[jj * 3 + 1], new_flat[jj * 3 + 2]);
            }
        }
    }
    surf.degree_v = new_degree;
    surf.control_points = new_poles;
    surf.weights = new_weights;
    // rebuild the flat V knot vector from the elevated (knots, mults)
    surf.knots_v = {
        let mut flat = Vec::with_capacity(new_mults.iter().map(|&m| m as usize).sum());
        for (k, m) in new_knots.iter().zip(new_mults.iter()) {
            for _ in 0..*m {
                flat.push(*k);
            }
        }
        flat
    };
}

fn flatten_column(column: &[DVec3], weights: &[f64]) -> (Vec<f64>, usize) {
    let mut flat = Vec::with_capacity(column.len() * 4);
    for (p, w) in column.iter().zip(weights.iter()) {
        flat.extend([p.x * w, p.y * w, p.z * w, *w]);
    }
    (flat, 4)
}

impl SectionLaw for NSections {
    /// OCCT D0 (L270-294) — the V-iso poles/weights.
    fn d0(&self, v: f64, poles: &mut [DVec3], weights: &mut [f64]) -> bool {
        let Some(surface) = &self.my_surface else {
            return false;
        };
        let iso = self.v_iso(v);
        let l = poles.len();
        for ii in 0..l {
            poles[ii] = iso.control_points[ii];
            weights[ii] = iso.weights[ii];
        }
        let _ = surface;
        true
    }

    /// OCCT D1 (L296-382).
    fn d1(
        &self,
        v: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        let Some(surface) = &self.my_surface else {
            return false;
        };
        let ok = self.d0(v, poles, weights);
        if !ok {
            return false;
        }
        let l = poles.len();
        let derivative_request = 1usize;
        let rational = surface.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
        let gap = if rational { 4 } else { 3 };
        let dim_result = surface.control_points.len() * gap;
        // OCCT IsVPeriodic branch: rcad BSplineSurface is non-periodic.
        let mut result = vec![0.0f64; dim_result];
        result_eval(surface, v, derivative_request, &mut result);
        let eps_w = 10.0 * P_CONFUSION;
        let mut null_weight = false;
        if !rational {
            for dw in dweights.iter_mut() {
                *dw = 0.0;
            }
        }
        let mut indice = 1usize;
        // recopie des poles du resultat sous forme de points 3D et de poids
        for ii in 0..l {
            if null_weight {
                break;
            }
            dpoles[ii] = DVec3::new(result[indice - 1], result[indice], result[indice + 1]);
            if rational {
                let ww = weights[ii];
                if ww < eps_w {
                    null_weight = true;
                } else {
                    dweights[ii] = result[indice + 2];
                    dpoles[ii] = (dpoles[ii] - dweights[ii] * poles[ii]) / ww;
                }
            }
            indice += gap;
        }
        !null_weight
    }

    /// OCCT D2 (L384-498).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        v: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        weights: &mut [f64],
        dweights: &mut [f64],
        d2weights: &mut [f64],
    ) -> bool {
        let Some(surface) = &self.my_surface else {
            return false;
        };
        // pb dans BSplCLib::Eval() pour les surfaces rationnelles de degre 1
        // si l'ordre de derivation est egal a 2.
        if surface.degree_v < 2 {
            return false;
        }
        let ok = self.d1(v, poles, dpoles, weights, dweights);
        if !ok {
            return false;
        }
        let l = poles.len();
        let derivative_request = 2usize;
        let rational = surface.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
        let gap = if rational { 4 } else { 3 };
        let dim_result = surface.control_points.len() * gap;
        let mut result = vec![0.0f64; dim_result];
        result_eval(surface, v, derivative_request, &mut result);
        let eps_w = 10.0 * P_CONFUSION;
        let mut null_weight = false;
        if !rational {
            for dw in dweights.iter_mut() {
                *dw = 0.0;
            }
            for d2w in d2weights.iter_mut() {
                *d2w = 0.0;
            }
        }
        let mut indice = 1usize;
        for ii in 0..l {
            if null_weight {
                break;
            }
            dpoles[ii] = DVec3::new(result[indice - 1], result[indice], result[indice + 1]);
            if rational {
                let ww = weights[ii];
                if ww < eps_w {
                    null_weight = true;
                } else {
                    dweights[ii] = result[indice + 2];
                    dpoles[ii] = (dpoles[ii] - dweights[ii] * poles[ii]) / ww;
                }
            }
            d2poles[ii] = DVec3::new(result[indice - 1], result[indice], result[indice + 1]);
            if rational {
                let ww = weights[ii];
                if ww < eps_w {
                    null_weight = true;
                } else {
                    d2weights[ii] = result[indice + 2];
                    d2poles[ii] = (d2poles[ii] - d2weights[ii] * poles[ii]
                        - 2.0 * dweights[ii] * dpoles[ii])
                        / ww;
                }
            }
            indice += gap;
        }
        !null_weight
    }

    /// OCCT BSplineSurface (hxx) — the computed surface.
    fn bspline_surface(&self) -> Option<&BSplineSurface> {
        self.my_surface.as_ref()
    }

    /// OCCT SectionShape (L619-629).
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize) {
        let Some(surface) = &self.my_surface else {
            return;
        };
        *nb_poles = surface.control_points.len();
        *nb_knots = knots_mults_of(&surface.knots_u).0.len();
        *degree = surface.degree_u;
    }

    /// OCCT Knots (L631-637).
    fn knots(&self, t_knots: &mut [f64]) {
        let Some(surface) = &self.my_surface else {
            return;
        };
        let (knots, _) = knots_mults_of(&surface.knots_u);
        t_knots[..knots.len()].copy_from_slice(&knots);
    }

    /// OCCT Mults (L639-645).
    fn mults(&self, t_mults: &mut [i32]) {
        let Some(surface) = &self.my_surface else {
            return;
        };
        let (_, mults) = knots_mults_of(&surface.knots_u);
        t_mults[..mults.len()].copy_from_slice(&mults);
    }

    /// OCCT IsRational (L647-656).
    fn is_rational(&self) -> bool {
        let Some(surface) = &self.my_surface else {
            return false;
        };
        surface
            .weights
            .iter()
            .any(|col| col.iter().any(|&w| w != 1.0))
    }

    /// OCCT IsUPeriodic (L658-666) — rcad BSplineSurface is non-periodic.
    fn is_u_periodic(&self) -> bool {
        false
    }

    /// OCCT IsVPeriodic (L668-676) — rcad BSplineSurface is non-periodic.
    fn is_v_periodic(&self) -> bool {
        false
    }

    /// OCCT NbIntervals (L678-690).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let Some(surface) = &self.my_surface else {
            return 0;
        };
        let (vk, vm) = knots_mults_of(&surface.knots_v);
        let arr = bspl_intervals(
            &vk,
            &vm,
            surface.degree_v,
            false,
            continuity_rank(s),
            surface.knots_v[surface.degree_v],
            surface.knots_v[surface.knots_v.len() - surface.degree_v - 1],
            P_CONFUSION,
        );
        arr.len() - 1
    }

    /// OCCT Intervals (L692-704).
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let Some(surface) = &self.my_surface else {
            return;
        };
        let (vk, vm) = knots_mults_of(&surface.knots_v);
        let arr = bspl_intervals(
            &vk,
            &vm,
            surface.degree_v,
            false,
            continuity_rank(s),
            surface.knots_v[surface.degree_v],
            surface.knots_v[surface.knots_v.len() - surface.degree_v - 1],
            P_CONFUSION,
        );
        t[..arr.len()].copy_from_slice(&arr);
    }

    /// OCCT SetInterval (L706-714) — "rien a faire : mySurface est supposee
    /// Cn en V".
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// OCCT GetInterval (L716-721).
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        *first = self.v_first;
        *last = self.v_last;
    }

    /// OCCT GetDomain (L723-728).
    fn get_domain(&self, first: &mut f64, last: &mut f64) {
        *first = self.v_first;
        *last = self.v_last;
    }

    /// OCCT GetTolerance (L730-742).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, _angle_tol: f64, tol3d: &mut [f64]) {
        for t in tol3d.iter_mut() {
            *t = surf_tol;
        }
        if bound_tol < surf_tol {
            tol3d[0] = bound_tol;
            tol3d[tol3d.len() - 1] = bound_tol;
        }
    }

    /// OCCT BarycentreOfSurf (L791-817).
    fn barycentre_of_surf(&self) -> DVec3 {
        let mut bary = DVec3::ZERO;
        let Some(surface) = &self.my_surface else {
            return bary;
        };
        let (_u0, u1, _v0, v1) = (
            surface.knots_u[surface.degree_u],
            surface.knots_u[surface.knots_u.len() - surface.degree_u - 1],
            surface.knots_v[surface.degree_v],
            surface.knots_v[surface.knots_v.len() - surface.degree_v - 1],
        );
        let u0 = surface.knots_u[surface.degree_u];
        let v0 = surface.knots_v[surface.degree_v];
        let delta_v = (v1 - v0) / 20.0;
        let delta_u = (u1 - u0) / 20.0;
        let mut v = v0;
        for _jj in 0..=20 {
            let mut u = u0;
            for _ii in 0..=20 {
                let p = surface.point_at(u, v);
                bary += p;
                u += delta_u;
            }
            v += delta_v;
        }
        bary /= (21.0 * 21.0);
        bary
    }

    /// OCCT MaximalSection (L819-836).
    fn maximal_section(&self) -> f64 {
        let mut lmax = 0.0f64;
        for section in &self.my_sections {
            let first = section.default_domain()[0];
            let last = section.default_domain()[1];
            let l = arc_length(section, first, last);
            if l > lmax {
                lmax = l;
            }
        }
        lmax
    }

    /// OCCT GetMinimalWeight (L838-878).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        let Some(surface) = &self.my_surface else {
            return;
        };
        let u_rational = surface
            .weights
            .iter()
            .any(|col| col.iter().any(|&w| w != 1.0));
        if u_rational {
            let nb_u = surface.control_points.len();
            let nb_v = surface.control_points[0].len();
            for i in 0..nb_u {
                let mut min = surface.weights[i][0];
                for j in 1..nb_v {
                    if min > surface.weights[i][j] {
                        min = surface.weights[i][j];
                    }
                }
                weights[i] = min;
            }
        } else {
            for w in weights.iter_mut() {
                *w = 1.0;
            }
        }
    }

    /// OCCT IsConstant (L880-948).
    fn is_constant(&self, error: &mut f64) -> bool {
        // on se limite a 2 sections
        let mut isconst = self.my_sections.len() == 2;
        *error = 0.0;
        if isconst {
            let c1 = &self.my_sections[0];
            let c2 = &self.my_sections[1];
            // les sections doivent avoir le meme type
            isconst = std::mem::discriminant(c1) == std::mem::discriminant(c2);
            if isconst {
                match (c1, c2) {
                    (Curve3::Circle(a), Curve3::Circle(b)) => {
                        let tol = 1.0e-7;
                        let samedir = axes_parallel(a.normal, b.normal, 1.0e-4);
                        let samerad = (a.radius - b.radius).abs() < tol;
                        let mut samepos = a.center.distance(b.center) < tol;
                        if !samepos {
                            let d = b.center - a.center;
                            samepos = axes_parallel(a.normal, d, 1.0e-4);
                        }
                        isconst = samedir && samerad && samepos;
                    }
                    (Curve3::Line(l1), Curve3::Line(l2)) => {
                        let tol = 1.0e-7;
                        let samedir = axes_parallel(l1.direction, l2.direction, 1.0e-4);
                        let p11 = l1.origin;
                        let p12 = l1.origin + l1.direction;
                        let p21 = l2.origin;
                        let p22 = l2.origin + l2.direction;
                        let samelength =
                            (p11.distance(p12) - p21.distance(p22)).abs() < tol;
                        // l'ecart entre les 2 sections ne compte pas
                        let samepos = (p11.distance(p21) < tol && p12.distance(p22) < tol)
                            || (p12.distance(p21) < tol && p11.distance(p22) < tol);
                        isconst = samedir && samelength && samepos;
                    }
                    _ => {
                        isconst = false;
                    }
                }
            }
        }
        isconst
    }

    /// OCCT ConstantSection (L950-957).
    fn constant_section(&self) -> Curve3 {
        self.my_sections[0].clone()
    }

    /// OCCT IsConicalLaw (L959-1020).
    fn is_conical_law(&self, error: &mut f64) -> bool {
        let mut isconic = self.my_sections.len() == 2;
        *error = 0.0;
        if isconic {
            let c1 = &self.my_sections[0];
            let c2 = &self.my_sections[1];
            isconic = matches!(c1, Curve3::Circle(_)) && matches!(c2, Curve3::Circle(_));
            if isconic {
                let mut circ1 = match c1 {
                    Curve3::Circle(c) => *c,
                    _ => unreachable!(),
                };
                if !self.my_trsfs.is_empty() {
                    // C1.Transform(myTrsfs(1).Inverted())
                    let inv = self.my_trsfs[0].to_daffine3().inverse();
                    circ1.center = inv.transform_point3(circ1.center);
                    circ1.normal = inv.transform_vector3(circ1.normal);
                    circ1.x_dir = inv.transform_vector3(circ1.x_dir);
                    circ1.y_dir = inv.transform_vector3(circ1.y_dir);
                }
                let mut circ2 = match c2 {
                    Curve3::Circle(c) => *c,
                    _ => unreachable!(),
                };
                if !self.my_trsfs.is_empty() {
                    let inv = self.my_trsfs[1].to_daffine3().inverse();
                    circ2.center = inv.transform_point3(circ2.center);
                    circ2.normal = inv.transform_vector3(circ2.normal);
                    circ2.x_dir = inv.transform_vector3(circ2.x_dir);
                    circ2.y_dir = inv.transform_vector3(circ2.y_dir);
                }
                let tol = 1.0e-7;
                isconic = axes_parallel(circ1.normal, circ2.normal, 1.0e-4);
                if isconic {
                    // gp_Lin Line1(C1.Axis()); Line1.Distance(C2.Location())
                    let diff = circ2.center - circ1.center;
                    let dist = diff.distance(circ1.normal * diff.dot(circ1.normal));
                    isconic = dist < tol;
                    if isconic {
                        //// Modified by jgv, 18.02.2009 for OCC20866 ////
                        let (f1, l1) = circle_bounds(&circ1);
                        let (f2, l2) = circle_bounds(&circ2);
                        isconic = (f1 - f2).abs() <= P_CONFUSION
                            && (l1 - l2).abs() <= P_CONFUSION;
                        //////////////////////////////////////////////////
                    }
                }
                let _ = tol;
            }
        }
        isconic
    }

    /// OCCT CirclSection (L1022-1047).
    fn circl_section(&self, v: f64) -> Curve3 {
        let mut err = 0.0f64;
        if !self.is_conical_law(&mut err) {
            panic!("StdFail_NotDone: The Law is not Conical!");
        }
        let c1 = match &self.my_sections[0] {
            Curve3::Circle(c) => *c,
            _ => unreachable!(),
        };
        let c2 = match &self.my_sections[self.my_sections.len() - 1] {
            Curve3::Circle(c) => *c,
            _ => unreachable!(),
        };
        let p1 = self.my_params[0];
        let p2 = self.my_params[self.my_params.len() - 1];
        let radius = (c2.radius - c1.radius) * (v - p1) / (p2 - p1) + c1.radius;
        let mut circ = c1;
        circ.radius = radius;
        let (a_par_f, a_par_l) = circle_bounds(&circ);
        let a_period = 2.0 * std::f64::consts::PI;
        if (a_par_l - a_par_f - a_period).abs() > P_CONFUSION {
            // OCCT: trimmed circle over [aParF, aParL].
            use rcad_kernel::geom::TrimmedCurve3;
            return Curve3::Trimmed(TrimmedCurve3::new(Curve3::Circle(circ), a_par_f, a_par_l));
        }
        Curve3::Circle(circ)
    }
}

fn circle_bounds(c: &Circle3) -> (f64, f64) {
    (0.0, 2.0 * std::f64::consts::PI)
}

/// OCCT gp_Ax1::IsParallel — angle <= tol or PI - angle <= tol.
/// OCCT gp_Ax1::IsParallel — Angle(Other) <= Tol or PI - Angle(Other) <= Tol.
fn axes_parallel(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    let ang = super::frenet::gp_vec_angle(a, b);
    ang <= angular_tolerance || std::f64::consts::PI - ang <= angular_tolerance
}
