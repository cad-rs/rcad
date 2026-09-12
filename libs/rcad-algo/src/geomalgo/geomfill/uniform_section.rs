//! OCCT GeomFill_UniformSection (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_UniformSection.hxx (members) + GeomFill_UniformSection.cxx
//! (whole file L35-295).
//!
//! Architecture differences:
//! - `Geom_BSplineCurve` maps to rcad `BSplineCurve3` (flat knot vector;
//!   the OCCT knots/multiplicities arrays are recovered through
//!   `knots_mults()`), `Geom_BSplineSurface` to `BSplineSurface`
//!   (no periodic flag on the rcad surface — the OCCT constructor periodic
//!   argument is carried by the knot layout).
//! - `GeomConvert::CurveToBSplineCurve(C, Convert_QuasiAngular)` maps to
//!   `rcad_kernel::base::convert::curve_to_bspline` (the rcad converter has
//!   no conversion-mode argument; the sampling conversion is its
//!   corresponding path).
//! - `GCPnts_AbscissaPoint::Length(GeomAdaptor_Curve)` maps to the rcad
//!   `gcpnts::arc_length` over the whole domain.

use glam::DVec3;

use rcad_kernel::base::convert::curve_to_bspline;
use rcad_kernel::base::gcpnts::abscissa_point::arc_length;
use rcad_kernel::geom::{BSplineCurve3, BSplineSurface, Curve3, CurveEval};
use rcad_kernel::math::GeomAbsShape;

use super::section_law::SectionLaw;

/// OCCT Precision::Confusion().
const CONFUSION: f64 = 1.0e-12;

/// OCCT GeomFill_UniformSection (GeomFill_UniformSection.hxx L80-85):
/// mySection / myCurve / First / Last.
#[derive(Debug, Clone)]
pub struct UniformSection {
    /// OCCT handle(Geom_Curve) mySection — the copied original section.
    my_section: Curve3,
    /// OCCT handle(Geom_BSplineCurve) myCurve.
    my_curve: BSplineCurve3,
    /// The OCCT BSplineSurface() result built once (architecture
    /// difference: the trait returns Option<&BSplineSurface>, so the OCCT
    /// on-demand construction happens here).
    my_bs: Option<BSplineSurface>,
    /// OCCT double First.
    first: f64,
    /// OCCT double Last.
    last: f64,
}

impl UniformSection {
    /// OCCT GeomFill_UniformSection::GeomFill_UniformSection (L35-53).
    pub fn new(c: &Curve3, first_parameter: f64, last_parameter: f64) -> Self {
        // OCCT: mySection = down_cast<Geom_Curve>(C->Copy()).
        let my_section = c.clone();
        // OCCT: myCurve = down_cast<Geom_BSplineCurve>(C); if null,
        // myCurve = GeomConvert::CurveToBSplineCurve(C,
        // Convert_QuasiAngular).
        let my_curve = match c {
            Curve3::BSpline(bs) => bs.clone(),
            _ => {
                let mut bs = curve_to_bspline(c, 0);
                if bs.is_periodic {
                    // OCCT: int M = myCurve->Degree() / 2 + 1;
                    // myCurve->RemoveKnot(1, M, Precision::Confusion()).
                    let m = (bs.degree as i32) / 2 + 1;
                    let (knots, mults) = bs.knots_mults();
                    let flat_poles = flatten_poles(&bs);
                    let mut new_poles = vec![0.0f64; flat_poles.len()];
                    let mut new_knots: Vec<f64> = Vec::new();
                    let mut new_mults: Vec<i32> = Vec::new();
                    let done = rcad_kernel::math::bspl_lib::remove_knot(
                        1,
                        m,
                        bs.degree,
                        bs.is_periodic,
                        3,
                        &flat_poles,
                        &knots,
                        &mults,
                        &mut new_poles,
                        &mut new_knots,
                        &mut new_mults,
                        CONFUSION,
                    );
                    if done {
                        let nb_poles = new_mults.iter().sum::<i32>() as usize;
                        let unflat = unflatten_poles(&new_poles, nb_poles);
                        bs = BSplineCurve3::from_knots_mults(bs.degree, new_knots, new_mults, unflat);
                        bs.is_periodic = true;
                    }
                }
                bs
            }
        };
        let my_bs = Self::build_bspline_surface(&my_curve, first_parameter, last_parameter);
        UniformSection {
            my_section,
            my_curve,
            my_bs,
            first: first_parameter,
            last: last_parameter,
        }
    }

    /// OCCT BSplineSurface (L110-144) — the on-demand surface construction.
    fn build_bspline_surface(
        my_curve: &BSplineCurve3,
        first: f64,
        last: f64,
    ) -> Option<BSplineSurface> {
        let (_nb_poles, cur_knots_mults) = (my_curve.control_points.len(), my_curve.knots_mults());
        let (cur_knots, cur_mults) = cur_knots_mults;

        // OCCT: Poles(1, NbPoles, 1, 2) with Poles(ii, 1) = Poles(ii, 2) =
        // myCurve->Pole(ii); VKnots(1) = First; VKnots(2) = Last;
        // VMults.Init(2).
        let control_points: Vec<Vec<DVec3>> = my_curve
            .control_points
            .iter()
            .map(|p| vec![*p, *p])
            .collect();
        let weights: Vec<Vec<f64>> = my_curve
            .weights
            .iter()
            .map(|w| vec![*w, *w])
            .collect();
        // OCCT: VKnots(1) = First; VKnots(2) = Last; VMults.Init(2).
        let knots_v = [first, last];
        let mults_v = [2i32, 2];
        let mut flat_knots_v = Vec::with_capacity(4);
        for (k, m) in knots_v.iter().zip(mults_v.iter()) {
            for _ in 0..*m {
                flat_knots_v.push(*k);
            }
        }
        let mut flat_knots_u = Vec::new();
        for (k, m) in cur_knots.iter().zip(cur_mults.iter()) {
            for _ in 0..*m {
                flat_knots_u.push(*k);
            }
        }
        // OCCT: new Geom_BSplineSurface(Poles, UKnots, VKnots, UMults,
        // VMults, myCurve->Degree(), 1, myCurve->IsPeriodic()).
        Some(BSplineSurface {
            degree_u: my_curve.degree,
            degree_v: 1,
            knots_u: flat_knots_u,
            knots_v: flat_knots_v,
            control_points,
            weights,
            is_periodic_u: my_curve.is_periodic,
            is_periodic_v: false,
        })
    }
}

/// Flatten the 3d poles for the BSplCLib-style remove_knot call (pure
/// math/tool helper).
fn flatten_poles(bs: &BSplineCurve3) -> Vec<f64> {
    let mut flat = Vec::with_capacity(bs.control_points.len() * 3);
    for p in &bs.control_points {
        flat.extend_from_slice(&[p.x, p.y, p.z]);
    }
    flat
}

/// Unflatten the 3d poles (pure math/tool helper).
fn unflatten_poles(flat: &[f64], nb_poles: usize) -> Vec<DVec3> {
    (0..nb_poles)
        .map(|i| DVec3::new(flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]))
        .collect()
}

impl SectionLaw for UniformSection {
    /// OCCT D0 (L57-68) — Poles = myCurve->Poles(); Weights =
    /// myCurve->WeightsArray().
    fn d0(&self, _param: f64, poles: &mut [DVec3], weights: &mut [f64]) -> bool {
        for (ii, p) in self.my_curve.control_points.iter().enumerate() {
            poles[ii] = *p;
        }
        for (ii, w) in self.my_curve.weights.iter().enumerate() {
            weights[ii] = *w;
        }
        true
    }

    /// OCCT D1 (L70-86).
    fn d1(
        &self,
        _param: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        weights: &mut [f64],
        dweights: &mut [f64],
    ) -> bool {
        for (ii, p) in self.my_curve.control_points.iter().enumerate() {
            poles[ii] = *p;
        }
        for (ii, w) in self.my_curve.weights.iter().enumerate() {
            weights[ii] = *w;
        }
        // OCCT: gp_Vec V0(0, 0, 0); DPoles.Init(V0); DWeights.Init(0).
        for dp in dpoles.iter_mut() {
            *dp = DVec3::ZERO;
        }
        for dw in dweights.iter_mut() {
            *dw = 0.0;
        }
        true
    }

    /// OCCT D2 (L88-108).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        _param: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        weights: &mut [f64],
        dweights: &mut [f64],
        d2weights: &mut [f64],
    ) -> bool {
        for (ii, p) in self.my_curve.control_points.iter().enumerate() {
            poles[ii] = *p;
        }
        for (ii, w) in self.my_curve.weights.iter().enumerate() {
            weights[ii] = *w;
        }
        // OCCT: DPoles.Init(V0); DWeights.Init(0); D2Poles.Init(V0);
        // D2Weights.Init(0).
        for dp in dpoles.iter_mut() {
            *dp = DVec3::ZERO;
        }
        for dw in dweights.iter_mut() {
            *dw = 0.0;
        }
        for dp in d2poles.iter_mut() {
            *dp = DVec3::ZERO;
        }
        for dw in d2weights.iter_mut() {
            *dw = 0.0;
        }
        true
    }

    /// OCCT BSplineSurface (L110-144) — the surface built by the ctor.
    fn bspline_surface(&self) -> Option<&BSplineSurface> {
        self.my_bs.as_ref()
    }

    /// OCCT SectionShape (L146-151).
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize) {
        *nb_poles = self.my_curve.control_points.len();
        let (knots, _) = self.my_curve.knots_mults();
        *nb_knots = knots.len();
        *degree = self.my_curve.degree;
    }

    /// OCCT Knots (L153-158).
    fn knots(&self, t_knots: &mut [f64]) {
        let (knots, _) = self.my_curve.knots_mults();
        t_knots[..knots.len()].copy_from_slice(&knots);
    }

    /// OCCT Mults (L161-167).
    fn mults(&self, t_mults: &mut [i32]) {
        let (_, mults) = self.my_curve.knots_mults();
        t_mults[..mults.len()].copy_from_slice(&mults);
    }

    /// OCCT IsRational (L169-175).
    fn is_rational(&self) -> bool {
        self.my_curve.is_rational()
    }

    /// OCCT IsUPeriodic (L177-183) — myCurve->IsPeriodic().
    fn is_u_periodic(&self) -> bool {
        self.my_curve.is_periodic
    }

    /// OCCT IsVPeriodic (L185-191).
    fn is_v_periodic(&self) -> bool {
        true
    }

    /// OCCT NbIntervals (L194-200).
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Intervals (L202-211).
    fn intervals(&self, t: &mut [f64], _s: GeomAbsShape) {
        let n = t.len();
        t[0] = self.first;
        t[n - 1] = self.last;
    }

    /// OCCT SetInterval (L213-219) — "Ne fait Rien".
    fn set_interval(&mut self, _first: f64, _last: f64) {}

    /// OCCT GetInterval (L221-227).
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        *first = self.first;
        *last = self.last;
    }

    /// OCCT GetDomain (L229-235).
    fn get_domain(&self, first: &mut f64, last: &mut f64) {
        *first = self.first;
        *last = self.last;
    }

    /// OCCT GetTolerance (L237-248).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, _angle_tol: f64, tol3d: &mut [f64]) {
        for value in tol3d.iter_mut() {
            *value = surf_tol;
        }
        let n = tol3d.len();
        if bound_tol < surf_tol {
            tol3d[0] = bound_tol;
            tol3d[n - 1] = bound_tol;
        }
    }

    /// OCCT BarycentreOfSurf (L250-266).
    fn barycentre_of_surf(&self) -> DVec3 {
        let curve = Curve3::BSpline(self.my_curve.clone());
        let mut u = curve.default_domain()[0];
        let mut bary = DVec3::ZERO;
        let delta = (self.my_curve.last_parameter() - u) / 20.0;
        for _ii in 0..=20 {
            bary += curve.point_at(u);
            u += delta;
        }
        bary /= 21.0;
        bary
    }

    /// OCCT MaximalSection (L268-272) — GCPnts_AbscissaPoint::Length over
    /// the whole section.
    fn maximal_section(&self) -> f64 {
        let curve = Curve3::BSpline(self.my_curve.clone());
        arc_length(&curve, curve.default_domain()[0], curve.default_domain()[1])
    }

    /// OCCT GetMinimalWeight (L274-277).
    fn get_minimal_weight(&self, weights: &mut [f64]) {
        for (ii, w) in self.my_curve.weights.iter().enumerate() {
            weights[ii] = *w;
        }
    }

    /// OCCT IsConstant (L279-283).
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        true
    }

    /// OCCT ConstantSection (L285-295).
    fn constant_section(&self) -> Curve3 {
        self.my_section.clone()
    }
}
