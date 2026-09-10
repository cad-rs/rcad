//! GAP carriers for the W1-4 batch (ShapeConstruct_ProjectCurveOnSurface).
//!
//! These stand-ins cover dependencies owned by *other, not-yet-landed* TKShHealing
//! batches.  Each carrier names its OCCT anchor and the owning batch, keeps
//! OCCT's failure path where the real algorithm is missing, and is replaced
//! wholesale when the owning batch lands:
//!
//! - `ShapeAnalysisSurface` — OCCT `ShapeAnalysis_Surface`
//!   (ShapeAnalysis_Surface.cxx, 1,946 LOC) — W2 docket row.
//! - `shape_analysis_curve_project` / `_next_project` — OCCT
//!   `ShapeAnalysis_Curve::Project/NextProject` (ShapeAnalysis_Curve.cxx) —
//!   W2 docket row.
//! - `int_curve_int_conic_conic_lin_lin` — OCCT `IntCurve_IntConicConic`
//!   (TKG2D IntCurve_IntConicConic.gxx) bridged onto the rcad IntAna2d port.
//! - `shape_analysis_adjust_by_period` / `_to_period` — exact re-hosts of the
//!   OCCT `ShapeAnalysis` statics (ShapeAnalysis.cxx L44-59 / L66-69), same
//!   pattern as `feat/loc_ope_wires_on_shape_b.rs` L165-183.
//! - `surface_u_period` / `surface_v_period` — pure-math re-hosts of the
//!   `Geom_Surface` period queries over `Surface3` (arch. diff.).
//! - `surface_is_cn_u` / `surface_is_cn_v` — OCCT
//!   `Geom_Surface::IsCNu/IsCNv` re-hosted over `Surface3` (arch. diff.).

use glam::{DVec2, DVec3};
use rcad_kernel::base::geom_api::project::{closest_point_on_curve_range, closest_point_on_surface};
use rcad_kernel::base::int_ana2d::AnaIntersection2d;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve_b::GeomSurfaceHandle;
use rcad_kernel::geom::{Curve3, Line2d, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox;
use std::cell::Cell;
use std::sync::Arc;

// OCCT ShapeAnalysis_Surface (ShapeAnalysis_Surface.cxx) — W2 docket row.
//
// GAP carrier reduced to the members consumed by
// ShapeConstruct_ProjectCurveOnSurface (OCCT
// ShapeConstruct_ProjectCurveOnSurface.cxx).  The projection members are
// re-hosted onto the rcad GeomAPI_ProjectPointOnSurf port
// (base/geom_api/project.rs); the singularity / degeneracy / iso machinery is
// reduced to OCCT's "nothing found" failure path.  The W2 1:1 translation of
// ShapeAnalysis_Surface replaces this carrier.
pub struct ShapeAnalysisSurface {
    /// OCCT mySurf.
    surf: Surface3,
    /// OCCT myGap — gap of the last projection.  OCCT mutates it from const
    /// methods; rcad uses interior mutability to keep the same signatures.
    my_gap: Cell<f64>,
    /// OCCT myBoxes[4] (GetBoxUF/UL/VF/VL), built from the surface corners.
    my_boxes: [BndBox; 4],
}

impl ShapeAnalysisSurface {
    /// OCCT ShapeAnalysis_Surface(Surf).
    pub fn new(the_surf: Surface3) -> Self {
        // OCCT builds the iso boxes lazily (BoxUF/BoxUL/BoxVF/BoxVL); the
        // carrier fills a conservative box from the four domain corners.
        let [uf, ul, vf, vl] = surface_bounds(&the_surf);
        let mut my_boxes = [BndBox::new(), BndBox::new(), BndBox::new(), BndBox::new()];
        let corners = [(uf, vf), (ul, vf), (uf, vl), (ul, vl)];
        for (bi, (u, v)) in corners.iter().enumerate() {
            my_boxes[bi].add_point(the_surf.point_at(*u, *v));
        }
        ShapeAnalysisSurface {
            surf: the_surf,
            my_gap: Cell::new(0.0),
            my_boxes,
        }
    }

    /// OCCT Surface() — the underlying Geom_Surface.
    pub fn surface(&self) -> &Surface3 {
        &self.surf
    }

    /// OCCT ValueOfUV(P3D, preci) (ShapeAnalysis_Surface.cxx L412-459):
    /// nearest-UV projection.  Carrier: the rcad closest-point port; OCCT's
    /// Extrema_ExtPS + Newton refinement (W2) may differ numerically.
    pub fn value_of_uv(&self, the_p3d: DVec3, _the_preci: f64) -> DVec2 {
        let proj = closest_point_on_surface(&self.surf, the_p3d, 64);
        self.my_gap.set(proj.distance);
        DVec2::new(proj.params.0, proj.params.1)
    }

    /// OCCT NextValueOfUV(uvPrev, P3D, preci, maxgap) (ShapeAnalysis_Surface.cxx
    /// L466-533): projection warm-started from `uv_prev`.  Carrier: full
    /// re-projection (the hint-based Newton restart is W2 scope).
    pub fn next_value_of_uv(
        &self,
        _the_uv_prev: DVec2,
        the_p3d: DVec3,
        _the_preci: f64,
        _the_maxgap: f64,
    ) -> DVec2 {
        self.value_of_uv(the_p3d, 0.0)
    }

    /// OCCT Gap() — the gap of the last projection.
    pub fn gap(&self) -> f64 {
        self.my_gap.get()
    }

    /// OCCT IsUClosed(preci) (ShapeAnalysis_Surface.cxx L107-127).
    /// Carrier: the plain Geom IsUClosed query (the BSpline tolerance-aware
    /// refinement is W2 scope).
    pub fn is_u_closed(&self, _the_preci: f64) -> bool {
        self.surf.is_u_closed()
    }

    /// OCCT IsVClosed(preci) (ShapeAnalysis_Surface.cxx L129-149).
    pub fn is_v_closed(&self, _the_preci: f64) -> bool {
        self.surf.is_v_closed()
    }

    /// OCCT Bounds(U1, U2, V1, V2) — surface parameter bounds.
    pub fn bounds(&self) -> [f64; 4] {
        surface_bounds(&self.surf)
    }

    /// OCCT UIso(U) (delegates to Geom_Surface::UIso).  GAP: iso-curve
    /// extraction is W2 scope; None keeps OCCT's null-iso failure path (the
    /// caller skips the isoline).
    pub fn u_iso(&self, _the_u: f64) -> Option<Curve3> {
        None
    }

    /// OCCT VIso(V) — see `u_iso`.
    pub fn v_iso(&self, _the_v: f64) -> Option<Curve3> {
        None
    }

    /// OCCT GetBoxUF().
    pub fn get_box_uf(&self) -> &BndBox {
        &self.my_boxes[0]
    }

    /// OCCT GetBoxUL().
    pub fn get_box_ul(&self) -> &BndBox {
        &self.my_boxes[1]
    }

    /// OCCT GetBoxVF().
    pub fn get_box_vf(&self) -> &BndBox {
        &self.my_boxes[2]
    }

    /// OCCT GetBoxVL().
    pub fn get_box_vl(&self) -> &BndBox {
        &self.my_boxes[3]
    }

    /// OCCT NbSingularities(preci) (ShapeAnalysis_Surface.cxx L537-...).
    /// GAP: singularity detection is W2 scope; 0 keeps OCCT's "no
    /// singularity" path.
    pub fn nb_singularities(&self, _the_preci: f64) -> i32 {
        0
    }

    /// OCCT Singularity(...) (ShapeAnalysis_Surface.cxx L560-...).  GAP:
    /// returns false (no more singularities) — the caller's enumeration loop
    /// breaks exactly as it does in OCCT when the list is exhausted.
    #[allow(clippy::too_many_arguments)]
    pub fn singularity(
        &self,
        _the_num: i32,
        _the_preci: &mut f64,
        _the_p3d: &mut DVec3,
        _the_first_p2d: &mut DVec2,
        _the_last_p2d: &mut DVec2,
        _the_first_par: &mut f64,
        _the_last_par: &mut f64,
        _the_is_iso_line: &mut bool,
    ) -> bool {
        false
    }

    /// OCCT IsDegenerated(P, preci).  GAP: false keeps OCCT's "not
    /// degenerated" path.
    pub fn is_degenerated(&self, _the_p: DVec3, _the_preci: f64) -> bool {
        false
    }

    /// OCCT ProjectDegenerated(nbPnt, Points, p2d, preci, direct)
    /// (ShapeAnalysis_Surface.cxx L2xx).  GAP: no-op keeps OCCT's path when
    /// no point lies in a degenerate region.
    pub fn project_degenerated(
        &self,
        _the_nbr: i32,
        _the_points: &[DVec3],
        _the_p2d: &mut [DVec2],
        _the_preci: f64,
        _the_direct: bool,
    ) {
    }

    /// OCCT Adaptor3d() — the surface as an Adaptor3d_Surface for ProjLib.
    pub fn adaptor3d(&self) -> GeomSurfaceHandle {
        Arc::new(GeomSurfaceAdaptor::new(self.surf.clone()))
    }
}

/// Surface parameter bounds `[u1, u2, v1, v2]` (Geom_Surface::Bounds over
/// Surface3's default domain).
fn surface_bounds(the_surf: &Surface3) -> [f64; 4] {
    the_surf.default_domain()
}

// OCCT ShapeAnalysis::AdjustByPeriod(Val, ToVal, Period) (ShapeAnalysis.cxx
// L44-59) — exact re-host (the feat/loc_ope_wires_on_shape_b.rs precedent).
pub fn shape_analysis_adjust_by_period(the_val: f64, to_val: f64, period: f64) -> f64 {
    let diff = the_val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

// OCCT ShapeAnalysis::AdjustToPeriod(Val, ValMin, ValMax) (ShapeAnalysis.cxx
// L66-69) — exact re-host; returns the shift to add to `the_val`.
pub fn shape_analysis_adjust_to_period(the_val: f64, val_min: f64, val_max: f64) -> f64 {
    shape_analysis_adjust_by_period(the_val, 0.5 * (val_min + val_max), val_max - val_min)
}

/// OCCT ShapeAnalysis_Curve::Project(C3D, P3d, preci, proj, param, Cf, Cl)
/// (ShapeAnalysis_Curve.cxx) — W2 docket row.
///
/// GAP carrier: projects through the rcad GeomAPI_ProjectPointOnCurve-range
/// port (`closest_point_on_curve_range`) instead of OCCT's Extrema-based
/// walk.  Returns `(proj, param, dist)`; the caller compares `dist` against
/// `preci` exactly as in OCCT.
pub fn shape_analysis_curve_project(
    the_c3d: &Curve3,
    the_p3d: DVec3,
    _the_preci: f64,
    the_cf: f64,
    the_cl: f64,
) -> (DVec3, f64, f64) {
    let proj = closest_point_on_curve_range(the_c3d, the_p3d, the_cf, the_cl, 64);
    (proj.point, proj.param, proj.distance)
}

/// OCCT ShapeAnalysis_Curve::NextProject(upar, C3D, P3d, preci, proj, param,
/// Cf, Cl, Adjust) (ShapeAnalysis_Curve.cxx) — W2 docket row.
///
/// GAP carrier: OCCT re-projects warm-started from `upar`; the carrier does a
/// full re-projection in `[Cf, Cl]` (the `Adjust` flag is the single call
/// site's `false`).
pub fn shape_analysis_curve_next_project(
    _the_upar: f64,
    the_c3d: &Curve3,
    the_p3d: DVec3,
    _the_preci: f64,
    the_cf: f64,
    the_cl: f64,
) -> (DVec3, f64, f64) {
    shape_analysis_curve_project(the_c3d, the_p3d, 0.0, the_cf, the_cl)
}

/// OCCT IntCurve_IntConicConic(Conic1, Dom1, Conic2, Dom2, Toler1, Toler2)
/// (IntCurve_IntConicConic.gxx), bridged onto the rcad IntAna2d port
/// (`int_ana2d::AnaIntersection2d`).
///
/// The single call site (correctExtremity) intersects two `gp_Lin2d`s with
/// default-constructed (unbounded) `IntRes2d_Domain`s, so the bridge covers
/// the line-line case only.  Returns `Some(value)` for the first intersection
/// point (`Intersector.Done() && !Empty() → Point(1).Value()`), else None.
pub fn int_curve_int_conic_conic_lin_lin(
    the_conic1: &Line2d,
    the_conic2: &Line2d,
) -> Option<DVec2> {
    let mut intersector = AnaIntersection2d::new();
    intersector.perform_lin_lin(the_conic1, the_conic2);
    if intersector.is_done() && !intersector.is_empty() {
        Some(intersector.point(1).value())
    } else {
        None
    }
}

/// OCCT Geom_BSplineSurface::IsCNu(N) recomputed from rcad's flat knot
/// vector: the largest interior-knot multiplicity bounds the continuity
/// (`continuity = degree - maxmult`); analytic surfaces are CN infinite.
/// Architecture bridge over `Surface3`'s flat-knot B-spline encoding.
pub fn surface_is_cn_u(the_surf: &Surface3, the_n: usize) -> bool {
    surface_is_cn(the_surf, the_n, true)
}

/// OCCT Geom_Surface::IsCNv(N) — see `surface_is_cn_u`.
pub fn surface_is_cn_v(the_surf: &Surface3, the_n: usize) -> bool {
    surface_is_cn(the_surf, the_n, false)
}

fn surface_is_cn(the_surf: &Surface3, the_n: usize, the_u: bool) -> bool {
    match the_surf {
        // Analytic surfaces: OCCT returns Standard_True for any N.
        Surface3::Plane(_)
        | Surface3::Cylinder(_)
        | Surface3::Sphere(_)
        | Surface3::Cone(_)
        | Surface3::Torus(_)
        | Surface3::Ellipsoid(_) => true,
        // Trimmed / offset surfaces delegate to the basis surface
        // (Geom_RectangularTrimmedSurface / Geom_OffsetSurface IsCN*
        // semantics, reduced to the basis query).
        Surface3::Trimmed(s) => surface_is_cn(&s.basis, the_n, the_u),
        Surface3::Offset(s) => surface_is_cn(&s.basis, the_n, the_u),
        Surface3::BSpline(s) => {
            // OCCT Geom_BSplineSurface::IsCNu(N): continuity in a direction is
            // degree - max interior multiplicity; rcad stores flat knots, so
            // interior multiplicities are counted as runs of equal knots.
            let knots = if the_u { &s.knots_u } else { &s.knots_v };
            let degree = if the_u { s.degree_u } else { s.degree_v };
            let mut max_mult = 0usize;
            let mut i = 1usize;
            while i + 1 < knots.len() {
                let mut mult = 1usize;
                while i + mult < knots.len() && knots[i + mult] == knots[i] {
                    mult += 1;
                }
                // A run that reaches the last flat knot is the end knot, not
                // an interior one.
                if i + mult < knots.len() {
                    max_mult = max_mult.max(mult);
                }
                i += mult;
            }
            degree >= the_n + max_mult
        }
        Surface3::Bezier(s) => {
            // Geom_BezierSurface degree = pole count - 1 per direction.
            let degree_u = s.control_points.len().saturating_sub(1);
            let degree_v = s.control_points.first().map(|r| r.len()).unwrap_or(1).saturating_sub(1);
            let degree = if the_u { degree_u } else { degree_v };
            degree >= the_n
        }
        // Composite / procedural surfaces: conservative false (OCCT's C0 path).
        _ => false,
    }
}

/// OCCT Geom_ConicSurface/Geom_SphericalSurface::UPeriod — the U period of
/// the elementary periodic surfaces (pure-math re-host, arch. diff.).
pub fn surface_u_period(surf: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match surf {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => TAU,
        _ => 0.0,
    }
}

/// OCCT Geom_SphericalSurface/ToroidalSurface::VPeriod (2*PI sphere,
/// 2*minorRadius torus; pure-math re-host, arch. diff.).
pub fn surface_v_period(surf: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match surf {
        Surface3::Sphere(_) => TAU,
        Surface3::Torus(t) => 2.0 * t.minor_radius,
        _ => 0.0,
    }
}
