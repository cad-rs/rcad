//! Surface-surface intersection (IntSS).
//!
//! Covers analytic pairs:
//!
//! | Pair | Result |
//! |------|--------|
//! | Plane 脳 Plane | Line or parallel/coincident |
//! | Plane 脳 Sphere | Circle |
//! | Plane 脳 Cylinder | Circle / Ellipse / Lines |
//! | Plane 脳 Cone | Circle / Ellipse / Lines |
//! | Sphere 脳 Sphere | Circle (intersection plane 鈯?line-of-centres) |
//! | Sphere 脳 Cylinder | Circle (axis 鈯?case) / Numeric |
//! | Cylinder 脳 Cylinder | Ellipse (axes 鈭?, numeric otherwise |
//! | Cylinder 脳 Cone | Circle (coaxial) / Numeric |
//! | Sphere 脳 Cone | Circle (apex-centred case) / Numeric |
//! | Cone 脳 Cone | Circle (coaxial) / Numeric |
//! | Everything else | Numeric polylines via marching |
//!
//! The numerical branch uses a geometric **tolerance band** (cell scale and optional
//! floor 鈮?[`crate::tolerance::TOLERANCE_ABS`] by default ([`numeric_intss_with_domains`] /
//! [`numeric_intss_with_density`] with `geom_tol_floor: None`) unless a tighter bound is injected
//! via [`intersect_surfaces_with_density_tol`]): XOR edge
//! topology, chord minima for grazing contacts, split checks when
//! both edge corners lie inside the band, and closest-approach seeding when crossings
//! are still empty. Seeds are tightened with [`project_onto_intersection_tol`](crate::inttools::marching::project_onto_intersection_tol).
//!
//! Analogous to OCCT `GeomAPI_IntSS`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    BSplineCurve3, Circle3, ConicalSurface, Curve2d, Curve3, CylindricalSurface, Ellipse3,
    Hyperbola3, Line3, Parabola3, Plane, SphericalSurface, Surface3, CurveEval, SurfaceEval,
};

use crate::inttools::{
    cone_cone::{ConeConeResult, intersect_cone_cone, intersect_cone_cone_with_tolerance},
    cylinder_cone::{CylinderConeResult, intersect_cylinder_cone},
    cylinder_cylinder::{CylinderCylinderResult, intersect_cylinder_cylinder_with_tolerance, sample_perpendicular_offset_curves},
    cylinder_torus::CylinderTorusResult as CylinderTorusResultAlias,
    marching::{project_onto_intersection_tol, surface_implicit},
    pcurve_derive::{
        circle_pcurve_on_cone, circle_pcurve_on_cylinder, circle_pcurve_on_plane,
        circle_pcurve_on_sphere, ellipse_pcurve_on_cone, ellipse_pcurve_on_plane,
        fallback_pcurve_by_projection, line_pcurve_on_cone, line_pcurve_on_cylinder,
        line_pcurve_on_plane, polyline_pcurve_by_projection, sampled_pcurve_on_cone,
    },
    plane_cone::{PlaneConicalResult, intersect_plane_cone},
    plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
    plane_plane::{PlanePlaneResult, intersect_plane_plane},
    plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
    plane_torus::{PlaneTorusResult, intersect_plane_torus, intersect_plane_torus_skew},
    sphere_cone::{SphereConeResult, intersect_sphere_cone_with_tolerance},
    torus_cone::{TorusConeResult, intersect_torus_cone_with_tolerance},
    torus_torus::{TorusTorusResult, intersect_torus_torus_with_tolerance},
};
use rcad_kernel::projection::closest_point_on_surface;
use crate::tolerance::*;

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Public result types
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// A single intersection component between two surfaces.
#[derive(Debug, Clone)]
pub enum SurfaceCurve {
    /// An exact analytic circle.
    Circle(Circle3),
    /// An exact analytic ellipse.
    Ellipse(Ellipse3),
    /// An exact analytic line (infinite).
    Line(Line3),
    /// An exact analytic parabola.
    Parabola(Parabola3),
    /// An exact analytic hyperbola (single branch representation).
    Hyperbola(Hyperbola3),
    /// A tangent point (zero-dimensional contact).
    Point(DVec3),
    /// Numerically sampled polyline (fallback for non-analytic pairs).
    ///
    /// Polylines from `numeric_intss_impl` are automatically converted to
    /// `BSplineCurve` via `polyline_to_bspline`.  Skew-quartic and other
    /// analytic-fallback paths still produce raw polylines;
    /// TODO: convert those paths as well.
    Polyline(Vec<DVec3>),
    /// BSpline approximation of a polyline intersection curve.
    ///
    /// 鉁?OCCT-aligned: matches `GeomInt_IntSS::MakeBSpline` output.
    /// Provides C2 continuity and exact parameter evaluation for BRep edge construction.
    /// Boxed to keep the enum size manageable.
    BSplineCurve(Box<BSplineCurve3>),
}

/// One intersection result: 3D curve plus optional PCurves on each surface.
#[derive(Debug, Clone)]
pub struct SurfaceIntersectionResult {
    pub curve_3d: SurfaceCurve,
    /// PCurve on surface A (populated in Task 3+).
    pub pcurve_on_a: Option<Curve2d>,
    /// PCurve on surface B (populated in Task 3+).
    pub pcurve_on_b: Option<Curve2d>,
}

/// All intersection curves / components found between two surfaces.
#[derive(Debug, Clone, Default)]
pub struct SurfaceSurfaceIntersection {
    pub curves: Vec<SurfaceIntersectionResult>,
}

impl SurfaceSurfaceIntersection {
    pub fn is_empty(&self) -> bool {
        self.curves.is_empty()
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Main dispatch
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the intersection between two `Surface3` values.
///
/// Returns exact analytic curves where possible; falls back to numerical
/// polylines for unsupported surface-type combinations.
///
/// Uses the same grid density as [`intersect_surfaces_with_density`] (`n = 48`).
/// On the **numerical** fallback, applies at least [`TOLERANCE_ABS`](crate::tolerance::TOLERANCE_ABS) as the
/// geometric tolerance floor (sign-change threshold and polyline chaining). Use
/// [`intersect_surfaces_with_density`] if you need the numerical branch without
/// an explicit floor.
pub fn intersect_surfaces(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    intersect_surfaces_with_density_tol(s1, s2, 48, TOLERANCE_ABS)
}

/// Like [`intersect_surfaces`] but allows fuzzy geometric tolerance routing for
/// selected analytic cases, and strengthens the **numeric** marching branch when
/// a floor is imposed (see [`intersect_surfaces_with_density_tol`]).
///
/// Applies to `Cone x Cone`, `Cylinder x Cylinder`, `Sphere x Cylinder`,
/// `Sphere x Cone`, `Torus x Cone`, and `Torus x Torus` near-degenerate routing.
/// All other pairs use the same analytic dispatch as strict mode; the numerical
/// fallback then uses `fuzzy_tol.max(TOLERANCE_ABS)` as the geometric floor and
/// enables **tolerance intersection** refinements (interior edge minima, localized
/// closest-approach recovery, intersection projection) for near-tangent / narrow
/// contacts.
pub fn intersect_surfaces_with_tolerance(
    s1: &Surface3,
    s2: &Surface3,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    if fuzzy_tol <= 0.0 {
        return intersect_surfaces(s1, s2);
    }

    let numeric_floor = fuzzy_tol.max(TOLERANCE_ABS);

    match (s1, s2) {
        (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
            cylinder_x_cylinder_with_tolerance(c1, c2, fuzzy_tol)
        }
        (Surface3::Cone(k1), Surface3::Cone(k2)) => {
            cone_x_cone_with_tolerance(k1, k2, fuzzy_tol)
        }
        (Surface3::Sphere(s), Surface3::Cylinder(c))
        | (Surface3::Cylinder(c), Surface3::Sphere(s)) => {
            sphere_x_cylinder_with_tolerance(s, c, fuzzy_tol)
        }
        (Surface3::Sphere(s), Surface3::Cone(k))
        | (Surface3::Cone(k), Surface3::Sphere(s)) => {
            sphere_x_cone_with_tolerance(s, k, fuzzy_tol)
        }
        (Surface3::Torus(t), Surface3::Cone(k))
        | (Surface3::Cone(k), Surface3::Torus(t)) => {
            torus_x_cone_with_tolerance(t, k, fuzzy_tol)
        }
        (Surface3::Torus(t1), Surface3::Torus(t2)) => {
            torus_x_torus_with_tolerance(t1, t2, fuzzy_tol)
        }
        _ => intersect_surfaces_with_density_tol(s1, s2, 48, numeric_floor),
    }
}

/// Like [`intersect_surfaces`] but lets the caller specify the grid density `n`
/// for the numerical fallback.  Analytic pairs always return exact results
/// regardless of `n`.  The numerical fallback uses an `n脳n` parameter-space
/// grid to find sign-change crossings.
///
/// Higher `n` gives more accurate intersection polylines at the cost of O(n虏)
/// work.  The default used by [`intersect_surfaces`] is `n = 48`.
///
/// When using [`intersect_surfaces_with_density`] (caller does not supply `geom_tol_floor`),
/// numerical marching resolves to at least [`TOLERANCE_ABS`](crate::tolerance::TOLERANCE_ABS)
/// inside [`numeric_intss_impl`] (aligned with [`intersect_surfaces`]' default density).
///
/// Prefer [`crate::tolerance::intss_geom_tol_floor`] (or pairwise [`crate::tolerance::combined_linear_tol_*`])
/// via [`intersect_surfaces_with_density_tol`] when BRep topology tolerances are known.
pub fn intersect_surfaces_with_density(
    s1: &Surface3,
    s2: &Surface3,
    grid_n: usize,
) -> SurfaceSurfaceIntersection {
    intersect_surfaces_with_density_impl(s1, s2, grid_n, None)
}

/// Like [`intersect_surfaces_with_density`] but supplies a **minimum geometric tolerance**
/// for the numerical fallback (sign-change threshold and polyline chaining).
/// Analytic branches ignore `geom_tol_floor`.
pub fn intersect_surfaces_with_density_tol(
    s1: &Surface3,
    s2: &Surface3,
    grid_n: usize,
    geom_tol_floor: f64,
) -> SurfaceSurfaceIntersection {
    let floor = if geom_tol_floor.is_finite() && geom_tol_floor > 0.0 {
        Some(geom_tol_floor)
    } else {
        None
    };
    intersect_surfaces_with_density_impl(s1, s2, grid_n, floor)
}

fn intersect_surfaces_with_density_impl(
    s1: &Surface3,
    s2: &Surface3,
    grid_n: usize,
    geom_tol_floor: Option<f64>,
) -> SurfaceSurfaceIntersection {
    use Surface3::*;
    match (s1, s2) {
        // 鈹€鈹€ Plane 脳 * 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Plane(p1), Plane(p2)) => plane_x_plane(p1, p2),
        (Plane(p), Sphere(s)) | (Sphere(s), Plane(p)) => plane_x_sphere(p, s),
        (Plane(p), Cylinder(c)) | (Cylinder(c), Plane(p)) => plane_x_cylinder(p, c),
        (Plane(p), Cone(c)) | (Cone(c), Plane(p)) => plane_x_cone(p, c),
        (Plane(p), Torus(t)) | (Torus(t), Plane(p)) => torus_x_plane(t, p),

        // 鈹€鈹€ Sphere 脳 * 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Sphere(s1), Sphere(s2)) => sphere_x_sphere(s1, s2),
        (Sphere(s), Cylinder(c)) | (Cylinder(c), Sphere(s)) => sphere_x_cylinder(s, c),
        (Sphere(s), Cone(c)) | (Cone(c), Sphere(s)) => sphere_x_cone(s, c),
        (Sphere(s), Torus(t)) | (Torus(t), Sphere(s)) => torus_x_sphere(t, s),

        // 鈹€鈹€ Cylinder 脳 Cylinder 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Cylinder(c1), Cylinder(c2)) => cylinder_x_cylinder(c1, c2),

        // 鈹€鈹€ Cylinder 脳 Cone 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Cylinder(c), Cone(k)) | (Cone(k), Cylinder(c)) => cylinder_x_cone(c, k),

        // 鈹€鈹€ Cone 脳 Cone 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Cone(k1), Cone(k2)) => cone_x_cone(k1, k2),

        // 鈹€鈹€ Torus 脳 Cylinder 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Torus(t), Cylinder(c)) | (Cylinder(c), Torus(t)) => torus_x_cylinder(t, c),

        // 鈹€鈹€ Torus 脳 Cone 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Torus(t), Cone(k)) | (Cone(k), Torus(t)) => torus_x_cone(t, k),

        // 鈹€鈹€ Torus 脳 Torus 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        (Torus(t1), Torus(t2)) => torus_x_torus(t1, t2),

        // 鈹€鈹€ All others 鈫?numeric marching 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        _ => numeric_intss_with_density(s1, s2, grid_n, geom_tol_floor),
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Plane 脳 Plane
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn plane_x_plane(p1: &Plane, p2: &Plane) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_plane(p1, p2) {
        PlanePlaneResult::Line(l) => {
            let pca = line_pcurve_on_plane(&l, p1);
            let pcb = line_pcurve_on_plane(&l, p2);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlanePlaneResult::Coincident => {} // surfaces identical 鈥?infinite overlap
        PlanePlaneResult::Parallel => {}   // no intersection
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Plane 脳 Sphere
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn plane_x_sphere(p: &Plane, s: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_sphere(p, s) {
        PlaneSphereResult::Circle(c) => {
            let pca = circle_pcurve_on_plane(&c, p);
            let pcb = circle_pcurve_on_sphere(&c, s);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(c),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneSphereResult::TangentPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneSphereResult::NoIntersection => {}
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Plane 脳 Cylinder
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn plane_x_cylinder(p: &Plane, c: &CylindricalSurface) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_cylinder(p, c) {
        PlaneCylinderResult::Circle(circ) => {
            let pca = circle_pcurve_on_plane(&circ, p);
            let pcb = circle_pcurve_on_cylinder(&circ, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::Ellipse(e) => {
            let pca = ellipse_pcurve_on_plane(&e, p);
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Ellipse(e),
                &[0.0, TAU],
                &Surface3::Cylinder(*c),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Ellipse(e),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::TangentLine(l) => {
            let pca = line_pcurve_on_plane(&l, p);
            let pcb = line_pcurve_on_cylinder(&l, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneCylinderResult::TwoLines(l1, l2) => {
            let pca1 = line_pcurve_on_plane(&l1, p);
            let pcb1 = line_pcurve_on_cylinder(&l1, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l1),
                pcurve_on_a: Some(pca1),
                pcurve_on_b: Some(pcb1),
            });
            let pca2 = line_pcurve_on_plane(&l2, p);
            let pcb2 = line_pcurve_on_cylinder(&l2, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l2),
                pcurve_on_a: Some(pca2),
                pcurve_on_b: Some(pcb2),
            });
        }
        PlaneCylinderResult::NoIntersection => {}
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Plane 脳 Cone
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn plane_x_cone(p: &Plane, c: &ConicalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_cone(p, c) {
        PlaneConicalResult::Circle(circ) => {
            let pca = circle_pcurve_on_plane(&circ, p);
            let pcb = circle_pcurve_on_cone(&circ, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Ellipse(e) => {
            let pca = ellipse_pcurve_on_plane(&e, p);
            let pcb = ellipse_pcurve_on_cone(&e, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Ellipse(e),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::SingleLine(l) => {
            let pca = line_pcurve_on_plane(&l, p);
            let pcb = line_pcurve_on_cone(&l, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::TwoLines(l1, l2) => {
            let pca1 = line_pcurve_on_plane(&l1, p);
            let pcb1 = line_pcurve_on_cone(&l1, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l1),
                pcurve_on_a: Some(pca1),
                pcurve_on_b: Some(pcb1),
            });
            let pca2 = line_pcurve_on_plane(&l2, p);
            let pcb2 = line_pcurve_on_cone(&l2, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(l2),
                pcurve_on_a: Some(pca2),
                pcurve_on_b: Some(pcb2),
            });
        }
        PlaneConicalResult::Parabola(par) => {
            // Sample over a reasonable bounded domain for the PCurves
            let pca = fallback_pcurve_by_projection(&Curve3::Parabola(par), &[-20.0, 20.0], &Surface3::Plane(*p));
            let pcb = sampled_pcurve_on_cone(&Curve3::Parabola(par), &[-20.0, 20.0], c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Parabola(par),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Hyperbola(hyp) => {
            // Each branch sampled separately; use the principal branch domain
            let pca = fallback_pcurve_by_projection(&Curve3::Hyperbola(hyp), &[-10.0, 10.0], &Surface3::Plane(*p));
            let pcb = sampled_pcurve_on_cone(&Curve3::Hyperbola(hyp), &[-10.0, 10.0], c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Hyperbola(hyp),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneConicalResult::Point(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneConicalResult::NoIntersection => {}
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Sphere 脳 Sphere
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Two spheres intersect in a circle (or are tangent/disjoint).
///
/// The intersection circle lies on the radical plane, whose normal is the
/// line-of-centres direction.
fn sphere_x_sphere(s1: &SphericalSurface, s2: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    let d_vec = s2.center - s1.center;
    let d = d_vec.length();

    // Concentric spheres: coincident (same r) or no intersection
    if d < TOLERANCE_ABS {
        return out; // treat as no intersection (or coincident, infinite)
    }

    let r1 = s1.radius;
    let r2 = s2.radius;

    // Disjoint or one contains the other
    if d > r1 + r2 + TOLERANCE_ABS || d < (r1 - r2).abs() - TOLERANCE_ABS {
        return out;
    }

    let axis = d_vec / d;

    // Distance from s1.center to the intersection plane (radical plane)
    let a = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);

    // Tangent case
    let r_sq = r1 * r1 - a * a;
    if r_sq < -TOLERANCE_ABS {
        return out;
    }
    let r_circle = r_sq.max(0.0).sqrt();
    let center = s1.center + axis * a;

    if r_circle < TOLERANCE_ABS {
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Point(center),
            pcurve_on_a: None,
            pcurve_on_b: None,
        });
    } else {
        let circle = Circle3 {
            center,
            normal: axis,
            radius: r_circle,
        };
        let pca = circle_pcurve_on_sphere(&circle, s1);
        let pcb = circle_pcurve_on_sphere(&circle, s2);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(circle),
            pcurve_on_a: Some(pca),
            pcurve_on_b: Some(pcb),
        });
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Sphere 脳 Cylinder
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Sphere-cylinder intersection.
///
/// Analytic case: cylinder axis passes through the sphere centre 鈫?
/// two parallel circles (or one if tangent).
/// All other cases fall back to numerical marching.
fn sphere_x_cylinder(s: &SphericalSurface, c: &CylindricalSurface) -> SurfaceSurfaceIntersection {
    sphere_x_cylinder_with_tolerance(s, c, 0.0)
}

fn sphere_x_cylinder_with_tolerance(
    s: &SphericalSurface,
    c: &CylindricalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);
    // Project sphere centre onto cylinder axis
    let t = (s.center - c.origin).dot(c.axis);
    let foot = c.origin + c.axis * t;
    let d_perp = (s.center - foot).length();

    // If the sphere centre is on the cylinder axis, section planes 鈯?axis give
    // circles at heights where r_sphere(z)虏 = R_cylinder虏
    // r_sphere(z)虏 = R虏 - (z - z_c)虏 where z_c is the axial position of sphere center
    // Solve: R虏 - z虏 = r_cyl虏  (in local frame with sphere center as origin along axis)
    if d_perp < tol {
        // Sphere centre on axis 鈥?analytic circles
        let dz_sq = s.radius * s.radius - c.radius * c.radius;
        if dz_sq < -tol {
            // Sphere smaller than cylinder 鈥?no intersection if dz_sq < 0
            // Actually: sphere radius < cylinder radius means sphere inside cyl,
            // could still intersect if large enough. Recheck:
            // Points on intersection: distance from axis = c.radius AND on sphere.
            // If s.radius < c.radius: sphere never reaches cylinder surface 鈫?no intersect.
            return SurfaceSurfaceIntersection::default();
        }
        let mut out = SurfaceSurfaceIntersection::default();
        if dz_sq.abs() < tol {
            // Tangent 鈥?single circle at sphere center height
            let circle = Circle3 {
                center: s.center,
                normal: c.axis,
                radius: c.radius,
            };
            let pca = circle_pcurve_on_sphere(&circle, s);
            let pcb = circle_pcurve_on_cylinder(&circle, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circle),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        } else {
            let dz = dz_sq.sqrt();
            for &sign in &[1.0f64, -1.0] {
                let center = s.center + c.axis * (sign * dz);
                let circle = Circle3 {
                    center,
                    normal: c.axis,
                    radius: c.radius,
                };
                let pca = circle_pcurve_on_sphere(&circle, s);
                let pcb = circle_pcurve_on_cylinder(&circle, c);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circle),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        return out;
    }

    // General / off-axis case: try analytic quartic solver first.
    use super::sphere_cylinder::{SphereCylinderResult, intersect_sphere_cylinder_with_tolerance};
    match intersect_sphere_cylinder_with_tolerance(s, c, fuzzy_tol) {
        SphereCylinderResult::SkewQuartic(branches) => {
            let mut out = SurfaceSurfaceIntersection::default();
            let s_sph = Surface3::Sphere(*s);
            let s_cyl = Surface3::Cylinder(*c);
            for branch in &branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(branch, &s_sph);
                let pcb = polyline_pcurve_by_projection(branch, &s_cyl);
                // TODO: Try polyline_to_bspline conversion here for compact 3D
                // representation when branch has >= 4 points.
                // This SkewQuartic solver path currently bypasses
                // numeric_intss_impl's auto-BSpline conversion.
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch.clone()),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
            out
        }
        _ => {
            // Fall back to numeric marching if the quartic solver returns nothing.
            numeric_intss(&Surface3::Sphere(*s), &Surface3::Cylinder(*c))
        }
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Sphere x Cone
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Sphere-cone intersection.
///
/// Analytic case: sphere centre on cone axis -> circles at intersecting heights.
/// General case -> numerical.
fn sphere_x_cone(s: &SphericalSurface, c: &ConicalSurface) -> SurfaceSurfaceIntersection {
    sphere_x_cone_with_tolerance(s, c, 0.0)
}

fn sphere_x_cone_with_tolerance(
    s: &SphericalSurface,
    c: &ConicalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_sphere_cone_with_tolerance(s, c, fuzzy_tol) {
        SphereConeResult::NoIntersection => {}
        SphereConeResult::SingleCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Sphere(*s),
            );
            let pcb = circle_pcurve_on_cone(&circ, c);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        SphereConeResult::TwoCircles(circ1, circ2) => {
            for circ in [circ1, circ2] {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Sphere(*s),
                );
                let pcb = circle_pcurve_on_cone(&circ, c);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        SphereConeResult::TangentPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        SphereConeResult::General => {
            return numeric_intss(&Surface3::Sphere(*s), &Surface3::Cone(*c));
        }
        SphereConeResult::Polyline(branches) => {
            let s1 = Surface3::Sphere(*s);
            let s2 = Surface3::Cone(*c);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s1);
                let pcb = polyline_pcurve_by_projection(&branch, &s2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Cylinder 脳 Cylinder
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Cylinder-cylinder intersection.
///
/// Analytic case: parallel axes 鈫?ellipses (or circles if same radius and
/// same orientation).  General case (skew/crossing axes) 鈫?numerical.
fn cylinder_x_cylinder(
    c1: &CylindricalSurface,
    c2: &CylindricalSurface,
) -> SurfaceSurfaceIntersection {
    use crate::inttools::cylinder_cylinder::intersect_cylinder_cylinder;
    cylinder_x_cylinder_from_result(c1, c2, intersect_cylinder_cylinder(c1, c2))
}

fn cylinder_x_cylinder_with_tolerance(
    c1: &CylindricalSurface,
    c2: &CylindricalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    cylinder_x_cylinder_from_result(
        c1,
        c2,
        intersect_cylinder_cylinder_with_tolerance(c1, c2, fuzzy_tol),
    )
}

fn cylinder_x_cylinder_from_result(
    c1: &CylindricalSurface,
    c2: &CylindricalSurface,
    cc: CylinderCylinderResult,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match cc {
        CylinderCylinderResult::NoIntersection | CylinderCylinderResult::Coaxial => {}
        CylinderCylinderResult::OneGeneratorLine(line) => {
            let pca = line_pcurve_on_cylinder(&line, c1);
            let pcb = line_pcurve_on_cylinder(&line, c2);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Line(line),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        CylinderCylinderResult::TwoGeneratorLines(line1, line2) => {
            for line in [line1, line2] {
                let pca = line_pcurve_on_cylinder(&line, c1);
                let pcb = line_pcurve_on_cylinder(&line, c2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Line(line),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        CylinderCylinderResult::TwoEllipses(e1, e2) => {
            for e in [e1, e2] {
                let pca = fallback_pcurve_by_projection(&Curve3::Ellipse(e), &[0.0, TAU], &Surface3::Cylinder(*c1));
                let pcb = fallback_pcurve_by_projection(&Curve3::Ellipse(e), &[0.0, TAU], &Surface3::Cylinder(*c2));
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Ellipse(e),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        CylinderCylinderResult::TwoCircles(circ1, circ2) => {
            for circ in [circ1, circ2] {
                let pca = fallback_pcurve_by_projection(&Curve3::Circle(circ), &[0.0, TAU], &Surface3::Cylinder(*c1));
                let pcb = fallback_pcurve_by_projection(&Curve3::Circle(circ), &[0.0, TAU], &Surface3::Cylinder(*c2));
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        CylinderCylinderResult::PerpendicularOffsetCurves { ref cyl1, ref cyl2, dist } => {
            let branches = sample_perpendicular_offset_curves(cyl1, cyl2, dist, 16);
            let s_c1 = Surface3::Cylinder(*c1);
            let s_c2 = Surface3::Cylinder(*c2);
            for branch in branches {
                if branch.len() < 2 { continue; }
                let pca = polyline_pcurve_by_projection(&branch, &s_c1);
                let pcb = polyline_pcurve_by_projection(&branch, &s_c2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
        CylinderCylinderResult::General => {
            return numeric_intss(&Surface3::Cylinder(*c1), &Surface3::Cylinder(*c2));
        }
        CylinderCylinderResult::SkewQuartic(branches) => {
            let s_c1 = Surface3::Cylinder(*c1);
            let s_c2 = Surface3::Cylinder(*c2);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_c1);
                let pcb = polyline_pcurve_by_projection(&branch, &s_c2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Cylinder 脳 Cone
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Cylinder-cone intersection.
///
/// Analytic case: coaxial axes 鈫?single circle.
/// General case 鈫?numerical marching.
fn cylinder_x_cone(
    cyl: &CylindricalSurface,
    cone: &ConicalSurface,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_cylinder_cone(cyl, cone) {
        CylinderConeResult::NoIntersection => {}
        CylinderConeResult::CoaxialCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Cylinder(*cyl),
            );
            let pcb = circle_pcurve_on_cone(&circ, cone);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        CylinderConeResult::CoaxialTwoCircles(c1, c2) => {
            for circ in [c1, c2] {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl),
                );
                let pcb = circle_pcurve_on_cone(&circ, cone);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        CylinderConeResult::General => {
            return numeric_intss(&Surface3::Cylinder(*cyl), &Surface3::Cone(*cone));
        }
        CylinderConeResult::ParallelOffsetPolyline(branches) => {
            let s_cyl = Surface3::Cylinder(*cyl);
            let s_cone = Surface3::Cone(*cone);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_cyl);
                let pcb = polyline_pcurve_by_projection(&branch, &s_cone);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
        CylinderConeResult::SkewQuartic(branches) => {
            let s_cyl = Surface3::Cylinder(*cyl);
            let s_cone = Surface3::Cone(*cone);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_cyl);
                let pcb = polyline_pcurve_by_projection(&branch, &s_cone);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Cone 脳 Cone
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Cone-cone intersection.
///
/// Analytic case: coaxial cones 鈫?circle (or point if touching at apex).
/// General case 鈫?numerical marching.
fn cone_x_cone(
    k1: &ConicalSurface,
    k2: &ConicalSurface,
) -> SurfaceSurfaceIntersection {
    cone_x_cone_from_result(k1, k2, intersect_cone_cone(k1, k2))
}

fn cone_x_cone_with_tolerance(
    k1: &ConicalSurface,
    k2: &ConicalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    cone_x_cone_from_result(
        k1,
        k2,
        intersect_cone_cone_with_tolerance(k1, k2, fuzzy_tol),
    )
}

fn cone_x_cone_from_result(
    k1: &ConicalSurface,
    k2: &ConicalSurface,
    cc: ConeConeResult,
) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match cc {
        ConeConeResult::NoIntersection => {}
        ConeConeResult::Coaxial => {} // identical cones 鈥?infinite overlap
        ConeConeResult::CoaxialCircle(circ) => {
            let pca = circle_pcurve_on_cone(&circ, k1);
            let pcb = circle_pcurve_on_cone(&circ, k2);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        ConeConeResult::CoaxialPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        ConeConeResult::General => {
            return numeric_intss(&Surface3::Cone(*k1), &Surface3::Cone(*k2));
        }
        ConeConeResult::SkewQuartic(branches) => {
            let s_k1 = Surface3::Cone(*k1);
            let s_k2 = Surface3::Cone(*k2);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_k1);
                let pcb = polyline_pcurve_by_projection(&branch, &s_k2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Torus 脳 Plane
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Torus-plane intersection.
///
/// **Analytic case 鈥?plane 鈯?torus axis**:
///   If `plane.normal 鈭?torus.axis` (the plane is perpendicular to the torus
///   axis), the intersection consists of up to two circles coaxial with the
///   torus.  Let `d` be the signed distance from the torus center to the
///   plane along the axis.  The plane cuts the torus tube at the radii
///   `sqrt((R 卤 sqrt(r虏-d虏))虏)` where R is the major radius and r the minor.
///   More simply: the intersection circles have radii `R + sqrt(r虏-d虏)` and
///   `R - sqrt(r虏-d虏)` (the latter exists only when it is positive).
///   This simplifies to: for d虏 鈮?r虏, the two circle radii are
///   `R 卤 sqrt(r虏 - d虏)`, both centered at the torus center projected onto
///   the plane.
///
/// **Analytic case 鈥?plane 鈭?torus axis**:
///   If the plane is parallel to the torus axis, the tube center circle (radius
///   `R`) intersects the plane at up to two points, producing up to two
///   circles of radius `r` (the tube radius) in the plane.
///
/// **All other planes** fall back to numerical marching.
fn torus_x_plane(
    torus: &rcad_kernel::geom::ToroidalSurface,
    plane: &rcad_kernel::geom::Plane,
) -> SurfaceSurfaceIntersection {
    let axis = torus.axis.normalize();
    let normal = plane.normal.normalize();

    let cos_angle = axis.dot(normal).abs();
    const TOL: f64 = TOLERANCE_MESH_LEGACY;

    // 鈹€鈹€ Perpendicular to axis 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if (cos_angle - 1.0).abs() < TOL {
        // Existing perpendicular analytic handling.
        return torus_x_plane_perp(torus, plane, axis, normal);
    }

    // 鈹€鈹€ Parallel to axis (delegate to plane_torus module) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if cos_angle < TOL {
        return torus_x_plane_parallel(torus, plane);
    }

    // 鈹€鈹€ Skew 鈫?u-parameterized analytic solver 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let branches = intersect_plane_torus_skew(plane, torus);
    if !branches.is_empty() {
        let mut out = SurfaceSurfaceIntersection::default();
        let s_torus = Surface3::Torus(*torus);
        let s_plane = Surface3::Plane(*plane);
        for branch in branches {
            if branch.len() < 2 { continue; }
            let pca = polyline_pcurve_by_projection(&branch, &s_torus);
            let pcb = polyline_pcurve_by_projection(&branch, &s_plane);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Polyline(branch),
                pcurve_on_a: pca,
                pcurve_on_b: pcb,
            });
        }
        return out;
    }

    // 鈹€鈹€ Fallback 鈫?numeric marching 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    numeric_intss(
        &Surface3::Torus(*torus),
        &Surface3::Plane(*plane),
    )
}

/// Handle plane 鈭?torus axis via the analytic `plane_torus` module.
fn torus_x_plane_parallel(
    torus: &rcad_kernel::geom::ToroidalSurface,
    plane: &rcad_kernel::geom::Plane,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_torus(plane, torus) {
        PlaneTorusResult::NoIntersection => {}
        PlaneTorusResult::TangentCircle(c) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(c), &[0.0, TAU], &Surface3::Torus(*torus));
            let pcb = circle_pcurve_on_plane(&c, plane);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(c),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneTorusResult::TwoCircles(c1, c2) => {
            for circ in [c1, c2] {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ), &[0.0, TAU], &Surface3::Torus(*torus));
                let pcb = circle_pcurve_on_plane(&circ, plane);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        PlaneTorusResult::SkewPolyline(branches) => {
            let s_torus = Surface3::Torus(*torus);
            let s_plane = Surface3::Plane(*plane);
            for branch in branches {
                if branch.len() < 2 { continue; }
                let pca = polyline_pcurve_by_projection(&branch, &s_torus);
                let pcb = polyline_pcurve_by_projection(&branch, &s_plane);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
        }
        PlaneTorusResult::General => {
            return numeric_intss(&Surface3::Torus(*torus), &Surface3::Plane(*plane));
        }
    }
    out
}

/// Handle plane 鉄?torus axis (two coaxial circles).
fn torus_x_plane_perp(
    torus: &rcad_kernel::geom::ToroidalSurface,
    plane: &rcad_kernel::geom::Plane,
    axis: DVec3,
    normal: DVec3,
) -> SurfaceSurfaceIntersection {
    // Signed distance from torus center to the plane along the axis.
    let d = (plane.origin - torus.center).dot(normal) * normal.dot(axis).signum();
    let d_sq = d * d;
    let r_sq = torus.minor_radius * torus.minor_radius;

    if d_sq > r_sq + TOLERANCE_ABS {
        // Plane misses the torus tube
        return SurfaceSurfaceIntersection::default();
    }

    // Two intersection circles (when d虏=r虏 they merge into one)
    let delta = (r_sq - d_sq).max(0.0).sqrt();
    let r1 = torus.major_radius + delta;
    let r2 = torus.major_radius - delta;

    // Center of intersection circles: projection of torus center onto plane.
    let center_proj = torus.center + axis * d;

    // Build circle normal (same as plane normal, oriented to match plane)
    let circle_normal = normal;

    let mut out = SurfaceSurfaceIntersection::default();

    // Outer circle (r1 > 0 always for valid torus)
    if r1 > TOLERANCE_ABS {
        let pcurve_a = pcurve_for_torus_circle(torus, center_proj, r1, plane);
        let pcurve_b = crate::inttools::pcurve_derive::circle_pcurve_on_plane(
            &rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r1,
            },
            plane,
        );
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r1,
            }),
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
        });
    }

    // Inner circle (r2 > 0 only if delta < major_radius)
    if r2 > TOLERANCE_ABS {
        let pcurve_a = pcurve_for_torus_circle(torus, center_proj, r2, plane);
        let pcurve_b = crate::inttools::pcurve_derive::circle_pcurve_on_plane(
            &rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r2,
            },
            plane,
        );
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(rcad_kernel::geom::Circle3 {
                center: center_proj,
                normal: circle_normal,
                radius: r2,
            }),
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
        });
    }

    out
}

/// Compute a UV PCurve for a circle of intersection on a torus.
///
/// The circle has given center and radius in 3D, lying on the plane
/// perpendicular to the torus axis. Each point on the circle corresponds to
/// a unique major angle u (azimuth around the torus axis) and a fixed minor
/// angle v (around the tube). Returns a BSpline2 approximation.
fn pcurve_for_torus_circle(
    torus: &rcad_kernel::geom::ToroidalSurface,
    circle_center: DVec3,
    circle_radius: f64,
    _plane: &rcad_kernel::geom::Plane,
) -> rcad_kernel::geom::Curve2d {
    use rcad_kernel::projection::closest_point_on_surface;

    // Sample the circle in 3D and project each point onto the torus UV domain.
    let n = 33_usize;
    let u_ax = rcad_kernel::geom::any_perpendicular(torus.axis);
    let v_ax = torus.axis.cross(u_ax).normalize();

    let mut uv_pts: Vec<glam::DVec2> = (0..n)
        .map(|i| {
            let theta = 2.0 * std::f64::consts::PI * i as f64 / (n - 1) as f64;
            let p = circle_center + (u_ax * theta.cos() + v_ax * theta.sin()) * circle_radius;
            let proj = closest_point_on_surface(&Surface3::Torus(*torus), p, 16);
            glam::DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    // Unwrap u discontinuities across the 2蟺 seam.
    for i in 1..uv_pts.len() {
        let du = uv_pts[i].x - uv_pts[i - 1].x;
        if du > std::f64::consts::PI {
            for p in &mut uv_pts[i..] { p.x -= std::f64::consts::TAU; }
        } else if du < -std::f64::consts::PI {
            for p in &mut uv_pts[i..] { p.x += std::f64::consts::TAU; }
        }
    }

    rcad_kernel::fit::interpolate_points_2d(&uv_pts)
        .map(rcad_kernel::geom::Curve2d::BSpline)
        .unwrap_or_else(|_| rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: uv_pts[0],
            direction: glam::DVec2::X,
        }))
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Torus 脳 Sphere
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Torus-sphere intersection.
///
/// **Analytic case 鈥?sphere centre on torus axis**:
///   By rotational symmetry, the intersection consists of circles at heights
///   where the sphere's cross-section radius equals the torus tube's cross-section.
///   Solve `(d_perp - R)虏 + h虏 = r虏` (torus) and `d_perp虏 + h虏 = R_s虏` (sphere)
///   via root-finding on `f(z) = torus_radius(z) - sphere_radius(z)`.
///
/// **All other cases** fall back to numerical marching.
fn torus_x_sphere(
    torus: &rcad_kernel::geom::ToroidalSurface,
    sphere: &SphericalSurface,
) -> SurfaceSurfaceIntersection {
    let axis = torus.axis.normalize();

    // Project sphere center onto torus axis
    let t = (sphere.center - torus.center).dot(axis);
    let foot = torus.center + axis * t;
    let d_perp = (sphere.center - foot).length();

    // On-axis: sphere center on torus axis 鈫?circles (existing path)
    if d_perp < TOLERANCE_ABS {
        return torus_x_sphere_on_axis(torus, sphere, axis);
    }

    // Off-axis: try analytic torus-parameterized solver before numeric
    let st_result =
        super::sphere_torus::intersect_skew_sphere_torus(sphere, torus);
    if !st_result.is_empty() {
        let mut out = SurfaceSurfaceIntersection::default();
        let s_torus = Surface3::Torus(*torus);
        let s_sphere = Surface3::Sphere(*sphere);
        for branch in st_result {
            if branch.len() < 2 {
                continue;
            }
            let pca = polyline_pcurve_by_projection(&branch, &s_torus);
            let pcb = polyline_pcurve_by_projection(&branch, &s_sphere);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Polyline(branch),
                pcurve_on_a: pca,
                pcurve_on_b: pcb,
            });
        }
        if !out.curves.is_empty() {
            return out;
        }
    }

    numeric_intss(&Surface3::Torus(*torus), &Surface3::Sphere(*sphere))
}

#[allow(non_snake_case)]
fn torus_x_sphere_on_axis(
    torus: &rcad_kernel::geom::ToroidalSurface,
    sphere: &SphericalSurface,
    axis: DVec3,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let R = torus.major_radius;
    let r = torus.minor_radius;
    let R_s = sphere.radius;

    // In the plane through the axis: torus is a circle of radius r centered at (R, 0),
    // sphere is a circle of radius R_s centered at (0, z_s) where z_s is sphere center's
    // axial offset from torus center.
    // Actually in local coords: torus tube centerline is at distance R from axis,
    // tube radius = r. Sphere center is on axis at height z_s from torus center.
    let z_s = (sphere.center - torus.center).dot(axis);

    // Find intersection of two circles in the (蟻, z) half-plane:
    // Torus tube circle: (蟻 - R)虏 + (z - 0)虏 = r虏  (tube center at (R, 0))
    // Sphere cross-section: 蟻虏 + (z - z_s)虏 = R_s虏  (sphere center at (0, z_s))
    //
    // We need to find (蟻, z) where both are satisfied, with 蟻 > 0.
    // From sphere: 蟻虏 = R_s虏 - (z - z_s)虏
    // Substitute into torus: (sqrt(R_s虏 - (z-z_s)虏) - R)虏 + z虏 = r虏
    //
    // Sample z and find sign changes of f(z) = torus_residual(z).
    let mut out = SurfaceSurfaceIntersection::default();
    let n = 128usize;
    let z_lo = z_s - R_s;
    let z_hi = z_s + R_s;
    let mut prev_f = f64::NAN;
    let mut prev_z = 0.0f64;

    for i in 0..=n {
        let z = z_lo + (z_hi - z_lo) * i as f64 / n as f64;
        let dz_sphere = z - z_s;
        let rho_s_sq = R_s * R_s - dz_sphere * dz_sphere;
        if rho_s_sq < 0.0 {
            prev_f = f64::NAN;
            prev_z = z;
            continue;
        }
        let rho_s = rho_s_sq.sqrt();
        // Residual: distance from (rho_s, z) to torus tube circle center (R, 0)
        let f = (rho_s - R).powi(2) + z * z - r * r;

        if !prev_f.is_nan() && prev_f * f < 0.0 {
            let mut lo = prev_z;
            let mut hi = z;
            for _ in 0..64 {
                let mid = (lo + hi) * 0.5;
                let dm = mid - z_s;
                let rm = (R_s * R_s - dm * dm).max(0.0).sqrt();
                let fm = (rm - R).powi(2) + mid * mid - r * r;
                let flo = ((R_s * R_s - (lo - z_s).powi(2)).max(0.0).sqrt() - R).powi(2) + lo * lo - r * r;
                if fm * flo < 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let z_sol = (lo + hi) * 0.5;
            let dz = z_sol - z_s;
            let rho_sol = (R_s * R_s - dz * dz).max(0.0).sqrt();
            if rho_sol > TOLERANCE_ABS {
                let center = torus.center + axis * z_sol;
                let circle = Circle3 {
                    center,
                    normal: axis,
                    radius: rho_sol,
                };
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[0.0, TAU],
                    &Surface3::Torus(*torus),
                );
                let pcb = circle_pcurve_on_sphere(&circle, sphere);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circle),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        prev_f = f;
        prev_z = z;
    }

    if out.curves.is_empty() {
        numeric_intss(&Surface3::Torus(*torus), &Surface3::Sphere(*sphere))
    } else {
        out
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Torus 脳 Cylinder
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Torus-cylinder intersection.
///
/// **Analytic case 鈥?cylinder axis = torus axis (coaxial)**:
///   Intersection consists of circles at heights where the torus tube
///   cross-section meets the cylinder radius. Solve `(R_cyl - R)虏 + h虏 = r虏`.
///
/// **All other cases** fall back to numerical marching.
fn torus_x_cylinder(
    torus: &rcad_kernel::geom::ToroidalSurface,
    cyl: &CylindricalSurface,
) -> SurfaceSurfaceIntersection {
    let t_axis = torus.axis.normalize();
    let c_axis = cyl.axis.normalize();
    let cross = t_axis.cross(c_axis);
    let sin_angle = cross.length();

    let delta = cyl.origin - torus.center;
    let d_perp = (delta - t_axis * delta.dot(t_axis)).length();

    // Coaxial: same axis line
    if sin_angle < TOLERANCE_ANG && d_perp < TOLERANCE_ABS {
        return torus_x_cylinder_coaxial(torus, cyl, t_axis);
    }

    // Skew: try quartic-based analytic solver
    let ct_result = super::cylinder_torus::intersect_cylinder_torus_with_tolerance(cyl, torus, 0.0);
    if let CylinderTorusResultAlias::SkewQuartic(branches) = &ct_result {
        let mut out = SurfaceSurfaceIntersection::default();
        let s_torus = Surface3::Torus(*torus);
        let s_cyl = Surface3::Cylinder(*cyl);
        for branch in branches {
            if branch.len() < 2 {
                continue;
            }
            let pca = polyline_pcurve_by_projection(branch, &s_torus);
            let pcb = polyline_pcurve_by_projection(branch, &s_cyl);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Polyline(branch.clone()),
                pcurve_on_a: pca,
                pcurve_on_b: pcb,
            });
        }
        return out;
    }

    numeric_intss(&Surface3::Torus(*torus), &Surface3::Cylinder(*cyl))
}

#[allow(non_snake_case)]
fn torus_x_cylinder_coaxial(
    torus: &rcad_kernel::geom::ToroidalSurface,
    cyl: &CylindricalSurface,
    axis: DVec3,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let R = torus.major_radius;
    let r = torus.minor_radius;
    let r_cyl = cyl.radius;

    // In the (蟻, z) plane: torus tube is circle of radius r at (R, 0).
    // Cylinder is vertical line at 蟻 = r_cyl.
    // Intersection: (r_cyl - R)虏 + h虏 = r虏  鉄? h = 卤sqrt(r虏 - (r_cyl - R)虏)
    let dr = r_cyl - R;
    let h_sq = r * r - dr * dr;

    let mut out = SurfaceSurfaceIntersection::default();
    if h_sq < -TOLERANCE_ABS {
        return out; // cylinder outside torus tube
    }

    let h = h_sq.max(0.0).sqrt();
    let heights = if h.abs() < TOLERANCE_ABS {
        vec![0.0f64]
    } else {
        vec![-h, h]
    };

    for &hz in &heights {
        let center = torus.center + axis * hz;
        let circle = Circle3 {
            center,
            normal: axis,
            radius: r_cyl,
        };
        let pca = fallback_pcurve_by_projection(
            &Curve3::Circle(circle),
            &[0.0, TAU],
            &Surface3::Torus(*torus),
        );
        let pcb = circle_pcurve_on_cylinder(&circle, cyl);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Circle(circle),
            pcurve_on_a: Some(pca),
            pcurve_on_b: Some(pcb),
        });
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Torus x Cone
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Torus-cone intersection.
///
/// **Analytic case -- cone apex on torus axis, cone axis = torus axis**:
///   By rotational symmetry, intersections are circles. Solve
///   torus tube equation vs cone surface in the (rho, z) half-plane.
///
/// **All other cases** fall back to numerical marching.
fn torus_x_cone(
    torus: &rcad_kernel::geom::ToroidalSurface,
    cone: &ConicalSurface,
) -> SurfaceSurfaceIntersection {
    torus_x_cone_with_tolerance(torus, cone, 0.0)
}

fn torus_x_cone_with_tolerance(
    torus: &rcad_kernel::geom::ToroidalSurface,
    cone: &ConicalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_torus_cone_with_tolerance(torus, cone, fuzzy_tol) {
        TorusConeResult::NoIntersection => {}
        TorusConeResult::SingleCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*torus),
            );
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Cone(*cone),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        TorusConeResult::TwoCircles(circ1, circ2) => {
            for circ in [circ1, circ2] {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Torus(*torus),
                );
                let pcb = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone),
                );
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        TorusConeResult::TangentCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*torus),
            );
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Cone(*cone),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        TorusConeResult::SkewQuartic(branches) => {
            let s_torus = Surface3::Torus(*torus);
            let s_cone = Surface3::Cone(*cone);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_torus);
                let pcb = polyline_pcurve_by_projection(&branch, &s_cone);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
            return out;
        }
        TorusConeResult::General => {
            return numeric_intss(&Surface3::Torus(*torus), &Surface3::Cone(*cone));
        }
    }
    out
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Torus x Torus
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Torus-torus intersection.
///
/// **Analytic case -- coaxial tori (same axis line)**:
///   By rotational symmetry, intersections are circles at heights where
///   the torus tube circles (in the rho-z half-plane) meet.
///
/// **All other cases** fall back to numerical marching.
fn torus_x_torus(
    t1: &rcad_kernel::geom::ToroidalSurface,
    t2: &rcad_kernel::geom::ToroidalSurface,
) -> SurfaceSurfaceIntersection {
    torus_x_torus_with_tolerance(t1, t2, 0.0)
}

fn torus_x_torus_with_tolerance(
    t1: &rcad_kernel::geom::ToroidalSurface,
    t2: &rcad_kernel::geom::ToroidalSurface,
    fuzzy_tol: f64,
) -> SurfaceSurfaceIntersection {
    use std::f64::consts::TAU;
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_torus_torus_with_tolerance(t1, t2, fuzzy_tol) {
        TorusTorusResult::NoIntersection => {}
        TorusTorusResult::SingleCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*t1),
            );
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*t2),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        TorusTorusResult::TwoCircles(circ1, circ2) => {
            for circ in [circ1, circ2] {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Torus(*t1),
                );
                let pcb = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &Surface3::Torus(*t2),
                );
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Circle(circ),
                    pcurve_on_a: Some(pca),
                    pcurve_on_b: Some(pcb),
                });
            }
        }
        TorusTorusResult::TangentCircle(circ) => {
            let pca = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*t1),
            );
            let pcb = fallback_pcurve_by_projection(
                &Curve3::Circle(circ),
                &[0.0, TAU],
                &Surface3::Torus(*t2),
            );
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(circ),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        TorusTorusResult::SkewQuartic(branches) => {
            let s_t1 = Surface3::Torus(*t1);
            let s_t2 = Surface3::Torus(*t2);
            for branch in branches {
                if branch.len() < 2 {
                    continue;
                }
                let pca = polyline_pcurve_by_projection(&branch, &s_t1);
                let pcb = polyline_pcurve_by_projection(&branch, &s_t2);
                out.curves.push(SurfaceIntersectionResult {
                    curve_3d: SurfaceCurve::Polyline(branch),
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                });
            }
            return out;
        }
        TorusTorusResult::Coaxial => {
            // Identical tori - infinite overlap, return empty
        }
        TorusTorusResult::General => {
            return numeric_intss(&Surface3::Torus(*t1), &Surface3::Torus(*t2));
        }
    }
    out
}



/// Minimum of `approx_dist` along segment `pa`鈥揱pb` for `t 鈭?[0,1]` (ternary search).
///
/// Used when endpoints stay outside the distance band but a **narrow tangent** trough
/// exists along the grid edge 鈥?the classic topology-only `(da<th) XOR (db<th)` test misses this.
fn min_approx_dist_on_segment(
    pa: DVec3,
    pb: DVec3,
    approx_dist: &impl Fn(DVec3) -> f64,
) -> (f64, f64) {
    let mut a = 0.0_f64;
    let mut b = 1.0_f64;
    let mut best_t = 0.5;
    let mut best_d = approx_dist(pa.lerp(pb, best_t));

    for _ in 0..18 {
        if (b - a) < TOLERANCE_FLOAT_DEDUP {
            break;
        }
        let t1 = a + (b - a) / 3.0;
        let t2 = a + 2.0 * (b - a) / 3.0;
        let d1 = approx_dist(pa.lerp(pb, t1));
        let d2 = approx_dist(pa.lerp(pb, t2));
        if d1 < best_d {
            best_d = d1;
            best_t = t1;
        }
        if d2 < best_d {
            best_d = d2;
            best_t = t2;
        }
        if d1 < d2 {
            b = t2;
        } else {
            a = t1;
        }
    }

    (best_t, best_d)
}

fn dedup_points_spatial(pts: Vec<DVec3>, tol: f64) -> Vec<DVec3> {
    if pts.is_empty() || !tol.is_finite() || tol <= 0.0 {
        return pts;
    }
    let tol_sq = tol * tol;
    let mut kept: Vec<DVec3> = Vec::with_capacity(pts.len());
    for p in pts {
        if !kept
            .iter()
            .any(|q| (*q - p).length_squared() <= tol_sq)
        {
            kept.push(p);
        }
    }
    kept
}



/// Numerical surface-surface intersection via sign-change edge marching.
///
/// **Algorithm**:
/// 1. Sample `s1` on an N脳N grid; for each sample, compute approximate distance to `s2`
///    using a pre-sampled `s2` grid.
/// 2. Detect edges (horizontal or vertical between adjacent grid cells) where the
///    distance changes sign (one end < threshold, other 鈮?threshold or vice versa).
///    Linearly interpolate each crossing 鈫?candidate intersection points.
/// 3. BFS-greedy sort: start from any unvisited point, repeatedly extend the chain
///    by picking the nearest unvisited neighbor. Repeat until all points are visited.
///    This produces ordered polylines suitable for UV splitting.
fn numeric_intss(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    numeric_intss_with_density(s1, s2, 48, None)
}

/// Same as `numeric_intss` but with configurable grid density N.
pub fn numeric_intss_with_density(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
    geom_tol_floor: Option<f64>,
) -> SurfaceSurfaceIntersection {
    numeric_intss_impl(s1, s2, n, None, None, geom_tol_floor)
}

/// Same as `numeric_intss_with_density` but uses caller-supplied UV domains
/// for s1 and s2 instead of `default_domain()`.  Pass `None` for either to
/// use the surface's own default domain (with infinite-domain clamping).
///
/// `geom_tol_floor` (when `Some`) lower-bounds the sign-change threshold
/// (`max(cell_scale脳2, floor)`) and the polyline chain-join tolerances so
/// coarse grids still admit intersections on models with larger face/edge tolerances.
pub fn numeric_intss_with_domains(
    s1: &Surface3,
    s2: &Surface3,
    n: usize,
    dom1_override: Option<[f64; 4]>,
    dom2_override: Option<[f64; 4]>,
    geom_tol_floor: Option<f64>,
) -> SurfaceSurfaceIntersection {
    numeric_intss_impl(s1, s2, n, dom1_override, dom2_override, geom_tol_floor)
}

include!("extra1.rs");
include!("extra2.rs");
#[cfg(test)]
mod tests {
    include!("tests_inc.rs");
}
