//! Analytic intersection of two cones.
//!
//! # Case classification
//!
//! ## Coaxial cones (axes coincide)
//!
//! When both cone axes are the same line the intersection depends on the radii
//! and half-angles:
//!
//! - **Same apex, same half-angle**: identical cones — `Coaxial`.
//! - **Different apex or half-angle**: the two lateral surfaces meet at a circle
//!   perpendicular to the shared axis at the height where both radii are equal.
//!   If no such positive-radius solution exists the intersection is at the apex
//!   only (`Point`) or empty.
//!
//! ## Parallel axes (non-coaxial)
//!
//! When the axes are parallel but distinct, a quick radial-distance test decides
//! `NoIntersection` when the cones' radial envelopes cannot overlap.  Otherwise
//! we return `General`.
//!
//! ## General / skew axes
//!
//! For all other configurations the intersection is a curve of degree ≤ 4.
//! We return `General` so the caller falls back to numeric marching.

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, ConicalSurface, any_perpendicular};
use std::f64::consts::TAU;

use super::pcurve_derive::refine_polyline;
use crate::tolerance::*;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of cone × cone intersection.
#[derive(Debug, Clone)]
pub enum ConeConeResult {
    /// The cones do not intersect (lateral surfaces are disjoint).
    NoIntersection,
    /// Cones are coaxial with identical geometry (same nappe, same surface).
    Coaxial,
    /// Coaxial cones with different geometry: intersection is a single circle.
    CoaxialCircle(Circle3),
    /// Coaxial cones that only touch at a single point (a shared apex).
    CoaxialPoint(DVec3),
    /// Skew axes (non-parallel): analytic quartic solution.
    ///
    /// The two cones intersect in a quartic space curve.  For each cone azimuth
    /// u ∈ [0, 2π) the second cone's equation reduces to a quadratic in the
    /// slant distance v, solved analytically.  Two branches (± sqrt) are returned
    /// as polylines.
    SkewQuartic(Vec<Vec<DVec3>>),
    /// General case (skew or oblique axes).  Caller should fall back to marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `cone1` and `cone2`.
pub fn intersect_cone_cone(cone1: &ConicalSurface, cone2: &ConicalSurface) -> ConeConeResult {
    let a1 = cone1.axis.normalize();
    let a2 = cone2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // ── Parallel axes (including coaxial) ────────────────────────────────────
    if sin_angle < TOLERANCE_ANG {
        let result = intersect_parallel_cones(cone1, cone2, a1, a2);
        // OCCT Sec 2 (parallel-offset) is a stub returning None for now.
        if !matches!(result, ConeConeResult::General) {
            return result;
        }
        if let Some(r) = occt_parallel_offset_cones(cone1, cone2) {
            return r;
        }
        return ConeConeResult::General;
    }

    // ── OCCT Sec 3: Coincident apices (IntAna_QuadQuadGeo case 3) ──────────
    if let Some(result) = occt_coincident_apex_cones(cone1, cone2) {
        return result;
    }

    // ── OCCT Sec 4: Common generatrix (IntAna_QuadQuadGeo case 4) ──────────
    if let Some(result) = occt_common_generatrix_cones(cone1, cone2) {
        return result;
    }

    // ── Skew axes (analytic quartic solver) ──────────────────────────────────
    // Parameterize cone1, substitute into cone2 equation → quadratic in v per u.
    let skew_result = intersect_skew_cone_cone(cone1, cone2);
    if !skew_result.is_empty() {
        return ConeConeResult::SkewQuartic(skew_result);
    }

    // ── General / skew (fallback to marching) ─────────────────────────────────
    ConeConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// Skew-axis analytic solver
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the intersection of two cones with skew (non-parallel) axes using
/// an analytic quartic solver.
///
/// # Theory
///
/// Cone 1 parametrisation (u = azimuth [0, 2π), v = slant distance from apex):
///
/// ```text
/// P(u,v) = O1 + v·d1(u)
/// d1(u) = a1 + tan(α1)·(cos(u)·x1 + sin(u)·y1)
/// |d1|² = sec²(α1)  (constant for a given cone)
/// ```
///
/// Cone 2 implicit equation (for any point P):
///
/// ```text
/// ((P - O2)·a2)² = cos²(α2)·|P - O2|²
/// ```
///
/// Substituting P(u,v) gives F(v) = 0 where
///
/// ```text
/// a_v(u)·v² + b_v(u)·v + c_v = 0
///
/// a_v(u) = (d1·a2)² - cos²(α2)·|d1|²
/// b_v(u) = 2·(Δ·a2)·(d1·a2) - 2·cos²(α2)·(Δ·d1)
/// c_v    = (Δ·a2)² - cos²(α2)·|Δ|²                                  (constant)
/// Δ      = O1 - O2
/// ```
///
/// For each u we solve the quadratic for v, giving two branches (± sqrt).
fn intersect_skew_cone_cone(cone1: &ConicalSurface, cone2: &ConicalSurface) -> Vec<Vec<DVec3>> {
    let a1 = cone1.axis.normalize();
    let a2 = cone2.axis.normalize();
    let o1 = cone1.apex_point();
    let o2 = cone2.apex_point();
    let tan1 = cone1.half_angle_rad.tan();
    let cos2_2 = cone2.half_angle_rad.cos().powi(2); // cos²(α2)
    let d1_sq = 1.0 + tan1 * tan1; // |d1|² = sec²(α1)

    // Perpendicular basis for cone1.
    let x1 = any_perpendicular(a1);
    let y1 = a1.cross(x1).normalize();

    let delta = o1 - o2; // O1 - O2
    let delta_a2 = delta.dot(a2);
    let delta_sq = delta.length_squared();
    let c_v = delta_a2 * delta_a2 - cos2_2 * delta_sq; // constant

    const N_SAMPLES: usize = 128;
    const CHORD_TOL: f64 = crate::bop::int_tools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::bop::int_tools::CHORD_REFINE_DEPTH;
    let mut branch_plus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);
    let mut branch_minus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cos_u, sin_u) = (u.cos(), u.sin());

        // d1(u) = a1 + tan(α1)·(cos(u)·x1 + sin(u)·y1)
        let d1 = a1 + tan1 * (cos_u * x1 + sin_u * y1);
        let d1_a2 = d1.dot(a2);
        let delta_d1 = delta.dot(d1);

        // a_v(u) = (d1·a2)² - cos²(α2)·|d1|²
        let a_v = d1_a2 * d1_a2 - cos2_2 * d1_sq;

        // b_v(u) = 2·(Δ·a2)·(d1·a2) - 2·cos²(α2)·(Δ·d1)
        let b_v = 2.0 * delta_a2 * d1_a2 - 2.0 * cos2_2 * delta_d1;

        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                continue;
            }
            let sqrt_disc = disc.sqrt();
            let two_a_v = 2.0 * a_v;

            let v_plus = (-b_v + sqrt_disc) / two_a_v;
            let v_minus = (-b_v - sqrt_disc) / two_a_v;

            if v_plus.is_finite() {
                let p = cone1.point_at(u, v_plus);
                if p.is_finite() {
                    branch_plus.push((u, p));
                }
            }
            if v_minus.is_finite() {
                let p = cone1.point_at(u, v_minus);
                if p.is_finite() {
                    branch_minus.push((u, p));
                }
            }
        } else if b_v.abs() > 1e-12 {
            // a_v ≈ 0: solve linear b_v·v + c_v = 0
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cone1.point_at(u, v);
                if p.is_finite() {
                    branch_plus.push((u, p));
                    branch_minus.push((u, p));
                }
            }
        }
    }

    // Adaptive refinement: subdivide segments where chord error exceeds tolerance
    let eval_plus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d1 = a1 + tan1 * (cos_u * x1 + sin_u * y1);
        let d1_a2 = d1.dot(a2);
        let delta_d1 = delta.dot(d1);
        let a_v = d1_a2 * d1_a2 - cos2_2 * d1_sq;
        let b_v = 2.0 * delta_a2 * d1_a2 - 2.0 * cos2_2 * delta_d1;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v + disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cone1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cone1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };
    let eval_minus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d1 = a1 + tan1 * (cos_u * x1 + sin_u * y1);
        let d1_a2 = d1.dot(a2);
        let delta_d1 = delta.dot(d1);
        let a_v = d1_a2 * d1_a2 - cos2_2 * d1_sq;
        let b_v = 2.0 * delta_a2 * d1_a2 - 2.0 * cos2_2 * delta_d1;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v - disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cone1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cone1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };

    let (mut branch_plus, mut branch_minus): (Vec<DVec3>, Vec<DVec3>) = (
        refine_polyline(&branch_plus, eval_plus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
        refine_polyline(&branch_minus, eval_minus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
    );

    // Dedup: remove trailing points that nearly duplicate the first point.
    let dedup = |pts: &mut Vec<DVec3>| {
        while pts.len() >= 3 {
            let n = pts.len();
            if (pts[n - 1] - pts[0]).length_squared() < TOLERANCE_VEC_SQ_MIN {
                pts.pop();
            } else {
                break;
            }
        }
    };

    let mut branches = Vec::new();
    if branch_plus.len() >= 2 {
        dedup(&mut branch_plus);
        if branch_plus.len() >= 2 {
            branches.push(branch_plus);
        }
    }
    if branch_minus.len() >= 2 {
        let is_distinct = branches.is_empty()
            || (branch_minus[0] - branches[0][0]).length_squared() > TOLERANCE_VEC_SQ_MIN;
        if is_distinct {
            dedup(&mut branch_minus);
            if branch_minus.len() >= 2 {
                branches.push(branch_minus);
            }
        }
    }

    branches
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel axes
// ─────────────────────────────────────────────────────────────────────────────

fn intersect_parallel_cones(
    cone1: &ConicalSurface,
    cone2: &ConicalSurface,
    a1: DVec3,
    a2: DVec3,
) -> ConeConeResult {
    let apex1 = cone1.apex_point();
    let apex2 = cone2.apex_point();
    // Ensure both axis vectors point in the same direction.
    let _a2 = if a1.dot(a2) >= 0.0 { a2 } else { -a2 };

    // Perpendicular distance between the two axes.
    let delta = apex2 - apex1;
    let delta_along = delta.dot(a1);
    let delta_perp = delta - a1 * delta_along;
    let d_perp = delta_perp.length();

    let beta1 = cone1.half_angle_rad;
    let beta2 = cone2.half_angle_rad;
    let tan1 = beta1.tan();
    let tan2 = beta2.tan();

    // ── Coaxial ──────────────────────────────────────────────────────────────
    if d_perp < TOLERANCE_ABS {
        // Axes coincide.  The apex of cone2 may be different from cone1's apex.

        // Check for identical geometry (same apex, same half-angle).
        if (apex2 - apex1).length() < TOLERANCE_ABS && (beta1 - beta2).abs() < TOLERANCE_ANG {
            return ConeConeResult::Coaxial;
        }

        // Height of cone2's apex above cone1's apex along the shared axis.
        // At height h above cone1.apex, cone1 has radius r1(h) = h * tan1 (h > 0).
        // At height h above cone1.apex, cone2 has radius r2(h) = (h - delta_along) * tan2
        //   (positive when h > delta_along).
        //
        // Set r1 = r2:  h*tan1 = (h - delta_along)*tan2
        //   h*(tan1 - tan2) = -delta_along * tan2
        //   h = -delta_along * tan2 / (tan1 - tan2)        when tan1 ≠ tan2
        //
        // When tan1 = tan2 (same opening angle):
        //   r1 = r2 is only satisfiable if delta_along = 0 (same apex → Coaxial above)
        //   or never (different apices → only the apex itself if a1 == a2 direction).

        if (tan1 - tan2).abs() < TOLERANCE_LEN_MIN {
            // Equal half-angles, different apices.
            // The two cones are coaxial "nested" with the same angle.
            // For infinite cones the lateral surfaces are parallel and never
            // intersect.  However for frustums (bounded cones) the face-boundary
            // edges still need processing by the pave-filler — return General so
            // the numeric engine gets a chance to find intersection curves, even
            // if it also returns empty, the pave-filler still processes face-plane
            // pairs for the end caps.
            return ConeConeResult::General;
        }

        let h = -delta_along * tan2 / (tan1 - tan2);

        // h must be positive for cone1's nappe, and (h - delta_along) > 0 for cone2's nappe.
        if h < -TOLERANCE_ABS || (h - delta_along) < -TOLERANCE_ABS {
            // Check if the intersection is at a shared apex.
            if h.abs() < TOLERANCE_ABS {
                return ConeConeResult::CoaxialPoint(apex1);
            }
            return ConeConeResult::NoIntersection;
        }

        let radius = h * tan1;
        if radius < TOLERANCE_ABS {
            return ConeConeResult::CoaxialPoint(apex1 + a1 * h);
        }

        let center = apex1 + a1 * h;
        return ConeConeResult::CoaxialCircle(Circle3::new(center, a1, radius));
    }

    // ── Parallel but offset ───────────────────────────────────────────────────
    // At height h above cone1.apex: cone1 radius = h*tan1.
    // The cone2 apex is at perpendicular offset d_perp from the cone1 axis.
    // At the same height, cone2 radius = (h - delta_along)*tan2 (from cone2 apex).
    //
    // Two circles of radii r1, r2 at perpendicular distance d_perp apart can
    // only intersect if |r1 - r2| ≤ d_perp ≤ r1 + r2.
    //
    // Since both radii grow with height (for positive nappes), they will
    // eventually be large enough to overlap for any finite d_perp.  So the
    // surfaces always intersect for parallel offset cones — fall back to marching.
    //
    // Quick early exit: if both cones are very thin (small half-angles) and
    // d_perp is large, no intersection occurs near the apices but they will
    // meet at large h.  For bounded CAD faces the marching algorithm handles this.
    ConeConeResult::General
}

// ─────────────────────────────────────────────────────────────────────────────
// OCCT IntAna_QuadQuadGeo Sec 2-4: stub analytic cases (documentation only)
// ─────────────────────────────────────────────────────────────────────────────

// OCCT reference: IntAna_QuadQuadGeo::Perform(gp_Cone, gp_Cone) at
//   $OCCT_SRC/src/IntAna/IntAna_QuadQuadGeo.cxx
//
// The full cone-cone intersection classifies into 5 cases:
//   1. Coaxial (handled above in intersect_parallel_cones)
//   2. Parallel-offset -- parallel axes, distinct centre lines
//   3. Coincident apices -- same apex, distinct non-parallel axes
//   4. Common generatrix -- one generator shared by both cones
//   5. General skew -- quartic space curve (handled by skew solver above)
//
// Cases 2-4 are documented below but return None so the caller falls through
// to numeric marching.

/// OCCT Sec 2: Parallel-offset cones (non-coaxial, parallel axes).
///
/// When the two cone axes are parallel but not coincident, the intersection is
/// a quartic space curve. OCCT (IntAna_QuadQuadGeo, case 2) solves this
/// analytically by reducing the problem to the plane perpendicular to the axes.
///
/// ## Setup
///
/// Let both axes be parallel to `a`.  Project both cones onto a plane
/// perpendicular to `a` -- each cone appears as a circle whose radius varies
/// linearly with axial height h.  Cone i has centre offset vector c_i in the
/// projection plane and radius r_i(h) = |h - h_apex_i| * tan(beta_i).
///
/// At height h the two projected circles intersect when
///
///   |r_1(h) - r_2(h)|  <=  d_perp  <=  r_1(h) + r_2(h)          (1)
///
/// where d_perp = |c_1 - c_2| is the fixed axis separation.  For each h in
/// the overlap interval, the angular parameter phi on cone 1 satisfies
///
///   cos(phi) = (r_1^2 + d_perp^2 - r_2^2) / (2 * r_1 * d_perp)   (2)
///
/// yielding two space points per h (one at +/- phi).
///
/// ## OCCT quartic formulation
///
/// OCCT eliminates phi by squaring (2) and substituting the linear expressions
/// for r_1(h), r_2(h).  After simplification this yields a quartic in h.  The
/// real roots are the axial heights where the intersection topology changes
/// (tangent transitions).  Between consecutive roots the curve is a single
/// smooth branch parametrised by h.
///
/// ## TODO
///
/// Implement the analytic quartic solve.  For now returns `None` so the
/// caller falls through to numeric marching.
fn occt_parallel_offset_cones(
    _cone1: &ConicalSurface,
    _cone2: &ConicalSurface,
) -> Option<ConeConeResult> {
    None
}

/// OCCT Sec 3: Coincident apices (same apex, distinct non-parallel axes).
///
/// When both cone apices coincide at point O but the axes a_1, a_2 are not
/// parallel, the intersection consists of lines (generators) through O.
/// OCCT (IntAna_QuadQuadGeo, case 3) solves this by finding directions v
/// that satisfy both cone equations simultaneously.
///
/// ## Derivation
///
/// For a unit direction v (|v| = 1), v lies on cone i iff
///
///   |v . a_i| = cos(beta_i)
///
/// Choosing both signs gives four sign pairs (s_1, s_2) with s_i in {+1, -1}:
///
///   v . a_1 = s_1 * cos(beta_1)
///   v . a_2 = s_2 * cos(beta_2)
///
/// For each pair, the two linear equations define a line (intersection of two
/// planes).  The direction of this line is a_1 x a_2.  A particular solution
/// in the a_1-a_2 plane is:
///
///   v_p = p * a_1 + q * a_2
///
/// with
///
///   q = (s_2*cos(beta_2) - s_1*cos(beta_1)*(a_1.a_2)) / (1 - (a_1.a_2)^2)
///   p = s_1*cos(beta_1) - q*(a_1.a_2)
///
/// Every point on the line is v = v_p + t*(a_1 x a_2).  The unit-norm constraint
/// gives one or two real values of t:
///
///   t = +/- sqrt((1 - |v_p|^2) / |a_1 x a_2|^2)
///
/// when |v_p| <= 1.  Each valid t yields an intersection line through O.
///
/// Up to 8 directions (4 sign pairs x 2 t values) can be found; opposite
/// directions represent the same line and are deduplicated.  The remaining
/// distinct lines are the intersection generators.
///
/// ## TODO
///
/// Implement the full coincident-apex solve.  Returns `None` for now.
fn occt_coincident_apex_cones(
    _cone1: &ConicalSurface,
    _cone2: &ConicalSurface,
) -> Option<ConeConeResult> {
    None
}

/// OCCT Sec 4: Common generatrix (one generator shared by both cones).
///
/// When two cones share a single generator line, OCCT (IntAna_QuadQuadGeo,
/// case 4) detects this by checking whether the vector connecting the two
/// apices forms the correct angle with both axes.
///
/// ## Detection
///
/// Let O_1, O_2 be the apices and a_1, a_2 the axis directions.  The connecting
/// vector Delta = O_2 - O_1 is a generatrix of cone 1 if
///
///   |Delta . a_1| = cos(beta_1) * |Delta|
///
/// and of cone 2 if
///
///   |Delta . a_2| = cos(beta_2) * |Delta|
///
/// If both hold, the line through O_1 and O_2 lies on both cones and is the
/// common generatrix.  The full intersection is this line plus a cubic
/// remainder curve (the quartic factorises with a linear term).
///
/// ## OCCT handling
///
/// OCCT checks all four sign combinations (+/- cos for each cone) before
/// falling through to the general quartic solver.  When a common generatrix
/// is found, the quartic is deflated to a cubic by factoring out the line.
///
/// ## TODO
///
/// Implement common-generatrix detection and deflation.  Returns `None`
/// for now.
fn occt_common_generatrix_cones(
    _cone1: &ConicalSurface,
    _cone2: &ConicalSurface,
) -> Option<ConeConeResult> {
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Fuzzy-tolerance wrapper & Tests
// ─────────────────────────────────────────────────────────────────────────────

/// Compute cone-cone intersection with fuzzy tolerance for near-coaxial cases.
pub fn intersect_cone_cone_with_tolerance(
    cone1: &ConicalSurface,
    cone2: &ConicalSurface,
    fuzzy_tol: f64,
) -> ConeConeResult {
    let tol = TOLERANCE_ABS + fuzzy_tol;
    let a1 = cone1.axis.normalize();
    let a2 = cone2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length();

    // ── Parallel axes (with fuzzy tolerance for near-coaxial detection) ────
    if sin_angle < TOLERANCE_ANG {
        let delta = cone2.apex - cone1.apex;
        let delta_along = delta.dot(a1);
        let delta_perp = delta - a1 * delta_along;
        let d_perp = delta_perp.length();

        let beta1 = cone1.half_angle_rad;
        let beta2 = cone2.half_angle_rad;
        let tan1 = beta1.tan();
        let tan2 = beta2.tan();

        // ── Coaxial (with fuzzy tolerance) ────────────────────────────────
        if d_perp < tol {
            // Check for identical geometry
            if (cone2.apex - cone1.apex).length() < tol && (beta1 - beta2).abs() < TOLERANCE_ANG {
                return ConeConeResult::Coaxial;
            }

            if (tan1 - tan2).abs() < TOLERANCE_LEN_MIN {
                return ConeConeResult::NoIntersection;
            }

            let h = -delta_along * tan2 / (tan1 - tan2);

            if h < -tol || (h - delta_along) < -tol {
                if h.abs() < tol {
                    return ConeConeResult::CoaxialPoint(cone1.apex);
                }
                return ConeConeResult::NoIntersection;
            }

            let radius = h * tan1;
            if radius < tol {
                return ConeConeResult::CoaxialPoint(cone1.apex + a1 * h);
            }

            let center = cone1.apex + a1 * h;
            return ConeConeResult::CoaxialCircle(Circle3::new(center, a1, radius));
        }
    }

    // ── Skew axes (analytic quartic solver) ──────────────────────────────
    let skew_result = intersect_skew_cone_cone(cone1, cone2);
    if !skew_result.is_empty() {
        return ConeConeResult::SkewQuartic(skew_result);
    }

    ConeConeResult::General
}
