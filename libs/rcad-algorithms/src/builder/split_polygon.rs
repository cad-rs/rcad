use std::collections::HashSet;
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*;
use crate::bopds::ds::*;
use crate::tolerance::*;
use crate::builder::types::{BooleanOpType, FaceSampleData, WireFace};
use crate::builder::split_polygon2::point_in_polygon_2d;
use crate::builder::wire_splitter::world_to_uv;

pub(crate) fn curved_subface_boundary_3d(
    uv_poly: &[DVec2],
    trim_polylines_uv: &[Vec<DVec2>],
    surface: &Surface3,
) -> Vec<DVec3> {
    // EDGE_SAMPLES must divide evenly into the 3D curve pre-sampling
    // density (128) used in split_planar_face so sphere and plane boundary
    // vertices share the same 3D positions along intersection curves.
    // Use fewer samples for high-vertex-count UV polygons (e.g. trims from
    // sphere_closed_trim_to_open_isolines with 65 vertices per meridian)
    // since each edge is already short.
    let edge_samples: usize = if matches!(&surface, rcad_kernel::geom::Surface3::Cylinder(_)) {
        // Cylinder: UV edges are straight lines, 2 samples per edge is sufficient.
        // OCCT BOPAlgo_BuilderFace uses exact edges, not sampled polylines.
        2
    } else if uv_poly.len() > 80 {
        4
    } else if uv_poly.len() > 30 {
        8
    } else if uv_poly.len() > 15 {
        16
    } else {
        32
    };

    let mut pts: Vec<DVec3> = Vec::new();

    // 1. Sample each UV edge and evaluate 3D positions
    let n = uv_poly.len();
    // Compute the u-span to detect winding polygons. When the UV polygon
    // spans > 蟺 in u, edges near the seam wrap the "long way" around the
    // sphere in 3D. We redirect such edges to go through the seam instead,
    // producing a compact 3D boundary.
    let pu_min = uv_poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let pu_max = uv_poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let is_winding = pu_max - pu_min > std::f64::consts::PI;
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];
        let du = b.x - a.x;
        if is_winding && du.abs() > std::f64::consts::PI {
            // Edge crosses the seam in a winding polygon.  Sample through
            // the seam (wrapping around) instead of the direct line, so the
            // 3D boundary goes the SHORT way around the sphere.
            let delta = if du > 0.0 { du - std::f64::consts::TAU } else { du + std::f64::consts::TAU };
            for k in 0..edge_samples {
                let t = k as f64 / edge_samples as f64;
                let u = a.x + t * delta;
                let v = a.y + t * (b.y - a.y);
                pts.push(surface.point_at(u, v));
            }
        } else {
            for k in 0..edge_samples {
                let t = k as f64 / edge_samples as f64;
                let uv = DVec2::new(a.x + t * du, a.y + t * (b.y - a.y));
                pts.push(surface.point_at(uv.x, uv.y));
            }
        }
    }

    // CHECKPOINT 5: after 3D point sampling
    if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
        let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        eprintln!("[SPHERE_SPLIT] checkpoint=5 uv_poly_nverts={} sampled_pts={} is_winding={} u_range=[{:.4},{:.4}] v_range=[{:.4},{:.4}]",
            uv_poly.len(), pts.len(), pu_max - pu_min > std::f64::consts::PI, pu_min, pu_max, v_min, v_max);
    }

    // 2. Consecutive deduplication 閳?collapse runs of pole/apex samples
    let mut deduped: Vec<DVec3> = Vec::new();
    for p in &pts {
        if deduped.is_empty() || (*p - deduped[deduped.len() - 1]).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS {
            deduped.push(*p);
        }
    }
    // Close the loop: remove last point if it equals the first
    if deduped.len() > 1 && (deduped[0] - deduped[deduped.len() - 1]).length_squared() < TOLERANCE_ABS * TOLERANCE_ABS {
        deduped.pop();
    }

    // 3. If still degenerate, supplement with trim polyline 3D points
    if deduped.len() < 3 {
        for trim_uv in trim_polylines_uv {
            if trim_uv.len() < 2 {
                continue;
            }
            for uv in trim_uv {
                let p3 = surface.point_at(uv.x, uv.y);
                if point_in_polygon_2d(uv_poly, *uv) || point_near_polygon_2d(uv_poly, *uv, 0.1) {
                    // Only add if not already in deduped
                    if deduped.iter().all(|q| (p3 - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
                        deduped.push(p3);
                    }
                }
            }
        }
    }

    // 4. Final global dedup
    let mut result: Vec<DVec3> = Vec::new();
    for p in &deduped {
        if result.iter().all(|q| (*p - *q).length_squared() > TOLERANCE_ABS * TOLERANCE_ABS) {
            result.push(*p);
        }
    }

    result
}

/// Check if a 2D point is within `margin` of any edge of a polygon.
pub(crate) fn point_near_polygon_2d(poly: &[DVec2], pt: DVec2, margin: f64) -> bool {
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = poly[i];
        let b = poly[j];
        let ab = b - a;
        let len_sq = ab.length_squared();
        let t = if len_sq < TOLERANCE_FLOAT_LOOSE { 0.0 } else { ((pt - a).dot(ab) / len_sq).clamp(0.0, 1.0) };
        let closest = a + t * ab;
        if (pt - closest).length() < margin {
            return true;
        }
    }
    false
}

/// Detect and handle UV seam crossings for periodic surfaces.
/// Returns a list of split polygons if the UV polygon crosses the seam.
pub(crate) fn handle_periodic_seam_crossing(
    uv_poly: &[DVec2],
    u_period: f64,
    seam_u: f64,
) -> Vec<Vec<DVec2>> {
    let n = uv_poly.len();
    if n < 3 || u_period <= 0.0 {
        return vec![uv_poly.to_vec()];
    }

    // Find all edges that cross the seam
    let mut seam_crossings: Vec<(usize, f64, DVec2)> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let u_i = uv_poly[i].x;
        let u_j = uv_poly[j].x;

        // Check for seam crossing (jump > period/2)
        let du = u_j - u_i;
        if du.abs() > u_period * 0.5 {
            // Compute intersection point with seam
            let t = if du > 0.0 {
                (seam_u + u_period - u_i) / du
            } else {
                (seam_u - u_i) / du
            };

            if t > 0.0 && t < 1.0 {
                let v_i = uv_poly[i].y;
                let v_j = uv_poly[j].y;
                let seam_v = v_i + t * (v_j - v_i);
                let seam_pt = DVec2::new(seam_u, seam_v);
                seam_crossings.push((i, t, seam_pt));
            }
        }
    }

    // If no crossings or odd number of crossings (invalid), return original
    if seam_crossings.is_empty() || !seam_crossings.len().is_multiple_of(2) {
        return vec![uv_poly.to_vec()];
    }

    // Sort crossings by edge index
    seam_crossings.sort_by_key(|&(idx, _, _)| idx);

    // For now, handle the simple case of exactly 2 crossings
    if seam_crossings.len() == 2 {
        let (idx1, _, pt1) = seam_crossings[0];
        let (idx2, _, pt2) = seam_crossings[1];

        // Build two sub-polygons
        let mut poly1: Vec<DVec2> = Vec::new();
        let mut poly2: Vec<DVec2> = Vec::new();

        // poly1: from crossing1 to crossing2 (wrapping the other way)
        poly1.push(pt1);
        for i in (idx1 + 1)..=idx2 {
            if i < n {
                poly1.push(uv_poly[i]);
            }
        }
        poly1.push(pt2);

        // poly2: from crossing2 back to crossing1
        poly2.push(pt2);
        for i in (idx2 + 1)..n {
            poly2.push(uv_poly[i]);
        }
        for i in 0..=idx1 {
            poly2.push(uv_poly[i]);
        }
        poly2.push(pt1);

        let mut result = Vec::new();
        if poly1.len() >= 3 {
            result.push(poly1);
        }
        if poly2.len() >= 3 {
            result.push(poly2);
        }

        if result.is_empty() {
            vec![uv_poly.to_vec()]
        } else {
            result
        }
    } else {
        // Multiple crossing pairs - complex case, return original for now
        vec![uv_poly.to_vec()]
    }
}

/// Split a polygon along a vertical u-isoline.
///
/// Used for sphere UV polygons whose u-span exceeds pi after normalisation.
/// Finds where the polygon crosses u=u_split and splits it into left and right
/// pieces, each bounded by the original polygon boundary on one side and the
/// isoline on the other.
pub(crate) fn split_polygon_at_u_isoline(poly: &[DVec2], u_split: f64) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }

    // Find all edges crossing u=u_split
    let mut crossings: Vec<(usize, DVec2)> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let u0 = poly[i].x;
        let u1 = poly[j].x;
        // Check if this edge crosses u_split
        if (u0 - u_split).abs() < TOLERANCE_COORD_SUB {
            // Vertex is on the isoline 鈥?use it directly
            if crossings.is_empty() || crossings.last().unwrap().0 != i {
                crossings.push((i, poly[i]));
            }
        } else if (u0 < u_split && u1 > u_split) || (u0 > u_split && u1 < u_split) {
            let t = (u_split - u0) / (u1 - u0);
            let v = poly[i].y + t * (poly[j].y - poly[i].y);
            crossings.push((i, DVec2::new(u_split, v)));
        }
    }

    if crossings.len() != 2 {
        return vec![poly.to_vec()];
    }

    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];

    // Build left polygon: edges from idx1+1 to idx2, plus pt1 and pt2
    let mut left: Vec<DVec2> = vec![pt1];
    for i in (idx1 + 1)..=idx2 {
        if i < n {
            left.push(poly[i]);
        }
    }
    left.push(pt2);

    // Build right polygon: edges from idx2+1 to n, then 0 to idx1, plus pt1 and pt2
    let mut right: Vec<DVec2> = vec![pt2];
    for i in (idx2 + 1)..n {
        right.push(poly[i]);
    }
    for i in 0..=idx1 {
        right.push(poly[i]);
    }
    right.push(pt1);

    let mut result = Vec::new();
    if left.len() >= 3 {
        result.push(left);
    }
    if right.len() >= 3 {
        result.push(right);
    }
    if result.is_empty() {
        vec![poly.to_vec()]
    } else {
        result
    }
}

struct BBox2 { u_min: f64, u_max: f64, v_min: f64, v_max: f64 }

pub(crate) fn bbox_of_poly(poly: &[DVec2]) -> BBox2 {
    let u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let v_min = poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    BBox2 { u_min, u_max, v_min, v_max }
}

/// Compute a tighter UV bounding box by sampling interior points of the polygon.
///
/// The polygon's boundary vertices may include trim curves that inflate the
/// bounding box far beyond the actual interior region (e.g. a trim curve that
/// wanders from u=-pi to u=pi but bounds a region that only occupies u=[0,pi]).
/// Sampling interior points and taking their min/max gives the true extent.
pub(crate) fn compute_interior_uv_bounds(
    poly: &[DVec2],
    bnd_u_min: f64,
    bnd_u_max: f64,
    bnd_v_min: f64,
    bnd_v_max: f64,
) -> (f64, f64, f64, f64) {
    const N_U: usize = 11;
    const N_V: usize = 11;
    let du = (bnd_u_max - bnd_u_min) / (N_U as f64 + 1.0);
    let dv = (bnd_v_max - bnd_v_min) / (N_V as f64 + 1.0);
    if du <= 0.0 || dv <= 0.0 {
        return (bnd_u_min, bnd_u_max, bnd_v_min, bnd_v_max);
    }

    let mut in_u_min = f64::INFINITY;
    let mut in_u_max = f64::NEG_INFINITY;
    let mut in_v_min = f64::INFINITY;
    let mut in_v_max = f64::NEG_INFINITY;
    let mut found = false;

    for iu in 1..=N_U {
        let u = bnd_u_min + du * iu as f64;
        for iv in 1..=N_V {
            let v = bnd_v_min + dv * iv as f64;
            if point_in_polygon_2d(poly, DVec2::new(u, v)) {
                in_u_min = in_u_min.min(u);
                in_u_max = in_u_max.max(u);
                in_v_min = in_v_min.min(v);
                in_v_max = in_v_max.max(v);
                found = true;
            }
        }
    }

    if found {
        // Expand slightly to account for sampling grid granularity
        let pad_u = du * 0.6;
        let pad_v = dv * 0.6;
        (
            (in_u_min - pad_u).max(bnd_u_min),
            (in_u_max + pad_u).min(bnd_u_max),
            (in_v_min - pad_v).max(bnd_v_min),
            (in_v_max + pad_v).min(bnd_v_max),
        )
    } else {
        (bnd_u_min, bnd_u_max, bnd_v_min, bnd_v_max)
    }
}

/// Detect degenerate points (poles, apex) and handle them in UV polygon.
/// Returns a modified 3D boundary that correctly handles surface singularities.
pub(crate) fn handle_degenerate_points(
    uv_poly: &[DVec2],
    surface: &Surface3,
) -> Vec<DVec3> {
    match surface {
        Surface3::Sphere(s) => {
            // Sphere has two poles at v=0 and v=锜?
            let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
            let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

            let mut boundary_3d = Vec::new();
            let pole_tol = 0.01; // Tolerance for detecting near-pole

            // Check if polygon touches the north pole (v 閳?0)
            let touches_north_pole = v_min < pole_tol;
            // Check if polygon touches the south pole (v 閳?锜?
            let touches_south_pole = v_max > std::f64::consts::PI - pole_tol;

            if touches_north_pole || touches_south_pole {
                // Sample the UV polygon edges more densely near poles

                // Detect winding polygon (UV spans > pi in u) so edges that
                // cross the seam go the SHORT way around the sphere in 3D.
                let pu_min = uv_poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                let pu_max = uv_poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                let is_winding = pu_max - pu_min > std::f64::consts::PI;

                // Sample UV edges
                let n = uv_poly.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    let a = uv_poly[i];
                    let b = uv_poly[j];

                    let du = b.x - a.x;

                    // More samples if edge is near pole
                    let near_pole = (a.y < pole_tol || a.y > std::f64::consts::PI - pole_tol)
                        || (b.y < pole_tol || b.y > std::f64::consts::PI - pole_tol);
                    let n_samples = if near_pole { 16 } else { 4 };

                    if is_winding && du.abs() > std::f64::consts::PI {
                        // Edge crosses the seam in a winding polygon. Sample
                        // through the seam (wrapping around) instead of the
                        // direct line, so the 3D boundary goes the SHORT way.
                        let delta = if du > 0.0 { du - std::f64::consts::TAU } else { du + std::f64::consts::TAU };
                        for k in 0..n_samples {
                            let t = k as f64 / n_samples as f64;
                            let u = a.x + t * delta;
                            let v = a.y + t * (b.y - a.y);
                            // CHECKPOINT 4A: before v-clamp in winding path
                            if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
                                eprintln!("[SPHERE_SPLIT] checkpoint=4a v_before_clamp={:.6}", v);
                            }
                            let v_clamped = v.clamp(0.001, std::f64::consts::PI - 0.001);
                            let pt = s.point_at(u, v_clamped);
                            boundary_3d.push(pt);
                        }
                    } else {
                        for k in 0..n_samples {
                            let t = k as f64 / n_samples as f64;
                            let uv = DVec2::new(
                                a.x + t * du,
                                a.y + t * (b.y - a.y),
                            );

                            // Clamp v to avoid pole singularity
                            // CHECKPOINT 4B: before v-clamp in non-winding path
                            if std::env::var("RCAD_DEBUG_SPHERE_SPLIT").is_ok() {
                                eprintln!("[SPHERE_SPLIT] checkpoint=4b v_before_clamp={:.6}", uv.y);
                            }
                            let v_clamped = uv.y.clamp(0.001, std::f64::consts::PI - 0.001);
                            let pt = s.point_at(uv.x, v_clamped);

                            boundary_3d.push(pt);
                        }
                    }
                }

                // NOTE: we do NOT add a separate pole-point vertex here.  The
                // clamped-v edge samples (v clamped to 0.001 / PI-0.001) already
                // span the full u-range of the face.  Adding a pole point at
                // (0.0, v=0|PI) creates a diagonal closing edge for faces whose
                // u-range doesn't include 0, deforming the UV polygon.
            } else {
                // No pole involvement - standard sampling
                for &uv in uv_poly {
                    boundary_3d.push(surface.point_at(uv.x, uv.y));
                }
            }

            // Deduplicate
            dedup_3d_points(&boundary_3d)
        }
        Surface3::Cone(c) => {
            // Cone has apex at v=0
            let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

            if v_min < 0.01 {
                // Near apex - need special handling
                let apex = c.apex_point();
                let mut boundary_3d = Vec::new();

                let n = uv_poly.len();
                for i in 0..n {
                    let j = (i + 1) % n;
                    let a = uv_poly[i];
                    let b = uv_poly[j];

                    // Check if edge crosses near apex
                    let near_apex = a.y < 0.1 || b.y < 0.1;
                    let n_samples = if near_apex { 16 } else { 4 };

                    for k in 0..n_samples {
                        let t = k as f64 / n_samples as f64;
                        let uv = DVec2::new(
                            a.x + t * (b.x - a.x),
                            a.y + t * (b.y - a.y),
                        );

                        // Clamp v to avoid apex singularity
                        let v_clamped = uv.y.max(0.001);
                        let pt = c.point_at(uv.x, v_clamped);

                        // Skip points very close to apex
                        if (pt - apex).length() > 0.01 {
                            boundary_3d.push(pt);
                        }
                    }
                }

                // Add apex if polygon contains it
                boundary_3d.push(apex);

                dedup_3d_points(&boundary_3d)
            } else {
                // Standard case
                uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
            }
        }
        _ => {
            // No degenerate points - standard mapping
            uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
        }
    }
}

/// Enhanced handling of degenerate UV polygons on surfaces with singularities.
///
/// This function handles UV polygons where vertices collapse at surface singularities:
/// - Sphere poles (v=0 or v=锜?
/// - Cone apex (v=0)
///
/// The function:
/// 1. Detects pole/apex proximity
/// 2. Handles triangulation specially for degenerate triangles
/// 3. Ensures edge PCurve tolerance near poles/apex
///
/// Returns a 3D boundary that correctly handles surface singularities.
pub fn handle_degenerate_uv_polygon(uv_poly: &[DVec2], surface: &Surface3) -> Vec<DVec3> {
    match surface {
        Surface3::Sphere(s) => {
            handle_sphere_degenerate_uv(uv_poly, s)
        }
        Surface3::Cone(c) => {
            handle_cone_degenerate_uv(uv_poly, c)
        }
        _ => {
            // No degenerate points - standard mapping
            uv_poly.iter().map(|uv| surface.point_at(uv.x, uv.y)).collect()
        }
    }
}

/// Handle degenerate UV polygons on sphere surfaces.
pub(crate) fn handle_sphere_degenerate_uv(uv_poly: &[DVec2], sphere: &SphericalSurface) -> Vec<DVec3> {
    let pole_tol = 0.01; // Tolerance for detecting near-pole

    // Find min/max v values to detect pole proximity
    let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max = uv_poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

    // Check if polygon touches either pole
    let touches_north_pole = v_min < pole_tol;
    let touches_south_pole = v_max > std::f64::consts::PI - pole_tol;

    if !touches_north_pole && !touches_south_pole {
        // No pole involvement - standard mapping
        return uv_poly.iter().map(|uv| sphere.point_at(uv.x, uv.y)).collect();
    }

    let mut boundary_3d = Vec::new();

    // Determine which pole(s) are involved
    let north_pole = sphere.center + sphere.axis * sphere.radius;
    let south_pole = sphere.center - sphere.axis * sphere.radius;

    // Sample UV polygon edges more densely near poles
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];

        // More samples if edge is near pole
        let near_pole = (a.y < pole_tol || a.y > std::f64::consts::PI - pole_tol)
            || (b.y < pole_tol || b.y > std::f64::consts::PI - pole_tol);
        let n_samples = if near_pole { 16 } else { 4 };

        for k in 0..n_samples {
            let t = k as f64 / n_samples as f64;
            let uv = DVec2::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
            );

            // Clamp v to avoid pole singularity
            let v_clamped = uv.y.clamp(0.001, std::f64::consts::PI - 0.001);
            let pt = sphere.point_at(uv.x, v_clamped);

            // Skip points very close to pole (will add pole point separately)
            let near_north = (pt - north_pole).length() < sphere.radius * 0.1;
            let near_south = (pt - south_pole).length() < sphere.radius * 0.1;
            if !near_north && !near_south {
                boundary_3d.push(pt);
            }
        }
    }

    // Add pole point(s) if polygon contains them
    if touches_north_pole {
        boundary_3d.push(north_pole);
    }
    if touches_south_pole {
        boundary_3d.push(south_pole);
    }

    dedup_3d_points(&boundary_3d)
}

/// Handle degenerate UV polygons on cone surfaces.
pub(crate) fn handle_cone_degenerate_uv(uv_poly: &[DVec2], cone: &ConicalSurface) -> Vec<DVec3> {
    let apex_tol = 0.01; // Tolerance for detecting near-apex

    // Find min v value to detect apex proximity
    let v_min = uv_poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

    if v_min >= apex_tol {
        // No apex involvement - standard mapping
        return uv_poly.iter().map(|uv| cone.point_at(uv.x, uv.y)).collect();
    }

    let mut boundary_3d = Vec::new();
    let apex = cone.apex_point();

    // Sample UV polygon edges more densely near apex
    let n = uv_poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let a = uv_poly[i];
        let b = uv_poly[j];

        // More samples if edge is near apex
        let near_apex = a.y < apex_tol * 10.0 || b.y < apex_tol * 10.0;
        let n_samples = if near_apex { 16 } else { 4 };

        for k in 0..n_samples {
            let t = k as f64 / n_samples as f64;
            let uv = DVec2::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
            );

            // Clamp v to avoid apex singularity
            let v_clamped = uv.y.max(0.001);
            let pt = cone.point_at(uv.x, v_clamped);

            // Skip points very close to apex
            if (pt - apex).length() > 0.01 {
                boundary_3d.push(pt);
            }
        }
    }

    // Add apex point
    boundary_3d.push(apex);

    dedup_3d_points(&boundary_3d)
}

/// Split an edge at a periodic seam if it crosses the U=0/2锜?boundary.
///
/// This function detects if an edge on a periodic surface (cylinder, sphere, torus)
/// crosses the seam and splits it at the crossing point.
///
/// Returns:
/// - `None` if the edge doesn't cross the seam
/// - `Some(vec![seg1, seg2])` where each segment is `[start_uv, end_uv]`
pub fn split_edge_at_periodic_seam(
    start_uv: DVec2,
    end_uv: DVec2,
    surface: &Surface3,
) -> Option<Vec<Vec<DVec2>>> {
    // Get the U period for this surface type
    let u_period = match surface {
        Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            std::f64::consts::TAU
        }
        Surface3::Cone(_) => {
            // Cone is also periodic in U
            std::f64::consts::TAU
        }
        _ => {
            // Non-periodic surface
            return None;
        }
    };

    let u1 = start_uv.x;
    let u2 = end_uv.x;
    let v1 = start_uv.y;
    let v2 = end_uv.y;
    let du = u2 - u1;

    // Check for seam crossing (jump > period/2)
    if du.abs() <= u_period * 0.5 {
        return None;
    }

    // Determine which way we're crossing
    let is_low_to_high = du < 0.0; // u1 is high, u2 is low

    // Calculate intersection point at seam
    let (t, seam_u) = if is_low_to_high {
        // u1 is near period, u2 is near 0
        // Find t where u = period
        let t = (u_period - u1) / ((u2 + u_period) - u1);
        (t, u_period)
    } else {
        // u1 is near 0, u2 is near period
        // Find t where u = 0
        let t = -u1 / ((u2 - u_period) - u1);
        (t, 0.0)
    };

    // Clamp t to [0, 1] for numerical stability
    let t = t.clamp(0.0, 1.0);
    let seam_v = v1 + t * (v2 - v1);

    // Build two segments
    let seam_point = DVec2::new(seam_u, seam_v);
    let opposite_seam_point = if seam_u < u_period * 0.5 {
        DVec2::new(u_period, seam_v)
    } else {
        DVec2::new(0.0, seam_v)
    };

    // First segment: from start to seam
    let seg1 = vec![start_uv, seam_point];
    // Second segment: from opposite seam to end
    let seg2 = vec![opposite_seam_point, end_uv];

    Some(vec![seg1, seg2])
}

/// Split a UV polygon at both U and V seams for torus double periodicity.
///
/// The torus has two periodic parameters:
/// - U period: 2锜?(around major circle)
/// - V period: 2锜?(around tube circle)
///
/// This function handles UV polygon splitting in both directions.
pub fn split_uv_polygon_torus_double(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // First, split at U seam
    let u_split = split_uv_polygon_at_seam(uv_polygon, period);

    // Then, split each result at V seam
    let mut result = Vec::new();
    for poly in u_split {
        let v_split = split_uv_polygon_at_v_seam(&poly, period);
        result.extend(v_split);
    }

    result
}

/// Split a UV polygon at the V periodic seam (V=0/period boundary).
///
/// This is similar to split_uv_polygon_at_seam but for the V parameter.
pub(crate) fn split_uv_polygon_at_v_seam(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // Find all edges crossing the V seam
    let mut crossings: Vec<(usize, f64, DVec2)> = Vec::new();

    for i in 0..uv_polygon.len() {
        let j = (i + 1) % uv_polygon.len();
        let v1 = uv_polygon[i].y;
        let v2 = uv_polygon[j].y;
        let dv = v2 - v1;

        // Check for seam crossing (jump > period/2)
        if dv.abs() > period * 0.5 {
            let u1 = uv_polygon[i].x;
            let u2 = uv_polygon[j].x;

            // Determine which way we're crossing
            let is_low_to_high = dv < 0.0; // v1 is high, v2 is low

            // Calculate intersection point
            let (t, seam_v) = if is_low_to_high {
                let t = (period - v1) / ((v2 + period) - v1);
                (t, period)
            } else {
                let t = -v1 / ((v2 - period) - v1);
                (t, 0.0)
            };

            let t = t.clamp(0.0, 1.0);
            let seam_u = u1 + t * (u2 - u1);

            crossings.push((i, t, DVec2::new(seam_u, seam_v)));
        }
    }

    if crossings.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    // For now, handle simple cases
    if crossings.len() != 2 {
        // Complex case - return original
        return vec![uv_polygon.to_vec()];
    }

    // Build two sub-polygons
    let (_idx1, _, _pt1) = crossings[0];
    let (_idx2, _, _pt2) = crossings[1];

    let mut low_polygon: Vec<DVec2> = Vec::new();
    let mut high_polygon: Vec<DVec2> = Vec::new();

    let is_low = |v: f64| v < period * 0.5;

    let n = uv_polygon.len();

    // Traverse polygon and assign vertices
    for i in 0..n {
        let curr = uv_polygon[i];

        // Add current vertex to appropriate polygon
        if is_low(curr.y) {
            low_polygon.push(curr);
        } else {
            high_polygon.push(curr);
        }

        // Check for crossing between i and i+1
        for (cross_idx, _, cross_pt) in &crossings {
            if *cross_idx == i {
                // Add seam points to both polygons
                let low_pt = DVec2::new(cross_pt.x, 0.0);
                let high_pt = DVec2::new(cross_pt.x, period);

                if is_low(curr.y) {
                    low_polygon.push(low_pt);
                    high_polygon.push(high_pt);
                } else {
                    high_polygon.push(high_pt);
                    low_polygon.push(low_pt);
                }
            }
        }
    }

    let mut result = Vec::new();
    if low_polygon.len() >= 3 {
        result.push(low_polygon);
    }
    if high_polygon.len() >= 3 {
        result.push(high_polygon);
    }

    if result.is_empty() {
        vec![uv_polygon.to_vec()]
    } else {
        result
    }
}

/// Deduplicate 3D points within tolerance.
pub(crate) fn dedup_3d_points(points: &[DVec3]) -> Vec<DVec3> {
    let mut result: Vec<DVec3> = Vec::new();
    let tol_sq = TOLERANCE_ABS * TOLERANCE_ABS;

    for &p in points {
        if result.iter().all(|q: &DVec3| (p - *q).length_squared() > tol_sq) {
            result.push(p);
        }
    }

    result
}

/// Check if a UV trim is a closed loop (first and last points coincide).
pub(crate) fn is_closed_uv_trim(trim: &[DVec2]) -> bool {
    if trim.len() < 3 {
        return false;
    }
    let d_sq = (trim[0] - trim[trim.len() - 1]).length_squared();
    d_sq < TOLERANCE_LINEAR_ULTRA_STRICT
}

/// Check if a UV polygon is valid (has sufficient area and no degenerate edges).
pub(crate) fn is_valid_uv_polygon(poly: &[DVec2]) -> bool {
    if poly.len() < 3 {
        return false;
    }

    // Check for sufficient area (shoelace formula)
    let mut area = 0.0;
    let n = poly.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly[i].x * poly[j].y;
        area -= poly[j].x * poly[i].y;
    }
    area = area.abs() * 0.5;

    // Area should be significant
    area > TOLERANCE_LINEAR_ULTRA_STRICT
}

/// Convert a closed loop trim on a sphere face to open boundary-to-boundary
/// meridian isolines.  Sphere great-circle PCurves often produce closed UV
/// loops because the UV parameterization has a singularity at the poles
/// (atan2(0,0)=0).  The min and max u-values of such a closed loop directly
/// give the two meridian positions of the great circle.
///
/// Returns one or two open isolines, or `None` if the trim is not a convertible
/// great-circle loop.
pub(crate) fn sphere_closed_trim_to_open_isolines(
    trim: &[DVec2],
    uv_boundary: &[DVec2],
) -> Option<Vec<Vec<DVec2>>> {
    if trim.len() < 4 {
        return None;
    }
    let first = trim[0];
    let last = trim[trim.len() - 1];
    if (first - last).length_squared() >= TOLERANCE_LEN_MIN {
        return None; // not closed
    }

    let bnd_u_min = uv_boundary.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let bnd_u_max = uv_boundary.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let bnd_v_min = uv_boundary.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let bnd_v_max = uv_boundary.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let bnd_u_span = bnd_u_max - bnd_u_min;
    let bnd_v_span = bnd_v_max - bnd_v_min;

    let trim_u_min = trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let trim_u_max = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let trim_v_min = trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let trim_v_max = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

    // Must cover most of the UV rectangle (great circle).
    let u_coverage = (trim_u_max - trim_u_min) / bnd_u_span.abs();
    let v_coverage = (trim_v_max - trim_v_min) / bnd_v_span.abs();
    if u_coverage < 0.35 || v_coverage < 0.75 {
        // Latitude (constant-v) great circles: v 閳?constant but u spans full range.
        // These are great circles like the equator that DON'T pass through the poles,
        // so they form a horizontal line in UV space, not a closed pole-to-pole loop.
        if (trim_v_max - trim_v_min).abs() <= TOLERANCE_COORD_SUB && u_coverage >= 0.9 {
            let v_level = (trim_v_min + trim_v_max) / 2.0;
            if (v_level - bnd_v_min).abs() > TOLERANCE_COORD_SUB
                && (v_level - bnd_v_max).abs() > TOLERANCE_COORD_SUB
            {
                return Some(vec![
                    vec![DVec2::new(bnd_u_min, v_level), DVec2::new(bnd_u_max, v_level)]
                ]);
            }
        }
        return None;
    }

    // The min and max u-values give the two meridian positions.
    let mut u_vals = vec![trim_u_min, trim_u_max];
    // Deduplicate: if the two values are within 5% of the period, they're the same meridian
    let period = bnd_u_span;
    let diff = (u_vals[1] - u_vals[0]).abs();
    if diff > period * 0.5 {
        // The values straddle the seam (e.g. 锜?and -锜? 閳?wrap to get the effective difference
        let wrapped = (u_vals[1] + period - u_vals[0]).abs();
        if wrapped < period * 0.05 {
            // Same point 閳?only one meridian
            u_vals.pop();
        }
    } else if diff < period * 0.05 {
        u_vals.pop();
    }

    let mut isolines: Vec<Vec<DVec2>> = Vec::new();
    for &u in &u_vals {
        // Skip if the meridian is ON the boundary edge (within 1% of period)
        let dist_to_left = (u - bnd_u_min).abs();
        let dist_to_right = (u - bnd_u_max).abs();
        let edge_tol = period * 0.01;
        if dist_to_left < edge_tol || dist_to_right < edge_tol {
            continue;
        }
        // Sample 64 intermediate points along the meridian so the 3D
        // boundary accurately follows the sphere surface (instead of a
        // straight chord between the two endpoints).
        const MERIDIAN_N: usize = 64;
        let mut line: Vec<DVec2> = Vec::with_capacity(MERIDIAN_N + 1);
        let v_step = (bnd_v_max - bnd_v_min) / MERIDIAN_N as f64;
        for i in 0..=MERIDIAN_N {
            let v = bnd_v_min + v_step * i as f64;
            line.push(DVec2::new(u, v));
        }
        isolines.push(line);
    }

    if isolines.is_empty() { None } else { Some(isolines) }
}

pub(crate) fn periodic_trim_to_open_isoline(poly: &[DVec2], trim: &[DVec2], u_period: f64) -> Option<Vec<DVec2>> {
    if poly.len() < 3 || trim.len() < 3 || u_period <= 0.0 {
        return None;
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];
    let close_sq = uv_polyline_trim_closed_len_sq_from_uv_poly(poly);
    let is_closed = (trim_start - trim_end).length_squared() < close_sq;
    if !is_closed {
        return None;
    }

    let u_min_trim = trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let u_max_trim = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let v_min_trim = trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let v_max_trim = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let u_span = u_max_trim - u_min_trim;
    let v_span = v_max_trim - v_min_trim;

    let poly_u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let poly_u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_min = poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let poly_v_max = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let poly_v_span = poly_v_max - poly_v_min;

    if u_span < 0.9 * u_period {
        return None;
    }
    if poly_v_span <= TOLERANCE_LEN_MIN || v_span > 0.1 * poly_v_span {
        return None;
    }

    let v_level = trim.iter().map(|p| p.y).sum::<f64>() / trim.len() as f64;
    if v_level <= poly_v_min + TOLERANCE_COORD_SUB || v_level >= poly_v_max - TOLERANCE_COORD_SUB {
        return None;
    }

    Some(vec![
        DVec2::new(poly_u_min, v_level),
        DVec2::new(poly_u_max, v_level),
    ])
}

/// Split a UV polygon at periodic seams (U=0/period boundary).
///
/// For periodic surfaces like cylinders, the U parameter wraps around.
/// When a polygon crosses the seam (U=0 or U=period), we need to split it
/// into separate polygons, each with consistent U coordinates.
///
/// Algorithm:
/// 1. Find edges that cross the seam (|du| > period * 0.5)
/// 2. For each crossing edge, compute the exact intersection point at U=0 or U=period
/// 3. Build output polygons by inserting intersection points
///
/// Returns one or more polygons that don't cross the seam.
pub fn split_uv_polygon_at_seam(uv_polygon: &[DVec2], period: f64) -> Vec<Vec<DVec2>> {
    if uv_polygon.len() < 3 {
        return vec![uv_polygon.to_vec()];
    }

    // Structure to hold information about seam crossings
    struct SeamCrossing {
        edge_idx: usize,
        intersection: DVec2,
        is_low_to_high: bool, // true: crossing from low u (near 0) to high u (near period)
    }

    // Find all edges crossing the seam and compute intersection points
    let mut crossings: Vec<SeamCrossing> = Vec::new();
    for i in 0..uv_polygon.len() {
        let j = (i + 1) % uv_polygon.len();
        let u1 = uv_polygon[i].x;
        let u2 = uv_polygon[j].x;
        let v1 = uv_polygon[i].y;
        let v2 = uv_polygon[j].y;
        let du = u2 - u1;

        // Large jump indicates seam crossing
        if du.abs() > period * 0.5 {
            // Determine which way we're crossing
            // du > 0: wrapping from low u to high u (crossing U=0 going backwards in unwrapped space)
            // du < 0: wrapping from high u to low u (crossing U=period going backwards in unwrapped space)
            let is_low_to_high = du < 0.0; // u1 is high, u2 is low

            // Calculate intersection point using linear interpolation
            // We need to find the V coordinate where the edge crosses the seam
            //
            // For an edge from (u1, v1) to (u2, v2) crossing the seam:
            // If u1 is near period and u2 is near 0: unwrap u2 to u2 + period, find where U = period
            // If u1 is near 0 and u2 is near period: unwrap u2 to u2 - period, find where U = 0
            let (t, seam_u) = if is_low_to_high {
                // u1 is near period, u2 is near 0
                // Unwrap u2: consider edge from (u1, v1) to (u2 + period, v2)
                // Find t where u = period
                let t = (period - u1) / ((u2 + period) - u1);
                (t, period)
            } else {
                // u1 is near 0, u2 is near period
                // Unwrap u2: consider edge from (u1, v1) to (u2 - period, v2)
                // Find t where u = 0 (which equals period in the unwrapped space)
                // Or equivalently: the edge goes from u1 to u2-period (negative)
                // We want u = 0, so t = (0 - u1) / ((u2 - period) - u1) = -u1 / (u2 - period - u1)
                let t = -u1 / ((u2 - period) - u1);
                (t, 0.0)
            };

            // Clamp t to [0, 1] to handle numerical edge cases
            let t = t.clamp(0.0, 1.0);
            let intersection_v = v1 + t * (v2 - v1);

            crossings.push(SeamCrossing {
                edge_idx: i,
                intersection: DVec2::new(seam_u, intersection_v),
                is_low_to_high,
            });        }
    }

    if crossings.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    // Build output polygons
    // We need to partition the vertices and insert intersection points
    // Each output polygon will have consistent U values (all low or all high)

    // Collect all vertices and their positions relative to the seam
    // "low" means u < period * 0.5, "high" means u >= period * 0.5
    let is_low = |u: f64| u < period * 0.5;

    // Build two polygons: one for low-u region, one for high-u region
    let mut low_polygon: Vec<DVec2> = Vec::new();
    let mut high_polygon: Vec<DVec2> = Vec::new();

    // Sort crossings by edge index for efficient lookup
    let crossing_map: std::collections::HashMap<usize, &SeamCrossing> = crossings
        .iter()
        .map(|c| (c.edge_idx, c))
        .collect();

    // Traverse the polygon and assign vertices to appropriate output polygons
    for i in 0..uv_polygon.len() {
        let curr = uv_polygon[i];
        let next_idx = (i + 1) % uv_polygon.len();
        let _next = uv_polygon[next_idx];

        // Add current vertex to appropriate polygon
        if is_low(curr.x) {
            low_polygon.push(curr);
        } else {
            high_polygon.push(curr);
        }

        // Check if edge (i, i+1) crosses the seam
        if let Some(crossing) = crossing_map.get(&i) {
            // Add intersection point to both polygons
            // The intersection point is at the seam (u = 0 or u = period)
            // For the low polygon, we want u = 0
            // For the high polygon, we want u = period
            let low_intersection = DVec2::new(0.0, crossing.intersection.y);
            let high_intersection = DVec2::new(period, crossing.intersection.y);

            if crossing.is_low_to_high {
                // Going from high u to low u
                // Add period-point to high polygon first, then 0-point to low polygon
                high_polygon.push(high_intersection);
                low_polygon.push(low_intersection);
            } else {
                // Going from low u to high u
                // Add 0-point to low polygon first, then period-point to high polygon
                low_polygon.push(low_intersection);
                high_polygon.push(high_intersection);
            }
        }
    }

    // Build result - only include valid polygons (at least 3 vertices)
    let mut result = Vec::new();

    if low_polygon.len() >= 3 {
        result.push(low_polygon);
    }
    if high_polygon.len() >= 3 {
        result.push(high_polygon);
    }

    // If we didn't get valid polygons, return the original
    if result.is_empty() {
        return vec![uv_polygon.to_vec()];
    }

    result
}
