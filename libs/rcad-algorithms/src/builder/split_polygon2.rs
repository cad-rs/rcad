use std::collections::HashSet;
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*;
use crate::bopds::ds::*;
use crate::tolerance::*;
use crate::builder::types::{BooleanOpType, FaceSampleData, WireFace};
use crate::builder::wire_splitter::world_to_uv;


/// Split a 2D UV polygon by a 2D trim polyline.
///
/// Algorithm:
/// 1. Find trim start/end's closest edge on the polygon boundary.
/// 2. Project trim endpoints onto boundary edges to find exact split points.
/// 3. Split polygon into two halves at those points, inserting the trim polyline
///    between them.
///
/// For closed trim polylines (start 閳?end), uses a closed-curve splitting
/// algorithm: the trim forms an interior polygon that divides the outer polygon
/// into "inside trim" and "outside trim" regions.
///
/// Returns 1 polygon if splitting is degenerate, or 2 sub-polygons otherwise.
pub(crate) fn split_uv_polygon_by_trim(poly: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 || trim.len() < 2 {
        return vec![poly.to_vec()];
    }

    let trim_start = trim[0];
    let trim_end = trim[trim.len() - 1];

    // Find closest point on each polygon edge for a query point.
    // Returns (edge_index, t_param, projected_point).
    let closest_on_boundary = |q: DVec2| -> (usize, f64, DVec2) {
        let mut best_edge = 0usize;
        let mut best_t = 0.0f64;
        let mut best_pt = poly[0];
        let mut best_dist = f64::INFINITY;
        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let len_sq = ab.dot(ab);
            let t = if len_sq < TOLERANCE_FLOAT_LOOSE {
                0.0
            } else {
                ((q - a).dot(ab) / len_sq).clamp(0.0, 1.0)
            };
            let proj = a + t * ab;
            let dist = (q - proj).length_squared();
            if dist < best_dist {
                best_dist = dist;
                best_edge = i;
                best_t = t;
                best_pt = proj;
            }
        }
        (best_edge, best_t, best_pt)
    };

    // Cast a 2D ray from `origin` along `dir` and return the first boundary edge
    // intersection with t > -eps (including slightly behind for on-boundary starts).
    // Returns None if no intersection is found within a reasonable range.
    let ray_to_boundary = |origin: DVec2, dir: DVec2| -> Option<(usize, DVec2)> {
        let dir_len = dir.length();
        if dir_len < TOLERANCE_LEN_MIN {
            return None;
        }
        let dir = dir / dir_len;
        let mut best_t = f64::INFINITY;
        let mut best_edge = 0usize;
        let mut best_pt = poly[0];
        for i in 0..n {
            let j = (i + 1) % n;
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            // Solve: origin + t*dir = a + s*ab
            // => t*(dir鑴砤b) = (a-origin)鑴砤b  (2D cross: x.x*y.y - x.y*y.x)
            let denom = dir.x * ab.y - dir.y * ab.x;
            if denom.abs() < TOLERANCE_FLOAT_LOOSE {
                continue; // parallel
            }
            let oa = a - origin;
            let t_ray = (oa.x * ab.y - oa.y * ab.x) / denom;
            let s_seg = (oa.x * dir.y - oa.y * dir.x) / denom;
            if t_ray > -TOLERANCE_COORD_SUB && (-TOLERANCE_COORD_SUB..=1.0 + TOLERANCE_COORD_SUB).contains(&s_seg) && t_ray < best_t {
                best_t = t_ray;
                best_edge = i;
                best_pt = a + s_seg.clamp(0.0, 1.0) * ab;
            }
        }
        if best_t.is_finite() {
            Some((best_edge, best_pt))
        } else {
            None
        }
    };

    // Compute UV polygon bounding box to compute a "near-boundary" threshold
    let u_span = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
        - poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let v_span = poly.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max)
        - poly.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let boundary_snap_tol = (u_span + v_span) * 0.05;

    // 閳光偓閳光偓 Closed trim detection 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
    // Detect truly-closed trim: start 閳?end in UV space (e.g. a small loop entirely
    // inside the face).  Wrapped-closed trims (start and end differ by ~2锜?in u,
    // representing a full-circle cut around a cylinder or sphere) are intentionally
    // NOT treated as closed loops here 閳?they are open trims whose endpoints lie on
    // opposite sides of the UV boundary seam and should split the face into two bands.
    let close_sq = uv_polyline_trim_closed_len_sq_from_uv_poly(poly);
    let is_closed_trim = (trim_start - trim_end).length_squared() < close_sq;
    if is_closed_trim {
        // 閳光偓閳光偓 INTERIOR CLOSED LOOP 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
        // The trim is a truly closed loop entirely inside the polygon.
        // Don't split by closed trims 閳?return the original polygon unchanged.
        // The closed trim will be detected as an inner wire (hole) during sub-face
        // construction below, avoiding overlapping UV polygons that would cause
        // double-counting in surface area computation.
        let trim_centroid = trim.iter().copied().sum::<DVec2>() / trim.len() as f64;
        if point_in_polygon_2d(poly, trim_centroid) {
            return vec![poly.to_vec()];
        }
        return vec![poly.to_vec()];
    }

    // For each trim endpoint: if it lies close to the boundary already, use closest_on_boundary.
    // Otherwise, extrapolate along the trim tangent to find the proper boundary edge.
    let locate_endpoint =
        |endpoint: DVec2, tangent_from: DVec2| -> (usize, DVec2) {
            let (_, _, proj) = closest_on_boundary(endpoint);
            let dist_to_bnd = (endpoint - proj).length();
            if dist_to_bnd <= boundary_snap_tol {
                // Already on/near boundary
                let (edge, _, pt) = closest_on_boundary(endpoint);
                (edge, pt)
            } else {
                // Interior endpoint 閳?cast ray along trim tangent toward boundary
                let tang = (endpoint - tangent_from).normalize_or_zero();
                if let Some((edge, pt)) = ray_to_boundary(endpoint, tang) {
                    (edge, pt)
                } else {
                    // Fallback to closest projection
                    let (edge, _, pt) = closest_on_boundary(endpoint);
                    (edge, pt)
                }
            }
        };

    let interior_from_start = if trim.len() >= 2 { trim[1] } else { trim_end };
    let interior_from_end = if trim.len() >= 2 { trim[trim.len() - 2] } else { trim_start };

    let (edge_s, pt_s) = locate_endpoint(trim_start, interior_from_start);
    let (edge_e, pt_e) = locate_endpoint(trim_end, interior_from_end);

    // Ensure ia <= ib for consistent polygon walking
    let (ia, ib, p_a, p_b, trim_forward) = if edge_s <= edge_e {
        (edge_s, edge_e, pt_s, pt_e, true)
    } else {
        (edge_e, edge_s, pt_e, pt_s, false)
    };

    eprintln!("[DBG_SPLIT] poly={:?} n={}", poly, poly.len());
    eprintln!("[DBG_SPLIT] trim_start={:?} trim_end={:?}", trim_start, trim_end);
    eprintln!("[DBG_SPLIT] edge_s={} edge_e={} ia={} ib={}", edge_s, edge_e, ia, ib);
    eprintln!("[DBG_SPLIT] p_a={:?} p_b={:?}", p_a, p_b);

    // If both endpoints project to the same edge, inserting them as polygon
    // vertices creates distinct sub-edges that the standard split can handle
    // without self-overlapping sub-polygons.
    if ia == ib {
        let edge_a = poly[ia];
        let edge_b = poly[(ia + 1) % n];
        let edge_vec = edge_b - edge_a;
        let edge_len_sq = edge_vec.dot(edge_vec);
        if edge_len_sq > TOLERANCE_FLOAT_LOOSE && (p_a - p_b).length_squared() > TOLERANCE_FLOAT_ULTRA {
            let t_a = ((p_a - edge_a).dot(edge_vec) / edge_len_sq).clamp(0.0, 1.0);
            let t_b = ((p_b - edge_a).dot(edge_vec) / edge_len_sq).clamp(0.0, 1.0);
            let (p_first, p_second) = if t_a <= t_b { (p_a, p_b) } else { (p_b, p_a) };
            let mut new_poly = poly[..=ia].to_vec();
            new_poly.push(p_first);
            new_poly.push(p_second);
            new_poly.extend_from_slice(&poly[ia + 1..]);
            return split_uv_polygon_by_trim(&new_poly, trim);
        }
        // Degenerate: endpoints are coincident 鈥?no split possible, return original.
        return vec![poly.to_vec()];
    }

    // Build the trim points in the correct order for each half
    let trim_pts: Vec<DVec2> = if trim_forward {
        trim.to_vec()
    } else {
        trim.iter().copied().rev().collect()
    };

    // Detect wrap-around: trim u-span significantly exceeds polygon u-span.
    // When a trim wraps around the periodic domain, including the full interior
    // in both sub-polygons makes them overlap.  We split the trim at the polygon's
    // u-midpoint (u=0 for [-Pi,Pi]) so Sub A gets the left portion and Sub B
    // gets the right portion, matching the boundary paths they already use.
    let poly_u_min = poly.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let poly_u_max = poly.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let poly_u_span = poly_u_max - poly_u_min;
    let trim_u_min = trim_pts.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
    let trim_u_max = trim_pts.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
    let is_wrap_around = poly_u_span > 0.0 && (trim_u_max - trim_u_min) > poly_u_span * 0.8;

    // Find the index where the trim crosses the polygon's u-midpoint.
    // For a monotonic wrap-around trim, there is exactly one crossing.
    let u_mid = (poly_u_min + poly_u_max) / 2.0;
    let mut split_idx: Option<usize> = None;
    if is_wrap_around {
        for i in 0..trim_pts.len().saturating_sub(1) {
            let u0 = trim_pts[i].x;
            let u1 = trim_pts[i + 1].x;
            if (u0 - u_mid).abs() <= TOLERANCE_COORD_SUB {
                split_idx = Some(i);
                break;
            }
            if (u0 < u_mid && u1 > u_mid) || (u0 > u_mid && u1 < u_mid) {
                // The crossing is between points i and i+1; use i+1 as the split
                split_idx = Some(i + 1);
                break;
            }
        }
    }

    // ✅ OCCT-aligned: 瀛愬杈瑰舰鍙寘鍚?trim 鐨勭鐐?宸叉姇褰卞埌杈圭晫),涓嶅寘鍚唴閮ㄧ偣銆?
    //    OCCT 鐨?BOPAlgo_BuilderFace 鐢?MakeBlocks 鐢熸垚鐨?section edge
    //    (姣忔潯杈逛笉鍒嗘)鐩存帴鏋勫缓闈㈢嚎妗嗐€俽cad 鐨?split_uv_polygon_by_trim
    //    濡傛灉鎶?trim 鍐呴儴鐐归兘澶嶅埗杩涘瓙澶氳竟褰?姣忎釜 trim 浼氳础鐚鏉¤竟(3鐐光啋2杈?
    //    65鐐光啋64杈?,鑰屼笉鏄?OCCT 鐨?1 section edge / 鏇茬嚎銆?
    //    Sub-polygon A: poly[0..=ia] + p_a + p_b + poly[ib+1..]
    let mut sub_a: Vec<DVec2> = poly[..=ia].to_vec();
    sub_a.push(p_a);
    if let Some(si) = split_idx {
        if si > 0 && si < trim_pts.len() {
            sub_a.push(trim_pts[si]); // split point shared with Sub B
        } else {
            sub_a.push(p_b);
        }
    } else {
        sub_a.push(p_b);
    }
    sub_a.push(p_b);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // ✅ OCCT-aligned: 瀛愬杈瑰舰 B 涓嶅惈 trim 鍐呴儴鐐广€?
    //    Sub-polygon B: p_a + poly[ia+1..=ib] + p_b
    let mut sub_b: Vec<DVec2> = vec![p_a];
    sub_b.extend_from_slice(&poly[ia + 1..=ib]);
    sub_b.push(p_b);

    // Deduplicate consecutive near-equal points
    let dedup_2d = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > TOLERANCE_FLOAT_ULTRA {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < TOLERANCE_FLOAT_ULTRA {
            result.pop();
        }
        result
    };

    let sub_a = dedup_2d(sub_a);
    let sub_b = dedup_2d(sub_b);

    eprintln!("[DBG_SPLIT] sub_a: {} pts, sub_b: {} pts", sub_a.len(), sub_b.len());
    if sub_a.len() >= 3 {
        let u_min = sub_a.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = sub_a.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let v_min = sub_a.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = sub_a.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        eprintln!("[DBG_SPLIT] sub_a: u=[{:.6}, {:.6}] v=[{:.6}, {:.6}]", u_min, u_max, v_min, v_max);
    }
    if sub_b.len() >= 3 {
        let u_min = sub_b.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = sub_b.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let v_min = sub_b.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = sub_b.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        eprintln!("[DBG_SPLIT] sub_b: u=[{:.6}, {:.6}] v=[{:.6}, {:.6}]", u_min, u_max, v_min, v_max);
    }

    let sub_a_deduped = dedup_2d(sub_a);
    let sub_b_deduped = dedup_2d(sub_b);

    // ✅ OCCT-aligned: 濡傛灉瀛愬杈瑰舰閫€鍖?<3椤剁偣),杩斿洖鍘熷澶氳竟褰€?
    //    鍙戠敓鍦╰rim涓庡杈瑰舰杈圭晫閲嶅悎鏃?濡傚懆鏈熸€ф煴闈=2蟺鐨勮竟)銆?
    let sub_a_valid = sub_a_deduped.len() >= 3;
    let sub_b_valid = sub_b_deduped.len() >= 3;

    if sub_a_valid && sub_b_valid {
        vec![sub_a_deduped, sub_b_deduped]
    } else if sub_a_valid {
        vec![sub_a_deduped]
    } else if sub_b_valid {
        vec![sub_b_deduped]
    } else {
        vec![poly.to_vec()]
    }
}

/// Split a 2D polygon by a circle boundary.
///
/// Vertices inside the circle (distance < radius) are on the "inside" group,
/// vertices outside (distance > radius) are on the "outside" group.
/// Returns up to 2 sub-polygons: the part inside and the part outside.
///
/// When the circle is fully inside the polygon (all vertices outside),
/// samples the circle at N_CIRCLE_SAMPLES points and returns both
/// the approximate circular cap and the annular region.
/// Find the point where a segment [a, b] crosses a circle boundary.
/// `a` should be outside (sd > 0) and `b` inside (sd < 0) or vice versa.
/// Returns Some(crossing_point) or None if no valid crossing is found.
pub(crate) fn find_circle_segment_crossing(a: DVec2, b: DVec2, center: DVec2, radius: f64, tol: f64) -> Option<DVec2> {
    let ab = b - a;
    let ac = a - center;
    let qa = ab.dot(ab);
    if qa < 1e-30 { return None; }
    let qb = 2.0 * ab.dot(ac);
    let qc = ac.dot(ac) - radius * radius;
    let disc = qb * qb - 4.0 * qa * qc;
    if disc < 0.0 { return None; }
    let sq = disc.sqrt();
    for &sign in &[-1.0_f64, 1.0_f64] {
        let t = (-qb + sign * sq) / (2.0 * qa);
        if t > tol && t < 1.0 - tol {
            return Some(a + t * ab);
        }
    }
    None
}

pub(crate) fn split_polygon_by_circle_2d(poly: &[DVec2], center: DVec2, radius: f64, op: Option<BooleanOpType>) -> (Vec<Vec<DVec2>>, Vec<Vec<DVec2>>) {
    const N_CIRCLE_SAMPLES: usize = 24;
    let n = poly.len();
    if n < 3 {
        return (vec![poly.to_vec()], vec![]);
    }

    let tol = TOLERANCE_ABS;
    // If the circle center coincides with a polygon vertex, distance-to-circle and arc angles
    // degenerate; nudge the center slightly toward the polygon centroid (inside the face for
    // typical box/sphere trims) so segment閳ユ彿ircle intersections and arc sampling stay stable.
    let mut center = center;
    for &p in poly {
        if (p - center).length() < tol * 50.0 {
            let c0 = poly.iter().copied().fold(DVec2::ZERO, |a, q| a + q) / (n as f64);
            let dir = (c0 - center).normalize_or_zero();
            if dir.length_squared() > TOLERANCE_FLOAT_ULTRA {
                center = center + dir * (tol * 200.0).max(TOLERANCE_MESH_LEGACY);
                break;
            }
        }
    }

    // Signed distance: negative = inside circle, positive = outside
    let signed_dist = |p: DVec2| -> f64 { (p - center).length() - radius };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();

    // Check if all vertices are on the same side
    let all_inside = dists.iter().all(|&d| d <= tol);
    let all_outside = dists.iter().all(|&d| d >= -tol);

    if all_inside {
        // All polygon vertices inside circle 閳?keep whole polygon
        return (vec![poly.to_vec()], vec![]);
    }

    if all_outside {
        let center_in_poly = point_in_polygon_2d(poly, center);
        if center_in_poly {
            let circle_poly: Vec<DVec2> = (0..N_CIRCLE_SAMPLES)
                .map(|i| {
                    let theta = std::f64::consts::TAU * i as f64 / N_CIRCLE_SAMPLES as f64;
                    center + DVec2::new(theta.cos(), theta.sin()) * radius
                })
                .collect();
            let circle_fully_inside = circle_poly.iter().all(|&p| point_in_polygon_2d(poly, p));
            if circle_fully_inside {
                match op {
                    Some(BooleanOpType::Union) | Some(BooleanOpType::Difference) => {
                        // Keep polygon, subtract circle as hole (inner_wire).
                        return (vec![poly.to_vec()], vec![circle_poly]);
                    }
                Some(BooleanOpType::Intersection) => {
                    // For Intersection A鈭〣, the inner_wire (hole) represents the
                // region of A outside B. The caller's crossing split
                // produces the non-overlapping circle region separately.
                        return (vec![circle_poly], vec![]);
                    }
                    _ => {} // Other ops: fall through to crossing-based split
                }
            }
            // Circle extends beyond polygon boundary 鈥?clip the circle to the
            // polygon and use the clipped region as an inner wire (hole).
            // This avoids the N-crossing (N > 2) case in the crossing-based split,
            // which only handles exactly 2 crossings correctly.
            let clipped = clip_polygon_by_convex_polygon(&circle_poly, poly);
            if clipped.len() >= 3 {
                match op {
                    Some(BooleanOpType::Union) | Some(BooleanOpType::Difference) => {
                        // Outer polygon with clipped circle as hole.
                        return (vec![poly.to_vec()], vec![clipped]);
                    }
                    Some(BooleanOpType::Intersection) => {
                        // For Intersection, the clipped circle IS the result.
                        return (vec![clipped], vec![]);
                    }
                    _ => {} // Fall through to crossing-based split as backup
                }
            }
            // If clipping failed (degenerate result), fall through to
            // crossing-based split with same-edge crossing detection.
        }
    }

    // Find crossings: edges where signed distance changes sign
    let mut crossings: Vec<(usize, DVec2)> = Vec::new();

    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];

        let on_i = di.abs() < tol;
        let on_j = dj.abs() < tol;

        // Both on circle 閳?edge lies on boundary, no crossing
        if on_i && on_j {
            continue;
        }
                // 鉁?OCCT_aligned: Handle one vertex on circle, the other not on circle
        //   BOPAlgo_BuilderFace uses Hatcher to split parametric domain with 2D pcurves,
        //   correctly handling vertices on cutting curves.
        //
        //   When the non-on-circle vertex is INSIDE (di/dj < -tol), the crossing IS at
        //   the on-circle vertex; record it directly.
        //   When the non-on-circle vertex is OUTSIDE (di/dj > tol), check edge midpoint:
        //     midpoint inside 鈫?edge pierces through circle interior, find interior crossing
        //     midpoint outside 鈫?crossing is at the on-circle vertex
        //  Before fix: INSIDE鈫扥N and ON鈫扞NSIDE edges missed crossings because the
        //    midpoint was inside the circle but both (mid, end) or (start, mid) were
        //    fully inside.
        if on_i && !on_j {
            if dj < -tol {
                // poly[j] is INSIDE the circle: crossing at on-circle vertex poly[i]
                crossings.push((i, poly[i]));
            } else if dj > tol {
                // poly[j] is OUTSIDE the circle: check for interior crossing
                let mid = (poly[i] + poly[j]) * 0.5;
                if signed_dist(mid) < -tol {
                    // Edge goes from on-circle INTO circle, then back out.
                    if let Some(pt) = find_circle_segment_crossing(mid, poly[j], center, radius, tol) {
                        crossings.push((i, pt));
                    }
                } else {
                    crossings.push((i, poly[i]));
                }
            }
            continue;
        }
        if !on_i && on_j {
            if di < -tol {
                // poly[i] is INSIDE the circle: crossing at on-circle vertex poly[j]
                crossings.push((i, poly[j]));
            } else if di > tol {
                // poly[i] is OUTSIDE the circle: check for interior crossing
                let mid = (poly[i] + poly[j]) * 0.5;
                if signed_dist(mid) < -tol {
                    // Edge goes from outside INTO circle, then reaches on-circle vertex.
                    if let Some(pt) = find_circle_segment_crossing(poly[i], mid, center, radius, tol) {
                        crossings.push((i, pt));
                    }
                } else {
                    crossings.push((i, poly[j]));
                }
            }
            continue;
        }

        if di * dj < 0.0 {
            // Edge crosses the circle boundary
            // Find exact crossing: solve |a + t*(b-a) - center|铏?= r铏?
            let a = poly[i];
            let b = poly[j];
            let ab = b - a;
            let ac = a - center;
            let qa = ab.dot(ab);
            let qb = 2.0 * ab.dot(ac);
            let qc = ac.dot(ac) - radius * radius;
            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 {
                continue;
            }
            let sq = disc.sqrt();
            for &sign in &[-1.0_f64, 1.0_f64] {
                let t = (-qb + sign * sq) / (2.0 * qa);
                if t > -tol && t < 1.0 + tol {
                    let t = t.clamp(0.0, 1.0);
                    let pt = a + t * ab;
                    crossings.push((i, pt));
                    break; // take the first valid crossing on this edge
                }
            }
        }
    }

    // Check for all_outside + center_inside (Union) case where endpoints
    // are all outside the circle but edges need crossing detection via midpoint.
    if crossings.len() < 2 && all_outside && point_in_polygon_2d(poly, center) {
        let mut ec: Vec<(usize, DVec2)> = Vec::new();
        for ei in 0..n {
            let ej = (ei + 1) % n;
            let mid = (poly[ei] + poly[ej]) * 0.5;
            if signed_dist(mid) < -tol {
                // Both endpoints are outside the circle, but the edge passes through
                // it (midpoint inside). Find BOTH crossings: entry (start鈫抦id) and
                // exit (mid鈫抏nd). This gives 2 crossings on the same edge.
                if let Some(pt) = find_circle_segment_crossing(poly[ei], mid, center, radius, tol) {
                    ec.push((ei, pt));
                }
                if let Some(pt) = find_circle_segment_crossing(mid, poly[ej], center, radius, tol) {
                    ec.push((ei, pt));
                }
            }
        }
        if ec.len() >= 2 {
            crossings = ec;
        }
    }

    if crossings.len() < 2 {
        // Can't split 閳?keep as-is
        return (vec![poly.to_vec()], vec![]);
    }

    // Sort crossings by edge index
    crossings.sort_by_key(|(idx, _)| *idx);

    // Deduplicate crossings at the same spatial position (degenerate on-circle vertices
    // can produce crossings on both adjacent edges at the same point).
    crossings.dedup_by(|a, b| (a.1 - b.1).length_squared() < tol * tol);

    if crossings.len() < 2 {
        return (vec![poly.to_vec()], vec![]);
    }

    // N > 2 crossings with all_outside + center_inside: the polygon completely
    // encircles the circle.  Build the inner wire (clipped region) from crossings:
    //   - polygon edge segments between crossings on the same edge (inside circle)
    //   - circle arcs between crossings on different edges (inside polygon)
    if crossings.len() > 2 && all_outside && point_in_polygon_2d(poly, center) {
        // Group crossings by edge index.
        let mut inner: Vec<DVec2> = Vec::new();
        for ci in 0..crossings.len() {
            let (e_i, pt_a) = crossings[ci];
            let (e_j, pt_b) = crossings[(ci + 1) % crossings.len()];
            if e_i == e_j {
                // Both crossings on the same polygon edge 鈥?the edge segment
                // between them is inside the circle. Add the polygon vertices
                // on this segment, starting at pt_a and ending at pt_b.
                inner.push(pt_a);
                let e_end = poly[(e_i + 1) % n];
                let e_start = poly[e_i];
                let evec = e_end - e_start;
                let elen2 = evec.length_squared();
                let t_a = if elen2 > 1e-30 { (pt_a - e_start).dot(evec) / elen2 } else { 0.0 };
                let t_b = if elen2 > 1e-30 { (pt_b - e_start).dot(evec) / elen2 } else { 0.0 };
                // Add polygon vertices between pt_a and pt_b (sorted by t).
                let (t_lo, t_hi, rev) = if t_a < t_b { (t_a, t_b, false) } else { (t_b, t_a, true) };
                let _mids: Vec<DVec2> = Vec::new();
                // Check the polygon edge for interior vertices (not the crossing points).
                let _vi = e_i;
                let _vj = (e_i + 1) % n;
                // Walk the polygon from vi to vj, collecting vertex parameters.
                let mut verts_on_edge: Vec<(f64, DVec2)> = Vec::new();
                // The endpoints are crossings 鈥?don't add them here.
                verts_on_edge.push((t_lo, pt_a));
                // Find interior vertices of the polygon edge.
                // The polygon edge is from poly[vi] to poly[vj].
                // Interior vertices would be at t between 0 and 1.
                // For a polygon edge, the only vertices are vi and vj (no interior vertices).
                verts_on_edge.push((t_hi, pt_b));
                verts_on_edge.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                if rev { verts_on_edge.reverse(); }
                for (_, p) in verts_on_edge.iter().skip(1) {
                    inner.push(*p);
                }
            } else {
                // Crossings on different edges 鈥?the circle arc between them is
                // inside the polygon.  Sample 12 points on this arc.
                let a1 = (pt_a - center).to_angle();
                let a2 = (pt_b - center).to_angle();
                let d_ccw = (a2 - a1 + std::f64::consts::TAU) % std::f64::consts::TAU;
                const N_ARC: usize = 12;
                for k in 1..=N_ARC {
                    let t = k as f64 / N_ARC as f64;
                    let ang = a1 + d_ccw * t;
                    inner.push(center + DVec2::new(ang.cos(), ang.sin()) * radius);
                }
            }
        }
        if inner.len() >= 3 {
            return (vec![poly.to_vec()], vec![inner]);
        }
    }

    // Take the first two crossings
    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];

    if idx1 == idx2 {
        // Both crossings on the same polygon edge. Create inside (circle-interior)
        // and outside (circle-exterior) sub-polygons.

        // Determine which crossing is closer to edge start (poly[idx1])
        // versus edge end (poly[(idx1+1)%n]).
        let e_end = poly[(idx1 + 1) % n];
        let e_start = poly[idx1];
        let evec = e_end - e_start;
        let elen2 = evec.length_squared();
        let t_pt1 = if elen2 > 1e-30 { (pt1 - e_start).dot(evec) / elen2 } else { 0.0 };
        let t_pt2 = if elen2 > 1e-30 { (pt2 - e_start).dot(evec) / elen2 } else { 0.0 };
        let (near_start, near_end) = if t_pt1 < t_pt2 { (pt1, pt2) } else { (pt2, pt1) };

        // Interior arc: near_start 鈫?near_end through inner_mid_theta (circle interior side).
        // The chord midpoint points from center toward the chord 鈥?the arc nearest the chord
        // is the interior (smaller) arc, which is the circle-interior side.
        let chord_mid = (near_start + near_end) * 0.5;
        let inner_mid_theta = (chord_mid - center).to_angle();
        let theta_start = (near_start - center).to_angle();
        let theta_end = (near_end - center).to_angle();
        let int_delta = {
            let mut d = theta_end - theta_start;
            let go_ccw = if theta_start < theta_end {
                inner_mid_theta > theta_start && inner_mid_theta < theta_end
            } else {
                inner_mid_theta > theta_start || inner_mid_theta < theta_end
            };
            if go_ccw {
                while d < 0.0 { d += std::f64::consts::TAU; }
                if d > std::f64::consts::TAU { d -= std::f64::consts::TAU; }
            } else {
                while d > 0.0 { d -= std::f64::consts::TAU; }
                if d < -std::f64::consts::TAU { d += std::f64::consts::TAU; }
            }
            d
        };
        let int_arc_n = ((N_CIRCLE_SAMPLES as f64 * int_delta.abs() / std::f64::consts::TAU)
            as usize).max(3);
        let interior_arc: Vec<DVec2> = (0..=int_arc_n)
            .map(|i| {
                let t = i as f64 / int_arc_n as f64;
                let theta = theta_start + int_delta * t;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect();

        // Inside sub-polygon (circular segment = chord + interior arc):
        // near_start 鈫?interior_arc 鈫?near_end (chord closes implicitly).
        let mut sub_inside: Vec<DVec2> = Vec::new();
        sub_inside.push(near_start);
        for &p in interior_arc.iter().skip(1) {
            let last = *sub_inside.last().unwrap();
            if (p - last).length_squared() > TOLERANCE_FLOAT_ULTRA {
                sub_inside.push(p);
            }
        }

        // Outside sub-polygon: near_start 鈫?backward polygon walk 鈫?near_end
        // 鈫?interior_arc_rev (closing through the large/exterior arc).
        let mut sub_outside: Vec<DVec2> = Vec::new();
        sub_outside.push(near_start);
        // Walk polygon vertices backward from idx1 (through idx1-1, idx1-2, ...,
        // wrapping around to idx1+1).  This is the long path from near_start
        // to near_end that stays outside the circle.
        for k in 0..n {
            let vi = (idx1 + n - k) % n;
            let v = poly[vi];
            let last = *sub_outside.last().unwrap();
            if (v - last).length_squared() > TOLERANCE_FLOAT_ULTRA {
                sub_outside.push(v);
            }
        }
        // Add near_end on edge idx1 (closer to poly[idx1+1]).
        {
            let last = *sub_outside.last().unwrap();
            if (near_end - last).length_squared() > TOLERANCE_FLOAT_ULTRA {
                sub_outside.push(near_end);
            }
        }
        // Add interior_arc reversed (near_end 鈫?... 鈫?near_start through the
        // large/exterior arc) to close the outside polygon.
        for &p in interior_arc.iter().rev() {
            let last = *sub_outside.last().unwrap();
            if (p - last).length_squared() > TOLERANCE_FLOAT_ULTRA {
                sub_outside.push(p);
            }
        }

        // Dedup consecutive near-coincident vertices and trailing-first match
        let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
            let mut result: Vec<DVec2> = Vec::new();
            for p in v {
                if result.is_empty() || (p - result[result.len() - 1]).length_squared() > TOLERANCE_FLOAT_ULTRA {
                    result.push(p);
                }
            }
            if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < TOLERANCE_FLOAT_ULTRA {
                result.pop();
            }
            result
        };
        let sub_inside = dedup(sub_inside);
        let sub_outside = dedup(sub_outside);

        let mut out = Vec::new();
        if sub_inside.len() >= 3 { out.push(sub_inside); }
        if sub_outside.len() >= 3 { out.push(sub_outside); }

        return if out.is_empty() { (vec![poly.to_vec()], vec![]) } else { (out, vec![]) };
    }

    // Sample the arc between pt1 and pt2 (going through the inside of the polygon)
    // Determine which arc (minor or major) connects pt1 to pt2 and stays inside the polygon
    let theta1 = (pt1 - center).to_angle();
    let theta2 = (pt2 - center).to_angle();

    // For the "inside" sub-polygon, we need the arc that passes through the inside of the polygon.
    // Try both arcs and pick the one whose midpoint is inside the polygon.
    let mid_theta_cw = (theta1 + theta2) * 0.5;
    let mid_theta_ccw = mid_theta_cw + std::f64::consts::PI;
    let mid_cw = center + DVec2::new(mid_theta_cw.cos(), mid_theta_cw.sin()) * radius;
    let _mid_ccw = center + DVec2::new(mid_theta_ccw.cos(), mid_theta_ccw.sin()) * radius;

    // The arc midpoint that is inside the polygon corresponds to the "inside" portion
    let arc_goes_cw_inside = point_in_polygon_2d(poly, mid_cw);
    let inner_mid_theta = if arc_goes_cw_inside {
        mid_theta_cw
    } else {
        mid_theta_ccw
    };

    // Determine angular span and direction for the inner arc
    let arc_n = ((N_CIRCLE_SAMPLES as f64 * (theta2 - theta1).abs() / std::f64::consts::TAU)
        as usize)
        .max(3);

    // Build arc points from pt1 to pt2 going through inner_mid_theta
    let inner_arc: Vec<DVec2> = {
        // Compute proper arc from theta1 through inner_mid_theta to theta2
        let delta = {
            let mut d = theta2 - theta1;
            // Adjust delta to go through inner_mid_theta.
            // inner_mid_theta is the arc waypoint inside the polygon.
            // The CCW arc from theta1 to theta2:
            //   if theta1 < theta2: spans [theta1, theta2]
            //   if theta1 > theta2: wraps around 閳?[theta1, 2锜? 閳?[0, theta2]
            let going_ccw = if theta1 < theta2 {
                inner_mid_theta > theta1 && inner_mid_theta < theta2
            } else {
                inner_mid_theta > theta1 || inner_mid_theta < theta2
            };
            if going_ccw {
                while d < 0.0 {
                    d += std::f64::consts::TAU;
                }
                if d > std::f64::consts::TAU {
                    d -= std::f64::consts::TAU;
                }
            } else {
                while d > 0.0 {
                    d -= std::f64::consts::TAU;
                }
                if d < -std::f64::consts::TAU {
                    d += std::f64::consts::TAU;
                }
            }
            d
        };
        (0..=arc_n)
            .map(|i| {
                let t = i as f64 / arc_n as f64;
                let theta = theta1 + delta * t;
                center + DVec2::new(theta.cos(), theta.sin()) * radius
            })
            .collect()
    };

    // Sub-polygon "inside" (circle side): pt1 閳?arc 閳?pt2 + polygon walk from idx2 to idx1
    // Actually: vertices of polygon that are INSIDE the circle + arc from pt1 to pt2
    let poly_inside_verts: Vec<DVec2> = poly[idx1 + 1..=idx2].to_vec();

    let mut sub_inside: Vec<DVec2> = vec![pt1];
    sub_inside.extend_from_slice(&poly_inside_verts);
    // Avoid duplicating pt2 when it's already the last element of poly_inside_verts
    // (happens when pt2 is at a polygon vertex, e.g. an on-circle vertex).
    if poly_inside_verts.last() != Some(&pt2) {
        sub_inside.push(pt2);
    }
    // Add arc back (reversed, so the boundary goes: inside polygon verts, then arc back to pt1)
    for &p in inner_arc.iter().skip(1).rev().skip(1) {
        sub_inside.push(p);
    }

    // Sub-polygon "outside" (non-circle side): pt2 閳?arc 閳?pt1 + polygon walk
    let poly_outside_verts_a: Vec<DVec2> = poly[..=idx1].to_vec();
    let poly_outside_verts_b: Vec<DVec2> = poly[idx2 + 1..].to_vec();

    let mut sub_outside: Vec<DVec2> = poly_outside_verts_a;
    // Avoid duplicating pt1 when it's already the last element of poly_outside_verts_a
    if sub_outside.last() != Some(&pt1) {
        sub_outside.push(pt1);
    }
    // Add inner arc forward (pt1 鈫?pt2) as the closing boundary.
    // The sub_inside polygon uses the arc REVERSED (pt2 鈫?pt1), so sub_outside
    // must use the FORWARD direction (pt1 鈫?pt2) to create a non-self-intersecting
    // boundary that correctly encloses the non-circle-side region.
    // Using the reversed arc here would cause self-intersecting sub_outside polygons
    // when the circle crossings are at corner vertices (e.g. sphere-plane cut at origin
    // corner of a box where the arc passes through two corners of the face).
    let n_arc = inner_arc.len();
    for &p in inner_arc.iter().skip(1).take(n_arc.saturating_sub(2)) {
        sub_outside.push(p);
    }
    // Avoid duplicating pt2 when it's already the last element added, or when
    // it would duplicate the first element of poly_outside_verts_b
    if sub_outside.last() != Some(&pt2) && poly_outside_verts_b.first() != Some(&pt2) {
        sub_outside.push(pt2);
    }
    sub_outside.extend(poly_outside_verts_b);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > TOLERANCE_FLOAT_ULTRA {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < TOLERANCE_FLOAT_ULTRA {
            result.pop();
        }
        result
    };

    let sub_inside = dedup(sub_inside);
    let sub_outside = dedup(sub_outside);

    let mut out = Vec::new();
    if sub_inside.len() >= 3 {
        out.push(sub_inside);
    }
    if sub_outside.len() >= 3 {
        out.push(sub_outside);
    }

    if out.is_empty() {
        (vec![poly.to_vec()], vec![])
    } else {
        (out, vec![])
    }
}


/// Clip a subject polygon against a convex clip polygon using Sutherland鈥揌odgman.
///
/// Both polygons are assumed to be in 2D, with vertices ordered CCW.
/// The result is the intersection of the two polygons (also CCW).
pub(crate) fn clip_polygon_by_convex_polygon(subject: &[DVec2], clip: &[DVec2]) -> Vec<DVec2> {
    if subject.len() < 3 || clip.len() < 3 {
        return Vec::new();
    }
    let tol = TOLERANCE_ABS;
    let mut result: Vec<DVec2> = subject.to_vec();
    let nclip = clip.len();
    for ci in 0..nclip {
        if result.is_empty() {
            return Vec::new();
        }
        let cj = (ci + 1) % nclip;
        let edge_start = clip[ci];
        let edge_end = clip[cj];
        let edge = edge_end - edge_start;

        let mut next_ring: Vec<DVec2> = Vec::new();
        let nsub = result.len();
        for si in 0..nsub {
            let sj = (si + 1) % nsub;
            let current = result[si];
            let next = result[sj];

            // Inside test: cross product (edge 脳 (P - edge_start)) >= 0
            // For a CCW clip polygon, interior is to the LEFT of each edge.
            let inside_curr = edge.perp_dot(current - edge_start) >= -tol;
            let inside_next = edge.perp_dot(next - edge_start) >= -tol;

            if inside_curr {
                next_ring.push(current);
            }
            if inside_curr != inside_next {
                // Edge crosses the clipping boundary 鈥?find intersection point
                let delta = next - current;
                let num = edge.perp_dot(current - edge_start);
                let den = edge.perp_dot(delta);
                if den.abs() > TOLERANCE_FLOAT_ULTRA {
                    let t = -num / den;
                    let t = t.clamp(0.0, 1.0);
                    next_ring.push(current + delta * t);
                }
            }
        }
        result = next_ring;
    }
    // Dedup near-coincident consecutive vertices
    let mut deduped: Vec<DVec2> = Vec::with_capacity(result.len());
    for p in &result {
        if deduped.is_empty()
            || (*p - *deduped.last().unwrap()).length_squared() > TOLERANCE_FLOAT_ULTRA
        {
            deduped.push(*p);
        }
    }
    if deduped.len() > 1
        && (deduped[0] - *deduped.last().unwrap()).length_squared() < TOLERANCE_FLOAT_ULTRA
    {
        deduped.pop();
    }
    deduped
}

/// Check if a 2D point is inside a 2D polygon using ray casting.
pub(crate) fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Insert imprint points that fall on (or very near) polygon edges so ResultBuilder wires share
/// vertices along coplanar seams instead of creating overlapping segments with T-junctions.
pub(crate) fn insert_points_on_polygon_edges(poly: &[DVec2], imprint: &[DVec2], tol: f64) -> Vec<DVec2> {
    let n = poly.len();
    if n < 3 {
        return poly.to_vec();
    }
    let mut out: Vec<DVec2> = Vec::with_capacity(n + imprint.len());
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        out.push(a);
        let mut splits: Vec<(f64, DVec2)> = Vec::new();
        for &p in imprint {
            if let Some(t) = segment_closest_param_2d(a, b, p, tol)
                && t > tol && t < 1.0 - tol {
                    splits.push((t, a + (b - a) * t));
                }
        }
        splits.sort_by(|u, v| u.0.partial_cmp(&v.0).unwrap_or(std::cmp::Ordering::Equal));
        for (_, q) in splits {
            if out
                .last()
                .map(|last| (*last - q).length() > tol * 0.5)
                .unwrap_or(true)
            {
                out.push(q);
            }
        }
    }
    dedup_consecutive_poly2d(&out, tol)
}

/// Closest-point parameter t in [0,1] on segment ab if p is within `tol` of the segment.
pub(crate) fn segment_closest_param_2d(a: DVec2, b: DVec2, p: DVec2, tol: f64) -> Option<f64> {
    let ab = b - a;
    let l2 = ab.length_squared();
    if l2 < tol * tol {
        return None;
    }
    let t = ((p - a).dot(ab) / l2).clamp(0.0, 1.0);
    let closest = a + ab * t;
    // Lenient perpendicular tolerance: imprint projections can sit slightly off the segment
    // after mixed plane lifts (union box test).
    if (p - closest).length() <= tol * 200.0 {
        Some(t)
    } else {
        None
    }
}

pub(crate) fn dedup_consecutive_poly2d(poly: &[DVec2], tol: f64) -> Vec<DVec2> {
    if poly.is_empty() {
        return vec![];
    }
    let mut v: Vec<DVec2> = Vec::with_capacity(poly.len());
    for &p in poly {
        if v.is_empty() || (p - v[v.len() - 1]).length() > tol * 0.5 {
            v.push(p);
        }
    }
    if v.len() > 2 && (v[0] - v[v.len() - 1]).length() <= tol * 0.5 {
        v.pop();
    }
    v
}

/// Split a 2D polygon by an infinite line through `point` with direction `dir`.
///
/// Vertices on the positive side (cross product > 0) form one group, negative side the other.
pub(crate) fn split_polygon_2d_by_line(poly: &[DVec2], point: DVec2, dir: DVec2) -> Vec<Vec<DVec2>> {
    let n = poly.len();
    if n < 3 {
        return vec![poly.to_vec()];
    }
    let tol = TOLERANCE_ABS;

    // Signed distance from line
    let signed_dist = |p: DVec2| -> f64 {
        let d = p - point;
        dir.x * d.y - dir.y * d.x // perpendicular component
    };

    let dists: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();
    // Vertices exactly on the line (|d| < tol) are neutral 閳?they don't count as
    // "all on one side".  Only strictly positive (> tol) or strictly negative (< -tol)
    // vertices determine whether the polygon crosses the line.
    let all_pos = dists.iter().all(|&d| d > tol);
    let all_neg = dists.iter().all(|&d| d < -tol);

    if all_pos || all_neg {
        return vec![poly.to_vec()];
    }

    let mut crossings: Vec<(usize, DVec2)> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let di = dists[i];
        let dj = dists[j];

        // When a vertex lies exactly on the split line (|d| < tol), the original
        // edge-crossing test would skip the edge entirely, losing the crossing.
        // This happens when a circular face with few boundary vertices is split
        // by a line passing through two vertices (e.g. an inscribed square on a
        // cylinder cap, where box edge passes through circle polygon vertices).
        //
        // Fix: two cases:
        // 1. *Current* vertex on line, *next* off: search backward for the first
        //    non-on-line vertex. If its sign opposes the next vertex's sign, the
        //    line crosses at this vertex.
        // 2. *Next* vertex on line, *current* off: search forward from the next
        //    for the first non-on-line vertex. If its sign opposes the current
        //    vertex's sign, the line crosses at the next vertex.
        if di.abs() < tol && dj.abs() >= tol {
            let mut pi = (i + n - 1) % n;
            while pi != i && dists[pi].abs() < tol {
                pi = (pi + n - 1) % n;
            }
            if pi != i && dists[pi] * dj < 0.0 {
                crossings.push((i, poly[i]));
                continue;
            }
        }
        if di.abs() >= tol && dj.abs() < tol {
            let mut nj = (j + 1) % n;
            while nj != j && dists[nj].abs() < tol {
                nj = (nj + 1) % n;
            }
            if nj != j && di * dists[nj] < 0.0 {
                crossings.push((j, poly[j]));
                continue;
            }
        }

        if di.abs() < tol || dj.abs() < tol {
            continue;
        }
        if di * dj < 0.0 {
            let t = di / (di - dj);
            let pt = poly[i] + t * (poly[j] - poly[i]);
            crossings.push((i, pt));
        }
    }

    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    // Deduplicate: the forward-search and backward-search may both detect
    // a crossing at the same vertex from adjacent edges (e.g. a diamond
    // polygon split by a line through two opposite vertices).
    crossings.sort_by_key(|(idx, _)| *idx);
    crossings.dedup_by(|a, b| a.0 == b.0);
    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    let (idx1, pt1) = crossings[0];
    let (idx2, pt2) = crossings[1];
    if idx1 == idx2 {
        return vec![poly.to_vec()];
    }

    let mut sub_a: Vec<DVec2> = poly[..=idx1].to_vec();
    sub_a.push(pt1);
    sub_a.push(pt2);
    sub_a.extend_from_slice(&poly[idx2 + 1..]);

    let mut sub_b: Vec<DVec2> = vec![pt1];
    sub_b.extend_from_slice(&poly[idx1 + 1..=idx2]);
    sub_b.push(pt2);

    let dedup = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut result: Vec<DVec2> = Vec::new();
        for p in v {
            if result.is_empty() || (p - result[result.len() - 1]).length_squared() > TOLERANCE_FLOAT_ULTRA {
                result.push(p);
            }
        }
        if result.len() > 1 && (result[0] - result[result.len() - 1]).length_squared() < TOLERANCE_FLOAT_ULTRA {
            result.pop();
        }
        result
    };

    let sub_a = dedup(sub_a);
    let sub_b = dedup(sub_b);
    let mut out = Vec::new();
    if sub_a.len() >= 3 {
        out.push(sub_a);
    }
    if sub_b.len() >= 3 {
        out.push(sub_b);
    }
    if out.is_empty() {
        vec![poly.to_vec()]
    } else {
        out
    }
}

/// Split a 2D polygon by a segment from `seg_start` to `seg_end`.
pub(crate) fn split_polygon_2d_by_segment(
    poly: &[DVec2],
    seg_start: DVec2,
    seg_end: DVec2,
) -> Vec<Vec<DVec2>> {
    let dir = seg_end - seg_start;
    if dir.length_squared() < TOLERANCE_FLOAT_ULTRA {
        return vec![poly.to_vec()];
    }
    split_polygon_2d_by_line(poly, seg_start, dir.normalize())
}

// ============================================================================
// Glue Path Enhancement Types and Functions
// ============================================================================