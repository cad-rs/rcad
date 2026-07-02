
/// Get a point on the rect perimeter at a given absolute position.
/// `t_abs` ranges from 0 to `total_perim = 2*(w+h)`.
fn rect_perimeter_point(bmin: DVec2, bmax: DVec2, t_abs: f64) -> DVec2 {
    let w = bmax.x - bmin.x;
    let h = bmax.y - bmin.y;
    let perim = 2.0 * (w + h);
    let t = t_abs.rem_euclid(perim);
    if t <= w {
        DVec2::new(bmin.x + t, bmin.y)
    } else if t <= w + h {
        DVec2::new(bmax.x, bmin.y + (t - w))
    } else if t <= 2.0 * w + h {
        DVec2::new(bmax.x - (t - w - h), bmax.y)
    } else {
        DVec2::new(bmin.x, bmax.y - (t - 2.0 * w - h))
    }
}

/// Build a triangulated BRep for `box - cone` using Z-slice tessellation.
///
/// The box is axis-aligned `[bmin, bmax]`. The cone is a Z-aligned conical frustum
/// with center at `(cx, cy)` in XY, extending from Z `z_lo` to `z_hi`, with bottom
/// radius `r_lo` and top radius `r_hi`.
fn build_box_minus_cone_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    z_lo: f64, z_hi: f64,
    r_lo: f64, r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r_lo < tol && r_hi < tol { return None; }

    // Adjust to box Z-range
    let z0 = z_lo.max(bmin.z);
    let z1 = z_hi.min(bmax.z);
    if z1 <= z0 + tol { return None; }

    let n_slices = 128usize;
    let n_boundary = 256usize;
    let dz = (z1 - z0) / n_slices as f64;
    let dr = (r_hi - r_lo) / (z_hi - z_lo);

    let bmin2 = DVec2::new(bmin.x, bmin.y);
    let bmax2 = DVec2::new(bmax.x, bmax.y);

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // ---- Generate Z-slice boundary polygons ----
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n_slices + 1);
    for i in 0..=n_slices {
        let z = z0 + dz * i as f64;
        let r = r_lo + dr * (z - z_lo);
        if r <= tol {
            slices.push(vec![]);
        } else {
            let poly = rect_minus_circle_boundary(bmin2, bmax2, cx, cy, r, n_boundary);
            slices.push(poly);
        }
    }

    // ---- Remap each boundary to n_boundary equally-spaced, aligned points ----
    // This ensures boundaries at all Z-levels have the same length and start
    // from the same physical reference point (bottom-left rect corner).
    let ref_pt = DVec2::new(bmin.x, bmin.y);
    for poly in &mut slices {
        if poly.len() < 3 { continue; }

        // Find the index of the point closest to the reference
        let mut best_idx = 0;
        let mut best_dist = (poly[0] - ref_pt).length_squared();
        for (idx, p) in poly.iter().enumerate() {
            let d = (*p - ref_pt).length_squared();
            if d < best_dist { best_dist = d; best_idx = idx; }
        }

        // Rotate the array to start from best_idx
        poly.rotate_left(best_idx);

        // Resample to exactly n_boundary equally-spaced points via arc-length
        let n_bnd = n_boundary.max(4);
        // Compute cumulative arc length
        let mut arc_len = vec![0.0_f64; poly.len() + 1];
        for i in 1..=poly.len() {
            let j = i % poly.len();
            let k = (i - 1) % poly.len();
            arc_len[i] = arc_len[i - 1] + (*poly)[k].distance((*poly)[j]);
        }
        let total = arc_len[poly.len()];
        if total <= tol { continue; }

        let mut new_poly = Vec::with_capacity(n_bnd);
        let mut src_idx = 0;
        for i in 0..n_bnd {
            let target = total * i as f64 / n_bnd as f64;
            while src_idx < poly.len() && arc_len[src_idx + 1] < target {
                src_idx += 1;
            }
            let t0 = arc_len[src_idx];
            let t1 = arc_len[(src_idx + 1) % (poly.len() + 1)];
            if (t1 - t0).abs() < 1e-15 {
                new_poly.push(poly[src_idx % poly.len()]);
            } else {
                let frac = (target - t0) / (t1 - t0);
                let a = src_idx % poly.len();
                let b = (src_idx + 1) % poly.len();
                new_poly.push(poly[a].lerp(poly[b], frac));
            }
        }
        *poly = new_poly;
    }

    // ---- Build wall faces between adjacent Z-slices ----
    for i in 0..n_slices {
        let bot = &slices[i];
        let top = &slices[i + 1];
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        // Both slices have boundaries �?build wall
        if !bot.is_empty() && !top.is_empty() {
            let n = bot.len().min(top.len());
            let mut idx = Vec::with_capacity(2 * n);
            for p in bot.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for p in top.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_top))); }

            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                let b0 = idx[j];
                let b1 = idx[k];
                let t0 = idx[n + j];
                let t1 = idx[n + k];
                tris.push([b0, b1, t1]);
                tris.push([b0, t1, t0]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        } else if !bot.is_empty() {
            // Top is empty (cone closed off at this Z) �?cap the top
            // Build a triangle fan from the last valid boundary to a center point
            let n = bot.len();
            let center_z = (z_bot + z_top) * 0.5;
            // Compute the center of the remaining region
            let mut center = DVec3::ZERO;
            for p in bot.iter() { center += DVec3::new(p.x, p.y, z_bot); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in bot.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        } else if !top.is_empty() {
            // Bottom is empty (cone opened at this Z) �?cap the bottom
            let n = top.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in top.iter() { center += DVec3::new(p.x, p.y, z_top); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in top.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_top))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        }
    }

    // ---- Build cap faces at z0 and z1 if boundary is non-empty ----
    // Bottom cap at z0 �?triangulated via ear-clipping
    if !slices[0].is_empty() && slices[0].len() >= 3 {
        let empty_wire = || Wire { edges: vec![] };
        let poly_3d: Vec<DVec3> = slices[0].iter()
            .map(|p| DVec3::new(p.x, p.y, z0))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, -DVec3::Z);

        // Remap vertex indices
        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    // Top cap at z1 �?triangulated via ear-clipping
    if !slices[n_slices].is_empty() && slices[n_slices].len() >= 3 {
        let empty_wire = || Wire { edges: vec![] };
        let poly_3d: Vec<DVec3> = slices[n_slices].iter()
            .map(|p| DVec3::new(p.x, p.y, z1))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path for `box �?cone` boolean difference.
///
/// Detects an axis-aligned box minus a Z-aligned conical frustum (possibly Z-rotated
/// and translated in XY). Builds the result via Z-slice tessellation.
pub fn try_difference_box_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    // Detect axis-aligned box (a)
    let [bmin, bmax] = try_as_axis_aligned_box(a)?;

    // Detect Z-aligned cone frustum (b)
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(b)?;

    // Compute Z overlap
    let z_lo = cz_lo.max(bmin.z);
    let z_hi = cz_hi.min(bmax.z);
    if z_hi <= z_lo + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    let cx = center_xy.x;
    let cy = center_xy.y;

    // Check if there's any XY overlap at all Z levels in the range
    // (quick check: if the cone is entirely outside the box XY at all Z)
    let _max_dist_xy = (bmax.x - cx).max(cx - bmin.x)
        .max((bmax.y - cy).max(cy - bmin.y));
    let min_r = if cr_lo < cr_hi { cr_lo } else { cr_hi };
    if min_r < TOLERANCE_LEN_MIN {
        // Sharp cone (radius near zero) �?can't form a proper void
        return None;
    }
    let r_at_zlo = cr_lo + (cr_hi - cr_lo) * (z_lo - cz_lo) / (cz_hi - cz_lo);
    let r_at_zhi = cr_lo + (cr_hi - cr_lo) * (z_hi - cz_lo) / (cz_hi - cz_lo);
    let min_overlap_r = if r_at_zlo < r_at_zhi { r_at_zlo } else { r_at_zhi };

    // If the cone is always outside the box XY, no void
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    if dist_center > box_half_diag + min_overlap_r + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    // Clamp radii to non-negative
    let r_lo = r_at_zlo.max(TOLERANCE_COORD_SUB);
    let r_hi = r_at_zhi.max(TOLERANCE_COORD_SUB);

    // Analytic fast path (extracted to cone_box_analytic module)
    if let Some(result) = crate::cone_box_analytic::build_box_minus_cone_analytic(a, b) {
        return Some(result);
    }

    build_box_minus_cone_tessellated(bmin, bmax, cx, cy, z_lo, z_hi, r_lo, r_hi)
}

/// Build a triangulated BRep for `cone - box` using Z-slice tessellation.
///
/// The box is axis-aligned `[bmin, bmax]`. The cone is a Z-aligned conical frustum
/// with center at `(cx, cy)` in XY, extending from Z `z_lo` to `z_hi`, with bottom
/// radius `r_lo` and top radius `r_hi`.  This is the inverse of
/// [`build_box_minus_cone_tessellated`]: the kept shape is the cone, and the box
/// cuts a channel through it, so the cross-section is `circle - rect`.
fn build_cone_minus_box_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    z_lo: f64, z_hi: f64,
    r_lo: f64, r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r_lo < tol && r_hi < tol { return None; }

    // Adjust to box Z-range
    let z0 = z_lo.max(bmin.z);
    let z1 = z_hi.min(bmax.z);
    if z1 <= z0 + tol { return None; }

    let n_slices = 128usize;
    let n_arc = 128usize;
    let dz = (z1 - z0) / n_slices as f64;
    let dr = (r_hi - r_lo) / (z_hi - z_lo);

    let bmin2 = DVec2::new(bmin.x, bmin.y);
    let bmax2 = DVec2::new(bmax.x, bmax.y);
    let tau = std::f64::consts::TAU;

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // Helper: test if a point is inside the rectangle
    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin2.x - tol && p.x <= bmax2.x + tol
            && p.y >= bmin2.y - tol && p.y <= bmax2.y + tol
    };

    // Find circle-rect intersections and generate boundary.
    // Only generates points for circle arcs OUTSIDE the rect.
    // Arc segments INSIDE the rect are skipped �?the polygon's closing
    // edge from last-to-first serves as the return along the rect perimeter.
    let gen_boundary = |r: f64| -> Vec<DVec2> {
        if r <= tol { return vec![]; }

        let edges = [
            (DVec2::new(bmin2.x, bmin2.y), DVec2::new(bmax2.x, bmin2.y)),
            (DVec2::new(bmax2.x, bmin2.y), DVec2::new(bmax2.x, bmax2.y)),
            (DVec2::new(bmax2.x, bmax2.y), DVec2::new(bmin2.x, bmax2.y)),
            (DVec2::new(bmin2.x, bmax2.y), DVec2::new(bmin2.x, bmin2.y)),
        ];

        struct Intersection { t: f64, edge: usize, pt: DVec2 }
        let mut ints: Vec<Intersection> = Vec::new();

        for (ei, (p0, p1)) in edges.iter().enumerate() {
            let d = *p1 - *p0;
            let a0 = *p0 - DVec2::new(cx, cy);
            let A = d.dot(d);
            if A < 1e-30 { continue; }
            let B = 2.0 * a0.dot(d);
            let C = a0.dot(a0) - r * r;
            let disc = B * B - 4.0 * A * C;
            if disc < 0.0 { continue; }
            let sqrt_disc = disc.sqrt();
            for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
                if t >= -tol && t <= 1.0 + tol {
                    let tc = t.clamp(0.0, 1.0);
                    let pt = *p0 + d * tc;
                    ints.push(Intersection { t: tc, edge: ei, pt });
                }
            }
        }

        ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
        ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);
        ints.sort_by(|a, b| a.pt.x.partial_cmp(&b.pt.x).unwrap()
            .then(a.pt.y.partial_cmp(&b.pt.y).unwrap()));
        ints.dedup_by(|a, b| (a.pt - b.pt).length_squared() < tol * tol);

        // No-intersection cases
        if ints.is_empty() {
            let center_inside = cx >= bmin2.x - tol && cx <= bmax2.x + tol
                && cy >= bmin2.y - tol && cy <= bmax2.y + tol;
            let corners = [
                DVec2::new(bmin2.x, bmin2.y), DVec2::new(bmax2.x, bmin2.y),
                DVec2::new(bmax2.x, bmax2.y), DVec2::new(bmin2.x, bmax2.y),
            ];
            let any_corner_outside = corners.iter().any(|p| {
                (*p - DVec2::new(cx, cy)).length_squared() > r * r + tol
            });
            if center_inside && !any_corner_outside {
                return vec![];
            }
            // Full circle boundary (no overlap or rect inside circle)
            let mut result = Vec::with_capacity(n_arc * 2);
            for i in 0..n_arc * 2 {
                let ang = tau * i as f64 / (n_arc * 2) as f64;
                let (s, c) = ang.sin_cos();
                result.push(DVec2::new(cx + r * c, cy + r * s));
            }
            return result;
        }

        // Sort by CCW circle angle
        ints.sort_by(|a, b| {
            f64::atan2(a.pt.y - cy, a.pt.x - cx)
                .partial_cmp(&f64::atan2(b.pt.y - cy, b.pt.x - cx))
                .unwrap()
        });

        let m = ints.len();
        let mut result = Vec::new();

        for i in 0..m {
            let j = (i + 1) % m;
            let ei = &ints[i];
            let pj = &ints[j];

            let v1 = ei.pt - DVec2::new(cx, cy);
            let v2 = pj.pt - DVec2::new(cx, cy);
            let a1 = f64::atan2(v1.y, v1.x);
            let a2 = f64::atan2(v2.y, v2.x);
            let da_ccw = (a2 - a1).rem_euclid(tau);

            // Midpoint of CCW circle arc
            let mid_ang = a1 + da_ccw * 0.5;
            let mid_pt = DVec2::new(cx + r * mid_ang.cos(), cy + r * mid_ang.sin());

            if !point_in_rect(mid_pt) {
                // Arc outside rect �?sample with n_arc points (includes both endpoints)
                for k in 0..=n_arc {
                    let frac = k as f64 / n_arc as f64;
                    let ang = a1 + da_ccw * frac;
                    let (s, c) = ang.sin_cos();
                    result.push(DVec2::new(cx + r * c, cy + r * s));
                }
            }
            // Arc inside rect �?skip. The closing edge from the last polygon
            // vertex back to the first serves as the rect perimeter return path.
        }

        result
    };

    // ---- Generate Z-slice boundary polygons ----
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n_slices + 1);
    for i in 0..=n_slices {
        let z = z0 + dz * i as f64;
        let r = r_lo + dr * (z - z_lo);
        slices.push(gen_boundary(r));
    }

    // ---- Build wall faces between adjacent Z-slices ----
    for i in 0..n_slices {
        let bot = &slices[i];
        let top = &slices[i + 1];
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        if !bot.is_empty() && !top.is_empty() {
            let n = bot.len().min(top.len());
            let mut idx = Vec::with_capacity(2 * n);
            for p in bot.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for p in top.iter() { idx.push(add_v(DVec3::new(p.x, p.y, z_top))); }

            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                let b0 = idx[j];
                let b1 = idx[k];
                let t0 = idx[n + j];
                let t1 = idx[n + k];
                tris.push([b0, b1, t1]);
                tris.push([b0, t1, t0]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        } else if !bot.is_empty() {
            // Top is empty (void closed off) �?cap the top with triangle fan
            let n = bot.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in bot.iter() { center += DVec3::new(p.x, p.y, z_bot); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in bot.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_bot))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        } else if !top.is_empty() {
            // Bottom is empty (void opened at this Z) �?cap the bottom with triangle fan
            let n = top.len();
            let center_z = (z_bot + z_top) * 0.5;
            let mut center = DVec3::ZERO;
            for p in top.iter() { center += DVec3::new(p.x, p.y, z_top); }
            center /= n as f64;
            center.z = center_z;

            let mut tris = Vec::new();
            let c_idx = add_v(center);
            let mut ring = Vec::with_capacity(n);
            for p in top.iter() { ring.push(add_v(DVec3::new(p.x, p.y, z_top))); }
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([c_idx, ring[j], ring[k]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        }
    }

    // ---- Build cap faces at z0 and z1 if boundary is non-empty ----
    // Bottom cap at z0
    if !slices[0].is_empty() && slices[0].len() >= 3 {
        let poly_3d: Vec<DVec3> = slices[0].iter()
            .map(|p| DVec3::new(p.x, p.y, z0))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, -DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    // Top cap at z1
    if !slices[n_slices].is_empty() && slices[n_slices].len() >= 3 {
        let poly_3d: Vec<DVec3> = slices[n_slices].iter()
            .map(|p| DVec3::new(p.x, p.y, z1))
            .collect();
        let tris = crate::triangulate::triangulate_polygon(&poly_3d, DVec3::Z);

        let mut remapped_tris = Vec::with_capacity(tris.len());
        let local_verts: Vec<usize> = poly_3d.iter().map(|p| add_v(*p)).collect();
        for t in &tris {
            remapped_tris.push([local_verts[t[0]], local_verts[t[1]], local_verts[t[2]]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: remapped_tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Build `cylinder �?box` by Z-slice tessellation.
///
/// Handles the case where clip-plane corners fall outside the cylinder �?/// the gap routing in `build_cylinder_box_clipped_brep` cannot create correct
/// per-clip-plane side faces when the corner is outside the cylinder.
///
/// The cross-section `circle �?rect` at every Z-level is identical (constant
/// cylinder radius, constant box size). This function computes the 2D boundary
/// as one or more closed polygons (disconnected components) and builds a
/// triangulated BRep by connecting Z-slices.
///
/// Parameters are in the box's UV frame: the box is `[-eu, eu] �?[-ev, ev]`,
/// the circle center is at `(cu, cv)` with radius `r`, and the cylinder extends
/// vertically from `z_lo` to `z_hi`. The world-space position of a point `(u, v, z)`
/// is `bc + u * u_ax + v * v_ax + z * DVec3::Z`.
fn build_cylinder_box_diff_tessellated(
    bc: DVec3,
    u_ax: DVec3,
    v_ax: DVec3,
    cu: f64,
    cv: f64,
    r: f64,
    _h: f64,
    eu: f64,
    ev: f64,
    z_lo: f64,
    z_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if z_hi <= z_lo + tol { return None; }
    if r < tol { return None; }

    let n_slices = 64usize;
    let n_arc = 128usize;
    let dz = (z_hi - z_lo) / n_slices as f64;
    let tau = std::f64::consts::TAU;

    let bmin = DVec2::new(-eu, -ev);
    let bmax = DVec2::new(eu, ev);
    let cx = cu;
    let cy = cv;

    let empty_wire = || Wire { edges: vec![] };

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    // Transform box UV coords to world space.
    // z is a world Z coordinate, so we must NOT add bc.z.
    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(bc.x + u_ax.x * u + v_ax.x * v, bc.y + u_ax.y * u + v_ax.y * v, z)
    };

    // ---- 1. Find circle-rect intersections ----
    // Edge order: bottom (L閳�?, right (B閳�?, top (R閳�?, left (T閳�?
    let rect_edges = [
        (DVec2::new(bmin.x, bmin.y), DVec2::new(bmax.x, bmin.y)),
        (DVec2::new(bmax.x, bmin.y), DVec2::new(bmax.x, bmax.y)),
        (DVec2::new(bmax.x, bmax.y), DVec2::new(bmin.x, bmax.y)),
        (DVec2::new(bmin.x, bmax.y), DVec2::new(bmin.x, bmin.y)),
    ];

    struct Intersection { t: f64, edge: usize, pt: DVec2 }
    let mut ints: Vec<Intersection> = Vec::new();

    for (ei, (p0, p1)) in rect_edges.iter().enumerate() {
        let d = *p1 - *p0;
        let a0 = *p0 - DVec2::new(cx, cy);
        let A = d.dot(d);
        if A < 1e-30 { continue; }
        let B = 2.0 * a0.dot(d);
        let C = a0.dot(a0) - r * r;
        let disc = B * B - 4.0 * A * C;
        if disc < 0.0 { continue; }
        let sqrt_disc = disc.sqrt();
        for t in [(-B - sqrt_disc) / (2.0 * A), (-B + sqrt_disc) / (2.0 * A)] {
            if t >= -tol && t <= 1.0 + tol {
                let tc = t.clamp(0.0, 1.0);
                let pt = *p0 + d * tc;
                ints.push(Intersection { t: tc, edge: ei, pt });
            }
        }
    }

    if ints.is_empty() {
        return None;
    }

    // Same-edge dedup by t parameter.
    // Keep corner duplicates (same point on 2 edges) �?they prevent
    // per-edge routing gaps and the zero-length arc they create is
    // handled transparently in run grouping.
    ints.sort_by(|a, b| a.edge.cmp(&b.edge).then(a.t.partial_cmp(&b.t).unwrap()));
    ints.dedup_by(|a, b| a.edge == b.edge && (a.t - b.t).abs() < tol);

    // Sort by CCW angle around circle center
    ints.sort_by(|a, b| {
        f64::atan2(a.pt.y - cy, a.pt.x - cx)
            .partial_cmp(&f64::atan2(b.pt.y - cy, b.pt.x - cx))
            .unwrap()
    });

    let m = ints.len();
    if m < 2 { return None; }

    // ---- 2. Classify each CCW arc as KEPT or SKIPPED ----
    let point_in_rect = |p: DVec2| -> bool {
        p.x >= bmin.x - tol && p.x <= bmax.x + tol
            && p.y >= bmin.y - tol && p.y <= bmax.y + tol
    };

    let mut is_kept = vec![false; m];
    let mut arc_zero = vec![false; m];
    for i in 0..m {
        let j = (i + 1) % m;
        let v1 = ints[i].pt - DVec2::new(cx, cy);
        let v2 = ints[j].pt - DVec2::new(cx, cy);
        let a1 = f64::atan2(v1.y, v1.x);
        let a2 = f64::atan2(v2.y, v2.x);
        let da_ccw = (a2 - a1).rem_euclid(tau);
        arc_zero[i] = da_ccw < 1e-12;
        let mid_ang = a1 + da_ccw * 0.5;
        let mid_pt = DVec2::new(cx + r * mid_ang.cos(), cy + r * mid_ang.sin());
        is_kept[i] = !point_in_rect(mid_pt);
    }

    // ---- 3. Group consecutive KEPT arcs into runs ----
    // Note: zero-length arcs (corner intersections) are NOT promoted �?they
    // separate distinct boundary components around the box.
    let is_kept_run = &is_kept;
    struct Run { start: usize, end: usize }
    let mut runs: Vec<Run> = Vec::new();

    // Handle wrap-around: if both first and last arcs are KEPT, they are one run.
    let first_kept = is_kept_run[0];
    let last_kept = is_kept_run[m - 1];
    let mut merged_wrap = false;

    if first_kept && last_kept {
        // Merge wrap: find the first SKIPPED arc, runs start after it.
        let mut split = 0;
        while split < m && is_kept_run[split] { split += 1; }
        if split < m {
            // Collect runs starting from `split`
            let mut i = split;
            while i < split + m {
                let ii = i % m;
                if is_kept_run[ii] {
                    let run_start = ii;
                    while i < split + m && is_kept_run[i % m] { i += 1; }
                    let run_end = (i - 1) % m;
                    runs.push(Run { start: run_start, end: run_end });
                    merged_wrap = true;
                } else {
                    i += 1;
                }
            }
        }
    }

    if !merged_wrap {
        let mut i = 0;
        while i < m {
            if is_kept_run[i] {
                let start = i;
                while i < m && is_kept_run[i] { i += 1; }
                runs.push(Run { start, end: i - 1 });
            } else {
                i += 1;
            }
        }
    }

    if runs.is_empty() { return None; }

    // ---- 4. For each run, build the boundary polygon ----
    // Map each rect edge to its two intersection indices (0..m).
    // Edge order CW: bottom(0), right(1), top(2), left(3) �?but CW order is
    // bottom, left, top, right. We'll build by edge index and traverse manually.
    let mut edge_idxs: [Vec<usize>; 4] = [
        Vec::new(), Vec::new(), Vec::new(), Vec::new()
    ];
    for (idx, inter) in ints.iter().enumerate() {
        if inter.edge < 4 {
            edge_idxs[inter.edge].push(idx);
        }
    }

    // CW edge traversal order for the rect perimeter.
    // Starting from a corner and going CW: bottom閳姰ight閳姲op閳姡eft is CCW!
    // CW is: bottom閳姡eft閳姲op閳姰ight (right-hand rule around +Z).
    // Edge 0 (bottom) in CW direction: right閳姡eft (decreasing x)
    // Edge 3 (left) in CW direction: bottom閳姲op (increasing y)
    // Edge 2 (top) in CW direction: left閳姰ight (increasing x)
    // Edge 1 (right) in CW direction: top閳妼ottom (decreasing y)
    let cw_edge_order = [0usize, 3, 2, 1];

    // Helper: get the two intersection indices on an edge, ordered by CW traversal
    // (first = encountered first when walking CW along the rect).
    let edge_cw_order = |edge: usize| -> Option<(usize, usize)> {
        let e = &edge_idxs[edge];
        if e.len() < 2 { return None; }
        let (a, b) = (e[0], e[1]);
        // Edge t parameter increases along the edge direction (as in rect_edges).
        // CW direction is OPPOSITE to edge direction for ALL edges:
        // - edge 0 (bottom, A閳�?left閳姰ight): CW = right閳姡eft = larger t first
        // - edge 1 (right, B閳�?bottom閳姲op): CW = top閳妼ottom = larger t first
        // - edge 2 (top,    C閳�?right閳姡eft): CW = left閳姰ight = larger t first
        // - edge 3 (left,   D閳�?top閳妼ottom): CW = bottom閳姲op = larger t first
        if ints[a].t > ints[b].t { Some((a, b)) } else { Some((b, a)) }
    };

    // Build polygons: one per run
    for run in &runs {
        let mut pts: Vec<DVec2> = Vec::new();

        // Walk CCW from run.start through the consecutive KEPT arcs to run.end.
        // Zero-length SKIPPED arcs (promoted to KEPT in is_kept_run) are
        // transparent �?skip their points but continue through them.
        let mut idx = run.start;
        loop {
            let j = (idx + 1) % m;
            if is_kept[idx] {
                let v1 = ints[idx].pt - DVec2::new(cx, cy);
                let v2 = ints[j].pt - DVec2::new(cx, cy);
                let a1 = f64::atan2(v1.y, v1.x);
                let a2 = f64::atan2(v2.y, v2.x);
                let da_ccw = (a2 - a1).rem_euclid(tau);
                for k in 0..=n_arc {
                    let frac = k as f64 / n_arc as f64;
                    let ang = a1 + da_ccw * frac;
                    let (s, c) = ang.sin_cos();
                    pts.push(DVec2::new(cx + r * c, cy + r * s));
                }
            }
            if idx == run.end || !is_kept_run[j] { break; }
            idx = j;
        }

        // Chord return path: walk CW along rect perimeter from the CCW arc
        // endpoint back to the arc startpoint.  The last CCW arc goes from
        // ints[run.end] to ints[(run.end+1)%m]; the polygon ends near the latter.
        let arc_end_idx = (run.end + 1) % m;
        let end_edge = ints[arc_end_idx].edge;
        let start_edge = ints[run.start].edge;

        // Find which edges to traverse in CW order from end to start.
        let end_pos = cw_edge_order.iter().position(|e| *e == end_edge).unwrap_or(0);
        let _start_pos = cw_edge_order.iter().position(|e| *e == start_edge).unwrap_or(0);

        // Traverse edges from end_edge CW until we've processed start_edge.
        let mut cur_pos = end_pos;
        let mut first_edge = true;
        loop {
            let edge = cw_edge_order[cur_pos];
            if let Some((cw_first, cw_second)) = edge_cw_order(edge) {
                if first_edge {
                    // The polygon ends near ints[arc_end_idx].pt on this edge.
                    // CW direction on this edge: cw_first �?cw_second.
                    // If arc_end_idx == cw_first: we're at the CW start, add cw_second.
                    // If arc_end_idx == cw_second: we're at the CW end, no points to add.
                    if arc_end_idx == cw_first {
                        let last_pt = pts.last().copied();
                        if last_pt.map_or(true, |lp| (lp - ints[cw_second].pt).length_squared() > tol * tol) {
                            pts.push(ints[cw_second].pt);
                        }
                        if cw_second == run.start { break; }
                    }
                    first_edge = false;
                } else {
                    // Add points in CW order, stopping at run.start on the
                    // start edge to avoid over-traversing into the next run's
                    // boundary segment (which would double-count wall area).
                    let last_pt = pts.last().copied();
                    if last_pt.map_or(true, |lp| (lp - ints[cw_first].pt).length_squared() > tol * tol) {
                        pts.push(ints[cw_first].pt);
                    }
                    if edge == start_edge && cw_first == run.start { break; }
                    let last_pt2 = pts.last().copied();
                    if last_pt2.map_or(true, |lp| (lp - ints[cw_second].pt).length_squared() > tol * tol) {
                        pts.push(ints[cw_second].pt);
                    }
                    if edge == start_edge && cw_second == run.start { break; }
                }
            }

            // Check if we've reached start_edge
            if edge == start_edge {
                break;
            }

            // Move to next edge in CW order
            cur_pos = (cur_pos + 1) % 4;
        }

        // Create boundary polygon
        if pts.len() >= 3 {
            // ---- 5. Build Z-slice wall and cap faces for this component ----
            for i in 0..n_slices {
                let z0 = z_lo + dz * i as f64;
                let z1 = z_lo + dz * (i + 1) as f64;
                let n = pts.len();

                let mut idx = Vec::with_capacity(2 * n);
                for p in &pts { idx.push(add_v(to_world(p.x, p.y, z0))); }
                for p in &pts { idx.push(add_v(to_world(p.x, p.y, z1))); }

                let mut tris = Vec::with_capacity(n * 2);
                for j in 0..n {
                    let k = (j + 1) % n;
                    tris.push([idx[j], idx[k], idx[n + k]]);
                    tris.push([idx[j], idx[n + k], idx[n + j]]);
                }

                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: tris,
                    sample_point: None, mesh_dirty: false,
                surface_idx: None,
                });
            }

            // Bottom cap
            let poly_lo: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z_lo)).collect();
            let tris_lo = crate::triangulate::triangulate_polygon(&poly_lo, -DVec3::Z);
            if !tris_lo.is_empty() {
                let mut remapped = Vec::with_capacity(tris_lo.len());
                let local: Vec<usize> = poly_lo.iter().map(|p| add_v(*p)).collect();
                for t in &tris_lo { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: remapped,
                    sample_point: None, mesh_dirty: false,
                surface_idx: None,
                });
            }

            // Top cap
            let poly_hi: Vec<DVec3> = pts.iter().map(|p| to_world(p.x, p.y, z_hi)).collect();
            let tris_hi = crate::triangulate::triangulate_polygon(&poly_hi, DVec3::Z);
            if !tris_hi.is_empty() {
                let mut remapped = Vec::with_capacity(tris_hi.len());
                let local: Vec<usize> = poly_hi.iter().map(|p| add_v(*p)).collect();
                for t in &tris_hi { remapped.push([local[t[0]], local[t[1]], local[t[2]]]); }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: remapped,
                    sample_point: None, mesh_dirty: false,
                surface_idx: None,
                });
            }
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

// 閳光偓閳光偓 Cylinder-Box Union Tessellation 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光�
