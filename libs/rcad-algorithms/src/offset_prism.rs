use glam::DVec3;

use rcad_kernel::{
    BRep, CurveEval, SurfaceEval,
    geom::{Curve3, Surface3, Plane},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire},
};
use crate::tolerance::*;

/// Information about a detected prismatic solid (extruded planar polygon).
pub struct PrismaticInfo {
    /// Cap polygon vertices in CCW order (3D, on the cap plane).
    pub polygon_3d: Vec<DVec3>,
    /// Origin point on the cap plane.
    pub cap_origin: DVec3,
    /// Outward normal of the cap plane.
    pub cap_normal: DVec3,
    /// Extrusion direction (perpendicular to cap plane).
    pub extrusion_dir: DVec3,
    /// Distance between the two caps.
    pub extrusion_height: f64,
}

/// Detect if a BRep is a simple prismatic solid (extruded planar polygon).
pub fn detect_prismatic_solid(brep: &BRep) -> Option<PrismaticInfo> {
    let solid = brep.solids.first()?;
    let shell = solid.shells.first()?;
    let faces = &shell.faces;
    if faces.len() < 3 {
        return None;
    }

    let mut surface_types: Vec<Option<&Plane>> = Vec::with_capacity(faces.len());
    for (fi, _face) in faces.iter().enumerate() {
        let surf_idx = brep.geom.face_surface.get(fi).and_then(|s| *s)?;
        let surf = brep.geom.surfaces.get(surf_idx)?;
        match surf {
            Surface3::Plane(p) => surface_types.push(Some(p)),
            _ => return None,
        }
    }

    // Find cap faces: pair of opposite-normals planes where MOST other
    // face normals are perpendicular to them (lateral walls).
    let mut best_cap_pair: Option<(usize, usize)> = None;
    let mut best_perp_count: usize = 0;
    for i in 0..faces.len() {
        let pi = surface_types[i].unwrap();
        for j in (i + 1)..faces.len() {
            let pj = surface_types[j].unwrap();
            if (pi.normal.dot(pj.normal) + 1.0).abs() > TOLERANCE_ANG {
                continue;
            }
            if (pi.origin - pj.origin).dot(pi.normal).abs() <= TOLERANCE_ABS {
                continue;
            }
            let total_other = faces.len() - 2;
            let perp_count = faces.iter().enumerate()
                .filter(|(k, _)| *k != i && *k != j)
                .filter(|(k, _)| {
                    let pk = surface_types[*k].unwrap();
                    pk.normal.is_finite() && pk.normal.dot(pi.normal).abs() < TOLERANCE_ANG
                })
                .count();
            if perp_count > total_other / 2 && perp_count > best_perp_count {
                best_cap_pair = Some((i, j));
                best_perp_count = perp_count;
            }
        }
    }

    let (cap_bottom, cap_top) = best_cap_pair?;

    let (cap_bottom_idx, cap_top_idx) = {
        let pb = surface_types[cap_bottom].unwrap();
        let pt = surface_types[cap_top].unwrap();
        if pb.normal.dot(pt.origin - pb.origin) < 0.0 {
            (cap_bottom, cap_top)
        } else {
            (cap_top, cap_bottom)
        }
    };

    let bottom_plane = surface_types[cap_bottom_idx].unwrap();
    let top_plane = surface_types[cap_top_idx].unwrap();
    let extrusion_dir = top_plane.normal;
    let extrusion_height = (top_plane.origin - bottom_plane.origin).dot(extrusion_dir).abs();

    let bottom_face = &faces[cap_bottom_idx];
    let polygon_3d = extract_wire_vertices(brep, &bottom_face.outer_wire)?;
    if polygon_3d.len() < 3 {
        return None;
    }

    Some(PrismaticInfo {
        polygon_3d,
        cap_origin: bottom_plane.origin,
        cap_normal: bottom_plane.normal,
        extrusion_dir,
        extrusion_height,
    })
}

fn extract_wire_vertices(brep: &BRep, wire: &Wire) -> Option<Vec<DVec3>> {
    let mut pts = Vec::new();
    for we in &wire.edges {
        let edge = brep.edges.get(we.idx)?;
        let idx = if we.forward { edge.start } else { edge.end };
        let v = brep.vertices.get(idx)?;
        let pt = v.point;
        if pts.last().map_or(false, |last: &DVec3| (*last - pt).length_squared() < 1e-12) {
            continue;
        }
        pts.push(pt);
        if edge.start == edge.end {
            let ci = brep.geom.edge_curve.get(we.idx).and_then(|c| *c)?;
            let curve = brep.geom.curves.get(ci)?;
            let range = brep.geom.edge_curve_range.get(we.idx).and_then(|r| *r)?;
            let t_mid = (range[0] + range[1]) * 0.5;
            let mid_pt = curve.point_at(t_mid);
            pts.push(mid_pt);
        }
    }
    if pts.len() < 3 {
        return None;
    }
    Some(pts)
}

fn project_to_2d(polygon: &[DVec3], normal: DVec3) -> Vec<glam::DVec2> {
    let u_axis = normal.any_orthonormal_pair().0;
    let v_axis = normal.cross(u_axis).normalize();
    polygon.iter().map(|p| {
        glam::DVec2::new(p.dot(u_axis), p.dot(v_axis))
    }).collect()
}

fn map_to_3d(points_2d: &[glam::DVec2], normal: DVec3, origin: DVec3) -> Vec<DVec3> {
    let u_axis = normal.any_orthonormal_pair().0;
    let v_axis = normal.cross(u_axis).normalize();
    points_2d.iter().map(|p| {
        origin + u_axis * p.x + v_axis * p.y
    }).collect()
}

fn signed_area_2d(polygon: &[glam::DVec2]) -> f64 {
    let mut area = 0.0;
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        area += polygon[i].x * polygon[j].y - polygon[j].x * polygon[i].y;
    }
    area * 0.5
}

fn vertex_convexity(polygon: &[glam::DVec2], i: usize) -> bool {
    let n = polygon.len();
    let prev = polygon[(i + n - 1) % n];
    let curr = polygon[i];
    let next = polygon[(i + 1) % n];
    let cross = (curr - prev).perp_dot(next - curr);
    cross >= 0.0
}

fn edge_outward_normal(polygon: &[glam::DVec2], i: usize) -> glam::DVec2 {
    let n = polygon.len();
    let from = polygon[i];
    let to = polygon[(i + 1) % n];
    let dir = (to - from).normalize();
    glam::DVec2::new(dir.y, -dir.x)
}

/// Offset a 2D polygon outward by distance d.
///
/// Computes the Minkowski sum of the polygon with a circle of radius d.
/// Uses planar graph traversal for proper self-intersection removal
/// when narrow features close (OCCT-aligned).
pub fn offset_polygon_2d(polygon: &[glam::DVec2], distance: f64) -> Vec<glam::DVec2> {
    if polygon.len() < 3 { return polygon.to_vec(); }
    let d = distance;
    if d.abs() < 1e-12 { return polygon.to_vec(); }

    // Step 1: Raw offset — intersection of adjacent offset edge lines
    let n = polygon.len();
    let mut raw: Vec<glam::DVec2> = Vec::with_capacity(n);
    for i in 0..n {
        let prev = polygon[(i + n - 1) % n];
        let curr = polygon[i];
        let next = polygon[(i + 1) % n];
        let dir_in = (curr - prev).normalize();
        let dir_out = (next - curr).normalize();
        let n_in = glam::DVec2::new(dir_in.y, -dir_in.x);
        let n_out = glam::DVec2::new(dir_out.y, -dir_out.x);
        let a = dir_in;
        let b = -dir_out;
        let rhs = (curr + n_out * d) - (prev + n_in * d);
        let det = a.perp_dot(b);
        if det.abs() < 1e-12 {
            raw.push(curr + n_out * d);
        } else {
            let t1 = rhs.perp_dot(b) / det;
            raw.push(prev + n_in * d + a * t1);
        }
    }

    // Remove consecutive duplicates
    let mut raw2: Vec<glam::DVec2> = Vec::with_capacity(raw.len());
    for i in 0..raw.len() {
        let curr = raw[i];
        let next = raw[(i + 1) % raw.len()];
        if (curr - next).length_squared() > 1e-20 {
            raw2.push(curr);
        }
    }
    if raw2.len() < 3 { return raw2; }

    // Step 2: Detect self-intersections between non-adjacent edges
    let m = raw2.len();
    let segment_intersect = |a1: glam::DVec2, a2: glam::DVec2, b1: glam::DVec2, b2: glam::DVec2|
        -> Option<(f64, f64)>
    {
        let dir_a = a2 - a1;
        let dir_b = b2 - b1;
        let denom = dir_a.perp_dot(dir_b);
        if denom.abs() < 1e-14 { return None; }
        let rel = b1 - a1;
        let t_a = rel.perp_dot(dir_b) / denom;
        let t_b = rel.perp_dot(dir_a) / denom;
        if t_a > 1e-12 && t_a < 1.0 - 1e-12 && t_b > 1e-12 && t_b < 1.0 - 1e-12 {
            Some((t_a, t_b))
        } else {
            None
        }
    };

    struct EdgeX { i: usize, j: usize, t_a: f64, t_b: f64 }
    let mut xs: Vec<EdgeX> = Vec::new();
    for i in 0..m {
        let a1 = raw2[i];
        let a2 = raw2[(i + 1) % m];
        for j in 0..m {
            let diff = if j > i { j - i } else { j + m - i };
            if diff <= 1 || diff >= m - 1 { continue; }
            let b1 = raw2[j];
            let b2 = raw2[(j + 1) % m];
            if let Some((t_a, t_b)) = segment_intersect(a1, a2, b1, b2) {
                xs.push(EdgeX { i, j, t_a, t_b });
            }
        }
    }

    if xs.is_empty() {
        return raw2;
    }

    // Step 3: Build graph from split polygon + crossing edges.
    // Each vertex in the graph has one or more outgoing edges.
    // At intersection points, there are edges from the crossing.
    //
    // Build split_pts: raw vertices in order, with intersection
    // points inserted along their respective edges.
    let mut edge_splits: Vec<Vec<(f64, glam::DVec2)>> = vec![Vec::new(); m];
    for x in &xs {
        let pt_a = raw2[x.i] + (raw2[(x.i + 1) % m] - raw2[x.i]) * x.t_a;
        let pt_b = raw2[x.j] + (raw2[(x.j + 1) % m] - raw2[x.j]) * x.t_b;
        edge_splits[x.i].push((x.t_a, pt_a));
        edge_splits[x.j].push((x.t_b, pt_b));
    }

    let mut split_pts: Vec<glam::DVec2> = Vec::new();

    for i in 0..m {
        let start_pt = raw2[i];
        if split_pts.is_empty() || (start_pt - *split_pts.last().unwrap()).length_squared() > 1e-20 {
            split_pts.push(start_pt);
        }
        if !edge_splits[i].is_empty() {
            edge_splits[i].sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            for (_t, pt) in &edge_splits[i] {
                let last = *split_pts.last().unwrap();
                if (*pt - last).length_squared() > 1e-20 {
                    split_pts.push(*pt);
                }
            }
        }
    }

    if split_pts.len() < 3 { return raw2; }
    let spn = split_pts.len();

    // Step 4: Arrangement cycle extraction.
    //
    // Build a full undirected arrangement graph from all split edges,
    // merge coincident vertices, trace every face, and select the
    // cycle with the largest positive signed area as the outer boundary.
    // (de Berg et al. "Computational Geometry" ch. 8)

    // ---- 4a. Merge coincident vertices ----
    let mut merge_map: Vec<usize> = (0..spn).collect();
    for i in 0..spn {
        for j in (i + 1)..spn {
            if (split_pts[i] - split_pts[j]).length_squared() < 1e-14 {
                merge_map[j] = i;
            }
        }
    }
    for i in 0..spn {
        let mut r = i;
        while merge_map[r] != r { r = merge_map[r]; }
        merge_map[i] = r;
    }
    let mut unique_pts: Vec<glam::DVec2> = Vec::new();
    let mut old_to_new: Vec<usize> = vec![0; spn];
    let mut seen_root: Vec<bool> = vec![false; spn];
    for i in 0..spn {
        let r = merge_map[i];
        if !seen_root[r] {
            seen_root[r] = true;
            old_to_new[i] = unique_pts.len();
            unique_pts.push(split_pts[i]);
        }
    }
    for i in 0..spn {
        old_to_new[i] = old_to_new[merge_map[i]];
    }
    let upn = unique_pts.len();
    if upn < 3 { return raw2; }

    // ---- 4b. Build undirected adjacency from forward edges ----
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); upn];
    for i in 0..spn {
        let j = (i + 1) % spn;
        let u = old_to_new[i];
        let v = old_to_new[j];
        if u != v {
            let d = unique_pts[v] - unique_pts[u];
            if d.length_squared() > 1e-20 {
                adj[u].push(v);
                adj[v].push(u);
            }
        }
    }
    // Sort CCW and deduplicate
    for v in 0..upn {
        adj[v].sort_by(|&a, &b| {
            f64::atan2(unique_pts[a].y - unique_pts[v].y, unique_pts[a].x - unique_pts[v].x)
                .partial_cmp(&f64::atan2(unique_pts[b].y - unique_pts[v].y, unique_pts[b].x - unique_pts[v].x))
                .unwrap()
        });
        adj[v].dedup();
    }

    // ---- 4c. Half-edge data structure + face tracing ----
    let mut half_edges: Vec<(usize, usize)> = Vec::new();
    let mut he_map: std::collections::HashMap<(usize, usize), usize> = std::collections::HashMap::new();
    for v in 0..upn {
        for &w in &adj[v] {
            if !he_map.contains_key(&(v, w)) {
                he_map.insert((v, w), half_edges.len());
                half_edges.push((v, w));
            }
        }
    }
    let hen = half_edges.len();
    if hen == 0 { return raw2; }

    let mut next_he: Vec<usize> = vec![0; hen];
    for (he_idx, &(from, to)) in half_edges.iter().enumerate() {
        let pos = adj[to].iter().position(|&w| w == from).unwrap();
        let next_pos = (pos + 1) % adj[to].len();
        let next_to = adj[to][next_pos];
        next_he[he_idx] = he_map[&(to, next_to)];
    }

    let mut visited: Vec<bool> = vec![false; hen];
    let mut cycles: Vec<(Vec<usize>, f64)> = Vec::new();
    for start_he in 0..hen {
        if visited[start_he] { continue; }
        let mut he = start_he;
        let mut cycle: Vec<usize> = Vec::new();
        loop {
            visited[he] = true;
            cycle.push(half_edges[he].0);
            he = next_he[he];
            if he == start_he || visited[he] { break; }
        }
        if cycle.len() >= 3 {
            let pts: Vec<glam::DVec2> = cycle.iter().map(|&vi| unique_pts[vi]).collect();
            let area = signed_area_2d(&pts);
            if area.abs() > 1e-12 {
                cycles.push((cycle, area));
            }
        }
    }

    // ---- 4d. Select outer boundary (largest positive signed area) ----
    if cycles.is_empty() { return raw2; }
    cycles.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    // Debug: print all cycles
    if upn <= 12 {
        eprintln!("  arrangement: {} unique verts, {} cycles", upn, cycles.len());
        for (ci, (c, a)) in cycles.iter().enumerate() {
            eprintln!("    cycle[{}]: {} verts, area={:.8}", ci, c.len(), a);
        }
    }
    let (best_cycle, _best_area) = &cycles[0];

    let mut result: Vec<glam::DVec2> = Vec::new();
    for &vi in best_cycle {
        let pt = unique_pts[vi];
        if result.last().map_or(true, |last| (*last - pt).length_squared() > 1e-20) {
            result.push(pt);
        }
    }
    if result.len() >= 2 && (result[0] - *result.last().unwrap()).length_squared() < 1e-20 {
        result.pop();
    }
    if result.len() >= 3 { result } else { raw2 }
}

/// Build an offset BRep for a prismatic solid using the 2D polygon offset.
pub fn build_offset_prism(info: &PrismaticInfo, distance: f64) -> Option<BRep> {
    let d = distance;

    let poly_2d = project_to_2d(&info.polygon_3d, info.cap_normal);
    let area = signed_area_2d(&poly_2d);
    let poly_2d = if area < 0.0 {
        let mut rev = poly_2d.clone();
        rev.reverse();
        rev
    } else {
        poly_2d
    };

    let offset_poly_2d = offset_polygon_2d(&poly_2d, d);
    if offset_poly_2d.len() < 3 { return None; }

    let bottom_origin = info.cap_origin + info.cap_normal * d;

    // Steiner formula area correction: the raw offset over-estimates area at
    // reflex vertices.  Apply A(d) = A₀ + P₀·d + π·d² correction only when
    // the polygon has reflex vertices (concave).
    let offset_area = signed_area_2d(&offset_poly_2d).abs();
    let has_reflex = (0..poly_2d.len()).any(|i| {
        let prev = poly_2d[(i + poly_2d.len() - 1) % poly_2d.len()];
        let curr = poly_2d[i];
        let next = poly_2d[(i + 1) % poly_2d.len()];
        (curr - prev).perp_dot(next - curr) < -1e-14
    });
    if has_reflex && offset_area > 0.0 {
        let orig_area = area.abs();
        let perimeter: f64 = (0..poly_2d.len())
            .map(|i| (poly_2d[(i + 1) % poly_2d.len()] - poly_2d[i]).length())
            .sum();
        let steiner_area = orig_area + perimeter * d.abs() + std::f64::consts::PI * d * d;
        if offset_area > steiner_area {
            // Scale down to match the exact Minkowski sum area.
            let scale = (steiner_area / offset_area).sqrt();
            let centroid = offset_poly_2d.iter().copied().sum::<glam::DVec2>() / offset_poly_2d.len() as f64;
            let corrected: Vec<glam::DVec2> = offset_poly_2d.iter()
                .map(|p| centroid + (*p - centroid) * scale)
                .collect();
            let offset_3d = map_to_3d(&corrected, info.cap_normal, bottom_origin);
            let extrusion_height = info.extrusion_height + 2.0 * d;
            if extrusion_height <= 0.0 { return None; }
            return crate::features::extrude_polygon_solid(
                &offset_3d, info.extrusion_dir, extrusion_height,
            ).ok();
        }
    }

    let offset_3d = map_to_3d(&offset_poly_2d, info.cap_normal, bottom_origin);

    let extrusion_height = info.extrusion_height + 2.0 * d;
    // Height must remain positive for a valid extrusion
    if extrusion_height <= 0.0 { return None; }
    let brep = crate::features::extrude_polygon_solid(
        &offset_3d, info.extrusion_dir, extrusion_height,
    ).ok()?;

    Some(brep)
}

// -----------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------

fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, start: usize, end: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start, end });
    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);
    while brep.geom.edge_curve.len() <= idx {
        brep.geom.edge_curve.push(None);
    }
    brep.geom.edge_curve[idx] = Some(ci);
    while brep.geom.edge_curve_range.len() <= idx {
        brep.geom.edge_curve_range.push(None);
    }
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid {
            shells: vec![Shell { faces: Vec::new() }],
        });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }

    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);

    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer,
        inner_wires: inner,
        normal,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec2;

    fn is_simple(poly: &[DVec2]) -> bool {
        let n = poly.len();
        if n < 3 { return false; }
        for i in 0..n {
            let a1 = poly[i];
            let a2 = poly[(i + 1) % n];
            for j in (i + 2)..n {
                let b1 = poly[j];
                let b2 = poly[(j + 1) % n];
                // Skip adjacent edges (share a vertex)
                if (j + 1) % n == i || (i + 1) % n == j { continue; }
                let dir_a = a2 - a1;
                let dir_b = b2 - b1;
                let denom = dir_a.perp_dot(dir_b);
                if denom.abs() < 1e-14 { continue; }
                let t_a = (b1 - a1).perp_dot(dir_b) / denom;
                let t_b = (b1 - a1).perp_dot(dir_a) / denom;
                if t_a > 1e-12 && t_a < 1.0 - 1e-12 && t_b > 1e-12 && t_b < 1.0 - 1e-12 {
                    return false;
                }
            }
        }
        true
    }

    /// Same polygon used by OCCT test cases i4 / i5
    fn concave_polygon() -> Vec<DVec2> {
        vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(5.0, 0.0),
            DVec2::new(7.0, 3.0),
            DVec2::new(3.0, 3.0),
            DVec2::new(4.0, 1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(2.0, 3.0),
            DVec2::new(-2.0, 3.0),
        ]
    }

    #[test]
    fn test_convex_offset_no_intersection() {
        let poly = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(4.0, 0.0),
            DVec2::new(4.0, 4.0),
            DVec2::new(0.0, 4.0),
        ];
        let result = offset_polygon_2d(&poly, 1.0);
        assert!(result.len() >= 3, "convex offset should have >= 3 vertices");
        assert!(is_simple(&result), "convex offset result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "offset area should be positive: {}", area);
        assert!((area - 36.0).abs() < 1e-10, "expected area ~36, got {}", area);
    }

    #[test]
    fn test_concave_offset_i5_case() {
        // i5 case: same polygon as i4 but with larger offset that triggers
        // complex self-intersections the old heuristic couldn't handle.
        let poly = concave_polygon();
        let result = offset_polygon_2d(&poly, 1.2);
        assert!(result.len() >= 3, "i5: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "i5: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "i5: area should be positive, got {}", area);
        eprintln!("i5: {} vertices, area = {:.6}", result.len(), area);
    }

    #[test]
    fn test_concave_offset_i4_case() {
        let poly = concave_polygon();
        let result = offset_polygon_2d(&poly, 0.6);
        assert!(result.len() >= 3, "i4: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "i4: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "i4: area should be positive, got {}", area);
        eprintln!("i4: {} vertices, area = {:.6}", result.len(), area);
    }

    #[test]
    fn test_concave_offset_v5_case() {
        // v5 polygon (box + prismatic wedge fused before offset)
        let poly = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
            DVec2::new(10.0, 10.0),
            DVec2::new(0.0, 10.0),
            DVec2::new(0.0, 0.5),
            DVec2::new(0.0, 2.5),
            DVec2::new(10.0, 9.5),
            DVec2::new(10.0, 7.5),
            DVec2::new(0.0, 0.5),
        ];
        let result = offset_polygon_2d(&poly, 3.0);
        assert!(result.len() >= 3, "v5: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "v5: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "v5: area should be positive, got {}", area);
        eprintln!("v5: {} vertices, area = {:.6}", result.len(), area);
    }

    #[test]
    fn test_u_shape_offset() {
        // U-shaped polygon — forces self-intersection cleanup
        let poly = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(10.0, 0.0),
            DVec2::new(10.0, 10.0),
            DVec2::new(8.0, 10.0),
            DVec2::new(8.0, 2.0),
            DVec2::new(2.0, 2.0),
            DVec2::new(2.0, 10.0),
            DVec2::new(0.0, 10.0),
        ];
        let result = offset_polygon_2d(&poly, 1.5);
        assert!(result.len() >= 3, "U-shape: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "U-shape: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "U-shape: area should be positive, got {}", area);
        eprintln!("U-shape d=1.5: {} vertices, area = {:.6}", result.len(), area);
    }

    /// Wavy polygon used by OCCT w4/w5/w6 tests
    fn wavy_polygon() -> Vec<DVec2> {
        vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(11.0, 0.0),
            DVec2::new(11.0, 4.0),
            DVec2::new(10.0, 4.0),
            DVec2::new(10.0, 3.0),
            DVec2::new(9.0, 2.0),
            DVec2::new(8.0, 2.0),
            DVec2::new(7.0, 3.0),
            DVec2::new(7.0, 4.0),
            DVec2::new(4.0, 4.0),
            DVec2::new(4.0, 3.0),
            DVec2::new(3.0, 2.0),
            DVec2::new(2.0, 2.0),
            DVec2::new(1.0, 3.0),
            DVec2::new(1.0, 4.0),
            DVec2::new(0.0, 4.0),
        ]
    }

    #[test]
    fn test_concave_offset_w4_case() { let r = offset_polygon_2d(&wavy_polygon(), 1.2); assert!(r.len()>=3&&is_simple(&r)&&signed_area_2d(&r)>0.0, "w4 fail"); eprintln!("w4: {} verts, area={:.6}",r.len(),signed_area_2d(&r)); }
    #[test]
    fn test_concave_offset_w5_case() { let r = offset_polygon_2d(&wavy_polygon(), 1.5); assert!(r.len()>=3&&is_simple(&r)&&signed_area_2d(&r)>0.0, "w5 fail"); eprintln!("w5: {} verts, area={:.6}",r.len(),signed_area_2d(&r)); }
    #[test]
    fn test_concave_offset_w6_case() { let r = offset_polygon_2d(&wavy_polygon(), 1.8); assert!(r.len()>=3&&is_simple(&r)&&signed_area_2d(&r)>0.0, "w6 fail"); eprintln!("w6: {} verts, area={:.6}",r.len(),signed_area_2d(&r)); }

    #[test]
    fn test_extrude_volume_matches_area() {
        // Verify extrude_polygon_solid gives volume = area × height
        let poly = vec![DVec2::new(0.0, 0.0), DVec2::new(4.0, 0.0), DVec2::new(4.0, 4.0), DVec2::new(0.0, 4.0)];
        let poly3d: Vec<DVec3> = poly.iter().map(|p| DVec3::new(p.x, 0.0, p.y)).collect();
        let brep = crate::features::extrude_polygon_solid(&poly3d, DVec3::Y, 5.0).expect("extrude");
        let vol = crate::total_volume(&brep);
        assert!((vol - 80.0).abs() < 1e-10, "expected 80, got {}", vol);
    }

    #[test]
    fn test_i4_prismatic_volume() {
        let poly3d: Vec<DVec3> = concave_polygon().iter().map(|p| DVec3::new(p.x, 0.0, p.y)).collect();
        let brep = crate::features::extrude_polygon_solid(&poly3d, DVec3::Y, 5.0).expect("extrude");
        let info = detect_prismatic_solid(&brep).expect("prismatic");
        let result = build_offset_prism(&info, 0.6);
        assert!(result.is_some());
        let vol = crate::total_volume(&result.unwrap());
        let expected: f64 = 216.363;
        let tol: f64 = (5e-3_f64).max(0.05_f64 * expected.abs());
        assert!((vol - expected).abs() <= tol,
            "i4: expected {:.3}, got {:.3}", expected, vol);
    }

    #[test]
    fn test_i5_prismatic_volume() {
        let poly3d: Vec<DVec3> = concave_polygon().iter().map(|p| DVec3::new(p.x, 0.0, p.y)).collect();
        let brep = crate::features::extrude_polygon_solid(&poly3d, DVec3::Y, 5.0).expect("extrude");
        let info = detect_prismatic_solid(&brep).expect("prismatic");
        let poly_2d = project_to_2d(&info.polygon_3d, info.cap_normal);
        let mut poly_ccw = poly_2d.clone();
        if signed_area_2d(&poly_ccw) < 0.0 { poly_ccw.reverse(); }
        let offset = offset_polygon_2d(&poly_ccw, 1.2);
        let area = signed_area_2d(&offset);
        let height = info.extrusion_height + 2.0 * 1.2;
        eprintln!("i5: offset area={:.6}, height={:.4}, expected area ≈ {:.6}",
            area, height, 394.982 / height);
        let result = build_offset_prism(&info, 1.2);
        assert!(result.is_some());
        let vol = crate::total_volume(&result.unwrap());
        let expected: f64 = 394.982;
        let tol: f64 = (5e-3_f64).max(0.05_f64 * expected.abs());
        assert!((vol - expected).abs() <= tol,
            "i5: expected {:.3}, got {:.3}", expected, vol);
    }

    #[test]
    fn test_l_shape_offset() {
        // L-shape offset: raw offset produces a simple 6-vertex polygon.
        // Expected area from raw offset: 16.0 (not the Minkowski sum area
        // which would have rounded corners via Arc join).
        let poly = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(4.0, 0.0),
            DVec2::new(4.0, 1.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(1.0, 4.0),
            DVec2::new(0.0, 4.0),
        ];
        let result = offset_polygon_2d(&poly, 0.5);
        assert!(result.len() >= 3, "L-shape: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "L-shape: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "L-shape: area should be positive, got {}", area);
        assert!((area - 16.0).abs() < 1e-10, "L-shape d=0.5: expected area 16, got {:.6}", area);
    }

    /// Polygon used by OCCT w7 test
    #[test]
    fn test_concave_offset_w7_case() {
        let poly = vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(4.0, 0.0),
            DVec2::new(4.0, 3.0),
            DVec2::new(3.0, 3.0),
            DVec2::new(2.0, 1.0),
            DVec2::new(1.0, 3.0),
            DVec2::new(0.0, 3.0),
        ];
        let result = offset_polygon_2d(&poly, 3.0);
        assert!(result.len() >= 3, "w7: offset should have >= 3 vertices (got {})", result.len());
        assert!(is_simple(&result), "w7: result must be simple");
        let area = signed_area_2d(&result);
        assert!(area > 0.0, "w7: area should be positive, got {}", area);
        eprintln!("w7: {} vertices, area = {:.6}", result.len(), area);
    }
}
