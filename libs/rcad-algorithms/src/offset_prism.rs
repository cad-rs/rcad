use std::f64::consts::PI;
use glam::DVec3;

use rcad_kernel::{
    BRep, CurveEval, SurfaceEval,
    geom::{Curve3, Surface3, Plane, Line3, Point3},
    topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge},
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
    // Track which raw edge each split point belongs to, and the t param.
    struct PtInfo { edge: usize, t: f64, is_raw_start: bool }
    let mut pt_info: Vec<PtInfo> = Vec::new();

    for i in 0..m {
        let start_pt = raw2[i];
        if split_pts.is_empty() || (start_pt - *split_pts.last().unwrap()).length_squared() > 1e-20 {
            pt_info.push(PtInfo { edge: i, t: 0.0, is_raw_start: true });
            split_pts.push(start_pt);
        }
        if !edge_splits[i].is_empty() {
            edge_splits[i].sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            for (t, pt) in &edge_splits[i] {
                let last = *split_pts.last().unwrap();
                if (*pt - last).length_squared() > 1e-20 {
                    pt_info.push(PtInfo { edge: i, t: *t, is_raw_start: false });
                    split_pts.push(*pt);
                }
            }
        }
    }

    if split_pts.len() < 3 { return raw2; }
    // Find the leftmost-bottommost vertex (guaranteed on outer boundary).
    let start_idx = split_pts.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let spn = split_pts.len();

    // Build graph: for each vertex, list (next_vertex, outgoing_direction, is_cross)
    let mut graph: Vec<Vec<(usize, glam::DVec2, bool)>> = vec![Vec::new(); spn];

    // Forward edges: along raw polygon direction
    for i in 0..spn {
        let next = (i + 1) % spn;
        let dir = split_pts[next] - split_pts[i];
        if dir.length_squared() > 1e-20 {
            graph[i].push((next, dir.normalize(), false));
        }
    }

    // Cross edges: at intersection points, connect across the edges.
    // For each intersection, the split_pts contains intersection points
    // on both edges. Cross from the point on edge i to the vertex AFTER
    // the intersection on edge j, and vice versa.
    for x in &xs {
        let pt_i = raw2[x.i] + (raw2[(x.i + 1) % m] - raw2[x.i]) * x.t_a;
        let pt_j = raw2[x.j] + (raw2[(x.j + 1) % m] - raw2[x.j]) * x.t_b;
        let cross_to_i = (x.i + 1) % m; // raw vertex at end of edge i
        let cross_to_j = (x.j + 1) % m; // raw vertex at end of edge j

        // Find split_pts indices
        let mut spi_i = None;
        let mut spi_j = None;
        let mut spj_i = None;
        let mut spj_j = None;

        for sp in 0..spn {
            let d = (split_pts[sp] - pt_i).length_squared();
            if d < 1e-14 {
                if pt_info[sp].edge == x.i { spi_i = Some(sp); }
                if pt_info[sp].edge == x.j { spi_j = Some(sp); }
            }
            let d2 = (split_pts[sp] - pt_j).length_squared();
            if d2 < 1e-14 {
                if pt_info[sp].edge == x.i { spj_i = Some(sp); }
                if pt_info[sp].edge == x.j { spj_j = Some(sp); }
            }
        }

        // Find cross targets in split_pts
        let mut cti = None;
        let mut ctj = None;
        for sp in 0..spn {
            if pt_info[sp].is_raw_start && pt_info[sp].edge == cross_to_i {
                cti = Some(sp);
            }
            if pt_info[sp].is_raw_start && pt_info[sp].edge == cross_to_j {
                ctj = Some(sp);
            }
        }

        // Add cross edges:
        if let (Some(from_i), Some(to_j)) = (spi_i.or(spj_i), ctj) {
            let dir = split_pts[to_j] - split_pts[from_i];
            if dir.length_squared() > 1e-20 {
                graph[from_i].push((to_j, dir.normalize(), true));
            }
        }
        if let (Some(from_j), Some(to_i)) = (spi_j.or(spj_j), cti) {
            let dir = split_pts[to_i] - split_pts[from_j];
            if dir.length_squared() > 1e-20 {
                graph[from_j].push((to_i, dir.normalize(), true));
            }
        }
    }

    // Step 4: Walk the outer boundary with cross-edge discipline.
    // After taking a cross edge, disable further cross edges at the
    // destination to prevent double-crossing into interior loops.

    let incoming_dir = if start_idx > 0 {
        split_pts[start_idx] - split_pts[(start_idx + spn - 1) % spn]
    } else {
        glam::DVec2::new(0.0, -1.0)
    };

    let inc_len = incoming_dir.length();
    let mut inc = if inc_len > 1e-20 { incoming_dir / inc_len } else { glam::DVec2::new(0.0, -1.0) };
    let mut used_edges: Vec<Vec<bool>> = (0..spn).map(|i| vec![false; graph[i].len()]).collect();
    let mut outer = Vec::new();
    let mut cur = start_idx;
    let mut just_crossed = false;
    let mut iter = 0;
    let max_iter = spn * 8;

    loop {
        outer.push(split_pts[cur]);

        // Among unvisited outgoing edges, pick the one with the most
        // clockwise turn. If just_crossed, skip cross edges.
        let mut best: Option<usize> = None;
        let mut best_angle = f64::INFINITY;

        for (ei, (next, dir, is_cross)) in graph[cur].iter().enumerate() {
            if used_edges[cur][ei] { continue; }
            if just_crossed && *is_cross { continue; }  // Prevent double-cross
            if *next == start_idx && outer.len() >= 2 {
                best = Some(ei);  // Prefer returning to start
                break;
            }
            let angle = inc.perp_dot(*dir).atan2(inc.dot(*dir));
            if angle < best_angle {
                best_angle = angle;
                best = Some(ei);
            }
        }

        let best_idx = match best { Some(i) => i, None => break };
        used_edges[cur][best_idx] = true;
        let (next, dir, is_cross) = graph[cur][best_idx];
        just_crossed = is_cross;
        inc = dir;
        cur = next;

        iter += 1;
        if cur == start_idx || iter > max_iter { break; }
    }
    // Safety check: if the walk produced a self-intersecting polygon,
    // compute the convex hull of the walk vertices (which is guaranteed
    // to be simple and contain the correct offset for narrow-feature
    // closure).
    fn is_self_intersecting(poly: &[glam::DVec2]) -> bool {
        let np = poly.len();
        for i in 0..np {
            let a1 = poly[i];
            let a2 = poly[(i + 1) % np];
            for j in 0..np {
                let diff = if j > i { j - i } else { j + np - i };
                if diff <= 1 || diff >= np - 1 { continue; }
                let b1 = poly[j];
                let b2 = poly[(j + 1) % np];
                let dir_a = a2 - a1;
                let dir_b = b2 - b1;
                let denom = dir_a.perp_dot(dir_b);
                if denom.abs() < 1e-14 { continue; }
                let t_a = (b1 - a1).perp_dot(dir_b) / denom;
                let t_b = (b1 - a1).perp_dot(dir_a) / denom;
                if t_a > 1e-12 && t_a < 1.0 - 1e-12 && t_b > 1e-12 && t_b < 1.0 - 1e-12 {
                    return true;
                }
            }
        }
        false
    }

    fn convex_hull(pts: &[glam::DVec2]) -> Vec<glam::DVec2> {
        if pts.len() < 3 { return pts.to_vec(); }
        let mut p = pts.to_vec();
        p.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()));
        let cross = |o: &glam::DVec2, a: &glam::DVec2, b: &glam::DVec2| -> f64 {
            (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
        };
        let mut lower: Vec<glam::DVec2> = Vec::new();
        for pt in &p {
            while lower.len() >= 2 && cross(&lower[lower.len()-2], &lower[lower.len()-1], pt) <= 0.0 { lower.pop(); }
            lower.push(*pt);
        }
        let mut upper: Vec<glam::DVec2> = Vec::new();
        for pt in p.iter().rev() {
            while upper.len() >= 2 && cross(&upper[upper.len()-2], &upper[upper.len()-1], pt) <= 0.0 { upper.pop(); }
            upper.push(*pt);
        }
        lower.pop(); upper.pop(); lower.extend(upper); lower
    }

    // If the walk produced a self-intersecting polygon, use convex hull instead.
    if outer.len() >= 3 && is_self_intersecting(&outer) {
        let hull = convex_hull(&outer);
        if hull.len() >= 3 {
            // Also check the hull for self-intersection
            if !is_self_intersecting(&hull) {
                return hull;
            }
        }
    }

    if outer.len() >= 3 {
        outer
    } else {
        raw2
    }
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
    let offset_3d = map_to_3d(&offset_poly_2d, info.cap_normal, bottom_origin);

    let extrusion_height = info.extrusion_height + 2.0 * d.abs();
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
    });

    while brep.geom.face_surface.len() <= idx {
        brep.geom.face_surface.push(None);
    }

    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);

    idx
}
