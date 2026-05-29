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
///
/// A prismatic solid has:
/// - All planar faces
/// - Two cap faces (parallel, opposite normals)
/// - N lateral wall faces (perpendicular to caps)
/// - All lateral faces form a closed chain
pub fn detect_prismatic_solid(brep: &BRep) -> Option<PrismaticInfo> {
    let solid = brep.solids.first()?;
    let shell = solid.shells.first()?;
    let faces = &shell.faces;
    if faces.len() < 3 {
        return None;
    }

    // Step 1: Classify faces by surface type
    let mut surface_types: Vec<Option<&Plane>> = Vec::with_capacity(faces.len());

    for (fi, _face) in faces.iter().enumerate() {
        let surf_idx = brep.geom.face_surface.get(fi).and_then(|s| *s)?;
        let surf = brep.geom.surfaces.get(surf_idx)?;
        match surf {
            Surface3::Plane(p) => surface_types.push(Some(p)),
            _ => return None, // All faces must be planar
        }
    }

    // Step 2: Find cap faces (parallel planes with opposite normals).
    // The CORRECT cap pair is the one whose normal is perpendicular to
    // most other faces (all lateral walls are perpendicular to the caps).
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

    // Determine which cap is "bottom" (normal points against extrusion direction)
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

    // Extrusion direction = from bottom to top
    let extrusion_dir = top_plane.normal;
    let extrusion_height = (top_plane.origin - bottom_plane.origin).dot(extrusion_dir).abs();

    // Step 3: Extract polygon vertices from the bottom cap's wire
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

/// Extract 3D vertices from a wire in order, skipping self-loop edges.
fn extract_wire_vertices(brep: &BRep, wire: &Wire) -> Option<Vec<DVec3>> {
    let mut pts = Vec::new();
    for we in &wire.edges {
        let edge = brep.edges.get(we.idx)?;
        let idx = if we.forward { edge.start } else { edge.end };
        let v = brep.vertices.get(idx)?;
        let pt = v.point;
        if pts.last().map_or(false, |last: &DVec3| (*last - pt).length_squared() < 1e-12) { continue; }
        pts.push(pt);
        // For self-loop edges (start == end), also push the midpoint
        // to ensure the polygon has all the boundary points.
        if edge.start == edge.end {
            // Get the edge curve and sample midpoint
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

/// Project 3D polygon vertices onto a 2D plane, returning 2D coordinates.
fn project_to_2d(polygon: &[DVec3], normal: DVec3) -> Vec<glam::DVec2> {
    let u_axis = normal.any_orthonormal_pair().0;
    let v_axis = normal.cross(u_axis).normalize();
    polygon.iter().map(|p| {
        glam::DVec2::new(p.dot(u_axis), p.dot(v_axis))
    }).collect()
}

/// Map 2D coordinates back to 3D using the plane basis.
fn map_to_3d(points_2d: &[glam::DVec2], normal: DVec3, origin: DVec3) -> Vec<DVec3> {
    let u_axis = normal.any_orthonormal_pair().0;
    let v_axis = normal.cross(u_axis).normalize();
    points_2d.iter().map(|p| {
        origin + u_axis * p.x + v_axis * p.y
    }).collect()
}

/// Compute the signed area of a 2D polygon (positive for CCW).
fn signed_area_2d(polygon: &[glam::DVec2]) -> f64 {
    let mut area = 0.0;
    for i in 0..polygon.len() {
        let j = (i + 1) % polygon.len();
        area += polygon[i].x * polygon[j].y - polygon[j].x * polygon[i].y;
    }
    area * 0.5
}

/// Determine if vertex i is convex (true) or concave (false).
/// Uses cross product of edge (i-1 鈫?i) and (i 鈫?i+1).
fn vertex_convexity(polygon: &[glam::DVec2], i: usize) -> bool {
    let n = polygon.len();
    let prev = polygon[(i + n - 1) % n];
    let curr = polygon[i];
    let next = polygon[(i + 1) % n];
    let cross = (curr - prev).perp_dot(next - curr);
    cross >= 0.0 // CCW polygon: positive = convex, negative = concave
}

/// Compute the outward normal of edge i in 2D (perpendicular to edge direction).
fn edge_outward_normal(polygon: &[glam::DVec2], i: usize) -> glam::DVec2 {
    let n = polygon.len();
    let from = polygon[i];
    let to = polygon[(i + 1) % n];
    let dir = (to - from).normalize();
    // For CCW polygon, outward normal is to the RIGHT of the edge direction
    // (right = dir rotated -90掳 = (dir.y, -dir.x))
    // Actually for CCW, the interior is to the LEFT:
    // left normal = (-dir.y, dir.x)
    // outward = -left = (dir.y, -dir.x)
    glam::DVec2::new(dir.y, -dir.x)
}

/// Offset a 2D polygon outward by distance d.
///
/// Computes the Minkowski sum of the polygon with a circle of radius d.
/// For JoinType::Intersection (sharp corners):
/// - Convex vertices: intersection of adjacent offset edges
/// - Concave vertices: offset edges are extended to intersect
///
/// Returns the offset polygon (possibly with fewer vertices if narrow features close).
/// Offset a 2D polygon outward by distance d.
///
/// Computes the Minkowski sum of the polygon with a circle of radius d.
/// Properly handles narrow-feature closure by detecting self-intersections
/// in the raw offset and removing interior loops via convex hull of the
/// split polygon (correct for JoinType::Intersection where narrow features
/// close completely).
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

    // Find all edge-edge intersection points
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
        return raw2; // No self-intersections
    }

    // Step 3: Build split polygon by inserting intersection points along edges.
    // For each edge, store its split points sorted by parameter t.
    let mut edge_splits: Vec<Vec<(f64, glam::DVec2)>> = vec![Vec::new(); m];
    for x in &xs {
        let pt_a = raw2[x.i] + (raw2[(x.i + 1) % m] - raw2[x.i]) * x.t_a;
        let pt_b = raw2[x.j] + (raw2[(x.j + 1) % m] - raw2[x.j]) * x.t_b;
        edge_splits[x.i].push((x.t_a, pt_a));
        edge_splits[x.j].push((x.t_b, pt_b));
    }

    // Build the split polygon by walking raw2 and inserting split points
    let mut split_pts: Vec<glam::DVec2> = Vec::new();
    for i in 0..m {
        // Start vertex of this edge
        let start_pt = raw2[i];
        if split_pts.is_empty() || (start_pt - *split_pts.last().unwrap()).length_squared() > 1e-20 {
            split_pts.push(start_pt);
        }
        // Insert split points along this edge in order of t
        if !edge_splits[i].is_empty() {
            edge_splits[i].sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            for (_, pt) in &edge_splits[i] {
                let last = *split_pts.last().unwrap();
                if (*pt - last).length_squared() > 1e-20 {
                    split_pts.push(*pt);
                }
            }
        }
    }

    if split_pts.len() < 3 {
        return raw2;
    }

/// Build an offset BRep for a prismatic solid using the 2D polygon offset.
    // Step 4: Walk the outer boundary at intersection points.
    // At each self-intersection, cross to the intersecting edge instead of
    // continuing along the current edge. This traces the outer boundary of
    // the Minkowski sum, removing interior loops from narrow-feature closure.

    // Build an adjacency map: for each intersection point P on raw edge E,
    // the cross target is the END vertex of the OTHER intersecting edge.
    use std::collections::HashMap;

    // Map: (intersection_point_index_in_split_pts) -> target_index_in_split_pts
    let mut cross_map: HashMap<usize, usize> = HashMap::new();

    for x in &xs {
        // Intersection points on edges i and j
        let pt_i = raw2[x.i] + (raw2[(x.i + 1) % m] - raw2[x.i]) * x.t_a;
        let pt_j = raw2[x.j] + (raw2[(x.j + 1) % m] - raw2[x.j]) * x.t_b;

        // The target raw vertices to cross to
        let target_i = (x.i + 1) % m;  // end of edge i
        let target_j = (x.j + 1) % m;  // end of edge j

        // Cross from edge i to end of edge j
        // Cross from edge j to end of edge i
        // Determine which edge each split point belongs to by checking
        // proximity to the end vertices of each edge.
        let edge_i_start = raw2[x.i];
        let edge_i_end = raw2[(x.i + 1) % m];
        let edge_j_start = raw2[x.j];
        let edge_j_end = raw2[(x.j + 1) % m];

        for spi in 0..split_pts.len() {
            let d_i = (split_pts[spi] - pt_i).length_squared();
            let d_j = (split_pts[spi] - pt_j).length_squared();
            if d_i < 1e-12 && d_j < 1e-12 {
                let on_edge_i = {
                    let dir = edge_i_end - edge_i_start;
                    let dir_len_sq = dir.length_squared();
                    if dir_len_sq < 1e-20 { false }
                    else {
                        let t = (split_pts[spi] - edge_i_start).dot(dir) / dir_len_sq;
                        t >= -0.01 && t <= 1.01
                    }
                };
                let on_edge_j = {
                    let dir = edge_j_end - edge_j_start;
                    let dir_len_sq = dir.length_squared();
                    if dir_len_sq < 1e-20 { false }
                    else {
                        let t = (split_pts[spi] - edge_j_start).dot(dir) / dir_len_sq;
                        t >= -0.01 && t <= 1.01
                    }
                };
                if on_edge_i && !on_edge_j {
                    for spj in 0..split_pts.len() {
                        if (split_pts[spj] - raw2[target_j]).length_squared() < 1e-12 {
                            cross_map.insert(spi, spj);
                            break;
                        }
                    }
                } else if on_edge_j && !on_edge_i {
                    for spj in 0..split_pts.len() {
                        if (split_pts[spj] - raw2[target_i]).length_squared() < 1e-12 {
                            cross_map.insert(spi, spj);
                            break;
                        }
                    }
                }
            }
        }
    }

    // Walk the outer boundary
    let start_idx = split_pts.iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let mut outer: Vec<glam::DVec2> = Vec::new();
    let mut current = start_idx;
    let mut visited = vec![false; split_pts.len()];
    let mut iter_count = 0;
    let max_iter = split_pts.len() * 4;

    loop {
        if visited[current] { break; }
        visited[current] = true;
        outer.push(split_pts[current]);

        // Check if we should cross to another edge here
        let next = if let Some(&cross_target) = cross_map.get(&current) {
            cross_target
        } else {
            (current + 1) % split_pts.len()
        };

        current = next;
        iter_count += 1;
        if current == start_idx || iter_count > max_iter { break; }
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

    // Project polygon to 2D
    let poly_2d = project_to_2d(&info.polygon_3d, info.cap_normal);

    // Ensure CCW winding
    let area = signed_area_2d(&poly_2d);
    let poly_2d = if area < 0.0 {
        let mut rev = poly_2d.clone();
        rev.reverse();
        rev
    } else {
        poly_2d
    };

    // Compute 2D offset
    let offset_poly_2d = offset_polygon_2d(&poly_2d, d);
    if offset_poly_2d.len() < 3 {
        return None;
    }

    // Map offset polygon to 3D on the bottom cap plane (shifted by d along cap normal)
    let bottom_origin = info.cap_origin + info.cap_normal * d;
    let offset_3d = map_to_3d(&offset_poly_2d, info.cap_normal, bottom_origin);

    // Extrude to create the offset solid
    let extrusion_height = info.extrusion_height + 2.0 * d.abs();
    let brep = crate::features::extrude_polygon_solid(
        &offset_3d,
        info.extrusion_dir,
        extrusion_height,
    ).ok()?;

    Some(brep)
}
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
