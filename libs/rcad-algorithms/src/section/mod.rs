//! Section: surface-solid intersection returning curves and wires.
//!
//! Analogous to OCCT `BRepAlgoAPI_Section`. Computes the intersection of a
//! cutting surface with the faces of a BRep, returning curves and wires.
//!
//! # Capabilities
//!
//! - Section by plane (original)
//! - Section by cylinder (cylindrical cut)
//! - Section by sphere (spherical cut)
//! - Section by arbitrary BRep surface
//! - Section by arbitrary analytic surface (cone, torus, etc.)
//!
//! - Returns exact analytic curves when possible (circle, ellipse, line)
//! - BSpline approximation for general cases
//! - Handles closed loops properly
//!
//! - Computes section properties: area, centroid, moments of inertia, perimeter
//!
//! - Multiple section support: parallel planes, cross-sections along a path

use glam::DVec3;
use glam::DVec2;
use rcad_kernel::{face_tolerance};
use rcad_kernel::topods::{BRep, TShape, ShapeRef};
pub use rcad_kernel::topods;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve3, CurveEval, CylindricalSurface, Ellipse3, Line3, Plane,
    SphericalSurface, Surface3, ToroidalSurface, any_perpendicular,
};
use std::f64::consts::PI;

use crate::inttools::{
    intersect_surfaces_with_density_tol, SurfaceCurve, SurfaceIntersectionResult,
};
use crate::inttools::plane_plane::PlanePlaneResult;
use crate::inttools::plane_cylinder::PlaneCylinderResult;
use crate::inttools::plane_sphere::PlaneSphereResult;
use crate::inttools::plane_cone::PlaneConicalResult;
use crate::tolerance::{
    intss_geom_tol_floor,
    max_face_tolerance_or_abs_pair,
    tessellation_merge_linear_from_brep,
    tessellation_merge_linear_from_two_breps,
    TOLERANCE_ABS,
    TOLERANCE_AREA_REL,
    TOLERANCE_CLAMP_MIN,
    TOLERANCE_COORD_SUB,
    TOLERANCE_LEN_MIN,
    TOLERANCE_MESH_LEGACY,
};
use crate::triangulate::{mesh_brep, triangulate_polygon, TessellationParams};

// = =  Helpers for extracting face data from topods::BRep  = = = = = = = = = =

/// Iterate all Solid ShapeRefs in a BRep.
fn iter_solid_refs(brep: &BRep) -> Vec<ShapeRef> {
    brep.tshapes.iter().enumerate().filter_map(|(i, ts)| {
        if matches!(ts.as_ref(), TShape::Solid(_)) {
            Some(ShapeRef::synthetic(i))
        } else {
            None
        }
    }).collect()
}

/// Iterate all shell ShapeRefs inside a solid.
fn iter_shell_refs(brep: &BRep, solid: ShapeRef) -> Vec<ShapeRef> {
    match &*brep.tshapes[solid.index] {
        TShape::Solid(sd) => sd.shells.iter().copied().collect(),
        _ => vec![],
    }
}

/// Iterate all face ShapeRefs inside a shell.
fn iter_face_refs(brep: &BRep, shell: ShapeRef) -> Vec<ShapeRef> {
    match &*brep.tshapes[shell.index] {
        TShape::Shell(sd) => sd.faces.iter().copied().collect(),
        _ => vec![],
    }
}

/// Get the surface of a face given its ShapeRef.
fn face_surface_from_ref(brep: &BRep, face: ShapeRef) -> Option<&Surface3> {
    match &*brep.tshapes[face.index] {
        TShape::Face(fd) => fd.surface.as_ref(),
        _ => None,
    }
}

/// Extract all boundary points of a face's outer wire by walking ShapeRef chain.
fn face_wire_points(brep: &BRep, face: ShapeRef) -> Vec<DVec3> {
    let outer_wire = match &*brep.tshapes[face.index] {
        TShape::Face(fd) => fd.outer_wire,
        _ => return vec![],
    };
    let edges = match &*brep.tshapes[outer_wire.index] {
        TShape::Wire(wd) => &wd.edges,
        _ => return vec![],
    };
    let mut pts = Vec::new();
    for edge_sr in edges {
        match &*brep.tshapes[edge_sr.index] {
            TShape::Edge(ed) => {
                // Use first vertex position (forward oriented)
                let loc = brep.get_location(edge_sr.location);
                let p = loc.transform_point3(brep.vertex(ed.first).point);
                pts.push(p);
            }
            _ => {}
        }
    }
    pts
}

/// Collect all face triangle vertices from a BRep's TShapes.
/// Walks Solid → Shell → Face → wire vertices, returns Vec<[DVec3; 3]>.
fn collect_face_triangles_tshape(brep: &BRep, face_ref: ShapeRef) -> Vec<[DVec3; 3]> {
    let fd = match &*brep.tshapes[face_ref.index] {
        TShape::Face(fd) => fd,
        _ => return vec![],
    };
    // Walk outer wire edges to get boundary vertices
    let wire_sr = fd.outer_wire;
    let edges = match &*brep.tshapes[wire_sr.index] {
        TShape::Wire(wd) => &wd.edges,
        _ => return vec![],
    };
    let wire_pts: Vec<DVec3> = edges.iter().map(|edge_sr| {
        match &*brep.tshapes[edge_sr.index] {
            TShape::Edge(ed) => {
                let loc = brep.get_location(edge_sr.location);
                loc.transform_point3(brep.vertex(ed.first).point)
            }
            _ => DVec3::ZERO,
        }
    }).collect();
    if wire_pts.len() < 3 {
        return vec![];
    }
    // Compute normal for fan triangulation
    let normal = (wire_pts[1] - wire_pts[0]).cross(wire_pts[2] - wire_pts[0]).normalize_or_zero();
    let tris = crate::triangulate::triangulate_polygon(&wire_pts, normal);
    tris.iter().filter_map(|&[i, j, k]| {
        let a = wire_pts.get(i)?;
        let b = wire_pts.get(j)?;
        let c = wire_pts.get(k)?;
        Some([*a, *b, *c])
    }).collect()
}

/// Find the global flat face index for a ShapeRef in the BRep's TShape tree.
/// Walks all Solid TShapes and their Shell→Face chains in order.
fn face_flat_index(brep: &BRep, target_face: ShapeRef) -> Option<usize> {
    let mut idx = 0usize;
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            for shell_sr in &sd.shells {
                if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if face_sr.is_same(&target_face) {
                            return Some(idx);
                        }
                        idx += 1;
                    }
                }
            }
        }
    }
    None
}

// = =  Internal helpers = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Flat BRep face index => geometric floor for numerical IntSS (phase C).
#[inline]
fn face_geom_floor(brep: &BRep, flat_face_idx: usize) -> f64 {
    intss_geom_tol_floor(TOLERANCE_ABS, face_tolerance(brep, flat_face_idx))
}

/// World-space merge / closed-loop slack for triangle-soup plane section (phase C).
#[inline]
fn plane_section_mesh_merge_eps(brep: &BRep) -> f64 {
    tessellation_merge_linear_from_brep(brep)
}

/// Signed distance from a point to the plane (positive on the normal side).
#[inline]
fn plane_dist(plane: &Plane, p: DVec3) -> f64 {
    plane.normal.dot(p - plane.origin)
}

/// Intersect a line segment (a, b) with the plane.
/// Returns the intersection point if the segment straddles the plane.
fn segment_plane_intersect(plane: &Plane, a: DVec3, b: DVec3) -> Option<DVec3> {
    let da = plane_dist(plane, a);
    let db = plane_dist(plane, b);
    if da.signum() == db.signum() || (da.abs() < TOLERANCE_ABS && db.abs() < TOLERANCE_ABS) {
        return None;
    }
    if da.abs() < TOLERANCE_ABS {
        return Some(a);
    }
    if db.abs() < TOLERANCE_ABS {
        return Some(b);
    }
    let t = da / (da - db);
    Some(a + t * (b - a))
}

/// Collect triangles for a face (pre-triangulated or fan-triangulated from wire).
/// Uses TShape-based access for the new BRep API.
fn face_triangles(brep: &BRep, face_ref: ShapeRef) -> Vec<[DVec3; 3]> {
    collect_face_triangles_tshape(brep, face_ref)
}

/// Collect all face triangles from a BRep (TShape-based).
///
/// Used for mesh-based intersection approximations (e.g. OCCT-style section).
pub fn brep_triangle_soup(brep: &BRep) -> Vec<[DVec3; 3]> {
    let mut out = Vec::new();
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            for shell_sr in &sd.shells {
                if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        out.extend(collect_face_triangles_tshape(brep, *face_sr));
                    }
                }
            }
        }
    }
    out
}

/// Intersect two triangle soups and chain intersection segments into polylines.
///
/// Uses [`TOLERANCE_MESH_LEGACY`] for pair/merge (historic default; unchanged for back-compat).
/// When both [`BRep`] operands are known, prefer [`intersect_triangle_soups_adaptive`] or [`intersect_triangle_soups_for_brep_tolerance`].
pub fn intersect_triangle_soups(tris_a: &[[DVec3; 3]], tris_b: &[[DVec3; 3]]) -> Vec<Vec<DVec3>> {
    intersect_triangle_soups_eps(tris_a, tris_b, TOLERANCE_MESH_LEGACY, TOLERANCE_MESH_LEGACY)
}

/// Triangle-triangle chaining with pair/merge epsilon from [`crate::tolerance::tessellation_merge_linear_from_two_breps`]
/// (**Relaxed adaptive** + **`TOLERANCE_MESH_LEGACY`** minimum + pairwise **`model_tolerance`**).
pub fn intersect_triangle_soups_adaptive(
    tris_a: &[[DVec3; 3]],
    tris_b: &[[DVec3; 3]],
    brep_a: &BRep,
    brep_b: &BRep,
) -> Vec<Vec<DVec3>> {
    let e = crate::tolerance::tessellation_merge_linear_from_two_breps(brep_a, brep_b);
    intersect_triangle_soups_eps(tris_a, tris_b, e, e)
}

/// Intersect two triangle soups; `pair_eps` feeds [`triangle_triangle_intersect_eps`],
/// `merge_eps` feeds [`chain_segments_eps`]. Both clamp up to [`TOLERANCE_ABS`].
pub fn intersect_triangle_soups_eps(
    tris_a: &[[DVec3; 3]],
    tris_b: &[[DVec3; 3]],
    pair_eps: f64,
    merge_eps: f64,
) -> Vec<Vec<DVec3>> {
    let pair = pair_eps.max(TOLERANCE_ABS);
    let merge = merge_eps.max(TOLERANCE_ABS);
    let mut segments: Vec<[DVec3; 2]> = Vec::new();
    for ta in tris_a {
        for tb in tris_b {
            if let Some(seg) = triangle_triangle_intersect_eps(ta, tb, pair) {
                segments.push(seg);
            }
        }
    }
    chain_segments_eps(segments, merge)
}

/// Intersect two triangle soups using [`crate::tolerance::max_face_tolerance_or_abs_pair`] as both
/// pair and merge epsilon (see [`intersect_triangle_soups_eps`]).
pub fn intersect_triangle_soups_for_brep_tolerance(
    tris_a: &[[DVec3; 3]],
    tris_b: &[[DVec3; 3]],
    brep_a: &BRep,
    brep_b: &BRep,
) -> Vec<Vec<DVec3>> {
    let e = max_face_tolerance_or_abs_pair(brep_a, brep_b);
    intersect_triangle_soups_eps(tris_a, tris_b, e, e)
}

/// Intersect a single triangle with the plane. Returns a segment [p0, p1] if
/// the triangle straddles the plane, or `None` otherwise.
fn triangle_section(plane: &Plane, tri: [DVec3; 3]) -> Option<[DVec3; 2]> {
    triangle_section_eps(plane, tri, TOLERANCE_ABS)
}

fn triangle_section_eps(plane: &Plane, tri: [DVec3; 3], dedup_eps: f64) -> Option<[DVec3; 2]> {
    let dedup = dedup_eps.max(TOLERANCE_CLAMP_MIN);
    let [a, b, c] = tri;
    let edges = [[a, b], [b, c], [c, a]];
    let mut pts = Vec::new();
    for [p, q] in edges {
        if let Some(hit) = segment_plane_intersect(plane, p, q) {
            // Deduplicate near-identical hits (e.g. at a vertex)
            if pts.iter().all(|&x: &DVec3| (x - hit).length() > dedup) {
                pts.push(hit);
            }
        }
    }
    if pts.len() >= 2 {
        Some([pts[0], pts[1]])
    } else {
        None
    }
}

/// Check if two points are close (within tolerance).
#[inline]
fn pts_close_eps(a: DVec3, b: DVec3, eps: f64) -> bool {
    (a - b).length() < eps
}

fn chain_segments_eps(segments: Vec<[DVec3; 2]>, merge_eps: f64) -> Vec<Vec<DVec3>> {
    if segments.is_empty() {
        return Vec::new();
    }

    let merge_eps = merge_eps.max(TOLERANCE_CLAMP_MIN);

    // Represent each segment as (start, end); build adjacency by proximity
    let mut remaining: Vec<[DVec3; 2]> = segments;
    let mut chains: Vec<Vec<DVec3>> = Vec::new();

    while !remaining.is_empty() {
        // Start a new chain with the first segment
        let first = remaining.remove(0);
        let mut chain = vec![first[0], first[1]];

        // Extend forward
        let mut extended = true;
        while extended {
            extended = false;
            let tail = *chain.last().expect("chain is non-empty (initialized with 2 points)");
            for i in 0..remaining.len() {
                if pts_close_eps(remaining[i][0], tail, merge_eps) {
                    chain.push(remaining[i][1]);
                    remaining.remove(i);
                    extended = true;
                    break;
                } else if pts_close_eps(remaining[i][1], tail, merge_eps) {
                    chain.push(remaining[i][0]);
                    remaining.remove(i);
                    extended = true;
                    break;
                }
            }
        }

        // Extend backward
        let mut extended = true;
        while extended {
            extended = false;
            let head = chain[0];
            for i in 0..remaining.len() {
                if pts_close_eps(remaining[i][1], head, merge_eps) {
                    chain.insert(0, remaining[i][0]);
                    remaining.remove(i);
                    extended = true;
                    break;
                } else if pts_close_eps(remaining[i][0], head, merge_eps) {
                    chain.insert(0, remaining[i][1]);
                    remaining.remove(i);
                    extended = true;
                    break;
                }
            }
        }

        chains.push(chain);
    }

    chains
}

// = =  Public API: Plane Section = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Compute the section of a BRep with a cutting plane.
///
/// Returns a new BRep containing only edges and wires (no faces/solids)
/// representing the section curves. Each closed loop is a separate wire.
///
/// Triangle-soup segment chaining uses [`crate::tolerance::tessellation_merge_linear_from_brep`]
/// (phase C), not a fixed [`TOLERANCE_MESH_LEGACY`] merge alone.
///
/// Analogous to OCCT `BRepAlgoAPI_Section`.
pub fn section(brep: &BRep, plane: &Plane) -> BRep {
    // Collect all section segments from all triangles
    let mut segments: Vec<[DVec3; 2]> = Vec::new();

    for solid_ref in iter_solid_refs(brep) {
        for shell_ref in iter_shell_refs(brep, solid_ref) {
            for face_ref in iter_face_refs(brep, shell_ref) {
                for tri in collect_face_triangles_tshape(brep, face_ref) {
                    if let Some(seg) = triangle_section(plane, tri) {
                        segments.push(seg);
                    }
                }
            }
        }
    }

    if segments.is_empty() {
        return BRep::new();
    }

    // Chain segments into loops
    let merge_eps = plane_section_mesh_merge_eps(brep);
    let loops = chain_segments_eps(segments, merge_eps);

    // Build result BRep using TShape API
    build_brep_from_polylines(&loops)
}

/// Pre-computed planar face info for fast 2D point-in-polygon.
struct PlanarFaceInfo {
    plane: Plane,
    x_axis: DVec3,
    y_axis: DVec3,
    wire_2d: Vec<DVec2>,
}

/// Extract planar face info if the face at `flat_idx` (tshape index) is a Plane surface.
fn extract_planar_face_info(brep: &BRep, flat_face_idx: usize) -> Option<PlanarFaceInfo> {
    // Find the face at flat_face_idx by walking the TShape tree
    let (face_ref, _) = find_face_by_flat_index(brep, flat_face_idx)?;
    let surface = face_surface_from_ref(brep, face_ref)?;
    let plane = match surface {
        Surface3::Plane(p) => *p,
        _ => return None,
    };
    let x_axis = plane.u_dir;
    let y_axis = plane.v_dir;

    let wire_pts = face_wire_points(brep, face_ref);
    let wire_2d: Vec<DVec2> = wire_pts.iter().map(|p| {
        let v = *p - plane.origin;
        DVec2::new(v.dot(x_axis), v.dot(y_axis))
    }).collect();

    if wire_2d.len() >= 3 {
        Some(PlanarFaceInfo { plane, x_axis, y_axis, wire_2d })
    } else {
        None
    }
}

/// Find a face ShapeRef and its flat count by flat index (walking Solid→Shell→Face in order).
fn find_face_by_flat_index(brep: &BRep, target: usize) -> Option<(ShapeRef, usize)> {
    let mut idx = 0usize;
    for ts in &brep.tshapes {
        if let TShape::Solid(sd) = ts.as_ref() {
            for shell_sr in &sd.shells {
                if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    for &face_sr in &shd.faces {
                        if idx == target {
                            return Some((face_sr, idx));
                        }
                        idx += 1;
                    }
                }
            }
        }
    }
    None
}

/// 2D ray-casting point-in-polygon test with tolerance margin.
fn point_in_polygon_2d(p: DVec2, poly: &[DVec2], tol: f64) -> bool {
    let mut inside = false;
    let mut j = poly.len() - 1;
    for i in 0..poly.len() {
        let yi = poly[i].y;
        let yj = poly[j].y;
        if ((yi + tol > p.y) != (yj + tol > p.y))
            || ((yi - tol > p.y) != (yj - tol > p.y))
        {
            let x_intersect = poly[j].x
                + (poly[i].x - poly[j].x) * (p.y - yj) / (yi - yj);
            if p.x < x_intersect + tol {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Check whether `point` (which lies in the plane of the planar face) is
/// inside the face's boundary polygon.
fn point_in_planar_face(
    point: DVec3,
    info: &PlanarFaceInfo,
    tol: f64,
) -> bool {
    let v = point - info.plane.origin;
    let p2 = DVec2::new(v.dot(info.x_axis), v.dot(info.y_axis));
    point_in_polygon_2d(p2, &info.wire_2d, tol)
}

/// Sample an infinite Line3 and keep the portion(s) within both planar faces.
/// Returns polylines (possibly multiple segments if the line re-enters).
fn sample_line_trimmed_to_planar_faces(
    line: &Line3,
    info_a: Option<&PlanarFaceInfo>,
    info_b: Option<&PlanarFaceInfo>,
    point_tol: f64,
) -> Vec<Vec<DVec3>> {
    // Sample over a generous range based on the line direction.
    let n = 2000usize;
    let t_range = 100.0;
    let pts: Vec<DVec3> = (0..n)
        .map(|i| line.point_at(-t_range + 2.0 * t_range * i as f64 / (n - 1) as f64))
        .collect();

    // Classify each point.
    let in_both: Vec<bool> = pts
        .iter()
        .map(|p| {
            info_a.map_or(true, |ia| point_in_planar_face(*p, ia, point_tol))
                && info_b.map_or(true, |ib| point_in_planar_face(*p, ib, point_tol))
        })
        .collect();

    // Extract contiguous runs.
    let mut result = Vec::new();
    let mut i = 0;
    while i < pts.len() {
        if in_both[i] {
            let start = i;
            while i < pts.len() && in_both[i] {
                i += 1;
            }
            if i - start >= 2 {
                result.push(pts[start..i].to_vec());
            }
        } else {
            i += 1;
        }
    }
    result
}

/// Sample a closed Curve3 and keep the portion(s) within a planar face.
fn sample_closed_curve_trimmed_to_planar_faces<C>(
    curve: &C,
    sample_fn: &dyn Fn(&C, usize) -> Vec<DVec3>,
    info_a: Option<&PlanarFaceInfo>,
    info_b: Option<&PlanarFaceInfo>,
    point_tol: f64,
) -> Vec<Vec<DVec3>> {
    let n = 128usize;
    let pts = sample_fn(curve, n);

    let in_both: Vec<bool> = pts
        .iter()
        .map(|p| {
            info_a.map_or(true, |ia| point_in_planar_face(*p, ia, point_tol))
                && info_b.map_or(true, |ib| point_in_planar_face(*p, ib, point_tol))
        })
        .collect();

    let mut result = Vec::new();
    let mut i = 0;
    while i < pts.len() {
        if in_both[i] {
            let start = i;
            while i < pts.len() && in_both[i] {
                i += 1;
            }
            if i - start >= 2 {
                result.push(pts[start..i].to_vec());
            }
        } else {
            i += 1;
        }
    }
    result
}

/// Helper to sample a Circle3.
fn sample_circle(c: &Circle3, n: usize) -> Vec<DVec3> {
    (0..n)
        .map(|i| c.point_at(2.0 * PI * i as f64 / (n - 1) as f64))
        .collect()
}

/// Helper to sample an Ellipse3.
fn sample_ellipse(e: &Ellipse3, n: usize) -> Vec<DVec3> {
    (0..n)
        .map(|i| e.point_at(2.0 * PI * i as f64 / (n - 1) as f64))
        .collect()
}

/// Push segments from a polyline into the segment list.
fn push_polyline_segments(polyline: &[DVec3], segments: &mut Vec<[DVec3; 2]>) {
    for w in polyline.windows(2) {
        segments.push([w[0], w[1]]);
    }
}

/// Try analytic intersection for a face pair, returning true if any segments
/// were produced.
fn try_analytic_face_pair(
    brep_a: &BRep,
    flat_a: usize,
    a_info: Option<&PlanarFaceInfo>,
    brep_b: &BRep,
    flat_b: usize,
    b_info: Option<&PlanarFaceInfo>,
    point_tol: f64,
    segments: &mut Vec<[DVec3; 2]>,
) -> bool {
    use Surface3::*;

    let a_surf = find_face_by_flat_index(brep_a, flat_a)
        .and_then(|(face_ref, _)| face_surface_from_ref(brep_a, face_ref));
    let b_surf = find_face_by_flat_index(brep_b, flat_b)
        .and_then(|(face_ref, _)| face_surface_from_ref(brep_b, face_ref));

    let (sa, sb) = match (a_surf, b_surf) {
        (Some(sa), Some(sb)) => (sa, sb),
        _ => return false,
    };

    match (sa, sb) {
        // = =  Plane vs Plane = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
        (Plane(pa), Plane(pb)) => {
            let result = crate::inttools::plane_plane::intersect_plane_plane(pa, pb);
            match result {
                PlanePlaneResult::Line(line) => {
                    let polylines = sample_line_trimmed_to_planar_faces(
                        &line,
                        a_info,
                        b_info,
                        point_tol,
                    );
                    for pl in &polylines {
                        push_polyline_segments(pl, segments);
                    }
                    !polylines.is_empty()
                }
                _ => false,
            }
        }

        // = =  Plane vs Cylinder = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
        (Plane(pa), Cylinder(cb)) => {
            intersect_plane_cylinder_pair(pa, cb, a_info, b_info, point_tol, segments)
        }
        (Cylinder(ca), Plane(pb)) => {
            intersect_plane_cylinder_pair(pb, ca, b_info, a_info, point_tol, segments)
        }

        // = =  Plane vs Sphere = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
        (Plane(pa), Sphere(sb)) => {
            intersect_plane_sphere_pair(pa, sb, a_info, b_info, point_tol, segments)
        }
        (Sphere(sa), Plane(pb)) => {
            intersect_plane_sphere_pair(pb, sa, b_info, a_info, point_tol, segments)
        }

        // = =  Plane vs Cone = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
        (Plane(pa), Cone(cb)) => {
            intersect_plane_cone_pair(pa, cb, a_info, b_info, point_tol, segments)
        }
        (Cone(ca), Plane(pb)) => {
            intersect_plane_cone_pair(pb, ca, b_info, a_info, point_tol, segments)
        }

        _ => false,
    }
}

fn intersect_plane_cylinder_pair(
    plane: &Plane,
    cyl: &CylindricalSurface,
    plane_info: Option<&PlanarFaceInfo>,
    cyl_info: Option<&PlanarFaceInfo>,
    point_tol: f64,
    segments: &mut Vec<[DVec3; 2]>,
) -> bool {
    let result = crate::inttools::plane_cylinder::intersect_plane_cylinder(plane, cyl);
    let mut found = false;
    match result {
        PlaneCylinderResult::TangentLine(line) => {
            let polylines =
                sample_line_trimmed_to_planar_faces(&line, plane_info, cyl_info, point_tol);
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneCylinderResult::TwoLines(l1, l2) => {
            for line in [l1, l2] {
                let polylines =
                    sample_line_trimmed_to_planar_faces(&line, plane_info, cyl_info, point_tol);
                for pl in &polylines {
                    push_polyline_segments(pl, segments);
                    found = true;
                }
            }
        }
        PlaneCylinderResult::Circle(c) => {
            let polylines = sample_closed_curve_trimmed_to_planar_faces(
                &c,
                &sample_circle,
                plane_info,
                cyl_info,
                point_tol,
            );
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneCylinderResult::Ellipse(e) => {
            let polylines = sample_closed_curve_trimmed_to_planar_faces(
                &e,
                &sample_ellipse,
                plane_info,
                cyl_info,
                point_tol,
            );
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneCylinderResult::NoIntersection => {}
    }
    found
}

fn intersect_plane_sphere_pair(
    plane: &Plane,
    sphere: &SphericalSurface,
    plane_info: Option<&PlanarFaceInfo>,
    sphere_info: Option<&PlanarFaceInfo>,
    point_tol: f64,
    segments: &mut Vec<[DVec3; 2]>,
) -> bool {
    let result = crate::inttools::plane_sphere::intersect_plane_sphere(plane, sphere);
    match result {
        PlaneSphereResult::Circle(c) => {
            let polylines = sample_closed_curve_trimmed_to_planar_faces(
                &c,
                &sample_circle,
                plane_info,
                sphere_info,
                point_tol,
            );
            for pl in &polylines {
                push_polyline_segments(pl, segments);
            }
            !polylines.is_empty()
        }
        _ => false,
    }
}

fn intersect_plane_cone_pair(
    plane: &Plane,
    cone: &ConicalSurface,
    plane_info: Option<&PlanarFaceInfo>,
    cone_info: Option<&PlanarFaceInfo>,
    point_tol: f64,
    segments: &mut Vec<[DVec3; 2]>,
) -> bool {
    let result = crate::inttools::plane_cone::intersect_plane_cone(plane, cone);
    let mut found = false;
    match result {
        PlaneConicalResult::Circle(c) => {
            let polylines = sample_closed_curve_trimmed_to_planar_faces(
                &c, &sample_circle, plane_info, cone_info, point_tol,
            );
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneConicalResult::Ellipse(e) => {
            let polylines = sample_closed_curve_trimmed_to_planar_faces(
                &e, &sample_ellipse, plane_info, cone_info, point_tol,
            );
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneConicalResult::Parabola(p) => {
            let pts: Vec<DVec3> = (0..128)
                .map(|i| p.point_at(-20.0 + 40.0 * i as f64 / 127.0))
                .collect();
            let in_both: Vec<bool> = pts.iter()
                .map(|p| {
                    plane_info.map_or(true, |ia| point_in_planar_face(*p, ia, point_tol))
                        && cone_info.map_or(true, |ib| point_in_planar_face(*p, ib, point_tol))
                })
                .collect();
            let mut i = 0;
            while i < pts.len() {
                if in_both[i] {
                    let start = i;
                    while i < pts.len() && in_both[i] { i += 1; }
                    if i - start >= 2 {
                        push_polyline_segments(&pts[start..i], segments);
                        found = true;
                    }
                } else { i += 1; }
            }
        }
        PlaneConicalResult::Hyperbola(h) => {
            let pts: Vec<DVec3> = (0..128)
                .map(|i| h.point_at(-10.0 + 20.0 * i as f64 / 127.0))
                .collect();
            let in_both: Vec<bool> = pts.iter()
                .map(|p| {
                    plane_info.map_or(true, |ia| point_in_planar_face(*p, ia, point_tol))
                        && cone_info.map_or(true, |ib| point_in_planar_face(*p, ib, point_tol))
                })
                .collect();
            let mut i = 0;
            while i < pts.len() {
                if in_both[i] {
                    let start = i;
                    while i < pts.len() && in_both[i] { i += 1; }
                    if i - start >= 2 {
                        push_polyline_segments(&pts[start..i], segments);
                        found = true;
                    }
                } else { i += 1; }
            }
        }
        PlaneConicalResult::SingleLine(l) | PlaneConicalResult::TwoLines(l, _) => {
            let polylines =
                sample_line_trimmed_to_planar_faces(&l, plane_info, cone_info, point_tol);
            for pl in &polylines {
                push_polyline_segments(pl, segments);
                found = true;
            }
        }
        PlaneConicalResult::Point(_) | PlaneConicalResult::NoIntersection => {}
    }
    found
}

/// Section curves between two BReps (all face-pair intersections).
///
/// Uses fast closed-form analytic intersection for analytic surface pairs
/// (plane-plane, plane-cylinder, plane-sphere, plane-cone) with 2D
/// point-in-polygon trimming to face boundaries. Falls back to triangle-soup
/// intersection for other surface types.
///
/// Like OCCT `BRepAlgoAPI_Section` applied to two whole shapes.
pub fn brep_section(a: &BRep, b: &BRep) -> BRep {
    // Tessellate both BReps so curved surfaces produce triangles for fallback.
    let mut a_tess = a.clone();
    let mut b_tess = b.clone();
    mesh_brep(&mut a_tess, &TessellationParams::standard());
    mesh_brep(&mut b_tess, &TessellationParams::analysis());

    let pair_eps = crate::tolerance::tessellation_merge_linear_from_two_breps(a, b)
        .max(TOLERANCE_ABS);
    let merge_eps = pair_eps;
    let point_tol = pair_eps.max(TOLERANCE_ABS);

    // Build flat face lists for both BReps (indexed by TShape tree order).
    let a_face_list: Vec<(ShapeRef, Option<PlanarFaceInfo>)> = {
        let mut list = Vec::new();
        let mut idx = 0usize;
        for ts in &a_tess.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                for shell_sr in &sd.shells {
                    if let TShape::Shell(shd) = &*a_tess.tshapes[shell_sr.index] {
                        for &face_sr in &shd.faces {
                            let info = extract_planar_face_info(&a_tess, idx);
                            list.push((face_sr, info));
                            idx += 1;
                        }
                    }
                }
            }
        }
        list
    };
    let b_face_list: Vec<(ShapeRef, Option<PlanarFaceInfo>)> = {
        let mut list = Vec::new();
        let mut idx = 0usize;
        for ts in &b_tess.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                for shell_sr in &sd.shells {
                    if let TShape::Shell(shd) = &*b_tess.tshapes[shell_sr.index] {
                        for &face_sr in &shd.faces {
                            let info = extract_planar_face_info(&b_tess, idx);
                            list.push((face_sr, info));
                            idx += 1;
                        }
                    }
                }
            }
        }
        list
    };

    // Collect segments from all non-coplanar face pairs.
    let mut segments: Vec<[DVec3; 2]> = Vec::new();
    {
        for (a_flat, (a_face_ref, a_info)) in a_face_list.iter().enumerate() {
            for (b_flat, (b_face_ref, b_info)) in b_face_list.iter().enumerate() {
                if !is_coplanar_face_pair(&a_tess, a_flat, &b_tess, b_flat) {
                    let used_analytic = try_analytic_face_pair(
                        &a_tess,
                        a_flat,
                        a_info.as_ref(),
                        &b_tess,
                        b_flat,
                        b_info.as_ref(),
                        point_tol,
                        &mut segments,
                    );

                    if !used_analytic {
                        // Fall back to triangle-triangle intersection.
                        let a_tris = collect_face_triangles_tshape(&a_tess, *a_face_ref);
                        let b_tris = collect_face_triangles_tshape(&b_tess, *b_face_ref);
                        for ta in &a_tris {
                            for tb in &b_tris {
                                if let Some(seg) =
                                    triangle_triangle_intersect_eps(ta, tb, pair_eps)
                                {
                                    segments.push(seg);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if segments.is_empty() {
        return BRep::new();
    }

    // Deduplicate colinear overlapping segments.
    let same_line_eps = merge_eps.max(0.01);
    dedup_colinear_segments(&mut segments, same_line_eps);

    let polylines = chain_segments_eps(segments, merge_eps);
    build_brep_from_polylines(&polylines)
}

/// Remove colinear overlapping segments, merging them so the chainer produces
/// clean non-oscillating polylines.
fn dedup_colinear_segments(segs: &mut Vec<[DVec3; 2]>, eps: f64) {
    let eps = eps.max(TOLERANCE_CLAMP_MIN);
    let mut i = 0;
    while i < segs.len() {
        let [a, b] = segs[i];
        let dir = (b - a).normalize();
        let len = (b - a).length();

        if len < eps {
            segs.swap_remove(i);
            continue;
        }

        let mut merged = false;
        let mut j = i + 1;
        while j < segs.len() {
            let [c, d] = segs[j];
            let c_dir = (d - c).normalize();
            let c_len = (d - c).length();

            if c_len < eps {
                segs.swap_remove(j);
                continue;
            }

            // Check colinearity and same offset
            if dir.dot(c_dir).abs() > 1.0 - eps {
                // Same line check
                let perp_dist = ((c - a) - (c - a).dot(dir) * dir).length();
                if perp_dist > eps {
                    j += 1;
                    continue;
                }
                // Project [c, d] onto the line defined by [a, b]
                let t_c = (c - a).dot(dir);
                let t_d = (d - a).dot(dir);
                let c_t0 = t_c.min(t_d);
                let c_t1 = t_c.max(t_d);

                // Check overlap
                if c_t1 >= -eps && c_t0 <= len + eps {
                    let new_t0 = 0.0f64.min(c_t0);
                    let new_t1 = len.max(c_t1);
                    segs[i] = [a + dir * new_t0, a + dir * new_t1];
                    segs.swap_remove(j);
                    merged = true;
                    break;
                }
            }
            j += 1;
        }

        if merged {
            continue;
        }
        i += 1;
    }
}

/// Check whether the two faces (identified by their flat indices) are coplanar
/// planes — i.e. both are [`Surface3::Plane`] with the same normal and origin
/// offset.
fn is_coplanar_face_pair(brep_a: &BRep, flat_a: usize, brep_b: &BRep, flat_b: usize) -> bool {
    let surf_a = find_face_by_flat_index(brep_a, flat_a)
        .and_then(|(face_ref, _)| face_surface_from_ref(brep_a, face_ref));
    let surf_b = find_face_by_flat_index(brep_b, flat_b)
        .and_then(|(face_ref, _)| face_surface_from_ref(brep_b, face_ref));

    match (surf_a, surf_b) {
        (Some(Surface3::Plane(pa)), Some(Surface3::Plane(pb))) => {
            // Same normal direction
            if (pa.normal.dot(pb.normal).abs() - 1.0).abs() > TOLERANCE_ABS {
                return false;
            }
            // Same distance from origin
            let dist = (pa.origin - pb.origin).dot(pa.normal).abs();
            dist < TOLERANCE_ABS
        }
        _ => false,
    }
}

/// Convenience: extract all section polylines as ordered lists of 3D points.
///
/// Each entry is one closed (or open) loop of points from the plane section.
/// Chaining tolerance follows [`crate::tolerance::tessellation_merge_linear_from_brep`].
pub fn section_polylines(brep: &BRep, plane: &Plane) -> Vec<Vec<DVec3>> {
    let mut segments: Vec<[DVec3; 2]> = Vec::new();

    for solid_ref in iter_solid_refs(brep) {
        for shell_ref in iter_shell_refs(brep, solid_ref) {
            for face_ref in iter_face_refs(brep, shell_ref) {
                for tri in collect_face_triangles_tshape(brep, face_ref) {
                    if let Some(seg) = triangle_section(plane, tri) {
                        segments.push(seg);
                    }
                }
            }
        }
    }

    let merge_eps = plane_section_mesh_merge_eps(brep);
    chain_segments_eps(segments, merge_eps)
}

// = =  Public API: Curved Surface Section = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Cutting surface for section operations.
///
/// Supports plane, cylinder, sphere, cone, torus, and arbitrary analytic surfaces.
#[derive(Debug, Clone)]
pub enum CuttingSurface {
    /// Planar cut (original behavior).
    Plane(Plane),
    /// Cylindrical cut.
    Cylinder(CylindricalSurface),
    /// Spherical cut.
    Sphere(SphericalSurface),
    /// Conical cut.
    Cone(ConicalSurface),
    /// Toroidal cut.
    Torus(ToroidalSurface),
    /// Arbitrary analytic surface.
    Surface(Surface3),
    /// Arbitrary BRep surface (uses face index).
    BRepSurface {
        /// Source BRep containing the cutting face.
        brep: Box<BRep>,
        /// Index of the face in the BRep to use as cutting surface.
        face_idx: usize,
    },
}

/// Result of a section operation with curves and properties.
#[derive(Debug, Clone)]
pub struct SectionResult {
    /// The section curves as a BRep (wires only).
    pub brep: BRep,
    /// Individual section curves (analytic or polyline).
    pub curves: Vec<SectionCurveResult>,
    /// Section properties (computed if section is planar and closed).
    pub properties: Option<SectionProperties>,
}

/// One curve from a section result.
#[derive(Debug, Clone)]
pub struct SectionCurveResult {
    /// The 3D curve (analytic or polyline approximation).
    pub curve: SectionCurveType,
    /// Whether this curve forms a closed loop.
    pub is_closed: bool,
    /// Parameter range for the curve.
    pub param_range: [f64; 2],
}

/// Type of section curve.
#[derive(Debug, Clone)]
pub enum SectionCurveType {
    /// Analytic line.
    Line(Line3),
    /// Analytic circle.
    Circle(Circle3),
    /// Analytic ellipse.
    Ellipse(Ellipse3),
    /// BSpline curve approximation.
    BSpline(rcad_kernel::geom::BSplineCurve3),
    /// Polyline (sampled points).
    Polyline(Vec<DVec3>),
}

impl SectionCurveType {
    /// Sample points on this curve for display or computation.
    pub fn sample_points(&self, n: usize) -> Vec<DVec3> {
        match self {
            SectionCurveType::Line(line) => {
                let [t0, t1] = [0.0, 100.0];
                (0..n)
                    .map(|i| line.point_at(t0 + (t1 - t0) * i as f64 / (n - 1).max(1) as f64))
                    .collect()
            }
            SectionCurveType::Circle(circle) => {
                (0..n)
                    .map(|i| circle.point_at(2.0 * PI * i as f64 / (n - 1).max(1) as f64))
                    .collect()
            }
            SectionCurveType::Ellipse(ellipse) => {
                (0..n)
                    .map(|i| ellipse.point_at(2.0 * PI * i as f64 / (n - 1).max(1) as f64))
                    .collect()
            }
            SectionCurveType::BSpline(bspline) => {
                let [t0, t1] = bspline.default_domain();
                (0..n)
                    .map(|i| bspline.point_at(t0 + (t1 - t0) * i as f64 / (n - 1).max(1) as f64))
                    .collect()
            }
            SectionCurveType::Polyline(pts) => pts.clone(),
        }
    }
}

/// Properties of a planar section.
#[derive(Debug, Clone, Copy)]
pub struct SectionProperties {
    /// Area of the section.
    pub area: f64,
    /// Centroid (center of mass) of the section.
    pub centroid: DVec3,
    /// Second moment of area about the centroidal X axis (Ixx).
    pub ixx: f64,
    /// Second moment of area about the centroidal Y axis (Iyy).
    pub iyy: f64,
    /// Product moment of area (Ixy).
    pub ixy: f64,
    /// Perimeter of the section.
    pub perimeter: f64,
}

impl SectionProperties {
    /// Compute the polar moment of inertia (J = Ixx + Iyy).
    pub fn polar_moment(&self) -> f64 {
        self.ixx + self.iyy
    }

    /// Compute principal moments and axes.
    ///
    /// Returns ((I1, I2), angle) where angle is the rotation from X axis to principal axis.
    pub fn principal_moments(&self) -> ((f64, f64), f64) {
        let avg = 0.5 * (self.ixx + self.iyy);
        let diff = 0.5 * (self.ixx - self.iyy);
        let rad = (diff * diff + self.ixy * self.ixy).sqrt();

        let i1 = avg + rad;
        let i2 = avg - rad;

        // Angle to principal axis (measured from X axis)
        let angle = 0.5 * (2.0 * self.ixy).atan2(self.ixx - self.iyy);

        ((i1, i2), angle)
    }
}

/// Compute the section of a BRep with a cutting surface.
///
/// This is the main entry point for section operations, supporting various
/// cutting surface types.
///
/// # Arguments
///
/// * `brep` - The BRep to section.
/// * `cutting_surface` - The surface to cut with.
///
/// # Returns
///
/// A `SectionResult` containing the section curves as a BRep, individual curve
/// data, and computed properties (if applicable).
pub fn section_with_surface(brep: &BRep, cutting_surface: &CuttingSurface) -> SectionResult {
    match cutting_surface {
        CuttingSurface::Plane(plane) => section_by_plane(brep, plane),
        CuttingSurface::Cylinder(cyl) => section_by_cylinder(brep, cyl),
        CuttingSurface::Sphere(sphere) => section_by_sphere(brep, sphere),
        CuttingSurface::Cone(cone) => section_by_cone(brep, cone),
        CuttingSurface::Torus(torus) => section_by_torus(brep, torus),
        CuttingSurface::Surface(surface) => section_by_analytic_surface(brep, surface),
        CuttingSurface::BRepSurface { brep: tool_brep, face_idx } => {
            section_by_brep_surface(brep, tool_brep, *face_idx)
        }
    }
}

/// Section by plane with full result.
fn section_by_plane(brep: &BRep, plane: &Plane) -> SectionResult {
    let merge_eps = plane_section_mesh_merge_eps(brep);
    let polylines = section_polylines(brep, plane);

    let mut curves = Vec::new();

    for polyline in &polylines {
        if polyline.len() < 2 {
            continue;
        }

        let is_closed = pts_close_eps(polyline[0], *polyline.last().unwrap(), merge_eps);

        // Try to fit a BSpline for smooth representation
        let curve = if polyline.len() >= 4 && !is_closed {
            match rcad_kernel::fit::interpolate_points(polyline) {
                Ok(bspline) => SectionCurveType::BSpline(bspline),
                Err(_) => SectionCurveType::Polyline(polyline.clone()),
            }
        } else {
            SectionCurveType::Polyline(polyline.clone())
        };

        curves.push(SectionCurveResult {
            curve,
            is_closed,
            param_range: [0.0, polyline.len() as f64],
        });
    }

    // Build BRep from polylines
    let result_brep = build_brep_from_polylines(&polylines);

    // Compute properties if section is planar
    let properties = compute_planar_section_properties(&polylines, plane);

    SectionResult {
        brep: result_brep,
        curves,
        properties,
    }
}

/// Section by cylinder.
fn section_by_cylinder(brep: &BRep, cyl: &CylindricalSurface) -> SectionResult {
    section_by_analytic_surface(brep, &Surface3::Cylinder(*cyl))
}

/// Section by sphere.
fn section_by_sphere(brep: &BRep, sphere: &SphericalSurface) -> SectionResult {
    section_by_analytic_surface(brep, &Surface3::Sphere(*sphere))
}

/// Section by cone.
fn section_by_cone(brep: &BRep, cone: &ConicalSurface) -> SectionResult {
    section_by_analytic_surface(brep, &Surface3::Cone(*cone))
}

/// Section by torus.
fn section_by_torus(brep: &BRep, torus: &ToroidalSurface) -> SectionResult {
    section_by_analytic_surface(brep, &Surface3::Torus(*torus))
}

/// Section by arbitrary analytic surface.
fn section_by_analytic_surface(brep: &BRep, cutting_surface: &Surface3) -> SectionResult {
    let mesh_merge_eps = plane_section_mesh_merge_eps(brep);
    let mut curves = Vec::new();
    let mut polylines = Vec::new();

    let a_face_list: Vec<(ShapeRef, Option<usize>)> = {
        let mut list = Vec::new();
        let mut idx = 0usize;
        for ts in &brep.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                for shell_sr in &sd.shells {
                    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                        for &face_sr in &shd.faces {
                            list.push((face_sr, Some(idx)));
                            idx += 1;
                        }
                    }
                }
            }
        }
        list
    };

    for (face_ref, flat_idx_opt) in &a_face_list {
        let flat_idx = flat_idx_opt.unwrap_or(0);
        let geom_floor = face_geom_floor(brep, flat_idx);

        // Get the analytic surface for this face from TShape data
        let face_surface = face_surface_from_ref(brep, *face_ref).cloned();

        if let Some(face_surf) = face_surface {
            let intersection = intersect_surfaces_with_density_tol(
                &face_surf,
                cutting_surface,
                48,
                geom_floor,
            );

            for curve_result in intersection.curves {
                let (curve, polyline, is_closed) =
                    convert_surface_curve(&curve_result, mesh_merge_eps);

                curves.push(SectionCurveResult {
                    curve,
                    is_closed,
                    param_range: [0.0, 1.0],
                });

                if let Some(pts) = polyline {
                    polylines.push(pts);
                }
            }
        } else {
            // Fall back to triangle-based section for non-analytic faces
            let face_polylines =
                section_face_by_surface_marching(brep, *face_ref, cutting_surface, geom_floor);
            for pts in &face_polylines {
                let is_closed = pts.len() > 2
                    && pts_close_eps(pts[0], *pts.last().unwrap(), geom_floor);
                curves.push(SectionCurveResult {
                    curve: SectionCurveType::Polyline(pts.clone()),
                    is_closed,
                    param_range: [0.0, pts.len() as f64],
                });
            }
            polylines.extend(face_polylines);
        }
    }

    let result_brep = build_brep_from_polylines(&polylines);

    SectionResult {
        brep: result_brep,
        curves,
        properties: None,
    }
}

/// Section by a face from another BRep.
fn section_by_brep_surface(brep: &BRep, tool_brep: &BRep, face_idx: usize) -> SectionResult {
    // Get the cutting surface from the tool BRep
    let cutting_surface = find_face_by_flat_index(tool_brep, face_idx)
        .and_then(|(face_ref, _)| face_surface_from_ref(tool_brep, face_ref))
        .cloned();

    match cutting_surface {
        Some(surface) => section_by_analytic_surface(brep, &surface),
        None => {
            // Fall back to triangle-based intersection
            let cutting_face_ref = find_face_by_flat_index(tool_brep, face_idx)
                .map(|(sr, _)| sr);
            let (polylines, merge_eps) =
                section_by_face_triangles(brep, tool_brep, face_idx, cutting_face_ref);

            let curves = polylines
                .iter()
                .map(|pts| SectionCurveResult {
                    curve: SectionCurveType::Polyline(pts.clone()),
                    is_closed: pts.len() > 2
                        && pts_close_eps(pts[0], *pts.last().unwrap(), merge_eps),
                    param_range: [0.0, pts.len() as f64],
                })
                .collect();

            let result_brep = build_brep_from_polylines(&polylines);

            SectionResult {
                brep: result_brep,
                curves,
                properties: None,
            }
        }
    }
}

/// Convert a SurfaceCurve from intersection to SectionCurveType.
fn convert_surface_curve(
    result: &SurfaceIntersectionResult,
    polyline_close_eps: f64,
) -> (SectionCurveType, Option<Vec<DVec3>>, bool) {
    match &result.curve_3d {
        SurfaceCurve::Line(line) => {
            let pts = (0..10)
                .map(|i| line.point_at(-50.0 + 100.0 * i as f64 / 9.0))
                .collect();
            (SectionCurveType::Line(*line), Some(pts), false)
        }
        SurfaceCurve::Circle(circle) => {
            let pts = (0..33)
                .map(|i| circle.point_at(2.0 * PI * i as f64 / 32.0))
                .collect();
            (SectionCurveType::Circle(*circle), Some(pts), true)
        }
        SurfaceCurve::Ellipse(ellipse) => {
            let pts = (0..33)
                .map(|i| ellipse.point_at(2.0 * PI * i as f64 / 32.0))
                .collect();
            (SectionCurveType::Ellipse(*ellipse), Some(pts), true)
        }
        SurfaceCurve::Parabola(parabola) => {
            let pts: Vec<DVec3> = (0..33)
                .map(|i| {
                    let t = -10.0 + 20.0 * i as f64 / 32.0;
                    parabola.point_at(t)
                })
                .collect();
            match rcad_kernel::fit::interpolate_points(&pts) {
                Ok(bspline) => (SectionCurveType::BSpline(bspline), Some(pts.clone()), false),
                Err(_) => (SectionCurveType::Polyline(pts.clone()), Some(pts), false),
            }
        }
        SurfaceCurve::Hyperbola(hyperbola) => {
            let pts: Vec<DVec3> = (0..33)
                .map(|i| {
                    let t = -5.0 + 10.0 * i as f64 / 32.0;
                    hyperbola.point_at(t)
                })
                .collect();
            match rcad_kernel::fit::interpolate_points(&pts) {
                Ok(bspline) => (SectionCurveType::BSpline(bspline), Some(pts.clone()), false),
                Err(_) => (SectionCurveType::Polyline(pts.clone()), Some(pts), false),
            }
        }
        SurfaceCurve::Point(_) => (SectionCurveType::Polyline(vec![]), None, false),
        SurfaceCurve::BSplineCurve(b) => {
            (SectionCurveType::BSpline((**b).clone()), None, false)
        }
        SurfaceCurve::Polyline(pts) => {
            let is_closed = pts.len() > 2
                && pts_close_eps(pts[0], *pts.last().unwrap(), polyline_close_eps);
            if pts.len() >= 4
                && let Ok(bspline) = rcad_kernel::fit::approximate_points(pts, (pts.len() / 2).max(4)) {
                return (SectionCurveType::BSpline(bspline), Some(pts.clone()), is_closed);
            }
            (SectionCurveType::Polyline(pts.clone()), Some(pts.clone()), is_closed)
        }
    }
}

/// Section a face by marching along a surface.
fn section_face_by_surface_marching(
    brep: &BRep,
    face_ref: ShapeRef,
    cutting_surface: &Surface3,
    geom_floor: f64,
) -> Vec<Vec<DVec3>> {
    let edge_eps = geom_floor.max(TOLERANCE_ABS);
    let bisect_tol =
        (edge_eps * TOLERANCE_AREA_REL).clamp(TOLERANCE_CLAMP_MIN, edge_eps.max(TOLERANCE_MESH_LEGACY));

    // Get triangles for the face
    let triangles = collect_face_triangles_tshape(brep, face_ref);

    // For each triangle, find intersection with cutting surface
    let mut segments: Vec<[DVec3; 2]> = Vec::new();

    for tri in triangles {
        if let Some(seg) =
            triangle_surface_intersect(&tri, cutting_surface, edge_eps, bisect_tol)
        {
            segments.push(seg);
        }
    }

    // Chain segments into polylines
    chain_segments_eps(segments, edge_eps)
}

/// Intersect a triangle with a surface.
fn triangle_surface_intersect(
    tri: &[DVec3; 3],
    surface: &Surface3,
    edge_eps: f64,
    bisect_tol: f64,
) -> Option<[DVec3; 2]> {
    // Sample points on triangle edges and find where surface distance changes sign
    let edges = [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]];

    let mut intersection_points = Vec::new();

    for [a, b] in edges {
        let n_samples = 10;
        let mut prev_dist = signed_distance_to_surface(a, surface);

        for i in 1..=n_samples {
            let t = i as f64 / n_samples as f64;
            let p = a.lerp(b, t);
            let dist = signed_distance_to_surface(p, surface);

            if prev_dist * dist < 0.0 || dist.abs() < edge_eps {
                let intersection = find_surface_intersection(a, b, surface, bisect_tol);
                if let Some(pt) = intersection {
                    if intersection_points
                        .iter()
                        .all(|&x: &DVec3| (x - pt).length() > edge_eps)
                    {
                        intersection_points.push(pt);
                    }
                }
            }

            prev_dist = dist;
        }
    }

    if intersection_points.len() >= 2 {
        Some([intersection_points[0], intersection_points[1]])
    } else {
        None
    }
}

/// Signed distance from a point to a surface.
///
/// Positive = outside, negative = inside (for closed surfaces).
fn signed_distance_to_surface(p: DVec3, surface: &Surface3) -> f64 {
    match surface {
        Surface3::Plane(plane) => {
            plane.normal.dot(p - plane.origin)
        }
        Surface3::Sphere(sphere) => {
            (p - sphere.center).length() - sphere.radius
        }
        Surface3::Cylinder(cyl) => {
            let axis = cyl.axis.normalize();
            let v = p - cyl.origin;
            let along = v.dot(axis);
            let perp = (v - axis * along).length();
            perp - cyl.radius
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            let apex = cone.apex_point();
            let v = p - apex;
            let along = v.dot(axis);
            let perp = (v - axis * along).length();
            let expected_radius = along * cone.half_angle_rad.tan();
            perp - expected_radius
        }
        Surface3::Torus(torus) => {
            let axis = torus.axis.normalize();
            let v = p - torus.center;
            let along = v.dot(axis);
            let perp_vec = v - axis * along;
            let _perp_dist = perp_vec.length();
            let major_circle_pt = torus.center + perp_vec.normalize_or_zero() * torus.major_radius;
            let dist_from_major = (p - major_circle_pt - axis * along).length();
            dist_from_major - torus.minor_radius
        }
        _ => {
            let proj = rcad_kernel::projection::closest_point_on_surface(surface, p, 8);
            (p - proj.point).length() * if proj.params.0.fract().abs() < 0.5 { 1.0 } else { -1.0 }
        }
    }
}

/// Find the intersection of a line segment with a surface using binary search.
fn find_surface_intersection(
    a: DVec3,
    b: DVec3,
    surface: &Surface3,
    bisect_tol: f64,
) -> Option<DVec3> {
    let dist_a = signed_distance_to_surface(a, surface);
    let dist_b = signed_distance_to_surface(b, surface);

    if dist_a * dist_b > 0.0 {
        return None;
    }

    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;

    for _ in 0..20 {
        let mid = 0.5 * (lo + hi);
        let p = a.lerp(b, mid);
        let dist_mid = signed_distance_to_surface(p, surface);

        if dist_mid.abs() < bisect_tol {
            return Some(p);
        }

        if dist_a * dist_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some(a.lerp(b, 0.5 * (lo + hi)))
}

/// Section by face triangles (for non-analytic surfaces).
fn section_by_face_triangles(
    brep: &BRep,
    tool_brep: &BRep,
    tool_face_idx: usize,
    cutting_face_ref: Option<ShapeRef>,
) -> (Vec<Vec<DVec3>>, f64) {
    let cutting_face_ref = match cutting_face_ref {
        Some(sr) => sr,
        None => return (Vec::new(), TOLERANCE_ABS),
    };

    let tool_floor = face_geom_floor(tool_brep, tool_face_idx);
    let cutting_triangles = collect_face_triangles_tshape(tool_brep, cutting_face_ref);

    let mut segments: Vec<[DVec3; 2]> = Vec::new();
    let mut merge_eps = tool_floor;

    // Build flat face list for brep
    let brep_face_list: Vec<ShapeRef> = {
        let mut list = Vec::new();
        for ts in &brep.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                for shell_sr in &sd.shells {
                    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                        for &face_sr in &shd.faces {
                            list.push(face_sr);
                        }
                    }
                }
            }
        }
        list
    };

    for (obj_idx, &face_ref) in brep_face_list.iter().enumerate() {
        let obj_floor = face_geom_floor(brep, obj_idx);
        merge_eps = merge_eps.max(obj_floor);
        let pair_eps = obj_floor.max(tool_floor);
        let brep_triangles = collect_face_triangles_tshape(brep, face_ref);

        for brep_tri in &brep_triangles {
            for cut_tri in &cutting_triangles {
                if let Some(seg) =
                    triangle_triangle_intersect_eps(brep_tri, cut_tri, pair_eps)
                {
                    segments.push(seg);
                }
            }
        }
    }

    merge_eps = merge_eps.max(tessellation_merge_linear_from_two_breps(brep, tool_brep));

    (chain_segments_eps(segments, merge_eps), merge_eps)
}

fn triangle_triangle_intersect_eps(
    tri1: &[DVec3; 3],
    tri2: &[DVec3; 3],
    pair_eps: f64,
) -> Option<[DVec3; 2]> {
    let pe = pair_eps.max(TOLERANCE_ABS);
    let degen_len = (TOLERANCE_LEN_MIN).max(pe * TOLERANCE_COORD_SUB);

    // Compute plane of tri2
    let normal2 = (tri2[1] - tri2[0]).cross(tri2[2] - tri2[0]);
    let len2 = normal2.length();
    if len2 < degen_len {
        return None;
    }
    let normal2 = normal2 / len2;
    let plane2 = Plane::new(tri2[0], normal2);

    // Find intersection of tri1 with plane of tri2
    let seg = triangle_section_eps(&plane2, *tri1, pe)?;

    // Clip segment to triangle 2 bounds
    clip_segment_to_triangle_eps(&seg, tri2, pe)
}

/// Clip a segment to a triangle's bounds.
fn clip_segment_to_triangle_eps(
    seg: &[DVec3; 2],
    tri: &[DVec3; 3],
    pair_eps: f64,
) -> Option<[DVec3; 2]> {
    let a_inside = point_in_triangle_eps(seg[0], tri, pair_eps);
    let b_inside = point_in_triangle_eps(seg[1], tri, pair_eps);

    if a_inside && b_inside {
        return Some(*seg);
    }

    if a_inside || b_inside {
        return Some(*seg);
    }

    None
}

/// Check if a point is inside a triangle (2D projection), with length-scale margin from `pair_eps`.
fn point_in_triangle_eps(p: DVec3, tri: &[DVec3; 3], pair_eps: f64) -> bool {
    let pe = pair_eps.max(TOLERANCE_ABS);
    let degen_len = (TOLERANCE_LEN_MIN).max(pe * TOLERANCE_COORD_SUB);

    let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
    let len = normal.length();
    if len < degen_len {
        return false;
    }
    let normal = normal / len;

    // Project to plane
    let v0 = tri[0];
    let v1 = tri[1] - tri[0];
    let v2 = tri[2] - tri[0];

    // Build local 2D basis
    let e1 = v1.normalize_or_zero();
    let e2 = normal.cross(e1).normalize_or_zero();

    let p_local = p - v0;
    let u = p_local.dot(e1);
    let v = p_local.dot(e2);

    let v1_local = DVec3::new(v1.length(), 0.0, 0.0);
    let v2_local = DVec3::new(v2.dot(e1), v2.dot(e2), 0.0);

    // Barycentric check
    let denom = v1_local.x * v2_local.y - v2_local.x * v1_local.y;
    if denom.abs() < degen_len {
        return false;
    }

    let s = (u * v2_local.y - v * v2_local.x) / denom;
    let t = (v * v1_local.x - u * v1_local.y) / denom;

    let e0 = tri[1] - tri[0];
    let e1_edge = tri[2] - tri[1];
    let e2_edge = tri[2] - tri[0];
    let scale = e0.length().max(e1_edge.length()).max(e2_edge.length()).max(TOLERANCE_LEN_MIN);
    let margin = (pe / scale).min(0.05).max(TOLERANCE_LEN_MIN);

    s >= -margin && t >= -margin && s + t <= 1.0 + margin
}

/// Build a BRep from polylines using the new TShape API.
fn build_brep_from_polylines(polylines: &[Vec<DVec3>]) -> BRep {
    let mut result = BRep::new();

    for polyline in polylines {
        if polyline.len() < 2 {
            continue;
        }

        let mut wire_edges = Vec::new();

        for i in 0..polyline.len().saturating_sub(1) {
            let a = polyline[i];
            let b = polyline[i + 1];

            let va = result.add_tvertex(a);
            let vb = result.add_tvertex(b);

            let len = (b - a).length();
            let dir = if len > TOLERANCE_ABS { (b - a) / len } else { DVec3::X };
            let curve = Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            });

            let edge = result.add_tedge(Some(curve), va, vb, [0.0, len]);
            wire_edges.push(edge);
        }

        if !wire_edges.is_empty() {
            let wire = result.add_twire(wire_edges);
            let face = result.add_tface(None, wire, vec![], None, None, vec![], true);
            let shell = result.add_tshell(vec![face]);
            result.add_tsolid(vec![shell]);
        }
    }

    result
}

// = =  Section Properties Computation = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Compute properties of a planar section.
///
/// Returns `None` if the section is not closed or not planar.

include!("extra.rs");
