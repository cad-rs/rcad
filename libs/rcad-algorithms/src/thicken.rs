//! Shell thickening — analogous to OCCT `BRepOffsetAPI_MakeThickSolid`.
//!
//! # Algorithm
//!
//! 1. Identify boundary wires (edges appearing in exactly one face).
//! 2. Offset each face along its normal by the given thickness.
//! 3. For each boundary edge, create a lateral ruled face connecting
//!    the original edge to the corresponding offset edge.
//! 4. Assemble offset faces + lateral faces into a closed solid.
//!
//! # Supported surfaces
//!
//! Plane, Sphere, Cylinder, Cone, Torus — each has a known parallel-surface
//! construction. B-spline and trimmed surfaces are skipped.

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

use crate::tolerance::TOLERANCE_ABS;
use crate::triangulate::{TessellationParams, mesh_brep};

/// Result of a thickening operation.
#[derive(Debug, Clone)]
pub struct ThickeningResult {
    /// The thickened solid as a new BRep.
    pub brep: BRep,
    /// Number of offset faces (one per input face).
    pub offset_faces: usize,
    /// Number of lateral faces connecting boundaries.
    pub lateral_faces: usize,
    /// Whether self-intersection was detected (thickness > half min face distance).
    pub self_intersection: bool,
}

// ── Inline BRep builder helpers (avoids rcad_modeling dependency) ────────────

fn add_vertex(brep: &mut BRep, point: DVec3) -> usize {
    let idx = brep.vertices.len();
    brep.vertices.push(Vertex { point });
    idx
}

fn add_edge(brep: &mut BRep, curve: Curve3, t0: f64, t1: f64, v0: usize, v1: usize) -> usize {
    let idx = brep.edges.len();
    brep.edges.push(Edge { start: v0, end: v1 });
    let ci = brep.geom.curves.len();
    brep.geom.curves.push(curve);
    while brep.geom.edge_curve.len() <= idx { brep.geom.edge_curve.push(None); }
    while brep.geom.edge_curve_range.len() <= idx { brep.geom.edge_curve_range.push(None); }
    while brep.geom.edge_degenerated.len() <= idx { brep.geom.edge_degenerated.push(false); }
    brep.geom.edge_curve[idx] = Some(ci);
    brep.geom.edge_curve_range[idx] = Some([t0, t1]);
    idx
}

fn add_face(brep: &mut BRep, surface: Surface3, outer: Wire, inner: Vec<Wire>) -> usize {
    if brep.solids.is_empty() {
        brep.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });
    }
    if brep.solids[0].shells.is_empty() {
        brep.solids[0].shells.push(Shell { faces: Vec::new() });
    }
    let idx = brep.solids[0].shells[0].faces.len();
    let normal = surface.normal_at(0.0, 0.0);
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: outer, inner_wires: inner, normal, triangles: Vec::new(),
    });
    while brep.geom.face_surface.len() <= idx { brep.geom.face_surface.push(None); }
    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surface);
    brep.geom.face_surface[idx] = Some(si);
    idx
}

// ── Surface offset ───────────────────────────────────────────────────────────

fn offset_surface(surf: &Surface3, d: f64) -> Option<Surface3> {
    use rcad_kernel::geom::*;
    match surf {
        Surface3::Plane(p) => Some(Surface3::Plane(Plane {
            origin: p.origin + p.normal * d,
            normal: p.normal,
        })),
        Surface3::Sphere(s) => {
            let r = s.radius + d;
            if r <= TOLERANCE_ABS { return None; }
            Some(Surface3::Sphere(SphericalSurface { center: s.center, axis: s.axis, radius: r }))
        }
        Surface3::Cylinder(c) => {
            let r = c.radius + d;
            if r <= TOLERANCE_ABS { return None; }
            Some(Surface3::Cylinder(CylindricalSurface { origin: c.origin, axis: c.axis, radius: r }))
        }
        Surface3::Cone(c) => {
            let sin_a = c.half_angle_rad.sin();
            let shift = if sin_a.abs() > 1e-10 { d / sin_a } else { d };
            let new_r = c.radius + d;
            if new_r <= TOLERANCE_ABS { return None; }
            Some(Surface3::Cone(ConicalSurface {
                apex: c.apex - c.axis * shift, axis: c.axis,
                radius: new_r, half_angle_rad: c.half_angle_rad,
            }))
        }
        Surface3::Torus(t) => {
            let r = t.minor_radius + d;
            if r <= TOLERANCE_ABS { return None; }
            Some(Surface3::Torus(ToroidalSurface {
                center: t.center, axis: t.axis,
                major_radius: t.major_radius, minor_radius: r,
            }))
        }
        _ => None,
    }
}

// ── Vertex normals ───────────────────────────────────────────────────────────

fn vertex_normal(shell: &Shell, brep: &BRep, vidx: usize) -> DVec3 {
    let mut n = DVec3::ZERO;
    let mut count = 0;
    for face in &shell.faces {
        let uses = face.outer_wire.edges.iter().any(|we| {
            let e = &brep.edges[we.idx];
            e.start == vidx || e.end == vidx
        });
        if uses { n += face.normal; count += 1; }
    }
    if count > 0 { (n / count as f64).normalize_or(DVec3::Z) } else { DVec3::Z }
}

// ── Edge chaining ────────────────────────────────────────────────────────────

fn chain_edges(edge_indices: &[usize], edges: &[Edge]) -> Vec<Vec<usize>> {
    if edge_indices.is_empty() { return vec![]; }
    let mut remaining: HashSet<usize> = edge_indices.iter().copied().collect();
    let mut loops = Vec::new();

    while let Some(&start_idx) = remaining.iter().next() {
        remaining.remove(&start_idx);
        let mut chain = vec![start_idx];
        let mut current_end = edges[start_idx].end;

        loop {
            let next = remaining.iter().find(|&&ei| {
                edges[ei].start == current_end || edges[ei].end == current_end
            }).copied();
            match next {
                Some(ei) => {
                    remaining.remove(&ei);
                    chain.push(ei);
                    let e = &edges[ei];
                    current_end = if e.start == current_end { e.end } else { e.start };
                }
                None => break,
            }
        }
        if chain.len() >= 2 { loops.push(chain); }
    }
    loops
}

// ── Public API ───────────────────────────────────────────────────────────────

/// Thicken a solid by removing specified faces, offsetting the remaining
/// faces, and building lateral ruled faces at the removed-face boundaries.
///
/// This is analogous to OCCT `BRepOffsetAPI_MakeThickSolid`.
///
/// - `brep`: input solid (must have at least one shell with geometry).
/// - `removed_face_indices`: indices of faces to remove (relative to
///   `brep.solids[0].shells[0].faces`).
/// - `thickness`: positive = inward (material removed), negative = outward.
///
/// Returns `None` if all faces are removed, thickness is zero, or the offset
/// would create degenerate surfaces.
pub fn thick_solid_with_removed_faces(
    brep: &BRep,
    removed_face_indices: &[usize],
    thickness: f64,
) -> Option<ThickeningResult> {
    if thickness.abs() < 1e-12 {
        return None;
    }

    let shell = brep.solids.first()?.shells.first()?;
    if shell.faces.is_empty() {
        return None;
    }

    let removed_set: HashSet<usize> = removed_face_indices.iter().copied().collect();
    if removed_set.len() >= shell.faces.len() {
        return None; // can't remove all faces
    }

    let d = thickness;

    // ── Step 1: build the "kept" shell (original faces minus removed) ──
    let kept_faces: Vec<(usize, &Face)> = shell
        .faces
        .iter()
        .enumerate()
        .filter(|(i, _)| !removed_set.contains(i))
        .collect();

    if kept_faces.is_empty() {
        return None;
    }

    // ── Step 2: find boundary edges of the kept shell ──────────────────
    // An edge is on the boundary if it appears in exactly one kept face.
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for (_, face) in &kept_faces {
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }
    let boundary_edges: Vec<usize> = edge_use
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(idx, _)| idx)
        .collect();

    // ── Step 3: compute offset vertex positions ────────────────────────
    // Use only kept faces for vertex normal computation.
    let kept_shell = Shell {
        faces: kept_faces.iter().map(|(_, f)| (*f).clone()).collect(),
    };
    let new_pts: Vec<DVec3> = brep.vertices.iter().enumerate().map(|(i, _)| {
        let n = vertex_normal(&kept_shell, brep, i);
        brep.vertices[i].point + n * d
    }).collect();

    // ── Step 4: build result BRep ──────────────────────────────────────
    let mut out = BRep::new();
    out.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    let mut orig_vidx: Vec<usize> = Vec::new();
    for v in &brep.vertices {
        orig_vidx.push(add_vertex(&mut out, v.point));
    }
    let mut off_vidx: Vec<usize> = Vec::new();
    for &p in &new_pts {
        off_vidx.push(add_vertex(&mut out, p));
    }

    // ── Step 5: offset kept faces ──────────────────────────────────────
    let mut offset_face_count = 0;
    for &(fi, face) in &kept_faces {
        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|o| *o) {
            Some(s) => s,
            None => continue,
        };
        let surf = &brep.geom.surfaces[surf_idx];
        let off_surf = match offset_surface(surf, d) {
            Some(s) => s,
            None => continue,
        };

        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = off_vidx[e.start];
            let ve = off_vidx[e.end];
            let dir = (out.vertices[ve].point - out.vertices[vs].point).normalize_or(DVec3::X);
            let len = (out.vertices[ve].point - out.vertices[vs].point).length();
            let curve = Curve3::Line(Line3 {
                origin: out.vertices[vs].point,
                direction: dir,
            });
            let eidx = add_edge(&mut out, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut out, off_surf, Wire { edges: wire_edges }, Vec::new());
        offset_face_count += 1;
    }

    if offset_face_count == 0 {
        return None;
    }

    // ── Step 6: lateral faces along boundary edges ─────────────────────
    let mut lateral_count = 0;
    let loops = chain_edges(&boundary_edges, &brep.edges);

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let e = &brep.edges[eidx];
            let o_vs = orig_vidx[e.start];
            let o_ve = orig_vidx[e.end];
            let f_vs = off_vidx[e.start];
            let f_ve = off_vidx[e.end];

            let p0 = out.vertices[o_vs].point;
            let p1 = out.vertices[o_ve].point;
            let p3 = out.vertices[f_vs].point;

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < 1e-10 {
                continue;
            }

            let surf = Surface3::Plane(rcad_kernel::geom::Plane {
                origin: p0,
                normal,
            });

            let vseq = [o_vs, o_ve, f_ve, f_vs];
            let mut edges = Vec::new();
            for i in 0..4 {
                let s = vseq[i];
                let en = vseq[(i + 1) % 4];
                let dir =
                    (out.vertices[en].point - out.vertices[s].point).normalize_or(DVec3::X);
                let len = (out.vertices[en].point - out.vertices[s].point).length();
                let curve = Curve3::Line(Line3 {
                    origin: out.vertices[s].point,
                    direction: dir,
                });
                edges.push(WireEdge::fwd(add_edge(&mut out, curve, 0.0, len, s, en)));
            }

            add_face(&mut out, surf, Wire { edges }, Vec::new());
            lateral_count += 1;
        }
    }

    // ── Step 7: triangulate ────────────────────────────────────────────
    mesh_brep(&mut out, &TessellationParams::default());

    // ── Step 8: self-intersection detection ────────────────────────────
    // For closed shells (no boundary edges), check if thickness exceeds
    // half the minimum distance between non-adjacent faces.
    let self_intersection = if boundary_edges.is_empty() && removed_face_indices.is_empty() {
        detect_self_intersection(brep, thickness)
    } else {
        false
    };

    Some(ThickeningResult {
        brep: out,
        offset_faces: offset_face_count,
        lateral_faces: lateral_count,
        self_intersection,
    })
}

/// Detect self-intersection for closed-shell inward offsetting.
///
/// Computes the minimum distance between non-adjacent face centroids.
/// If `thickness > min_distance / 2`, the offset faces will self-intersect.
fn detect_self_intersection(brep: &BRep, thickness: f64) -> bool {
    let shell = brep.solids.first().and_then(|s| s.shells.first());
    let shell = match shell {
        Some(s) => s,
        None => return false,
    };

    // Compute face centroids
    let centroids: Vec<DVec3> = shell.faces.iter().map(|face| {
        let mut sum = DVec3::ZERO;
        let mut count = 0;
        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            sum += brep.vertices[e.start].point;
            count += 1;
        }
        if count > 0 { sum / count as f64 } else { DVec3::ZERO }
    }).collect();

    // Find minimum distance between non-adjacent faces
    let mut min_dist = f64::MAX;
    for i in 0..centroids.len() {
        for j in (i + 1)..centroids.len() {
            // Check if faces share an edge (adjacent)
            let share_edge = shell.faces[i].outer_wire.edges.iter()
                .any(|we_i| shell.faces[j].outer_wire.edges.iter().any(|we_j| we_i.idx == we_j.idx));
            if share_edge {
                continue;
            }
            let dist = (centroids[i] - centroids[j]).length();
            if dist < min_dist {
                min_dist = dist;
            }
        }
    }

    if min_dist == f64::MAX {
        return false; // no non-adjacent faces
    }

    thickness.abs() > min_dist * 0.5
}

/// Thicken an open shell by offsetting faces along their normals and
/// filling the gaps with lateral ruled faces.
///
/// The input BRep must have at least one face with populated surface data
/// (e.g. created via `make_box_brep` which populates analytic surfaces).
///
/// `thickness` > 0 offsets outward, < 0 offsets inward.
/// Returns `None` if the shell is closed, has no geometry, or the offset
/// would create degenerate surfaces.
pub fn thicken_shell(brep: &BRep, thickness: f64) -> Option<ThickeningResult> {
    if thickness.abs() < 1e-12 { return None; }

    let shell = brep.solids.first()?.shells.first()?;
    if shell.faces.is_empty() { return None; }

    let d = thickness;

    // ── Step 1: find boundary edges ──────────────────────────────────────
    let mut edge_use: HashMap<usize, usize> = HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_use.entry(we.idx).or_insert(0) += 1;
        }
    }
    let boundary_edges: Vec<usize> = edge_use.into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(idx, _)| idx)
        .collect();

    // ── Step 2: compute offset vertex positions ──────────────────────────
    let new_pts: Vec<DVec3> = brep.vertices.iter().enumerate().map(|(i, _)| {
        let n = vertex_normal(shell, brep, i);
        brep.vertices[i].point + n * d
    }).collect();

    // ── Step 3: build result BRep with original + offset vertices ────────
    let mut out = BRep::new();
    out.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    let mut orig_vidx: Vec<usize> = Vec::new();
    for v in &brep.vertices {
        orig_vidx.push(add_vertex(&mut out, v.point));
    }
    let mut off_vidx: Vec<usize> = Vec::new();
    for &p in &new_pts {
        off_vidx.push(add_vertex(&mut out, p));
    }

    // ── Step 4: offset faces ─────────────────────────────────────────────
    let mut offset_face_count = 0;
    for (fi, face) in shell.faces.iter().enumerate() {
        let surf_idx = match brep.geom.face_surface.get(fi).and_then(|o| *o) {
            Some(s) => s, None => continue,
        };
        let surf = &brep.geom.surfaces[surf_idx];
        let off_surf = match offset_surface(surf, d) {
            Some(s) => s, None => continue,
        };

        // Build wire from offset vertices
        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let e = &brep.edges[we.idx];
            let vs = off_vidx[e.start];
            let ve = off_vidx[e.end];
            let dir = (out.vertices[ve].point - out.vertices[vs].point).normalize_or(DVec3::X);
            let len = (out.vertices[ve].point - out.vertices[vs].point).length();
            let curve = Curve3::Line(Line3 { origin: out.vertices[vs].point, direction: dir });
            let eidx = add_edge(&mut out, curve, 0.0, len, vs, ve);
            wire_edges.push(WireEdge::fwd(eidx));
        }

        add_face(&mut out, off_surf, Wire { edges: wire_edges }, Vec::new());
        offset_face_count += 1;
    }

    if offset_face_count == 0 { return None; }

    // ── Step 5: lateral faces along boundary edges ───────────────────────
    let mut lateral_count = 0;
    let loops = chain_edges(&boundary_edges, &brep.edges);

    for loop_edges in &loops {
        for &eidx in loop_edges {
            let e = &brep.edges[eidx];
            let o_vs = orig_vidx[e.start];
            let o_ve = orig_vidx[e.end];
            let f_vs = off_vidx[e.start];
            let f_ve = off_vidx[e.end];

            let p0 = out.vertices[o_vs].point;
            let p1 = out.vertices[o_ve].point;
            let p3 = out.vertices[f_vs].point;

            let normal = (p1 - p0).cross(p3 - p0).normalize_or(DVec3::Z);
            if normal.length() < 1e-10 { continue; }

            let surf = Surface3::Plane(rcad_kernel::geom::Plane { origin: p0, normal });

            // Quad: orig_start → orig_end → off_end → off_start
            let vseq = [o_vs, o_ve, f_ve, f_vs];
            let mut edges = Vec::new();
            for i in 0..4 {
                let s = vseq[i];
                let en = vseq[(i + 1) % 4];
                let dir = (out.vertices[en].point - out.vertices[s].point).normalize_or(DVec3::X);
                let len = (out.vertices[en].point - out.vertices[s].point).length();
                let curve = Curve3::Line(Line3 { origin: out.vertices[s].point, direction: dir });
                edges.push(WireEdge::fwd(add_edge(&mut out, curve, 0.0, len, s, en)));
            }

            add_face(&mut out, surf, Wire { edges }, Vec::new());
            lateral_count += 1;
        }
    }

    // ── Step 6: triangulate ──────────────────────────────────────────────
    mesh_brep(&mut out, &TessellationParams::default());

    Some(ThickeningResult { brep: out, offset_faces: offset_face_count, lateral_faces: lateral_count, self_intersection: false })
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_modeling::make_box_brep;

    #[test]
    fn offset_plane_translates() {
        let plane = Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO, normal: DVec3::Z,
        });
        let off = offset_surface(&plane, 0.5).unwrap();
        if let Surface3::Plane(p) = off {
            assert!((p.origin.z - 0.5).abs() < 1e-9);
        } else { panic!("expected Plane"); }
    }

    #[test]
    fn offset_sphere_grows() {
        let s = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: DVec3::ZERO, axis: DVec3::Z, radius: 2.0,
        });
        let off = offset_surface(&s, 0.5).unwrap();
        if let Surface3::Sphere(s) = off {
            assert!((s.radius - 2.5).abs() < 1e-9);
        } else { panic!("expected Sphere"); }
    }

    #[test]
    fn offset_cylinder_grows() {
        let c = Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: DVec3::ZERO, axis: DVec3::Z, radius: 1.0,
        });
        let off = offset_surface(&c, 0.3).unwrap();
        if let Surface3::Cylinder(c) = off {
            assert!((c.radius - 1.3).abs() < 1e-9);
        } else { panic!("expected Cylinder"); }
    }

    #[test]
    fn thicken_closed_box_no_lateral_faces() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);
        eprintln!("box faces={}, verts={}, edges={}", box_brep.solids.first().map(|s| s.shells.first().map(|sh| sh.faces.len()).unwrap_or(0)).unwrap_or(0), box_brep.vertices.len(), box_brep.edges.len());
        eprintln!("geom surfaces={}, face_surface={:?}", box_brep.geom.surfaces.len(), box_brep.geom.face_surface);
        let result = thicken_shell(&box_brep, 0.1);
        // Closed shell → 6 offset faces, no lateral faces
        let r = result.expect("closed shell should still offset faces");
        assert_eq!(r.offset_faces, 6);
        assert_eq!(r.lateral_faces, 0);
    }

    #[test]
    fn thicken_open_box_produces_lateral_faces() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        // Remove top face → open shell
        let mut open_brep = box_brep.clone();
        if let Some(s) = open_brep.solids.first_mut() {
            if let Some(sh) = s.shells.first_mut() {
                if sh.faces.len() > 1 { sh.faces.pop(); }
            }
        }

        let result = thicken_shell(&open_brep, 0.1);
        assert!(result.is_some(), "open shell thickening should succeed");
        let r = result.unwrap();
        assert_eq!(r.offset_faces, 5, "should offset 5 faces");
        assert!(r.lateral_faces > 0, "should create lateral faces for the open boundary");
    }

    #[test]
    fn thicken_negative_thickness_inwards() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);
        let mut open_brep = box_brep.clone();
        if let Some(s) = open_brep.solids.first_mut() {
            if let Some(sh) = s.shells.first_mut() {
                if sh.faces.len() > 1 { sh.faces.pop(); }
            }
        }

        let result = thicken_shell(&open_brep, -0.1);
        assert!(result.is_some(), "negative thickness should work (inward offset)");
    }

    #[test]
    fn thicken_zero_returns_none() {
        let box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        assert!(thicken_shell(&box_brep, 0.0).is_none());
    }

    #[test]
    fn thick_solid_remove_one_face() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        // Remove top face (index 5) and thicken inward
        let result = thick_solid_with_removed_faces(&box_brep, &[5], 0.1);
        assert!(result.is_some(), "should succeed with one face removed");
        let r = result.unwrap();
        assert_eq!(r.offset_faces, 5, "should offset 5 kept faces");
        assert!(r.lateral_faces > 0, "should create lateral faces at the removed-face boundary");
    }

    #[test]
    fn thick_solid_remove_multiple_faces() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        // Remove top (5) and bottom (0) faces
        let result = thick_solid_with_removed_faces(&box_brep, &[0, 5], 0.1);
        assert!(result.is_some(), "should succeed with two faces removed");
        let r = result.unwrap();
        assert_eq!(r.offset_faces, 4, "should offset 4 kept faces");
        assert!(r.lateral_faces > 0, "should create lateral faces at both removed boundaries");
    }

    #[test]
    fn thick_solid_remove_all_faces_returns_none() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        let result = thick_solid_with_removed_faces(&box_brep, &[0, 1, 2, 3, 4, 5], 0.1);
        assert!(result.is_none(), "removing all faces should return None");
    }

    #[test]
    fn thick_solid_zero_thickness_returns_none() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        let result = thick_solid_with_removed_faces(&box_brep, &[5], 0.0);
        assert!(result.is_none(), "zero thickness should return None");
    }

    #[test]
    fn thick_solid_closed_box_detects_self_intersection() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        // A 1x1x1 box has minimum face-to-face distance of 1.0.
        // Inward offset with thickness > 0.5 would self-intersect.
        // The function should still produce a result but warn about self-intersection.
        let result = thick_solid_with_removed_faces(&box_brep, &[], 0.6);
        assert!(result.is_some(), "should produce a result even with self-intersection");
        let r = result.unwrap();
        assert!(
            r.self_intersection,
            "should detect self-intersection for thickness > half min dimension"
        );
    }

    #[test]
    fn thick_solid_no_self_intersection_small_thickness() {
        let mut box_brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut box_brep);

        // A 2x2x2 box: min face-to-face distance is 2.0, so thickness 0.5 is safe.
        let result = thick_solid_with_removed_faces(&box_brep, &[], 0.5);
        assert!(result.is_some());
        let r = result.unwrap();
        assert!(
            !r.self_intersection,
            "should not self-intersect for small thickness"
        );
    }
}
