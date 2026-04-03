//! Chamfer and fillet operations on convex BRep edges.
//!
//! Both operations are limited to:
//! - Convex edges shared by exactly two planar faces
//! - The first solid's first shell (`solids[0].shells[0]`)
//! - One edge at a time (no corner blending)
//!
//! Both functions return a new BRep rather than modifying in place.

use glam::DVec3;
use rcad_kernel::geom::{CylindricalSurface, Line3, Plane, Surface3};
use rcad_kernel::topology::{Face, Vertex, WireEdge};
use rcad_kernel::BRep;

use crate::builder::brep_builder::{make_edge, make_face, make_wire};
use crate::builder::BuildError;

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Find the two faces in `solids[0].shells[0]` that share `edge_idx`.
/// Returns `None` if not exactly two faces reference the edge.
fn find_adjacent_faces(brep: &BRep, edge_idx: usize) -> Option<(usize, usize)> {
    let shell = brep.solids.first()?.shells.first()?;
    let mut found = Vec::new();
    for (fi, face) in shell.faces.iter().enumerate() {
        if face.outer_wire.edges.iter().any(|we| we.idx == edge_idx) {
            found.push(fi);
        }
        if found.len() > 2 {
            return None;
        }
    }
    if found.len() == 2 {
        Some((found[0], found[1]))
    } else {
        None
    }
}

/// Compute the inward direction from `edge_idx` into face `face_idx`.
///
/// Finds a vertex in the face that is NOT an endpoint of the edge, then
/// projects the vector from an edge endpoint to that vertex onto the plane
/// perpendicular to the edge direction. Returns the normalized result.
fn setback_direction(brep: &BRep, face_idx: usize, edge_idx: usize) -> DVec3 {
    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;
    let edge_dir = (p1 - p0).normalize_or_zero();

    let face = &brep.solids[0].shells[0].faces[face_idx];
    // Collect all vertex indices referenced by this face's outer wire.
    let edge_verts = [edge.start, edge.end];
    for we in &face.outer_wire.edges {
        if let Some(e) = brep.edges.get(we.idx) {
            for &vi in &[e.start, e.end] {
                if !edge_verts.contains(&vi) {
                    let v_other = brep.vertices[vi].point;
                    let diff = v_other - p0;
                    let along = edge_dir * edge_dir.dot(diff);
                    let perp = diff - along;
                    let n = perp.normalize_or_zero();
                    if n.length_squared() > 1e-20 {
                        return n;
                    }
                }
            }
        }
    }
    DVec3::ZERO
}

/// Return the face normal for face `face_idx` in `solids[0].shells[0]`.
/// Falls back to the stored `face.normal`.
fn face_normal(brep: &BRep, face_idx: usize) -> DVec3 {
    brep.solids[0].shells[0].faces[face_idx].normal
}

// ── Rebuild helpers ───────────────────────────────────────────────────────────

/// Get the start vertex of a WireEdge (respecting orientation).
#[allow(dead_code)]
fn wire_edge_start(brep: &BRep, we: &WireEdge) -> usize {
    let e = &brep.edges[we.idx];
    if we.forward { e.start } else { e.end }
}

/// Get the end vertex of a WireEdge (respecting orientation).
#[allow(dead_code)]
fn wire_edge_end(brep: &BRep, we: &WireEdge) -> usize {
    let e = &brep.edges[we.idx];
    if we.forward { e.end } else { e.start }
}

/// Add a straight line edge between two vertices in `dst`.
fn add_line_edge(dst: &mut BRep, va: usize, vb: usize) -> Result<usize, BuildError> {
    let pa = dst.vertices[va].point;
    let pb = dst.vertices[vb].point;
    let dir = (pb - pa).normalize_or_zero();
    let len = (pb - pa).length();
    let curve = rcad_kernel::geom::Curve3::Line(Line3 { origin: pa, direction: dir });
    make_edge(dst, curve, 0.0, len, va, vb)
}

/// Add a triangular closing face with the given three vertex indices.
fn add_closing_triangle(dst: &mut BRep, v0: usize, v1: usize, v2: usize) -> Result<(), BuildError> {
    let p0 = dst.vertices[v0].point;
    let p1 = dst.vertices[v1].point;
    let p2 = dst.vertices[v2].point;
    let n = (p1 - p0).cross(p2 - p0).normalize_or_zero();
    let surf = Surface3::Plane(Plane { origin: p0, normal: n });

    let e01 = add_line_edge(dst, v0, v1)?;
    let e12 = add_line_edge(dst, v1, v2)?;
    let e20 = add_line_edge(dst, v2, v0)?;
    let wire = make_wire(vec![WireEdge::fwd(e01), WireEdge::fwd(e12), WireEdge::fwd(e20)]);
    make_face(dst, surf, wire, vec![])?;
    Ok(())
}

/// Copy all vertices and geom data from `src` into `dst` (start of rebuild).
/// After this call, vertex indices in `dst` match those in `src`.
fn copy_vertices_from(dst: &mut BRep, src: &BRep) {
    for v in &src.vertices {
        dst.vertices.push(Vertex { point: v.point });
    }
}

/// Push a vertex into `dst` and return its index.
fn push_vertex(dst: &mut BRep, point: DVec3) -> usize {
    let idx = dst.vertices.len();
    dst.vertices.push(Vertex { point });
    idx
}

/// Rebuild one face from `src` into `dst`, applying vertex remapping.
///
/// `vi_remap[old_vi]` gives the new vertex index for `old_vi`. If `vi_remap[i] == i`
/// no remapping occurs (identity). Edges are recreated as straight lines.
fn copy_face_remapped(
    dst: &mut BRep,
    src: &BRep,
    face: &Face,
    vi_remap: &[usize],
) -> Result<(), BuildError> {
    // Build a new wire by walking the outer_wire and remapping vertices.
    let mut new_wire_edges = Vec::new();

    for we in &face.outer_wire.edges {
        let src_edge = &src.edges[we.idx];
        let old_va = if we.forward { src_edge.start } else { src_edge.end };
        let old_vb = if we.forward { src_edge.end } else { src_edge.start };
        let new_va = vi_remap[old_va];
        let new_vb = vi_remap[old_vb];

        // Add straight-line edge in dst
        let ei = add_line_edge(dst, new_va, new_vb)?;
        new_wire_edges.push(WireEdge::fwd(ei));
    }

    let wire = make_wire(new_wire_edges);

    // Compute face normal from first 3 vertices of the wire.
    let normal = {
        let pts: Vec<DVec3> = face.outer_wire.edges.iter()
            .take(3)
            .filter_map(|we| {
                let src_e = &src.edges[we.idx];
                let old_v = if we.forward { src_e.start } else { src_e.end };
                Some(dst.vertices[vi_remap[old_v]].point)
            })
            .collect();
        if pts.len() >= 3 {
            (pts[1] - pts[0]).cross(pts[2] - pts[0]).normalize_or_zero()
        } else {
            face.normal
        }
    };

    let surf = Surface3::Plane(Plane { origin: dst.vertices[vi_remap[
        src.edges[face.outer_wire.edges[0].idx].start
    ]].point, normal });

    make_face(dst, surf, wire, vec![])?;
    Ok(())
}

// ── Chamfer ───────────────────────────────────────────────────────────────────

/// Chamfer (straight bevel) the edge at `edge_idx` in `brep` by `dist` on each side.
///
/// Produces a new BRep with:
/// - The two adjacent faces shortened by `dist` on each side of the edge
/// - A new planar chamfer face replacing the edge
/// - Two triangular closing faces at the edge endpoints
///
/// # Errors
/// Returns `BuildError::InvalidIndex` if `edge_idx` is out of bounds.
/// Returns `BuildError::DegenerateGeometry` if the edge is not shared by exactly two faces,
/// or if the setback direction cannot be computed.
/// Returns `BuildError::NonPositiveValue` if `dist <= 0`.
pub fn chamfer_edge(brep: &BRep, edge_idx: usize, dist: f64) -> Result<BRep, BuildError> {
    if dist <= 0.0 {
        return Err(BuildError::NonPositiveValue("dist"));
    }
    if edge_idx >= brep.edges.len() {
        return Err(BuildError::InvalidIndex(edge_idx));
    }

    let (f0, f1) = find_adjacent_faces(brep, edge_idx)
        .ok_or(BuildError::DegenerateGeometry("edge must be shared by exactly 2 faces"))?;

    let s0 = setback_direction(brep, f0, edge_idx);
    let s1 = setback_direction(brep, f1, edge_idx);
    if s0.length_squared() < 1e-20 || s1.length_squared() < 1e-20 {
        return Err(BuildError::DegenerateGeometry("cannot compute setback direction"));
    }

    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    // 4 new setback vertices
    let nv0a = p0 + dist * s0; // face f0, at v0
    let nv1a = p1 + dist * s0; // face f0, at v1
    let nv0b = p0 + dist * s1; // face f1, at v0
    let nv1b = p1 + dist * s1; // face f1, at v1

    rebuild_with_chamfer_verts(brep, edge_idx, f0, f1, nv0a, nv1a, nv0b, nv1b)
}

// ── Fillet ────────────────────────────────────────────────────────────────────

/// Fillet (rounded) the edge at `edge_idx` in `brep` with `radius`.
///
/// Produces a new BRep with:
/// - The two adjacent faces shortened by `setback = radius / tan(beta/2)` on each side
/// - A new cylindrical fillet face replacing the edge
/// - Two triangular closing faces at the edge endpoints
///
/// # Errors
/// Returns `BuildError::InvalidIndex` if `edge_idx` is out of bounds.
/// Returns `BuildError::DegenerateGeometry` if the edge is not shared by exactly two faces.
/// Returns `BuildError::NonPositiveValue` if `radius <= 0`.
pub fn fillet_edge(brep: &BRep, edge_idx: usize, radius: f64) -> Result<BRep, BuildError> {
    if radius <= 0.0 {
        return Err(BuildError::NonPositiveValue("radius"));
    }
    if edge_idx >= brep.edges.len() {
        return Err(BuildError::InvalidIndex(edge_idx));
    }

    let (f0, f1) = find_adjacent_faces(brep, edge_idx)
        .ok_or(BuildError::DegenerateGeometry("edge must be shared by exactly 2 faces"))?;

    let s0 = setback_direction(brep, f0, edge_idx);
    let s1 = setback_direction(brep, f1, edge_idx);
    if s0.length_squared() < 1e-20 || s1.length_squared() < 1e-20 {
        return Err(BuildError::DegenerateGeometry("cannot compute setback direction"));
    }

    // Compute dihedral angle from face normals (exterior angle beta).
    // For a convex edge, n0·n1 = cos(π - beta) → beta = π - acos(n0·n1)
    let n0 = face_normal(brep, f0);
    let n1 = face_normal(brep, f1);
    let cos_angle = n0.dot(n1).clamp(-1.0, 1.0);
    // exterior dihedral angle between the faces
    let beta = (std::f64::consts::PI - cos_angle.acos()).abs();
    let half_beta = beta * 0.5;
    let setback = if half_beta.tan().abs() < 1e-10 {
        radius // fallback for 180° (flat) case
    } else {
        radius / half_beta.tan()
    };

    let edge = &brep.edges[edge_idx];
    let p0 = brep.vertices[edge.start].point;
    let p1 = brep.vertices[edge.end].point;

    let nv0a = p0 + setback * s0;
    let nv1a = p1 + setback * s0;
    let nv0b = p0 + setback * s1;
    let nv1b = p1 + setback * s1;

    rebuild_with_fillet_verts(brep, edge_idx, f0, f1, nv0a, nv1a, nv0b, nv1b, radius)
}

// ── Shared rebuild core ───────────────────────────────────────────────────────

/// Rebuild the BRep replacing `edge_idx` with a chamfer quad face.
fn rebuild_with_chamfer_verts(
    brep: &BRep,
    edge_idx: usize,
    f0: usize,
    f1: usize,
    nv0a: DVec3,
    nv1a: DVec3,
    nv0b: DVec3,
    nv1b: DVec3,
) -> Result<BRep, BuildError> {
    let orig_edge = &brep.edges[edge_idx];
    let v0_orig = orig_edge.start;
    let v1_orig = orig_edge.end;

    let mut dst = BRep::new();
    copy_vertices_from(&mut dst, brep);

    // New vertices
    let nv0a_idx = push_vertex(&mut dst, nv0a);
    let nv1a_idx = push_vertex(&mut dst, nv1a);
    let nv0b_idx = push_vertex(&mut dst, nv0b);
    let nv1b_idx = push_vertex(&mut dst, nv1b);

    let shell = &brep.solids[0].shells[0];
    let n_faces = shell.faces.len();

    for fi in 0..n_faces {
        let face = &shell.faces[fi];
        if fi == f0 {
            let remap = build_remap(brep, face, v0_orig, v1_orig, nv0a_idx, nv1a_idx);
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        } else if fi == f1 {
            let remap = build_remap(brep, face, v0_orig, v1_orig, nv0b_idx, nv1b_idx);
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        } else {
            // Copy unchanged — identity remap
            let remap: Vec<usize> = (0..brep.vertices.len()).collect();
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        }
    }

    // Chamfer face: nv0a → nv1a → nv1b → nv0b (quad)
    {
        let pa = nv0a;
        let pb = nv1a;
        let pc = nv1b;
        let _pd = nv0b;
        let n = (pb - pa).cross(pc - pa).normalize_or_zero();
        let surf = Surface3::Plane(Plane { origin: pa, normal: n });

        let ea = add_line_edge(&mut dst, nv0a_idx, nv1a_idx)?;
        let eb = add_line_edge(&mut dst, nv1a_idx, nv1b_idx)?;
        let ec = add_line_edge(&mut dst, nv1b_idx, nv0b_idx)?;
        let ed = add_line_edge(&mut dst, nv0b_idx, nv0a_idx)?;
        let wire = make_wire(vec![WireEdge::fwd(ea), WireEdge::fwd(eb), WireEdge::fwd(ec), WireEdge::fwd(ed)]);
        make_face(&mut dst, surf, wire, vec![])?;
    }

    // Closing triangles at v0 and v1
    add_closing_triangle(&mut dst, v0_orig, nv0a_idx, nv0b_idx)?;
    add_closing_triangle(&mut dst, v1_orig, nv1b_idx, nv1a_idx)?;

    Ok(dst)
}

/// Rebuild the BRep replacing `edge_idx` with a cylindrical fillet face.
fn rebuild_with_fillet_verts(
    brep: &BRep,
    edge_idx: usize,
    f0: usize,
    f1: usize,
    nv0a: DVec3,
    nv1a: DVec3,
    nv0b: DVec3,
    nv1b: DVec3,
    radius: f64,
) -> Result<BRep, BuildError> {
    let orig_edge = &brep.edges[edge_idx];
    let v0_orig = orig_edge.start;
    let v1_orig = orig_edge.end;

    let mut dst = BRep::new();
    copy_vertices_from(&mut dst, brep);

    let nv0a_idx = push_vertex(&mut dst, nv0a);
    let nv1a_idx = push_vertex(&mut dst, nv1a);
    let nv0b_idx = push_vertex(&mut dst, nv0b);
    let nv1b_idx = push_vertex(&mut dst, nv1b);

    let shell = &brep.solids[0].shells[0];
    let n_faces = shell.faces.len();

    for fi in 0..n_faces {
        let face = &shell.faces[fi];
        if fi == f0 {
            let remap = build_remap(brep, face, v0_orig, v1_orig, nv0a_idx, nv1a_idx);
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        } else if fi == f1 {
            let remap = build_remap(brep, face, v0_orig, v1_orig, nv0b_idx, nv1b_idx);
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        } else {
            let remap: Vec<usize> = (0..brep.vertices.len()).collect();
            copy_face_remapped(&mut dst, brep, face, &remap)?;
        }
    }

    // Fillet face: cylindrical surface along the original edge direction
    {
        let p_orig0 = brep.vertices[v0_orig].point;
        let p_orig1 = brep.vertices[v1_orig].point;
        let edge_dir = (p_orig1 - p_orig0).normalize_or_zero();

        // Use CylindricalSurface as the fillet face geometry
        let surf = Surface3::Cylinder(CylindricalSurface {
            origin: p_orig0,
            axis: edge_dir,
            radius,
        });

        let ea = add_line_edge(&mut dst, nv0a_idx, nv1a_idx)?;
        let eb = add_line_edge(&mut dst, nv1a_idx, nv1b_idx)?;
        let ec = add_line_edge(&mut dst, nv1b_idx, nv0b_idx)?;
        let ed = add_line_edge(&mut dst, nv0b_idx, nv0a_idx)?;
        let wire = make_wire(vec![WireEdge::fwd(ea), WireEdge::fwd(eb), WireEdge::fwd(ec), WireEdge::fwd(ed)]);
        make_face(&mut dst, surf, wire, vec![])?;
    }

    // Closing triangles
    add_closing_triangle(&mut dst, v0_orig, nv0a_idx, nv0b_idx)?;
    add_closing_triangle(&mut dst, v1_orig, nv1b_idx, nv1a_idx)?;

    Ok(dst)
}

/// Build a vertex remap array that replaces `v0_orig` with `new_v0` and
/// `v1_orig` with `new_v1` everywhere in a face.
fn build_remap(
    brep: &BRep,
    face: &Face,
    v0_orig: usize,
    v1_orig: usize,
    new_v0: usize,
    new_v1: usize,
) -> Vec<usize> {
    let mut remap: Vec<usize> = (0..brep.vertices.len()).collect();
    // Check if this face actually references v0_orig / v1_orig
    let face_verts: Vec<usize> = face.outer_wire.edges.iter()
        .filter_map(|we| brep.edges.get(we.idx))
        .flat_map(|e| [e.start, e.end])
        .collect();
    if face_verts.contains(&v0_orig) {
        remap[v0_orig] = new_v0;
    }
    if face_verts.contains(&v1_orig) {
        remap[v1_orig] = new_v1;
    }
    remap
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::PrimitiveSolid;

    fn box_brep_2x2x2() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 })
    }

    #[test]
    fn chamfer_edge_produces_more_faces() {
        let brep = box_brep_2x2x2();
        let result = chamfer_edge(&brep, 0, 0.2).unwrap();
        let n_faces = result.solids[0].shells[0].faces.len();
        // 6 original - 2 trimmed + 1 chamfer + 2 closing triangles = 7 total?
        // Actually: original 6 faces all kept (2 trimmed in-place), +1 chamfer +2 closing = 9
        // But our rebuild copies all 6 faces (2 with remap, 4 unchanged), + 1 chamfer + 2 closing = 9
        assert_eq!(n_faces, 9, "expected 9 faces after chamfer (6 + 1 chamfer + 2 closing)");
    }

    #[test]
    fn fillet_edge_produces_more_faces() {
        let brep = box_brep_2x2x2();
        let result = fillet_edge(&brep, 0, 0.2).unwrap();
        let n_faces = result.solids[0].shells[0].faces.len();
        assert_eq!(n_faces, 9, "expected 9 faces after fillet (6 + 1 fillet + 2 closing)");
    }

    #[test]
    fn chamfer_rejects_zero_dist() {
        let brep = box_brep_2x2x2();
        assert!(chamfer_edge(&brep, 0, 0.0).is_err());
    }

    #[test]
    fn fillet_rejects_zero_radius() {
        let brep = box_brep_2x2x2();
        assert!(fillet_edge(&brep, 0, 0.0).is_err());
    }

    #[test]
    fn chamfer_invalid_edge_index() {
        let brep = box_brep_2x2x2();
        assert!(chamfer_edge(&brep, 999, 0.2).is_err());
    }
}
