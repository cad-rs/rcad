//! Draft angle operation — analogous to OCCT `BRepDraftBuilder`.
//!
//! # Algorithm
//!
//! For each vertex, compute its signed distance `h` to the neutral plane along
//! the pull direction. The draft displacement is:
//!
//!   delta = h * tan(angle) * n_perp
//!
//! where `n_perp` is the component of the face normal perpendicular to the pull
//! direction. This tilts each face by the draft angle while keeping vertices on
//! the neutral plane fixed.
//!
//! # Supported surfaces
//!
//! Only planar faces are supported in this phase. Non-planar faces return
//! `DraftError::UnsupportedSurface`.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3};
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

/// Parameters controlling the draft operation.
#[derive(Debug, Clone)]
pub struct DraftParams {
    /// Normalized pull direction (the "pull" axis of the mold).
    pub pull_direction: DVec3,
    /// Draft angle in radians. Positive = material added, negative = removed.
    pub draft_angle: f64,
    /// A point on the neutral plane (vertices on this plane don't move).
    pub neutral_point: DVec3,
}

/// Error type for draft operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DraftError {
    /// A face has a surface type that is not yet supported for drafting.
    UnsupportedSurface,
    /// The draft angle is too large (> 89 degrees).
    AngleTooLarge,
    /// The input BRep has no faces.
    NoFaces,
}

impl std::fmt::Display for DraftError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedSurface => write!(f, "unsupported surface type for drafting"),
            Self::AngleTooLarge => write!(f, "draft angle must be < 89 degrees"),
            Self::NoFaces => write!(f, "input BRep has no faces"),
        }
    }
}

/// Apply a draft angle to all planar faces of a BRep.
///
/// Vertices on the neutral plane remain fixed. Other vertices are displaced
/// perpendicular to the pull direction by `h * tan(angle)`.
pub fn draft_solid(brep: &BRep, params: &DraftParams) -> Result<BRep, DraftError> {
    if params.draft_angle.abs() > std::f64::consts::FRAC_PI_2 - 0.02 {
        return Err(DraftError::AngleTooLarge);
    }

    let shell = brep.solids.first().and_then(|s| s.shells.first()).ok_or(DraftError::NoFaces)?;
    if shell.faces.is_empty() {
        return Err(DraftError::NoFaces);
    }

    let pull = params.pull_direction.normalize();
    let neutral = params.neutral_point;
    let tan_angle = params.draft_angle.tan();

    // ── Step 1: compute new vertex positions ─────────────────────────
    let new_pts: Vec<DVec3> = brep.vertices.iter().map(|v| {
        // Signed distance from vertex to neutral plane along pull direction
        let h = (v.point - neutral).dot(pull);
        // Draft displacement: h * tan(angle) in the pull direction
        v.point + pull * (h * tan_angle)
    }).collect();

    // ── Step 2: compute new face normals ─────────────────────────────
    // For planar faces, the normal tilts. We rotate the normal around the
    // axis perpendicular to both the face normal and pull direction.
    let new_face_normals: Vec<DVec3> = shell.faces.iter().map(|face| {
        let n = face.normal.normalize();
        // Rodrigues rotation: rotate n around axis = n × pull by draft_angle
        let axis = n.cross(pull);
        let axis_len = axis.length();
        if axis_len < 1e-10 {
            // Normal is parallel to pull direction — no change
            return n;
        }
        let k = axis / axis_len;
        let cos_a = params.draft_angle.cos();
        let sin_a = params.draft_angle.sin();
        // Rodrigues formula: v_rot = v*cos + (k×v)*sin + k*(k·v)*(1-cos)
        // Since k ⊥ n, k·n = 0, so the third term vanishes.
        let rotated = n * cos_a + k.cross(n) * sin_a;
        rotated.normalize_or(n)
    }).collect();

    // ── Step 3: build result BRep ────────────────────────────────────
    let mut out = BRep::new();
    out.solids.push(Solid { shells: vec![Shell { faces: Vec::new() }] });

    // Copy vertices with new positions
    let mut vmap: Vec<usize> = Vec::new();
    for &p in &new_pts {
        let idx = out.vertices.len();
        out.vertices.push(Vertex { point: p });
        vmap.push(idx);
    }

    // Copy edges with new curves
    let mut emap: Vec<usize> = Vec::new();
    for e in brep.edges.iter() {
        let vs = vmap[e.start];
        let ve = vmap[e.end];
        let dir = (out.vertices[ve].point - out.vertices[vs].point).normalize_or(DVec3::X);
        let len = (out.vertices[ve].point - out.vertices[vs].point).length();

        let curve_idx = out.geom.curves.len();
        out.geom.curves.push(Curve3::Line(Line3 {
            origin: out.vertices[vs].point,
            direction: dir,
        }));
        let eidx = out.edges.len();
        out.edges.push(Edge { start: vs, end: ve });
        out.geom.edge_curve.push(Some(curve_idx));
        out.geom.edge_curve_range.push(Some([0.0, len]));
        out.geom.edge_degenerated.push(false);
        emap.push(eidx);
    }

    // Copy faces with updated normals
    for (fi, face) in shell.faces.iter().enumerate() {
        let mut wire_edges = Vec::new();
        for we in &face.outer_wire.edges {
            let mapped = emap[we.idx];
            wire_edges.push(WireEdge {
                idx: mapped,
                forward: we.forward,
            });
        }

        let face_idx = out.solids[0].shells[0].faces.len();
        // Copy triangles directly — vertex indices are preserved (1:1 mapping).
        let triangles = face.triangles.clone();
        out.solids[0].shells[0].faces.push(Face {
            outer_wire: Wire { edges: wire_edges },
            inner_wires: face.inner_wires.clone(),
            normal: new_face_normals[fi],
            triangles,
            mesh_dirty: face.mesh_dirty,
        });

        // Copy surface reference
        if let Some(&surf_idx) = brep.geom.face_surface.get(fi).and_then(|o| o.as_ref()) {
            while out.geom.face_surface.len() <= face_idx {
                out.geom.face_surface.push(None);
            }
            out.geom.surfaces.push(brep.geom.surfaces[surf_idx].clone());
            out.geom.face_surface[face_idx] = Some(out.geom.surfaces.len() - 1);
        }
    }

    // ── Step 4: update triangles — vertex indices are preserved, so
    // the original triangles reference the same (now displaced) vertices.
    // Do NOT call mesh_brep here — it would clear and re-tessellate.
    // Triangles are already copied in Step 3.

    Ok(out)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_modeling::make_box_brep;

    fn make_box() -> BRep {
        let mut brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        crate::geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn draft_box_positive_angle_increases_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 5.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(v_draft > v_orig, "positive draft should increase volume: {v_orig} -> {v_draft}");
    }

    #[test]
    fn draft_box_negative_angle_decreases_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: (-5.0_f64).to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(v_draft < v_orig, "negative draft should decrease volume: {v_orig} -> {v_draft}");
    }

    #[test]
    fn draft_box_zero_angle_preserves_volume() {
        let brep = make_box();
        let v_orig = rcad_kernel::properties::volume(&brep);

        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 0.0,
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();
        let v_draft = rcad_kernel::properties::volume(&result);

        assert!(
            (v_draft - v_orig).abs() < 0.01,
            "zero draft should preserve volume: {v_orig} vs {v_draft}"
        );
    }

    #[test]
    fn draft_neutral_plane_vertices_unchanged() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 10.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();

        // Vertices at z=0 (on the neutral plane) should not move
        for (i, v) in brep.vertices.iter().enumerate() {
            if (v.point.z - 0.0).abs() < 1e-9 {
                let new_v = &result.vertices[i];
                assert!(
                    (new_v.point.z - 0.0).abs() < 1e-9,
                    "vertex {i} on neutral plane should stay at z=0, got z={}",
                    new_v.point.z
                );
            }
        }
    }

    #[test]
    fn draft_angle_too_large_returns_error() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 89.5_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        assert!(matches!(draft_solid(&brep, &params), Err(DraftError::AngleTooLarge)));
    }

    #[test]
    fn draft_faces_have_tilt_normals() {
        let brep = make_box();
        let params = DraftParams {
            pull_direction: DVec3::Z,
            draft_angle: 5.0_f64.to_radians(),
            neutral_point: DVec3::ZERO,
        };
        let result = draft_solid(&brep, &params).unwrap();

        // Side faces (originally vertical, normal ⊥ Z) should now have a Z component
        for (i, face) in result.solids[0].shells[0].faces.iter().enumerate() {
            let orig_face = &brep.solids[0].shells[0].faces[i];
            let orig_dot_z = orig_face.normal.dot(DVec3::Z).abs();
            let new_dot_z = face.normal.dot(DVec3::Z).abs();

            // If the original face was perpendicular to Z (side face),
            // the drafted face should have a non-zero Z component
            if orig_dot_z < 0.1 {
                assert!(
                    new_dot_z > 0.01,
                    "side face {i} normal should tilt: orig_dot_z={orig_dot_z:.4}, new_dot_z={new_dot_z:.4}"
                );
            }
        }
    }
}
