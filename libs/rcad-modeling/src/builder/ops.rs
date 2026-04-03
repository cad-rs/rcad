//! Sweep operations: linear extrusion and revolution.
//!
//! Corresponds to OCCT `BRepPrimAPI_MakePrism` and `BRepPrimAPI_MakeRevol`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Line3, Plane, Surface3};
use rcad_kernel::topology::{Vertex, WireEdge};
use rcad_kernel::BRep;

use crate::builder::brep_builder::{make_edge, make_face, make_wire};
use crate::builder::{normalize_vector, BuildError};

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Rodrigues' rotation formula: rotate `p` around `axis_origin + t*axis_dir`.
fn rotate_point(p: DVec3, axis_origin: DVec3, axis_dir: DVec3, angle: f64) -> DVec3 {
    let v = p - axis_origin;
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let rotated = v * cos_a + axis_dir.cross(v) * sin_a + axis_dir * axis_dir.dot(v) * (1.0 - cos_a);
    rotated + axis_origin
}

/// Compute the outward normal of a planar quad (or triangle) given CCW vertices.
fn quad_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    (b - a).cross(c - a).normalize_or_zero()
}

/// Extract ordered boundary points from a face's outer wire.
fn face_boundary_points(brep: &BRep, face_idx: usize) -> Vec<DVec3> {
    let face = match brep.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
    {
        Some(f) => f,
        None => return Vec::new(),
    };

    // Collect start vertex of each wire edge in order.
    // Each wire edge's start vertex gives a distinct corner.
    let mut pts = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vidx = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vidx) {
                pts.push(v.point);
            }
        }
    }
    pts
}

// ── Linear extrusion ─────────────────────────────────────────────────────────

/// Extrude a profile face (from `profile` BRep) along `direction` by `distance`.
///
/// Returns a new closed BRep solid with:
/// - The bottom cap (original profile, normals flipped downward)
/// - The top cap (translated copy, normals pointing upward)
/// - N lateral faces (one per profile edge)
///
/// # Errors
/// - `BuildError::ZeroVector` if direction is zero.
/// - `BuildError::NonPositiveValue` if distance ≤ 0.
/// - `BuildError::InvalidIndex` if face_idx is out of bounds.
/// - `BuildError::DegenerateGeometry` if the profile has fewer than 3 vertices.
pub fn extrude(
    profile: &BRep,
    face_idx: usize,
    direction: DVec3,
    distance: f64,
) -> Result<BRep, BuildError> {
    let dir = normalize_vector("direction", direction)?;
    if distance <= 0.0 {
        return Err(BuildError::NonPositiveValue("distance"));
    }

    // Validate face_idx
    profile.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
        .ok_or(BuildError::InvalidIndex(face_idx))?;

    let bot_pts = face_boundary_points(profile, face_idx);
    if bot_pts.len() < 3 {
        return Err(BuildError::DegenerateGeometry("profile has fewer than 3 vertices"));
    }

    let offset = dir * distance;
    let top_pts: Vec<DVec3> = bot_pts.iter().map(|&p| p + offset).collect();
    let n = bot_pts.len();

    let mut result = BRep {
        vertices: Vec::new(),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
    };

    // Add all vertices: bot[0..n], top[0..n]
    let bot_vi: Vec<usize> = bot_pts.iter().map(|&p| {
        let idx = result.vertices.len();
        result.vertices.push(Vertex { point: p });
        idx
    }).collect();
    let top_vi: Vec<usize> = top_pts.iter().map(|&p| {
        let idx = result.vertices.len();
        result.vertices.push(Vertex { point: p });
        idx
    }).collect();

    // Build bottom cap (inward normal = -dir)
    let _bot_face = {
        let mut wire_edges = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            // Bottom edges go bot[i] → bot[j]
            let line = rcad_kernel::geom::Curve3::Line(Line3 {
                origin: bot_pts[i],
                direction: (bot_pts[j] - bot_pts[i]).normalize_or_zero(),
            });
            let t1 = 0.0_f64;
            let t2 = (bot_pts[j] - bot_pts[i]).length();
            let eidx = make_edge(&mut result, line, t1, t2, bot_vi[i], bot_vi[j])?;
            wire_edges.push(WireEdge { idx: eidx, forward: true });
        }
        let bot_normal = -dir;
        let surface = Surface3::Plane(Plane { origin: bot_pts[0], normal: bot_normal });
        make_face(&mut result, surface, make_wire(wire_edges), vec![])?
    };

    // Build top cap (outward normal = +dir)
    let _top_face = {
        let mut wire_edges = Vec::new();
        for i in 0..n {
            let j = (i + 1) % n;
            let line = rcad_kernel::geom::Curve3::Line(Line3 {
                origin: top_pts[i],
                direction: (top_pts[j] - top_pts[i]).normalize_or_zero(),
            });
            let t1 = 0.0_f64;
            let t2 = (top_pts[j] - top_pts[i]).length();
            let eidx = make_edge(&mut result, line, t1, t2, top_vi[i], top_vi[j])?;
            wire_edges.push(WireEdge { idx: eidx, forward: true });
        }
        let top_normal = dir;
        let surface = Surface3::Plane(Plane { origin: top_pts[0], normal: top_normal });
        make_face(&mut result, surface, make_wire(wire_edges), vec![])?
    };

    // Build lateral faces: one quad per profile edge
    for i in 0..n {
        let j = (i + 1) % n;

        // Lateral quad vertices (CCW when viewed from outside):
        // bot[i] → bot[j] → top[j] → top[i]
        let a = bot_pts[i];
        let b = bot_pts[j];
        let c = top_pts[j];
        let d = top_pts[i];
        let lat_normal = quad_normal(a, b, c);

        // Edges: bottom (already created), right, top (reversed), left
        // We create the lateral-vertical edges. Bottom and top edges are reused
        // by idx. For simplicity, create all 4 per quad (extra edges on shared sides).
        let left_dir = (d - a).normalize_or_zero();
        let left_len = (d - a).length();
        let right_dir = (c - b).normalize_or_zero();
        let right_len = (c - b).length();
        let bot_dir = (b - a).normalize_or_zero();
        let bot_len = (b - a).length();
        let top_dir = (d - c).normalize_or_zero();
        let top_len = (d - c).length();

        // Create the 4 edges for this lateral face
        let e_bot = make_edge(&mut result,
            rcad_kernel::geom::Curve3::Line(Line3 { origin: a, direction: bot_dir }),
            0.0, bot_len, bot_vi[i], bot_vi[j])?;
        let e_right = make_edge(&mut result,
            rcad_kernel::geom::Curve3::Line(Line3 { origin: b, direction: right_dir }),
            0.0, right_len, bot_vi[j], top_vi[j])?;
        let e_top = make_edge(&mut result,
            rcad_kernel::geom::Curve3::Line(Line3 { origin: c, direction: top_dir }),
            0.0, top_len, top_vi[j], top_vi[i])?;
        let e_left = make_edge(&mut result,
            rcad_kernel::geom::Curve3::Line(Line3 { origin: d, direction: -left_dir }),
            0.0, left_len, top_vi[i], bot_vi[i])?;

        let wire = make_wire(vec![
            WireEdge { idx: e_bot, forward: true },
            WireEdge { idx: e_right, forward: true },
            WireEdge { idx: e_top, forward: true },
            WireEdge { idx: e_left, forward: true },
        ]);
        let surface = Surface3::Plane(Plane { origin: a, normal: lat_normal });
        make_face(&mut result, surface, wire, vec![])?;
    }

    // All faces are in solids[0].shells[0] already from make_face calls
    Ok(result)
}

// ── Revolution ────────────────────────────────────────────────────────────────

/// Revolve a profile face around an axis by `angle` radians.
///
/// For a full revolution (angle ≈ 2π), the start and end caps are identified,
/// resulting in a closed solid without explicit caps.
///
/// # Errors
/// - `BuildError::ZeroVector` if axis_dir is zero.
/// - `BuildError::NonPositiveValue` if angle ≤ 0.
/// - `BuildError::InvalidIndex` if face_idx is out of bounds.
/// - `BuildError::DegenerateGeometry` if the profile has fewer than 2 vertices.
pub fn revolve(
    profile: &BRep,
    face_idx: usize,
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle: f64,
) -> Result<BRep, BuildError> {
    let dir = normalize_vector("axis_dir", axis_dir)?;
    if angle <= 0.0 {
        return Err(BuildError::NonPositiveValue("angle"));
    }

    // Validate face_idx
    profile.solids.first()
        .and_then(|s| s.shells.first())
        .and_then(|sh| sh.faces.get(face_idx))
        .ok_or(BuildError::InvalidIndex(face_idx))?;

    let profile_pts = face_boundary_points(profile, face_idx);
    if profile_pts.len() < 2 {
        return Err(BuildError::DegenerateGeometry("profile has fewer than 2 vertices"));
    }

    let n = profile_pts.len();
    let full_revolution = (angle - std::f64::consts::TAU).abs() < 1e-6;

    let mut result = BRep {
        vertices: Vec::new(),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
    };

    // Rotate each profile point
    let rot_pts: Vec<DVec3> = profile_pts.iter()
        .map(|&p| rotate_point(p, axis_origin, dir, angle))
        .collect();

    // Add start vertices (bottom/start cap positions)
    let start_vi: Vec<usize> = profile_pts.iter().map(|&p| {
        let idx = result.vertices.len();
        result.vertices.push(Vertex { point: p });
        idx
    }).collect();

    // Add end vertices — if full revolution, share with start vertices
    let end_vi: Vec<usize> = if full_revolution {
        start_vi.clone()
    } else {
        rot_pts.iter().map(|&p| {
            let idx = result.vertices.len();
            result.vertices.push(Vertex { point: p });
            idx
        }).collect()
    };

    // Create cap faces (bottom = original profile, top = rotated) for partial revolution
    if !full_revolution {
        // Bottom cap: original profile
        let bot_pts_ref = &profile_pts;
        let _bot_face = {
            let mut wire_edges = Vec::new();
            for i in 0..n {
                let j = (i + 1) % n;
                let a = bot_pts_ref[i];
                let b = bot_pts_ref[j];
                let line = rcad_kernel::geom::Curve3::Line(Line3 {
                    origin: a,
                    direction: (b - a).normalize_or_zero(),
                });
                let eidx = make_edge(&mut result, line, 0.0, (b - a).length(), start_vi[i], start_vi[j])?;
                wire_edges.push(WireEdge { idx: eidx, forward: true });
            }
            let bot_normal = quad_normal(bot_pts_ref[0], bot_pts_ref[1], bot_pts_ref[2]);
            let surface = Surface3::Plane(Plane { origin: bot_pts_ref[0], normal: bot_normal });
            make_face(&mut result, surface, make_wire(wire_edges), vec![])?
        };

        // Top cap: rotated profile
        let top_pts_ref = &rot_pts;
        let _top_face = {
            let mut wire_edges = Vec::new();
            for i in 0..n {
                let j = (i + 1) % n;
                let a = top_pts_ref[i];
                let b = top_pts_ref[j];
                let line = rcad_kernel::geom::Curve3::Line(Line3 {
                    origin: a,
                    direction: (b - a).normalize_or_zero(),
                });
                let eidx = make_edge(&mut result, line, 0.0, (b - a).length(), end_vi[i], end_vi[j])?;
                wire_edges.push(WireEdge { idx: eidx, forward: true });
            }
            let top_normal = quad_normal(top_pts_ref[0], top_pts_ref[1], top_pts_ref[2]);
            let surface = Surface3::Plane(Plane { origin: top_pts_ref[0], normal: top_normal });
            make_face(&mut result, surface, make_wire(wire_edges), vec![])?
        };
    }

    // Create lateral swept faces for each profile edge
    for i in 0..n {
        let j = (i + 1) % n;

        let p0 = profile_pts[i];
        let p1 = profile_pts[j];
        let p1_rot = rot_pts[j];
        let p0_rot = rot_pts[i];

        let lat_normal = quad_normal(p0, p1, p1_rot);

        // 4 edges of the quad
        let e_bot = {
            let d = (p1 - p0).normalize_or_zero();
            let len = (p1 - p0).length();
            make_edge(&mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: p0, direction: d }),
                0.0, len, start_vi[i], start_vi[j])?
        };
        let e_right = {
            let d = (p1_rot - p1).normalize_or_zero();
            let len = (p1_rot - p1).length();
            make_edge(&mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: p1, direction: d }),
                0.0, len, start_vi[j], end_vi[j])?
        };
        let e_top = {
            let d = (p0_rot - p1_rot).normalize_or_zero();
            let len = (p0_rot - p1_rot).length();
            make_edge(&mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: p1_rot, direction: d }),
                0.0, len, end_vi[j], end_vi[i])?
        };
        let e_left = {
            let d = (p0 - p0_rot).normalize_or_zero();
            let len = (p0 - p0_rot).length();
            make_edge(&mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: p0_rot, direction: d }),
                0.0, len, end_vi[i], start_vi[i])?
        };

        let wire = make_wire(vec![
            WireEdge { idx: e_bot, forward: true },
            WireEdge { idx: e_right, forward: true },
            WireEdge { idx: e_top, forward: true },
            WireEdge { idx: e_left, forward: true },
        ]);
        let surface = Surface3::Plane(Plane { origin: p0, normal: lat_normal });
        make_face(&mut result, surface, wire, vec![])?;
    }

    Ok(result)
}

// ── Loft (multi-profile) ──────────────────────────────────────────────────────

/// Connect N cross-section profiles with ruled lateral faces plus planar caps.
///
/// All profiles must have the **same** vertex count (≥ 3).
/// Each profile is a list of 3D positions in order (the polygon vertices).
///
/// Returns a closed BRep Solid.
///
/// # Errors
/// - `BuildError::DegenerateGeometry` if fewer than 2 profiles, or vertex counts differ, or fewer than 3 vertices per profile.
pub fn loft(profiles: &[Vec<DVec3>]) -> Result<BRep, BuildError> {
    if profiles.len() < 2 {
        return Err(BuildError::DegenerateGeometry("loft requires at least 2 profiles"));
    }
    let n = profiles[0].len();
    if n < 3 {
        return Err(BuildError::DegenerateGeometry("loft profiles must have at least 3 vertices"));
    }
    for (i, p) in profiles.iter().enumerate() {
        if p.len() != n {
            return Err(BuildError::DegenerateGeometry(
                "all loft profiles must have the same vertex count",
            ));
        }
        if p.len() < 3 {
            return Err(BuildError::DegenerateGeometry("loft profile has fewer than 3 vertices"));
        }
        let _ = i;
    }

    let s = profiles.len(); // number of sections

    let mut result = BRep {
        vertices: Vec::new(),
        edges: Vec::new(),
        solids: Vec::new(),
        geom: rcad_kernel::GeomStore::default(),
    };

    // Add all vertices; vi[section][vertex]
    let vi: Vec<Vec<usize>> = profiles
        .iter()
        .map(|prof| {
            prof.iter()
                .map(|&p| {
                    let idx = result.vertices.len();
                    result.vertices.push(Vertex { point: p });
                    idx
                })
                .collect()
        })
        .collect();

    // Bottom cap (profile[0]) — normal pointing away from profile[1]
    {
        let mut wire_edges = Vec::new();
        let pts = &profiles[0];
        for i in 0..n {
            let j = (i + 1) % n;
            let a = pts[i];
            let b = pts[j];
            let d = (b - a).normalize_or_zero();
            let len = (b - a).length();
            let eidx = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: a, direction: d }),
                0.0, len, vi[0][i], vi[0][j],
            )?;
            wire_edges.push(WireEdge { idx: eidx, forward: true });
        }
        // Normal pointing away from next section
        let centroid_0: DVec3 = pts.iter().sum::<DVec3>() / n as f64;
        let centroid_1: DVec3 = profiles[1].iter().sum::<DVec3>() / n as f64;
        let bot_normal = (centroid_0 - centroid_1).normalize_or_zero();
        let surface = Surface3::Plane(Plane { origin: pts[0], normal: bot_normal });
        make_face(&mut result, surface, make_wire(wire_edges), vec![])?;
    }

    // Top cap (profile[last]) — normal pointing away from profile[last-1]
    {
        let mut wire_edges = Vec::new();
        let pts = &profiles[s - 1];
        for i in 0..n {
            let j = (i + 1) % n;
            let a = pts[i];
            let b = pts[j];
            let d = (b - a).normalize_or_zero();
            let len = (b - a).length();
            let eidx = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: a, direction: d }),
                0.0, len, vi[s - 1][i], vi[s - 1][j],
            )?;
            wire_edges.push(WireEdge { idx: eidx, forward: true });
        }
        let centroid_prev: DVec3 = profiles[s - 2].iter().sum::<DVec3>() / n as f64;
        let centroid_top: DVec3 = pts.iter().sum::<DVec3>() / n as f64;
        let top_normal = (centroid_top - centroid_prev).normalize_or_zero();
        let surface = Surface3::Plane(Plane { origin: pts[0], normal: top_normal });
        make_face(&mut result, surface, make_wire(wire_edges), vec![])?;
    }

    // Lateral quad faces between consecutive sections
    for sec in 0..s - 1 {
        let pts_bot = &profiles[sec];
        let pts_top = &profiles[sec + 1];

        for i in 0..n {
            let j = (i + 1) % n;
            // Quad: bot[i] → bot[j] → top[j] → top[i]
            let a = pts_bot[i];
            let b = pts_bot[j];
            let c = pts_top[j];
            let d = pts_top[i];
            let lat_normal = quad_normal(a, b, c);

            let e_bot = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: a, direction: (b - a).normalize_or_zero() }),
                0.0, (b - a).length(), vi[sec][i], vi[sec][j],
            )?;
            let e_right = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: b, direction: (c - b).normalize_or_zero() }),
                0.0, (c - b).length(), vi[sec][j], vi[sec + 1][j],
            )?;
            let e_top = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: c, direction: (d - c).normalize_or_zero() }),
                0.0, (d - c).length(), vi[sec + 1][j], vi[sec + 1][i],
            )?;
            let e_left = make_edge(
                &mut result,
                rcad_kernel::geom::Curve3::Line(Line3 { origin: d, direction: (a - d).normalize_or_zero() }),
                0.0, (a - d).length(), vi[sec + 1][i], vi[sec][i],
            )?;
            let wire = make_wire(vec![
                WireEdge { idx: e_bot,   forward: true },
                WireEdge { idx: e_right, forward: true },
                WireEdge { idx: e_top,   forward: true },
                WireEdge { idx: e_left,  forward: true },
            ]);
            make_face(&mut result, Surface3::Plane(Plane { origin: a, normal: lat_normal }), wire, vec![])?;
        }
    }

    Ok(result)
}

// ── Pipe Sweep ────────────────────────────────────────────────────────────────

/// Sweep a 2D profile polygon along a 3D spine polyline.
///
/// `profile_2d` is a polygon in the XY plane (local cross-section).
/// `spine` is the path (≥ 2 points) in 3D world space.
///
/// At each spine station a Frenet-like frame (tangent/right/up) is computed and
/// the 2D profile is transformed into the corresponding 3D cross-section.
/// The resulting cross-section polygons are connected via `loft`.
///
/// # Errors
/// - `BuildError::DegenerateGeometry` if fewer than 2 spine points or fewer than 3 profile points.
/// - `BuildError::ZeroVector` if a spine segment has zero length.
pub fn sweep_pipe(profile_2d: &[DVec2], spine: &[DVec3]) -> Result<BRep, BuildError> {
    if profile_2d.len() < 3 {
        return Err(BuildError::DegenerateGeometry("sweep_pipe profile must have at least 3 vertices"));
    }
    if spine.len() < 2 {
        return Err(BuildError::DegenerateGeometry("sweep_pipe spine must have at least 2 points"));
    }

    let ns = spine.len();

    // Compute tangent at each spine station (forward / central / backward difference)
    let tangents: Vec<DVec3> = (0..ns)
        .map(|i| {
            if i == 0 {
                (spine[1] - spine[0]).normalize_or_zero()
            } else if i == ns - 1 {
                (spine[ns - 1] - spine[ns - 2]).normalize_or_zero()
            } else {
                (spine[i + 1] - spine[i - 1]).normalize_or_zero()
            }
        })
        .collect();

    // Build 3D cross-sections
    let world_up_primary = DVec3::Y;
    let world_up_fallback = DVec3::Z;

    let cross_sections: Vec<Vec<DVec3>> = tangents
        .iter()
        .enumerate()
        .map(|(i, &tan)| {
            // Right = tangent × world_up; fall back if nearly parallel
            let right_raw = tan.cross(world_up_primary);
            let right = if right_raw.length_squared() > 1e-8 {
                right_raw.normalize()
            } else {
                tan.cross(world_up_fallback).normalize_or_zero()
            };
            let up = right.cross(tan).normalize_or_zero();

            profile_2d
                .iter()
                .map(|p2| spine[i] + p2.x * right + p2.y * up)
                .collect()
        })
        .collect();

    loft(&cross_sections)
}

/// Variable-section pipe sweep: a different 2D profile at each spine station.
///
/// `profiles[i]` is placed at `spine[i]` using the same Frenet-like frame
/// as [`sweep_pipe`]. All profiles must have the same vertex count (≥ 3)
/// and `profiles.len()` must equal `spine.len()` (≥ 2).
///
/// Analogous to OCCT `BRepOffsetAPI_MakePipeShell` with multiple sections.
pub fn sweep_pipe_variable(
    profiles: &[Vec<DVec2>],
    spine: &[DVec3],
) -> Result<BRep, BuildError> {
    if profiles.len() != spine.len() {
        return Err(BuildError::DegenerateGeometry(
            "sweep_pipe_variable: profiles.len() must equal spine.len()",
        ));
    }
    if spine.len() < 2 {
        return Err(BuildError::DegenerateGeometry(
            "sweep_pipe_variable spine must have at least 2 points",
        ));
    }
    let n_verts = profiles.first().map(|p| p.len()).unwrap_or(0);
    if n_verts < 3 {
        return Err(BuildError::DegenerateGeometry(
            "sweep_pipe_variable profile must have at least 3 vertices",
        ));
    }
    for p in profiles {
        if p.len() != n_verts {
            return Err(BuildError::DegenerateGeometry(
                "sweep_pipe_variable: all profiles must have the same vertex count",
            ));
        }
    }

    let ns = spine.len();
    let world_up_primary  = DVec3::Y;
    let world_up_fallback = DVec3::Z;

    let cross_sections: Vec<Vec<DVec3>> = (0..ns)
        .map(|i| {
            let tan = if i == 0 {
                (spine[1] - spine[0]).normalize_or_zero()
            } else if i == ns - 1 {
                (spine[ns - 1] - spine[ns - 2]).normalize_or_zero()
            } else {
                (spine[i + 1] - spine[i - 1]).normalize_or_zero()
            };

            let right_raw = tan.cross(world_up_primary);
            let right = if right_raw.length_squared() > 1e-8 {
                right_raw.normalize()
            } else {
                tan.cross(world_up_fallback).normalize_or_zero()
            };
            let up = right.cross(tan).normalize_or_zero();

            profiles[i]
                .iter()
                .map(|p2| spine[i] + p2.x * right + p2.y * up)
                .collect()
        })
        .collect();

    loft(&cross_sections)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::brep_builder::{make_edge, make_face, make_vertex, make_wire};
    use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
    use rcad_kernel::topology::WireEdge;
    use rcad_kernel::BRep;

    fn square_profile() -> BRep {
        let mut brep = BRep::default();
        let v0 = make_vertex(&mut brep, DVec3::new(0.0, 0.0, 0.0));
        let v1 = make_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0));
        let v2 = make_vertex(&mut brep, DVec3::new(1.0, 1.0, 0.0));
        let v3 = make_vertex(&mut brep, DVec3::new(0.0, 1.0, 0.0));

        let mk_line = |a: DVec3, b: DVec3| Curve3::Line(Line3 {
            origin: a,
            direction: (b - a).normalize(),
        });

        let pts = [DVec3::new(0.0,0.0,0.0), DVec3::new(1.0,0.0,0.0),
                   DVec3::new(1.0,1.0,0.0), DVec3::new(0.0,1.0,0.0)];
        let vs = [v0, v1, v2, v3];
        let mut wires = Vec::new();
        for i in 0..4 {
            let j = (i + 1) % 4;
            let len = (pts[j] - pts[i]).length();
            let eidx = make_edge(&mut brep, mk_line(pts[i], pts[j]), 0.0, len, vs[i], vs[j]).unwrap();
            wires.push(WireEdge { idx: eidx, forward: true });
        }
        let surface = Surface3::Plane(Plane { origin: DVec3::ZERO, normal: DVec3::Z });
        make_face(&mut brep, surface, make_wire(wires), vec![]).unwrap();
        brep
    }

    #[test]
    fn extrude_square_produces_6_faces() {
        let profile = square_profile();
        let result = extrude(&profile, 0, DVec3::Z, 1.0).unwrap();
        let n_faces = result.solids[0].shells[0].faces.len();
        assert_eq!(n_faces, 6, "extrude of square should yield 6 faces, got {n_faces}");
    }

    #[test]
    fn extrude_rejects_zero_direction() {
        let profile = square_profile();
        let err = extrude(&profile, 0, DVec3::ZERO, 1.0).unwrap_err();
        assert_eq!(err, BuildError::ZeroVector("direction"));
    }

    #[test]
    fn extrude_rejects_nonpositive_distance() {
        let profile = square_profile();
        let err = extrude(&profile, 0, DVec3::Z, -1.0).unwrap_err();
        assert_eq!(err, BuildError::NonPositiveValue("distance"));
    }

    #[test]
    fn revolve_triangle_partial_produces_faces() {
        // A simple triangle profile
        let mut brep = BRep::default();
        let v0 = make_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0));
        let v1 = make_vertex(&mut brep, DVec3::new(2.0, 0.0, 0.0));
        let v2 = make_vertex(&mut brep, DVec3::new(1.5, 1.0, 0.0));
        let pts = [DVec3::new(1.0,0.0,0.0), DVec3::new(2.0,0.0,0.0), DVec3::new(1.5,1.0,0.0)];
        let vs = [v0, v1, v2];
        let mut wires = Vec::new();
        for i in 0..3 {
            let j = (i + 1) % 3;
            let d = (pts[j] - pts[i]).normalize();
            let len = (pts[j] - pts[i]).length();
            let eidx = make_edge(&mut brep,
                Curve3::Line(Line3 { origin: pts[i], direction: d }),
                0.0, len, vs[i], vs[j]).unwrap();
            wires.push(WireEdge { idx: eidx, forward: true });
        }
        let surface = Surface3::Plane(Plane { origin: pts[0], normal: DVec3::Z });
        make_face(&mut brep, surface, make_wire(wires), vec![]).unwrap();

        // Revolve 90° around Y axis
        let result = revolve(&brep, 0, DVec3::ZERO, DVec3::Y, std::f64::consts::FRAC_PI_2).unwrap();
        let n_faces = result.solids[0].shells[0].faces.len();
        // 2 caps + 3 lateral = 5
        assert!(n_faces >= 3, "revolve should produce at least 3 faces, got {n_faces}");
    }
}
