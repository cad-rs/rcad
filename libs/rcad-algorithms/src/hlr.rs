//! Hidden-Line Removal (HLR).
//!
//! Projects a BRep's edges onto a view plane and classifies each edge segment
//! as **visible** or **hidden** by testing against the silhouette of all faces.
//!
//! Analogous to OCCT `HLRBRep_Algo` / `HLRBRep_HLRToShape`.
//!
//! # Algorithm
//!
//! For each edge:
//! 1. Project both endpoints onto the screen plane.
//! 2. Sample `N` points along the edge in 3D.
//! 3. For each sample, cast a ray from that point toward the camera.
//! 4. If any face triangle blocks the ray **closer** to the camera than the
//!    edge sample, that sample is hidden.
//! 5. Classify runs of consecutive samples → visible/hidden segments.
//!
//! The result is a set of `HlrSegment`s — 2D projected line segments labeled
//! visible or hidden.

use glam::{DMat4, DVec2, DVec3, DVec4};
use rcad_kernel::BRep;

// ── Public types ──────────────────────────────────────────────────────────────

/// A projected edge segment labeled as visible or hidden.
#[derive(Debug, Clone, PartialEq)]
pub struct HlrSegment {
    /// Start point in 2D screen space.
    pub start: DVec2,
    /// End point in 2D screen space.
    pub end: DVec2,
    /// Whether this segment is visible from the camera.
    pub visible: bool,
}

/// Output of an HLR computation.
#[derive(Debug, Clone, Default)]
pub struct HlrResult {
    pub segments: Vec<HlrSegment>,
}

impl HlrResult {
    /// Return only visible segments.
    pub fn visible(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| s.visible)
    }

    /// Return only hidden segments.
    pub fn hidden(&self) -> impl Iterator<Item = &HlrSegment> {
        self.segments.iter().filter(|s| !s.visible)
    }
}

/// Camera / view specification for HLR.
#[derive(Debug, Clone)]
pub struct HlrCamera {
    /// Camera position in world space.
    pub eye: DVec3,
    /// Target point (look-at).
    pub target: DVec3,
    /// Up direction.
    pub up: DVec3,
}

impl HlrCamera {
    pub fn new(eye: DVec3, target: DVec3) -> Self {
        Self { eye, target, up: DVec3::Y }
    }

    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up = up;
        self
    }

    /// Isometric-style view from the +X+Y+Z octant.
    pub fn isometric(distance: f64) -> Self {
        let d = distance / 3.0_f64.sqrt();
        Self::new(DVec3::splat(d), DVec3::ZERO)
    }

    /// Front view (looking along +Y, up = +Z).
    pub fn front(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, -distance, 0.0), DVec3::ZERO)
            .with_up(DVec3::Z)
    }

    /// Top view (looking down -Z).
    pub fn top(distance: f64) -> Self {
        Self::new(DVec3::new(0.0, 0.0, distance), DVec3::ZERO)
            .with_up(DVec3::Y)
    }

    /// Right-side view (looking along -X, up = +Z).
    pub fn right(distance: f64) -> Self {
        Self::new(DVec3::new(distance, 0.0, 0.0), DVec3::ZERO)
            .with_up(DVec3::Z)
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

/// Build a right-handed view matrix (world → camera space).
fn look_at(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let forward = (target - eye).normalize_or_zero();
    let right = forward.cross(up).normalize_or_zero();
    let up = right.cross(forward);

    DMat4::from_cols(
        DVec4::new(right.x,   right.y,   right.z,   -right.dot(eye)),
        DVec4::new(up.x,      up.y,      up.z,      -up.dot(eye)),
        DVec4::new(-forward.x,-forward.y,-forward.z, forward.dot(eye)),
        DVec4::new(0.0,       0.0,       0.0,        1.0),
    ).transpose()
}

/// Project a world-space point to 2D screen space using the view matrix.
/// Returns (x, y) in camera space (z is depth; ignored for 2D output).
fn project(p: DVec3, view: &DMat4) -> (DVec2, f64) {
    let hp = view.mul_vec4(DVec4::new(p.x, p.y, p.z, 1.0));
    (DVec2::new(hp.x, hp.y), hp.z)
}

/// Collect all triangles from a BRep (fan-triangulate faces without pre-triangulated data).
fn collect_triangles(brep: &BRep) -> Vec<[DVec3; 3]> {
    let mut tris = Vec::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if !face.triangles.is_empty() {
                    for &[i, j, k] in &face.triangles {
                        if let (Some(a), Some(b), Some(c)) = (
                            brep.vertices.get(i),
                            brep.vertices.get(j),
                            brep.vertices.get(k),
                        ) {
                            tris.push([a.point, b.point, c.point]);
                        }
                    }
                } else {
                    // Fan-triangulate from wire
                    let pts: Vec<DVec3> = face.outer_wire.edges.iter()
                        .filter_map(|we| {
                            let edge = brep.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            brep.vertices.get(vi).map(|v| v.point)
                        })
                        .collect();
                    if pts.len() >= 3 {
                        let origin = pts[0];
                        for i in 1..pts.len() - 1 {
                            tris.push([origin, pts[i], pts[i + 1]]);
                        }
                    }
                }
            }
        }
    }
    tris
}

/// Ray-triangle intersection (Möller–Trumbore). Returns `Some(t)` if the ray
/// `origin + t*dir` hits the triangle (t > epsilon, front-face only).
fn ray_triangle_intersect(origin: DVec3, dir: DVec3, tri: &[DVec3; 3]) -> Option<f64> {
    const EPS: f64 = 1e-8;
    let edge1 = tri[1] - tri[0];
    let edge2 = tri[2] - tri[0];
    let h = dir.cross(edge2);
    let a = edge1.dot(h);
    if a.abs() < EPS { return None; }
    let f = 1.0 / a;
    let s = origin - tri[0];
    let u = f * s.dot(h);
    if !(0.0..=1.0).contains(&u) { return None; }
    let q = s.cross(edge1);
    let v = f * dir.dot(q);
    if v < 0.0 || u + v > 1.0 { return None; }
    let t = f * edge2.dot(q);
    if t > EPS { Some(t) } else { None }
}

/// Test if a world-space point is occluded by any triangle when viewed from `eye`.
fn is_occluded(point: DVec3, eye: DVec3, triangles: &[[DVec3; 3]], dist_to_eye: f64) -> bool {
    let dir = (eye - point).normalize_or_zero();
    let origin = point + dir * 1e-5; // push off surface
    for tri in triangles {
        if let Some(t) = ray_triangle_intersect(origin, dir, tri) {
            if t < dist_to_eye - 1e-4 {
                return true;
            }
        }
    }
    false
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Perform hidden-line removal on a BRep from the given camera position.
///
/// Returns 2D projected segments labeled visible/hidden.
/// `samples` controls how finely each edge is subdivided for occlusion testing
/// (higher = more accurate but slower; 8 is a reasonable default).
pub fn hlr(brep: &BRep, camera: &HlrCamera, samples: usize) -> HlrResult {
    let view = look_at(camera.eye, camera.target, camera.up);
    let triangles = collect_triangles(brep);
    let samples = samples.max(2);
    let mut result = HlrResult::default();

    // Collect all unique edges from all faces + standalone edges
    let mut edge_indices: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for we in &face.outer_wire.edges {
                    edge_indices.insert(we.idx);
                }
                for inner in &face.inner_wires {
                    for we in &inner.edges {
                        edge_indices.insert(we.idx);
                    }
                }
            }
        }
    }
    for i in 0..brep.edges.len() {
        edge_indices.insert(i);
    }

    for &edge_idx in &edge_indices {
        let Some(edge) = brep.edges.get(edge_idx) else { continue };
        let Some(v_start) = brep.vertices.get(edge.start) else { continue };
        let Some(v_end) = brep.vertices.get(edge.end) else { continue };

        let p0 = v_start.point;
        let p1 = v_end.point;
        if (p1 - p0).length_squared() < 1e-12 { continue; }

        // Sample the edge
        let sample_vis: Vec<bool> = (0..samples).map(|i| {
            let t = i as f64 / (samples - 1) as f64;
            let world_pt = p0 + (p1 - p0) * t;
            let dist = (camera.eye - world_pt).length();
            !is_occluded(world_pt, camera.eye, &triangles, dist)
        }).collect();

        // Build projected 2D points
        let screen_pts: Vec<DVec2> = (0..samples).map(|i| {
            let t = i as f64 / (samples - 1) as f64;
            let world_pt = p0 + (p1 - p0) * t;
            project(world_pt, &view).0
        }).collect();

        // Merge consecutive same-visibility samples into segments
        let mut seg_start = 0usize;
        for i in 1..samples {
            let changed = sample_vis[i] != sample_vis[seg_start];
            let last = i == samples - 1;
            if changed || last {
                let end_idx = if last && !changed { i } else { i - 1 };
                let seg = HlrSegment {
                    start: screen_pts[seg_start],
                    end: screen_pts[end_idx],
                    visible: sample_vis[seg_start],
                };
                // Skip degenerate (zero-length) segments
                if (seg.end - seg.start).length_squared() > 1e-16 {
                    result.segments.push(seg);
                }
                if changed {
                    seg_start = i;
                }
            }
        }
    }

    result
}

/// Render HLR result as a simple SVG string.
///
/// Visible edges are drawn solid black; hidden edges are dashed gray.
/// `scale` controls pixel size per unit.
pub fn hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String {
    if result.segments.is_empty() {
        return "<svg xmlns=\"http://www.w3.org/2000/svg\"/>".to_string();
    }

    // Compute bounding box
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for seg in &result.segments {
        for p in [seg.start, seg.end] {
            min_x = min_x.min(p.x);
            min_y = min_y.min(p.y);
            max_x = max_x.max(p.x);
            max_y = max_y.max(p.y);
        }
    }

    // Flip Y (SVG Y grows downward, camera Y grows upward)
    let transform = |p: DVec2| -> (f64, f64) {
        let x = (p.x - min_x) * scale + margin;
        let y = (max_y - p.y) * scale + margin;
        (x, y)
    };

    let w = (max_x - min_x) * scale + 2.0 * margin;
    let h = (max_y - min_y) * scale + 2.0 * margin;

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">\n",
        w, h, w, h
    );
    svg.push_str("  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");

    for seg in &result.segments {
        let (x1, y1) = transform(seg.start);
        let (x2, y2) = transform(seg.end);
        if seg.visible {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"black\" stroke-width=\"1.5\"/>\n",
                x1, y1, x2, y2
            ));
        } else {
            svg.push_str(&format!(
                "  <line x1=\"{:.3}\" y1=\"{:.3}\" x2=\"{:.3}\" y2=\"{:.3}\" stroke=\"#999\" stroke-width=\"0.8\" stroke-dasharray=\"4,3\"/>\n",
                x1, y1, x2, y2
            ));
        }
    }
    svg.push_str("</svg>\n");
    svg
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_hlr_produces_segments() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        assert!(!result.segments.is_empty(), "HLR should produce segments for a box");
        assert!(result.visible().count() > 0, "some segments should be visible");
    }

    #[test]
    fn hlr_svg_is_valid_xml() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::isometric(5.0);
        let result = hlr(&brep, &camera, 8);
        let svg = hlr_to_svg(&result, 100.0, 20.0);
        assert!(svg.contains("<svg"), "output should be SVG");
        assert!(svg.contains("</svg>"), "SVG should close properly");
        assert!(svg.contains("<line"), "SVG should contain lines");
    }

    #[test]
    fn top_view_box_has_visible_top_edges() {
        let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0, height: 1.0, depth: 1.0,
        });
        let camera = HlrCamera::top(5.0);
        let result = hlr(&brep, &camera, 8);
        let vis = result.visible().count();
        let hid = result.hidden().count();
        assert!(vis > 0, "top view should have visible edges");
        assert!(hid > 0, "top view should have hidden (bottom) edges");
    }
}
