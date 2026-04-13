//! Section: plane-solid intersection returning Wire loops.
//!
//! Analogous to OCCT `BRepAlgoAPI_Section`. Computes the intersection of a
//! cutting plane with the faces of a BRep, returning a list of closed `Wire`s
//! (polyline approximations of the section curves).
//!
//! Also provides [`section_curves`] which returns exact analytic curves
//! (`Circle`, `Ellipse`, `Line`) when the face has an analytic surface.
//!
//! # Algorithm
//!
//! For each triangulated face in the BRep:
//! 1. Intersect each triangle edge with the cutting plane.
//! 2. Collect intersection segments (one per triangle that straddles the plane).
//! 3. Chain segments into closed or open polyline loops.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve3, Line3, Plane};
use rcad_kernel::topology::{Edge, Shell, Solid, Vertex, Wire, WireEdge};

use crate::triangulate::triangulate_polygon;

// ── Internal helpers ──────────────────────────────────────────────────────────

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
    if da.signum() == db.signum() || (da.abs() < 1e-10 && db.abs() < 1e-10) {
        return None;
    }
    if da.abs() < 1e-10 {
        return Some(a);
    }
    if db.abs() < 1e-10 {
        return Some(b);
    }
    let t = da / (da - db);
    Some(a + t * (b - a))
}

/// Collect triangles for a face (pre-triangulated or fan-triangulated from wire).
fn face_triangles(brep: &BRep, face: &rcad_kernel::Face) -> Vec<[DVec3; 3]> {
    if !face.triangles.is_empty() {
        return face
            .triangles
            .iter()
            .filter_map(|&[i, j, k]| {
                let a = brep.vertices.get(i)?.point;
                let b = brep.vertices.get(j)?.point;
                let c = brep.vertices.get(k)?.point;
                Some([a, b, c])
            })
            .collect();
    }

    // Fan-triangulate from wire vertices
    let wire_pts: Vec<DVec3> = face
        .outer_wire
        .edges
        .iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vidx = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vidx).map(|v| v.point)
        })
        .collect();

    if wire_pts.len() < 3 {
        return Vec::new();
    }

    // Need a normal to triangulate
    let normal = face.normal;
    let tris = triangulate_polygon(&wire_pts, normal);
    tris.iter()
        .filter_map(|&[i, j, k]| {
            let a = wire_pts.get(i)?;
            let b = wire_pts.get(j)?;
            let c = wire_pts.get(k)?;
            Some([*a, *b, *c])
        })
        .collect()
}

/// Intersect a single triangle with the plane. Returns a segment [p0, p1] if
/// the triangle straddles the plane, or `None` otherwise.
fn triangle_section(plane: &Plane, tri: [DVec3; 3]) -> Option<[DVec3; 2]> {
    let [a, b, c] = tri;
    let edges = [[a, b], [b, c], [c, a]];
    let mut pts = Vec::new();
    for [p, q] in edges {
        if let Some(hit) = segment_plane_intersect(plane, p, q) {
            // Deduplicate near-identical hits (e.g. at a vertex)
            if pts.iter().all(|&x: &DVec3| (x - hit).length() > 1e-8) {
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
fn pts_close(a: DVec3, b: DVec3) -> bool {
    (a - b).length() < 1e-6
}

/// Chain a set of unordered segments into ordered polylines.
///
/// Returns a list of loops (each loop is an ordered list of DVec3 points).
/// Attempts to close loops; open chains are also returned as-is.
fn chain_segments(segments: Vec<[DVec3; 2]>) -> Vec<Vec<DVec3>> {
    if segments.is_empty() {
        return Vec::new();
    }

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
                if pts_close(remaining[i][0], tail) {
                    chain.push(remaining[i][1]);
                    remaining.remove(i);
                    extended = true;
                    break;
                } else if pts_close(remaining[i][1], tail) {
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
                if pts_close(remaining[i][1], head) {
                    chain.insert(0, remaining[i][0]);
                    remaining.remove(i);
                    extended = true;
                    break;
                } else if pts_close(remaining[i][0], head) {
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

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute the section of a BRep with a cutting plane.
///
/// Returns a new BRep containing only edges and wires (no faces/solids)
/// representing the section curves. Each closed loop is a separate wire.
///
/// For rendering, callers can extract vertices from the returned BRep's wires.
///
/// Analogous to OCCT `BRepAlgoAPI_Section`.
pub fn section(brep: &BRep, plane: &Plane) -> BRep {
    // Collect all section segments from all triangles
    let mut segments: Vec<[DVec3; 2]> = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for tri in face_triangles(brep, face) {
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
    let loops = chain_segments(segments);

    // Build result BRep
    let mut result = BRep::new();
    let mut wire_list: Vec<Wire> = Vec::new();

    for loop_pts in loops {
        if loop_pts.len() < 2 {
            continue;
        }

        let mut wire_edges = Vec::new();

        // Add vertices and edges for the loop
        for i in 0..loop_pts.len().saturating_sub(1) {
            let a = loop_pts[i];
            let b = loop_pts[i + 1];

            let vi_a = result.vertices.len();
            result.vertices.push(Vertex { point: a });
            let vi_b = result.vertices.len();
            result.vertices.push(Vertex { point: b });

            let edge_idx = result.edges.len();
            result.edges.push(Edge {
                start: vi_a,
                end: vi_b,
            });

            // Register curve in geom
            let len = (b - a).length();
            let dir = if len > 1e-10 { (b - a) / len } else { DVec3::X };
            let curve_idx = result.geom.curves.len();
            result.geom.curves.push(Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }));

            while result.geom.edge_curve.len() <= edge_idx {
                result.geom.edge_curve.push(None);
            }
            while result.geom.edge_curve_range.len() <= edge_idx {
                result.geom.edge_curve_range.push(None);
            }
            while result.geom.edge_degenerated.len() <= edge_idx {
                result.geom.edge_degenerated.push(false);
            }
            result.geom.edge_curve[edge_idx] = Some(curve_idx);
            result.geom.edge_curve_range[edge_idx] = Some([0.0, len]);

            wire_edges.push(WireEdge::fwd(edge_idx));
        }

        wire_list.push(Wire { edges: wire_edges });
    }

    // Pack wires into a single shell/solid so callers can iterate normally
    // Each loop becomes a face-less shell entry. We use a minimal Solid with
    // one "open shell" per wire.
    // For simplicity, pack all wires as a flat list in a degenerate solid.
    if !wire_list.is_empty() {
        // Store wires as faces with no surface (open section wires, not closed faces)
        use rcad_kernel::topology::Face;
        let faces: Vec<_> = wire_list
            .into_iter()
            .map(|w| Face {
                outer_wire: w,
                inner_wires: vec![],
                normal: plane.normal,
                triangles: vec![],
                mesh_dirty: true,
            })
            .collect();
        result.solids.push(Solid {
            shells: vec![Shell { faces }],
        });
    }

    result
}

/// Convenience: extract all section polylines as ordered lists of 3D points.
///
/// Each entry is one closed (or open) loop of points from the plane section.
pub fn section_polylines(brep: &BRep, plane: &Plane) -> Vec<Vec<DVec3>> {
    // Collect segments directly without building full BRep
    let mut segments: Vec<[DVec3; 2]> = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for tri in face_triangles(brep, face) {
                    if let Some(seg) = triangle_section(plane, tri) {
                        segments.push(seg);
                    }
                }
            }
        }
    }

    chain_segments(segments)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

// ── Analytic section curves ───────────────────────────────────────────────────

/// One result curve from [`section_curves`].
#[derive(Debug, Clone)]
pub enum SectionCurve {
    /// Exact analytic curve returned when the face has a recognized analytic surface.
    Analytic(Curve3),
    /// Polyline fallback for parametric surfaces (BSpline, Bezier, Offset, Torus, …).
    Polyline(Vec<DVec3>),
}

/// Section a BRep with a plane, returning analytic curves where possible.
///
/// For faces backed by `Plane`, `Sphere`, `Cylinder`, or `Cone` surfaces the
/// function dispatches to the exact analytical intersection tools and returns
/// `SectionCurve::Analytic`. For all other surfaces it falls back to the
/// triangle-mesh polyline method and returns `SectionCurve::Polyline`.
///
/// Curves that do not intersect the given plane are silently omitted.
///
/// Analogous to OCCT `BRepAlgoAPI_Section` returning proper edge geometry.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::{BRep, geom::{Plane, PrimitiveSolid}};
/// use rcad_algorithms::section_curves;
///
/// let sphere = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 2.0 });
/// let plane = Plane { origin: DVec3::ZERO, normal: DVec3::Z };
/// let curves = section_curves(&sphere, &plane);
/// // Equatorial section of a sphere → one Circle
/// assert!(!curves.is_empty());
/// ```
pub fn section_curves(brep: &BRep, plane: &Plane) -> Vec<SectionCurve> {
    use crate::inttools::{
        plane_cone::{PlaneConicalResult, intersect_plane_cone},
        plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
        plane_plane::{PlanePlaneResult, intersect_plane_plane},
        plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
    };
    use rcad_kernel::geom::Surface3;

    let mut results: Vec<SectionCurve> = Vec::new();

    if brep.solids.is_empty() {
        return results;
    }

    let mut face_global_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Look up the analytic surface for this face
                let surf_opt = brep
                    .geom
                    .face_surface
                    .get(face_global_idx)
                    .and_then(|o| *o)
                    .and_then(|si| brep.geom.surfaces.get(si));

                if let Some(surface) = surf_opt {
                    let analytic = match surface {
                        Surface3::Plane(face_plane) => {
                            match intersect_plane_plane(plane, face_plane) {
                                PlanePlaneResult::Line(line) => Some(Curve3::Line(line)),
                                _ => None,
                            }
                        }
                        Surface3::Sphere(sph) => match intersect_plane_sphere(plane, sph) {
                            PlaneSphereResult::Circle(c) => Some(Curve3::Circle(c)),
                            PlaneSphereResult::TangentPoint(_) => None,
                            PlaneSphereResult::NoIntersection => None,
                        },
                        Surface3::Cylinder(cyl) => match intersect_plane_cylinder(plane, cyl) {
                            PlaneCylinderResult::Circle(c) => Some(Curve3::Circle(c)),
                            PlaneCylinderResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
                            PlaneCylinderResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
                            PlaneCylinderResult::TangentLine(_) => None,
                            PlaneCylinderResult::NoIntersection => None,
                        },
                        Surface3::Cone(cone) => match intersect_plane_cone(plane, cone) {
                            PlaneConicalResult::Circle(c) => Some(Curve3::Circle(c)),
                            PlaneConicalResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
                            PlaneConicalResult::Parabola(par) => Some(Curve3::Parabola(par)),
                            PlaneConicalResult::Hyperbola(hyp) => Some(Curve3::Hyperbola(hyp)),
                            PlaneConicalResult::SingleLine(l) => Some(Curve3::Line(l)),
                            PlaneConicalResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
                            PlaneConicalResult::Point(_) => None,
                            PlaneConicalResult::NoIntersection => None,
                        },
                        // All other surfaces: use polyline fallback
                        _ => {
                            let segs: Vec<[DVec3; 2]> = face_triangles(brep, face)
                                .into_iter()
                                .filter_map(|tri| triangle_section(plane, tri))
                                .collect();
                            if !segs.is_empty() {
                                let chains = chain_segments(segs);
                                for chain in chains {
                                    if chain.len() >= 2 {
                                        results.push(SectionCurve::Polyline(chain));
                                    }
                                }
                            }
                            face_global_idx += 1;
                            continue;
                        }
                    };

                    if let Some(curve) = analytic {
                        results.push(SectionCurve::Analytic(curve));
                    }
                } else {
                    // No analytic surface: triangle fallback
                    let segs: Vec<[DVec3; 2]> = face_triangles(brep, face)
                        .into_iter()
                        .filter_map(|tri| triangle_section(plane, tri))
                        .collect();
                    if !segs.is_empty() {
                        let chains = chain_segments(segs);
                        for chain in chains {
                            if chain.len() >= 2 {
                                results.push(SectionCurve::Polyline(chain));
                            }
                        }
                    }
                }

                face_global_idx += 1;
            }
        }
    }

    results
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn section_of_unit_box_at_midplane_z() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 0.5),
            normal: DVec3::Z,
        };

        let polylines = section_polylines(&brep, &plane);
        assert!(
            !polylines.is_empty(),
            "section of unit box should yield at least one loop"
        );

        // All points should be at z ≈ 0.5
        for poly in &polylines {
            for &p in poly {
                assert!(
                    (p.z - 0.5).abs() < 1e-5,
                    "section point z should be 0.5, got {}",
                    p.z
                );
            }
        }
    }

    #[test]
    fn section_misses_when_plane_outside() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let plane = Plane {
            origin: DVec3::new(0.0, 0.0, 5.0),
            normal: DVec3::Z,
        };

        let polylines = section_polylines(&brep, &plane);
        assert!(polylines.is_empty(), "section outside box should be empty");
    }

    #[test]
    fn section_points_within_box_bounds() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        let plane = Plane {
            origin: DVec3::new(0.0, 1.5, 0.0),
            normal: DVec3::Y,
        };

        let polylines = section_polylines(&brep, &plane);
        assert!(!polylines.is_empty());

        for poly in &polylines {
            for &p in poly {
                assert!(p.x >= -1e-5 && p.x <= 2.0 + 1e-5);
                assert!(p.z >= -1e-5 && p.z <= 4.0 + 1e-5);
            }
        }
    }
}
