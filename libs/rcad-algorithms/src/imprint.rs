//! Phase R.B — Face imprinting and gap/overlap detection.
//!
//! **Face imprinting** (`imprint_brep`): splits each face of `target` wherever
//! the boundary of `tool` crosses it, without performing a boolean classification.
//! The result is a new BRep whose faces share edges with the tool boundary — a
//! prerequisite for conformal meshing (FEM/FDTD).
//!
//! **Gap/overlap detection** (`detect_gaps_overlaps`): reports pairs of faces
//! from two BReps that are either too close (gap) or interpenetrating (overlap),
//! using face bounding-box pre-filtering and `closest_point_on_surface`.

use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::projection::closest_point_on_surface;
use rcad_kernel::topology::*;
use rcad_kernel::{BRep};

use crate::bopds::ds::{DS, ShapeOrigin};
use crate::builder::SubFace;
use crate::pave_filler::PaveFiller;
use crate::triangulate::triangulate_polygon;

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Result of imprinting `tool` geometry onto `target`.
#[derive(Debug)]
pub struct ImprintResult {
    /// Modified target BRep whose faces are split wherever the tool boundary crosses.
    pub brep: BRep,
    /// Pairs of (target face index in result, source tool face index) that share
    /// a seam (imprinted edge).
    pub seam_edges: Vec<(usize, usize)>,
}

/// Detected gap between two faces.
#[derive(Debug, Clone)]
pub struct Gap {
    /// Face index in BRep A.
    pub face_a: usize,
    /// Face index in BRep B.
    pub face_b: usize,
    /// Maximum gap distance found between the two faces.
    pub max_gap: f64,
    /// A world-space point on face A that is closest to face B.
    pub sample_point: DVec3,
}

/// Detected overlap (interpenetration) between two faces.
#[derive(Debug, Clone)]
pub struct Overlap {
    /// Face index in BRep A.
    pub face_a: usize,
    /// Face index in BRep B.
    pub face_b: usize,
    /// Estimated penetration depth (positive = overlapping).
    pub penetration_depth: f64,
}

/// Report from gap/overlap detection.
#[derive(Debug, Default)]
pub struct GapOverlapReport {
    pub gaps: Vec<Gap>,
    pub overlaps: Vec<Overlap>,
    /// Pairs of faces (a, b) that are perfectly coincident and coplanar.
    pub shared_faces: Vec<(usize, usize)>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Face imprinting
// ─────────────────────────────────────────────────────────────────────────────

/// Imprint the boundary of `tool` onto the faces of `target`.
///
/// This runs the PaveFiller intersection pass between the two BReps, then splits
/// each target face by the intersection curves recorded in its `FaceInfo`.
/// No boolean classification is performed — all faces of `target` are preserved,
/// but split where the tool boundary crosses them.
///
/// Analogy: OCCT `BRepAlgoAPI_Splitter` (lightweight variant — keeps all target faces).
pub fn imprint_brep(target: &BRep, tool: &BRep) -> ImprintResult {
    // Run PaveFiller to compute intersections
    let mut ds = DS::new(target, tool);
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();

    // Identify which DS faces came from target (ShapeA)
    let target_face_indices: Vec<usize> = ds
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.origin == ShapeOrigin::ShapeA)
        .map(|(i, _)| i)
        .collect();

    // Identify tool face indices per DS face (for seam tracking)
    let tool_face_indices: Vec<usize> = ds
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.origin == ShapeOrigin::ShapeB)
        .map(|(i, _)| i)
        .collect();

    let mut result_faces: Vec<Face> = Vec::new();
    let mut seam_edges: Vec<(usize, usize)> = Vec::new();

    for &dfi in &target_face_indices {
        let sub_faces = split_face_by_curves(&ds, dfi);

        let has_intersection = !ds.faces[dfi].face_info.curves_in.is_empty();
        let result_face_start = result_faces.len();

        for sf in sub_faces {
            let triangles = triangulate_polygon(&sf.boundary, sf.normal);
            result_faces.push(Face {
                outer_wire: Wire {
                    edges: sf
                        .boundary
                        .windows(2)
                        .enumerate()
                        .map(|(i, _)| WireEdge {
                            idx: i,
                            forward: true,
                        })
                        .collect(),
                },
                inner_wires: vec![],
                normal: sf.normal,
                triangles,
            });
        }

        // Record seam edges for every tool face that has curves on this target face
        if has_intersection {
            for &tfi in &tool_face_indices {
                // Check if any curve on this target face came from a FF interference with this tool face
                let shares_curve = ds.interferences.iter().any(|iv| {
                    if let crate::bopds::ds::Interference::FaceFace { f1, f2, curves, .. } = iv {
                        let (ta, tb) = if *f1 == dfi { (*f1, *f2) } else { (*f2, *f1) };
                        ta == dfi && tb == tfi && !curves.is_empty()
                    } else {
                        false
                    }
                });
                if shares_curve {
                    for ri in result_face_start..result_faces.len() {
                        seam_edges.push((ri, ds.faces[tfi].source_face_idx));
                    }
                }
            }
        }
    }

    // Assemble result BRep from split faces
    let brep = BRep {
        vertices: target.vertices.clone(),
        edges: target.edges.clone(),
        solids: vec![Solid {
            shells: vec![Shell {
                faces: result_faces,
            }],
        }],
        geom: target.geom.clone(),
    };

    ImprintResult { brep, seam_edges }
}

/// Split a single DS face by its intersection curves.
/// Shared with builder logic — produces a list of SubFace.
fn split_face_by_curves(ds: &DS, face_idx: usize) -> Vec<SubFace> {
    let face = &ds.faces[face_idx];
    let fi = &face.face_info;

    if fi.curves_in.is_empty() {
        let boundary = face
            .boundary_verts
            .iter()
            .map(|&vi| ds.vertices[vi].point)
            .collect();
        return vec![SubFace {
            boundary,
            surface: face.surface.clone(),
            normal: face.normal,
        }];
    }

    match &face.surface.clone() {
        Surface3::Plane(plane) => split_planar_face_simple(ds, face_idx, plane),
        _ => {
            // For curved surfaces: return whole face (splitting would require
            // the full BooleanBuilder machinery)
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| ds.vertices[vi].point)
                .collect();
            vec![SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
            }]
        }
    }
}

fn split_planar_face_simple(ds: &DS, face_idx: usize, plane: &Plane) -> Vec<SubFace> {
    use crate::inttools::edge_face::plane_local_basis;

    let face = &ds.faces[face_idx];
    let boundary_3d: Vec<DVec3> = face
        .boundary_verts
        .iter()
        .map(|&vi| ds.vertices[vi].point)
        .collect();

    let mut segments: Vec<(DVec3, DVec3)> = Vec::new();
    for &ci in &face.face_info.curves_in {
        let ic = &ds.intersection_curves[ci];
        let p0 = ds.vertices[ic.start_vertex].point;
        let p1 = ds.vertices[ic.end_vertex].point;
        if (p1 - p0).length_squared()
            > crate::tolerance::TOLERANCE_ABS * crate::tolerance::TOLERANCE_ABS
        {
            segments.push((p0, p1));
        }
    }

    if segments.is_empty() {
        return vec![SubFace {
            boundary: boundary_3d,
            surface: face.surface.clone(),
            normal: face.normal,
        }];
    }

    let (u_axis, v_axis) = plane_local_basis(plane);

    let project = |p: DVec3| -> [f64; 2] {
        let d = p - plane.origin;
        [d.dot(u_axis), d.dot(v_axis)]
    };
    let unproject = |uv: [f64; 2]| -> DVec3 { plane.origin + u_axis * uv[0] + v_axis * uv[1] };

    let mut polygons_2d: Vec<Vec<[f64; 2]>> =
        vec![boundary_3d.iter().map(|&p| project(p)).collect()];

    for (seg_a, seg_b) in &segments {
        let sa = project(*seg_a);
        let sb = project(*seg_b);
        let mut next: Vec<Vec<[f64; 2]>> = Vec::new();
        for poly in polygons_2d.drain(..) {
            let split = split_poly_2d(&poly, sa, sb);
            next.extend(split);
        }
        polygons_2d = next;
    }

    polygons_2d
        .into_iter()
        .filter(|p| p.len() >= 3)
        .map(|poly_2d| {
            let boundary: Vec<DVec3> = poly_2d.iter().map(|&uv| unproject(uv)).collect();
            SubFace {
                boundary,
                surface: face.surface.clone(),
                normal: face.normal,
            }
        })
        .collect()
}

/// Split a 2D polygon by a directed segment. Returns 1 polygon if no split, 2 if split.
fn split_poly_2d(poly: &[[f64; 2]], sa: [f64; 2], sb: [f64; 2]) -> Vec<Vec<[f64; 2]>> {
    let n = poly.len();
    let seg_dir = [sb[0] - sa[0], sb[1] - sa[1]];

    let signed_dist =
        |p: [f64; 2]| -> f64 { seg_dir[0] * (p[1] - sa[1]) - seg_dir[1] * (p[0] - sa[0]) };

    let sides: Vec<f64> = poly.iter().map(|&p| signed_dist(p)).collect();

    // Find crossings
    let mut crossings: Vec<(usize, [f64; 2])> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        let di = sides[i];
        let dj = sides[j];
        if di * dj < 0.0 {
            let t = di / (di - dj);
            let cx = poly[i][0] + t * (poly[j][0] - poly[i][0]);
            let cy = poly[i][1] + t * (poly[j][1] - poly[i][1]);
            crossings.push((i, [cx, cy]));
        }
    }

    if crossings.len() < 2 {
        return vec![poly.to_vec()];
    }

    // Use first two crossings
    let (i0, c0) = crossings[0];
    let (i1, c1) = crossings[1];

    let (ia, ib, ca, cb) = if i0 <= i1 {
        (i0, i1, c0, c1)
    } else {
        (i1, i0, c1, c0)
    };

    // Sub-poly A: [0..=ia] + ca + cb + [ib+1..]
    let mut sub_a = poly[..=ia].to_vec();
    sub_a.push(ca);
    sub_a.push(cb);
    sub_a.extend_from_slice(&poly[ib + 1..]);

    // Sub-poly B: [ia+1..=ib] + cb + ca
    let mut sub_b = poly[ia + 1..=ib].to_vec();
    sub_b.push(cb);
    sub_b.push(ca);

    let mut result = Vec::new();
    if sub_a.len() >= 3 {
        result.push(sub_a);
    }
    if sub_b.len() >= 3 {
        result.push(sub_b);
    }
    if result.is_empty() {
        result.push(poly.to_vec());
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// Gap / overlap detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect gaps and overlaps between two BReps.
///
/// For each pair of faces (one from each BRep) that are within `tolerance` of
/// each other, samples points on face A and measures the distance to surface B.
///
/// - Distance ∈ (0, tolerance]: **Gap**
/// - Distance ≈ 0 and normals anti-parallel: **SharedFace**
/// - Distance < 0 (interpenetration, estimated): **Overlap**
pub fn detect_gaps_overlaps(a: &BRep, b: &BRep, tolerance: f64) -> GapOverlapReport {
    let mut report = GapOverlapReport::default();

    // Flatten faces with their surface indices
    let faces_a = collect_faces_with_surfaces(a);
    let faces_b = collect_faces_with_surfaces(b);

    for (fa_idx, fa_pts, _fa_surf, fa_normal) in &faces_a {
        // Bounding box of face A
        let (a_min, a_max) = aabb(fa_pts);

        for (fb_idx, _fb_pts, fb_surf, fb_normal) in &faces_b {
            let (b_min, b_max) = {
                let fb_pts2 = collect_face_points(b, *fb_idx);
                aabb(&fb_pts2)
            };

            // AABB pre-filter: skip if clearly too far
            let gap_max = (b_min - a_max).max(a_min - b_max).max_element();
            if gap_max > tolerance * 2.0 + 1.0 {
                continue;
            }

            // Sample up to 5 points on face A and measure distance to surface B
            let samples = sample_face_points(fa_pts, 5);
            let mut max_dist: f64 = f64::NEG_INFINITY;
            let mut min_dist: f64 = f64::INFINITY;
            let mut closest_sample = fa_pts[0];

            for &sp in &samples {
                let proj = closest_point_on_surface(fb_surf, sp, 8);
                let d = proj.distance;
                if d < min_dist {
                    min_dist = d;
                    closest_sample = sp;
                }
                if d > max_dist {
                    max_dist = d;
                }
            }

            // Classify
            let normals_antiparallel = fa_normal.dot(*fb_normal) < -0.9;

            if min_dist.abs() < tolerance * 0.1 && normals_antiparallel {
                report.shared_faces.push((*fa_idx, *fb_idx));
            } else if min_dist > 0.0 && min_dist <= tolerance {
                report.gaps.push(Gap {
                    face_a: *fa_idx,
                    face_b: *fb_idx,
                    max_gap: max_dist,
                    sample_point: closest_sample,
                });
            } else if min_dist < -tolerance * 0.1 {
                report.overlaps.push(Overlap {
                    face_a: *fa_idx,
                    face_b: *fb_idx,
                    penetration_depth: -min_dist,
                });
            }
        }
    }

    report
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Collect (flat_face_idx, vertex_points, surface, normal) for every face in brep.
fn collect_faces_with_surfaces(brep: &BRep) -> Vec<(usize, Vec<DVec3>, Surface3, DVec3)> {
    let mut result = Vec::new();
    let mut flat_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let pts: Vec<DVec3> = face
                    .outer_wire
                    .edges
                    .iter()
                    .filter_map(|we| brep.edges.get(we.idx))
                    .map(|e| brep.vertices[e.start].point)
                    .collect();

                let surface = brep
                    .geom
                    .face_surface
                    .get(flat_idx)
                    .and_then(|&si| si)
                    .map(|si| brep.geom.surfaces[si].clone())
                    .unwrap_or_else(|| {
                        Surface3::Plane(Plane {
                            origin: DVec3::ZERO,
                            normal: face.normal,
                        })
                    });

                result.push((flat_idx, pts, surface, face.normal));
                flat_idx += 1;
            }
        }
    }
    result
}

fn collect_face_points(brep: &BRep, flat_idx: usize) -> Vec<DVec3> {
    let mut idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if idx == flat_idx {
                    return face
                        .outer_wire
                        .edges
                        .iter()
                        .filter_map(|we| brep.edges.get(we.idx))
                        .map(|e| brep.vertices[e.start].point)
                        .collect();
                }
                idx += 1;
            }
        }
    }
    vec![]
}

fn aabb(pts: &[DVec3]) -> (DVec3, DVec3) {
    if pts.is_empty() {
        return (DVec3::ZERO, DVec3::ZERO);
    }
    let mut mn = pts[0];
    let mut mx = pts[0];
    for &p in pts.iter().skip(1) {
        mn = mn.min(p);
        mx = mx.max(p);
    }
    (mn, mx)
}

/// Pick up to `n` evenly spaced sample points from a face boundary.
fn sample_face_points(pts: &[DVec3], n: usize) -> Vec<DVec3> {
    if pts.is_empty() {
        return vec![];
    }
    let step = (pts.len() as f64 / n as f64).ceil() as usize;
    let step = step.max(1);
    pts.iter().step_by(step).copied().collect()
}
