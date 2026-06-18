use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, *};
use rcad_kernel::{BRep, CurveEval};

use super::face_info::FaceInfo;
use super::pave::{Pave, PaveBlock};
use crate::tolerance::*;

/// Identifies which input shape a sub-shape came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeOrigin {
    ShapeA,
    ShapeB,
}

/// Information about shared topology between the two input shapes.
///
/// This is used by the glue path to skip interference detection for
/// sub-shapes that are already coincident between the two inputs.
#[derive(Debug, Clone, Default)]
pub struct SharedTopologyInfo {
    /// Pairs of vertex indices (v_a, v_b) that are coincident.
    /// v_a is from ShapeA (index < a_vertex_count), v_b from ShapeB.
    pub shared_vertices: Vec<(usize, usize)>,
    /// Pairs of edge indices (e_a, e_b) that share the same geometry.
    /// e_a is from ShapeA (index < a_edge_count), e_b from ShapeB.
    pub shared_edges: Vec<(usize, usize)>,
    /// Pairs of face indices (f_a, f_b) that have shared topology.
    /// This includes both fully-overlapping faces and faces with partial overlap.
    pub shared_faces: Vec<(usize, usize)>,
    /// Face pairs with full boundary overlap (can be skipped entirely).
    pub fully_glued_faces: Vec<(usize, usize)>,
    /// Face pairs with partial edge sharing.
    pub partially_glued_faces: Vec<(usize, usize)>,
}

/// A vertex in the DS pool.
#[derive(Debug, Clone)]
pub struct DSVertex {
    pub point: DVec3,
    /// None for vertices created at intersections.
    pub origin: Option<ShapeOrigin>,
    /// Model tolerance at this vertex (`vertex_tolerance` from source BRep when loaded;
    /// [`TOLERANCE_ABS`](crate::tolerance::TOLERANCE_ABS) for vertices added by the DS).
    pub geom_tol: f64,
}

/// An edge in the DS pool with curve reference.
#[derive(Debug, Clone)]
pub struct DSEdge {
    /// Index into DS.vertices.
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub curve: Curve3,
    /// Parametric range `[t_start, t_end]` on the curve.
    pub t_range: [f64; 2],
    pub origin: ShapeOrigin,
    /// Model tolerance on the source edge (`edge_tolerance` from BRep when loaded).
    pub geom_tol: f64,
    /// Paves inserted on this edge by intersection passes (unsorted until build_split_edges).
    pub paves: Vec<Pave>,
    /// After `build_split_edges`, the edge is represented by these sub-segments.
    pub pave_blocks: Vec<PaveBlock>,
}

/// A face in the DS pool with surface reference.
#[derive(Debug, Clone)]
pub struct DSFace {
    pub surface: Surface3,
    /// Boundary vertex indices (ordered, into DS.vertices) — outer wire.
    pub boundary_verts: Vec<usize>,
    /// Boundary edge indices (into DS.edges) — outer wire.
    pub boundary_edges: Vec<usize>,
    pub normal: DVec3,
    pub origin: ShapeOrigin,
    pub face_info: FaceInfo,
    /// Original face index within the source BRep's flattened face list.
    pub source_face_idx: usize,
    /// Model tolerance on the source face (`face_tolerance` from BRep when loaded).
    pub geom_tol: f64,
    /// UV-space boundary polygon on this face's surface (populated in Task 3+).
    pub uv_boundary: Option<Vec<DVec2>>,
}

/// Record of an intersection between two sub-shapes.
#[derive(Debug, Clone)]
pub enum Interference {
    VertexVertex {
        v1: usize,
        v2: usize,
        merged_vertex: usize,
    },
    VertexEdge {
        vertex: usize,
        edge: usize,
        param: f64,
    },
    EdgeEdge {
        e1: usize,
        e2: usize,
        point: DVec3,
        param1: f64,
        param2: f64,
        new_vertex: usize,
    },
    VertexFace {
        vertex: usize,
        face: usize,
    },
    EdgeFace {
        edge: usize,
        face: usize,
        point: DVec3,
        edge_param: f64,
        new_vertex: usize,
    },
    FaceFace {
        f1: usize,
        f2: usize,
        /// Intersection curve indices (into DS.intersection_curves).
        curves: Vec<usize>,
        /// Tangent touch point vertices.
        points: Vec<usize>,
    },
}

/// An intersection curve from F-F intersection, bounded by vertices.
#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    /// Sampled points from numerical marching (non-empty for marched curves).
    /// When non-empty this takes priority over `curve` for face splitting.
    pub polyline: Vec<DVec3>,
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub t_range: [f64; 2],
    /// PCurve (2D parametric curve) of this intersection on surface A (populated in Task 3+).
    pub pcurve_on_a: Option<Curve2d>,
    /// PCurve (2D parametric curve) of this intersection on surface B (populated in Task 3+).
    pub pcurve_on_b: Option<Curve2d>,
}

/// Information about detected extreme geometry conditions.
///
/// This stores the results of pre-analysis for near-tangent and near-coincident
/// geometry, enabling automatic tolerance adjustment during boolean operations.
#[derive(Debug, Clone, Default)]
pub struct ExtremeGeometryInfo {
    /// Near-tangent face pairs detected during pre-analysis.
    pub near_tangent_faces: Vec<NearTangentFacePair>,
    /// Near-coincident face pairs detected during pre-analysis.
    pub near_coincident_faces: Vec<NearCoincidentFacePair>,
    /// Recommended fuzzy tolerance adjustment.
    pub recommended_fuzzy_adjustment: f64,
    /// Whether extreme geometry was detected that requires special handling.
    pub has_extreme_geometry: bool,
}

/// A near-tangent face pair with detailed information.
#[derive(Debug, Clone)]
pub struct NearTangentFacePair {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Distance between faces at closest point.
    pub distance: f64,
    /// Type of near-tangency.
    pub tangent_type: NearTangentType,
    /// Suggested fuzzy tolerance for this pair.
    pub suggested_fuzzy: f64,
}

/// Type of near-tangency between faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearTangentType {
    /// Planes that are nearly parallel.
    PlaneParallel,
    /// Cylinder tangent to plane.
    CylinderPlane,
    /// Sphere tangent to plane.
    SpherePlane,
    /// Two cylinders tangent along a generator.
    CylinderCylinder,
    /// Cone tangent to plane.
    ConePlane,
    /// General surface tangency.
    General,
}

/// A near-coincident face pair with detailed information.
#[derive(Debug, Clone)]
pub struct NearCoincidentFacePair {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Maximum distance between faces in overlap region.
    pub max_distance: f64,
    /// Overlap ratio (0.0 to 1.0).
    pub overlap_ratio: f64,
    /// Suggested fuzzy tolerance for this pair.
    pub suggested_fuzzy: f64,
}

/// Central data structure (OCCT: BOPDS_DS).
#[derive(Debug)]
pub struct DS {
    pub vertices: Vec<DSVertex>,
    pub edges: Vec<DSEdge>,
    pub faces: Vec<DSFace>,
    pub interferences: Vec<Interference>,
    pub intersection_curves: Vec<IntersectionCurve>,
    /// Fuzzy tolerance used during interference detection.
    ///
    /// Vertices/edges within this distance are considered coincident.
    /// When set to a value larger than `TOLERANCE_ABS`, approximate
    /// near-miss intersections (analogous to OCCT `BOPAlgo_Options::SetFuzzyValue`).
    pub fuzzy_tol: f64,
    /// Number of vertices loaded from shape A (first shape). Shape A DS vertex indices are 0..a_vertex_count.
    pub a_vertex_count: usize,
    /// Number of edges loaded from shape A. Shape A DS edge indices are 0..a_edge_count.
    pub a_edge_count: usize,
    /// Number of faces loaded from shape A. Shape A DS face indices are 0..a_face_count.
    pub a_face_count: usize,
    /// Shared topology information for glue path optimization.
    pub shared_topology: SharedTopologyInfo,
    /// Extreme geometry analysis results.
    pub extreme_geometry: ExtremeGeometryInfo,
    /// Pre-computed overlap polygons for same-domain (coplanar) face pairs.
    /// Each entry is (face_a_index, face_b_index, overlap_boundary_in_3d).
    /// Populated during PaveFiller's coplanar analysis, consumed by Builder.
    pub same_domain_overlaps: Vec<(usize, usize, Vec<DVec3>)>,

    /// Edge image mapping (OCCT: BOPAlgo_Builder::myImages).
    /// Indexed by original edge index, each entry lists sub-edge indices
    /// created by `build_edge_images()`.
    pub my_images: Vec<Vec<usize>>,
    /// Edge origin mapping (OCCT: BOPAlgo_Builder::myOrigins).
    /// Indexed by sub-edge index, value is the original edge index.
    pub my_origins: Vec<usize>,
}

impl DS {
    /// ✅ OCCT-aligned: BRep_Tool::Degenerated(edge) equivalent.
    pub fn is_edge_degenerated(&self, edge_idx: usize) -> bool {
        self.edges[edge_idx].start_vertex == self.edges[edge_idx].end_vertex
    }

    /// Build DS from two BReps using the default absolute tolerance.
    pub fn new(a: &BRep, b: &BRep) -> Self {
        Self::new_with_fuzzy(a, b, crate::tolerance::TOLERANCE_ABS)
    }

    /// Build DS with a caller-supplied fuzzy tolerance.
    ///
    /// `fuzzy_tol` must be ≥ `TOLERANCE_ABS`; smaller values are clamped up.
    pub fn new_with_fuzzy(a: &BRep, b: &BRep, fuzzy_tol: f64) -> Self {
        let tol = fuzzy_tol.max(crate::tolerance::TOLERANCE_ABS);
        let mut ds = DS {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            interferences: Vec::new(),
            intersection_curves: Vec::new(),
            fuzzy_tol: tol,
            a_vertex_count: 0,
            a_edge_count: 0,
            a_face_count: 0,
            shared_topology: SharedTopologyInfo::default(),
            extreme_geometry: ExtremeGeometryInfo::default(),
            same_domain_overlaps: Vec::new(),
            my_images: Vec::new(),
            my_origins: Vec::new(),
        };

        ds.load_brep(a, ShapeOrigin::ShapeA);
        ds.a_vertex_count = ds.vertices.len();
        ds.a_edge_count = ds.edges.len();
        ds.a_face_count = ds.faces.len();
        ds.load_brep(b, ShapeOrigin::ShapeB);
        ds.compute_uv_boundaries();

        ds
    }

    /// Compute the characteristic scale of the model from all vertices.
    /// Returns the diagonal of the bounding box, or 1.0 if empty.
    pub fn model_scale(&self) -> f64 {
        use glam::DVec3;
        let mut min_pt = DVec3::splat(f64::INFINITY);
        let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
        let mut has_vertices = false;

        for v in &self.vertices {
            min_pt = min_pt.min(v.point);
            max_pt = max_pt.max(v.point);
            has_vertices = true;
        }

        if !has_vertices {
            return 1.0;
        }

        let diagonal = (max_pt - min_pt).length();
        diagonal.max(TOLERANCE_MODEL_SCALE_MIN)
    }

    /// Compute UV boundary for all curved faces by projecting 3D boundary
    /// points onto the face surface's parameter domain.
    ///
    /// For each boundary edge, we sample `N_SAMPLES` evenly-spaced points along
    /// the edge curve so that the resulting UV polygon is well-defined even when
    /// the wire has very few vertices (e.g. a sphere with only 2 poles).
    pub fn compute_uv_boundaries(&mut self) {
        use std::f64::consts::PI;
        const N_SAMPLES: usize = 8;

        for fi in 0..self.faces.len() {
            let surface = self.faces[fi].surface.clone();

            // For sphere and cylinder, the UV boundary is the full parameter
            // domain rectangle. The topological boundary (seam edge) maps to a
            // degenerate line in UV space and cannot be used as a polygon.
            match &surface {
                Surface3::Sphere(_) => {
                    // Sphere param from projection: u = longitude [-π, π] (atan2 range),
                    // v = colatitude [0, π]. Use the full domain as UV boundary.
                    let uv = vec![
                        DVec2::new(-PI, 0.0),
                        DVec2::new(PI, 0.0),
                        DVec2::new(PI, PI),
                        DVec2::new(-PI, PI),
                    ];
                    self.faces[fi].uv_boundary = Some(uv);
                    continue;
                }
                Surface3::Cylinder(cyl) => {
                    // Cylinder param: u = azimuth [0, 2π] (matches CylindricalSurface::point_at),
                    // v = height along axis.  Estimate height range from boundary edge samples.
                    let boundary_edges = self.faces[fi].boundary_edges.clone();
                    let mut h_min = f64::INFINITY;
                    let mut h_max = f64::NEG_INFINITY;
                    let axis = cyl.axis.normalize();
                    let origin = cyl.origin;
                    // Seam edges on a cylinder are lines parallel to the axis.  When a face has
                    // multiple seam edges (e.g. a cylinder with explicit front/back seams), the
                    // u-range is bounded by them rather than the full [0, 2π].
                    let u_ax = any_perpendicular(axis);
                    let v_ax = axis.cross(u_ax).normalize();
                    let mut seam_u_vals: Vec<f64> = Vec::new();
                    for ei in &boundary_edges {
                        let edge = &self.edges[*ei];
                        let [t0, t1] = edge.t_range;
                        // Detect seam edges: lines whose direction is parallel to the cylinder axis
                        if let Curve3::Line(line) = &edge.curve {
                            let dir = line.direction.normalize();
                            if dir.dot(axis).abs() > 1.0 - TOLERANCE_ABS {
                                let mid = edge.curve.point_at(0.5 * (t0 + t1));
                                let v_comp = (mid - origin).dot(axis);
                                let radial = mid - origin - v_comp * axis;
                                let u = radial.dot(v_ax).atan2(radial.dot(u_ax));
                                seam_u_vals.push(u);
                            }
                        }
                        // v-range sampling (same as before)
                        for k in 0..=N_SAMPLES {
                            let t = t0 + (t1 - t0) * k as f64 / N_SAMPLES as f64;
                            let p = edge.curve.point_at(t);
                            let h = (p - origin).dot(axis);
                            h_min = h_min.min(h);
                            h_max = h_max.max(h);
                        }
                    }
                    if !h_min.is_finite() || !h_max.is_finite() {
                        h_min = -1.0;
                        h_max = 1.0;
                    }
                    // Deduplicate seam u-values (same edge may appear fwd+rev in the wire)
                    seam_u_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    seam_u_vals.dedup_by(|a, b| (*a - *b).abs() < 1e-6);
                    let (u_lo, u_hi) = if seam_u_vals.len() >= 2 {
                        (seam_u_vals[0], seam_u_vals[seam_u_vals.len() - 1])
                    } else {
                        (0.0, 2.0 * PI)
                    };
                    // Add a small margin to the v-range so that intersection polyline
                    // endpoints near v=0 or v=h are not clipped by
                    // extend_trim_to_uv_boundary.  A 1% margin matches the cone case
                    // below and has worked well in practice.
                    let v_range = h_max - h_min;
                    let margin = if v_range > TOLERANCE_COORD_SUB {
                        v_range * 0.01 + TOLERANCE_COORD_SUB
                    } else {
                        TOLERANCE_COORD_SUB
                    };
                    let uv = vec![
                        DVec2::new(u_lo, h_min - margin),
                        DVec2::new(u_hi, h_min - margin),
                        DVec2::new(u_hi, h_max + margin),
                        DVec2::new(u_lo, h_max + margin),
                    ];
                    self.faces[fi].uv_boundary = Some(uv);
                    continue;
                }
                Surface3::Cone(cone) => {
                    // Cone param: u = azimuth [0, 2π], v = slant distance from the
                    // reference circle centered at `cone.apex`.
                    // Estimate the full slant range from boundary edge samples so
                    // reference-circle cones keep the correct UV window.
                    let boundary_edges = self.faces[fi].boundary_edges.clone();
                    let mut v_min = f64::INFINITY;
                    let mut v_max = f64::NEG_INFINITY;
                    let ref_point = cone.apex;
                    let axis = cone.axis_dir();
                    for ei in &boundary_edges {
                        let edge = &self.edges[*ei];
                        let [t0, t1] = edge.t_range;
                        for k in 0..=N_SAMPLES {
                            let t = t0 + (t1 - t0) * k as f64 / N_SAMPLES as f64;
                            let p = edge.curve.point_at(t);
                            let local = p - ref_point;
                            let along = local.dot(axis);
                            let slant = cone.slant_from_axial(along);
                            v_min = v_min.min(slant);
                            v_max = v_max.max(slant);
                        }
                    }
                    if !v_min.is_finite() || !v_max.is_finite() {
                        v_min = 0.0;
                        v_max = 1.0;
                    }
                    if (v_max - v_min).abs() < TOLERANCE_COORD_SUB {
                        v_min -= 0.5;
                        v_max += 0.5;
                    }
                    let margin = (v_max - v_min) * 0.01 + TOLERANCE_COORD_SUB;
                    let uv = vec![
                        DVec2::new(0.0, v_min - margin),
                        DVec2::new(2.0 * PI, v_min - margin),
                        DVec2::new(2.0 * PI, v_max + margin),
                        DVec2::new(0.0, v_max + margin),
                    ];
                    self.faces[fi].uv_boundary = Some(uv);
                    continue;
                }
                Surface3::Torus(_) => {
                    // Torus param: u = major angle [0, 2π], v = minor angle [0, 2π].
                    // Full parameter domain is always the UV boundary.
                    let uv = vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(2.0 * PI, 0.0),
                        DVec2::new(2.0 * PI, 2.0 * PI),
                        DVec2::new(0.0, 2.0 * PI),
                    ];
                    self.faces[fi].uv_boundary = Some(uv);
                    continue;
                }
                _ => {}
            }

            let boundary_edges = self.faces[fi].boundary_edges.clone();

            if boundary_edges.is_empty() {
                continue;
            }

            let mut pts_3d: Vec<DVec3> = Vec::new();
            for ei in &boundary_edges {
                let edge = &self.edges[*ei];
                let [t0, t1] = edge.t_range;
                for k in 0..N_SAMPLES {
                    let t = t0 + (t1 - t0) * (k as f64) / (N_SAMPLES as f64);
                    pts_3d.push(edge.curve.point_at(t));
                }
            }

            if pts_3d.is_empty() {
                continue;
            }

            let uv_pts: Vec<DVec2> = pts_3d
                .iter()
                .map(|&p| {
                    let proj = rcad_kernel::projection::closest_point_on_surface(&surface, p, 16);
                    DVec2::new(proj.params.0, proj.params.1)
                })
                .collect();

            // For planar surfaces (Plane and planar BSpline), decimate colinear
            // UV points. OCCT's BOPAlgo_BuilderFace uses the face's topological
            // edges directly (each edge has exactly 2 vertices), not sampled UV
            // polylines.  The 8-sample-per-edge UV boundary creates edge fragments.
            // Tolerance 1e-6 * edge_length accounts for projection noise from
            // closest_point_on_surface on BSpline surfaces (Newton iteration).
            let is_planar = matches!(&surface, Surface3::Plane(_))
                || (if let Surface3::BSpline(ref bsp) = surface { rcad_kernel::geom::bspline_is_planar(bsp, 1e-7) } else { false });
            let decimated = if is_planar {
                let n = uv_pts.len();
                if n > 2 {
                    let mut kept: Vec<DVec2> = Vec::with_capacity(n);
                    for i in 0..n {
                        let prev = uv_pts[(i + n - 1) % n];
                        let curr = uv_pts[i];
                        let next = uv_pts[(i + 1) % n];
                        let d1 = curr - prev;
                        let d2 = next - curr;
                        let cross = (d1.x * d2.y - d1.y * d2.x).abs();
                        let len1 = d1.length_squared();
                        let len2 = d2.length_squared();
                        if cross < 1e-6 * len1.max(len2).max(f64::MIN_POSITIVE)
                            && len1 > f64::MIN_POSITIVE && len2 > f64::MIN_POSITIVE
                        {
                            continue;
                        }
                        kept.push(curr);
                    }
                    if kept.len() >= 3 && kept.len() < uv_pts.len() {
                        kept
                    } else { uv_pts }
                } else { uv_pts }
            } else {
                uv_pts
            };

            self.faces[fi].uv_boundary = Some(decimated);
        }
    }

    fn load_brep(&mut self, brep: &BRep, origin: ShapeOrigin) {
        let vert_offset = self.vertices.len();
        let edge_offset = self.edges.len();

        // Vertices
        for (local_i, v) in brep.vertices.iter().enumerate() {
            self.vertices.push(DSVertex {
                point: v.point,
                origin: Some(origin),
                geom_tol: rcad_kernel::vertex_tolerance(brep, local_i),
            });
        }

        // Edges
        for (i, edge) in brep.edges.iter().enumerate() {
            let start = edge.start + vert_offset;
            let end = edge.end + vert_offset;

            let curve = brep
                .geom
                .edge_curve
                .get(i)
                .and_then(|c| *c)
                .map(|ci| brep.geom.curves[ci].clone())
                .unwrap_or_else(|| {
                    // Fallback: synthesize line from vertices
                    let p0 = brep.vertices[edge.start].point;
                    let p1 = brep.vertices[edge.end].point;
                    let dir = (p1 - p0).normalize();
                    Curve3::Line(Line3 {
                        origin: p0,
                        direction: dir,
                    })
                });

            // Compute parametric range
            let t_range = match &curve {
                Curve3::Line(line) => {
                    let p0 = brep.vertices[edge.start].point;
                    let p1 = brep.vertices[edge.end].point;
                    let t0 = (p0 - line.origin).dot(line.direction);
                    let t1 = (p1 - line.origin).dot(line.direction);
                    [t0, t1]
                }
                _ => brep
                    .geom
                    .edge_curve_range
                    .get(i)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| curve.default_domain()),
            };

            self.edges.push(DSEdge {
                start_vertex: start,
                end_vertex: end,
                curve,
                t_range,
                origin,
                geom_tol: rcad_kernel::edge_tolerance(brep, i),
                paves: Vec::new(),
                pave_blocks: Vec::new(),
            });
        }

        // Faces
        let mut face_idx = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let surface = brep
                        .geom
                        .face_surface
                        .get(face_idx)
                        .and_then(|s| *s)
                        .map(|si| brep.geom.surfaces[si].clone())
                        .unwrap_or_else(|| {
                            // Fallback: synthesize plane from face normal
                            // Use first vertex from outer wire, or origin if no wire
                            let origin = if !face.triangles.is_empty() {
                                brep.vertices[face.triangles[0][0]].point
                            } else if !face.outer_wire.edges.is_empty() {
                                let first_edge = &brep.edges[face.outer_wire.edges[0].idx];
                                brep.vertices[first_edge.start].point
                            } else {
                                DVec3::ZERO
                            };
                            Surface3::Plane(Plane {
                                origin,
                                normal: face.normal,
                            })
                        });

                    // Collect boundary vertices from wire edges
                    let boundary_edges: Vec<usize> = face
                        .outer_wire
                        .edges
                        .iter()
                        .map(|we| we.idx + edge_offset)
                        .collect();

                    // Trace the wire edges to get ordered boundary vertices.
                    // Wire edges are not necessarily in traversal order;
                    // we must find shared vertices between consecutive edges.
                    let boundary_verts: Vec<usize> = {
                        let edges_in_wire = &face.outer_wire.edges;
                        if edges_in_wire.is_empty() {
                            Vec::new()
                        } else if edges_in_wire.len() == 1 {
                            let e = &brep.edges[edges_in_wire[0].idx];
                            vec![e.start + vert_offset, e.end + vert_offset]
                        } else {
                            // For each consecutive pair of wire edges, find the
                            // shared vertex → the other vertex of the first edge
                            // is the boundary vertex contributed by that edge.
                            let mut verts = Vec::with_capacity(edges_in_wire.len());
                            for i in 0..edges_in_wire.len() {
                                let next_i = (i + 1) % edges_in_wire.len();
                                let e = &brep.edges[edges_in_wire[i].idx];
                                let en = &brep.edges[edges_in_wire[next_i].idx];

                                // The shared vertex between e and en
                                let shared = if e.start == en.start || e.start == en.end {
                                    e.start
                                } else {
                                    e.end
                                };

                                // The non-shared vertex of e is the boundary vertex
                                let non_shared = if shared == e.start { e.end } else { e.start };
                                verts.push(non_shared + vert_offset);
                            }
                            verts
                        }
                    };

                    self.faces.push(DSFace {
                        surface,
                        boundary_verts,
                        boundary_edges,
                        normal: face.normal,
                        origin,
                        face_info: FaceInfo::default(),
                        source_face_idx: face_idx,
                        geom_tol: rcad_kernel::face_tolerance(brep, face_idx),
                        uv_boundary: None,
                    });

                    face_idx += 1;
                }
            }
        }
    }

    /// ✅ OCCT-aligned: find existing vertex within tolerance (PutPaveOnCurve equivalent).
    /// OCCT's IsVertexOnLine checks if a boundary vertex lies on the intersection
    /// curve, then places the EXISTING vertex index on the curve's pave block,
    /// ensuring the section edge reuses the same TopoDS_Vertex.  This tolerance-
    /// based scan achieves the same sharing for rcad's flat vertex array.
    pub fn find_vertex_near(&self, point: DVec3, tol: f64) -> Option<usize> {
        let tol2 = tol * tol;
        self.vertices.iter().position(|v| (v.point - point).length_squared() <= tol2)
    }

    /// Add a vertex, deduplicating against existing vertices.
    ///
    /// Coincidence uses `max(fuzzy_tol, TOLERANCE_ABS, each vertex's geom_tol)` so
    /// imported vertex tolerances and pave fuzzy both widen merging.
    pub fn add_vertex(&mut self, point: DVec3) -> usize {
        let new_base = self.fuzzy_tol.max(TOLERANCE_ABS);
        for (i, v) in self.vertices.iter().enumerate() {
            let merge_tol = new_base.max(v.geom_tol);
            if (v.point - point).length() <= merge_tol {
                return i;
            }
        }
        let idx = self.vertices.len();
        self.vertices.push(DSVertex {
            point,
            origin: None,
            geom_tol: new_base,
        });
        idx
    }

    /// Collect 3D boundary points for a face.
    ///
    /// When the topological wire produces a degenerate polygon (< 3 unique points),
    /// falls back to the UV boundary if available and the face is planar. This ensures
    /// e.g. cylinder caps (2 boundary vertices on a diameter line) produce a proper
    /// circular polygon for face-face intersection clipping.
    pub fn face_boundary_points(&self, face_idx: usize) -> Vec<DVec3> {
        let face = &self.faces[face_idx];
        let pts: Vec<DVec3> = face
            .boundary_verts
            .iter()
            .map(|&vi| self.vertices[vi].point)
            .collect();
        if pts.len() < 3 {
            if let Some(uv_bnd) = &face.uv_boundary {
                if uv_bnd.len() >= 3 && matches!(face.surface, Surface3::Plane(_)) {
                    return uv_bnd
                        .iter()
                        .map(|uv| face.surface.point_at(uv.x, uv.y))
                        .collect();
                }
            }
        }
        pts
    }

    /// Detect shared topology between ShapeA and ShapeB.
    ///
    /// This method populates `self.shared_topology` with information about
    /// coincident vertices, edges, and faces. It should be called after
    /// the DS is fully constructed but before interference detection.
    ///
    /// # Arguments
    /// * `tolerance` - Base distance for glue-style coincidence; combined per pair with
    ///   relevant `geom_tol` on vertices, edges, or faces (`max` of all).
    ///
    /// # Returns
    /// A reference to the populated `SharedTopologyInfo`.
    pub fn detect_shared_topology(&mut self, tolerance: f64) -> &SharedTopologyInfo {
        let tol = tolerance.max(TOLERANCE_ABS);

        // Clear any previous data
        self.shared_topology = SharedTopologyInfo::default();

        // Detect shared vertices
        for vi_a in 0..self.a_vertex_count {
            for vi_b in self.a_vertex_count..self.vertices.len() {
                let p_a = self.vertices[vi_a].point;
                let p_b = self.vertices[vi_b].point;
                let pair_tol = tol
                    .max(self.vertices[vi_a].geom_tol)
                    .max(self.vertices[vi_b].geom_tol);
                let pair_tol_sq = pair_tol * pair_tol;
                if (p_a - p_b).length_squared() <= pair_tol_sq {
                    self.shared_topology.shared_vertices.push((vi_a, vi_b));
                }
            }
        }

        // Detect shared edges
        for ei_a in 0..self.a_edge_count {
            for ei_b in self.a_edge_count..self.edges.len() {
                let edge_tol = tol
                    .max(self.edges[ei_a].geom_tol)
                    .max(self.edges[ei_b].geom_tol);
                if self.edges_geometry_compatible(ei_a, ei_b, edge_tol) {
                    self.shared_topology.shared_edges.push((ei_a, ei_b));
                }
            }
        }

        // Detect shared faces (full and partial)
        for fi_a in 0..self.a_face_count {
            for fi_b in self.a_face_count..self.faces.len() {
                if self.faces[fi_a].origin == self.faces[fi_b].origin {
                    continue; // Same shape, skip
                }

                let face_tol = tol
                    .max(self.faces[fi_a].geom_tol)
                    .max(self.faces[fi_b].geom_tol);
                let full_overlap = self.faces_boundary_fully_overlap(fi_a, fi_b, face_tol);
                let partial_overlap =
                    !full_overlap && self.faces_share_edges(fi_a, fi_b, face_tol);

                if full_overlap {
                    self.shared_topology.fully_glued_faces.push((fi_a, fi_b));
                    self.shared_topology.shared_faces.push((fi_a, fi_b));
                } else if partial_overlap {
                    self.shared_topology.partially_glued_faces.push((fi_a, fi_b));
                    self.shared_topology.shared_faces.push((fi_a, fi_b));
                }
            }
        }

        &self.shared_topology
    }

    /// Check if two edges have compatible geometry.
    fn edges_geometry_compatible(&self, e1: usize, e2: usize, tol: f64) -> bool {
        let edge1 = &self.edges[e1];
        let edge2 = &self.edges[e2];

        // Check curve compatibility
        match (&edge1.curve, &edge2.curve) {
            (Curve3::Line(l1), Curve3::Line(l2)) => {
                // Lines must be collinear
                let d1 = l1.direction.normalize_or_zero();
                let d2 = l2.direction.normalize_or_zero();
                if d1.dot(d2).abs() < 0.999 {
                    return false;
                }
                // Origins must be on the same line
                let v = l2.origin - l1.origin;
                let perp = v - d1 * v.dot(d1);
                if perp.length() > tol {
                    return false;
                }
                // Check parameter ranges overlap
                let _p1_start = l1.origin + d1 * edge1.t_range[0];
                let _p1_end = l1.origin + d1 * edge1.t_range[1];
                let p2_start = l2.origin + d2 * edge2.t_range[0];
                let p2_end = l2.origin + d2 * edge2.t_range[1];

                // Project edge2 endpoints onto edge1 line
                let t2_start = (p2_start - l1.origin).dot(d1);
                let t2_end = (p2_end - l1.origin).dot(d1);
                let (t2_min, t2_max) = if t2_start < t2_end {
                    (t2_start, t2_end)
                } else {
                    (t2_end, t2_start)
                };

                // Check for overlap
                t2_min <= edge1.t_range[1] + tol && t2_max >= edge1.t_range[0] - tol
            }
            (Curve3::Circle(c1), Curve3::Circle(c2)) => {
                // Circles must be same radius and coplanar
                (c1.center - c2.center).length() <= tol
                    && c1.normal.dot(c2.normal).abs() >= 0.999
                    && (c1.radius - c2.radius).abs() <= tol
            }
            (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
                // Ellipses must be same
                (e1.center - e2.center).length() <= tol
                    && e1.normal.dot(e2.normal).abs() >= 0.999
                    && (e1.major_radius - e2.major_radius).abs() <= tol
                    && (e1.minor_radius - e2.minor_radius).abs() <= tol
            }
            _ => false,
        }
    }

    /// Check if two faces have fully overlapping boundaries.
    fn faces_boundary_fully_overlap(&self, f1: usize, f2: usize, tol: f64) -> bool {
        let pts1 = self.face_boundary_points(f1);
        let pts2 = self.face_boundary_points(f2);

        if pts1.len() < 3 || pts2.len() < 3 {
            return false;
        }

        // Each point in pts1 must have a matching point in pts2
        let tol_sq = tol * tol;
        let mut used = vec![false; pts2.len()];

        for p1 in &pts1 {
            let mut found = false;
            for (j, p2) in pts2.iter().enumerate() {
                if used[j] {
                    continue;
                }
                if (*p1 - *p2).length_squared() <= tol_sq {
                    used[j] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
        }

        true
    }

    /// Check if two faces share any edges.
    fn faces_share_edges(&self, f1: usize, f2: usize, tol: f64) -> bool {
        let edges1: std::collections::HashSet<usize> =
            self.faces[f1].boundary_edges.iter().copied().collect();
        let edges2: std::collections::HashSet<usize> =
            self.faces[f2].boundary_edges.iter().copied().collect();

        // Check for geometry-compatible edges
        for &e1 in &edges1 {
            for &e2 in &edges2 {
                if self.edges_geometry_compatible(e1, e2, tol) {
                    return true;
                }
            }
        }

        false
    }

    /// Check if a face pair is fully glued (can skip intersection entirely).
    pub fn is_fully_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
        self.shared_topology
            .fully_glued_faces
            .iter()
            .any(|&(a, b)| (a == f1 && b == f2) || (a == f2 && b == f1))
    }

    /// Check if a face pair is partially glued (has shared edges).
    pub fn is_partially_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
        self.shared_topology
            .partially_glued_faces
            .iter()
            .any(|&(a, b)| (a == f1 && b == f2) || (a == f2 && b == f1))
    }

    /// Get shared vertices for a face pair.
    pub fn get_shared_vertices_for_faces(&self, f1: usize, f2: usize) -> Vec<(usize, usize)> {
        let boundary1: std::collections::HashSet<usize> =
            self.faces[f1].boundary_verts.iter().copied().collect();
        let boundary2: std::collections::HashSet<usize> =
            self.faces[f2].boundary_verts.iter().copied().collect();

        self.shared_topology
            .shared_vertices
            .iter()
            .filter(|(v1, v2)| {
                (boundary1.contains(v1) && boundary2.contains(v2))
                    || (boundary1.contains(v2) && boundary2.contains(v1))
            })
            .copied()
            .collect()
    }

    /// ✅ OCCT-aligned: Build edge images from pave blocks (BOPAlgo_Builder::FillImagesEdges).
    ///
    /// For edges that were split by PaveFiller (pave_blocks.len() > 1), create sub-edges
    /// in `self.edges` and populate `my_images` / `my_origins` mappings.
    ///
    /// This must be called after `build_split_edges()` (end of `make_blocks`).
    pub fn build_edge_images(&mut self) {
        let n_edges = self.edges.len();
        self.my_images = vec![Vec::new(); n_edges];
        self.my_origins = Vec::new();

        // Pre-collect edge data to avoid borrow conflict with self.edges.push
        struct EdgeData {
            curve: Curve3,
            origin: ShapeOrigin,
            geom_tol: f64,
            blocks: Vec<(usize, usize, f64, f64)>, // (sv, ev, t_start, t_end)
        }
        let edge_data: Vec<EdgeData> = self.edges.iter().map(|e| {
            let blocks = e.pave_blocks.iter()
                .filter(|pb| pb.pave1.vertex_idx != pb.pave2.vertex_idx)
                .map(|pb| {
                    (pb.pave1.vertex_idx, pb.pave2.vertex_idx, pb.pave1.param, pb.pave2.param)
                })
                .collect();
            EdgeData {
                curve: e.curve.clone(),
                origin: e.origin,
                geom_tol: e.geom_tol,
                blocks,
            }
        }).collect();

        for ei in 0..n_edges {
            if edge_data[ei].blocks.is_empty() {
                continue;
            }
            let data = &edge_data[ei];
            for &(sv, ev, t_start, t_end) in &data.blocks {
                let sub_ei = self.edges.len();
                self.edges.push(DSEdge {
                    start_vertex: sv,
                    end_vertex: ev,
                    curve: data.curve.clone(),
                    t_range: [t_start, t_end],
                    origin: data.origin,
                    geom_tol: data.geom_tol,
                    paves: Vec::new(),
                    pave_blocks: Vec::new(),
                });
                self.my_images[ei].push(sub_ei);
                self.my_origins.push(ei);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom_populate::populate_box_geom;
    use rcad_kernel::tolerance::{set_edge_tolerance, set_face_tolerance, set_vertex_tolerance};
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn ds_load_brep_copies_geom_tolerances_into_pool() {
        let mut a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        populate_box_geom(&mut a);
        set_vertex_tolerance(&mut a, 0, 2.0 * TOLERANCE_ADAPTIVE_MAX);
        set_edge_tolerance(&mut a, 0, 3e-3);
        set_face_tolerance(&mut a, 0, 4e-3);

        let b = BRep::default();
        let ds = DS::new(&a, &b);

        assert!((ds.vertices[0].geom_tol - 2.0 * TOLERANCE_ADAPTIVE_MAX).abs() < TOLERANCE_LEN_MIN);
        assert!((ds.edges[0].geom_tol - 3e-3).abs() < TOLERANCE_LEN_MIN);
        assert!((ds.faces[0].geom_tol - 4e-3).abs() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn ds_add_vertex_dedup_respects_fuzzy_tol() {
        let empty = BRep::default();
        let mut ds = DS::new_with_fuzzy(&empty, &empty, TOLERANCE_RETRY_LADDER_COARSE);
        let a = ds.add_vertex(DVec3::ZERO);
        let b = ds.add_vertex(DVec3::new(5e-5, 0.0, 0.0));
        assert_eq!(a, b);
        assert!((ds.vertices[a].geom_tol - TOLERANCE_RETRY_LADDER_COARSE).abs() < TOLERANCE_FLOAT_ULTRA);
    }

    #[test]
    fn ds_from_two_boxes() {
        let mut a = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let mut b = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        populate_box_geom(&mut a);
        populate_box_geom(&mut b);

        let ds = DS::new(&a, &b);
        assert_eq!(ds.vertices.len(), 16); // 8 + 8
        assert_eq!(ds.edges.len(), 24); // 12 + 12
        assert_eq!(ds.faces.len(), 12); // 6 + 6

        // Check origin tags
        assert!(ds.vertices[0].origin == Some(ShapeOrigin::ShapeA));
        assert!(ds.vertices[8].origin == Some(ShapeOrigin::ShapeB));
        assert!(ds.edges[0].origin == ShapeOrigin::ShapeA);
        assert!(ds.edges[12].origin == ShapeOrigin::ShapeB);
    }

    #[test]
    fn ds_sphere_has_uv_boundary() {
        use rcad_modeling::make_sphere_brep;

        let a = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
        let ds = DS::new(&a, &b);

        // Sphere faces should have uv_boundary computed
        let sphere_faces: Vec<_> = ds
            .faces
            .iter()
            .filter(|f| matches!(f.surface, Surface3::Sphere(_)))
            .collect();
        assert!(!sphere_faces.is_empty(), "should have sphere faces");
        for f in &sphere_faces {
            assert!(
                f.uv_boundary.is_some(),
                "sphere face should have uv_boundary"
            );
            let uv = f.uv_boundary.as_ref().unwrap();
            assert!(uv.len() >= 3, "uv boundary should have at least 3 points");
        }
    }

    #[test]
    fn ds_cone_uv_boundary_uses_reference_circle_slant_range() {
        use rcad_modeling::make_cone_brep;

        let a = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).unwrap();
        let b = BRep::default();
        let ds = DS::new(&a, &b);

        let cone_face = ds
            .faces
            .iter()
            .find(|face| matches!(face.surface, Surface3::Cone(_)))
            .expect("should have a cone face");
        let uv = cone_face
            .uv_boundary
            .as_ref()
            .expect("cone face should have uv_boundary");

        let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        assert!(v_min < 0.0, "expected apex-side slant range below the reference circle, got {v_min}");
        assert!(v_max > 0.0, "expected base-side slant range above the reference circle, got {v_max}");
    }
}
