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

/// A vertex in the DS pool.
#[derive(Debug, Clone)]
pub struct DSVertex {
    pub point: DVec3,
    /// None for vertices created at intersections.
    pub origin: Option<ShapeOrigin>,
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

/// Central data structure (OCCT: BOPDS_DS).
#[derive(Debug)]
pub struct DS {
    pub vertices: Vec<DSVertex>,
    pub edges: Vec<DSEdge>,
    pub faces: Vec<DSFace>,
    pub interferences: Vec<Interference>,
    pub intersection_curves: Vec<IntersectionCurve>,
}

impl DS {
    /// Build DS from two BReps. Both must have populated GeomStore.
    pub fn new(a: &BRep, b: &BRep) -> Self {
        let mut ds = DS {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            interferences: Vec::new(),
            intersection_curves: Vec::new(),
        };

        ds.load_brep(a, ShapeOrigin::ShapeA);
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
        diagonal.max(1e-10)
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
            if matches!(self.faces[fi].surface, Surface3::Plane(_)) {
                continue; // Planar faces use existing 2D projection logic
            }

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
                    for ei in &boundary_edges {
                        let edge = &self.edges[*ei];
                        let [t0, t1] = edge.t_range;
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
                    // Add small margin
                    let margin = (h_max - h_min) * 0.01 + 1e-9;
                    // Use [0, 2π] to match CylindricalSurface::point_at parameterisation.
                    // circle_pcurve_on_cylinder also uses u ∈ [0, 2π], so the trim polyline
                    // will lie entirely inside this UV boundary.
                    let uv = vec![
                        DVec2::new(0.0, h_min - margin),
                        DVec2::new(2.0 * PI, h_min - margin),
                        DVec2::new(2.0 * PI, h_max + margin),
                        DVec2::new(0.0, h_max + margin),
                    ];
                    self.faces[fi].uv_boundary = Some(uv);
                    continue;
                }
                Surface3::Cone(cone) => {
                    // Cone param: u = azimuth [0, 2π], v = slant distance from apex (v ≥ 0).
                    // Estimate v range from boundary edge samples.
                    let boundary_edges = self.faces[fi].boundary_edges.clone();
                    let mut v_max = 0.0_f64;
                    let apex = cone.apex;
                    let axis = cone.axis.normalize();
                    let tan_h = cone.half_angle_rad.tan();
                    for ei in &boundary_edges {
                        let edge = &self.edges[*ei];
                        let [t0, t1] = edge.t_range;
                        for k in 0..=N_SAMPLES {
                            let t = t0 + (t1 - t0) * k as f64 / N_SAMPLES as f64;
                            let p = edge.curve.point_at(t);
                            let local = p - apex;
                            let along = local.dot(axis);
                            // slant distance s satisfies: z = s, r = s*tan(half)
                            // → s = along / (1 + tan²) + radial/tan ... approximate: s ≈ along
                            let radial_len = (local - axis * along).length();
                            let s = if tan_h > 1e-14 {
                                (along + radial_len / tan_h) * 0.5
                            } else {
                                along
                            };
                            v_max = v_max.max(s.max(0.0));
                        }
                    }
                    if v_max < 1e-9 {
                        v_max = 1.0;
                    }
                    let margin = v_max * 0.01 + 1e-9;
                    let uv = vec![
                        DVec2::new(0.0, 0.0),
                        DVec2::new(2.0 * PI, 0.0),
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

            self.faces[fi].uv_boundary = Some(uv_pts);
        }
    }

    fn load_brep(&mut self, brep: &BRep, origin: ShapeOrigin) {
        let vert_offset = self.vertices.len();
        let edge_offset = self.edges.len();

        // Vertices
        for v in &brep.vertices {
            self.vertices.push(DSVertex {
                point: v.point,
                origin: Some(origin),
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
                            // Fallback: synthesize plane from face normal and first triangle
                            let origin = brep.vertices[face.triangles[0][0]].point;
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
                        uv_boundary: None,
                    });

                    face_idx += 1;
                }
            }
        }
    }

    /// Add a vertex, deduplicating against existing vertices.
    pub fn add_vertex(&mut self, point: DVec3) -> usize {
        for (i, v) in self.vertices.iter().enumerate() {
            if points_coincide(v.point, point) {
                return i;
            }
        }
        let idx = self.vertices.len();
        self.vertices.push(DSVertex {
            point,
            origin: None,
        });
        idx
    }

    /// Collect 3D boundary points for a face.
    pub fn face_boundary_points(&self, face_idx: usize) -> Vec<DVec3> {
        self.faces[face_idx]
            .boundary_verts
            .iter()
            .map(|&vi| self.vertices[vi].point)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom_populate::populate_box_geom;
    use rcad_kernel::PrimitiveSolid;

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
}
