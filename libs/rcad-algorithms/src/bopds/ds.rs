use glam::DVec3;
use rcad_kernel::geom::*;
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
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub t_range: [f64; 2],
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

        ds
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
                                let non_shared = if shared == e.start {
                                    e.end
                                } else {
                                    e.start
                                };
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
}
