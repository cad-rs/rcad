use glam::DVec3;
use serde::{Deserialize, Serialize};

/// Geometric (analytic) model types: position, curve, surface, primitive descriptors.
///
/// This module describes *what shape is*.
pub mod geom;

/// Topology model types: vertex/edge/face/shell/solid incidence relationships.
///
/// This module describes *how things are connected*.
pub mod topology;

pub use geom::PrimitiveSolid;
pub use geom::{Curve3, Surface3};
pub use topology::{Edge, Face, Shell, Solid, Vertex, Wire};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GeomStore {
    pub curves: Vec<Curve3>,
    pub surfaces: Vec<Surface3>,
    /// Indexed by `BRep.edges` index; value is index into `curves`.
    pub edge_curve: Vec<Option<usize>>,
    /// Flattened face order across solids/shells; value is index into `surfaces`.
    pub face_surface: Vec<Option<usize>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub solids: Vec<Solid>,
    #[serde(default)]
    pub geom: GeomStore,
}

impl Default for BRep {
    fn default() -> Self {
        Self::new()
    }
}

impl BRep {
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            solids: Vec::new(),
            geom: GeomStore::default(),
        }
    }

    /// Creates a unit box B-Rep.
    ///
    /// Vertex layout:
    ///   0:(0,0,0)  1:(w,0,0)  2:(w,h,0)  3:(0,h,0)   <- front face (z=0)
    ///   4:(0,0,d)  5:(w,0,d)  6:(w,h,d)  7:(0,h,d)   <- back face  (z=d)
    pub fn create_box(width: f64, height: f64, depth: f64) -> Self {
        let (w, h, d) = (width, height, depth);

        let vertices = vec![
            Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: DVec3::new(w,   0.0, 0.0) }, // 1
            Vertex { point: DVec3::new(w,   h,   0.0) }, // 2
            Vertex { point: DVec3::new(0.0, h,   0.0) }, // 3
            Vertex { point: DVec3::new(0.0, 0.0, d  ) }, // 4
            Vertex { point: DVec3::new(w,   0.0, d  ) }, // 5
            Vertex { point: DVec3::new(w,   h,   d  ) }, // 6
            Vertex { point: DVec3::new(0.0, h,   d  ) }, // 7
        ];

        // 12 edges: 4 front + 4 back + 4 lateral
        let edges = vec![
            Edge { start: 0, end: 1 }, // 0  front-bottom
            Edge { start: 1, end: 2 }, // 1  front-right
            Edge { start: 2, end: 3 }, // 2  front-top
            Edge { start: 3, end: 0 }, // 3  front-left
            Edge { start: 4, end: 5 }, // 4  back-bottom
            Edge { start: 5, end: 6 }, // 5  back-right
            Edge { start: 6, end: 7 }, // 6  back-top
            Edge { start: 7, end: 4 }, // 7  back-left
            Edge { start: 0, end: 4 }, // 8  lateral-bl
            Edge { start: 1, end: 5 }, // 9  lateral-br
            Edge { start: 2, end: 6 }, // 10 lateral-tr
            Edge { start: 3, end: 7 }, // 11 lateral-tl
        ];

        let faces = vec![
            // Front  (z=0, normal -Z)
            Face { outer_wire: Wire { edges: vec![0,1,2,3] }, inner_wires: vec![], normal: DVec3::new(0.0, 0.0, -1.0), triangles: vec![[0,1,2],[0,2,3]] },
            // Back   (z=d, normal +Z)
            Face { outer_wire: Wire { edges: vec![4,5,6,7] }, inner_wires: vec![], normal: DVec3::new(0.0, 0.0,  1.0), triangles: vec![[5,4,7],[5,7,6]] },
            // Bottom (y=0, normal -Y)
            Face { outer_wire: Wire { edges: vec![0,9,4,8] }, inner_wires: vec![], normal: DVec3::new(0.0,-1.0, 0.0), triangles: vec![[0,1,5],[0,5,4]] },
            // Top    (y=h, normal +Y)
            Face { outer_wire: Wire { edges: vec![2,10,6,11] }, inner_wires: vec![], normal: DVec3::new(0.0, 1.0, 0.0), triangles: vec![[3,2,6],[3,6,7]] },
            // Left   (x=0, normal -X)
            Face { outer_wire: Wire { edges: vec![3,11,7,8] }, inner_wires: vec![], normal: DVec3::new(-1.0,0.0, 0.0), triangles: vec![[0,3,7],[0,7,4]] },
            // Right  (x=w, normal +X)
            Face { outer_wire: Wire { edges: vec![1,10,5,9] }, inner_wires: vec![], normal: DVec3::new( 1.0,0.0, 0.0), triangles: vec![[1,2,6],[1,6,5]] },
        ];

        BRep {
            vertices,
            edges,
            solids: vec![Solid { shells: vec![Shell { faces }] }],
            geom: GeomStore::default(),
        }
    }

    /// Creates a triangulated UV sphere centered at origin.
    pub fn create_sphere(radius: f64, u_segments: usize, v_segments: usize) -> Self {
        let u = u_segments.max(3);
        let v = v_segments.max(2);

        let mut points = Vec::with_capacity((u + 1) * (v + 1));
        for y in 0..=v {
            let phi = std::f64::consts::PI * y as f64 / v as f64;
            let sin_phi = phi.sin();
            let cos_phi = phi.cos();

            for x in 0..=u {
                let theta = 2.0 * std::f64::consts::PI * x as f64 / u as f64;
                points.push(DVec3::new(
                    radius * theta.cos() * sin_phi,
                    radius * cos_phi,
                    radius * theta.sin() * sin_phi,
                ));
            }
        }

        let mut triangles = Vec::with_capacity(u * v * 2);
        let stride = u + 1;
        for y in 0..v {
            for x in 0..u {
                let i0 = y * stride + x;
                let i1 = i0 + 1;
                let i2 = i0 + stride;
                let i3 = i2 + 1;

                if y != 0 {
                    triangles.push([i0, i2, i1]);
                }
                if y != v - 1 {
                    triangles.push([i1, i2, i3]);
                }
            }
        }

        Self::from_triangle_soup(points, triangles)
    }

    /// Creates a triangulated cylinder along Y axis, centered at origin.
    pub fn create_cylinder(radius: f64, height: f64, segments: usize) -> Self {
        let seg = segments.max(3);
        let half_h = height * 0.5;

        let mut points = Vec::with_capacity(seg * 2 + 2);
        for i in 0..seg {
            let t = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            let x = radius * t.cos();
            let z = radius * t.sin();
            points.push(DVec3::new(x, -half_h, z));
            points.push(DVec3::new(x, half_h, z));
        }
        let bottom_center = points.len();
        points.push(DVec3::new(0.0, -half_h, 0.0));
        let top_center = points.len();
        points.push(DVec3::new(0.0, half_h, 0.0));

        let mut triangles = Vec::with_capacity(seg * 4);
        for i in 0..seg {
            let next = (i + 1) % seg;
            let b0 = i * 2;
            let t0 = b0 + 1;
            let b1 = next * 2;
            let t1 = b1 + 1;

            triangles.push([b0, b1, t0]);
            triangles.push([t0, b1, t1]);
            triangles.push([bottom_center, b1, b0]);
            triangles.push([top_center, t0, t1]);
        }

        Self::from_triangle_soup(points, triangles)
    }

    /// Creates a triangulated cone along Y axis, apex at +Y.
    pub fn create_cone(base_radius: f64, height: f64, segments: usize) -> Self {
        let seg = segments.max(3);
        let half_h = height * 0.5;

        let mut points = Vec::with_capacity(seg + 2);
        for i in 0..seg {
            let t = 2.0 * std::f64::consts::PI * i as f64 / seg as f64;
            points.push(DVec3::new(base_radius * t.cos(), -half_h, base_radius * t.sin()));
        }
        let apex = points.len();
        points.push(DVec3::new(0.0, half_h, 0.0));
        let base_center = points.len();
        points.push(DVec3::new(0.0, -half_h, 0.0));

        let mut triangles = Vec::with_capacity(seg * 2);
        for i in 0..seg {
            let next = (i + 1) % seg;
            triangles.push([i, next, apex]);
            triangles.push([base_center, next, i]);
        }

        Self::from_triangle_soup(points, triangles)
    }

    /// Creates a triangulated torus around Y axis, centered at origin.
    pub fn create_torus(
        major_radius: f64,
        minor_radius: f64,
        major_segments: usize,
        minor_segments: usize,
    ) -> Self {
        let major = major_segments.max(3);
        let minor = minor_segments.max(3);

        let mut points = Vec::with_capacity(major * minor);
        for i in 0..major {
            let u = 2.0 * std::f64::consts::PI * i as f64 / major as f64;
            let cu = u.cos();
            let su = u.sin();

            for j in 0..minor {
                let v = 2.0 * std::f64::consts::PI * j as f64 / minor as f64;
                let cv = v.cos();
                let sv = v.sin();

                let ring = major_radius + minor_radius * cv;
                points.push(DVec3::new(ring * cu, minor_radius * sv, ring * su));
            }
        }

        let mut triangles = Vec::with_capacity(major * minor * 2);
        for i in 0..major {
            let ni = (i + 1) % major;
            for j in 0..minor {
                let nj = (j + 1) % minor;

                let a = i * minor + j;
                let b = ni * minor + j;
                let c = ni * minor + nj;
                let d = i * minor + nj;

                triangles.push([a, b, d]);
                triangles.push([b, c, d]);
            }
        }

        Self::from_triangle_soup(points, triangles)
    }

    pub fn from_primitive(primitive: PrimitiveSolid) -> Self {
        match primitive {
            PrimitiveSolid::Box {
                width,
                height,
                depth,
            } => Self::create_box(width, height, depth),
            PrimitiveSolid::Sphere {
                radius,
                u_segments,
                v_segments,
            } => Self::create_sphere(radius, u_segments, v_segments),
            PrimitiveSolid::Cylinder {
                radius,
                height,
                segments,
            } => Self::create_cylinder(radius, height, segments),
            PrimitiveSolid::Cone {
                base_radius,
                height,
                segments,
            } => Self::create_cone(base_radius, height, segments),
            PrimitiveSolid::Torus {
                major_radius,
                minor_radius,
                major_segments,
                minor_segments,
            } => Self::create_torus(major_radius, minor_radius, major_segments, minor_segments),
        }
    }

    pub fn center(&self) -> DVec3 {
        if self.vertices.is_empty() {
            return DVec3::ZERO;
        }
        let mut sum = DVec3::ZERO;
        for v in &self.vertices {
            sum += v.point;
        }
        sum / self.vertices.len() as f64
    }

    fn from_triangle_soup(points: Vec<DVec3>, triangles: Vec<[usize; 3]>) -> Self {
        let vertices = points.into_iter().map(|point| Vertex { point }).collect();
        let face = Face {
            outer_wire: Wire { edges: Vec::new() },
            inner_wires: Vec::new(),
            normal: DVec3::Z,
            triangles,
        };

        Self {
            vertices,
            edges: Vec::new(),
            solids: vec![Solid {
                shells: vec![Shell { faces: vec![face] }],
            }],
            geom: GeomStore::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle_count(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .map(|f| f.triangles.len())
            .sum()
    }

    #[test]
    fn creates_sphere_with_triangles() {
        let brep = BRep::create_sphere(1.0, 24, 12);
        assert!(!brep.vertices.is_empty());
        assert!(triangle_count(&brep) > 0);
    }

    #[test]
    fn creates_cylinder_with_triangles() {
        let brep = BRep::create_cylinder(1.0, 2.0, 24);
        assert!(!brep.vertices.is_empty());
        assert!(triangle_count(&brep) > 0);
    }

    #[test]
    fn creates_cone_with_triangles() {
        let brep = BRep::create_cone(1.0, 2.0, 24);
        assert!(!brep.vertices.is_empty());
        assert!(triangle_count(&brep) > 0);
    }

    #[test]
    fn creates_torus_with_triangles() {
        let brep = BRep::create_torus(1.0, 0.3, 24, 16);
        assert!(!brep.vertices.is_empty());
        assert!(triangle_count(&brep) > 0);
    }
}
