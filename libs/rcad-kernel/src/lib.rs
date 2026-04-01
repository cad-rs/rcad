use glam::DVec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Vertex {
    pub point: DVec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Wire {
    pub edges: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    pub outer_wire: Wire,
    pub inner_wires: Vec<Wire>,
    pub normal: DVec3,
    /// Pre-triangulated vertex index triples (into BRep.vertices)
    pub triangles: Vec<[usize; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    pub faces: Vec<Face>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    pub shells: Vec<Shell>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub vertices: Vec<Vertex>,
    pub edges: Vec<Edge>,
    pub solids: Vec<Solid>,
}

impl Default for BRep {
    fn default() -> Self {
        Self::new()
    }
}

impl BRep {
    pub fn new() -> Self {
        Self { vertices: Vec::new(), edges: Vec::new(), solids: Vec::new() }
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
}
