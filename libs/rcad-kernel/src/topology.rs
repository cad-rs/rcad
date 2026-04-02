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
