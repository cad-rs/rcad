//! OCCT Shape + TShape architecture: flat per-type storage.
//!
//! Shape = { kind, index, location, orientation }  -- 24 bytes, cache-friendly.
//! ShapeStore = per-type Vecs for compact, cache-efficient iteration.
//!
//! Unlike TShape enum (each element = max variant size, wastes cache),
//! per-type Vecs keep vertex data contiguous, edge data contiguous, etc.

use glam::DVec3;
use crate::geom::{Curve2d, Curve3, Surface3};
use crate::topods::{
    TVertexData, TEdgeData, TWireData, TFaceData, TShellData, TSolidData,
    Orientation, ShapeType, ShapeRef,
};

use std::sync::Arc;

// ======================================================================
// Shape -- compact TopoDS_Shape equivalent (24 bytes)
// ======================================================================

#[derive(Debug, Clone, Copy)]
pub struct Shape {
    pub kind: ShapeType,
    pub index: usize,
    pub location: u32,
    pub orientation: Orientation,
}

impl Shape {
    pub fn shape_type(&self) -> ShapeType { self.kind }

    /// Access typed data through ShapeStore.
    pub fn vertex<'a>(&self, s: &'a ShapeStore) -> Option<&'a TVertexData> {
        if self.kind != ShapeType::Vertex { return None; }
        s.vertices.get(self.index)
    }
    pub fn edge<'a>(&self, s: &'a ShapeStore) -> Option<&'a TEdgeData> {
        if self.kind != ShapeType::Edge { return None; }
        s.edges.get(self.index)
    }
    pub fn wire<'a>(&self, s: &'a ShapeStore) -> Option<&'a TWireData> {
        if self.kind != ShapeType::Wire { return None; }
        s.wires.get(self.index)
    }
    pub fn face<'a>(&self, s: &'a ShapeStore) -> Option<&'a TFaceData> {
        if self.kind != ShapeType::Face { return None; }
        s.faces.get(self.index)
    }
    pub fn shell<'a>(&self, s: &'a ShapeStore) -> Option<&'a TShellData> {
        if self.kind != ShapeType::Shell { return None; }
        s.shells.get(self.index)
    }
    pub fn solid<'a>(&self, s: &'a ShapeStore) -> Option<&'a TSolidData> {
        if self.kind != ShapeType::Solid { return None; }
        s.solids.get(self.index)
    }
}

// Identity by (ptr_id from Arc pointer, location, orientation).
impl PartialEq for Shape {
    fn eq(&self, other: &Self) -> bool {
        self.kind == other.kind
            && self.index == other.index
            && self.location == other.location
            && self.orientation as u8 == other.orientation as u8
    }
}
impl Eq for Shape {}

// ======================================================================
// ShapeStore -- per-type flat arrays
// ======================================================================

#[derive(Debug, Clone)]
pub struct ShapeStore {
    pub vertices: Vec<TVertexData>,
    pub edges: Vec<TEdgeData>,
    pub wires: Vec<TWireData>,
    pub faces: Vec<TFaceData>,
    pub shells: Vec<TShellData>,
    pub solids: Vec<TSolidData>,
    pub compsolids: Vec<Vec<ShapeRef>>,
    pub compounds: Vec<Vec<ShapeRef>>,
    pub locations: Vec<glam::DAffine3>,
}

impl ShapeStore {
    pub fn new() -> Self {
        ShapeStore {
            vertices: Vec::new(), edges: Vec::new(), wires: Vec::new(),
            faces: Vec::new(), shells: Vec::new(), solids: Vec::new(),
            compsolids: Vec::new(), compounds: Vec::new(),
            locations: Vec::new(),
        }
    }

    pub fn add_vertex(&mut self, v: TVertexData) -> Shape {
        let idx = self.vertices.len();
        self.vertices.push(v);
        Shape { kind: ShapeType::Vertex, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_edge(&mut self, e: TEdgeData) -> Shape {
        let idx = self.edges.len();
        self.edges.push(e);
        Shape { kind: ShapeType::Edge, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_face(&mut self, f: TFaceData) -> Shape {
        let idx = self.faces.len();
        self.faces.push(f);
        Shape { kind: ShapeType::Face, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_wire(&mut self, w: TWireData) -> Shape {
        let idx = self.wires.len();
        self.wires.push(w);
        Shape { kind: ShapeType::Wire, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_shell(&mut self, s: TShellData) -> Shape {
        let idx = self.shells.len();
        self.shells.push(s);
        Shape { kind: ShapeType::Shell, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_solid(&mut self, s: TSolidData) -> Shape {
        let idx = self.solids.len();
        self.solids.push(s);
        Shape { kind: ShapeType::Solid, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_compsolid(&mut self, c: Vec<ShapeRef>) -> Shape {
        let idx = self.compsolids.len();
        self.compsolids.push(c);
        Shape { kind: ShapeType::CompSolid, index: idx, location: 0, orientation: Orientation::Forward }
    }
    pub fn add_compound(&mut self, c: Vec<ShapeRef>) -> Shape {
        let idx = self.compounds.len();
        self.compounds.push(c);
        Shape { kind: ShapeType::Compound, index: idx, location: 0, orientation: Orientation::Forward }
    }

    pub fn total(&self) -> usize {
        self.vertices.len() + self.edges.len() + self.wires.len()
            + self.faces.len() + self.shells.len() + self.solids.len()
            + self.compsolids.len() + self.compounds.len()
    }
}

impl Default for ShapeStore { fn default() -> Self { Self::new() } }
