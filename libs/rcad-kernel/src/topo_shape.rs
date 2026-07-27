//! OCCT TopoDS_Shape: Arc<TShape> handle + Location + Orientation.
//! Self-contained -- no external ShapeStore needed.
//! Shape follows Arc<TShape> (like OCCT's TShape* handle).

use std::sync::Arc;
use glam::DVec3;
use crate::topods::{
    TVertexData, TEdgeData, TWireData, TFaceData, TShellData, TSolidData,
    Orientation, ShapeType, TShape,
};

// ---- Typed Shape wrappers (OCCT: TopoDS_Vertex, TopoDS_Edge, ...) ----

#[derive(Debug, Clone)]
pub struct Vertex(pub Shape);
impl Vertex {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Vertex); Vertex(s) }
}

#[derive(Debug, Clone)]
pub struct Edge(pub Shape);
impl Edge {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Edge); Edge(s) }
    pub fn tedge_data(&self) -> &TEdgeData {
        if let TShape::Edge(ref ed) = *self.0.data { ed } else { panic!("not an Edge") }
    }
}

#[derive(Debug, Clone)]
pub struct Wire(pub Shape);
impl Wire {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Wire); Wire(s) }
}

#[derive(Debug, Clone)]
pub struct Face(pub Shape);
impl Face {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Face); Face(s) }
    pub fn tface_data(&self) -> &TFaceData {
        if let TShape::Face(ref fd) = *self.0.data { fd } else { panic!("not a Face") }
    }
}

#[derive(Debug, Clone)]
pub struct Shell(pub Shape);
impl Shell {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Shell); Shell(s) }
}

#[derive(Debug, Clone)]
pub struct Solid(pub Shape);
impl Solid {
    pub fn shape(&self) -> &Shape { &self.0 }
    pub fn into_shape(self) -> Shape { self.0 }
    pub fn new(s: Shape) -> Self { assert!(s.shape_type() == ShapeType::Solid); Solid(s) }
}

#[derive(Debug, Clone)]
pub struct Shape {
    /// Arc<TShape> handle (OCCT: Handle(TShape) / TShape*).
    pub data: Arc<TShape>,
    /// TopLoc_Location index; 0 = identity.
    pub location: u32,
    /// TopAbs_Orientation.
    pub orientation: Orientation,
}

impl Shape {
    pub fn new(data: Arc<TShape>, location: u32, orientation: Orientation) -> Self {
        Shape { data, location, orientation }
    }

    pub fn shape_type(&self) -> ShapeType { self.data.shape_type() }

    // Type-safe accessors (like OCCT's TopoDS::Vertex(s)).
    pub fn as_vertex(&self) -> Option<&TVertexData> {
        if let TShape::Vertex(ref vd) = *self.data { Some(vd) } else { None }
    }
    pub fn as_edge(&self) -> Option<&TEdgeData> {
        if let TShape::Edge(ref ed) = *self.data { Some(ed) } else { None }
    }
    pub fn as_wire(&self) -> Option<&TWireData> {
        if let TShape::Wire(ref wd) = *self.data { Some(wd) } else { None }
    }
    pub fn as_face(&self) -> Option<&TFaceData> {
        if let TShape::Face(ref fd) = *self.data { Some(fd) } else { None }
    }
    pub fn as_shell(&self) -> Option<&TShellData> {
        if let TShape::Shell(ref sd) = *self.data { Some(sd) } else { None }
    }
    pub fn as_solid(&self) -> Option<&TSolidData> {
        if let TShape::Solid(ref sd) = *self.data { Some(sd) } else { None }
    }

    /// ptr_id for Hash/Eq (Arc pointer identity, like OCCT's TShape*).
    pub fn ptr_id(&self) -> u64 { Arc::as_ptr(&self.data) as u64 }

    pub fn is_null(&self) -> bool { false }
    pub fn is_vertex(&self) -> bool { self.shape_type() == ShapeType::Vertex }
    pub fn is_edge(&self) -> bool { self.shape_type() == ShapeType::Edge }
    pub fn is_wire(&self) -> bool { self.shape_type() == ShapeType::Wire }
    pub fn is_face(&self) -> bool { self.shape_type() == ShapeType::Face }
    pub fn is_shell(&self) -> bool { self.shape_type() == ShapeType::Shell }
    pub fn is_solid(&self) -> bool { self.shape_type() == ShapeType::Solid }
}

impl PartialEq for Shape {
    fn eq(&self, other: &Self) -> bool {
        Arc::as_ptr(&self.data) == Arc::as_ptr(&other.data)
            && self.location == other.location
            && self.orientation as u8 == other.orientation as u8
    }
}
impl Eq for Shape {}
impl std::hash::Hash for Shape {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        (Arc::as_ptr(&self.data) as u64).hash(state);
    }
}
