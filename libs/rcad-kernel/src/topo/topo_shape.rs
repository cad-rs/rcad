//! OCCT TopoDS_Shape: Arc<TShape> handle + Location + Orientation.
//! Self-contained -- no external ShapeStore needed.
//! Shape follows Arc<TShape> (like OCCT's TShape* handle).
//!
//! This version replaces the old ShapeRef: Shape now carries an `index` field
//! for O(1) BRep.tshapes[] access.

use std::sync::Arc;
use glam::DVec3;
use crate::topo::topods::{
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

/// TopoDS_Shape equivalent: Arc<TShape> pointer identity + Orientation + Index.
///
/// Replaces the old ShapeRef:
/// - `data`: Arc<TShape> (real identity, OCCT TShape* handle)
/// - `index`: O(1) access into BRep.tshapes[], or usize::MAX for synthetic/null shapes
/// - `location`: TopLoc_Location index; 0 = identity
/// - `orientation`: TopAbs_Orientation
#[derive(Debug, Clone)]
pub struct Shape {
    /// Arc<TShape> handle (OCCT: Handle(TShape) / TShape*).
    pub data: Arc<TShape>,
    /// Index into BRep.tshapes[] for O(1) access.
    /// usize::MAX for null/synthetic shapes without a real TShape in BRep.
    pub index: usize,
    /// TopLoc_Location index; 0 = identity.
    pub location: u32,
    /// TopAbs_Orientation.
    pub orientation: Orientation,
}

impl Shape {
    /// Create a Shape with no BRep index (index = usize::MAX).
    pub fn new(data: Arc<TShape>, location: u32, orientation: Orientation) -> Self {
        Shape { data, location, orientation, index: usize::MAX }
    }

    /// Create a Shape with a known BRep index.
    pub fn from_parts(data: Arc<TShape>, index: usize, location: u32, orientation: Orientation) -> Self {
        Shape { data, index, location, orientation }
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

    /// OCCT TopoDS_Shape::IsNull -- true if this is a null/uninitialized shape.
    pub fn is_null(&self) -> bool { self.index == usize::MAX && self.ptr_id() == 0 }

    /// OCCT TopoDS_Shape::IsSame -- same TShape (ignores Location and Orientation).
    pub fn is_same(&self, other: &Shape) -> bool {
        self.ptr_id() == other.ptr_id()
    }

    /// OCCT TopoDS_Shape::IsPartner -- same TShape AND same Location.
    pub fn is_partner(&self, other: &Shape) -> bool {
        self.ptr_id() == other.ptr_id() && self.location == other.location
    }

    /// OCCT TopoDS_Shape::IsEqual -- same TShape, Location, AND Orientation.
    pub fn is_equal(&self, other: &Shape) -> bool {
        self.ptr_id() == other.ptr_id()
            && self.location == other.location
            && self.orientation as u8 == other.orientation as u8
    }

    /// Returns a copy with a different TopLoc_Location index.
    pub fn with_location(self, location: u32) -> Self {
        Shape { location, ..self }
    }

    /// Null/uninitialized shape (OCCT: TopoDS_Shape() default constructor).
    pub fn null() -> Self {
        use crate::topo::topods::tshape_flags;
        Shape {
            data: Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE | tshape_flags::CLOSED | tshape_flags::CONVEX,
                point: DVec3::ZERO,
                tolerance: 0.0,
                points: Vec::new(),
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// Synthetic shape from a flat index (for DS adaptor code).
    /// The index is stored in the `index` field.
    pub fn synthetic(index: usize, orientation: Orientation) -> Self {
        use crate::topo::topods::tshape_flags;
        Shape {
            data: Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE | tshape_flags::CLOSED | tshape_flags::CONVEX,
                point: DVec3::ZERO,
                tolerance: 0.0,
                points: Vec::new(),
            })),
            index,
            location: 0,
            orientation,
        }
    }

    /// Synthetic shape from a flat index with location (for DS adaptor code).
    pub fn synthetic_with_location(index: usize, orientation: Orientation, location: u32) -> Self {
        use crate::topo::topods::tshape_flags;
        Shape {
            data: Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE | tshape_flags::CLOSED | tshape_flags::CONVEX,
                point: DVec3::ZERO,
                tolerance: 0.0,
                points: Vec::new(),
            })),
            index,
            location,
            orientation,
        }
    }

    /// OCCT TopoDS_Shape::ShapeType -- returns shape type from BRep.
    pub fn shape_type_from_brep(&self, brep: &crate::topo::topods::BRep) -> ShapeType {
        if self.is_null() {
            return ShapeType::Shape;
        }
        brep.tshapes
            .get(self.index)
            .map_or(ShapeType::Shape, |ts| ts.shape_type())
    }

    pub fn is_vertex(&self) -> bool { self.shape_type() == ShapeType::Vertex }
    pub fn is_edge(&self) -> bool { self.shape_type() == ShapeType::Edge }
    pub fn is_wire(&self) -> bool { self.shape_type() == ShapeType::Wire }
    pub fn is_face(&self) -> bool { self.shape_type() == ShapeType::Face }
    pub fn is_shell(&self) -> bool { self.shape_type() == ShapeType::Shell }
    pub fn is_solid(&self) -> bool { self.shape_type() == ShapeType::Solid }
}

impl PartialEq for Shape {
    fn eq(&self, other: &Self) -> bool {
        if self.index != usize::MAX && other.index != usize::MAX {
            self.index == other.index
                && self.location == other.location
                && self.orientation as u8 == other.orientation as u8
        } else {
            Arc::as_ptr(&self.data) == Arc::as_ptr(&other.data)
                && self.location == other.location
                && self.orientation as u8 == other.orientation as u8
        }
    }
}
impl Eq for Shape {}
impl std::hash::Hash for Shape {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        if self.index != usize::MAX {
            self.index.hash(state);
        } else {
            (Arc::as_ptr(&self.data) as u64).hash(state);
        }
    }
}

// ---- Custom Serialize: only index/orientation/location ----
// The `data` Arc is NOT serialized (it is reconstructed from the BRep on deserialization).

use serde::{Serialize, Deserialize, Serializer, Deserializer};
use serde::ser::SerializeStruct;

impl Serialize for Shape {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut s = serializer.serialize_struct("Shape", 3)?;
        s.serialize_field("index", &self.index)?;
        s.serialize_field("orientation", &self.orientation)?;
        s.serialize_field("location", &self.location)?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for Shape {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        use crate::topo::topods::tshape_flags;
        #[derive(Deserialize)]
        struct Temp {
            index: usize,
            orientation: Orientation,
            location: u32,
        }
        let t = Temp::deserialize(deserializer)?;
        Ok(Shape {
            data: Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE | tshape_flags::CLOSED | tshape_flags::CONVEX,
                point: DVec3::ZERO,
                tolerance: 0.0,
                points: Vec::new(),
            })),
            index: t.index,
            orientation: t.orientation,
            location: t.location,
        })
    }
}
