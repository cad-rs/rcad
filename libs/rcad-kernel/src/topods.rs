use std::sync::Arc;
use std::collections::HashMap;
use glam::DVec3;
use serde::{Deserialize, Serialize};
use crate::geom::{Curve2d, Curve3, Surface3};

/// Quantized vertex position for identity-based sharing.
/// Two geometrically coincident points at TOLERANCE_ABS scale produce the same key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VertexKey(u64);

impl VertexKey {
    fn from(p: DVec3) -> Self {
        const S: f64 = 1.0 / 1e-7;
        let q = |c: f64| (c * S).round() as i64;
        VertexKey(
            (q(p.x) as u64).wrapping_mul(0xbf58476d1ce4e5b9)
            ^ (q(p.y) as u64).wrapping_mul(0x94d049bb133111eb)
            ^ (q(p.z) as u64).wrapping_mul(0x9e3779b97f4a7c15)
        )
    }
}

/// OCCT TopAbs_Orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reversed,
    Internal,
    External,
}

impl Orientation {
    pub const fn is_forward(self) -> bool {
        matches!(self, Orientation::Forward)
    }
}

// ── TopoDS_TShape::myState flags (OCCT: BitLayout enum) ────────────────────
// These mirror TopoDS_TShape.hxx L69-82. Bits 0-3 are the shape type;
// bits 4-11 are boolean flags. Only the flag masks are exposed here.
pub mod tshape_flags {
    pub const FREE: u16 = 0x0010;       // TopoDS_TShape::Bit_Free
    pub const MODIFIED: u16 = 0x0020;   // TopoDS_TShape::Bit_Modified
    pub const CHECKED: u16 = 0x0040;    // TopoDS_TShape::Bit_Checked
    pub const ORIENTABLE: u16 = 0x0080; // TopoDS_TShape::Bit_Orientable
    pub const CLOSED: u16 = 0x0100;     // TopoDS_TShape::Bit_Closed
    pub const INFINITE: u16 = 0x0200;   // TopoDS_TShape::Bit_Infinite
    pub const CONVEX: u16 = 0x0400;     // TopoDS_TShape::Bit_Convex
    pub const LOCKED: u16 = 0x0800;     // TopoDS_TShape::Bit_Locked
    /// Default flags for a new TShape: Free | Modified | Orientable
    pub const DEFAULT: u16 = FREE | MODIFIED | ORIENTABLE;
}

/// OCCT TopAbs_ShapeEnum
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ShapeType {
    /// TopAbs_SHAPE — generic shape (null/unknown).
    Shape,
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
}

/// Sentinel ptr_id value for ShapeRef::synthetic(usize) — marks synthetic (non-Arc) identity.
/// High bit patterns avoid collision with real heap addresses.
const SYNTH_PTR_ID: u64 = 0xFFFFFFFF_DEAD0000;

/// TopoDS_Shape equivalent: Arc<TShape> pointer identity + Orientation.
/// This is a value type (Clone, portable) — identity is by Arc pointer, not array position.
/// `index` field is retained for O(1) access into BRep.tshapes[].
#[derive(Debug, Clone, Copy)]
pub struct ShapeRef {
    /// Arc pointer identity (Arc::as_ptr() as u64). Hash/Eq/Ord based on this.
    /// 0 = NULL; 0xFFFFFFFF_DEAD0000 | idx = synthetic (from ShapeRef::new).
    pub ptr_id: u64,
    /// Index into BRep.tshapes[] — fast accessor, NOT identity.
    pub index: usize,
    pub orientation: Orientation,
    /// TopLoc_Location index into BRep.locations[]; 0 = identity.
    pub location: u32,
}

/// Custom PartialEq: identity by ptr_id only.
impl PartialEq for ShapeRef {
    fn eq(&self, other: &Self) -> bool {
        self.ptr_id == other.ptr_id
    }
}
impl Eq for ShapeRef {}

/// Custom Hash: identity by ptr_id only.
impl std::hash::Hash for ShapeRef {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ptr_id.hash(state);
    }
}

/// Custom Ord/PartialOrd: by ptr_id.
impl PartialOrd for ShapeRef {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.ptr_id.cmp(&other.ptr_id))
    }
}
impl Ord for ShapeRef {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.ptr_id.cmp(&other.ptr_id)
    }
}

/// Custom Serialize: only index/orientation/location — ptr_id is runtime-only.
impl Serialize for ShapeRef {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeStruct;
        let mut s = serializer.serialize_struct("ShapeRef", 3)?;
        s.serialize_field("index", &self.index)?;
        s.serialize_field("orientation", &self.orientation)?;
        s.serialize_field("location", &self.location)?;
        s.end()
    }
}

/// Custom Deserialize: ptr_id defaults to 0 (caller must ensure index is valid).
impl<'de> Deserialize<'de> for ShapeRef {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        struct Temp {
            index: usize,
            orientation: Orientation,
            location: u32,
        }
        let t = Temp::deserialize(deserializer)?;
        Ok(ShapeRef { ptr_id: 0, index: t.index, orientation: t.orientation, location: t.location })
    }
}

impl ShapeRef {
    /// OCCT TopoDS_Shape::IsNull — null/uninitialized shape.
    pub const NULL: ShapeRef = ShapeRef { ptr_id: 0, index: usize::MAX, orientation: Orientation::Forward, location: 0 };

    /// Construct a synthetic ShapeRef with a computed sentinel ptr_id.
    /// Use for sentinel keys (shells/solids/compounds) and DS-adaptor code
    /// where no real TShape Arc pointer is available.
    pub const fn synthetic(index: usize) -> Self {
        Self { ptr_id: SYNTH_PTR_ID | (index as u64), index, orientation: Orientation::Forward, location: 0 }
    }
    pub const fn synthetic_with_orientation(index: usize, orientation: Orientation) -> Self {
        Self { ptr_id: SYNTH_PTR_ID | (index as u64), index, orientation, location: 0 }
    }
    pub const fn synthetic_with_location(index: usize, orientation: Orientation, location: u32) -> Self {
        Self { ptr_id: SYNTH_PTR_ID | (index as u64), index, orientation, location }
    }
    /// Construct a ShapeRef from an Arc<TShape> pointer — real identity.
    pub fn from_arc(tshape: &Arc<TShape>, orientation: Orientation, location: u32) -> Self {
        let index = 0; // caller must set index separately for BRep access
        let ptr_id = Arc::as_ptr(tshape) as u64;
        Self { ptr_id, index, orientation, location }
    }
    /// OCCT-aligned: return a copy with a different TopLoc_Location index.
    pub const fn with_location(self, location: u32) -> Self {
        Self { location, ..self }
    }
    /// OCCT TopoDS_Shape::IsNull — true if this ShapeRef is null/uninitialized.
    pub const fn is_null(&self) -> bool {
        self.index == usize::MAX
    }
    /// OCCT TopoDS_Shape::IsSame — same TShape (ignores Location and Orientation).
    pub const fn is_same(&self, other: &ShapeRef) -> bool {
        self.ptr_id == other.ptr_id
    }
    /// OCCT TopoDS_Shape::IsPartner — same TShape AND same Location.
    pub const fn is_partner(&self, other: &ShapeRef) -> bool {
        self.ptr_id == other.ptr_id && self.location == other.location
    }
    /// OCCT TopoDS_Shape::IsEqual — same TShape, Location, AND Orientation.
    pub const fn is_equal(&self, other: &ShapeRef) -> bool {
        self.ptr_id == other.ptr_id
            && self.location == other.location
            && self.orientation as u8 == other.orientation as u8
    }
    /// OCCT TopoDS_Shape::ShapeType — returns shape type from BRep.
    pub fn shape_type(&self, brep: &BRep) -> ShapeType {
        if self.is_null() { return ShapeType::Shape; }
        brep.tshapes.get(self.index).map_or(ShapeType::Shape, |ts| ts.shape_type())
    }
    /// OCCT TopoDS_Shape::NbChildren — direct sub-shape count.
    pub fn nb_children(&self, brep: &BRep) -> usize {
        brep.nb_children(*self)
    }
    /// OCCT TopoDS_Shape::EmptyCopy — create new TShape of same type, no sub-shapes.
    pub fn empty_copy(&self, brep: &mut BRep) -> ShapeRef {
        brep.empty_copy(*self)
    }
}

/// TShape — shared geometric/topological data (analogous to TopoDS_TShape + subclasses).
/// Stored in Arc<TShape> within BRep so multiple ShapeRefs share the same data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TShape {
    Vertex(TVertexData),
    Edge(TEdgeData),
    Wire(TWireData),
    Face(TFaceData),
    Shell(TShellData),
    Solid(TSolidData),
    CompSolid(Vec<ShapeRef>),
    Compound(Vec<ShapeRef>),
}

/// OCCT BRep_PointRepresentation — stores vertex parameter on a curve or surface.
/// Used for SameParameter tolerance propagation and history tracking.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "t")]
pub enum PointRepresentation {
    /// BRep_PointOnCurve — parameter on a 3D curve.
    #[serde(rename = "c")]
    PointOnCurve {
        curve: usize,
        parameter: f64,
        tolerance: f64,
    },
    /// BRep_PointOnSurface — UV parameters on a surface.
    #[serde(rename = "s")]
    PointOnSurface {
        face: usize,
        u: f64,
        v: f64,
        tolerance: f64,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TVertexData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub point: DVec3,
    #[serde(default)]
    pub tolerance: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub points: Vec<PointRepresentation>,
}

/// OCCT BRep_CurveRepresentation — how an edge lies on a face or in 3D.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CurveRepresentation {
    /// BRep_GCurve — 3D curve.
    Curve3D {
        curve: usize,
        location: u32,
    },
    /// BRep_CurveOnSurface — pcurve on a face.
    CurveOnSurface {
        face: usize,
        pcurve: Curve2d,
        range: [f64; 2],
    },
    /// BRep_CurveOnClosedSurface — two pcurves for periodic surfaces.
    CurveOnClosedSurface {
        face: usize,
        pcurve1: Curve2d,
        pcurve2: Curve2d,
        range: [f64; 2],
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TEdgeData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub curve: Option<Curve3>,
    pub first: ShapeRef,
    pub last: ShapeRef,
    pub range: [f64; 2],
    #[serde(default)]
    pub degenerated: bool,
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub pcurves: HashMap<usize, (Curve2d, f64, f64)>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub representations: Vec<CurveRepresentation>,
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub vertex_params: HashMap<usize, f64>,
    #[serde(default)]
    pub tolerance: f64,
    #[serde(default)]
    pub same_parameter: bool,
    #[serde(default)]
    pub same_range: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TWireData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub edges: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TFaceData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub surface: Option<Surface3>,
    /// TopLoc_Location index for the face surface; 0 = identity.
    /// OCCT BRep_TFace::Location — transforms surface to world coordinates.
    #[serde(default)]
    pub surface_location: u32,
    pub outer_wire: ShapeRef,
    pub inner_wires: Vec<ShapeRef>,
    pub sample_point: Option<DVec3>,
    /// UV domain [umin, umax, vmin, vmax] — used by surface area calculation.
    pub uv_domain: Option<[f64; 4]>,
    /// INTERNAL vertices (OCCT: TopAbs_INTERNAL sub-shapes, BRep_Builder.Add(aF, aV)).
    pub internal_vertices: Vec<ShapeRef>,
    /// BRep_Tool::Tolerance(aF) equivalent.
    #[serde(default)]
    pub tolerance: f64,
    /// BRep_Tool::NaturalRestriction equivalent — true when the face surface
    /// has natural boundaries (full untrimmed sphere, cylinder, cone, etc.).
    #[serde(default)]
    pub natural_restriction: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TShellData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub faces: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TSolidData {
    pub my_shapes: Vec<ShapeRef>,
    pub flags: u16,
    pub shells: Vec<ShapeRef>,
    pub internal_vertices: Vec<ShapeRef>,
    pub internal_edges: Vec<ShapeRef>,
}

/// BRep top-level shape container — all TShapes in a single pool with shared Arc ownership.
/// Analogous to OCCT's Doc/assembly structure where all TShapes live in a shared scope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub tshapes: Vec<Arc<TShape>>,
    /// 3D transformations (TopLoc_Location equivalent). Index 0 = identity.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub locations: Vec<glam::DAffine3>,
    /// OCCT-aligned: vertex identity cache — quantized position → ShapeRef.
    /// Same geometric point → same TShape::Vertex across all code paths.
    #[serde(skip)]
    pub vert_by_pos: HashMap<VertexKey, ShapeRef>,
    /// OCCT-aligned: face identity cache — wire key → ShapeRef.
    /// Two faces with the same wire structure share the same TShape::Face.
    /// Key: (outer_wire.ptr_id as usize, sorted inner_wire.ptr_ids as usize).
    #[serde(skip)]
    pub face_by_key: HashMap<(usize, Vec<usize>), ShapeRef>,
    /// OCCT-aligned: edge identity cache — (first, last, degenerated) → ShapeRef.
    #[serde(skip)]
    pub edge_by_key: HashMap<(u64, u64, bool), ShapeRef>,
}

impl Default for BRep {
    fn default() -> Self { Self::new() }
}

impl BRep {
    pub fn new() -> Self {
        Self {
            tshapes: Vec::new(),
            locations: Vec::new(),
            vert_by_pos: HashMap::new(),
            face_by_key: HashMap::new(),
            edge_by_key: HashMap::new(),
        }
    }

    /// Add a location transformation and return its index.
    /// If the transformation is identity, returns 0.
    pub fn add_location(&mut self, loc: glam::DAffine3) -> u32 {
        if loc == glam::DAffine3::IDENTITY {
            return 0;
        }
        let idx = self.locations.len();
        self.locations.push(loc);
        // Index 0 is identity — shift by 1 so 0 = identity, 1 = first real location
        (idx + 1) as u32
    }

    /// Get a location by index (0 = identity).
    pub fn get_location(&self, idx: u32) -> glam::DAffine3 {
        if idx == 0 { glam::DAffine3::IDENTITY } else {
            self.locations.get((idx - 1) as usize).copied().unwrap_or(glam::DAffine3::IDENTITY)
        }
    }

    pub fn add_tvertex(&mut self, point: DVec3) -> ShapeRef {
        // OCCT-aligned: identity-based sharing — same position → same TShape::Vertex.
        let key = VertexKey::from(point);
        if let Some(&sr) = self.vert_by_pos.get(&key) {
            return sr;
        }
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Vertex(TVertexData { my_shapes: Vec::new(), flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE | tshape_flags::CLOSED | tshape_flags::CONVEX, point, tolerance: 0.0, points: Vec::new() }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        let sr = ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 };
        self.vert_by_pos.insert(key, sr);
        sr
    }

    pub fn add_tedge(&mut self, curve: Option<Curve3>, first: ShapeRef, last: ShapeRef, range: [f64; 2]) -> ShapeRef {
        // OCCT-aligned: edge identity by (first, last) — curve carried on TShape directly.
        let ekey = (first.ptr_id, last.ptr_id, false);
        if let Some(&sr) = self.edge_by_key.get(&ekey) {
            return sr;
        }
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Edge(TEdgeData { my_shapes: vec![first, last], flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE, curve, first, last, range, degenerated: false, pcurves: HashMap::new(), representations: Vec::new(), vertex_params: HashMap::new(), tolerance: 0.0, same_parameter: true, same_range: true }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        let sr = ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 };
        self.edge_by_key.insert(ekey, sr);
        self.tshapes.push(tshape);
        sr
    }

    pub fn add_twire(&mut self, edges: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Wire(TWireData { my_shapes: edges.clone(), flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE, edges }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    pub fn add_tface(&mut self, surface: Option<Surface3>, outer_wire: ShapeRef, inner_wires: Vec<ShapeRef>, sample_point: Option<DVec3>, uv_domain: Option<[f64; 4]>, internal_vertices: Vec<ShapeRef>, natural_restriction: bool) -> ShapeRef {
        // OCCT-aligned: identity-based sharing — same wire structure → same TShape::Face.
        let mut inners: Vec<usize> = inner_wires.iter().map(|w| w.ptr_id as usize).collect();
        inners.sort_unstable();
        let key = (outer_wire.ptr_id as usize, inners);
        if let Some(&sr) = self.face_by_key.get(&key) {
            return sr;
        }
        let index = self.tshapes.len();
        // OCCT: myShapes = [outer_wire, inner_wire_1, ..., inner_wire_n].
        // Internal vertices are stored in separate list per OCCT (TopAbs_INTERNAL).
        let mut face_shapes = Vec::with_capacity(1 + inner_wires.len());
        face_shapes.push(outer_wire);
        face_shapes.extend_from_slice(&inner_wires);
        let tshape = Arc::new(TShape::Face(TFaceData { my_shapes: face_shapes, flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE, surface, surface_location: 0, outer_wire, inner_wires, sample_point, uv_domain, internal_vertices, tolerance: 0.0, natural_restriction }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        let sr = ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 };
        self.face_by_key.insert(key, sr);
        sr
    }

    pub fn add_tshell(&mut self, faces: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Shell(TShellData { my_shapes: faces.clone(), flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE, faces }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    pub fn add_tsolid(&mut self, shells: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Solid(TSolidData { my_shapes: shells.clone(), flags: tshape_flags::FREE | tshape_flags::MODIFIED, shells, internal_vertices: Vec::new(), internal_edges: Vec::new() }));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    pub fn add_tcompsolid(&mut self, solids: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::CompSolid(solids));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    pub fn add_tcompound(&mut self, shapes: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Compound(shapes));
        let ptr_id = Arc::as_ptr(&tshape) as u64;
        self.tshapes.push(tshape);
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    /// Remove all Solid TShapes from the BRep and return their indices.
    /// OCCT BuildRC removes unwanted solids from myShape after BuildResult(SOLID) added them.
    pub fn clear_solids(&mut self) -> usize {
        let before = self.tshapes.len();
        self.tshapes.retain(|ts| !matches!(&**ts, TShape::Solid(_)));
        before - self.tshapes.len()
    }

    /// Apply an affine transform to all vertex positions, edge curves, and
    /// face surfaces in-place.  Equivalent to `rcad_kernel::BRep::apply_transform`.
    pub fn apply_transform(&mut self, mat: glam::DAffine3) {
        use glam::DAffine3;
        use crate::geom::{Curve3, Surface3};

        fn xf_curve(c: &mut Curve3, mat: DAffine3) {
            match c {
                Curve3::Line(l) => {
                    l.origin = mat.transform_point3(l.origin);
                    l.direction = mat.transform_vector3(l.direction).normalize_or_zero();
                }
                Curve3::Circle(c3) => {
                    c3.center = mat.transform_point3(c3.center);
                    c3.normal = mat.transform_vector3(c3.normal).normalize_or_zero();
                }
                Curve3::Ellipse(e) => {
                    e.center = mat.transform_point3(e.center);
                    e.normal = mat.transform_vector3(e.normal).normalize_or_zero();
                    e.major_dir = mat.transform_vector3(e.major_dir).normalize_or_zero();
                }
                Curve3::BSpline(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                Curve3::Bezier(b) => {
                    for p in &mut b.control_points {
                        *p = mat.transform_point3(*p);
                    }
                }
                Curve3::Offset(o) => {
                    xf_curve(&mut o.basis, mat);
                    o.offset_dir = mat.transform_vector3(o.offset_dir).normalize_or_zero();
                }
                Curve3::Hyperbola(h) => {
                    h.center = mat.transform_point3(h.center);
                    h.normal = mat.transform_vector3(h.normal).normalize_or_zero();
                    h.major_dir = mat.transform_vector3(h.major_dir).normalize_or_zero();
                }
                Curve3::Parabola(p) => {
                    p.vertex = mat.transform_point3(p.vertex);
                    p.normal = mat.transform_vector3(p.normal).normalize_or_zero();
                    p.axis_dir = mat.transform_vector3(p.axis_dir).normalize_or_zero();
                }
                Curve3::CircularHelix(h) => {
                    h.origin = mat.transform_point3(h.origin);
                    h.axis = mat.transform_vector3(h.axis).normalize_or_zero();
                    h.ref_dir = mat.transform_vector3(h.ref_dir).normalize_or_zero();
                }
                Curve3::SineWave(s) => {
                    s.origin = mat.transform_point3(s.origin);
                    s.baseline_dir = mat.transform_vector3(s.baseline_dir).normalize_or_zero();
                    s.amplitude_dir = mat.transform_vector3(s.amplitude_dir).normalize_or_zero();
                }
            }
        }

        fn xf_surface(s: &mut Surface3, mat: DAffine3) {
            match s {
                Surface3::Plane(p) => {
                    p.origin = mat.transform_point3(p.origin);
                    p.normal = mat.transform_vector3(p.normal).normalize_or_zero();
                }
                Surface3::Cylinder(c) => {
                    c.origin = mat.transform_point3(c.origin);
                    c.axis = mat.transform_vector3(c.axis).normalize_or_zero();
                    c.ref_dir = mat.transform_vector3(c.ref_dir).normalize_or_zero();
                }
                Surface3::Sphere(s) => {
                    s.center = mat.transform_point3(s.center);
                    s.axis = mat.transform_vector3(s.axis).normalize_or_zero();
                    s.ref_dir = mat.transform_vector3(s.ref_dir).normalize_or_zero();
                }
                Surface3::Cone(c) => {
                    c.apex = mat.transform_point3(c.apex);
                    c.axis = mat.transform_vector3(c.axis).normalize_or_zero();
                }
                Surface3::Torus(t) => {
                    t.center = mat.transform_point3(t.center);
                    t.axis = mat.transform_vector3(t.axis).normalize_or_zero();
                }
                Surface3::Ellipsoid(e) => {
                    e.center = mat.transform_point3(e.center);
                    e.axis = mat.transform_vector3(e.axis).normalize_or_zero();
                    e.ref_dir = mat.transform_vector3(e.ref_dir).normalize_or_zero();
                }
                Surface3::Helicoid(h) => {
                    h.origin = mat.transform_point3(h.origin);
                    h.axis = mat.transform_vector3(h.axis).normalize_or_zero();
                    h.ref_dir = mat.transform_vector3(h.ref_dir).normalize_or_zero();
                }
                Surface3::Pipe(p) => {
                    xf_curve(&mut p.spine, mat);
                    p.ref_dir = mat.transform_vector3(p.ref_dir).normalize_or_zero();
                }
                Surface3::BSpline(b) => {
                    for row in &mut b.control_points {
                        for p in row {
                            *p = mat.transform_point3(*p);
                        }
                    }
                }
                Surface3::Bezier(b) => {
                    for row in &mut b.control_points {
                        for p in row {
                            *p = mat.transform_point3(*p);
                        }
                    }
                }
                Surface3::Coons(c) => {
                    xf_curve(&mut c.south, mat);
                    xf_curve(&mut c.north, mat);
                    xf_curve(&mut c.west, mat);
                    xf_curve(&mut c.east, mat);
                }
                Surface3::Offset(o) => {
                    xf_surface(&mut o.basis, mat);
                }
                Surface3::Revolution(r) => {
                    xf_curve(&mut r.profile, mat);
                    r.axis_origin = mat.transform_point3(r.axis_origin);
                    r.axis_dir = mat.transform_vector3(r.axis_dir).normalize_or_zero();
                }
                Surface3::LinearExtrusion(e) => {
                    xf_curve(&mut e.profile, mat);
                    e.direction = mat.transform_vector3(e.direction).normalize_or_zero();
                }
                Surface3::Ruled(r) => {
                    xf_curve(&mut r.start, mat);
                    xf_curve(&mut r.end, mat);
                }
                Surface3::Trimmed(t) => {
                    xf_surface(&mut t.basis, mat);
                }
                Surface3::TriBezier(_) => {}
            }
        }

        for ts in &mut self.tshapes {
            match &mut *Arc::make_mut(ts) {
                TShape::Vertex(vd) => {
                    vd.point = mat.transform_point3(vd.point);
                }
                TShape::Edge(ed) => {
                    if let Some(ref mut curve) = ed.curve {
                        xf_curve(curve, mat);
                    }
                }
                TShape::Face(fd) => {
                    if let Some(ref mut surface) = fd.surface {
                        xf_surface(surface, mat);
                    }
                }
                _ => {}
            }
        }
    }

    /// OCCT TopoDS_TShape::EmptyCopy — create a new TShape of the same type
    /// with no sub-shapes. Preserves flags. Returns the new TShape index.
    pub fn empty_copy(&mut self, r: ShapeRef) -> ShapeRef {
        let ts = &self.tshapes[r.index];
        let new = match &**ts {
            TShape::Vertex(vd) => Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: vd.flags,
                point: vd.point,
                tolerance: vd.tolerance,
                points: Vec::new(),
            })),
            TShape::Edge(ed) => Arc::new(TShape::Edge(TEdgeData {
                my_shapes: Vec::new(),
                flags: ed.flags,
                curve: ed.curve.clone(),
                first: ShapeRef::NULL,
                last: ShapeRef::NULL,
                range: ed.range,
                degenerated: ed.degenerated,
                pcurves: HashMap::new(),
                representations: Vec::new(),
                vertex_params: HashMap::new(),
                tolerance: ed.tolerance,
                same_parameter: ed.same_parameter,
                same_range: ed.same_range,
            })),
            TShape::Wire(wd) => Arc::new(TShape::Wire(TWireData {
                my_shapes: Vec::new(),
                flags: wd.flags,
                edges: Vec::new(),
            })),
            TShape::Face(fd) => Arc::new(TShape::Face(TFaceData {
                my_shapes: Vec::new(),
                flags: fd.flags,
                surface: fd.surface.clone(),
                surface_location: fd.surface_location,
                outer_wire: ShapeRef::NULL,
                inner_wires: Vec::new(),
                sample_point: fd.sample_point,
                uv_domain: fd.uv_domain,
                internal_vertices: Vec::new(),
                tolerance: fd.tolerance,
                natural_restriction: fd.natural_restriction,
            })),
            TShape::Shell(sd) => Arc::new(TShape::Shell(TShellData {
                my_shapes: Vec::new(),
                flags: sd.flags,
                faces: Vec::new(),
            })),
            TShape::Solid(sd) => Arc::new(TShape::Solid(TSolidData {
                my_shapes: Vec::new(),
                flags: sd.flags,
                shells: Vec::new(),
                internal_vertices: Vec::new(),
                internal_edges: Vec::new(),
            })),
            TShape::CompSolid(_) => Arc::new(TShape::CompSolid(Vec::new())),
            TShape::Compound(_) => Arc::new(TShape::Compound(Vec::new())),
        };
        let index = self.tshapes.len();
        self.tshapes.push(new);
        let ptr_id = Arc::as_ptr(&self.tshapes[index]) as u64;
        ShapeRef { ptr_id, index, orientation: Orientation::Forward, location: 0 }
    }

    /// Count face TShapes in this BRep (for shifted-key pcurve lookup).
    pub fn nb_faces(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Face(_))).count()
    }

    /// OCCT TopoDS_TShape::NbChildren — number of direct sub-shapes.
    pub fn nb_children(&self, r: ShapeRef) -> usize {
        match &*self.tshapes[r.index] {
            TShape::Vertex(_) => 0,
            TShape::Edge(_) => 2,
            TShape::Wire(wd) => wd.edges.len(),
            TShape::Face(fd) => 1 + fd.inner_wires.len() + fd.internal_vertices.len(),
            TShape::Shell(sd) => sd.faces.len(),
            TShape::Solid(sd) => sd.shells.len() + sd.internal_vertices.len() + sd.internal_edges.len(),
            TShape::CompSolid(cd) => cd.len(),
            TShape::Compound(cd) => cd.len(),
        }
    }

    /// OCCT-aligned: check if a shape has a given flag set.
    pub fn has_flag(&self, r: ShapeRef, flag: u16) -> bool {
        let flags = match &*self.tshapes[r.index] {
            TShape::Vertex(vd) => vd.flags,
            TShape::Edge(ed) => ed.flags,
            TShape::Wire(wd) => wd.flags,
            TShape::Face(fd) => fd.flags,
            TShape::Shell(sd) => sd.flags,
            TShape::Solid(sd) => sd.flags,
            TShape::CompSolid(_) => tshape_flags::FREE,
            TShape::Compound(_) => tshape_flags::FREE,
        };
        (flags & flag) != 0
    }

    /// OCCT-aligned: set a single flag on a shape.
    pub fn set_flag(&mut self, r: ShapeRef, flag: u16, on: bool) {
        let arc = Arc::get_mut(&mut self.tshapes[r.index]).expect("set_flag: Arc still shared");
        let flags = match arc {
            TShape::Vertex(vd) => &mut vd.flags,
            TShape::Edge(ed) => &mut ed.flags,
            TShape::Wire(wd) => &mut wd.flags,
            TShape::Face(fd) => &mut fd.flags,
            TShape::Shell(sd) => &mut sd.flags,
            TShape::Solid(sd) => &mut sd.flags,
            _ => return,
        };
        if on { *flags |= flag; } else { *flags &= !flag; }
    }

    pub fn vertex(&self, r: ShapeRef) -> &TVertexData {
        match &*self.tshapes[r.index] {
            TShape::Vertex(v) => v,
            _ => panic!("ShapeRef {} is not a Vertex", r.index),
        }
    }

    pub fn edge(&self, r: ShapeRef) -> &TEdgeData {
        match &*self.tshapes[r.index] {
            TShape::Edge(e) => e,
            _ => panic!("ShapeRef {} is not an Edge", r.index),
        }
    }

    pub fn wire(&self, r: ShapeRef) -> &TWireData {
        match &*self.tshapes[r.index] {
            TShape::Wire(w) => w,
            _ => panic!("ShapeRef {} is not a Wire", r.index),
        }
    }

    pub fn face(&self, r: ShapeRef) -> &TFaceData {
        match &*self.tshapes[r.index] {
            TShape::Face(f) => f,
            _ => panic!("ShapeRef {} is not a Face", r.index),
        }
    }

    pub fn shell(&self, r: ShapeRef) -> &TShellData {
        match &*self.tshapes[r.index] {
            TShape::Shell(s) => s,
            _ => panic!("ShapeRef {} is not a Shell", r.index),
        }
    }

    pub fn solid(&self, r: ShapeRef) -> &TSolidData {
        match &*self.tshapes[r.index] {
            TShape::Solid(s) => s,
            _ => panic!("ShapeRef {} is not a Solid", r.index),
        }
    }

    /// Build a minimal Cube (for testing).
    pub fn build_unit_cube() -> (Self, ShapeRef) {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        bld.build_unit_cube(&mut brep);
        let root = bld.root.expect("cube should have root solid");
        (brep, root)
    }

    // -----------------------------------------------------------------------
    // Mutable accessors (OCCT BRep_Builder pattern — mutate after creation)
    // -----------------------------------------------------------------------

    /// Mutate a vertex's data (panics if vertex index is out of range or Arc is shared).
    pub fn vertex_mut(&mut self, r: ShapeRef) -> &mut TVertexData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("vertex_mut: Arc still shared") {
            TShape::Vertex(v) => v,
            _ => panic!("vertex_mut: ShapeRef {} is not a Vertex", r.index),
        }
    }

    /// Mutate an edge's data.
    pub fn edge_mut(&mut self, r: ShapeRef) -> &mut TEdgeData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("edge_mut: Arc still shared") {
            TShape::Edge(e) => e,
            _ => panic!("edge_mut: ShapeRef {} is not an Edge", r.index),
        }
    }

    /// Mutate a wire's data.
    pub fn wire_mut(&mut self, r: ShapeRef) -> &mut TWireData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("wire_mut: Arc still shared") {
            TShape::Wire(w) => w,
            _ => panic!("wire_mut: ShapeRef {} is not a Wire", r.index),
        }
    }

    /// Mutate a face's data.
    pub fn face_mut(&mut self, r: ShapeRef) -> &mut TFaceData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("face_mut: Arc still shared") {
            TShape::Face(f) => f,
            _ => panic!("face_mut: ShapeRef {} is not a Face", r.index),
        }
    }

    /// Mutate a shell's data.
    pub fn shell_mut(&mut self, r: ShapeRef) -> &mut TShellData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("shell_mut: Arc still shared") {
            TShape::Shell(s) => s,
            _ => panic!("shell_mut: ShapeRef {} is not a Shell", r.index),
        }
    }

    /// Mutate a solid's data.
    pub fn solid_mut(&mut self, r: ShapeRef) -> &mut TSolidData {
        match Arc::get_mut(&mut self.tshapes[r.index]).expect("solid_mut: Arc still shared") {
            TShape::Solid(s) => s,
            _ => panic!("solid_mut: ShapeRef {} is not a Solid", r.index),
        }
    }

    /// Count vertex TShapes (includes orphan vertices not referenced by any edge).
    pub fn nb_vertices(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Vertex(_))).count()
    }

    /// Count edge TShapes.
    pub fn nb_edges(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Edge(_))).count()
    }

    /// Count shell TShapes.
    pub fn nb_shells(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Shell(_))).count()
    }

    /// Count solid TShapes.
    pub fn nb_solids(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(ts.as_ref(), TShape::Solid(_))).count()
    }

    /// Axis-aligned bounding box computed from all vertex positions.
    /// Returns `None` if the BRep has no vertices.
    pub fn bounding_box(&self) -> Option<[DVec3; 2]> {
        let mut mn = DVec3::splat(f64::INFINITY);
        let mut mx = DVec3::splat(f64::NEG_INFINITY);
        for ts in &self.tshapes {
            if let TShape::Vertex(v) = ts.as_ref() {
                mn = mn.min(v.point);
                mx = mx.max(v.point);
            }
        }
        if mn.x.is_infinite() { None } else { Some([mn, mx]) }
    }

    /// Centroid of all vertex positions (simple average).
    pub fn center(&self) -> DVec3 {
        let mut sum = DVec3::ZERO;
        let mut count = 0usize;
        for ts in &self.tshapes {
            if let TShape::Vertex(v) = ts.as_ref() {
                sum += v.point;
                count += 1;
            }
        }
        if count == 0 { DVec3::ZERO } else { sum / count as f64 }
    }

    /// Build a compound from multiple topods::BRep shapes.
    pub fn compound_from_shapes(shapes: &[BRep]) -> BRep {
        let mut t = BRep::new();
        let mut refs = Vec::new();
        for s in shapes {
            let mut solid_refs = Vec::new();
            for (ti, ts) in s.tshapes.iter().enumerate() {
                if matches!(ts.as_ref(), TShape::Solid(_)) {
                    solid_refs.push(t.add_tsolid(Vec::new()));
                }
            }
            refs.extend(solid_refs);
        }
        if !refs.is_empty() {
            t.add_tcompound(refs);
        }
        t
    }
}

// ---------------------------------------------------------------------------
// BRepTool trait — OCCT BRep_Tool free-function equivalents
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool equivalent: parameter/tolerance/pcurve queries on a BRep.
///
/// In OCCT these are free functions (`BRep_Tool::Parameter(aV, aE, aF)` etc.).
/// Here they are methods on a `BRepTool` trait so the boolean pipeline can
/// be generic over the data source (real BRep or DS adaptor).
pub trait BRepTool {
    /// BRep_Tool::Pnt(aV) — 3D position of a vertex.
    fn vertex_position(&self, v: ShapeRef) -> DVec3;
    /// BRep_Tool::Tolerance(aV) — 3D tolerance of a vertex.
    fn vertex_tolerance(&self, v: ShapeRef) -> f64;
    /// BRep_Tool::Degenerated(aE).
    fn is_edge_degenerated(&self, e: ShapeRef) -> bool;
    /// TopExp: given one vertex of an edge, return the other vertex.
    fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef) -> ShapeRef;
    /// First vertex of an edge (FORWARD orientation, OCCT TopExp::FirstVertex).
    fn first_vertex(&self, edge: ShapeRef) -> ShapeRef;
    /// Last vertex of an edge (FORWARD orientation, OCCT TopExp::LastVertex).
    fn last_vertex(&self, edge: ShapeRef) -> ShapeRef;
    /// TopExp::FirstVertex on an oriented edge — canonical first for FORWARD,
    /// canonical last for REVERSED.  Matches OCCT's orientation-aware topology.
    fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation) -> ShapeRef;
    /// BRep_Tool::Parameter(aV, aE, aF) — vertex parameter on edge's pcurve.
    fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, face: ShapeRef) -> Option<f64>;
    /// BRep_Tool::CurveOnSurface(aE, aF) — pcurve of edge on face.
    fn curve_on_surface(&self, edge: ShapeRef, face: ShapeRef) -> Option<&(Curve2d, f64, f64)>;
    /// BRep_Tool::Surface(aF) — face surface (local coordinates, no Location applied).
    fn face_surface(&self, face: ShapeRef) -> Option<&Surface3>;
    /// BRep_Tool::Surface(aF) with Location applied — returns world-coordinate surface.
    fn face_surface_world(&self, face: ShapeRef) -> Option<Surface3>;
    /// BRep_Tool::Curve(aE) with Location applied — returns 3D curve and range in world coordinates.
    fn edge_curve_world(&self, edge: ShapeRef) -> Option<(Curve3, [f64; 2])>;
    /// UResolution: parameter tolerance in U direction (OCCT: BRepAdaptor_Surface::UResolution).
    fn u_resolution(&self, face: ShapeRef, tol3d: f64) -> f64;
    /// VResolution: parameter tolerance in V direction.
    fn v_resolution(&self, face: ShapeRef, tol3d: f64) -> f64;
    /// OCCT L204-207: vertex orientation (TopAbs_INTERNAL for split-edge interior vertices).
    /// Default: Forward (non-INTERNAL). Override when INTERNAL vertex data is available.
    fn vertex_orientation(&self, _v: ShapeRef) -> Orientation { Orientation::Forward }
    /// BRep_Tool::IsClosed(aE, aF) — true when the edge appears twice on the face
    /// (periodic surface seam).  Checks for CurveOnClosedSurface representation.
    fn is_edge_closed_on_face(&self, edge: ShapeRef, face: ShapeRef) -> bool {
        self.curve_on_surface(edge, face).is_some()
            && self.curve_on_surface_second(edge, face).is_some()
    }
    /// Retrieve the second pcurve for periodic surfaces (CurveOnClosedSurface).
    /// Returns None for non-periodic faces or if only one pcurve exists.
    fn curve_on_surface_second(&self, edge: ShapeRef, face: ShapeRef) -> Option<&(Curve2d, f64, f64)> {
        let _ = (edge, face);
        None
    }

    // ── OCCT BRep_Tool / TopoDS_Shape convenience queries ──

    /// OCCT TopoDS_Shape::Closed — checks CLOSED flag on the TShape.
    /// For Shell, this is a simple flag check; use `is_shell_closed` for
    /// the full edge-count verification (BRepCheck_Shell).
    fn is_closed(&self, s: ShapeRef) -> bool {
        self.has_flag(s, tshape_flags::CLOSED)
    }

    /// BRep_Tool::SameParameter(aE) — true when edge pcurves match 3D curve.
    fn edge_same_parameter(&self, e: ShapeRef) -> bool {
        self.edge_data(e).map(|ed| ed.same_parameter).unwrap_or(true)
    }

    /// BRep_Tool::SameRange(aE) — true when edge pcurve ranges match 3D range.
    fn edge_same_range(&self, e: ShapeRef) -> bool {
        self.edge_data(e).map(|ed| ed.same_range).unwrap_or(true)
    }

    /// BRep_Tool::NaturalRestriction(aF) — true when face surface bounds are
    /// determined by the underlying surface's natural domain.
    fn face_natural_restriction(&self, f: ShapeRef) -> bool {
        self.face_data(f).map(|fd| fd.natural_restriction).unwrap_or(true)
    }

    /// BRep_Tool::Curve(aE) — raw 3D curve reference (no Location applied).
    /// Returns the curve data directly (OCCT-aligned: geometry on TShape).
    fn edge_curve_data(&self, e: ShapeRef) -> Option<Curve3> {
        self.edge_data(e).and_then(|ed| ed.curve.clone())
    }

    /// BRep_Tool::Range(aE) — 3D curve parameter range.
    fn edge_range(&self, e: ShapeRef) -> [f64; 2] {
        self.edge_data(e).map(|ed| ed.range).unwrap_or([0.0, 0.0])
    }

    /// BRep_Tool::Tolerance(s) — geometric tolerance for any shape type.
    fn tolerance(&self, s: ShapeRef) -> f64;

    // ── Extension helpers (not in OCCT's BRep_Tool) ──

    /// TopoDS_Shape::ShapeType — returns the shape type.
    fn shape_type(&self, s: ShapeRef) -> ShapeType;

    /// Check CLOSED flag directly (bypasses trait default).
    fn has_flag(&self, s: ShapeRef, flag: u16) -> bool;

    /// Access Edge data (for implementing default methods).
    fn edge_data(&self, e: ShapeRef) -> Option<&TEdgeData>;

    /// Access Face data (for implementing default methods).
    fn face_data(&self, f: ShapeRef) -> Option<&TFaceData>;
}

impl BRepTool for BRep {
    fn vertex_position(&self, v: ShapeRef) -> DVec3 {
        let pt = self.vertex(v).point;
        self.get_location(v.location).transform_point3(pt)
    }

    fn vertex_tolerance(&self, v: ShapeRef) -> f64 {
        self.vertex(v).tolerance
    }

    fn is_edge_degenerated(&self, e: ShapeRef) -> bool {
        self.edge(e).degenerated
    }

    fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef) -> ShapeRef {
        let ed = self.edge(edge);
        if ed.first.index == v.index { ed.last } else { ed.first }
    }

    fn first_vertex(&self, edge: ShapeRef) -> ShapeRef {
        self.edge(edge).first
    }

    fn last_vertex(&self, edge: ShapeRef) -> ShapeRef {
        self.edge(edge).last
    }

    fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation) -> ShapeRef {
        if orientation == Orientation::Reversed {
            self.last_vertex(edge)
        } else {
            self.first_vertex(edge)
        }
    }

    fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, _face: ShapeRef) -> Option<f64> {
        self.edge(edge).vertex_params.get(&vertex.index).copied()
    }

    fn curve_on_surface(&self, edge: ShapeRef, face: ShapeRef) -> Option<&(Curve2d, f64, f64)> {
        self.edge(edge).pcurves.get(&face.index)
    }

    fn is_edge_closed_on_face(&self, edge: ShapeRef, face: ShapeRef) -> bool {
        let ed = self.edge(edge);
        ed.representations.iter().any(|r| matches!(r, CurveRepresentation::CurveOnClosedSurface { face: f, .. } if *f == face.index))
            || ed.pcurves.contains_key(&face.index)
                && ed.pcurves.contains_key(&(face.index + self.nb_faces()))
    }

    fn curve_on_surface_second(&self, edge: ShapeRef, face: ShapeRef) -> Option<&(Curve2d, f64, f64)> {
        let shifted = face.index + self.nb_faces();
        self.edge(edge).pcurves.get(&shifted)
    }

    fn face_surface(&self, face: ShapeRef) -> Option<&Surface3> {
        self.face(face).surface.as_ref()
    }

    fn face_surface_world(&self, face: ShapeRef) -> Option<Surface3> {
        let fd = self.face(face);
        let surface = fd.surface.as_ref()?.clone();
        let loc = self.get_location(face.location);
        if loc == glam::DAffine3::IDENTITY {
            Some(surface)
        } else {
            Some(crate::geom::transform_surface(&surface, &loc))
        }
    }

    fn edge_curve_world(&self, edge: ShapeRef) -> Option<(Curve3, [f64; 2])> {
        let ed = self.edge(edge);
        let crv = ed.curve.as_ref()?.clone();
        let loc = self.get_location(edge.location);
        if loc == glam::DAffine3::IDENTITY {
            Some((crv, ed.range))
        } else {
            Some((crate::geom::transform_curve(&crv, &loc), ed.range))
        }
    }

    fn u_resolution(&self, face: ShapeRef, tol3d: f64) -> f64 {
        match self.face(face).surface.as_ref() {
            Some(surf) => u_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }

    fn v_resolution(&self, face: ShapeRef, tol3d: f64) -> f64 {
        match self.face(face).surface.as_ref() {
            Some(surf) => v_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }

    fn tolerance(&self, s: ShapeRef) -> f64 {
        match &*self.tshapes[s.index] {
            TShape::Vertex(vd) => vd.tolerance,
            TShape::Edge(ed) => ed.tolerance,
            TShape::Face(fd) => fd.tolerance,
            _ => 0.0,
        }
    }

    fn shape_type(&self, s: ShapeRef) -> ShapeType {
        s.shape_type(self)
    }

    fn has_flag(&self, s: ShapeRef, flag: u16) -> bool {
        self.has_flag(s, flag)
    }

    fn edge_data(&self, e: ShapeRef) -> Option<&TEdgeData> {
        match &*self.tshapes[e.index] {
            TShape::Edge(ed) => Some(ed),
            _ => None,
        }
    }

    fn face_data(&self, f: ShapeRef) -> Option<&TFaceData> {
        match &*self.tshapes[f.index] {
            TShape::Face(fd) => Some(fd),
            _ => None,
        }
    }
}

/// Surface-aware UResolution (OCCT: BRepAdaptor_Surface::UResolution).
pub fn u_resolution_for_surface(surf: &Surface3, tol3d: f64) -> f64 {
    match surf {
        Surface3::Sphere(s) => tol3d / s.radius.max(1e-15),
        Surface3::Cylinder(c) => tol3d / c.radius.max(1e-15),
        Surface3::Cone(_) => tol3d * 1e-3,
        Surface3::Torus(t) => tol3d / t.major_radius.max(1e-15),
        _ => tol3d,
    }
}

/// Surface-aware VResolution (OCCT: BRepAdaptor_Surface::VResolution).
pub fn v_resolution_for_surface(surf: &Surface3, tol3d: f64) -> f64 {
    match surf {
        Surface3::Sphere(s) => tol3d / s.radius.max(1e-15),
        Surface3::Cylinder(_) => tol3d,
        Surface3::Cone(_) => tol3d,
        Surface3::Torus(t) => tol3d / t.minor_radius.max(1e-15),
        _ => tol3d,
    }
}

// ---------------------------------------------------------------------------
// BRepBuilder — OCCT BRep_Builder equivalent for incrementally constructing BRep
// ---------------------------------------------------------------------------

// ── Backward-compat flat-index access methods ──────────────────────
// These provide the old brep.vertices/edges/solids/geom patterns
// for modules that haven't been migrated yet.
impl BRep {
    /// Collect vertex points in tshape order (like old brep.vertices).
    pub fn flat_vertices(&self) -> Vec<DVec3> {
        self.tshapes.iter().filter_map(|ts| {
            if let TShape::Vertex(vd) = &**ts { Some(vd.point) } else { None }
        }).collect()
    }

    /// Collect edge endpoint pairs in tshape order (like old brep.edges).
    pub fn flat_edges(&self) -> Vec<(usize, usize)> {
        self.tshapes.iter().filter_map(|ts| {
            if let TShape::Edge(ed) = &**ts { Some((ed.first.index, ed.last.index)) } else { None }
        }).collect()
    }

    /// Count vertex TShapes.
    pub fn vertex_count(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(&**ts, TShape::Vertex(_))).count()
    }

    /// Count edge TShapes.
    pub fn edge_count(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(&**ts, TShape::Edge(_))).count()
    }

    /// Count face TShapes.
    pub fn face_count(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(&**ts, TShape::Face(_))).count()
    }

    /// Count solid TShapes.
    pub fn solid_count(&self) -> usize {
        self.tshapes.iter().filter(|ts| matches!(&**ts, TShape::Solid(_))).count()
    }

    /// Check if any solid TShapes exist.
    pub fn has_solids(&self) -> bool {
        self.tshapes.iter().any(|ts| matches!(&**ts, TShape::Solid(_)))
    }

    /// Get vertex point by tshape index.
    pub fn vertex_point(&self, idx: usize) -> Option<DVec3> {
        self.tshapes.get(idx).and_then(|ts| {
            if let TShape::Vertex(vd) = &**ts { Some(vd.point) } else { None }
        })
    }

    /// Add a new edge from flat indices (creates ShapeRef for each vertex).
    /// Returns the tshape index of the new edge.
    pub fn add_edge_flat(&mut self, start_idx: usize, end_idx: usize, curve: Option<Curve3>, range: [f64; 2]) -> usize {
        let first = self.tshapes.get(start_idx).map(|ts| ShapeRef {
            ptr_id: Arc::as_ptr(ts) as u64,
            index: start_idx,
            orientation: Orientation::Forward,
            location: 0,
        }).unwrap_or(ShapeRef::NULL);
        let last = self.tshapes.get(end_idx).map(|ts| ShapeRef {
            ptr_id: Arc::as_ptr(ts) as u64,
            index: end_idx,
            orientation: Orientation::Forward,
            location: 0,
        }).unwrap_or(ShapeRef::NULL);
        let sr = self.add_tedge(curve, first, last, range);
        sr.index
    }
}

pub struct BRepBuilder {
    pub root: Option<ShapeRef>,
    vertex_cache: Vec<[f64; 3]>,
}

impl BRepBuilder {
    pub fn new() -> Self {
        Self { root: None, vertex_cache: Vec::new() }
    }

    fn find_or_add_vertex(&mut self, brep: &mut BRep, pt: DVec3) -> ShapeRef {
        for (i, &cached) in self.vertex_cache.iter().enumerate() {
            let dp = DVec3::new(cached[0], cached[1], cached[2]) - pt;
            if dp.length_squared() < 1e-30 {
                return ShapeRef::synthetic(i);
            }
        }
        let r = brep.add_tvertex(pt);
        self.vertex_cache.push([pt.x, pt.y, pt.z]);
        r
    }

    pub fn build_unit_cube(&mut self, brep: &mut BRep) {
        // 8 vertices
        let v: Vec<ShapeRef> = [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
        ].iter().map(|&p| self.find_or_add_vertex(brep, DVec3::new(p[0], p[1], p[2]))).collect();

        // 12 edges (unit box: 4 around each of 3 axis-aligned rings)
        let faces_verts: [[usize; 4]; 6] = [
            [0, 1, 2, 3], // Z- (bottom)
            [4, 5, 6, 7], // Z+ (top)
            [0, 1, 5, 4], // Y- (front)
            [2, 3, 7, 6], // Y+ (back)
            [0, 3, 7, 4], // X- (left)
            [1, 2, 6, 5], // X+ (right)
        ];

        // Build 12 edges, reuse shared vertex references between faces
        let mut edge_map = std::collections::HashMap::new();
        let mut edge_for_face = Vec::new();

        for &fv in &faces_verts {
            let mut face_edges = Vec::new();
            for j in 0..4 {
                let a = fv[j];
                let b = fv[(j + 1) % 4];
                let key = if a < b { (a, b, true) } else { (b, a, false) };
                let e_idx = *edge_map.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    let va = v[a];
                    let vb = v[b];
                    brep.add_tedge(None, va, ShapeRef::synthetic_with_orientation(vb.index, Orientation::Reversed), [0.0, 1.0]).index
                });
                let orient = if a < b { Orientation::Forward } else { Orientation::Reversed };
                face_edges.push(ShapeRef::synthetic_with_orientation(e_idx, orient));
            }
            edge_for_face.push(face_edges);
        }

        // Build 6 faces, collecting their refs for shell building
        let mut face_refs = Vec::new();
        for (i, face_edges) in edge_for_face.into_iter().enumerate() {
            let wire = brep.add_twire(face_edges);
            let face = brep.add_tface(None, wire, vec![], None, None, vec![], true);
            let face_ref = face;

            // Outer normal: orient face based on the face index
            // For a unit cube at origin, faces 0,2,4 have inward normals, orient REVERSED
            // faces 1,3,5 have outward normals, orient FORWARD
            let orient = match i {
                0 => Orientation::Reversed,  // Z-
                1 => Orientation::Forward,   // Z+
                2 => Orientation::Reversed,  // Y-
                3 => Orientation::Forward,   // Y+
                4 => Orientation::Reversed,  // X-
                5 => Orientation::Forward,   // X+
                _ => unreachable!(),
            };
            face_refs.push(ShapeRef::synthetic_with_orientation(face.index, orient));
        }

        // Build shell from the 6 oriented faces
        let shell = brep.add_tshell(face_refs);
        // Build solid from the shell
        let solid = brep.add_tsolid(vec![shell]);
        self.root = Some(solid);
    }

    // -----------------------------------------------------------------------
    // Extended BRepBuilder API — OCCT BRep_Builder equivalent
    // -----------------------------------------------------------------------

    /// Add a vertex with tolerance.
    pub fn add_vertex(&mut self, brep: &mut BRep, pt: DVec3, tol: f64) -> ShapeRef {
        let r = brep.add_tvertex(pt);
        brep.vertex_mut(r).tolerance = tol;
        r
    }

    /// Update vertex tolerance (max with existing).
    pub fn update_vertex_tolerance(&mut self, brep: &mut BRep, v: ShapeRef, tol: f64) {
        let vd = brep.vertex_mut(v);
        vd.tolerance = vd.tolerance.max(tol);
    }

    /// Add an edge with curve, vertices, and range.
    pub fn add_edge(&mut self, brep: &mut BRep,
        curve: Option<Curve3>, v1: ShapeRef, v2: ShapeRef, range: [f64; 2]) -> ShapeRef {
        brep.add_tedge(curve, v1, v2, range)
    }

    /// Add a pcurve to an edge for a specific face.
    pub fn add_pcurve(&mut self, brep: &mut BRep,
        edge: ShapeRef, face: ShapeRef, pc: Curve2d, t1: f64, t2: f64) {
        brep.edge_mut(edge).pcurves.insert(face.index, (pc, t1, t2));
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, theTol) — update edge tolerance.
    pub fn update_edge_tolerance(&mut self, brep: &mut BRep, edge: ShapeRef, tol: f64) {
        let ed = brep.edge_mut(edge);
        ed.tolerance = ed.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, aC2d, aF, theTol) — set pcurve on face.
    pub fn update_edge_pcurve(&mut self, brep: &mut BRep,
        edge: ShapeRef, pcurve: Curve2d, face: ShapeRef, tol: f64)
    {
        let ed = brep.edge_mut(edge);
        let (ta, tb) = pc_parameter_range(&pcurve);
        ed.pcurves.insert(face.index, (pcurve.clone(), ta, tb));
        ed.representations.push(CurveRepresentation::CurveOnSurface {
            face: face.index,
            pcurve,
            range: [ta, tb],
        });
        ed.tolerance = ed.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, aC3d) — set 3D curve.
    pub fn update_edge_curve3d(&mut self, brep: &mut BRep,
        edge: ShapeRef, curve: usize, location: u32) {
        let ed = brep.edge_mut(edge);
        ed.representations.push(CurveRepresentation::Curve3D { curve, location });
    }

    /// Set vertex parameter on an edge's pcurve.
    pub fn set_vertex_param(&mut self, brep: &mut BRep,
        edge: ShapeRef, vertex: ShapeRef, param: f64) {
        brep.edge_mut(edge).vertex_params.insert(vertex.index, param);
    }

    /// Set degenerated flag on an edge.
    pub fn set_edge_degenerated(&mut self, brep: &mut BRep, edge: ShapeRef, flag: bool) {
        brep.edge_mut(edge).degenerated = flag;
    }

    /// OCCT BRep_Builder::Add(aW).Closed(aW) — mark wire as closed.
    pub fn close_wire(&mut self, brep: &mut BRep, wire: ShapeRef) {
        brep.wire_mut(wire).flags |= tshape_flags::CLOSED;
    }

    /// OCCT BRep_Builder::Add(aShell).Closed(aShell) — mark shell as closed via Closed flag.
    pub fn close_shell(&mut self, brep: &mut BRep, shell: ShapeRef) {
        brep.shell_mut(shell).flags |= tshape_flags::CLOSED;
    }

    /// Make a wire (empty container).
    pub fn make_wire(&mut self, brep: &mut BRep) -> ShapeRef {
        brep.add_twire(vec![])
    }

    /// Add an edge to a wire.
    pub fn add_to_wire(&mut self, brep: &mut BRep, wire: ShapeRef, edge: ShapeRef) {
        let wd = brep.wire_mut(wire);
        wd.edges.push(edge);
        wd.my_shapes.push(edge);
    }

    /// Build a wire from edges.
    pub fn build_wire(&mut self, brep: &mut BRep, edges: Vec<ShapeRef>) -> ShapeRef {
        brep.add_twire(edges)
    }

    /// Make a face from a surface and outer wire.
    pub fn make_face(&mut self, brep: &mut BRep,
        surface: Option<Surface3>, outer_wire: ShapeRef) -> ShapeRef {
        brep.add_tface(surface, outer_wire, vec![], None, None, vec![], true)
    }

    /// Add an inner wire to a face.
    pub fn add_to_face(&mut self, brep: &mut BRep, face: ShapeRef, inner_wire: ShapeRef) {
        let fd = brep.face_mut(face);
        fd.inner_wires.push(inner_wire);
        fd.my_shapes.push(inner_wire);
    }

    /// Add an internal vertex to a face.
    pub fn add_internal_vertex(&mut self, brep: &mut BRep, face: ShapeRef, v: ShapeRef) {
        let fd = brep.face_mut(face);
        fd.internal_vertices.push(v);
        fd.my_shapes.push(v);
    }

    /// Add an edge with section-curve semantics (MakeSectEdge equivalent).
    /// Creates an edge with pcurves for both faces.
    pub fn add_section_edge(&mut self, brep: &mut BRep,
        curve: Option<Curve3>, v1: ShapeRef, v2: ShapeRef, range: [f64; 2],
        pc_a: Option<&Curve2d>, face_a: Option<ShapeRef>,
        pc_b: Option<&Curve2d>, face_b: Option<ShapeRef>,
    ) -> ShapeRef {
        let e = brep.add_tedge(curve, v1, v2, range);
        if let (Some(pc), Some(fa)) = (pc_a, face_a) {
            let (t1, t2) = pc_parameter_range(pc);
            brep.edge_mut(e).pcurves.insert(fa.index, (pc.clone(), t1, t2));
        }
        if let (Some(pc), Some(fb)) = (pc_b, face_b) {
            let (t1, t2) = pc_parameter_range(pc);
            brep.edge_mut(e).pcurves.insert(fb.index, (pc.clone(), t1, t2));
        }
        e
    }

    /// Make a shell (empty container).
    pub fn make_shell(&mut self, brep: &mut BRep) -> ShapeRef {
        brep.add_tshell(vec![])
    }

    /// Add a face to a shell.
    pub fn add_to_shell(&mut self, brep: &mut BRep, shell: ShapeRef, face: ShapeRef) {
        let sd = brep.shell_mut(shell);
        sd.faces.push(face);
        sd.my_shapes.push(face);
    }

    /// Make a solid from shells.
    pub fn make_solid(&mut self, brep: &mut BRep, shells: Vec<ShapeRef>) -> ShapeRef {
        brep.add_tsolid(shells)
    }

    /// Make a compsolid from solids.
    pub fn make_compsolid(&mut self, brep: &mut BRep, solids: Vec<ShapeRef>) -> ShapeRef {
        brep.add_tcompsolid(solids)
    }

    /// Make a compound from shapes.
    pub fn make_compound(&mut self, brep: &mut BRep, shapes: Vec<ShapeRef>) -> ShapeRef {
        brep.add_tcompound(shapes)
    }

    /// Add a shape to an existing compound.
    pub fn add_to_compound(&self, brep: &mut BRep, compound: ShapeRef, shape: ShapeRef) {
        let ts = Arc::get_mut(&mut brep.tshapes[compound.index])
            .expect("add_to_compound: unique ownership required");
        match ts {
            TShape::Compound(shapes) => shapes.push(shape),
            _ => panic!("add_to_compound: shape is not a Compound"),
        }
    }
}

/// Get the parameter range for a Curve2d (Trimmed → stored range, Circle → [0, 2π]).
fn pc_parameter_range(curve: &Curve2d) -> (f64, f64) {
    match curve {
        Curve2d::Trimmed(tc) => (tc.t_min, tc.t_max),
        Curve2d::Circle(_) => (0.0, std::f64::consts::TAU),
        _ => (0.0, 1.0),
    }
}

// ---------------------------------------------------------------------------
// ShapeType helpers
// ---------------------------------------------------------------------------

impl TShape {
    pub fn shape_type(&self) -> ShapeType {
        match self {
            TShape::Vertex(_) => ShapeType::Vertex,
            TShape::Edge(_) => ShapeType::Edge,
            TShape::Wire(_) => ShapeType::Wire,
            TShape::Face(_) => ShapeType::Face,
            TShape::Shell(_) => ShapeType::Shell,
            TShape::Solid(_) => ShapeType::Solid,
            TShape::CompSolid(_) => ShapeType::CompSolid,
            TShape::Compound(_) => ShapeType::Compound,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_values() {
        assert!(Orientation::Forward.is_forward());
        assert!(!Orientation::Reversed.is_forward());
        assert!(!Orientation::Internal.is_forward());
        assert!(!Orientation::External.is_forward());
    }

    #[test]
    fn test_shape_ref_construction() {
        let r = ShapeRef::synthetic(5);
        assert_eq!(r.index, 5);
        assert_eq!(r.orientation, Orientation::Forward);

        let r2 = ShapeRef::synthetic_with_orientation(3, Orientation::Reversed);
        assert_eq!(r2.index, 3);
        assert_eq!(r2.orientation, Orientation::Reversed);
    }

    #[test]
    fn test_empty_brep() {
        let brep = BRep::new();
        assert!(brep.tshapes.is_empty());
    }

    #[test]
    fn test_vertex_creation() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(brep.tshapes.len(), 1);
        assert_eq!(v.index, 0);
        assert_eq!(brep.vertex(v).point, DVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_vertex_sharing() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::ZERO);
        // Two separate calls create two distinct TShapes (OCCT: each call adds a unique TShape)
        assert_eq!(brep.tshapes.len(), 2);
        assert_eq!(brep.vertex(v0).point, brep.vertex(v1).point);
        // Same Arc — points are equal but TShape identity differs (OCCT behavior)
        assert!(!Arc::ptr_eq(&brep.tshapes[v0.index], &brep.tshapes[v1.index]));
    }

    #[test]
    fn test_edge_creation() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
        assert_eq!(brep.tshapes.len(), 3);
        let ed = brep.edge(e);
        assert_eq!(ed.first.index, v0.index);
        assert_eq!(ed.last.index, v1.index);
    }

    #[test]
    fn test_edge_shares_vertex_tshape() {
        // Two edges sharing the same vertex TShape
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let v2 = brep.add_tvertex(DVec3::new(1.0, 1.0, 0.0));

        let e0 = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v1, v2, [0.0, 1.0]);

        // Both edges reference v1 at index 1
        assert_eq!(brep.edge(e0).last.index, v1.index);
        assert_eq!(brep.edge(e1).first.index, v1.index);
        // Same TShape identity (v1)
        assert!(Arc::ptr_eq(&brep.tshapes[brep.edge(e0).last.index], &brep.tshapes[brep.edge(e1).first.index]));
    }

    #[test]
    fn test_wire_and_face() {
        let mut brep = BRep::new();
        let v = (0..4).map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0))).collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0], v[1], [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v[1], v[2], [0.0, 1.0]);
        let e2 = brep.add_tedge(None, v[2], v[3], [0.0, 1.0]);
        let wire = brep.add_twire(vec![e0, e1, e2]);
        let face = brep.add_tface(None, wire, vec![], None, None, vec![], true);
        let fd = brep.face(face);
        assert_eq!(brep.tshapes.len(), 9); // 4V + 3E + 1W + 1F
        assert!(fd.inner_wires.is_empty());
        assert!(fd.sample_point.is_none());
    }

    #[test]
    fn test_shape_type_discrimination() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        assert_eq!(brep.tshapes[v.index].shape_type(), ShapeType::Vertex);

        let e = brep.add_tedge(None, v, v, [0.0, 1.0]);
        assert_eq!(brep.tshapes[e.index].shape_type(), ShapeType::Edge);

        let w = brep.add_twire(vec![e]);
        assert_eq!(brep.tshapes[w.index].shape_type(), ShapeType::Wire);

        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        assert_eq!(brep.tshapes[f.index].shape_type(), ShapeType::Face);

        let sh = brep.add_tshell(vec![f]);
        assert_eq!(brep.tshapes[sh.index].shape_type(), ShapeType::Shell);

        let so = brep.add_tsolid(vec![sh]);
        assert_eq!(brep.tshapes[so.index].shape_type(), ShapeType::Solid);
    }

    #[test]
    fn test_orientation_on_shape_reference() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        // Edge references v1 with REVERSED orientation (last vertex traversed backward)
        let e = brep.add_tedge(None, v0, ShapeRef::synthetic_with_orientation(v1.index, Orientation::Reversed), [0.0, 1.0]);
        let ed = brep.edge(e);
        assert_eq!(ed.last.orientation, Orientation::Reversed);
    }

    #[test]
    fn test_clone_preserves_shape_identity() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let _v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(None, v0, v0, [0.0, 1.0]);

        let cloned = brep.clone();
        // Same TShape identity in clone (Arc::ptr_eq across clone)
        assert_eq!(cloned.tshapes.len(), brep.tshapes.len());
        assert!(Arc::ptr_eq(&cloned.tshapes[e.index], &brep.tshapes[e.index]));
    }

    #[test]
    fn test_serde_roundtrip() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        let e = brep.add_tedge(None, v, v, [0.0, 1.0]);

        let json = serde_json::to_string(&brep).unwrap();
        let restored: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.tshapes.len(), 2);
        assert_eq!(restored.vertex(ShapeRef::synthetic(0)).point, DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(restored.edge(ShapeRef::synthetic(1)).range, [0.0, 1.0]);
    }

    #[test]
    fn test_brep_builder_unit_cube() {
        let (_brep, root) = BRep::build_unit_cube();
        assert_eq!(root.orientation, Orientation::Forward);
    }

    #[test]
    fn test_wire_edge_orientation() {
        let mut brep = BRep::new();
        let v = (0..4).map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0))).collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0], v[1], [0.0, 1.0]);

        // Wire with forward edge
        let wire_fwd = brep.add_twire(vec![e0]);
        assert_eq!(brep.tshapes[wire_fwd.index].shape_type(), ShapeType::Wire);

        // Wire with reversed edge
        let e0_rev = ShapeRef::synthetic_with_orientation(e0.index, Orientation::Reversed);
        let wire_rev = brep.add_twire(vec![e0_rev]);
        if let TShape::Wire(ref wd) = *brep.tshapes[wire_rev.index] {
            assert_eq!(wd.edges[0].orientation, Orientation::Reversed);
        } else {
            panic!("expected Wire");
        }

        // Same TShape for the edge, different orientation on the reference
        assert!(Arc::ptr_eq(&brep.tshapes[e0.index], &brep.tshapes[e0_rev.index]));
    }

    #[test]
    fn test_to_topods_roundtrip_sa() {
        // Build a simple box BRep, convert to topods and back, check SA
        let mut orig = crate::BRep::new();
        // 8 vertices of a unit cube
        let v: Vec<_> = [
            [0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0],
            [0.0,0.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0],[0.0,1.0,1.0],
        ].iter().map(|&p| orig.vertices.push(crate::topology::Vertex { point: DVec3::new(p[0],p[1],p[2]) })).collect::<Vec<_>>();
        // Simple box with 6 faces, 12 edges
        let face_vert_idxs: [[usize;4]; 6] = [
            [0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[0,3,7,4],[1,2,6,5],
        ];
        for &fv in &face_vert_idxs {
            for j in 0..4 {
                let a = fv[j]; let b = fv[(j+1)%4];
                // Check if edge already exists
                let exists = orig.edges.iter().position(|e| (e.start==a&&e.end==b)||(e.start==b&&e.end==a));
                if exists.is_none() {
                    orig.edges.push(crate::topology::Edge { start: a, end: b });
                }
            }
        }
        // Convert to topods and back
        let t = orig.to_topods();
        let back = crate::BRep::from_topods(&t);
        assert_eq!(back.vertices.len(), 8, "vertex count");
        assert_eq!(back.edges.len(), 12, "edge count");
        // Check surface area (unit cube SA = 6.0)
        let sa = crate::surface_area(&back);
        assert!((sa - 6.0).abs() < 0.01, "unit cube SA: expected 6.0, got {}", sa);
    }

    #[test]
    fn test_face_natural_restriction_default_true() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        let w = brep.add_twire(vec![]);
        // default (no explicit nr) → true via add_tface with natural_restriction=true
        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        let fd = brep.face(f);
        assert!(fd.natural_restriction);
    }

    #[test]
    fn test_face_natural_restriction_false() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], false);
        let fd = brep.face(f);
        assert!(!fd.natural_restriction);
    }

    #[test]
    fn test_face_natural_restriction_serialize_roundtrip() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], false);
        let json = serde_json::to_string(&brep).unwrap();
        let restored: BRep = serde_json::from_str(&json).unwrap();
        let fd = restored.face(f);
        assert!(!fd.natural_restriction);
    }
}
