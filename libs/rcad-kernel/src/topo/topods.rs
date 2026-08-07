use crate::geom::{Curve2d, Curve3, Surface3, SurfaceEval};
use crate::core::precision::{CONFUSION, parametric_default};
use crate::math::bspl::{
    bezier_curve_resolution, bezier_surface_resolution, bspline_curve_resolution,
    bspline_surface_resolution,
};
use std::f64::consts::PI;
pub use crate::topo::topo_shape::Shape;
use glam::DVec3;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;

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
                ^ (q(p.z) as u64).wrapping_mul(0x9e3779b97f4a7c15),
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

    /// OCCT TopAbs::Compose (TopAbs.hxx L69-78) — the cumulated orientation of
    /// a sub-shape with orientation `other` inside a parent with orientation
    /// `self`. External dominates, then Internal, then Forward/Reversed xor.
    pub const fn compose(self, other: Orientation) -> Orientation {
        match self {
            Orientation::Forward => other,
            Orientation::Reversed => match other {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
                o => o,
            },
            Orientation::Internal => Orientation::Internal,
            Orientation::External => Orientation::External,
        }
    }
}

/// OCCT TopAbs_State — classification result (TopAbs_State.hxx L27-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum State {
    In,
    Out,
    On,
    Unknown,
}

// 鈹€鈹€ TopoDS_TShape::myState flags (OCCT: BitLayout enum) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// These mirror TopoDS_TShape.hxx L69-82. Bits 0-3 are the shape type;
// bits 4-11 are boolean flags. Only the flag masks are exposed here.
pub mod tshape_flags {
    pub const FREE: u16 = 0x0010; // TopoDS_TShape::Bit_Free
    pub const MODIFIED: u16 = 0x0020; // TopoDS_TShape::Bit_Modified
    pub const CHECKED: u16 = 0x0040; // TopoDS_TShape::Bit_Checked
    pub const ORIENTABLE: u16 = 0x0080; // TopoDS_TShape::Bit_Orientable
    pub const CLOSED: u16 = 0x0100; // TopoDS_TShape::Bit_Closed
    pub const INFINITE: u16 = 0x0200; // TopoDS_TShape::Bit_Infinite
    pub const CONVEX: u16 = 0x0400; // TopoDS_TShape::Bit_Convex
    pub const LOCKED: u16 = 0x0800; // TopoDS_TShape::Bit_Locked
    /// Default flags for a new TShape: Free | Modified | Orientable
    pub const DEFAULT: u16 = FREE | MODIFIED | ORIENTABLE;
}

/// OCCT TopAbs_ShapeEnum
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ShapeType {
    /// TopAbs_SHAPE �?generic shape (null/unknown).
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

/// Sentinel ptr_id value for Shape::synthetic(usize) �?marks synthetic (non-Arc) identity.
/// High bit patterns avoid collision with real heap addresses.
const SYNTH_PTR_ID: u64 = 0xFFFFFFFF_DEAD0000;



/// TShape �?shared geometric/topological data (analogous to TopoDS_TShape + subclasses).
/// Stored in Arc<TShape> within BRep so multiple ShapeRefs share the same data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TShape {
    Vertex(TVertexData),
    Edge(TEdgeData),
    Wire(TWireData),
    Face(TFaceData),
    Shell(TShellData),
    Solid(TSolidData),
    CompSolid(Vec<Shape>),
    Compound(Vec<Shape>),
}

/// OCCT BRep_PointRepresentation �?stores vertex parameter on a curve or surface.
/// Used for SameParameter tolerance propagation and history tracking.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "t")]
pub enum PointRepresentation {
    /// BRep_PointOnCurve �?parameter on a 3D curve.
    #[serde(rename = "c")]
    PointOnCurve {
        curve: usize,
        parameter: f64,
        tolerance: f64,
    },
    /// BRep_PointOnSurface �?UV parameters on a surface.
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
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub point: DVec3,
    #[serde(default)]
    pub tolerance: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub points: Vec<PointRepresentation>,
}

/// OCCT BRep_CurveRepresentation �?how an edge lies on a face or in 3D.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CurveRepresentation {
    /// BRep_GCurve �?3D curve.
    Curve3D { curve: usize, location: u32 },
    /// BRep_CurveOnSurface �?pcurve on a face.
    CurveOnSurface {
        face: usize,
        pcurve: Curve2d,
        range: [f64; 2],
    },
    /// BRep_CurveOnClosedSurface �?two pcurves for periodic surfaces.
    CurveOnClosedSurface {
        face: usize,
        pcurve1: Curve2d,
        pcurve2: Curve2d,
        range: [f64; 2],
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TEdgeData {
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub curve: Option<Curve3>,
    pub first: Shape,
    pub last: Shape,
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
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub edges: Vec<Shape>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TFaceData {
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub surface: Option<Surface3>,
    /// TopLoc_Location index for the face surface; 0 = identity.
    /// OCCT BRep_TFace::Location �?transforms surface to world coordinates.
    #[serde(default)]
    pub surface_location: u32,
    pub outer_wire: Shape,
    pub inner_wires: Vec<Shape>,
    pub sample_point: Option<DVec3>,
    /// UV domain [umin, umax, vmin, vmax] �?used by surface area calculation.
    pub uv_domain: Option<[f64; 4]>,
    /// INTERNAL vertices (OCCT: TopAbs_INTERNAL sub-shapes, BRep_Builder.Add(aF, aV)).
    pub internal_vertices: Vec<Shape>,
    /// BRep_Tool::Tolerance(aF) equivalent.
    #[serde(default)]
    pub tolerance: f64,
    /// BRep_Tool::NaturalRestriction equivalent �?true when the face surface
    /// has natural boundaries (full untrimmed sphere, cylinder, cone, etc.).
    #[serde(default)]
    pub natural_restriction: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TShellData {
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub faces: Vec<Shape>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TSolidData {
    pub my_shapes: Vec<Shape>,
    pub flags: u16,
    pub shells: Vec<Shape>,
    pub internal_vertices: Vec<Shape>,
    pub internal_edges: Vec<Shape>,
}

/// BRep top-level shape container �?all TShapes in a single pool with shared Arc ownership.
/// Analogous to OCCT's Doc/assembly structure where all TShapes live in a shared scope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub tshapes: Vec<Arc<TShape>>,
    /// 3D transformations (TopLoc_Location equivalent). Index 0 = identity.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub locations: Vec<glam::DAffine3>,
    /// OCCT-aligned: vertex identity cache �?quantized position �?Shape.
    /// Same geometric point �?same TShape::Vertex across all code paths.
    #[serde(skip)]
    pub vert_by_pos: HashMap<VertexKey, Shape>,
    /// OCCT-aligned: face identity cache �?wire key �?Shape.
    /// Two faces with the same wire structure share the same TShape::Face.
    /// Key: (outer_wire.ptr_id() as usize, sorted inner_wire.ptr_ids as usize).
    #[serde(skip)]
    pub face_by_key: HashMap<(usize, Vec<usize>), Shape>,
    /// OCCT-aligned: edge identity cache �?(first, last, degenerated) �?Shape.
    #[serde(skip)]
    pub edge_by_key: HashMap<(u64, u64, bool), Shape>,
}

impl Default for BRep {
    fn default() -> Self {
        Self::new()
    }
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
        // Index 0 is identity �?shift by 1 so 0 = identity, 1 = first real location
        (idx + 1) as u32
    }

    /// Get a location by index (0 = identity).
    pub fn get_location(&self, idx: u32) -> glam::DAffine3 {
        if idx == 0 {
            glam::DAffine3::IDENTITY
        } else {
            self.locations
                .get((idx - 1) as usize)
                .copied()
                .unwrap_or(glam::DAffine3::IDENTITY)
        }
    }

    pub fn add_tvertex(&mut self, point: DVec3) -> Shape {
        // OCCT-aligned: identity-based sharing �?same position �?same TShape::Vertex.
        let key = VertexKey::from(point);
        if let Some(sr) = self.vert_by_pos.get(&key) {
            return sr.clone();
        }
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point,
            tolerance: CONFUSION,
            points: Vec::new(),
        }));

        self.tshapes.push(tshape);
        let sr = Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        };
        self.vert_by_pos.insert(key, sr.clone());
        sr
    }

    pub fn add_tedge(
        &mut self,
        curve: Option<Curve3>,
        first: Shape,
        last: Shape,
        range: [f64; 2],
    ) -> Shape {
        // OCCT-aligned: edge identity by (first, last) �?curve carried on TShape directly.
        let ekey = (first.ptr_id(), last.ptr_id(), false);
        if let Some(sr) = self.edge_by_key.get(&ekey) {
            return sr.clone();
        }
        let index = self.tshapes.len();
        // OCCT: BRep_Tool::Parameter stores vertex→param mapping on edge creation.
        // Compute by matching vertex positions to curve range endpoints.
        let vertex_params = {
            let mut vp = HashMap::new();
            if let Some(ref c) = curve {
                use crate::geom::CurveEval;
                let p0 = c.point_at(range[0]);
                let p1 = c.point_at(range[1]);
                if let Some(vd) = first.as_vertex() {
                    let d0 = (p0 - vd.point).length_squared();
                    let d1 = (p1 - vd.point).length_squared();
                    vp.insert(first.index, if d0 <= d1 { range[0] } else { range[1] });
                }
                if let Some(vd) = last.as_vertex() {
                    let d0 = (p0 - vd.point).length_squared();
                    let d1 = (p1 - vd.point).length_squared();
                    vp.insert(last.index, if d0 <= d1 { range[0] } else { range[1] });
                }
            }
            vp
        };
        // OCCT BRep_Tool::Degenerated: an edge whose two vertices coincide and
        // which carries no geometric curve (e.g. the sphere apex edges) is
        // degenerated. The DS then skips it in FillShrunkData (no pave blocks).
        let is_degenerated = curve.is_none() && first.ptr_id() == last.ptr_id();
        let tshape = Arc::new(TShape::Edge(TEdgeData {
            my_shapes: vec![first.clone(), last.clone()],
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            curve,
            first,
            last,
            range,
            degenerated: is_degenerated,
            pcurves: HashMap::new(),
            representations: Vec::new(),
            vertex_params,
            tolerance: CONFUSION,
            same_parameter: true,
            same_range: true,
        }));

        self.tshapes.push(tshape);
        let sr = Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        };
        self.edge_by_key.insert(ekey, sr.clone());
        sr
    }

    pub fn add_twire(&mut self, edges: Vec<Shape>) -> Shape {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Wire(TWireData {
            my_shapes: edges.clone(),
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            edges,
        }));

        self.tshapes.push(tshape);
        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    pub fn add_tface(
        &mut self,
        surface: Option<Surface3>,
        outer_wire: Shape,
        inner_wires: Vec<Shape>,
        sample_point: Option<DVec3>,
        uv_domain: Option<[f64; 4]>,
        internal_vertices: Vec<Shape>,
        natural_restriction: bool,
    ) -> Shape {
        // OCCT-aligned: identity-based sharing �?same wire structure �?same TShape::Face.
        let mut inners: Vec<usize> = inner_wires.iter().map(|w| w.ptr_id() as usize).collect();
        inners.sort_unstable();
        let key = (outer_wire.ptr_id() as usize, inners);
        if let Some(sr) = self.face_by_key.get(&key) {
            return sr.clone();
        }
        let index = self.tshapes.len();
        // OCCT: myShapes = [outer_wire, inner_wire_1, ..., inner_wire_n].
        // Internal vertices are stored in separate list per OCCT (TopAbs_INTERNAL).
        let mut face_shapes = Vec::with_capacity(1 + inner_wires.len());
        face_shapes.push(outer_wire.clone());
        face_shapes.extend_from_slice(&inner_wires);
        let tshape = Arc::new(TShape::Face(TFaceData {
            my_shapes: face_shapes,
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            surface,
            surface_location: 0,
            outer_wire,
            inner_wires,
            sample_point,
            uv_domain,
            internal_vertices,
            // OCCT BRep_Tool::Tolerance(face) — a standard BRep face carries the
            // Precision::Confusion() tolerance (1e-7), not zero.
            tolerance: CONFUSION,
            natural_restriction,
        }));

        self.tshapes.push(tshape);
        let sr = Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        };
        self.face_by_key.insert(key, sr.clone());
        sr
    }

    /// Ensure a Vertex TShape exists at the given flat index.
    /// Grows tshapes vec if needed (fills with dummy vertices).
    /// Returns Shape with the correct flat index.
    pub fn ensure_vertex_at(&mut self, idx: usize, point: DVec3) -> Shape {
        if idx < self.tshapes.len() {
            if let TShape::Vertex(_) = &*self.tshapes[idx] {

                return Shape {
                    data: self.tshapes[idx].clone(),
                    index: idx,
                    orientation: Orientation::Forward,
                    location: 0,
                };
            }
        }
        let tshape = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point,
            tolerance: CONFUSION,
            points: Vec::new(),
        }));

        let dummy = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point: DVec3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        }));
        while self.tshapes.len() <= idx {
            self.tshapes.push(dummy.clone());
        }
        self.tshapes[idx] = tshape;
        Shape {
            data: self.tshapes[idx].clone(),
            index: idx,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// Ensure an Edge TShape exists at the given flat index.
    /// Grows tshapes vec if needed. The vertex references (first/last) use
    /// their own indices (caller must ensure they are set at correct flat indices).
    pub fn ensure_edge_at(
        &mut self,
        idx: usize,
        curve: Option<Curve3>,
        first: Shape,
        last: Shape,
        range: [f64; 2],
    ) -> Shape {
        if idx < self.tshapes.len() {
            if let TShape::Edge(_) = &*self.tshapes[idx] {

                return Shape {
                    data: self.tshapes[idx].clone(),
                    index: idx,
                    orientation: Orientation::Forward,
                    location: 0,
                };
            }
        }
        let tshape = Arc::new(TShape::Edge(TEdgeData {
            my_shapes: vec![first.clone(), last.clone()],
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            curve,
            first,
            last,
            range,
            degenerated: false,
            pcurves: HashMap::new(),
            representations: Vec::new(),
            vertex_params: HashMap::new(),
            tolerance: 0.0,
            same_parameter: true,
            same_range: true,
        }));

        let dummy = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point: DVec3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        }));
        while self.tshapes.len() <= idx {
            self.tshapes.push(dummy.clone());
        }
        self.tshapes[idx] = tshape;
        Shape {
            data: self.tshapes[idx].clone(),
            index: idx,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// Ensure a Wire TShape exists at the given flat index.
    /// Grows tshapes vec if needed. The edges list contains ShapeRefs to edge TShapes.
    pub fn ensure_wire_at(&mut self, idx: usize, edges: Vec<Shape>) -> Shape {
        if idx < self.tshapes.len() {
            if let TShape::Wire(_) = &*self.tshapes[idx] {

                return Shape {
                    data: self.tshapes[idx].clone(),
                    index: idx,
                    orientation: Orientation::Forward,
                    location: 0,
                };
            }
        }
        let tshape = Arc::new(TShape::Wire(TWireData {
            my_shapes: edges.clone(),
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            edges,
        }));

        let dummy = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point: DVec3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        }));
        while self.tshapes.len() <= idx {
            self.tshapes.push(dummy.clone());
        }
        self.tshapes[idx] = tshape;
        Shape {
            data: self.tshapes[idx].clone(),
            index: idx,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// Ensure a Face TShape exists at the given flat index.
    /// Grows tshapes vec if needed.
    pub fn ensure_face_at(
        &mut self,
        idx: usize,
        surface: Option<Surface3>,
        outer_wire: Shape,
        inner_wires: Vec<Shape>,
        sample_point: Option<DVec3>,
        uv_domain: Option<[f64; 4]>,
        internal_vertices: Vec<Shape>,
        natural_restriction: bool,
    ) -> Shape {
        if idx < self.tshapes.len() {
            if let TShape::Face(_) = &*self.tshapes[idx] {

                return Shape {
                    data: self.tshapes[idx].clone(),
                    index: idx,
                    orientation: Orientation::Forward,
                    location: 0,
                };
            }
        }
        let mut face_shapes = Vec::with_capacity(1 + inner_wires.len());
        face_shapes.push(outer_wire.clone());
        face_shapes.extend_from_slice(&inner_wires);
        let tshape = Arc::new(TShape::Face(TFaceData {
            my_shapes: face_shapes,
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            surface,
            surface_location: 0,
            outer_wire,
            inner_wires,
            sample_point,
            uv_domain,
            internal_vertices,
            // OCCT BRep_Tool::Tolerance(face) — a standard BRep face carries the
            // Precision::Confusion() tolerance (1e-7), not zero.
            tolerance: CONFUSION,
            natural_restriction,
        }));

        let dummy = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::FREE
                | tshape_flags::MODIFIED
                | tshape_flags::ORIENTABLE
                | tshape_flags::CLOSED
                | tshape_flags::CONVEX,
            point: DVec3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        }));
        while self.tshapes.len() <= idx {
            self.tshapes.push(dummy.clone());
        }
        self.tshapes[idx] = tshape;

        Shape {
            data: self.tshapes[idx].clone(),
            index: idx,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    pub fn add_tshell(&mut self, faces: Vec<Shape>) -> Shape {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Shell(TShellData {
            my_shapes: faces.clone(),
            flags: tshape_flags::FREE | tshape_flags::MODIFIED | tshape_flags::ORIENTABLE,
            faces,
        }));

        self.tshapes.push(tshape);
        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    pub fn add_tsolid(&mut self, shells: Vec<Shape>) -> Shape {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Solid(TSolidData {
            my_shapes: shells.clone(),
            flags: tshape_flags::FREE | tshape_flags::MODIFIED,
            shells,
            internal_vertices: Vec::new(),
            internal_edges: Vec::new(),
        }));

        self.tshapes.push(tshape);
        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    pub fn add_tcompsolid(&mut self, solids: Vec<Shape>) -> Shape {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::CompSolid(solids));

        self.tshapes.push(tshape);
        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    pub fn add_tcompound(&mut self, shapes: Vec<Shape>) -> Shape {
        let index = self.tshapes.len();
        let tshape = Arc::new(TShape::Compound(shapes));

        self.tshapes.push(tshape);
        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// Remove all Solid TShapes from the BRep and return their indices.
    /// OCCT BuildRC removes unwanted solids from myShape after BuildResult(SOLID) added them.
    pub fn clear_solids(&mut self) -> usize {
        let before = self.tshapes.len();
        self.tshapes
            .retain(|ts| !matches!(ts.as_ref(), &TShape::Solid(_)));
        before - self.tshapes.len()
    }

    /// Apply an affine transform to all vertex positions, edge curves, and
    /// face surfaces in-place.  Equivalent to `rcad_kernel::BRep::apply_transform`.
    pub fn apply_transform(&mut self, mat: glam::DAffine3) {
        use crate::geom::{Curve3, Surface3};
        use glam::DAffine3;

        fn xf_curve(c: &mut Curve3, mat: DAffine3) {
            match c {
                Curve3::Line(l) => {
                    l.origin = mat.transform_point3(l.origin);
                    l.direction = mat.transform_vector3(l.direction).normalize_or_zero();
                }
                Curve3::Circle(c3) => {
                    c3.center = mat.transform_point3(c3.center);
                    c3.normal = mat.transform_vector3(c3.normal).normalize_or_zero();
                    // OCCT gp_Circ::Transform updates all frame axes.
                    c3.x_dir = mat.transform_vector3(c3.x_dir).normalize_or_zero();
                    c3.y_dir = c3.normal.cross(c3.x_dir).normalize();
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
                Curve3::Trimmed(tc) => xf_curve(&mut tc.curve, mat),
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

    /// OCCT TopoDS_TShape::EmptyCopy �?create a new TShape of the same type
    /// with no sub-shapes. Preserves flags. Returns the new TShape index.
    pub fn empty_copy(&mut self, r: Shape) -> Shape {
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
                first: Shape::null(),
                last: Shape::null(),
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
                outer_wire: Shape::null(),
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

        Shape {
            data: self.tshapes[index].clone(),
            index,
            orientation: Orientation::Forward,
            location: 0,
        }
    }

    /// Count face TShapes in this BRep (for shifted-key pcurve lookup).

    /// Create a Shape from a flat tshape index (Forward orientation, no location).
    pub fn shape_at(&self, idx: usize) -> Shape {
        Shape::from_parts(self.tshapes[idx].clone(), idx, 0, Orientation::Forward)
    }

    /// Create a Shape from a flat tshape index with orientation and location.
    pub fn shape_at_full(&self, idx: usize, orientation: Orientation, location: u32) -> Shape {
        Shape::from_parts(self.tshapes[idx].clone(), idx, location, orientation)
    }

    /// Count face TShapes in this BRep (for shifted-key pcurve lookup).
pub fn nb_faces(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
            .count()
    }

    /// OCCT TopoDS_TShape::NbChildren �?number of direct sub-shapes.
    pub fn nb_children(&self, r: Shape) -> usize {
        match &*self.tshapes[r.index] {
            TShape::Vertex(_) => 0,
            TShape::Edge(_) => 2,
            TShape::Wire(wd) => wd.edges.len(),
            TShape::Face(fd) => 1 + fd.inner_wires.len() + fd.internal_vertices.len(),
            TShape::Shell(sd) => sd.faces.len(),
            TShape::Solid(sd) => {
                sd.shells.len() + sd.internal_vertices.len() + sd.internal_edges.len()
            }
            TShape::CompSolid(cd) => cd.len(),
            TShape::Compound(cd) => cd.len(),
        }
    }

    /// OCCT-aligned: check if a shape has a given flag set.
    pub fn has_flag(&self, r: Shape, flag: u16) -> bool {
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
    pub fn set_flag(&mut self, r: Shape, flag: u16, on: bool) {
        let arc = Arc::make_mut(&mut self.tshapes[r.index]);
        let flags = match arc {
            TShape::Vertex(vd) => &mut vd.flags,
            TShape::Edge(ed) => &mut ed.flags,
            TShape::Wire(wd) => &mut wd.flags,
            TShape::Face(fd) => &mut fd.flags,
            TShape::Shell(sd) => &mut sd.flags,
            TShape::Solid(sd) => &mut sd.flags,
            _ => return,
        };
        if on {
            *flags |= flag;
        } else {
            *flags &= !flag;
        }
    }

    pub fn vertex(&self, r: Shape) -> &TVertexData {
        match &*self.tshapes[r.index] {
            TShape::Vertex(v) => v,
            _ => panic!("Shape {} is not a Vertex", r.index),
        }
    }

    pub fn edge(&self, r: Shape) -> &TEdgeData {
        match &*self.tshapes[r.index] {
            TShape::Edge(e) => e,
            _ => panic!("Shape {} is not an Edge", r.index),
        }
    }

    pub fn wire(&self, r: Shape) -> &TWireData {
        match &*self.tshapes[r.index] {
            TShape::Wire(w) => w,
            _ => panic!("Shape {} is not a Wire", r.index),
        }
    }

    pub fn face(&self, r: Shape) -> &TFaceData {
        match &*self.tshapes[r.index] {
            TShape::Face(f) => f,
            _ => panic!("Shape {} is not a Face", r.index),
        }
    }

    pub fn shell(&self, r: Shape) -> &TShellData {
        match &*self.tshapes[r.index] {
            TShape::Shell(s) => s,
            _ => panic!("Shape {} is not a Shell", r.index),
        }
    }

    pub fn solid(&self, r: Shape) -> &TSolidData {
        match &*self.tshapes[r.index] {
            TShape::Solid(s) => s,
            _ => panic!("Shape {} is not a Solid", r.index),
        }
    }

    /// Build a minimal Cube (for testing).
    pub fn build_unit_cube() -> (Self, Shape) {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        bld.build_unit_cube(&mut brep);
        let root = bld.root.expect("cube should have root solid");
        (brep, root)
    }

    // -----------------------------------------------------------------------
    // Mutable accessors (OCCT BRep_Builder pattern �?mutate after creation)
    // -----------------------------------------------------------------------

    /// Mutate a vertex's data (panics if vertex index is out of range or Arc is shared).
    pub fn vertex_mut(&mut self, r: Shape) -> &mut TVertexData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Vertex(v) => v,
            _ => panic!("vertex_mut: Shape {} is not a Vertex", r.index),
        }
    }

    /// Mutate an edge's data.
    pub fn edge_mut(&mut self, r: Shape) -> &mut TEdgeData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Edge(e) => e,
            _ => panic!("edge_mut: Shape {} is not an Edge", r.index),
        }
    }

    /// Mutate an edge's TShape in place, preserving Arc identity. OCCT
    /// `BRep_Builder::UpdateEdge` edits the edge TShape in place (BRepPrim adds
    /// pcurves to an edge that is already referenced by a wire); `Arc::make_mut`
    /// would clone-on-write and split the identity, leaving the wire edges on
    /// the old data. Safe because the primitive is being built sequentially.
    pub fn edge_mut_inplace(&mut self, r: Shape) -> &mut TEdgeData {
        // SAFETY: the caller holds &mut BRep (exclusive borrow of the tshape
        // slot) and is building the shape sequentially; every other reference
        // (face wires) observes the change, matching OCCT BRep_Builder::UpdateEdge
        // editing the edge TShape in place.
        let ptr = Arc::as_ptr(&self.tshapes[r.index]) as *mut TShape;
        unsafe {
            match &mut *ptr {
                TShape::Edge(e) => e,
                _ => panic!("edge_mut_inplace: Shape {} is not an Edge", r.index),
            }
        }
    }

    /// Mutate a wire's data.
    pub fn wire_mut(&mut self, r: Shape) -> &mut TWireData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Wire(w) => w,
            _ => panic!("wire_mut: Shape {} is not a Wire", r.index),
        }
    }

    /// Mutate a face's data.
    pub fn face_mut(&mut self, r: Shape) -> &mut TFaceData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Face(f) => f,
            _ => panic!("face_mut: Shape {} is not a Face", r.index),
        }
    }

    /// Mutate a shell's data.
    pub fn shell_mut(&mut self, r: Shape) -> &mut TShellData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Shell(s) => s,
            _ => panic!("shell_mut: Shape {} is not a Shell", r.index),
        }
    }

    /// Mutate a solid's data.
    pub fn solid_mut(&mut self, r: Shape) -> &mut TSolidData {
        match Arc::make_mut(&mut self.tshapes[r.index]) {
            TShape::Solid(s) => s,
            _ => panic!("solid_mut: Shape {} is not a Solid", r.index),
        }
    }

    /// Count vertex TShapes (includes orphan vertices not referenced by any edge).
    pub fn nb_vertices(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Vertex(_)))
            .count()
    }

    /// Count edge TShapes.
    pub fn nb_edges(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Edge(_)))
            .count()
    }

    /// Count shell TShapes.
    pub fn nb_shells(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Shell(_)))
            .count()
    }

    /// Count solid TShapes.
    pub fn nb_solids(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Solid(_)))
            .count()
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
        if mn.x.is_infinite() {
            None
        } else {
            Some([mn, mx])
        }
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
        if count == 0 {
            DVec3::ZERO
        } else {
            sum / count as f64
        }
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
// BRepTool trait �?OCCT BRep_Tool free-function equivalents
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool equivalent: parameter/tolerance/pcurve queries on a BRep.
///
/// In OCCT these are free functions (`BRep_Tool::Parameter(aV, aE, aF)` etc.).
/// Here they are methods on a `BRepTool` trait so the boolean pipeline can
/// be generic over the data source (real BRep or DS adaptor).
pub trait BRepTool {
    /// BRep_Tool::Pnt(aV) �?3D position of a vertex.
    fn vertex_position(&self, v: &Shape) -> DVec3;
    /// BRep_Tool::Tolerance(aV) �?3D tolerance of a vertex.
    fn vertex_tolerance(&self, v: &Shape) -> f64;
    /// BRep_Tool::Degenerated(aE).
    fn is_edge_degenerated(&self, e: &Shape) -> bool;
    /// TopExp: given one vertex of an edge, return the other vertex.
    fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape;
    /// First vertex of an edge (FORWARD orientation, OCCT TopExp::FirstVertex).
    fn first_vertex(&self, edge: &Shape) -> Shape;
    /// Last vertex of an edge (FORWARD orientation, OCCT TopExp::LastVertex).
    fn last_vertex(&self, edge: &Shape) -> Shape;
    /// TopExp::FirstVertex on an oriented edge �?canonical first for FORWARD,
    /// canonical last for REVERSED.  Matches OCCT's orientation-aware topology.
    fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape;
    /// BRep_Tool::Parameter(aV, aE, aF) �?vertex parameter on edge's pcurve.
    fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, face: &Shape) -> Option<f64>;
    /// BRep_Tool::CurveOnSurface(aE, aF) �?pcurve of edge on face.
    fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)>;
    /// BRep_Tool::Surface(aF) �?face surface (local coordinates, no Location applied).
    fn face_surface(&self, face: &Shape) -> Option<&Surface3>;
    /// BRep_Tool::Surface(aF) with Location applied �?returns world-coordinate surface.
    fn face_surface_world(&self, face: &Shape) -> Option<Surface3>;
    /// BRep_Tool::Curve(aE) with Location applied �?returns 3D curve and range in world coordinates.
    fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])>;
    /// UResolution: parameter tolerance in U direction (OCCT: BRepAdaptor_Surface::UResolution).
    fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64;
    /// VResolution: parameter tolerance in V direction.
    fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64;
    /// OCCT L204-207: vertex orientation (TopAbs_INTERNAL for split-edge interior vertices).
    /// Default: Forward (non-INTERNAL). Override when INTERNAL vertex data is available.
    fn vertex_orientation(&self, _v: &Shape) -> Orientation {
        Orientation::Forward
    }
    /// BRep_Tool::IsClosed(aE, aF) �?true when the edge appears twice on the face
    /// (periodic surface seam).  Checks for CurveOnClosedSurface representation.
    fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape) -> bool {
        self.curve_on_surface(edge, face).is_some()
            && self.curve_on_surface_second(edge, face).is_some()
    }
    /// Retrieve the second pcurve for periodic surfaces (CurveOnClosedSurface).
    /// Returns None for non-periodic faces or if only one pcurve exists.
    fn curve_on_surface_second(
        &self,
        edge: &Shape,
        face: &Shape,
    ) -> Option<&(Curve2d, f64, f64)> {
        None
    }

    // 鈹€鈹€ OCCT BRep_Tool / TopoDS_Shape convenience queries 鈹€鈹€

    /// OCCT TopoDS_Shape::Closed �?checks CLOSED flag on the TShape.
    /// For Shell, this is a simple flag check; use `is_shell_closed` for
    /// the full edge-count verification (BRepCheck_Shell).
    fn is_closed(&self, s: &Shape) -> bool {
        self.has_flag(s, tshape_flags::CLOSED)
    }

    /// BRep_Tool::SameParameter(aE) �?true when edge pcurves match 3D curve.
    fn edge_same_parameter(&self, e: &Shape) -> bool {
        self.edge_data(e)
            .map(|ed| ed.same_parameter)
            .unwrap_or(true)
    }

    /// BRep_Tool::SameRange(aE) �?true when edge pcurve ranges match 3D range.
    fn edge_same_range(&self, e: &Shape) -> bool {
        self.edge_data(e).map(|ed| ed.same_range).unwrap_or(true)
    }

    /// BRep_Tool::NaturalRestriction(aF) �?true when face surface bounds are
    /// determined by the underlying surface's natural domain.
    fn face_natural_restriction(&self, f: &Shape) -> bool {
        self.face_data(f)
            .map(|fd| fd.natural_restriction)
            .unwrap_or(true)
    }

    /// BRep_Tool::Curve(aE) �?raw 3D curve reference (no Location applied).
    /// Returns the curve data directly (OCCT-aligned: geometry on TShape).
    fn edge_curve_data(&self, e: &Shape) -> Option<Curve3> {
        self.edge_data(e).and_then(|ed| ed.curve.clone())
    }

    /// BRep_Tool::Range(aE) �?3D curve parameter range.
    fn edge_range(&self, e: &Shape) -> [f64; 2] {
        self.edge_data(e).map(|ed| ed.range).unwrap_or([0.0, 0.0])
    }

    /// BRep_Tool::Tolerance(s) �?geometric tolerance for any shape type.
    fn tolerance(&self, s: &Shape) -> f64;

    // 鈹€鈹€ Extension helpers (not in OCCT's BRep_Tool) 鈹€鈹€

    /// TopoDS_Shape::ShapeType �?returns the shape type.
    fn shape_type(&self, s: &Shape) -> ShapeType;

    /// Check CLOSED flag directly (bypasses trait default).
    fn has_flag(&self, s: &Shape, flag: u16) -> bool;

    /// Access Edge data (for implementing default methods).
    fn edge_data(&self, e: &Shape) -> Option<&TEdgeData>;

    /// Access Face data (for implementing default methods).
    fn face_data(&self, f: &Shape) -> Option<&TFaceData>;
}

impl BRepTool for BRep {
    fn vertex_position(&self, v: &Shape) -> DVec3 {
        let pt = self.vertex(v.clone()).point;
        self.get_location(v.location).transform_point3(pt)
    }

    fn vertex_tolerance(&self, v: &Shape) -> f64 {
        self.vertex(v.clone()).tolerance
    }

    fn is_edge_degenerated(&self, e: &Shape) -> bool {
        self.edge(e.clone()).degenerated
    }

    fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {
        let ed = self.edge(edge.clone());
        if ed.first.clone().index == v.index {
            ed.last.clone()
        } else {
            ed.first.clone()
        }
    }

    fn first_vertex(&self, edge: &Shape) -> Shape {
        self.edge(edge.clone()).first.clone()
    }

    fn last_vertex(&self, edge: &Shape) -> Shape {
        self.edge(edge.clone()).last.clone()
    }

    fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape {
        if orientation == Orientation::Reversed {
            self.last_vertex(edge)
        } else {
            self.first_vertex(edge)
        }
    }

    fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {
        self.edge(edge.clone()).vertex_params.get(&vertex.index).copied()
    }

    fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)> {
        self.edge(edge.clone()).pcurves.get(&face.index)
    }

    fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape) -> bool {
        let ed = self.edge(edge.clone());
        ed.representations.iter().any(|r| matches!(r, CurveRepresentation::CurveOnClosedSurface { face: f, .. } if *f == face.index))
            || ed.pcurves.contains_key(&face.index)
                && ed.pcurves.contains_key(&(face.index + self.nb_faces()))
    }

    fn curve_on_surface_second(
        &self,
        edge: &Shape,
        face: &Shape,
    ) -> Option<&(Curve2d, f64, f64)> {
        let shifted = face.index + self.nb_faces();
        self.edge(edge.clone()).pcurves.get(&shifted)
    }

    fn face_surface(&self, face: &Shape) -> Option<&Surface3> {
        self.face(face.clone()).surface.as_ref()
    }

    fn face_surface_world(&self, face: &Shape) -> Option<Surface3> {
        let fd = self.face(face.clone());
        let surface = fd.surface.as_ref()?.clone();
        let loc = self.get_location(face.location);
        if loc == glam::DAffine3::IDENTITY {
            Some(surface)
        } else {
            Some(crate::geom::transform_surface(&surface, &loc))
        }
    }

    fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {
        let ed = self.edge(edge.clone());
        let crv = ed.curve.as_ref()?.clone();
        let loc = self.get_location(edge.location);
        if loc == glam::DAffine3::IDENTITY {
            Some((crv, ed.range))
        } else {
            Some((crate::geom::transform_curve(&crv, &loc), ed.range))
        }
    }

    fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(face.clone()).surface.as_ref() {
            Some(surf) => u_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }

    fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(face.clone()).surface.as_ref() {
            Some(surf) => v_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }

    fn tolerance(&self, s: &Shape) -> f64 {
        match &*self.tshapes[s.index] {
            TShape::Vertex(vd) => vd.tolerance,
            TShape::Edge(ed) => ed.tolerance,
            TShape::Face(fd) => fd.tolerance,
            _ => 0.0,
        }
    }

    fn shape_type(&self, s: &Shape) -> ShapeType {
        s.shape_type_from_brep(self)
    }

    fn has_flag(&self, s: &Shape, flag: u16) -> bool {
        self.has_flag(s.clone(), flag)
    }

    fn edge_data(&self, e: &Shape) -> Option<&TEdgeData> {
        match &*self.tshapes[e.index] {
            TShape::Edge(ed) => Some(ed),
            _ => None,
        }
    }

    fn face_data(&self, f: &Shape) -> Option<&TFaceData> {
        match &*self.tshapes[f.index] {
            TShape::Face(fd) => Some(fd),
            _ => None,
        }
    }
}

/// OCCT GeomAdaptor_Surface::Load(S) — unwraps a Geom_RectangularTrimmedSurface
/// to its basis surface while keeping the trimmed parameter bounds
/// (GeomAdaptor_Surface.cxx L417-430). Returns the basis surface and the
/// parameter range (u1,u2,v1,v2) the adaptor reports.
pub fn surface_adaptor_basis_and_bounds(surf: &Surface3) -> (&Surface3, [f64; 4]) {
    match surf {
        Surface3::Trimmed(ts) => (ts.basis.as_ref(), ts.trim),
        s => (s, s.default_domain()),
    }
}

/// OCCT-aligned: GeomAdaptor_Curve::Resolution (GeomAdaptor_Curve.cxx
/// L1116-1148). Line -> R3d; Circle -> 2*asin(R3d/(2R)) when R > R3d/2 else
/// 2*PI; Ellipse -> R3d / MajorRadius; Bezier/BSpline -> the curve's
/// Resolution (Geom_BezierCurve/Geom_BSplineCurve, see math/bspl.rs);
/// TrimmedCurve unwraps to its basis (GeomAdaptor_Curve::load L252-255); all
/// remaining curve types (OffsetCurve/Hyperbola/Parabola/OtherCurve) ->
/// Precision::Parametric(R3d) = R3d * 0.01.
pub fn curve_resolution(curve: &Curve3, r3d: f64) -> f64 {
    match curve {
        Curve3::Line(_) => r3d,
        Curve3::Circle(c) => {
            let r = c.radius;
            if r > r3d / 2.0 {
                2.0 * (r3d / (2.0 * r)).asin()
            } else {
                2.0 * PI
            }
        }
        Curve3::Ellipse(e) => r3d / e.major_radius,
        Curve3::Bezier(b) => bezier_curve_resolution(b, r3d),
        Curve3::BSpline(b) => bspline_curve_resolution(b, r3d),
        Curve3::Trimmed(t) => curve_resolution(&t.curve, r3d),
        // GeomAbs_OffsetCurve / Hyperbola / Parabola / OtherCurve -> default.
        Curve3::Offset(_)
        | Curve3::Hyperbola(_)
        | Curve3::Parabola(_)
        | Curve3::CircularHelix(_)
        | Curve3::SineWave(_) => parametric_default(r3d),
    }
}

/// OCCT-aligned: GeomAdaptor_Surface::UResolution (GeomAdaptor_Surface.cxx
/// L1819-1896). Analytic branches (Torus/Sphere/Cylinder/Cone/Plane) are 1:1;
/// the adaptor wraps the stored surface and its parameter bounds come from
/// Geom_Surface::Bounds() (a Geom_RectangularTrimmedSurface reports its trim).
/// Bezier/BSpline call Geom_Surface::Resolution(Tol, Ures, Vres) = Tol *
/// UMaxDerivInv (BSplSLib::Resolution, see math/bspl.rs) and Offset recurses
/// the basis adaptor. The SurfaceOfExtrusion branch calls BasisCurve->
/// Resolution (GeomAdaptor_Curve::Resolution -> curve_resolution, the 1D
/// curve chain in math/bspl.rs).
pub fn u_resolution_for_surface(surf: &Surface3, tol3d: f64) -> f64 {
    let (basis, [_, _, v1, v2]) = surface_adaptor_basis_and_bounds(surf);
    let mut res = match basis {
        Surface3::Torus(t) => {
            let r = t.major_radius + t.minor_radius;
            if r > CONFUSION {
                tol3d / (2.0 * r)
            } else {
                0.0
            }
        }
        Surface3::Sphere(s) => {
            if s.radius > CONFUSION {
                tol3d / (2.0 * s.radius)
            } else {
                0.0
            }
        }
        Surface3::Cylinder(c) => {
            if c.radius > CONFUSION {
                tol3d / (2.0 * c.radius)
            } else {
                0.0
            }
        }
        Surface3::Cone(c) => {
            // OCCT L1856-1867: unbounded cone (V range > 1e10) -> unknown
            // resolution, return Precision::Parametric(R3d). Otherwise the
            // resolution is R3d / R where R is the larger of the VFirst/VLast
            // iso-circle radii.
            if v2 - v1 > 1e10 {
                return parametric_default(tol3d);
            }
            let rayon1 = c.radius_at_slant(v2);
            let rayon2 = c.radius_at_slant(v1);
            let r = rayon1.max(rayon2);
            return if r > CONFUSION { tol3d / r } else { 0.0 };
        }
        Surface3::Plane(_) => return tol3d,
        // SurfaceOfExtrusion L1824-1827: BasisCurve->Resolution(R3d) — the 1D
        // curve chain (curve_resolution, GeomAdaptor_Curve::Resolution).
        Surface3::LinearExtrusion(e) => return curve_resolution(&e.profile, tol3d),
        // BezierSurface/BSplineSurface L1872-1881: surface->Resolution(R3d,
        // Ures, Vres) -> Ures = R3d * UMaxDerivInv (BSplSLib::Resolution).
        Surface3::Bezier(b) => {
            let (ures, _) = bezier_surface_resolution(b, tol3d);
            return ures;
        }
        Surface3::BSpline(b) => {
            let (ures, _) = bspline_surface_resolution(b, tol3d);
            return ures;
        }
        // OffsetSurface L1882-1885: recurse BasisAdaptor->UResolution(R3d).
        Surface3::Offset(o) => return u_resolution_for_surface(o.basis.as_ref(), tol3d),
        // default L1886-1888: Precision::Parametric(R3d).
        _ => return parametric_default(tol3d),
    };
    // OCCT L1890-1895.
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * PI
    }
}

/// OCCT-aligned: GeomAdaptor_Surface::VResolution (GeomAdaptor_Surface.cxx
/// L1900-1962). Analytic branches (Torus/Sphere/Cylinder/Cone/Plane, plus
/// Extrusion returning R3d) are 1:1; Bezier/BSpline call Geom_Surface::
/// Resolution(Tol, Ures, Vres) = Tol * VMaxDerivInv (BSplSLib::Resolution,
/// see math/bspl.rs) and Offset recurses the basis adaptor. The
/// SurfaceOfRevolution branch calls BasisCurve->Resolution (GeomAdaptor_Curve::
/// Resolution -> curve_resolution, the 1D curve chain in math/bspl.rs).
pub fn v_resolution_for_surface(surf: &Surface3, tol3d: f64) -> f64 {
    let (basis, _) = surface_adaptor_basis_and_bounds(surf);
    let mut res = match basis {
        Surface3::Torus(t) => {
            let r = t.minor_radius;
            if r > CONFUSION {
                tol3d / (2.0 * r)
            } else {
                0.0
            }
        }
        Surface3::Sphere(s) => {
            if s.radius > CONFUSION {
                tol3d / (2.0 * s.radius)
            } else {
                0.0
            }
        }
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Plane(_) => return tol3d,
        // SurfaceOfRevolution L1906-1909: BasisCurve->Resolution(R3d) — the 1D
        // curve chain (curve_resolution, GeomAdaptor_Curve::Resolution).
        Surface3::Revolution(r) => return curve_resolution(&r.profile, tol3d),
        // SurfaceOfExtrusion L1928-1933: return R3d.
        Surface3::LinearExtrusion(_) => return tol3d,
        // BezierSurface/BSplineSurface L1934-1943: surface->Resolution(R3d,
        // Ures, Vres) -> Vres = R3d * VMaxDerivInv (BSplSLib::Resolution).
        Surface3::Bezier(b) => {
            let (_, vres) = bezier_surface_resolution(b, tol3d);
            return vres;
        }
        Surface3::BSpline(b) => {
            let (_, vres) = bspline_surface_resolution(b, tol3d);
            return vres;
        }
        // OffsetSurface L1944-1947: recurse BasisAdaptor->VResolution(R3d).
        Surface3::Offset(o) => return v_resolution_for_surface(o.basis.as_ref(), tol3d),
        // default L1948-1950: Precision::Parametric(R3d).
        _ => return parametric_default(tol3d),
    };
    // OCCT L1952-1962.
    if res <= 1.0 {
        2.0 * res.asin()
    } else {
        2.0 * PI
    }
}

// ---------------------------------------------------------------------------
// BRepBuilder �?OCCT BRep_Builder equivalent for incrementally constructing BRep
// ---------------------------------------------------------------------------

// 鈹€鈹€ Backward-compat flat-index access methods 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// These provide the old brep.vertices/edges/solids/geom patterns
// for modules that haven't been migrated yet.
impl BRep {
    /// Collect vertex points in tshape order (like old brep.vertices).
    pub fn flat_vertices(&self) -> Vec<DVec3> {
        self.tshapes
            .iter()
            .filter_map(|ts| {
                if let TShape::Vertex(vd) = &**ts {
                    Some(vd.point)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Collect edge endpoint pairs in tshape order (like old brep.edges).
    pub fn flat_edges(&self) -> Vec<(usize, usize)> {
        self.tshapes
            .iter()
            .filter_map(|ts| {
                if let TShape::Edge(ed) = &**ts {
                    Some((ed.first.clone().index, ed.last.clone().index))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Count vertex TShapes.
    pub fn vertex_count(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!((&**ts).as_ref(), TShape::Vertex(_)))
            .count()
    }

    /// Count edge TShapes.
    pub fn edge_count(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!((&**ts).as_ref(), TShape::Edge(_)))
            .count()
    }

    /// Count face TShapes.
    pub fn face_count(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!((&**ts).as_ref(), TShape::Face(_)))
            .count()
    }

    /// Count solid TShapes.
    pub fn solid_count(&self) -> usize {
        self.tshapes
            .iter()
            .filter(|ts| matches!((&**ts).as_ref(), TShape::Solid(_)))
            .count()
    }

    /// Check if any solid TShapes exist.
    pub fn has_solids(&self) -> bool {
        self.tshapes
            .iter()
            .any(|ts| matches!(ts.as_ref(), &TShape::Solid(_)))
    }

    /// Get vertex point by tshape index.
    pub fn vertex_point(&self, idx: usize) -> Option<DVec3> {
        self.tshapes.get(idx).and_then(|ts| {
            if let TShape::Vertex(vd) = &**ts {
                Some(vd.point)
            } else {
                None
            }
        })
    }

    /// Add a new edge from flat indices (creates Shape for each vertex).
    /// Returns the tshape index of the new edge.
    pub fn add_edge_flat(
        &mut self,
        start_idx: usize,
        end_idx: usize,
        curve: Option<Curve3>,
        range: [f64; 2],
    ) -> usize {
        let first = self
            .tshapes
            .get(start_idx)
            .map(|ts| Shape {
                data: ts.clone(),
                index: start_idx,
                orientation: Orientation::Forward,
                location: 0,
            })
            .unwrap_or(Shape::null());
        let last = self
            .tshapes
            .get(end_idx)
            .map(|ts| Shape {
                data: ts.clone(),
                index: end_idx,
                orientation: Orientation::Forward,
                location: 0,
            })
            .unwrap_or(Shape::null());
        let sr = self.add_tedge(curve, first, last, range);
        sr.index
    }
}

pub struct BRepBuilder {
    pub root: Option<Shape>,
    vertex_cache: Vec<[f64; 3]>,
}

// --- Backward-compat flat-index methods ---
// These emulate the old brep.solids, brep.edges, brep.vertices, brep.geom access
// for modules not yet migrated to the tshape-based API.
impl BRep {
    /// Return flat list of vertex points in tshape order.
    pub fn vertices(&self) -> Vec<crate::Vertex> {
        self.tshapes
            .iter()
            .filter_map(|ts| {
                if let TShape::Vertex(vd) = ts.as_ref() {
                    Some(crate::Vertex { point: vd.point })
                } else {
                    None
                }
            })
            .collect()
    }

    /// Return flat list of edge topologies in tshape order (start/end are old-style vertex indices).
    pub fn edges(&self) -> Vec<crate::topo::topology::Edge> {
        self.tshapes
            .iter()
            .filter_map(|ts| {
                if let TShape::Edge(ed) = ts.as_ref() {
                    Some(crate::topo::topology::Edge {
                        start: ed.first.clone().index,
                        end: ed.last.clone().index,
                    })
                } else {
                    None
                }
            })
            .collect()
    }

    /// Return flat list of solid topologies (materialized from tshapes).
    pub fn solids(&self) -> Vec<crate::topo::topology::Solid> {
        let mut out = Vec::new();
        for ts in &self.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                let mut shells = Vec::new();
                for shell_sr in &sd.shells {
                    if shell_sr.index >= self.tshapes.len() {
                        continue;
                    }
                    if let TShape::Shell(shd) = &*self.tshapes[shell_sr.index] {
                        let mut faces = Vec::new();
                        for face_sr in &shd.faces {
                            if face_sr.index >= self.tshapes.len() {
                                continue;
                            }
                            if let TShape::Face(fd) = &*self.tshapes[face_sr.index] {
                                let outer_wire = {
                                    let mut edges = Vec::new();
                                    if fd.outer_wire.index < self.tshapes.len() {
                                        if let TShape::Wire(wd) =
                                            &*self.tshapes[fd.outer_wire.index]
                                        {
                                            for esr in &wd.edges {
                                                edges.push(crate::topo::topology::WireEdge {
                                                    idx: esr.index,
                                                    forward: true,
                                                });
                                            }
                                        }
                                    }
                                    crate::topo::topology::Wire { edges }
                                };
                                let inner_wires = fd
                                    .inner_wires
                                    .iter()
                                    .map(|iw_sr| {
                                        let mut edges = Vec::new();
                                        if iw_sr.index < self.tshapes.len() {
                                            if let TShape::Wire(wd) = &*self.tshapes[iw_sr.index] {
                                                for esr in &wd.edges {
                                                    edges.push(crate::topo::topology::WireEdge {
                                                        idx: esr.index,
                                                        forward: true,
                                                    });
                                                }
                                            }
                                        }
                                        crate::topo::topology::Wire { edges }
                                    })
                                    .collect();
                                let normal = fd
                                    .surface
                                    .as_ref()
                                    .map(|s| crate::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
                                    .unwrap_or_default();
                                faces.push(crate::topo::topology::Face {
                                    outer_wire,
                                    inner_wires,
                                    normal,
                                    triangles: Vec::new(),
                                    sample_point: fd.sample_point,
                                    mesh_dirty: true,
                                    surface_idx: None,
                                });
                            }
                        }
                        shells.push(crate::topo::topology::Shell { faces });
                    }
                }
                out.push(crate::topo::topology::Solid { shells });
            }
        }
        out
    }

    /// Return old-style geom store (empty �?geom now embedded in TShapes).
    pub fn geom(&self) -> crate::GeomStore {
        crate::GeomStore::new()
    }

    /// OCCT BRepLib::SameParameter.
    /// For each edge: compute pcurves on adjacent faces (if missing) and update
    /// edge tolerance to max(sampled_3D_2D_deviation, CONFUSION).
    pub fn same_parameter(&mut self) {
        use crate::geom::{
            Curve2d, Curve2dEval, CurveEval, Line2d, Surface3, SurfaceEval, any_perpendicular,
        };
        use glam::DVec2;
        // Build edge �?faces map
        let mut edge_faces: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (ti, ts) in self.tshapes.iter().enumerate() {
            if let TShape::Face(fd) = &**ts {
                let process_wire =
                    |sr: &Shape, map: &mut std::collections::HashMap<usize, Vec<usize>>| {
                        if sr.index < self.tshapes.len() {
                            if let TShape::Wire(wd) = &*self.tshapes[sr.index] {
                                for esr in &wd.edges {
                                    if esr.index < self.tshapes.len() {
                                        map.entry(esr.index).or_default().push(ti);
                                    }
                                }
                            }
                        }
                    };
                process_wire(&fd.outer_wire, &mut edge_faces);
                for iw_sr in &fd.inner_wires {
                    process_wire(iw_sr, &mut edge_faces);
                }
            }
        }
        // Process each edge
        for (ei, face_indices) in &edge_faces {
            if *ei >= self.tshapes.len() {
                continue;
            }
            let edge_data = match &*self.tshapes[*ei] {
                TShape::Edge(ed) => ed.clone(),
                _ => continue,
            };
            let Some(ref curve3) = edge_data.curve else {
                continue;
            };
            let t_range = edge_data.range;
            let mut max_dev = 0.0f64;
            for &fi in face_indices {
                let face_data = match &*self.tshapes[fi] {
                    TShape::Face(fd) => fd.clone(),
                    _ => continue,
                };
                let Some(ref surf) = face_data.surface else {
                    continue;
                };
                // Compute pcurve if missing (planar surfaces only)
                let has_pcurve = edge_data.pcurves.contains_key(&fi);
                if !has_pcurve {
                    if let Surface3::Plane(p) = surf {
                        let p0 = curve3.point_at(t_range[0]);
                        let p1 = curve3.point_at(t_range[1]);
                        // Project 3D points to UV on plane
                        let u_axis = crate::geom::any_perpendicular(p.normal);
                        let v_axis = p.normal.cross(u_axis);
                        let to_uv = |pt: DVec3| -> DVec2 {
                            let local = pt - p.origin;
                            DVec2::new(local.dot(u_axis), local.dot(v_axis))
                        };
                        let uv0 = to_uv(p0);
                        let uv1 = to_uv(p1);
                        let du = uv1.x - uv0.x;
                        let dv = uv1.y - uv0.y;
                        let duv_len = (du * du + dv * dv).sqrt();
                        if duv_len > 1e-15 {
                            let dir2d = glam::DVec2::new(du, dv) / duv_len;
                            let pc = Curve2d::Line(Line2d {
                                origin: uv0,
                                direction: dir2d,
                            });
                            // Insert pcurve into edge
                            if let TShape::Edge(ed) = Arc::make_mut(&mut self.tshapes[*ei]) {
                                ed.pcurves.insert(fi, (pc.clone(), t_range[0], t_range[1]));
                            }
                            // Sample deviation
                            let n_samples = 7;
                            for si in 0..=n_samples {
                                let t = t_range[0]
                                    + (t_range[1] - t_range[0]) * si as f64 / n_samples as f64;
                                let p3d = curve3.point_at(t);
                                let uv = pc.point_at(t);
                                if !uv.is_finite() {
                                    continue;
                                }
                                let p_surf = surf.point_at(uv.x, uv.y);
                                let dev = (p3d - p_surf).length();
                                if dev > max_dev {
                                    max_dev = dev;
                                }
                            }
                        }
                    }
                }
            }
            // Update edge tolerance
            let new_tol = max_dev.max(crate::core::precision::CONFUSION);
            if let TShape::Edge(ed) = Arc::make_mut(&mut self.tshapes[*ei]) {
                ed.tolerance = ed.tolerance.max(new_tol);
                ed.same_parameter = true;
            }
        }
    }
}

impl BRepBuilder {
    pub fn new() -> Self {
        Self {
            root: None,
            vertex_cache: Vec::new(),
        }
    }

    fn find_or_add_vertex(&mut self, brep: &mut BRep, pt: DVec3) -> Shape {
        for (i, &cached) in self.vertex_cache.iter().enumerate() {
            let dp = DVec3::new(cached[0], cached[1], cached[2]) - pt;
            if dp.length_squared() < 1e-30 {
                return Shape::synthetic(i, Orientation::Forward);
            }
        }
        let r = brep.add_tvertex(pt);
        self.vertex_cache.push([pt.x, pt.y, pt.z]);
        r
    }

    pub fn build_unit_cube(&mut self, brep: &mut BRep) {
        // 8 vertices
        let v: Vec<Shape> = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        .iter()
        .map(|&p| self.find_or_add_vertex(brep, DVec3::new(p[0], p[1], p[2])))
        .collect();

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
                    let va = v[a].clone();
                    let vb = v[b].clone();
                    brep.add_tedge(
                        None,
                        va,
                        Shape::synthetic(vb.index, Orientation::Reversed),
                        [0.0, 1.0],
                    )
                    .index
                });
                let orient = if a < b {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
                face_edges.push(Shape::synthetic(e_idx, orient));
            }
            edge_for_face.push(face_edges);
        }

        // Build 6 faces, collecting their refs for shell building
        let mut face_refs = Vec::new();
        for (i, face_edges) in edge_for_face.into_iter().enumerate() {
            let wire = brep.add_twire(face_edges);
            let face = brep.add_tface(None, wire, vec![], None, None, vec![], true);
            // Outer normal: orient face based on the face index
            // For a unit cube at origin, faces 0,2,4 have inward normals, orient REVERSED
            // faces 1,3,5 have outward normals, orient FORWARD
            let orient = match i {
                0 => Orientation::Reversed, // Z-
                1 => Orientation::Forward,  // Z+
                2 => Orientation::Reversed, // Y-
                3 => Orientation::Forward,  // Y+
                4 => Orientation::Reversed, // X-
                5 => Orientation::Forward,  // X+
                _ => unreachable!(),
            };
            face_refs.push(Shape::synthetic(face.index, orient));
        }

        // Build shell from the 6 oriented faces
        let shell = brep.add_tshell(face_refs);
        // Build solid from the shell
        let solid = brep.add_tsolid(vec![shell]);
        self.root = Some(solid);
    }

    // -----------------------------------------------------------------------
    // Extended BRepBuilder API �?OCCT BRep_Builder equivalent
    // -----------------------------------------------------------------------

    /// Add a vertex with tolerance.
    pub fn add_vertex(&mut self, brep: &mut BRep, pt: DVec3, tol: f64) -> Shape {
        let r = brep.add_tvertex(pt);
        brep.vertex_mut(r.clone()).tolerance = tol;
        r
    }

    /// Update vertex tolerance (max with existing).
    pub fn update_vertex_tolerance(&mut self, brep: &mut BRep, v: Shape, tol: f64) {
        let vd = brep.vertex_mut(v);
        vd.tolerance = vd.tolerance.max(tol);
    }

    /// Add an edge with curve, vertices, and range.
    pub fn add_edge(
        &mut self,
        brep: &mut BRep,
        curve: Option<Curve3>,
        v1: Shape,
        v2: Shape,
        range: [f64; 2],
    ) -> Shape {
        brep.add_tedge(curve, v1, v2, range)
    }

    /// Add a pcurve to an edge for a specific face.
    pub fn add_pcurve(
        &mut self,
        brep: &mut BRep,
        edge: Shape,
        face: Shape,
        pc: Curve2d,
        t1: f64,
        t2: f64,
    ) {
        brep.edge_mut(edge).pcurves.insert(face.index, (pc, t1, t2));
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, theTol) �?update edge tolerance.
    pub fn update_edge_tolerance(&mut self, brep: &mut BRep, edge: Shape, tol: f64) {
        let ed = brep.edge_mut(edge);
        ed.tolerance = ed.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, aC2d, aF, theTol) �?set pcurve on face.
    pub fn update_edge_pcurve(
        &mut self,
        brep: &mut BRep,
        edge: Shape,
        pcurve: Curve2d,
        face: Shape,
        tol: f64,
    ) {
        let ed = brep.edge_mut(edge);
        let (ta, tb) = pc_parameter_range(&pcurve);
        ed.pcurves.insert(face.index, (pcurve.clone(), ta, tb));
        ed.representations
            .push(CurveRepresentation::CurveOnSurface {
                face: face.index,
                pcurve,
                range: [ta, tb],
            });
        ed.tolerance = ed.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, aC1, aC2, aF, theTol) — set two
    /// pcurves of an edge on the same face. The edge lies on the closing curve
    /// (seam) of a closed surface and carries a BRep_CurveOnClosedSurface
    /// representation (BRepPrim_OneAxis::LateralFace L434-438).
    ///
    /// `aFirst`/`aLast` are the edge's first and last parameter (BRep_Tool::Range);
    /// both pcurves are evaluated over the same parameter interval as the edge
    /// (same-parameter edge). The second pcurve is stored under the shifted key
    /// `face.index + nb_faces` (rcad convention read by
    /// [`BRepTool::curve_on_surface_second`]).
    pub fn update_edge_pcurve_closed(
        &mut self,
        brep: &mut BRep,
        edge: Shape,
        pcurve1: Curve2d,
        pcurve2: Curve2d,
        face: Shape,
        a_first: f64,
        a_last: f64,
        tol: f64,
    ) {
        let nb_faces = brep.nb_faces();
        let ed = brep.edge_mut(edge);
        ed.pcurves
            .insert(face.index, (pcurve1.clone(), a_first, a_last));
        ed.pcurves.insert(
            face.index + nb_faces,
            (pcurve2.clone(), a_first, a_last),
        );
        ed.representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: face.index,
                pcurve1,
                pcurve2,
                range: [a_first, a_last],
            });
        ed.tolerance = ed.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(aE, aC3d) �?set 3D curve.
    pub fn update_edge_curve3d(
        &mut self,
        brep: &mut BRep,
        edge: Shape,
        curve: usize,
        location: u32,
    ) {
        let ed = brep.edge_mut(edge);
        ed.representations
            .push(CurveRepresentation::Curve3D { curve, location });
    }

    /// Set vertex parameter on an edge's pcurve.
    pub fn set_vertex_param(
        &mut self,
        brep: &mut BRep,
        edge: Shape,
        vertex: Shape,
        param: f64,
    ) {
        brep.edge_mut(edge)
            .vertex_params
            .insert(vertex.index, param);
    }

    /// Set degenerated flag on an edge.
    pub fn set_edge_degenerated(&mut self, brep: &mut BRep, edge: Shape, flag: bool) {
        brep.edge_mut(edge).degenerated = flag;
    }

    /// OCCT BRep_Builder::Add(aW).Closed(aW) �?mark wire as closed.
    pub fn close_wire(&mut self, brep: &mut BRep, wire: Shape) {
        brep.wire_mut(wire).flags |= tshape_flags::CLOSED;
    }

    /// OCCT BRep_Builder::Add(aShell).Closed(aShell) �?mark shell as closed via Closed flag.
    pub fn close_shell(&mut self, brep: &mut BRep, shell: Shape) {
        brep.shell_mut(shell).flags |= tshape_flags::CLOSED;
    }

    /// Make a wire (empty container).
    pub fn make_wire(&mut self, brep: &mut BRep) -> Shape {
        brep.add_twire(vec![])
    }

    /// Add an edge to a wire.
    pub fn add_to_wire(&mut self, brep: &mut BRep, wire: Shape, edge: Shape) {
        let wd = brep.wire_mut(wire);
        wd.edges.push(edge.clone());
        wd.my_shapes.push(edge);
    }

    /// Build a wire from edges.
    pub fn build_wire(&mut self, brep: &mut BRep, edges: Vec<Shape>) -> Shape {
        brep.add_twire(edges)
    }

    /// Make a face from a surface and outer wire.
    pub fn make_face(
        &mut self,
        brep: &mut BRep,
        surface: Option<Surface3>,
        outer_wire: Shape,
    ) -> Shape {
        brep.add_tface(surface, outer_wire, vec![], None, None, vec![], true)
    }

    /// Add an inner wire to a face.
    pub fn add_to_face(&mut self, brep: &mut BRep, face: Shape, inner_wire: Shape) {
        let fd = brep.face_mut(face);
        fd.inner_wires.push(inner_wire.clone());
        fd.my_shapes.push(inner_wire);
    }

    /// Add an internal vertex to a face.
    pub fn add_internal_vertex(&mut self, brep: &mut BRep, face: Shape, v: Shape) {
        let fd = brep.face_mut(face);
        fd.internal_vertices.push(v.clone());
        fd.my_shapes.push(v);
    }

    /// Add an edge with section-curve semantics (MakeSectEdge equivalent).
    /// Creates an edge with pcurves for both faces.
    pub fn add_section_edge(
        &mut self,
        brep: &mut BRep,
        curve: Option<Curve3>,
        v1: Shape,
        v2: Shape,
        range: [f64; 2],
        pc_a: Option<&Curve2d>,
        face_a: Option<Shape>,
        pc_b: Option<&Curve2d>,
        face_b: Option<Shape>,
    ) -> Shape {
        let e = brep.add_tedge(curve, v1, v2, range);
        if let (Some(pc), Some(fa)) = (pc_a, face_a) {
            let (t1, t2) = pc_parameter_range(pc);
            brep.edge_mut(e.clone())
                .pcurves
                .insert(fa.index, (pc.clone(), t1, t2));
        }
        if let (Some(pc), Some(fb)) = (pc_b, face_b) {
            let (t1, t2) = pc_parameter_range(pc);
            brep.edge_mut(e.clone())
                .pcurves
                .insert(fb.index, (pc.clone(), t1, t2));
        }
        e
    }

    /// Make a shell (empty container).
    pub fn make_shell(&mut self, brep: &mut BRep) -> Shape {
        brep.add_tshell(vec![])
    }

    /// Add a face to a shell.
    pub fn add_to_shell(&mut self, brep: &mut BRep, shell: Shape, face: Shape) {
        let sd = brep.shell_mut(shell);
        sd.faces.push(face.clone());
        sd.my_shapes.push(face);
    }

    /// Make a solid from shells.
    pub fn make_solid(&mut self, brep: &mut BRep, shells: Vec<Shape>) -> Shape {
        brep.add_tsolid(shells)
    }

    /// Make a compsolid from solids.
    pub fn make_compsolid(&mut self, brep: &mut BRep, solids: Vec<Shape>) -> Shape {
        brep.add_tcompsolid(solids)
    }

    /// Make a compound from shapes.
    pub fn make_compound(&mut self, brep: &mut BRep, shapes: Vec<Shape>) -> Shape {
        brep.add_tcompound(shapes)
    }

    /// Add a shape to an existing compound.
    pub fn add_to_compound(&self, brep: &mut BRep, compound: Shape, shape: Shape) {
        let ts = Arc::make_mut(&mut brep.tshapes[compound.index]);
        match ts {
            TShape::Compound(shapes) => shapes.push(shape),
            _ => panic!("add_to_compound: shape is not a Compound"),
        }
    }

    // ======================================================================
    // Additional BRep_Builder methods (completing coverage)
    // ======================================================================

    /// OCCT BRep_Builder::Range(aE, First, Last) �?set edge parametric range.
    /// Updates the edge-level range field.
    pub fn set_edge_range(&mut self, brep: &mut BRep, edge: Shape, first: f64, last: f64) {
        brep.edge_mut(edge).range = [first, last];
    }

    /// OCCT BRep_Builder::SameParameter(aE, theFlag) �?set same-parameter flag.
    pub fn set_edge_same_parameter(&mut self, brep: &mut BRep, edge: Shape, flag: bool) {
        brep.edge_mut(edge).same_parameter = flag;
    }

    /// OCCT BRep_Builder::SameRange(aE, theFlag) �?set same-range flag.
    pub fn set_edge_same_range(&mut self, brep: &mut BRep, edge: Shape, flag: bool) {
        brep.edge_mut(edge).same_range = flag;
    }

    /// OCCT BRep_Builder::NaturalRestriction(aF, theFlag) �?set natural-restriction flag.
    pub fn set_face_natural_restriction(&mut self, brep: &mut BRep, face: Shape, flag: bool) {
        brep.face_mut(face).natural_restriction = flag;
    }

    /// OCCT BRep_Builder::UpdateFace(aF, theTol) �?update face tolerance.
    pub fn update_face_tolerance(&mut self, brep: &mut BRep, face: Shape, tol: f64) {
        let fd = brep.face_mut(face);
        fd.tolerance = fd.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::UpdateVertex(aV, aP, theTol) �?set vertex point and tolerance.
    pub fn update_vertex_point(&mut self, brep: &mut BRep, vertex: Shape, pt: DVec3, tol: f64) {
        let vd = brep.vertex_mut(vertex);
        vd.point = pt;
        vd.tolerance = vd.tolerance.max(tol);
    }

    /// Remove an edge from a wire (OCCT BRep_Builder::Remove).
    pub fn remove_from_wire(&mut self, brep: &mut BRep, wire: Shape, edge: Shape) {
        let wd = brep.wire_mut(wire);
        wd.edges
            .retain(|e| e.index != edge.index || e.orientation != edge.orientation);
        wd.my_shapes
            .retain(|s| s.index != edge.index || s.orientation != edge.orientation);
    }

    /// Remove an inner wire from a face (OCCT BRep_Builder::Remove).
    pub fn remove_from_face(&mut self, brep: &mut BRep, face: Shape, inner_wire: Shape) {
        let fd = brep.face_mut(face);
        fd.inner_wires.retain(|w| w.index != inner_wire.index);
        fd.my_shapes.retain(|s| s.index != inner_wire.index);
    }

    /// Remove a face from a shell (OCCT BRep_Builder::Remove).
    pub fn remove_from_shell(&mut self, brep: &mut BRep, shell: Shape, face: Shape) {
        let sd = brep.shell_mut(shell);
        sd.faces.retain(|f| f.index != face.index);
        sd.my_shapes.retain(|s| s.index != face.index);
    }

    /// OCCT BRep_Builder::Transfert(aEin, aEout) �?copy 3D curve from one edge to another.
    /// Copies the Curve3D representation (first one found) from edge_in to edge_out.
    pub fn transfert_edge_curve(&mut self, brep: &mut BRep, edge_in: Shape, edge_out: Shape) {
        let curve_clone = brep.edge(edge_in.clone()).curve.clone();
        let pcurves_clone = brep.edge(edge_in).pcurves.clone();
        if let Some(curve) = curve_clone {
            brep.edge_mut(edge_out.clone()).curve = Some(curve);
        }
        if !pcurves_clone.is_empty() {
            brep.edge_mut(edge_out).pcurves.extend(pcurves_clone);
        }
    }

    /// OCCT BRep_Builder::Transfert(aEin, aEout, aVin, aVout) �?copy vertex parameter
    /// from vertex_in on edge_in to vertex_out on edge_out.
    pub fn transfert_vertex_param(
        &mut self,
        brep: &mut BRep,
        edge_in: Shape,
        vertex_in: Shape,
        edge_out: Shape,
        vertex_out: Shape,
    ) {
        let param = brep
            .edge(edge_in)
            .vertex_params
            .get(&vertex_in.index)
            .copied();
        if let Some(param) = param {
            brep.edge_mut(edge_out)
                .vertex_params
                .insert(vertex_out.index, param);
        }
    }

    /// OCCT BRep_Builder::Degenerated(aE, true) �?set degenerated flag AND clear 3D curve.
    /// OCCT removes the 3D curve when marking an edge as degenerated.
    pub fn set_edge_degenerated_with_clear(&mut self, brep: &mut BRep, edge: Shape, flag: bool) {
        let ed = brep.edge_mut(edge);
        ed.degenerated = flag;
        if flag {
            ed.curve = None;
        }
    }
}

/// Get the parameter range for a Curve2d (Trimmed �?stored range, Circle �?[0, 2蟺]).
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
    /// OCCT TopoDS_TShape::Locked — returns true if the LOCKED flag is set.
    pub fn locked(&self) -> bool {
        let flags = match self {
            TShape::Vertex(v) => v.flags,
            TShape::Edge(e) => e.flags,
            TShape::Wire(w) => w.flags,
            TShape::Face(f) => f.flags,
            TShape::Shell(s) => s.flags,
            TShape::Solid(s) => s.flags,
            TShape::CompSolid(_) => 0,
            TShape::Compound(_) => 0,
        };
        flags & tshape_flags::LOCKED != 0
    }

    /// OCCT TopoDS_TShape::Locked(theFlag) — sets or clears the LOCKED flag.
    pub fn set_locked(&mut self, the_flag: bool) {
        let flags = match self {
            TShape::Vertex(v) => &mut v.flags,
            TShape::Edge(e) => &mut e.flags,
            TShape::Wire(w) => &mut w.flags,
            TShape::Face(f) => &mut f.flags,
            TShape::Shell(s) => &mut s.flags,
            TShape::Solid(s) => &mut s.flags,
            TShape::CompSolid(_) | TShape::Compound(_) => return,
        };
        if the_flag {
            *flags |= tshape_flags::LOCKED;
        } else {
            *flags &= !tshape_flags::LOCKED;
        }
    }

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
impl ShapeType {
    /// OCCT BOPDS_Tools::HasBRep — true for Vertex/Edge/Face (shapes with geometric data).
    pub fn has_brep(&self) -> bool {
        matches!(self, ShapeType::Vertex | ShapeType::Edge | ShapeType::Face)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::*;

    #[test]
    fn test_orientation_values() {
        assert!(Orientation::Forward.is_forward());
        assert!(!Orientation::Reversed.is_forward());
        assert!(!Orientation::Internal.is_forward());
        assert!(!Orientation::External.is_forward());
    }

    #[test]
    fn test_shape_ref_construction() {
        let r = Shape::synthetic(5, Orientation::Forward);
        assert_eq!(r.index, 5);
        assert_eq!(r.orientation, Orientation::Forward);

        let r2 = Shape::synthetic(3, Orientation::Reversed);
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
        assert_eq!(brep.vertex(v.clone()).point, DVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_vertex_sharing() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::ZERO);
        // OCCT-aligned: same position returns the same TShape (identity sharing)
        assert_eq!(brep.tshapes.len(), 1);
        assert_eq!(v0.index, v1.index);
    }

    #[test]
    fn test_edge_creation() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let v0_idx = v0.index;
        let v1_idx = v1.index;
        let e = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
        assert_eq!(brep.tshapes.len(), 3);
        let ed = brep.edge(e.clone());
        assert_eq!(ed.first.clone().index, v0_idx);
        assert_eq!(ed.last.clone().index, v1_idx);
    }

    #[test]
    fn test_edge_shares_vertex_tshape() {
        // Two edges sharing the same vertex TShape
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let v2 = brep.add_tvertex(DVec3::new(1.0, 1.0, 0.0));
        let v1_idx = v1.index;

        let e0 = brep.add_tedge(None, v0, v1.clone(), [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v1, v2, [0.0, 1.0]);

        // Both edges reference v1 at index 1
        assert_eq!(brep.edge(e0.clone()).last.index, v1_idx);
        assert_eq!(brep.edge(e1.clone()).first.index, v1_idx);
        // Same TShape identity (v1)
        assert!(Arc::ptr_eq(
            &brep.tshapes[brep.edge(e0.clone()).last.index],
            &brep.tshapes[brep.edge(e1.clone()).first.index]
        ));
    }

    #[test]
    fn test_wire_and_face() {
        let mut brep = BRep::new();
        let v = (0..4)
            .map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0)))
            .collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0].clone(), v[1].clone(), [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v[1].clone(), v[2].clone(), [0.0, 1.0]);
        let e2 = brep.add_tedge(None, v[2].clone(), v[3].clone(), [0.0, 1.0]);
        let wire = brep.add_twire(vec![e0, e1, e2]);
        let face = brep.add_tface(None, wire, vec![], None, None, vec![], true);
        let fd = brep.face(face.clone());
        assert_eq!(brep.tshapes.len(), 9); // 4V + 3E + 1W + 1F
        assert!(fd.inner_wires.is_empty());
        assert!(fd.sample_point.is_none());
    }

    #[test]
    fn test_shape_type_discrimination() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        assert_eq!(brep.tshapes[v.index].shape_type(), ShapeType::Vertex);

        let e = brep.add_tedge(None, v.clone(), v, [0.0, 1.0]);
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
        let e = brep.add_tedge(
            None,
            v0,
            Shape::synthetic(v1.index, Orientation::Reversed),
            [0.0, 1.0],
        );
        let ed = brep.edge(e.clone());
        assert_eq!(ed.last.clone().orientation, Orientation::Reversed);
    }

    #[test]
    fn test_clone_preserves_shape_identity() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let _v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(None, v0.clone(), v0, [0.0, 1.0]);

        let cloned = brep.clone();
        // Same TShape identity in clone (Arc::ptr_eq across clone)
        assert_eq!(cloned.tshapes.len(), brep.tshapes.len());
        assert!(Arc::ptr_eq(
            &cloned.tshapes[e.index],
            &brep.tshapes[e.index]
        ));
    }

    #[test]
    fn test_serde_roundtrip() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        let e = brep.add_tedge(None, v.clone(), v, [0.0, 1.0]);

        let json = serde_json::to_string(&brep).unwrap();
        let restored: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.tshapes.len(), 2);
        assert_eq!(
            restored.vertex(Shape::synthetic(0, Orientation::Forward)).point,
            DVec3::new(1.0, 2.0, 3.0)
        );
        assert_eq!(
            restored.edge(Shape::synthetic(1, Orientation::Forward)).range,
            [0.0, 1.0]
        );
    }

    #[test]
    fn test_brep_builder_unit_cube() {
        let (_brep, root) = BRep::build_unit_cube();
        assert_eq!(root.orientation, Orientation::Forward);
    }

    #[test]
    fn test_wire_edge_orientation() {
        let mut brep = BRep::new();
        let v = (0..4)
            .map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0)))
            .collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0].clone(), v[1].clone(), [0.0, 1.0]);
        let e0_idx = e0.index;

        // Wire with forward edge
        let wire_fwd = brep.add_twire(vec![e0]);
        assert_eq!(brep.tshapes[wire_fwd.index].shape_type(), ShapeType::Wire);

        // Wire with reversed edge
        let e0_rev = Shape::synthetic(e0_idx, Orientation::Reversed);
        let e0_rev_idx = e0_rev.index;
        let wire_rev = brep.add_twire(vec![e0_rev]);
        if let TShape::Wire(ref wd) = *brep.tshapes[wire_rev.index] {
            assert_eq!(wd.edges[0].orientation, Orientation::Reversed);
        } else {
            panic!("expected Wire");
        }

        // Same TShape for the edge, different orientation on the reference
        assert!(Arc::ptr_eq(
            &brep.tshapes[e0_idx],
            &brep.tshapes[e0_rev_idx]
        ));
    }

    #[test]
    fn test_topods_roundtrip_sa() {
        // Build a unit cube via BRepBuilder and verify basic structure counts.
        // Note: build_unit_cube creates topological-only BRep (no face surfaces),
        // so surface_area computation is not applicable here.
        let mut brep = crate::BRep::new();
        let mut builder = BRepBuilder::new();
        builder.build_unit_cube(&mut brep);
        // Sanity-check: cube has 8 vertices, 12 edges, 6 faces
        let n_verts = brep
            .tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Vertex(_)))
            .count();
        let n_edges = brep
            .tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Edge(_)))
            .count();
        let n_faces = brep
            .tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
            .count();
        assert_eq!(n_verts, 8, "cube: 8 vertices");
        assert_eq!(n_edges, 12, "cube: 12 edges");
        assert_eq!(n_faces, 6, "cube: 6 faces");
    }

    #[test]
    fn test_face_natural_restriction_default_true() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        let w = brep.add_twire(vec![]);
        // default (no explicit nr) �?true via add_tface with natural_restriction=true
        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        let fd = brep.face(f.clone());
        assert!(fd.natural_restriction);
    }

    #[test]
    fn test_face_natural_restriction_false() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], false);
        let fd = brep.face(f.clone());
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

    // ---- BRepBuilder extension method tests ----

    #[test]
    fn test_builder_set_edge_range() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = bld.add_edge(&mut brep, None, v0, v1, [0.0, 1.0]);
        assert_eq!(brep.edge(e.clone()).range, [0.0, 1.0]);
        bld.set_edge_range(&mut brep, e.clone(), 0.5, 2.5);
        assert_eq!(brep.edge(e.clone()).range, [0.5, 2.5]);
    }

    #[test]
    fn test_builder_set_edge_same_parameter() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = bld.add_edge(&mut brep, None, v0, v1, [0.0, 1.0]);
        assert!(brep.edge(e.clone()).same_parameter); // default is true in add_tedge
        bld.set_edge_same_parameter(&mut brep, e.clone(), false);
        assert!(!brep.edge(e.clone()).same_parameter);
        bld.set_edge_same_parameter(&mut brep, e.clone(), true);
        assert!(brep.edge(e.clone()).same_parameter);
    }

    #[test]
    fn test_builder_set_edge_same_range() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = bld.add_edge(&mut brep, None, v0, v1, [0.0, 1.0]);
        assert!(brep.edge(e.clone()).same_range); // default is true in add_tedge
        bld.set_edge_same_range(&mut brep, e.clone(), false);
        assert!(!brep.edge(e.clone()).same_range);
    }

    #[test]
    fn test_builder_set_face_natural_restriction() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        assert!(brep.face(f.clone()).natural_restriction);
        bld.set_face_natural_restriction(&mut brep, f.clone(), false);
        assert!(!brep.face(f.clone()).natural_restriction);
    }

    #[test]
    fn test_builder_update_face_tolerance() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        assert!((brep.face(f.clone()).tolerance - 0.0).abs() < 1e-15);
        bld.update_face_tolerance(&mut brep, f.clone(), 1.5);
        assert!((brep.face(f.clone()).tolerance - 1.5).abs() < 1e-15);
        bld.update_face_tolerance(&mut brep, f.clone(), 0.5); // smaller max keeps 1.5
        assert!((brep.face(f.clone()).tolerance - 1.5).abs() < 1e-15);
    }

    #[test]
    fn test_builder_update_vertex_point() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        assert_eq!(brep.vertex(v.clone()).point, DVec3::ZERO);
        bld.update_vertex_point(&mut brep, v.clone(), DVec3::new(5.0, 0.0, 0.0), 0.1);
        assert_eq!(brep.vertex(v.clone()).point, DVec3::new(5.0, 0.0, 0.0));
        assert!((brep.vertex(v.clone()).tolerance - 0.1).abs() < 1e-15);
    }

    #[test]
    fn test_builder_remove_from_wire() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = bld.add_edge(&mut brep, None, v0, v1, [0.0, 1.0]);
        let wire = bld.build_wire(&mut brep, vec![e.clone()]);
        assert_eq!(brep.wire(wire.clone()).edges.len(), 1);
        bld.remove_from_wire(&mut brep, wire.clone(), e);
        assert!(brep.wire(wire.clone()).edges.is_empty());
    }

    #[test]
    fn test_builder_remove_from_face() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let outer = brep.add_twire(vec![]);
        let inner = brep.add_twire(vec![]);
        let f = brep.add_tface(None, outer, vec![inner.clone()], None, None, vec![], true);
        assert_eq!(brep.face(f.clone()).inner_wires.len(), 1);
        bld.remove_from_face(&mut brep, f.clone(), inner);
        assert!(brep.face(f.clone()).inner_wires.is_empty());
    }

    #[test]
    fn test_builder_remove_from_shell() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let w = brep.add_twire(vec![]);
        let f = brep.add_tface(None, w, vec![], None, None, vec![], true);
        let shell = bld.make_shell(&mut brep);
        bld.add_to_shell(&mut brep, shell.clone(), f.clone());
        assert_eq!(brep.shell(shell.clone()).faces.len(), 1);
        bld.remove_from_shell(&mut brep, shell.clone(), f.clone());
        assert!(brep.shell(shell.clone()).faces.is_empty());
    }

    #[test]
    fn test_builder_transfert_edge_curve() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let curve = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let e_in = bld.add_edge(&mut brep, Some(curve.clone()), v0, v1, [0.0, 1.0]);

        // Use different vertices so add_tedge creates a distinct edge (dedup by vertex pair)
        let v2 = brep.add_tvertex(DVec3::Y);
        let v3 = brep.add_tvertex(DVec3::new(1.0, 1.0, 0.0));
        let e_out = bld.add_edge(&mut brep, None, v2, v3, [0.0, 1.0]);

        assert!(brep.edge(e_out.clone()).curve.is_none());
        bld.transfert_edge_curve(&mut brep, e_in, e_out.clone());
        assert!(
            brep.edge(e_out.clone()).curve.is_some(),
            "Curve should be transferred"
        );
    }

    #[test]
    fn test_builder_transfert_vertex_param() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = bld.add_edge(&mut brep, None, v0.clone(), v1.clone(), [0.0, 1.0]);
        bld.set_vertex_param(&mut brep, e.clone(), v0.clone(), 0.5);
        assert_eq!(brep.edge(e.clone()).vertex_params.get(&v0.index), Some(&0.5));

        // Transfer from v0 on e to v_new on a new edge
        let v_new = brep.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let e_new = bld.add_edge(&mut brep, None, v_new.clone(), v1.clone(), [0.0, 1.0]);
        assert!(brep.edge(e_new.clone()).vertex_params.get(&v_new.index).is_none());
        bld.transfert_vertex_param(&mut brep, e.clone(), v0.clone(), e_new.clone(), v_new.clone());
        assert_eq!(brep.edge(e_new.clone()).vertex_params.get(&v_new.index), Some(&0.5));
    }

    #[test]
    fn test_builder_set_edge_degenerated_with_clear() {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let curve = Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        });
        let e = bld.add_edge(&mut brep, Some(curve), v0, v1, [0.0, 1.0]);
        assert!(!brep.edge(e.clone()).degenerated);
        assert!(brep.edge(e.clone()).curve.is_some());
        bld.set_edge_degenerated_with_clear(&mut brep, e.clone(), true);
        assert!(brep.edge(e.clone()).degenerated);
        assert!(
            brep.edge(e.clone()).curve.is_none(),
            "3D curve should be cleared when setting degenerated"
        );
    }
}
