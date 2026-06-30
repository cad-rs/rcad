use std::sync::Arc;
use std::collections::HashMap;
use glam::DVec3;
use serde::{Deserialize, Serialize};
use crate::geom::{Curve2d, Curve3, Surface3};

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

/// OCCT TopAbs_ShapeEnum
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ShapeType {
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
}

/// TopoDS_Shape equivalent: index into the TShape pool + Orientation.
/// This is a value type (Copy, small) — the geometric data lives in the TShape Arc.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ShapeRef {
    /// Index into BRep.tshapes[]
    pub index: usize,
    pub orientation: Orientation,
}

impl ShapeRef {
    pub const fn new(index: usize) -> Self {
        Self { index, orientation: Orientation::Forward }
    }
    pub const fn with_orientation(index: usize, orientation: Orientation) -> Self {
        Self { index, orientation }
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TVertexData {
    pub point: DVec3,
    /// BRep_Tool::Tolerance(aV) equivalent — vertex tolerance.
    #[serde(default)]
    pub tolerance: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TEdgeData {
    pub curve: Option<usize>,
    pub first: ShapeRef,
    pub last: ShapeRef,
    pub range: [f64; 2],
    /// BRep_Tool::Degenerated(aE) equivalent.
    #[serde(default)]
    pub degenerated: bool,
    /// Per-face pcurves: face ShapeRef index → (curve2d, t_first, t_last).
    /// OCCT: BRep_Tool::CurveOnSurface(aE, aF) → Handle(Geom2d_Curve).
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub pcurves: HashMap<usize, (Curve2d, f64, f64)>,
    /// Per-vertex parameter on this edge: vertex ShapeRef index → param.
    /// OCCT: BRep_Tool::Parameter(aV, aE, aF).
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub vertex_params: HashMap<usize, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TWireData {
    pub edges: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TFaceData {
    pub surface: Option<usize>,
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
    pub faces: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TSolidData {
    pub shells: Vec<ShapeRef>,
    /// INTERNAL vertices embedded in the solid (OCCT: sub-shapes of TopoDS_Solid).
    pub internal_vertices: Vec<ShapeRef>,
    /// INTERNAL edges embedded in the solid (OCCT: sub-shapes of TopoDS_Solid).
    pub internal_edges: Vec<ShapeRef>,
}

/// BRep top-level shape container — all TShapes in a single pool with shared Arc ownership.
/// Analogous to OCCT's Doc/assembly structure where all TShapes live in a shared scope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub tshapes: Vec<Arc<TShape>>,
    pub curves: Vec<Curve3>,
    pub surfaces: Vec<Surface3>,
    pub curve2ds: Vec<crate::geom::Curve2d>,
}

impl Default for BRep {
    fn default() -> Self { Self::new() }
}

impl BRep {
    pub fn new() -> Self {
        Self { tshapes: Vec::new(), curves: Vec::new(), surfaces: Vec::new(), curve2ds: Vec::new() }
    }

    pub fn add_tvertex(&mut self, point: DVec3) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Vertex(TVertexData { point, tolerance: 0.0 })));
        ShapeRef::new(index)
    }

    pub fn add_tedge(&mut self, curve: Option<usize>, first: ShapeRef, last: ShapeRef, range: [f64; 2]) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Edge(TEdgeData { curve, first, last, range, degenerated: false, pcurves: HashMap::new(), vertex_params: HashMap::new() })));
        ShapeRef::new(index)
    }

    pub fn add_twire(&mut self, edges: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Wire(TWireData { edges })));
        ShapeRef::new(index)
    }

    pub fn add_tface(&mut self, surface: Option<usize>, outer_wire: ShapeRef, inner_wires: Vec<ShapeRef>, sample_point: Option<DVec3>, uv_domain: Option<[f64; 4]>, internal_vertices: Vec<ShapeRef>, natural_restriction: bool) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Face(TFaceData { surface, outer_wire, inner_wires, sample_point, uv_domain, internal_vertices, tolerance: 0.0, natural_restriction })));
        ShapeRef::new(index)
    }

    pub fn add_tshell(&mut self, faces: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Shell(TShellData { faces })));
        ShapeRef::new(index)
    }

    pub fn add_tsolid(&mut self, shells: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Solid(TSolidData { shells, internal_vertices: Vec::new(), internal_edges: Vec::new() })));
        ShapeRef::new(index)
    }

    pub fn add_tcompsolid(&mut self, solids: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::CompSolid(solids)));
        ShapeRef::new(index)
    }

    pub fn add_tcompound(&mut self, shapes: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Compound(shapes)));
        ShapeRef::new(index)
    }

    /// Remove all Solid TShapes from the BRep and return their indices.
    /// OCCT BuildRC removes unwanted solids from myShape after BuildResult(SOLID) added them.
    pub fn clear_solids(&mut self) -> usize {
        let before = self.tshapes.len();
        self.tshapes.retain(|ts| !matches!(&**ts, TShape::Solid(_)));
        before - self.tshapes.len()
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
    /// BRep_Tool::Surface(aF) — face surface.
    fn face_surface(&self, face: ShapeRef) -> Option<&Surface3>;
    /// UResolution: parameter tolerance in U direction (OCCT: BRepAdaptor_Surface::UResolution).
    fn u_resolution(&self, face: ShapeRef, tol3d: f64) -> f64;
    /// VResolution: parameter tolerance in V direction.
    fn v_resolution(&self, face: ShapeRef, tol3d: f64) -> f64;
    /// OCCT L204-207: vertex orientation (TopAbs_INTERNAL for split-edge interior vertices).
    /// Default: Forward (non-INTERNAL). Override when INTERNAL vertex data is available.
    fn vertex_orientation(&self, _v: ShapeRef) -> Orientation { Orientation::Forward }
}

impl BRepTool for BRep {
    fn vertex_position(&self, v: ShapeRef) -> DVec3 {
        self.vertex(v).point
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

    fn face_surface(&self, face: ShapeRef) -> Option<&Surface3> {
        let fi = face.index;
        self.face(face).surface.map(|si| &self.surfaces[si])
    }

    fn u_resolution(&self, face: ShapeRef, tol3d: f64) -> f64 {
        let surf_idx = match self.face(face).surface {
            Some(si) => si,
            None => return tol3d,
        };
        u_resolution_for_surface(&self.surfaces[surf_idx], tol3d)
    }

    fn v_resolution(&self, face: ShapeRef, tol3d: f64) -> f64 {
        let surf_idx = match self.face(face).surface {
            Some(si) => si,
            None => return tol3d,
        };
        v_resolution_for_surface(&self.surfaces[surf_idx], tol3d)
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
                return ShapeRef::new(i);
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
                    brep.add_tedge(None, va, ShapeRef::with_orientation(vb.index, Orientation::Reversed), [0.0, 1.0]).index
                });
                let orient = if a < b { Orientation::Forward } else { Orientation::Reversed };
                face_edges.push(ShapeRef::with_orientation(e_idx, orient));
            }
            edge_for_face.push(face_edges);
        }

        // Build 6 faces, collecting their refs for shell building
        let mut face_refs = Vec::new();
        for (i, face_edges) in edge_for_face.into_iter().enumerate() {
            let wire = brep.add_twire(face_edges);
            let face = brep.add_tface(None, wire, vec![], None, None, vec![], true);            let face_ref = ShapeRef::new(face.index);

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
            face_refs.push(ShapeRef::with_orientation(face.index, orient));
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
        curve: Option<usize>, v1: ShapeRef, v2: ShapeRef, range: [f64; 2]) -> ShapeRef {
        brep.add_tedge(curve, v1, v2, range)
    }

    /// Add a pcurve to an edge for a specific face.
    pub fn add_pcurve(&mut self, brep: &mut BRep,
        edge: ShapeRef, face: ShapeRef, pc: Curve2d, t1: f64, t2: f64) {
        brep.edge_mut(edge).pcurves.insert(face.index, (pc, t1, t2));
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

    /// Make a wire (empty container).
    pub fn make_wire(&mut self, brep: &mut BRep) -> ShapeRef {
        brep.add_twire(vec![])
    }

    /// Add an edge to a wire.
    pub fn add_to_wire(&mut self, brep: &mut BRep, wire: ShapeRef, edge: ShapeRef) {
        let wd = brep.wire_mut(wire);
        wd.edges.push(edge);
    }

    /// Build a wire from edges.
    pub fn build_wire(&mut self, brep: &mut BRep, edges: Vec<ShapeRef>) -> ShapeRef {
        brep.add_twire(edges)
    }

    /// Make a face from a surface and outer wire.
    pub fn make_face(&mut self, brep: &mut BRep,
        surface: Option<usize>, outer_wire: ShapeRef) -> ShapeRef {
        brep.add_tface(surface, outer_wire, vec![], None, None, vec![], true)
    }

    /// Add an inner wire to a face.
    pub fn add_to_face(&mut self, brep: &mut BRep, face: ShapeRef, inner_wire: ShapeRef) {
        let fd = brep.face_mut(face);
        fd.inner_wires.push(inner_wire);
    }

    /// Add an internal vertex to a face.
    pub fn add_internal_vertex(&mut self, brep: &mut BRep, face: ShapeRef, v: ShapeRef) {
        let fd = brep.face_mut(face);
        fd.internal_vertices.push(v);
    }

    /// Add an edge with section-curve semantics (MakeSectEdge equivalent).
    /// Creates an edge with pcurves for both faces.
    pub fn add_section_edge(&mut self, brep: &mut BRep,
        curve: Option<usize>, v1: ShapeRef, v2: ShapeRef, range: [f64; 2],
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
        let r = ShapeRef::new(5);
        assert_eq!(r.index, 5);
        assert_eq!(r.orientation, Orientation::Forward);

        let r2 = ShapeRef::with_orientation(3, Orientation::Reversed);
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
        let e = brep.add_tedge(None, v0, ShapeRef::with_orientation(v1.index, Orientation::Reversed), [0.0, 1.0]);
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
    fn test_serialize_roundtrip() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        let e = brep.add_tedge(None, v, v, [0.0, 1.0]);

        let json = serde_json::to_string(&brep).unwrap();
        let restored: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.tshapes.len(), 2);
        assert_eq!(restored.vertex(ShapeRef::new(0)).point, DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(restored.edge(ShapeRef::new(1)).range, [0.0, 1.0]);
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
        let e0_rev = ShapeRef::with_orientation(e0.index, Orientation::Reversed);
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
